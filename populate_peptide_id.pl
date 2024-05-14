#!/usr/bin/env perl
#Populate peptide_id field in DNA locus corresponding to a peptide locus.
#Written by Keith Jolley 2024
#Version 20240513
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	HOST             => undef,                  #Use values in config.xml
	PORT             => undef,                  #But you can override here.
	USER             => undef,
	PASSWORD         => undef
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use Getopt::Long qw(:config no_ignore_case);
use Bio::Seq;
use BIGSdb::Offline::Script;
use BIGSdb::Constants qw(LOG_TO_SCREEN);

#Direct all library logging calls to screen
my $log_conf = LOG_TO_SCREEN;
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my %opts;
GetOptions(
	'c|check_existing'     => \$opts{'check'},
	'd|database=s'         => \$opts{'database'},
	'n|nucleotide_locus=s' => \$opts{'nucleotide_locus'},
	'peptide_id_field=s'   => \$opts{'peptide_id_field'},
	'p|protein_locus=s'    => \$opts{'protein_locus'},
	'q|quiet'              => \$opts{'quiet'},
) or die("Error in command line arguments\n");
my $script = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		instance         => $opts{'database'},
		logger           => $logger
	}
);
die "Database not defined.\n" if !$opts{'database'};
die "Script initialization failed - check logs (authentication problems or server too busy?).\n"
  if !defined $script->{'db'};
die "This script can only be run against a seqdef database.\n"
  if ( $script->{'system'}->{'dbtype'} // '' ) ne 'sequences';
$opts{'peptide_id_field'} //= 'peptide_id';
check_loci();
main();
undef $script;

sub check_loci {
	die "Nucleotide locus not set.\n" if !defined $opts{'nucleotide_locus'};
	my $nuc_locus = $script->{'datastore'}->get_locus_info( $opts{'nucleotide_locus'} );
	die "Nucleotide locus $opts{'nucleotide_locus'} does not exist.\n"            if !defined $nuc_locus;
	die "Nucleotide locus $opts{'nucleotide_locus'} is not set as a DNA locus.\n" if $nuc_locus->{'data_type'} ne 'DNA';
	die "Protein locus not set.\n"                                                if !defined $opts{'protein_locus'};
	my $protein_locus = $script->{'datastore'}->get_locus_info( $opts{'protein_locus'} );
	die "Protein locus $opts{'protein_locus'} does not exist.\n" if !defined $protein_locus;
	die "Protein locus $opts{'protein_locus'} is not set as a protein locus.\n"
	  if $protein_locus->{'data_type'} ne 'peptide';
	return;
}

sub main {
	if ( $opts{'check'} ) {
		check();
	}
}

sub get_nucleotide_sequences {
	my $nuc_locus = $script->{'datastore'}->get_locus_info( $opts{'nucleotide_locus'} );
	my $order_by  = $nuc_locus->{'allele_id_format'} eq 'text' ? 's.allele_id' : 'CAST(s.allele_id AS int)';
	my $alleles   = $script->{'datastore'}->run_query(
		"SELECT s.allele_id,s.sequence,se.value AS peptide_id FROM sequences s LEFT JOIN "
		  . "sequence_extended_attributes se ON "
		  . "(s.locus,s.allele_id,se.field)=(se.locus,se.allele_id,?) WHERE s.locus=? AND "
		  . "s.allele_id NOT IN ('N','0','P') ORDER BY $order_by",
		[ $opts{'peptide_id_field'}, $opts{'nucleotide_locus'} ],
		{ fetch => 'all_arrayref', slice => {} }
	);
	return $alleles;
}

sub get_protein_sequences {
	my $protein_locus = $script->{'datastore'}->get_locus_info( $opts{'protein_locus'} );
	my $peptides =
	  $script->{'datastore'}
	  ->run_query( "SELECT allele_id,sequence FROM sequences WHERE locus=? AND allele_id NOT IN ('N','0','P')",
		$opts{'protein_locus'}, { fetch => 'all_arrayref', slice => {} } );
	my %seqs = map { $_->{'allele_id'} => $_->{'sequence'} } @$peptides;
	return \%seqs;
}

sub check {
	my $alleles     = get_nucleotide_sequences();
	my $peptides    = get_protein_sequences();
	my $codon_table = $script->{'system'}->{'codon_table'} // 11;
	my $nuc_locus   = $script->{'datastore'}->get_locus_info( $opts{'nucleotide_locus'} );
	my $orf         = $nuc_locus->{'orf'} // 1;
	foreach my $allele (@$alleles) {
		next if !defined $allele->{ $opts{'peptide_id_field'} };
		my $translate = translate_sequence( $allele->{'sequence'}, $codon_table, $orf );
		my $id        = "$opts{'nucleotide_locus'}-$allele->{'allele_id'} [peptide_id: $allele->{'peptide_id'}]";
		if ( $allele->{'peptide_id'} ) {
			my $peptide_seq = $peptides->{ $allele->{'peptide_id'} };
			if ( !defined $peptide_seq ) {
				say "$id: Peptide $allele->{'peptide_id'} does not exist";
				next;
			}
			if ( $translate->{'error'} ) {
				say "$id: $translate->{'error'}";
				next;
			}
			if ( $translate->{'peptide'} ne $peptide_seq ) {
				say "$id: Translated seq does not match peptide record";
				next;
			}
		} else {
			if ( !$translate->{'error'} ) {
				say "$id: No peptide_id set";
				next;
			}
		}
	}
}

sub translate_sequence {
	my ( $sequence, $codon_table, $orf ) = @_;
	my $reverse;
	if ( $orf > 3 ) {
		$reverse = 1;
		$orf     = $orf - 3;
	}
	my $seq = $reverse ? BIGSdb::Utils::reverse_complement($sequence) : $sequence;
	if ( $orf > 1 && $orf <= 3 ) {
		$seq = substr( $seq, $orf - 1 );
	}
	my $seq_obj = Bio::Seq->new( -seq => $seq, -alphabet => 'dna' );
	my $peptide = $seq_obj->translate( -codontable_id => $codon_table )->seq;
	$peptide =~ s/\*$//x;    #Remove terminal stop codon.
	my $error;
	$error = 'Internal stop codon' if $peptide =~ /\*/gx;
	return {
		peptide => $peptide,
		error   => $error
	};
}
