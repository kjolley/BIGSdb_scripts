#!/usr/bin/env perl
#Update A. baumannii database with
#Written by Julia Moreno & Keith Jolley, 2022.
#Version:20221107
use strict;
use warnings;
use 5.010;
use Bio::Seq;
use Log::Log4perl qw(get_logger);
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	SEQDEF_DATABASE  => 'pubmlst_abaumannii_seqdef',
	ISOLATE_DATABASE => 'pubmlst_abaumannii_isolates',
};
use constant OXA_LOCI => qw(ACIN07895 ACIN20010 ACIN20011 ACIN20012 ACIN20013 ACIN20014);
use constant MBL_LOCI => qw(ACIN20004 ACIN20005 ACIN20006 ACIN20007 ACIN20008 ACIN20009);
use constant CLASS    => (
	ACIN07895 => 'Oxacilliniases (Class D)',
	ACIN20010 => 'Oxacilliniases (Class D)',
	ACIN20011 => 'Oxacilliniases (Class D)',
	ACIN20012 => 'Oxacilliniases (Class D)',
	ACIN20013 => 'Oxacilliniases (Class D)',
	ACIN20014 => 'Oxacilliniases (Class D)',
	ACIN20004 => 'Metallo-beta-lactamases (Class B)',
	ACIN20005 => 'Metallo-beta-lactamases (Class B)',
	ACIN20006 => 'Metallo-beta-lactamases (Class B)',
	ACIN20007 => 'Metallo-beta-lactamases (Class B)',
	ACIN20008 => 'Metallo-beta-lactamases (Class B)',
	ACIN20009 => 'Metallo-beta-lactamases (Class B)'
);
use constant FAMILY => (
	ACIN07895 => 'blaOXA-51-like',
	ACIN20010 => 'blaOXA-23-like',
	ACIN20011 => 'blaOXA-24/40-like',
	ACIN20012 => 'blaOXA-58-like',
	ACIN20013 => 'blaOXA-143-like',
	ACIN20014 => 'blaOXA-134-like',
	ACIN20004 => 'blaIMP',
	ACIN20005 => 'blaNDM',
	ACIN20006 => 'blaVIM',
	ACIN20007 => 'blaSIM',
	ACIN20008 => 'blaGIM',
	ACIN20009 => 'blaSPM'
);
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Constants qw(LOG_TO_SCREEN);
use Data::Dumper;

#Direct all library logging calls to screen
my $log_conf = LOG_TO_SCREEN;
Log::Log4perl->init( \$log_conf );
my $logger      = Log::Log4perl::get_logger('BIGSdb.Script');
my $seqdef_obj  = get_script_obj(SEQDEF_DATABASE);
my $isolate_obj = get_script_obj(ISOLATE_DATABASE);
main();
undef $seqdef_obj;
undef $isolate_obj;

sub main {
	my $peptide_sequences = get_peptide_data();
	my $peptide_lookup    = {};
	foreach my $peptide (@$peptide_sequences) {
		$peptide_lookup->{ $peptide->{'sequence'} } = $peptide;
		delete $peptide_lookup->{ $peptide->{'sequence'} }->{'sequence'};
	}


	my %family = FAMILY;
	my %class  = CLASS;

	#	say Dumper $peptide_sequences;
	my $allele_sequences = get_allele_data();
	foreach my $allele (@$allele_sequences) {
		if ( $family{ $allele->{'locus'} } ) {
			$allele->{'family'} = $family{ $allele->{'locus'} };
		} else {
			die "$allele->{'locus'} does not have a family defined!\n";
		}
		if ( $class{ $allele->{'locus'} } ) {
			$allele->{'class'} = $class{ $allele->{'locus'} };
		} else {
			die "$allele->{'locus'} does not have a class defined!\n";
		}
		if ( $allele->{'sequence'} =~ /([^ACGT])/ ) {
			die "$allele->{'locus'} $allele->{'allele_id'} has an invalid char $1.\n";
		}
		my $bioseq_obj = Bio::Seq->new(
			-display_id => 'allele',
			-seq        => $allele->{'sequence'}
		);
		my $peptide_seq = $bioseq_obj->translate( -codontable_id => 11 )->seq();
		$peptide_seq =~ s/\*$//;
		if ( defined $peptide_lookup->{$peptide_seq}->{'b-lactamase'} ) {
			$allele->{'b-lactamase'} = $peptide_lookup->{$peptide_seq}->{'b-lactamase'};
		}
	}
	say Dumper $allele_sequences;
}

sub get_peptide_data {
	my $peptides = $seqdef_obj->{'datastore'}->run_query(
		'SELECT locus,allele_id,sequence FROM sequences WHERE locus IN (?,?)',
		[ 'OXA_peptide', 'MBL_peptide' ],
		{ fetch => 'all_arrayref', slice => {} }
	);
	my $blactamases_list = $seqdef_obj->{'datastore'}->run_query(
		'SELECT locus,allele_id,value FROM sequence_extended_attributes WHERE locus IN (?,?) AND field=?',
		[ 'OXA_peptide', 'MBL_peptide', 'beta-lactamase' ],
		{ fetch => 'all_arrayref', slice => {} }
	);
	my $blactamases = {};
	foreach my $blactamase (@$blactamases_list) {
		$blactamases->{ $blactamase->{'locus'} }->{ $blactamase->{'allele_id'} } = $blactamase->{'value'};
	}
	foreach my $peptide (@$peptides) {
		if ( defined $blactamases->{ $peptide->{'locus'} }->{ $peptide->{'allele_id'} } ) {
			$peptide->{'b-lactamase'} = $blactamases->{ $peptide->{'locus'} }->{ $peptide->{'allele_id'} };
		}
	}
	return $peptides;
}

sub get_allele_data {
	my @loci = ( OXA_LOCI, MBL_LOCI );
	my @placeholders = ('?') x @loci;
	local $" = ',';
	my $alleles =
	  $seqdef_obj->{'datastore'}
	  ->run_query( "SELECT locus,allele_id,sequence FROM sequences WHERE locus IN (@placeholders) AND allele_id != 'N'",
		[@loci], { fetch => 'all_arrayref', slice => {} } );
	return $alleles;
}

sub get_script_obj {
	my ($instance) = @_;
	return BIGSdb::Offline::Script->new(
		{
			config_dir       => CONFIG_DIR,
			lib_dir          => LIB_DIR,
			dbase_config_dir => DBASE_CONFIG_DIR,
			instance         => $instance,
			logger           => $logger
		}
	);
}
