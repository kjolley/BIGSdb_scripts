#!/usr/bin/env perl
#Update A. baumannii database with betalactamase information.
#Written by Julia Moreno & Keith Jolley, 2022.
#Version:20221108
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
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Constants qw(LOG_TO_SCREEN);
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
use Data::Dumper;
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
my %opts;
GetOptions(
	'help'         => \$opts{'h'},
	'i|isolates=s' => \$opts{'i'},
	'q|quiet'      => \$opts{'quiet'},
	'x|min=i'      => \$opts{'x'},
	'y|max=i'      => \$opts{'y'},
) or die("Error in command line arguments\n");

if ( $opts{'h'} ) {
	show_help();
	exit;
}

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
	my $allele_sequences = get_allele_data();
	link_b_lactamase_to_alleles( $allele_sequences, $peptide_lookup );
	my $allele_sequence_hash = {};
	foreach my $as (@$allele_sequences) {
		$allele_sequence_hash->{ $as->{'locus'} }->{ $as->{'allele_id'} } = $as;
	}
	my $isolates = get_isolates();
	foreach my $isolate_id (@$isolates) {
		process_isolate( $isolate_id, $allele_sequence_hash );
	}
	return;
}

sub process_isolate {
	my ( $isolate_id, $allele_sequence_hash ) = @_;
	my @loci = ( OXA_LOCI, MBL_LOCI );
	my @placeholders = ('?') x OXA_LOCI;
	local $" = ',';
	my ( %oxa_family, %oxa_class, %oxa_bl, %mbl_family, %mbl_class, %mbl_bl, %bl );
	my $allele_designations = $isolate_obj->{'datastore'}->run_query(
		"SELECT locus,allele_id FROM allele_designations WHERE isolate_id=? AND locus IN (@placeholders)",
		[ $isolate_id, OXA_LOCI ],
		{ fetch => 'all_arrayref', slice => {} }
	);
	foreach my $ad (@$allele_designations) {
		if ( $allele_sequence_hash->{ $ad->{'locus'} }->{ $ad->{'allele_id'} } ) {
			my $as = $allele_sequence_hash->{ $ad->{'locus'} }->{ $ad->{'allele_id'} };
			$oxa_family{ $as->{'family'} } = 1;
			$oxa_class{ $as->{'class'} }   = 1;
			if ( defined $as->{'b-lactamase'} ) {
				$oxa_bl{"$as->{'family'} [$as->{'b-lactamase'}]"} = 1;
				$bl{ $as->{'b-lactamase'} } = 1;
			}
		}
	}
	@placeholders        = ('?') x MBL_LOCI;
	$allele_designations = $isolate_obj->{'datastore'}->run_query(
		"SELECT locus,allele_id FROM allele_designations WHERE isolate_id=? AND locus IN (@placeholders)",
		[ $isolate_id, MBL_LOCI ],
		{ fetch => 'all_arrayref', slice => {} }
	);
	foreach my $ad (@$allele_designations) {
		if ( $allele_sequence_hash->{ $ad->{'locus'} }->{ $ad->{'allele_id'} } ) {
			my $as = $allele_sequence_hash->{ $ad->{'locus'} }->{ $ad->{'allele_id'} };
			$mbl_family{ $as->{'family'} } = 1;
			$mbl_class{ $as->{'class'} }   = 1;
			if ( defined $as->{'b-lactamase'} ) {
				$mbl_bl{"$as->{'family'} [$as->{'b-lactamase'}]"} = 1;
				$bl{"$as->{'b-lactamase'} "}                      = 1;
			}
		}
	}
	my @oxa_family = sort keys %oxa_family;
	my @oxa_class  = sort keys %oxa_class;
	my @oxa_bl     = sort keys %oxa_bl;
	my @mbl_family = sort keys %mbl_family;
	my @mbl_class  = sort keys %mbl_class;
	my @mbl_bl     = sort keys %mbl_bl;
	my @bl         = sort keys %bl;
	local $" = ';';
	return if !@bl;
	say "id:$isolate_id; OXA_family:@oxa_family; OXA_class:@oxa_class; OXA_bl:@oxa_bl; "
	  . "MBL_family:@mbl_family; MBL_class:@mbl_class; MBL_bl:@mbl_bl; bl:@bl"
	  if !$opts{'quiet'};
	eval {
		$isolate_obj->{'db'}->do(
			'UPDATE isolates SET (oxa_family,oxa_class,oxa_betalactamase,mbl_family,mbl_class,'
			  . 'mbl_betalactamase,class_b_d_betalactamase)=(?,?,?,?,?,?,?) WHERE id=?',
			undef,
			@oxa_family ? BIGSdb::Utils::get_pg_array( \@oxa_family ) : undef,
			@oxa_class  ? BIGSdb::Utils::get_pg_array( \@oxa_class )  : undef,
			@oxa_bl     ? BIGSdb::Utils::get_pg_array( \@oxa_bl )     : undef,
			@mbl_family ? BIGSdb::Utils::get_pg_array( \@mbl_family ) : undef,
			@mbl_class  ? BIGSdb::Utils::get_pg_array( \@mbl_class )  : undef,
			@mbl_bl     ? BIGSdb::Utils::get_pg_array( \@mbl_bl )     : undef,
			@bl         ? BIGSdb::Utils::get_pg_array( \@bl )         : undef,
			$isolate_id
		);
	};

	if ($@) {
		$logger->error($@);
		$isolate_obj->{'db'}->rollback;
		exit;
	}
	$isolate_obj->{'db'}->commit;
}

sub get_isolates {
	my $list =
	  $isolate_obj->{'datastore'}
	  ->run_query( "SELECT id FROM $isolate_obj->{'system'}->{'view'} WHERE class_b_d_betalactamase IS NULL",
		undef, { fetch => 'col_arrayref' } );
	my $isolates = $isolate_obj->get_isolates_with_linked_seqs( { size => 1_000_000 } );
	my %no_betalactamase_set = map { $_ => 1 } @$list;
	my $filtered_list = [];
	foreach my $id (@$isolates) {
		push @$filtered_list, $id if $no_betalactamase_set{$id};
	}
	$filtered_list = $isolate_obj->filter_and_sort_isolates($filtered_list);
	return $filtered_list;
}

sub link_b_lactamase_to_alleles {
	my ( $allele_sequences, $peptide_lookup ) = @_;
	my %family = FAMILY;
	my %class  = CLASS;
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
	return;
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
	my $options = $instance eq ISOLATE_DATABASE ? \%opts : {};
	return BIGSdb::Offline::Script->new(
		{
			config_dir       => CONFIG_DIR,
			lib_dir          => LIB_DIR,
			dbase_config_dir => DBASE_CONFIG_DIR,
			instance         => $instance,
			options          => $options,
			logger           => $logger
		}
	);
}

sub show_help {
	my $termios = POSIX::Termios->new;
	$termios->getattr;
	my $ospeed = $termios->getospeed;
	my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
	my ( $norm, $bold, $under ) = map { $t->Tputs( $_, 1 ) } qw(me md us);
	say << "HELP";
${bold}NAME$norm
    ${bold}AB_blactamase.pl$norm - Populate beta_lactamase fields in A. baumannii database

${bold}SYNOPSIS$norm
    ${bold}AB_blactamase.pl$norm [${under}options$norm]

${bold}OPTIONS$norm

${bold}--help$norm
    This help page.
    
${bold}--isolates$norm ${under}LIST$norm  
    Comma-separated list of isolate ids to scan.
      
${bold}--quiet$norm
    Suppress output, only showing errors.
       
${bold}-x, --min$norm ${under}ID$norm
    Minimum isolate id.

${bold}-y, --max$norm ${under}ID$norm
    Maximum isolate id.   
HELP
	return;
}
