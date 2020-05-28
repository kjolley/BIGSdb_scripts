#!/usr/bin/env perl
#Predict cross-reactivity to Bexsero and Trumenba vaccine components
#Written by Keith Jolley
#Copyright (c) 2018-2020, University of Oxford
#E-mail: keith.jolley@zoo.ox.ac.uk
#Version: 20200522
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
use List::MoreUtils qw(uniq);
use BIGSdb::Offline::Script;

#Direct all library logging calls to screen
my $log_conf =
    qq(log4perl.category.BIGSdb.Script        = INFO, Screen\n)
  . qq(log4perl.category.BIGSdb.Dataconnector = WARN, Screen\n)
  . qq(log4perl.category.BIGSdb.Datastore     = WARN, Screen\n)
  . qq(log4perl.category.BIGSdb.Scheme        = WARN, Screen\n)
  . qq(log4perl.appender.Screen               = Log::Log4perl::Appender::Screen\n)
  . qq(log4perl.appender.Screen.stderr        = 1\n)
  . qq(log4perl.appender.Screen.layout        = Log::Log4perl::Layout::SimpleLayout\n);
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my %opts;
GetOptions(
	'd|database=s'        => \$opts{'database'},
	'l|loci_designated=i' => \$opts{'loci_designated'},
	'q|quiet'             => \$opts{'quiet'},
	'r|refresh'           => \$opts{'refresh'}
) or die("Error in command line arguments\n");
if ( !$opts{'database'} ) {
	say "\nUsage: vaccine_protection_prediction.pl --database <NAME> \n";
	exit;
}
$opts{'loci_designated'} //= 1000;
my $script = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		host             => HOST,
		port             => PORT,
		user             => USER,
		password         => PASSWORD,
		options          => \%opts,
		instance         => $opts{'database'},
		logger           => $logger
	}
);
die "This script can only be run against an isolate database.\n"
  if ( $script->{'system'}->{'dbtype'} // '' ) ne 'isolates';
main();
undef $script;

sub main {
	my $qry =
	    qq(SELECT id FROM $script->{'system'}->{'view'} WHERE species='Neisseria meningitidis' AND id )
	  . q(IN (SELECT isolate_id FROM seqbin_stats WHERE total_length>1500000));
	$qry .= q(AND id NOT IN (SELECT isolate_id FROM eav_text WHERE field='Bexsero_reactivity')) if !$opts{'refresh'};
	$qry .= q(ORDER BY id);
	my $bexsero_id_list = $script->{'datastore'}->run_query( $qry, undef, { fetch => 'col_arrayref' } );
	foreach my $isolate_id (@$bexsero_id_list) {
		next if !contains_bexsero_antigens($isolate_id);
		my $designated_count =
		  $script->{'datastore'}->run_query( 'SELECT COUNT(*) FROM allele_designations WHERE isolate_id=?',
			$isolate_id, { cache => 'designation_count' } );
		next if $designated_count < $opts{'loci_designated'};
		my $has_exact_bexsero_match = has_exact_bexsero_match($isolate_id);
		if ( $has_exact_bexsero_match->{'result'} ) {
			say qq($isolate_id Bexsero exact match) if !$opts{'quiet'};
			set_value( $isolate_id, 'Bexsero_reactivity', 'exact match' );
			local $" = q(; );
			set_value( $isolate_id, 'Bexsero_notes', qq(@{ $has_exact_bexsero_match->{'notes'} } ) );
			next;
		}
		my $cross_reacts_with_bexsero = cross_reacts_with_bexsero($isolate_id);
		if ( $cross_reacts_with_bexsero->{'result'} ) {
			say qq($isolate_id Bexsero cross-reactive) if !$opts{'quiet'};
			set_value( $isolate_id, 'Bexsero_reactivity', 'cross-reactive' );
			local $" = q(; );
			set_value( $isolate_id, 'Bexsero_notes', qq(@{ $cross_reacts_with_bexsero->{'notes'} } ) );
			next;
		}
		my $no_reactivity_with_bexsero = no_reactivity_with_bexsero($isolate_id);
		if ( $no_reactivity_with_bexsero->{'result'} ) {
			say qq($isolate_id Bexsero none) if !$opts{'quiet'};
			set_value( $isolate_id, 'Bexsero_reactivity', 'none' );
			local $" = q(; );
			set_value( $isolate_id, 'Bexsero_notes', qq(@{ $no_reactivity_with_bexsero->{'notes'} } ) );
			next;
		}
		say qq($isolate_id Bexsero insufficient data) if !$opts{'quiet'};
		set_value( $isolate_id, 'Bexsero_reactivity', 'insufficient data' );
	}
	$qry =
	    qq(SELECT id FROM $script->{'system'}->{'view'} WHERE species='Neisseria meningitidis' AND id )
	  . q(IN (SELECT isolate_id FROM seqbin_stats WHERE total_length>1500000));
	$qry .= q(AND id NOT IN (SELECT isolate_id FROM eav_text WHERE field='Trumenba_reactivity')) if !$opts{'refresh'};
	$qry .= q(ORDER BY id);
	my $trumenba_id_list = $script->{'datastore'}->run_query( $qry, undef, { fetch => 'col_arrayref' } );
	foreach my $isolate_id (@$trumenba_id_list) {
		next if !contains_trumenba_antigens($isolate_id);
		my $designated_count =
		  $script->{'datastore'}->run_query( 'SELECT COUNT(*) FROM allele_designations WHERE isolate_id=?',
			$isolate_id, { cache => 'designation_count' } );
		next if $designated_count < $opts{'loci_designated'};
		my $has_exact_trumenba_match = has_exact_trumenba_match($isolate_id);
		if ( $has_exact_trumenba_match->{'result'} ) {
			say qq($isolate_id Trumenba exact match) if !$opts{'quiet'};
			set_value( $isolate_id, 'Trumenba_reactivity', 'exact match' );
			set_value( $isolate_id, 'Trumenba_notes',      @{ $has_exact_trumenba_match->{'notes'} } );
			next;
		}
		my $cross_reacts_with_trumenba = cross_reacts_with_trumenba($isolate_id);
		if ( $cross_reacts_with_trumenba->{'result'} ) {
			say qq($isolate_id Trumenba cross-reactive) if !$opts{'quiet'};
			set_value( $isolate_id, 'Trumenba_reactivity', 'cross-reactive' );
			local $" = q(; );
			set_value( $isolate_id, 'Trumenba_notes', qq(@{ $cross_reacts_with_trumenba->{'notes'} } ) );
			next;
		}
		if (
			$script->{'datastore'}->run_query(
				'SELECT EXISTS(SELECT * FROM allele_designations WHERE (isolate_id,locus,allele_id)=(?,?,?))',
				, [ $isolate_id, 'fHbp_peptide', '0' ]
			)
		  )
		{
			say qq($isolate_id Trumenba none) if !$opts{'quiet'};
			set_value( $isolate_id, 'Trumenba_reactivity', 'none' );
			set_value( $isolate_id, 'Trumenba_notes',      'fHbp_peptide is missing' );
			next;
		}
		say qq($isolate_id Trumenba insufficient data) if !$opts{'quiet'};
		set_value( $isolate_id, 'Trumenba_reactivity', 'insufficient data' );
	}
	return;
}

sub contains_bexsero_antigens {
	my ($isolate_id) = @_;
	return $script->{'datastore'}->run_query(
		q[SELECT EXISTS(SELECT * FROM allele_designations WHERE locus ]
		  . q[IN ('fHbp_peptide','NHBA_peptide','NadA_peptide','PorA_VR2') AND isolate_id=?)],
		$isolate_id,
		{ cache => 'contains_bexsero_antigens' }
	);
}

sub has_exact_bexsero_match {
	my ($isolate_id) = @_;
	my $components = [
		{ locus => 'fHbp_peptide', variant => '1', references => { 27521232 => 'peptide sequence match' } },
		{ locus => 'NHBA_peptide', variant => '2', references => { 27521232 => 'peptide sequence match' } },
		{ locus => 'NadA_peptide', variant => '8', references => { 27521232 => 'peptide sequence match' } },
		{
			locus      => 'PorA_VR2',
			variant    => '4',
			references => { 27521232 => 'peptide sequence match' }
		}
	];
	return format_match_results( $isolate_id, $components, 'is exact match to vaccine variant' );
}

sub format_match_results {
	my ( $isolate_id, $components, $reason ) = @_;
	my $result;
	my $notes = [];
	my %assays = map { $_ => 1 } qw(MATS SBA MEASURE);
	foreach my $component (@$components) {
		if ( allele_designated( $isolate_id, $component->{'locus'}, $component->{'variant'} ) ) {
			$result = 1;
			my @note_evidence;
			my @evidence_list = uniq( sort { $a cmp $b } values %{ $component->{'references'} } );
			my $evidence_count = 0;
			foreach my $evidence (@evidence_list) {
				$evidence_count++;
				my @refs;
				foreach my $ref ( sort { $a <=> $b } keys %{ $component->{'references'} } ) {
					push @refs, $ref if $component->{'references'}->{$ref} eq $evidence;
				}
				if (@refs) {
					local $" = q(, PMID:);
					if ( $assays{$evidence} ) {
						if ( $evidence_count == 1 ) {
							$evidence = qq(data derived from $evidence assays);
						} elsif ( $evidence_count == @evidence_list ) {
							$evidence = qq(and $evidence assays);
						} else {
							$evidence = qq(, $evidence assays);
						}
					}
					push @note_evidence, qq($evidence (PMID:@refs));
				}
			}
			if (@note_evidence) {
				local $" = q(, );
				push @$notes, qq($component->{'locus'}: $component->{'variant'} $reason - @note_evidence);
			}
		}
	}
	return { result => $result, notes => $notes };
}

sub cross_reacts_with_bexsero {
	my ($isolate_id) = @_;
	my $components = [
		{
			locus      => 'fHbp_peptide',
			variant    => '4',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30135218 => 'MATS',
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '10',
			references => {
				25630407 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS',
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '12',
			references => {
				29152576 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '14',
			references => {
				23414709 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '15',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '37',
			references => {
				23414709 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '110',
			references => {
				26950303 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '144',
			references => {
				26686998 => 'MATS',
				26950303 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '215',
			references => {
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '232',
			references => {
				23414709 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '1',
			references => {
				27355628 => 'MATS',
				29152576 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '5',
			references => {
				29152576 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '10',
			references => {
				23414709 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '113',
			references => {
				23414709 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '243',
			references => {
				26686998 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '607',
			references => {
				30135218 => 'MATS'
			}
		},
		{
			locus      => 'NadA_peptide',
			variant    => '3',
			references => {
				23588089 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				29950334 => 'MATS'
			}
		},
		{
			locus      => 'NadA_peptide',
			variant    => '6',
			references => {
				29950334 => 'MATS'
			}
		}
	];
	return format_match_results( $isolate_id, $components, 'is cross-reactive to vaccine variant' );
}

sub no_reactivity_with_bexsero {
	my ($isolate_id) = @_;
	my $notes        = [];
	my $components   = [
		{
			locus      => 'fHbp_peptide',
			variant    => '16',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '19',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				25630407 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '21',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '22',
			references => {
				23414709 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '24',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '25',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '29',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS',
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '30',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '31',
			references => {
				26950303 => 'MATS',
				27083425 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '45',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '47',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '59',
			references => {
				23414709 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '76',
			references => {
				29152576 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '109',
			references => {
				23414709 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '119',
			references => {
				23414709 => 'MATS',
				28366725 => 'MATS'
			}
		},
	];

	#FHbp
	my $matched = is_not_cross_reactive( $isolate_id, $components, $notes );
	if ( allele_designated( $isolate_id, 'fHbp_peptide', '0' ) ) {
		$matched = 1;
		push @$notes, qq(fHbp_peptide is missing);
	}
	return { result => 0 } if !$matched;
	$components = [
		{
			locus      => 'NHBA_peptide',
			variant    => '6',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '9',
			references => {
				26950303 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '17',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				26950303 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				30135218 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '18',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '25',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '30',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '31',
			references => {
				23414709 => 'MATS',
				27083425 => 'MATS',
				28366725 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '43',
			references => {
				23414709 => 'MATS',
				27355628 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '47',
			references => {
				23414709 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '63',
			references => {
				26686998 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '112',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '120',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '160',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '187',
			references => {
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '197',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS'
			}
		},
	];

	#NHBA
	$matched = is_not_cross_reactive( $isolate_id, $components, $notes );
	if ( allele_designated( $isolate_id, 'NHBA_peptide', '0' ) ) {
		$matched = 1;
		push @$notes, qq(NHBA_peptide is missing);
	}
	return { result => 0 } if !$matched;
	$components = [
		{
			locus      => 'NadA_peptide',
			variant    => '1',
			references => {
				26950303 => 'MATS',
				28366725 => 'MATS',
				29950334 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NadA_peptide',
			variant    => '21',
			references => {
				26950303 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NadA_peptide',
			variant    => '100',
			references => {
				29950334 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS'
			}
		},
	];

	#NadA
	$matched = is_not_cross_reactive( $isolate_id, $components, $notes );
	if ( allele_designated( $isolate_id, 'NadA_peptide', '0' ) ) {
		$matched = 1;
		push @$notes, qq(NadA_peptide is missing);
	}
	return { result => 0 } if !$matched;

	#PorA
	my $PorA_VR2_missing =
	  $script->{'datastore'}
	  ->run_query( 'SELECT EXISTS(SELECT * FROM allele_designations WHERE (isolate_id,locus,allele_id)=(?,?,?))',
		, [ $isolate_id, 'PorA_VR2', '0' ] );
	if ($PorA_VR2_missing) {
		push @$notes, qq(PorA_VR2 is missing);
		return { result => 1, notes => $notes };
	}
	my $contains_PorA_VR2_4 =
	  $script->{'datastore'}
	  ->run_query( 'SELECT EXISTS(SELECT * FROM allele_designations WHERE (isolate_id,locus,allele_id)=(?,?,?))',
		[ $isolate_id, 'PorA_VR2', '4' ] );
	if ( !$contains_PorA_VR2_4 ) {
		push @$notes, qq(PorA_VR2 is not variant 4);
	}
	return { result => 1, notes => $notes };
}

sub is_not_cross_reactive {
	my ( $isolate_id, $components, $notes ) = @_;
	my $matched = 0;
	foreach my $component (@$components) {
		if ( allele_designated( $isolate_id, $component->{'locus'}, $component->{'variant'} ) ) {
			$matched = 1;
			my @note_evidence;
			my @evidence_list = uniq( sort { $a cmp $b } values %{ $component->{'references'} } );
			my $evidence_count = 0;
			foreach my $evidence (@evidence_list) {
				$evidence_count++;
				my @refs;
				foreach my $ref ( sort { $a <=> $b } keys %{ $component->{'references'} } ) {
					push @refs, $ref if $component->{'references'}->{$ref} eq $evidence;
				}
				if (@refs) {
					if ( $evidence_count == 1 ) {
						$evidence = qq(data derived from $evidence assays);
					} elsif ( $evidence_count == @evidence_list ) {
						$evidence = qq(and $evidence assays);
					} else {
						$evidence = qq(, $evidence assays);
					}
					local $" = q(, PMID:);
					push @note_evidence, qq($evidence (PMID:@refs));
				}
			}
			if (@note_evidence) {
				local $" = q(; );
				push @$notes,
				  qq($component->{'locus'}: $component->{'variant'} is not )
				  . qq(cross-reactive with vaccine variant - @note_evidence);
			}
		}
	}
	return $matched;
}

sub allele_designated {
	my ( $isolate_id, $locus, $allele_id ) = @_;
	return $script->{'datastore'}->run_query(
		q(SELECT EXISTS(SELECT * FROM allele_designations WHERE (isolate_id,locus,allele_id)=(?,?,?))),
		[ $isolate_id, $locus, $allele_id ],
		{ cache => 'exact_match' }
	);
}

sub contains_trumenba_antigens {
	my ($isolate_id) = @_;
	return $script->{'datastore'}->run_query(
		q[SELECT EXISTS(SELECT * FROM allele_designations WHERE (isolate_id,locus)=(?,?))],
		[ $isolate_id, 'fHbp_peptide' ],
		{ cache => 'contains_trumenba_antigens' }
	);
}

sub has_exact_trumenba_match {
	my ($isolate_id) = @_;
	my $components = [
		{ locus => 'fHbp_peptide', variant => '45', references => { 20619376 => 'peptide sequence match' } },
		{ locus => 'fHbp_peptide', variant => '55', references => { 20619376 => 'peptide sequence match' } },
	];
	return format_match_results( $isolate_id, $components, 'is exact match to vaccine variant' );
}

sub cross_reacts_with_trumenba {
	my ($isolate_id) = @_;
	my $components = [
		{
			locus      => 'fHbp_peptide',
			variant    => '1',
			references => {
				22569484 => 'SBA',
				22718089 => 'SBA',
				22871351 => 'SBA',
				23114369 => 'SBA',
				23352429 => 'SBA',
				26407272 => 'SBA',
				26707218 => 'SBA',
				26803328 => 'SBA',
				26835974 => 'SBA',
				26974889 => 'SBA',
				27745812 => 'SBA',
				27846061 => 'SBA',
				28196734 => 'SBA',
				28566335 => 'SBA',
				29236639 => 'SBA',
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '4',
			references => {
				27846061 => 'SBA',
				29236639 => 'SBA',
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '13',
			references => {
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '14',
			references => {
				22569484 => 'SBA',
				22718089 => 'SBA',
				27846061 => 'SBA',
				28196734 => 'SBA',
				29236639 => 'SBA',
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '15',
			references => {
				22569484 => 'SBA',
				22718089 => 'SBA',
				23352429 => 'SBA',
				26407272 => 'SBA',
				26707218 => 'SBA',
				26803328 => 'SBA',
				26835974 => 'SBA',
				27745812 => 'SBA',
				27846061 => 'SBA',
				28196734 => 'SBA',
				28566335 => 'SBA',
				29236639 => 'SBA',
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '16',
			references => {
				27846061 => 'SBA',
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '19',
			references => {
				22569484 => 'SBA',
				22718089 => 'SBA',
				22871351 => 'SBA',
				23114369 => 'SBA',
				23352429 => 'SBA',
				26407272 => 'SBA',
				26707218 => 'SBA',
				26803328 => 'SBA',
				26835974 => 'SBA',
				26974889 => 'SBA',
				27745812 => 'SBA',
				27846061 => 'SBA',
				27846061 => 'SBA',
				28196734 => 'SBA',
				28566335 => 'SBA',
				29236639 => 'SBA',
				29236639 => 'SBA',
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '21',
			references => {
				27846061 => 'SBA',
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '23',
			references => {
				28566335 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '25',
			references => {
				28566335 => 'SBA',
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '30',
			references => {
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '30',
			references => {
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '47',
			references => {
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '49',
			references => {
				22871351 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '76',
			references => {
				28566335 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '87',
			references => {
				22569484 => 'SBA',
				22718089 => 'SBA',
				22871351 => 'SBA',
				23114369 => 'SBA',
				23352429 => 'SBA',
				27846061 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '121',
			references => {
				28196734 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '180',
			references => {
				22569484 => 'SBA',
				27846061 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '187',
			references => {
				22569484 => 'SBA',
				26407272 => 'SBA',
				26707218 => 'SBA',
				26803328 => 'SBA',
				26835974 => 'SBA',
				27745812 => 'SBA',
				27846061 => 'SBA',
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '252',
			references => {
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '276',
			references => {
				27846061 => 'SBA',
				28566335 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '510',
			references => {
				28566335 => 'SBA'
			}
		}
	];
	return format_match_results( $isolate_id, $components, 'is cross-reactive to vaccine variant' );
}

sub set_value {
	my ( $isolate_id, $field, $value ) = @_;
	my $value_set = $script->{'datastore'}->run_query(
		'SELECT EXISTS(SELECT * FROM eav_text WHERE (isolate_id,field)=(?,?))',
		[ $isolate_id, $field ],
		{ cache => 'set_value' }
	);
	if ($value_set) {
		eval {
			$script->{'db'}
			  ->do( 'UPDATE eav_text SET value=? WHERE (isolate_id,field)=(?,?)', undef, $value, $isolate_id, $field );
		};
		if ($@) {
			$script->{'db'}->rollback;
			die "$@\n";
		}
		$script->{'db'}->commit;
	} else {
		eval {
			$script->{'db'}->do( 'INSERT INTO eav_text (isolate_id,field,value) VALUES (?,?,?)',
				undef, $isolate_id, $field, $value );
		};
		if ($@) {
			$script->{'db'}->rollback;
			die "$@\n";
		}
		$script->{'db'}->commit;
	}
	return;
}
