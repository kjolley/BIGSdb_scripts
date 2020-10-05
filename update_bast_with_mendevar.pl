#!/usr/bin/env perl
#Predict cross-reactivity to Bexsero and Trumenba vaccine components
#and update BAST scheme
#Written by Keith Jolley
#Copyright (c) 2020, University of Oxford
#E-mail: keith.jolley@zoo.ox.ac.uk
#Version: 20200930
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use FindBin;
use lib "$FindBin::Bin/lib";
use Getopt::Long qw(:config no_ignore_case);
use List::MoreUtils qw(uniq);
use BIGSdb::Offline::Script;
use MenVaccine;
use constant LOCI => qw(fHbp_peptide NadA_peptide NHBA_peptide PorA_VR2);

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
	'd|database=s' => \$opts{'database'},
	'q|quiet'      => \$opts{'quiet'},
	's|scheme=i'   => \$opts{'scheme'}
) or die("Error in command line arguments\n");
$opts{'database'} //= 'pubmlst_neisseria_seqdef';
$opts{'scheme'}   //= 53;
my $script = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		options          => \%opts,
		instance         => $opts{'database'},
		logger           => $logger
	}
);
die "This script can only be run against a sequence definition database.\n"
  if ( $script->{'system'}->{'dbtype'} // '' ) ne 'sequences';
main();
undef $script;

sub main {
	my $indices =
	  $script->{'datastore'}->run_query( 'SELECT locus,index FROM scheme_warehouse_indices WHERE scheme_id=?',
		$opts{'scheme'}, { fetch => 'all_arrayref', slice => {} } );
	my %locus_index = map { $_->{'locus'} => $_->{'index'} } @$indices;
	my $bast_profiles =
	  $script->{'datastore'}->run_query( "SELECT * FROM mv_scheme_$opts{'scheme'} ORDER BY CAST(BAST AS int)",
		undef, { fetch => 'all_arrayref', slice => {} } );
	eval {
		$script->{'db'}->do( 'DELETE FROM profile_fields WHERE scheme_id=? AND scheme_field=?',
			undef, $opts{'scheme'}, 'MenDeVAR_Bexsero_reactivity' );
		$script->{'db'}->do( 'DELETE FROM profile_fields WHERE scheme_id=? AND scheme_field=?',
			undef, $opts{'scheme'}, 'MenDeVAR_Trumenba_reactivity' );
		foreach my $profile (@$bast_profiles) {
			my $variants = {};
			foreach my $locus (LOCI) {
				$variants->{$locus} = $profile->{'profile'}->[ $locus_index{$locus} - 1 ];
			}
			my $bexsero;
			if ( has_exact_bexsero_match($variants) ) {
				$bexsero = 'exact match';
			} elsif ( cross_reacts_with_bexsero($variants) ) {
				$bexsero = 'cross-reactive';
			} elsif ( no_reactivity_with_bexsero($variants) ) {
				$bexsero = 'none';
			} else {
				$bexsero = 'insufficient data';
			}
			my $trumenba;
			if ( has_exact_trumenba_match($variants) ) {
				$trumenba = 'exact match';
			} elsif ( cross_reacts_with_trumenba($variants) ) {
				$trumenba = 'cross-reactive';
			} elsif ( $variants->{'fHbp_peptide'} eq '0' ) {
				$trumenba = 'none';
			} else {
				$trumenba = 'insufficient data';
			}
			$script->{'db'}->do(
				'INSERT INTO profile_fields (scheme_id,scheme_field,profile_id,value,curator,datestamp) '
				  . 'VALUES (?,?,?,?,?,?)', undef, $opts{'scheme'}, 'MenDeVAR_Bexsero_reactivity', $profile->{'bast'}, $bexsero,
				0, 'now'
			);
			$script->{'db'}->do(
				'INSERT INTO profile_fields (scheme_id,scheme_field,profile_id,value,curator,datestamp) '
				  . 'VALUES (?,?,?,?,?,?)', undef, $opts{'scheme'}, 'MenDeVAR_Trumenba_reactivity', $profile->{'bast'}, $trumenba,
				0, 'now'
			);
			say qq($profile->{'bast'}\t$bexsero\t$trumenba) if !$opts{'quiet'};
		}
	};
	if ($@) {
		$script->{'db'}->rollback;
		die "$@\n";
	}
	$script->{'db'}->commit;
}

sub has_exact_bexsero_match {
	my ($variants) = @_;
	my $components = MenVaccine::get_exact_bexsero_match();
	foreach my $component (@$components) {
		my $locus = $component->{'locus'};
		return 1 if $variants->{$locus} eq $component->{'variant'};
	}
	return;
}

sub cross_reacts_with_bexsero {
	my ($variants) = @_;
	my $components = MenVaccine::get_cross_reacts_with_bexsero();
	foreach my $component (@$components) {
		my $locus = $component->{'locus'};
		return 1 if $variants->{$locus} eq $component->{'variant'};
	}
	return;
}

sub no_reactivity_with_bexsero {
	my ($variants) = @_;

	#FHbp
	my $components = MenVaccine::get_fhbp_no_reactivity_with_bexsero();
	my $matched    = 0;
	foreach my $component (@$components) {
		$matched = 1
		  if $variants->{ $component->{'locus'} } eq $component->{'variant'}
		  || $variants->{ $component->{'locus'} } eq '0';
	}
	return if !$matched;

	#NHBA
	$components = MenVaccine::get_nhba_no_reactivity_with_bexsero();
	$matched    = 0;
	foreach my $component (@$components) {
		$matched = 1
		  if $variants->{ $component->{'locus'} } eq $component->{'variant'}
		  || $variants->{ $component->{'locus'} } eq '0';
	}
	return if !$matched;

	#NadA
	$components = MenVaccine::get_nadA_no_reactivity_with_bexsero();
	$matched    = 0;
	foreach my $component (@$components) {
		$matched = 1
		  if $variants->{ $component->{'locus'} } eq $component->{'variant'}
		  || $variants->{ $component->{'locus'} } eq '0';
	}
	return if !$matched;

	#PorA
	$matched = 0;
	$matched = 1 if $variants->{'PorA_VR2'} ne '4' && $variants->{'PorA_VR2'} ne '0';
	return if !$matched;
	return 1;
}

sub has_exact_trumenba_match {
	my ($variants) = @_;
	my $components = MenVaccine::get_exact_trumenba_match();
	foreach my $component (@$components) {
		my $locus = $component->{'locus'};
		return 1 if $variants->{$locus} eq $component->{'variant'};
	}
	return;
}

sub cross_reacts_with_trumenba {
	my ($variants) = @_;
	my $components = MenVaccine::get_cross_reacts_with_trumenba();
	foreach my $component (@$components) {
		my $locus = $component->{'locus'};
		return 1 if $variants->{$locus} eq $component->{'variant'};
	}
	return;
}

sub no_reactivity_with_trumenba {
}
