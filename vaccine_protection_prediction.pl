#!/usr/bin/perl -T
#Predict cross-reactivity to Bexsero and Trumenba vaccine components
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
		if ( has_exact_bexsero_match($isolate_id) ) {
			say qq($isolate_id Bexsero exact match) if !$opts{'quiet'};
			set_value( $isolate_id, 'Bexsero_reactivity', 'exact match' );
			next;
		}
		if ( cross_reacts_with_bexsero($isolate_id) ) {
			say qq($isolate_id Bexsero cross-reactive) if !$opts{'quiet'};
			set_value( $isolate_id, 'Bexsero_reactivity', 'cross-reactive' );
			next;
		}
		say qq($isolate_id Bexsero none) if !$opts{'quiet'};
		set_value( $isolate_id, 'Bexsero_reactivity', 'none' );
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
		if ( has_exact_trumenba_match($isolate_id) ) {
			say qq($isolate_id Trumenba exact match) if !$opts{'quiet'};
			set_value( $isolate_id, 'Trumenba_reactivity', 'exact match' );
			next;
		}
		if ( cross_reacts_with_trumenba($isolate_id) ) {
			say qq($isolate_id Trumenba cross-reactive) if !$opts{'quiet'};
			set_value( $isolate_id, 'Trumenba_reactivity', 'cross-reactive' );
			next;
		}
		say qq($isolate_id Trumenba none) if !$opts{'quiet'};
		set_value( $isolate_id, 'Trumenba_reactivity', 'none' );
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
		{ locus => 'fHbp_peptide', variant => '1' },
		{ locus => 'NHBA_peptide', variant => '2' },
		{ locus => 'NadA_peptide', variant => '8' },
		{ locus => 'PorA_VR2',     variant => '4' }
	];
	foreach my $component (@$components) {
		return 1
		  if $script->{'datastore'}->run_query(
			q(SELECT EXISTS(SELECT * FROM allele_designations WHERE (isolate_id,locus,allele_id)=(?,?,?))),
			[ $isolate_id, $component->{'locus'}, $component->{'variant'} ],
			{ cache => 'exact_bexsero_match' }
		  );
	}
	return;
}

sub cross_reacts_with_bexsero {
	my ($isolate_id) = @_;
	my $designations = $script->{'datastore'}->run_query(
		q[SELECT locus,allele_id FROM allele_designations WHERE isolate_id=? AND locus IN ]
		  . q[('fHbp_peptide','NadA_peptide')],
		$isolate_id,
		{ fetch => 'all_arrayref', slice => {}, cache => 'cross_reacts_with_bexsero' }
	);
	foreach my $designation (@$designations) {
		if ( $designation->{'locus'} eq 'fHbp_peptide' ) {
			my %fhbp_match = map { $_ => 1 } qw(4 13 14 15 37 232);
			return 1 if $fhbp_match{ $designation->{'allele_id'} };
		}
		if ( $designation->{'locus'} eq 'NadA_peptide' ) {
			my %nadA_match = map { $_ => 1 } qw(1 2 3 4 5 6 7 9 25 26 28 58 85 86 91 101 103 105
			  106 108 110 112 113 115 118 119 120 121 127 128 129 130 131 132 133 134 135 136 137
			  141 142 143 145 146 147 148 149 150 152 153 154 155 156 157 162 163 164 165 168 169
			  170 173 174 175 178 180 181
			);
			return 1 if $nadA_match{ $designation->{'allele_id'} };
		}
	}
	return;
}
sub contains_trumenba_antigens {
	my ($isolate_id) = @_;
	return $script->{'datastore'}->run_query(
		q[SELECT EXISTS(SELECT * FROM allele_designations WHERE (isolate_id,locus)=(?,?))],
		[$isolate_id,'fHbp_peptide'],
		{ cache => 'contains_trumenba_antigens' }
	);
}


sub has_exact_trumenba_match {
	my ($isolate_id) = @_;
	return $script->{'datastore'}->run_query(
		q[SELECT EXISTS(SELECT * FROM allele_designations WHERE (isolate_id,locus)=(?,?) ]
		  . q[AND allele_id IN ('45','55'))],
		[ $isolate_id, 'fHbp_peptide' ],
		{ cache => 'has_exact_trumenba_match' }
	);
}

sub cross_reacts_with_trumenba {
	my ($isolate_id) = @_;
	return $script->{'datastore'}->run_query(
		    q[SELECT EXISTS(SELECT * FROM allele_designations WHERE (isolate_id,locus)=(?,?) ]
		  . q[AND allele_id IN ('1','4','13','14','15','16','19','21','23','24','25','30',]
		  . q['47','76','87','180','187','252','276','510'))],
		[ $isolate_id, 'fHbp_peptide' ],
		{ cache => 'cross_reacts_with_trumenba' }
	);
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
			$script->{'db'}
			  ->do( 'INSERT INTO eav_text (isolate_id,field,value) VALUES (?,?,?)', undef, $isolate_id, $field,
				$value );
		};
		if ($@) {
			$script->{'db'}->rollback;
			die "$@\n";
		}
		$script->{'db'}->commit;
	}
	return;
}
