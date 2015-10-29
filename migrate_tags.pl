#!/usr/bin/perl -T
#Move sequences, designations and tags between BIGSdb databases
#Corresponding isolate records need to exist in both databases prior to migration.
#Written by Keith Jolley 2011-2015
use strict;
use warnings;
use 5.010;
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	HOST             => 'localhost',
	PORT             => 5432,
	USER             => 'apache',
	PASSWORD         => undef
};
#######End Local configuration###############################
use lib (LIB_DIR);
use List::MoreUtils qw(any);
use Getopt::Std;
use BIGSdb::Utils;
use BIGSdb_Scripts::Migrate;
use Log::Log4perl qw(get_logger);

#Direct all library logging calls to screen
my $log_conf = q(
	log4perl.category.BIGSdb.Script        = INFO, Screen
	log4perl.category.BIGSdb.Dataconnector = WARN, Screen
	log4perl.category.BIGSdb.Datastore     = WARN, Screen
	log4perl.appender.Screen               = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr        = 1
    log4perl.appender.Screen.layout        = Log::Log4perl::Layout::SimpleLayout
);
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my %opts;
getopts( 'a:b:i:j:hqn', \%opts );
if ( $opts{'h'} ) {
	show_help();
	exit;
}
if ( any { !$opts{$_} } qw (a b i j) ) {
	say "Usage: migrate_tags.pl -a <source database config> -b <destination database config> -i <id> -j <id>\n";
	say "Help: migrate_tags.pl -h";
	exit;
}
$opts{'throw_busy_exception'} = 0;
my $script = BIGSdb_Scripts::Migrate->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		host             => HOST,
		port             => PORT,
		user             => USER,
		password         => PASSWORD,
		options          => \%opts,
		instance         => $opts{'a'},
		logger           => $logger
	}
);
my $i = $opts{'i'};    #source isolate id
my $j = $opts{'j'};    #destination isolate id

#Do some checks before actually doing anything
die "This script should only be called on isolate databases.\n" if any { $_ ne 'isolates' } $script->get_db_types;
die "Source isolate id must be an integer.\n"                   if !BIGSdb::Utils::is_int($i);
die "Destination isolate id must be an integer.\n"              if !BIGSdb::Utils::is_int($j);
die "Isolate id does not exist in source database.\n" if !$script->isolate_exists_in_source($i);
die "Isolate id does not exist in destination database.\n" if !$script->isolate_exists_in_destination($j);
my $seqs_exist_in_destination =
  $script->{'datastore'}
  ->run_query( "SELECT EXISTS(SELECT * FROM sequence_bin WHERE isolate_id=?)", $j, { db => $script->{'db2'}->{ $opts{'b'} } } );
die "Isolate in destination database already has sequences associated with record.\n" if $seqs_exist_in_destination;
my $designations_exist_in_destination =
  $script->{'datastore'}
  ->run_query( "SELECT EXISTS(SELECT * FROM allele_designations WHERE isolate_id=?)", $j, { db => $script->{'db2'}->{ $opts{'b'} } } );
die "Isolate in destination database already has allele designations set.\n" if $designations_exist_in_destination && !$opts{'n'};
my $missing_users = $script->get_missing_designation_users_in_destination($i);

if ( keys %$missing_users ) {
	print "Missing sender/curator in destination database:\n";
	print "$missing_users->{$_}\n" foreach keys %$missing_users;
	exit;
}
if ( !$opts{'n'} ) {
	my $designations =
	  $script->{'datastore'}
	  ->run_query( "SELECT * FROM allele_designations WHERE isolate_id=? ORDER BY locus", $i, { fetch => 'all_arrayref', slice => {} } );
	my $sql =
	  $script->{'db2'}->{ $opts{'b'} }->prepare( "INSERT INTO allele_designations (isolate_id,locus,allele_id,sender,status,method,"
		  . "curator,date_entered,datestamp,comments) VALUES (?,?,?,?,?,?,?,?,?,?)" );
	foreach (@$designations) {
		if ( $script->is_locus_in_destination( $_->{'locus'} ) ) {
			eval {
				$sql->execute(
					$j, $_->{'locus'}, $_->{'allele_id'}, $script->{'user_map'}->{ $_->{'sender'} },
					$_->{'status'}, $_->{'method'}, $script->{'user_map'}->{ $_->{'curator'} },
					$_->{'date_entered'}, $_->{'datestamp'}, $_->{'comments'}
				);
			};
			if ($@) {
				$script->{'db2'}->{ $opts{'b'} }->rollback;
				die "$@\n" if $@;
			}
		} else {
			print "Locus $_->{'locus'} not in destination - skipping.\n" if !$opts{'q'};
		}
	}
}
my $sequences =
  $script->{'datastore'}
  ->run_query( "SELECT * FROM sequence_bin WHERE isolate_id=? ORDER BY id", $i, { fetch => 'all_arrayref', slice => {} } );
my $sql =
  $script->{'db2'}->{ $opts{'b'} }->prepare( "INSERT INTO sequence_bin (id,isolate_id,sequence,method,original_designation,comments,"
	  . "sender,curator,date_entered,datestamp) VALUES (?,?,?,?,?,?,?,?,?,?)" );
my %sequence_map;
my %allele_sequence_map;
my $last_one;
foreach (@$sequences) {
	my $next = $script->get_next_id( 'sequence_bin', $last_one );
	$sequence_map{ $_->{'id'} } = $next;
	eval {
		$sql->execute(
			$next, $j, $_->{'sequence'}, $_->{'method'}, $_->{'original_designation'},
			$_->{'comments'},
			$script->{'user_map'}->{ $_->{'sender'} },
			$script->{'user_map'}->{ $_->{'curator'} },
			$_->{'date_entered'}, $_->{'datestamp'}
		);
	};
	if ($@) {
		$script->{'db2'}->{ $opts{'b'} }->rollback;
		die "$@\n" if $@;
	}
	$last_one = $next;
}
if ( !$opts{'n'} ) {
	my $allele_sequences =
	  $script->{'datastore'}
	  ->run_query( "SELECT allele_sequences.* FROM allele_sequences WHERE isolate_id=?", $i, { fetch => 'all_arrayref', slice => {} } );
	$sql =
	  $script->{'db2'}->{ $opts{'b'} }->prepare( "INSERT INTO allele_sequences (seqbin_id,locus,start_pos,end_pos,reverse,complete,curator,"
		  . "datestamp) VALUES (?,?,?,?,?,?,?,?)" );
	my $allele_seq_sql =
	  $script->{'db2'}->{ $opts{'b'} }->prepare("SELECT id FROM allele_sequences WHERE (seqbin_id,locus,start_pos,end_pos)=(?,?,?,?)");
	foreach (@$allele_sequences) {
		if ( $script->is_locus_in_destination( $_->{'locus'} ) ) {
			eval {
				$sql->execute(
					$sequence_map{ $_->{'seqbin_id'} },
					$_->{'locus'}, $_->{'start_pos'}, $_->{'end_pos'}, $_->{'reverse'}, $_->{'complete'},
					$script->{'user_map'}->{ $_->{'curator'} },
					$_->{'datestamp'}
				);
				$allele_seq_sql->execute( $sequence_map{ $_->{'seqbin_id'} }, $_->{'locus'}, $_->{'start_pos'}, $_->{'end_pos'} );
				my $new_id = $allele_seq_sql->fetchrow_array;
				if ($new_id) {
					$allele_sequence_map{ $_->{'id'} } = $new_id;
				}
			};
			if ($@) {
				$script->{'db2'}->{ $opts{'b'} }->rollback;
				die "$@\n" if $@;
			}
		} else {
			print "Locus $_->{'locus'} not in destination - skipping.\n" if !$opts{'q'};
		}
	}
	my $flags = $script->{'datastore'}->run_query(
		"SELECT allele_sequences.*,sequence_flags.* FROM sequence_flags LEFT JOIN allele_sequences ON sequence_flags.id = "
		  . "allele_sequences.id WHERE isolate_id=?",
		$i,
		{ fetch => 'all_arrayref', slice => {} }
	);
	$sql = $script->{'db2'}->{ $opts{'b'} }->prepare("INSERT INTO sequence_flags (id, flag, curator, datestamp) VALUES (?,?,?,?)");
	my $sql2 =
	  $script->{'db2'}->{ $opts{'b'} }->prepare("SELECT id FROM allele_sequences WHERE (seqbin_id,locus,start_pos,end_pos) = (?,?,?,?)");
	foreach (@$flags) {
		if ( $script->is_locus_in_destination( $_->{'locus'} ) ) {
			eval {
				$sql->execute(
					$allele_sequence_map{ $_->{'id'} },
					$_->{'flag'}, $script->{'user_map'}->{ $_->{'curator'} },
					$_->{'datestamp'}
				);
			};
			if ($@) {
				$script->{'db2'}->{ $opts{'b'} }->rollback;
				die "$@\n" if $@;
			}
		}
	}
}
eval {
	$script->{'db2'}->{ $opts{'b'} }->do(q(SELECT setval('sequence_bin_id_seq',(SELECT max(id) FROM sequence_bin))));
};
if ($@){
	die qq(Grant update permission to sequence_bin_id_seq: \n)
	  . qq(GRANT USAGE,SELECT,UPDATE ON SEQUENCE sequence_bin_id_seq TO apache;\n);
}
$script->{'db2'}->{ $opts{'b'} }->commit;

sub show_help {
	print << "HELP";

Usage migrate_tags.pl -a <source database config> -b <destination database config> -i <id> -j <id>

Options
-------
-a <name>  Source database configuration name.
-b <name>  Destination database configuration name.
-h         This help page.
-i <id>    Source database isolate id.
-j <id>    Destination database isolate id.
-n         Don't transfer tags
-q         Quiet


HELP
	return;
}
