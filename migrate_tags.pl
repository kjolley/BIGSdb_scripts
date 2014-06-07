#!/usr/bin/perl -T
#Move sequences, designations and tags between BIGSdb databases
#Corresponding isolate records need to exist in both databases prior to migration.
#Written by Keith Jolley 2011
use strict;
use warnings;
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	HOST             => 'localhost',
	PORT             => 5432,
	USER             => 'apache',
	PASSWORD         => ''
};
#######End Local configuration###############################
use lib (LIB_DIR);
use List::MoreUtils qw(any);
use Getopt::Std;
use BIGSdb::Utils;
use BIGSdb_Scripts::Migrate;
my %opts;
getopts( 'a:b:i:j:hqn', \%opts );

if ( $opts{'h'} ) {
	show_help();
	exit;
}
if ( any { !$opts{$_} } qw (a b i j) ) {
	print "\nUsage: migrate_tags.pl -a <source database config> -b <destination database config> -i <id> -j <id>\n\n";
	print "Help: migrate_tags.pl -h\n";
	exit;
}
$opts{'throw_busy_exception'} = 1;
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
		writable         => 1
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
my $seqs_exist_in_destination = $script->run_simple_query( $opts{'b'}, "SELECT COUNT(*) FROM sequence_bin WHERE isolate_id=?", $j )->[0];
die "Isolate in destination database already has sequences associated with record.\n" if $seqs_exist_in_destination;
my $designations_exist_in_destination =
  $script->run_simple_query( $opts{'b'}, "SELECT COUNT(*) FROM allele_designations WHERE isolate_id=?", $j )->[0];
die "Isolate in destination database already has allele designations set.\n" if $designations_exist_in_destination && !$opts{'n'};
my $missing_users = $script->get_missing_designation_users_in_destination($i);

if ( keys %$missing_users ) {
	print "Missing sender/curator in destination database:\n";
	print "$missing_users->{$_}\n" foreach keys %$missing_users;
	exit;
}
if ( !$opts{'n'} ) {
	my $designations =
	  $script->{'datastore'}->run_list_query_hashref( "SELECT * FROM allele_designations WHERE isolate_id=? ORDER BY locus", $i );
	my $sql =
	  $script->{'db2'}->{ $opts{'b'} }->prepare(
"INSERT INTO allele_designations (isolate_id,locus,allele_id,sender,status,method,curator,date_entered,datestamp,comments) VALUES (?,?,?,?,?,?,?,?,?,?)"
	  );
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
#	my $pending =
#	  $script->{'datastore'}->run_list_query_hashref( "SELECT * FROM pending_allele_designations WHERE isolate_id=? ORDER BY locus", $i );
#	$sql =
#	  $script->{'db2'}->{ $opts{'b'} }->prepare(
#"INSERT INTO pending_allele_designations (isolate_id,locus,allele_id,sender,method,curator,date_entered,datestamp,comments) VALUES (?,?,?,?,?,?,?,?,?)"
#	  );
#	foreach (@$pending) {
#		if ( $script->is_locus_in_destination( $_->{'locus'} ) ) {
#			eval {
#				$sql->execute(
#					$j,                                        $_->{'locus'},     $_->{'allele_id'},
#					$script->{'user_map'}->{ $_->{'sender'} }, $_->{'method'},    $script->{'user_map'}->{ $_->{'curator'} },
#					$_->{'date_entered'},                      $_->{'datestamp'}, $_->{'comments'}
#				);
#			};
#			if ($@) {
#				$script->{'db2'}->{ $opts{'b'} }->rollback;
#				die "$@\n" if $@;
#			}
#		}
#	}
}
my $sequences = $script->{'datastore'}->run_list_query_hashref( "SELECT * FROM sequence_bin WHERE isolate_id=? ORDER BY id", $i );
my $sql =
  $script->{'db2'}->{ $opts{'b'} }->prepare(
"INSERT INTO sequence_bin (id, isolate_id, sequence, method, original_designation, comments, sender, curator, date_entered, datestamp) VALUES (?,?,?,?,?,?,?,?,?,?)"
  );
my %sequence_map;
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
	  $script->{'datastore'}->run_list_query_hashref( "SELECT allele_sequences.* FROM allele_sequences WHERE isolate_id=?", $i );
	$sql =
	  $script->{'db2'}->{ $opts{'b'} }->prepare(
"INSERT INTO allele_sequences (seqbin_id, locus, start_pos, end_pos, reverse, complete, curator, datestamp) VALUES (?,?,?,?,?,?,?,?)"
	  );
	foreach (@$allele_sequences) {
		if ( $script->is_locus_in_destination( $_->{'locus'} ) ) {
			eval {
				$sql->execute(
					$sequence_map{ $_->{'seqbin_id'} },
					$_->{'locus'}, $_->{'start_pos'}, $_->{'end_pos'}, $_->{'reverse'}, $_->{'complete'},
					$script->{'user_map'}->{ $_->{'curator'} },
					$_->{'datestamp'}
				);
			};
			if ($@) {
				$script->{'db2'}->{ $opts{'b'} }->rollback;
				die "$@\n" if $@;
				$opts{'throw_busy_exception'} = 1;
			}
		} else {
			print "Locus $_->{'locus'} not in destination - skipping.\n" if !$opts{'q'};
		}
	}
	my $flags =
	  $script->{'datastore'}->run_list_query_hashref(
		"SELECT sequence_flags.* FROM sequence_flags LEFT JOIN sequence_bin ON seqbin_id = sequence_bin.id WHERE isolate_id=?", $i );
	$sql =
	  $script->{'db2'}->{ $opts{'b'} }
	  ->prepare("INSERT INTO sequence_flags (seqbin_id, locus, start_pos, end_pos, flag, curator, datestamp) VALUES (?,?,?,?,?,?,?)");
	foreach (@$flags) {
		if ( $script->is_locus_in_destination( $_->{'locus'} ) ) {
			eval {
				$sql->execute(
					$sequence_map{ $_->{'seqbin_id'} },
					$_->{'locus'}, $_->{'start_pos'}, $_->{'end_pos'}, $_->{'flag'}, $script->{'user_map'}->{ $_->{'curator'} },
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
