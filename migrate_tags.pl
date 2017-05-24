#!/usr/bin/perl -T
#Move sequences, designations and tags between BIGSdb databases
#Corresponding isolate records need to exist in both databases prior to migration.
#Written by Keith Jolley 2011-2017
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
my $log_conf = <<"LOG";
	log4perl.category.BIGSdb.Script        = INFO, Screen
	log4perl.category.BIGSdb.Dataconnector = WARN, Screen
	log4perl.category.BIGSdb.Datastore     = WARN, Screen
	log4perl.appender.Screen               = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr        = 1
    log4perl.appender.Screen.layout        = Log::Log4perl::Layout::SimpleLayout
LOG
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my %opts;
getopts( 'a:b:i:j:u:hqn', \%opts );
if ( $opts{'h'} ) {
	show_help();
	exit;
}
if ( any { !$opts{$_} } qw (a b i j u) ) {
	say
qq(Usage: migrate_tags.pl -a <source database config> -b <destination database config> -i <id> -j <id> -u <user>\n);
	say q(Help: migrate_tags.pl -h);
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
my $u = $opts{'u'};    #User id

#Do some checks before actually doing anything
die "This script should only be called on isolate databases.\n" if any { $_ ne 'isolates' } $script->get_db_types;
die "Source isolate id must be an integer.\n"              if !BIGSdb::Utils::is_int($i);
die "Destination isolate id must be an integer.\n"         if !BIGSdb::Utils::is_int($j);
die "User id must be an integer.\n"                        if !BIGSdb::Utils::is_int($u);
die "Isolate id does not exist in source database.\n"      if !$script->isolate_exists_in_source($i);
die "Isolate id does not exist in destination database.\n" if !$script->isolate_exists_in_destination($j);
main();

sub main {
	say "Source id: $i; Destination id: $j";
	my $seqs_exist_in_destination =
	  $script->{'datastore'}->run_query( 'SELECT EXISTS(SELECT * FROM sequence_bin WHERE isolate_id=?)',
		$j, { db => $script->{'db2'}->{ $opts{'b'} } } );
	die "Isolate in destination database already has sequences associated with record.\n" if $seqs_exist_in_destination;
	my $designations_exist_in_destination =
	  $script->{'datastore'}->run_query( 'SELECT EXISTS(SELECT * FROM allele_designations WHERE isolate_id=?)',
		$j, { db => $script->{'db2'}->{ $opts{'b'} } } );
	die "Isolate in destination database already has allele designations set.\n"
	  if $designations_exist_in_destination && !$opts{'n'};
	if ( !$opts{'n'} ) {
		my $designations =
		  $script->{'datastore'}->run_query( 'SELECT * FROM allele_designations WHERE isolate_id=? ORDER BY locus',
			$i, { fetch => 'all_arrayref', slice => {} } );
		my $sql =
		  $script->{'db2'}->{ $opts{'b'} }
		  ->prepare( 'INSERT INTO allele_designations (isolate_id,locus,allele_id,sender,status,method,'
			  . 'curator,date_entered,datestamp,comments) VALUES (?,?,?,?,?,?,?,?,?,?)' );
		foreach my $des (@$designations) {
			if ( $script->is_locus_in_destination( $des->{'locus'} ) ) {
				eval {
					$sql->execute(
						$j,                  $des->{'locus'},  $des->{'allele_id'}, $u,
						$des->{'status'},    $des->{'method'}, $u,                  $des->{'date_entered'},
						$des->{'datestamp'}, $des->{'comments'}
					);
				};
				if ($@) {
					$script->{'db2'}->{ $opts{'b'} }->rollback;
					die "$@\n" if $@;
				}
			} else {
				print "Locus $des->{'locus'} not in destination - skipping.\n" if !$opts{'q'};
			}
		}
	}
	my $sequences = $script->{'datastore'}->run_query( 'SELECT * FROM sequence_bin WHERE isolate_id=? ORDER BY id',
		$i, { fetch => 'all_arrayref', slice => {} } );
	my $sql =
	  $script->{'db2'}->{ $opts{'b'} }
	  ->prepare( 'INSERT INTO sequence_bin (isolate_id,sequence,method,original_designation,comments,'
		  . 'sender,curator,date_entered,datestamp) VALUES (?,?,?,?,?,?,?,?,?) RETURNING id' );
	my %sequence_map;
	my %allele_sequence_map;
	foreach my $seq (@$sequences) {
		eval {
			$sql->execute(
				$j,                             $seq->{'sequence'},     $seq->{'method'},
				$seq->{'original_designation'}, $seq->{'comments'},     $u,
				$u,                             $seq->{'date_entered'}, $seq->{'datestamp'}
			);
		};
		my ($last_seqbin_id) = $sql->fetchrow_array;
		$sequence_map{ $seq->{'id'} } = $last_seqbin_id;
		if ($@) {
			$script->{'db2'}->{ $opts{'b'} }->rollback;
			die "$@\n" if $@;
		}
	}
	if ( !$opts{'n'} ) {
		my $allele_sequences =
		  $script->{'datastore'}->run_query( 'SELECT allele_sequences.* FROM allele_sequences WHERE isolate_id=?',
			$i, { fetch => 'all_arrayref', slice => {} } );
		$sql =
		  $script->{'db2'}->{ $opts{'b'} }
		  ->prepare( 'INSERT INTO allele_sequences (seqbin_id,locus,start_pos,end_pos,reverse,complete,curator,'
			  . 'datestamp) VALUES (?,?,?,?,?,?,?,?) RETURNING id' );
		foreach my $allele_seq (@$allele_sequences) {
			if ( $script->is_locus_in_destination( $allele_seq->{'locus'} ) ) {
				eval {
					$sql->execute(
						$sequence_map{ $allele_seq->{'seqbin_id'} }, $allele_seq->{'locus'},
						$allele_seq->{'start_pos'},                  $allele_seq->{'end_pos'},
						$allele_seq->{'reverse'},                    $allele_seq->{'complete'},
						$u,                                          $allele_seq->{'datestamp'}
					);
					my ($new_id) = $sql->fetchrow_array;
					if ($new_id) {
						$allele_sequence_map{ $allele_seq->{'id'} } = $new_id;
					}
				};
				if ($@) {
					$script->{'db2'}->{ $opts{'b'} }->rollback;
					die "$@\n" if $@;
				}
			} else {
				print "Locus $allele_seq->{'locus'} not in destination - skipping.\n" if !$opts{'q'};
			}
		}
		my $flags = $script->{'datastore'}->run_query(
			'SELECT allele_sequences.*,sequence_flags.* FROM sequence_flags LEFT JOIN '
			  . 'allele_sequences ON sequence_flags.id=allele_sequences.id WHERE isolate_id=?',
			$i,
			{ fetch => 'all_arrayref', slice => {} }
		);
		$sql =
		  $script->{'db2'}->{ $opts{'b'} }
		  ->prepare('INSERT INTO sequence_flags (id, flag, curator, datestamp) VALUES (?,?,?,?)');
		my $sql2 =
		  $script->{'db2'}->{ $opts{'b'} }
		  ->prepare('SELECT id FROM allele_sequences WHERE (seqbin_id,locus,start_pos,end_pos) = (?,?,?,?)');
		foreach my $flag (@$flags) {
			if ( $script->is_locus_in_destination( $flag->{'locus'} ) ) {
				eval {
					$sql->execute( $allele_sequence_map{ $flag->{'id'} }, $flag->{'flag'}, $u, $flag->{'datestamp'} );
				};
				if ($@) {
					$script->{'db2'}->{ $opts{'b'} }->rollback;
					die "$@\n" if $@;
				}
			}
		}
	}
	eval {
		$script->{'db2'}->{ $opts{'b'} }
		  ->do(q(SELECT setval('sequence_bin_id_seq',(SELECT max(id) FROM sequence_bin))));
	};
	if ($@) {
		die qq(Grant update permission to sequence_bin_id_seq: \n)
		  . qq(GRANT USAGE,SELECT,UPDATE ON SEQUENCE sequence_bin_id_seq TO apache;\n);
	}
	$script->{'db2'}->{ $opts{'b'} }->commit;
	return;
}

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
-u <id>    User id for sender/curator

HELP
	return;
}
