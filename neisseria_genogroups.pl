#!/usr/bin/env perl
#Run capsule characterization script
#(https://github.com/ntopaz/characterize_neisseria_capsule) against genome
#records in BIGSdb Neisseria database
#Written by Keith Jolley
#Copyright (c) 2019, University of Oxford
#E-mail: keith.jolley@zoo.ox.ac.uk
#
#This is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#BIGSdb is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this software.  If not, see <http://www.gnu.org/licenses/>.
#
#Version: 20191004
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	TMP_DIR          => '/var/tmp',
	DBASE            => 'pubmlst_neisseria_isolates_all',
	SCRIPT_DIR       => '/usr/local/characterize_neisseria_capsule',
	REST_URL         => 'http://rest.pubmlst.org',
	USER_ID          => -3
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Constants qw(LOG_TO_SCREEN);
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
use File::Path qw(make_path remove_tree);

#Direct all library logging calls to screen
my $log_conf = LOG_TO_SCREEN;
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my %opts;
GetOptions(
	'database=s'   => \$opts{'d'},
	'help'         => \$opts{'h'},
	'i|isolates=s' => \$opts{'i'},
	'q|quiet'      => \$opts{'q'},
	'size=i'       => \$opts{'size'},
	'threads=i'    => \$opts{'threads'},
	'x|min=i'      => \$opts{'x'},
	'y|max=i'      => \$opts{'y'},
) or die("Error in command line arguments\n");
if ( $opts{'h'} ) {
	show_help();
	exit;
}
$opts{'threads'} //= 4;
my $script = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		options          => \%opts,
		instance         => $opts{'d'} // DBASE,
		logger           => $logger
	}
);
die "Script initialization failed - server too busy?\n"
  if !defined $script->{'db'};
die "This script can only be run against an isolate database.\n"
  if ( $script->{'system'}->{'dbtype'} // '' ) ne 'isolates';
main();
undef $script;

sub main {
	my $isolates = get_ids_with_no_genogroup();
	$isolates = $script->filter_and_sort_isolates($isolates);
	return if !@$isolates;
	my $dir        = create_blast_database();
	my $fasta_dir  = get_isolate_fasta_files($isolates);
	my $out_dir    = TMP_DIR . '/' . BIGSdb::Utils::get_random();
	my $script_dir = SCRIPT_DIR;
	my $quiet      = $opts{'q'} ? q(>/dev/null 2>&1) : q();
	system("python3 characterize_neisseria_capsule.py -d $fasta_dir -t $opts{'threads'} -o $out_dir $quiet");
	chdir TMP_DIR;
	remove_tree($dir);
	remove_tree($fasta_dir);
	my ($output_file) = glob("$out_dir/serogroup/*.tab");
	return if !$output_file;
	process_output($output_file);
	remove_tree($out_dir);
}

sub process_output {
	my ($output_file) = @_;
	my $allowed_genogroups = $script->{'xmlHandler'}->get_field_option_list('genogroup');
	my %allowed_genogroups = map { $_ => 1 } @$allowed_genogroups;
	open( my $fh, '<', $output_file ) or die "Cannot open $output_file for reading.\n";
	while ( my $line = <$fh> ) {
		next if $line =~ /^Query/x;
		my ( $id, $genogroup, undef, $notes ) = split /\t/x, $line;
		next if !$id;
		$genogroup = 'cnl' if $genogroup eq 'cnl_1';
		if ( !$allowed_genogroups{$genogroup} ) {
			say qq(id: $id; Invalid genogroup: $genogroup);
			next;
		}
		say qq($id\t$genogroup\t$notes);
		eval {
			if ($notes) {
				$notes =~ s/\s+$//gx;
				$notes .= '. ';
			}
			$notes .= 'Prediction code: https://github.com/ntopaz/characterize_neisseria_capsule.';
			$script->{'db'}->do( 'UPDATE isolates SET (genogroup,genogroup_notes)=(?,?)', undef, $genogroup, $notes );
			$script->{'db'}->do( 'INSERT INTO history (isolate_id,timestamp,action,curator) VALUES (?,?,?,?)',
				undef, $id, 'now', "genogroup: '' -> '$genogroup'", USER_ID );
		};
		if ($@) {
			$script->{'db'}->rollback;
			die "$@\n";
		}
		$script->{'db'}->commit;
	}
	close $fh;
}

sub get_ids_with_no_genogroup {
	my $no_genogroup =
	  $script->{'datastore'}->run_query(
		"SELECT id FROM $script->{'system'}->{'view'} WHERE species='Neisseria meningitidis' AND genogroup IS NULL",
		undef, { fetch => 'col_arrayref' } );
	my %no_genogroup = map { $_ => 1 } @$no_genogroup;
	my $isolates = $script->get_isolates_with_linked_seqs( { size => $opts{'size'} // 2_000_000 } );
	my $filtered_list = [];
	foreach my $id (@$isolates) {
		push @$filtered_list, $id if $no_genogroup{$id};
	}
	return $filtered_list;
}

sub get_isolate_fasta_files {
	my ($isolates) = @_;
	my $dir = TMP_DIR . '/' . BIGSdb::Utils::get_random();
	make_path($dir);
	foreach my $id (@$isolates) {
		my $contig_ids =
		  $script->{'datastore'}
		  ->run_query( 'SELECT id FROM sequence_bin WHERE isolate_id=?', $id, { fetch => 'col_arrayref' } );
		my $contigs  = $script->{'contigManager'}->get_contigs_by_list($contig_ids);
		my $filename = "$dir/$id.fas";
		open( my $fh, '>', $filename ) || die "Cannot open $filename for writing.\n";
		foreach my $contig_id ( keys %$contigs ) {
			say $fh qq(>$contig_id);
			say $fh qq($contigs->{$contig_id});
		}
		close $fh;
	}
	return $dir;
}

sub create_blast_database {
	my $dir        = TMP_DIR . '/' . BIGSdb::Utils::get_random();
	my $script_dir = SCRIPT_DIR;
	make_path($dir);
	system "cp -R $script_dir/* $dir";
	chdir $dir;
	my $quiet = $opts{'q'} ? q(>/dev/null) : q();
	system("python3 $dir/build_neisseria_dbs.py$quiet");
	return $dir;
}

sub show_help {
	my $termios = POSIX::Termios->new;
	$termios->getattr;
	my $ospeed = $termios->getospeed;
	my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
	my ( $norm, $bold, $under ) = map { $t->Tputs( $_, 1 ) } qw(me md us);
	say << "HELP";
${bold}NAME$norm
    ${bold}neisseria_genogroups.pl$norm - Run capsule characterization script against 
    Neisseria database

${bold}SYNOPSIS$norm
    ${bold}neisseria_genogroups.pl $norm [${under}options$norm]

${bold}OPTIONS$norm

${bold}--database$norm ${under}NAME$norm
    Database configuration name.
      
${bold}--help$norm
    This help page.
    
${bold}--isolates$norm ${under}LIST$norm  
    Comma-separated list of isolate ids to scan.
    
${bold}--size$norm ${under}SIZE$norm
    Limit isolates to a sequence bin size >= value set (default 2,000,000 bp).
    
${bold}--threads$norm ${under}THREADS$norm
    Number of threads to use for BLAST query.
    
${bold}-x, --min$norm ${under}ID$norm
    Minimum isolate id.

${bold}-y, --max$norm ${under}ID$norm
    Maximum isolate id.
    
HELP
	return;
}
