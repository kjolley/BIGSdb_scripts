#!/usr/bin/env perl
#Run capsule characterization script
#(https://github.com/scwatts/hicap) against genome
#records in BIGSdb Haemophilus database
#Written by Keith Jolley
#Copyright (c) 2023, University of Oxford
#E-mail: keith.jolley@biology.ox.ac.uk
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
#Version: 20231010
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	DBASE            => 'pubmlst_hinfluenzae_isolates',
	ENV_SETUP        => 'source /home/bigsdb/miniconda3/bin/activate hicap',
	SCRIPT_DIR       => '/home/bigsdb/miniconda3/envs/hicap/bin'
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Constants qw(LOG_TO_SCREEN);
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
use JSON;

#Direct all library logging calls to screen
my $log_conf = LOG_TO_SCREEN;
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my %opts;
GetOptions(
	'database=s'   => \$opts{'d'},
	'help'         => \$opts{'h'},
	'i|isolates=s' => \$opts{'i'},
	'l|limit=i'    => \$opts{'limit'},
	'q|quiet'      => \$opts{'q'},
	'size=i'       => \$opts{'size'},
	'x|min=i'      => \$opts{'x'},
	'y|max=i'      => \$opts{'y'},
) or die("Error in command line arguments\n");
if ( $opts{'h'} ) {
	show_help();
	exit;
}
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
	my $isolates = get_ids_with_no_hicap_analysis();
	return if !@$isolates;
	if ( defined $opts{'limit'} && $opts{'limit'} > 0 ) {
		my $upper_bound = $opts{'limit'} > @$isolates ? @$isolates - 1 : $opts{'limit'} - 1;
		@$isolates = @{$isolates}[ 0 .. $upper_bound ];
	}
	my $prefix = BIGSdb::Utils::get_random();
	my $quiet  = $opts{'q'} ? q(>/dev/null 2>&1) : q();
	foreach my $isolate_id (@$isolates) {
		my $fasta = create_fasta_file( $prefix, $isolate_id );
		my $version;
		eval {
			system (ENV_SETUP);
			system( SCRIPT_DIR
				  . "/hicap --query_fp $fasta --output_dir $script->{'config'}->{'secure_tmp_dir'} $quiet" );
			my $version_cmd = SCRIPT_DIR . '/hicap --version';
			$version = `$version_cmd`;
			chomp $version;
		};
		my $tsv = "$script->{'config'}->{'secure_tmp_dir'}/${prefix}_$isolate_id.tsv";
		my $data;
		
		if ( -e $tsv ) {
			$data = read_tsv($tsv);
		} else {
			$data = {predicted_serotype => 'NT'};
		}
		$data->{'version'} = $version;
		my $svg = "$script->{'config'}->{'secure_tmp_dir'}/${prefix}_$isolate_id.svg";
		if ( -e $svg ) {
			my $contents_ref = BIGSdb::Utils::slurp($svg);
			$data->{'svg'} = $$contents_ref;
		}
		$script->delete_temp_files("$script->{'config'}->{'secure_tmp_dir'}/$prefix*");
		eval {
			$script->{'db'}->do('INSERT INTO analysis_results (name,isolate_id,datestamp,results) VALUES (?,?,?,?)',undef,'hicap',$isolate_id,
			'now',encode_json($data));
		};
		if ($@){
			$logger->error($@);
			$script->{'db'}->rollback;
		} else {
			$script->{'db'}->commit;
		}
	}
}

sub read_tsv {
	my ($tsv_file) = @_;
	open( my $fh, '<', $tsv_file ) || die "Cannot read TSV file $tsv_file.\n";
	my $header = <$fh>;
	chomp $header;
	my @fields = split /\t/x, $header;
	my $values = <$fh>;
	chomp $values;
	my @values = split /\t/x, $values;
	close $fh;
	my $data = {};

	for ( my $i = 0 ; $i < @fields ; $i++ ) {
		next if $fields[$i] eq '#isolate';
		$data->{ $fields[$i] } = $values[$i];
	}
	return $data;
}

sub create_fasta_file {
	my ( $prefix, $isolate_id ) = @_;
	my $contig_ids =
	  $script->{'datastore'}
	  ->run_query( 'SELECT id FROM sequence_bin WHERE isolate_id=?', $isolate_id, { fetch => 'col_arrayref' } );
	my $contigs  = $script->{'contigManager'}->get_contigs_by_list($contig_ids);
	my $filename = "$script->{'config'}->{'secure_tmp_dir'}/${prefix}_$isolate_id.fas";
	open( my $fh, '>', $filename ) || die "Cannot open $filename for writing.\n";
	foreach my $contig_id ( keys %$contigs ) {
		say $fh qq(>$contig_id);
		say $fh qq($contigs->{$contig_id});
	}
	close $fh;
	return $filename;
}

sub get_ids_with_no_hicap_analysis {

	#Use script method as this includes filtering options we may wish to use.
	my $isolates = $script->get_isolates_with_linked_seqs( { size => $opts{'size'} // 1_600_000 } );
	my $hicap    = $script->{'datastore'}
	  ->run_query( 'SELECT isolate_id FROM analysis_results WHERE name=?', 'hicap', { fetch => 'col_arrayref' } );
	my %hicap = map { $_ => 1 } @$hicap;
	my $hi    = $script->{'datastore'}->run_query(
		"SELECT id FROM $script->{'system'}->{'view'} WHERE species=?",
		'Haemophilus influenzae',
		{ fetch => 'col_arrayref' }
	);
	my %hi            = map { $_ => 1 } @$hi;
	my $filtered_list = [];
	foreach my $id (@$isolates) {
		next if defined $opts{'x'} && $id < $opts{'x'};
		next if defined $opts{'y'} && $id > $opts{'y'};
		push @$filtered_list, $id if $hi{$id} && !$hicap{$id};
	}
	return $filtered_list;
}

sub show_help {
	my $termios = POSIX::Termios->new;
	$termios->getattr;
	my $ospeed = $termios->getospeed;
	my $t      = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
	my ( $norm, $bold, $under ) = map { $t->Tputs( $_, 1 ) } qw(me md us);
	say << "HELP";
${bold}NAME$norm
    ${bold}hinfluenzae_hicap.pl$norm - Run capsule characterization script against 
    Haemophilus database

${bold}SYNOPSIS$norm
    ${bold}hinfluenzae_hicap.pl $norm [${under}options$norm]

${bold}OPTIONS$norm

${bold}--database$norm ${under}NAME$norm
    Database configuration name.
      
${bold}--help$norm
    This help page.
    
${bold}--isolates$norm ${under}LIST$norm  
    Comma-separated list of isolate ids to scan.
    
${bold}--limit$norm ${under}LIMIT$norm
    Stop after processing the number of isolates.
    
${bold}--quiet$norm
    Suppress output from underlying script.
    
${bold}--size$norm ${under}SIZE$norm
    Limit isolates to a sequence bin size >= value set (default 1,600,000 bp).
      
${bold}-x, --min$norm ${under}ID$norm
    Minimum isolate id.

${bold}-y, --max$norm ${under}ID$norm
    Maximum isolate id.
    
HELP
	return;
}
