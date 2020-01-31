#!/usr/bin/env perl
#Assemble short read data for BIGSdb records with an ENA_run_accession
#and no contigs in sequence bin
#Written by Keith Jolley
#Copyright (c) 2020, University of Oxford
#E-mail: keith.jolley@zoo.ox.ac.uk
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with BIGSdb.  If not, see <http://www.gnu.org/licenses/>.
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	HOST             => undef,                                      #Use values in config.xml
	PORT             => undef,                                      #But you can override here.
	USER             => undef,
	PASSWORD         => undef,
	SHOVILL_PATH     => '/home/linuxbrew/.linuxbrew/bin/shovill',
	BREW_PATH        => '/home/linuxbrew/.linuxbrew/bin',
	TMP_DIR          => '/var/tmp',
	ASSEMBLER_USER   => -2
};
#######End Local configuration#############################################
use constant ENA_URI => 'https://www.ebi.ac.uk/ena';
use constant PASSIVE => 1;
use LWP::UserAgent;
use Net::FTP;
use File::Path qw(make_path remove_tree);
use Digest::MD5::File qw(file_md5_hex);
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
use POSIX;
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Constants qw(LOG_TO_SCREEN);

#Direct all library logging calls to screen
my $log_conf = LOG_TO_SCREEN;
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my %opts;
GetOptions(
	'd|database=s'         => \$opts{'database'},
	'debug'                => \$opts{'debug'},
	'help'                 => \$opts{'help'},
	'i|isolates=s'         => \$opts{'i'},
	'isolate_list_file=s'  => \$opts{'isolate_list_file'},
	'I|exclude_isolates=s' => \$opts{'I'},
	'overwrite'            => \$opts{'overwrite'},
	'min_length'           => \$opts{'min_length'},
	'x|min=i'              => \$opts{'x'},
	'y|max=i'              => \$opts{'y'},
	'p|projects=s'         => \$opts{'p'},
	'P|exclude_projects=s' => \$opts{'P'},
	'quiet'                => \$opts{'quiet'},
	'ram=i'                => \$opts{'ram'},
	'tmp_dir=s'            => \$opts{'tmp_dir'},
	'user_id=i'            => \$opts{'user_id'}
) or die("Error in command line arguments\n");
if ( $opts{'help'} ) {
	show_help();
	exit;
}
if ( !$opts{'database'} ) {
	say "\nUsage: assemble_and_upload.pl -d <database configuration>\n";
	say 'Help: assemble_and_upload.pl -h';
	exit;
}
$ENV{'PATH'} .= ':' . BREW_PATH;
$opts{'tmp_dir'} //= TMP_DIR;
$opts{'tmp_dir'} =~ s/\/$//x;    #Remove trailing slash
$opts{'min_length'} //= 100;
$opts{'ram'}        //= 32;
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
exit;

sub main {
	return if required_fields_not_exist();
	return if assembler_user_not_exist();
	my $isolates = $script->get_isolates;
	$isolates = $script->filter_and_sort_isolates($isolates);
	$isolates = filter_out_those_with_seqbin($isolates);
	$isolates = filter_out_those_without_run_accession($isolates);
	foreach my $id ( sort { $a <=> $b } @$isolates ) {
		say qq(Processing id-$id...) if !$opts{'quiet'};
		my $accession = $script->{'datastore'}->run_query( 'SELECT ENA_run_accession FROM isolates WHERE id=?', $id );
		if ( invalid_accession($accession) ) {
			say qq(id-$id has an invalid run accession: $accession);
			next;
		}
		my $fastq_uris = get_fastq_uris($accession);
		eval {
			download_fastqs( $accession, $fastq_uris );
			run_shovill( $accession, $fastq_uris );
		};
		say $@ if $@;
	}
}

sub required_fields_not_exist {
	my $fields = $script->{'xmlHandler'}->get_field_list;
	my %fields = map { $_ => 1 } @$fields;
	foreach my $field (qw(ENA_run_accession assembly_details)) {
		if ( !$fields{$field} ) {
			say qq($field field does not exist.);
			return 1;
		}
	}
	return;
}

sub assembler_user_not_exist {
	my $user_id = $opts{'user_id'} // ASSEMBLER_USER;
	if ( !$script->{'datastore'}
		->run_query( 'SELECT EXISTS(SELECT * FROM users WHERE (id,user_name)=(?,?))', [ $user_id, 'assembler' ] ) )
	{
		say qq(Assembler user not defined - should be id:$user_id; user_name:assembler);
		return 1;
	}
	return;
}

sub filter_out_those_with_seqbin {
	my ($isolates) = @_;
	my $ids_with_seqbin =
	  $script->{'datastore'}->run_query( 'SELECT isolate_id FROM seqbin_stats', undef, { fetch => 'col_arrayref' } );
	my %ids_with_seqbin = map { $_ => 1 } @$ids_with_seqbin;
	my $list = [];
	foreach my $id (@$isolates) {
		push @$list, $id if !$ids_with_seqbin{$id};
	}
	return $list;
}

sub filter_out_those_without_run_accession {
	my ($isolates) = @_;
	my $ids_with_accession =
	  $script->{'datastore'}
	  ->run_query( 'SELECT id FROM isolates WHERE ena_run_accession IS NOT NULL', undef, { fetch => 'col_arrayref' } );
	my %ids_with_accession = map { $_ => 1 } @$ids_with_accession;
	my $list = [];
	foreach my $id (@$isolates) {
		push @$list, $id if $ids_with_accession{$id};
	}
	return $list;
}

sub invalid_accession {
	my ($accession) = @_;
	return 1 if $accession !~ /^(E|D|S)RR[0-9]{6,}$/x;
	return;
}

sub get_fastq_uris {
	my ($accession) = @_;
	my $uri = ENA_URI . qq(/data/warehouse/filereport?accession=$accession&result=read_run&fields=fastq_ftp,fastq_md5);
	my $ua  = LWP::UserAgent->new;
	my $response = $ua->get($uri);
	my $uris     = [];
	my $md5      = [];
	if ( $response->is_success ) {
		my $content = $response->decoded_content;
		my $pairs   = 0;
		foreach my $line ( split /\n/x, $content ) {
			chomp $line;
			next if !$line;
			next if $line =~ /^fastq_ftp/x;
			if ( $line =~ /\t/x ) {
				my ( $uri_string, $md5_string ) = split /\t/x, $line;
				if ( $uri_string =~ /;/x && $md5_string =~ /;/x ) {
					my @uris = split /;/x, $uri_string;
					my @md5s = split /;/x, $md5_string;
					foreach my $i ( 0 .. @uris - 1 ) {
						push @$uris,
						  {
							uri => $uris[$i],
							md5 => $md5s[$i]
						  };
					}
					$pairs++;
				}
			}
		}
		if ( !$pairs ) {
			die qq(No FASTQ files associated with this accession.\n);
		} elsif ( $pairs > 1 ) {
			die qq(Multiple paired FASTQ files associated with this accession.\n);
		}
		if ( @$uris != 2 ) {
			die scalar @$uris . q( links to FASTQ files. There should be 2.);
		}
		return $uris;
	} else {
		die $response->status_line;
	}
}

sub download_fastqs {
	my ( $accession, $fastq_uris ) = @_;
	my $tmp_dir = qq($opts{'tmp_dir'}/$accession);
	if ( $opts{'overwrite'} ) {
		remove_tree($tmp_dir);
	}
	if ( -d $tmp_dir ) {
		die qq(Temp directory $tmp_dir already exists. Aborting.\n);
	}
	make_path($tmp_dir);
	foreach my $fastq_uri (@$fastq_uris) {
		my ( $host, $path, $filename );
		if ( $fastq_uri->{'uri'} =~ /^([^\/]+)\/(.*)\/([^\/]+)$/x ) {
			( $host, $path, $filename ) = ( $1, $2, $3 );
		} else {
			die qq(Cannot extract host and path from URI: $fastq_uri->{'uri'}\n);
		}
		my $ftp = Net::FTP->new( $host, Passive => PASSIVE );
		$ftp->login;
		$ftp->cwd($path);
		$ftp->binary;
		$ftp->get( $filename, qq($tmp_dir/$filename) ) || die $ftp->message;
		$ftp->quit;
		my $md5_hash = file_md5_hex("$tmp_dir/$filename");

		if ( $md5_hash ne $fastq_uri->{'md5'} ) {
			die qq(MD5 hash of downloaded file $filename does not match archive MD5. Aborting.\n);
		}
		say qq(Downloaded $filename to $tmp_dir. MD5 matches.) if !$opts{'quiet'};
		$fastq_uri->{'filename'} = qq($tmp_dir/$filename);
	}
	return;
}

sub run_shovill {
	my ( $accession, $fastq_uris ) = @_;
	my $out_dir = "$opts{'tmp_dir'}/$accession/shovill";
	my @params  = (
		'--outdir' => $out_dir,
		'--R1'     => $fastq_uris->[0]->{'filename'},
		'--R2'     => $fastq_uris->[1]->{'filename'},
		'--minlen' => $opts{'min_length'},
		'--ram'    => $opts{'ram'}
	);
	say q(Running Shovill...) if !$opts{'quiet'};
	my $quiet = $opts{'debug'} ? q() : q(> /dev/null 2>&1);
	system( SHOVILL_PATH . qq( @params $quiet) );
	my $assembly_details = {};
	if ( -e "$out_dir/shovill.log" ) {
		open( my $fh, '<', "$out_dir/shovill.log" ) || die "Cannot open shovill.log.\n";
		while ( my $line = <$fh> ) {
			$assembly_details->{'shovill_parameters'} = 'default';
			if ( $line =~ /\[shovill\]\ This\ is\ shovill\ (.*)$/x ) {
				$assembly_details->{'shovill_version'} = $1;
			}
			if ( $line =~ /\[shovill\]\ Using\ kmers:\ (.*)$/x ) {
				$assembly_details->{'kmers'} = $1;
			}
			if ( $line =~ /^\[spades\]\s+SPAdes\ version:\ (.*)$/x ) {
				$assembly_details->{'spades_version'} = $1;
			}
		}
		close $fh;
	} else {
		die "No shovill.log file.\n";
	}
	use Data::Dumper;
	say Dumper $assembly_details;
}

sub show_help {
	my $termios = POSIX::Termios->new;
	$termios->getattr;
	my $ospeed = $termios->getospeed;
	my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
	my ( $norm, $bold, $under ) = map { $t->Tputs( $_, 1 ) } qw/me md us/;
	say << "HELP";
${bold}NAME$norm
    ${bold}assemble_and_upload.pl$norm - Assemble short read data for BIGSdb 
    records with an ENA_run_accession and no contigs in sequence bin.
    
${bold}SYNOPSIS$norm
    ${bold}assemble_and_upload.pl --database ${under}NAME$norm [${under}options$norm]

${bold}OPTIONS$norm
    
${bold}--database$norm ${under}CONFIG$norm
    Database configuration name
    
${bold}--debug$norm
    Output progress messages from Shovill.

${bold}--help$norm
    This help page.
    
${bold}-i, --isolates$norm ${under}LIST$norm  
    Comma-separated list of isolate ids to scan (ignored if -p used).
    
${bold}--isolate_list_file$norm ${under}FILE$norm  
    File containing list of isolate ids (ignored if -i or -p used).
           
${bold}-I, --exclude_isolates$norm ${under}LIST$norm
    Comma-separated list of isolate ids to ignore.
    
${bold}--min_length$norm ${under}LENGTH$norm
    Minimum length of contig
    
${bold}--overwrite$norm
    Overwrite download directory.
    
${bold}-p, --projects$norm ${under}LIST$norm
    Comma-separated list of project isolates to scan.

${bold}-P, --exclude_projects$norm ${under}LIST$norm
    Comma-separated list of projects whose isolates will be excluded.
    
${bold}--quiet$norm
    Suppress progress messages.
    
${bold}--ram$norm ${under}RAM (GB)$norm
    Amount of RAM to use (GB). Default is 32GB.
    
${bold}--tmp_dir$norm ${under}DIR$norm
    Temporary directory.
    
${bold}--user_id$norm ${under}ID$norm
    Assembler user id.
    
${bold}-x, --min$norm ${under}ID$norm
    Minimum isolate id.

${bold}-y, --max$norm ${under}ID$norm
    Maximum isolate id.
    
HELP
	return;
}
