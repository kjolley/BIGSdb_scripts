#!/usr/bin/env perl
#Update PubMLST Neisseria database with Ng azithromycin resistance predictions.
#Written by Keith Jolley
#Copyright (c) 2026, University of Oxford
#E-mail: keith.jolley@biology.ox.ac.uk
#Version: 20260312
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	HOST             => undef,                          #Use values in config.xml
	PORT             => undef,                          #But you can override here.
	USER             => undef,
	PASSWORD         => undef,
	SCRIPT           => 'neisseria_azi_resistance.py'
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use FindBin;
use lib "$FindBin::Bin/lib";
use JSON;
use IPC::Open3;
use Symbol       qw(gensym);
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
	'bin_dir=s'      => \$opts{'bin_dir'},
	'd|database=s'   => \$opts{'database'},
	'dir=s'          => \$opts{'dir'},
	'python=s'       => \$opts{'python'},
	'q|quiet'        => \$opts{'quiet'},
	'r|refresh'      => \$opts{'refresh'},
	'rRNA_23S=s'     => \$opts{'rRNA_23S'},
	'pro_NEIS1635=s' => \$opts{'pro_NEIS1635'},
	'NEIS1633=s'     => \$opts{'NEIS1633'},
	'tmp_dir=s'      => \$opts{'tmp_dir'}
) or die("Error in command line arguments\n");
if ( !$opts{'database'} ) {
	say "\nUsage: populate_neisseria_azi_resistance.pl --database <NAME> \n";
	exit;
}

$opts{'python'}  //= '/usr/bin/python3';
$opts{'bin_dir'} //= '/usr/local/bin/';
$opts{'bin_dir'} .= '/' if $opts{'bin_dir'} !~ /\/$/x;
$opts{'dir'} //= './';
$opts{'dir'} .= '/' if $opts{'dir'} !~ /\/$/x;
$opts{'tmp_dir'} //= '/var/tmp/';
$opts{'tmp_dir'} .= '/' if $opts{'tmp_dir'} !~ /\/$/x;

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
check_files();
main();
undef $script;

sub check_files {
	if ( !-e $opts{'bin_dir'} . SCRIPT ) {
		die $opts{'bin_dir'} . SCRIPT . " does not exist.\n";
	}
	if ( !-e $opts{'python'} ) {
		die "$opts{'python'} does not exist.\n";
	}
	if ( !-x $opts{'python'} ) {
		die "$opts{'python'} is not executable.\n";
	}
	foreach my $arg (qw (rRNA_23S pro_NEIS1635 NEIS1633)) {
		if ( !defined $opts{$arg} ) {
			die "--$arg not passed.\n";
		}
		my $file = "$opts{'dir'}$opts{$arg}";
		if ( !-e $file ) {
			die "$file does not exist.\n";
		}
		if ( !-r $file ) {
			die "$file exists but is not readable.\n";
		}
	}
}

sub main {
	my $qry =
		qq(SELECT id FROM $script->{'system'}->{'view'} WHERE species='Neisseria gonorrhoeae' AND id )
	  . q(IN (SELECT isolate_id FROM seqbin_stats WHERE total_length>1500000) );
	$qry .= q(AND id NOT IN (SELECT isolate_id FROM eav_text WHERE field='azithromycin'))
	  if !$opts{'refresh'};
	$qry .= q(ORDER BY id);
	my $ids = $script->{'datastore'}->run_query( $qry, undef, { fetch => 'col_arrayref' } );
	return if !@$ids;
	my $data = [];
	local $| = 1;
	my $count  = @$ids;
	my $plural = $count == 1 ? q() : q(s);
	print qq(Extracting allele data from database ($count isolate$plural) ...) if !$opts{'quiet'};

	foreach my $isolate_id (@$ids) {
		my $alleles = $script->{'datastore'}->run_query(
			q[SELECT locus,allele_id FROM allele_designations WHERE isolate_id=? AND locus IN ]
			  . q[('23S_rRNA','pro_NEIS1635','NEIS1633')],
			$isolate_id,
			{ fetch => 'all_arrayref', slice => {}, cache => 'get_alleles' }
		);
		next if @$alleles != 3;
		my %alleles;
		foreach my $allele (@$alleles) {
			$alleles{ $allele->{'locus'} } = $allele->{'allele_id'};
		}
		next if keys %alleles != 3;    #Just in case there are 2 designations for 1 locus and none for another.
		push @$data,
		  {
			id           => $isolate_id,
			'23S_rRNA'   => $alleles{'23S_rRNA'},
			NEIS1633     => $alleles{'NEIS1633'},
			pro_NEIS1635 => $alleles{'pro_NEIS1635'}
		  };

	}
	my $json         = encode_json($data);
	my $isolate_file = "$opts{'tmp_dir'}isolate_input_$$.json";
	open( my $fh, '>', $isolate_file ) or die "Cannot open $isolate_file for writing";
	say $fh $json;
	close $fh;
	say q(done. )                          if !$opts{'quiet'};
	print q(Running prediction script ...) if !$opts{'quiet'};

	my %params = (
		'--rRNA_23S'     => "$opts{'dir'}$opts{'rRNA_23S'}",
		'--pro_NEIS1635' => "$opts{'dir'}$opts{'pro_NEIS1635'}",
		'--NEIS1633'     => "$opts{'dir'}$opts{'NEIS1633'}",
		'--isolates'     => $isolate_file
	);
	my $err_fh = gensym;
	my $results_json;
	eval {
		my $pid = open3( my $in_fh, my $out_fh, $err_fh, $opts{'python'}, $opts{'bin_dir'} . SCRIPT, %params );
		close $in_fh;
		local $/ = undef;
		my $out = <$out_fh>;
		my $err = <$err_fh>;
		waitpid( $pid, 0 );
		my $exit_code = $? >> 8;
		if ($exit_code) {
			say "error!\n" if !$opts{'quiet'};
			say $err       if $err;
			say $out       if $out;
			exit(1);
		}
		$results_json = $out;
	};
	die "$@\n"  if $@;
	say 'done.' if !$opts{'quiet'};
	my $results;
	eval { $results = decode_json($results_json); };
	if ($@) {
		die "Invalid results file. $@\n";
	}
	print q(Updating database ...) if !$opts{'quiet'};
	eval {
		foreach my $result (@$results) {
			$script->{'db'}->do( 'INSERT INTO eav_text (isolate_id,field,value) VALUES (?,?,?)',
				undef, $result->{'id'}, 'azithromycin', $result->{'prediction'} );
		}
	};
	if ($@) {
		say 'failed!' if !$opts{'quiet'};
		$script->{'db'}->rollback;
		die "$@\n";
	}
	$script->{'db'}->commit;
	say 'done.' if !$opts{'quiet'};

	unlink $isolate_file;
}
