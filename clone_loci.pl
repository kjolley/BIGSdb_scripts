#!/usr/bin/perl -T
#Copy locus configurations between BIGSdb isolate databases
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
	PASSWORD         => undef
};
#######End Local configuration###############################

use lib (LIB_DIR);
use List::MoreUtils qw(any);
use Getopt::Std;
use BIGSdb_Scripts::Migrate;

my %opts;
getopts( 'a:b:l:h', \%opts );

if ($opts{'h'}){
	show_help();
	exit;
}

if (any {!$opts{$_}} qw (a b l)){
	print "\nUsage: clone_loci.pl -a <source database config> -b <destination database config> -l <loci>\n\n";
	print "Help: clone_loci.pl -h\n";
	exit;
}

my $script = BIGSdb_Scripts::Migrate->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		host             => HOST,
		port             => PORT,
		user             => USER,
		password         => PASSWORD,
		options			 => \%opts,
		instance		 => $opts{'a'},
		writable         => 1
	}
);

my @loci = split /,/, $opts{'l'};

#Do some checks before actually doing anything
die "This script should only be called on isolate databases.\n" if any {$_ ne 'isolates'} $script->get_db_types;
foreach my $locus(@loci){
	die "Locus $locus does not exist in database $opts{'a'}.\n" if !$script->locus_exists_in_source($locus);	
	die "Locus $locus already exists in database $opts{'b'}.\n" if $script->locus_exists_in_destination($locus);
}

foreach my $locus(@loci){
	$script->clone_locus($locus);
}

sub show_help {
	print << "HELP";

Usage clone_loci.pl -a <source database config> -b <destination database config>

Options
-------
-a <name>  Source database configuration name.
-b <name>  Destination database configuration name.
-h         This help page.
-l <list>  Loci - comma-separated list of loci to migrate.

HELP
	return;
}