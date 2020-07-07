#!/usr/bin/env perl
#Move loci between BIGSdb databases
#Written by Keith Jolley 2011-2020
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
use FindBin;
use lib "$FindBin::Bin/lib";
use List::MoreUtils qw(any);
use Getopt::Std;
use Migrate;

my %opts;
getopts( 'a:b:c:l:h', \%opts );

if ($opts{'h'}){
	show_help();
	exit;
}

if (any {!$opts{$_}} qw (a b l)){
	print "\nUsage: migrate_loci.pl -a <source database config> -b <destination database config> -l <loci>\n\n";
	print "Help: migrate_loci.pl -h\n";
	exit;
}
$opts{'throw_busy_exception'} = 1;
my $script = Migrate->new(
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
die "This script should only be called on seqdef databases.\n" if any {$_ eq 'isolates'} $script->get_db_types;
foreach my $locus(@loci){
	die "Locus $locus does not exist in database $opts{'a'}.\n" if !$script->locus_exists_in_source($locus);	
	die "Locus $locus already exists in database $opts{'b'}.\n" if $script->locus_exists_in_destination($locus);
	die "Locus $locus is a member of a scheme.\n" if $script->is_locus_in_scheme($locus);
}

foreach my $locus(@loci){
	$script->copy_locus($locus);
	$script->copy_alleles($locus);
	$script->update_locus_fields_in_clients($locus);
	$script->delete_locus_from_source($locus);
}

sub show_help {
	print << "HELP";

Usage migrate_loci.pl -a <source database config> -b <destination database config> [-c <client databases>]

Options
-------
-a <name>  Source database configuration name.
-b <name>  Destination database configuration name.
-c <list>  Client databases - comma-separated list of database configuration
           names of isolate databases that refer to loci.
-h         This help page.
-l <list>  Loci - comma-separated list of loci to migrate.

HELP
	return;
}

