#!/usr/bin/perl -T
#Move loci between BIGSdb databases
#Written by Keith Jolley 2011-2012
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
use BIGSdb_Scripts::Migrate;

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
die "This script should only be called on seqdef databases.\n" if any {$_ eq 'isolates'} $script->get_db_types;
foreach my $locus(@loci){
	die "Locus $locus does not exist in database $opts{'a'}.\n" if !$script->locus_exists_in_source($locus);	
	die "Locus $locus already exists in database $opts{'b'}.\n" if $script->locus_exists_in_destination($locus);
	die "Locus $locus is a member of a scheme.\n" if $script->is_locus_in_scheme($locus);
	my $missing_users = $script->get_missing_allele_seq_users_in_destination($locus);
	my $missing_locus_curators = $script->_get_missing_curator_in_locus_tables($locus);
	local $" = "\n";
	if (keys %$missing_users){
		print "Missing users in destination database:\n";
		print "$missing_users->{$_}\n" foreach keys %$missing_users;
		exit;
	}
	if (keys %$missing_locus_curators){
		print "Missing locus curators in destination database:\n";
		print "$missing_locus_curators->{$_}\n" foreach keys %$missing_locus_curators;
		exit;
	}
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

