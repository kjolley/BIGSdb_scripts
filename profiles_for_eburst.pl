#!/usr/bin/perl -T
#Generate tab-delimited text file of MLST profiles including a
#a column for ST
#Written by Keith Jolley 2014-2015
use strict;
use warnings;
use 5.010;
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	HOST             => 'zoo-oban',
	PORT             => 5432,
	USER             => 'apache',
	PASSWORD         => undef
};
#######End Local configuration###############################
use lib (LIB_DIR);
use List::MoreUtils qw(any);
use Getopt::Std;
use BIGSdb::Utils;
use BIGSdb::Offline::Script;
my %opts;
getopts( 'd:s:h', \%opts );

if ( $opts{'h'} ) {
	show_help();
	exit;
}
if ( !$opts{'d'} || !$opts{'s'}) {
	say "Usage: profiles_for_eburst.pl -d <isolate db config> -s <scheme id>";
	say "Help: profiles_for_eburst.pl -h";
	exit;
}
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
		instance         => $opts{'d'},
	}
);
main();
undef $script;

sub main {
	my ($self) = @_;
	my $scheme_id = BIGSdb::Utils::is_int( $opts{'s'} ) ? $opts{'s'} : 1;
	print_header($scheme_id);
	print_scheme($scheme_id);
	return;
}

sub print_header {
	my ($scheme_id) = @_;
	my $scheme_info = $script->{'datastore'}->get_scheme_info( $scheme_id, { get_pk => 1 } );
	if (!$scheme_info){
		$script->{'logger'}->error("Scheme does not exist.");
		undef $script;
		exit;
	}
	if ( !defined $scheme_info->{'primary_key'} ) {
		$script->{'logger'}->error("$scheme_info->{'description'}: No primary key set.");
		undef $script;
		exit;
	}
	my $loci = $script->{'datastore'}->get_scheme_loci($scheme_id);
	if ( $opts{'p'} ) {
		s/^$opts{'p'}// foreach @$loci;
	}
	local $" = "\t";
	say "$scheme_info->{'primary_key'}\t@$loci";
	return;
}

sub print_scheme {
	my ($scheme_id) = @_;
	my $scheme_info   = $script->{'datastore'}->get_scheme_info( $scheme_id, { get_pk => 1 } );
	my $loci          = $script->{'datastore'}->get_scheme_loci($scheme_id);
	my $scheme_fields = $script->{'datastore'}->get_scheme_fields($scheme_id);
	my $ids = $script->{'datastore'}->run_query("SELECT id FROM $script->{'system'}->{'view'} ORDER BY id", undef, {fetch=>'col_arrayref'});
	foreach my $id (@$ids){
		my $field_values = $script->{'datastore'}->get_scheme_field_values_by_isolate_id( $id, $scheme_id );
		foreach my $pk (keys %{$field_values->{lc $scheme_info->{'primary_key'}}} ){
			my $profile = $script->{'datastore'}->get_profile_by_primary_key($scheme_id, $pk);
			if ($profile){
				local $" = "\t";
				say "$pk\t@$profile";	
			}
			
		}
	}
	return;
}

sub show_help {
	print << "HELP";

Usage profiles_for_eburst.pl -a <isolate database config> -s <scheme id>

Options
-------
-a <name>  Isolate database configuration name.
-h         This help page.
-s <id>    Scheme id.

HELP
	return;
}
