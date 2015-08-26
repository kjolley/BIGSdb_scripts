#!/usr/bin/perl -T
#Generate tab-delimited text file of allele lengths for a set of
#isolates
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
use Log::Log4perl qw(get_logger);
my %opts;
getopts( 'd:i:p:s:h', \%opts );

if ( $opts{'h'} ) {
	show_help();
	exit;
}
if ( !$opts{'d'} || !$opts{'s'} ) {
	say "Usage: isolate_allele_lengths.pl -d <isolate db config> [-i <isolate ids>] [-p <project list>] -s <scheme_id>";
	say "Help: isolate_allele_lengths.pl -h";
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
my $scheme_id  = $opts{'s'};
my $project_id = $opts{'p'};
my $loci       = $script->{'datastore'}->get_scheme_loci($scheme_id);
local $" = "\t";
say "id\tisolate\t@$loci";
my $isolates = $script->get_isolates_with_linked_seqs;
my %locus_cache;

foreach my $isolate_id (@$isolates) {
	my $isolate    = get_isolate_name($isolate_id);
	my $allele_ids = $script->{'datastore'}->get_all_allele_ids($isolate_id);
	print "$isolate_id\t$isolate";
	foreach my $locus (@$loci) {
		if ( !$locus_cache{$locus} ) {
			$locus_cache{$locus} = $script->{'datastore'}->get_locus($locus);
		}
		my $allele_id;
		my $seq;
		if ( ref $allele_ids->{$locus} eq 'ARRAY' && @{ $allele_ids->{$locus} } ) {
			$allele_id = $allele_ids->{$locus}->[0];
			$seq       = $locus_cache{$locus}->get_allele_sequence($allele_id);
		}
		if ($$seq) {
			print "\t" . ( length $$seq );
		} else {
			print "\t";
		}
	}
	print "\n";
}
get_isolate_name( undef, { finish => 1 } );                                        #Finish statement handle
undef $script;

sub get_isolate_name {
	my ( $isolate_id, $options ) = @_;
	$options = {} if ref $options ne 'HASH';
	state $sql = $script->{'db'}->prepare("SELECT isolate FROM isolates WHERE id=?");
	if ( $options->{'finish'} ) {
		undef $sql;
		return;
	}
	$sql->execute($isolate_id);
	my $isolate = $sql->fetchrow_array;
	return $isolate;
}

sub show_help {
	print << "HELP";

Usage isolate_allele_lengths.pl -d <isolate database config> -s <scheme id>

Options
-------
-d <name>  Isolate database configuration name.
-h         This help page.
-i         Comma-separated list of isolate ids.
-p <id>    Project id.
-s <id>    Scheme id.

HELP
	return;
}
