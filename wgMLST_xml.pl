#!/usr/bin/perl -T
#Generate XML file of loci for wgMLST
#Written by Keith Jolley 2014
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
	PASSWORD         => ''
};
#######End Local configuration###############################
use lib (LIB_DIR);
use Getopt::Long;
use BIGSdb::Offline::Script;
use List::MoreUtils qw(all);
use constant DOMAIN => 'pubmlst.org';
my %opts;
GetOptions(
	'd|database=s'     => \$opts{'d'},
	'l|loci=s'         => \$opts{'l'},
	'L|exclude_loci=s' => \$opts{'L'},
	'n|name=s'         => \$opts{'n'},
	'p|species=s'      => \$opts{'p'},
	'R|locus_regex=s'  => \$opts{'R'},
	's|schemes=s'      => \$opts{'s'},
	'h|help'           => \$opts{'h'},
) or die("Error in command line arguments\n");

if ( $opts{'h'} ) {
	show_help();
	exit;
}
if ( !( all { $opts{$_} } qw(d n p) ) ) {
	say "Usage: wgMLST_xml.pl -d <seqdef db config> -n <scheme name> -p <species/genus>";
	say "Help: profiles_with_species.pl -h";
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
die "No connection to database (check logs).\n" if !defined $script->{'db'};
die "This script can only be run against a seqdef database.\n" if ( $script->{'system'}->{'dbtype'} // '' ) ne 'sequences';
my $loci = $script->get_selected_loci;
my $locus_desc_sql = $script->{'db'}->prepare("SELECT * FROM locus_descriptions WHERE locus=?");
say "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>";
say "<data>";
say "  <species>$opts{'p'}";
say "    <$opts{'n'}>";
say "      <loci>";

foreach my $locus (@$loci) {
	say "        <locus>$locus";
	say "          <url>http://" . DOMAIN . "/cgi-bin/bigsdb/bigsdb.pl?db=$opts{'d'}&amp;page=downloadAlleles&amp;locus=$locus</url>";
	$locus_desc_sql->execute($locus);
	my $desc =  $locus_desc_sql->fetchrow_hashref;
	if ($desc){
		foreach (qw (full_name product)){
			if (defined $desc->{$_}){ 
				chomp $desc->{$_};
				$desc->{$_} =~ s/&/\&amp;/g;
				say "          <$_>$desc->{$_}</$_>";
			}
		}
	}
	say "        </locus>";
}
say "      </loci>";
say "    </$opts{'n'}>";
say "  </species>";
say "</data>";

sub show_help {
	print << "HELP";

Usage profiles_with_species.pl -a <seqdef database config> -s <scheme id>

Options
-------
--d <name>                Seqdef database configuration name.
--database	

-h                        This help page.
--help

-l <list>                 Comma-separated list of loci to scan (ignored if -s
--loci <list>             used).

-L <list>                 Comma-separated list of loci to exclude
--exclude_loci <list>

-n <name>                 Name of scheme to use in XML
--name <name>

-p <species>              Species/genus name
--species <name>

-R <regex>                Regex for locus names
--locus_regex <regex>

-s <list>                 Comma-separated list of scheme loci to scan.
--schemes <list>

HELP
	return;
}
