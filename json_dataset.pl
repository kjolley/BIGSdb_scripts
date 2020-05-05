#!/usr/bin/env perl
#Extract JSON dataset for visualisation from BIGSdb database
#records in BIGSdb Neisseria database
#Written by Keith Jolley
#Copyright (c) 2020, University of Oxford
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
#Version: 20200501
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
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
	'database=s' => \$opts{'d'},
	'help'       => \$opts{'h'},
	'query=s'    => \$opts{'query'},
) or die("Error in command line arguments\n");
if ( $opts{'h'} ) {
	show_help();
	exit;
}
foreach my $required (qw(d query)) {
	if ( !$opts{$required} ) {
		show_help();
		exit;
	}
}
my $script = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		options          => \%opts,
		instance         => $opts{'d'},
		logger           => $logger
	}
);
die "Script initialization failed - server too busy?\n"
  if !defined $script->{'db'};
die "This script can only be run against an isolate database.\n"
  if ( $script->{'system'}->{'dbtype'} // '' ) ne 'isolates';
$opts{'fields'} //= '*';
$opts{'view'}   //= 'isolates';
main();
undef $script;

sub main {
	my $dataset = $script->{'datastore'}->run_query( $opts{'query'}, undef, { fetch => 'all_arrayref', slice => {} } )
	  ;
	say encode_json($dataset);
	return;
}

sub show_help {
	my $termios = POSIX::Termios->new;
	$termios->getattr;
	my $ospeed = $termios->getospeed;
	my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
	my ( $norm, $bold, $under ) = map { $t->Tputs( $_, 1 ) } qw(me md us);
	say << "HELP";
${bold}NAME$norm
    ${bold}json_dataset.pl$norm - Return JSON dataset for isolates defined
    in a project

${bold}SYNOPSIS$norm
    ${bold}json_dataset.pl --database [${under}DATABASE$norm${bold}] --query [${under}SQL$norm${bold}]$norm

${bold}OPTIONS$norm

${bold}--database$norm ${under}NAME$norm
    Database configuration name.
    
${bold}--query$norm ${under}SQL$norm
      
${bold}--help$norm
    This help page.
    
HELP
	return;
}
