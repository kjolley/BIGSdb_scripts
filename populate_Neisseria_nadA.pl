#!/usr/bin/perl -T
#Populate NadA_peptide with confirmed '0' if currently undefined and
#nadA allele (NEIS1969) is defined with flag 'internal stop codon' or 
#'contains IS element'.
#Written by Keith Jolley 2016-2019
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	HOST             => undef,                  #Use values in config.xml
	PORT             => undef,                  #But you can override here.
	USER             => undef,
	PASSWORD         => undef
};
#######End Local configuration#############################################
use constant AUTOTAGGER    => -1;
use constant ALLELE_LOCUS  => 'NEIS1969';
use constant PEPTIDE_LOCUS => 'NadA_peptide';
use lib (LIB_DIR);
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
GetOptions( 'd|database=s' => \$opts{'database'}, ) or die("Error in command line arguments\n");
if ( !$opts{'database'} ) {
	say "\nUsage: populate_Neisseria_NadA.pl --database <NAME> \n";
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
		instance         => $opts{'database'},
		logger           => $logger
	}
);
die "This script can only be run against an isolate database.\n"
  if ( $script->{'system'}->{'dbtype'} // '' ) ne 'isolates';
my $allele_locus = $script->{'datastore'}->get_locus(ALLELE_LOCUS);
my $isolates     = $script->get_isolates_with_linked_seqs;
ISOLATES: foreach my $id (@$isolates) {
	my $alleles = $script->{'datastore'}->get_allele_designations( $id, ALLELE_LOCUS );
	next ISOLATES if !@$alleles;
	my $peptides = $script->{'datastore'}->get_allele_designations( $id, PEPTIDE_LOCUS );
	next ISOLATES if @$peptides;
	foreach my $allele (@$alleles) {
		next ISOLATES if $allele->{'allele_id'} eq '0';
		my $flags = $allele_locus->get_flags( $allele->{'allele_id'} );
		my %flags = map { $_ => 1 } @$flags;
		if ( $flags{'internal stop codon'} || $flags{'frameshift'} || $flags{'contains IS element'} ) {
			say "Marking id: $id NadA_peptide as missing.";
			eval {
				$script->{'db'}->do(
					'INSERT INTO allele_designations (isolate_id,locus,allele_id,'
					  . 'sender,status,method,curator,date_entered,datestamp) VALUES (?,?,?,?,?,?,?,?,?)',
					undef, $id, PEPTIDE_LOCUS, '0', AUTOTAGGER, 'confirmed', 'automatic', AUTOTAGGER, 'now', 'now'
				);
			};
			$script->{'db'}->do( 'INSERT INTO history (isolate_id,timestamp,action,curator) VALUES (?,?,?,?)',
				undef, $id, 'now',
				PEPTIDE_LOCUS . q(: new designation '0' (allele with frameshift, internal stop codon or IS element)), AUTOTAGGER );
			if ($@) {
				$script->{'db'}->rollback;
				say $@;
				exit(1);
			} else {
				$script->{'db'}->commit;
			}
			next ISOLATES;
		}
	}
}
