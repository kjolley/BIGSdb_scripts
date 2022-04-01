#!/usr/bin/env perl
#Populate A. baumannii biofilm_production field based on OD values
#Written by Julia Moreno Manjon and Keith Jolley 2022
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	DBASE_CONFIG     => 'pubmlst_abaumannii_isolates',
	CURATOR_ID       => -2
};
#######End Local configuration#############################################
use lib (LIB_DIR);
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
my $script = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		instance         => DBASE_CONFIG,
		logger           => $logger
	}
);
main();

sub main {
	my $ids = get_isolate_ids();
	foreach my $id (@$ids) {
		my $OD_neg_control = $script->{'datastore'}
		  ->run_query( 'SELECT value FROM eav_float WHERE (isolate_id,field)=(?,?)', [ $id, 'OD_neg_control' ] );
		my $OD_neg_control_SD = $script->{'datastore'}
		  ->run_query( 'SELECT value FROM eav_float WHERE (isolate_id,field)=(?,?)', [ $id, 'OD_neg_control_SD' ] );
		my $OD_sample = $script->{'datastore'}
		  ->run_query( 'SELECT value FROM eav_float WHERE (isolate_id,field)=(?,?)', [ $id, 'OD_sample' ] );
		my $biofilm_production = biofilm_production( $OD_neg_control, $OD_neg_control_SD, $OD_sample );
		say "$id	$biofilm_production";
		eval {
			$script->{'db'}->do( 'INSERT INTO eav_text (isolate_id,field,value) VALUES (?,?,?)',
				undef, $id, 'biofilm_production', $biofilm_production );
		};
		if ($@) {
			say $@;
			$script->{'db'}->rollback;
			exit;
		}
		$script->{'db'}->commit;
	}
}

#Get isolate ids where biofilm_production field is NULL and OD values are set.
sub get_isolate_ids {
	return $script->{'datastore'}->run_query(
		"SELECT id FROM isolates WHERE id NOT IN (SELECT isolate_id FROM eav_text WHERE field=?) AND "
		  . 'id IN (SELECT isolate_id FROM eav_float WHERE field=?) AND id IN (SELECT isolate_id FROM eav_float WHERE '
		  . 'field=?) AND id IN (SELECT isolate_id FROM eav_float WHERE field=?)',
		[ 'biofilm_production', 'OD_neg_control_SD', 'OD_neg_control', 'OD_sample' ],
		{ fetch => 'col_arrayref' }
	);
}

sub biofilm_production {
	my ( $OD_neg_control, $OD_neg_control_SD, $OD_sample ) = @_;
	my $ODc = $OD_neg_control + ( 3 * $OD_neg_control_SD );
	if ( $OD_sample <= $ODc ) {
		return "none";
	}
	if ( $ODc < $OD_sample and $OD_sample <= 2 * $ODc ) {
		return "low";
	}
	if ( 2 * $ODc < $OD_sample and $OD_sample <= 4 * $ODc ) {
		return "medium";
	}
	if ( 4 * $ODc < $OD_sample ) {
		return "high";
	} else {
		return "null";
	}
}
