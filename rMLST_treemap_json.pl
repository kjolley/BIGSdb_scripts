#!/usr/bin/perl
#Generate treemap data from rMLST database
#Written by Keith Jolley, 2019
use strict;
use warnings;
use 5.010;
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	ISOLATE_DB       => 'pubmlst_rmlst_isolates',
};
#######End Local configuration###############################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use JSON;
use POSIX qw{strftime};
use constant RANKS => qw(genus family order class phylum);

#Direct all library logging calls to screen
my $log_conf =
    qq(log4perl.category.BIGSdb.Script        = INFO, Screen\n)
  . qq(log4perl.category.BIGSdb.Dataconnector = WARN, Screen\n)
  . qq(log4perl.category.BIGSdb.Datastore     = WARN, Screen\n)
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
		instance         => ISOLATE_DB,
		options          => { always_run => 1 }
	}
);
main();
undef $script;

sub main {
	my ($view) = @ARGV;
	$view //= 'public';
	my $qry =
	    q(select species,COALESCE(a1.value,'[undefined phylum]') AS phylum,)
	  . q(COALESCE(a2.value,'[undefined class]') AS class,)
	  . q(COALESCE(a3.value,'[undefined order]') AS "order",)
	  . q(COALESCE(a4.value,'[undefined family]')AS family,)
	  . qq(COALESCE(a5.value,'[undefined genus]') AS genus,COUNT(*) AS isolates FROM $view i FULL OUTER JOIN )
	  . q(isolate_value_extended_attributes a1 ON (a1.field_value,a1.attribute)=(i.species,'phylum') FULL OUTER JOIN )
	  . q(isolate_value_extended_attributes a2 ON (a2.field_value,a2.attribute)=(i.species,'class') FULL OUTER JOIN )
	  . q(isolate_value_extended_attributes a3 ON (a3.field_value,a3.attribute)=(i.species,'order') FULL OUTER JOIN )
	  . q(isolate_value_extended_attributes a4 ON (a4.field_value,a4.attribute)=(i.species,'family') FULL OUTER JOIN )
	  . q(isolate_value_extended_attributes a5 ON (a5.field_value,a5.attribute)=(i.species,'genus') )
	  . q(WHERE species IS NOT NULL GROUP BY species,phylum,class,"order",family,genus);
	my $taxonomy = $script->{'datastore'}->run_query( $qry, undef, { fetch => 'all_arrayref', slice => {} } );
	my $hierarchy = {};
	foreach my $taxon (@$taxonomy) {
		$hierarchy->{ $taxon->{'phylum'} }->{ $taxon->{'class'} }->{ $taxon->{'order'} }->{ $taxon->{'family'} }
		  ->{ $taxon->{'genus'} }->{ $taxon->{'species'} } = $taxon->{'isolates'};
	}
	my $data = [];
	foreach my $phylum ( keys %$hierarchy ) {
		my $class_list = [];
		foreach my $class ( keys %{ $hierarchy->{$phylum} } ) {
			my $order_list = [];
			foreach my $order ( keys %{ $hierarchy->{$phylum}->{$class} } ) {
				my $family_list = [];
				foreach my $family ( keys %{ $hierarchy->{$phylum}->{$class}->{$order} } ) {
					my $genus_list = [];
					foreach my $genus ( keys %{ $hierarchy->{$phylum}->{$class}->{$order}->{$family} } ) {
						my $species_list = [];
						foreach my $species ( keys %{ $hierarchy->{$phylum}->{$class}->{$order}->{$family}->{$genus} } )
						{
							push @$species_list,
							  {
								name => $species,
								size => $hierarchy->{$phylum}->{$class}->{$order}->{$family}->{$genus}->{$species}
							  };
						}
						push @$genus_list, { name => $genus, children => $species_list };
					}
					push @$family_list, { name => $family, children => $genus_list };
				}
				push @$order_list, { name => $order, children => $family_list };
			}
			push @$class_list, { name => $class, children => $order_list };
		}
		push @$data, { name => $phylum, children => $class_list };
	}
	say encode_json( { name => 'Bacteria', children => $data } );
}
