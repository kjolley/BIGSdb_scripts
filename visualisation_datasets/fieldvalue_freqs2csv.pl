#!/usr/bin/env perl
#Written by Keith Jolley
#Copyright (c) 2021, University of Oxford
#E-mail: keith.jolley@zoo.ox.ac.uk
#
#Outputs CSV file for visualisations.
#GPL3.
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	DBASE_CONFIG     => 'pubmlst_neisseria_isolates',
	PROJECT_ID       => 3
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Constants qw(LOG_TO_SCREEN);
use Getopt::Long qw(:config no_ignore_case);
my %opts;
GetOptions(
	'fields=s'   => \$opts{'fields'},
	'headings=s' => \$opts{'headings'}
);

#Direct all library logging calls to screen
my $log_conf = LOG_TO_SCREEN;
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my $script = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		instance         => DBASE_CONFIG,
		logger           => $logger,
		options          => { always_run => 1 },
	}
);
die "Script initialization failed.\n" if !defined $script->{'db'};
$opts{'fields'}   //= q(age_range);
$opts{'headings'} //= $opts{'fields'};
my $fields   = [ split /\s*,\s*/x, $opts{'fields'} ];
my $headings = [ split /\s*,\s*/x, $opts{'headings'} ];
check_fields_and_headings( $fields, $headings );
main();
undef $script;

sub check_fields_and_headings {
	my ( $fields, $headings ) = @_;
	if ( @$fields != @$headings ) {
		die "Unequal number of fields and headings.\n";
	}
	my $prov_fields = $script->{'xmlHandler'}->get_field_list;
	my %allowed_fields = map { $_ => 1 } ( @$prov_fields, 'clonal_complex' );
	foreach my $field (@$fields) {
		die "Invalid field '$field'.\n" if !$allowed_fields{$field};
	}
	return;
}

sub main {
	local $" = q(,);
	my $isolates = $script->{'datastore'}->run_query(
		"SELECT @$fields FROM isolates i LEFT JOIN "
		  . 'temp_isolates_scheme_fields_1 s ON i.id=s.id WHERE i.id IN '
		  . '(SELECT isolate_id FROM project_members WHERE project_id=?)',
		PROJECT_ID,
		{ fetch => 'all_arrayref', slice => {} }
	);
	my $dataset = {};
	my %combinations;
	foreach my $isolate (@$isolates) {
		my $values = rewrite_labels($isolate);
		local $" = q(|);
		$combinations{qq(@$values)}++;
	}
	say qq(@$headings,value);
	foreach my $combination ( keys %combinations ) {
		my @values = split /\|/x, $combination;
		local $" = q(,);
		say qq(@values,$combinations{$combination});
	}
	return;
}

sub rewrite_labels {
	my ($isolate) = @_;
	state $age_labels = get_age_labels();
	my $values = [];
	foreach my $field (@$fields) {
		my $value;
		if ( $field eq 'age_range' ) {
			$isolate->{$field} //= 'null';
			$value = $age_labels->{ $isolate->{$field} };
			$value = qq(@$value) if ref $value;
		} else {
			$value = $isolate->{$field};
			$value = qq(@$value) if ref $value;
		}
		push @$values, $value // 'none';
	}
	return $values;
}

sub get_age_labels {
	my $values = $script->{'xmlHandler'}->get_field_option_list('age_range');
	my $i      = 1;
	my $labels = {};
	foreach my $value (@$values) {
		$labels->{$value} = "$i: $value";
		$i++;
	}
	$labels->{'null'} = "$i: unknown";
	return $labels;
}
