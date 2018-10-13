#!/usr/bin/env perl
#Written by Keith Jolley
#Copyright (c) 2018, University of Oxford
#E-mail: keith.jolley@zoo.ox.ac.uk
#
#Generates report for EMERT dataset.
#
#BIGSdb is free software: you can redistribute it and/or modify
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
#along with BIGSdb.  If not, see <http://www.gnu.org/licenses/>.
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	IMAGE_DIR        => '.'
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Utils;
use perlchartdir;
use Date::Manip;
my $db_config = 'pubmlst_neisseria_emert';

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
		instance         => $db_config,
		logger           => $logger,
		options          => { always_run => 1 },
	}
);
my $view = $script->{'system'}->{'view'};
die "Script initialization failed.\n" if !defined $script->{'db'};
main();
undef $script;

sub main {
	print_head();
	say q(<h1>EMERT statistics</h1>);
	my $datestamp = BIGSdb::Utils::get_datestamp();
	say qq(<h2>Report generated: $datestamp</h2>);
	say q(<div class="box resultstable">);
	print_totals();
	print_group_with_time();
	print_straintype_with_time();
	print_cc_with_time();
	say q(<div style="clear:both"></div>);
	say q(</html>);
	return;
}

sub print_head {
	say << 'HTML';
<!DOCTYPE html>
<head>
<meta name="viewport" content="width=device-width" />
<title>EMERT II - PubMLST.org</title>
<link rel="stylesheet" type="text/css" href="https://pubmlst.org/css/jquery-ui.css" media="Screen"/>
<link rel="stylesheet" type="text/css" href="https://pubmlst.org/css/fontawesome-all.css" media="Screen"/>
<link rel="stylesheet" type="text/css" href="https://pubmlst.org/css/bigsdb.css" media="Screen"/>
<script src="https://pubmlst.org/javascript/jquery.js?" type="text/Javascript"></script>
<script src="https://pubmlst.org/javascript/jquery.tablesorter.js" type="text/Javascript"></script>
<script src="https://pubmlst.org/javascript/jquery.slimbox2.js?" type="text/Javascript"></script>
<script type="text/Javascript">
$(function () {
	$("table#summary_table").tablesorter({widgets:['zebra']});
});
</script>
</head>	
<body>
<div class="content" style="min-height: 400px">
HTML
	return;
}

sub print_totals {
	say q(<h2>Summary totals</h2>);
	my $countries = $script->{'datastore'}
	  ->run_query( "SELECT DISTINCT country FROM $view ORDER BY country", undef, { fetch => 'col_arrayref' } );
	say q(<div class="scrollable">);
	say q(<table class="resultstable tablesorter" id="summary_table"><thead><tr><th rowspan="2">Country</th>)
	  . q(<th colspan="4" data-sorter="false">Isolates submitted</th></tr>);
	say q(<tr><th>past 6 months</th><th>past year</th><th>past 2 years</th><th>Total</th></tr></thead><tbody>);
	my $td = 1;
	foreach my $country (@$countries) {
		my @values;
		foreach my $time ( '6 months', '1 year', '2 years' ) {
			push @values,
			  $script->{'datastore'}->run_query(
				qq(SELECT COUNT(*) FROM $view WHERE LEAST(date_received,date_entered) > )
				  . qq((NOW()-interval '$time') AND country=?),
				$country
			  );
		}
		push @values, $script->{'datastore'}->run_query( qq(SELECT COUNT(*) FROM $view WHERE country=?), $country );
		local $" = q(</td><td>);
		say qq(<tr class="td$td"><td>$country</td><td>@values</td></tr>);
		$td = $td == 1 ? 2 : 1;
	}
	my @values;
	foreach my $time ( '6 months', '1 year', '2 years' ) {
		push @values,
		  $script->{'datastore'}->run_query(
			qq(SELECT COUNT(*) FROM $view WHERE LEAST(date_received,date_entered) > ) . qq((NOW()-interval '$time') ),
		  );
	}
	say q(</tbody><tfoot>);
	push @values, $script->{'datastore'}->run_query(qq(SELECT COUNT(*) FROM $view));
	local $" = q(</th><th>);
	say qq(<tr class="td$td tablesorter-ignoreRow"><th>Total</th><th>@values</th></tr>);
	say q(</tfoot></table>);
	say q(</div>);
	return;
}

sub print_group_with_time {
	say q(<h2>Capsule groups</h2>);

	#Cumulative totals
	my $dates  = get_dates( undef, { yearly => 1 } );
	my $days   = @$dates;
	my $values = get_groups($dates);
	process_dates( $days, $dates );
	my $file_name = 'capsule_cumulative.png';
	my $full_path = IMAGE_DIR . qq(/$file_name);
	print_chart( $dates, $values, $full_path );
	print_img_fieldset( 'Cumulative', 'Capsule groups (cumulative totals)', $file_name );

	#Past 2 years
	$days   = 365 * 2;
	$dates  = get_dates($days);
	$values = get_groups($dates);
	process_dates( $days, $dates );
	$file_name = 'capsule_2y.png';
	$full_path = IMAGE_DIR . qq(/$file_name);
	print_chart( $dates, $values, $full_path, 10 );
	print_img_fieldset( 'Past 2 years', 'Capsule groups (past 2 years)', $file_name );

	#Past year
	$days   = 365;
	$dates  = get_dates($days);
	$values = get_groups($dates);
	process_dates( $days, $dates );
	$file_name = 'capsule_1y.png';
	$full_path = IMAGE_DIR . qq(/$file_name);
	print_chart( $dates, $values, $full_path, 10 );
	print_img_fieldset( 'Past year', 'Capsule groups (past year)', $file_name );
	say q(<div style="clear:both"></div>);
	return;
}

sub get_groups {
	my ($dates) = @_;
	my $values = [];
	foreach my $date (@$dates) {
		my $data = $script->{'datastore'}->run_query(
			"SELECT capsule_group,count(*) AS count FROM $view "
			  . 'WHERE LEAST(date_received,date_entered) <= ? AND LEAST(date_received,date_entered) >=? '
			  . 'AND capsule_group IS NOT NULL GROUP BY '
			  . 'capsule_group',
			[ $date, $dates->[0] ],
			{ fetch => 'all_arrayref', slice => {}, cache => 'get_cc' }
		);
		my %groups = map { $_->{'capsule_group'} => $_->{'count'} } @$data;
		push @$values, \%groups;
	}
	return $values;
}

sub print_straintype_with_time {
	say q(<h2>Strain types</h2>);
	say q(<p>Only records with no missing data for ST, capsule group, PorA and FetA are included.</p>);

	#Cumulative totals
	my $dates  = get_dates( undef, { yearly => 1 } );
	my $days   = @$dates;
	my $values = get_strains($dates);
	process_dates( $days, $dates );
	my $file_name = 'strain_cumulative.png';
	my $full_path = IMAGE_DIR . qq(/$file_name);
	print_chart( $dates, $values, $full_path, 8 );
	print_img_fieldset( 'Cumulative', 'Strain types (cumulative totals)', $file_name );

	#Past 2 years
	$days   = 365 * 2;
	$dates  = get_dates($days);
	$values = get_strains($dates);
	process_dates( $days, $dates );
	$file_name = 'strain_2y.png';
	$full_path = IMAGE_DIR . qq(/$file_name);
	print_chart( $dates, $values, $full_path, 8 );
	print_img_fieldset( 'Past 2 years', 'Strain types (past 2 years)', $file_name );

	#Past year
	$days   = 365;
	$dates  = get_dates($days);
	$values = get_strains($dates);
	process_dates( $days, $dates );
	$file_name = 'strain_1y.png';
	$full_path = IMAGE_DIR . qq(/$file_name);
	print_chart( $dates, $values, $full_path, 8 );
	print_img_fieldset( 'Past year', 'Strain types (past year)', $file_name );
	say q(<div style="clear:both"></div>);
	return;
}

sub get_strains {
	my ($dates)     = @_;
	my $values      = [];
	my $prov_fields = $script->{'datastore'}
	  ->run_query( "SELECT id,capsule_group FROM $view", undef, { fetch => 'all_arrayref', slice => {} } );
	my %prov_data = map { $_->{'id'} => { capsule_group => $_->{'capsule_group'} } } @$prov_fields;
	foreach my $date (@$dates) {
		my $ids = $script->{'datastore'}->run_query(
			"SELECT id FROM $view WHERE LEAST(date_received,date_entered)<=? AND "
			  . 'LEAST(date_received,date_entered)>=?',
			[ $date, $dates->[0] ],
			{ fetch => 'col_arrayref' }
		);
		my $composites = {};
		foreach my $id (@$ids) {
			state %cache;
			if ( !$cache{$id} ) {
				my $value =
				  $script->{'datastore'}
				  ->get_composite_value( $id, 'strain_designation', $prov_data{$id}, { no_format => 1 } );
				$cache{$id} = $value;
			}
			next if $cache{$id} =~ /ND/x;
			$composites->{ $cache{$id} }++;
		}
		push @$values, $composites;
	}
	return $values;
}

sub print_cc_with_time {
	say q(<h2>Clonal complexes</h2>);

	#Cumulative totals
	my $dates  = get_dates( undef, { yearly => 1 } );
	my $days   = @$dates;
	my $values = get_ccs($dates);
	process_dates( $days, $dates );
	my $file_name = 'cc_cumulative.png';
	my $full_path = IMAGE_DIR . qq(/$file_name);
	print_chart( $dates, $values, $full_path, 10 );
	print_img_fieldset( 'Cumulative', 'Clonal complexes (cumulative totals)', $file_name );

	#Past 2 years
	$days   = 365 * 2;
	$dates  = get_dates($days);
	$values = get_ccs($dates);
	process_dates( $days, $dates );
	$file_name = 'cc_2y.png';
	$full_path = IMAGE_DIR . qq(/$file_name);
	print_chart( $dates, $values, $full_path, 10 );
	print_img_fieldset( 'Past 2 years', 'Clonal complexes (past 2 years)', $file_name );

	#Past year
	$days   = 365;
	$dates  = get_dates($days);
	$values = get_ccs($dates);
	process_dates( $days, $dates );
	$file_name = 'cc_1y.png';
	$full_path = IMAGE_DIR . qq(/$file_name);
	print_chart( $dates, $values, $full_path, 10 );
	print_img_fieldset( 'Past year', 'Clonal complexes (past year)', $file_name );
	say q(<div style="clear:both"></div>);
	return;
}

sub get_ccs {
	my ($dates) = @_;
	my $values = [];
	foreach my $date (@$dates) {
		my $data = $script->{'datastore'}->run_query(
			"SELECT clonal_complex,count(*) AS count FROM $view v JOIN temp_isolates_scheme_fields_1 t ON "
			  . 'v.id=t.id WHERE LEAST(date_received,date_entered) <= ? AND LEAST(date_received,date_entered) >=? AND clonal_complex IS NOT NULL GROUP BY '
			  . 'clonal_complex',
			[ $date, $dates->[0] ],
			{ fetch => 'all_arrayref', slice => {}, cache => 'get_cc' }
		);
		my %ccs = map { $_->{'clonal_complex'} => $_->{'count'} } @$data;
		push @$values, \%ccs;
	}
	return $values;
}

sub print_img_fieldset {
	my ( $title, $caption, $path ) = @_;
	say qq(<fieldset style="float:left"><legend>$title</legend>);
	say qq(<a href="$path" data-rel="lightbox-1" class="lightbox" title="$caption">)
	  . qq(<img src="$path" alt="$caption" style="width:200px;border:1px dashed black" />)
	  . q(</a></fieldset>);
	return;
}

sub get_dates {
	my ( $elapsed_days, $options ) = @_;
	my $mindate = defined $elapsed_days ? ParseDate("$elapsed_days days ago") : undef;
	my $today = ParseDate( BIGSdb::Utils::get_datestamp() );
	if ( !$mindate ) {
		$mindate = ParseDate($today);
		my $mindate2 = $script->{'datastore'}->run_query("SELECT LEAST(date_received,date_entered) FROM $view");
		my $date1    = $mindate;
		my $date2    = ParseDate($mindate2);
		my $flag     = Date_Cmp( $date1, $date2 );
		if ( $flag > 0 ) {
			$mindate = ParseDate($mindate2);
		}
	}
	$mindate = ParseDate('2007-01-01') if Date_Cmp( '2007-01-01', $mindate ) > 0;
	my $minyear = UnixDate( $mindate, '%Y' );
	my $dates   = [];
	my $days    = &Delta_Format( &DateCalc( $mindate, $today ), 0, '%0dds' );
	for my $i ( 0 .. $days ) {
		my $date = &UnixDate( &DateCalc( $mindate, "+ $i days" ), '%Y-%m-%d' );
		if ( $options->{'monthly'} ) {
			if ( $date =~ /-01$/x ) {
				push @$dates, $date;
			}
		} elsif ( $options->{'yearly'} ) {
			if ( $date =~ /01-01$/x ) {
				push @$dates, $date;
			}
		} else {
			push @$dates, $date;
		}
	}
	return $dates;
}

sub process_dates {
	my ( $days, $dates ) = @_;
	my $mindate = $dates->[0];
	if ( $days <= 365 ) {

		#label beginning of each month
		for my $i ( 0 .. @$dates - 1 ) {
			if ( $dates->[$i] !~ /-01$/x ) {
				$dates->[$i] = q();
			} else {
				$dates->[$i] =~ s/-01$//x;    #Remove day part
			}
		}
	} elsif ( $days <= 365 * 2 ) {

		#label beginning of each month
		for my $i ( 0 .. @$dates - 1 ) {
			if ( $dates->[$i] !~ /\d[02468]-01$/x ) {
				$dates->[$i] = q();
			} else {
				$dates->[$i] =~ s/-01$//x;    #Remove day part
			}
		}
	} else {

		#label beginning of each year
		for my $i ( 0 .. @$dates - 1 ) {
			if ( $dates->[$i] !~ /01-01$/x ) {
				$dates->[$i] = q();
			} else {
				$dates->[$i] =~ s/-01-01$//x;    #Remove month and day part
			}
		}

		#if mindate is more than two months from beginning of the next year, label it
		my $too_near;
		for my $i ( 0 .. 60 ) {
			if ( $dates->[$i] ne q() ) {
				$too_near = 1;
			}
		}
		if ( !$too_near ) {
			$dates->[0] = &UnixDate( $mindate, '%Y-%m-%d' );
		}
	}
	return;
}

sub print_chart {
	my ( $dates, $values, $file_name, $top_hits ) = @_;
	my $size     = 1.5;
	my $bgcolour = 0xffffff;
	my $c        = XYChart->new( 500 * $size, 500 * $size, $bgcolour, 0x0, 1 );
	$c->setPlotArea( 55 * $size, 80 * $size, 420 * $size, 345 * $size, 0xffffff );
	$c->addLegend( 55 * $size, 20 * $size, 0, '', 8 * $size )->setBackground($perlchartdir::Transparent);
	$c->xAxis()->setWidth( 2 * $size );
	$c->xAxis()->setTitle( 'Date received', '', 8 * $size );
	$c->xAxis()->setLabels($dates);
	$c->xAxis()->setLabelStyle( '', 6 * $size, 0x000000, 45 );
	$c->yAxis()->setWidth( 2 * $size );
	$c->yAxis()->setTitle( 'Cumulative frequency', '', 8 * $size );
	$c->yAxis()->setLabelStyle( '', 6 * $size, 0x000000 );
	my $layer = $c->addLineLayer();
	$layer->setLineWidth( 2 * $size );
	my @keys = sort { $values->[-1]->{$b} <=> $values->[-1]->{$a} } keys %{ $values->[-1] };
	my $count = 0;

	foreach my $key (@keys) {
		my $types = [];
		foreach my $value (@$values) {
			push @$types, $value->{$key};
		}
		$layer->addDataSet( $types, -1, $key );
		$count++;
		if ( $top_hits && $count == $top_hits ) {
			last;
		}
	}
	$c->makeChart($file_name);
	return;
}
