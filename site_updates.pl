#!/usr/bin/env perl
#Written by Keith Jolley
#Copyright (c) 2018, University of Oxford
#E-mail: keith.jolley@zoo.ox.ac.uk
#
#Generates HTML report of recent updates on PubMLST.
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
use Getopt::Long qw(:config no_ignore_case);
use REST::Client;
use JSON;
use Date::Manip;
use Parallel::ForkManager;
use constant BASE_URI => 'http://rest.pubmlst.org';

#Term::Cap and POSIX are just used for formatting help page
use Term::Cap;
use POSIX;
my %opts;
GetOptions(
	'check=i'   => \$opts{'check'},
	'help'      => \$opts{'help'},
	'show=i'    => \$opts{'show'},
	'threads=i' => \$opts{'threads'}
) or die("Error in command line arguments\n");
if ( $opts{'help'} ) {
	show_help();
	exit;
}
$opts{'check'} //= 5;
$opts{'check'} = 1 if $opts{'check'} < 1;
$opts{'show'} //= $opts{'check'};
my $client = REST::Client->new();
main();

sub main {
	my $last_updates = get_updates();
	my $dates        = get_dates();
	my $buffer;
	my $shown = 0;
	foreach my $date (@$dates) {
		my @date_buffer;
		foreach my $species ( sort keys %$last_updates ) {
			my @species_buffer;
			if ( $last_updates->{$species}->{'sequences'}->{$date} ) {
				my $url =
				    qq(/bigsdb?db=$last_updates->{$species}->{'seqdef_db'}&amp;page=tableQuery&amp;table=sequences&amp;)
				  . qq(s1=datestamp&amp;y1==&amp;t1=$date&amp;submit=1);
				my $plural = $last_updates->{$species}->{'sequences'}->{$date} == 1 ? q() : q(s);
				push @species_buffer,
				  qq(<a href="$url">$last_updates->{$species}->{'sequences'}->{$date} allele$plural</a>);
			}
			if ( keys %{ $last_updates->{$species}->{'schemes'} } ) {
				foreach my $scheme ( keys %{ $last_updates->{$species}->{'schemes'} } ) {
					if ( $last_updates->{$species}->{'schemes'}->{$scheme}->{'dates'}->{$date} ) {
						my $url = qq(/bigsdb?db=$last_updates->{$species}->{'seqdef_db'}&amp;page=query&amp;)
						  . qq(scheme_id=$scheme&amp;s1=datestamp&amp;y1==&amp;t1=$date&amp;submit=1);
						my $plural =
						  $last_updates->{$species}->{'schemes'}->{$scheme}->{'dates'}->{$date} == 1 ? q() : q(s);
						push @species_buffer,
						  qq(<a href="$url">$last_updates->{$species}->{'schemes'}->{$scheme}->{'dates'}->{$date} )
						  . qq($last_updates->{$species}->{'schemes'}->{$scheme}->{'description'} profile$plural</a>);
					}
				}
			}
			if ( $last_updates->{$species}->{'isolates'}->{$date} ) {
				my $url = qq(/bigsdb?db=$last_updates->{$species}->{'isolate_db'}&amp;page=query&amp;)
				  . qq(prov_field1=datestamp&amp;prov_operator1==&amp;prov_value1=$date&amp;submit=1);
				my $plural = $last_updates->{$species}->{'isolates'}->{$date} == 1 ? q() : q(s);
				push @species_buffer,
				  qq(<a href="$url">$last_updates->{$species}->{'isolates'}->{$date} isolate$plural</a>);
			}
			local $" = q(</li><li>);
			push @date_buffer, qq($species:<ul><li>@species_buffer</li></ul>) if @species_buffer;
		}
		if (@date_buffer) {
			$buffer .= qq(<p><b>$date</b>);
			local $" = qq(</li><li>\n);
			$buffer .= qq(<ul><li>@date_buffer</li></ul>\n);
			$shown++;
			last if $shown == $opts{'show'};
		}
	}
	say $buffer;
	return;
}

sub get_dates {
	my $dates = [];
	for my $i ( 0 .. $opts{'check'} - 1 ) {
		push @$dates, UnixDate( DateCalc( &datestamp, "- $i days" ), '%Y-%m-%d' );
	}
	return $dates;
}

sub get_updates {
	my $all_updates = {};
	$client->request( 'GET', BASE_URI );
	my $resources = from_json( $client->responseContent );
	my $pm        = Parallel::ForkManager->new( $opts{'threads'} );
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data ) = @_;
			my $name    = $data->{'name'};
			my $updates = $data->{'updates'};
			return if !$updates->{'isolates'} && !$updates->{'sequences'} && !$updates->{'schemes'};
			foreach my $key ( keys %{$updates} ) {
				$all_updates->{$name}->{$key} = $updates->{$key};
			}
		}
	);
	my @databases;
	foreach my $resource (@$resources) {
		next if !$resource->{'databases'};
		foreach my $db ( @{ $resource->{'databases'} } ) {
			push @databases, $db;
		}
	}
	foreach my $db (@databases) {
		$pm->start and next;
		my $name;
		my $updates = {};
		$client->request( 'GET', $db->{'href'} );
		if ( $db->{'description'} =~ /(.+)\ sequence\/profile\ definitions/x ) {
			$name = $1;
			my $seqdef = from_json( $client->responseContent );
			if ( $seqdef->{'sequences'} ) {
				my $sequences = get_sequence_updates( $seqdef->{'sequences'} );
				$updates->{'sequences'} = $sequences if keys %$sequences;
			}
			if ( $seqdef->{'schemes'} ) {
				my $schemes = get_scheme_updates( $seqdef->{'schemes'} );
				$updates->{'schemes'} = $schemes if keys %$schemes;
			}
			if ( keys %$updates && $db->{'href'} =~ /db\/(.+)$/x ) {
				$updates->{'seqdef_db'} = $1;
			}
		} elsif ( $db->{'description'} =~ /(.+)\ isolates/x ) {
			$name = $1;
			my $isolate_db = from_json( $client->responseContent );
			if ( $isolate_db->{'isolates'} ) {
				my $isolates = get_isolate_updates( $isolate_db->{'isolates'} );
				$updates->{'isolates'} = $isolates if keys %$isolates;
			}
			if ( keys %$updates && $db->{'href'} =~ /db\/(.+)$/x ) {
				$updates->{'isolate_db'} = $1;
			}
		}
		$pm->finish( 0, { name => $name, updates => $updates } );
	}
	$pm->wait_all_children;
	return $all_updates;
}

sub get_sequence_updates {
	my ($url) = @_;
	my $data  = {};
	my $dates = get_dates();
	$client->request( 'GET', $url );
	my $sequences = from_json( $client->responseContent );
	if ( $sequences->{'last_updated'} && $sequences->{'last_updated'} ge $dates->[-1] ) {
		foreach my $date (@$dates) {
			$client->request( 'GET', "$url?updated_on=$date" );
			$sequences = from_json( $client->responseContent );
			$data->{$date} = $sequences->{'records'} if $sequences->{'records'};
		}
	}
	return $data;
}

sub get_scheme_updates {
	my ($url) = @_;
	my $data  = {};
	my $dates = get_dates();
	$client->request( 'GET', $url );
	my $schemes = from_json( $client->responseContent );
	foreach my $scheme ( @{ $schemes->{'schemes'} } ) {
		$client->request( 'GET', $scheme->{'scheme'} );
		my $this_scheme = from_json( $client->responseContent );
		if ( $this_scheme->{'last_updated'} && $this_scheme->{'last_updated'} ge $dates->[-1] ) {
			foreach my $date (@$dates) {
				$client->request( 'GET', "$scheme->{'scheme'}?updated_on=$date" );
				my $profiles = from_json( $client->responseContent );
				$data->{ $this_scheme->{'id'} }->{'dates'}->{$date} = $profiles->{'records'} if $profiles->{'records'};
			}
			$data->{ $this_scheme->{'id'} }->{'description'} = $this_scheme->{'description'};
		}
	}
	return $data;
}

sub get_isolate_updates {
	my ($url) = @_;
	my $data  = {};
	my $dates = get_dates();
	$client->request( 'GET', $url );
	my $isolates = from_json( $client->responseContent );
	if ( $isolates->{'last_updated'} && $isolates->{'last_updated'} ge $dates->[-1] ) {
		foreach my $date (@$dates) {
			$client->request( 'GET', "$url?updated_on=$date" );
			$isolates = from_json( $client->responseContent );
			$data->{$date} = $isolates->{'records'} if $isolates->{'records'};
		}
	}
	return $data;
}

sub datestamp {
	my @date = localtime;
	my $year = 1900 + $date[5];
	my $mon  = $date[4] + 1;
	my $day  = $date[3];
	return ( sprintf( '%d-%02d-%02d', $year, $mon, $day ) );
}

sub show_help {
	my $termios = POSIX::Termios->new;
	$termios->getattr;
	my $ospeed = $termios->getospeed;
	my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
	my ( $norm, $bold, $under ) = map { $t->Tputs( $_, 1 ) } qw/me md us/;
	say << "HELP";
${bold}NAME$norm
${bold}site_updates.pl$norm - Generate HTML report of recent updates on PubMLST

${bold}SYNOPSIS$norm
    ${bold}site_updates.pl$norm [${under}options$norm]

${bold}OPTIONS$norm

${bold}--check$norm ${under}DAYS$norm
    Number of days to check (default 5).

${bold}--help$norm
    This help page.

${bold}--show$norm ${under}DAYS$norm
    Number of days to show (default = days to check). 
    
${bold}--threads$norm [${under}THREADS$norm]
    Threads to use when querying API. Do not set too high or you may overload
    the remote server (and get banned if you are running this as a third 
    party).

HELP
	return;
}
