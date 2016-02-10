#!/usr/bin/perl 
#Written by Keith Jolley
#Copyright (c) 2016, University of Oxford
#E-mail: keith.jolley@zoo.ox.ac.uk
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with BIGSdb.  If not, see <http://www.gnu.org/licenses/>.
#
#Read API token from ./enterobase_token
use strict;
use warnings;
use 5.010;
use POSIX;
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
use LWP::UserAgent;
use MIME::Base64;
use JSON;
use Data::Dumper qw(Dumper);
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	UPDATE_USER      => 3
};
my $SERVER_ADDRESS = 'https://enterobase.warwick.ac.uk/api/v1.0';
my $TOKEN_FILE     = "$ENV{'HOME'}/.enterobase_token";
#######End Local configuration################################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
my %opts;
GetOptions(
	'd|database=s'      => \$opts{'d'},
	'e|enterobase_db=s' => \$opts{'e'},
	'l|loci'            => \$opts{'l'},
	'h|help'            => \$opts{'h'},
	'r|route=s'         => \$opts{'r'},
	's|scheme=s'        => \$opts{'s'},
	'u|update'          => \$opts{'u'}
) or die("Error in command line arguments\n");

if ( $opts{'h'} ) {
	show_help();
	exit;
}
my $ua           = LWP::UserAgent->new;
my $token        = get_api_token();
my $base64string = encode_base64("$token:");
$ua->default_header( Authorization => "Basic $base64string" );
main();

sub main {
	my %methods = (
		l => sub {
			my $loci = get_loci();
			local $" = qq(\n);
			say qq(@$loci);
		},
		r => sub {
			get_route();
		},
		u => sub {
			update();
		}
	);
	foreach my $param ( sort keys %opts ) {
		if ( $opts{$param} && $methods{$param} ) {
			$methods{$param}->();
			return;
		}
	}
	get_route();
	return;
}

sub initiate_script_object {
	my $script = BIGSdb::Offline::Script->new(
		{
			config_dir       => CONFIG_DIR,
			lib_dir          => LIB_DIR,
			dbase_config_dir => DBASE_CONFIG_DIR,
			options          => \%opts,
			instance         => $opts{'d'},
		}
	);
	die "Fatal error - check log file.\n" if !$script->{'db'};
	return $script;
}

sub update {
	check_options(qw(d e s));
	my $script = initiate_script_object;
	my $loci   = get_loci();
#	my $locus_count = 0;
	foreach my $locus (@$loci) {
		my $existing_alleles =
		  $script->{'datastore'}->run_query( "SELECT allele_id,sequence FROM sequences WHERE locus=?",
			$locus, { fetch => 'all_arrayref', cache => 'get_all_alleles' } );
		my %existing = map { $_->[0] => $_->[1] } @$existing_alleles;
		my $url      = "$SERVER_ADDRESS/$opts{'e'}/$opts{'s'}/alleles?locus=$locus&limit=1000";
		say $url;
		while (1){
			my $resp     = $ua->get($url);
			if ( $resp->is_success ) {
				my $data        = decode_json( $resp->decoded_content );
				my $alleles = $data->{'alleles'};
				foreach my $allele (@$alleles){
					say ">$locus-$allele->{'allele_id'}";
					say $allele->{'seq'};
				}
				if ($data->{'paging'}->{'next'}){
					$url = $data->{'paging'}->{'next'};
					say Dumper $data->{'paging'};
				} else {
					last;
				}
			} else {
				say $resp->status_line;
				say $resp->decoded_content;
				last;
			}
		}
#		$locus_count++;
#		last if $locus_count == 10;
		last;
	}
	return;
}

sub get_route {
	my $url = $SERVER_ADDRESS;

	$url .= "/$opts{'e'}" if $opts{'e'};
	$url .= "/$opts{'r'}" if $opts{'r'};
	say $url;
	my $resp = $ua->get($url);
	if ( $resp->is_success ) {
		say $resp->status_line;
		my $data = decode_json( $resp->decoded_content );
		say Dumper $data;
	} else {
		say $resp->status_line;
		say $resp->decoded_content;
	}
	return;
}

sub check_options {
	my @options = @_;
	my %option_names = ( d => 'Database configuration', e => 'Enterobase database', s => 'scheme' );
	foreach my $option (@options) {
		die "$option_names{$option} (-$option) not provided.\n" if !$opts{$option};
	}
	my %enterobase_dbs = map { $_ => 1 } qw(senterica ecoli mcatarrhalis yersinia);
	if ( $opts{'e'} && !$enterobase_dbs{ $opts{'e'} } ) {
		die "$opts{'e'} is not a recognized Enterobase database.\n";
	}
	return;
}

sub get_loci {
	die "No scheme selected.\n" if !$opts{'s'};
	check_options(qw(e s));
	my $url  = "$SERVER_ADDRESS/$opts{'e'}/$opts{'s'}/loci?limit=10000";
	my $resp = $ua->get($url);
	if ( $resp->is_success ) {
		my $data        = decode_json( $resp->decoded_content );
		my $loci        = $data->{'loci'};
		my $locus_names = [];
		foreach my $locus (@$loci) {
			push @$locus_names, $locus->{'locus'};
		}
		return $locus_names;
	} else {
		say $resp->status_line;
		say $resp->decoded_content;
	}
	return [];
}

sub get_api_token {
	if ( !-e $TOKEN_FILE ) {
		die "$TOKEN_FILE does not exist.\n";
	}
	open( my $fh, '<:raw', $TOKEN_FILE ) || die "Cannot open $TOKEN_FILE for reading.\n";
	my $contents = do { local $/ = undef; <$fh> };
	close $fh;
	return $contents;
}

sub show_help {
	my $termios = POSIX::Termios->new;
	$termios->getattr;
	my $ospeed = $termios->getospeed;
	my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
	my ( $norm, $bold, $under ) = map { $t->Tputs( $_, 1 ) } qw/me md us/;
	say << "HELP";
${bold}NAME$norm
    ${bold}get_defs_from_enterobase.pl$norm - Download definitions from Enterobase

${bold}SYNOPSIS$norm
    ${bold}get_defs_from_enterobase.pl [${under}options$norm]

${bold}OPTIONS$norm

${bold}-d, --database$norm ${under}NAME$norm
    Database configuration name.

${bold}-e, --enterobase_db$norm ${under}DBASE$norm
    Enterobase database name.

${bold}-h, --help$norm.
    This help page.
    
${bold}-l, --loci$norm
    Retrieve list of loci.
    
${bold}-r, --route$norm
    Relative route.
    
${bold}-s, --scheme$norm ${under}SCHEME NAME$norm
    Scheme name.
HELP
	return;
}
