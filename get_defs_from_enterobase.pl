#!/usr/bin/perl 
#Written by Keith Jolley
#Copyright (c) 2016-2017, University of Oxford
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
use Digest::MD5;
use Data::Dumper qw(Dumper);
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	UPDATE_USER      => 3
};
my $SERVER_ADDRESS = 'https://enterobase.warwick.ac.uk/api/v2.0';
my $TOKEN_FILE     = "$ENV{'HOME'}/.enterobase_token";
#######End Local configuration################################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Utils;
my %opts;
GetOptions(
	'a|update_alleles'  => \$opts{'a'},
	'check_alleles'     => \$opts{'check_alleles'},
	'd|database=s'      => \$opts{'d'},
	'e|enterobase_db=s' => \$opts{'e'},
	'l|loci'            => \$opts{'l'},
	'limit=i'           => \$opts{'limit'},
	'locus_regex=s'     => \$opts{'locus_regex'},
	'h|help'            => \$opts{'h'},
	'p|update_profiles' => \$opts{'p'},
	'r|route=s'         => \$opts{'r'},
	'n|new_loci'        => \$opts{'n'},
	's|scheme=s'        => \$opts{'s'},
	'scheme_id=i'       => \$opts{'scheme_id'}
) or die("Error in command line arguments\n");

if ( $opts{'h'} ) {
	show_help();
	exit;
}
my $ua           = LWP::UserAgent->new;
my $token        = get_api_token();
my $base64string = encode_base64("$token:");
$ua->default_header( Authorization => "Basic $base64string" );
my $script;
if ( $opts{'d'} ) {
	$script = initiate_script_object();
}
local $| = 1;
main();
undef $script;

sub main {
	my %methods = (
		check_alleles => sub { check_alleles(); },
		a             => sub {
			update_alleles();
		},
		l => sub {
			my $loci = get_loci();
			local $" = qq(\n);
			say qq(@$loci);
		},
		r => sub {
			get_route();
		},
		p => sub {
			update_profiles();
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
	my $script_object = BIGSdb::Offline::Script->new(
		{
			config_dir       => CONFIG_DIR,
			lib_dir          => LIB_DIR,
			dbase_config_dir => DBASE_CONFIG_DIR,
			options          => \%opts,
			instance         => $opts{'d'},
		}
	);
	die "Fatal error - check log file.\n" if !$script_object->{'db'};
	return $script_object;
}

sub check_alleles {
	check_options(qw(e s));
	my $loci  = get_loci();
	my $first = 1;
	my $infinity = 9**99;
  LOCUS: foreach my $locus (@$loci) {
		my ( $allele_count, $complete_cds, $min_length, $max_length ) = ( 0, 0, $infinity, 0 );
		next LOCUS if $opts{'locus_regex'} && $locus !~ /$opts{'locus_regex'}/x;
		my $url = "$SERVER_ADDRESS/$opts{'e'}/$opts{'s'}/alleles?locus=$locus";
		$url .= "&limit=$opts{'limit'}" if $opts{'limit'};
		my %already_received;
	  PAGE: while (1) {
			my $resp = $ua->get($url);
			if ( !$resp->is_success ) {
				say $url;
				say $resp->status_line;
				say $resp->decoded_content;
				last PAGE;
			}
			my $data    = decode_json( $resp->decoded_content );
			my $alleles = $data->{'alleles'};
			foreach my $allele (@$alleles) {
				next if $already_received{ $allele->{'allele_id'} };
				$already_received{ $allele->{'allele_id'} } = 1;
				$allele_count++;
				my $cds = BIGSdb::Utils::is_complete_cds( $allele->{'seq'} );
				$complete_cds++ if $cds->{'cds'};
				my $length = length $allele->{'seq'};
				$min_length = $length if $length < $min_length;
				$max_length = $length if $length > $max_length;
			}
			if ( @$alleles && $data->{'paging'}->{'next'} ) {
				$url = $data->{'paging'}->{'next'};
			} else {
				last PAGE;
			}
		}
		say qq(Locus\tAlleles\tCDS\t%CDS\tMin length\tMax length) if $first;
		$first = 0;
		print qq($locus\t$allele_count);
		if ( $allele_count > 0 ) {
			my $pc_cds = BIGSdb::Utils::decimal_place( $complete_cds / $allele_count * 100, 1 );
			print qq(\t$complete_cds\t$pc_cds\t$min_length\t$max_length);
		}
		print qq(\n);
	}
	return;
}

#BIGSdb doesn't allows locus names to begin with an underscore or a digit.
sub get_mapped_loci {
	return {
		'93282_00335'  => '_93282_00335',
		'93282_01680'  => '_93282_01680',
		'93282_01824'  => '_93282_01824',
		'93282_01825'  => '_93282_01825',
		'93282_02692'  => '_93282_02692',
		'93282_02869'  => '_93282_02869',
		'93282_02880'  => '_93282_02880',
		'93282_02906'  => '_93282_02906',
		'93282_02907'  => '_93282_02907',
		'93282_03486'  => '_93282_03486',
		'93282_03487'  => '_93282_03487',
		'93282_03488'  => '_93282_03488',
		'EP-D108_gp14' => 'EP_D108_gp14',
		'EP-D108_gp15' => 'EP_D108_gp15',
		'EP-D108_gp19' => 'EP_D108_gp19',
		'EP-D108_gp31' => 'EP_D108_gp31',
		'EP-D108_gp38' => 'EP_D108_gp38',
		'EP-D108_gp41' => 'EP_D108_gp41',
		'EP-D108_gp42' => 'EP_D108_gp42',
		'EP-D108_gp57' => 'EP_D108_gp57',
		'O96-HN_00746' => 'O96_HN_00746',
		'O96-HN_00747' => 'O96_HN_00747',
		'O96-HN_00748' => 'O96_HN_00748',
		'O96-HN_01860' => 'O96_HN_01860',
		'O96-HN_01861' => 'O96_HN_01861',
		'O96-HN_02062' => 'O96_HN_02062',
		'O96-HN_02063' => 'O96_HN_02063',
		'O96-HN_02064' => 'O96_HN_02064',
		'O96-HN_03653' => 'O96_HN_03653',
		'O96-HN_03657' => 'O96_HN_03657',
	};
}

sub update_alleles {
	check_options(qw(d e s));
	my $loci        = get_loci();
	my $mapped_loci = get_mapped_loci();
  LOCUS: foreach my $locus (@$loci) {
		next LOCUS if $opts{'locus_regex'} && $locus !~ /$opts{'locus_regex'}/x;
		my $locus_name = $mapped_loci->{$locus} // $locus;
		my $existing_alleles =
		  $script->{'datastore'}->run_query( 'SELECT allele_id,sequence FROM sequences WHERE locus=?',
			$locus_name, { fetch => 'all_arrayref', cache => 'get_all_alleles' } );
		next if @$existing_alleles && $opts{'n'};
		my %existing = map { $_->[0] => $_->[1] } @$existing_alleles;
		my %existing_seqs = map { Digest::MD5::md5_hex( $_->[1] ) => $_->[0] } @$existing_alleles;
		my $url = "$SERVER_ADDRESS/$opts{'e'}/$opts{'s'}/alleles?locus=$locus";
		$url .= "&limit=$opts{'limit'}" if $opts{'limit'};
		my %already_received;
	  PAGE: while (1) {
			my $resp = $ua->get($url);
			if ( !$resp->is_success ) {
				say $url;
				say $resp->status_line;
				say $resp->decoded_content;
				last PAGE;
			}
			my $data    = decode_json( $resp->decoded_content );
			my $alleles = $data->{'alleles'};
			foreach my $allele (@$alleles) {

				#Don't trust Enterobase API not to serve up bad alleles.
				if ( $allele->{'allele_id'} < 1 ) {
					say "$locus-$allele->{'allele_id'} has a -ve allele id!";
					next;
				}
				if ( $already_received{ $allele->{'allele_id'} } ) {
					say "$locus-$allele->{'allele_id'} has already been received in this download!";
					next;
				}
				( my $new_seq = uc( $allele->{'seq'} ) ) =~ s/\s//gx;
				if ( $new_seq =~ /N/x ) {
					say "$locus-$allele->{'allele_id'} contains Ns!";
					next;
				}
				$already_received{ $allele->{'allele_id'} } = 1;
				if ( $existing{ $allele->{'allele_id'} } ) {
					next if $existing{ $allele->{'allele_id'} } eq $new_seq;
					say "$locus-$allele->{'allele_id'} has changed!";
				} elsif ( defined $existing_seqs{ Digest::MD5::md5_hex( $allele->{'seq'} ) } ) {
					say "$locus-$allele->{'allele_id'} is already defined as "
					  . "$locus-$existing_seqs{Digest::MD5::md5_hex($allele->{'seq'})}!";
				} else {
					eval {
						say "Inserting $locus-$allele->{'allele_id'}";
						$script->{'db'}->do(
							'INSERT INTO sequences (locus,allele_id,sequence,status,date_entered,'
							  . 'datestamp,sender,curator) VALUES (?,?,?,?,?,?,?,?)',
							undef,
							$locus_name,
							$allele->{'allele_id'},
							$new_seq,
							'unchecked',
							'now',
							'now',
							UPDATE_USER,
							UPDATE_USER
						);
					};
					if ($@) {
						say $@;
						$script->{'db'}->rollback;
						last LOCUS;
					}
					$existing_seqs{ Digest::MD5::md5_hex( $allele->{'seq'} ) } = $allele->{'allele_id'};
				}
			}
			if ( @$alleles && $data->{'paging'}->{'next'} ) {
				$url = $data->{'paging'}->{'next'};
			} else {
				last PAGE;
			}
		}
		$script->{'db'}->commit;
	}
	return;
}

sub update_profiles {
	check_options(qw(d e s scheme_id));
	my $loci = get_loci();
	local $" = q(,);
	my $existing_alleles = {};
	foreach my $locus (@$loci) {
		my $allele_ids = $script->{'datastore'}->run_query( 'SELECT allele_id FROM sequences WHERE locus=?',
			$locus, { fetch => 'col_arrayref', cache => 'get_all_allele_ids' } );
		%{ $existing_alleles->{$locus} } = map { $_ => 1 } @$allele_ids;
	}
	my $existing_profiles =
	  $script->{'datastore'}
	  ->run_query( qq(SELECT ST,profile FROM mv_scheme_$opts{'scheme_id'}), undef, { fetch => 'all_arrayref' } );
	my $positions = $script->{'datastore'}->get_scheme_locus_indices( $opts{'scheme_id'} );
	my $existing  = {};
	foreach my $st_profile (@$existing_profiles) {
		my ( $st, $profile ) = @$st_profile;
		foreach my $locus (@$loci) {
			$existing->{$st}->{$locus} = $profile->[ $positions->{$locus} ];
		}
	}
	my $url = "$SERVER_ADDRESS/$opts{'e'}/$opts{'s'}/sts?show_alleles=true";
	$url .= "&limit=$opts{'limit'}" if $opts{'limit'};
	my %already_received;
  PAGE: while (1) {
		my $resp = $ua->get($url);
		if ( !$resp->is_success ) {
			say $url;
			say $resp->status_line;
			say $resp->decoded_content;
			last PAGE;
		}
		my $data = decode_json( $resp->decoded_content );
		
		my $profiles = $data->{'STs'};
	  PROFILE: foreach my $profile (@$profiles) {
			my $st = $profile->{'ST_id'};
			if ( $already_received{$st} ) {
				say "ST-$st has already been received in this download!";
				next;
			}
			$already_received{$st} = 1;
			my %alleles;
			my $alleles = $profile->{'alleles'};
			foreach my $allele (@$alleles) {
				$alleles{ $allele->{'locus'} } = $allele->{'allele_id'};
			}
			next PROFILE if keys %alleles != @$loci;
			if ( $existing->{$st} ) {
				my $changed;
				foreach my $locus (@$loci) {

					#					say "$locus $alleles{$locus} $existing->{$st}->{$locus}";
					$changed = 1 if $alleles{$locus} ne $existing->{$st}->{$locus};
				}
				if ($changed) {
					say "ST-$st has changed!";
					exit;
				}
				next PROFILE;
			} else {
				foreach my $locus (@$loci) {
					if ( !$existing_alleles->{$locus}->{ $alleles{$locus} } ) {
						say "Cannot insert ST-$st - $locus-$alleles{$locus} is not defined!";
						next PROFILE;
					}
				}
				say "Inserting ST-$st.";
				eval {
					$script->{'db'}->do(
						'INSERT INTO profiles (scheme_id,profile_id,sender,curator,'
						  . 'date_entered,datestamp) VALUES (?,?,?,?,?,?)',
						undef, $opts{'scheme_id'}, $st, UPDATE_USER, UPDATE_USER, 'now', 'now'
					);
					$script->{'db'}->do(
						'INSERT INTO profile_fields (scheme_id,profile_id,scheme_field,value,curator,'
						  . 'datestamp) VALUES (?,?,?,?,?,?)',
						undef, $opts{'scheme_id'}, $st, 'ST', $st, UPDATE_USER, 'now'
					);
					foreach my $locus (@$loci) {
						$script->{'db'}->do(
							'INSERT INTO profile_members (scheme_id,locus,profile_id,allele_id,'
							  . 'curator,datestamp) VALUES (?,?,?,?,?,?)',
							undef, $opts{'scheme_id'}, $locus, $st, $alleles{$locus}, UPDATE_USER, 'now'
						);
					}
				};
				if ($@) {
					say $@;
					$script->{'db'}->rollback;
					last PAGE;
				}
			}
		}
		if ( @$profiles && $data->{'paging'}->{'next'} ) {
			$url = $data->{'paging'}->{'next'};
		} else {
			last PAGE;
		}
	}
	$script->{'db'}->commit;
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
	my %option_names =
	  ( d => 'Database configuration', e => 'Enterobase database', s => 'Scheme', scheme_id => 'Scheme id' );
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
	my $url  = "$SERVER_ADDRESS/$opts{'e'}/$opts{'s'}/loci?limit=50000";
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
		say $url;
		say $resp->status_line;
		say $resp->decoded_content;
		exit;
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

${bold}-a, --update_alleles$norm
    Update allele sequences.
    
${bold}--check_alleles$norm
    Report count of alleles and number that are complete coding sequences.
    This can be used to check if there are problems with allele definitions
    in cgMLST schemes as most alleles should be complete CDS.
    Cannot be used with --update_alleles.
    
${bold}-d, --database$norm ${under}NAME$norm
    Database configuration name.

${bold}-e, --enterobase_db$norm ${under}DBASE$norm
    Enterobase database name.

${bold}-h, --help$norm.
    This help page.
    
${bold}-l, --loci$norm
    Retrieve list of loci.
    
${bold}--limit$norm ${under}LIMIT$norm
    Request LIMIT number of alleles or profiles per page. 
    
${bold}--locus_regex$norm ${under}REGEX$norm
    Restrict updated loci to those that match regular expression. 
    
${bold}-n, --new_loci$norm
    Only download new loci (those without any existing alleles defined).
    
${bold}-p, --update_profiles$norm
    Update allelic profiles.
    
${bold}-r, --route$norm
    Relative route.
    
${bold}-s, --scheme$norm ${under}SCHEME NAME$norm
    Scheme name.
    
${bold}--scheme_id$norm ${under}SCHEME ID$norm
    Scheme id number in BIGSdb database.
HELP
	return;
}
