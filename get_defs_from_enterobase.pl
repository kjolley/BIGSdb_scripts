#!/usr/bin/perl 
#Written by Keith Jolley
#Copyright (c) 2016-2018, University of Oxford
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
use Time::HiRes qw(usleep);
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
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
	'commit'            => \$opts{'commit'},
	'd|database=s'      => \$opts{'d'},
	'e|enterobase_db=s' => \$opts{'e'},
	'l|loci'            => \$opts{'l'},
	'limit=i'           => \$opts{'limit'},
	'locus_regex=s'     => \$opts{'locus_regex'},
	'h|help'            => \$opts{'h'},
	'p|update_profiles' => \$opts{'p'},
	'r|route=s'         => \$opts{'r'},
	'n|new_loci'        => \$opts{'n'},
	'm|allow_missing'   => \$opts{'allow_missing'},
	'no_errors'         => \$opts{'no_errors'},
	'reldate=i'         => \$opts{'reldate'},
	's|scheme=s'        => \$opts{'s'},
	'scheme_id=i'       => \$opts{'scheme_id'},
	'user_id=i'         => \$opts{'user_id'}
) or die("Error in command line arguments\n");

if ( $opts{'h'} ) {
	show_help();
	exit;
}
my $ua = LWP::UserAgent->new;
$ua->timeout(600);
my $token        = get_api_token();
my $base64string = encode_base64("$token:");
$ua->default_header( Authorization => "Basic $base64string" );
my $script;
if ( $opts{'d'} ) {
	$script = initiate_script_object();
}
local $| = 1;
$opts{'user_id'} //= -10;
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
	my $loci     = get_loci();
	my $first    = 1;
	my $infinity = 9**99;
  LOCUS: foreach my $locus (@$loci) {
		my ( $allele_count, $complete_cds, $min_length, $max_length ) = ( 0, 0, $infinity, 0 );
		next LOCUS if $opts{'locus_regex'} && $locus !~ /$opts{'locus_regex'}/x;
		my $url = "$SERVER_ADDRESS/$opts{'e'}/$opts{'s'}/alleles?locus=$locus";
		$url .= "&limit=$opts{'limit'}"     if $opts{'limit'};
		$url .= "&reldate=$opts{'reldate'}" if $opts{'reldate'};
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
			usleep(500_000);    #Rate-limiting
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

sub get_mapped_loci {
	if ( $opts{'e'} eq 'yersinia' ) {
		return {
			adk  => 'Yersinia_pseudotuberculosis_MLST_adk',
			argA => 'Yersinia_pseudotuberculosis_MLST_argA',
			aroA => 'Yersinia_pseudotuberculosis_MLST_aroA',
			glnA => 'Yersinia_pseudotuberculosis_MLST_glnA',
			thrA => 'Yersinia_pseudotuberculosis_MLST_thrA',
			tmk  => 'Yersinia_pseudotuberculosis_MLST_tmk',
			trpE => 'Yersinia_pseudotuberculosis_MLST_trpE'
		};
	} elsif ( $opts{'e'} eq 'mcatarrhalis' ) {
		return {
			abcZ    => 'Moraxella_catarrhalis_MLST_abcZ',
			adk     => 'Moraxella_catarrhalis_MLST_adk',
			efp     => 'Moraxella_catarrhalis_MLST_efp',
			fumC    => 'Moraxella_catarrhalis_MLST_fumC',
			glyBeta => 'Moraxella_catarrhalis_MLST_glyBeta',
			mutY    => 'Moraxella_catarrhalis_MLST_mutY',
			ppa     => 'Moraxella_catarrhalis_MLST_ppa',
			trpE    => 'Moraxella_catarrhalis_MLST_trpE'
		};
	}
	return {};
}

sub update_alleles {
	check_options(qw(d e s));
	my $loci        = get_loci();
	my $mapped_loci = get_mapped_loci();
	usleep(500_000);    #Rate-limiting
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
		$url .= "&limit=$opts{'limit'}"     if $opts{'limit'};
		$url .= "&reldate=$opts{'reldate'}" if $opts{'reldate'};
		my %already_received;
		local $| = 1;
	  PAGE: while (1) {
			usleep(500_000);    #Rate-limiting
			my $resp;
		  ATTEMPT: for my $attempt ( 1 .. 10 ) {
				$resp = $ua->get($url);
				if ( !$resp->is_success ) {
					if ( $attempt < 10 ) {
						say $url;
						say $resp->status_line;
						say $resp->decoded_content;
						my $delay = $attempt * 10;
						say "Retrying in $delay seconds...";
						sleep $delay;
						next ATTEMPT;
					}
					say "Giving up on locus $locus.";
					last PAGE;
				} else {
					last ATTEMPT;
				}
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
							$opts{'user_id'},
							$opts{'user_id'}
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
			usleep(500_000);    #Rate-limiting
			if ( @$alleles && $data->{'paging'}->{'next'} ) {
				$url = $data->{'paging'}->{'next'};
			} else {
				$script->{'db'}->commit if $opts{'commit'};
				last PAGE;
			}
		}
		$script->{'db'}->commit;
	}
	return;
}

sub update_profiles {
	check_options(qw(d e s scheme_id));
	my $loci        = get_loci();
	my $mapped_loci = get_mapped_loci();
	usleep(500_000);    #Rate-limiting
	local $" = q(,);
	my $existing_alleles = {};
	foreach my $locus (@$loci) {
		my $locus_name = $mapped_loci->{$locus} // $locus;
		my $allele_ids = $script->{'datastore'}->run_query( 'SELECT allele_id FROM sequences WHERE locus=?',
			$locus_name, { fetch => 'col_arrayref', cache => 'get_all_allele_ids' } );
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
			my $locus_name = $mapped_loci->{$locus} // $locus;
			$existing->{$st}->{$locus} = $profile->[ $positions->{$locus_name} ];
		}
	}
	my $url = "$SERVER_ADDRESS/$opts{'e'}/$opts{'s'}/sts?show_alleles=true";
	$url .= "&limit=$opts{'limit'}"     if $opts{'limit'};
	$url .= "&reldate=$opts{'reldate'}" if $opts{'reldate'};
	my %already_received;
  PAGE: while (1) {
		my $resp;
	  ATTEMPT: for my $attempt ( 1 .. 10 ) {
			$resp = $ua->get($url);
			if ( !$resp->is_success ) {
				if ( $attempt < 10 ) {
					say $url;
					say $resp->status_line;
					say $resp->decoded_content;
					my $delay = $attempt * 10;
					say "Retrying in $delay seconds...";
					sleep $delay;
					next ATTEMPT;
				}
				say "Giving up.";
				last PAGE;
			} else {
				last ATTEMPT;
			}
		}
		my $data     = decode_json( $resp->decoded_content );
		my $profiles = $data->{'STs'};
	  PROFILE: foreach my $profile (@$profiles) {
			my $st = $profile->{'ST_id'};
			my $cc = $profile->{'info'}->{'st_complex'};
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
			if ( keys %alleles != @$loci ) {
				if ( $opts{'allow_missing'} ) {
					foreach my $locus (@$loci) {
						$alleles{$locus} //= 'N';
					}
				} else {
					next PROFILE;
				}
			}
			if ( $existing->{$st} ) {
				my $changed;
				foreach my $locus (@$loci) {
					my $locus_name = $mapped_loci->{$locus} // $locus;

					#					say "$locus $alleles{$locus} $existing->{$st}->{$locus}";
					$changed = 1 if $alleles{$locus} ne $existing->{$st}->{$locus};
				}
				if ($changed) {
					say "ST-$st has changed! Deleting old version";
					$script->{'db'}
					  ->do( 'DELETE FROM profiles WHERE (scheme_id,profile_id)=(?,?)', undef, $opts{'scheme_id'}, $st );
				}
				next PROFILE;
			} else {
				foreach my $locus (@$loci) {
					my $locus_name = $mapped_loci->{$locus} // $locus;
					if ( !$existing_alleles->{$locus}->{ $alleles{$locus} } ) {
						say "Cannot insert ST-$st - $locus-$alleles{$locus} is not defined!" if !$opts{'no_errors'};
						next PROFILE;
					}
				}
				say "Inserting ST-$st.";
				eval {
					$script->{'db'}->do(
						'INSERT INTO profiles (scheme_id,profile_id,sender,curator,'
						  . 'date_entered,datestamp) VALUES (?,?,?,?,?,?)',
						undef, $opts{'scheme_id'}, $st, $opts{'user_id'}, $opts{'user_id'}, 'now', 'now'
					);
					$script->{'db'}->do(
						'INSERT INTO profile_fields (scheme_id,profile_id,scheme_field,value,curator,'
						  . 'datestamp) VALUES (?,?,?,?,?,?)',
						undef, $opts{'scheme_id'}, $st, 'ST', $st, $opts{'user_id'}, 'now'
					);
					if ( defined $cc ) {
						$script->{'db'}->do(
							'INSERT INTO profile_fields (scheme_id,profile_id,scheme_field,value,curator,'
							  . 'datestamp) VALUES (?,?,?,?,?,?)',
							undef, $opts{'scheme_id'}, $st, 'clonal_complex', $cc, $opts{'user_id'}, 'now'
						);
					}
					foreach my $locus (@$loci) {
						my $locus_name = $mapped_loci->{$locus} // $locus;
						$alleles{$locus} = 'N' if $alleles{$locus} eq '0';
						$script->{'db'}->do(
							'INSERT INTO profile_members (scheme_id,locus,profile_id,allele_id,'
							  . 'curator,datestamp) VALUES (?,?,?,?,?,?)',
							undef, $opts{'scheme_id'}, $locus_name, $st, $alleles{$locus}, $opts{'user_id'}, 'now'
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
		usleep(500_000);    #Rate-limiting
		if ( @$profiles && $data->{'links'}->{'paging'}->{'next'} ) {
			$url = $data->{'links'}->{'paging'}->{'next'};
			$url .= q(&show_alleles=true);
			$script->{'db'}->commit if $opts{'commit'};
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
	my $url         = "$SERVER_ADDRESS/$opts{'e'}/$opts{'s'}/loci?limit=2000";
	my $locus_names = [];
	while (1) {
		my $resp = $ua->get($url);
		if ( $resp->is_success ) {
			my $data = decode_json( $resp->decoded_content );
			my $loci = $data->{'loci'};
			foreach my $locus (@$loci) {
				push @$locus_names, $locus->{'locus'};
			}
			if ( @$loci && $data->{'links'}->{'paging'}->{'next'} ) {
				$url = $data->{'links'}->{'paging'}->{'next'};
				usleep(500_000);
			} else {
				return $locus_names;
			}
		} else {
			say $url;
			say $resp->status_line;
			say $resp->decoded_content;
			exit;
		}
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
    
${bold}--commit$norm
    Commit after every download page.
    
${bold}-d, --database$norm ${under}NAME$norm
    Database configuration name.

${bold}-e, --enterobase_db$norm ${under}DBASE$norm
    Enterobase database name.

${bold}-h, --help$norm
    This help page.
    
${bold}-l, --loci$norm
    Retrieve list of loci.
    
${bold}--limit$norm ${under}LIMIT$norm
    Request LIMIT number of alleles or profiles per page. 
    
${bold}--locus_regex$norm ${under}REGEX$norm
    Restrict updated loci to those that match regular expression. 
    
${bold}-m, --allow_missing$norm
    Allow profile definitions with missing loci - an N will be used for the
    missing alleles.
    
${bold}-n, --new_loci$norm
    Only download new loci (those without any existing alleles defined).
    
${bold}--no_errors$norm
    Silently ignore insertion errors caused by profiles being defined on
    Enterobase for which no alleles are defined.
    
${bold}-p, --update_profiles$norm
    Update allelic profiles.
    
${bold}-r, --route$norm
    Relative route.
    
${bold}--reldate$norm ${under}DAYS$norm
    Only return data from added in the last set number of days.
    
${bold}-s, --scheme$norm ${under}SCHEME NAME$norm
    Scheme name.
    
${bold}--scheme_id$norm ${under}SCHEME ID$norm
    Scheme id number in BIGSdb database.
    
${bold}--user_id$norm ${under}USER_ID$norm
    User id number for sender/curator
HELP
	return;
}
