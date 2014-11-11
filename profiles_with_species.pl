#!/usr/bin/perl -T
#Generate tab-delimited text file of MLST profiles including a
#a column for species
#Written by Keith Jolley 2014
use strict;
use warnings;
use 5.010;
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	HOST             => 'zoo-oban',
	PORT             => 5432,
	USER             => 'apache',
	PASSWORD         => ''
};
#######End Local configuration###############################
use lib (LIB_DIR);
use List::MoreUtils qw(any);
use Getopt::Std;
use BIGSdb::Utils;
use BIGSdb::Offline::Script;
my %opts;
getopts( 'a:b:g:n:N:p:s:ht', \%opts );

if ( $opts{'h'} ) {
	show_help();
	exit;
}
if ( !$opts{'a'} ) {
	say "Usage: profiles_with_species.pl -a <seqdef db config>";
	say "Help: profiles_with_species.pl -h";
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
		instance         => $opts{'a'},
	}
);
main();
undef $script;

sub main {
	my ($self) = @_;
	my $scheme_id = BIGSdb::Utils::is_int( $opts{'s'} ) ? $opts{'s'} : 1;
	initiate_isolate_db();
	print_header($scheme_id);
	print_scheme($scheme_id);
	return;
}

sub print_header {
	my ($scheme_id) = @_;
	my $scheme_info = $script->{'datastore'}->get_scheme_info( $scheme_id, { get_pk => 1 } );
	if ( !defined $scheme_info->{'primary_key'} ) {
		$script->{'logger'}->error("No primary key set.");
		undef $script;
		exit;
	}
	my $loci = $script->{'datastore'}->get_scheme_loci($scheme_id);
	if ( $opts{'p'} ) {
		s/^$opts{'p'}// foreach @$loci;
	}
	my $scheme_fields = $script->{'datastore'}->get_scheme_fields($scheme_id);
	local $" = "\t";
	print "$scheme_info->{'primary_key'}\t@$loci";
	foreach my $field (@$scheme_fields) {
		next if $opts{'t'} && $field eq 'species';
		print "\t$field" if $field ne $scheme_info->{'primary_key'};
	}
	say "\tspecies";
	return;
}

sub print_scheme {
	my ($scheme_id) = @_;
	my $scheme_info   = $script->{'datastore'}->get_scheme_info( $scheme_id, { get_pk => 1 } );
	my $loci          = $script->{'datastore'}->get_scheme_loci($scheme_id);
	my $scheme_fields = $script->{'datastore'}->get_scheme_fields($scheme_id);
	my $matview =
	  $script->{'datastore'}->run_simple_query( "SELECT EXISTS(SELECT 1 FROM matviews WHERE v_name=?)", "scheme_$scheme_id" )->[0];
	my $view = $matview ? "mv_scheme_$scheme_id" : "scheme_$scheme_id";
	my $profiles = $script->{'datastore'}->run_query( "SELECT * FROM $view ORDER BY CAST($scheme_info->{'primary_key'} AS INT)",
		undef, { fetch => 'all_arrayref', slice => {} } );
	my ( $isolates_sql, $species_sql );

	if ( $opts{'b'} ) {

		#Check species field exists
		my $check_sql =
		  $script->{'dbi'}
		  ->prepare( "SELECT EXISTS(SELECT * FROM information_schema.columns WHERE " . "table_name='isolates' AND column_name='species')" );
		eval { $check_sql->execute };
		$script->{'logger'}->error($@) if $@;
		my ($field_exists) = $check_sql->fetchrow_array;
		if ( !$field_exists ) {
			delete $opts{'b'};
		} else {
			my $locus_count = @$loci;
			my $qry =
			  "SELECT isolates.id FROM isolates LEFT JOIN allele_designations ON " . "isolates.id=allele_designations.isolate_id WHERE ";
			my @locus_clause;
			foreach my $locus (@$loci) {
				my $cleaned_locus = $locus;
				if ( $opts{'p'} ) {
					$cleaned_locus =~ s/^$opts{'p'}//;
				}
				push @locus_clause, "(locus='$cleaned_locus' AND allele_id=?)";
			}
			local $" = ' OR ';
			$qry .= "@locus_clause GROUP BY isolates.id HAVING COUNT(isolates.id)=$locus_count";
			$isolates_sql = $script->{'dbi'}->prepare($qry);
			$species_sql  = $script->{'dbi'}->prepare("SELECT species FROM isolates WHERE id=?");
		}
	}
	my %ignore_species;
	if ( $opts{'N'} ) {
		my @ignore = split ',', $opts{'N'};
		%ignore_species = map { $_ => 1 } @ignore;
	}
	foreach my $profile (@$profiles) {
		print $profile->{ lc( $scheme_info->{'primary_key'} ) };
		my @allelic_profile;
		push @allelic_profile, $profile->{ lc($_) } foreach @$loci;
		local $" = "\t";
		print "\t@allelic_profile";
		foreach my $field (@$scheme_fields) {
			next if $field eq $scheme_info->{'primary_key'};
			next if $opts{'t'} && $field eq 'species';
			my $value = $profile->{ lc($field) } // '';
			print "\t$value";
		}
		my $species_name;
		if ( $opts{'b'} ) {
			my %species_count;
			eval { $isolates_sql->execute(@allelic_profile) };
			$script->{'logger'}->error($@) if $@;
			while ( ( my $id ) = $isolates_sql->fetchrow_array ) {
				eval { $species_sql->execute($id) };
				$script->{'logger'}->error($@) if $@;
				my ($species) = $species_sql->fetchrow_array // '';
				if ( $opts{'g'} ) {
					my $initial = substr( $opts{'g'}, 0, 1 );
					$species =~ s/^$initial\./$opts{'g'}/;
				}
				$species =~ s/ or .*$//;    #Special case of alternative names in Aeromonas
				$species_count{$species}++ if $species && !$ignore_species{$species};
			}
			my @species_names = sort { $species_count{$b} <=> $species_count{$a} } keys %species_count;
			local $" = ',';
			$species_name = "@species_names";
		}
		$species_name ||= $opts{'n'} // '';
		print "\t$species_name\n";
	}
}

sub initiate_isolate_db {
	return if !$opts{'b'};
	my $xml_handler = BIGSdb::Parser->new;
	my $parser      = XML::Parser::PerlSAX->new( Handler => $xml_handler );
	my $full_path   = DBASE_CONFIG_DIR . "/$opts{'b'}/config.xml";
	eval { $parser->parse( Source => { SystemId => $full_path } ); };
	if ($@) {
		$script->{'logger'}->fatal("Invalid XML description: $@");
		say "Invalid database $opts{'b'}";
		undef $script;
		exit;
	}
	my $system_hash = $xml_handler->get_system_hash;
	$script->{'dbi'} = $script->{'dataConnector'}->get_connection( { dbase_name => $system_hash->{'db'} } );
}

sub show_help {
	print << "HELP";

Usage profiles_with_species.pl -a <seqdef database config> -s <scheme id>

Options
-------
-a <name>  Seqdef database configuration name.
-b <name>  Isolate database configuration name.
-g <genus> Expand genus name when single initial used in species name.
-h         This help page.
-n         Default species name (use if not defined in isolate database)
-N         Ignore species name (comma-separated list)
-p         Locus prefix (this will be stripped out)
-s <id>    Scheme id.
-t         Ignore species field in profile definition

HELP
	return;
}
