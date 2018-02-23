#!/usr/bin/perl
#Wrapper script for rMLST database species identifier code.
#Written by Keith Jolley, 2018
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	RMLST_DATABASE   => 'pubmlst_rmlst_isolates',
	HOST             => undef,                                           #Use values in config.xml
	PORT             => undef,                                           #But you can override here.
	USER             => undef,
	PASSWORD         => undef,
	ANALYSIS_SCRIPT  => '/usr/local/TaxonomyFinder/TaxonomyFinder.pl',
	TOTAL_LOCI       => 53
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use JSON;
use BIGSdb::Offline::Script;
use BIGSdb::Utils;
Log::Log4perl->init( CONFIG_DIR . '/rest_logging.conf' );
my $logger = Log::Log4perl::get_logger('BIGSdb.Rest');
my $script = get_script_obj();
main();
undef $script;

sub main {
	my $json_results = $ARGV[0];
	die "No results file passed.\n" if !$json_results;
	die "Results file $json_results does not exist.\n" if !-e $json_results;
	my $json_ref = BIGSdb::Utils::slurp($json_results);
	my $prefix   = BIGSdb::Utils::get_random();
	if ($json_ref) {
		my $results       = decode_json($$json_ref);
		my $taxonomy_file = "$script->{'config'}->{'secure_tmp_dir'}/${prefix}_taxonomy.txt";
		my $scan_file     = "$script->{'config'}->{'secure_tmp_dir'}/${prefix}_scan.txt";
		my $out_file      = "$script->{'config'}->{'secure_tmp_dir'}/${prefix}_analysis.txt";
		make_taxonomy_file( $taxonomy_file, $results );
		make_scan_file( $scan_file, $results );
		my ( $unlinked_matches, $total_matches ) = count_unlinked_matches($results);
		my @params = (
			-in       => $scan_file,
			-taxonomy => $taxonomy_file,
			-out      => $out_file
		);
		local $" = q( );
		my $command = ANALYSIS_SCRIPT . qq( @params >/dev/null 2>&1);
		eval { system($command); };
		$logger->error($@) if $@;

		if ( -e $out_file ) {
			my $analysis = BIGSdb::Utils::slurp($out_file);
			my $output = parse_analysis( $analysis, $total_matches );
			say encode_json($output);
		}
		unlink $taxonomy_file, $scan_file, $out_file, "${out_file}.info";
		return;
	}
	say q(Script failed!);
	return;
}

sub count_unlinked_matches {
	my ($results) = @_;
	my @loci      = keys %{ $results->{'exact_matches'} };
	my $unlinked  = 0;
	my $matches   = 0;
	foreach my $locus (@loci) {
		my $locus_matches = $results->{'exact_matches'}->{$locus};
		next if !ref $locus_matches;
		foreach my $match (@$locus_matches) {
			my $allele_species = $match->{'linked_data'}->{'rMLST genome database'}->{'species'};
			$unlinked++ if !$allele_species;
			$matches++;
		}
	}
	return ( $unlinked, $matches );
}

sub parse_analysis {
	my ( $analysis, $total_matches ) = @_;
	my @ranks = qw(phylum class order family genus species);
	my @lines = split /\n/x, $$analysis;
	my @matches;
	foreach my $line (@lines) {
		last if $line =~ /^\#/x;
		my @record = split /\t/x, $line;
		my $taxonomy = $record[13];
		$taxonomy =~ s/^\|ROOT\|//x;
		$taxonomy =~ s/\|$//x;
		my @rank_values = split /\|/x, $taxonomy;
		my $tax_string;
		foreach my $i ( 0 .. 5 ) {
			last if !$rank_values[$i];
			$tax_string .= ' > ' if $i;
			$rank_values[$i] = '[unclassified]' if $rank_values[$i] eq 'NULL';
			$tax_string .= qq($rank_values[$i]);
		}
		my $percent_support_across_all = $record[8];
		my $allele_matches             = $record[9];
		my $support;
		if ($total_matches) {
			my $percent_linked_allele_matches = ( 100 * $allele_matches ) / $total_matches;
			$support = ( $percent_support_across_all * $percent_linked_allele_matches ) / 100;
		} else {
			$support = 0;
		}
		push @matches,
		  {
			rank     => $record[3],
			taxon    => $record[4],
			support  => int($support),
			taxonomy => $tax_string
		  };
		@matches = sort { $b->{'support'} <=> $a->{'support'} || $a->{'taxon'} cmp $b->{'taxon'} } @matches;
	}
	return \@matches;
}

sub make_scan_file {
	my ( $filename, $results ) = @_;
	my @loci = keys %{ $results->{'exact_matches'} };
	open( my $fh, '>', $filename ) || $logger->error("Cannot open $filename for writing.");
	foreach my $locus (@loci) {
		my $matches = $results->{'exact_matches'}->{$locus};
		next if !ref $matches;
		foreach my $match (@$matches) {
			my $allele_species = $match->{'linked_data'}->{'rMLST genome database'}->{'species'};
			next if !$allele_species;
			my $allele_id = $match->{'allele_id'};
			foreach my $species (@$allele_species) {
				say $fh qq($locus\t$allele_id\t$species->{'value'});
			}
		}
	}
	close $fh;
	return;
}

sub make_taxonomy_file {
	my ( $filename, $results ) = @_;
	my @loci = keys %{ $results->{'exact_matches'} };
	my %species;
	foreach my $locus (@loci) {
		foreach my $match ( @{ $results->{'exact_matches'}->{$locus} } ) {
			my $allele_id      = $match->{'allele_id'};
			my $allele_species = $match->{'linked_data'}->{'rMLST genome database'}->{'species'};
			next if !$allele_species;
			$species{ $_->{'value'} } = 1 foreach @$allele_species;
		}
	}
	my $temp_table = $script->{'datastore'}->create_temp_list_table_from_array( 'text', [ keys %species ] );
	my $data = $script->{'datastore'}->run_query(
		q{SELECT a.isolate_field,a.attribute,a.field_value,a.value FROM isolate_value_extended_attributes a }
		  . qq{JOIN $temp_table t ON field_value=t.value WHERE isolate_field='species'},
		undef,
		{ fetch => 'all_arrayref', slice => {} }
	);
	open( my $fh, '>', $filename ) || $logger->error("Cannot open $filename for writing.");
	say $fh qq(isolate_field\tattribute\tfield_value\tvalue);
	local $" = qq(\t);
	foreach my $row (@$data) {
		say $fh qq(@{$row}{qw(isolate_field attribute field_value value)});
	}
	close $fh;
	return;
}

sub get_script_obj {
	my $script_obj = BIGSdb::Offline::Script->new(
		{
			config_dir       => CONFIG_DIR,
			lib_dir          => LIB_DIR,
			dbase_config_dir => DBASE_CONFIG_DIR,
			host             => HOST,
			port             => PORT,
			user             => USER,
			password         => PASSWORD,
			options          => { always_run => 1, logger => $logger },
			instance         => RMLST_DATABASE,
		}
	);
	return $script_obj;
}
