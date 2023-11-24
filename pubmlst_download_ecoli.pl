#!/usr/bin/env perl
#Demo script to download basic metadata and contigs
#from PubMLST Escherichia database
#Written by Keith Jolley 2023
#License: GPL3
#Version 20231124
use REST::Client;
use LWP::Simple;
use JSON;
use strict;
use warnings;
use 5.010;
use constant MAX_ATTEMPTS  => 10;
use constant API_URI       => 'https://rest.pubmlst.org';
use constant DB            => 'pubmlst_escherichia_isolates';
use constant OUTPUT_DIR    => './';
use constant METADATA_FILE => 'metadata.tsv';
my $client = REST::Client->new();
main();

sub main {
	my $url = API_URI . '/db/' . DB;

	#Get list of genome URLs;
	my $response    = call("$url/genomes?return_all=1");
	my $genome_urls = $response->{'isolates'};
	my $count       = @$genome_urls;
	say "Found $count genomes.";
	my @fields = qw(id isolate aliases species country year ST source phylogroup pathotype
	  pubmed_ids bioproject_accession biosample_accession run_accession contig_count total_length);
	my %is_provenance_field =
	  map { $_ => 1 }
	  qw(id isolate species country year source phylogroup pathotype bioproject_accession
	  biosample_accession run_accession);
	my %is_contig_field = map { $_ => 1 } qw(contig_count total_length);
	local $" = qq(\t);              #Value separator (set to a tab)
	write_metadata(qq(@fields));    #Print header line

	foreach my $genome_url (@$genome_urls) {
		my $record = call($genome_url);
		my @values;
		foreach my $field (@fields) {
			if ( $is_provenance_field{$field} ) {
				my $field_value = $record->{'provenance'}->{$field};

				#If field has multiple values then convert to semi-colon separated list.
				if ( ref $field_value eq 'ARRAY' ) {
					local $" = q(;);
					$field_value = qq(@$field_value);
				}
				$field_value //= q();    #Set value to empty string if undefined.
				push @values, $field_value;
			}
			if ( $is_contig_field{$field} ) {
				my $field_value = $record->{'sequence_bin'}->{$field} // q();    #Set to empty string if not defined
				push @values, $field_value;
			}
			if ( $field eq 'aliases' ) {
				local $" = q(;);                                                 #Output as semi-colon separated list.
				my $aliases = $record->{'aliases'} // [];                        #Set to empty list if not defined
				push @values, qq(@$aliases);
			}
			if ( $field eq 'ST' ) {

				#We need to read through each of the scheme results
				my $schemes = $record->{'schemes'};
				my $value   = q();
				foreach my $scheme (@$schemes) {
					if ( $scheme->{'description'} eq 'MLST (Achtman)' ) {
						$value = $scheme->{'fields'}->{'ST'} // q();
						last;
					}
				}
				push @values, $value;
			}
			if ( $field eq 'pubmed_ids' ) {
				my @pubmed_ids;
				my $publications = $record->{'publications'} // [];
				foreach my $pub (@$publications) {
					push @pubmed_ids, $pub->{'pubmed_id'};
				}
				local $" = q(;);    #Output as semi-colon separated list.
				push @values, qq(@pubmed_ids);
			}
		}
		write_metadata( qq(@values), { append => 1 } );
		my $contigs_url       = $record->{'sequence_bin'}->{'contigs_fasta'};
		my $assembly_filename = OUTPUT_DIR . "id-$record->{'provenance'}->{'id'}.fasta";
		getstore( "$contigs_url?header=original_designation", $assembly_filename );
	}
}

sub write_metadata {
	my ( $line, $options ) = @_;
	say $line;
	my $file_open_method = $options->{'append'} ? '>>' : '>';
	my $path             = OUTPUT_DIR . METADATA_FILE;
	open( my $fh, "$file_open_method:encoding(utf8)", $path ) || die "Cannot write to $path.\n";
	say $fh $line;
	close $fh;
	return;
}

sub call {
	my ($url)        = @_;
	my %server_error = map { $_ => 1 } ( 500, 502, 503, 504 );
	my $attempt      = 0;
	my $response;
	while (1) {
		$attempt++;
		$response = $client->request( 'GET', $url );
		my $code = $response->responseCode;
		if ( $server_error{$code} ) {
			if ( $attempt > MAX_ATTEMPTS ) {
				die "API call failed multiple times on $url. Giving up!\n";
			}
			say "$code error from API. Retrying in 10s...";
			sleep 10;
			next;
		}
		last;
	}
	my $content = from_json( $response->responseContent );
	return $content;
}
