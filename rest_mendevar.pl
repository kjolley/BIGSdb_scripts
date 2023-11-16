#!/usr/bin/perl
#Run MendeVar from REST API
#Written by Keith Jolley, 2022-2023
#Version: 20231116
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	DATABASE         => 'pubmlst_neisseria_seqdef',
	MIN_ALIGNMENT    => 50,
	MIN_IDENTITY     => 70
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use JSON;
use BIGSdb::Offline::Script;
use BIGSdb::Utils;
use List::MoreUtils qw(uniq);
use FindBin;
use lib "$FindBin::Bin/lib";
use MenVaccine;
Log::Log4perl->init( CONFIG_DIR . '/rest_logging.conf' );
my $logger   = Log::Log4perl::get_logger('BIGSdb.Rest');
my $script   = get_script_obj();
my @loci     = qw(fHbp_peptide NHBA_peptide NadA_peptide PorA_VR2);
my %nuc_loci = (
	fHbp_peptide => 'NEIS0349',
	NHBA_peptide => 'NEIS2109',
	NadA_peptide => 'NEIS1969'
);
main();
undef $script;

sub main {
	my $json_results = $ARGV[0];
	die "No results file passed.\n"                    if !$json_results;
	die "Results file $json_results does not exist.\n" if !-e $json_results;
	my $json_ref = BIGSdb::Utils::slurp($json_results);
	my $prefix   = BIGSdb::Utils::get_random();
	if ( !$json_ref ) {
		say q(Script failed!);
		return;
	}
	my $response;
	my $results = extract_results($json_ref);
	$response->{'processed_matches'} = $results;
	$response->{'bexsero'}           = process_bexsero_outcome($results);
	$response->{'trumenba'}          = process_trumenba_outcome($results);
	say encode_json($response);
	return;
}

sub extract_results {
	my ($json_ref) = @_;
	my $processed  = {};
	my $results    = decode_json($$json_ref);
  LOCUS: foreach my $locus (@loci) {
		if ( $results->{'exact_matches'}->{$locus} && @{ $results->{'exact_matches'}->{$locus} } ) {
			$processed->{$locus} = { variants => $results->{'exact_matches'}->{$locus} };
			next LOCUS;
		}
		if ( defined $nuc_loci{$locus} ) {
			if ( $results->{'exact_matches'}->{ $nuc_loci{$locus} }
				&& @{ $results->{'exact_matches'}->{ $nuc_loci{$locus} } } )
			{
			  ALLELE: foreach my $match ( @{ $results->{'exact_matches'}->{ $nuc_loci{$locus} } } ) {
					my $allele_id = ref $match ? $match->{'allele_id'} : 0;
					my $flags     = $script->{'datastore'}->run_query(
						'SELECT flag FROM allele_flags WHERE (locus,allele_id)=(?,?)',
						[ $nuc_loci{$locus}, $allele_id ],
						{ fetch => 'col_arrayref' }
					);
					my %flags = map { $_ => 1 } @$flags;
					if ( $flags{'internal stop codon'} || $flags{'frameshift'} || $flags{'contains IS element'} ) {
						$processed->{$locus} = {
							variants => [0],
							text     => 'CDS has frameshift, internal stop codon or IS element'
						};
					} else {
						$processed->{$locus} = {
							variants => [0],
							text => "allele identified ($nuc_loci{$locus}:$allele_id), but peptide variant not defined"
						};
					}
					last ALLELE;
				}
				next LOCUS;
			}
			if ( $results->{'partial_matches'}->{ $nuc_loci{$locus} }->{'allele'} ) {
				my $match = $results->{'partial_matches'}->{ $nuc_loci{$locus} };
				if ( ( 100 * $match->{'alignment'} / $match->{'length'} ) < MIN_ALIGNMENT
					|| $match->{'identity'} < MIN_IDENTITY )
				{
					my ( $alignment, $identity ) = ( MIN_ALIGNMENT, MIN_IDENTITY );
					$processed->{$locus} = { variants => [0] };
				} else {
					$processed->{$locus} = { variants => ['possible new variant'] };
				}
			} else {
				$processed->{$locus} = { variants => [0] };
			}
		} else {
			$processed->{$locus} = { variants => [0] };
		}
	}
	return $processed;
}

sub process_bexsero_outcome {
	my ($results) = @_;
	if ( !contains_bexsero_antigens($results) ) {
		return { result => 'none', notes => 'none of the vaccine antigens found' };
	}
	local $" = q(; );
	my $has_exact_bexsero_match = has_exact_bexsero_match($results);
	if ( $has_exact_bexsero_match->{'result'} ) {
		return { result => 'exact match', notes => qq(@{$has_exact_bexsero_match->{'notes'}}) };
	}
	my $cross_reacts_with_bexsero = cross_reacts_with_bexsero($results);
	if ( $cross_reacts_with_bexsero->{'result'} ) {
		return { result => 'cross-reactive', notes => qq(@{$cross_reacts_with_bexsero->{'notes'}}) };
	}
	my $no_reactivity_with_bexsero = no_reactivity_with_bexsero($results);
	if ( $no_reactivity_with_bexsero->{'result'} ) {
		return { result => 'none', notes => qq(@{$no_reactivity_with_bexsero->{'notes'}}) };
	}
	return {
		result => 'insufficient data',
		notes  => 'insufficient data about cross-reactivity of antigens to make an assessment'
	};
}

sub process_trumenba_outcome {
	my ($results) = @_;
	if ( !contains_trumenba_antigens($results) ) {
		return { result => 'none', notes => 'fHbp_peptide is missing' };
	}
	local $" = q(; );
	my $has_exact_trumenba_match = has_exact_trumenba_match($results);
	if ( $has_exact_trumenba_match->{'result'} ) {
		return { result => 'exact match', notes => qq(@{$has_exact_trumenba_match->{'notes'}}) };
	}
	my $cross_reacts_with_trumenba = cross_reacts_with_trumenba($results);
	if ( $cross_reacts_with_trumenba->{'result'} ) {
		return { result => 'cross-reactive', notes => qq(@{$cross_reacts_with_trumenba->{'notes'}}) };
	}
	return {
		result => 'insufficient data',
		notes  => 'insufficient data about cross-reactivity of antigens to make an assessment'
	};
}

sub contains_bexsero_antigens {
	my ($results) = @_;
	foreach my $locus (@loci) {
		return 1 if $results->{$locus}->{'variants'} && $results->{$locus}->{'variants'}->[0] ne '0';
	}
	return;
}

sub contains_trumenba_antigens {
	my ($results) = @_;
	return 1 if $results->{'fHbp_peptide'}->{'variants'} && $results->{'fHbp_peptide'}->{'variants'}->[0] ne '0';
	return;
}

sub has_exact_bexsero_match {
	my ($results) = @_;
	my $components = MenVaccine::get_exact_bexsero_match();
	return format_match_results( $results, $components, 'is exact match to vaccine variant' );
}

sub has_exact_trumenba_match {
	my ($results) = @_;
	my $components = MenVaccine::get_exact_trumenba_match();
	return format_match_results( $results, $components, 'is exact match to vaccine variant' );
}

sub cross_reacts_with_bexsero {
	my ($results) = @_;
	my $components = MenVaccine::get_cross_reacts_with_bexsero();
	return format_match_results( $results, $components, 'is cross-reactive to vaccine variant' );
}

sub cross_reacts_with_trumenba {
	my ($results) = @_;
	my $components = MenVaccine::get_cross_reacts_with_trumenba();
	return format_match_results( $results, $components, 'is cross-reactive to vaccine variant' );
}

sub no_reactivity_with_bexsero {
	my ($results) = @_;
	my $notes = [];

	#FHbp
	my $components = MenVaccine::get_fhbp_no_reactivity_with_bexsero();
	my $matched    = is_not_cross_reactive( $results, $components, $notes );
	if ( ( $results->{'fHbp_peptide'}->{'variants'}->[0] // '0' ) eq '0' ) {
		$matched = 1;
		push @$notes, qq(fHbp_peptide is missing);
	}
	return { result => 0 } if !$matched;
	$components = MenVaccine::get_nhba_no_reactivity_with_bexsero();

	#NHBA
	$matched = is_not_cross_reactive( $results, $components, $notes );
	if ( ( $results->{'NHBA_peptide'}->{'variants'}->[0] // '0' ) eq '0' ) {
		$matched = 1;
		push @$notes, qq(NHBA_peptide is missing);
	}
	return { result => 0 } if !$matched;

	#NadA
	$components = MenVaccine::get_nadA_no_reactivity_with_bexsero();
	$matched    = is_not_cross_reactive( $results, $components, $notes );
	if ( ( $results->{'NadA_peptide'}->{'variants'}->[0] // '0' ) eq '0' ) {
		$matched = 1;
		push @$notes, qq(NadA_peptide is missing);
	}
	return { result => 0 } if !$matched;

	#PorA
	if ( ( $results->{'PorA_VR2'}->{'variants'}->[0] // '0' ) eq '0' ) {
		push @$notes, qq(PorA_VR2 is missing);
		return { result => 1, notes => $notes };
	}
	if ( ( $results->{'PorA_VR2'}->{'variants'}->[0] // '0' ) ne '4' ) {
		push @$notes, qq(PorA_VR2 is not variant 4);
	}
	return { result => 1, notes => $notes };
}

sub is_not_cross_reactive {
	my ( $results, $components, $notes ) = @_;
	my $matched = 0;
	foreach my $component (@$components) {
		my $variants =
		  ref $results->{ $component->{'locus'} }->{'variants'}
		  ? $results->{ $component->{'locus'} }->{'variants'}
		  : [];
		my %variant_found;
		foreach my $variant (@$variants) {
			next if $variant eq '0';
			$variant_found{ $variant->{'allele_id'} } = 1;
		}
		if ( $variant_found{ $component->{'variant'} } ) {
			$matched = 1;
			my @note_evidence;
			my @evidence_list  = uniq( sort { $a cmp $b } values %{ $component->{'references'} } );
			my $evidence_count = 0;
			foreach my $evidence (@evidence_list) {
				$evidence_count++;
				my @refs;
				foreach my $ref ( sort { $a <=> $b } keys %{ $component->{'references'} } ) {
					push @refs, $ref if $component->{'references'}->{$ref} eq $evidence;
				}
				if (@refs) {
					if ( $evidence_count == 1 ) {
						$evidence = qq(data derived from $evidence assays);
					} elsif ( $evidence_count == @evidence_list ) {
						$evidence = qq(and $evidence assays);
					} else {
						$evidence = qq(, $evidence assays);
					}
					local $" = q(, PMID:);
					push @note_evidence, qq($evidence (PMID:@refs));
				}
			}
			if (@note_evidence) {
				local $" = q(; );
				push @$notes,
				  qq($component->{'locus'}: $component->{'variant'} is not )
				  . qq(cross-reactive with vaccine variant - @note_evidence);
			}
		}
	}
	return $matched;
}

sub format_match_results {
	my ( $results, $components, $reason ) = @_;
	my $result;
	my $notes  = [];
	my %assays = map { $_ => 1 } qw(MATS SBA MEASURE);
	foreach my $component (@$components) {
		my $match_variants =
		  ref $results->{ $component->{'locus'} }->{'variants'}
		  ? $results->{ $component->{'locus'} }->{'variants'}
		  : [];
		my %match_variants;
		foreach my $variant (@$match_variants) {
			if ( ref $variant ) {
				$match_variants{ $variant->{'allele_id'} } = 1;
			} elsif ( $variant eq '0' ) {
				$match_variants{'0'} = 1;
			}
		}
		if ( $match_variants{ $component->{'variant'} } ) {
			$result = 1;
			my @note_evidence;
			my @evidence_list  = uniq( sort { $a cmp $b } values %{ $component->{'references'} } );
			my $evidence_count = 0;
			foreach my $evidence (@evidence_list) {
				$evidence_count++;
				my @refs;
				foreach my $ref ( sort { $a <=> $b } keys %{ $component->{'references'} } ) {
					push @refs, $ref if $component->{'references'}->{$ref} eq $evidence;
				}
				if (@refs) {
					local $" = q(, PMID:);
					if ( $assays{$evidence} ) {
						if ( $evidence_count == 1 ) {
							$evidence = qq(data derived from $evidence assays);
						} elsif ( $evidence_count == @evidence_list ) {
							$evidence = qq(and $evidence assays);
						} else {
							$evidence = qq(, $evidence assays);
						}
					}
					push @note_evidence, qq($evidence (PMID:@refs));
				}
			}
			if (@note_evidence) {
				local $" = q(, );
				push @$notes, qq($component->{'locus'}: $component->{'variant'} $reason - @note_evidence);
			}
		}
	}
	return { result => $result, notes => $notes };
}

sub get_script_obj {
	my $script_obj = BIGSdb::Offline::Script->new(
		{
			config_dir       => CONFIG_DIR,
			lib_dir          => LIB_DIR,
			dbase_config_dir => DBASE_CONFIG_DIR,
			options          => { always_run => 1, logger => $logger },
			instance         => DATABASE,
		}
	);
	return $script_obj;
}
