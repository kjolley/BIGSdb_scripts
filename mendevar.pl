#!/usr/bin/env perl
#Perform Neisseria vaccine reactivity prediction.
#Written by Keith Jolley, 2020.
#Version:20200709
use strict;
use warnings;
use 5.010;
use JSON;
use Log::Log4perl qw(get_logger);
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	SEQDEF_DATABASE  => 'pubmlst_neisseria_seqdef',
	MIN_ALIGNMENT    => 50,
	MIN_IDENTITY     => 70
};
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Utils;
use FindBin;
use lib "$FindBin::Bin/lib";
use MenVaccine;
use List::MoreUtils qw(uniq);
Log::Log4perl->init( CONFIG_DIR . '/logging.conf' );
my $logger   = Log::Log4perl::get_logger('BIGSdb.Page');
my $script   = get_script_obj();
my @loci     = qw(fHbp_peptide NHBA_peptide NadA_peptide PorA_VR2);
my %nuc_loci = (
	fHbp_peptide => 'NEIS0349',
	NHBA_peptide => 'NEIS2109',
	NadA_peptide => 'NEIS1969'
);
main();

sub main {
	my $json_results = $ARGV[0];
	die "No results file passed.\n" if !$json_results;
	die "Results file $json_results does not exist.\n" if !-e $json_results;
	my $json_ref = BIGSdb::Utils::slurp($json_results);
	if ( !$json_ref ) {
		say q(Script failed!);
		return;
	}
	my $results = extract_results($json_ref);
	say q(<div class="box" id="resultspanel"><div class="scrollable">);
	say q(<h2>Matches</h2>);
	say q(<dl class="data">);
	foreach my $locus (@loci) {
		local $" = q(, );
		say qq(<dt>$locus</dt>);
		my $variants = defined $results->{$locus}->{'variants'} ? qq(@{$results->{$locus}->{'variants'}}) : 'missing';
		if ( $variants eq '0' ) {
			$variants = 'missing';
			$variants .= " - $results->{$locus}->{'text'}" if $results->{$locus}->{'text'};
		}
		say qq(<dd>$variants</dd>);
	}
	say q(</dl>);
	say q(<h2>Vaccine cross-reactivity</h2>);
	my $outcome = process_bexsero_outcome($results);
	my $notes   = get_bexsero_notes();
	say q(<h3>Bexsero<sup>&reg;</sup></h3>);
	print_results( $outcome, $notes );
	$outcome = process_trumenba_outcome($results);
	$notes   = get_trumenba_notes();
	say q(<h3>Trumenba<sup>&reg;</sup></h3>);
	print_results( $outcome, $notes );
	return;
}

sub print_results {
	my ( $outcome, $info ) = @_;
	state $div = 1;
	$outcome->{'notes'} =~ s/PMID:(\d+)/PMID:<a href="https:\/\/pubmed.ncbi.nlm.nih.gov\/$1">$1<\/a>/gx;
	say q(<div style="display:flex">);
	say q(<div style="padding:0.5em 0.5em 0 0">);
	my $traffic_light = get_traffic_light( $outcome->{'result'} );
	say qq(<p style="font-size:3em;text-align:center">$traffic_light</p>);
	say q(<p style="background:#88a;color:#fff;font-size:1.5em;padding:0.2em 0.5em;text-align:center;)
	  . qq(margin-top:-0.5em;min-width:6em">$outcome->{'result'}</p>);
	say q(</div>);
	say qq(<div style="margin:-3em 1em 1em 1em">);
	say qq(<div id="notes_$div" class="expandable_retracted_large">);
	say qq(<h4>Explanation:</h4><p>$outcome->{'notes'}</p>);
	say q(<h4>Notes:</h4>);
	say qq(<div class="scrollable">$info</div></div>);
	say qq(<div class="expand_link" id="expand_notes_$div"><span class="fas fa-chevron-down"></span></div>);
	say q(</div>);
	say q(</div>);
	$div++;
}

sub process_bexsero_outcome {
	my ($results) = @_;
	if ( !contains_bexsero_antigens($results) ) {
		return { result => 'none', notes => '<ul><li>none of the vaccine antigens found</li></ul>' };
	}
	local $" = q(</li><li>);
	my $has_exact_bexsero_match = has_exact_bexsero_match($results);
	if ( $has_exact_bexsero_match->{'result'} ) {
		return { result => 'exact match', notes => qq(<ul><li>@{$has_exact_bexsero_match->{'notes'}}</li></ul>) };
	}
	my $cross_reacts_with_bexsero = cross_reacts_with_bexsero($results);
	if ( $cross_reacts_with_bexsero->{'result'} ) {
		return { result => 'cross-reactive', notes => qq(<ul><li>@{$cross_reacts_with_bexsero->{'notes'}}</li></ul>) };
	}
	my $no_reactivity_with_bexsero = no_reactivity_with_bexsero($results);
	if ( $no_reactivity_with_bexsero->{'result'} ) {
		return { result => 'none', notes => qq(<ul><li>@{$no_reactivity_with_bexsero->{'notes'}}</li></ul>) };
	}
	return {
		result => 'insufficient data',
		notes  => '<ul><li>insufficient data about cross-reactivity of antigens to make an assessment</li></ul>'
	};
}

sub get_bexsero_notes {
	return << 'NOTES';
<p style="color:#606060"><strong>Bexsero<sup>&reg;</sup> (4CMenB) is a multicomponent vaccine.</strong></p>
<ul>
<li>Protein-based meningococcal vaccines contain surface proteins as vaccine antigens, these proteins 
demonstrate nucleotide and amino acid sequence diversity.</li>
<li>Peptide sequence diversity can be analysed using the Bexsero Antigen Sequence Typing (BAST) 
scheme<sup>1</sup>.</li>
<li>Bexsero<sup>&reg;</sup> contains: fHbp peptide 1; NHBA peptide 2; NadA peptide 8; PorA VR2 4.</li>
</ul>
<p>The Deduced Vaccine Antigen Reactivity (MenDeVAR) Index was developed to combine multiple, complex data
that inform the reactivity of each vaccine against specific antigenic variants.</p>
<p style="color:#606060"><strong>The MenDeVAR Index:</strong></p>
<ul>
<li><span class="fa-stack" style="margin-left:0.5em;font-size:0.8em">
<span class="fas fa-circle fa-stack-2x" style="color:#009800"></span>
<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span> isolate contains &ge;1 exact 
sequence match to antigenic variants found in the vaccine.</li>
<li><span class="fa-stack" style="margin-left:0.5em;font-size:0.8em">
<span class="fas fa-circle fa-stack-2x" style="color:#ef780f"></span>
<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span> isolate contains &ge;1 antigenic 
variant deemed cross-reactive to vaccine variants through experimental studies.</li>
<li><span class="fa-stack" style="margin-left:0.5em;font-size:0.8em">
<span class="fas fa-circle fa-stack-2x" style="color:#d90013"></span>
<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span> all the isolate's antigenic 
variants have been deemed not cross-reactive to vaccine variants through experimental studies.</li>
<li><span class="fa-stack" style="margin-left:0.5em;font-size:0.8em">
<span class="fas fa-circle fa-stack-2x" style="color:#d3d3d3"></span>
<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span> isolate contains antigens 
for which there is insufficient data from or are yet to be tested in experimental studies.</li>
</ul>

<p style="color:#606060"><strong>It is important to understand the caveats to interpreting the MenDeVAR Index:</strong></p>
<p><strong>Source of data</strong> - These data combine multiple sources of information including: 
peptide sequence identity through whole genome sequencing; experimental assays developed as indirect 
measures of the breadth of vaccine protection against diverse meningococci; and assays developed to 
assess immunogenicity. The Meningococcal Antigen Typing System (MATS)<sup>2</sup> assay was used for 
Bexsero<sup>&reg;</sup>.</p>  <p><strong>Cross-reactivity definition</strong> - An antigenic variant 
was considered cross-reactive if it had been tested in &ge;5 isolates/subjects and was above the 
accepted threshold in &ge;75% of those isolates. This was established through combined analysis of 
published experimental studies (PMID provided for each variant), not from genomic data. These assays 
were based on serogroup B disease isolates.</p>  <p><strong>Protein expression</strong> - We have not 
inferred from genomic data, therefore there may be isolates that possess genes but do no express the 
protein <i>in vivo</i>.</p>

<p><strong>Age of vaccinees</strong> - For MATS assay development<sup>2</sup>, Bexsero<sup>&reg;</sup> 
vaccine recipients were infants who had received 3 doses of vaccine and then a booster at 12 months. 
The pooled sera used for the MATS assay were taken from the toddlers at 13 months of age.</p>

<ol>
<li>Brehony C, Rodrigues CMC, Borrow R, <i>et al.</i> Distribution of Bexsero<sup>&reg;</sup> Antigen 
Sequence Types (BASTs) in invasive meningococcal disease isolates: Implications for immunisation. 
<a href="https://www.ncbi.nlm.nih.gov/pubmed/27521232" target="_blank">Vaccine 2016 34(39):4690-7</a></li>
<li>Donnelly J, Medini D, Boccadifuoco G, <i>et al</i>. Qualitative and quantitative assessment of 
meningococcal antigens to evaluate the potential strain coverage of protein-based vaccines. 
<a href="https://www.ncbi.nlm.nih.gov/pubmed/20962280" target="_blank">Proc Natl Acad Sci USA 2010;107(45):19490-19495</a>
</li>
</ol>
<p><em>MenDeVAR is described in <a href="https://doi.org/10.1101/2020.08.18.256834">Rodrigues <i>et al.</i> 2020, 
bioRxiv 2020.08.18.256834</a>.  Please contact <a href="mailto:charlene.rodrigues@zoo.ox.ac.uk">us</a> 
if you have queries.</em></p>
NOTES
}

sub process_trumenba_outcome {
	my ($results) = @_;
	if ( !contains_trumenba_antigens($results) ) {
		return { result => 'none', notes => '<ul><li>fHbp_peptide is missing</li></ul>' };
	}
	local $" = q(</li><li>);
	my $has_exact_trumenba_match = has_exact_trumenba_match($results);
	if ( $has_exact_trumenba_match->{'result'} ) {
		return { result => 'exact match', notes => qq(<ul><li>@{$has_exact_trumenba_match->{'notes'}}</li></ul>) };
	}
	my $cross_reacts_with_trumenba = cross_reacts_with_trumenba($results);
	if ( $cross_reacts_with_trumenba->{'result'} ) {
		return { result => 'cross-reactive', notes => qq(<ul><li>@{$cross_reacts_with_trumenba->{'notes'}}</li></ul>) };
	}
	return {
		result => 'insufficient data',
		notes  => '<ul><li>insufficient data about cross-reactivity of antigens to make an assessment</li></ul>'
	};
}

sub get_trumenba_notes {
	return << 'NOTES';
<p style="color:#606060"><strong>Trumenba<sup>&reg;</sup> (rLP2086) is a bivalent fHbp-containing vaccine.</strong></p>
<ul>
<li>Protein-based meningococcal vaccines contain surface proteins as vaccine antigens, these proteins demonstrate 
nucleotide and amino acid sequence diversity.</li>
<li>Peptide sequence diversity can be analysed using the fhbp peptide locus.</li>
<li>Trumenba<sup>&reg;</sup> vaccine contains fHbp peptides 45 and 551.</li>
</ul>
<p>The Deduced Vaccine Antigen Reactivity (MenDeVAR ) Index was developed to combine multiple, complex data that 
inform the reactivity of each vaccine against specific antigenic variants.</p>

<p style="color:#606060"><strong>The MenDeVAR Index:</strong></p>
<ul>
<li><span class="fa-stack" style="margin-left:0.5em;font-size:0.8em">
<span class="fas fa-circle fa-stack-2x" style="color:#009800"></span>
<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span> isolate contains &ge;1 exact sequence 
match to antigenic variants found in the vaccine.</li>
<li><span class="fa-stack" style="margin-left:0.5em;font-size:0.8em">
<span class="fas fa-circle fa-stack-2x" style="color:#ef780f"></span>
<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span> isolate contains &ge;1 antigenic 
variant deemed cross-reactive to vaccine variants through experimental studies.</li>
<li><span class="fa-stack" style="margin-left:0.5em;font-size:0.8em">
<span class="fas fa-circle fa-stack-2x" style="color:#d90013"></span>
<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span> all the isolate's antigenic 
variants have been deemed not cross-reactive to vaccine variants through experimental studies.</li>
<li><span class="fa-stack" style="margin-left:0.5em;font-size:0.8em">
<span class="fas fa-circle fa-stack-2x" style="color:#d3d3d3"></span>
<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span> isolate contains antigens 
for which there is insufficient data from or are yet to be tested in experimental studies.</li>
</ul>

<p style="color:#606060"><strong>It is important to understand the caveats to interpreting the 
MenDeVAR Index:</strong></p>
<p><strong>Source of data</strong> - These data combine multiple sources of information including: 
peptide sequence identity through whole genome sequencing; experimental assays developed as indirect 
measures of the breadth of vaccine protection against diverse meningococci; and assays developed to 
assess immunogenicity. The meningococcal antigen surface expression (MEASURE)<sup>2</sup> and serum 
bactericidal activity (SBA) assays were used for Trumenba<sup>&reg;</sup>.</p>    
<p><strong>Cross-reactivity definition</strong> - An antigenic variant was considered cross-reactive 
if it had been tested in &ge;5 isolates/subjects and was above the accepted threshold in &ge;75% of 
those isolates. This was established through combined analysis of published experimental studies 
(PMID provided for each variant), not from genomic data. These assays were based on serogroup B disease 
isolates.</p>
<p><strong>Age of vaccinees</strong> - The age of vaccine recipients in the experimental studies varies 
widely, ranging from toddlers to adults, and needs to be taken into consideration when interpreting results. 
Vaccine studies used different schedules and doses of vaccines.    
<ol>
<li>Jiang HQ, Hoiseth SK, Harris SL, <i>et al.</i> Broad vaccine coverage predicted for a bivalent 
recombinant factor H binding protein based vaccine to prevent serogroup B meningococcal disease. 
<a href="https://www.ncbi.nlm.nih.gov/pubmed/20619376" target="_blank">Vaccine 2010;28(37):6086-93</a></li>
<li>McNeil LK, Donald RGK, Gribenko A, <i>et al.</i> Predicting the Susceptibility of Meningococcal Serogroup 
B Isolates to Bactericidal Antibodies Elicited by Bivalent rLP2086, a Novel Prophylactic Vaccine. 
<a href="https://www.ncbi.nlm.nih.gov/pubmed/29535195" target="_blank">mBio 2018;9(2):e00036-18</a></li>
</ol>
<p><em>MenDeVAR is described in <a href="https://doi.org/10.1101/2020.08.18.256834">Rodrigues <i>et al.</i> 2020, 
bioRxiv 2020.08.18.256834</a>.  Please contact <a href="mailto:charlene.rodrigues@zoo.ox.ac.uk">us</a> 
if you have queries.</em></p>
NOTES
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
			  ALLELE: foreach my $allele_id ( @{ $results->{'exact_matches'}->{ $nuc_loci{$locus} } } ) {
					my $flags = $script->{'datastore'}->run_query(
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

sub get_script_obj {
	my $script_obj = BIGSdb::Offline::Script->new(
		{
			config_dir       => CONFIG_DIR,
			lib_dir          => LIB_DIR,
			dbase_config_dir => DBASE_CONFIG_DIR,
			options          => { always_run => 1, logger => $logger },
			instance         => SEQDEF_DATABASE,
		}
	);
	return $script_obj;
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
	my $matched = is_not_cross_reactive( $results, $components, $notes );
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
	$matched = is_not_cross_reactive( $results, $components, $notes );
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
		my %variants = map { $_ => 1 } @$variants;
		if ( $variants{ $component->{'variant'} } ) {
			$matched = 1;
			my @note_evidence;
			my @evidence_list = uniq( sort { $a cmp $b } values %{ $component->{'references'} } );
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
	my $notes = [];
	my %assays = map { $_ => 1 } qw(MATS SBA MEASURE);
	foreach my $component (@$components) {
		my $match_variants =
		  ref $results->{ $component->{'locus'} }->{'variants'}
		  ? $results->{ $component->{'locus'} }->{'variants'}
		  : [];
		my %match_variants = map { $_ => 1 } @$match_variants;
		if ( $match_variants{ $component->{'variant'} } ) {
			$result = 1;
			my @note_evidence;
			my @evidence_list = uniq( sort { $a cmp $b } values %{ $component->{'references'} } );
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

sub get_traffic_light {
	my ($value) = @_;
	my $lights = {
		'exact match' => q(<a title="Exact match to at least one component antigen">)
		  . q(<span class="fa-stack" style="font-size:0.9em">)
		  . q(<span class="fas fa-circle fa-stack-2x" style="color:#009800"></span>)
		  . q(<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span></a>),
		'cross-reactive' => q(<a title="Predicted cross-reactive to at least one component antigen">)
		  . q(<span class="fa-stack" style="font-size:0.9em">)
		  . q(<span class="fas fa-circle fa-stack-2x" style="color:#ef780f">)
		  . q(</span><span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span></a>),
		'none' => q(<a title="All antigens have been shown to not be cross-reactive to vaccine variants">)
		  . q(<span class="fa-stack" style="font-size:0.9em">)
		  . q(<span class="fas fa-circle fa-stack-2x" style="color:#d90013"></span>)
		  . q(<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span></a>),
		'insufficient data' => q(<a title="Insufficient data about cross-reactivity of antigens to make assessment">)
		  . q(<span class="fa-stack" style="font-size:0.9em">)
		  . q(<span class="fas fa-circle fa-stack-2x" style="color:#d3d3d3"></span>)
		  . q(<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span></a>)
	};
	return $lights->{$value} // q();
}
