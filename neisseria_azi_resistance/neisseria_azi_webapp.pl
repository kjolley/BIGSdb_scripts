#!/usr/bin/env perl
#Perform Neisseria azithromycin resistance prediction.
#Written by Keith Jolley, 2026.
#Version:20260317
use strict;
use warnings;
use 5.010;
use JSON;
use Log::Log4perl qw(get_logger);
use constant {
	CONFIG_DIR        => '/etc/bigsdb',
	LIB_DIR           => '/usr/local/lib',
	MIN_ALIGNMENT     => 50,
	MIN_IDENTITY      => 70,
	PYTHON            => '/home/bigsdb/neisseria_azi_resistance/.venv/bin/python3',
	ASSIGNMENT_SCRIPT => '/home/bigsdb/neisseria_azi_resistance/neisseria_azi_resistance.py',
	DIR               => '/home/bigsdb/neisseria_azi_resistance/',
	TMP               => '/var/bigsdb/'
};
use lib (LIB_DIR);
use IPC::Open3;
use Symbol qw(gensym);
Log::Log4perl->init( CONFIG_DIR . '/logging.conf' );
my $logger = Log::Log4perl::get_logger('BIGSdb.Page');
my @loci   = qw(23S_rRNA pro_NEIS1635 NEIS1633);

main();

sub main {
	my $scan_json_results = $ARGV[0];
	$logger->logdie('No results file passed.')                    if !$scan_json_results;
	$logger->logdie("Results file $scan_json_results does not exist.") if !-e $scan_json_results;
	my $scan_json_ref = slurp($scan_json_results);
	if ( !$scan_json_ref ) {
		$logger->logdie('Script failed!');
		return;
	}
		say << "CSS";
<style>
p.result_message {
	line-height: 1.2em;
	font-size: 3em;
}
p.no_result,p.no_value {
	color: #66a;
	background: -webkit-linear-gradient(#aaa, #66a);
	-webkit-background-clip: text;
	background-clip: text;
	-webkit-text-fill-color: transparent	
}
p.sensitive {
	color: green;
	background: -webkit-linear-gradient(#aaa, #6a6);
	-webkit-background-clip: text;
	background-clip: text;
	-webkit-text-fill-color: transparent
}
p.resistant {
	color: red;
	background: -webkit-linear-gradient(#aaa, #a66);
	-webkit-background-clip: text;
	background-clip: text;
	-webkit-text-fill-color: transparent
}
</style>
CSS
	
	my $scan_results      = extract_results($scan_json_ref);
	my @missing;
	foreach my $locus (@loci){
		push @missing,$locus if !defined $scan_results->{$locus};
	}
	if (keys %$scan_results != 4){

		say q(<div class="box" id="resultspanel"><div class="scrollable">);
		
		my $traffic_light = get_traffic_light('insufficient data');
		
		say qq(<p style="float:left;margin-right:2em">$traffic_light</p>);
		local $" = q(, );
		say qq(<p class="result_message no_result">No result.</p>);
		say qq(<p>No matches for: @missing.</p>);
		say q(</div></div>);
		return;
	}
	my $json_upload  = encode_json([$scan_results]);
	my $isolate_file = make_isolate_file($json_upload);

	my %params = (
		'--rRNA_23S'     => DIR . '23S_rRNA.json',
		'--pro_NEIS1635' => DIR . 'pro_NEIS1635.json',
		'--NEIS1633'     => DIR . 'NEIS1633.json',
		'--isolates'     => $isolate_file
	);
	my $err_fh = gensym;
	my $results;
	eval {
		my $pid = open3( my $in_fh, my $out_fh, $err_fh, PYTHON, ASSIGNMENT_SCRIPT, %params );
		close $in_fh;
		local $/ = undef;
		my $out = <$out_fh>;
		my $err = <$err_fh>;
		waitpid( $pid, 0 );
		my $exit_code = $? >> 8;
		if ($exit_code) {
			$logger->error($err) if $err;
			exit(1);
		}
		$results = decode_json($out);
	};
	$logger->logdie($@) if $@;
	

	say q(<div class="box" id="resultspanel"><div class="scrollable">);
	if (!$results || !$results->[0] || !$results->[0]->{'prediction'}){
		my $traffic_light = get_traffic_light('insufficient data');
		say qq(<p style="float:left;margin-right:2em">$traffic_light</p>);
		local $" = q(, );
		say qq(<p class="result_message no_result">No result.</p>);
	} else {
		my $traffic_light = get_traffic_light($results->[0]->{'prediction'});
		say qq(<p style="float:left;margin-right:2em">$traffic_light</p>);
		if ($results->[0]->{'prediction'} eq 'R'){
			 say qq(<p class="result_message resistant">Resistant</p>);
		} elsif ($results->[0]->{'prediction'} eq 'S'){
			say qq(<p class="result_message sensitive">Sensitive</p>);
		} elsif ($results->[0]->{'prediction'} eq 'insufficient data'){
			say qq(<p class="result_message no_result">Insufficient data.</p>);
		}
	}
	say q(</div></div>);
	return;
}

sub make_isolate_file {
	my ($json) = @_;
	my $filename = TMP . 'BIGSdb_isolates_' . ( int( rand(9_000_000_000) ) + 1_000_000_000 ) . '.json';
	open( my $fh, '>', $filename ) or $logger->logdie( "Cannot open $filename for writing.");
	say $fh $json;
	close $fh;
	return $filename;
}

sub extract_results {
	my ($json_ref) = @_;
	my $processed  = {};
	my $results    = decode_json($$json_ref);
	my $matches    = $results->{'exact_matches'};
	if ($matches) {
		$processed->{'id'} = 1;
		foreach my $locus (@loci) {
			if ( $matches->{$locus} && defined $matches->{$locus}->[0] ) {
				$processed->{$locus} = $matches->{$locus}->[0];
			}
		}
	}
	return $processed;
}

sub get_traffic_light {
	my ($value) = @_;
	my $lights = {
		'S' => q(<a title="Predicted to be sensitive to azithromycin">)
		  . q(<span class="fa-stack" style="font-size:3em">)
		  . q(<span class="fas fa-circle fa-stack-2x" style="color:#009800"></span>)
		  . q(<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span></a>),
		'R' => q(<a title="Predicted to be resistant to azithromycin">)
		  . q(<span class="fa-stack" style="font-size:3em">)
		  . q(<span class="fas fa-circle fa-stack-2x" style="color:#d90013"></span>)
		  . q(<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span></a>),
		'insufficient data' => q(<a title="Insufficient data to make assessment">)
		  . q(<span class="fa-stack" style="font-size:3em">)
		  . q(<span class="fas fa-circle fa-stack-2x" style="color:#d3d3d3"></span>)
		  . q(<span class="fas fa-traffic-light fa-stack-1x fa-inverse"></span></span></a>)
	};
	return $lights->{$value} // q();
}

sub slurp {
	my ($file_path) = @_;
	open( my $fh, '<:raw', $file_path )
	  || $logger->logdie("Cannot open $file_path for reading.");
	my $contents = do { local $/ = undef; <$fh> };
	return \$contents;
}
