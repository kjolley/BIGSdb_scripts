#!/usr/bin/env perl
#Written by Keith Jolley
#Generate image file for Zooniverse from results of BIGSdb scan.
#Version: 20200218
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	WKHTMLIMAGE_PATH => '/usr/bin/wkhtmltoimage'
};
###########################################################################
use lib (LIB_DIR);
use BIGSdb::Offline::Scan;
use BIGSdb::ExtractedSequencePage;
use BIGSdb::Constants qw(LOG_TO_SCREEN);
use constant ALIGN_WIDTH => 100;
use CGI;
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;

#Direct all library logging calls to screen
my $log_conf = LOG_TO_SCREEN;
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my %opts;
GetOptions(
	'cds'                  => \$opts{'cds'},
	'd|database=s'         => \$opts{'d'},
	'dir=s'                => \$opts{'dir'},
	'f|flanking=i'         => \$opts{'flanking'},
	'help'                 => \$opts{'help'},
	'i|isolates=s'         => \$opts{'i'},
	'isolate_list_file=s'  => \$opts{'isolate_list_file'},
	'I|exclude_isolates=s' => \$opts{'I'},
	'l|loci=s'             => \$opts{'l'},
	'L|exclude_loci=s'     => \$opts{'L'},
	'm|min_size=i'         => \$opts{'m'},
	'p|projects=s'         => \$opts{'p'},
	'P|exclude_projects=s' => \$opts{'P'},
	'quiet'                => \$opts{'quiet'},
	'R|locus_regex=s'      => \$opts{'R'},
	's|schemes=s'          => \$opts{'s'},
	'w|word_size=i'        => \$opts{'w'},
	'x|min=i'              => \$opts{'x'},
	'y|max=i'              => \$opts{'y'},
	'v|view=s'             => \$opts{'v'}
) or die("Error in command line arguments\n");
if ( $opts{'help'} ) {
	show_help();
	exit;
}
if ( !$opts{'d'} ) {
	say "\nUsage: scan_image.pl --database <NAME>\n";
	say 'Help: scan_image.pl --help';
	exit;
}
$opts{'flanking'} //= 50;
$opts{'dir'}      //= '.';
$opts{'dir'} =~ s/\/$//x;
my $script = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		options          => { always_run => 1, %opts },
		instance         => $opts{'d'},
		logger           => $logger
	}
);
if ( ( $script->{'system'}->{'dbtype'} // q() ) ne 'isolates' ) {
	say q(Script can only be run against isolate databases.);
	undef $script;
	exit;
}
main();
undef $script;

sub main {
	my $isolates       = $script->get_isolates_with_linked_seqs;
	my $isolate_list   = $script->filter_and_sort_isolates($isolates);
	my $loci           = $script->get_loci_with_ref_db;
	my $isolate_prefix = BIGSdb::Utils::get_random();
	my $locus_prefix   = BIGSdb::Utils::get_random();
	my $params         = {
		identity  => $opts{'identity'}  // 70,
		alignment => $opts{'alignment'} // 50,
		word_size => $opts{'word_size'} // 20
	};
	my $scan = BIGSdb::Offline::Scan->new(
		{
			config_dir       => CONFIG_DIR,
			lib_dir          => LIB_DIR,
			dbase_config_dir => DBASE_CONFIG_DIR,
			options          => { always_run => 1, query_only => 1, %opts },
			instance         => $opts{'d'},
			logger           => $logger
		}
	);
	my $page = BIGSdb::ExtractedSequencePage->new(
		cgi           => CGI->new,
		system        => $script->{'system'},
		config        => $scan->{'config'},
		contigManager => $scan->{'contigManager'},
	);

	foreach my $isolate_id (@$isolates) {
		foreach my $locus (@$loci) {
			my ( $exact_matches, $partial_matches ) =
			  $scan->blast( $params, $locus, $isolate_id, "${isolate_prefix}_$isolate_id", $locus_prefix );
			next if ref $exact_matches && @$exact_matches;
			next if !ref $partial_matches || !@$partial_matches;
		  MATCH: foreach my $match (@$partial_matches) {
				next if $match->{'predicted_start'} < 1;
				my $seqbin_length = $script->{'datastore'}->run_query(
					'SELECT GREATEST(r.length,length(s.sequence)) FROM sequence_bin s LEFT JOIN '
					  . 'remote_contigs r ON s.id=r.seqbin_id WHERE s.id=?',
					$match->{'seqbin_id'},
					{ cache => 'contig_length' }
				);
				next if $match->{'predicted_end'} > $seqbin_length;
				my $seq_features = $page->get_seq_features(
					{
						seqbin_id => $match->{'seqbin_id'},
						reverse   => $match->{'reverse'},
						start     => $match->{'predicted_start'},
						end       => $match->{'predicted_end'},
						flanking  => 0
					}
				);
				foreach my $feature (@$seq_features) {
					if ( $feature->{'feature'} eq 'allele_seq' ) {
						my $seq = $feature->{'sequence'};
						next MATCH if $scan->is_complete_gene($seq) && !$opts{'cds'};
					}
				}
				my $html = get_html( $page, $match );
				my $html_file = create_html_file( $isolate_id, $locus, $html );
				( my $png_file = $html_file ) =~ s/\.html$/.png/x;
				my $program = WKHTMLIMAGE_PATH;
				system "xvfb-run -a $program --format jpg --quality 50 $html_file $png_file";
			}
			delete_temp_files($locus_prefix);
		}
		delete_temp_files($isolate_prefix);
	}
}

sub create_html_file {
	my ( $isolate_id, $locus, $html ) = @_;
	my $filename = "$opts{'dir'}/${isolate_id}_$locus.html";
	open( my $fh, '>', $filename ) || die "Cannot open $filename for writing.\n";
	say $fh '<!DOCTYPE html>';
	say $fh '<head>';
	say $fh '<style>';
	say $fh '.coding {background:#ccf}';
	say $fh '.startcodon {background:#cfc}';
	say $fh '.stopcodon {background:#fcc}';
	say $fh '</style>';
	say $fh '</head>';
	say $fh '<body>';
	say $fh '<pre>';
	say $fh $html;
	say $fh '</pre>';
	say $fh '</body>';
	close $fh;
	return $filename;
}

sub get_html {
	my ( $page, $match ) = @_;
	my $seq_features = $page->get_seq_features(
		{
			seqbin_id => $match->{'seqbin_id'},
			reverse   => $match->{'reverse'},
			start     => $match->{'predicted_start'},
			end       => $match->{'predicted_end'},
			flanking  => $opts{'flanking'}
		}
	);
	$page->{'prefs'}->{'alignwidth'} = ALIGN_WIDTH;
	return $page->get_sixpack_display($seq_features);
}

sub delete_temp_files {
	my ($prefix) = @_;
	my @files = glob("$script->{'config'}->{'secure_tmp_dir'}/*$prefix*");
	foreach (@files) { unlink $1 if /^(.*BIGSdb.*)$/x }
	return;
}

sub show_help {
	my $termios = POSIX::Termios->new;
	$termios->getattr;
	my $ospeed = $termios->getospeed;
	my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
	my ( $norm, $bold, $under ) = map { $t->Tputs( $_, 1 ) } qw(me md us);
	say << "HELP";
${bold}NAME$norm
    ${bold}scan_image.pl$norm - Generate sixpack image following BIGSdb scan.
    
${bold}SYNOPSIS$norm
    ${bold}scan_image.pl --database ${under}NAME$norm [${under}options$norm]

${bold}OPTIONS$norm

${bold}--cds$norm
    Produce output for predicted complete coding sequences. By default, output
    is only generated for sequences that the scan tool does not identify as
    complete CDS sequences.

${bold}--database$norm ${under}NAME$norm
    Database configuration name.
    
${bold}--dir$norm ${under}DIR$norm
    Output directory.
  
${bold}-f, --flanking$norm ${under}LENGTH$norm
    Length of flanking sequence to include (default: 50bp).  
    
${bold}--help$norm
    This help page.

${bold}-i, --isolates$norm ${under}LIST$norm  
    Comma-separated list of isolate ids to scan (ignored if -p used).
    
${bold}--isolate_list_file$norm ${under}FILE$norm  
    File containing list of isolate ids (ignored if -i or -p used).
           
${bold}-I, --exclude_isolates$norm ${under}LIST$norm
    Comma-separated list of isolate ids to ignore.

${bold}-l, --loci$norm ${under}LIST$norm
    Comma-separated list of loci to scan (ignored if -s used).

${bold}-L, --exclude_loci$norm ${under}LIST$norm
    Comma-separated list of loci to exclude

${bold}-m, --min_size$norm ${under}SIZE$norm
    Minimum size of seqbin (bp) - limit search to isolates with at least this
    much sequence.

${bold}-p, --projects$norm ${under}LIST$norm
    Comma-separated list of project isolates to scan.

${bold}-P, --exclude_projects$norm ${under}LIST$norm
    Comma-separated list of projects whose isolates will be excluded.
    
${bold}--quiet$norm
    Suppress normal output.
    
${bold}-R, --locus_regex$norm ${under}REGEX$norm
    Regex for locus names.

${bold}-s, --schemes$norm ${under}LIST$norm
    Comma-separated list of scheme loci to scan.
    
${bold}-v, --view$norm ${under}VIEW$norm
    Isolate database view (overrides value set in config.xml).

${bold}-w, --word_size$norm ${under}SIZE$norm
    BLASTN word size.

${bold}-x, --min$norm ${under}ID$norm
    Minimum isolate id.

${bold}-y, --max$norm ${under}ID$norm
    Maximum isolate id.
    
HELP
	return;
}
