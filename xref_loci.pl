#!/usr/bin/perl -T
#Cross-reference loci in sequence definition database with those
#in annotated reference
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
	PASSWORD         => '',
	BLAST_PATH       => '/usr/bin'
};
#######End Local configuration###############################
use lib (LIB_DIR);
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
use POSIX;
use Bio::SeqIO;
use Bio::DB::GenBank;
use BIGSdb::Offline::Script;
use BIGSdb::Utils;
use Log::Log4perl qw(get_logger);

#Direct all library logging calls to screen
my $log_conf = q(
	log4perl.category.BIGSdb.Script        = INFO, Screen
	log4perl.category.BIGSdb.Dataconnector = WARN, Screen
	log4perl.category.BIGSdb.Datastore     = WARN, Screen
	log4perl.appender.Screen               = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr        = 1
    log4perl.appender.Screen.layout        = Log::Log4perl::Layout::SimpleLayout
);
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
( my $script_name = $0 ) =~ s/.*\///;
say "Command: $script_name @ARGV\n";
my %opts;
GetOptions(
	'A|alignment=i'       => \$opts{'A'},
	'B|identity=i'        => \$opts{'B'},
	'a|accession=s'       => \$opts{'a'},
	'd|database=s'        => \$opts{'d'},
	'h|help'              => \$opts{'h'},
	'f|file=s'            => \$opts{'f'},
	'l|loci=s'            => \$opts{'l'},
	'L|exclude_loci=s'    => \$opts{'L'},
	'R|locus_regex=s'     => \$opts{'R'},
	's|schemes=s'         => \$opts{'s'},
	't|threads=i'         => \$opts{'t'},
	'u|local_accession=s' => \$opts{'u'},
	'w|word_size=i'       => \$opts{'w'}
) or die("Error in command line arguments\n");

if ( $opts{'h'} ) {
	show_help();
	exit;
}
if ( !$opts{'d'} ) {
	say "Usage: xref_loci.pl --database <seqdef db config>";
	say "Help: xref_loci.pl --help";
	exit;
}
if ( !$opts{'a'} && !$opts{'u'} ) {
	say "You must specify an accession using --accession or --upload_accession.";
	exit(1);
}
my $EXIT = 0;
local @SIG{qw (INT TERM HUP)} = ( sub { $EXIT = 1 } ) x 3;    #Allow temp files to be cleaned on kill signals
$opts{'A'} //= 99;
$opts{'B'} //= 99;
my $accession_name = $opts{'a'} // $opts{'u'};
$accession_name =~ s/\..+$//;
my $seq_obj;

if ( $opts{'u'} ) {
	if ( !-e $opts{'u'} ) {
		die "File $opts{'u'} doesn't exist.\n";
	}
	eval {
		my $seqio_obj = Bio::SeqIO->new( -file => $opts{'u'} );
		$seq_obj = $seqio_obj->next_seq;
	};
	if ($@) {
		die "Invalid data in local annotation $opts{'u'}.\n";
	}
} else {
	my $seq_db = Bio::DB::GenBank->new;
	eval { $seq_obj = $seq_db->get_Seq_by_acc( $opts{'a'} ) };
	if ( !$seq_obj ) {

		#Why is this printed twice??
		die "No data returned for accession number $opts{'a'}.\n";
	}
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
		instance         => $opts{'d'},
		logger           => $logger
	}
);
main();

sub main {
	my $fasta = make_fasta_from_accession();
	system( BLAST_PATH . '/makeblastdb', ( -in => $fasta, -logfile => '/dev/null', -dbtype => 'nucl' ) );
	my $loci = $script->get_selected_loci;
	print_heading("Locally defined loci in $accession_name");
	say "Locus\t$accession_name";
	my $out_prefix = BIGSdb::Utils::get_random();
	my $outfile    = "$out_prefix.BLAST.out";
	my %args       = (
		-num_threads => $opts{'t'} // 1,
		-max_target_seqs => 10,
		-word_size       => $opts{'w'} // 15,
		-out             => $outfile,
		-outfmt          => 6,
		-dust            => 'no'
	);

	foreach my $locus (@$loci) {
		my $locus_fasta = make_fasta_from_locus($locus);
		system( BLAST_PATH . '/blastn', %args, -db => $fasta, -query => $locus_fasta );
		last if $EXIT;
		my $match_locus = parse_blast_local( $locus, $outfile );
		unlink $locus_fasta;
		say "$locus\t$match_locus";
	}
	$script->delete_temp_files("$fasta*");
	unlink $outfile;
	exit(1) if $EXIT;
	say "\n-- ";
	$fasta = make_fasta_from_database($loci);
	system( BLAST_PATH . '/makeblastdb', ( -in => $fasta, -logfile => '/dev/null', -dbtype => 'nucl' ) );
	print_heading("Loci from $accession_name defined locally");
	say "$accession_name\tLocus";
	my $cds = get_cds($seq_obj);

	foreach my $cds (@$cds) {
		my $locus_fasta = make_locus_fasta_from_accession($cds);
		system( BLAST_PATH . '/blastn', %args, -db => $fasta, -query => $locus_fasta );
		last if $EXIT;
		my $match_locus = parse_blast_accession( $cds, $outfile );
		unlink $locus_fasta;
		say "$cds->{'id'}\t$match_locus";
	}
	$script->delete_temp_files("$fasta*");
	unlink $outfile;
	return;
}

sub print_heading {
	my ($heading) = @_;
	say $heading;
	say '-' x length $heading;
	return;
}

sub parse_blast_local {
	my ( $locus, $outfile ) = @_;
	open( my $fh, '<', $outfile ) || $logger->error("Can't open $outfile file");
	my @matching_loci;
	my %seen;
	while ( my $line = <$fh> ) {
		my @data = split /\t/, $line;
		my ( $allele_id, $match_locus, $identity, $alignment ) = @data[ 0 .. 3 ];
		my $allele_seq = $script->{'datastore'}->get_sequence( $locus, $allele_id );
		my $allele_length = length $$allele_seq;
		next if !$allele_length;
		my $percent_alignment = ( $alignment / $allele_length ) * 100;
		if ( $percent_alignment >= $opts{'A'} && $identity >= $opts{'B'} ) {
			push @matching_loci, $match_locus if !$seen{$match_locus};
			$seen{$match_locus} = 1;
		}
	}
	close $fh;
	local $" = ',';
	return "@matching_loci";
}

sub parse_blast_accession {
	my ( $cds, $outfile ) = @_;
	open( my $fh, '<', $outfile ) || $logger->error("Can't open $outfile file");
	my @matching_loci;
	my %seen;
	while ( my $line = <$fh> ) {
		my @data = split /\t/, $line;
		my ( $match_allele, $identity, $alignment ) = @data[ 1 .. 3 ];
		my $allele_length     = length $cds->{'seq'};
		my $percent_alignment = ( $alignment / $allele_length ) * 100;
		( my $match_locus = $match_allele ) =~ s/::.*$//;
		if ( $percent_alignment >= $opts{'A'} && $identity >= $opts{'B'} ) {
			push @matching_loci, $match_locus if !$seen{$match_locus};
			$seen{$match_locus} = 1;
		}
	}
	close $fh;
	local $" = ',';
	return "@matching_loci";
}

sub make_fasta_from_database {
	my ($loci)   = @_;
	my $prefix   = BIGSdb::Utils::get_random();
	my $filename = "$prefix.fas";
	open( my $fh, '>', $filename ) || die "Can't create temporary file $filename";
	foreach my $locus (@$loci) {
		my $sequences =
		  $script->{'datastore'}
		  ->run_query( "SELECT allele_id,sequence FROM sequences WHERE locus=? AND locus IN (SELECT id FROM loci WHERE data_type='DNA')",
			$locus, { fetch => 'all_arrayref', slice => {}, cache => 'make_fasta_from_database' } );
		foreach (@$sequences) {
			say $fh ">$locus" . '::' . $_->{'allele_id'};
			say $fh $_->{'sequence'};
		}
	}
	close $fh;
	return $filename;
}

sub make_locus_fasta_from_accession {
	my ($cds)    = @_;
	my $prefix   = BIGSdb::Utils::get_random();
	my $filename = "$prefix . fas";
	open( my $fh, '>', $filename ) || die " Can't create temporary file $filename";
	say $fh ">$cds->{'id'}";
	say $fh $cds->{'seq'};
	close $fh;
	return $filename;
}

sub make_fasta_from_accession {
	my $cds      = get_cds($seq_obj);
	my $prefix   = BIGSdb::Utils::get_random();
	my $filename = "$prefix.fas";
	open( my $fh, '>', $filename ) || die " Can't create temporary file $filename";
	foreach (@$cds) {
		say $fh ">$_->{'id'}";
		say $fh $_->{'seq'};
	}
	close $fh;
	return $filename;
}

sub make_fasta_from_locus {
	my ($locus) = @_;
	my $sequences = $script->{'datastore'}->run_query( "SELECT allele_id,sequence FROM sequences WHERE locus = ? ORDER BY allele_id",
		$locus, { fetch => 'all_arrayref', slice => {} } );
	my $filename = "$locus.fas";
	open( my $fh, '>', $filename ) || die "Can't create $filename";
	foreach (@$sequences) {
		say $fh ">$_->{'allele_id'}";
		say $fh $_->{'sequence'};
	}
	close $fh;
	return $filename;
}

sub get_cds {
	my ($seq_obj) = @_;
	my $seqs = [];
	foreach my $cds ( $seq_obj->get_SeqFeatures ) {
		if ( $cds->primary_tag eq 'CDS' ) {
			my $locus;
			foreach (qw (locus_tag gene gene_synonym old_locus_tag)) {
				my @values = $cds->has_tag($_) ? $cds->get_tag_values($_) : ();
				foreach my $value (@values) {
					if ( !defined $locus ) {
						$locus = $value;
					}
				}
			}
			local $" = '|';
			push @$seqs, { id => $locus, seq => $cds->seq->seq };
		}
	}
	return $seqs;
}

sub show_help {
	my $termios = POSIX::Termios->new;
	$termios->getattr;
	my $ospeed = $termios->getospeed;
	my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
	my ( $norm, $bold, $under ) = map { $t->Tputs( $_, 1 ) } qw/me md us/;
	say << "HELP";
${bold}NAME$norm
    ${bold}xref_loci.pl$norm - Cross-reference loci in seqdef database against 
    annotated reference genome
    
${bold}SYNOPSIS$norm
    ${bold}xref_loci.pl --database$norm ${under}DATABASE$norm ${bold}--accession$norm ${under}ACCESSION$norm [${under}options$norm]
    
    Output is directed to STDOUT.  
    Direct to output file: 
    
    ${bold}xref_loci.pl --database$norm ${under}DATABASE$norm ${bold}--accession$norm ${under}ACCESSION$norm >  ${bold}${under}OUTPUT_FILE$norm

${bold}OPTIONS$norm

${bold}-a, --accession$norm ${under}ACCESSION$norm  
    Comma-separated list of Genbank accession ids.
    
${bold}-A, --alignment$norm ${under}INT$norm
    Percentage alignment (default: 99).

${bold}-B, --identity$norm ${under}INT$norm
    Percentage identity (default: 99).

${bold}-d, --database$norm ${under}DATABASE$norm  
	Sequence definition database configuration.
	
${bold}-l, --loci$norm ${under}LIST$norm
    Comma-separated list of loci to scan (ignored if -s used).

${bold}-L, --exclude_loci$norm ${under}LIST$norm
    Comma-separated list of loci to exclude

${bold}-R, --locus_regex$norm ${under}REGEX$norm
    Regex for locus names.
    
${bold}-s, --schemes$norm ${under}LIST$norm
    Comma-separated list of scheme loci to scan.
    
${bold}-t, --threads$norm ${under}INT$norm
    Number of BLASTN threads to use.
    
${bold}-u, --local_accession$norm ${under}FILE$norm
    Filename of local accession file.  
    
${bold}-w, --word_size$norm ${under}INT$norm
    BLASTN word size.   

HELP
	return;
}
