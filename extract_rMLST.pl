#!/usr/bin/perl -T
#Extract rMLST alleles for genomes with Genbank accessions
#Written by Keith Jolley, 2012-2015.
use strict;
use warnings;
use 5.010;
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	HOST             => 'zoo-oban.zoo.ox.ac.uk',
	PORT             => 5432,
	USER             => 'apache',
	PASSWORD         => 'remote',
	CACHE_DIR        => '/var/tmp/rMLST'
};
#######End Local configuration###############################
use lib (LIB_DIR);
use Getopt::Std;
use Error qw(:try);
use Bio::DB::GenBank;
use Bio::SeqIO;
use BIGSdb::Offline::Script;
use BIGSdb::Utils;
my %opts;
getopts( 'd:f:l:x:y:achs', \%opts );

if ( $opts{'h'} ) {
	show_help();
	exit;
}
if ( !$opts{'d'} || !$opts{'f'} ) {
	say "Usage: extract_rMLST.pl -d <database config> -f <data file>\n";
	say "Help: extract_rMLST.pl -h";
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
		instance         => $opts{'d'}
	}
);
die "Database $opts{'d'} does not exist!\n" if !$script->{'db'};
die "Must be a seqdef database!\n" if ( $script->{'system'}->{'dbtype'} // '' ) ne 'sequences';
my $data = read_data_file( $opts{'f'} );
my $loci =
  $script->{'datastore'}->run_query( "SELECT id FROM loci WHERE id IN (SELECT locus FROM scheme_members WHERE scheme_id=?) ORDER BY id",
	1, { fetch => 'col_arrayref' } );
if ( $opts{'s'} ) {
	say @$loci . " loci";
	say @$data . " records";
	exit;
}
my $retrieved = 0;
my $prefix    = BIGSdb::Utils::get_random();
my $i         = 0;
foreach my $record (@$data) {
	$i++;
	next if BIGSdb::Utils::is_int( $opts{'x'} ) && $i < $opts{'x'};
	my $accession_list = $record->{'accessions'};
	my $isolate_id     = ( $record->{'isolate_id'} =~ /(.*)/ ) ? $1 : undef;    #untaint
	my $seq_file       = CACHE_DIR . "/$isolate_id";
	if ( !$opts{'a'} && ( !-e $seq_file || !-s $seq_file ) ) {
		last if $opts{'l'} && BIGSdb::Utils::is_int( $opts{'l'} ) && $opts{'l'} <= $retrieved;
		sleep 60;                                                               #Let's not annoy NCBI.
		foreach my $accession (@$accession_list) {
			my $seq_obj;
			my $seq_db = Bio::DB::GenBank->new;
			$seq_db->verbose(2);                                                #convert warn to exception
			try {
				$seq_obj = $seq_db->get_Seq_by_acc($accession);
			}
			catch Bio::Root::Exception with {
				my $err = shift;
				die "No data returned for accession number $accession. $err\n";
			};
			open( my $acc_fh, '>>', $seq_file ) || die "Cannot open $seq_file for appending";
			say $acc_fh ">" . $seq_obj->id;
			say $acc_fh $seq_obj->seq;
			close $acc_fh;
			sleep 10;
		}
		$retrieved++;
	}
	next if $opts{'c'};
	next if $opts{'a'} && !-e $seq_file;
	foreach my $locus (@$loci) {
		my $ref_seq_file = create_locus_FASTA_db( $locus, $prefix );
		my $out_file = "$script->{'config'}->{'secure_tmp_dir'}/$prefix\_isolate_$isolate_id\_outfile.txt";
		blast( 15, $isolate_id, $ref_seq_file, $out_file );
		my $match = parse_blast( $locus, $out_file );
		my $extracted_seq = extract_sequence( $isolate_id, $match ) // '';
		$script->{'logger'}->error("$isolate_id $locus missing.") if !$extracted_seq;
		say "$isolate_id\t$locus\t$extracted_seq";
	}
	last if BIGSdb::Utils::is_int( $opts{'y'} ) && $i >= $opts{'y'};
}

sub read_data_file {
	my ($datafile) = @_;
	die "Data file $datafile does not exist!\n" if !-e $datafile;
	my @data;
	open( my $fh, '<', $datafile ) || die "Cannot open datafile $datafile";
	<$fh>;    #ignore header line
	while (<$fh>) {
		next if /^\s*$/;
		my ( $isolate_id, $acc, $genus_species ) = split /\t/, $_;
		my @accessions = split /\s+/, $acc;
		push @data, { isolate_id => $isolate_id, accessions => \@accessions, genus_species => $genus_species };
	}
	close $fh;
	return \@data;
}

sub parse_blast {

	#return best match
	my ( $locus, $blast_file ) = @_;
	my $identity  = 70;
	my $alignment = 50;
	my $match;
	my $quality = 0;    #simple metric of alignment length x percentage identity
	my $ref_seq_sql = $script->{'db'}->prepare("SELECT length(reference_sequence) FROM loci WHERE id=?");
	my %lengths;
	my @blast;
	open( my $blast_fh, '<', $blast_file ) || die "Can't open BLAST output file $blast_file. $!";
	push @blast, $_ foreach <$blast_fh>;    #slurp file to prevent file handle remaining open during database queries.
	close $blast_fh;

	foreach my $line (@blast) {
		next if !$line || $line =~ /^#/;
		my @record = split /\s+/, $line;
		if ( !$lengths{ $record[1] } ) {
			if ( $record[1] eq 'ref' ) {
				eval {
					$ref_seq_sql->execute($locus);
					( $lengths{ $record[1] } ) = $ref_seq_sql->fetchrow_array;
				};
				die "$@\n" if $@;
			} else {
				my $seq =
				  $script->{'datastore'}
				  ->run_query( "SELECT sequence FROM sequences WHERE locus=? AND allele_id=?", [$locus, $record[1]] );
				$lengths{ $record[1] } = length($seq);
			}
		}
		my $length       = $lengths{ $record[1] };
		my $this_quality = $record[3] * $record[2];
		if (   ( !$match->{'exact'} && $record[2] == 100 && $record[3] == $length )
			|| ( $this_quality > $quality && $record[3] > $alignment * 0.01 * $length && $record[2] >= $identity ) )
		{
			#Always score exact match higher than a longer partial match
			next if $match->{'exact'} && !( $record[2] == 100 && $record[3] == $length );
			$quality              = $this_quality;
			$match->{'seq_id'}    = $record[0];
			$match->{'allele'}    = $record[1];
			$match->{'identity'}  = $record[2];
			$match->{'length'}    = $length;
			$match->{'alignment'} = $record[3];
			$match->{'start'}     = $record[6];
			$match->{'end'}       = $record[7];

			if ( ( $record[8] > $record[9] && $record[7] > $record[6] ) || ( $record[8] < $record[9] && $record[7] < $record[6] ) ) {
				$match->{'reverse'} = 1;
			} else {
				$match->{'reverse'} = 0;
			}
			$match->{'exact'} = 1 if $match->{'identity'} == 100 && $match->{'alignment'} == $length;
			if ( $length > $match->{'alignment'} ) {
				if ( $match->{'reverse'} ) {
					if ( $record[8] < $record[9] ) {
						$match->{'predicted_start'} = $match->{'start'} - $length + $record[9];
						$match->{'predicted_end'}   = $match->{'end'} + $record[8] - 1;
					} else {
						$match->{'predicted_start'} = $match->{'start'} - $length + $record[8];
						$match->{'predicted_end'}   = $match->{'end'} + $record[9] - 1;
					}
				} else {
					if ( $record[8] < $record[9] ) {
						$match->{'predicted_start'} = $match->{'start'} - $record[8] + 1;
						$match->{'predicted_end'}   = $match->{'end'} + $length - $record[9];
					} else {
						$match->{'predicted_start'} = $match->{'start'} - $record[9] + 1;
						$match->{'predicted_end'}   = $match->{'end'} + $length - $record[8];
					}
				}
			} else {
				$match->{'predicted_start'} = $match->{'start'};
				$match->{'predicted_end'}   = $match->{'end'};
			}
		}
	}
	return $match;
}

sub blast {
	my ( $word_size, $isolate_id, $fasta_file, $out_file ) = @_;
	my $in_file = CACHE_DIR . "/$isolate_id";
	die "Isolate FASTA file does not exist for $isolate_id.\n" if !-e $fasta_file;
	system(
"$script->{'config'}->{'blast+_path'}/blastn -max_target_seqs 10 -parse_deflines -word_size $word_size -db $fasta_file -query $in_file -out $out_file -outfmt 6 -dust no"
	);
	return;
}

sub create_locus_FASTA_db {
	my ( $locus, $prefix ) = @_;
	my $clean_locus = $locus;
	$clean_locus =~ s/\W/_/g;
	$clean_locus = $1 if $locus =~ /(\w*)/;    #avoid taint check
	my $temp_fastafile = "$script->{'config'}->{'secure_tmp_dir'}/$prefix\_fastafile_$clean_locus.txt";
	$temp_fastafile =~ s/\\/\\\\/g;
	$temp_fastafile =~ s/'/__prime__/g;
	if ( !-e $temp_fastafile ) {
		my $locus_info = $script->{'datastore'}->get_locus_info($locus);
		my $file_buffer;
		my $seqs_ref =
		  $script->{'datastore'}
		  ->run_query( "SELECT allele_id,sequence FROM sequences WHERE locus=?", $locus, { fetch => 'all_arrayref', slice => {} } );
		foreach (@$seqs_ref) {
			next if !length $_->{'sequence'};
			$file_buffer .= ">$_->{'allele_id'}\n$_->{'sequence'}\n";
		}
		open( my $fasta_fh, '>', $temp_fastafile ) || die("Can't open $temp_fastafile for writing");
		print $fasta_fh $file_buffer if $file_buffer;
		close $fasta_fh;
		my $dbtype = $locus_info->{'data_type'} eq 'DNA' ? 'nucl' : 'prot';
		system("$script->{'config'}->{'blast+_path'}/makeblastdb -in $temp_fastafile -logfile /dev/null -parse_seqids -dbtype $dbtype");
	}
	return $temp_fastafile;
}

sub extract_sequence {
	my ( $isolate_id, $match ) = @_;
	my $start = $match->{'predicted_start'};
	my $end   = $match->{'predicted_end'};
	return if !defined $start || !defined $end;
	my $length = abs( $end - $start ) + 1;
	if ( $end < $start ) {
		$start = $end;
	}
	my $fasta_file = CACHE_DIR . "/$isolate_id";
	my $in = Bio::SeqIO->new( -file => $fasta_file, -format => 'fasta' );
	while ( my $seq_obj = $in->next_seq ) {
		my $id = $seq_obj->primary_id;
		next if $id ne $match->{'seq_id'};
		my $seq = $seq_obj->primary_seq->seq;
		my $extracted = substr( $seq, $start - 1, $length );
		$extracted = BIGSdb::Utils::reverse_complement($extracted) if $match->{'reverse'};
		return $extracted;
	}
	return '';
}

sub show_help {
	print << "HELP";

Usage extract_rMLST.pl -d <database config> -f <data file>

Options
-------
-a             Analyse only - only extract from genomes already downloaded
-c             Prepare genome cache
-d <name>      Database configuration name.
-f <filename>  Data filename.
-h             This help page.
-l <limit>     Number of genomes to retrieve to cache.
-s             Get basic stats.
-x             Start at row x
-y             Finish at row y

HELP
	return;
}
