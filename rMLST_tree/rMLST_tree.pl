#!/usr/bin/perl
#Generate trees from rMLST profile data.
#Written by Keith Jolley, 2017-2020
#Version 20201120
use strict;
use warnings;
use 5.010;
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	ISOLATE_DB       => 'pubmlst_rmlst_isolates',
	SEQDEF_DB        => 'pubmlst_rmlst_seqdef',
	RMLST_SCHEME_ID  => 1,
	TMP_DIR          => '/var/tmp',
	RAPIDNJ_PATH     => '/usr/local/bin/rapidnj'
};
#######End Local configuration###############################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Utils;
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
use File::Path qw(make_path);
use File::Copy;
use Bio::SeqIO;
use JSON;
use POSIX qw{strftime};
use constant RANKS => qw(genus family order class phylum);

#Direct all library logging calls to screen
my $log_conf =
    qq(log4perl.category.BIGSdb.Script        = INFO, Screen\n)
  . qq(log4perl.category.BIGSdb.Dataconnector = WARN, Screen\n)
  . qq(log4perl.category.BIGSdb.Datastore     = WARN, Screen\n)
  . qq(log4perl.appender.Screen               = Log::Log4perl::Appender::Screen\n)
  . qq(log4perl.appender.Screen.stderr        = 1\n)
  . qq(log4perl.appender.Screen.layout        = Log::Log4perl::Layout::SimpleLayout\n);
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my %opts;
GetOptions(
	'depth=i'           => \$opts{'depth'},
	'dir=s'             => \$opts{'dir'},
	'format=s'          => \$opts{'format'},
	'help'              => \$opts{'help'},
	'hyperlinks'        => \$opts{'hyperlinks'},
	'include_top_level' => \$opts{'include_top_level'},
	'isolate_count'     => \$opts{'isolate_count'},
	'method=s'          => \$opts{'method'},
	'new_only'          => \$opts{'new_only'},
	'public'            => \$opts{'public'},
	'quiet'             => \$opts{'quiet'},
	'rank=s'            => \$opts{'rank'},
	'rst_count'         => \$opts{'rst_count'},
	'taxon=s'           => \$opts{'taxon'},
	'threads=i'         => \$opts{'threads'},
	'trees'             => \$opts{'trees'}
) or die("Error in command line arguments.\n");
if ( $opts{'help'} ) {
	show_help();
	exit;
}
die "No method specified.\n" if !$opts{'method'};
check_valid_rank();
my $seqdef_db = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		instance         => SEQDEF_DB,
		options          => { always_run => 1 }
	}
);
my $isolate_db = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		instance         => ISOLATE_DB,
		options          => { always_run => 1 }
	}
);
$opts{'dir'}   //= '/var/tmp/taxonomy';
$opts{'depth'} //= 7;
$opts{'depth'}++;
my $methods = {
	list_rsts     => \&list_rsts,
	list_taxonomy => \&list_taxonomy,
	hierarchy     => \&show_hierarchy,
	trees         => \&trees
};
my %valid_method = map { $_ => 1 } keys %$methods;
die "Invalid method specified.\n" if !$valid_method{ $opts{'method'} };
$methods->{ $opts{'method'} }->();
undef $seqdef_db;
undef $isolate_db;

sub list_taxonomy {
	my $data  = get_taxonomy();
	my @ranks = RANKS;
	local $" = qq(\t);
	say qq(species\t@ranks);
	foreach my $species ( sort keys %$data ) {
		my @values = @{ $data->{$species} }{@ranks};
		$_ //= q() foreach @values;
		say qq($species\t@values);
	}
	return;
}

sub list_rsts {
	my $rsts = get_rsts();
	say $_ foreach @$rsts;
	return;
}

sub get_rsts {
	my ( $rank, $taxon ) = @_;
	my $data         = get_taxonomy( $rank, $taxon );
	my @species      = sort keys %$data;
	my $temp_table   = $seqdef_db->{'datastore'}->create_temp_list_table_from_array( 'text', \@species );
	my $scheme_cache = 'mv_scheme_' . RMLST_SCHEME_ID;
	my $list =
	  $seqdef_db->{'datastore'}->run_query(
		"SELECT rST FROM $scheme_cache c JOIN $temp_table l ON c.species=l.value ORDER BY CAST(c.rST AS integer)",
		undef, { fetch => 'col_arrayref' } );
	$seqdef_db->{'db'}->do("DROP TABLE $temp_table");
	return $list;
}

sub check_valid_rank {
	if ( $opts{'rank'} ) {
		my %defined_ranks = map { $_ => 1 } ( RANKS, 'domain', 'species' );
		die "$opts{'rank'} is not a valid rank.\n" if !$defined_ranks{ $opts{'rank'} };
		$opts{'taxon'} //= 'Bacteria' if $opts{'rank'} eq 'domain';
		die "No taxon provided.\n" if !$opts{'taxon'};
	}
	return;
}

sub get_taxonomy {
	my ( $rank, $taxon ) = @_;
	$rank  //= $opts{'rank'};
	$taxon //= $opts{'taxon'};
	my $qry =
	  q(SELECT attribute,field_value,value FROM isolate_value_extended_attributes WHERE isolate_field='species');
	if ( $opts{'public'} ) {
		$qry .= q( AND field_value IN (SELECT species FROM public));
	}
	my $data = $isolate_db->{'datastore'}->run_query( $qry, undef, { fetch => 'all_arrayref', slice => {} } );
	$isolate_db->{'db'}->commit;    #Prevent idle in transaction locks
	my $taxonomy = {};
	map { $taxonomy->{ $_->{'field_value'} }->{ $_->{'attribute'} } = $_->{'value'} } @$data;
	$rank //= 'species' if $taxon;
	if ($rank) {
		my $filtered = {};
		foreach my $species ( keys %$taxonomy ) {
			$taxonomy->{$species}->{'species'} = $species;
			next if ( $taxonomy->{$species}->{$rank} // q() ) ne $taxon && $rank ne 'domain';
			$filtered->{$species} = $taxonomy->{$species};
		}
		return $filtered;
	}
	return $taxonomy;
}

sub get_sub_rank {
	my ($rank) = @_;
	my $sub_rank = {
		domain => 'phylum',
		phylum => 'class',
		class  => 'order',
		order  => 'family',
		family => 'genus',
		genus  => 'species',
	};
	return $sub_rank->{$rank};
}

sub get_rst_counts {
	my ($taxonomy) = @_;
	$taxonomy //= get_taxonomy();
	my $rst_counts = {};
	if ( $opts{'rst_count'} ) {
		my $table = 'mv_scheme_' . RMLST_SCHEME_ID;
		my $rst_data =
		  $seqdef_db->{'datastore'}
		  ->run_query( "SELECT species,count(*) AS count FROM $table WHERE species IS NOT NULL GROUP BY species",
			undef, { fetch => 'all_arrayref', slice => {} } );
		$seqdef_db->{'db'}->commit;
		foreach my $values (@$rst_data) {
			foreach my $rank ( RANKS, 'species' ) {
				next if !$taxonomy->{ $values->{'species'} }->{$rank};
				$rst_counts->{$rank}->{ $taxonomy->{ $values->{'species'} }->{$rank} } += $values->{'count'} // 0;
			}
		}
	}
	return $rst_counts;
}

sub get_isolate_counts {
	my ($taxonomy) = @_;
	$taxonomy //= get_taxonomy();
	my $isolate_counts = {};
	if ( $opts{'isolate_count'} ) {
		my $rst_data = $isolate_db->{'datastore'}->run_query(
			"SELECT species,count(*) AS count FROM $isolate_db->{'system'}->{'view'} "
			  . 'WHERE species IS NOT NULL GROUP BY species',
			undef,
			{ fetch => 'all_arrayref', slice => {} }
		);
		$isolate_db->{'db'}->commit;
		foreach my $values (@$rst_data) {
			foreach my $rank ( RANKS, 'species' ) {
				next if !$taxonomy->{ $values->{'species'} }->{$rank};
				$isolate_counts->{$rank}->{ $taxonomy->{ $values->{'species'} }->{$rank} } += $values->{'count'} // 0;
			}
		}
	}
	return $isolate_counts;
}

sub hierarchy {
	my ( $selected_rank, $selected_taxon, $depth, $taxonomy ) = @_;
	$depth-- if $depth;
	$taxonomy //= get_taxonomy();
	my $hierarchy = {};
	$selected_rank  //= 'domain';
	$selected_taxon //= 'Bacteria';
	return $hierarchy if $selected_rank eq 'species';
	my $sub_rank = get_sub_rank($selected_rank);
	my %taxa;
	my $rst_counts     = get_rst_counts($taxonomy);
	my $isolate_counts = get_isolate_counts($taxonomy);

	if ($depth) {
		foreach my $species ( keys %$taxonomy ) {
			next if !$sub_rank;
			$taxonomy->{$species}->{'domain'} = 'Bacteria';
			next if !$taxonomy->{$species}->{$selected_rank};
			if ( $taxonomy->{$species}->{$selected_rank} eq $selected_taxon ) {
				if ( $taxonomy->{$species}->{$sub_rank} ) {
					$taxa{ $taxonomy->{$species}->{$sub_rank} } = 1;
				}
			}
		}
		foreach my $taxon ( sort keys %taxa ) {
			my $taxa = hierarchy( $sub_rank, $taxon, $depth, $taxonomy );
			$hierarchy->{$taxon}->{'rank'} = $sub_rank;
			if ( $opts{'rst_count'} ) {
				$hierarchy->{$taxon}->{'rSTs'} = $rst_counts->{$sub_rank}->{$taxon} // 0;
			}
			if ( $opts{'isolate_count'} ) {
				$hierarchy->{$taxon}->{'isolates'} = $isolate_counts->{$sub_rank}->{$taxon} // 0;
			}
			$hierarchy->{$taxon}->{'children'} = $taxa if keys %$taxa;
		}
	}
	return $hierarchy;
}

sub format_hierarchy {
	my ($hierarchy) = @_;
	$opts{'format'} //= 'text';
	my %allowed_formats = map { $_ => 1 } qw(text html json);
	if ( !$allowed_formats{ $opts{'format'} } ) {
		die "Invalid format selected.\n";
	}
	format_hierarchy_text($hierarchy) if $opts{'format'} eq 'text';
	format_hierarchy_html($hierarchy) if $opts{'format'} eq 'html';
	format_hierarchy_json($hierarchy) if $opts{'format'} eq 'json';
	return;
}

sub format_hierarchy_text {
	my ( $hierarchy, $depth ) = @_;
	$depth //= 0;
	foreach my $taxon ( sort keys %$hierarchy ) {
		print q( ) x ( 4 * $depth );
		my @values;
		push @values, qq(isolates:$hierarchy->{$taxon}->{'isolates'}) if $hierarchy->{$taxon}->{'isolates'};
		push @values, qq(rSTs:$hierarchy->{$taxon}->{'rSTs'})         if $hierarchy->{$taxon}->{'rSTs'};
		my $term = $taxon;
		local $" = q(; );
		$term .= qq( (@values)) if @values;
		say $term;

		if ( $hierarchy->{$taxon}->{'children'} ) {
			format_hierarchy_text( $hierarchy->{$taxon}->{'children'}, $depth + 1 );
		}
	}
	return;
}

sub format_hierarchy_html {
	my ( $hierarchy, $depth, $tree_file_path ) = @_;
	my $top_level = defined $depth ? 0 : 1;
	say q(<ul>) if $top_level;
	$depth          //= 0;
	$tree_file_path //= "$opts{'dir'}/trees";
	foreach my $taxon ( sort keys %$hierarchy ) {
		my $indent = 4 * $depth + 4;
		print q( ) x $indent;
		my @values;
		if ( $hierarchy->{$taxon}->{'isolates'} ) {
			if ( $opts{'hyperlinks'} ) {
				push @values, qq(isolates:<a data-i="1">$hierarchy->{$taxon}->{'isolates'}</a>);
			} else {
				push @values, qq(isolates:$hierarchy->{$taxon}->{'isolates'});
			}
		}
		my $local_tree_file_path = $tree_file_path;
		( my $cleaned_taxon = $taxon ) =~ s/\ /_/gx;
		$local_tree_file_path .= "/$cleaned_taxon";
		my $local_tree_file = "$local_tree_file_path.nwk";
		if ( $hierarchy->{$taxon}->{'rSTs'} ) {
			if ( -e $local_tree_file ) {
				push @values, qq(rSTs:$hierarchy->{$taxon}->{'rSTs'} <a data-t="1">tree</a>);
			} else {
				push @values, qq(rSTs:$hierarchy->{$taxon}->{'rSTs'});
			}
		}
		my $term = qq($taxon);
		local $" = q(; );
		$term .= qq( (@values)) if @values;
		print qq(<li data-rank="$hierarchy->{$taxon}->{'rank'}" data-taxon="$taxon">)
		  . qq(<span title="$hierarchy->{$taxon}->{'rank'}">$term</span>);
		my $closing_on_new_line = 0;
		if ( $hierarchy->{$taxon}->{'children'} ) {
			say q(<ul>);
			format_hierarchy_html( $hierarchy->{$taxon}->{'children'}, $depth + 1, $local_tree_file_path );
			print q( ) x ( $indent + 2 );
			say q(</ul>);
			$closing_on_new_line = 1;
		}
		print q( ) x $indent if $closing_on_new_line;
		say q(</li>);
	}
	say q(</ul>) if $top_level;
	return;
}

sub format_hierarchy_json {
	my ( $hierarchy, $depth, $tree_file_path ) = @_;
	my $top_level = defined $depth ? 0 : 1;
	say q([) if $top_level;
	$depth          //= 0;
	$tree_file_path //= "$opts{'dir'}/trees";
	my $first = 1;
	foreach my $taxon ( sort keys %$hierarchy ) {
		say q(,) if !$first;
		$first = 0;
		my $indent = 4 * $depth + 4;
		print q( ) x $indent;
		my @values = (
			qq("data":{"\$angularWidth":$hierarchy->{$taxon}->{'isolates'},"\$color":"#B2ABF4"}),
			qq("id":"$hierarchy->{$taxon}->{'rank'}/$taxon"),
			qq("name":"$taxon")
		);
		if ( $hierarchy->{$taxon}->{'isolates'} ) {
			push @values, qq("isolates":$hierarchy->{$taxon}->{'isolates'});
		}
		my $local_tree_file_path = $tree_file_path;
		( my $cleaned_taxon = $taxon ) =~ s/\ /_/gx;
		$local_tree_file_path .= "/$cleaned_taxon";
		my $local_tree_file = "$local_tree_file_path.nwk";
		if ( $hierarchy->{$taxon}->{'rSTs'} ) {
			push @values, qq("rSTs":$hierarchy->{$taxon}->{'rSTs'});
		}

		#		my $term = qq($taxon);
		local $" = q(,);
		print qq({@values);

		#		print qq(<li data-rank="$hierarchy->{$taxon}->{'rank'}" data-taxon="$taxon">)
		#		  . qq(<span title="$hierarchy->{$taxon}->{'rank'}">$term</span>);
		my $closing_on_new_line = 0;
		if ( $hierarchy->{$taxon}->{'children'} ) {
			say q(,"children":[);
			format_hierarchy_json( $hierarchy->{$taxon}->{'children'}, $depth + 1, $local_tree_file_path );
			print q( ) x ( $indent + 2 );
			say q(]);
			$closing_on_new_line = 1;
		}
		print q( ) x $indent if $closing_on_new_line;
		print q(});
	}
	say q(];) if $top_level;
	return;
}

sub show_hierarchy {
	my $hierarchy = hierarchy( $opts{'rank'}, $opts{'taxon'}, $opts{'depth'} );
	format_hierarchy($hierarchy);
	return;
}

sub trees {
	$opts{'rst_count'} = 1;
	undef $opts{'isolate_count'};
	my $hierarchy = hierarchy( $opts{'rank'}, $opts{'taxon'}, $opts{'depth'} );
	generate_trees( $hierarchy, $opts{'dir'} );
	return;
}

sub generate_trees {
	my ( $hierarchy, $dir, $depth ) = @_;
	$depth //= 0;
	if ( !-d $dir ) {
		make_path $dir or die "Failed to create path: $dir.\n";
	}
	make_tree( "$dir/$opts{'taxon'}.nwk", $opts{'rank'}, $opts{'taxon'} ) if $opts{'include_top_level'};
	my $parent_dir = $dir;
	foreach my $taxon ( sort keys %$hierarchy ) {
		next if ( $hierarchy->{$taxon}->{'rSTs'} // 0 ) < 3;
		my $filename = "$parent_dir/$taxon.nwk";
		make_tree( $filename, $hierarchy->{$taxon}->{'rank'}, $taxon );
		if ( $hierarchy->{$taxon}->{'children'} ) {
			$dir = "$parent_dir/$taxon";
			$dir =~ s/\ /_/gx;
			generate_trees( $hierarchy->{$taxon}->{'children'}, $dir, $depth + 1 );
		}
	}
	return;
}

sub make_tree {
	my ( $tree_file, $rank, $taxon ) = @_;
	$tree_file =~ s/\ /_/gx;
	print "Tree $tree_file..." if !$opts{'quiet'};
	if ( -e $tree_file && $opts{'new_only'} ) {
		say 'skipping.' if !$opts{'quiet'};
		return;
	}
	my $start_time   = time;
	my $rsts         = get_rsts( $rank, $taxon );
	my $loci         = $seqdef_db->{'datastore'}->get_scheme_loci(RMLST_SCHEME_ID);
	my $scheme_table = 'mv_scheme_' . RMLST_SCHEME_ID;
	my $start        = 1;
	my $end;
	my $no_output = 1;
	my $job_id    = BIGSdb::Utils::get_random();
	my $xmfa_file = TMP_DIR . "/${job_id}.xmfa";
	open( my $fh, '>', $xmfa_file )
	  or $logger->error("Cannot open output file $xmfa_file for writing");
	$isolate_db->{'dataConnector'}->drop_all_connections;
	$seqdef_db->{'dataConnector'}->drop_all_connections;
	foreach my $locus_name (@$loci) {
		my %no_seq;
		$seqdef_db->reconnect;
		my $locus_info   = $seqdef_db->{'datastore'}->get_locus_info($locus_name);
		my $locus        = $seqdef_db->{'datastore'}->get_locus($locus_name);
		my $temp         = BIGSdb::Utils::get_random();
		my $temp_file    = "$seqdef_db->{'config'}->{secure_tmp_dir}/$temp.txt";
		my $aligned_file = "$seqdef_db->{'config'}->{secure_tmp_dir}/$temp.aligned";
		open( my $fh_unaligned, '>', $temp_file ) or die("Could not open temp file $temp_file.\n");
		my $count = 0;

		foreach my $id (@$rsts) {
			$count++;
			my $profile_data = $seqdef_db->{'datastore'}->run_query( "SELECT * FROM $scheme_table WHERE rST=?",
				$id, { fetch => 'row_hashref', cache => 'profile_data' } );
			my $profile_id = $profile_data->{'rst'};
			my $header;
			if ( defined $profile_id ) {
				( my $species = $profile_data->{'species'} ) =~ s/\ /_/gx;
				$header = ">$profile_id|$species";
				my $allele_id =
				  $seqdef_db->{'datastore'}->get_profile_allele_designation( RMLST_SCHEME_ID, $id, $locus_name )
				  ->{'allele_id'};
				my $allele_seq_ref = $seqdef_db->{'datastore'}->get_sequence( $locus_name, $allele_id );
				say $fh_unaligned $header;
				if ( $allele_id eq '0' || $allele_id eq 'N' ) {
					say $fh_unaligned 'NNN';
					$no_seq{$id} = 1;
				} else {
					say $fh_unaligned $$allele_seq_ref;
				}
			} else {
				next;
			}
		}
		close $fh_unaligned;
		$seqdef_db->{'dataConnector'}->drop_all_connections;
		append_sequences(
			{
				fh                => $fh,
				output_locus_name => $locus_name,
				aligned_file      => $aligned_file,
				temp_file         => $temp_file,
				start             => \$start,
				end               => \$end,
				no_output_ref     => \$no_output,
				no_seq            => \%no_seq
			}
		);
	}
	close $fh;
	$seqdef_db->{'db'}->commit;    #Prevent idle in transaction locks
	if ( !-e $xmfa_file ) {
		say 'failed.' if !$opts{'quiet'};
		return;
	}
	my $fasta_file = BIGSdb::Utils::xmfa2fasta($xmfa_file);
	unlink $xmfa_file;
	my $output_tree_file = TMP_DIR . "/$job_id.tree";
	my $threads = $opts{'threads'} // 1;
	my $cmd =
	    RAPIDNJ_PATH
	  . " $fasta_file --input-format fa --cores $threads -x $output_tree_file "
	  . '--alignment-type d > /dev/null 2>&1';
	system $cmd;
	if ( !-e $output_tree_file ) {
		say 'failed (no tree produced).';
		return;
	}
	move( $output_tree_file, $tree_file ) || die "Copy failed.\n";
	my $duration = get_nice_duration( time - $start_time );
	say "done ($duration)." if !$opts{'quiet'};
	unlink $fasta_file;
	$isolate_db->reconnect;
	$seqdef_db->reconnect;
	return;
}

sub get_nice_duration {
	my ($seconds) = @_;
	return strftime( '%H:%M:%S', gmtime($seconds) );
}

sub append_sequences {
	my ($args) = @_;
	my ( $fh, $output_locus_name, $aligned_file, $temp_file, $start, $end, $no_output_ref, $no_seq ) =
	  @{$args}{qw(fh output_locus_name aligned_file temp_file start end no_output_ref no_seq)};
	my $output_file;
	if ( -e $temp_file && -s $temp_file ) {
		my $threads = $opts{'threads'} // 1;
		system(
			"$seqdef_db->{'config'}->{'mafft_path'} --thread $threads --quiet --preservecase $temp_file > $aligned_file"
		);
		$output_file = $aligned_file;
	} else {
		$output_file = $temp_file;
	}
	if ( -e $output_file && !-z $output_file ) {
		$$no_output_ref = 0;
		my $seq_in = Bio::SeqIO->new( -format => 'fasta', -file => $output_file );
		while ( my $seq = $seq_in->next_seq ) {
			my $length = $seq->length;
			$$end = $$start + $length - 1;
			say $fh '>' . $seq->id . ":$$start-$$end + $output_locus_name";
			my $sequence = BIGSdb::Utils::break_line( $seq->seq, 60 );
			( my $id = $seq->id ) =~ s/\|.*$//x;
			$sequence =~ s/N/-/g if $no_seq->{$id};
			say $fh $sequence;
		}
		$$start = ( $$end // 0 ) + 1;
		say $fh q(=);
	}
	unlink $output_file;
	unlink $temp_file;
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
    ${bold}rMLST_tree.pl$norm - Generate trees from rMLST profile data

${bold}SYNOPSIS$norm
    ${bold}rMLST_tree.pl$norm [${under}options$norm]

${bold}OPTIONS$norm

${bold}--depth$norm [${under}depth$norm]
    Number of sub-ranks to traverse.
    
${bold}--dir$norm [${under}directory$norm]
    Top-level output directory for tree files. Default '/var/tmp/taxonomy'.

${bold}--help$norm
    This help page.
    
${bold}--format$norm [${under}format$norm]
    Either text, html or json. Default 'text'.

${bold}--hyperlinks$norm
    Include hyperlinks to the isolate database when generating HTML taxonomic 
    output. 
    
${bold}--include_top_level$norm
    Include top-level tree (whole bacterial domain by default).
    
${bold}--isolate_count$norm
    Include isolate count in output.
    
${bold}--method$norm [${under}method$norm]
    hierarchy: Output hierarchical list - combine with --format, --rank and 
       --taxon.
    list_taxonomy: Display defined species taxonomic values.
    list_rsts: Display list of rSTs that match specified rank/taxon.
    trees: Generate trees at each level below specified rank/taxon.
    
${bold}--new_only$norm
    Only create a tree file if one does not already exist.
      
${bold}--public$norm
    Only include species found in the public view.
    
${bold}--quiet$norm
    Don't display progress messages. These are only created when generating
    trees.
      
${bold}--rank$norm [${under}taxonomic rank$norm]
    Either phylum, class, order, family, genus, or species. Default 'species'.
    This is used in combination with --taxon to specify a group of organisms.
    
${bold}--rst_count$norm
    Include rST count in output.
    
${bold}--taxon$norm [${under}taxonomic group$norm]
    Taxonomic group (of rank defined by --rank) or species.
    
${bold}--threads$norm [${under}threads$norm]
    Number of threads to use for MAFFT alignment.    
    
${bold}--trees$norm
    Make tree files. Trees are generated from concatenated sequences for rSTs
    only. There must >2 records in a category for a tree to be created.
        
HELP
	return;
}
