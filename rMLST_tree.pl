#!/usr/bin/perl -T
#Generate trees from rMLST profile data.
#Written by Keith Jolley, 2017
use strict;
use warnings;
use 5.010;
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
###########Local configuration################################
use constant {
	CONFIG_DIR         => '/etc/bigsdb',
	LIB_DIR            => '/usr/local/lib',
	DBASE_CONFIG_DIR   => '/etc/bigsdb/dbases',
	ISOLATE_DB         => 'pubmlst_rmlst_isolates',
	SEQDEF_DB          => 'pubmlst_rmlst_seqdef',
	RMLST_SCHEME_CACHE => 'mv_scheme_1',
};
#######End Local configuration###############################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
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
	'depth=i'       => \$opts{'depth'},
	'format=s'      => \$opts{'format'},
	'help'          => \$opts{'help'},
	'hyperlinks'    => \$opts{'hyperlinks'},
	'isolate_count' => \$opts{'isolate_count'},
	'method=s'      => \$opts{'method'},
	'public'        => \$opts{'public'},
	'rank=s'        => \$opts{'rank'},
	'rst_count'     => \$opts{'rst_count'},
	'taxon=s'       => \$opts{'taxon'}
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
$opts{'depth'} //= 7;
$opts{'depth'}++;
my $methods = {
	list_rsts     => \&list_rsts,
	list_taxonomy => \&list_taxonomy,
	hierarchy     => sub {
		my $hierarchy = hierarchy( $opts{'rank'}, $opts{'taxon'}, $opts{'depth'} );
		format_hierarchy($hierarchy);
	  }
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
	my $data         = get_taxonomy();
	my @species      = sort keys %$data;
	my $temp_table   = $seqdef_db->{'datastore'}->create_temp_list_table_from_array( 'text', \@species );
	my $scheme_cache = RMLST_SCHEME_CACHE;
	return $seqdef_db->{'datastore'}->run_query(
		"SELECT rST FROM $scheme_cache c JOIN $temp_table l ON c.species=l.value ORDER BY CAST(c.rST AS integer)",
		undef, { fetch => 'col_arrayref' } );
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
	my $qry =
	  q(SELECT attribute,field_value,value FROM isolate_value_extended_attributes WHERE isolate_field='species');
	if ( $opts{'public'} ) {
		$qry .= q( AND field_value IN (SELECT species FROM public));
	}
	my $args     = [];
	my $data     = $isolate_db->{'datastore'}->run_query( $qry, $args, { fetch => 'all_arrayref', slice => {} } );
	my $taxonomy = {};
	map { $taxonomy->{ $_->{'field_value'} }->{ $_->{'attribute'} } = $_->{'value'} } @$data;
	$opts{'rank'} //= 'species' if $opts{'taxon'};
	if ( $opts{'rank'} ) {
		my $filtered = {};
		foreach my $species ( keys %$taxonomy ) {
			$taxonomy->{$species}->{'species'} = $species;
			next if ( $taxonomy->{$species}->{ $opts{'rank'} } // q() ) ne $opts{'taxon'} && $opts{'rank'} ne 'domain';
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
		my $table = RMLST_SCHEME_CACHE;
		my $rst_data =
		  $seqdef_db->{'datastore'}
		  ->run_query( "SELECT species,count(*) AS count FROM $table WHERE species IS NOT NULL GROUP BY species",
			undef, { fetch => 'all_arrayref', slice => {} } );
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
	my %allowed_formats = map { $_ => 1 } qw(text html);
	if ( !$allowed_formats{ $opts{'format'} } ) {
		die "Invalid format selected.\n";
	}
	format_hierarchy_text($hierarchy) if $opts{'format'} eq 'text';
	format_hierarchy_html($hierarchy) if $opts{'format'} eq 'html';
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
	my ( $hierarchy, $depth ) = @_;
	my $top_level = defined $depth ? 0 : 1;
	say q(<ul>) if $top_level;
	$depth //= 0;
	foreach my $taxon ( sort keys %$hierarchy ) {
		my $indent = 4 * $depth + 4;
		print q( ) x $indent;
		my @values;
		if ( $hierarchy->{$taxon}->{'isolates'} ) {
			if ( $opts{'hyperlinks'} ) {
				push @values, qq(isolates:<a data-rank="$hierarchy->{$taxon}->{'rank'}" data-taxon="$taxon">$hierarchy->{$taxon}->{'isolates'}</a>);
			} else {
				push @values, qq(isolates:$hierarchy->{$taxon}->{'isolates'});
			}
		}
		push @values, qq(rSTs:$hierarchy->{$taxon}->{'rSTs'}) if $hierarchy->{$taxon}->{'rSTs'};
		my $term = qq($taxon);
		local $" = q(; );
		$term .= qq( (@values)) if @values;
		print qq(<li><span title="$hierarchy->{$taxon}->{'rank'}">$term</span>);
		my $closing_on_new_line = 0;
		if ( $hierarchy->{$taxon}->{'children'} ) {
			say q(<ul>);
			format_hierarchy_html( $hierarchy->{$taxon}->{'children'}, $depth + 1 );
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

${bold}--depth [${under}depth$norm]
    Number of sub-ranks to show.

${bold}--help$norm
    This help page.
    
${bold}--format [${under}format$norm]
    Either text, html. Default 'text'.

${bold}--hyperlinks$norm
    Include hyperlinks to the isolate database when generating HTML taxonomic 
    output. 
    
${bold}--isolate_count
	Include isolate count in output.
      
${bold}--rank$norm [${under}taxonomic rank$norm]
    Either phylum, class, order, family, genus, or species. Default 'species'.
    This is used in combination with --taxon to specify a group of organisms.
    
${bold}--rst_count
	Include rST count in output.
    
${bold}--taxon$norm [${under}taxonomic group$norm]
    Taxonomic group (of rank defined by --rank) or species.
    
${bold}--method$norm [${under}method$norm]
    hierarchy: Output hierarchical list - combine with --format, --rank and 
       --taxon.
    list_taxonomy: Display defined species taxonomic values.
    list_rsts: Display list of rSTs that match specified rank/taxon.
    
${bold}--public$norm
    Only include species found in the public view.
        
HELP
	return;
}
