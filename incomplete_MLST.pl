#!/usr/bin/perl -T
#List isolates with incomplete profiles and predict clonal complex
#Written by Keith Jolley 2011-2015
use strict;
use warnings;
###########Local configuration################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	HOST             => 'localhost',
	PORT             => 5432,
	USER             => 'apache',
	PASSWORD         => undef
};
#######End Local configuration###############################
use lib (LIB_DIR);
use List::MoreUtils qw(any);
use Getopt::Std;
use BIGSdb::Offline::Script;
my %opts;
getopts( 'c:d:s:h', \%opts );

if ( $opts{'h'} ) {
	show_help();
	exit;
}
if ( !$opts{'d'} || !$opts{'c'} ) {
	print "\nUsage: incomplete_MLST.pl -d <database config> -c <complex_file> [-s <MLST scheme id>]\n\n";
	print "Help: incomplete_MLST.pl -h\n";
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
		instance         => $opts{'d'},
		writable         => 1
	}
);
my ( $complexes, $complex_names ) = read_complex_file( $opts{'c'} );
my $scheme_id = $opts{'s'} || 1;
my $loci = $script->{'datastore'}->get_scheme_loci($scheme_id );
my $sql = $script->{'db'}->prepare("SELECT id, $script->{'system'}->{'labelfield'} FROM $script->{'system'}->{'view'} ORDER BY id");
foreach my $i ( 1 .. @$loci ) {
	my $buffer;
	eval { $sql->execute };
	die $@ if $@;
	while ( my ( $id, $isolate ) = $sql->fetchrow_array ) {
		my $allele_ids = $script->{'datastore'}->get_all_allele_ids($id);
		my @profile;
		my $missing = 0;
		foreach (@$loci) {
			if (defined $allele_ids->{$_} && @{$allele_ids->{$_}} ){
				push @profile, $allele_ids->{$_}->[0];
			} else {
				push @profile,'-';
				$missing++;
			}
		}
		if ( $missing == $i ) {
			my $complex = '';
			if ( $i < 4 ) {
				$complex = predict_complex($script, $scheme_id, \@profile, $complexes, $complex_names);
			}
			local $" = "\t";
			$buffer .= "$id\t$isolate\t@profile\t$complex\n";
		}
	}
	if ($buffer) {
		local $" = "\t";
		print "$i " . ( $i == 1 ? 'locus' : 'loci' ) . " missing:\n\n";
		print "id\tisolate\t@$loci\tpredicted complex\n";
		print "$buffer\n\n";
	}
}

sub predict_complex {
	my ( $script, $scheme_id, $profile, $complexes, $complex_names ) = @_;
	my $scheme    = $script->{'datastore'}->get_scheme($scheme_id);
	my $max_match = 0;
	my $complex = '';
	foreach my $central_ST (@$complexes) {
		my $central_profile = $scheme->get_profile_by_primary_keys( [$central_ST] );
		my $match = 0;
		foreach my $i ( 0 .. @$loci ) {
			if ( defined $profile->[$i] && $profile->[$i] eq $central_profile->[$i] ) {
				$match++;
			}
		}
		if ( $match > 3 && $match > $max_match ) {
			$max_match = $match;
			$complex   = $complex_names->{$central_ST};
			chomp $complex;
		}
	}
	return $complex;
}

sub read_complex_file {
	my ($complex_file) = @_;
	my @complexes;
	my %complex_names;
	if ( -e $complex_file ) {
		open( my $complex_fh, '<', $complex_file );
		while ( my $line = <$complex_fh> ) {
			my ( $central_ST, $name ) = split /\t/, $line;
			$central_ST =~ s/\D//g;
			if ($central_ST) {
				push @complexes, $central_ST;
			}
			$complex_names{$central_ST} = $name;
		}
		close $complex_fh;
	} else {
		die "Complex list file '$complex_file' does not exist.\n";
	}
	return \@complexes, \%complex_names;
}

sub show_help {
	print << "HELP";

Usage incomplete_MLST.pl -d <database config> -c <complex file> [-s <MLST scheme id>]

Options
-------
-c <file>  Complex file
-d <name>  Database configuration name.
-h         This help page.
-s <id>    Scheme id

HELP
	return;
}
