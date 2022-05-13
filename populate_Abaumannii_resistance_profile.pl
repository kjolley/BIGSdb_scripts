#!/usr/bin/env perl
#Populate A. baumannii resistance profile field based on SIR values of
#indivual antibiotic fields
#Written by Julia Moreno Manjon and Keith Jolley 2022
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
	DBASE_CONFIG     => 'pubmlst_abaumannii_isolates',
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Constants qw(LOG_TO_SCREEN);

#Direct all library logging calls to screen
my $log_conf = LOG_TO_SCREEN;
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my $script = BIGSdb::Offline::Script->new(
	{
		config_dir       => CONFIG_DIR,
		lib_dir          => LIB_DIR,
		dbase_config_dir => DBASE_CONFIG_DIR,
		instance         => DBASE_CONFIG,
		logger           => $logger
	}
);
main();
undef $script;

sub main {
	my $records =
	  $script->{'datastore'}->run_query( 'SELECT * FROM isolates WHERE resistance_profile IS NULL ORDER BY id',
		undef, { fetch => 'all_arrayref', slice => {} } );
	foreach my $record (@$records) {
		my $antib_panel_count = antib_panel_count_formula($record);
		next if !$antib_panel_count;
		my $resistance_profile = get_resistance_profile($record);
		say "$record->{'id'}\t$resistance_profile\t$antib_panel_count";
		eval {
			$script->{'db'}->do( 'UPDATE isolates SET resistance_profile=? WHERE id=?',
				undef, $resistance_profile, $record->{'id'} );
		};
		if ($@) {
			say $@;
			$script->{'db'}->rollback;
			exit;
		}
		$script->{'db'}->commit;
	}
}

sub get_resistance_profile {
	my ($args)            = @_;
	my $antib_panel_count = antib_panel_count_formula($args);
	my $family_count      = family_count_formula($args);
	my $RI_count          = antib_RI_formula($args);
	my $S_count           = antib_S_formula($args);
	my $RIfam_count       = antib_RIfam_formula($args);
	my $SIR_count         = $RI_count + $S_count;
	if ( $antib_panel_count != $SIR_count ) {
		return
		  "Error: antib_panel_count ($antib_panel_count) different from SIR_count ($S_count+$RI_count = $SIR_count)";
	}
	if ( $antib_panel_count == 0 ) {
		return "Error: No antibiotic panel";
	}
	if ( $antib_panel_count != 0 ) {
		if ( $RI_count == 9 && $S_count == 0 ) {    #PDR: non-susceptible to all antimicrobial agents listed
			return "PDR";
		}
		if ( $RI_count == 0 && $S_count == 9 ) {    #S: susceptible to all antimicrobial agents listed
			return "sensitive";
		}
		if ( $RIfam_count < 3 ) {
			return "not MDR";
		}
		if ( $RIfam_count == 3 ) {                  #MDR: non-susceptible to >=1 agent in >=3 antimicrobial categories.
			return "MDR";
		}
		if ( $family_count - $RIfam_count <= 2 ) {    #XDR: non-susceptible to >=1 agent in all but <=2 categories
			return "XDR";
		} else {
			return "Error: unclassified";
		}
	} else {
		return "Error $!";
	}
}

# Antibiotic panel count
sub antib_panel_count_formula {
	my ($args) = @_;
	my (
		$amikacin_SIR,     $gentamicin_SIR, $cefotaxime_SIR, $cefepime_SIR, $levofloxacin_SIR,
		$tetracycline_SIR, $imipenem_SIR,   $meropenem_SIR,  $colistin_SIR
	  )
	  = @{$args}{
		qw(amikacin_sir gentamicin_sir cefotaxime_sir cefepime_sir levofloxacin_sir tetracycline_sir
		  imipenem_sir meropenem_sir colistin_sir)
	  };
	my $antib_panel_count = 0;
	foreach my $SIR (
		$amikacin_SIR,     $gentamicin_SIR, $cefotaxime_SIR, $cefepime_SIR, $levofloxacin_SIR,
		$tetracycline_SIR, $imipenem_SIR,   $meropenem_SIR,  $colistin_SIR
	  )
	{
		$antib_panel_count++ if defined $SIR;
	}
	return $antib_panel_count;
}

# Antibiotic family count
sub family_count_formula {
	my ($args) = @_;
	my (
		$amikacin_SIR,     $gentamicin_SIR, $cefotaxime_SIR, $cefepime_SIR, $levofloxacin_SIR,
		$tetracycline_SIR, $imipenem_SIR,   $meropenem_SIR,  $colistin_SIR
	  )
	  = @{$args}{
		qw(amikacin_sir gentamicin_sir cefotaxime_sir cefepime_sir levofloxacin_sir tetracycline_sir
		  imipenem_sir meropenem_sir colistin_sir)
	  };

	#Define variable as empty string if undefined.
	$_ //= '' foreach (
		$amikacin_SIR,     $gentamicin_SIR, $cefotaxime_SIR, $cefepime_SIR, $levofloxacin_SIR,
		$tetracycline_SIR, $imipenem_SIR,   $meropenem_SIR,  $colistin_SIR
	);
	my $family_count = 0;
	my %allowed = map { $_ => 1 } qw(S I R);
	if ( $allowed{$amikacin_SIR} or $allowed{$gentamicin_SIR} ) {    #Aminoglycosides
		$family_count++;
	}
	if ( $allowed{$cefotaxime_SIR} or $allowed{$cefepime_SIR} ) {    #Extended-spectrum cephalosporins
		$family_count++;
	}
	if ( $allowed{$levofloxacin_SIR} ) {                             #Antipseudomonal fluoroquinolones
		$family_count++;
	}
	if ( $allowed{$tetracycline_SIR} ) {                             #Tetracyclines
		$family_count++;
	}
	if ( $allowed{$imipenem_SIR} or $allowed{$meropenem_SIR} ) {     #Antipseudomonal carbapenems
		$family_count++;
	}
	if ( $allowed{$colistin_SIR} ) {                                 #Polymyxins
		$family_count++;
	}
	return $family_count;
}

# Count resistance (R/I)
sub antib_RI_formula {
	my ($args) = @_;
	my (
		$amikacin_SIR,     $gentamicin_SIR, $cefotaxime_SIR, $cefepime_SIR, $levofloxacin_SIR,
		$tetracycline_SIR, $imipenem_SIR,   $meropenem_SIR,  $colistin_SIR
	  )
	  = @{$args}{
		qw(amikacin_sir gentamicin_sir cefotaxime_sir cefepime_sir levofloxacin_sir tetracycline_sir
		  imipenem_sir meropenem_sir colistin_sir)
	  };

	#Define variable as empty string if undefined.
	$_ //= '' foreach (
		$amikacin_SIR,     $gentamicin_SIR, $cefotaxime_SIR, $cefepime_SIR, $levofloxacin_SIR,
		$tetracycline_SIR, $imipenem_SIR,   $meropenem_SIR,  $colistin_SIR
	);
	my $RI_count = 0;
	my %allowed = ( 'R' => 1, 'I' => 1 );
	if ( $allowed{$amikacin_SIR} ) {
		$RI_count++;
	}
	if ( $allowed{$gentamicin_SIR} ) {
		$RI_count++;
	}
	if ( $allowed{$cefotaxime_SIR} ) {
		$RI_count++;
	}
	if ( $allowed{$cefepime_SIR} ) {
		$RI_count++;
	}
	if ( $allowed{$levofloxacin_SIR} ) {
		$RI_count++;
	}
	if ( $allowed{$tetracycline_SIR} ) {
		$RI_count++;
	}
	if ( $allowed{$imipenem_SIR} ) {
		$RI_count++;
	}
	if ( $allowed{$meropenem_SIR} ) {
		$RI_count++;
	}
	if ( $allowed{$colistin_SIR} ) {
		$RI_count++;
	}
	return $RI_count;
}

# Count sensitive (S)
sub antib_S_formula {
	my ($args) = @_;
	my (
		$amikacin_SIR,     $gentamicin_SIR, $cefotaxime_SIR, $cefepime_SIR, $levofloxacin_SIR,
		$tetracycline_SIR, $imipenem_SIR,   $meropenem_SIR,  $colistin_SIR
	  )
	  = @{$args}{
		qw(amikacin_sir gentamicin_sir cefotaxime_sir cefepime_sir levofloxacin_sir tetracycline_sir
		  imipenem_sir meropenem_sir colistin_sir)
	  };

	#Define variable as empty string if undefined.
	$_ //= '' foreach (
		$amikacin_SIR,     $gentamicin_SIR, $cefotaxime_SIR, $cefepime_SIR, $levofloxacin_SIR,
		$tetracycline_SIR, $imipenem_SIR,   $meropenem_SIR,  $colistin_SIR
	);
	my $S_count = 0;
	my %allowed = ( 'S' => 1 );
	if ( $allowed{$amikacin_SIR} ) {
		$S_count++;
	}
	if ( $allowed{$gentamicin_SIR} ) {
		$S_count++;
	}
	if ( $allowed{$cefotaxime_SIR} ) {
		$S_count++;
	}
	if ( $allowed{$cefepime_SIR} ) {
		$S_count++;
	}
	if ( $allowed{$levofloxacin_SIR} ) {
		$S_count++;
	}
	if ( $allowed{$tetracycline_SIR} ) {
		$S_count++;
	}
	if ( $allowed{$imipenem_SIR} ) {
		$S_count++;
	}
	if ( $allowed{$meropenem_SIR} ) {
		$S_count++;
	}
	if ( $allowed{$colistin_SIR} ) {
		$S_count++;
	}
	return $S_count;
}

# Count resistant families (R/I)
sub antib_RIfam_formula {
	my ($args) = @_;
	my (
		$amikacin_SIR,     $gentamicin_SIR, $cefotaxime_SIR, $cefepime_SIR, $levofloxacin_SIR,
		$tetracycline_SIR, $imipenem_SIR,   $meropenem_SIR,  $colistin_SIR
	  )
	  = @{$args}{
		qw(amikacin_sir gentamicin_sir cefotaxime_sir cefepime_sir levofloxacin_sir tetracycline_sir
		  imipenem_sir meropenem_sir colistin_sir)
	  };

	#Define variable as empty string if undefined.
	$_ //= '' foreach (
		$amikacin_SIR,     $gentamicin_SIR, $cefotaxime_SIR, $cefepime_SIR, $levofloxacin_SIR,
		$tetracycline_SIR, $imipenem_SIR,   $meropenem_SIR,  $colistin_SIR
	);
	my $RIfam_count = 0;
	my %allowed = ( 'R' => 1, 'I' => 1 );
	if ( $allowed{$amikacin_SIR} or $allowed{$gentamicin_SIR} ) {    #Aminoglycosides
		$RIfam_count++;
	}
	if ( $allowed{$cefotaxime_SIR} or $allowed{$cefepime_SIR} ) {    #Extended-spectrum cephalosporins
		$RIfam_count++;
	}
	if ( $allowed{$levofloxacin_SIR} ) {                             #Antipseudomonal fluoroquinolones
		$RIfam_count++;
	}
	if ( $allowed{$tetracycline_SIR} ) {                             #Tetracyclines
		$RIfam_count++;
	}
	if ( $allowed{$imipenem_SIR} or $allowed{$meropenem_SIR} ) {     #Antipseudomonal carbapenems
		$RIfam_count++;
	}
	if ( $allowed{$colistin_SIR} ) {                                 #Polymyxins
		$RIfam_count++;
	}
	return $RIfam_count;
}
