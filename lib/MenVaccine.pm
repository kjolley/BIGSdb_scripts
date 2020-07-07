#Predict cross-reactivity to Bexsero and Trumenba vaccine components
#Written by Keith Jolley
#Copyright (c) 2020, University of Oxford
package MenVaccine;
use strict;
use warnings;
use 5.010;

sub get_exact_bexsero_match {
	return [
		{ locus => 'fHbp_peptide', variant => '1', references => { 27521232 => 'peptide sequence match' } },
		{ locus => 'NHBA_peptide', variant => '2', references => { 27521232 => 'peptide sequence match' } },
		{ locus => 'NadA_peptide', variant => '8', references => { 27521232 => 'peptide sequence match' } },
		{
			locus      => 'PorA_VR2',
			variant    => '4',
			references => { 27521232 => 'peptide sequence match' }
		}
	];
}

sub get_exact_trumenba_match {
	return [
		{ locus => 'fHbp_peptide', variant => '45', references => { 20619376 => 'peptide sequence match' } },
		{ locus => 'fHbp_peptide', variant => '55', references => { 20619376 => 'peptide sequence match' } },
	];
}

sub get_cross_reacts_with_bexsero {
	return [
		{
			locus      => 'fHbp_peptide',
			variant    => '4',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30135218 => 'MATS',
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '10',
			references => {
				25630407 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS',
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '12',
			references => {
				29152576 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '14',
			references => {
				23414709 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '15',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '37',
			references => {
				23414709 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '110',
			references => {
				26950303 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '144',
			references => {
				26686998 => 'MATS',
				26950303 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '215',
			references => {
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '232',
			references => {
				23414709 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '1',
			references => {
				27355628 => 'MATS',
				29152576 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '5',
			references => {
				29152576 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '10',
			references => {
				23414709 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '113',
			references => {
				23414709 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '243',
			references => {
				26686998 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '607',
			references => {
				30135218 => 'MATS'
			}
		},
		{
			locus      => 'NadA_peptide',
			variant    => '3',
			references => {
				23588089 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				29950334 => 'MATS'
			}
		},
		{
			locus      => 'NadA_peptide',
			variant    => '6',
			references => {
				29950334 => 'MATS'
			}
		}
	];
}

sub get_cross_reacts_with_trumenba {
	return [
		{
			locus      => 'fHbp_peptide',
			variant    => '1',
			references => {
				22569484 => 'SBA',
				22718089 => 'SBA',
				22871351 => 'SBA',
				23114369 => 'SBA',
				23352429 => 'SBA',
				26407272 => 'SBA',
				26707218 => 'SBA',
				26803328 => 'SBA',
				26835974 => 'SBA',
				26974889 => 'SBA',
				27745812 => 'SBA',
				27846061 => 'SBA',
				28196734 => 'SBA',
				28566335 => 'SBA',
				29236639 => 'SBA',
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '4',
			references => {
				27846061 => 'SBA',
				29236639 => 'SBA',
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '13',
			references => {
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '14',
			references => {
				22569484 => 'SBA',
				22718089 => 'SBA',
				27846061 => 'SBA',
				28196734 => 'SBA',
				29236639 => 'SBA',
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '15',
			references => {
				22569484 => 'SBA',
				22718089 => 'SBA',
				23352429 => 'SBA',
				26407272 => 'SBA',
				26707218 => 'SBA',
				26803328 => 'SBA',
				26835974 => 'SBA',
				27745812 => 'SBA',
				27846061 => 'SBA',
				28196734 => 'SBA',
				28566335 => 'SBA',
				29236639 => 'SBA',
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '16',
			references => {
				27846061 => 'SBA',
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '19',
			references => {
				22569484 => 'SBA',
				22718089 => 'SBA',
				22871351 => 'SBA',
				23114369 => 'SBA',
				23352429 => 'SBA',
				26407272 => 'SBA',
				26707218 => 'SBA',
				26803328 => 'SBA',
				26835974 => 'SBA',
				26974889 => 'SBA',
				27745812 => 'SBA',
				27846061 => 'SBA',
				27846061 => 'SBA',
				28196734 => 'SBA',
				28566335 => 'SBA',
				29236639 => 'SBA',
				29236639 => 'SBA',
				29535195 => 'MEASURE'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '21',
			references => {
				27846061 => 'SBA',
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '23',
			references => {
				28566335 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '25',
			references => {
				28566335 => 'SBA',
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '30',
			references => {
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '30',
			references => {
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '47',
			references => {
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '49',
			references => {
				22871351 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '76',
			references => {
				28566335 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '87',
			references => {
				22569484 => 'SBA',
				22718089 => 'SBA',
				22871351 => 'SBA',
				23114369 => 'SBA',
				23352429 => 'SBA',
				27846061 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '121',
			references => {
				28196734 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '180',
			references => {
				22569484 => 'SBA',
				27846061 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '187',
			references => {
				22569484 => 'SBA',
				26407272 => 'SBA',
				26707218 => 'SBA',
				26803328 => 'SBA',
				26835974 => 'SBA',
				27745812 => 'SBA',
				27846061 => 'SBA',
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '252',
			references => {
				29236639 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '276',
			references => {
				27846061 => 'SBA',
				28566335 => 'SBA'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '510',
			references => {
				28566335 => 'SBA'
			}
		}
	];
}

sub get_fhbp_no_reactivity_with_bexsero {
	[
		{
			locus      => 'fHbp_peptide',
			variant    => '16',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '19',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				25630407 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '21',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '22',
			references => {
				23414709 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '24',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				26950303 => 'MATS',
				27083425 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '25',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '29',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS',
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '30',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '31',
			references => {
				26950303 => 'MATS',
				27083425 => 'MATS',
				28366725 => 'MATS',
				29152576 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '45',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '47',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '59',
			references => {
				23414709 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '76',
			references => {
				29152576 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '109',
			references => {
				23414709 => 'MATS'
			}
		},
		{
			locus      => 'fHbp_peptide',
			variant    => '119',
			references => {
				23414709 => 'MATS',
				28366725 => 'MATS'
			}
		},
	];
}

sub get_nhba_no_reactivity_with_bexsero {
	return [
		{
			locus      => 'NHBA_peptide',
			variant    => '6',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '9',
			references => {
				26950303 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '17',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				26950303 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				30135218 => 'MATS',
				31770063 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '18',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				27355628 => 'MATS',
				28366725 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '25',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '30',
			references => {
				23414709 => 'MATS',
				25630407 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '31',
			references => {
				23414709 => 'MATS',
				27083425 => 'MATS',
				28366725 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '43',
			references => {
				23414709 => 'MATS',
				27355628 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '47',
			references => {
				23414709 => 'MATS',
				26686998 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '63',
			references => {
				26686998 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '112',
			references => {
				23414709 => 'MATS',
				23588089 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '120',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '160',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '187',
			references => {
				26950303 => 'MATS',
				28366725 => 'MATS'
			}
		},
		{
			locus      => 'NHBA_peptide',
			variant    => '197',
			references => {
				23414709 => 'MATS',
				26950303 => 'MATS'
			}
		},
	];
}

sub get_nadA_no_reactivity_with_bexsero {
	return [
		{
			locus      => 'NadA_peptide',
			variant    => '1',
			references => {
				26950303 => 'MATS',
				28366725 => 'MATS',
				29950334 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NadA_peptide',
			variant    => '21',
			references => {
				26950303 => 'MATS',
				30592763 => 'MATS'
			}
		},
		{
			locus      => 'NadA_peptide',
			variant    => '100',
			references => {
				29950334 => 'MATS',
				30135218 => 'MATS',
				30592763 => 'MATS'
			}
		},
	];
}
1;
