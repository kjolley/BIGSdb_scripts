#!/usr/bin/env perl
#Extract JSON dataset for visualisation from BIGSdb database
#records in BIGSdb Neisseria database
#Written by Keith Jolley
#Copyright (c) 2020-2022, University of Oxford
#E-mail: keith.jolley@zoo.ox.ac.uk
#
#This is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#BIGSdb is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this software.  If not, see <http://www.gnu.org/licenses/>.
#
#Version: 20220119
use strict;
use warnings;
use 5.010;
###########Local configuration#############################################
use constant {
	CONFIG_DIR       => '/etc/bigsdb',
	LIB_DIR          => '/usr/local/lib',
	DBASE_CONFIG_DIR => '/etc/bigsdb/dbases',
};
#######End Local configuration#############################################
use lib (LIB_DIR);
use BIGSdb::Offline::Script;
use BIGSdb::Constants qw(LOG_TO_SCREEN);
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
use JSON;
use constant COUNTRIES => {
	q(Afghanistan)                                  => { iso2 => q(AF), iso3 => q(AFG), continent => q(Asia) },
	q(Åland Islands)                               => { iso2 => q(AX), iso3 => q(ALA), continent => q(Europe) },
	q(Albania)                                      => { iso2 => q(AL), iso3 => q(ALB), continent => q(Europe) },
	q(Algeria)                                      => { iso2 => q(DZ), iso3 => q(DZA), continent => q(Africa) },
	q(American Samoa)                               => { iso2 => q(AS), iso3 => q(ASM), continent => q(Oceania) },
	q(Andorra)                                      => { iso2 => q(AD), iso3 => q(AND), continent => q(Europe) },
	q(Angola)                                       => { iso2 => q(AO), iso3 => q(AGO), continent => q(Africa) },
	q(Anguilla)                                     => { iso2 => q(AI), iso3 => q(AIA), continent => q(North America) },
	q(Antarctica)                                   => { iso2 => q(AQ), iso3 => q(ATA), continent => q(Antarctica) },
	q(Antigua and Barbuda)                          => { iso2 => q(AG), iso3 => q(ATG), continent => q(North America) },
	q(Argentina)                                    => { iso2 => q(AR), iso3 => q(ARG), continent => q(South America) },
	q(Armenia)                                      => { iso2 => q(AM), iso3 => q(ARM), continent => q(Asia) },
	q(Aruba)                                        => { iso2 => q(AW), iso3 => q(ABW), continent => q(North America) },
	q(Australia)                                    => { iso2 => q(AU), iso3 => q(AUS), continent => q(Oceania) },
	q(Austria)                                      => { iso2 => q(AT), iso3 => q(AUT), continent => q(Europe) },
	q(Azerbaijan)                                   => { iso2 => q(AZ), iso3 => q(AZE), continent => q(Europe) },
	q(Bahamas)                                      => { iso2 => q(BS), iso3 => q(BHS), continent => q(North America) },
	q(Bahrain)                                      => { iso2 => q(BH), iso3 => q(BHR), continent => q(Asia) },
	q(Bangladesh)                                   => { iso2 => q(BD), iso3 => q(BGD), continent => q(Asia) },
	q(Barbados)                                     => { iso2 => q(BB), iso3 => q(BRB), continent => q(North America) },
	q(Belarus)                                      => { iso2 => q(BY), iso3 => q(BLR), continent => q(Europe) },
	q(Belgium)                                      => { iso2 => q(BE), iso3 => q(BEL), continent => q(Europe) },
	q(Belize)                                       => { iso2 => q(BZ), iso3 => q(BLZ), continent => q(North America) },
	q(Benin)                                        => { iso2 => q(BJ), iso3 => q(BEN), continent => q(Africa) },
	q(Bermuda)                                      => { iso2 => q(BM), iso3 => q(BMU), continent => q(North America) },
	q(Bhutan)                                       => { iso2 => q(BT), iso3 => q(BTN), continent => q(Asia) },
	q(Bolivia)                                      => { iso2 => q(BO), iso3 => q(BOL), continent => q(South America) },
	q(Bonaire, Sint Eustatius and Saba)             => { iso2 => q(BQ), iso3 => q(BES), continent => q(North America) },
	q(Bosnia and Herzegovina)                       => { iso2 => q(BA), iso3 => q(BIH), continent => q(Europe) },
	q(Botswana)                                     => { iso2 => q(BW), iso3 => q(BWA), continent => q(Africa) },
	q(Bouvet Island)                                => { iso2 => q(BV), iso3 => q(BVT), continent => q(Antarctica) },
	q(Brazil)                                       => { iso2 => q(BR), iso3 => q(BRA), continent => q(South America) },
	q(British Indian Ocean Territory)               => { iso2 => q(IO), iso3 => q(IOT), continent => q(Asia) },
	q(British Virgin Islands)                       => { iso2 => q(VG), iso3 => q(VGB), continent => q(North America) },
	q(Brunei)                                       => { iso2 => q(BN), iso3 => q(BRN), continent => q(Asia) },
	q(Bulgaria)                                     => { iso2 => q(BG), iso3 => q(BGR), continent => q(Europe) },
	q(Burkina Faso)                                 => { iso2 => q(BF), iso3 => q(BFA), continent => q(Africa) },
	q(Burundi)                                      => { iso2 => q(BI), iso3 => q(BDI), continent => q(Africa) },
	q(Cambodia)                                     => { iso2 => q(KH), iso3 => q(KHM), continent => q(Asia) },
	q(Cameroon)                                     => { iso2 => q(CM), iso3 => q(CMR), continent => q(Africa) },
	q(Canada)                                       => { iso2 => q(CA), iso3 => q(CAN), continent => q(North America) },
	q(Cape Verde)                                   => { iso2 => q(CV), iso3 => q(CPV), continent => q(Africa) },
	q(Cayman Islands)                               => { iso2 => q(KY), iso3 => q(CYM), continent => q(North America) },
	q(Central African Republic)                     => { iso2 => q(CF), iso3 => q(CAF), continent => q(Africa) },
	q(Chad)                                         => { iso2 => q(TD), iso3 => q(TCD), continent => q(Africa) },
	q(Chile)                                        => { iso2 => q(CL), iso3 => q(CHL), continent => q(South America) },
	q(China)                                        => { iso2 => q(CN), iso3 => q(CHN), continent => q(Asia) },
	q(China [Hong Kong])                            => { iso2 => q(HK), iso3 => q(HKG), continent => q(Asia) },
	q(China [Macao])                                => { iso2 => q(MO), iso3 => q(MAC), continent => q(Asia) },
	q(Christmas Island)                             => { iso2 => q(CX), iso3 => q(CXR), continent => q(Asia) },
	q(Cocos (Keeling) Islands)                      => { iso2 => q(CC), iso3 => q(CCK), continent => q(Asia) },
	q(Colombia)                                     => { iso2 => q(CO), iso3 => q(COL), continent => q(South America) },
	q(Comoros)                                      => { iso2 => q(KM), iso3 => q(COM), continent => q(Africa) },
	q(Congo [DRC])                                  => { iso2 => q(CD), iso3 => q(COD), continent => q(Africa) },
	q(Congo [Republic])                             => { iso2 => q(CG), iso3 => q(COG), continent => q(Africa) },
	q(Cook Islands)                                 => { iso2 => q(CK), iso3 => q(COK), continent => q(Oceania) },
	q(Costa Rica)                                   => { iso2 => q(CR), iso3 => q(CRI), continent => q(North America) },
	q(Croatia)                                      => { iso2 => q(HR), iso3 => q(HRV), continent => q(Europe) },
	q(Cuba)                                         => { iso2 => q(CU), iso3 => q(CUB), continent => q(North America) },
	q(Curaçao)                                     => { iso2 => q(CW), iso3 => q(CUW), continent => q(North America) },
	q(Cyprus)                                       => { iso2 => q(CY), iso3 => q(CYP), continent => q(Europe) },
	q(Czech Republic)                               => { iso2 => q(CZ), iso3 => q(CZE), continent => q(Europe) },
	q(Denmark)                                      => { iso2 => q(DK), iso3 => q(DNK), continent => q(Europe) },
	q(Djibouti)                                     => { iso2 => q(DJ), iso3 => q(DJI), continent => q(Africa) },
	q(Dominica)                                     => { iso2 => q(DM), iso3 => q(DMA), continent => q(North America) },
	q(Dominican Republic)                           => { iso2 => q(DO), iso3 => q(DOM), continent => q(North America) },
	q(East Timor)                                   => { iso2 => q(TL), iso3 => q(TLS), continent => q(Asia) },
	q(Ecuador)                                      => { iso2 => q(EC), iso3 => q(ECU), continent => q(South America) },
	q(Egypt)                                        => { iso2 => q(EG), iso3 => q(EGY), continent => q(Africa) },
	q(El Salvador)                                  => { iso2 => q(SV), iso3 => q(SLV), continent => q(North America) },
	q(Equatorial Guinea)                            => { iso2 => q(GQ), iso3 => q(GNQ), continent => q(Africa) },
	q(Eritrea)                                      => { iso2 => q(ER), iso3 => q(ERI), continent => q(Africa) },
	q(Estonia)                                      => { iso2 => q(EE), iso3 => q(EST), continent => q(Europe) },
	q(Ethiopia)                                     => { iso2 => q(ET), iso3 => q(ETH), continent => q(Africa) },
	q(Falkland Islands (Malvinas))                  => { iso2 => q(FK), iso3 => q(FLK), continent => q(South America) },
	q(Faroe Islands)                                => { iso2 => q(FO), iso3 => q(FRO), continent => q(Europe) },
	q(Fiji)                                         => { iso2 => q(FJ), iso3 => q(FJI), continent => q(Oceania) },
	q(Finland)                                      => { iso2 => q(FI), iso3 => q(FIN), continent => q(Europe) },
	q(France)                                       => { iso2 => q(FR), iso3 => q(FRA), continent => q(Europe) },
	q(French Guiana)                                => { iso2 => q(GF), iso3 => q(GUF), continent => q(South America) },
	q(French Polynesia)                             => { iso2 => q(PF), iso3 => q(PYF), continent => q(Oceania) },
	q(French Southern Territories)                  => { iso2 => q(TF), iso3 => q(ATF), continent => q(Antarctica) },
	q(Gabon)                                        => { iso2 => q(GA), iso3 => q(GAB), continent => q(Africa) },
	q(Georgia)                                      => { iso2 => q(GE), iso3 => q(GEO), continent => q(Europe) },
	q(Germany)                                      => { iso2 => q(DE), iso3 => q(DEU), continent => q(Europe) },
	q(Ghana)                                        => { iso2 => q(GH), iso3 => q(GHA), continent => q(Africa) },
	q(Gibraltar)                                    => { iso2 => q(GI), iso3 => q(GIB), continent => q(Europe) },
	q(Greece)                                       => { iso2 => q(GR), iso3 => q(GRC), continent => q(Europe) },
	q(Greenland)                                    => { iso2 => q(GL), iso3 => q(GRL), continent => q(North America) },
	q(Grenada)                                      => { iso2 => q(GD), iso3 => q(GRD), continent => q(North America) },
	q(Guadeloupe)                                   => { iso2 => q(GP), iso3 => q(GLP), continent => q(North America) },
	q(Guam)                                         => { iso2 => q(GU), iso3 => q(GUM), continent => q(Oceania) },
	q(Guatemala)                                    => { iso2 => q(GT), iso3 => q(GTM), continent => q(North America) },
	q(Guernsey)                                     => { iso2 => q(GG), iso3 => q(GGY), continent => q(Europe) },
	q(Guinea)                                       => { iso2 => q(GN), iso3 => q(GIN), continent => q(Africa) },
	q(Guinea-Bissau)                                => { iso2 => q(GW), iso3 => q(GNB), continent => q(Africa) },
	q(Guyana)                                       => { iso2 => q(GY), iso3 => q(GUY), continent => q(South America) },
	q(Haiti)                                        => { iso2 => q(HT), iso3 => q(HTI), continent => q(North America) },
	q(Heard Island and McDonald Islands)            => { iso2 => q(HM), iso3 => q(HMD), continent => q(Antarctica) },
	q(Holy See)                                     => { iso2 => q(VA), iso3 => q(VAT), continent => q(Europe) },
	q(Honduras)                                     => { iso2 => q(HN), iso3 => q(HND), continent => q(North America) },
	q(Hungary)                                      => { iso2 => q(HU), iso3 => q(HUN), continent => q(Europe) },
	q(Iceland)                                      => { iso2 => q(IS), iso3 => q(ISL), continent => q(Europe) },
	q(India)                                        => { iso2 => q(IN), iso3 => q(IND), continent => q(Asia) },
	q(Indonesia)                                    => { iso2 => q(ID), iso3 => q(IDN), continent => q(Asia) },
	q(Iran)                                         => { iso2 => q(IR), iso3 => q(IRN), continent => q(Asia) },
	q(Iraq)                                         => { iso2 => q(IQ), iso3 => q(IRQ), continent => q(Asia) },
	q(Ireland)                                      => { iso2 => q(IE), iso3 => q(IRL), continent => q(Europe) },
	q(Isle of Man)                                  => { iso2 => q(IM), iso3 => q(IMN), continent => q(Europe) },
	q(Israel)                                       => { iso2 => q(IL), iso3 => q(ISR), continent => q(Asia) },
	q(Italy)                                        => { iso2 => q(IT), iso3 => q(ITA), continent => q(Europe) },
	q(Ivory Coast)                                  => { iso2 => q(CI), iso3 => q(CIV), continent => q(Africa) },
	q(Jamaica)                                      => { iso2 => q(JM), iso3 => q(JAM), continent => q(North America) },
	q(Japan)                                        => { iso2 => q(JP), iso3 => q(JPN), continent => q(Asia) },
	q(Jersey)                                       => { iso2 => q(JE), iso3 => q(JEY), continent => q(Europe) },
	q(Jordan)                                       => { iso2 => q(JO), iso3 => q(JOR), continent => q(Asia) },
	q(Kazakhstan)                                   => { iso2 => q(KZ), iso3 => q(KAZ), continent => q(Asia) },
	q(Kenya)                                        => { iso2 => q(KE), iso3 => q(KEN), continent => q(Africa) },
	q(Kiribati)                                     => { iso2 => q(KI), iso3 => q(KIR), continent => q(Oceania) },
	q(Kuwait)                                       => { iso2 => q(KW), iso3 => q(KWT), continent => q(Asia) },
	q(Kyrgyzstan)                                   => { iso2 => q(KG), iso3 => q(KGZ), continent => q(Asia) },
	q(Laos)                                         => { iso2 => q(LA), iso3 => q(LAO), continent => q(Asia) },
	q(Latvia)                                       => { iso2 => q(LV), iso3 => q(LVA), continent => q(Europe) },
	q(Lebanon)                                      => { iso2 => q(LB), iso3 => q(LBN), continent => q(Asia) },
	q(Lesotho)                                      => { iso2 => q(LS), iso3 => q(LSO), continent => q(Africa) },
	q(Liberia)                                      => { iso2 => q(LR), iso3 => q(LBR), continent => q(Africa) },
	q(Libya)                                        => { iso2 => q(LY), iso3 => q(LBY), continent => q(Africa) },
	q(Liechtenstein)                                => { iso2 => q(LI), iso3 => q(LIE), continent => q(Europe) },
	q(Lithuania)                                    => { iso2 => q(LT), iso3 => q(LTU), continent => q(Europe) },
	q(Luxembourg)                                   => { iso2 => q(LU), iso3 => q(LUX), continent => q(Europe) },
	q(Madagascar)                                   => { iso2 => q(MG), iso3 => q(MDG), continent => q(Africa) },
	q(Malawi)                                       => { iso2 => q(MW), iso3 => q(MWI), continent => q(Africa) },
	q(Malaysia)                                     => { iso2 => q(MY), iso3 => q(MYS), continent => q(Asia) },
	q(Maldives)                                     => { iso2 => q(MV), iso3 => q(MDV), continent => q(Asia) },
	q(Mali)                                         => { iso2 => q(ML), iso3 => q(MLI), continent => q(Africa) },
	q(Malta)                                        => { iso2 => q(MT), iso3 => q(MLT), continent => q(Europe) },
	q(Marshall Islands)                             => { iso2 => q(MH), iso3 => q(MHL), continent => q(Oceania) },
	q(Martinique)                                   => { iso2 => q(MQ), iso3 => q(MTQ), continent => q(North America) },
	q(Mauritania)                                   => { iso2 => q(MR), iso3 => q(MRT), continent => q(Africa) },
	q(Mauritius)                                    => { iso2 => q(MU), iso3 => q(MUS), continent => q(Africa) },
	q(Mayotte)                                      => { iso2 => q(YT), iso3 => q(MYT), continent => q(Africa) },
	q(Mexico)                                       => { iso2 => q(MX), iso3 => q(MEX), continent => q(North America) },
	q(Micronesia)                                   => { iso2 => q(FM), iso3 => q(FSM), continent => q(Oceania) },
	q(Moldova)                                      => { iso2 => q(MD), iso3 => q(MDA), continent => q(Europe) },
	q(Monaco)                                       => { iso2 => q(MC), iso3 => q(MCO), continent => q(Europe) },
	q(Mongolia)                                     => { iso2 => q(MN), iso3 => q(MNG), continent => q(Asia) },
	q(Montenegro)                                   => { iso2 => q(ME), iso3 => q(MNE), continent => q(Europe) },
	q(Montserrat)                                   => { iso2 => q(MS), iso3 => q(MSR), continent => q(North America) },
	q(Morocco)                                      => { iso2 => q(MA), iso3 => q(MAR), continent => q(Africa) },
	q(Mozambique)                                   => { iso2 => q(MZ), iso3 => q(MOZ), continent => q(Africa) },
	q(Myanmar)                                      => { iso2 => q(MM), iso3 => q(MMR), continent => q(Asia) },
	q(Namibia)                                      => { iso2 => q(NA), iso3 => q(NAM), continent => q(Africa) },
	q(Nauru)                                        => { iso2 => q(NR), iso3 => q(NRU), continent => q(Oceania) },
	q(Nepal)                                        => { iso2 => q(NP), iso3 => q(NPL), continent => q(Asia) },
	q(New Caledonia)                                => { iso2 => q(NC), iso3 => q(NCL), continent => q(Oceania) },
	q(New Zealand)                                  => { iso2 => q(NZ), iso3 => q(NZL), continent => q(Oceania) },
	q(Nicaragua)                                    => { iso2 => q(NI), iso3 => q(NIC), continent => q(North America) },
	q(Niger)                                        => { iso2 => q(NE), iso3 => q(NER), continent => q(Africa) },
	q(Nigeria)                                      => { iso2 => q(NG), iso3 => q(NGA), continent => q(Africa) },
	q(Niue)                                         => { iso2 => q(NU), iso3 => q(NIU), continent => q(Oceania) },
	q(Norfolk Island)                               => { iso2 => q(NF), iso3 => q(NFK), continent => q(Oceania) },
	q(North Korea)                                  => { iso2 => q(KP), iso3 => q(PRK), continent => q(Asia) },
	q(North Macedonia)                              => { iso2 => q(MK), iso3 => q(MKD), continent => q(Europe) },
	q(Northern Mariana Islands)                     => { iso2 => q(MP), iso3 => q(MNP), continent => q(Oceania) },
	q(Norway)                                       => { iso2 => q(NO), iso3 => q(NOR), continent => q(Europe) },
	q(Oman)                                         => { iso2 => q(OM), iso3 => q(OMN), continent => q(Asia) },
	q(Pakistan)                                     => { iso2 => q(PK), iso3 => q(PAK), continent => q(Asia) },
	q(Palau)                                        => { iso2 => q(PW), iso3 => q(PLW), continent => q(Oceania) },
	q(Palestinian territories)                      => { iso2 => q(PS), iso3 => q(PSE), continent => q(Asia) },
	q(Panama)                                       => { iso2 => q(PA), iso3 => q(PAN), continent => q(North America) },
	q(Papua New Guinea)                             => { iso2 => q(PG), iso3 => q(PNG), continent => q(Oceania) },
	q(Paraguay)                                     => { iso2 => q(PY), iso3 => q(PRY), continent => q(South America) },
	q(Peru)                                         => { iso2 => q(PE), iso3 => q(PER), continent => q(South America) },
	q(Philippines)                                  => { iso2 => q(PH), iso3 => q(PHL), continent => q(Asia) },
	q(Pitcairn)                                     => { iso2 => q(PN), iso3 => q(PCN), continent => q(Oceania) },
	q(Poland)                                       => { iso2 => q(PL), iso3 => q(POL), continent => q(Europe) },
	q(Portugal)                                     => { iso2 => q(PT), iso3 => q(PRT), continent => q(Europe) },
	q(Puerto Rico)                                  => { iso2 => q(PR), iso3 => q(PRI), continent => q(North America) },
	q(Qatar)                                        => { iso2 => q(QA), iso3 => q(QAT), continent => q(Asia) },
	q(Réunion)                                     => { iso2 => q(RE), iso3 => q(REU), continent => q(Africa) },
	q(Romania)                                      => { iso2 => q(RO), iso3 => q(ROU), continent => q(Europe) },
	q(Russia)                                       => { iso2 => q(RU), iso3 => q(RUS), continent => q(Asia) },
	q(Rwanda)                                       => { iso2 => q(RW), iso3 => q(RWA), continent => q(Africa) },
	q(Saint Barthélemy)                            => { iso2 => q(BL), iso3 => q(BLM), continent => q(North America) },
	q(Saint Helena)                                 => { iso2 => q(SH), iso3 => q(SHN), continent => q(Africa) },
	q(Saint Kitts and Nevis)                        => { iso2 => q(KN), iso3 => q(KNA), continent => q(North America) },
	q(Saint Lucia)                                  => { iso2 => q(LC), iso3 => q(LCA), continent => q(North America) },
	q(Saint Martin (French Part))                   => { iso2 => q(MF), iso3 => q(MAF), continent => q(North America) },
	q(Saint Pierre and Miquelon)                    => { iso2 => q(PM), iso3 => q(SPM), continent => q(North America) },
	q(Saint Vincent and the Grenadines)             => { iso2 => q(VC), iso3 => q(VCT), continent => q(North America) },
	q(Samoa)                                        => { iso2 => q(WS), iso3 => q(WSM), continent => q(Oceania) },
	q(San Marino)                                   => { iso2 => q(SM), iso3 => q(SMR), continent => q(Europe) },
	q(São Tomé and Príncipe)                     => { iso2 => q(ST), iso3 => q(STP), continent => q(Africa) },
	q(Sark)                                         => { iso2 => q(CQ), iso3 => q(),    continent => q(Europe) },
	q(Saudi Arabia)                                 => { iso2 => q(SA), iso3 => q(SAU), continent => q(Asia) },
	q(Senegal)                                      => { iso2 => q(SN), iso3 => q(SEN), continent => q(Africa) },
	q(Serbia)                                       => { iso2 => q(RS), iso3 => q(SRB), continent => q(Europe) },
	q(Seychelles)                                   => { iso2 => q(SC), iso3 => q(SYC), continent => q(Africa) },
	q(Sierra Leone)                                 => { iso2 => q(SL), iso3 => q(SLE), continent => q(Africa) },
	q(Singapore)                                    => { iso2 => q(SG), iso3 => q(SGP), continent => q(Asia) },
	q(Sint Maarten (Dutch part))                    => { iso2 => q(SX), iso3 => q(SXM), continent => q(North America) },
	q(Slovakia)                                     => { iso2 => q(SK), iso3 => q(SVK), continent => q(Europe) },
	q(Slovenia)                                     => { iso2 => q(SI), iso3 => q(SVN), continent => q(Europe) },
	q(Solomon Islands)                              => { iso2 => q(SB), iso3 => q(SLB), continent => q(Oceania) },
	q(Somalia)                                      => { iso2 => q(SO), iso3 => q(SOM), continent => q(Africa) },
	q(South Africa)                                 => { iso2 => q(ZA), iso3 => q(ZAF), continent => q(Africa) },
	q(South Georgia and the South Sandwich Islands) => { iso2 => q(GS), iso3 => q(SGS), continent => q(Antarctica) },
	q(South Korea)                                  => { iso2 => q(KR), iso3 => q(KOR), continent => q(Asia) },
	q(South Sudan)                                  => { iso2 => q(SS), iso3 => q(SSD), continent => q(Africa) },
	q(Spain)                                        => { iso2 => q(ES), iso3 => q(ESP), continent => q(Europe) },
	q(Sri Lanka)                                    => { iso2 => q(LK), iso3 => q(LKA), continent => q(Asia) },
	q(Sudan)                                        => { iso2 => q(SD), iso3 => q(SDN), continent => q(Africa) },
	q(Suriname)                                     => { iso2 => q(SR), iso3 => q(SUR), continent => q(South America) },
	q(Svalbard and Jan Mayen Islands)               => { iso2 => q(SJ), iso3 => q(SJM), continent => q(Europe) },
	q(Swaziland)                                    => { iso2 => q(SZ), iso3 => q(SWZ), continent => q(Africa) },
	q(Sweden)                                       => { iso2 => q(SE), iso3 => q(SWE), continent => q(Europe) },
	q(Switzerland)                                  => { iso2 => q(CH), iso3 => q(CHE), continent => q(Europe) },
	q(Syria)                                        => { iso2 => q(SY), iso3 => q(SYR), continent => q(Asia) },
	q(Taiwan)                                       => { iso2 => q(TW), iso3 => q(TWN), continent => q(Asia) },
	q(Tajikistan)                                   => { iso2 => q(TJ), iso3 => q(TJK), continent => q(Asia) },
	q(Tanzania)                                     => { iso2 => q(TZ), iso3 => q(TZA), continent => q(Africa) },
	q(Thailand)                                     => { iso2 => q(TH), iso3 => q(THA), continent => q(Asia) },
	q(The Gambia)                                   => { iso2 => q(GM), iso3 => q(GMB), continent => q(Africa) },
	q(The Netherlands)                              => { iso2 => q(NL), iso3 => q(NLD), continent => q(Europe) },
	q(Togo)                                         => { iso2 => q(TG), iso3 => q(TGO), continent => q(Africa) },
	q(Tokelau)                                      => { iso2 => q(TK), iso3 => q(TKL), continent => q(Oceania) },
	q(Tonga)                                        => { iso2 => q(TO), iso3 => q(TON), continent => q(Oceania) },
	q(Trinidad and Tobago)                          => { iso2 => q(TT), iso3 => q(TTO), continent => q(North America) },
	q(Tunisia)                                      => { iso2 => q(TN), iso3 => q(TUN), continent => q(Africa) },
	q(Turkey)                                       => { iso2 => q(TR), iso3 => q(TUR), continent => q(Asia) },
	q(Turkmenistan)                                 => { iso2 => q(TM), iso3 => q(TKM), continent => q(Asia) },
	q(Turks and Caicos Islands)                     => { iso2 => q(TC), iso3 => q(TCA), continent => q(North America) },
	q(Tuvalu)                                       => { iso2 => q(TV), iso3 => q(TUV), continent => q(Oceania) },
	q(Uganda)                                       => { iso2 => q(UG), iso3 => q(UGA), continent => q(Africa) },
	q(UK)                                           => { iso2 => q(GB), iso3 => q(GBR), continent => q(Europe) },
	q(UK [England])          => { iso2 => q(GB-ENG), iso3 => q(GBR), continent => q(Europe) },
	q(UK [Northern Ireland]) => { iso2 => q(GB-NIR), iso3 => q(GBR), continent => q(Europe) },
	q(UK [Scotland])         => { iso2 => q(GB-SCT), iso3 => q(GBR), continent => q(Europe) },
	q(UK [Wales])            => { iso2 => q(GB-WLS), iso3 => q(GBR), continent => q(Europe) },
	q(Ukraine)               => { iso2 => q(UA),     iso3 => q(UKR), continent => q(Europe) },
	q(United Arab Emirates)  => { iso2 => q(AE),     iso3 => q(ARE), continent => q(Asia) },
	q(Uruguay)               => { iso2 => q(UY),     iso3 => q(URY), continent => q(South America) },
	q(US Minor Outlying Islands) => { iso2 => q(UM), iso3 => q(UMI), continent => q(Oceania) },
	q(US Virgin Islands)         => { iso2 => q(VI), iso3 => q(VIR), continent => q(North America) },
	q(USA)                       => { iso2 => q(US), iso3 => q(USA), continent => q(North America) },
	q(Uzbekistan)                => { iso2 => q(UZ), iso3 => q(UZB), continent => q(Asia) },
	q(Vanuatu)                   => { iso2 => q(VU), iso3 => q(VUT), continent => q(Oceania) },
	q(Venezuela)                 => { iso2 => q(VE), iso3 => q(VEN), continent => q(South America) },
	q(Vietnam)                   => { iso2 => q(VN), iso3 => q(VNM), continent => q(Asia) },
	q(Wallis and Futuna Islands) => { iso2 => q(WF), iso3 => q(WLF), continent => q(Oceania) },
	q(Western Sahara)            => { iso2 => q(EH), iso3 => q(ESH), continent => q(Africa) },
	q(Yemen)                     => { iso2 => q(YE), iso3 => q(YEM), continent => q(Asia) },
	q(Zambia)                    => { iso2 => q(ZM), iso3 => q(ZMB), continent => q(Africa) },
	q(Zimbabwe)                  => { iso2 => q(ZW), iso3 => q(ZWE), continent => q(Africa) },
};

#Direct all library logging calls to screen
my $log_conf = LOG_TO_SCREEN;
Log::Log4perl->init( \$log_conf );
my $logger = Log::Log4perl::get_logger('BIGSdb.Script');
my %opts;
GetOptions(
	'database=s' => \$opts{'d'},
	'help'       => \$opts{'h'},
	'query=s'    => \$opts{'query'},
) or die("Error in command line arguments\n");
if ( $opts{'h'} ) {
	show_help();
	exit;
}
foreach my $required (qw(d query)) {
	if ( !$opts{$required} ) {
		show_help();
		exit;
	}
}
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
die "Script initialization failed - server too busy?\n"
  if !defined $script->{'db'};
die "This script can only be run against an isolate database.\n"
  if ( $script->{'system'}->{'dbtype'} // '' ) ne 'isolates';
main();
undef $script;

sub main {
	my $countries = COUNTRIES;
	my $dataset = $script->{'datastore'}->run_query( $opts{'query'}, undef, { fetch => 'all_arrayref', slice => {} } );
	foreach my $record (@$dataset) {
		$record->{'iso2'} = $countries->{ $record->{'country'} }->{'iso2'} // 'XX';
	}
	say encode_json($dataset);
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
    ${bold}json_dataset.pl$norm - Return JSON dataset for isolates defined
    in a project

${bold}SYNOPSIS$norm
    ${bold}json_dataset.pl --database [${under}DATABASE$norm${bold}] --query [${under}SQL$norm${bold}]$norm

${bold}OPTIONS$norm

${bold}--database$norm ${under}NAME$norm
    Database configuration name.
    
${bold}--query$norm ${under}SQL$norm
      
${bold}--help$norm
    This help page.
    
HELP
	return;
}
