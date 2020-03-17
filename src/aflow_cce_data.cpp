// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2020              *
// *                                                                         *
// ***************************************************************************
// Written by Rico Friedrich, Corey Oses, and Marco Esters
// rico.friedrich@duke.edu

#ifndef _AFLOW_CCE_CPP_
#define _AFLOW_CCE_CPP_

#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_cce.h"

namespace cce {
  //CCE_get_corrections_line////////////////////////////////////////////////////////
  // function to get corrections
  string CCE_get_corrections_line(const string& cor_identifier) {
  //CCE corrections per bond for DFT formation energies for polar materials from binary data using AFLOW (PBE, LDA, SCAN, and PBE+U_ICSD) 
  //with PAW data sets for VASP 5.4.4 and measured experimental values according to the CCE paper 
  //Friedrich et al., Coordination corrected ab initio formation enthalpies, npj Comput. Mater. 5, 59 (2019). 
  //The extensions after the binary formula indicate the source from which the experimental value was taken: 
  //NJ=NIST_JANAF 1998; Ba=Barin 1995; if no extension is given, the value from Kubaschewski et al. 1993 was used. 
  //For DFT+U corrections will only be included for the PBE(+U) calculations used to generate the AFLOW ICSD data
  //DFT+U calculations differ not only by the U value but also by other specifications of the implementation, see:
  //Kick, Reuter, and Oberhofer; Intricacies of DFT+U, Not Only in a Numeric Atom Centered Orbital Framework, JCTC (2019); DOI: 10.1021/acs.jctc.8b01211
  //The formation enthalpies per bond from the experimental data (Exp.) need to be multiplied by the number of M-O bonds 
  //and added (sign change compared to usual CCE corrections to be substracted from DFT formation energies) 
  //for all cations to get an estimate for the formation enthalpy only based on experimental data; 
  //positive values for the peroxide and superoxide O-O bond formation enthalpy indicate that taking the cation-O2- 
  //correction for the cation-O bonds is not a very good approximation for per- and superoxides, 
  //i.e. the oxidation state dependence with respect to the oxidation state of O is relevant!)
  //The number of M-O bonds per formula unit in the last column are for the compounds used to get the corrections (at the start of the return string).
  //They can be used to recalculate the used experimental values by multiplying with the exp. corrections per bond in the previous column.
  //For per- and superoxides the number given is the number of O-O bonds per formula unit of the compound used to get the correction.
                                       //cation species            ox. state  PBE (for 298.15K)  PBE (for 0K)  LDA (for 298.15K)  LDA (for 0K)  SCAN (for 298.15K)  SCAN (for 0K)  PBE(+U)_ICSD (for 298.15K)  PBE(+U)_ICSD (for 0K)  Exp. (for 298.15K)  M-X bonds
                                       //                                     (eV/bond)          (eV/bond)     (eV/bond)          (eV/bond)     (eV/bond)           (eV/bond)      (eV/bond)                   (eV/bond)              (eV/bond)           per f. u.
  //OXIDES
  if (cor_identifier=="Ag_+1_O")  {return "Ag_+1_Ag2O             +1         -0.00825           -0.00700      -0.05400           -0.05375      -0.06525            -0.06450       -0.06150                    -0.06000               -0.08050            4  ";}
  if (cor_identifier=="Al_+3_O")  {return "Al_+3_Al2O3            +3          0.18692            0.17783      -0.00733           -0.01683      -0.01242            -0.02217        0.18667                     0.17758               -1.44725            12 ";}
  if (cor_identifier=="As_+5_O")  {return "As_+5_As2O5            +5          0.20220            0.19190      -0.06360           -0.07520       0.02120             0.00920        0.19770                     0.18760               -0.95360            10 ";}
  if (cor_identifier=="B_+3_O")   {return "B_+3_B2O3              +3          0.19517            0.18250      -0.06933           -0.08350      -0.04717            -0.06117        0.19317                     0.18067               -2.19983            6  ";}
  if (cor_identifier=="Ba_+2_O")  {return "Ba_+2_BaO              +2          0.11833            0.11650       0.00983            0.00750       0.00417             0.00200        0.11817                     0.11633               -0.94683            6  ";}
  if (cor_identifier=="Be_+2_O")  {return "Be_+2_BeO              +2          0.19525            0.18750       0.00825            0.00000       0.00600            -0.00225        0.19475                     0.18700               -1.57900            4  ";}
  if (cor_identifier=="Bi_+3_O")  {return "Bi_+3_Bi2O3            +3         -0.02580           -0.02760      -0.17520           -0.17780      -0.03560            -0.03790       -0.02830                    -0.03010               -0.59150            10 ";}
  if (cor_identifier=="Ca_+2_O")  {return "Ca_+2_CaO              +2          0.10567            0.10017      -0.03317           -0.03950      -0.02217            -0.02800        0.10567                     0.10017               -1.09667            6  ";}
  if (cor_identifier=="Cd_+2_O")  {return "Cd_+2_CdO              +2          0.10417            0.10133       0.01617            0.01300       0.00533             0.00233        0.07417                     0.07117               -0.44633            6  ";}
  if (cor_identifier=="Co_+2_O")  {return "Co_+2_CoO              +2          0.23967            0.23733       0.16617            0.16450       0.12300             0.11933        0.00300                    -0.00017               -0.41067            6  ";}
  if (cor_identifier=="Cr_+3_O")  {return "Cr_+3_Cr2O3            +3          0.15283            0.14733       0.04542            0.03908      -0.01892            -0.02467       -0.18417                    -0.18958               -0.98000            12 ";}
  if (cor_identifier=="Cr_+6_O")  {return "Cr_+6_CrO3             +6         -0.13225           -0.14425      -0.30525           -0.32100      -0.28125            -0.29675       -0.06775                    -0.07950               -1.52100            4  ";}
  if (cor_identifier=="Cs_+1_O")  {return "Cs_+1_Cs2O             +1          0.10083            0.10083      -0.05667           -0.05833      -0.00500            -0.00600        0.10117                     0.10117               -0.59767            6  ";}
  if (cor_identifier=="Cu_+1_O")  {return "Cu_+1_Cu2O_NJ          +1          0.13100            0.12950       0.03400            0.03175       0.06350             0.06175        0.06800                     0.06625               -0.44225            4  ";}
  if (cor_identifier=="Cu_+2_O")  {return "Cu_+2_CuO_NJ           +2          0.10150            0.09725      -0.02575           -0.03075       0.00675             0.00200        0.06525                     0.06075               -0.40425            4  ";}
  if (cor_identifier=="Fe_+2_O")  {return "Fe_+2_FeO_NJ           +2          0.17633            0.17283       0.13100            0.12867       0.01883             0.01433       -0.00900                    -0.01300               -0.47000            6  ";}
  if (cor_identifier=="Fe_+3_O")  {return "Fe_+3_Fe2O3            +3          0.16325            0.15858       0.01300            0.00550      -0.06550            -0.07175       -0.03417                    -0.04517               -0.71117            12 ";}
  if (cor_identifier=="Ga_+3_O")  {return "Ga_+3_Ga2O3            +3          0.20090            0.19250       0.00680           -0.00220       0.03860             0.02890        0.16220                     0.15370               -1.12880            10 ";}
  if (cor_identifier=="Ge_+4_O")  {return "Ge_+4_GeO2             +4          0.19917            0.18950      -0.03567           -0.04617       0.03967             0.02900        0.19883                     0.18917               -1.00183            6  ";}
  if (cor_identifier=="Hf_+4_O")  {return "Hf_+4_HfO2_Ba          +4          0.16171            0.15657      -0.02957           -0.03529      -0.02686            -0.03257        0.16286                     0.15786               -1.69486            7  ";}
  if (cor_identifier=="Hg_+2_O")  {return "Hg_+2_HgO              +2          0.17000            0.15250      -0.06500           -0.08650       0.05800             0.03800        0.16050                     0.14300               -0.47050            2  ";}
  if (cor_identifier=="In_+3_O")  {return "In_+3_In2O3            +3          0.13525            0.13025      -0.01667           -0.02217      -0.01300            -0.01858        0.09750                     0.09258               -0.79967            12 ";}
  if (cor_identifier=="Ir_+4_O")  {return "Ir_+4_IrO2_Ba          +4          0.02017            0.01567      -0.18133           -0.18683       0.01800             0.01350       -0.04283                    -0.04717               -0.41917            6  ";}
  if (cor_identifier=="K_+1_O")   {return "K_+1_K2O               +1          0.08300            0.08025      -0.02625           -0.03013      -0.00350            -0.00712        0.08262                     0.07987               -0.47050            8  ";}
  if (cor_identifier=="Li_+1_O")  {return "Li_+1_Li2O             +1          0.07663            0.07038      -0.01538           -0.02225      -0.01175            -0.01863        0.08400                     0.07775               -0.77463            8  ";}
  if (cor_identifier=="Mg_+2_O")  {return "Mg_+2_MgO              +2          0.13350            0.12717       0.00250           -0.00417      -0.00233            -0.00900        0.13183                     0.12550               -1.03917            6  ";}
  if (cor_identifier=="Mn_+2_O")  {return "Mn_+2_MnO              +2          0.25133            0.25133       0.27000            0.26933      -0.03250            -0.03400       -0.01650                    -0.02033               -0.66483            6  ";}
  if (cor_identifier=="Mn_+4_O")  {return "Mn_+4_MnO2             +4          0.06000            0.05233      -0.09433           -0.10300      -0.15817            -0.16667        0.11017                     0.10233               -0.89983            6  ";}
  if (cor_identifier=="Mo_+4_O")  {return "Mo_+4_MoO2             +4          0.02917            0.02150      -0.18433           -0.19267      -0.10533            -0.11367        0.06400                     0.05650               -1.01550            6  ";}
  if (cor_identifier=="Mo_+6_O")  {return "Mo_+6_MoO3             +6         -0.04700           -0.06025      -0.34100           -0.35750      -0.25575            -0.27175        0.07800                     0.06375               -1.93075            4  ";}
  if (cor_identifier=="Na_+1_O")  {return "Na_+1_Na2O_NJ          +1          0.08225            0.07762      -0.00425           -0.00962      -0.01125            -0.01675        0.08200                     0.07737               -0.54150            8  ";}
  if (cor_identifier=="Nb_+2_O")  {return "Nb_+2_NbO              +2          0.05925            0.05325      -0.12625           -0.13275      -0.08450            -0.09100        0.03975                     0.03375               -1.08750            4  ";}
  if (cor_identifier=="Ni_+2_O")  {return "Ni_+2_NiO              +2          0.25583            0.25367       0.15517            0.15117       0.19800             0.19550       -0.07683                    -0.08017               -0.41400            6  ";}
  if (cor_identifier=="Os_+4_O")  {return "Os_+4_OsO2             +4          0.06167            0.05700      -0.14433           -0.14983       0.00117            -0.00400       -0.01267                    -0.01717               -0.50883            6  ";}
  if (cor_identifier=="Os_+8_O")  {return "Os_+8_OsO4             +8         -0.22250           -0.22950      -0.37925           -0.39200      -0.27725            -0.28800       -0.22450                    -0.23175               -1.02000            4  ";}
  if (cor_identifier=="Pb_+2_O")  {return "Pb_+2_PbO              +2          0.00275            0.00325      -0.10875           -0.10925      -0.05475            -0.05500       -0.00075                    -0.00025               -0.56850            4  ";}
  if (cor_identifier=="Pb_+4_O")  {return "Pb_+4_PbO2             +4          0.05683            0.05450      -0.12500           -0.12850      -0.02500            -0.02833        0.04683                     0.04467               -0.47417            6  ";}
  if (cor_identifier=="Pd_+2_O")  {return "Pd_+2_PdO              +2          0.05775            0.05475      -0.07925           -0.08300      -0.02450            -0.02800       -0.02650                    -0.02950               -0.29925            4  ";}
  if (cor_identifier=="Rb_+1_O")  {return "Rb_+1_Rb2O             +1          0.09500            0.09400      -0.02150           -0.02350       0.00488             0.00313        0.09450                     0.09350               -0.43900            8  ";}
  if (cor_identifier=="Re_+4_O")  {return "Re_+4_ReO2_Ba          +4          0.08967            0.08450      -0.12417           -0.13017       0.02117             0.01550        0.05067                     0.04650               -0.77550            6  ";}
  if (cor_identifier=="Re_+6_O")  {return "Re_+6_ReO3_Ba          +6         -0.07267           -0.08033      -0.30400           -0.31250      -0.15967            -0.16817       -0.08233                    -0.08967               -1.01767            6  ";}
  if (cor_identifier=="Rh_+3_O")  {return "Rh_+3_Rh2O3            +3          0.01408            0.00650      -0.13633           -0.14150      -0.05833            -0.06308       -0.11975                    -0.12783               -0.30717            12 ";}
  if (cor_identifier=="Ru_+4_O")  {return "Ru_+4_RuO2             +4         -0.00483           -0.01150      -0.20567           -0.21333      -0.11183            -0.11917       -0.05283                    -0.05800               -0.52683            6  ";}
  if (cor_identifier=="Sb_+3_O")  {return "Sb_+3_Sb2O3            +3          0.12067            0.11533      -0.11350           -0.12117      -0.01983            -0.02667        0.11733                     0.11217               -1.23700            6  ";}
  if (cor_identifier=="Sb_+5_O")  {return "Sb_+5_Sb2O5_Ba         +5          0.10517            0.09700      -0.13233           -0.14175      -0.05725            -0.06692        0.09825                     0.09008               -0.83942            12 ";}
  if (cor_identifier=="Sc_+3_O")  {return "Sc_+3_Sc2O3            +3          0.16175            0.15408      -0.00833           -0.01658      -0.02567            -0.03383        0.16608                     0.15867               -1.64817            12 ";}
  if (cor_identifier=="Se_+4_O")  {return "Se_+4_SeO2             +4          0.07500            0.06567      -0.22767           -0.23967      -0.09633            -0.10833        0.07300                     0.06367               -0.77767            3  ";}
  if (cor_identifier=="Si_+4_O")  {return "Si_+4_SiO2(al-quartz)  +4          0.25300            0.23800      -0.02325           -0.03900      -0.02075            -0.03675        0.25275                     0.23775               -2.36025            4  ";}
  if (cor_identifier=="Sn_+2_O")  {return "Sn_+2_SnO_Ba           +2          0.06700            0.06500      -0.06375           -0.06650      -0.01300            -0.01575        0.12025                     0.11825               -0.74050            4  ";}
  if (cor_identifier=="Sn_+4_O")  {return "Sn_+4_SnO2_Ba          +4          0.15050            0.14333      -0.04600           -0.05400      -0.01533            -0.02367        0.28650                     0.28000               -1.00333            6  ";}
  if (cor_identifier=="Sr_+2_O")  {return "Sr_+2_SrO              +2          0.10817            0.10483      -0.01917           -0.02317      -0.01550            -0.01933        0.10767                     0.10433               -1.02267            6  ";}
  if (cor_identifier=="Te_+4_O")  {return "Te_+4_TeO2_Ba          +4          0.06300            0.05575      -0.20325           -0.21225      -0.08850            -0.09700        0.06275                     0.05575               -0.83800            4  ";}
  if (cor_identifier=="Ti_+2_O")  {return "Ti_+2_TiO              +2          0.11313            0.10667      -0.06667           -0.07375      -0.02646            -0.03312        0.10479                     0.09958               -1.17188            4.8";}
  if (cor_identifier=="Ti_+3_O")  {return "Ti_+3_Ti2O3            +3          0.09800            0.09050      -0.08550           -0.09358      -0.06875            -0.07667        0.07933                     0.07233               -1.31358            12 ";}
  if (cor_identifier=="Ti_+4_O")  {return "Ti_+4_TiO2(rutile)     +4          0.10717            0.09717      -0.09650           -0.10717      -0.12367            -0.13450        0.05217                     0.04233               -1.63067            6  ";}
  if (cor_identifier=="Tl_+1_O")  {return "Tl_+1_Tl2O             +1         -0.00650           -0.00533      -0.06650           -0.06600      -0.06117            -0.06050       -0.00783                    -0.00667               -0.28917            6  ";}
  if (cor_identifier=="Tl_+3_O")  {return "Tl_+3_Tl2O3            +3          0.05358            0.05175      -0.09358           -0.09617      -0.01400            -0.01658        0.04492                     0.04300               -0.33717            12 ";}
  if (cor_identifier=="V_+2_O")   {return "V_+2_VO                +2          0.26367            0.26200       0.12033            0.11517       0.15683             0.15467        0.24633                     0.24167               -0.74583            6  ";}
  if (cor_identifier=="V_+3_O")   {return "V_+3_V2O3              +3          0.09825            0.09183      -0.06608           -0.07342      -0.05350            -0.06000       -0.04833                    -0.05458               -1.05267            12 ";}
  if (cor_identifier=="V_+4_O")   {return "V_+4_VO2               +4          0.04667            0.03750      -0.14967           -0.15983      -0.15367            -0.16367       -0.02633                    -0.03600               -1.23300            6  ";}
  if (cor_identifier=="V_+5_O")   {return "V_+5_V2O5              +5         -0.00820           -0.01890      -0.21180           -0.22480      -0.21790            -0.23070        0.00730                    -0.00370               -1.60670            10 ";}
  if (cor_identifier=="W_+4_O")   {return "W_+4_WO2               +4          0.05667            0.05117      -0.15867           -0.16483      -0.04833            -0.05433        0.05417                     0.04883               -1.01833            6  ";}
  if (cor_identifier=="W_+6_O")   {return "W_+6_WO3               +6          0.00500           -0.00250      -0.20783           -0.21650      -0.13617            -0.14483        0.00867                     0.00117               -1.45567            6  ";}
  if (cor_identifier=="Y_+3_O")   {return "Y_+3_Y2O3              +3          0.13583            0.13042      -0.02192           -0.02800      -0.04825            -0.05425        0.15800                     0.15250               -1.64533            12 ";}
  if (cor_identifier=="Zn_+2_O")  {return "Zn_+2_ZnO              +2          0.18525            0.18025       0.04225            0.03675       0.04550             0.03975        0.00300                    -0.00250               -0.90825            4  ";}
  if (cor_identifier=="Zr_+4_O")  {return "Zr_+4_ZrO2_NJ          +4          0.13929            0.13200      -0.04514           -0.05300      -0.06314            -0.07086        0.18086                     0.17371               -1.62486            7  ";}
  // corrections for per- and superoxides                                                                                                                                                                                                                 
  if (cor_identifier=="O2_-2_O")  {return "O2_-2_Li2O2            -1         -0.09300           -0.08556      -0.11600           -0.11100       0.24000             0.24756       -0.12700                    -0.12000                2.72556            1  ";}
  if (cor_identifier=="O2_-1_O")  {return "O2_-1_KO2              -0.5       -0.54500           -0.54350      -0.26900           -0.26970      -0.04900            -0.04680       -0.54320                    -0.54170                1.75600            1  ";}
  //NITRIDES
  //For GaN Lany, PRB 78, 245207 2008 mentions "that the long-time tabulated value Hf =−1.10 eV (Refs. 10–12) for GaN has recently been corrected to −1.62 eV (Ref. 21)". 
  //However when I compare the values with my calculated LDA, PBE and SCAN data, the -1.10 eV/formula unit fits better to the overall fact that LDA and SCAN overestimate the formation enthalpy (too negative value) and PBE does the opposite. 
  //With the newer value of -1.62 eV/formula unit this would not be the case. Therefore I took the old tabulated value for deducing the corrections (January 07 2020).
  if (cor_identifier=="Al_+3_N")  {return "Al_+3_AlN              +3          0.11875            0.10850      -0.04575           -0.05575      -0.03850            -0.04875        0.11800                     0.10825               -0.82500            4  ";}
  if (cor_identifier=="B_+3_N")   {return "B_+3_BN_NJ             +3          0.00100           -0.00867      -0.15933           -0.17000      -0.11000            -0.12067        0.00167                    -0.00767               -0.86700            3  ";}
  if (cor_identifier=="Ca_+2_N")  {return "Ca_+2_Ca3N2_Ba         +2          0.03967            0.03400      -0.09225           -0.09792      -0.04392            -0.04942        0.03958                     0.03442               -0.37225            12 ";}
  if (cor_identifier=="Li_+1_N")  {return "Li_+1_Li3N             +1          0.01663            0.00912      -0.06850           -0.07562      -0.03588            -0.04338        0.02725                     0.02013               -0.21350            8  ";}
  if (cor_identifier=="Zn_+2_N")  {return "Zn_+2_Zn3N2            +2          0.06733            0.06225      -0.03492           -0.03908       0.00417            -0.00025       -0.06792                    -0.07242               -0.01950            12 ";}
  if (cor_identifier=="Be_+2_N")  {return "Be_+2_Be3N2            +2          0.05908            0.05292      -0.06500           -0.07133      -0.03433            -0.04075        0.05875                     0.05258               -0.50917            12 ";}
  if (cor_identifier=="Cr_+3_N")  {return "Cr_+3_CrN              +3          0.08417            0.08167      -0.02100           -0.02583      -0.01867            -0.02233       -0.22883                    -0.23283               -0.20250            6  ";}
  if (cor_identifier=="Ga_+3_N")  {return "Ga_+3_GaN_!!!          +3          0.04725            0.03950      -0.11725           -0.12550      -0.04200            -0.05075        0.01550                     0.00750               -0.28400            4  ";}
  if (cor_identifier=="Hf_+3_N")  {return "Hf_+3_HfN              +3          0.05983            0.05633      -0.09283           -0.09683      -0.01433            -0.01817        0.05933                     0.05583               -0.64533            6  ";}
  if (cor_identifier=="In_+3_N")  {return "In_+3_InN_Ba           +3          0.08825            0.08350      -0.05700           -0.06225      -0.01125            -0.01650        0.04800                     0.04325               -0.04450            4  ";}
  if (cor_identifier=="La_+3_N")  {return "La_+3_LaN_Ba           +3          0.08300            0.07983      -0.00100           -0.00433      -0.01600            -0.01950        0.13183                     0.12883               -0.52400            6  ";}
  if (cor_identifier=="Mg_+2_N")  {return "Mg_+2_Mg3N2            +2          0.08300            0.07642      -0.03475           -0.04167      -0.01242            -0.01950        0.08083                     0.07425               -0.39858            12 ";}
  if (cor_identifier=="Nb_+3_N")  {return "Nb_+3_NbN              +3          0.07467            0.06967      -0.06550           -0.07083       0.02667             0.02150        0.07400                     0.06900               -0.40833            6  ";}
  if (cor_identifier=="Sc_+3_N")  {return "Sc_+3_ScN              +3         -0.09850           -0.10500      -0.23150           -0.23833      -0.20667            -0.21350       -0.07883                    -0.08517               -0.54200            6  ";}
  if (cor_identifier=="Si_+4_N")  {return "Si_+4_Si3N4            +4          0.00000           -0.01300      -0.22408           -0.23733      -0.15183            -0.16508       -0.00058                    -0.01358               -0.64325            12 ";}
  if (cor_identifier=="Ta_+3_N")  {return "Ta_+3_TaN              +3          0.15283            0.14933       0.00617            0.00233       0.13100             0.12733        0.15650                     0.15317               -0.43583            6  ";}
  if (cor_identifier=="Ti_+3_N")  {return "Ti_+3_TiN              +3          0.00367           -0.00283      -0.14550           -0.15233      -0.08733            -0.09400        0.03083                     0.02467               -0.58400            6  ";}
  if (cor_identifier=="V_+3_N")   {return "V_+3_VN                +3          0.04500            0.03917      -0.09483           -0.10100      -0.01683            -0.02250        0.03283                     0.03100               -0.37650            6  ";}
  if (cor_identifier=="Y_+3_N")   {return "Y_+3_YN                +3         -0.06283           -0.06750      -0.18267           -0.18783      -0.16950            -0.17450       -0.04933                    -0.05400               -0.51683            6  ";}
  if (cor_identifier=="Zr_+3_N")  {return "Zr_+3_ZrN_NJ           +3          0.04250            0.03733      -0.09933           -0.10500      -0.04067            -0.04617        0.06783                     0.06267               -0.63100            6  ";}
  else {return "";}                                                                                                      
  }

  //CCE_get_Bader_templates////////////////////////////////////////////////////////
  // function to get Bader charges from the binaries used to determine the corrections
  string CCE_get_Bader_templates(const string& element) {
  //Bader charges of the elements obtained from binary oxides
  //for DFT+U corrections will only be included for the PBE(+U) calculations used to generate the AFLOW ICSD data
  //DFT+U calculations differ not only by the U value but also by other specifications of the implementation, see:
  //Kick, Reuter, and Oberhofer; Intricacies of DFT+U, Not Only in a Numeric Atom Centered Orbital Framework, JCTC (2019); DOI: 10.1021/acs.jctc.8b01211
  //cation species           #ox.  ox       Bader charge (e charges)
  //                         stat  nr     PBE      LDA     SCAN    PBE(+U)_ICSD 
  if (element=="Ag")  {return "1   +1   0.4794   0.4659   0.4992   0.4513";}
  if (element=="Al")  {return "1   +3   2.4800   2.4533   2.5188   2.4706";}
  if (element=="As")  {return "1   +5   2.5884   2.6167   2.7515   2.5858";}
  if (element=="B")   {return "1   +3   2.3417   2.3087   2.4076   2.3261";}
  if (element=="Ba")  {return "1   +2   1.4375   1.4110   1.4636   1.4534";}
  if (element=="Be")  {return "1   +2   1.7043   1.6823   1.7276   1.7159";}
  if (element=="Bi")  {return "1   +3   1.7693   1.7505   1.7962   1.7707";}
  if (element=="Ca")  {return "1   +2   1.4772   1.4152   1.5117   1.4828";}
  if (element=="Cd")  {return "1   +2   1.1479   1.0778   1.2298   1.1791";}
  if (element=="Co")  {return "1   +2   1.1974   1.1351   1.2440   1.2665";}
  if (element=="Cr")  {return "2   +3   1.7009   1.6240   1.7560   1.7671   +6   2.0198   1.9856   2.1249   2.0498";}
  if (element=="Cs")  {return "1   +1   0.6818   0.6747   0.7132   0.6845";}
  if (element=="Cu")  {return "2   +1   0.5287   0.5236   0.5464   0.5105   +2   0.9651   0.9614   1.0033   1.0437";}
  if (element=="Fe")  {return "2   +2   1.2937   1.2688   1.3440   1.3784   +3   1.6233   1.3824   1.7882   1.8411";}
  if (element=="Ga")  {return "1   +3   1.8610   1.8334   1.9507   1.8600";}
  if (element=="Ge")  {return "1   +4   2.3971   2.3777   2.5471   2.4253";}
  if (element=="Hf")  {return "1   +4   2.5659   2.4994   2.6054   2.5587";}
  if (element=="Hg")  {return "1   +2   0.8813   0.8887   0.9201   0.8884";}
  if (element=="In")  {return "1   +3   1.8191   1.7843   1.9155   1.8147";}
  if (element=="Ir")  {return "1   +4   1.6523   1.6448   1.6469   1.6083";}
  if (element=="K")   {return "1   +1   0.7458   0.7101   0.7660   0.7265";}
  if (element=="Li")  {return "1   +1   0.8162   0.7974   0.8223   0.8372";}
  if (element=="Mg")  {return "1   +2   1.6973   1.6729   1.7147   1.6821";}
  if (element=="Mn")  {return "2   +2   1.3286   1.2956   1.4272   1.3775   +4   1.8460   1.7701   1.9347   1.8942";}
  if (element=="Mo")  {return "2   +4   2.0898   2.0799   2.1917   2.0983   +6   2.6270   2.6499   2.7532   2.5918";}
  if (element=="Na")  {return "1   +1   0.7868   0.7676   0.8073   0.7888";}
  if (element=="Nb")  {return "1   +2   1.3151   1.3781   1.3329   1.2806";}
  if (element=="Ni")  {return "1   +2   1.1075   1.0731   1.1578   1.2148";}
  if (element=="Os")  {return "2   +4   1.8518   1.8489   1.9032   1.8615   +8   2.6008   2.5884   2.6683   2.5591";}
  if (element=="Pb")  {return "2   +2   1.1513   1.1301   1.1698   1.1493   +4   2.0173   2.0324   2.1245   2.0469";}
  if (element=="Pd")  {return "1   +2   0.8312   0.8326   0.8630   0.8242";}
  if (element=="Rb")  {return "1   +1   0.7401   0.7124   0.7634   0.7192";}
  if (element=="Re")  {return "2   +4   2.0264   2.0190   2.0771   2.0627   +6   2.9344   2.9106   3.0224   2.8782";}
  if (element=="Rh")  {return "1   +3   1.2622   1.2465   1.3048   1.2867";}
  if (element=="Ru")  {return "1   +4   1.7022   1.6972   1.7991   1.8402";}
  if (element=="Sb")  {return "2   +3   1.7754   1.7743   1.8537   1.7742   +5   2.7273   2.7153   2.8873   2.7311";}
  if (element=="Sc")  {return "1   +3   2.0251   1.9722   2.0868   2.1215";}
  if (element=="Se")  {return "1   +4   1.9132   1.9535   2.0425   1.8766";}
  if (element=="Si")  {return "1   +4   3.2043   3.1655   3.2758   3.1716";}
  if (element=="Sn")  {return "2   +2   1.1884   1.1704   1.2280   1.1716   +4   2.2780   2.2608   2.4336   2.2097";}
  if (element=="Sr")  {return "1   +2   1.4895   1.4274   1.5335   1.4767";}
  if (element=="Te")  {return "1   +4   2.2551   2.2661   2.3516   2.2452";}
  if (element=="Ti")  {return "3   +2   1.2335   1.2183   1.2637   1.2658   +3   1.9228   1.8736   1.9955   2.0175   +4   2.2327   2.1844   2.3356   2.3956";}
  if (element=="Tl")  {return "2   +1   0.5628   0.5488   0.5779   0.5640   +3   1.4739   1.4718   1.5568   1.4781";}
  if (element=="V")   {return "4   +2   1.4271   1.4437   1.4306   1.5312   +3   1.7739   1.7286   1.8327   1.8783   +4   2.0728   2.0077   2.1737   2.1641   +5   2.1742   2.1430   2.2781   2.1982";}
  if (element=="W")   {return "2   +4   2.2270   2.1974   2.2933   2.2501   +6   2.9727   2.9284   3.0812   2.9838";}
  if (element=="Y")   {return "1   +3   2.1451   2.1027   2.2107   2.1473";}
  if (element=="Zn")  {return "1   +2   1.2134   1.2091   1.2764   1.1739";}
  if (element=="Zr")  {return "1   +4   2.5731   2.5281   2.6697   2.5711";}
  // corrections for per- (Li2O2) and superoxides (KO2)
  if (element=="O")   {return "2   -1   -0.8505  -0.8275  -0.8534  -0.8765  -0.5 -0.4334  -0.4157  -0.4338  -0.4313";}
  else {return "";}                                                                                                      
  }
} // namespace cce

#endif // _AFLOW_CCE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2020              *
// *                                                                         *
// ***************************************************************************
