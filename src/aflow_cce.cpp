// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2019              *
// *                                                                         *
// ***************************************************************************
// Written by Rico Friedrich and Corey Oses
// rico.friedrich@duke.edu

#ifndef _AFLOW_CCE_CPP_
#define _AFLOW_CCE_CPP_

#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_cce.h"

namespace pflow {
  // function to get Bader charges from the binaries used to determine the corrections
  string get_Bader_templates(const string& element) {
  //Bader charges fo the elements for obtained from binary oxides
  //cation species           #ox.  ox    Bader charge (e charges)
  //         stat  nr     PBE      LDA     SCAN 
  if (element=="Ag")  {return "1   +1   0.4794   0.4659   0.4992";}
  if (element=="Al")  {return "1   +3   2.4800   2.4533   2.5188";}
  if (element=="As")  {return "1   +5   2.5884   2.6167   2.7515";}
  if (element=="B")   {return "1   +3   2.3417   2.3087   2.4076";}
  if (element=="Ba")  {return "1   +2   1.4375   1.4110   1.4636";}
  if (element=="Be")  {return "1   +2   1.7043   1.6823   1.7276";}
  if (element=="Bi")  {return "1   +3   1.7693   1.7505   1.7962";}
  if (element=="Ca")  {return "1   +2   1.4772   1.4152   1.5117";}
  if (element=="Cd")  {return "1   +2   1.1479   1.0778   1.2298";}
  if (element=="Co")  {return "1   +2   1.1974   1.1351   1.2440";}
  if (element=="Cr")  {return "2   +3   1.7009   1.6240   1.7560   +6   2.0198   1.9856   2.1249";}
  if (element=="Cs")  {return "1   +1   0.6818   0.6747   0.7132";}
  if (element=="Cu")  {return "2   +1   0.5287   0.5236   0.5464   +2   0.9651   0.9614   1.0033";}
  if (element=="Fe")  {return "2   +2   1.2937   1.2688   1.3440   +3   1.6233   1.3824   1.7882";}
  if (element=="Ga")  {return "1   +3   1.8610   1.8334   1.9507";}
  if (element=="Ge")  {return "1   +4   2.3971   2.3777   2.5471";}
  if (element=="Hf")  {return "1   +4   2.5659   2.4994   2.6054";}
  if (element=="Hg")  {return "1   +2   0.8813   0.8887   0.9201";}
  if (element=="In")  {return "1   +3   1.8191   1.7843   1.9155";}
  if (element=="Ir")  {return "1   +4   1.6523   1.6448   1.6469";}
  if (element=="K")   {return "1   +1   0.7458   0.7101   0.7660";}
  if (element=="Li")  {return "1   +1   0.8162   0.7974   0.8223";}
  if (element=="Mg")  {return "1   +2   1.6973   1.6729   1.7147";}
  if (element=="Mn")  {return "2   +2   1.3286   1.2956   1.4272   +4   1.8460   1.7701   1.9347";}
  if (element=="Mo")  {return "2   +4   2.0898   2.0799   2.1917   +6   2.6270   2.6499   2.7532";}
  if (element=="Na")  {return "1   +1   0.7868   0.7676   0.8073";}
  if (element=="Nb")  {return "1   +2   1.3151   1.3781   1.3329";}
  if (element=="Ni")  {return "1   +2   1.1075   1.0731   1.1578";}
  if (element=="Os")  {return "2   +4   1.8518   1.8489   1.9032   +8   2.6008   2.5884   2.6683";}
  if (element=="Pb")  {return "2   +2   1.1513   1.1301   1.1698   +4   2.0173   2.0324   2.1245";}
  if (element=="Pd")  {return "1   +2   0.8312   0.8326   0.8630";}
  if (element=="Rb")  {return "1   +1   0.7401   0.7124   0.7634";}
  if (element=="Re")  {return "2   +4   2.0264   2.0190   2.0771   +6   2.9344   2.9106   3.0224";}
  if (element=="Rh")  {return "1   +3   1.2622   1.2465   1.3048";}
  if (element=="Ru")  {return "1   +4   1.7022   1.6972   1.7991";}
  if (element=="Sb")  {return "2   +3   1.7754   1.7743   1.8537   +5   2.7273   2.7153   2.8873";}
  if (element=="Sc")  {return "1   +3   2.0251   1.9722   2.0868";}
  if (element=="Se")  {return "1   +4   1.9132   1.9535   2.0425";}
  if (element=="Si")  {return "1   +4   3.2043   3.1655   3.2758";}
  if (element=="Sn")  {return "2   +2   1.1884   1.1704   1.2280   +4   2.2780   2.2608   2.4336";}
  if (element=="Sr")  {return "1   +2   1.4895   1.4274   1.5335";}
  if (element=="Te")  {return "1   +4   2.2551   2.2661   2.3516";}
  if (element=="Ti")  {return "3   +2   1.2335   1.2183   1.2637   +3   1.9228   1.8736   1.9955   +4   2.2327   2.1844   2.3356";}
  if (element=="Tl")  {return "2   +1   0.5628   0.5488   0.5779   +3   1.4739   1.4718   1.5568";}
  if (element=="V")   {return "4   +2   1.4271   1.4437   1.4306   +3   1.7739   1.7286   1.8327   +4   2.0728   2.0077   2.1737   +5   2.1742   2.1430   2.2781";}
  if (element=="W")   {return "2   +4   2.2270   2.1974   2.2933   +6   2.9727   2.9284   3.0812";}
  if (element=="Y")   {return "1   +3   2.1451   2.1027   2.2107";}
  if (element=="Zn")  {return "1   +2   1.2134   1.2091   1.2764";}
  if (element=="Zr")  {return "1   +4   2.5731   2.5281   2.6697";}
  // corrections for per- (Li2O2) and superoxides (KO2)
  if (element=="O")   {return "2   -1   -0.8505  -0.8275  -0.8534  -0.5 -0.4334  -0.4157  -0.4338";}
    else {return "";}                                                                                                      
  }

  // function to get corrections
  string get_corrections(const string& cor_identifier) {
  //CCE corrections per bond for oxide DFT formation energies from binary oxide data using AFLOW (PBE, LDA and SCAN) with PAW data sets for VASP 5.4.4 and measured experimental values according to the CCE paper of Rico Friedrich et al. (2019). The extensions after the oxide formula indicate the source from which the experimental value was taken: NJ=NIST_JANAF 1998; Ba=Barin 1995; if no extension is given, the value from Kubaschewski et al. 1993 was used. The Bader charges are kept for convenience.
               //cation species            ox. state  Bader PBE    PBE (for 298.15K)  PBE (for 0K)  Bader LDA    LDA (for 298.15K)  LDA (for 0K)  Bader SCAN   SCAN (for 298.15K)  SCAN (for 0K)
                                         //                                     (e charges)  (eV/bond)          (eV/bond)     (e charges)  (eV/bond)          (eV/bond)     (e charges)  (eV/bond)           (eV/bond)     
  if (cor_identifier=="Ag_+1_ox")  {return "Ag_+1_Ag2O             +1          0.4794      -0.00825           -0.00700       0.4659      -0.05400           -0.05375       0.4992      -0.06525            -0.06450";}
  if (cor_identifier=="Al_+3_ox")  {return "Al_+3_Al2O3            +3          2.4800       0.18692            0.17783       2.4533      -0.00733           -0.01683       2.5188      -0.01242            -0.02217";}
  if (cor_identifier=="As_+5_ox")  {return "As_+5_As2O5            +5          2.5884       0.20220            0.19190       2.6167      -0.06360           -0.07520       2.7515       0.02120             0.00920";}
  if (cor_identifier=="B_+3_ox")   {return "B_+3_B2O3              +3          2.3417       0.19517            0.18250       2.3087      -0.06933           -0.08350       2.4076      -0.04717            -0.06117";}
  if (cor_identifier=="Ba_+2_ox")  {return "Ba_+2_BaO              +2          1.4375       0.11833            0.11650       1.4110       0.00983            0.00750       1.4636       0.00417             0.00200";}
  if (cor_identifier=="Be_+2_ox")  {return "Be_+2_BeO              +2          1.7043       0.19525            0.18750       1.6823       0.00825            0.00000       1.7276       0.00600            -0.00225";}
  if (cor_identifier=="Bi_+3_ox")  {return "Bi_+3_Bi2O3            +3          1.7693      -0.02580           -0.02760       1.7505      -0.17520           -0.17780       1.7962      -0.03560            -0.03790";}
  if (cor_identifier=="Ca_+2_ox")  {return "Ca_+2_CaO              +2          1.4772       0.10567            0.10017       1.4152      -0.03317           -0.03950       1.5117      -0.02217            -0.02800";}
  if (cor_identifier=="Cd_+2_ox")  {return "Cd_+2_CdO              +2          1.1479       0.10417            0.10133       1.0778       0.01617            0.01300       1.2298       0.00533             0.00233";}
  if (cor_identifier=="Co_+2_ox")  {return "Co_+2_CoO              +2          1.1974       0.23967            0.23733       1.1351       0.16617            0.16450       1.2440       0.12300             0.11933";}
  if (cor_identifier=="Cr_+3_ox")  {return "Cr_+3_Cr2O3            +3          1.7009       0.15283            0.14733       1.6240       0.04542            0.03908       1.7560      -0.01892            -0.02467";}
  if (cor_identifier=="Cr_+6_ox")  {return "Cr_+6_CrO3             +6          2.0198      -0.13225           -0.14425       1.9856      -0.30525           -0.32100       2.1249      -0.28125            -0.29675";}
  if (cor_identifier=="Cs_+1_ox")  {return "Cs_+1_Cs2O             +1          0.6818       0.10083            0.10083       0.6747      -0.05667           -0.05833       0.7132      -0.00500            -0.00600";}
  if (cor_identifier=="Cu_+1_ox")  {return "Cu_+1_Cu2O_NJ          +1          0.5287       0.13100            0.12950       0.5236       0.03400            0.03175       0.5464       0.06350             0.06175";}
  if (cor_identifier=="Cu_+2_ox")  {return "Cu_+2_CuO_NJ           +2          0.9651       0.10150            0.09725       0.9614      -0.02575           -0.03075       1.0033       0.00675             0.00200";}
  if (cor_identifier=="Fe_+2_ox")  {return "Fe_+2_FeO_NJ           +2          1.2937       0.17633            0.17283       1.2688       0.13100            0.12867       1.3440       0.01883             0.01433";}
  if (cor_identifier=="Fe_+3_ox")  {return "Fe_+3_Fe2O3            +3          1.6233       0.16325            0.15858       1.3824       0.01300            0.00550       1.7882      -0.06550            -0.07175";}
  if (cor_identifier=="Ga_+3_ox")  {return "Ga_+3_Ga2O3            +3          1.8610       0.20090            0.19250       1.8334       0.00680           -0.00220       1.9507       0.03860             0.02890";}
  if (cor_identifier=="Ge_+4_ox")  {return "Ge_+4_GeO2             +4          2.3971       0.19917            0.18950       2.3777      -0.03567           -0.04617       2.5471       0.03967             0.02900";}
  if (cor_identifier=="Hf_+4_ox")  {return "Hf_+4_HfO2_Ba          +4          2.5659       0.16171            0.15657       2.4994      -0.02957           -0.03529       2.6054      -0.02686            -0.03257";}
  if (cor_identifier=="Hg_+2_ox")  {return "Hg_+2_HgO              +2          0.8813       0.17000            0.15250       0.8887      -0.06500           -0.08650       0.9201       0.05800             0.03800";}
  if (cor_identifier=="In_+3_ox")  {return "In_+3_In2O3            +3          1.8191       0.13525            0.13025       1.7843      -0.01667           -0.02217       1.9155      -0.01300            -0.01858";}
  if (cor_identifier=="Ir_+4_ox")  {return "Ir_+4_IrO2_Ba          +4          1.6523       0.02017            0.01567       1.6448      -0.18133           -0.18683       1.6469       0.01800             0.01350";}
  if (cor_identifier=="K_+1_ox")   {return "K_+1_K2O               +1          0.7458       0.08300            0.08025       0.7101      -0.02625           -0.03013       0.7660      -0.00350            -0.00712";}
  if (cor_identifier=="Li_+1_ox")  {return "Li_+1_Li2O             +1          0.8162       0.07663            0.07038       0.7974      -0.01538           -0.02225       0.8223      -0.01175            -0.01863";}
  if (cor_identifier=="Mg_+2_ox")  {return "Mg_+2_MgO              +2          1.6973       0.13350            0.12717       1.6729       0.00250           -0.00417       1.7147      -0.00233            -0.00900";}
  if (cor_identifier=="Mn_+2_ox")  {return "Mn_+2_MnO              +2          1.3286       0.25133            0.25133       1.2956       0.27000            0.26933       1.4272      -0.03250            -0.03400";}
  if (cor_identifier=="Mn_+4_ox")  {return "Mn_+4_MnO2             +4          1.8460       0.06000            0.05233       1.7701      -0.09433           -0.10300       1.9347      -0.15817            -0.16667";}
  if (cor_identifier=="Mo_+4_ox")  {return "Mo_+4_MoO2             +4          2.0898       0.02917            0.02150       2.0799      -0.18433           -0.19267       2.1917      -0.10533            -0.11367";}
  if (cor_identifier=="Mo_+6_ox")  {return "Mo_+6_MoO3             +6          2.6270      -0.04700           -0.06025       2.6499      -0.34100           -0.35750       2.7532      -0.25575            -0.27175";}
  if (cor_identifier=="Na_+1_ox")  {return "Na_+1_Na2O_NJ          +1          0.7868       0.08225            0.07762       0.7676      -0.00425           -0.00962       0.8073      -0.01125            -0.01675";}
  if (cor_identifier=="Nb_+2_ox")  {return "Nb_+2_NbO              +2          1.3151       0.05925            0.05325       1.3781      -0.12625           -0.13275       1.3329      -0.08450            -0.09100";}
  if (cor_identifier=="Ni_+2_ox")  {return "Ni_+2_NiO              +2          1.1075       0.25583            0.25367       1.0731       0.15517            0.15117       1.1578       0.19800             0.19550";}
  if (cor_identifier=="Os_+4_ox")  {return "Os_+4_OsO2             +4          1.8518       0.06167            0.05700       1.8489      -0.14433           -0.14983       1.9032       0.00117            -0.00400";}
  if (cor_identifier=="Os_+8_ox")  {return "Os_+8_OsO4             +8          2.6008      -0.22250           -0.22950       2.5884      -0.37925           -0.39200       2.6683      -0.27725            -0.28800";}
  if (cor_identifier=="Pb_+2_ox")  {return "Pb_+2_PbO              +2          1.1513       0.00275            0.00325       1.1301      -0.10875           -0.10925       1.1698      -0.05475            -0.05500";}
  if (cor_identifier=="Pb_+4_ox")  {return "Pb_+4_PbO2             +4          2.0173       0.05683            0.05450       2.0324      -0.12500           -0.12850       2.1245      -0.02500            -0.02833";}
  if (cor_identifier=="Pd_+2_ox")  {return "Pd_+2_PdO              +2          0.8312       0.05775            0.05475       0.8326      -0.07925           -0.08300       0.8630      -0.02450            -0.02800";}
  if (cor_identifier=="Rb_+1_ox")  {return "Rb_+1_Rb2O             +1          0.7401       0.09500            0.09400       0.7124      -0.02150           -0.02350       0.7634       0.00488             0.00313";}
  if (cor_identifier=="Re_+4_ox")  {return "Re_+4_ReO2_Ba          +4          2.0264       0.08967            0.08450       2.0190      -0.12417           -0.13017       2.0771       0.02117             0.01550";}
  if (cor_identifier=="Re_+6_ox")  {return "Re_+6_ReO3_Ba          +6          2.9344      -0.07267           -0.08033       2.9106      -0.30400           -0.31250       3.0224      -0.15967            -0.16817";}
  if (cor_identifier=="Rh_+3_ox")  {return "Rh_+3_Rh2O3            +3          1.2622       0.01408            0.00650       1.2465      -0.13633           -0.14150       1.3048      -0.05833            -0.06308";}
  if (cor_identifier=="Ru_+4_ox")  {return "Ru_+4_RuO2             +4          1.7022      -0.00483           -0.01150       1.6972      -0.20567           -0.21333       1.7991      -0.11183            -0.11917";}
  if (cor_identifier=="Sb_+3_ox")  {return "Sb_+3_Sb2O3            +3          1.7754       0.12067            0.11533       1.7743      -0.11350           -0.12117       1.8537      -0.01983            -0.02667";}
  if (cor_identifier=="Sb_+5_ox")  {return "Sb_+5_Sb2O5_Ba         +5          2.7273       0.10517            0.09700       2.7153      -0.13233           -0.14175       2.8873      -0.05725            -0.06692";}
  if (cor_identifier=="Sc_+3_ox")  {return "Sc_+3_Sc2O3            +3          2.0251       0.16175            0.15408       1.9722      -0.00833           -0.01658       2.0868      -0.02567            -0.03383";}
  if (cor_identifier=="Se_+4_ox")  {return "Se_+4_SeO2             +4          1.9132       0.07500            0.06567       1.9535      -0.22767           -0.23967       2.0425      -0.09633            -0.10833";}
  if (cor_identifier=="Si_+4_ox")  {return "Si_+4_SiO2(al-quartz)  +4          3.2043       0.25300            0.23800       3.1655      -0.02325           -0.03900       3.2758      -0.02075            -0.03675";}
  if (cor_identifier=="Sn_+2_ox")  {return "Sn_+2_SnO_Ba           +2          1.1884       0.06700            0.06500       1.1704      -0.06375           -0.06650       1.2280      -0.01300            -0.01575";}
  if (cor_identifier=="Sn_+4_ox")  {return "Sn_+4_SnO2_Ba          +4          2.2780       0.15050            0.14333       2.2608      -0.04600           -0.05400       2.4336      -0.01533            -0.02367";}
  if (cor_identifier=="Sr_+2_ox")  {return "Sr_+2_SrO              +2          1.4895       0.10817            0.10483       1.4274      -0.01917           -0.02317       1.5335      -0.01550            -0.01933";}
  if (cor_identifier=="Te_+4_ox")  {return "Te_+4_TeO2_Ba          +4          2.2551       0.06300            0.05575       2.2661      -0.20325           -0.21225       2.3516      -0.08850            -0.09700";}
  if (cor_identifier=="Ti_+2_ox")  {return "Ti_+2_TiO              +2          1.2335       0.11313            0.10667       1.2183      -0.06667           -0.07375       1.2637      -0.02646            -0.03312";}
  if (cor_identifier=="Ti_+3_ox")  {return "Ti_+3_Ti2O3            +3          1.9228       0.09800            0.09050       1.8736      -0.08550           -0.09358       1.9955      -0.06875            -0.07667";}
  if (cor_identifier=="Ti_+4_ox")  {return "Ti_+4_TiO2(rutile)     +4          2.2327       0.10717            0.09717       2.1844      -0.09650           -0.10717       2.3356      -0.12367            -0.13450";}
  if (cor_identifier=="Tl_+1_ox")  {return "Tl_+1_Tl2O             +1          0.5628      -0.00650           -0.00533       0.5488      -0.06650           -0.06600       0.5779      -0.06117            -0.06050";}
  if (cor_identifier=="Tl_+3_ox")  {return "Tl_+3_Tl2O3            +3          1.4739       0.05358            0.05175       1.4718      -0.09358           -0.09617       1.5568      -0.01400            -0.01658";}
  if (cor_identifier=="V_+2_ox")   {return "V_+2_VO                +2          1.4271       0.26367            0.26200       1.4437       0.12033            0.11517       1.4306       0.15683             0.15467";}
  if (cor_identifier=="V_+3_ox")   {return "V_+3_V2O3              +3          1.7739       0.09825            0.09183       1.7286      -0.06608           -0.07342       1.8327      -0.05350            -0.06000";}
  if (cor_identifier=="V_+4_ox")   {return "V_+4_VO2               +4          2.0728       0.04667            0.03750       2.0077      -0.14967           -0.15983       2.1737      -0.15367            -0.16367";}
  if (cor_identifier=="V_+5_ox")   {return "V_+5_V2O5              +5          2.1742      -0.00820           -0.01890       2.1430      -0.21180           -0.22480       2.2781      -0.21790            -0.23070";}
  if (cor_identifier=="W_+4_ox")   {return "W_+4_WO2               +4          2.2270       0.05667            0.05117       2.1974      -0.15867           -0.16483       2.2933      -0.04833            -0.05433";}
  if (cor_identifier=="W_+6_ox")   {return "W_+6_WO3               +6          2.9727       0.00500           -0.00250       2.9284      -0.20783           -0.21650       3.0812      -0.13617            -0.14483";}
  if (cor_identifier=="Y_+3_ox")   {return "Y_+3_Y2O3              +3          2.1451       0.13583            0.13042       2.1027      -0.02192           -0.02800       2.2107      -0.04825            -0.05425";}
  if (cor_identifier=="Zn_+2_ox")  {return "Zn_+2_ZnO              +2          1.2134       0.18525            0.18025       1.2091       0.04225            0.03675       1.2764       0.04550             0.03975";}
  if (cor_identifier=="Zr_+4_ox")  {return "Zr_+4_ZrO2_NJ          +4          2.5731       0.13929            0.13200       2.5281      -0.04514           -0.05300       2.6697      -0.06314            -0.07086";}
  // corrections for per- and superoxides
  if (cor_identifier=="O2_-2_ox")  {return "O2_-2_Li2O2            -1         -0.8505      -0.09300           -0.08556      -0.8275      -0.11600           -0.11100      -0.8534       0.24000             0.24756";}
  if (cor_identifier=="O2_-1_ox")  {return "O2_-1_KO2              -0.5       -0.4334      -0.54500           -0.54350      -0.4157      -0.26900           -0.26970      -0.4338      -0.04900            -0.04680";}
    else {return "";}                                                                                                      
  }

  // main CCE function for analyzing the structure, determining the oxidation numbers, assigning corrections, calculating total corrections and writing output
  void CCE(aurostd::xoption& flags){
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string soliloquy="pflow::CCE_CORRECTION():";
  stringstream message;

  // option to print user instructions and exit upon completion
  if(flags.flag("CCE_CORRECTION::USAGE")){
    cout << endl;
    cout << "Written by Rico Friedrich and Corey Oses, 2019" << endl;
    cout << endl;
    cout << "USER INSTRUCTIONS:" << endl;
    cout << endl;
    cout << "(i) GENERAL INFORMATION:" << endl;
    cout << "Implementation to obtain corrected DFT formation enthalpies based on the coordination corrected" << endl; 
    cout << "enthalpies (CE) methodology described in:" << endl; 
    cout << "Friedrich et al., Coordination corrected ab initio formation enthalpies, npj Comput. Mater. 5, 59 (2019);" << endl;
    cout << "https://doi.org/10.1038/s41524-019-0192-1." << endl;
    cout << "Please cite this article when using this method and/or this implementation." << endl;
    cout << "The corrections depend on the number of cation-anion bonds and on the cation oxidation state." << endl;
    cout << "More general information and a list of elements and oxidation states for which corrections are available" << endl; 
    cout << "can be found in README_AFLOW_CCE.TXT." << endl;
    cout << endl;
    cout << "(ii) AVAILABLE OPTIONS:" << endl;
    cout << "--cce                       Prints these user instructions." << endl;
    cout << "--cce=POSCAR_FILE_PATH      Provide the path to the structure file in VASP5 POSCAR format." << endl; 
    cout << "--oxidation_numbers=        Provide as a comma separated list the oxidation numbers. It is" << endl;
    cout << "                            assumed that: (i) one is provided for each atom of the structure and" << endl; 
    cout << "                            (ii) they are in the same sequence as the corresponding atoms in the" << endl;
    cout << "                            provided POSCAR file." << endl;
    cout << "--dft_formation_energies=   Provide a comma separated list of precalculated DFT formation energies," << endl; 
    cout << "                            they are assumed to be: (i) negative for compounds lower in energy" << endl; 
    cout << "                            than the elements, (ii) in eV/cell. Currently, corrections are available" << endl; 
    cout << "                            for PBE, LDA and SCAN (no DFT+U yet)." << endl;
    cout << "--functionals=              Provide a comma separated list with the functionals in the same sequence" << endl; 
    cout << "                            as the DFT formation energies they correspond to. These can be:" << endl; 
    cout << "                            (i) PBE, (ii) LDA or (iii) SCAN. Default: PBE (if only one DFT formation" << endl; 
    cout << "                            energy is provided)." << endl;
    cout << endl;
    cout << "(iii) EXAMPLE INPUT STRUCTURE FOR ROCKSALT MgO:" << endl;
    cout << "Mg1O1   [FCC,FCC,cF8] (STD_PRIM doi:10.1  [FCC,FCC,cF8] (STD_PRIM doi:10.1016/j.commatsci.2010.05.010)" << endl;
    cout << "1.224745" << endl;
      cout << "   0.00000000000000   1.73568248770103   1.73568248770103" << endl;
      cout << "   1.73568248770103   0.00000000000000   1.73568248770103" << endl;
      cout << "   1.73568248770103   1.73568248770103   0.00000000000000" << endl;
    cout << "Mg O" << endl;
    cout << "1 1" << endl;
    cout << "Direct(2) [A1B1]" << endl;
      cout << "   0.00000000000000   0.00000000000000   0.00000000000000  Mg" << endl;
      cout << "   0.50000000000000   0.50000000000000   0.50000000000000  O" << endl;
    cout << endl;
    cout << "(iv) EXAMPLE COMMANDS:" << endl;
    cout << "Assuming that AFLOW is in your PATH and you saved the above example structure file for MgO" << endl; 
    cout << "in the current directory as POSCAR, the following commands can be executed:" << endl;
    cout << endl;
    cout << "aflow --cce=POSCAR --dft_formation_energies=-5.434,-6.220,-6.249 --functionals=PBE,LDA,SCAN --oxidation_numbers=2,-2" << endl;
    cout << "This will give you the CCE corrections and CCE formation enthalpies for PBE, LDA and SCAN for MgO." << endl;
    cout << endl;
    cout << "aflow --cce=POSCAR --dft_formation_energies=-6.220 --functionals=LDA --oxidation_numbers=2,-2" << endl;
    cout << "This gives you only the CCE corrections and CCE formation enthalpies for LDA." << endl;
    cout << endl;
    cout << "aflow --cce=POSCAR --dft_formation_energies=-5.434 --oxidation_numbers=2,-2" << endl;
    cout << "This gives you the CCE corrections and CCE formation enthalpies for PBE with a warning that" << endl; 
    cout << "PBE is assumed as functional." << endl;
    cout << endl;
    cout << "aflow --cce=POSCAR --oxidation_numbers=2,-2" << endl;
    cout << "This gives you only the CCE corrections for PBE, LDA and SCAN." << endl;
    cout << endl;
    cout << "aflow --cce=POSCAR" << endl;
    cout << "This gives you the CCE corrections for PBE, LDA and SCAN if an aflow.in and a Bader charges file" << endl; 
    cout << "in AFLOW format are present and the oxidation numbers can be determined from them." << endl;
    cout << endl;
    //[CO190620 - do NOT exit]exit(EXIT_FAILURE);
    return;
  }

  cout << endl;

  //read structural data
        xstructure a(flags.getattachedscheme("CCE_CORRECTION::POSCAR_PATH"),IOAFLOW_AUTO);
  a.ReScale(1.0); // rescales scaling factor in second line of POSCAR to 1, needed for correct distances
  //let the program spit out what it thinks (input structure)
  if(LDEBUG){
    cout << "INPUT STRUCTURE:" << endl;
    cerr << soliloquy << " input structure:" << endl;
    cerr << a << endl;
  }
  // if species of atoms are not known as in VASP4 format, throw error
  if (a.atoms[0].name == ""){
    message << "BAD NEWS: It seems you are providing a POSCAR without species information as input. This implementation requires a POSCAR in VASP5 format with the species information included. Please adjust the structure file and rerun.";
      throw aurostd::xerror(soliloquy,message,_INPUT_ILLEGAL_);
  }
  // if there is only one specie it must be an elemental phase and is hence not correctable
  if (a.species.size() == 1){
    message << "BAD NEWS: There is only one species in this system. Hence it is an elemental phase whose enthalpy cannot be corrected based on the CCE methodology.";
      throw aurostd::xerror(soliloquy,message,_INPUT_ILLEGAL_);
  }

  //read precalculated DFT formation energies if provided
  string dft_energies_input_str=flags.getattachedscheme("CCE_CORRECTION::DFT_FORMATION_ENERGIES"); // option to provide precalculated DFT formation energies
  vector<double> dft_energies;
  if(!dft_energies_input_str.empty()){aurostd::string2tokens<double>(dft_energies_input_str,dft_energies,",");} //if the DFT energies input string is not empty, convert it to vector of doubles
  //read corresponding functionals
  string functionals;
  functionals=flags.getattachedscheme("CCE_CORRECTION::FUNCTIONALS"); // option to provide functional names corresponding to precalculated DFT formation energies
  vector<string> vfunctionals;
  if(!functionals.empty()){aurostd::string2tokens(functionals,vfunctionals,",");} // if functionals input string is not empty, convert it to vector of strings
  //use PBE as default if only 1 DFT formation energy is provided
      ostream& oss = cout;
  ofstream FileMESSAGE;
  _aflags aflags;aflags.Directory=".";
  bool pbe=0;
  bool lda=0;
  bool scan=0;
  if(functionals.empty() && dft_energies.size() == 1){
    message << "Setting functionals=PBE since only 1 DFT formation energy is provided and PBE is the default functional!";
              pflow::logger(soliloquy, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
    functionals="PBE";
    pbe=1;
    aurostd::string2tokens(functionals,vfunctionals,","); // if functionals input string is not empty, convert it to vector of strings
  // otherwise if sizes of provided DFT formation energies and functionals do not match, throw error
  } else if(dft_energies.size()!=vfunctionals.size()){ // checking only whether the sizes are equal should suffice since then none of them alone can be empty
    message << "BAD NEWS: The number of provided precalculated DFT formation energies and functionals must match.";
      throw aurostd::xerror(soliloquy,message,_INPUT_ILLEGAL_);
  }
  //let the program spit out what it thinks
  if(LDEBUG){
    bool roff=false;
    cout << "INPUT DFT FORMATION ENERGIES & FUNCTIONALS:" << endl;
    cerr << soliloquy << " input dft formation energies=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(dft_energies,6,roff),",") << " (assumed to be in eV/cell)" << endl;
    cerr << soliloquy << " input functionals=" << functionals << endl;
  }
  // determine functional from input to select right corrections to be used
  bool correctable=1; // first assuming that formation engery of system IS correctable; will be set to zero (not correctable) if, for any atom, no correction can be identified
  for(uint k=0,ksize=vfunctionals.size();k<ksize;k++){
    if(LDEBUG){
      cout << "vfunctionals[" << k << "]: " << vfunctionals[k] << endl;
    }
    if (vfunctionals[k]!="PBE" && vfunctionals[k]!="LDA" && vfunctionals[k]!="SCAN"){
      cout << "It seems you are not using the correct naming for the functionals. Please check and rerun." << endl;
      correctable=0;
    } else if (vfunctionals[k] == "PBE"){
      pbe=1;
    } else if (vfunctionals[k] == "LDA"){
      lda=1;
    } else if (vfunctionals[k] == "SCAN"){
      scan=1;
    }
  }
  // if only structure (and oxidation numbers) are provided, corrections should be given for all functionals
  if(functionals.empty() && dft_energies_input_str.empty()){
    pbe=1; lda=1; scan=1;
  }
  if(LDEBUG){
    cout << "PBE: " << pbe << endl;
    cout << "LDA: " << lda << endl;
    cout << "SCAN: " << scan << endl;
    cout << endl; // to get bit clearer output
  }
  // save functionals determined from input in another variable. If the functionals change based on the Bader analysis to determine the oxidation numbers, corrections and corrected formation enthalpies should still be given for all three functionals, but a WARNING should be thrown that the oxidation numbers were only determined from a calculations with a specific functional.
  bool pbe_0=pbe;
  bool lda_0=lda;
  bool scan_0=scan;

  //implement option to provide oxidation numbers by user
  double oxidation_states[a.atoms.size()]; // to check whether sum over all atoms equals zero; otherwise correction assignment might be wrong; double because for superoxides O ox. number is -0.5
  string oxidation_numbers_input_str=""; // option to provide precalculated DFT formation energies
  oxidation_numbers_input_str=flags.getattachedscheme("CCE_CORRECTION::OXIDATION_NUMBERS"); // option to provide precalculated DFT formation energies
  bool ox_nums_provided=0;
  if(oxidation_numbers_input_str != ""){
    ox_nums_provided=1;
    vector<double> oxidation_states_vec;
    if(!oxidation_numbers_input_str.empty()){aurostd::string2tokens<double>(oxidation_numbers_input_str,oxidation_states_vec,",");} //if the DFT energies input string is not empty, convert it to vector of doubles
    //sizes of oxidation numbers and atoms must match
    if(!oxidation_states_vec.empty() && oxidation_states_vec.size()!=a.atoms.size()){ 
      message << "BAD NEWS: The number of provided oxidation numbers does not match the number of atoms in the structure! Please correct and rerun.";
          throw aurostd::xerror(soliloquy,message,_INPUT_ILLEGAL_);
    }
    for(uint k=0,ksize=a.atoms.size();k<ksize;k++){
      oxidation_states[k]=oxidation_states_vec[k];
    }
    // print oxidation numbers
    cout << "INPUT OXIDATION NUMBERS:" << endl;
    for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
      cout << "Oxidation state of " << a.atoms[k].name << " (atom[" << k << "]): " << oxidation_states[k] << endl;
    }
    // calculate sum of oxidation numbers
    double oxidation_sum=0; // double because for superoxides O ox. number is -0.5
    for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
      oxidation_sum+= oxidation_states[k];
    }
    cout << "Sum over all oxidation numbers is (should be ZERO): " << oxidation_sum << endl;
    cout << endl;
    // system should not be regarded correctable if sum over oxidation states is not zero
    if (oxidation_sum != 0){
      message << "BAD NEWS: The formation energy of this system is not correctable! The oxidation numbers that you provided do not add up to zero! Please correct and rerun.";
          throw aurostd::xerror(soliloquy,message,_INPUT_ILLEGAL_);
    }
  } else {
    ox_nums_provided=0;
  }


  // DETERMINE NUMBER OF NEAREST O NEIGHBORS FOR EACH CATION: STRUCTURAL PART OF COORECTION ##############################################################################################################################################################
  // determine which species is oxygen
  uint number_O;
  if(LDEBUG){
    cout << "STRUCTURAL ANALYSIS:" << endl;
  }
  for(uint k=0,ksize=a.species.size();k<ksize;k++){
    //cout << "a.species[k] " << a.species[k] << endl;
    if ( a.species[k] == "O") {
      number_O=k+1; // for counting starting from 1, needed for distsij matrix later
    }
  }
  if(LDEBUG){
    cout << "species number of O is: " << number_O << endl;
  }

  //double nndist=NearestNeighbour(a);
  //vector<double> vbondxx=GetNBONDXX(a); //seems to give/store the nearest neighbour bond lengths between all species pairs, i. e. Mg-Mg, Mg-O, O-O; dist matrix considered to be more useful (intuitive) and hence get bond vector was not used

  // determine species selective nearest neighbor cation-O distances and then cutoffs accordingly
  double tolerance=0.5; // 0.5 Ang tolerance between shortest and longest bonds for each cation-anion pair; works best up to know; in future maybe bonding could be explicitly determined via Bader analysis
  double cutoffs[a.species.size()];
  double cutoffs_max;
  xmatrix<double> distsij=GetDistMatrix(a); //seems to give/store the matrix with the nearest neighbor distances between all species pairs with species running over rows and columns; example MgO: first row Mg, second O, first column Mg, second O
  for(uint i=1,isize=distsij.rows+1;i<isize;i++){
    for(uint j=i+1,jsize=distsij.cols+1;j<jsize;j++){  // need to consider only matrix elements above main diagonal j>i
      if (i==number_O){ // need to know the other index for species selective cutoff determination and assignment
        if(LDEBUG){
          cout << "nearest neighbor cation-O distance for species " << j << " is " << distsij(i,j) << " Ang." << endl;
          cout << "cutoff for the cation-O distance for species " << j << " is " << distsij(i,j)+tolerance << " Ang." << endl; 
        }
        cutoffs[j-1]=distsij(i,j)+tolerance; // -1 since counting of array elements starts from zero, NOT 1
      } else if (j==number_O){ // need to know the other index for species selective cutoff determination and assignment
        if(LDEBUG){
          cout << "nearest neighbor cation-O distance for species " << i << " is " << distsij(i,j) << " Ang." << endl;
          cout << "cutoff for the cation-O distance for species " << i << " is " << distsij(i,j)+tolerance << " Ang." << endl;
        }
        cutoffs[i-1]=distsij(i,j)+tolerance; // -1 since counting of array elements starts from zero, NOT 1
      }
    }
  }
  cutoffs[number_O-1]=0; // set cutoff for O species to zero
  cutoffs_max=0;
  // determine the maximum cutoff between all species and oxygen
  for (uint i=0,isize=a.species.size();i<isize;i++){
    //cout << "cutoffs[" << i << "]= " << cutoffs[i] << endl;
    if (cutoffs[i]>cutoffs_max){
      cutoffs_max=cutoffs[i];
    }
  }
  if(LDEBUG){
    cout << "cutoffs_max= " << cutoffs_max << " Ang." << endl;
  }

  // apply species selective cutoffs to determine/output only nearest neighbors within respective cutoff
  int num_neighbors[a.atoms.size()];
  double perox_cutoff=1.6; // O-O bonds in peroxides for the studied examples are all shorter than 1.6 Ang
  double superox_cutoff=1.4; // O-O bonds in superoxides for the studied examples are all shorter than 1.4 Ang
  int num_perox_bonds;
  int num_superox_bonds;
  int perox_count=0;
  int superox_count=0;
  int perox_indices[a.atoms.size()]; // array in which elements will be 1 for peroxide O atoms and 0 otherwise; needed for correct setting of oxidation numbers below
  int superox_indices[a.atoms.size()]; // array in which elements will be 1 for superoxide O atoms and 0 otherwise; needed for correct setting of oxidation numbers below
  for ( uint i = 0; i < a.atoms.size(); i++ ) { // initialize elements of arrays to 0 
    perox_indices[i] = 0; 
    superox_indices[i] = 0; 
  }
  deque<deque<_atom> > neigh_mat;
  a.GetStrNeighData(cutoffs_max,neigh_mat);
  for(uint i=0,isize=neigh_mat.size();i<isize;i++){ //same size as a.atoms.size(); number of atoms in the structure (not determined by cutoff (or cutoffs_max))
    if (a.atoms[i].name != "O"){
      //cout << "ATOM[" << i << "].cpos=" << a.atoms[i].cpos << "; ATOM[" << i << "].type " << a.atoms[i].type << "; cutoff[species " << a.atoms[i].type << "] " << cutoffs[a.atoms[i].type] << endl;
      int neighbors_count=0;
      bool warning=0;
      for(uint j=0,jsize=neigh_mat[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by the cutoff (cutoffs_max)
        const _atom& atom=neigh_mat[i][j];
        if (0.1 < AtomDist(a.atoms[i],atom) && AtomDist(a.atoms[i],atom) <= cutoffs[a.atoms[i].type] ){ // distance must be larger than 0.1 Ang to savely exclude the cation itself having distance zero to itself
          //cout << atom.name << " cpos=" << atom.cpos << " ijk=" << aurostd::joinWDelimiter(atom.ijk,",") << " dist=" << AtomDist(a,a.atoms[i],atom) <<  " dist_check=" << aurostd::modulus(a.atoms[i].cpos-atom.cpos) << endl; // the first argument to AtomDist function here seems to be unneccesary and might even lead to wrong results since rescaling was performed before; therefore it is omitted in the next line
          //cout << atom.name << " number=" << atom.number << " type=" << atom.type << " cpos=" << atom.cpos << " ijk=" << aurostd::joinWDelimiter(atom.ijk,",") << " dist=" << AtomDist(a.atoms[i],atom) <<  " dist_check=" << aurostd::modulus(a.atoms[i].cpos-atom.cpos) << endl; // ijk should yield the unit cell within which the neighbor is found with respect to the central unit cell
          if (atom.name == "O"){
            neighbors_count+=1;
          } else if (atom.name != "O"){  
            // implement check whether each nearest neighbor is O, otherwise throw warning; 
            warning=1;
            cout << "WARNING: Not all nearest neighbors within the distance tolerance are oxygen, there is also " << atom.name << endl;
          }
        }
      }
      if (warning){
        cout << "WARNING: Not all nearest neighbors within the distance tolerance are oxygen!" << endl;
      }
      num_neighbors[i]=neighbors_count; // zero-based counting as for cutoffs array above
      if(LDEBUG){
        cout << "number of oxygen nearest neighbors within " << tolerance << " Ang tolerance of " << a.atoms[i].name << " (ATOM[" << i << "]): " << num_neighbors[i] << endl;
        cout << endl;
      }
    } else if (a.atoms[i].name == "O"){ // identify per- and superoxides by O-O bond length
      num_neighbors[i]=0; // set number of neighbors for O atoms in structure to zero since per- and superoxide O-O bonds will be stored in the following variables
      for(uint j=0,jsize=neigh_mat[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by the cutoff (cutoffs_max)
        const _atom& atom=neigh_mat[i][j];
        if (atom.name == "O"){
          if (0.1 < AtomDist(a.atoms[i],atom) && AtomDist(a.atoms[i],atom) <= 1.2 ){ // distance must be larger than 0.1 Ang to savely exclude the oxygen anion itself having distance zero to itself; if O-O bond is shorter than in O2 molecule (approx. 1.21 Ang) the result of the structural relaxation is most likely wrong 
            cout << "THE DETERMINED OXYGEN-OXYGEN BOND LENGTH IS SHORTER THAN IN THE O2 MOLECULE; CHECK YOUR STRUCTURE! THE O-O BOND LENGTH IS: " << AtomDist(a.atoms[i],atom) << " Ang." << endl;
            correctable=0;
          } else if (1.2 < AtomDist(a.atoms[i],atom) && AtomDist(a.atoms[i],atom) <= superox_cutoff ){
            if(LDEBUG){
              cout << "WARNING: This should be a superoxide; the O-O bond length is: " << AtomDist(a.atoms[i],atom) << " Ang." << endl;
            }
            superox_count+=1;
            superox_indices[i]=1;
          } else if (superox_cutoff < AtomDist(a.atoms[i],atom) && AtomDist(a.atoms[i],atom) <= perox_cutoff ){
            if(LDEBUG){
              cout << "WARNING: This should be a peroxide; the O-O bond length is: " << AtomDist(a.atoms[i],atom) << " Ang." << endl;
            }
            perox_count+=1;
            perox_indices[i]=1;
          }
        }
      }
    }
  }
  num_perox_bonds=perox_count/2; // needs to be divided by two due to double counting of O-O bonds when looping over all atoms of the structure
  num_superox_bonds=superox_count/2;
  // print out the number of per- & superoxide O-O bonds;
  if (num_perox_bonds > 0){
    cout << "WARNING: This should be a peroxide!" << endl;
    cout << "Number of peroxide O-O bonds in cell: " << num_perox_bonds << endl;
    cout << endl;
  } else if (num_superox_bonds > 0) {
    cout << "WARNING: This should be a superoxide!" << endl;
    cout << "Number of superoxide O-O bonds in cell: " << num_superox_bonds << endl;
    cout << endl;
  }

  // DETERMINE CORRECTION FOR EACH CATION: OXIDATION STATE DEPENDENT PART OF COORECTION ##############################################################################################################################################################
  // try to assign oxidation numbers from Bader charges only if no input oxidation numbers are provided
  if(!ox_nums_provided){ 
    if(LDEBUG){
      cout << "DETERMINATION OF OXIDATION NUMBERS FROM BADER CHARGES:" << endl;
    }
    // for every cation in the structure assign oxidation numbers
    // in this part  functional dependent settings will be primarily evaluated with "if, else if" conditions since the calculation (especially for getting the Bader charges) could have only be done with one functional and not more than one
    // oxidation number should be maybe not functional dependent since when results for a specific calculation are corrected, this should be either a PBE or LDA or SCAN calculation and not more than one; if only a structure (POSCAR) is provided as input to the correction scheme the Bader charge based determination of oxidation numbers might anyway be skipped and the oxidation numbers might even be supplied by the user; hence there should be basically no case where corrections for multiple functionals are needed based on the Bader charges (of a specific functional); for a given structure there should also be only one correct assignment of oxidation numbers
    //double oxidation_states[a.atoms.size()]; // to check whether sum over all atoms equals zero; otherwise correction assignment might be wrong; double because for superoxides O ox. number is -0.5
    string system_name;
    // check whether aflow.in exists
    ifstream aflow_in_check("aflow.in"); // Create an input file stream.
    if (aflow_in_check){
      // determine oxidation state of each (cat-)ion and store in array Bader_charges
      // get system name to identify Bader charges file and functional to distinguish corrections to be loaded from aflow.in
      ifstream aflow_in; // Create an input file stream.
      aflow_in.open("aflow.in"); // Use it to read from a file
      string line_a;
      pbe=1;
      lda=0;
      scan=0;
      while (getline( aflow_in, line_a)){
        if (line_a.find("[AFLOW]SYSTEM=") != string::npos){ // string::npos is returned if string is not found
          system_name= line_a.substr(14,line_a.size()-14);
          // remove from system name potential carriage return characters that might be inserted when sending aflow.in files as email attachments; this characters lead to problems with the later string addition causing weird string substitution
          system_name=aurostd::RemoveSubString(system_name,"\r");
        }
        if (line_a.find("=potpaw_LDA") != string::npos){ // string::npos is returned if string is not found
          lda=1; pbe=0; scan=0;
        } else if (line_a.find("METAGGA=SCAN") != string::npos){
          scan=1; pbe=0; lda=0;
        }
      }
      if(LDEBUG){
        cout << "PBE: " << pbe << endl;
        cout << "LDA: " << lda << endl;
        cout << "SCAN: " << scan << endl;
      }
      aflow_in.close();
      // if functional determined from aflow.in is different from the ones given by the input options, throw warning that oxidation numbers are only determined on the basis of a specific functional
      //if(!functionals.empty()){
        if(pbe != pbe_0 || lda != lda_0 || scan != scan_0){
          if(pbe){
            cout << "WARNING: The oxidation numbers are only determined on the basis of a PBE calculation." << endl;
          } else if(lda){
            cout << "WARNING: The oxidation numbers are only determined on the basis of an LDA calculation." << endl;
          } else if(scan){
            cout << "WARNING: The oxidation numbers are only determined on the basis of a SCAN calculation." << endl;
          }
        }
      //}
                        
      // check whether Bader file exists
      // next five lines to get system_name characters in hexadecimal used to debug the problem when carriage return characters werre added when sending aflow.in file as email attachment 
      //for (uint i =0; i < system_name.size(); i++) {
        //  long long int x = system_name[i];
        //  std::cout << "0x" << std::hex << std::setw(2) << std::setfill('0') << x << " ";
      //}
      //std::cout << std::endl;
      static string Bader_file_name;
      string checkname1=system_name + "_abader.out";
      string checkname2=system_name + "_abader.out.xz";
      string checkname3=system_name + "_abader.out.bz2";
      string checkname4=system_name + "_abader.out.gz";
      ifstream Bader_file_check1(checkname1.c_str()); // Create an input file stream.
      ifstream Bader_file_check2(checkname2.c_str()); // Create an input file stream.
      ifstream Bader_file_check3(checkname3.c_str()); // Create an input file stream.
      ifstream Bader_file_check4(checkname4.c_str()); // Create an input file stream.
      if (Bader_file_check1 || Bader_file_check2 || Bader_file_check3 || Bader_file_check4 ){
        // determine Bader charges from Bader file; O Bader charges will be included although they might not be used
        double Bader_charges[a.atoms.size()];
        string zipping_method;
        if (Bader_file_check2){
          std::system("xz -d *_abader.out.xz"); // execute UNIX command to unzip Bader charge file
          zipping_method="xz";
        } else if (Bader_file_check3){
          std::system("bzip2 -d *_abader.out.bz2"); // execute UNIX command to unzip Bader charge file
          zipping_method="bzip2";
        } else if (Bader_file_check4){
          std::system("gunzip *_abader.out.gz"); // execute UNIX command to unzip Bader charge file
          zipping_method="gzip";
        }
        string Bader_file_extension="_abader.out";
        Bader_file_name= system_name + Bader_file_extension;
        ifstream Bader_file; 
        Bader_file.open(Bader_file_name.c_str()); // c_str() gives a char* representation (string as a C "string", more correctly a pointer to the first character) rather than a char
        string line;
        uint i=0;
        while (getline( Bader_file, line )){
          if (line[38] == '-' ){
            Bader_charges[i]= strtof((line.substr(38,7)).c_str(),0);
            i++;
          } else if (line[38] == ' ') {
            Bader_charges[i]= strtof((line.substr(39,6)).c_str(),0);
            i++;
          }
        }
        Bader_file.close();
        for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
          if(LDEBUG){
            cout << "Bader_charges[" << k << "]: " << Bader_charges[k] << endl;
          }
        }
        if(zipping_method != ""){
          string zipping_command=zipping_method + " *_abader.out";
          std::system(zipping_command.c_str()); // execute UNIX command to rezip Bader charge file
        }
                                
        //correctable=1; // first assuming that formation engery of system IS correctable; will be set to zero (not correctable) if, for any atom, no correction can be identified
        //new implementation determining oxidation number from Bader charges first
        for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (a.atoms[i].name != "O"){
            string Bader_templ_line;
            if ( get_Bader_templates(a.atoms[i].name) == "") {
              correctable=0; 
              cout << "VERY BAD NEWS: There is no correction for " << a.atoms[i].name << " (ATOM[" << i << "])" << " since this species was not included in the set for deducing corrections!"  << endl;
            } else {
              Bader_templ_line=get_Bader_templates(a.atoms[i].name);
              if(LDEBUG){
                cout << "Bader templates: " << Bader_templ_line << endl;
              }
              uint num_ox_states;
              num_ox_states= strtof((Bader_templ_line.substr(0,1)).c_str(),0); // the Bader charges for cations are positive and it is hence not necessary to check for minus signs as for the corrections
              if(LDEBUG){
                cout << "number of oxidation states for which Bader charges are available: " << num_ox_states << endl;
              }
              double Bader_template=0; // Bader cation charge from the binary oxide from which the correction was deduced
              double Bader_tolerance=0.26; // tolerance with which the Bader charge of each atom is allowed to deviate from any of the template values to still assign the correction (and oxidation number)
              double Bader_deviation=0.5; // deviation of the Bader charge of the cation from the template value; initial value safely larger than Bader tolerance for identifying oxidation number so that if value is not changed no correction is made
              double Bader_deviation_0=Bader_deviation; // initial Bader deviation for later check whether any corrections were found, i. e. Bader deviation was changed
              double Bader_deviation_min=17; // initialize with high value; used later to output smallest Bader deviation if no correction according to the Bader tolerance was identified
              for(uint n=0;n<num_ox_states;n++){ //loop over all oxidation number for which Bader charges are available and read the Bader charges from the respective positions
                if (pbe){
                  Bader_template= strtof((Bader_templ_line.substr(9+n*32,6)).c_str(),0); // the Bader charges for cations are positive and it is hence not necessary to check for minus signs as for the corrections
                  if(LDEBUG){
                    cout << "Bader_template PBE: " << Bader_template << endl;
                  }
                } else if (lda) {
                  Bader_template= strtof((Bader_templ_line.substr(18+n*32,6)).c_str(),0); // the Bader charges for cations are positive and it is hence not necessary to check for minus signs as for the corrections
                  if(LDEBUG){
                    cout << "Bader_template LDA: " << Bader_template << endl;
                  }
                } else if (scan) {
                  Bader_template= strtof((Bader_templ_line.substr(27+n*32,6)).c_str(),0); // the Bader charges for cations are positive and it is hence not necessary to check for minus signs as for the corrections
                  if(LDEBUG){
                    cout << "Bader_template SCAN: " << Bader_template << endl;
                  }
                }
                if ( abs(Bader_template-Bader_charges[i]) < Bader_tolerance && abs(Bader_template-Bader_charges[i]) < Bader_deviation ){ // must be compared to Bader_deviation to load correction for Bader charge closest to template and not overloaded by other corrections for later found cases
                  Bader_deviation= abs(Bader_template-Bader_charges[i]);
                  if(LDEBUG){
                    cout << "Bader_deviation: " << Bader_deviation << endl;
                  }
                  oxidation_states[i]= strtof((Bader_templ_line.substr(4+n*32,2)).c_str(),0);
                } else { // only if element is found but no corrections can be assigned because above conditions are not met, update Bader_deviation_min
                  if ( abs(Bader_template-Bader_charges[i]) < Bader_deviation_min ) {
                    Bader_deviation_min= abs(Bader_template-Bader_charges[i]); 
                  }
                }
              }
              if ( Bader_deviation == Bader_deviation_0 ){ // only if no oxidation state was found, i. e. Bader deviation was not changed but there are corrections (Bader charges) for this species, which was tested above by whether the function get_Bader_templates returns an empty string. The following warning (after handling the special cases) will be displayed  
                // FIRST point to deal with special cases known from (ternary) oxides; other cases will be dealt with below when looking at the oxidation numbers
                // since it can happen that the Bader charge for W in the compound (especially for W+6) is too far away from the Bader charge template, here W+6 is assumed and later the sum over the oxidation states will still be checked
                if (a.atoms[i].name == "W") {
                  oxidation_states[i]=+6;
                // since it can happen that the Bader charge for Pb in the compound (especially for Pb+2) is too far away from the Bader charge template, here Pb+2 is assumed and later the sum over the oxidation states will still be checked
                } else if (a.atoms[i].name == "Pb") {
                  oxidation_states[i]=+2;
                // since it can happen that the Bader charge for Ag in the compound (especially for Ag+1) is too far away from the Bader charge template, here Ag+1 is assumed and later the sum over the oxidation states will still be checked
                } else if (a.atoms[i].name == "Ag") {
                  oxidation_states[i]=+1;
                } else {
                  correctable=0; 
                  cout << "BAD NEWS: The correction for " << a.atoms[i].name << " (ATOM[" << i << "])" << " cannot be identified from the Bader charges!"  << endl;
                  cout << "The deviation of the Bader charge from the closest tested template value is: " << Bader_deviation_min << " electrons. This is larger than the tolerance of: " << Bader_tolerance  << " electrons." << endl;
                  // list all oxidation states of the element for which corrections are available
                  string ox_nums_avail="";
                  string separator=", ";
                  for(uint n=0;n<num_ox_states;n++){ //loop over all oxidation number for which Bader charges are available and read the Bader charges from the respective positions
                    if (n<num_ox_states-1){
                      ox_nums_avail+= Bader_templ_line.substr(4+n*32,2) + separator ; // the Bader charges for cations are positive and it is hence not necessary to check for minus signs as for the corrections
                    } else if (n==num_ox_states-1){
                      ox_nums_avail+= Bader_templ_line.substr(4+n*32,2); // the Bader charges for cations are positive and it is hence not necessary to check for minus signs as for the corrections
                    }
                  }
                  cout << "Available oxidation states for " << a.atoms[i].name << " coordinated by oxygen are: " << ox_nums_avail << endl;
                  cout << "If a correction for the desired oxidation state is available but it is just not correctly determined from the Bader charges," << endl;
                  cout << "you might want to consider supplying the oxidation numbers manually by using the option --oxidation_numbers=." << endl;
                  cout << endl;
                }
              }
            }
          } else if (a.atoms[i].name == "O") {
            oxidation_states[i]=-2; // oxidation numbers for O are assumed to be -2 and are corrected below if it is a per- or superoxide O atom as identified from the structural analysis above
            if (num_perox_bonds > 0){
              for (uint j=0,jsize=a.atoms.size();j<jsize;j++){
                if (perox_indices[j]==1 && j == i){
                  oxidation_states[i]=-1;
                }
              }
            }
            if (num_superox_bonds > 0){
              for (uint j=0,jsize=a.atoms.size();j<jsize;j++){
                if (superox_indices[j]==1 && j == i){
                  oxidation_states[i]=-0.5;
                }
              }
            }
          }
        }
        // check whether sum of oxidation numbers is equal zero and try to correct for knwon special cases in case it is not
        // SECOND point to deal with special cases known from (binary+ternary) oxides
        if (correctable){
          // print oxidation numbers
          cout << "OXIDATION NUMBERS:" << endl;
          for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
            cout << "Oxidation state of " << a.atoms[k].name << " (atom[" << k << "]): " << oxidation_states[k] << endl;
          }
          // calculate sum of oxidation numbers
          double oxidation_sum=0; // double because for superoxides O ox. number is -0.5
          for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
            oxidation_sum+= oxidation_states[k];
          }
          cout << "CHECK whether this is what you are expecting!" << endl;
          cout << "Sum over all oxidation numbers is (should be ZERO): " << oxidation_sum << endl;
          // Ti-O Magneli phases need to be treated specially since oxdiation numbers are not recognized appropriately; there are always 2xTi+3 per formula unit and the rest is Ti+4; fortunately, both Ti+3 and Ti+4 have 6 Ti-O bonds, hence one only needs to know how many ions of the respective oxidation state there are, not which one is which
          if ( a.species.size() == 2 && a.species[0] == "O" && a.species[1] == "Ti" ) {
            if(LDEBUG){
              cout << "Ti-O system, Magneli for Ti_nO_(2n-1), i.e. Ti-O ratio= " << 3.0/5 << ", " << 4.0/7 << ", " << 5.0/9 << ", " << 6.0/11 << ", " << 7.0/13 << ", " << 8.0/15 << ", " << 9.0/17 << ", " << 10.0/19 << "..." << endl;
            }
            double amount_O=0;
            double amount_Ti=0;
            for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
              if (a.atoms[i].name == "O"){
                amount_O+=1;
              } else if (a.atoms[i].name == "Ti"){
                amount_Ti+=1;
              }
            }
            if(LDEBUG){
              cout << "number of O ions= " << amount_O << endl;
              cout << "number of Ti ions= " << amount_Ti << endl;
            }
            double Ti_O_ratio=amount_Ti/amount_O;
            if(LDEBUG){
              cout << "ratio of Ti/O= " << Ti_O_ratio << endl;
            }
            double num_formula_units_in_cell;
            // check for Magneli composition Ti_(n)O_(2n-1)
            double n;
            bool magneli=0;
            for(n=3;n<101;n++){
              //cout << "n/(2*n-1)= " << n/(2*n-1) << endl;
              if ( Ti_O_ratio == n/(2*n-1) ){
                cout << "n= " << n << " Magneli composition Ti_nO_(2n-1)" << endl;
                magneli=1;
                num_formula_units_in_cell=amount_Ti/n;
                cout << "number of formula units in cell: " << num_formula_units_in_cell << endl;
              }
            }
            if ( magneli == 0){
              if(LDEBUG){
                cout << "Not a Magneli composition." << endl;
              }
            }
            if (magneli){
              for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
                if (a.atoms[i].name == "Ti"){
                  // taking the first Ti as +3 (without knowing whether they are actually +3) works only because for the Magneli phases both Ti+3 and Ti+4 have both always 6 Ti-O bonds
                  if ( i < amount_O+num_formula_units_in_cell*2 ){
                    oxidation_states[i]=+3;
                    cout << "setting oxidation state to Ti+3 " << endl;
                  } else {
                    oxidation_states[i]=+4;
                    cout << "setting oxidation state to Ti+4 " << endl;
                  }
                }
              }
              // print oxidation numbers
              for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                cout << "Oxidation state of " << a.atoms[k].name << " (atom[" << k << "]): " << oxidation_states[k] << endl;
              }
              // calculate sum of oxidation numbers
              oxidation_sum=0; // double because for superoxides O ox. number is -0.5
              for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                oxidation_sum+= oxidation_states[k];
              }
              cout << "Sum over all oxidation numbers is (should be ZERO): " << oxidation_sum << endl;
              // system should not be regarded correctable if sum over oxidation states is not zero
              if (oxidation_sum != 0){
                correctable =0;
                cout << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
              }
            }
          }
          // for Fe3O4 in inverse spinel structure the oxidation states are not identified properly via the Bader charges. According to Wikipedia the Fe2+ ions are octahedrally coordinated and the Fe3+ ions are evenly distributed between the octahedral and tetrahedral sites and one hence needs per formula unit 6x the Fe2+ correction and 4+6=10x the Fe3+ correction, see:
          // https://en.wikipedia.org/wiki/Iron(II,III)_oxide
          // for Co3O4 and Mn3O4 this is different; in the normal spinel structure Co2+/Mn2+ occupies only tetrahedral sites while Co3+/Mn3+ occupies octahedral sites; the correction hence needs to be different but at the present stage there are no corrections for Co3+ and Mn3+, see:
          // https://en.wikipedia.org/wiki/Cobalt(II,III)_oxide
          // https://en.wikipedia.org/wiki/Manganese(II,III)_oxide
          if ( a.species.size() == 2 && a.species[0] == "Fe" && a.species[1] == "O" ) {
            double amount_Fe=0;
            double amount_O=0;
            for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
              if (a.atoms[i].name == "Fe"){
                amount_Fe+=1;
              } else if (a.atoms[i].name == "O"){
                amount_O+=1;
              }
            }
            if(LDEBUG){
              cout << "number of Fe ions= " << amount_Fe << endl;
              cout << "number of O ions= " << amount_O << endl;
            }
            double Fe_O_ratio=amount_Fe/amount_O;
            if ( Fe_O_ratio == 0.75 ){
              cout << "Fe3O4 with ratio of Fe/O= " << Fe_O_ratio << endl;
              double num_formula_units_in_cell;
              num_formula_units_in_cell=a.atoms.size()/7; // 7 for Fe3O4
              cout << "number of formula units in cell: " << num_formula_units_in_cell << endl;
              for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
                if (a.atoms[i].name == "Fe"){
                  // taking the first Fe as +2 (without knowing whether they are actually +2) works only because the number of bonds for them will be adjusted to 6 (octahedral) as needed for Fe3O4 disregarding the actual number of Fe-O bonds (this is a hack since I know how many bonds there should be for each ion type)
                  if ( i < 1*num_formula_units_in_cell ){  // 1 Fe2+ ions per formula unit * number of formula units
                    oxidation_states[i]=+2;
                    cout << "setting oxidation state to Fe+2 " << endl;
                    num_neighbors[i]=6; // for Fe3O4 the Fe2+ ions are 6-fold coordinated by oxygen according to Wikipedia
                    cout << "Modified number of neighbors of Fe[" << i << "] taken as Fe+2: " << num_neighbors[i] << endl;
                  } else {
                    oxidation_states[i]=+3;
                    cout << "setting oxidation state to Fe+3 " << endl;
                    num_neighbors[i]=5; // for Fe3O4 the Fe3+ ions are evenly 6- and 4-fold, so on average 5-fold (set here as a hack), coordinated by oxygen according to Wikipedia
                    cout << "Modified number of neighbors of Fe[" << i << "] taken as Fe+3 (average between even 6- and 4-fold coordination): " << num_neighbors[i] << endl;
                  }
                }
              }
              // print oxidation numbers
              for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                cout << "Oxidation state of " << a.atoms[k].name << " (atom[" << k << "]): " << oxidation_states[k] << endl;
              }
              // calculate sum of oxidation numbers
              oxidation_sum=0; // double because for superoxides O ox. number is -0.5
              for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                oxidation_sum+= oxidation_states[k];
              }
              cout << "Sum over all oxidation numbers is (should be ZERO): " << oxidation_sum << endl;
              // system should not be regarded correctable if sum over oxidation states is not zero
              if (oxidation_sum != 0){
                correctable =0;
                cout << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
              }
            }
          }
          // for MnMoO4 the oxidation numbers are determined to be +4 for both Mn and Mo from the Bader charges; it should be Mn+2 & Mo+6; however since the sum of the falsely determined oxidation numbers is accidentally 0 it needs to be corrected individually
          if ( a.species.size() == 3 && a.species[0] == "Mn" && a.species[1] == "Mo" && a.species[2] == "O" ) {
            double amount_Mn=0;
            double amount_Mo=0;
            double amount_O=0;
            for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
              if (a.atoms[i].name == "Mn"){
                amount_Mn+=1;
              } else if (a.atoms[i].name == "Mo"){
                amount_Mo+=1;
              } else if (a.atoms[i].name == "O"){
                amount_O+=1;
              }
            }
            if(LDEBUG){
              cout << "number of Mn ions= " << amount_Mn << endl;
              cout << "number of Mo ions= " << amount_Mo << endl;
              cout << "number of O ions= " << amount_O << endl;
            }
            if (amount_Mn/amount_Mo == 1 && amount_Mn/amount_O == 0.25 && amount_Mo/amount_O == 0.25) {
              cout << "MnMoO4 special treatment since sum over oxdiation states is zero but individual oxidation numbers are wrong!!!" << endl;
              for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
                if (a.atoms[i].name == "Mn"){
                  oxidation_states[i]=+2;
                  cout << "setting oxidation state to Mn+2 " << endl;
                }
                if (a.atoms[i].name == "Mo"){
                  oxidation_states[i]=+6;
                  cout << "setting oxidation state to Mo+6 " << endl;
                }
              }
              // print oxidation numbers
              for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                cout << "Oxidation state of " << a.atoms[k].name << " (atom[" << k << "]): " << oxidation_states[k] << endl;
              }
              // calculate sum of oxidation numbers
              oxidation_sum=0; // double because for superoxides O ox. number is -0.5
              for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                oxidation_sum+= oxidation_states[k];
              }
              cout << "Sum over all oxidation numbers is (should be ZERO): " << oxidation_sum << endl;
              // system should not be regarded correctable if sum over oxidation states is not zero
              if (oxidation_sum != 0){
                correctable =0;
                cout << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
              }
            }
          }
          if (lda) {
            // for Ca2Fe2O5 and CaFe2O4 for LDA the oxidation numbers of Fe are not correctly determined to be Fe+3 but are partly Fe+2 which will be corrected here
            if ( a.species.size() == 3 && a.species[0] == "Ca" && a.species[1] == "Fe" && a.species[2] == "O" ) {
              double amount_Ca=0;
              double amount_Fe=0;
              double amount_O=0;
              for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
                if (a.atoms[i].name == "Ca"){
                  amount_Ca+=1;
                } else if (a.atoms[i].name == "Fe"){
                  amount_Fe+=1;
                } else if (a.atoms[i].name == "O"){
                  amount_O+=1;
                }
              }
              if(LDEBUG){
                cout << "number of Ca ions= " << amount_Ca << endl;
                cout << "number of Fe ions= " << amount_Fe << endl;
                cout << "number of O ions= " << amount_O << endl;
              }
              //making sure it is Ca2Fe2O5
              if (amount_Ca/amount_Fe == 1 && amount_Ca/amount_O == 0.4 && amount_Fe/amount_O == 0.4) {
                cout << "Ca2Fe2O5 special treatment for LDA since oxidation numbers for Fe, which should be Fe+3, are not correctly determined from Bader charges for all Fe!!!" << endl;
                for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
                  if (a.atoms[i].name == "Fe"){
                    oxidation_states[i]=+3;
                    cout << "setting oxidation state to Fe+3 " << endl;
                  }
                }
                // print oxidation numbers
                for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                  cout << "Oxidation state of " << a.atoms[k].name << " (atom[" << k << "]): " << oxidation_states[k] << endl;
                }
                // calculate sum of oxidation numbers
                oxidation_sum=0; // double because for superoxides O ox. number is -0.5
                for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                  oxidation_sum+= oxidation_states[k];
                }
                cout << "Sum over all oxidation numbers is (should be ZERO): " << oxidation_sum << endl;
                // system should not be regarded correctable if sum over oxidation states is not zero
                if (oxidation_sum != 0){
                  correctable =0;
                  cout << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
                }
              //making sure it is CaFe2O4
              } else if (amount_Ca/amount_Fe == 0.5 && amount_Ca/amount_O == 0.25 && amount_Fe/amount_O == 0.5) {
                cout << "CaFe2O4 special treatment for LDA since oxidation numbers for Fe, which should be Fe+3, are not correctly determined from Bader charges for all Fe!!!" << endl;
                for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
                  if (a.atoms[i].name == "Fe"){
                    oxidation_states[i]=+3;
                    cout << "setting oxidation state to Fe+3 " << endl;
                  }
                }
                // print oxidation numbers
                for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                  cout << "Oxidation state of " << a.atoms[k].name << " (atom[" << k << "]): " << oxidation_states[k] << endl;
                }
                // calculate sum of oxidation numbers
                oxidation_sum=0; // double because for superoxides O ox. number is -0.5
                for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                  oxidation_sum+= oxidation_states[k];
                }
                cout << "Sum over all oxidation numbers is (should be ZERO): " << oxidation_sum << endl;
                // system should not be regarded correctable if sum over oxidation states is not zero
                if (oxidation_sum != 0){
                  correctable =0;
                  cout << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
                }
              }
            }
            // for FeTiO3 for LDA the oxidation numbers of Ti are not correctly determined to be Ti+4 and using below fixes would modify both the Ti AND the Fe oxidation numbers resulting again in non-zero oxidation number sum which is fixed here
            if ( a.species.size() == 3 && a.species[0] == "Fe" && a.species[1] == "O" && a.species[2] == "Ti" ) {
              double amount_Fe=0;
              double amount_O=0;
              double amount_Ti=0;
              for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
                if (a.atoms[i].name == "Fe"){
                  amount_Fe+=1;
                } else if (a.atoms[i].name == "O"){
                  amount_O+=1;
                } else if (a.atoms[i].name == "Ti"){
                  amount_Ti+=1;
                }
              }
              if(LDEBUG){
                cout << "number of Fe ions= " << amount_Fe << endl;
                cout << "number of O ions= " << amount_O << endl;
                cout << "number of Ti ions= " << amount_Ti << endl;
              }
              //making sure it is FeTiO3
              if (amount_Fe/amount_Ti == 1 && amount_Fe/amount_O == 1.0/3 && amount_Ti/amount_O == 1.0/3) {
                cout << "FeTiO3 special treatment for LDA since oxidation numbers for Ti, which should be Ti+4, are not correctly determined from Bader charges and using other fixing would also change Fe oxidation numbers!!!" << endl;
                for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
                  if (a.atoms[i].name == "Ti"){
                    oxidation_states[i]=+4;
                    cout << "setting oxidation state to Ti+4 " << endl;
                  }
                }
                // print oxidation numbers
                for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                  cout << "Oxidation state of " << a.atoms[k].name << " (atom[" << k << "]): " << oxidation_states[k] << endl;
                }
                // calculate sum of oxidation numbers
                oxidation_sum=0; // double because for superoxides O ox. number is -0.5
                for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
                  oxidation_sum+= oxidation_states[k];
                }
                cout << "Sum over all oxidation numbers is (should be ZERO): " << oxidation_sum << endl;
                // system should not be regarded correctable if sum over oxidation states is not zero
                if (oxidation_sum != 0){
                  correctable =0;
                  cout << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
                }
              }
            }
          }
          // try to repair wrong oxidation states for special cases known (from oxides)
          if (oxidation_sum != 0){
            cout << "The sum over all oxidation numbers for all atoms of the system is NOT zero, trying to repair that based on known problematic cases (Ti, V, Fe). This might or might not work." << endl;
            // repairing only by considering non mixed valence oxides first
            for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
              if (a.atoms[i].name == "Ti"){
                if ( oxidation_states[i] == 3){
                  oxidation_states[i]=+4;
                  cout << "setting oxidation state to Ti+4 " << endl;
                }
              }
              if (a.atoms[i].name == "Fe"){
                if ( oxidation_states[i] == 2){
                  oxidation_states[i]=+3;
                  cout << "setting oxidation state to Fe+3 " << endl;
                } else if ( oxidation_states[i] == 3){
                  oxidation_states[i]=+2;
                  cout << "setting oxidation state to Fe+2 " << endl;
                }
              }
              if (a.atoms[i].name == "V"){
                if ( oxidation_states[i] == 4){
                  oxidation_states[i]=+5;
                  cout << "setting oxidation state to V+5 " << endl;
                }
              }
            }
            // print oxidation numbers
            cout << "OXIDATION NUMBERS:" << endl;
            for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
              cout << "Oxidation state of " << a.atoms[k].name << " (atom[" << k << "]): " << oxidation_states[k] << endl;
            }
            // calculate sum of oxidation numbers
            oxidation_sum=0; // double because for superoxides O ox. number is -0.5
            for (uint k=0,ksize=a.atoms.size();k<ksize;k++){
              oxidation_sum+= oxidation_states[k];
            }
            cout << "CHECK whether this is what you are expecting!" << endl;
            cout << "Sum over all oxidation numbers is (should be ZERO): " << oxidation_sum << endl;
            // system should not be regarded correctable if sum over oxidation states is not zero
            if (oxidation_sum != 0){
              correctable =0;
              cout << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
            }
          }
          cout << endl; // empty line to separate formation enthalpy values from the rest
        }
      } else {
        correctable=0;
        cout << "Bader file " << checkname1 << " (or xz, bz2, gz version) not found. A Bader file is required to determine the oxidation numbers." << endl;
        cout << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=" << endl;
      }
    } else {
      correctable=0;
      cout << "aflow.in file not found. An aflow.in file is required to identify the functional and to find the Bader charges file to determine the oxidation numbers." << endl;
      cout << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=" << endl;
    }
  }

  //new implementation assigning corrections after oxidation numbers are determined
  //here for every functional there will be a separate "if" evaluation, since when only a structure (+ oxidation numbers) are provided as input one might want to get the corrections for more than one functional
  double corrections_298_PBE[a.atoms.size()]; // for predicting formation enthalpies at 298.15K
  double corrections_0_PBE[a.atoms.size()]; // for predicting formation enthalpies at 0K
  double corrections_298_LDA[a.atoms.size()]; // for predicting formation enthalpies at 298.15K
  double corrections_0_LDA[a.atoms.size()]; // for predicting formation enthalpies at 0K
  double corrections_298_SCAN[a.atoms.size()]; // for predicting formation enthalpies at 298.15K
  double corrections_0_SCAN[a.atoms.size()]; // for predicting formation enthalpies at 0K
  double perox_cor_298_PBE;
  double superox_cor_298_PBE;
  double perox_cor_0_PBE;
  double superox_cor_0_PBE;
  double perox_cor_298_LDA;
  double superox_cor_298_LDA;
  double perox_cor_0_LDA;
  double superox_cor_0_LDA;
  double perox_cor_298_SCAN;
  double superox_cor_298_SCAN;
  double perox_cor_0_SCAN;
  double superox_cor_0_SCAN;
  if (correctable){
    if(LDEBUG){
      cout << "ASSIGNMENT OF CORRECTIONS:" << endl;
    }
    for(uint i=0,isize=a.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (a.atoms[i].name != "O"){
        string corrections_line;
        std::ostringstream ox_strs;
        ox_strs << oxidation_states[i];
        std::string ox_str = ox_strs.str();
        if ( get_corrections(a.atoms[i].name + "_+" + ox_str + "_ox") == "") {
          correctable=0; 
          cout << "BAD NEWS: No correction available for " << a.atoms[i].name << " (ATOM[" << i << "])" << " in oxidation state " << oxidation_states[i] << " when coordinated by oxygen." << endl;
          string Bader_templ_line;
          if ( get_Bader_templates(a.atoms[i].name) == "") {
            cout << "Currently no corrections available for " << a.atoms[i].name << " when coordinated by oxygen." << endl;
          } else {
            // list all oxidation states of the element for which corrections are available
            Bader_templ_line=get_Bader_templates(a.atoms[i].name);
            uint num_ox_states;
            string ox_nums_avail="";
            string separator=", ";
            num_ox_states= strtof((Bader_templ_line.substr(0,1)).c_str(),0); // the Bader charges for cations are positive and it is hence not necessary to check for minus signs as for the corrections
            for(uint n=0;n<num_ox_states;n++){ //loop over all oxidation number for which Bader charges are available and read the Bader charges from the respective positions
              if (n<num_ox_states-1){
                ox_nums_avail+= Bader_templ_line.substr(4+n*32,2) + separator ; // the Bader charges for cations are positive and it is hence not necessary to check for minus signs as for the corrections
              } else if (n==num_ox_states-1){
                ox_nums_avail+= Bader_templ_line.substr(4+n*32,2); // the Bader charges for cations are positive and it is hence not necessary to check for minus signs as for the corrections
              }
            }
            cout << "Available oxidation states for " << a.atoms[i].name << " coordinated by oxygen are: " << ox_nums_avail << endl;
          }
          cout << endl;
        } else {
          corrections_line=get_corrections(a.atoms[i].name + "_+" + ox_str + "_ox");
          if(LDEBUG){
            cout << "corrections line: " << corrections_line << endl;
          }
          if (pbe || pbe_0){
          // load PBE corrections
            if (corrections_line[47] == '-'){ // check whether correction is negative
              corrections_298_PBE[i]= strtof((corrections_line.substr(47,8)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 298.15K: " << corrections_298_PBE[i] << " eV/bond" << endl;
              }
            } else if (corrections_line[47] == ' '){
              corrections_298_PBE[i]= strtof((corrections_line.substr(48,7)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 298.15K: " << corrections_298_PBE[i] << " eV/bond" << endl;
              }
            }
            if (corrections_line[66] == '-'){
              corrections_0_PBE[i]= strtof((corrections_line.substr(66,8)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 0K: " << corrections_0_PBE[i] << " eV/bond" << endl;
              }
            } else if (corrections_line[66] == ' '){
              corrections_0_PBE[i]= strtof((corrections_line.substr(67,7)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 0K: " << corrections_0_PBE[i] << " eV/bond" << endl;
              }
            }
          } 
          if (lda || lda_0){
          // load LDA corrections
            if (corrections_line[93] == '-'){
              corrections_298_LDA[i]= strtof((corrections_line.substr(93,8)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 298.15K: " << corrections_298_LDA[i] << " eV/bond" << endl;
              }
            } else if (corrections_line[93] == ' '){
              corrections_298_LDA[i]= strtof((corrections_line.substr(94,7)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 298.15K: " << corrections_298_LDA[i] << " eV/bond" << endl;
              }
            }
            if (corrections_line[112] == '-'){
              corrections_0_LDA[i]= strtof((corrections_line.substr(112,8)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 0K: " << corrections_0_LDA[i] << " eV/bond" << endl;
              }
            } else if (corrections_line[112] == ' '){
              corrections_0_LDA[i]= strtof((corrections_line.substr(113,7)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 0K: " << corrections_0_LDA[i] << " eV/bond" << endl;
              }
            }
          } 
          if (scan || scan_0){
          // load SCAN corrections
            if (corrections_line[139] == '-'){
              corrections_298_SCAN[i]= strtof((corrections_line.substr(139,8)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 298.15K: " << corrections_298_SCAN[i] << " eV/bond" << endl;
              }
            } else if (corrections_line[139] == ' '){
              corrections_298_SCAN[i]= strtof((corrections_line.substr(140,7)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 298.15K: " << corrections_298_SCAN[i] << " eV/bond" << endl;
              }
            }
            if (corrections_line[159] == '-'){
              corrections_0_SCAN[i]= strtof((corrections_line.substr(159,8)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 0K: " << corrections_0_SCAN[i] << " eV/bond" << endl;
              }
            } else if (corrections_line[159] == ' '){
              corrections_0_SCAN[i]= strtof((corrections_line.substr(160,7)).c_str(),0);
              if(LDEBUG){
                cout << "correction for " << a.atoms[i].name << " for 0K: " << corrections_0_SCAN[i] << " eV/bond" << endl;
              }
            }
          }
        }
      } else if (a.atoms[i].name == "O") {
        // set corrections for O atoms to zero since only for cations number of bonds with O are needed; per- & superoxides will be dealt with differently
        if (pbe || pbe_0){
          corrections_298_PBE[i]=0;
          if(LDEBUG){
            cout << "correction for oxygen (atom[" << i << "]) for 298.15K: " << corrections_298_PBE[i] << " eV/bond" << endl;
          }
          corrections_0_PBE[i]=0;
          if(LDEBUG){
            cout << "correction for oxygen (atom[" << i << "]) for 0K: " << corrections_0_PBE[i] << " eV/bond" << endl;
          }
        } 
        if (lda || lda_0){
          corrections_298_LDA[i]=0;
          if(LDEBUG){
            cout << "correction for oxygen (atom[" << i << "]) for 298.15K: " << corrections_298_LDA[i] << " eV/bond" << endl;
          }
          corrections_0_LDA[i]=0;
          if(LDEBUG){
            cout << "correction for oxygen (atom[" << i << "]) for 0K: " << corrections_0_LDA[i] << " eV/bond" << endl;
          }
        } 
        if (scan || scan_0){
          corrections_298_SCAN[i]=0;
          if(LDEBUG){
            cout << "correction for oxygen (atom[" << i << "]) for 298.15K: " << corrections_298_SCAN[i] << " eV/bond" << endl;
          }
          corrections_0_SCAN[i]=0;
          if(LDEBUG){
            cout << "correction for oxygen (atom[" << i << "]) for 0K: " << corrections_0_SCAN[i] << " eV/bond" << endl;
          }
        }
      }
    }
    // load per- & superox. corrections if needed;
    if (num_perox_bonds > 0){
      string corrections_line;
      corrections_line=get_corrections("O2_-2_ox");
      if(LDEBUG){
        cout << "corrections line peroxides: " << corrections_line << endl;
      }
      if (pbe || pbe_0){
        if (corrections_line[47] == '-'){ 
          perox_cor_298_PBE= strtof((corrections_line.substr(47,8)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 298.15K: " << perox_cor_298_PBE << " eV/bond" << endl;
          }
        } else if (corrections_line[47] == ' '){
          perox_cor_298_PBE= strtof((corrections_line.substr(48,7)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 298.15K: " << perox_cor_298_PBE << " eV/bond" << endl;
          }
        }
        if (corrections_line[66] == '-'){ 
          perox_cor_0_PBE= strtof((corrections_line.substr(66,8)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 0K: " << perox_cor_0_PBE << " eV/bond" << endl;
          }
        } else if (corrections_line[66] == ' '){
          perox_cor_0_PBE= strtof((corrections_line.substr(67,7)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 0K: " << perox_cor_0_PBE << " eV/bond" << endl;
          }
        }
      } 
      if (lda || lda_0){
        if (corrections_line[93] == '-'){
          perox_cor_298_LDA= strtof((corrections_line.substr(93,8)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 298.15K: " << perox_cor_298_LDA << " eV/bond" << endl;
          }
        } else if (corrections_line[93] == ' '){
          perox_cor_298_LDA= strtof((corrections_line.substr(94,7)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 298.15K: " << perox_cor_298_LDA << " eV/bond" << endl;
          }
        }
        if (corrections_line[112] == '-'){
          perox_cor_0_LDA= strtof((corrections_line.substr(112,8)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 0K: " << perox_cor_0_LDA << " eV/bond" << endl;
          }
        } else if (corrections_line[112] == ' '){
          perox_cor_0_LDA= strtof((corrections_line.substr(113,7)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 0K: " << perox_cor_0_LDA << " eV/bond" << endl;
          }
        }
      } 
      if (scan || scan_0){
        if (corrections_line[139] == '-'){
          perox_cor_298_SCAN= strtof((corrections_line.substr(139,8)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 298.15K: " << perox_cor_298_SCAN << " eV/bond" << endl;
          }
        } else if (corrections_line[139] == ' '){
          perox_cor_298_SCAN= strtof((corrections_line.substr(140,7)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 298.15K: " << perox_cor_298_SCAN << " eV/bond" << endl;
          }
        }
        if (corrections_line[159] == '-'){
          perox_cor_0_SCAN= strtof((corrections_line.substr(159,8)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 0K: " << perox_cor_0_SCAN << " eV/bond" << endl;
          }
        } else if (corrections_line[159] == ' '){
          perox_cor_0_SCAN= strtof((corrections_line.substr(160,7)).c_str(),0);
          if(LDEBUG){
            cout << "Peroxide correction for 0K: " << perox_cor_0_SCAN << " eV/bond" << endl;
          }
        }
      }
    } 
    if (num_superox_bonds > 0) {
      string corrections_line;
      corrections_line=get_corrections("O2_-1_ox");
      if(LDEBUG){
        cout << "corrections line superoxides: " << corrections_line << endl;
      }
      if (pbe || pbe_0){
        if (corrections_line[47] == '-'){
          superox_cor_298_PBE= strtof((corrections_line.substr(47,8)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 298.15K: " << superox_cor_298_PBE << " eV/bond" << endl;
          }
        } else if (corrections_line[47] == ' '){
          superox_cor_298_PBE= strtof((corrections_line.substr(48,7)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 298.15K: " << superox_cor_298_PBE << " eV/bond" << endl;
          }
        }
        if (corrections_line[66] == '-'){
          superox_cor_0_PBE= strtof((corrections_line.substr(66,8)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 0K: " << superox_cor_0_PBE << " eV/bond" << endl;
          }
        } else if (corrections_line[66] == ' '){
          superox_cor_0_PBE= strtof((corrections_line.substr(67,7)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 0K: " << superox_cor_0_PBE << " eV/bond" << endl;
          }
        }
      } 
      if (lda || lda_0){
        if (corrections_line[93] == '-'){
          superox_cor_298_LDA= strtof((corrections_line.substr(93,8)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 298.15K: " << superox_cor_298_LDA << " eV/bond" << endl;
          }
        } else if (corrections_line[93] == ' '){
          superox_cor_298_LDA= strtof((corrections_line.substr(94,7)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 298.15K: " << superox_cor_298_LDA << " eV/bond" << endl;
          }
        }
        if (corrections_line[112] == '-'){
          superox_cor_0_LDA= strtof((corrections_line.substr(112,8)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 0K: " << superox_cor_0_LDA << " eV/bond" << endl;
          }
        } else if (corrections_line[112] == ' '){
          superox_cor_0_LDA= strtof((corrections_line.substr(113,7)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 0K: " << superox_cor_0_LDA << " eV/bond" << endl;
          }
        }
      } 
      if (scan || scan_0){
        if (corrections_line[139] == '-'){
          superox_cor_298_SCAN= strtof((corrections_line.substr(139,8)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 298.15K: " << superox_cor_298_SCAN << " eV/bond" << endl;
          }
        } else if (corrections_line[139] == ' '){
          superox_cor_298_SCAN= strtof((corrections_line.substr(140,7)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 298.15K: " << superox_cor_298_SCAN << " eV/bond" << endl;
          }
        }
        if (corrections_line[159] == '-'){
          superox_cor_0_SCAN= strtof((corrections_line.substr(159,8)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 0K: " << superox_cor_0_SCAN << " eV/bond" << endl;
          }
        } else if (corrections_line[159] == ' '){
          superox_cor_0_SCAN= strtof((corrections_line.substr(160,7)).c_str(),0);
          if(LDEBUG){
            cout << "Superoxide correction for 0K: " << superox_cor_0_SCAN << " eV/bond" << endl;
          }
        }
      }
    }
  }
  if(LDEBUG){
    cout << endl;
  }


  // CALCULATE PLAIN DFT AND CCE CORRECTED FORMATION ENTHALPIES AT 298.15 AND 0K ##############################################################################################################################################################
  // calculate CCE correction if the system is correctable
  if (correctable){
    //CCE correction should be for entire cell and normalized to per atom
    double cce_correction_part_298_PBE=0;
    double cce_correction_part_0_PBE=0;
    double cce_correction_part_298_LDA=0;
    double cce_correction_part_0_LDA=0;
    double cce_correction_part_298_SCAN=0;
    double cce_correction_part_0_SCAN=0;
    double cce_form_energy_cell_298_PBE;
    double cce_form_energy_cell_0_PBE;
    double cce_form_energy_atom_298_PBE;
    double cce_form_energy_atom_0_PBE;
    double cce_form_energy_cell_298_LDA;
    double cce_form_energy_cell_0_LDA;
    double cce_form_energy_atom_298_LDA;
    double cce_form_energy_atom_0_LDA;
    double cce_form_energy_cell_298_SCAN;
    double cce_form_energy_cell_0_SCAN;
    double cce_form_energy_atom_298_SCAN;
    double cce_form_energy_atom_0_SCAN;
    // do correction only using cations although the number of neighbors and corrections for O atoms should be set to zero; the last factor is for the normalization per formula unit, if one wants to have it per atom directly, one should divide by a.atoms.size() instead
    for(uint i=0,isize=a.atoms.size();i<isize;i++){
      if (a.atoms[i].name != "O"){
        if (pbe || pbe_0){
          cce_correction_part_298_PBE += (num_neighbors[i]*corrections_298_PBE[i]) ;
          cce_correction_part_0_PBE += (num_neighbors[i]*corrections_0_PBE[i]) ;
        } 
        if (lda || lda_0){
          cce_correction_part_298_LDA += (num_neighbors[i]*corrections_298_LDA[i]) ;
          cce_correction_part_0_LDA += (num_neighbors[i]*corrections_0_LDA[i]) ;
        } 
        if (scan || scan_0){
          cce_correction_part_298_SCAN += (num_neighbors[i]*corrections_298_SCAN[i]) ;
          cce_correction_part_0_SCAN += (num_neighbors[i]*corrections_0_SCAN[i]) ;
        }
      }
    }
    // add per- and superoxide correction if needed; the last factor is for the normalization per formula unit, if one wants to have it per atom directly, one should divide by a.atoms.size() instead
    if (num_perox_bonds > 0){
      if (pbe || pbe_0){
        cce_correction_part_298_PBE += (num_perox_bonds * perox_cor_298_PBE) ;
        cout << "Peroxide correction for 298.15K per cell: " << (num_perox_bonds * perox_cor_298_PBE) << " eV" << endl;
        cce_correction_part_0_PBE += (num_perox_bonds * perox_cor_0_PBE) ;
        cout << "Peroxide correction for 0K per cell: " << (num_perox_bonds * perox_cor_0_PBE) << " eV" << endl;
        cout << endl;
      } 
      if (lda || lda_0){
        cce_correction_part_298_LDA += (num_perox_bonds * perox_cor_298_LDA) ;
        cout << "Peroxide correction for 298.15K per cell: " << (num_perox_bonds * perox_cor_298_LDA) << " eV" << endl;
        cce_correction_part_0_LDA += (num_perox_bonds * perox_cor_0_LDA) ;
        cout << "Peroxide correction for 0K per cell: " << (num_perox_bonds * perox_cor_0_LDA) << " eV" << endl;
        cout << endl;
      } 
      if (scan || scan_0){
        cce_correction_part_298_SCAN += (num_perox_bonds * perox_cor_298_SCAN) ;
        cout << "Peroxide correction for 298.15K per cell: " << (num_perox_bonds * perox_cor_298_SCAN) << " eV" << endl;
        cce_correction_part_0_SCAN += (num_perox_bonds * perox_cor_0_SCAN) ;
        cout << "Peroxide correction for 0K per cell: " << (num_perox_bonds * perox_cor_0_SCAN) << " eV" << endl;
        cout << endl;
      }
    }
    if (num_superox_bonds > 0){
      if (pbe || pbe_0){
        cce_correction_part_298_PBE += (num_superox_bonds * superox_cor_298_PBE) ;
        cout << "Superoxide correction for 298.15K per cell: " << (num_superox_bonds * superox_cor_298_PBE) << " eV" << endl;
        cce_correction_part_0_PBE += (num_superox_bonds * superox_cor_0_PBE) ;
        cout << "Superoxide correction for 0K per cell: " << (num_superox_bonds * superox_cor_0_PBE) << " eV" << endl ;
        cout << endl;
      } 
      if (lda || lda_0){
        cce_correction_part_298_LDA += (num_superox_bonds * superox_cor_298_LDA) ;
        cout << "Superoxide correction for 298.15K per cell: " << (num_superox_bonds * superox_cor_298_LDA) << " eV" << endl;
        cce_correction_part_0_LDA += (num_superox_bonds * superox_cor_0_LDA) ;
        cout << "Superoxide correction for 0K per cell: " << (num_superox_bonds * superox_cor_0_LDA) << " eV" << endl ;
        cout << endl;
      } 
      if (scan || scan_0){
        cce_correction_part_298_SCAN += (num_superox_bonds * superox_cor_298_SCAN) ;
        cout << "Superoxide correction for 298.15K per cell: " << (num_superox_bonds * superox_cor_298_SCAN) << " eV" << endl;
        cce_correction_part_0_SCAN += (num_superox_bonds * superox_cor_0_SCAN) ;
        cout << "Superoxide correction for 0K per cell: " << (num_superox_bonds * superox_cor_0_SCAN) << " eV" << endl ;
        cout << endl;
      }
    }

    // print out CCE correction for functionals selected
    cout << "CCE CORRECTIONS:" << endl;
    if (pbe || pbe_0){
      cout << "CCE@PBE correction (to be subtracted from PBE formation energy) for 298.15K in eV per cell: " << std::setprecision(3) << std::fixed << std::setw(7) << cce_correction_part_298_PBE << " eV" << endl;
      cout << "CCE@PBE correction (to be subtracted from PBE formation energy) for 0K in eV per cell: " << std::setprecision(3) << std::fixed << std::setw(12) << cce_correction_part_0_PBE << " eV" << endl;
      cout << "CCE@PBE correction (to be subtracted from PBE formation energy) for 298.15K in eV per atom: " << std::setprecision(3) << std::fixed << std::setw(7) << cce_correction_part_298_PBE/a.atoms.size() << " eV" << endl;
      cout << "CCE@PBE correction (to be subtracted from PBE formation energy) for 0K in eV per atom: " << std::setprecision(3) << std::fixed << std::setw(12) << cce_correction_part_0_PBE/a.atoms.size() << " eV" << endl;
      cout << endl;
    } 
    if (lda || lda_0){
      cout << "CCE@LDA correction (to be subtracted from LDA formation energy) for 298.15K in eV per cell: " << std::setprecision(3) << std::fixed << std::setw(7) << cce_correction_part_298_LDA << " eV" << endl;
      cout << "CCE@LDA correction (to be subtracted from LDA formation energy) for 0K in eV per cell: " << std::setprecision(3) << std::fixed << std::setw(12) << cce_correction_part_0_LDA << " eV" << endl;
      cout << "CCE@LDA correction (to be subtracted from LDA formation energy) for 298.15K in eV per atom: " << std::setprecision(3) << std::fixed << std::setw(7) << cce_correction_part_298_LDA/a.atoms.size() << " eV" << endl;
      cout << "CCE@LDA correction (to be subtracted from LDA formation energy) for 0K in eV per atom: " << std::setprecision(3) << std::fixed << std::setw(12) << cce_correction_part_0_LDA/a.atoms.size() << " eV" << endl;
      cout << endl;
    } 
    if (scan || scan_0){
      cout << "CCE@SCAN correction (to be subtracted from SCAN formation energy) for 298.15K in eV per cell: " << std::setprecision(3) << std::fixed << std::setw(7) << cce_correction_part_298_SCAN << " eV" << endl;
      cout << "CCE@SCAN correction (to be subtracted from SCAN formation energy) for 0K in eV per cell: " << std::setprecision(3) << std::fixed << std::setw(12) << cce_correction_part_0_SCAN << " eV" << endl;
      cout << "CCE@SCAN correction (to be subtracted from SCAN formation energy) for 298.15K in eV per atom: " << std::setprecision(3) << std::fixed << std::setw(7) << cce_correction_part_298_SCAN/a.atoms.size() << " eV" << endl;
      cout << "CCE@SCAN correction (to be subtracted from SCAN formation energy) for 0K in eV per atom: " << std::setprecision(3) << std::fixed << std::setw(12) << cce_correction_part_0_SCAN/a.atoms.size() << " eV" << endl;
      cout << endl;
    }

    //calculated corrected DFT formation enthalpies if precalculated DFT formation energies are provided
    if(dft_energies.size()!=0){ // only when precalculated DFT formatio energies are provided calculated corrected values
      cout << "CCE FORMATION ENTHALPIES:" << endl;
      for(uint k=0,ksize=vfunctionals.size();k<ksize;k++){
        if (vfunctionals[k] == "PBE"){
          cce_form_energy_cell_298_PBE = dft_energies[k] - cce_correction_part_298_PBE ;
          cce_form_energy_cell_0_PBE = dft_energies[k] - cce_correction_part_0_PBE ;
          cce_form_energy_atom_298_PBE = cce_form_energy_cell_298_PBE/a.atoms.size();
          cce_form_energy_atom_0_PBE = cce_form_energy_cell_0_PBE/a.atoms.size();
          cout << "CCE@PBE formation enthalpy of the system at 298.15K per cell: " << std::setw(8) << cce_form_energy_cell_298_PBE << " eV" << endl;
          cout << "CCE@PBE formation enthalpy of the system at 0K per cell: " << std::setw(13)  << cce_form_energy_cell_0_PBE << " eV" << endl;
          cout << "CCE@PBE formation enthalpy of the system at 298.15K per atom: " << std::setw(8) << cce_form_energy_atom_298_PBE << " eV" << endl;
          cout << "CCE@PBE formation enthalpy of the system at 0K per atom: " << std::setw(13) << cce_form_energy_atom_0_PBE << " eV" << endl;
          cout << endl;
        } else if (vfunctionals[k] == "LDA"){
          cce_form_energy_cell_298_LDA = dft_energies[k] - cce_correction_part_298_LDA ;
          cce_form_energy_cell_0_LDA = dft_energies[k] - cce_correction_part_0_LDA ;
          cce_form_energy_atom_298_LDA = cce_form_energy_cell_298_LDA/a.atoms.size();
          cce_form_energy_atom_0_LDA = cce_form_energy_cell_0_LDA/a.atoms.size();
          cout << "CCE@LDA formation enthalpy of the system at 298.15K per cell: " << std::setw(8) << cce_form_energy_cell_298_LDA << " eV" << endl;
          cout << "CCE@LDA formation enthalpy of the system at 0K per cell: " << std::setw(13) << cce_form_energy_cell_0_LDA << " eV" << endl;
          cout << "CCE@LDA formation enthalpy of the system at 298.15K per atom: " << std::setw(8) << cce_form_energy_atom_298_LDA << " eV" << endl;
          cout << "CCE@LDA formation enthalpy of the system at 0K per atom: " << std::setw(13) << cce_form_energy_atom_0_LDA << " eV" << endl;
          cout << endl;
        } else if (vfunctionals[k] == "SCAN"){
          cce_form_energy_cell_298_SCAN = dft_energies[k] - cce_correction_part_298_SCAN ;
          cce_form_energy_cell_0_SCAN = dft_energies[k] - cce_correction_part_0_SCAN ;
          cce_form_energy_atom_298_SCAN = cce_form_energy_cell_298_SCAN/a.atoms.size();
          cce_form_energy_atom_0_SCAN = cce_form_energy_cell_0_SCAN/a.atoms.size();
          cout << "CCE@SCAN formation enthalpy of the system at 298.15K per cell: " << std::setw(8) << cce_form_energy_cell_298_SCAN << " eV" << endl;
          cout << "CCE@SCAN formation enthalpy of the system at 0K per cell: " << std::setw(13) << cce_form_energy_cell_0_SCAN << " eV" << endl;
          cout << "CCE@SCAN formation enthalpy of the system at 298.15K per atom: " << std::setw(8) << cce_form_energy_atom_298_SCAN << " eV" << endl;
          cout << "CCE@SCAN formation enthalpy of the system at 0K per atom: " << std::setw(13) << cce_form_energy_atom_0_SCAN << " eV" << endl;
          cout << endl;
        }
      }
    }
  cout << "#########################################################################################################" << endl;
  cout << "When you obtain results using the CCE methodology and/or this implementation," << endl; 
  cout << "please cite the following article:" << endl;
  cout << "Friedrich et al., Coordination corrected ab initio formation enthalpies, npj Comput. Mater. 5, 59 (2019);" << endl; 
  cout << "https://doi.org/10.1038/s41524-019-0192-1" << endl;
  cout << "#########################################################################################################" << endl;
  cout << endl;
  }
  }
} // namespace pflow

#endif // _AFLOW_CCE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2019              *
// *                                                                         *
// ***************************************************************************
