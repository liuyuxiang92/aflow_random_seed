// ***************************************************************************
// *                                                                         *
// *              Aflow COREY OSES - Duke University 2003-2020               *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses 2020

#ifndef _AUROSTD_XPARSER_CPP_
#define _AUROSTD_XPARSER_CPP_

#ifndef _AUROSTD_XPARSER_H_
#include "aurostd_xparser.h"
#endif

namespace aurostd {

  void VASP_PseudoPotential_CleanName_InPlace(string& species,bool capital_letters_only) { //CO20190712
    //WARNING: to anyone adding to this list, BE CAREFUL to avoid adding entries that contain capital letters
    //they must be added to CAPITAL_LETTERS_PP_LIST in aurostd.h
    //these pp suffixes cause problems when parsing compounds (capital letters)

    vector<string> vCAPITAL_LETTERS_PP;
    aurostd::string2tokens(CAPITAL_LETTERS_PP_LIST,vCAPITAL_LETTERS_PP,",");
    for(uint i=0;i<vCAPITAL_LETTERS_PP.size();i++){//capital letter ones to watch out for when parsing compounds
      aurostd::RemoveSubStringInPlace(species,vCAPITAL_LETTERS_PP[i]);
    }

    if(capital_letters_only==false){
      aurostd::RemoveSubStringInPlace(species,"_old");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Si_h_old
      aurostd::RemoveSubStringInPlace(species,".old");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Mg_pv.old
      aurostd::RemoveSubStringInPlace(species,"_vnew");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Pd_vnew
      aurostd::RemoveSubStringInPlace(species,"_new2");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Ti_sv_new2
      aurostd::RemoveSubStringInPlace(species,"_new");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Au_new

      aurostd::RemoveSubStringInPlace(species,"_pvf");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Cu_pvf
      aurostd::RemoveSubStringInPlace(species,"_rel");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Pb_d_rel
      aurostd::RemoveSubStringInPlace(species,"_ref");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Ge_d_GW_ref
      aurostd::RemoveSubStringInPlace(species,"_local");  //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/C_local
      aurostd::RemoveSubStringInPlace(species,"_nopc");  //CO20190712 - potpaw_LDA/potpaw_PBE.20100505/Si_nopc
      aurostd::RemoveSubStringInPlace(species,".nrel");  //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Ga_pv_GW.nrel
      aurostd::RemoveSubStringInPlace(species,"_nr");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/C_h_nr
      aurostd::RemoveSubStringInPlace(species,"_nc");  //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/H_nc_GW
      aurostd::RemoveSubStringInPlace(species,"_n");  //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/As_GW_n
      aurostd::RemoveSubStringInPlace(species,"_parsv");  //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Mg_pv_parsv_GW
      aurostd::RemoveSubStringInPlace(species,"_sv2");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Li_sv2
      aurostd::RemoveSubStringInPlace(species,"_sv");
      aurostd::RemoveSubStringInPlace(species,"_vs"); //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/N_vs
      aurostd::RemoveSubStringInPlace(species,"_pv");
      aurostd::RemoveSubStringInPlace(species,"_dr");  //CO20190712 - BEFORE _d //potpaw_LDA/potpaw_LDA.20100505/Pb_dr
      aurostd::RemoveSubStringInPlace(species,"_d3");  //CO20190712 - BEFORE _d //potpaw_PBE/potpaw_PBE.20100506/Ge_d3
      aurostd::RemoveSubStringInPlace(species,"_d2");  //CO20190712 - BEFORE _d //potpaw_LDA/potpaw_LDA.05May2010/As_d2_GW
      aurostd::RemoveSubStringInPlace(species,"_d");
      aurostd::RemoveSubStringInPlace(species,"_soft");  //CO20190712 - BEFORE _s
      aurostd::RemoveSubStringInPlace(species,"_s");
      //[CO20190712 - OBSOLETE really _n and _2]aurostd::RemoveSubStringInPlace(species,"_2_n");
      aurostd::RemoveSubStringInPlace(species,"_h");
      aurostd::RemoveSubStringInPlace(species,"_f");  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Cu_f
      aurostd::RemoveSubStringInPlace(species,"_af"); //CO20191110 - SHACHAR aflow pp 

      aurostd::RemoveSubStringInPlace(species,"_1");
      aurostd::RemoveSubStringInPlace(species,"_2");
      aurostd::RemoveSubStringInPlace(species,"_3");

      aurostd::RemoveSubStringInPlace(species,"1.75"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H1.75
      aurostd::RemoveSubStringInPlace(species,"1.66"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H1.66
      aurostd::RemoveSubStringInPlace(species,"1.33"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H1.33
      aurostd::RemoveSubStringInPlace(species,"1.25"); //CO20190712 - before all other decimal numbers
      aurostd::RemoveSubStringInPlace(species,"1.5"); //CO20190712 - potpaw_PBE/potpaw_PBE.06May2010/H1.5
      aurostd::RemoveSubStringInPlace(species,".75");  //CO20190712 - before 0.5
      aurostd::RemoveSubStringInPlace(species,".25");  //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.25
      aurostd::RemoveSubStringInPlace(species,".66"); //CO20190712 - potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.66
      aurostd::RemoveSubStringInPlace(species,".33"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H.33
      aurostd::RemoveSubStringInPlace(species,".42"); //CO20190712 - potpaw_PBE.54/potpaw_PBE.54.04Sep2015/H.42
      aurostd::RemoveSubStringInPlace(species,".58"); //CO20190712 - before 0.5 //potpaw_LDA.52/potpaw_LDA.52.19Apr2012/H.58
      aurostd::RemoveSubStringInPlace(species,".5");

      aurostd::RemoveSubStringInPlace(species,"+1");
      aurostd::RemoveSubStringInPlace(species,"+3");
      aurostd::RemoveSubStringInPlace(species,"+5");
      aurostd::RemoveSubStringInPlace(species,"+7");
      aurostd::RemoveSubStringInPlace(species,"-1");
      aurostd::RemoveSubStringInPlace(species,"-3");
      aurostd::RemoveSubStringInPlace(species,"-5");
      aurostd::RemoveSubStringInPlace(species,"-7");

      //from AFLOW.org database
      aurostd::RemoveSubStringInPlace(species,"pot_LDA/");
      aurostd::RemoveSubStringInPlace(species,"pot_GGA/");
      aurostd::RemoveSubStringInPlace(species,"pot_PBE/");
      aurostd::RemoveSubStringInPlace(species,"potpaw_LDA/");
      aurostd::RemoveSubStringInPlace(species,"potpaw_GGA/");
      aurostd::RemoveSubStringInPlace(species,"potpaw_PBE/");
      aurostd::RemoveSubStringInPlace(species,"potpaw_LDA.54/");
      aurostd::RemoveSubStringInPlace(species,"potpaw_PBE.54/");

      //general database
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POT_LDA+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POT_GGA+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POT_PBE+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_GGA+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_LDA_KIN+"/");
      aurostd::RemoveSubStringInPlace(species,DEFAULT_VASP_POTCAR_DIR_POTPAW_PBE_KIN+"/");

      aurostd::RemoveSubStringInPlace(species,"__"); //CO20190712 - BEFORE _ - potpaw_LDA/potpaw_LDA.05May2010/Si_sv_GW__
      aurostd::RemoveSubStringInPlace(species,"_");  //CO20190712  //potpaw_LDA/potpaw_LDA.05May2010/Si_sv_GW_
    }
  }

} // namespace aurostd

#endif // _AUROSTD_XPARSER_CPP_

// **************************************************************************
// *                                                                        *
// *              Aflow COREY OSES - Duke University 2003-2020              *
// *                                                                        *
// **************************************************************************
