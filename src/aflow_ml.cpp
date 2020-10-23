// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu

#ifndef _AFLOW_ML_CPP_
#define _AFLOW_ML_CPP_

#include "aflow.h"
#include "aflow_cce.h"

namespace aflowML {
  void writeCCECSV() {
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy=XPID+"aflowML::writeCCECSV():";

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    uint ielement=0,ioxidation=0,i=0,ipp=0;
    string correction_line="",input_pre="",input="";
    xelement::xelement xel;
    //cce::get_corrections_line_O("Al_+3_N");
    vector<vector<string> > vlines;
    vector<string> vitems_pre,vitems,vtokens;
    vlines.resize(2); //N,O
    
    vector<string> vheaders;
    vheaders.push_back("symbol_cation");
    vheaders.push_back("period_cation");
    vheaders.push_back("group_cation");
    vheaders.push_back("electronegativity_Allen_cation");
    vheaders.push_back("oxidation_cation");
    vheaders.push_back("Z_cation");
    vheaders.push_back("groundstate_energy_PBE_cation");
    vheaders.push_back("EATOM_PBE_cation");
    //
    vheaders.push_back("symbol_anion");
    vheaders.push_back("period_anion");
    vheaders.push_back("group_anion");
    vheaders.push_back("electronegativity_Allen_anion");
    vheaders.push_back("Z_anion");
    vheaders.push_back("groundstate_energy_PBE_anion");
    vheaders.push_back("EATOM_PBE_anion");
    //
    vheaders.push_back("PBE_298.15K");
    vheaders.push_back("PBE_0K");
    vheaders.push_back("LDA_298.15K");
    vheaders.push_back("LDA_0K");
    vheaders.push_back("SCAN_298.15K");
    vheaders.push_back("SCAN_0K");
    vheaders.push_back("PBE+U_298.15K");
    vheaders.push_back("PBE+U_0K");
    vheaders.push_back("exp_298.15K");
    vheaders.push_back("M-X_bonds");

    for(i=0;i<vlines.size();i++){
      vlines[i].push_back(aurostd::joinWDelimiter(vheaders,","));
    }
    
    int precision=5;
    string species_pp="";
    bool found_pp=false;
    for(ielement=0;ielement<100;ielement++){
      for(ioxidation=0;ioxidation<10;ioxidation++){
        xel.populate(ielement);
        input_pre=xel.symbol;
        input_pre+="_+"+aurostd::utype2string(ioxidation);
        if(LDEBUG){cerr << soliloquy << " input=" << input_pre << endl;}
        vitems_pre.clear();
        //cation
        vitems_pre.push_back(xel.symbol); //symbol
        //[get from xpseudo instead]vitems_pre.push_back(aurostd::utype2string(xel.Z,precision)); //Z
        vitems_pre.push_back(aurostd::utype2string(xel.Period,precision)); //Period
        vitems_pre.push_back(aurostd::utype2string(xel.Group,precision)); //Group
        vitems_pre.push_back(aurostd::utype2string(xel.electronegativity_Allen,precision)); //electronegativity_Allen
        vitems_pre.push_back(aurostd::utype2string(ioxidation,precision)); //ioxidation
        try{species_pp=AVASP_Get_PseudoPotential_PAW_PBE(xel.symbol);}
        catch(aurostd::xerror& excpt){continue;}
        found_pp=false;
        for(ipp=0;ipp<vxpseudopotential.size()&&found_pp==false;ipp++) {
          if(vxpseudopotential[ipp].species_pp_type[0]=="PAW_PBE" && vxpseudopotential[ipp].species_pp[0]==species_pp){
            found_pp=true;
            vitems_pre.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_Z[0],precision)); //Z
            vitems_pre.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_pp_groundstate_energy[0],precision)); //groundstate_energy
            vitems_pre.push_back(aurostd::utype2string(vxpseudopotential[ipp].vEATOM[0],precision)); //EATOM
          }
        }

        //N
        input=input_pre+"_N";
        correction_line=cce::get_corrections_line_N(input);
        if(!correction_line.empty()){
          vitems.clear();for(i=0;i<vitems_pre.size();i++){vitems.push_back(vitems_pre[i]);}
          //anion
          xel.populate("N");
          vitems.push_back(xel.symbol); //symbol
          //[get from xpseudo instead]vitems.push_back(aurostd::utype2string(xel.Z,precision)); //Z
          vitems.push_back(aurostd::utype2string(xel.Period,precision)); //Period
          vitems.push_back(aurostd::utype2string(xel.Group,precision)); //Group
          vitems.push_back(aurostd::utype2string(xel.electronegativity_Allen,precision)); //electronegativity_Allen

          try{species_pp=AVASP_Get_PseudoPotential_PAW_PBE(xel.symbol);}
          catch(aurostd::xerror& excpt){continue;}
          found_pp=false;
          for(ipp=0;ipp<vxpseudopotential.size()&&found_pp==false;ipp++) {
            if(vxpseudopotential[ipp].species_pp_type[0]=="PAW_PBE" && vxpseudopotential[ipp].species_pp[0]==species_pp){
              found_pp=true;
              if(vxpseudopotential[ipp].species_pp_groundstate_structure[0]!="diatom"){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"pp_groundstate_structure for "+xel.symbol+" is NOT diatom",_INPUT_ILLEGAL_);}
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_Z[0],precision)); //Z
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_pp_groundstate_energy[0],precision)); //groundstate_energy
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].vEATOM[0],precision)); //EATOM
            }
          }

          //
          if(LDEBUG){cerr << soliloquy << " correction_line=\"" << correction_line << "\"" << endl;}
          aurostd::string2tokens(correction_line,vtokens," ");
          if(vtokens.size()<12){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vtokens.size()<11",_FILE_CORRUPT_);}
          for(i=2;i<=11;i++){vitems.push_back(vtokens[i]);}
          vlines[0].push_back(aurostd::joinWDelimiter(vitems,","));
        }

        //O
        input=input_pre+"_O";
        correction_line=cce::get_corrections_line_O(input);
        if(!correction_line.empty()){
          vitems.clear();for(i=0;i<vitems_pre.size();i++){vitems.push_back(vitems_pre[i]);}
          //anion
          xel.populate("O");
          vitems.push_back(xel.symbol); //symbol
          //[get from xpseudo instead]vitems.push_back(aurostd::utype2string(xel.Z,precision)); //Z
          vitems.push_back(aurostd::utype2string(xel.Period,precision)); //Period
          vitems.push_back(aurostd::utype2string(xel.Group,precision)); //Group
          vitems.push_back(aurostd::utype2string(xel.electronegativity_Allen,precision)); //electronegativity_Allen

          try{species_pp=AVASP_Get_PseudoPotential_PAW_PBE(xel.symbol);}
          catch(aurostd::xerror& excpt){continue;}
          found_pp=false;
          for(ipp=0;ipp<vxpseudopotential.size()&&found_pp==false;ipp++) {
            if(vxpseudopotential[ipp].species_pp_type[0]=="PAW_PBE" && vxpseudopotential[ipp].species_pp[0]==species_pp){
              found_pp=true;
              if(vxpseudopotential[ipp].species_pp_groundstate_structure[0]!="diatom"){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"pp_groundstate_structure for "+xel.symbol+" is NOT diatom",_INPUT_ILLEGAL_);}
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_Z[0],precision)); //Z
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].species_pp_groundstate_energy[0],precision)); //groundstate_energy
              vitems.push_back(aurostd::utype2string(vxpseudopotential[ipp].vEATOM[0],precision)); //EATOM
            }
          }

          //
          if(LDEBUG){cerr << soliloquy << " correction_line=\"" << correction_line << "\"" << endl;}
          aurostd::string2tokens(correction_line,vtokens," ");
          if(vtokens.size()<12){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"vtokens.size()<11",_FILE_CORRUPT_);}
          for(i=2;i<=11;i++){vitems.push_back(vtokens[i]);}
          vlines[1].push_back(aurostd::joinWDelimiter(vitems,","));
        }
      }
    }

    stringstream file;
    //N
    aurostd::StringstreamClean(file);
    for(uint i=0;i<vlines[0].size();i++){file << vlines[0][i] << endl;}
    aurostd::stringstream2file(file,"cce_data_N.csv");
    //O
    aurostd::StringstreamClean(file);
    for(uint i=0;i<vlines[1].size();i++){file << vlines[1][i] << endl;}
    aurostd::stringstream2file(file,"cce_data_O.csv");
    //total
    aurostd::StringstreamClean(file);
    for(uint i=0;i<vlines[0].size();i++){file << vlines[0][i] << endl;}
    for(uint i=1;i<vlines[1].size();i++){file << vlines[1][i] << endl;} //skip header
    aurostd::stringstream2file(file,"cce_data_NO.csv");

  }
}

#endif  // _AFLOW_ML_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
