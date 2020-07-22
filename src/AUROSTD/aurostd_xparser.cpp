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

  //use only as a supplement for getElements(), do NOT use outside
  //this assumes a very simple Mn2Pd5, non-stoich is ok (e.g., Mn2.5Pd5)
  //no pseudo potential specification
  //no junk at the end (_ICSD_, :LDAU2, :PAW_PBE, .OLD, etc.), pre-process before
  //this is FASTER than getElements(), but not as robust for general input (specialized)
  void elementsFromCompositionString(const string& input,vector<string>& velements){vector<double> vcomposition;return elementsFromCompositionString(input,velements,vcomposition);}  //CO20190712
  template<class utype> void elementsFromCompositionString(const string& input,vector<string>& velements,vector<utype>& vcomposition){ //CO20190712
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "aurostd::getElementsFromCompositionString():";
    velements.clear();
    vcomposition.clear();  //ME20190628

    //////////////////////////////////////////////////////////////////////////////
    // START Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << soliloquy << " original input=" << input << endl;}

    //CO20180409 - running through input twice, no need, simply check at the end
    //uint numberOfElements = 0;
    //for (uint i = 0; i < input.size(); i++) {
    //  if(isupper(input[i])) {
    //    numberOfElements++;
    //  }
    //}
    //if(numberOfElements == 0) {
    //  pflow::logger(_AFLOW_FILE_NAME_, soliloquy, "Elements must be properly capitalized", FileMESSAGE, oss, _LOGGER_ERROR_);
    //  return velements;
    //}

    //////////////////////////////////////////////////////////////////////////////
    // END Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Parsing input
    //////////////////////////////////////////////////////////////////////////////

    //CO20180316 - fixed this function to be simpler, too complicated before
    string auxstr;
    for (uint i = 0; i < input.size(); i++) {
      if(isupper(input[i])) {
        auxstr.clear();
        auxstr+=input[i++];
        while (((i < input.size()) && (input[i]>='a' && input[i]<='z') )){auxstr+=input[i++];}
        i--;  //very important since we increase at the top of the loop (equivalent to i+j-1)
        //while (((i < input.size()) && isalpha(input[i]) && !isupper(input[i]))){auxstr+=input[i++];}
        //while(!(clean && !isalpha(input[i]))) //(input[i]=='_' || input[i]==':' || input[i]=='.' || isdigit(input[i]))))
        //isalpha() saves us from all issues with VASP_PseudoPotential_CleanName() except, e.g., potpaw_PBE/Na, we took care of that above
        //if(clean)
        //{ //CO20200106 - patching for auto-indenting
        //  auxstr = KBIN::VASP_PseudoPotential_CleanName(auxstr);  //fix vasp pp
        //  //CO20180409 - again, no need to run through essentially a third time, we already cut at these characters
        //  //look for bad characters and cut the string
        //  //for(uint j=1;j<auxstr.size();j++){
        //  //  if(auxstr[j]=='_' || auxstr[j]==':' || isdigit(auxstr[j])){auxstr=auxstr.substr(0,j);break;}  //fix aflow stuff like ':'
        //  //}
        //}
        if(LDEBUG) {cerr << soliloquy << " element found: " << auxstr << endl;}
        velements.push_back(auxstr);
        //ME20190628 - get composition, too
      } else if ( (input[i]>='0' && input[i]<='9') || (input[i] == '.')) {  //CO20190712 - just in case we have H.25 (not good form but try to catch anyway, never produced by aflow automatically)
        auxstr.clear();
        auxstr += input[i++];
        while ((i < input.size()) && ( (input[i]>='0' && input[i]<='9') || (input[i] == '.'))) {auxstr += input[i++];}
        i--;
        if (LDEBUG) {
          std::cerr << soliloquy << " found element count: " << auxstr << " of element " << (velements.size() - 1) << ".";
          if (vcomposition.size() != velements.size()) {
            std::cerr << " Will add ones to elements " << vcomposition.size() << " to " << (velements.size() - 2) << ".";
          }
          std::cerr << std::endl;
        }
        // Add implicit ones
        for (uint i = vcomposition.size(); i < velements.size() - 1; i++){vcomposition.push_back((utype)1.0);}
        vcomposition.push_back(aurostd::string2utype<utype>(auxstr));
      }
    }
    // Add implicit ones
    for (uint i = vcomposition.size(); i < velements.size(); i++) vcomposition.push_back((utype)1.0);
  }

  //use only as a supplement for getElements(), do NOT use outside
  //this assumes Mn_pvPt
  //no composition information
  //no junk at the end (_ICSD_, :LDAU2, :PAW_PBE, .OLD, etc.), pre-process before
  void elementsFromPPString(const string& input,vector<string>& velements,bool keep_pp){ //CO20190712
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "aurostd::getElementsFromPPString():";
    velements=getElements(input);
    if(LDEBUG){cerr << soliloquy << " velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}
    if(keep_pp==false){return;}

    //copy info into vspecies and clear velements
    vector<string> vspecies;
    for(uint i=0;i<velements.size();i++){vspecies.push_back(velements[i]);}
    velements.clear();

    //simply parse string around these elements
    string::size_type loc1=0,loc2=string::npos;
    vector<string> vCAPITAL_LETTERS_PP;
    aurostd::string2tokens(CAPITAL_LETTERS_PP_LIST,vCAPITAL_LETTERS_PP,",");
    bool found_CAPITAL_LETTERS_PP=false;
    bool found_CAPITAL_LETTERS=false;
    for(uint i=0;i<vspecies.size();i++){
      if((i+1)>=vspecies.size()){loc2=string::npos;}
      else{loc2=input.find(vspecies[i+1],loc1);}
      while(loc2!=string::npos){
        found_CAPITAL_LETTERS_PP=false;
        for(uint j=0;j<vCAPITAL_LETTERS_PP.size()&&found_CAPITAL_LETTERS_PP==false;j++){
          if((loc2-(vCAPITAL_LETTERS_PP[j].size()-1))<input.size()){continue;}
          found_CAPITAL_LETTERS=true;
          for(uint k=0;k<vCAPITAL_LETTERS_PP[j].size()&&found_CAPITAL_LETTERS==true;k++){
            if(input[loc2-k]!=vCAPITAL_LETTERS_PP[j][vCAPITAL_LETTERS_PP[j].size()-k-1]){found_CAPITAL_LETTERS=false;}
          }
          if(found_CAPITAL_LETTERS){found_CAPITAL_LETTERS_PP=true;}
        }
        if(found_CAPITAL_LETTERS_PP==false){break;}
        //[OBSOLETE] Do not pick W from _GW (tungsten)
        //[OBSOLETE]if (!( (loc2-2)<input.size() && (input[loc2] == 'W') && (input[loc2-1] == 'G') && (input[loc2-2] == '_') )){break;} //(loc2-2)<input.size() because loc2 is utype, it's always >0, loc2-2 can wrap around to a big number though
        loc2=input.find(vspecies[i],loc2+1);
      }
      if(LDEBUG){cerr << soliloquy << " loc1=" << loc1 << ", loc2=" << loc2 << endl;}
      velements.push_back(input.substr(loc1,loc2-loc1));  //loc2-loc1 because it is the distance
      loc1=loc2;
    }

    if(LDEBUG){cerr << soliloquy << " velements=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(velements,"\""),",") << endl;}

  }

  // ***************************************************************************
  // aurostd::getElements(string input,ostream&
  // oss,ofstream& FileMESSAGE)
  // ***************************************************************************
  // returns UNSORTED vector<string> from string
  vector<string> getElements(const string& input){ //CO20190712 //borrowed from XATOM_SplitAlloySpecies() //slow since we create many strings, but definitely works
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "aurostd::getElements():";
    if(LDEBUG){cerr << soliloquy << " original input=\"" << input << "\"" << endl;}
    string alloy=input;
    //[CO20190712 - no need for multiple passes anymore]for(uint i=1;i<=2;i++){alloy=KBIN::VASP_PseudoPotential_CleanName(alloy);} //be certain you clean everything, especially _GW (worst offender)
    aurostd::VASP_PseudoPotential_CleanName_InPlace(alloy); //be certain you clean everything, especially _GW (worst offender)
    aurostd::RemoveNumbersInPlace(alloy);              // remove composition
    if(LDEBUG){cerr << soliloquy << " cleaned input=\"" << alloy << "\"" << endl;}
    vector<string> vspecies;
    for(uint i=0;i<alloy.length();i++) {
      if(alloy[i]>='A' && alloy[i]<='Z') vspecies.push_back("");
      vspecies.back()+=alloy[i];
    }
    if(LDEBUG){cerr << soliloquy << " vspecies pre ASCII clean=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies,"\""),",") << endl;}
    for(uint i=0;i<vspecies.size();i++){aurostd::CleanStringASCII_InPlace(vspecies[i]);}
    if(LDEBUG){cerr << soliloquy << " vspecies post ASCII clean=" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(vspecies,"\""),",") << endl;}
    return vspecies;
  }
  vector<string> getElements(const string& input,elements_string_type e_str_type,bool clean,bool sort_elements,bool keep_pp,ostream& oss) {  // overload
    ofstream FileMESSAGE;
    return getElements(input,e_str_type,FileMESSAGE,clean,sort_elements,keep_pp,oss);
  }
  //ME20190628 - added variant that also determines the composition
  template<class utype> vector<string> getElements(const string& input,vector<utype>& vcomposition,bool clean,bool sort_elements,bool keep_pp,ostream& oss) {
    ofstream FileMESSAGE;
    return getElements(input,vcomposition,composition_string,FileMESSAGE,clean,sort_elements,keep_pp,oss);  //this gets composition_string by default, pp_string has no composition
  }
  //cannot deduce utype from this construction
  vector<string> getElements(const string& input,elements_string_type e_str_type,ofstream& FileMESSAGE,bool clean,bool sort_elements,bool keep_pp,ostream& oss) {  // overload
    vector<double> vcomposition;
    return getElements(input,vcomposition,e_str_type,FileMESSAGE,clean,sort_elements,keep_pp,oss);
  }
  template<class utype> vector<string> getElements(const string& _input,vector<utype>& vcomposition,elements_string_type e_str_type,ofstream& FileMESSAGE,bool clean,bool sort_elements,bool keep_pp,ostream& oss) { // main function
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "aurostd::getElements():";
    vector<string> velements;
    vcomposition.clear();  //ME20190628

    //////////////////////////////////////////////////////////////////////////////
    // START Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    if(LDEBUG) {cerr << soliloquy << " original input=" << _input << endl;}

    if(_input.empty()) {
      pflow::logger(_AFLOW_FILE_NAME_, soliloquy, "Empty input", FileMESSAGE, oss, _LOGGER_ERROR_);
      return velements;
    }

    string input=_input;

    if(clean && (e_str_type==composition_string || (e_str_type==pp_string && keep_pp==false))){aurostd::VASP_PseudoPotential_CleanName_InPlace(input);}  //in case we run into potpaw_PBE/Na, but only works for single elements, must be before check for isupper(input[0])

    if(!isupper(input[0])) {
      pflow::logger(_AFLOW_FILE_NAME_, soliloquy, "Elements must be properly capitalized (input="+input+")", FileMESSAGE, oss, _LOGGER_ERROR_);
      return velements;
    }

    //we have a LIB1 problem... grab first everything before :
    //this is safe, as aflow generally introduces : in prototype, e.g., :LDAU2
    //this is safe anyway because elements would be BEFORE :
    if(clean){
      //FAST
      string::size_type loc;
      //:
      loc=input.find(':');input=input.substr(0,loc);
      //_ICSD_
      loc=input.find("_ICSD_");input=input.substr(0,loc);
      //SLOW
      //vector<string> tokens;
      //aurostd::string2tokens(input,tokens,":");
      //input=tokens[0];
    }

    if(LDEBUG) {cerr << soliloquy << " checking input=" << input << endl;}

    //CO20180409 - running through input twice, no need, simply check at the end
    //uint numberOfElements = 0;
    //for (uint i = 0; i < input.size(); i++) {
    //  if(isupper(input[i])) {
    //    numberOfElements++;
    //  }
    //}
    //if(numberOfElements == 0) {
    //  pflow::logger(_AFLOW_FILE_NAME_, soliloquy, "Elements must be properly capitalized", FileMESSAGE, oss, _LOGGER_ERROR_);
    //  return velements;
    //}

    //////////////////////////////////////////////////////////////////////////////
    // END Checks for correct input by counting number of uppercase letters
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // START Parsing input
    //////////////////////////////////////////////////////////////////////////////

    if(e_str_type==composition_string){elementsFromCompositionString(input,velements,vcomposition);}
    else if(e_str_type==pp_string){elementsFromPPString(input,velements,keep_pp);}
    else{
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown compound designation",_INPUT_ILLEGAL_);
    }

    if(clean){
      for(uint i=0;i<velements.size();i++){aurostd::CleanStringASCII_InPlace(velements[i]);}  //CO20190712 - extra cleaning from XATOM_SplitAlloySpecies
    }

    // Add implicit ones
    for (uint i = vcomposition.size(); i < velements.size(); i++) vcomposition.push_back((utype)1.0);

    //////////////////////////////////////////////////////////////////////////////
    // END Parsing input
    //////////////////////////////////////////////////////////////////////////////

    if(velements.size()==0){pflow::logger(_AFLOW_FILE_NAME_, soliloquy, "No elements found", FileMESSAGE, oss, _LOGGER_ERROR_);}

    if(sort_elements && velements.size()>1){
      string etmp="";
      utype ctmp=(utype)0.0;
      for(uint i=0;i<velements.size()-1;i++){
        for(uint j=i+1;j<velements.size();j++){
          if(velements[i]>velements[j]){
            etmp=velements[j];  //fix old j
            velements[j]=velements[i];  //swap
            velements[i]=etmp;  //set i to old j
            ctmp=vcomposition[j]; //fix old j
            vcomposition[j]=vcomposition[i];  //swap
            vcomposition[i]=ctmp; //set i to old j
          }
        }
      }
    }

    return velements;
  }
} // namespace aurostd

#endif // _AUROSTD_XPARSER_CPP_

// **************************************************************************
// *                                                                        *
// *              Aflow COREY OSES - Duke University 2003-2020              *
// *                                                                        *
// **************************************************************************
