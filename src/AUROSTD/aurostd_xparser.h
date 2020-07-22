// ***************************************************************************
// *                                                                         *
// *              Aflow COREY OSES - Duke University 2003-2020               *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses 2020

#ifndef _AUROSTD_XPARSER_H_
#define _AUROSTD_XPARSER_H_

#include "aurostd.h"

//compound specification is how a compound is specified
//composition (Mn2Pt3) is ORTHOGONAL to pseudopotential string (Mn_pvPt)
//for instance, H1.25 can be a pseudopotential and NOT a composition
enum compound_designation {
  composition_string,
  pp_string,
};

//CO20190712 - see VASP_PseudoPotential_CleanName_InPlace() in aflow_ivasp.cpp
const string CAPITAL_LETTERS_PP_LIST="_GW2"    //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Li_AE_GW2
",_GW"    //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/As_GW
",_ZORA"  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Pt_ZORA
",_LDApU" //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Zn_sv_LDApU
",_AE"    //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Li_AE_GW2
",_NC2"   //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/As_NC2
",_200eV"
"";

namespace aurostd {
  void VASP_PseudoPotential_CleanName_InPlace(string& species,bool capital_letters_only=false); //CO20190712
  ////////////////////////////////////////////////////////////////////////////////
  vector<string> getElements(const string& _input); //CO20190712
  void elementsFromCompositionString(const string& input);  //CO20190712
  void elementsFromCompositionString(const string& input,vector<string>& velements,vector<double>& vcomposition); //CO20190712
  void elementsFromPPString(const string& input,vector<string>& velements,bool keep_pp=false); //CO20190712
  ////////////////////////////////////////////////////////////////////////////////
  // returns UNSORTED vector<string> from string
  vector<string> stringElements2VectorElements(const string& input, bool clean=true, bool sort_elements=false, compound_designation c_desig=composition_string, bool keep_pp=false, ostream& oss=cout);
  vector<string> stringElements2VectorElements(const string& input, vector<double>& vcomposition, bool clean=true, bool sort_elements=false, compound_designation c_desig=composition_string, bool keep_pp=false, ostream& oss=cout);  //ME20190628
  vector<string> stringElements2VectorElements(const string& input, ofstream& FileMESSAGE, bool clean=true, bool sort_elements=false, compound_designation c_desig=composition_string, bool keep_pp=false, ostream& oss=cout);
  vector<string> stringElements2VectorElements(const string& input, vector<double>& vcomposition, ofstream& FileMESSAGE, bool clean=true, bool sort_elements=false, compound_designation c_desig=composition_string, bool keep_pp=false, ostream& oss=cout);  //ME20190628
} // namespace aurostd

#endif // _AUROSTD_XPARSER_H_

// **************************************************************************
// *                                                                        *
// *              Aflow COREY OSES - Duke University 2003-2020              *
// *                                                                        *
// **************************************************************************
