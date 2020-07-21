// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo

#ifndef _AFLOWLIB_LIBRARIES_CPP_
#define _AFLOWLIB_LIBRARIES_CPP_

#include "aflow.h"
#include "aflowlib.h"
#include "aflow_pflow.h"
#include "aflow_bader.h"
#include "aflow_cce.h" //CO20200624
#include "aflow_pocc.h" //CO20200624
#include "aflow_matlab_funcs.cpp" //CO20200508
#include "aflow_gnuplot_funcs.cpp" //CO20200508

using std::vector;
using std::deque;
using std::string;

vector<string> vLibrary_ALL;vector<vector<string> > vLibrary_ALL_tokens;
vector<string> vLibrary_ICSD;vector<vector<string> > vLibrary_ICSD_tokens;
vector<string> vLibrary_LIB0;vector<vector<string> > vLibrary_LIB0_tokens;
vector<string> vLibrary_LIB1;vector<vector<string> > vLibrary_LIB1_tokens;
vector<string> vLibrary_LIB2;vector<vector<string> > vLibrary_LIB2_tokens;
vector<string> vLibrary_LIB3;vector<vector<string> > vLibrary_LIB3_tokens;
vector<string> vLibrary_LIB4;vector<vector<string> > vLibrary_LIB4_tokens;
vector<string> vLibrary_LIB5;vector<vector<string> > vLibrary_LIB5_tokens;
vector<string> vLibrary_LIB6;vector<vector<string> > vLibrary_LIB6_tokens;
vector<string> vLibrary_LIB7;vector<vector<string> > vLibrary_LIB7_tokens;
vector<string> vLibrary_LIB8;vector<vector<string> > vLibrary_LIB8_tokens;
vector<string> vLibrary_LIB9;vector<vector<string> > vLibrary_LIB9_tokens;

//[moved to aflow.h - CO20180705, list also looks a bit obsolete from avasp.cpp]#define SPECIE_RAW_LIB3 string("Ag,Al,As,Au,B_h,Bi_d,Cd,Co,Cr_pv,Cu_pv,Fe_pv,Ga_h,Ge_h,Hf_pv,Hg,In_d,Ir,La,Mg_pv,Mn_pv,Mo_pv,Nb_sv,Ni_pv,Os_pv,P,Pb_d,Pd_pv,Pt,Re_pv,Rh_pv,Ru_pv,Sb,Sc_sv,Se,Si,Sn,Ta_pv,Te,Tc_pv,Ti_sv,V_sv,W_pv,Y_sv,Zn,Zr_sv")

bool AFLOWLIB_VERBOSE=TRUE; // FALSE;
#define _EPSILON_COMPOSITION_ 0.001

#define USE_AFLOW_SG
//#define USE_PLATON_SG

#define AFLOW_MAX_ARGV 1024
#define AFLOW_CORE_TEMPERATURE_LIB2RAW 46.0

// ******************************************************************************************************************************************************
using aurostd::FileExist;

// ***************************************************************************
// aflowlib::TokenPresentAFLOWLIB
// ***************************************************************************
namespace aflowlib {
  bool TokenPresentAFLOWLIB(string line,string query) { return aurostd::substring2bool(line,query,TRUE); }
}

// ***************************************************************************
// aflowlib::TokenExtractAFLOWLIB
// ***************************************************************************
namespace aflowlib {
  string TokenExtractAFLOWLIB(string line,string query) {
    if(!TokenPresentAFLOWLIB(line,query)) return "";
    vector<string> tokens,tk;aurostd::string2tokens(line,tokens,"|");
    string qquery=query;aurostd::StringSubst(qquery,"=","");aurostd::StringSubst(qquery," ","");
    for(uint i=0;i<tokens.size();i++) {
      aurostd::StringSubst(tokens.at(i)," ","");
      aurostd::string2tokens(tokens.at(i),tk,"=");
      if(tk.size()==2) { //   cerr << XPID << tk.at(0) << endl;
        if(tk.at(0)==qquery) { return tk.at(1); }
        // return aurostd::substring2string(tokens.at(i),query,TRUE);
      }
    }
    return "";
  }
}

namespace aflowlib {
  bool TokenExtractAFLOWLIB(string line,string query,string &value) {
    // cerr << XPID << "HERE=[TokenExtractAFLOWLIB(string line,string query,string &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    value=TokenExtractAFLOWLIB(line,query);return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,int &value) {
    // cerr << XPID << "HERE=[TokenExtractAFLOWLIB(string line,string query,int &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    value=aurostd::string2utype<int>(TokenExtractAFLOWLIB(line,query));return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,uint &value) {
    // cerr << XPID << "HERE=[TokenExtractAFLOWLIB(string line,string query,uint &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    value=aurostd::string2utype<uint>(TokenExtractAFLOWLIB(line,query));return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,float &value) {
    // cerr << XPID << "HERE=[TokenExtractAFLOWLIB(string line,string query,float &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    value=aurostd::string2utype<float>(TokenExtractAFLOWLIB(line,query));return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,double &value) {
    // cerr << XPID << "HERE=[TokenExtractAFLOWLIB(string line,string query,double &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    value=aurostd::string2utype<double>(TokenExtractAFLOWLIB(line,query));return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,vector<string> &value) {
    // cerr << XPID << "HERE=[TokenExtractAFLOWLIB(string line,string query,vector<string> &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    aurostd::string2tokens(TokenExtractAFLOWLIB(line,query),value,",");return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,vector<int> &value) {
    // cerr << XPID << "HERE=[TokenExtractAFLOWLIB(string line,string query,vector<int> &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    vector<string> vvalue;aurostd::string2tokens(TokenExtractAFLOWLIB(line,query),vvalue,",");value=aurostd::vectorstring2vectorint(vvalue);return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,vector<uint> &value) {
    // cerr << XPID << "HERE=[TokenExtractAFLOWLIB(string line,string query,vector<uint> &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    vector<string> vvalue;aurostd::string2tokens(TokenExtractAFLOWLIB(line,query),vvalue,",");value=aurostd::vectorstring2vectoruint(vvalue);return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,vector<float> &value) {
    // cerr << XPID << "HERE=[TokenExtractAFLOWLIB(string line,string query,vector<float> &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    vector<string> vvalue;aurostd::string2tokens(TokenExtractAFLOWLIB(line,query),vvalue,",");value=aurostd::vectorstring2vectorfloat(vvalue);return TokenPresentAFLOWLIB(line,query); }
  bool TokenExtractAFLOWLIB(string line,string query,vector<double> &value) {
    // cerr << XPID << "HERE=[TokenExtractAFLOWLIB(string line,string query,vector<double> &value)]" << TokenExtractAFLOWLIB(line,query) << "<br>" << endl;
    vector<string> vvalue;aurostd::string2tokens(TokenExtractAFLOWLIB(line,query),vvalue,",");value=aurostd::vectorstring2vectordouble(vvalue);return TokenPresentAFLOWLIB(line,query); }
}

// ******************************************************************************************************************************************************
// aflowlib::GREP_Species_ALL
// ***************************************************************************
namespace aflowlib {
  uint GREP_Species_ALL(vector<string> vspecies,                         // IN  [0,nspecies[ the species Ag,Cd,...   nspecies=number of these items    nspecies=naries
      vector<string>& vspecies_pp,                     // IN  [0,nspecies[ the pseudopotentials Ag_pv, Cd_sv
      vector<vector<string> >& vList,                  // OUT [0,naries[*[0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
      vector<vector<vector<string> > > &vList_tokens,  // OUT [0,naries[*[0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
      vector<vector<vector<string> > > &vList_species, // OUT [0,naries[*[0,vList.size()[*nspecies  returns the species present
      vector<double> &vList_Hmin,                      // OUT [0,nspecies[ returns the min enthalpy for reference
      vector<string> &vList_Pmin,                      // OUT [0,nspecies[ returns the prototype for reference
      vector<uint> &vList_Imin,                        // OUT [0,nspecies[ returns the line index of vList.at(ispecies), in which we have the min enthalpy for reference
      vector<vector<vector<double> > > &vList_concs,   // OUT [0,naries[*[0,vList.size()[*[0.nspecies[ the concentrations AxAyCz... where x+y+z=1 and it contains also ZEROS so that 0 0.25 0.75 is allowed
      vector<vector<double> > &vList_Hf) {             // OUT [0,naries[*[0,vList.size()[ returns the formation enthalpy of the list
    bool LVERBOSE=FALSE;
    bool LIBVERBOSE=TRUE;
    if(LVERBOSE) cerr << XPID << "GREP_Species_ALL: START" << endl;
    if(vspecies.size()==1) LOAD_Library_ALL("no_aflowlib_lib3,no_aflowlib_lib4,no_aflowlib_lib5,no_aflowlib_lib6,no_aflowlib_lib7,no_aflowlib_lib8,no_aflowlib_lib9,no_aflowlib_lib2",LIBVERBOSE);
    else LOAD_Library_ALL(LIBVERBOSE); // might not have been loaded, yet
    // clean stuff

    // start
    vector<string> tokens,tokens_species,tokens_species_pp;
    vector<uint> tokens_composition;
    uint natoms,nspecies=vspecies.size(),naries=vspecies.size();
    string species;
    if(nspecies==0) return 0;
    if(vspecies_pp.size()==0 || vspecies.size()!=vspecies_pp.size())  { // generate some default ones
      vspecies_pp.clear(); for(uint i=0;i<nspecies;i++)vspecies_pp.push_back(AVASP_Get_PseudoPotential_PAW_PBE(vspecies.at(i)));
    }
    for(uint i=0;i<nspecies;i++)
      if(vspecies_pp.at(i).empty()) vspecies_pp.at(i)=AVASP_Get_PseudoPotential_PAW_PBE(vspecies.at(i)); // in the event that I pass an empty string

    // generate space
    vector<string> dummy_vstring;
    vList.clear(); for(uint i=0;i<naries;i++) vList.push_back(dummy_vstring);                  // create space  [0,naries[*[0,vList.size()[  <strings>
    vector<vector<string> > dummy_vvstring;
    vList_tokens.clear(); for(uint i=0;i<naries;i++) vList_tokens.push_back(dummy_vvstring);   // create space  [0,naries[*[0,vList.size()[*[0,size_tokens[  <strings>
    vList_Hmin.clear(); for(uint i=0;i<nspecies;i++) vList_Hmin.push_back(1E6);                // create space  [0,nspecies[ <double>
    vList_Pmin.clear(); for(uint i=0;i<nspecies;i++) vList_Pmin.push_back("");                 // create space  [0,nspecies[ <string>
    vList_Imin.clear(); for(uint i=0;i<nspecies;i++) vList_Imin.push_back(0);                  // create space  [0,nspecies[ <uint>
    vector<double> dummy_vdouble;
    vList_Hf.clear(); for(uint i=0;i<naries;i++) vList_Hf.push_back(dummy_vdouble);            // create space  [0,naries[*[0,vList.size()[  <double>
    vector<vector<double> > dummy_vvdouble;
    vList_concs.clear(); for(uint i=0;i<naries;i++) vList_concs.push_back(dummy_vvdouble);     // create space  [0,naries[*[0,vList.size()[*[0.nspecies[   <double>
    vList_species.clear(); for(uint i=0;i<naries;i++) vList_species.push_back(dummy_vvstring); // create space  [0,naries[*[0,vList.size()[*[0.nspecies[   <string>

    uint num_species=0;string species_pp;double enthalpy_atom,enthalpy;
    // for(uint ispecies=1;ispecies<=nspecies;ispecies++)
    for(uint iall=0;iall<vLibrary_ALL.size();iall++)
    { //CO20200106 - patching for auto-indenting
      bool found_some_species=FALSE;
      for(uint ispecies=0;((ispecies<nspecies) && (!found_some_species));ispecies++)
        if(aurostd::substring2bool(vLibrary_ALL.at(iall),vspecies_pp.at(ispecies),TRUE)) found_some_species=TRUE;
      if(found_some_species) {
        TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"nspecies=",num_species); // search for num_species
        for(uint ispecies=1;ispecies<=nspecies;ispecies++) {
          bool found_nspecies=FALSE,found_species_pp=FALSE;
          if(num_species==ispecies) found_nspecies=TRUE;
          if(found_nspecies && num_species==1) {
            if(!TokenPresentAFLOWLIB(vLibrary_ALL.at(iall),"LIB0")) found_nspecies=FALSE; // patch to use only pure
            if(!TokenPresentAFLOWLIB(vLibrary_ALL.at(iall),"LIB1")) found_nspecies=FALSE; // patch to use only pure
            if(found_nspecies) {
              if(TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=")=="Rh_pv" && TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=")=="A6") found_nspecies=FALSE;  // patches for incorrect relaxation
              if(TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=")=="Ni_pv" && TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=")=="A6") found_nspecies=FALSE;  // patches for incorrect relaxation
              // if(TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=")=="La" && TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=")=="A1") found_nspecies=FALSE;  // patches for incorrect relaxation
              if(TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=")=="Au" && TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=")=="A7") found_nspecies=FALSE;  // patches for incorrect relaxation
              if(TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=")=="Au" && TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=")=="A8") found_nspecies=FALSE;  // patches for incorrect relaxation
              // fix La Mn and Hg
            }
          }
          // search for species_pp
          if(found_nspecies) {
            TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=",species_pp);
            TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"species_pp=",tokens_species);
            uint species_pp_matches=0;
            for(uint k1=0;k1<vspecies_pp.size();k1++)
              for(uint k2=0;k2<tokens_species.size();k2++)
                if(vspecies_pp.at(k1)==tokens_species.at(k2)) species_pp_matches++;
            if(species_pp_matches==ispecies) found_species_pp=TRUE;
          }
          if(found_species_pp) {
            // cerr << XPID << ispecies << " " << vLibrary_ALL.at(iall) << endl;
            string prototype="";
            TokenExtractAFLOWLIB(vLibrary_ALL.at(iall),"prototype=",prototype);
            if(prototype=="64" || prototype=="65" || prototype=="549" || prototype=="550") {  // check for broken calcs
              if(LVERBOSE) cerr << XPID << "GREP_Species_ALL: Removing prototype=" << prototype << endl;
            } else {
              vList.at(ispecies-1).push_back(vLibrary_ALL.at(iall));       // add lines
              aurostd::string2tokens(vLibrary_ALL.at(iall),tokens,"|");    // add tokens
              vList_tokens.at(ispecies-1).push_back(tokens);               // add tokens
              if(LVERBOSE) if(!mod((int) vList.at(ispecies-5).at(vList.at(ispecies-1).size()-1).size(),1)) cerr << XPID << vLibrary_ALL_tokens.at(iall).at(0) << " " << species_pp << endl;  // for DEBUG printout
            }
          }
        }
      }
    }
    if(LVERBOSE) for(uint ispecies=0;ispecies<nspecies;ispecies++) cerr << XPID << ispecies << " " << vList.at(ispecies).size() << " " << vList_tokens.at(ispecies).size() << endl;
    // LIST PREPARED
    if(LVERBOSE) cerr << XPID << "GREP_Species_ALL: LIST PREPARED" << endl;
    // creating reference enthalpies
    if(LVERBOSE) cerr << XPID << "GREP_Species_ALL: creating reference enthalpies" << endl;
    for(uint ispecies=0;ispecies<nspecies;ispecies++) {
      if(LVERBOSE) cerr << XPID << "GREP_Species_ALL: getting [" << vspecies.at(ispecies) << "]" << endl;
      for(uint j=0;j<vList.at(0).size();j++) { // only pure
        if(LVERBOSE) cerr << XPID << j << " " << vList_tokens.at(0).at(j).at(0) << endl;
        TokenExtractAFLOWLIB(vList.at(0).at(j),"species=",species);
        TokenExtractAFLOWLIB(vList.at(0).at(j),"enthalpy_atom=",enthalpy_atom);
        if(species==vspecies.at(ispecies)) {
          if(LVERBOSE) cerr << XPID << j << " " << vList_tokens.at(0).at(j).at(0) << endl;
          if(enthalpy_atom<vList_Hmin.at(ispecies)) {
            if(LVERBOSE) cerr << XPID << "GREP_Species_ALL: ***** FOUND " << j << " " << vList_tokens.at(0).at(j).at(0) << endl;
            vList_Hmin.at(ispecies)=enthalpy_atom;
            vList_Pmin.at(ispecies)=TokenExtractAFLOWLIB(vList.at(0).at(j),"prototype=");
            vList_Imin.at(ispecies)=j;
          }
        }
      }
    }

    //JX THIS IS THE WAY TO READ vList_Hmin.at(ispecies) vList_Pmin.at(ispecies) vList_Imin.at(ispecies) right use of vList_tokens.at(0,1,2 A,AB,ABC).at(vList_Imin.at(ispecies)).at(TOKENS INDEX THAT YOU WANT)
    // HOWEVER YOU CAN EXTRACT INFORMATION WITH TokenExtractAFLOWLIB(vList..at(j),"species=",species);
    for(uint ispecies=0;ispecies<nspecies;ispecies++) {
      if(LVERBOSE) cerr << XPID << vspecies.at(ispecies) << " " << vList_Hmin.at(ispecies) << " " << vList_Pmin.at(ispecies) << " " << vList_Imin.at(ispecies) << " " << vList_tokens.at(0).at(vList_Imin.at(ispecies)).at(0) << endl;
    }
    for(uint iary=0;iary<naries;iary++) {
      if(LVERBOSE) cerr << XPID << vList_tokens.at(iary).size() << endl;
    }

    // creating formation enthalpies
    for(uint iaries=0;iaries<naries;iaries++) {
      for(uint j=0;j<vList.at(iaries).size();j++) { // only pure
        vector<string> vspecies_local;
        vector<double> vconcs_local;
        TokenExtractAFLOWLIB(vList.at(iaries).at(j),"natoms=",natoms); // extract natoms
        TokenExtractAFLOWLIB(vList.at(iaries).at(j),"enthalpy_cell=",enthalpy); // extract enthalpy
        TokenExtractAFLOWLIB(vList.at(iaries).at(j),"species=",tokens_species);
        TokenExtractAFLOWLIB(vList.at(iaries).at(j),"species_pp=",tokens_species_pp);
        TokenExtractAFLOWLIB(vList.at(iaries).at(j),"composition=",tokens_composition);

        if(tokens_species.size()!=tokens_composition.size()) { cerr << XPID << "GREP_Species_ALL: tokens_species.size()!=tokens_composition.size()" << endl;exit(0); }
        // cerr << XPID << "GREP_Species_ALL: E=" << enthalpy << " ";
        for(uint k=0;k<nspecies;k++) // order by species, no entries
          for(uint j=0;j<tokens_species.size();j++)
            // if(tokens_species.at(j)==vspecies.at(k)) // found specie k associated with position k
            if(tokens_species_pp.at(j)==vspecies_pp.at(k)) // found specie_pp k associated with position k
            { //CO20200106 - patching for auto-indenting
              enthalpy=enthalpy-(double) vList_Hmin.at(k)*tokens_composition.at(j); // remove formation
              vspecies_local.push_back(vspecies_pp.at(k));
              vconcs_local.push_back((double) tokens_composition.at(j)/(double) natoms);
            }
        enthalpy=enthalpy/(double) natoms; // normalize per atom
        if(aurostd::abs(enthalpy)<0.0001) enthalpy=0.0; // <0.1meV equal zero
        // cerr << enthalpy << endl;

        vList_Hf.at(iaries).push_back(enthalpy);              // plug scalar enthalpy_formation_atom
        vList_species.at(iaries).push_back(vspecies_local); // plug VECTOR species
        vList_concs.at(iaries).push_back(vconcs_local);     // plug VECTOR concentrations
      }
    }

    for(uint iaries=0;iaries<naries;iaries++) {
      for(uint jentry=0;jentry<vList.at(iaries).size();jentry++) { // only pure
        if(LVERBOSE) cerr << XPID << vList_tokens.at(iaries).at(jentry).at(0) << "  " << vList_tokens.at(iaries).at(jentry).at(1) << "   Hf = " << vList_Hf.at(iaries).at(jentry) << "  Compound = ";
        for(uint k=0;k<vList_species.at(iaries).at(jentry).size();k++) if(LVERBOSE) cerr << XPID << vList_species.at(iaries).at(jentry).at(k) << "_" << vList_concs.at(iaries).at(jentry).at(k) << " ";
        if(LVERBOSE) cerr << endl;
      }
    }

    // RETURN
    uint numentries=0;
    for(uint iaries=0;iaries<naries;iaries++) numentries+=vList.at(iaries).size();
    return numentries;
  }
}

// ***************************************************************************
// aflowlib::GREP_Species_ALL
// ***************************************************************************
namespace aflowlib {
  uint GREP_Species_ALL(string species,                                   // IN   the species Ag,...   nspecies=number of these items
      string& species_pp,                               // IN   the pseudopotentials Ag_pv,
      vector<string> & vList,                           // OUT  [0,nspecies[*[0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
      vector<vector<string> > &vList_tokens,            // OUT [0,nspecies[*[0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
      double& List_Hmin,                                // OUT  returns the min enthalpy for reference
      string& List_Pmin,                                // OUT  returns the prototype for reference
      uint& List_Imin) {                                // OUT  returns the line index of vList, in which we have the min enthalpy for reference
    vector<string> vspecies;vspecies.push_back(species);
    vector<string> vspecies_pp;vspecies_pp.push_back(species_pp);

    vector<vector<string> > _vList;
    vector<vector<vector<string> > > _vList_tokens;
    vector<vector<vector<string> > > vList_species;
    vector<double> vList_Hmin;
    vector<string> vList_Pmin;
    vector<uint> vList_Imin;
    vector<vector<vector<double> > > vList_concs;
    vector<vector<double> > vList_Hf;

    uint out=GREP_Species_ALL(vspecies,vspecies_pp,_vList,_vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imin,vList_concs,vList_Hf);
    species_pp=vspecies_pp.at(0);
    List_Hmin=vList_Hmin.at(0);  // they are for ONE SPECIE =
    List_Pmin=vList_Pmin.at(0);  // they are for ONE SPECIE
    List_Imin=vList_Imin.at(0);  // they are for ONE SPECIE
    for(uint i=0;i<_vList.at(0).size();i++) vList.push_back(_vList.at(0).at(i));
    for(uint i=0;i<_vList_tokens.at(0).size();i++) vList_tokens.push_back(_vList_tokens.at(0).at(i));

    return out;
  }
}

// ***************************************************************************
// aflowlib::GREP_Species_ALL
// ***************************************************************************
namespace aflowlib {
  uint GREP_Species_ALL(vector<string> vspecies,vector<string>& vspecies_pp,vector<vector<string> >& vList,vector<vector<vector<string> > > &vList_tokens) {
    vector<vector<vector<string> > > vList_species;
    vector<vector<vector<double> > > vList_concs;
    vector<vector<double> > vList_Hf;
    vector<double> vList_Hmin;
    vector<string> vList_Pmin;
    vector<uint> vList_Imin;

    return GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imin,vList_concs,vList_Hf);
  }
}

// ***************************************************************************
// aflowlib::GREP_Species_ALL
// ***************************************************************************
namespace aflowlib {
  uint GREP_Species_ALL(vector<string> vspecies,vector<string>& vspecies_pp,vector<vector<string> >& vList) {
    bool LVERBOSE=FALSE;
    bool LIBVERBOSE=TRUE;
    if(LVERBOSE) cerr << XPID << "GREP_Species_ALL: START" << endl;
    if(vspecies.size()==1) LOAD_Library_ALL("no_aflowlib_lib3,no_aflowlib_lib4,no_aflowlib_lib5,no_aflowlib_lib6,no_aflowlib_lib7,no_aflowlib_lib8,no_aflowlib_lib9,no_aflowlib_lib2",LIBVERBOSE);
    else LOAD_Library_ALL(LIBVERBOSE); // might not have been loaded, yet
    // clean stuff

    // start //
    vector<string> tokens,tokens_species,tokens_species_pp;
    uint nspecies=vspecies.size(),naries=vspecies.size();
    string species;
    if(nspecies==0) return 0;
    if(vspecies_pp.size()==0 || vspecies.size()!=vspecies_pp.size())  { // generate some default ones
      vspecies_pp.clear();
      for(uint i=0;i<nspecies;i++) {
        vspecies_pp.push_back(AVASP_Get_PseudoPotential_PAW_PBE(vspecies.at(i)));
      }
    }
    for(uint i=0;i<nspecies;i++) {
      if(vspecies_pp.at(i).empty()) {
        vspecies_pp.at(i)=AVASP_Get_PseudoPotential_PAW_PBE(vspecies.at(i)); // in the event that I pass an empty string
      }
    }
    // generate space
    vector<string> dummy_vstring;
    vList.clear(); for(uint i=0;i<naries;i++) vList.push_back(dummy_vstring);                  // create space  [0,naries[*[0,vList.size()[  <strings>

    cerr << XPID << "START GREP" << endl;

    uint num_species=0;string species_pp;
    // for(uint ispecies=1;ispecies<=nspecies;ispecies++) { //[CO20200106 - close bracket for indenting]}
    string vLibrary_ALL_at_iall;
    for(uint iall=0;iall<vLibrary_ALL.size();iall++) {
      vLibrary_ALL_at_iall=vLibrary_ALL.at(iall);
      bool found_some_species=FALSE;
      for(uint ispecies=0;((ispecies<nspecies) && (!found_some_species));ispecies++)
        if(aurostd::substring2bool(vLibrary_ALL_at_iall,vspecies_pp.at(ispecies),TRUE)) found_some_species=TRUE;
      if(found_some_species) {
        for(uint ispecies=1;ispecies<=nspecies;ispecies++) {
          if(aurostd::substring2bool(vLibrary_ALL_at_iall,string("nspecies="+aurostd::utype2string(ispecies)))) {
            TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"nspecies=",num_species); // search for num_species
            // cerr << XPID << " [nspecies=" << aurostd::utype2string(ispecies) << "] " << endl;
            bool found_nspecies=FALSE,found_species_pp=FALSE;
            if(num_species==ispecies) found_nspecies=TRUE;
            if(found_nspecies && num_species==1) {
              if(!TokenPresentAFLOWLIB(vLibrary_ALL_at_iall,"LIB0")) found_nspecies=FALSE; // patch to use only pure
              if(!TokenPresentAFLOWLIB(vLibrary_ALL_at_iall,"LIB1")) found_nspecies=FALSE; // patch to use only pure
              if(found_nspecies) {
                if(TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=")=="Rh_pv" && TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=")=="A6") found_nspecies=FALSE;  // patches for incorrect relaxation
                if(TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=")=="Ni_pv" && TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=")=="A6") found_nspecies=FALSE;  // patches for incorrect relaxation
                // if(TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=")=="La" && TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=")=="A1") found_nspecies=FALSE;  // patches for incorrect relaxation
                if(TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=")=="Au" && TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=")=="A7") found_nspecies=FALSE;  // patches for incorrect relaxation
                if(TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=")=="Au" && TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=")=="A8") found_nspecies=FALSE;  // patches for incorrect relaxation
                // fix La Mn and Hg
              }
            }
            // search for species_pp
            if(found_nspecies) {
              TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=",species_pp);
              TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"species_pp=",tokens_species);
              uint species_pp_matches=0;
              for(uint k1=0;k1<vspecies_pp.size();k1++)
                for(uint k2=0;k2<tokens_species.size();k2++)
                  if(vspecies_pp.at(k1)==tokens_species.at(k2)) species_pp_matches++;
              if(species_pp_matches==ispecies) found_species_pp=TRUE;
            }
            if(found_species_pp) {
              // cerr << XPID << ispecies << " " << vLibrary_ALL_at_iall << endl;
              string prototype="";
              TokenExtractAFLOWLIB(vLibrary_ALL_at_iall,"prototype=",prototype);
              if(prototype=="64" || prototype=="65" || prototype=="549" || prototype=="550") {  // check for broken calcs
                if(LVERBOSE) cerr << XPID << "GREP_Species_ALL: Removing prototype=" << prototype << endl;
              } else {
                vList.at(ispecies-1).push_back(vLibrary_ALL_at_iall);       // add lines
                aurostd::string2tokens(vLibrary_ALL_at_iall,tokens,"|");    // add tokens
              }
            }
          }
        }
      }
    }

    cerr << XPID << "END GREP" << endl;
    // RETURN
    uint numentries=0;
    for(uint iaries=0;iaries<naries;iaries++) numentries+=vList.at(iaries).size();
    return numentries;
  }
}

// ***************************************************************************
// aflowlib::GREP_Species_ALL
// ***************************************************************************
namespace aflowlib {
  uint GREP_Species_ALL(string species,                                   // IN   the species Ag,...   nspecies=number of these items
      string& species_pp,                               // IN   the pseudopotentials Ag_pv,
      double& List_Hmin,                                // OUT  returns the min enthalpy for reference
      string& List_Pmin) {                              // OUT  returns the prototype for reference
    vector<string> vList;
    vector<vector<string> > vList_tokens;
    uint List_Imin;
    return GREP_Species_ALL(species,species_pp,vList,vList_tokens,List_Hmin,List_Pmin,List_Imin);
  }
}

// ***************************************************************************
// aflowlib::LIBS_EFormation
// ***************************************************************************
namespace aflowlib {
  bool LIBS_EFormation(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::ALIBRARIES: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    bool LVERBOSE=FALSE;
    aflowlib::LOAD_Library_LIBRARY("pure","",LVERBOSE);

    string species,species_pp;
    // chose somehow the species
    if(tokens.size()>=1) { species=tokens.at(0); } else { species="Nb"; }

    //  cerr << XPID << "vLibrary_ALL.size()=" << vLibrary_ALL.size() << " vLibrary_ALL_tokens.size()=" << vLibrary_ALL_tokens.size() << endl;

    vector<string> vspecies,vspecies_pp;
    aurostd::string2tokens(SPECIE_RAW_LIB3,vspecies,",");
    for(uint i=0;i<vspecies.size();i++) {
      vspecies.at(i)=KBIN::VASP_PseudoPotential_CleanName(vspecies.at(i));
      vspecies_pp.push_back("");//AVASP_Get_PseudoPotential_PAW_PBE(vspecies.at(i)));

      vector<string> vList;
      vector<vector<string> > vList_tokens;
      double List_Hmin=0.0;
      string List_Pmin="";
      uint List_Imin=0;

      species_pp=AVASP_Get_PseudoPotential_PAW_PBE(species);

      GREP_Species_ALL(vspecies.at(i),          // IN    the species Ag,...   nspecies=number of these items
          vspecies_pp.at(i),       // IN    the pseudopotentials Ag_pv,
          vList,            // OUT   [0,vList.size()[ returns the lines of the library containing A,B,C,AB,AC,BC,ABC....
          vList_tokens,     // OUT   [0,vList.size()[*[0,size_tokens[ returns the tokens for each line of the vList
          List_Hmin,        // OUT   returns the min enthalpy for reference
          List_Pmin,        // OUT   returns the prototype for reference
          List_Imin);
      cerr << XPID << vspecies.at(i) << " " << vspecies_pp.at(i) << " " << List_Hmin << " " << List_Pmin << " " << List_Imin << " " << vList_tokens.at(List_Imin).at(0) << endl;
    }
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::ALIBRARIES
// ***************************************************************************
namespace aflowlib {
  bool ALIBRARIES(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::ALIBRARIES: BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");

    if(tokens.size()>10) {
      init::ErrorOption(options,"aflowlib::ALIBRARIES","aflow --qhull[=Nb[,Pt[,Rh[,...]]]]");
      exit(0);
    }

    bool LVERBOSE=TRUE;

    LOAD_Library_ALL(LVERBOSE);
    //  LOAD_Library_ALL();
    cerr << XPID << "vLibrary_ALL.size()=" << vLibrary_ALL.size() << " vLibrary_ALL_tokens.size()=" << vLibrary_ALL_tokens.size() << endl;

    vector<string> vspecies,vspecies_pp;
    vector<vector<string> > vList;
    vector<vector<vector<string> > > vList_tokens;
    vector<double> vList_Hmin;
    vector<string> vList_Pmin;
    vector<uint> vList_Imin;
    vector<vector<vector<string> > > vList_species;
    vector<vector<vector<double> > > vList_concs;
    vector<vector<double> > vList_Hf;

    // uint nentries=GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imin,vList_concs,vList_Hf);  // complete
    // uint nentries=GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imi);

    // [OBSOLETE] double List_Hmin=1E6;
    string List_Pmin="";
    // [OBSOLETE] uint List_Imin=0;

    // chose somehow the species
    if(tokens.size()==0) {
      vspecies.push_back("Nb");vspecies.push_back("Pt");vspecies.push_back("Rh");  // DO IN ALPHABETIC ORDER OTHERWISE COMPLAIN
    } else {
      for(uint i=0;i<tokens.size();i++)
        vspecies.push_back(tokens.at(i));
    }

    //  uint nentries=GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens);
    uint nentries=GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens,vList_species,vList_Hmin,vList_Pmin,vList_Imin,vList_concs,vList_Hf);
    // uint nentries=GREP_Species_ALL(vspecies,vspecies_pp,vList,vList_tokens);
    //  uint nentries=GREP_Species_ALL(vspecies.at(0),"",vList,vList_tokens,List_Hmin,List_Pmin,List_Imin);
    for(uint iaries=0;iaries<vList.size();iaries++) {
      //  cerr << XPID << iaries << " " << vList.at(iaries).size() << endl;
      //    for(uint j=0;j<vList.at(iaries).size();j++) cerr << vList.at(iaries).at(j) << endl;
    }

    for(uint iaries=0;iaries<vList.size();iaries++) {
      if(LVERBOSE) cerr << XPID << "***************[" << iaries+1 << "]*************" << endl;
      for(uint jentry=0;jentry<vList.at(iaries).size();jentry++) { // only pure
        if(LVERBOSE) cerr << XPID << vList_tokens.at(iaries).at(jentry).at(0) << "  " << vList_tokens.at(iaries).at(jentry).at(1) << "   Hf = " << vList_Hf.at(iaries).at(jentry) << "  Compound = ";
        for(uint k=0;k<vList_species.at(iaries).at(jentry).size();k++) if(LVERBOSE) cerr << vList_species.at(iaries).at(jentry).at(k) << "_" << vList_concs.at(iaries).at(jentry).at(k) << " ";
        if(LVERBOSE) cerr << endl;
      }
    }

    if(LVERBOSE) cerr << XPID << nentries << endl;
    cerr << XPID << "vList.at(iaries).size()="; for(uint iaries=0;iaries<vList.size();iaries++) cerr << vList.at(iaries).size() << " "; cerr << endl;
    cerr << XPID << "vList_tokens.at(iaries).size()="; for(uint iaries=0;iaries<vList_tokens.size();iaries++) cerr << vList_tokens.at(iaries).size() << " "; cerr << endl;
    // cerr << XPID << "List_Hmin=" << List_Hmin << endl;
    // cerr << XPID << "List_Pmin=" << List_Pmin << endl;
    // cerr << XPID << "List_Imin=" << List_Imin << endl;

    // exit(0);
    if(LDEBUG) cerr << XPID << "aflowlib::ALIBRARIES: END" << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LOAD_Library_LIBRARY
// ***************************************************************************
namespace aflowlib {
  uint LOAD_Library_LIBRARY(string file,string grep,bool LVERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // LDEBUG=TRUE;
    string FileLibrary="";
    vector<string> tokens;
    bool found=FALSE;
    vector<string>* pvLibrary=NULL; vector<vector<string> >* pvLibrary_tokens=NULL;
    // ----------------------------------------------------------------------
    // ----------------------------------------------------------------------
    if(!found && (file=="aflowlib_lib0" ||
          file=="aflowlib_lib1" ||
          file=="aflowlib_lib2" ||
          file=="aflowlib_lib3" ||
          file=="aflowlib_lib4" ||
          file=="aflowlib_lib5" ||
          file=="aflowlib_lib6" ||
          file=="aflowlib_lib7" ||
          file=="aflowlib_lib8" ||
          file=="aflowlib_lib9" ||
          file=="aflowlib_icsd" )) { // from aflow_data
      if(file=="aflowlib_icsd") { // ICSD IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_ICSD;pvLibrary_tokens=&vLibrary_ICSD_tokens;
        init::InitGlobalObject("aflowlib_icsd",grep,LVERBOSE);
        if(XHOST_aflowlib_icsd.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_icsd,(*pvLibrary));
        found=TRUE; }
      if(file=="aflowlib_lib0") { // LIB0 IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_LIB0;pvLibrary_tokens=&vLibrary_LIB0_tokens;
        init::InitGlobalObject("aflowlib_lib0",grep,LVERBOSE);
        if(XHOST_aflowlib_lib0.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib0,(*pvLibrary));
        found=TRUE; }
      if(file=="aflowlib_lib1") { // LIB1 IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_LIB1;pvLibrary_tokens=&vLibrary_LIB1_tokens;
        init::InitGlobalObject("aflowlib_lib1",grep,LVERBOSE);
        if(XHOST_aflowlib_lib1.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib1,(*pvLibrary));
        found=TRUE; }
      if(file=="aflowlib_lib2") { // LIB2 IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_LIB2;pvLibrary_tokens=&vLibrary_LIB2_tokens;
        init::InitGlobalObject("aflowlib_lib2",grep,LVERBOSE);
        if(XHOST_aflowlib_lib2.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib2,(*pvLibrary));
        found=TRUE; }
      if(file=="aflowlib_lib3") { // LIB3 IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_LIB3;pvLibrary_tokens=&vLibrary_LIB3_tokens;
        init::InitGlobalObject("aflowlib_lib3",grep,LVERBOSE);
        if(XHOST_aflowlib_lib3.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib3,(*pvLibrary));
        found=TRUE; }
      if(file=="aflowlib_lib4") { // LIB4 IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_LIB4;pvLibrary_tokens=&vLibrary_LIB4_tokens;
        init::InitGlobalObject("aflowlib_lib4",grep,LVERBOSE);
        if(XHOST_aflowlib_lib4.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib4,(*pvLibrary));
        found=TRUE; }
      if(file=="aflowlib_lib5") { // LIB5 IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_LIB5;pvLibrary_tokens=&vLibrary_LIB5_tokens;
        init::InitGlobalObject("aflowlib_lib5",grep,LVERBOSE);
        if(XHOST_aflowlib_lib5.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib5,(*pvLibrary));
        found=TRUE; }
      if(file=="aflowlib_lib6") { // LIB6 IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_LIB6;pvLibrary_tokens=&vLibrary_LIB6_tokens;
        init::InitGlobalObject("aflowlib_lib6",grep,LVERBOSE);
        if(XHOST_aflowlib_lib6.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib6,(*pvLibrary));
        found=TRUE; }
      if(file=="aflowlib_lib7") { // LIB7 IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_LIB7;pvLibrary_tokens=&vLibrary_LIB7_tokens;
        init::InitGlobalObject("aflowlib_lib7",grep,LVERBOSE);
        if(XHOST_aflowlib_lib7.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib7,(*pvLibrary));
        found=TRUE; }
      if(file=="aflowlib_lib8") { // LIB8 IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_LIB8;pvLibrary_tokens=&vLibrary_LIB8_tokens;
        init::InitGlobalObject("aflowlib_lib8",grep,LVERBOSE);
        if(XHOST_aflowlib_lib8.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib8,(*pvLibrary));
        found=TRUE; }
      if(file=="aflowlib_lib9") { // LIB9 IS INSIDE init::InitGlobalObject
        pvLibrary=&vLibrary_LIB9;pvLibrary_tokens=&vLibrary_LIB9_tokens;
        init::InitGlobalObject("aflowlib_lib9",grep,LVERBOSE);
        if(XHOST_aflowlib_lib9.length()!=0) aurostd::string2vectorstring(XHOST_aflowlib_lib9,(*pvLibrary));
        found=TRUE; }
      // cerr << "(*pvLibrary_tokens).size()=" << (uint) (*pvLibrary_tokens).size() << endl;
      // cerr << "(*pvLibrary).size()=" << (*pvLibrary).size() << endl;
      if(found==TRUE) {
        if(LDEBUG || LVERBOSE) cerr << " ... size=" << (*pvLibrary).size() << " ..  ";// << endl;
      } else {
        cerr << "WWWWW  AFLOW_LIBRARY not found! " << endl;exit(0);
        return FALSE;
      }
      if((*pvLibrary_tokens).size()==0) {
        if(LDEBUG || LVERBOSE) cerr << "processing...   ";
        for(uint i=0;i<(*pvLibrary).size();i++) {
          tokens.clear();
          aurostd::string2tokens((*pvLibrary).at(i),tokens,"|");
          (*pvLibrary_tokens).push_back(tokens);
        }
        if(LDEBUG || LVERBOSE) cerr << "done." << endl;
      }
      // ----------------------------------------------------------------------
    }
    return (*pvLibrary).size();
  }
}

// ***************************************************************************
// aflowlib::LOAD_Library_ALL
// ***************************************************************************
namespace aflowlib {
  uint LOAD_Library_ALL(string options,string grep,bool LVERBOSE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(!aurostd::substring2bool(options,"no_aflowlib_icsd",TRUE) && vLibrary_ICSD.size()==0) LOAD_Library_LIBRARY("aflowlib_icsd",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib0",TRUE) && vLibrary_LIB0.size()==0) LOAD_Library_LIBRARY("aflowlib_lib0",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib1",TRUE) && vLibrary_LIB1.size()==0) LOAD_Library_LIBRARY("aflowlib_lib1",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib2",TRUE) && vLibrary_LIB2.size()==0) LOAD_Library_LIBRARY("aflowlib_lib2",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib3",TRUE) && vLibrary_LIB3.size()==0) LOAD_Library_LIBRARY("aflowlib_lib3",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib4",TRUE) && vLibrary_LIB4.size()==0) LOAD_Library_LIBRARY("aflowlib_lib4",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib5",TRUE) && vLibrary_LIB5.size()==0) LOAD_Library_LIBRARY("aflowlib_lib5",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib6",TRUE) && vLibrary_LIB6.size()==0) LOAD_Library_LIBRARY("aflowlib_lib6",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib7",TRUE) && vLibrary_LIB7.size()==0) LOAD_Library_LIBRARY("aflowlib_lib7",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib8",TRUE) && vLibrary_LIB8.size()==0) LOAD_Library_LIBRARY("aflowlib_lib8",grep,LVERBOSE);
    if(!aurostd::substring2bool(options,"no_aflowlib_lib9",TRUE) && vLibrary_LIB9.size()==0) LOAD_Library_LIBRARY("aflowlib_lib9",grep,LVERBOSE);
    if(LDEBUG) cerr << XPID << "vLibrary_ICSD.size()=" << vLibrary_ICSD.size() << " vLibrary_ICSD_tokens.size()=" << vLibrary_ICSD_tokens.size() << endl;
    if(LDEBUG) cerr << XPID << "vLibrary_LIB0.size()=" << vLibrary_LIB0.size() << " vLibrary_LIB0_tokens.size()=" << vLibrary_LIB0_tokens.size() << endl;
    if(LDEBUG) cerr << XPID << "vLibrary_LIB1.size()=" << vLibrary_LIB1.size() << " vLibrary_LIB1_tokens.size()=" << vLibrary_LIB1_tokens.size() << endl;
    if(LDEBUG) cerr << XPID << "vLibrary_LIB2.size()=" << vLibrary_LIB2.size() << " vLibrary_LIB2_tokens.size()=" << vLibrary_LIB2_tokens.size() << endl;
    if(LDEBUG) cerr << XPID << "vLibrary_LIB3.size()=" << vLibrary_LIB3.size() << " vLibrary_LIB3_tokens.size()=" << vLibrary_LIB3_tokens.size() << endl;
    if(LDEBUG) cerr << XPID << "vLibrary_LIB4.size()=" << vLibrary_LIB4.size() << " vLibrary_LIB4_tokens.size()=" << vLibrary_LIB4_tokens.size() << endl;
    if(LDEBUG) cerr << XPID << "vLibrary_LIB5.size()=" << vLibrary_LIB5.size() << " vLibrary_LIB5_tokens.size()=" << vLibrary_LIB5_tokens.size() << endl;
    if(LDEBUG) cerr << XPID << "vLibrary_LIB6.size()=" << vLibrary_LIB6.size() << " vLibrary_LIB6_tokens.size()=" << vLibrary_LIB6_tokens.size() << endl;
    if(LDEBUG) cerr << XPID << "vLibrary_LIB7.size()=" << vLibrary_LIB7.size() << " vLibrary_LIB7_tokens.size()=" << vLibrary_LIB7_tokens.size() << endl;
    if(LDEBUG) cerr << XPID << "vLibrary_LIB8.size()=" << vLibrary_LIB8.size() << " vLibrary_LIB8_tokens.size()=" << vLibrary_LIB8_tokens.size() << endl;
    if(LDEBUG) cerr << XPID << "vLibrary_LIB9.size()=" << vLibrary_LIB9.size() << " vLibrary_LIB9_tokens.size()=" << vLibrary_LIB9_tokens.size() << endl;
    // now make vLibrary_ALL
    if(vLibrary_ALL.size()==0) {
      vLibrary_ALL.clear();
      vector<string> tokens;
      if(LDEBUG || LVERBOSE) cerr << XPID << "00000  MESSAGE AFLOW LIBRARY  ALL  ";
      if(LDEBUG || LVERBOSE) cerr << "loading... ";
      for(uint i=0;i<vLibrary_ICSD.size();i++) vLibrary_ALL.push_back(vLibrary_ICSD.at(i));
      for(uint i=0;i<vLibrary_LIB0.size();i++) vLibrary_ALL.push_back(vLibrary_LIB0.at(i));
      for(uint i=0;i<vLibrary_LIB1.size();i++) vLibrary_ALL.push_back(vLibrary_LIB1.at(i));
      for(uint i=0;i<vLibrary_LIB2.size();i++) vLibrary_ALL.push_back(vLibrary_LIB2.at(i));
      for(uint i=0;i<vLibrary_LIB3.size();i++) vLibrary_ALL.push_back(vLibrary_LIB3.at(i));
      for(uint i=0;i<vLibrary_LIB4.size();i++) vLibrary_ALL.push_back(vLibrary_LIB4.at(i));
      for(uint i=0;i<vLibrary_LIB5.size();i++) vLibrary_ALL.push_back(vLibrary_LIB5.at(i));
      for(uint i=0;i<vLibrary_LIB6.size();i++) vLibrary_ALL.push_back(vLibrary_LIB6.at(i));
      for(uint i=0;i<vLibrary_LIB7.size();i++) vLibrary_ALL.push_back(vLibrary_LIB7.at(i));
      for(uint i=0;i<vLibrary_LIB8.size();i++) vLibrary_ALL.push_back(vLibrary_LIB8.at(i));
      for(uint i=0;i<vLibrary_LIB9.size();i++) vLibrary_ALL.push_back(vLibrary_LIB9.at(i));
      if(LDEBUG || LVERBOSE) cerr << "processing...   ";
      for(uint i=0;i<vLibrary_ALL.size();i++) {
        tokens.clear();
        aurostd::string2tokens(vLibrary_ALL.at(i),tokens,"|");
        vLibrary_ALL_tokens.push_back(tokens);
      }
      if(LDEBUG || LVERBOSE) cerr << "done." << endl;
    }
    if(LDEBUG) cerr << XPID << "vLibrary_ALL.size()=" << vLibrary_ALL.size() << " vLibrary_ALL_tokens.size()=" << vLibrary_ALL_tokens.size() << endl;
    return vLibrary_ALL.size();
  }
}

namespace aflowlib {
  uint LOAD_Library_ALL(string options,bool LVERBOSE) {
    return LOAD_Library_ALL(options,"",LVERBOSE);
  }
}

namespace aflowlib {
  uint LOAD_Library_ALL(bool LVERBOSE) {
    return LOAD_Library_ALL("","",LVERBOSE);
  }
}

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

// ***************************************************************************
// aflowlib::XLIB2RAW_CheckProjectFromDirectory
// ***************************************************************************
namespace aflowlib {
  string LIB2RAW_CheckProjectFromDirectory(const string& directory) {
    string soliloquy=XPID+"aflowlib::LIB2RAW_CheckProjectFromDirectory():";
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::LIB2RAW_CheckProjectFromDirectory" << endl;  
    CheckMaterialServer("aflowlib::LIB2RAW_CheckProjectFromDirectory"); // must be in AFLOW_MATERIALS_SERVER
    // find from PWD
    string PROJECT_LIBRARY="NOTHING",directory_pwd=aurostd::getPWD();aurostd::StringSubst(directory_pwd,"\n","");  //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")
    CleanDirectoryLIB(directory_pwd);
    //[CO20200624 - do inside CleanDirectoryLIB()]if(directory_pwd=="/common/GNDSTATE") directory_pwd="/common/LIB2"; // [HISTORIC]
    //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory_pwd,"common/SCINT","common/ICSD"); // [HISTORIC]
    //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory_pwd,"ELPASOLITES","AURO"); // PATCH  //CO20200624 - MOVED FROM BELOW, PROBABLY OBSOLETE
    //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory_pwd,"common/ELPASOLITES","common/AURO"); // [HISTORIC]
    //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory_pwd,"SCINT","ICSD"); // PATCH  //CO20200624 - MOVED FROM BELOW, PROBABLY OBSOLETE
    //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory_pwd,"LIB2/RAW","LIB2/LIB"); // [HISTORIC]

    // if not found by PWD switch to directory
    for(uint i=0;i<vAFLOW_PROJECTS_DIRECTORIES.size();i++){
      if(LDEBUG){cerr << soliloquy << " looking for " << vAFLOW_PROJECTS_DIRECTORIES.at(i) << " in " << directory << endl;}  //CO20200624
      if(aurostd::substring2bool(directory,vAFLOW_PROJECTS_DIRECTORIES.at(i))) PROJECT_LIBRARY=vAFLOW_PROJECTS_DIRECTORIES.at(i);
    }
    if(PROJECT_LIBRARY!="NOTHING") {
      if(LDEBUG) cerr << XPID << "aflowlib::LIB2RAW_CheckProjectFromDirectory: FOUND from directory: " << PROJECT_LIBRARY << endl;
      return PROJECT_LIBRARY;
    }

    // cerr << XPID << directory_pwd << endl;
    for(uint i=0;i<vAFLOW_PROJECTS_DIRECTORIES.size();i++) {
      if(LDEBUG) cerr << XPID << vAFLOW_PROJECTS_DIRECTORIES.at(i) << endl;
      if(aurostd::substring2bool(directory_pwd,vAFLOW_PROJECTS_DIRECTORIES.at(i))) PROJECT_LIBRARY=vAFLOW_PROJECTS_DIRECTORIES.at(i);
    }
    if(PROJECT_LIBRARY!="NOTHING") {
      if(LDEBUG) cerr << XPID << "aflowlib::LIB2RAW_CheckProjectFromDirectory: FOUND from pwd: " << PROJECT_LIBRARY << endl;
      return PROJECT_LIBRARY;
    }

    if(PROJECT_LIBRARY=="NOTHING") {
      cerr << XPID << "aflowlib::LIB2RAW_CheckProjectFromDirectory: Nothing found from pwd or directory ... exiting... [directory=" << directory << "]   [pwd=" << directory_pwd << "]" << endl;
      exit(0);
    }
    return XHOST.tmpfs;
  }
}

// ***************************************************************************
// aflowlib::XFIX_LIBRARY_ALL
// ***************************************************************************
namespace aflowlib {
  void XFIX_LIBRARY_ALL(string LIBRARY_IN,vector<string> argv) {
    CheckMaterialServer("aflowlib::XFIX_LIBRARY_ALL"); // must be in AFLOW_MATERIALS_SERVER
    stringstream aus_exec;
    string system=argv.at(2),structure=argv.at(3);   // if 3 inputs need to fix with fewer inputs
    string systemstructure=system+"/"+structure; // if 3 inputs need to fix with fewer inputs
    cout << "Fixing " << LIBRARY_IN << "/LIB/" << systemstructure << endl;
    aurostd::RemoveFile(string("_aflow_xfix_LIBRARY_"));
    aus_exec << "#!/bin/csh" << endl;
    aus_exec << "# Cleaning up " << LIBRARY_IN << " system/structure " << endl;
    aus_exec << "if !(-e " << LIBRARY_IN << "/FIX/" << system << ") then" << endl;
    aus_exec << "   mkdir " << LIBRARY_IN << "/FIX/" << system << endl;
    aus_exec << "endif" << endl;
    aus_exec << "mv " << LIBRARY_IN << "/LIB/" << systemstructure << " " << LIBRARY_IN << "/FIX/" << system << " " << endl;
    aus_exec << "mv " << LIBRARY_IN << "/RAW/" << systemstructure << " " << XHOST.tmpfs << "/" << system << " " << endl;
    aus_exec << "rm -rf " << LIBRARY_IN << "/RAW/" << systemstructure << endl;
    aurostd::stringstream2file(aus_exec,string("_aflow_xfix_LIBRARY_"));
    aurostd::ChmodFile("755",string("_aflow_xfix_LIBRARY_"));
    // aurostd::execute((char*) "cat ./_aflow_xfix_LIBRARY_");
    aurostd::execute((char*) "./_aflow_xfix_LIBRARY_");
    aurostd::RemoveFile(string("_aflow_xfix_LIBRARY_"));
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_ALL
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_ALL(string options,bool flag_FORCE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::LIB2RAW_ALL: BEGIN" << endl;
    if(LDEBUG) cerr << XPID << "aflowlib::LIB2RAW_ALL: options=" << options << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,","); // QUICK FIX

    if(LDEBUG)  cerr << XPID << "aflowlib::LIB2RAW_ALL: tokens.size()=" << tokens.size() << endl;
    if(tokens.size()>2) {
      init::ErrorOption(options,"aflowlib::LIB2RAW_ALL",aurostd::liststring2string("aflow --lib2raw=directory","aflow --lib2raw=all[,dir]"));
      exit(0);
    }
    if(tokens.size()>=1) {
      if(tokens.at(0)!="all")  {
        init::ErrorOption(options,"aflowlib::LIB2RAW_ALL",aurostd::liststring2string("aflow --lib2raw=directory","aflow --lib2raw=all[,dir]"));
        exit(0);
      }
    }

    CheckMaterialServer("aflowlib::LIB2RAW_ALL"); // must be in AFLOW_MATERIALS_SERVER
    string PROJECT_LIBRARY;
    if(tokens.size()==2) PROJECT_LIBRARY=aflowlib::LIB2RAW_CheckProjectFromDirectory(tokens.at(1));
    else PROJECT_LIBRARY=aflowlib::LIB2RAW_CheckProjectFromDirectory(aurostd::getPWD()); //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")
    cerr << XPID << "aflowlib::LIB2RAW_ALL FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;

    int multi_sh_value=XHOST.CPU_Cores;
    //ME20181109 - Handle NCPUS=MAX
    if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS") && !(XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX")))
      multi_sh_value=aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS"));

    cerr << XPID << "multi_sh_value=" << multi_sh_value << endl;
    //    exit(0);
    deque<string> dcmds;

    if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("ICSD") ||
        PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB0") ||
        PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB1") ||
        PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB2") ||
        PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB3") ||
        PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB4") ||
        PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB5") ||
        PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB6") ||
        PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB7") ||
        PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB8") ||
        PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB9")) {
      string command="find \""+PROJECT_LIBRARY+"/LIB/\" -name \"" + _AFLOWIN_ + "\" | sort ";
      string directory_list=aurostd::execute2string(command);
      vector<string> tokens;
      aurostd::string2tokens(directory_list,tokens,"\n");
      for(uint i=0;i<tokens.size();i++) {
        aurostd::StringSubst(tokens.at(i),"/"+_AFLOWIN_,"");
        aurostd::StringSubst(tokens.at(i),"/"+_AFLOWLOCK_,"");
        if(tokens.at(i).size()>=5 && tokens.at(i)[tokens.at(i).size()-5]=='/' && tokens.at(i)[tokens.at(i).size()-4]=='c' && tokens.at(i)[tokens.at(i).size()-3]=='o' && tokens.at(i)[tokens.at(i).size()-2]=='r' && tokens.at(i)[tokens.at(i).size()-1]=='e'){tokens.at(i)=tokens.at(i).substr(0,tokens.at(i).size()-5);} //aurostd::StringSubst(tokens.at(i),"/core",""); //CO20200624 - prevent /home/corey -> /homey
        string cmd="aflow";
        if(XHOST.vflag_control.flag("BEEP")) cmd+=" --beep";
        if(flag_FORCE) cmd+=" --force";
        cmd+=" --lib2raw="+ tokens.at(i);
        dcmds.push_back(cmd);
      };
      if(multi_sh_value==0 || multi_sh_value==1) {
        long double delta_seconds,eta_seconds,reference_seconds=aurostd::get_seconds();
        vector<long double> vdelta_seconds;
        for(uint i=0;i<dcmds.size();i++) {
          //  cout << delta_seconds << " " << dcmds.at(i) << endl;
          aflowlib::LIB2RAW(tokens.at(i),flag_FORCE);
          delta_seconds=aurostd::get_delta_seconds(reference_seconds);
          if(delta_seconds>1.0) {
            vdelta_seconds.push_back(delta_seconds);
            eta_seconds=aurostd::mean(vdelta_seconds)*(dcmds.size()-vdelta_seconds.size());
            cout << XPID << "aflowlib::LIB2RAW_ALL: [STEP]"
              << "  DONE= " << vdelta_seconds.size() << " / " << dcmds.size()-vdelta_seconds.size() << " (" << 100*vdelta_seconds.size()/dcmds.size()   << ")"
              << "  iSEC=" << vdelta_seconds.at(vdelta_seconds.size()-1)
              << "  aSEC=" << aurostd::mean(vdelta_seconds)
              << "  TO_DO=" << dcmds.size()-vdelta_seconds.size()
              //       << "  ETA(secs)=" << eta_seconds
              //       << "  ETA(mins)=" << eta_seconds/(60.0)
              //       << "  ETA(hours)=" << eta_seconds/(60*60)
              << "  ETA(days)=" << eta_seconds/(60*60*24)
              //       << "  ETA(weeks)=" << eta_seconds/(60*60*24*7)
              //       << "  ETA(years)=" << eta_seconds/(60*60*24*365)
              << endl;
          }
        }
      } else {
        // DO MULTI
        aurostd::multithread_execute(dcmds,multi_sh_value,FALSE);
      }
      return TRUE;
    }
    cerr << XPID << "aflowlib::LIB2RAW_ALL Project Not Found" << endl;exit(0);
    return FALSE;
  }
}

// ***************************************************************************
// aflowlib::GetSpeciesDirectory
// ***************************************************************************
namespace aflowlib {
  uint GetSpeciesDirectory(const string& directory,vector<string>& vspecies) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"aflowlib::GetSpeciesDirectory():";
    if(LDEBUG) cerr << soliloquy << " BEGIN" << endl;
    vspecies.clear();vector<string> vs,tokens;
    stringstream oss;
    //[CO20200624 - OBSOLETE]string temp_file=aurostd::TmpFileCreate("getspecies");

    if(XHOST.vext.size()!=XHOST.vcat.size()) { cerr << XPID << "ERROR - aflowlib::GetSpeciesDirectory: XHOST.vext.size()!=XHOST.vcat.size()" << endl;exit(0); }

    if(LDEBUG) cerr << soliloquy << " [1]" << endl;
    for(uint iext=0;iext<XHOST.vext.size();iext++) { // check _AFLOWIN_.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/" + _AFLOWIN_ + XHOST.vext.at(iext))) {
        aurostd::string2vectorstring(aurostd::execute2string(XHOST.vcat.at(iext)+" \""+directory+"/" + _AFLOWIN_ + XHOST.vext.at(iext)+"\""+" | grep VASP_POTCAR_FILE"),vspecies);
        for(uint i=0;i<vspecies.size();i++) aurostd::StringSubst(vspecies.at(i),"[VASP_POTCAR_FILE]","");
        for(uint i=0;i<vspecies.size();i++) KBIN::VASP_PseudoPotential_CleanName(vspecies.at(i));
      }
    }
    if(LDEBUG) cerr << soliloquy << " [2]" << endl;
    for(uint iext=0;iext<XHOST.vext.size();iext++) { // check POSCAR.orig.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/POSCAR.orig"+XHOST.vext.at(iext))) {
        oss.str(aurostd::execute2string(XHOST.vcat.at(iext)+" \""+directory+"/POSCAR.orig"+XHOST.vext.at(iext)+"\""));
        xstructure xstr(oss,IOVASP_POSCAR);
        if(xstr.species.size()>0) if(xstr.species.at(0)!="") for(uint i=0;i<xstr.species.size();i++) vspecies.push_back(xstr.species.at(i)); // dont change order
      }
    }
    if(LDEBUG) cerr << soliloquy << " [3]" << endl;
    for(uint iext=0;iext<XHOST.vext.size();iext++) {  // check POSCAR.relax1.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/POSCAR.relax1"+XHOST.vext.at(iext))) {
        oss.str(aurostd::execute2string(XHOST.vcat.at(iext)+" \""+directory+"/POSCAR.relax1"+XHOST.vext.at(iext)+"\""));
        xstructure xstr(oss,IOVASP_POSCAR);
        if(xstr.species.size()>0) if(xstr.species.at(0)!="") for(uint i=0;i<xstr.species.size();i++) vspecies.push_back(xstr.species.at(i)); // dont change order
      }
    }
    if(LDEBUG) cerr << soliloquy << " [4]" << endl;
    for(uint iext=0;iext<XHOST.vext.size();iext++) { // check POSCAR.bands.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/POSCAR.bands"+XHOST.vext.at(iext))) {
        oss.str(aurostd::execute2string(XHOST.vcat.at(iext)+" \""+directory+"/POSCAR.bands"+XHOST.vext.at(iext)+"\""));
        xstructure xstr(oss,IOVASP_POSCAR);
        if(xstr.species.size()>0) if(xstr.species.at(0)!="") for(uint i=0;i<xstr.species.size();i++) vspecies.push_back(xstr.species.at(i)); // dont change order
      }
    }
    if(LDEBUG) cerr << soliloquy << " [5]" << endl;
    for(uint iext=0;iext<XHOST.vext.size();iext++) { // check OUTCAR.relax1.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/OUTCAR.relax1"+XHOST.vext.at(iext))) {
        aurostd::string2vectorstring(aurostd::execute2string(XHOST.vcat.at(iext)+" \""+directory+"/OUTCAR.relax1"+XHOST.vext.at(iext)+"\""+" | grep TITEL"),vs);
        for(uint i=0;i<vs.size();i++) { aurostd::string2tokens(vs.at(i),tokens," ");vspecies.push_back(KBIN::VASP_PseudoPotential_CleanName(tokens.at(3))); }
      }
    }
    if(LDEBUG) cerr << soliloquy << " [6]" << endl;
    for(uint iext=0;iext<XHOST.vext.size();iext++) { // check OUTCAR.static.EXT
      if(!vspecies.size() && aurostd::FileExist(directory+"/OUTCAR.static"+XHOST.vext.at(iext))) {
        aurostd::string2vectorstring(aurostd::execute2string(XHOST.vcat.at(iext)+" \""+directory+"/OUTCAR.static"+XHOST.vext.at(iext)+"\""+" | grep TITEL"),vs);
        for(uint i=0;i<vs.size();i++) { aurostd::string2tokens(vs.at(i),tokens," ");vspecies.push_back(KBIN::VASP_PseudoPotential_CleanName(tokens.at(3))); }
      }
    }
    //  for(uint i=0;i<vspecies.size();i++) cerr << XPID << i << " - " << vspecies.at(i) << endl; exit(0);
    if(LDEBUG) cerr << soliloquy << " END" << endl;
    return vspecies.size();
  }
}

namespace aflowlib {
  void CleanDirectoryLIB(string& directory){  //CO20200624
    if(directory.empty()){return;}
    directory=aurostd::CleanFileName(directory);
    //[CO20200624 - do below]if(directory.at(directory.size()-1)=='/')  directory=directory.substr(0,directory.size()-1);
    aurostd::StringSubst(directory,"/"+_AFLOWIN_,"");
    aurostd::StringSubst(directory,"/"+_AFLOWLOCK_,"");
    if(directory.size()>=5 && directory[directory.size()-5]=='/' && directory[directory.size()-4]=='c' && directory[directory.size()-3]=='o' && directory[directory.size()-2]=='r' && directory[directory.size()-1]=='e'){directory=directory.substr(0,directory.size()-5);} //aurostd::StringSubst(directory,"/core",""); //CO20200624 - prevent /home/corey -> /homey
    aurostd::StringSubst(directory,"/common/GNDSTATE","/common/LIB2");  // [HISTORIC]
    aurostd::StringSubst(directory,"common/SCINT","common/ICSD");       // [HISTORIC]
    aurostd::StringSubst(directory,"SCINT","ICSD");                     // [HISTORIC]
    aurostd::StringSubst(directory,"common/ELPASOLITES","common/AURO"); // [HISTORIC]
    aurostd::StringSubst(directory,"ELPASOLITES","AURO");               // [HISTORIC]
    aurostd::StringSubst(directory,"LIBRARYX/RAW","LIBRARYX/LIB");      // [HISTORIC]
    aurostd::StringSubst(directory,"LIB2/RAW","LIB2/LIB");              // [HISTORIC]
    aurostd::StringSubst(directory,"/RAW/","/LIB/");                    // to help 
    while(directory.size()>0 && directory.at(directory.size()-1)=='/'){directory=directory.substr(0,directory.size()-1);} //CO20200624 - MOVED FROM BELOW
    if(directory=="." || directory.empty()) {directory=aurostd::getPWD();} //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW(const string& options,bool flag_FORCE,bool LOCAL) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"aflowlib::LIB2RAW():";
    long double seconds_begin=aurostd::get_seconds();
    if(LDEBUG) cerr << soliloquy << " BEGIN" << endl;
    if(LDEBUG) cerr << soliloquy << " options=" << options << endl;
    
    // Call to LIB2LIB function to run postprocessing - added by CT20181212
    bool perform_LIB2LIB=false; //CO20200624 - good for debugging the rest of lib2raw, lib2lib can be slow
    if(perform_LIB2LIB && !LIB2LIB(options, flag_FORCE, LOCAL)) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"aflowlib::LIB2LIB() failed",_RUNTIME_ERROR_);
      //[CO20200624 - OBSOLETE]cerr << soliloquy << " LIB2LIB failed" << endl;
      //[CO20200624 - OBSOLETE]exit(0);
    }
    
    vector<string> tokens;
    if(!aurostd::substring2bool(options,"POCC")) {
      aurostd::string2tokens(options,tokens,","); // QUICK FIX
    } else {
      tokens.clear();
      tokens.push_back(options);
    }

    if(LDEBUG) cerr << soliloquy << " tokens.size()=" << tokens.size() << endl;
    if(tokens.size()==0 || tokens.size()>2) {
      init::ErrorOption(options,soliloquy,aurostd::liststring2string("aflow --lib2raw=directory","aflow --lib2raw=all[,dir]"));
      return false;  //CO20200624
    }

    if(tokens.size()>=1) {
      if(tokens.at(0)=="all") {
        XHOST.sensors_allowed=FALSE;
        XHOST.vflag_pflow.flag("MULTI=SH",FALSE);
        return LIB2RAW_ALL(options,flag_FORCE);
        XHOST.sensors_allowed=TRUE;
      }
    }

    string directory_LIB="",directory_RAW="",directory_WEB="";
    bool flag_WEB=FALSE;
    bool flag_files_LIB=FALSE,flag_files_RAW=FALSE,flag_files_WEB=FALSE;
    string PROJECT_LIBRARY;
    if(LOCAL) {
      flag_FORCE=true;
      string directory=options;
      CleanDirectoryLIB(directory);
      //[CO20200624 - do inside CleanDirectoryLIB()]string directory=aurostd::CleanFileName(options);
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"./","");
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"/RAW/","/LIB/"); // to help 
      //[CO20200624 - do inside CleanDirectoryLIB()]if(directory=="." || directory.empty()) { directory=aurostd::getPWD(); } //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")
      flag_WEB=FALSE;
      flag_files_LIB=FALSE,flag_files_RAW=FALSE,flag_files_WEB=FALSE;
      directory_LIB=directory;
      directory_RAW=directory+"/RAW";
      directory_WEB=directory+"/WEB";
      PROJECT_LIBRARY=directory_LIB;
    } else {    //normal run
      //  cout << soliloquy << " AFLOW (" << VERSION << ")" << endl;
      CheckMaterialServer(soliloquy); // must be in AFLOW_MATERIALS_SERVER
      
      string directory=options;
      CleanDirectoryLIB(directory);
      //[CO20200624 - do inside CleanDirectoryLIB()]if(directory.at(directory.size()-1)=='/')  directory=directory.substr(0,directory.size()-1);
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"/"+_AFLOWIN_,"");
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"/"+_AFLOWLOCK_,"");
      //[CO20200624 - do inside CleanDirectoryLIB()]if(directory.size()>=5 && directory[directory.size()-5]=='/' && directory[directory.size()-4]=='c' && directory[directory.size()-3]=='o' && directory[directory.size()-2]=='r' && directory[directory.size()-1]=='e'){directory=directory.substr(0,directory.size()-5);} //aurostd::StringSubst(directory,"/core",""); //CO20200624 - prevent /home/corey -> /homey
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"common/SCINT","common/ICSD");   // [HISTORIC]
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"common/ELPASOLITES","common/AURO");   // [HISTORIC]
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"LIBRARYX/RAW","LIBRARYX/LIB");   // [HISTORIC]
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"LIB2/RAW","LIB2/LIB");

      PROJECT_LIBRARY=aflowlib::LIB2RAW_CheckProjectFromDirectory(directory);
      // cout << soliloquy << " FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
      cout << soliloquy << " directory=" << directory << endl;

      //  bool flag_LIB=FALSE;
      if(LDEBUG) cerr << soliloquy << " scan libraries BEGIN [1]" << endl;
      if(LDEBUG) cerr << soliloquy << " PROJECT_LIBRARY=" << PROJECT_LIBRARY << endl;
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_LIB0=" << XHOST_LIBRARY_LIB0 << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB0")) { flag_WEB=FALSE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_LIB1=" << XHOST_LIBRARY_LIB1 << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB1")) { flag_WEB=FALSE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_LIB2=" << XHOST_LIBRARY_LIB2 << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB2")) { flag_WEB=FALSE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_LIB3=" << XHOST_LIBRARY_LIB3 << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB3")) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_LIB4=" << XHOST_LIBRARY_LIB4 << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB4")) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_LIB5=" << XHOST_LIBRARY_LIB5 << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB5")) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_LIB6=" << XHOST_LIBRARY_LIB6 << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB6")) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_LIB7=" << XHOST_LIBRARY_LIB7 << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB7")) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_LIB8=" << XHOST_LIBRARY_LIB8 << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB8")) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_LIB9=" << XHOST_LIBRARY_LIB9 << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB9")) { flag_WEB=TRUE;flag_files_RAW=TRUE; }
      if(LDEBUG) cerr << soliloquy << " XHOST_LIBRARY_ICSD=" << XHOST_LIBRARY_ICSD << endl;
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("ICSD")) { flag_WEB=TRUE;flag_files_WEB=TRUE; }
      if(LDEBUG) cerr << soliloquy << " scan libraries END [2]"   << endl;

      stringstream aus_exec;
      //   vector<string> tokens;
      aurostd::string2tokens(directory,tokens,"/");
      directory_LIB.clear();directory_RAW.clear();
      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("ICSD")) {
        for(uint i=tokens.size()-1;i>0;i--)
          if(aurostd::substring2bool(tokens.at(i),"_ICSD_")) {
            if(i>=1) {
              directory_LIB=tokens.at(i-1)+"/"+tokens.at(i);
              //    directory_WEB=tokens.at(i-1);  removed to make things consistent
              directory_WEB=directory;
            }
          }
        if(directory_LIB.length()==0) {
          cerr << soliloquy << " FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
          cerr << soliloquy << " you must specify the directory including the whole lattice type" << endl;
          cerr << " such as  aflow --lib2raw=FCC/La1Se1_ICSD_27104  " << endl;
          return false;  //CO20200624
        }
      }

      // strip the directory_LIB of everything else
      directory_LIB=directory; // somewhere to start
      if(aurostd::substring2bool(directory_LIB,"LIB/")) directory_LIB=aurostd::substring2string(directory_LIB,"LIB/",FALSE);
      directory_RAW=directory_LIB;
      directory_WEB=directory_LIB;
      directory_LIB=aurostd::CleanFileName(PROJECT_LIBRARY+"/LIB/"+directory_LIB);
      while(directory_LIB.size()>0 && directory_LIB.at(directory_LIB.size()-1)=='/'){directory_LIB=directory_LIB.substr(0,directory_LIB.size()-1);} //CO20200624 - REPETITA IUVANT

      //CO20200624 - BELOW HERE DIRECTORY_LIB IS FIXED!!!!!!!!!
      
      directory_RAW=aurostd::CleanFileName(PROJECT_LIBRARY+"/RAW/"+directory_RAW);aurostd::StringSubst(directory_RAW,"RAW/LIB","RAW");
      directory_WEB=aurostd::CleanFileName(PROJECT_LIBRARY+"/WEB/"+directory_WEB);aurostd::StringSubst(directory_WEB,"WEB/LIB","WEB");
      if(flag_WEB==FALSE) directory_WEB=aurostd::CleanFileName(PROJECT_LIBRARY+"/WEB/");

      if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB0") ||
          PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB1") ||
          PROJECT_LIBRARY==init::AFLOW_Projects_Directories("LIB2")) {
        //      directory_WEB=directory_RAW;
      }
      cout << soliloquy << " PROJECT_LIBRARY=" << PROJECT_LIBRARY << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("LIB0")=" << init::AFLOW_Projects_Directories("LIB0") << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("LIB1")=" << init::AFLOW_Projects_Directories("LIB1") << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("LIB2")=" << init::AFLOW_Projects_Directories("LIB2") << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("LIB3")=" << init::AFLOW_Projects_Directories("LIB3") << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("LIB4")=" << init::AFLOW_Projects_Directories("LIB4") << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("LIB5")=" << init::AFLOW_Projects_Directories("LIB5") << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("LIB6")=" << init::AFLOW_Projects_Directories("LIB6") << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("LIB7")=" << init::AFLOW_Projects_Directories("LIB7") << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("LIB8")=" << init::AFLOW_Projects_Directories("LIB8") << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("LIB9")=" << init::AFLOW_Projects_Directories("LIB9") << endl;
      // cout << soliloquy << " init::AFLOW_Projects_Directories("ICSD")=" << init::AFLOW_Projects_Directories("ICSD") << endl;
    }
    cout << soliloquy << " directory_LIB=" << directory_LIB << endl;
    cout << soliloquy << " directory_RAW=" << directory_RAW << endl;
    cout << soliloquy << " directory_WEB=" << directory_WEB << endl;

    if(!directory_LIB.empty()) {XHOST.vflag_control.flag("DIRECTORY_CLEAN",TRUE);XHOST.vflag_control.push_attached("DIRECTORY_CLEAN",directory_LIB);} //CO20200624 - fix error messages to point to directory_LIB

    bool perform_LOCK=FALSE;  //CO20200624 - turning off in general, check below
    bool perform_BANDS=FALSE,perform_BADER=FALSE,perform_MAGNETIC=FALSE,perform_THERMODYNAMICS=FALSE;
    bool perform_AGL=FALSE,perform_AEL=FALSE;
    bool perform_POCC=FALSE;  //CO20200624
    bool perform_PATCH=FALSE; // to inject updates while LIB2RAW  //CO20200624 - turning off in general, check below

    //CO20200624 - do NOT use EFileExist(), we want to EXPLICITLY catch those that are zipped
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      if(aurostd::FileExist(directory_LIB+"/"+_AFLOWLOCK_+XHOST.vext.at(iext))) perform_LOCK=TRUE;  //CO20200624
      if(aurostd::FileExist(directory_LIB+"/OUTCAR.relax1"+XHOST.vext.at(iext)) ||
          aurostd::FileExist(directory_LIB+"/OUTCAR.relax2"+XHOST.vext.at(iext)) ||
          aurostd::FileExist(directory_LIB+"/OUTCAR.relax3"+XHOST.vext.at(iext)) || //CO20180827 
          aurostd::FileExist(directory_LIB+"/OUTCAR.static"+XHOST.vext.at(iext)) || //CO20180827
          FALSE) perform_THERMODYNAMICS=TRUE; //CO20180827
      if(aurostd::FileExist(directory_LIB+"/OUTCAR.static"+XHOST.vext.at(iext))) perform_MAGNETIC=TRUE;
      if(aurostd::FileExist(directory_LIB+"/OUTCAR.bands"+XHOST.vext.at(iext))) perform_BANDS=TRUE;
      if(aurostd::FileExist(directory_LIB+"/AECCAR0.static"+XHOST.vext.at(iext)) && aurostd::FileExist(directory_LIB+"/AECCAR2.static"+XHOST.vext.at(iext))) perform_BADER=TRUE;
      if(aurostd::FileExist(directory_LIB+"/aflow.agl.out"+XHOST.vext.at(iext))) perform_AGL=TRUE;
      if(aurostd::FileExist(directory_LIB+"/aflow.ael.out"+XHOST.vext.at(iext))) perform_AEL=TRUE;
      if(aurostd::FileExist(directory_LIB+"/"+POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE+XHOST.vext.at(iext))) perform_POCC=TRUE; //CO20200624
    }
    if((perform_THERMODYNAMICS || perform_BANDS || perform_MAGNETIC)){perform_PATCH=true;} //CO20200624

    // override for the web
    //    if(PROJECT_LIBRARY==init::AFLOW_Projects_Directories("ICSD")) directory_WEB=aurostd::CleanFileName(PROJECT_AFLOWLIB_WEB);

    // cerr << XPID << directory_WEB << endl; cerr << XPID << directory_LIB << endl; cerr << XPID << directory_RAW << endl; exit(0);
    if(!aurostd::FileExist(directory_LIB+"/"+_AFLOWIN_)) {
      cout << soliloquy << " FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
      cout << soliloquy << " directory does not exist: " << directory_LIB << endl;
      return FALSE;
      // exit(0);
    }
    if(flag_FORCE==FALSE) {
      //CO20200624 - checking files in RAW from previous lib2raw run
      if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) {  // directory_RAW+"/"+_AFLOWIN_)
        _aflags aflags;
        aurostd::ZIP2ZIP(directory_LIB,"bz2","xz",FALSE,XHOST.sPID); // PATCH FOR REFRESH (to be removed)
        aurostd::ZIP2ZIP(directory_LIB,"gz","xz",FALSE,XHOST.sPID); // PATCH FOR REFRESH (to be removed)

        cout << soliloquy << " ALREADY CALCULATED = " << directory_RAW << "   END_DATE - [v=" << string(AFLOW_VERSION) << "] -" << Message(aflags,"time",_AFLOW_FILE_NAME_) << endl;
        return FALSE;
      }
      if(perform_BANDS) {
        if(aurostd::FileExist(directory_RAW+"/EIGENVAL.bands") || aurostd::EFileExist(directory_RAW+"/EIGENVAL.bands")) { //CO20200624 - this check does not make sense, we don't rely on EIGENVAL.bands, we need OUTCAR.bands!
          // return FALSE;
          cout << soliloquy << " directory is skipped because of BANDS: " << directory_RAW << endl;
          return FALSE;
          //  exit(0);
        }
      }
    }
    // directory_LIB exists and directory_RAW does not exist can move on:
    if(aurostd::RemoveWhiteSpaces(directory_RAW).empty()) {  //imagine this empty, then do this: "rm -f /"+directory_RAW+"/*" DEAD
      cout << soliloquy << " directory_RAW is empty" << endl;
      return FALSE;
    }
    //if(LOCAL){aurostd::execute("rm -rf \""+directory_RAW)+"\"";} //CO20190321 - safer NOT to -rf
    //{ aurostd::DirectoryMake(directory_RAW);aurostd::execute("rm -f \"/"+directory_RAW+"/*\""); } //CO20190321 - the prepending of '/' is BAD, only works for SC but can kill a local users whole computer if "/"+directory_RAW exists
    {
      aurostd::DirectoryMake(directory_RAW);
      aurostd::RemoveFile(directory_RAW+"/*");
      // [OBSOLETE] aurostd::execute("rm -f \""+directory_RAW+"/*\"");
    }
    if(!aurostd::IsDirectory(directory_RAW)) {
      cout << soliloquy << " directory is skipped because directory_RAW cannot be created: " << directory_RAW << endl;
      return FALSE;
    }
    if(flag_WEB) {
      if(aurostd::RemoveWhiteSpaces(directory_WEB).empty()) {  //imagine this empty, then do this: "rm -f /"+directory_WEB+"/*" DEAD
        cout << soliloquy << " directory_WEB is empty" << endl;
        return FALSE;
      }
      //if(LOCAL){aurostd::execute("rm -rf \""+directory_WEB+"\"");} //CO20190321 - safer NOT to -rf
      //aurostd::DirectoryMake(directory_WEB);aurostd::execute("rm -f \"/"+directory_WEB+"/*\""); //CO20190321 - the prepending of '/' is BAD, only works for SC but can kill a local users whole computer if "/"+directory_RAW exists
      aurostd::DirectoryMake(directory_WEB);
      aurostd::RemoveFile(directory_WEB+"/*");
      // [OBSOLETE] aurostd::execute("rm -f \""+directory_WEB+"/*\"");
      if(!aurostd::IsDirectory(directory_WEB)) {
        cout << soliloquy << " directory is skipped because directory_WEB cannot be created: " << directory_WEB << endl;
        return FALSE;
      }
    }

    cout << soliloquy << " FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
    
    //[CO20200624 - OBSOLETE]if((perform_THERMODYNAMICS || perform_BANDS ||  perform_MAGNETIC))
    _aflags aflags;
    cout << soliloquy << " dir=" << directory_LIB << " BEGIN_DATE = " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    aurostd::ZIP2ZIP(directory_LIB,"bz2","xz",TRUE,XHOST.sPID); 
    aurostd::ZIP2ZIP(directory_LIB,"gz","xz",TRUE,XHOST.sPID);

    if(LDEBUG) cerr << soliloquy << " [1]" << endl;
    //DX20190516 [OBSOLETE - overwrites file everytime, messes up timestamp] aurostd::RemoveBinaryCharactersFromFile(directory_LIB,_AFLOWIN_); //DX20190211 - Ensure no binary characters in aflow.in
    //DX20190516 [ADD THIS LINE ONCE TESTED WITH REAL-WORLD CASE] if(aurostd::RemoveControlCodeCharactersFromFile(directory_LIB,_AFLOWIN_)) { //DX20190516 - Ensure no control characters in aflow.in; only modifies if control characters are detected 
    //DX20190516 [ADD THIS LINE ONCE TESTED WITH REAL-WORLD CASE]   cout << soliloquy << " Control characters detected in aflow.in. The invalid characters have been removed, and the file has been overwritten (the original file has moved to aflow.in_old)." << endl; // signals file has been cleaned
    //DX20190516 [ADD THIS LINE ONCE TESTED WITH REAL-WORLD CASE] }

    vector<string> vfile;   // the needed files
    aflowlib::_aflowlib_entry aflowlib_data;
    //[CO20200624 - OBSOLETE]vector<string> vspecies;   // the species
    GetSpeciesDirectory(directory_LIB,aflowlib_data.vspecies);
    aflowlib_data.species=aurostd::joinWDelimiter(aflowlib_data.vspecies,",");

    if(LDEBUG) cerr << soliloquy << " [2]" << endl;

    //[CO20200624 - OBSOLETE]for(uint i=0;i<aflowlib_data.vspecies.size();i++) {
    //[CO20200624 - OBSOLETE]  aflowlib_data.species+=aflowlib_data.vspecies.at(i);
    //[CO20200624 - OBSOLETE]  if(i<vspecies.size()-1)  aflowlib_data.species+=",";
    //[CO20200624 - OBSOLETE]}    // plug vspecies into aflowlib_data
    

    string system_name=KBIN::ExtractSystemName(directory_LIB);
    
    if(LDEBUG) cerr << soliloquy << " [3]" << endl;

    if(LOCAL) {
      aflowlib_data.aurl=aflowlib_data.auid=directory_LIB; //dummy
    } else {
      // build aflowlib_data.aurl
      aflowlib_data.aurl="aflowlib.duke.edu:"+directory_LIB;
      
      //CO REMOVE ME - DEBUGGING ONLY!!!!!!
      aurostd::StringSubst(aflowlib_data.aurl,"/mnt/MAIN/STAGING/qrats_finished_runs/LIB6/carbide","common"); //CO REMOVE ME - DEBUGGING ONLY!!!!!!
      aurostd::StringSubst(aflowlib_data.aurl,"/home/corey/common","common"); //CO REMOVE ME - DEBUGGING ONLY!!!!!!
      //CO REMOVE ME - DEBUGGING ONLY!!!!!!

      if(LDEBUG){cerr << soliloquy << " aurl(PRE )=" << aflowlib_data.aurl << endl;}

      if(aurostd::substring2bool(aflowlib_data.aurl,"LIBRARYX")) { aurostd::StringSubst(aflowlib_data.aurl,"common/GNDSTATE/LIBRARYX/LIB","AFLOWDATA/LIB2_RAW");aflowlib_data.catalog="LIBRARYX"; } // [HISTORIC]
      if(aurostd::substring2bool(aflowlib_data.aurl,"ICSD")) { aurostd::StringSubst(aflowlib_data.aurl,"common/ICSD/LIB","AFLOWDATA/ICSD_WEB");aflowlib_data.catalog="ICSD"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB0")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB0/LIB","AFLOWDATA/LIB0_RAW");aflowlib_data.catalog="LIB0"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB1")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB1/LIB","AFLOWDATA/LIB1_RAW");aflowlib_data.catalog="LIB1"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB2")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB2/LIB","AFLOWDATA/LIB2_RAW");aflowlib_data.catalog="LIB2"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB3")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB3/LIB","AFLOWDATA/LIB3_RAW");aflowlib_data.catalog="LIB3"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB4")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB4/LIB","AFLOWDATA/LIB4_RAW");aflowlib_data.catalog="LIB4"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB5")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB5/LIB","AFLOWDATA/LIB5_RAW");aflowlib_data.catalog="LIB5"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB6")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB6/LIB","AFLOWDATA/LIB6_RAW");aflowlib_data.catalog="LIB6"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB7")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB7/LIB","AFLOWDATA/LIB7_RAW");aflowlib_data.catalog="LIB7"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB8")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB8/LIB","AFLOWDATA/LIB8_RAW");aflowlib_data.catalog="LIB8"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB9")) { aurostd::StringSubst(aflowlib_data.aurl,"common/LIB9/LIB","AFLOWDATA/LIB9_RAW");aflowlib_data.catalog="LIB9"; }
      if(aurostd::substring2bool(aflowlib_data.aurl,"AURO")) { aurostd::StringSubst(aflowlib_data.aurl,"common/AURO/LIB","AFLOWDATA/AURO_RAW");aflowlib_data.catalog="AURO"; }
      aurostd::StringSubst(aflowlib_data.aurl,":/AFLOWDATA",":AFLOWDATA");
      
      // cout << soliloquy << " AURL = " << aurostd::PaddedPOST(aflowlib_data.aurl,60) << endl;//"   " << directory_LIB << endl;  //CO20181226

      if(LDEBUG){cerr << soliloquy << " aurl(POST)=" << aflowlib_data.aurl << endl;}

      // build aflowlib_data.auid
      
      if(LDEBUG) cerr << soliloquy << " [AUID=0] directory_LIB=" << directory_LIB << endl;
      if(perform_POCC){ //CO20200624
        string metadata_auid_json=aflowlib_data.directory2MetadataAUIDjsonfile(directory_LIB);
        aurostd::string2file(metadata_auid_json,aurostd::CleanFileName(directory_RAW+"/"+"metadata_auid.json"));
      }else{
        // aflowlib::directory2auid(directory_LIB,aflowlib_data.aurl,aflowlib_data.auid,aflowlib_data.vauid); // OLD STYLE
        aflowlib_data.directory2auid(directory_LIB); // NEW STLYE
      }
      //  if(AFLOWLIB_VERBOSE)
      if(LDEBUG) cerr << soliloquy << " [AUID=0b] directory_LIB=" << directory_LIB << endl;
      cout << soliloquy << " AURL  = " << aurostd::PaddedPOST(aflowlib_data.aurl,60) << endl;//"   " << directory_LIB << endl; 
      cout << soliloquy << " AUID  = " << aurostd::PaddedPOST(aflowlib_data.auid,60) << endl;//"   " << directory_LIB << endl;
      cout << soliloquy << " VAUID = " << aflowlib::auid2directory(aflowlib_data.auid) << endl;
      //ME20200207 - the system name is the canonical title
      aflowlib_data.title = system_name;  //ME20200207  //CO20200624

      //[OBSOLETE] //ME20190125 BEGIN - Build the title of the calculation
      //[OBSOLETE] vector<string> tokens;
      //[OBSOLETE] aurostd::string2tokens(aflowlib_data.aurl, tokens, "/");
      //[OBSOLETE] if (aurostd::substring2bool(aflowlib_data.aurl, "ICSD")) {
      //[OBSOLETE]   aflowlib_data.title = tokens.back();
      //[OBSOLETE] } else {
      //[OBSOLETE]   //ME20190821 - Made more general
      //[OBSOLETE]   // [OBSOLETE] aflowlib_data.title = tokens[tokens.size() - 2] + "." + tokens[tokens.size() - 1];
      //[OBSOLETE]   uint i = 0;
      //[OBSOLETE]   vector<string> vtitle;
      //[OBSOLETE]   for (i = 0; i < tokens.size(); i++) {
      //[OBSOLETE]     if (tokens[i] == "AFLOWDATA") break;
      //[OBSOLETE]   }
      //[OBSOLETE]   for (uint j = i + 2; j < tokens.size(); j++) {  // i + 2 to skip AFLOWDATA/LIBX_RAW/
      //[OBSOLETE]     vtitle.push_back(tokens[j]);
      //[OBSOLETE]   }
      //[OBSOLETE]   aflowlib_data.title = aurostd::joinWDelimiter(vtitle, ".");
      //[OBSOLETE] }
      //[OBSOLETE] //ME20190125 END
      if(LDEBUG) cerr << soliloquy << " [AUID=2]" << endl;
      cout << soliloquy << " CATALOG = " << aurostd::PaddedPOST(aflowlib_data.catalog,60) << endl;//"   " << directory_LIB << endl;
      if(LDEBUG) cerr << soliloquy << " [AUID=3]" << endl;
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the THERMNODYNAMYCS
    if(perform_THERMODYNAMICS) {
      cout << soliloquy << " THERMODYNAMIC LOOP ---------------------------------------------------------------------------------" << endl;//perform_BANDS=FALSE;
      aflowlib::LIB2RAW_Loop_Thermodynamics(directory_LIB,directory_RAW,vfile,aflowlib_data,soliloquy+" (thermodynamics):",LOCAL); // identifier inside
    }
    if(LDEBUG) cerr << soliloquy << " [4]" << endl;
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the BANDS
    if(perform_BANDS) {
      cout << soliloquy << " BANDS LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_Bands(directory_LIB,directory_RAW,vfile,aflowlib_data,soliloquy+" (bands):");
      // MOVE/LINK PICS data
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the MAGNETIC
    if((perform_MAGNETIC || perform_BANDS)) { //JX
      cout << soliloquy << " MAGNETIC LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_Magnetic(directory_LIB,directory_RAW,vfile,aflowlib_data,soliloquy+" (magnetic):"); 
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the BADER
    //      if(aurostd::substring2bool(aflowlib_data.aurl,"LIB6")) perform_BADER=FALSE; // hack
    // [OBSOLETE]    if(0) // DEBUG NOW
    if(perform_BADER) {
      cout << soliloquy << " BADER LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_Bader(directory_LIB,directory_RAW,vfile,aflowlib_data,soliloquy+" (bader):");
      // MOVE/LINK PICS data
      // [OBSOLETE] [MOVED DOWN] if(flag_WEB) {
      // [OBSOLETE] [MOVED DOWN] aurostd::LinkFile(directory_RAW+"/*_abader.out",directory_WEB);      // LINK
      // [OBSOLETE] [MOVED DOWN] aurostd::LinkFile(directory_RAW+"/*jvxl",directory_WEB);            // LINK
      // [OBSOLETE] [MOVED DOWN] }
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the AGL
    if(perform_AGL) {
      cout << soliloquy << " AGL LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_AGL(directory_LIB,directory_RAW,vfile,aflowlib_data,soliloquy+" (agl):");
      if(flag_WEB) {
        if(aurostd::FileExist(directory_RAW+"/aflow.agl.out")) aurostd::LinkFile(directory_RAW+"/aflow.agl.out",directory_WEB);    // LINK  //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AGL.out")) aurostd::LinkFile(directory_RAW+"/AGL.out",directory_WEB);    // LINK  //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AGL_energies_temperature.out")) aurostd::LinkFile(directory_RAW+"/AGL_energies_temperature.out",directory_WEB);    // LINK  //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AGL_thermal_properties_temperature.out")) aurostd::LinkFile(directory_RAW+"/AGL_thermal_properties_temperature.out",directory_WEB);    // LINK  //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AGL_edos_gap_pressure.out")) aurostd::LinkFile(directory_RAW+"/AGL_edos_gap_pressure.out",directory_WEB);    // LINK //CT20181212 //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AGL_edos_gap_pressure.json")) aurostd::LinkFile(directory_RAW+"/AGL_edos_gap_pressure.json",directory_WEB);    // LINK //CT20181212 //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AGL_energy.json")) aurostd::LinkFile(directory_RAW+"/AGL_energy.json",directory_WEB);    // LINK //CT20181212 //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AGL_energy_structures.json")) aurostd::LinkFile(directory_RAW+"/AGL_energy_structures.json",directory_WEB);    // LINK //CT20181212 //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AGL_energy_volume.out")) aurostd::LinkFile(directory_RAW+"/AGL_energy_volume.out",directory_WEB);    // LINK //CT20181212 //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AGL_gibbs_energy_pT.out")) aurostd::LinkFile(directory_RAW+"/AGL_gibbs_energy_pT.out",directory_WEB);    // LINK //CT20181212 //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AGL_Hugoniot.out")) aurostd::LinkFile(directory_RAW+"/AGL_Hugoniot.out",directory_WEB);    // LINK //CT20181212 //CO20200624 - adding FileExist() check
      }
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the AEL
    if(perform_AEL) {
      cout << soliloquy << " AEL LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_AEL(directory_LIB,directory_RAW,vfile,aflowlib_data,soliloquy+" (ael):");
      if(flag_WEB) {
        if(aurostd::FileExist(directory_RAW+"/aflow.ael.out")) aurostd::LinkFile(directory_RAW+"/aflow.ael.out",directory_WEB);    // LINK  //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AEL_Elastic_constants.out")) aurostd::LinkFile(directory_RAW+"/AEL_Elastic_constants.out",directory_WEB);    // LINK  //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AEL_Compliance_tensor.out")) aurostd::LinkFile(directory_RAW+"/AEL_Compliance_tensor.out",directory_WEB);    // LINK  //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AEL_elastic_tensor.json")) aurostd::LinkFile(directory_RAW+"/AEL_elastic_tensor.json",directory_WEB);    // LINK //CT20181212 //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/AEL_energy_structures.json")) aurostd::LinkFile(directory_RAW+"/AEL_energy_structures.json",directory_WEB);    // LINK //CT20181212 //CO20200624 - adding FileExist() check
      }
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the POCC
    if(perform_POCC) {
      cout << soliloquy << " POCC LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_POCC(directory_LIB,directory_RAW,vfile,aflowlib_data,soliloquy+" (POCC):"); 
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the LOCK
    if(perform_LOCK) {
      cout << soliloquy << " LOCK LOOP ---------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_LOCK(directory_LIB,directory_RAW,vfile,aflowlib_data,soliloquy+" (LOCK):"); 
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // do the PATCH
    if(perform_PATCH) {
      cout << soliloquy << " PATCH LOOP --------------------------------------------------------------------------------" << endl;
      aflowlib::LIB2RAW_Loop_PATCH(directory_LIB,directory_RAW,vfile,aflowlib_data,soliloquy+" (PATCH):"); 
    }
    // ---------------------------------------------------------------------------------------------------------------------------------
    // write DOS + BANDS JSON //CO20171025
    if((aurostd::FileExist(directory_LIB+"/DOSCAR.static") || aurostd::EFileExist(directory_LIB+"/DOSCAR.static")) &&
        (aurostd::FileExist(directory_LIB+"/POSCAR.static") || aurostd::EFileExist(directory_LIB+"/POSCAR.static"))) {
      cout << soliloquy << " DOS + BANDS JSON  ---------------------------------------------------------------------------------" << endl;
      stringstream _json_file,json_file;
      aurostd::xoption vpflow;  //dummy
      if(estructure::DOSDATA_JSON(vpflow,directory_LIB,_json_file,true)) {
        json_file << _json_file.str() << endl; _json_file.str("");
        cout << soliloquy << " compressing file: " << string(directory_RAW+"/"+system_name+"_dosdata.json") << endl; cout.flush();
        aurostd::stringstream2compressfile(DEFAULT_KZIP_BIN,json_file,directory_RAW+"/"+system_name+"_dosdata.json");
        json_file.str("");
        if(aurostd::EFileExist(directory_LIB+"/EIGENVAL.bands") && aurostd::EFileExist(directory_LIB+"/KPOINTS.bands")) {
          if(estructure::BANDSDATA_JSON(vpflow,directory_LIB,_json_file,true)) {
            json_file << _json_file.str() << endl; _json_file.str("");
            cout << soliloquy << " compressing file: " << string(directory_RAW+"/"+system_name+"_bandsdata.json") << endl; cout.flush();
            aurostd::stringstream2compressfile(DEFAULT_KZIP_BIN,json_file,directory_RAW+"/"+system_name+"_bandsdata.json");
            json_file.str("");
          }
        }
      }
    }
    

    // ---------------------------------------------------------------------------------------------------------------------------------
    // DO THE COMPRESSING OF VASP FILES
    //      cout << "COMPRESSING" << endl;
    cout << soliloquy << " COMPRESSING REMOVING LINKING --------------------------------------------------------------" << endl;

    // generic for compressing and linking
    deque<string> vtype; aurostd::string2tokens(".orig,.relax,.relax1,.relax2,.relax3,.relax4,.static,.bands",vtype,",");
    deque<string> vout; aurostd::string2tokens(".out,.json",vout,",");

    //CO20200624 - removed COMPRESSED eps
    for(uint iext=1;iext<XHOST.vext.size();iext++) {
      aurostd::RemoveFile(directory_RAW+"/*.eps"+XHOST.vext.at(iext));
    }

    // [OBSOLETE] for(uint iext=0;iext<XHOST.vext.size();iext++) {  //[CO20200106 - close bracket for indenting]}
    // [OBSOLETE] 	//  aurostd::RemoveFile(directory_RAW+"/core*");
    // [OBSOLETE] 	  for(uint ifile=0;ifile<vfile.size();ifile++) {
    // [OBSOLETE] 	  cout << vfile.at(ifile) << endl;
    // [OBSOLETE]    if(!aurostd::substring2bool(vfile.at(ifile),XHOST.vext.at(iext)) && aurostd::FileExist(directory_RAW+"/"+vfile.at(ifile))) {
    // [OBSOLETE] 	    if(LOCAL) { //CO20171025
    // [OBSOLETE] 	      if(aurostd::EFileExist(directory_LIB+"/"+vfile.at(ifile))) { aurostd::RemoveFile(directory_RAW+"/"+vfile.at(ifile)); } //it's all sitting in the directory above, wasteful
    // [OBSOLETE] 	      if(vfile.at(ifile)=="KPOINTS.relax") { aurostd::RemoveFile(directory_RAW+"/"+vfile.at(ifile)); }
    // [OBSOLETE] 	      if(vfile.at(ifile)=="EIGENVAL.bands.old") { aurostd::RemoveFile(directory_RAW+"/"+vfile.at(ifile)); }
    // [OBSOLETE] 	      if(vfile.at(ifile)=="OUTCAR.relax") { aurostd::RemoveFile(directory_RAW+"/"+vfile.at(ifile)); }
    // [OBSOLETE] 	    } 
    // [OBSOLETE] 	} // vfile
    // [OBSOLETE]  } // iext

    // FILES to remove
    deque<string> vfile2remove; aurostd::string2tokens("KPOINTS.bands.old,EIGENVAL.bands.old,OUTCAR.relax1",vfile2remove,","); //CO20200404 - removed OUTCAR.bands from this list, we need for Egap. EIGENVAL.bands does not have occupancies
    vfile2remove.push_back("aflow.pgroupk_xtal.out"); // comes from nowhere DX - we have variants for each VASP run (relax, static, bands), so this is NOW redundant
    for(uint iremove=0;iremove<vfile2remove.size();iremove++) {
      if(aurostd::FileExist(directory_RAW+"/"+vfile2remove.at(iremove))) { // need to be present
        //if(LDEBUG)
        cout << soliloquy << " removing file: " << string(directory_RAW+"/"+vfile2remove.at(iremove)) << endl;
        aurostd::RemoveFile(directory_RAW+"/"+vfile2remove.at(iremove));
      } // FileExist
    } // iremove
    // FILES to compress if not compressed already (linked), in such case they will be deleted.
    deque<string> vfile2compress0; aurostd::string2tokens("OUTCAR.relax",vfile2compress0,",");
    for(uint icompress=0;icompress<vfile2compress0.size();icompress++) {
      if(aurostd::FileExist(directory_RAW+"/"+vfile2compress0.at(icompress)+"."+DEFAULT_KZIP_BIN)) { // test if there is one already compressed (possibly linked)
        if(1||LDEBUG) { cout << soliloquy << " found compressed file: " << string(directory_RAW+"/"+vfile2compress0.at(icompress)+"."+DEFAULT_KZIP_BIN) << endl;  cout.flush(); }
        if(1||LDEBUG) { cout << soliloquy << " removing file: " << string(directory_RAW+"/"+vfile2compress0.at(icompress)) << endl;  cout.flush(); }
        aurostd::RemoveFile(directory_RAW+"/"+vfile2compress0.at(icompress));
      }
      if(aurostd::FileExist(directory_RAW+"/"+vfile2compress0.at(icompress))) { // need to be present
        if(1||LDEBUG) { cout << soliloquy << " compressing file: " << string(directory_RAW+"/"+vfile2compress0.at(icompress)) << endl;  cout.flush(); }
        aurostd::CompressFile(directory_RAW+"/"+vfile2compress0.at(icompress),DEFAULT_KZIP_BIN);
      } // FILES to compress
    } // icompress
    // [OBSOLETE] no compress jvxl     if(perform_BADER) {
    // [OBSOLETE] no compress jvxl   if(LDEBUG) { cout << soliloquy << " compressing file: " << string(directory_RAW+"/"+"*.jvxl") << endl;  cout.flush(); }
    // [OBSOLETE] no compress jvxl   aurostd::CompressFile(directory_RAW+"/"+"*.jvxl",DEFAULT_KZIP_BIN);
    // [OBSOLETE] no compress jvxl     }
    // FILES to compress
    deque<string> vfile2compress1; aurostd::string2tokens("aflow.pgroup,aflow.pgroup_xtal,aflow.pgroupk,aflow.pgroupk_xtal,aflow.pgroupk_Patterson,aflow.fgroup,aflow.iatoms,aflow.agroup",vfile2compress1,","); //DX20200520 - added Patterson
    for(uint ilink=0;ilink<vfile2compress1.size();ilink++) {
      for(uint iout=0;iout<vout.size();iout++) {
        for(uint itype=0;itype<vtype.size();itype++) {   
          if(aurostd::FileExist(directory_RAW+"/"+vfile2compress1.at(ilink)+vtype.at(itype)+vout.at(iout))) {
            if(1||LDEBUG) { cout << soliloquy << " compressing file (" << ilink << "," << iout << "," << itype << "): " << string(directory_RAW+"/"+vfile2compress1.at(ilink)+vtype.at(itype)+vout.at(iout)) << endl; cout.flush(); }
            aurostd::CompressFile(directory_RAW+"/"+vfile2compress1.at(ilink)+vtype.at(itype)+vout.at(iout),DEFAULT_KZIP_BIN);
          } // FileExist
        } // itype
      } // iout
    } // ilink

    // FILES to link LIB to RAW
    for(uint ifile=0;ifile<vfile.size();ifile++) {
      if(aurostd::FileExist(directory_RAW+"/"+vfile.at(ifile))) { // need to be present
        deque<string> vfile2link0; aurostd::string2tokens("DOSCAR.static,OUTCAR.static,CHGCAR.static,AECCAR0.static,AECCAR2.static,EIGENVAL.bands,OSZICAR.bands",vfile2link0,",");
        for(uint ilink=0;ilink<vfile2link0.size();ilink++) {
          if(vfile.at(ifile)==vfile2link0.at(ilink)) {
            // cout << soliloquy << " linking file RAW->LIB: " << vfile2link0.at(ilink) << endl;       
            if(aurostd::FileExist(directory_LIB+"/"+vfile.at(ifile)) || aurostd::EFileExist(directory_LIB+"/"+vfile.at(ifile))) { // need to be present in LIB also  
              aurostd::RemoveFile(directory_RAW+"/"+vfile.at(ifile)); // remove RAW original
              if(aurostd::FileExist(directory_LIB+"/"+vfile.at(ifile))) { 
                cout << soliloquy << " linking file RAW->LIB: " << string(directory_LIB+"/"+vfile.at(ifile)) << endl; cout.flush();   
                aurostd::LinkFile(directory_LIB+"/"+vfile.at(ifile),directory_RAW); // link LIB to RAW (save space)
              }
              for(uint iext=0;iext<XHOST.vext.size();iext++) {
                if(aurostd::FileExist(directory_LIB+"/"+vfile.at(ifile)+XHOST.vext.at(iext))) {
                  cout << soliloquy << " linking file RAW->LIB: " << string(directory_LIB+"/"+vfile.at(ifile)+XHOST.vext.at(iext)) << endl;  cout.flush();          
                  aurostd::LinkFile(directory_LIB+"/"+vfile.at(ifile)+XHOST.vext.at(iext),directory_RAW); // link LIB to RAW (save space)
                }
              }
            }
          }
        } // files to link LIB to RAW
      } // File Exist
    } // ifile

    //[CO20200624 - OBSOLETE]// FILES to leave as is
    //[CO20200624 - OBSOLETE]for(uint ifile=0;ifile<vfile.size();ifile++) {
    //[CO20200624 - OBSOLETE]  if(aurostd::FileExist(directory_RAW+"/"+vfile.at(ifile))) { // need to be present
        // [NO COMPRESS] if(vfile.at(ifile)=="KPOINTS.relax") aurostd::execute(DEFAULT_KZIP_BIN+" -9f \""+directory_RAW+"/"+vfile.at(ifile)+"\""); // seems to work, although I do not like if(SC-0914)
        // [NO COMPRESS] if(vfile.at(ifile)=="KPOINTS.static") aurostd::execute(DEFAULT_KZIP_BIN+" -9f \""+directory_RAW+"/"+vfile.at(ifile)+"\""); // seems to work, although I do not like if(SC-0914)
        // [NO COMPRESS] if(vfile.at(ifile)=="KPOINTS.bands") aurostd::execute(DEFAULT_KZIP_BIN+" -9f \""+directory_RAW+"/"+vfile.at(ifile)+"\""); // seems to work, although I do not like if(SC-0914)
        // [NO COMPRESS] if(vfile.at(ifile)=="INCAR.relax") aurostd::execute(DEFAULT_KZIP_BIN+" -9f \""+directory_RAW+"/"+vfile.at(ifile)+"\""); // seems to work, although I do not like if(SC-0914)
        // [NO COMPRESS] if(vfile.at(ifile)=="INCAR.static") aurostd::execute(DEFAULT_KZIP_BIN+" -9f \""+directory_RAW+"/"+vfile.at(ifile)+"\""); // seems to work, although I do not like if(SC-0914)
        // [NO COMPRESS] if(vfile.at(ifile)=="INCAR.bands") aurostd::execute(DEFAULT_KZIP_BIN+" -9f \""+directory_RAW+"/"+vfile.at(ifile)+"\""); // seems to work, although I do not like if(SC-0914)
        // [NO COMPRESS] if(vfile.at(ifile)=="OUTCAR.relax") aurostd::execute(DEFAULT_KZIP_BIN+" -9f \""+directory_RAW+"/"+vfile.at(ifile)+"\"");  // seems to work, although I do not like if(SC-0512)
    //[CO20200624 - OBSOLETE]  } // File Exist
    //[CO20200624 - OBSOLETE]} // ifile

    // DO THE FINISH LINK/COPY FOR WEB
    cout << soliloquy << " linking stuff flag_WEB=" << flag_WEB << endl;
    // MOVE/LINK PICS data

    if(flag_WEB) {
      cout << soliloquy << " linking SYSTEM=" << system_name << endl;
      if(aurostd::FileExist(directory_RAW+"/"+system_name+".png") ||
          aurostd::FileExist(directory_RAW+"/"+system_name+"_banddos.png") ||  //ME20190621 - new file name convention
          FALSE){
        aurostd::LinkFile(directory_RAW+"/*png",directory_WEB);            // LINK
      }
      if(aurostd::FileExist(directory_RAW+"/"+system_name+".cif")) aurostd::LinkFile(directory_RAW+"/*cif",directory_WEB);            // LINK //CO20200624 - adding FileExist() check

      if(aurostd::FileExist("/www/AFLOWDATA/api_index.php")) aurostd::LinkFile("/www/AFLOWDATA/api_index.php",directory_RAW+"/index.php");                      // LINK //CO20200624 - adding FileExist() check
      if(aurostd::FileExist("/www/AFLOWDATA/api_index.php")) aurostd::LinkFile("/www/AFLOWDATA/api_index.php",directory_WEB+"/index.php");                      // LINK //CO20200624 - adding FileExist() check
      //CO20200624 - these are PRE-LINKS: these files have NOT been written yet
      //if we check if they exist first, then the links will not be created
      //we want to do this NOW so that it can be captured by vfiles
      //[CO20200624 - THIS IS A PRE-LINK!]if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) 
        aurostd::LinkFile(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,directory_WEB);    // LINK  //CO20200624 - adding FileExist() check
      //[CO20200624 - THIS IS A PRE-LINK!]if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON)) 
        aurostd::LinkFile(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON,directory_WEB);   // LINK //CO20200624 - adding FileExist() check
      //[CO20200624 - THIS IS A PRE-LINK!]if(aurostd::FileExist(directory_RAW+"/"+aflowlib_data.auid+".out")) 
        aurostd::LinkFile(directory_RAW+"/"+aflowlib_data.auid+".out",directory_WEB);     // LINK for AUID  //CO20200624 - adding FileExist() check
      //[CO20200624 - THIS IS A PRE-LINK!]if(aurostd::FileExist(directory_RAW+"/"+aflowlib_data.auid+".json")) 
        aurostd::LinkFile(directory_RAW+"/"+aflowlib_data.auid+".json",directory_WEB);     // LINK for AUID  //CO20200624 - adding FileExist() check

      deque<string> vfile2link1; aurostd::string2tokens("aflow.pgroup,aflow.pgroup_xtal,aflow.pgroupk,aflow.pgroupk_xtal,aflow.fgroup,aflow.iatoms,aflow.agroup,edata,data",vfile2link1,",");
      for(uint ilink=0;ilink<vfile2link1.size();ilink++) {
        for(uint iout=0;iout<vout.size();iout++) {
          for(uint itype=0;itype<vtype.size();itype++) {   
            if(aurostd::FileExist(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout))) { // no compression
              if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout)) << endl; cout.flush(); }              
              aurostd::LinkFile(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout),directory_WEB);
            }  // FileExist
            for(uint iext=0;iext<XHOST.vext.size();iext++) {
              if(aurostd::FileExist(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout)+XHOST.vext.at(iext))) { // with compression
                if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout)+XHOST.vext.at(iext)) << endl; cout.flush(); }
                aurostd::LinkFile(directory_RAW+"/"+vfile2link1.at(ilink)+vtype.at(itype)+vout.at(iout)+XHOST.vext.at(iext),directory_WEB);
              } // FileExist
            } // iext
          } //itype
        } // iout
      } // ilink
      // PREVIOUSLY IN NETFLIX // cheat for void string
      for(uint iout=0;iout<vout.size();iout++) {
        if(aurostd::FileExist(directory_RAW+"/"+system_name+"_structure_relax"+vout.at(iout))) {
          if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_structure_relax"+vout.at(iout)) << endl; cout.flush(); }
          aurostd::LinkFile(directory_RAW+"/"+system_name+"_structure_relax"+vout.at(iout),directory_WEB);  //CO20171024
        }  // FileExist
      } // iout
      for(uint iout=0;iout<vout.size();iout++) {
        if(aurostd::FileExist(directory_RAW+"/"+system_name+"_structure_relax1"+vout.at(iout))) {
          if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_structure_relax1"+vout.at(iout)) << endl; cout.flush(); }
          aurostd::LinkFile(directory_RAW+"/"+system_name+"_structure_relax1"+vout.at(iout),directory_WEB);  //CO20171024
        } // FileExist
      } // iout
      for(uint iout=0;iout<vout.size();iout++) {
        if(aurostd::FileExist(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout)))  {  // NO EXTENSION
          if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout)) << endl; cout.flush(); }
          aurostd::LinkFile(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout),directory_WEB);       //CO20171024
        }
        for(uint iext=0;iext<XHOST.vext.size();iext++) {
          if(aurostd::FileExist(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout)+XHOST.vext.at(iext)))  {
            if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout)+XHOST.vext.at(iext)) << endl; cout.flush(); }
            aurostd::LinkFile(directory_RAW+"/"+system_name+"_dosdata"+vout.at(iout)+XHOST.vext.at(iext),directory_WEB);       //CO20171024
          } // FileExist
        } // iext
      } // iout
      for(uint iout=0;iout<vout.size();iout++) {
        if(aurostd::FileExist(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout))) { // NO EXTENSION
          if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout)) << endl; cout.flush(); }
          aurostd::LinkFile(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout),directory_WEB);     //CO20171024
        } // FileExist
        for(uint iext=0;iext<XHOST.vext.size();iext++) {
          if(aurostd::FileExist(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout)+XHOST.vext.at(iext))) {
            if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout)+XHOST.vext.at(iext)) << endl; cout.flush(); }
            aurostd::LinkFile(directory_RAW+"/"+system_name+"_bandsdata"+vout.at(iout)+XHOST.vext.at(iext),directory_WEB);     //CO20171024
          } // FileExist
        } // iext
      } // iout
      deque<string> vfile2link2; aurostd::string2tokens("EIGENVAL.bands,DOSCAR.static,OUTCAR.static,CONTCAR.relax,CONTCAR.relax1,POSCAR.bands,CONTCAR.relax.vasp,CONTCAR.relax.qe,CONTCAR.relax.abinit,CONTCAR.relax.aims,KPOINTS.relax,KPOINTS.static,KPOINTS.bands,INCAR.bands,AECCAR0.static,AECCAR2.static,CHGCAR.static",vfile2link2,",");
      for(uint ilink=0;ilink<vfile2link2.size();ilink++) {
        if(aurostd::FileExist(directory_RAW+"/"+vfile2link2.at(ilink))) { // NO EXTENSION
          if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+vfile2link2.at(ilink)) << endl; cout.flush(); }
          aurostd::LinkFile(directory_RAW+"/"+vfile2link2.at(ilink),directory_WEB);         // LINK
        }
        for(uint iext=0;iext<XHOST.vext.size();iext++) {
          if(aurostd::FileExist(directory_RAW+"/"+vfile2link2.at(ilink)+XHOST.vext.at(iext))) {
            if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+vfile2link2.at(ilink)+XHOST.vext.at(iext)) << endl; cout.flush(); }
            aurostd::LinkFile(directory_RAW+"/"+vfile2link2.at(ilink)+XHOST.vext.at(iext),directory_WEB);         // LINK
          }
        } // iext
      } // ilink

      if(perform_BADER) {
        //CO20200624 - we don't need to worry about zipped variants, this analysis is done in RAW only (no LIB compression)
        //H15Th4_ICSD_638495_Bader_50_H.jvxl - 50 will always exist
        if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+"*jvxl") << endl; cout.flush(); }
        if(aflowlib_data.vspecies.size()>0 && aurostd::FileExist(directory_RAW+"/"+system_name+"_Bader_50_"+aflowlib_data.vspecies[0]+".jvxl")){  //CO20200624 - adding FileExist() check
          aurostd::LinkFile(directory_RAW+"/"+"*jvxl*",directory_WEB);         // LINK
        }
        //H15Th4_ICSD_638495_abader.out
        if(LDEBUG) { cout << soliloquy << " linking file WEB->RAW: " << string(directory_RAW+"/"+system_name+"_abader.out") << endl; cout.flush(); }  //CO20200624  //CO20200624 - adding FileExist() check
        if(aurostd::FileExist(directory_RAW+"/"+system_name+"_abader.out")) aurostd::LinkFile(directory_RAW+"/"+system_name+"_abader.out",directory_WEB); // LINK //CO20200624
      }
    } // flag_WEB
    
    // DONE
    // write files if necessary
    vector<string> vdirectory;
    // do the directories
    if(flag_files_LIB) {
      aurostd::DirectoryLS(directory_LIB,vdirectory);
      for(uint i=0;i<vdirectory.size();i++)
        aflowlib_data.vfiles.push_back(vdirectory.at(i));
    }
    if(flag_files_RAW) {
      aurostd::DirectoryLS(directory_RAW,vdirectory);
      for(uint i=0;i<vdirectory.size();i++)
        aflowlib_data.vfiles.push_back(vdirectory.at(i));
    }
    if(flag_files_WEB) {
      aurostd::DirectoryLS(directory_WEB,vdirectory);
      for(uint i=0;i<vdirectory.size();i++)
        aflowlib_data.vfiles.push_back(vdirectory.at(i));
    }
    // DO THE FINAL WRITING
    
    if(aflowlib_data.vaflowlib_date.size()!=1){ //CO20200624 - this means we didn't get the LOCK date, spit warning
      pflow::logger(_AFLOW_FILE_NAME_,soliloquy,"LOCK date NOT found",_LOGGER_WARNING_);
    }
    aflowlib_data.vaflowlib_date.push_back(aurostd::get_datetime()+"_GMT-5"); //CO20200624 - adding LOCK date

    //     cout << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << ": " << aflowlib_data.aflowlib2file(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT);
    //      aurostd::LinkFile("../../"+_XENTRY_","directory_RAW+"/"+_XENTRY_);
    if(!LOCAL) { //CO20171025
      if(aurostd::FileExist("/www/AFLOWDATA/api_index.php")) aurostd::LinkFile("/www/AFLOWDATA/api_index.php",directory_RAW+"/"+_XENTRY_);  //CO20200624 - adding FileExist() check
    }

    // write aflowlib.out
    cout << soliloquy << " writing " << string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT) << endl; cout.flush(); 
    cout << XPID << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << ": " << aflowlib_data.aflowlib2file(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,"out");
    cout << soliloquy << " linking file RAW->RAW: " << string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT) << " -> " << string(directory_RAW+"/"+aflowlib_data.auid+".out") << endl; cout.flush(); 
    if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT)) aurostd::LinkFile(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT,directory_RAW+"/"+aflowlib_data.auid+".out");         // LINK //CO20200624 - adding FileExist() check
    // cout << DEFAULT_FILE_AFLOWLIB_ENTRY_OUT << ": " << aflowlib_data.aflowlib2file(directory_WEB+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_OUT);

    // write aflowlib.json
    cout << soliloquy << " writing " << string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON) << endl; cout.flush(); 
    cout << XPID << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << ": " << aflowlib_data.aflowlib2file(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON,"json");
    cout << soliloquy << " linking file RAW->RAW: " << string(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON) << " -> " << string(directory_RAW+"/"+aflowlib_data.auid+".json") << endl; cout.flush(); 
    if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON)) aurostd::LinkFile(directory_RAW+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON,directory_RAW+"/"+aflowlib_data.auid+".json");         // LINK  //CO20200624 - adding FileExist() check
    // cout << DEFAULT_FILE_AFLOWLIB_ENTRY_JSON << ": " << aflowlib_data.aflowlib2file(directory_WEB+"/"+DEFAULT_FILE_AFLOWLIB_ENTRY_JSON);

    if(flag_WEB) {
      if(!LOCAL) { //CO20171025
        if(aurostd::FileExist("/www/AFLOWDATA/api_index.php")) aurostd::LinkFile("/www/AFLOWDATA/api_index.php",directory_WEB+"/"+_XENTRY_);  //CO20200624 - adding FileExist() check
      }
    }

    // CHECK FOR AUID-WEB-LINKS
    //      LDEBUG=TRUE;
    if (!LOCAL) {
      if(LDEBUG) cout << soliloquy << " flag_WEB=" << flag_WEB << endl;
      if(LDEBUG) cout << soliloquy << " directory_LIB=" << directory_LIB << endl;
      if(LDEBUG) cout << soliloquy << " directory_RAW=" << directory_RAW << endl;
      if(LDEBUG) cout << soliloquy << " directory_WEB=" << directory_WEB << endl;
      //if(LDEBUG)
      cout << soliloquy << " aflowlib_data.auid=" << aflowlib_data.auid << endl;
      // **----------------------------------------------------------------------------
      if(XHOST.hostname=="nietzsche.mems.duke.edu" && (XHOST.user=="auro"||XHOST.user=="common"||XHOST.user=="stefano")){  //CO20200624 - cannot create these links otherwise
        // **----------------------------------------------------------------------------
        // NEW BUT STILL DOING
        string directory_AUID="";
        string directory_AUID_LIB="",directory_AUID_RAW="",directory_AUID_WEB="";
        directory_AUID=init::AFLOW_Projects_Directories("AUID")+"/"+aflowlib::auid2directory(aflowlib_data.auid);
        aurostd::DirectoryMake(directory_AUID);
        directory_AUID_LIB=directory_AUID+"/LIB"; 
        directory_AUID_RAW=directory_AUID+"/RAW";
        directory_AUID_WEB=directory_AUID+"/WEB";
        aurostd::RemoveFile(directory_AUID_LIB); // to avoid auto-linking SC20181205
        aurostd::RemoveFile(directory_AUID_RAW); // to avoid auto-linking SC20181205
        aurostd::RemoveFile(directory_AUID_WEB); // to avoid auto-linking SC20181205

        // directory_AUID_LIB
        if(LDEBUG) cout << soliloquy << " (AUID_NEW) directory_AUID_LIB=" << directory_AUID_LIB << " -> " << directory_LIB << endl;
        cout << soliloquy << " (AUID_NEW) linking file AUID_LIB->LIB: " << directory_AUID_LIB << " -> " << directory_LIB << endl; cout.flush();
        if(aurostd::IsDirectory(directory_LIB)) aurostd::LinkFile(directory_LIB,directory_AUID_LIB);         // LINK  //CO20200624 - adding FileExist() check
        // directory_AUID_RAW
        if(LDEBUG) cout << soliloquy << " (AUID_NEW) directory_AUID_RAW=" << directory_AUID_RAW << " -> " << directory_RAW << endl;
        cout << soliloquy << " (AUID_NEW) linking file AUID_RAW->LIB: " << directory_AUID_RAW << " -> " << directory_RAW << endl; cout.flush();
        if(aurostd::IsDirectory(directory_RAW)) aurostd::LinkFile(directory_RAW,directory_AUID_RAW);         // LINK  //CO20200624 - adding FileExist() check

        if(flag_WEB) {
          if(LDEBUG) cout << soliloquy << " (AUID_NEW) directory_AUID_WEB=" << directory_AUID_WEB << " -> " << directory_WEB << endl;
          cout << soliloquy << " (AUID_NEW) linking file AUID_WEB->LIB: " << directory_AUID_WEB << " -> " << directory_WEB << endl; cout.flush();
          if(aurostd::IsDirectory(directory_WEB)) aurostd::LinkFile(directory_WEB,directory_AUID_WEB);         // LINK  //CO20200624 - adding FileExist() check
        } else {
          if(LDEBUG) cout << soliloquy << " (AUID_NEW) directory_AUID_WEB=" << directory_AUID_WEB << " -> " << directory_RAW << endl;
          cout << soliloquy << " (AUID_NEW) linking file AUID_WEB->LIB: " << directory_AUID_WEB << " -> " << directory_RAW << endl; cout.flush();
          if(aurostd::IsDirectory(directory_RAW)) aurostd::LinkFile(directory_RAW,directory_AUID_WEB);         // LINK  //CO20200624 - adding FileExist() check
        }

        // ICSD2LINK
        if(aflowlib_data.catalog=="ICSD") {
          aurostd::string2tokens(directory_LIB,tokens,"_");
          if(tokens.size()>2) {
            if(tokens.at(tokens.size()-2)=="ICSD") {
              string directory_ICSD2LINK=init::AFLOW_Projects_Directories("AUID")+"/icsd:/"+tokens.at(tokens.size()-1);
              aurostd::DirectoryMake(directory_ICSD2LINK);	      
              cout << soliloquy << " (ICSD2LINK) making ICSD2LINK: " << directory_ICSD2LINK << endl; cout.flush();
              // LIB
              aurostd::RemoveFile(directory_ICSD2LINK+"/LIB");  // to avoid auto-linking SC20190830
              if(aurostd::IsDirectory(directory_LIB)) aurostd::LinkFile(directory_LIB,directory_ICSD2LINK+"/LIB"); // LINK  //CO20200624 - adding FileExist() check
              cout << soliloquy << " (ICSD2LINK) linking file LIB->ICSD2LINK/LIB: " << directory_LIB << " -> " << directory_ICSD2LINK << "/LIB" << endl; cout.flush();
              // RAW
              aurostd::RemoveFile(directory_ICSD2LINK+"/RAW");  // to avoid auto-linking SC20190830
              if(aurostd::IsDirectory(directory_RAW)) aurostd::LinkFile(directory_RAW,directory_ICSD2LINK+"/RAW"); // LINK  //CO20200624 - adding FileExist() check
              cout << soliloquy << " (ICSD2LINK) linking file RAW->ICSD2LINK/RAW: " << directory_RAW << " -> " << directory_ICSD2LINK << "/RAW" << endl; cout.flush();
              // WEB
              aurostd::RemoveFile(directory_ICSD2LINK+"/WEB");  // to avoid auto-linking SC20190830
              if(aurostd::IsDirectory(directory_WEB)) aurostd::LinkFile(directory_WEB,directory_ICSD2LINK+"/WEB"); // LINK  //CO20200624 - adding FileExist() check
              cout << soliloquy << " (ICSD2LINK) linking file WEB->ICSD2LINK/WEB: " << directory_WEB << " -> " << directory_ICSD2LINK << "/WEB" << endl; cout.flush();
            }
          }
        }
      }


      // for(uint i=0;i<aflowlib_data.vauid.size();i++) { cout << soliloquy << " DEBUG  aflowlib_data.vauid.at(" << i << ")=" << aflowlib_data.vauid.at(i) << endl; }
      // DONE
      //      cout << soliloquy << " dir=" << directory_LIB << "   END_DATE - [v=" << string(AFLOW_VERSION) << "] -" << " [time=" << aurostd::utype2string(aurostd::get_seconds(seconds_begin),2) << "] " << Message(aflags,"time",_AFLOW_FILE_NAME_) << endl;
      cout << soliloquy << " dir=" << directory_LIB << "   END_DATE - [v=" << string(AFLOW_VERSION) << "] -" << Message(aflags,"time",_AFLOW_FILE_NAME_) << " [time=" << aurostd::utype2string(aurostd::get_seconds(seconds_begin),2,FIXED_STREAM) << "] " << endl; //CO20200624 - added FIXED_STREAM
      if(XHOST.vflag_control.flag("BEEP")) aurostd::beep(aurostd::min(6000,aurostd::abs(int(1*aflowlib_data.aflowlib2string().length()-2000))),50);
    }
    
    // COMPRESS
    /*
    //  aurostd::execute(DEFAULT_KZIP_BIN+" -9f \""+directory_RAW+"/POSCAR.bands\"");
    //  aurostd::execute(DEFAULT_KZIP_BIN+" -9f \""+directory_RAW+"/plotbz.sh\"");
    */
    // DELETE STUFF

    // CHANGE PERMISSIONS
    // changing order of permission editing, if LOCAL, then order matters, otherwise NOT REALLY
    // files first, since we do /* (just in case directory_RAW is inside directory_LIB)
    // then directories
    // FILES

    //[CO20200624 - leave as what it is]LDEBUG=FALSE; //CO20180321
    //CO20180216 - more robust for any type of directory/file setup
    vector<string> Chmod_Files;
    if(LDEBUG) { cerr << soliloquy << " pwd=" << aurostd::getPWD() << endl;}  //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")

    //[CO20190321 - bust if find grabs nothing]aurostd::execute("chmod 755 `find \""+directory_RAW+"\" -type d`");
    aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find")+" \""+directory_RAW+"\" -type d"),Chmod_Files);
    for(uint i=0;i<Chmod_Files.size();i++) {
      if(LDEBUG) { cerr << soliloquy << " \"" << XHOST.command("chmod")+" 755 \""+Chmod_Files[i] << "\"\"" << endl;}
      aurostd::ChmodFile("755",Chmod_Files[i]);
    }

    //[CO20190321 - bust if find grabs nothing]aurostd::execute("chmod 644 `find \""+directory_RAW+"\" -type f`");
    aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find")+" \""+directory_RAW+"\" -type f"),Chmod_Files);
    for(uint i=0;i<Chmod_Files.size();i++) {
      if(LDEBUG) { cerr << soliloquy << " \"" << XHOST.command("chmod")+" 644 \""+Chmod_Files[i] << "\"\"" << endl;}
      aurostd::ChmodFile("644",Chmod_Files[i]);
    }

    if(CHMODWEB) {
      //[CO20190321 - bust if find grabs nothing]aurostd::execute("chmod 755 `find \""+directory_WEB+"\" -type d`");
      aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find")+" \""+directory_WEB+"\" -type d"),Chmod_Files);
      for(uint i=0;i<Chmod_Files.size();i++) {
        if(LDEBUG) { cerr << soliloquy << " \"" << XHOST.command("chmod")+" 755 \""+Chmod_Files[i] << "\"\"" << endl;}
        aurostd::ChmodFile("755",Chmod_Files[i]);
      }

      //[CO20190321 - bust if find grabs nothing]aurostd::execute("chmod 644 `find \""+directory_WEB+"\" -type f`");
      aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find")+" \""+directory_WEB+"\" -type f"),Chmod_Files);
      for(uint i=0;i<Chmod_Files.size();i++) {
        if(LDEBUG) { cerr << soliloquy << " \"" << XHOST.command("chmod")+" 644 \""+Chmod_Files[i] << "\"\"" << endl;}
        aurostd::ChmodFile("644",Chmod_Files[i]);
      }
    }

    //[CO20190321 - bust if find grabs nothing]aurostd::execute("chmod 755 `find \""+directory_LIB+"\" -type d`");
    aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find")+" \""+directory_LIB+"\" -type d"),Chmod_Files);
    for(uint i=0;i<Chmod_Files.size();i++) {
      if(LDEBUG) { cerr << soliloquy << " \"" << XHOST.command("chmod")+" 755 \""+Chmod_Files[i] << "\"\"" << endl;}
      aurostd::ChmodFile("755",Chmod_Files[i]);
    }

    //[CO20190321 - bust if find grabs nothing]aurostd::execute("chmod 644 `find \""+directory_LIB+"\" -type f`");
    aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find")+" \""+directory_LIB+"\" -type f"),Chmod_Files);
    for(uint i=0;i<Chmod_Files.size();i++) {
      if(LDEBUG) { cerr << soliloquy << " \"" << XHOST.command("chmod")+" 644 \""+Chmod_Files[i] << "\"\"" << endl;}
      aurostd::ChmodFile("644",Chmod_Files[i]);
    }

    // ALTERNATIVE but bombs xchmod {}
    // aurostd::execute("find \""+directory_LIB+"\" -type d ! -perm 755 -print -exec xchmod 755 {} \\;");
    // aurostd::execute("find \""+directory_LIB+\"" -type f ! -perm 644 -print -exec xchmod 644 {} \\;");
    // aurostd::execute("find \""+directory_RAW+"\" -type d ! -perm 755 -print -exec xchmod 755 {} \\;");
    // aurostd::execute("find \""+directory_RAW+"\" -type f ! -perm 644 -print -exec xchmod 644 {} \\;");
    // if(CHMODWEB) aurostd::execute("find \""+directory_WEB+"\" -type d ! -perm 755 -print -exec xchmod 755 {} \\;");
    // if(CHMODWEB) aurostd::execute("find \""+directory_WEB+"\" -type f ! -perm 644 -print -exec xchmod 644 {} \\;");

    // [OBSOLETE CO20180216] aurostd::ChmodFile("755",directory_LIB);
    // [OBSOLETE CO20180216] aurostd::ChmodFile("644",directory_LIB+"/*");
    // [OBSOLETE CO20180216] aurostd::ChmodFile("755",directory_LIB+"/ARUN*");
    // [OBSOLETE CO20180216] aurostd::ChmodFile("755",directory_RAW);
    // [OBSOLETE CO20180216] aurostd::ChmodFile("644",directory_RAW+"/*");
    // [OBSOLETE CO20180216] aurostd::ChmodFile("755",directory_RAW+"/ARUN*");
    // [OBSOLETE CO20180216] if(CHMODWEB) if(flag_WEB) { aurostd::ChmodFile("755",directory_WEB); }
    // [OBSOLETE CO20180216] if(CHMODWEB) if(flag_WEB) { aurostd::ChmodFile("644",directory_WEB+"/*.png"); }
    // [OBSOLETE CO20180216] if(CHMODWEB) if(flag_WEB) { aurostd::ChmodFile("644",directory_WEB+"/*.html"); }
    // [OBSOLETE CO20180216] if(CHMODWEB) if(flag_WEB) { aurostd::ChmodFile("644",directory_WEB+"/*.jpg"); }
    // [OBSOLETE CO20180216] if(CHMODWEB) if(flag_WEB) { aurostd::ChmodFile("755",directory_WEB+"/ARUN*"); }
    // NOW WRITE DOWN THE FILE FOR THE LIBRARY

    // done
    //    cout << soliloquy << " [1]" << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_FileNeeded
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_FileNeeded(string directory_LIB,string fileLIB,string directory_RAW,string fileRAW,vector<string> &vfile,const string& MESSAGE) {
    //  bool LDEBUG=(FALSE || XHOST.DEBUG);

    if(XHOST.vext.size()!=XHOST.vzip.size()) { cerr << XPID << "ERROR - aflowlib::LIB2RAW_FileNeeded: XHOST.vext.size()!=XHOST.vzip.size()" << endl;exit(0); }

    string file_LIB=directory_LIB+"/"+fileLIB;
    string file_RAW=directory_RAW+"/"+fileRAW;
    string file_LIB_nocompress=directory_LIB+"/"+fileLIB;
    for(uint iext=1;iext<XHOST.vext.size();iext++) aurostd::StringSubst(file_LIB_nocompress,XHOST.vext.at(iext),"");  // SKIP uncompressed
    string file_RAW_nocompress=directory_RAW+"/"+fileRAW;
    for(uint iext=1;iext<XHOST.vext.size();iext++) aurostd::StringSubst(file_RAW_nocompress,XHOST.vext.at(iext),"");  // SKIP uncompressed
    if(aurostd::FileExist(file_RAW)) return TRUE;   // already there
    if(aurostd::FileExist(file_RAW_nocompress)) return TRUE; // already there

    if(!aurostd::FileExist(file_LIB) && !aurostd::FileExist(file_LIB_nocompress) && !aurostd::EFileExist(file_LIB_nocompress)) {
      if(!aurostd::FileExist(file_LIB)) {
        cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_FileNeeded: file not found " << file_LIB << endl;exit(0); }
      if(!aurostd::FileExist(file_LIB_nocompress)) {
        cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_FileNeeded: file not found " << file_LIB_nocompress << endl;exit(0); }
      if(!aurostd::EFileExist(file_LIB_nocompress)) {
        cout << MESSAGE << " ERROR - aflowlib::LIB2RAW_FileNeeded: file not found " << file_LIB_nocompress << ".EXT" << endl;exit(0); }
    }
    if(aurostd::FileExist(file_LIB)) aurostd::CopyFile(file_LIB,file_RAW);
    if(aurostd::FileExist(file_LIB_nocompress)) aurostd::CopyFile(file_LIB_nocompress,file_RAW_nocompress);
    for(uint iext=1;iext<XHOST.vext.size();iext++)  { // SKIP uncompressed
      if(aurostd::FileExist(file_LIB_nocompress+XHOST.vext.at(iext))) {
        aurostd::CopyFile(file_LIB_nocompress+XHOST.vext.at(iext),file_RAW_nocompress+XHOST.vext.at(iext));
      }
    }
    for(uint iext=1;iext<XHOST.vext.size();iext++) {   // SKIP uncompressed
      if(aurostd::FileExist(file_RAW) && aurostd::substring2bool(file_RAW,XHOST.vext.at(iext))) aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+file_RAW+"\"");
      if(aurostd::FileExist(file_RAW_nocompress+XHOST.vext.at(iext))) aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+file_RAW_nocompress+XHOST.vext.at(iext)+"\"");
    }
    string file2add=fileRAW;
    for(uint iext=1;iext<XHOST.vext.size();iext++) {   // SKIP uncompressed
      aurostd::StringSubst(file2add,XHOST.vext.at(iext),""); 
    }
    vfile.push_back(file2add);
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Bands
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_Bands(const string& directory_LIB,const string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,const string& MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // LDEBUG=TRUE;
    if(LDEBUG) cerr << XPID << "aflowlib::LIB2RAW_Loop_Bands [1]" << endl;
    // Stefano Curtarolo 2009-2010-2011-2012
    //[CO20200624 - OBSOLETE]vector<string> vspecies;aurostd::string2tokens(data.species,vspecies,",");

    cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Bands: begin " << directory_LIB << endl;
    cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Bands: species = " << data.vspecies.size() << endl;
    data.vloop.push_back("bands");

    stringstream command;command.clear();command.str(std::string());
    stringstream aus_exec;
    // directories must exist already
    bool flag_DATA_BANDS_=FALSE;
    bool flag_use_MATLAB=FALSE,flag_use_GNUPLOT=!flag_use_MATLAB;  //KY
    // bool flag_use_MATLAB=TRUE,flag_use_GNUPLOT=!flag_use_MATLAB;   //WSETYAWAN
    // [OBSOLETE]  bool flag_ORIG=FALSE;

    //[CO20200624 - do inside CleanDirectoryLIB(), PROBABLY OBSOLETE]aurostd::StringSubst(directory_LIB,"ELPASOLITES","AURO"); // PATCH
    //[CO20200624 - do inside CleanDirectoryLIB(), PROBABLY OBSOLETE]aurostd::StringSubst(directory_LIB,"SCINT","ICSD"); // PATCH

    // copy _AFLOWIN_ LOCK DOSCAR.static.EXT EIGENVAL.bands.EXT KPOINTS.bands.EXT POSCAR.bands.EXT
    string file_LIB,file_RAW;
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWIN_,directory_RAW,_AFLOWIN_,vfile,MESSAGE);  // _AFLOWIN_
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWLOCK_,directory_RAW,_AFLOWLOCK_,vfile,MESSAGE);  // LOCK
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"DOSCAR.static",directory_RAW,"DOSCAR.static",vfile,MESSAGE);  // DOSCAR.static
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.static",directory_RAW,"OUTCAR.static",vfile,MESSAGE);  // OUTCAR.static
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OSZICAR.static",directory_RAW,"OSZICAR.static",vfile,MESSAGE);  // OSZICAR.static
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OSZICAR.bands",directory_RAW,"OSZICAR.bands",vfile,MESSAGE);  // OSZICAR.bands
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"EIGENVAL.bands",directory_RAW,"EIGENVAL.bands",vfile,MESSAGE);  // EIGENVAL.bands
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"KPOINTS.static",directory_RAW,"KPOINTS.static",vfile,MESSAGE);  // KPOINTS.static  // not needed but good for show off SC 0914
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"KPOINTS.bands",directory_RAW,"KPOINTS.bands",vfile,MESSAGE);  // KPOINTS.bands
    //  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"INCAR.static",directory_RAW,"INCAR.static",vfile,MESSAGE);  // INCAR.static  // not needed but good for show off SC 0914
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,"INCAR.bands",directory_RAW,"INCAR.bands",vfile,MESSAGE);  // INCAR.bands  // not needed but good for show off SC 0914
    //  aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.relax1",directory_RAW,"POSCAR.relax1",vfile,MESSAGE);  // POSCAR.relax1, Get Atomic Species Name


    if(TRUE || flag_DATA_BANDS_) { // POSCAR.bands.EXT
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.bands",directory_RAW,"POSCAR.bands",vfile,MESSAGE);  // POSCAR.bands
    }

    xKPOINTS kpoints_static;
    kpoints_static.GetPropertiesFile(directory_RAW+"/KPOINTS.static");
    data.kpoints_nnn_static=kpoints_static.nnn_kpoints;
    data.kpoints+=";"+aurostd::utype2string(kpoints_static.nnn_kpoints[1])+","+aurostd::utype2string(kpoints_static.nnn_kpoints[2])+","+aurostd::utype2string(kpoints_static.nnn_kpoints[3]);
    xKPOINTS kpoints_bands;
    kpoints_bands.GetPropertiesFile(directory_RAW+"/KPOINTS.bands");
    //get pairs
    data.kpoints_pairs.clear();
    if(kpoints_bands.vpath.size()%2==0) {  //if even
      for(uint i=0;i<kpoints_bands.vpath.size();i+=2) {
        data.kpoints_pairs.push_back(kpoints_bands.vpath.at(i)+"-"+kpoints_bands.vpath.at(i+1));
      }
    }
    data.kpoints_bands_path_grid=kpoints_bands.path_grid;
    data.kpoints+=";"+kpoints_bands.path+";"+aurostd::utype2string(kpoints_bands.path_grid);
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " KPOINTS = " << data.kpoints << endl;

    if(flag_use_MATLAB) { // MATLAB STUFF  OLD WSETYAWAN+SC
      // PERFORM THE MATLAB STEP
      cout << MESSAGE << " MATLAB start: " << directory_RAW << endl;
      //WRITE plotbz.sh
      stringstream gnuplot_plotbz;
      gnuplot_plotbz.clear();gnuplot_plotbz.str(std::string());
      // [OBSOLETE]    gnuplot_plotbz << GNUPLOT_FUNCS_plotbz("./","") << endl;
      gnuplot_plotbz << GNUPLOT_FUNCS_plotbz() << endl;
      aurostd::stringstream2file(gnuplot_plotbz,string(directory_RAW+"/plotbz.sh"));

      // WRITE PARAM.M
      stringstream matlab_param;
      matlab_param.clear();matlab_param.str(std::string());
      matlab_param << MATLAB_FUNCS_param() << endl;
      // aurostd::stringstream2file(matlab_param,string(directory_RAW+"/param.m")); // it seems that param is not needed anymore
      // WRITE PLOTBAND.M
      stringstream matlab_plotband;
      matlab_plotband.clear();matlab_plotband.str(std::string());
      // [OBSOLETE]   matlab_plotband << MATLAB_FUNCS_plotband("./","normal") << endl; // normal OR log
      matlab_plotband << MATLAB_FUNCS_plotband() << endl; // normal OR log
      matlab_plotband << "exit;"   << endl;
      aurostd::stringstream2file(matlab_plotband,string(directory_RAW+"/plotband.m"));

      // NEW STUFF
      command.clear();command.str(std::string());
      command << "cd " << directory_RAW << endl;
      command << "mv KPOINTS.bands KPOINTS.bands.old" << endl; vfile.push_back("KPOINTS.bands.old"); // so it is compressed
      command << "mv EIGENVAL.bands EIGENVAL.bands.old" << endl; vfile.push_back("EIGENVAL.bands.old"); // so it is compressed
      aurostd::execute(command);
      command.clear();command.str(std::string());
      _aflags aflags;
      aflags.Directory=directory_RAW;

      if(pflow::FIXBANDS(aflags,"POSCAR.bands,KPOINTS.bands.old,EIGENVAL.bands.old,KPOINTS.bands,EIGENVAL.bands")==FALSE) {
        cout << "ERROR_RERUN " << directory_LIB << endl;
        return FALSE;
      }

      // EXECUTE PLOTBZ.SH using ksh
      command.clear();command.str(std::string());
      command << "cd \"" << directory_RAW << "\"" << endl;
      command << "ksh plotbz.sh" << endl;
      aurostd::execute(command);

      // EXECUTE MATLAB
      command.clear();command.str(std::string());
      command << "cd \"" << directory_RAW << "\"" << endl;
      command << "export DISPLAY=:0.0" << endl;
      aurostd::CommandRequired(DEFAULT_KBIN_MATLAB_BIN); // MATLAB MUST BE AVAILABLE
      command << DEFAULT_KBIN_MATLAB_BIN << " -r " << string("plotband") << endl;
      aurostd::execute(command);
    }

    if(flag_use_GNUPLOT) { // GNUPLOT STUFF NEW KY+SC
      //ME20190614 BEGIN
      // This has to come first because FIXBANDS messes up the EIGENVAL files
      aurostd::xoption opts, plotoptions;
      string dosscale = aurostd::utype2string<double>(DEFAULT_DOS_SCALE);
      opts.push_attached("PLOT_DOS", directory_RAW + ",,," + dosscale);
      opts.push_attached("PLOT_BANDDOS", directory_RAW + ",,," + dosscale);
      opts.push_attached("PLOT_PDOS", directory_RAW + ",-1,,," + dosscale);
      opts.push_attached("PLOTTER::PRINT", "png");
      plotoptions = plotter::getPlotOptionsEStructure(opts, "PLOT_DOS");
      plotter::PLOT_DOS(plotoptions);
      plotoptions = plotter::getPlotOptionsEStructure(opts, "PLOT_BANDDOS");
      plotter::PLOT_BANDDOS(plotoptions);
      plotoptions = plotter::getPlotOptionsEStructure(opts, "PLOT_PDOS", true);
      plotter::PLOT_PDOS(plotoptions);
      //ME20190614 END

      //KY WRITE THE CODE HERE
      cout << MESSAGE << " GNUPLOT start: " << directory_RAW << endl;
      // WRITE plotbz.sh
      stringstream gnuplot_plotbz;
      gnuplot_plotbz.clear();gnuplot_plotbz.str(std::string());
      // [OBSOLETE]    gnuplot_plotbz << GNUPLOT_FUNCS_plotbz("./","") << endl;
      gnuplot_plotbz << GNUPLOT_FUNCS_plotbz() << endl;
      aurostd::stringstream2file(gnuplot_plotbz,string(directory_RAW+"/plotbz.sh"));

      // NEW STUFF
      command.clear();command.str(std::string());
      command << "cd \"" << directory_RAW << "\"" << endl;
      command << "mv KPOINTS.bands KPOINTS.bands.old" << endl; vfile.push_back("KPOINTS.bands.old"); // so it is compressed
      command << "mv EIGENVAL.bands EIGENVAL.bands.old" << endl; vfile.push_back("EIGENVAL.bands.old"); // so it is compressed
      aurostd::execute(command);

      command.clear();command.str(std::string());
      _aflags aflags;
      aflags.Directory=directory_RAW;
      if(pflow::FIXBANDS(aflags,"POSCAR.bands,KPOINTS.bands.old,EIGENVAL.bands.old,KPOINTS.bands,EIGENVAL.bands")==FALSE) {
        cout << "ERROR_RERUN " << directory_LIB << endl;
        return FALSE;
      }

      // EXECUTE PLOTBZ.SH using ksh
      command.clear();command.str(std::string());
      command << "cd \"" << directory_RAW << "\"" << endl;
      command << "ksh plotbz.sh" << endl;
      aurostd::execute(command);

      // EXECUTE PLOTBAND
      //[CO20191112 - OBSOLETE]char work_dir[1024];
      //[CO20191112 - OBSOLETE]string cdir; //, wdir;
      //[CO20191112 - OBSOLETE]cdir = getcwd(work_dir, 1024);  //Get the working directory
      //[CO20191112 - OBSOLETE back and forth directory change]string work_dir=aurostd::getPWD();    //Get the working directory //CO20191112

      //[CO20191112 - OBSOLETE]char raw_dir[1024];
      //[CO20191112 - OBSOLETE]strcpy(raw_dir, directory_RAW.c_str());
      //[CO20191112 - OBSOLETE]chdir(raw_dir);               //Change into the RAW directory
      //[CO20191112 - OBSOLETE back and forth directory change]chdir(directory_RAW.c_str());                //Change into the RAW directory  //CO20191112

      // [OBSOLETE]  vector<string> directory;
      // [OBSOLETE]  directory.push_back(" ");
      // [OBSOLETE]  directory.push_back(" ");
      // [OBSOLETE]  directory.push_back("./");
      // [OBSOLETE]  estructure::PLOT_BANDDOS(directory);
      // [OBSOLETE]  estructure::PLOT_PEDOSALL_AFLOWLIB(directory, aflags);

      // [OBSOLETE - ME20190614]  estructure::PLOT_BANDDOS("./");
      // [OBSOLETE - ME20190614]  estructure::PLOT_PEDOSALL_AFLOWLIB("./", aflags);

      //[CO20191112 - OBSOLETE back and forth directory change]chdir(work_dir.c_str());  //Go to the working directory //CO20191112
    }

    //[CO20191112] - NO PDF JPG detected in aflow_gnuplot_funcs.cpp, so this is useless anyway
    //[CO20191112 - dangerous command, will delete ANY pdf jpg created before this routine]//KY adds it
    //[CO20191112 - dangerous command, will delete ANY pdf jpg created before this routine]command << "cd \"" << directory_RAW << "\"" << endl;
    //[CO20191112 - dangerous command, will delete ANY pdf jpg created before this routine]command << "rm -f *pdf *jpg " << endl;
    //[CO20191112 - dangerous command, will delete ANY pdf jpg created before this routine]aurostd::execute(command);

    // DONE
    cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Bands: end " << directory_LIB << endl;
    return TRUE;
  }
}
// ***************************************************************************
// aflowlib::LIB2RAW_Loop_DATA
// ***************************************************************************
// namespace aflowlib {
//   bool LIB2RAW_Loop_DATA(const string& directory_LIB,const string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,const string& MESSAGE) {
//     bool LDEBUG=(FALSE || XHOST.DEBUG);
//     // LDEBUG=TRUE;
//     if(LDEBUG) cerr << XPID << "aflowlib::LIB2RAW_Loop_DATAs [1]" << endl;
//     // Stefano Curtarolo 2009-2010-2011-2012
//     vector<string> vspecies;aurostd::string2tokens(data.species,vspecies,",");

//     cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATAs: begin " << directory_LIB << endl;
//     cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATAs: species = " << vspecies.size() << endl;
//     data.vloop.push_back("data");

//     uint relax_max=10;
//     bool flag_EDATA_ORIG_=FALSE,flag_EDATA_RELAX_=FALSE;
//     bool flag_DATA_ORIG_=FALSE,flag_DATA_RELAX_=FALSE,flag_DATA_BANDS_=FALSE;


//     aurostd::StringSubst(directory_LIB,"ELPASOLITES","AURO"); // PATCH
//     aurostd::StringSubst(directory_LIB,"SCINT","ICSD"); // PATCH
//     if(aurostd::substring2bool(directory_LIB,"ICSD"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB0"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB1"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB2"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB3"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB4"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB5"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB6"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB7"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB8"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIB9"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"AURO"))     { flag_EDATA_ORIG_=TRUE; }
//     if(aurostd::substring2bool(directory_LIB,"LIBRARYX")) { flag_EDATA_ORIG_=TRUE; } // [HISTORIC]

//     string file_LIB,file_RAW;
//     aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWIN_,directory_RAW,_AFLOWIN_,vfile,MESSAGE);  // _AFLOWIN_
//     aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWLOCK_,directory_RAW,_AFLOWLOCK_,vfile,MESSAGE);  // LOCK

//     if(flag_DATA_ORIG_ || flag_EDATA_ORIG_) {  // POSCAR.orig.EXT
//       // if(flag_ORIG==FALSE) {
//       bool found=FALSE;
//       file_LIB=directory_LIB+"/POSCAR.orig"+EXT;
//       if(!found && (found=aurostd::FileExist(file_LIB))) {
// 	cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATAs: building POSCAR.orig from POSCAR.orig"+EXT << endl;
// 	aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.orig",directory_RAW,"POSCAR.orig",vfile,MESSAGE);  // POSCAR.orig
//       }
//       file_LIB=directory_LIB+"/POSCAR.relax1"+EXT;
//       if(!found && (found=aurostd::FileExist(file_LIB))) {
// 	cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATAs: building POSCAR.orig from POSCAR.relax1"+EXT << endl;
// 	aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.relax1",directory_RAW,"POSCAR.orig",vfile,MESSAGE);  // POSCAR.orig
//       }
//       if(!found) {
// 	found=TRUE;
// 	cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATAs: building POSCAR.orig from " << _AFLOWIN_ << "" << endl;
// 	aurostd::execute("cat \""+directory_RAW+"/" + _AFLOWIN_ + "\" | aflow --justbetween=\"[VASP_POSCAR_MODE_EXPLICIT]START\",\"[VASP_POSCAR_MODE_EXPLICIT]STOP\" > "+directory_RAW+"/POSCAR.orig");
// 	//  ExtractToStringEXPLICIT(Library_ICSD,Library_ICSD0,"[README_LIBRARY_ICSD1.TXT]START","[README_LIBRARY_ICSD1.TXT]STOP");
//       }
//     }

//     if(flag_DATA_RELAX_ || flag_EDATA_RELAX_) {  // CONTCAR.relax.EXT
//       for(uint i=1;i<=relax_max;i++) {
// 	file_LIB=directory_LIB+"/CONTCAR.relax"+aurostd::utype2string(i)+EXT;file_RAW=directory_RAW+"/CONTCAR.relax"+EXT;
// 	if(aurostd::FileExist(file_LIB)) { aurostd::CopyFile(file_LIB,file_RAW);aurostd::execute(DEFAULT_KZIP_BIN+" -9fd \""+file_RAW+"\"");vfile.push_back("CONTCAR.relax"); }
//       }
//       if(aurostd::FileExist(file_RAW)) { cout << MESSAGE << ": ERROR - aflowlib::LIB2RAW_Loop_DATAs:[1] - file not prepared " << file_LIB << endl;exit(0); }
//     }

//     if(!flag_DATA_BANDS_) {
//       file_LIB=directory_LIB+"/POSCAR.bands"+EXT;
//       if(aurostd::FileExist(file_LIB)) { aurostd::CopyFile(file_LIB,file_RAW);aurostd::execute(DEFAULT_KZIP_BIN+" -9fd \""+file_RAW+"\"");vfile.push_back("POSCAR.bands"); }
//     }

//     xstructure str,str_sp,str_sc;
//     vector<xstructure> vcif;

//     // PERFORM EDATA STEP
//     if(flag_EDATA_ORIG_ || flag_EDATA_RELAX_) cout << MESSAGE << " EDATA start: " << directory_RAW << endl;
//     if(flag_EDATA_ORIG_) { // ORIG
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.orig")) {
// 	cout << MESSAGE << " EDATA doing orig (POSCAR.orig): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"EDATA"); // 1=EDATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT);
// 	vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
//       }
//     }
//     if(flag_EDATA_RELAX_) { // RELAX
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT) && aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
// 	cout << MESSAGE << " EDATA doing relax (CONTCAR.relax): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"EDATA"); // EDATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT);
// 	vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
//       }
//     }

//     //  if(flag_EDATA_BANDS_)
//     { // BANDS IF AVAILABLE DO IT
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_BANDS_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.bands")) {
// 	cout << MESSAGE << " EDATA doing bands (POSCAR.bands): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/POSCAR.bands",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"EDATA"); // EDATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_BANDS_OUT);
// 	vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
//       }
//     }

//     // PERFORM DATA STEP
//     if(flag_DATA_ORIG_ || flag_DATA_RELAX_ || flag_DATA_BANDS_) cout << MESSAGE << " DATA start: " << directory_RAW << endl;
//     if(flag_DATA_ORIG_) { // ORIG
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_ORIG_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.orig")) {
// 	cout << MESSAGE << " DATA doing orig (POSCAR.orig): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"DATA"); // DATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_ORIG_OUT);
//       }
//     }
//     if(flag_DATA_RELAX_) { // RELAX
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_RELAX_OUT) && aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
// 	cout << MESSAGE << " DATA doing relax (CONTCAR.relax): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"DATA"); // DATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_RELAX_OUT);
//       }
//     }
//     if(flag_DATA_BANDS_) { // BANDS
//       if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_BANDS_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.bands")) {
// 	cout << MESSAGE << " DATA doing bands (POSCAR.bands): " << directory_RAW << endl;
// 	str=xstructure(directory_RAW+"/POSCAR.bands",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
// 	stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
// 	pflow::PrintData(str,str_sp,str_sc,sss,"DATA"); // DATA
// 	aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_BANDS_OUT);
//       }
//     }

//     // NOW DO THE CIFS
//     aflowlib::LIB2RAW_FileNeeded(directory_LIB,"CONTCAR.relax2",directory_RAW,"CONTCAR.relax2",vfile,MESSAGE);  // POSCAR.bands

//     //      cout << MESSAGE << " CIF creation: " << directory_LIB << endl;
//     if(vcif.size()==0) vcif.push_back(xstructure(directory_RAW+"/POSCAR.bands",IOVASP_AUTO));
//     xvector<double> nvec(3);nvec(1)=1;nvec(2)=1;nvec(3)=1;
//     double angle=45;

//     for (uint j=0;j<vcif.size();j++) {
//       vcif.at(j)=GetLTFVCell(nvec,angle,vcif.at(j));
//       if(j==0) cout << MESSAGE << " CIF creation: data.spacegroup_relax=" << data.spacegroup_relax << endl;
//       if(j==0) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing normal" << endl;
//       if(j==1) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sprim" << endl;
//       if(j==2) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sconv" << endl;
//       stringstream oss;
//       //    for(uint i=0;i<vspecies.size();i++) cerr << XPID << vspecies.at(i) << endl;
//       vcif.at(j).species.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species.push_back(vspecies.at(i));
//       vcif.at(j).species_pp.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp.push_back(vspecies.at(i));
//       vcif.at(j).species_pp_type.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_type.push_back("");
//       vcif.at(j).species_pp_version.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_version.push_back("");
//       vcif.at(j).species_pp_ZVAL.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_ZVAL.push_back(0.0);
//       vcif.at(j).species_pp_vLDAU.clear();for(uint i=0;i<vspecies.size();i++) vcif.at(j).species_pp_vLDAU.push_back(deque<double>());
//       for(uint i=0;i<vcif.at(j).atoms.size();i++) vcif.at(j).atoms.at(i).name=vcif.at(j).species.at(vcif.at(j).atoms.at(i).type);
//       for(uint i=0;i<vcif.at(j).atoms.size();i++) vcif.at(j).atoms.at(i).cleanname=vcif.at(j).species.at(vcif.at(j).atoms.at(i).type);
//       pflow::PrintCIF(oss,vcif.at(j),1);//aurostd::string2utype<int>(data.spacegroup_relax));
//       if(j==0) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+".cif");
//       if(j==1) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sprim.cif");
//       if(j==2) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sconv.cif");
//       vcif.at(j).AddCorners();
//       oss.clear();oss.str("");
//       pflow::PrintCIF(oss,vcif.at(j),1);//aurostd::string2utype<int>(data.spacegroup_relax));
//       if(j==0) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing normal_corner" << endl;
//       if(j==1) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sprim_corner" << endl;
//       if(j==2) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sconv_corner" << endl;
//       if(j==0) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_corner.cif");
//       if(j==1) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sprim_corner.cif");
//       if(j==2) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sconv_corner.cif");

//       // [OBSOLETE] aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+".cif");
//     }

//     // DONE
//     cout << MESSAGE << " aflowlib::LIB2RAW_Loop_DATA: end " << directory_LIB << endl;
//     return TRUE;
//   }
// }


// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Thermodynamics
// ***************************************************************************
namespace aflowlib {
  bool ExtractOUT_from_VASP_OUTCAR(string _file,const double& data_natoms,xOUTCAR& xOUT) {
    bool flag=xOUT.GetPropertiesFile(_file);
    if(aurostd::abs(data_natoms-(double) xOUT.natoms)>0.1) { cerr << XPID << "ERROR ExtractOUT_from_VASP_OUTCAR: data_natoms(" << data_natoms << ")!= (int) xOUT.natoms(" << xOUT.natoms << ") ..." << endl;exit(0); }
    return flag;
  }
}

namespace aflowlib {
  bool AddFileNameBeforeExtension(string _file,string addendum,string& out_file) { //CO20171025
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string file=aurostd::CleanFileName(_file);
    out_file=file;
    if(file[file.size()-1]=='/' || file[file.size()-1]=='\\') { return false; } //not doing directories
    //if(aurostd::IsDirectory(file)) { return false; }
    //if(!aurostd::FileExist(file)) { return false; }
    //grab actually file
    string path="",rfile="";
    std::size_t pos=file.find_last_of("/\\");
    if(pos!=string::npos) { path=file.substr(0,pos); } //keep path empty if there's no directory info (just relative file name)
    rfile=file.substr(pos+1);
    if(LDEBUG) {
      cerr << XPID << "aflowlib::AddFileNameBeforeExtension:: path=" << path << endl;
      cerr << XPID << "aflowlib::AddFileNameBeforeExtension:: file=" << rfile << endl;
    }
    string possible_extensions;
    possible_extensions="out,txt,json,png,pdf,jpg,jpeg";  //normal file types first
    possible_extensions+=",tar,tbz,bz2,zip,gz";           //zip extensions last
    vector<string> vpossible_extensions;
    aurostd::string2tokens(possible_extensions,vpossible_extensions,",");
    //now split file by "."
    vector<string> parts;
    aurostd::string2tokens(rfile,parts,".");
    bool found=false;
    string rfile_new="";
    for(uint i=0;i<parts.size();i++) {
      for(uint j=0;j<vpossible_extensions.size()&&!found;j++) {
        if(parts[i]==vpossible_extensions[j]) {
          found=true;
          rfile_new+=(rfile_new.empty()?"":".")+addendum;
        }
      }
      rfile_new+=(rfile_new.empty()?"":".")+parts[i];
    }
    if(!found) { rfile_new+=(rfile_new.empty()?"":".")+addendum; }
    out_file=aurostd::CleanFileName((path.empty()?"":path)+(pos==string::npos?"":string(file.substr(pos,1)))+rfile_new);
    return true;
  }
}

namespace aflowlib {
  bool LIB2RAW_Loop_Thermodynamics(const string& directory_LIB,const string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,const string& MESSAGE,bool LOCAL) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"aflowlib::LIB2RAW_Loop_Thermodynamics():";
    stringstream message; //dummy stringstream, can output to cout, but not used right now

    if(LDEBUG) cerr << soliloquy << " [-1]" << endl;
    // ZIP-AGNOSTIC
    if(XHOST.vext.size()!=XHOST.vzip.size()) { cerr << XPID << "ERROR - " << soliloquy << " XHOST.vext.size()!=XHOST.vzip.size()" << endl;exit(0); }
    //CO+DX START 20170713 - adding symmetry output to RAW
    _aflags aflags;
    aflags.Directory=directory_RAW;
    ofstream FileMESSAGE; //dummy ofstream, not really used
    if(LDEBUG) cerr << soliloquy << " [-2]" << endl;
    string system_name=KBIN::ExtractSystemName(directory_LIB);
    if(LDEBUG) cerr << soliloquy << " [-3]" << endl;
    //CO+DX STOP 20170713 - adding symmetry output to RAW
    //    if(LDEBUG) cerr << soliloquy << " data.aurl=" << data.aurl << endl;
    // aurostd::StringSubst(directory_LIB,"/./","/");
    // aurostd::StringSubst(directory_RAW,"/./","/");
    //[CO20200624 - OBSOLETE]vector<string> vspecies;aurostd::string2tokens(data.species,vspecies,",");
    deque<string> deq_species;aurostd::string2tokens(data.species,deq_species,","); //DX20190620
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << soliloquy << " - begin " << directory_LIB << endl;
    if(LDEBUG) cerr << XPID << "directory_LIB=\"" << directory_LIB << "\"" << endl;
    if(LDEBUG) cerr << XPID << "directory_RAW=\"" << directory_RAW << "\"" << endl;
    //[CO20200624 - MOVING UP]if(directory_LIB.at(directory_LIB.size()-1)=='/')  directory_LIB=directory_LIB.substr(0,directory_LIB.size()-1);
    //[CO20200624 - MOVING UP]if(directory_RAW.at(directory_RAW.size()-1)=='/')  directory_RAW=directory_RAW.substr(0,directory_RAW.size()-1);
    if(LDEBUG) cerr << XPID << "directory_LIB=\"" << directory_LIB << "\"" << endl;
    if(LDEBUG) cerr << XPID << "directory_RAW=\"" << directory_RAW << "\"" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << soliloquy << " - species = " << data.vspecies.size() << endl;
    uint relax_max=10;
    stringstream command;command.clear();command.str(std::string());
    stringstream aus_exec;
    // directories must exist already
    bool flag_EDATA_ORIG_=FALSE,flag_EDATA_RELAX_=FALSE,flag_EDATA_BANDS_=FALSE;
    bool flag_DATA_ORIG_=FALSE,flag_DATA_RELAX_=FALSE,flag_DATA_BANDS_=FALSE;
    bool flag_TIMING=FALSE;
    bool flag_ENERGY1=FALSE;
    bool flag_SG1=FALSE;
    bool flag_SG2=FALSE;
    bool flag_VOLDISTPARAMS=FALSE;
    bool flag_VOLDISTEVOLUTION=FALSE;
    bool flag_ICSD=FALSE,flag_MAGNETIC=FALSE,flag_LIB0=FALSE,flag_LIB1=FALSE,flag_LIB2=FALSE;
    bool flag_ERROR=FALSE;

    string fileA_LIB,fileA_RAW,fileE_LIB,fileE_RAW,fileX_LIB,fileX_RAW,fileK_LIB,fileK_RAW,fileI_LIB,fileI_RAW,fileJ_LIB,fileJ_RAW,cmd;
    string FileName_OUTCAR_relax="";

    if(LOCAL) {
      flag_EDATA_ORIG_=flag_EDATA_RELAX_=flag_EDATA_BANDS_=TRUE;
      flag_DATA_ORIG_=flag_DATA_RELAX_=flag_DATA_BANDS_=TRUE;
      flag_TIMING=TRUE;
      flag_ENERGY1=FALSE;
      flag_SG1=TRUE;
      flag_SG2=TRUE;
      flag_VOLDISTPARAMS=TRUE;
      flag_VOLDISTEVOLUTION=TRUE;
      flag_ICSD=flag_LIB0=flag_LIB1=flag_LIB2=FALSE;
      flag_MAGNETIC=TRUE;
      flag_ERROR=FALSE;
    } else {
      //[CO20200624 - do inside CleanDirectoryLIB(), PROBABLY OBSOLETE]aurostd::StringSubst(directory_LIB,"ELPASOLITES","AURO"); // PATCH
      //[CO20200624 - do inside CleanDirectoryLIB(), PROBABLY OBSOLETE]aurostd::StringSubst(directory_LIB,"SCINT","ICSD"); // PATCH
      if(aurostd::substring2bool(directory_LIB,"ICSD"))     { flag_ICSD=TRUE;    flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB0"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=FALSE;flag_EDATA_RELAX_=FALSE;flag_TIMING=TRUE;flag_SG1=FALSE;flag_SG2=FALSE;flag_VOLDISTPARAMS=FALSE;flag_VOLDISTEVOLUTION=FALSE; }
      if(aurostd::substring2bool(directory_LIB,"LIB1"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB2"))     { flag_LIB2=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB3"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB4"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB5"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB6"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB7"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB8"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIB9"))     { flag_MAGNETIC=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; }
      if(aurostd::substring2bool(directory_LIB,"LIBRARYX")) { flag_LIB2=TRUE;flag_EDATA_ORIG_=TRUE;flag_EDATA_RELAX_=TRUE;flag_TIMING=TRUE;flag_SG1=TRUE;flag_SG2=TRUE;flag_VOLDISTPARAMS=TRUE;flag_VOLDISTEVOLUTION=TRUE; } // [HISTORIC]
    }
    flag_VOLDISTEVOLUTION=false;  //CO20200624 - simply churns, does nothing

    flag_LIB0=aurostd::substring2bool(directory_LIB,"LIB0");

    // check for flag_EDATA_BANDS_
    if(flag_EDATA_ORIG_ && flag_EDATA_RELAX_) {
      for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
        fileA_LIB=directory_LIB+"/POSCAR.bands"+XHOST.vext.at(iext);
        fileA_RAW=directory_RAW+"/POSCAR.bands"+XHOST.vext.at(iext);
        if(aurostd::FileExist(fileA_LIB)) {
          flag_EDATA_BANDS_=TRUE;
          aurostd::CopyFile(fileA_LIB,fileA_RAW);
          aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileA_RAW+"\"");
          vfile.push_back("POSCAR.bands");
        }
      }
    }
    // check for flag_DATA_BANDS_
    if(flag_DATA_ORIG_ && flag_DATA_RELAX_) {
      for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
        fileA_LIB=directory_LIB+"/POSCAR.bands"+XHOST.vext.at(iext);
        fileA_RAW=directory_RAW+"/POSCAR.bands"+XHOST.vext.at(iext);
        if(aurostd::FileExist(fileA_LIB)) {
          flag_DATA_BANDS_=TRUE;
          aurostd::CopyFile(fileA_LIB,fileA_RAW);
          aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileA_RAW+"\"");
          vfile.push_back("POSCAR.bands");
        }
      }
    }
    if(flag_ICSD) { ; } // dummy load
    if(flag_MAGNETIC) { ; } // dummy load
    if(flag_LIB0) { ; } // dummy load
    if(flag_LIB1) { ; } // dummy load

    xstructure str_orig,str_relax1,str_relax;

    //DX START
    //DX20180526 [OBSOLETE] str_orig.directory = str_relax1.directory = str_relax.directory = aflags.Directory;
    //DX END

    vector<string> tokens;

    double data1_energy_cell=0.0,data1_energy_atom=0.0;
    string data_sg1_pre="",data_sg1_mid="",data_sg1_post="",data_sg2_pre="",data_sg2_mid="",data_sg2_post="";
    string data_v_atom="",data_ucelld="";
    xvector<double> data_abcabc;
    vector<double> data_vcomposition;
    vector<xvector<double> > data_vforces;                      // QM FORCES calculation
    vector<xvector<double> > data_vpositions_cartesian;         // QM POSITIONS_CARTESIAN calculation

    xOUTCAR xOUT;
    // xDOSCAR xDOS;
    // xEIGENVAL xEIG;
    xKPOINTS kpoints;

    if(LDEBUG) cerr << soliloquy << " [2]" << endl;
    // copy _AFLOWIN_ LOCK
    // _AFLOWIN_
    data.vloop.push_back("thermodynamics");
    // star

    aurostd::string2tokens(directory_LIB,tokens,"/");
    data.prototype=tokens.at(tokens.size()-1);
    aurostd::StringSubst(data.prototype,":LDAU2","");
    aurostd::StringSubst(data.prototype,"\n","");aurostd::StringSubst(data.prototype," ","");aurostd::StringSubst(data.prototype," ","");

    // FILES
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWIN_,directory_RAW,_AFLOWIN_,vfile,MESSAGE);  // _AFLOWIN_
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWLOCK_,directory_RAW,_AFLOWLOCK_,vfile,MESSAGE);  // LOCK

    if(LDEBUG) cerr << soliloquy << " [3]" << endl;
    if(TRUE || flag_DATA_ORIG_ || flag_EDATA_ORIG_ || flag_SG1 || flag_SG2) {  // POSCAR.orig.EXT
      if(LDEBUG) cerr << soliloquy << " [3.1]" << endl;
      // if(flag_ORIG==FALSE) { //[CO20200106 - close bracket for indenting]}
      bool found=FALSE;
      for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
        if(LDEBUG) cerr << soliloquy << ": BUILDING POSCAR.orig from POSCAR.orig.EXT" << endl;
        if(!found && aurostd::FileExist(directory_LIB+"/POSCAR.orig"+XHOST.vext.at(iext))) {
          found=TRUE;
          if(LDEBUG) cerr << soliloquy << ": BUILDING POSCAR.orig from POSCAR.orig" << XHOST.vext.at(iext) << endl;
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.orig",directory_RAW,"POSCAR.orig",vfile,MESSAGE);  // POSCAR.orig
        }
      }
      for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
        if(LDEBUG) cerr << soliloquy << ": BUILDING POSCAR.orig from POSCAR.relax1.EXT" << endl;
        if(!found && aurostd::FileExist(directory_LIB+"/POSCAR.relax1"+XHOST.vext.at(iext))) {
          found=TRUE;
          if(LDEBUG) cerr << soliloquy << ": BUILDING POSCAR.orig from POSCAR.relax1" << XHOST.vext.at(iext) << endl;
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.relax1",directory_RAW,"POSCAR.orig",vfile,MESSAGE);  // POSCAR.orig
        }
      }
      if(!found) {
        found=TRUE;
        if(LDEBUG) cerr << soliloquy << ": BUILDING POSCAR.orig from " << _AFLOWIN_ << "" << endl;
        aurostd::execute(string("cat ")+"\""+directory_RAW+"/"+_AFLOWIN_+"\""+" | aflow --justbetween=\"[VASP_POSCAR_MODE_EXPLICIT]START\",\"[VASP_POSCAR_MODE_EXPLICIT]STOP\" > \""+directory_RAW+"/POSCAR.orig"+"\"");
        //    ExtractToStringEXPLICIT(Library_ICSD,Library_ICSD0,"[README_LIBRARY_ICSD1.TXT]START","[README_LIBRARY_ICSD1.TXT]STOP");
      }
    }
    if(LDEBUG) cerr << soliloquy << " [4]" << endl;
    if(!flag_LIB0) { // no LIB0
      if(TRUE || flag_DATA_RELAX_ || flag_EDATA_RELAX_ || TRUE || flag_SG1 || flag_SG2) {  // CONTCAR.relax.EXT
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          for(uint i=1;i<=relax_max;i++) {
            fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.relax"+aurostd::utype2string(i)+XHOST.vext.at(iext));
            fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+XHOST.vext.at(iext));
            fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.relax"+aurostd::utype2string(i)+XHOST.vext.at(iext));
            fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+XHOST.vext.at(iext));
            fileK_LIB=aurostd::CleanFileName(directory_LIB+"/KPOINTS.relax"+aurostd::utype2string(i)+XHOST.vext.at(iext));
            fileK_RAW=aurostd::CleanFileName(directory_RAW+"/KPOINTS.relax"+XHOST.vext.at(iext));
            fileI_LIB=aurostd::CleanFileName(directory_LIB+"/INCAR.relax"+aurostd::utype2string(i)+XHOST.vext.at(iext));
            fileI_RAW=aurostd::CleanFileName(directory_RAW+"/INCAR.relax"+XHOST.vext.at(iext));
            if(aurostd::FileExist(fileX_LIB)) {
              aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileX_RAW+"\"");vfile.push_back("CONTCAR.relax");
            }
            if(aurostd::FileExist(fileE_LIB)) {
              aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileE_RAW+"\"");vfile.push_back("OUTCAR.relax");
              FileName_OUTCAR_relax=fileE_LIB;
            }
            // if(aurostd::FileExist(fileK_LIB)) {
            //   aurostd::CopyFile(fileK_LIB,fileK_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileK_RAW+"\"");vfile.push_back("KPOINTS.relax");
            // }
            // if(aurostd::FileExist(fileI_LIB)) {
            //   aurostd::CopyFile(fileI_LIB,fileI_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileI_RAW+"\"");vfile.push_back("INCAR.relax");
            // }
          }
        }
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          if(!aurostd::FileExist(directory_RAW+"/CONTCAR.relax") && aurostd::FileExist(directory_LIB+"/CONTCAR.static"+XHOST.vext.at(iext)) &&
              !aurostd::FileExist(directory_RAW+"/OUTCAR.relax") && aurostd::FileExist(directory_LIB+"/OUTCAR.static"+XHOST.vext.at(iext)) &&
              !aurostd::FileExist(directory_RAW+"/KPOINTS.relax") && aurostd::FileExist(directory_LIB+"/KPOINTS.static"+XHOST.vext.at(iext))) {
            fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.static"+XHOST.vext.at(iext));
            fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+XHOST.vext.at(iext));
            fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.static"+XHOST.vext.at(iext));
            fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+XHOST.vext.at(iext));
            fileK_LIB=aurostd::CleanFileName(directory_LIB+"/KPOINTS.static"+XHOST.vext.at(iext));
            fileK_RAW=aurostd::CleanFileName(directory_RAW+"/KPOINTS.relax"+XHOST.vext.at(iext));
            //  if(AFLOWLIB_VERBOSE)
            cout << MESSAGE << " WARNING - PATCHING CONTCAR.relax with CONTCAR.static " << fileX_LIB << endl;
            //	  if(AFLOWLIB_VERBOSE)
            cout << MESSAGE << " WARNING - PATCHING OUTCAR.relax with OUTCAR.static " << fileE_LIB << endl;
            //	  if(AFLOWLIB_VERBOSE)
            cout << MESSAGE << " WARNING - PATCHING KPOINTS.relax with KPOINTS.static " << fileE_LIB << endl;
            if(aurostd::FileExist(fileX_LIB)) {
              aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileX_RAW+"\"");vfile.push_back("CONTCAR.relax");
            }
            if(aurostd::FileExist(fileE_LIB)) {
              aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileE_RAW+"\"");vfile.push_back("OUTCAR.relax");
              FileName_OUTCAR_relax=fileE_LIB;
            }
            if(aurostd::FileExist(fileK_LIB)) {
              aurostd::CopyFile(fileK_LIB,fileK_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileK_RAW+"\"");vfile.push_back("KPOINTS.relax");
            }
          }
        }
        if(aurostd::FileExist(fileX_RAW)) { cout << MESSAGE << " ERROR - " << soliloquy << " [1] - file not prepared " << fileX_LIB << endl;exit(0); }
        if(aurostd::FileExist(fileE_RAW)) { cout << MESSAGE << " ERROR - " << soliloquy << " [2] - file not prepared " << fileE_LIB << endl;exit(0); }
        if(aurostd::FileExist(fileK_RAW)) { cout << MESSAGE << " ERROR - " << soliloquy << " [3] - file not prepared " << fileK_LIB << endl;exit(0); }
        if(aurostd::FileExist(fileI_RAW)) { cout << MESSAGE << " ERROR - " << soliloquy << " [3] - file not prepared " << fileI_LIB << endl;exit(0); }
      }
    } // no LIB0

    if(flag_LIB0) {
      for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
        if(aurostd::FileExist(directory_LIB+"/CONTCAR.static"+XHOST.vext.at(iext)) &&
            aurostd::FileExist(directory_LIB+"/OUTCAR.static"+XHOST.vext.at(iext))) {
          fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.static"+XHOST.vext.at(iext));
          fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+XHOST.vext.at(iext));
          fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.static"+XHOST.vext.at(iext));
          fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+XHOST.vext.at(iext));
          fileK_LIB=aurostd::CleanFileName(directory_LIB+"/KPOINTS.static"+XHOST.vext.at(iext));
          fileK_RAW=aurostd::CleanFileName(directory_RAW+"/KPOINTS.relax"+XHOST.vext.at(iext));
          // if(AFLOWLIB_VERBOSE)
          cout << MESSAGE << " WARNING - PATCHING CONTCAR.relax with CONTCAR.static " << fileX_LIB << endl;
          // if(AFLOWLIB_VERBOSE)
          cout << MESSAGE << " WARNING - PATCHING OUTCAR.relax with OUTCAR.static " << fileE_LIB << endl;
          // if(AFLOWLIB_VERBOSE)
          cout << MESSAGE << " WARNING - PATCHING KPOINTS.relax with KPOINTS.static " << fileK_LIB << endl;
          if(aurostd::FileExist(fileX_LIB)) {
            aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileX_RAW+"\"");vfile.push_back("CONTCAR.relax");
          }
          if(aurostd::FileExist(fileE_LIB)) {
            aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileE_RAW+"\"");vfile.push_back("OUTCAR.relax");
            FileName_OUTCAR_relax=fileE_LIB;
          }
          if(aurostd::FileExist(fileK_LIB)) {
            aurostd::CopyFile(fileK_LIB,fileK_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileK_RAW+"\"");vfile.push_back("KPOINTS.relax");
          }
        }
      }
      if(aurostd::FileExist(fileX_RAW)) { cout << MESSAGE << " ERROR - " << soliloquy << " [1b] - file not prepared " << fileX_LIB << endl;exit(0); }
      if(aurostd::FileExist(fileE_RAW)) { cout << MESSAGE << " ERROR - " << soliloquy << " [2b] - file not prepared " << fileE_LIB << endl;exit(0); }
      if(aurostd::FileExist(fileK_RAW)) { cout << MESSAGE << " ERROR - " << soliloquy << " [3b] - file not prepared " << fileK_LIB << endl;exit(0); }
      if(aurostd::FileExist(fileI_RAW)) { cout << MESSAGE << " ERROR - " << soliloquy << " [3b] - file not prepared " << fileI_LIB << endl;exit(0); }
    }

    if(LDEBUG) cerr << soliloquy << " [5]" << endl;
    if(flag_SG1 || flag_SG2) {  // CONTCAR.relax1
      for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
        fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.relax1"+XHOST.vext.at(iext));
        fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax1"+XHOST.vext.at(iext));
        fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.relax1"+XHOST.vext.at(iext));
        fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax1"+XHOST.vext.at(iext));
        if(aurostd::FileExist(fileX_LIB)) {
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,"CONTCAR.relax1",directory_RAW,"CONTCAR.relax1",vfile,MESSAGE);
        }
        if(aurostd::FileExist(fileE_LIB)) {
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.relax1",directory_RAW,"OUTCAR.relax1",vfile,MESSAGE);
        }
        if(!aurostd::FileExist(directory_RAW+"/CONTCAR.relax1") && aurostd::FileExist(directory_LIB+"/CONTCAR.static"+XHOST.vext.at(iext)) &&
            !aurostd::FileExist(directory_RAW+"/OUTCAR.relax1") && aurostd::FileExist(directory_LIB+"/OUTCAR.static"+XHOST.vext.at(iext))) {
          fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.static"+XHOST.vext.at(iext));
          fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax1"+XHOST.vext.at(iext));
          fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.static"+XHOST.vext.at(iext));
          fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax1"+XHOST.vext.at(iext));
          //  if(AFLOWLIB_VERBOSE)
          cout << MESSAGE << " WARNING - PATCHING CONTCAR.relax1 with CONTCAR.static " << fileX_LIB << endl;
          //	  if(AFLOWLIB_VERBOSE)
          cout << MESSAGE << " WARNING - PATCHING OUTCAR.relax1 with OUTCAR.static " << fileE_LIB << endl;
          if(aurostd::FileExist(fileX_LIB)) {
            aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileX_RAW+"\"");vfile.push_back("CONTCAR.relax1");
          }
          if(aurostd::FileExist(fileE_LIB)) {
            aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileE_RAW+"\"");vfile.push_back("OUTCAR.relax1");
          }
        }
        if(aurostd::FileExist(fileX_RAW)) {
          cout << MESSAGE << " ERROR - " << soliloquy << " [4] - file not prepared " << fileX_LIB << endl;exit(0); }
        if(aurostd::FileExist(fileE_RAW)) {
          cout << MESSAGE << " ERROR - " << soliloquy << " [5] - file not prepared " << fileE_LIB << endl;exit(0); }
      }
    }
    if(LDEBUG) cerr << soliloquy << " [6]" << endl;
    if(flag_DATA_RELAX_ || flag_EDATA_RELAX_ || TRUE) {  // OSZICAR.relax.EXT
      for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
        // trying to get an OUTCAR
        if(aurostd::FileExist(directory_LIB+"/OUTCAR.relax1"+XHOST.vext.at(iext)))
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.relax1",directory_RAW,"OUTCAR.relax1",vfile,MESSAGE);  // _AFLOWIN_
        if(aurostd::FileExist(directory_LIB+"/CONTCAR.relax1"+XHOST.vext.at(iext)))
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,"CONTCAR.relax1",directory_RAW,"CONTCAR.relax1",vfile,MESSAGE);  // _AFLOWIN_
        // if(aurostd::FileExist(directory_LIB+"/KPOINTS.relax1"+XHOST.vext.at(iext)))
        // aflowlib::LIB2RAW_FileNeeded(directory_LIB,"KPOINTS.relax1",directory_RAW,"KPOINTS.relax1",vfile,MESSAGE);  // _AFLOWIN_
        fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.static"+XHOST.vext.at(iext));
        fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+XHOST.vext.at(iext));
        fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.static"+XHOST.vext.at(iext));
        fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+XHOST.vext.at(iext));
        fileK_LIB=aurostd::CleanFileName(directory_LIB+"/KPOINTS.static"+XHOST.vext.at(iext));
        fileK_RAW=aurostd::CleanFileName(directory_RAW+"/KPOINTS.relax"+XHOST.vext.at(iext));
        //SC STATIC=>RELAX
        //   if(aurostd::FileExist(fileE_LIB) && aurostd::FileExist(fileX_LIB) && aurostd::FileExist(fileK_LIB)) {
        //// a bug in old aflow the fileK_LIB was directory_LIB+"/KPOINTS.relax"+aurostd::utype2string(i)+XHOST.vext.at(iext) with i=3 or more, so it was never picked
        //// and this part was never done.
        //aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileE_RAW+"\"");vfile.push_back("OUTCAR.relax");
        //FileName_OUTCAR_relax=fileE_LIB;
        //// aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileX_RAW+"\"");vfile.push_back("CONTCAR.relax");
        //aurostd::CopyFile(fileK_LIB,fileK_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileK_RAW+"\"");vfile.push_back("KPOINTS.relax");
        //} else

        {
          fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.relax2"+XHOST.vext.at(iext));
          fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+XHOST.vext.at(iext));
          fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.relax2"+XHOST.vext.at(iext));
          fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+XHOST.vext.at(iext));
          fileK_LIB=aurostd::CleanFileName(directory_LIB+"/KPOINTS.relax2"+XHOST.vext.at(iext));
          fileK_RAW=aurostd::CleanFileName(directory_RAW+"/KPOINTS.relax"+XHOST.vext.at(iext));
          fileI_LIB=aurostd::CleanFileName(directory_LIB+"/INCAR.relax2"+XHOST.vext.at(iext));
          fileI_RAW=aurostd::CleanFileName(directory_RAW+"/INCAR.relax"+XHOST.vext.at(iext));
          if(aurostd::FileExist(fileE_LIB) && aurostd::FileExist(fileX_LIB) && aurostd::FileExist(fileK_LIB) && aurostd::FileExist(fileI_LIB)) {
            aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileE_RAW+"\"");vfile.push_back("OUTCAR.relax");
            FileName_OUTCAR_relax=fileE_LIB;
            aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileX_RAW+"\"");vfile.push_back("CONTCAR.relax");
            aurostd::CopyFile(fileK_LIB,fileK_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileK_RAW+"\"");vfile.push_back("KPOINTS.relax");
            // nocopy aurostd::CopyFile(fileI_LIB,fileI_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileI_RAW+"\"");vfile.push_back("INCAR.relax");
          } else {
            fileE_LIB=aurostd::CleanFileName(directory_LIB+"/OUTCAR.relax1"+XHOST.vext.at(iext));
            fileE_RAW=aurostd::CleanFileName(directory_RAW+"/OUTCAR.relax"+XHOST.vext.at(iext));
            fileX_LIB=aurostd::CleanFileName(directory_LIB+"/CONTCAR.relax1"+XHOST.vext.at(iext));
            fileX_RAW=aurostd::CleanFileName(directory_RAW+"/CONTCAR.relax"+XHOST.vext.at(iext));
            fileK_LIB=aurostd::CleanFileName(directory_LIB+"/KPOINTS.relax1"+XHOST.vext.at(iext));
            fileK_RAW=aurostd::CleanFileName(directory_RAW+"/KPOINTS.relax"+XHOST.vext.at(iext));
            fileI_LIB=aurostd::CleanFileName(directory_LIB+"/INCAR.relax1"+XHOST.vext.at(iext));
            fileI_RAW=aurostd::CleanFileName(directory_RAW+"/INCAR.relax"+XHOST.vext.at(iext));
            if(aurostd::FileExist(fileE_LIB) && aurostd::FileExist(fileX_LIB) && aurostd::FileExist(fileK_LIB) && aurostd::FileExist(fileI_LIB)) {
              aurostd::CopyFile(fileE_LIB,fileE_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileE_RAW+"\"");vfile.push_back("OUTCAR.relax");
              FileName_OUTCAR_relax=fileE_LIB;
              aurostd::CopyFile(fileX_LIB,fileX_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileX_RAW+"\"");vfile.push_back("CONTCAR.relax");
              aurostd::CopyFile(fileK_LIB,fileK_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileK_RAW+"\"");vfile.push_back("KPOINTS.relax");
              // nocopy aurostd::CopyFile(fileI_LIB,fileI_RAW);aurostd::execute(XHOST.vzip.at(iext)+" -dqf \""+fileI_RAW+"\"");vfile.push_back("INCAR.relax");
            }
          }
        }
      }
    }

    if(LDEBUG) cerr << soliloquy << " [6.9]" << endl;
    // get code
    data.code="nan";
    aurostd::string2tokens(aurostd::execute2string("cat \""+directory_RAW+"/OUTCAR.relax"+"\" | grep vasp | head -1"),tokens," ");
    if(tokens.size()>1) data.code=tokens.at(0);
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CODE = " << data.code << endl;

    // create structures
    //CO20171021 - and write xstr_json
    stringstream xstr_js;
    if(LDEBUG) cerr << soliloquy << " [7]" << endl;

    if(str_orig.num_each_type.size()==0 && aurostd::FileExist(directory_RAW+"/POSCAR.orig")) {
      xstructure _str_orig(directory_RAW+"/POSCAR.orig",IOVASP_AUTO);
      str_orig=_str_orig;
      str_orig.SetSpecies(deq_species); //DX20190620 - add species to xstructure
      str_orig.ReScale(1.0);
      //xstr_js.str("");xstr_js << xstructure2json(str_orig);aurostd::stringstream2file(xstr_js,directory_RAW+"/"+system_name+"_structure_orig.json");
    } //CO20171025

    if(LDEBUG) cerr << soliloquy << " [7.1]" << endl;

    if(str_relax.num_each_type.size()==0 && aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
      xstructure _str_relax(directory_RAW+"/CONTCAR.relax",IOVASP_AUTO);
      str_relax=_str_relax;
      str_relax.SetSpecies(deq_species); //DX20190620 - add species to xstructure
      str_relax.ReScale(1.0);
      xstr_js.str("");
      xstr_js << xstructure2json(str_relax); //DX20190620 - bug fix; orig->relax
      aurostd::stringstream2file(xstr_js,directory_RAW+"/"+system_name+"_structure_relax.json");
    } //CO20171025

    if(LDEBUG) cerr << soliloquy << " [7.2]" << endl;

    if(str_relax.num_each_type.size()==0 && aurostd::FileExist(directory_RAW+"/CONTCAR.static")) {
      xstructure _str_relax(directory_RAW+"/CONTCAR.static",IOVASP_AUTO);
      str_relax=_str_relax;
      str_relax.SetSpecies(deq_species); //DX20190620 - add species to xstructure
      str_relax.ReScale(1.0);
      xstr_js.str("");
      xstr_js << xstructure2json(str_relax); //DX20190620 - bug fix; orig->relax
      aurostd::stringstream2file(xstr_js,directory_RAW+"/"+system_name+"_structure_relax.json");
    } //CO20171025

    if(LDEBUG) cerr << soliloquy << " [7.3]" << endl;

    if(aurostd::FileExist(directory_RAW+"/CONTCAR.relax1")) {
      xstructure _str_relax1(directory_RAW+"/CONTCAR.relax1",IOVASP_AUTO);
      str_relax1=_str_relax1;
      str_relax1.SetSpecies(deq_species); //DX20190620 - add species to xstructure
      str_relax1.ReScale(1.0);
      xstr_js.str("");
      xstr_js << xstructure2json(str_relax1); //DX20190620 - bug fix; orig->relax1
      aurostd::stringstream2file(xstr_js,directory_RAW+"/"+system_name+"_structure_relax1.json");
    } //CO20171025
    // do the extractions

    if(LDEBUG) cerr << soliloquy << " [8]" << endl;

    // LOAD STRUCTURES
    data.nspecies=str_relax.num_each_type.size();
    data.natoms=str_relax.atoms.size();
    data.volume_cell=str_relax.GetVolume();
    data.volume_atom=str_relax.GetVolume()/(double) data.natoms;
    data_abcabc=Getabc_angles(str_relax.lattice,DEGREES);
    data.vgeometry.clear(); for(uint i=1;i<=6;i++) data.vgeometry.push_back(data_abcabc(i));

    //DX20190124 - add original crystal info - START
    // LOAD ORIGINAL STRUCTURE
    data.natoms_orig=str_orig.atoms.size();
    data.volume_cell_orig=str_orig.GetVolume();
    data.volume_atom_orig=str_orig.GetVolume()/(double) data.natoms_orig;
    data_abcabc=Getabc_angles(str_orig.lattice,DEGREES);
    data.vgeometry_orig.clear(); for(uint i=1;i<=6;i++) data.vgeometry_orig.push_back(data_abcabc(i));

    //DX20190124 - add original crystal info - END

    //CO, get fpos now, cpos comes from outcar later (either way works)
    vector<string> fpos_strings;
    vector<string> fpos_strings_combined;
    data.vpositions_fractional.clear();
    for(uint i=0;i<str_relax.atoms.size();i++) {
      data.vpositions_fractional.push_back(str_relax.atoms.at(i).fpos);
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " POSITIONS_FRACTIONAL = " << data.vpositions_fractional.at(i)[1] << "," << data.vpositions_fractional.at(i)[2] << "," << data.vpositions_fractional.at(i)[3] << "," << endl;
      //prepare for string variant
      for(uint j=1;j<(uint)str_relax.atoms.at(i).fpos.rows+1;j++) {
        fpos_strings.push_back(aurostd::utype2string(str_relax.atoms.at(i).fpos[j],8));
      }
      fpos_strings_combined.push_back(aurostd::joinWDelimiter(fpos_strings,","));
      fpos_strings.clear();
    }
    data.positions_fractional=aurostd::joinWDelimiter(fpos_strings_combined,";");

    data.geometry="";
    if(data.vgeometry.size()) {
      for(uint i=0;i<data.vgeometry.size();i++) {
        data.geometry+=aurostd::utype2string(data.vgeometry.at(i),7)+(i<data.vgeometry.size()-1?",":"");
      }
    }

    //DX20190124 - add original crystal info - START
    data.geometry_orig="";
    if(data.vgeometry_orig.size()) {
      for(uint i=0;i<data.vgeometry_orig.size();i++) {
        data.geometry_orig+=aurostd::utype2string(data.vgeometry_orig.at(i),7)+(i<data.vgeometry_orig.size()-1?",":"");
      }
    }
    //DX20190124 - add original crystal info - END

    data_vcomposition.clear();
    for(uint i=0;i<str_relax.num_each_type.size();i++) data_vcomposition.push_back(double(str_relax.num_each_type.at(i))/double(str_relax.atoms.size()));
    data.vstoichiometry=data_vcomposition;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NSPECIES = " << data.nspecies << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NATOMS = " << data.natoms << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " VOLUME (A^3) = " << data.volume_cell << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " VOLUME_ATOM (A^3) = " << data.volume_atom << "   " << directory_LIB << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " GEOMETRY (A,A,A,deg,deg,deg) = " << data.geometry << endl;
    // deque<int> num_each_type;              // WARNING: we use starting from 0
    // std::deque<_atom> atoms;               // WARNING: we use starting from 0
    // deque<string> species,species_pp;      // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    // deque<double> species_volume;          // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    // deque<double> species_mass;            // WARNING: we use starting from 0 // CAN BE THE ONES OF VASP5
    // cerr << XPID << str_relax.atoms.size() << " " << str_relax.atoms.at(0) << endl;
    // cerr << XPID << str_relax.num_each_type.size() << " " << str_relax.num_each_type.at(0) << endl;
    // cerr << XPID << str_relax.species.size() << " " << str_relax.species.size() << endl;
    // cerr << XPID << str_relax.species_pp.size() << " " << str_relax.species_pp.size() << endl;
    // cerr << XPID << str_relax.species_pp_type.size() << " " << str_relax.species_pp_type.size() << endl;
    // cerr << XPID << str_relax.species_pp_version.size() << " " << str_relax.species_pp_version.size() << endl;
    // cerr << XPID << str_relax.species_pp_ZVAL.size() << " " << str_relax.species_pp_ZVAL.size() << endl;
    // cerr << XPID << str_relax.species_pp_vLDAU.size() << " " << str_relax.species_pp_vLDAU.size() << endl;
    // cerr << XPID << str_relax.species_volume.size() << " " << str_relax.species_volume.size() << endl;
    // cerr << XPID << str_relax.species_mass.size() << " " << str_relax.species_mass.size() << endl;

    // LOAD ENERGY DATA1
    if(flag_ENERGY1) {
      xOUT.GetPropertiesFile(directory_RAW+"/OUTCAR.relax1",data.natoms,TRUE);
      //   ExtractDataOSZICAR(directory_RAW+"/OSZICAR.relax1",data.natoms,data1_dE,data1_dEN,data1_mag,data1_mag_atom);
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENERGY1 total E0 (eV) = " << (data1_energy_cell=xOUT.energy_cell) << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENERGY1 per atom E0/N (eV) = " << (data1_energy_atom=xOUT.energy_atom) << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY1 total E0 (eV) = " << xOUT.enthalpy_cell << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY1 per atom E0/N (eV) = " << xOUT.enthalpy_atom << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPIN1 mag (\\mu) = " << xOUT.mag_cell << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPIN1 per atom mag/N (\\mu) = " << xOUT.mag_atom << endl;
    }
    // LOAD ENERGY DATA
    xOUT.GetPropertiesFile(directory_RAW+"/OUTCAR.relax",data.natoms,TRUE);
    kpoints.GetPropertiesFile(directory_RAW+"/KPOINTS.relax",TRUE);

    data.pressure=xOUT.pressure;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PRESSURE (kB) = " << data.pressure << endl;
    data.vstress_tensor.clear();
    data.vstress_tensor.push_back(xOUT.stress(1,1));data.vstress_tensor.push_back(xOUT.stress(1,2));data.vstress_tensor.push_back(xOUT.stress(1,3));
    data.vstress_tensor.push_back(xOUT.stress(2,1));data.vstress_tensor.push_back(xOUT.stress(2,2));data.vstress_tensor.push_back(xOUT.stress(2,3));
    data.vstress_tensor.push_back(xOUT.stress(3,1));data.vstress_tensor.push_back(xOUT.stress(3,2));data.vstress_tensor.push_back(xOUT.stress(3,3));
    data.stress_tensor="";
    if(data.vstress_tensor.size())
      for(uint i=0;i<data.vstress_tensor.size();i++)
        data.stress_tensor+=aurostd::utype2string(data.vstress_tensor.at(i),7)+(i<data.vstress_tensor.size()-1?",":"");
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " STRESS_TENSOR (kB) = " << data.stress_tensor << endl;
    data.pressure_residual=xOUT.pressure_residual;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PRESSURE_RESIDUAL (kB) = " << data.pressure_residual << endl;
    data.Pulay_stress=xOUT.Pulay_stress;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PULAY_STRESS (kB) = " << data.Pulay_stress << endl;
    data.energy_cell=xOUT.energy_cell;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENERGY total E0 (eV) = " << data.energy_cell << endl;
    data.energy_atom=xOUT.energy_atom;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENERGY per atom E0/N (eV) = " << data.energy_atom << endl;
    data.enthalpy_cell=xOUT.enthalpy_cell;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY total E0 (eV) = " << data.enthalpy_cell << endl;
    data.enthalpy_atom=xOUT.enthalpy_atom;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY per atom E0/N (eV) = " << data.enthalpy_atom << "   " << directory_LIB << endl;
    data.eentropy_cell=xOUT.eentropy_cell;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " E-ENTROPY total E0 (eV) = " << data.eentropy_cell << endl;
    data.eentropy_atom=xOUT.eentropy_atom;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " E-ENTROPY per atom E0/N (eV) = " << data.eentropy_atom << endl;
    data.PV_cell=xOUT.PV_cell;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PV total E0 (eV) = " << data.PV_cell << endl;
    data.PV_atom=xOUT.PV_atom;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PV per atom E0/N (eV) = " << data.PV_atom << endl;
    data.spin_cell=xOUT.mag_cell;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPIN mag (\\mu) = " << data.spin_cell << endl;
    data.spin_atom=xOUT.mag_atom;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPIN per atom mag/N (\\mu) = " << data.spin_atom << "   " << directory_LIB << endl;
    //CO20180130 START
    //moving FROM magnetic loop so we keep spin/cell, spin/atom, and spinD all together
    data.spinD="";
    data.vspinD.clear();
    if(xOUT.vmag.size()) {
      for(uint i=0;i<(uint) xOUT.vmag.size();i++) {
        data.spinD+=aurostd::utype2string<double>(xOUT.vmag.at(i),5)+(i<xOUT.vmag.size()-1?",":"");
        data.vspinD.push_back(xOUT.vmag.at(i));
      }
    } else {
      for(uint i=0;i<xOUT.natoms;i++) {  //use xOUT.natoms as there can be a primitivization between relax and static
        data.spinD+=aurostd::utype2string<double>(0)+(i<xOUT.natoms-1?",":"");
        data.vspinD.push_back(0.0);
      }
    }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPIND (\\mu) = " << data.spinD << "   " << directory_LIB << endl;
    //CO20180130 STOP
    data.energy_cutoff=xOUT.ENCUT;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENERGY_CUTOFF (eV) = " << data.energy_cutoff << endl;
    data.delta_electronic_energy_convergence=xOUT.total_energy_change;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DELTA_ELECTRONIC_ENERGY_CONVERGENCE (eV) = " << data.delta_electronic_energy_convergence << endl; // CORMAC
    data.delta_electronic_energy_threshold=xOUT.EDIFF;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DELTA_ELECTRONIC_ENERGY_THRESHOLD (eV) = " << data.delta_electronic_energy_threshold << endl; // CORMAC
    data.nkpoints=kpoints.nkpoints;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NKPOINTS (from KPOINTS) = " << data.nkpoints << endl;
    data.kpoints_nnn_relax=kpoints.nnn_kpoints;
    data.kpoints=aurostd::utype2string(kpoints.nnn_kpoints[1])+","+aurostd::utype2string(kpoints.nnn_kpoints[2])+","+aurostd::utype2string(kpoints.nnn_kpoints[3]);
    data.nkpoints_irreducible=xOUT.nkpoints_irreducible;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NKPOINTS_IRREDUCIBLE (from OUTCAR) = " << data.nkpoints_irreducible << endl;
    data.kppra=data.natoms*kpoints.nkpoints;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NKPPRA = " << data.kppra << endl;
    data_vforces.clear(); for(uint i=0;i<xOUT.vforces.size();i++) { data_vforces.push_back(xOUT.vforces.at(i)); }
    data.vforces=data_vforces;
    data_vpositions_cartesian.clear(); for(uint i=0;i<xOUT.vpositions_cartesian.size();i++)  data_vpositions_cartesian.push_back(xOUT.vpositions_cartesian.at(i));
    data.vpositions_cartesian=data_vpositions_cartesian;
    int precfp=8; //DX20190320 - changed from uint to int, otherwise breaks
    data.forces="";
    data.positions_cartesian="";
    for(uint i=0;i<(uint) data.natoms;i++) {
      data.positions_cartesian+=aurostd::utype2string(data_vpositions_cartesian.at(i)[1],precfp)+","+aurostd::utype2string(data_vpositions_cartesian.at(i)[2],precfp)+","+aurostd::utype2string(data_vpositions_cartesian.at(i)[3],precfp);
      if(i<data.natoms-1) data.positions_cartesian+=";";
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " POSITIONS_CARTESIAN = " << data_vpositions_cartesian.at(i)[1] << "," << data_vpositions_cartesian.at(i)[2] << "," << data_vpositions_cartesian.at(i)[3] << "," << endl;
    }
    for(uint i=0;i<(uint) data.natoms;i++) {
      data.forces+=aurostd::utype2string(data_vforces.at(i)[1],precfp)+","+aurostd::utype2string(data_vforces.at(i)[2],precfp)+","+aurostd::utype2string(data_vforces.at(i)[3],precfp);
      if(i<data.natoms-1) data.forces+=";";
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " FORCES(eV/Angst) = " << data_vforces.at(i)[1] << "," << data_vforces.at(i)[2] << "," << data_vforces.at(i)[3] << "," << endl;
    }

    // LOAD SPECIES
    if(LDEBUG) cerr << soliloquy << ": " << data.nspecies << " " << str_relax.species.size() << " " << str_relax.species_pp.size() << " " << str_relax.species_pp_type.size() << " " << str_relax.species_pp_version.size() << " " << str_relax.species_pp_ZVAL.size() << " " << str_relax.species_volume.size() << " " << str_relax.species_mass.size() << endl;

    str_relax.species.clear();str_relax.species_pp.clear();str_relax.species_pp_type.clear();str_relax.species_pp_version.clear();str_relax.species_pp_ZVAL.clear();
    str_relax.species_pp_vLDAU.clear();str_relax.species_volume.clear();str_relax.species_mass.clear();

    // try OUTCARs
    data.vspecies.clear();
    data.vspecies_pp.clear();
    data.vspecies_pp_version.clear();
    data.vspecies_pp_ZVAL.clear();
    data.vspecies_pp_AUID.clear();
    data.METAGGA="";

    for(uint itry=1;itry<=5&&str_relax.species.size()==0;itry++) {  // ONLY OUTCAR as POTCARS ARE OBSOLETE
      for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
        string stry="";
        if(itry==1) stry="OUTCAR.relax1"+XHOST.vext.at(iext);
        if(itry==2) stry="OUTCAR.relax2"+XHOST.vext.at(iext);
        if(itry==3) stry="OUTCAR.relax3"+XHOST.vext.at(iext);
        if(itry==4) stry="OUTCAR.static"+XHOST.vext.at(iext);
        if(itry==5) stry="OUTCAR.bands"+XHOST.vext.at(iext);
        fileE_LIB=directory_LIB+"/"+stry;
        if(aurostd::FileExist(fileE_LIB)) {
          if(LDEBUG) cerr << soliloquy << ": fileE_LIB=" << fileE_LIB << endl;
          // [OBSOLETE]	stringstream stream_outcar(aurostd::execute2string("bzcat"+" \""+fileE_LIB+"\""));
          // [OBSOLETE]	xOUT.GetProperties(stream_outcar);
          xOUT.GetPropertiesFile(fileE_LIB);
          // DEBUG	for(uint i=0;i<xOUT.species.size();i++) cerr << XPID << "xOUT.species.at(i)=" << xOUT.species.at(i) << endl;
          str_relax.species.clear(); for(uint i=0;i<xOUT.species.size();i++) { str_relax.species.push_back(xOUT.species.at(i));data.vspecies.push_back(xOUT.species.at(i)); } // for aflowlib_libraries.cpp
          if(LDEBUG) cerr << soliloquy << ": xOUT.species.size()=" << xOUT.species.size() << endl;
          str_relax.species_pp.clear(); for(uint i=0;i<xOUT.species_pp.size();i++) { str_relax.species_pp.push_back(xOUT.species_pp.at(i));data.vspecies_pp.push_back(xOUT.species_pp.at(i)); } // for aflowlib_libraries.cpp
          str_relax.species_pp_type.clear(); for(uint i=0;i<xOUT.species_pp_type.size();i++) { str_relax.species_pp_type.push_back(xOUT.species_pp_type.at(i));/*data.vspecies_pp_type.push_back(xOUT.vspecies_pp_type.at(i));*/} // for aflowlib_libraries.cpp
          str_relax.species_pp_version.clear(); for(uint i=0;i<xOUT.species_pp_version.size();i++) { str_relax.species_pp_version.push_back(xOUT.species_pp_version.at(i));data.vspecies_pp_version.push_back(xOUT.species_pp_version.at(i)); } // for aflowlib_libraries.cpp
          str_relax.species_pp_ZVAL.clear(); for(uint i=0;i<xOUT.vZVAL.size();i++) { str_relax.species_pp_ZVAL.push_back(xOUT.vZVAL.at(i));data.vspecies_pp_ZVAL.push_back(xOUT.vZVAL.at(i)); } // for aflowlib_libraries.cpp
          data.vspecies_pp_AUID.clear(); for(uint i=0;i<xOUT.species_pp_AUID.size();i++) { data.vspecies_pp_AUID.push_back(xOUT.species_pp_AUID.at(i)); } // for aflowlib_libraries.cpp
          data.dft_type=xOUT.pp_type;
          data.vdft_type.push_back(xOUT.pp_type);  //CO, this is technically a vector (RESTAPI paper)
          str_relax.species_pp_vLDAU.clear(); for(uint i=0;i<xOUT.species_pp_vLDAU.size();i++) str_relax.species_pp_vLDAU.push_back(xOUT.species_pp_vLDAU.at(i));  // for aflowlib_libraries.cpp
          data.ldau_TLUJ=xOUT.string_LDAU;	  	  
          data.METAGGA=xOUT.METAGGA;	  	  

          //ME20190124 BEGIN - Store LDAU information individually
          // Note that the vector here has the species in the columns, not the
          // rows because this is closer to the format in the out and json files.
          data.vLDAU.resize(4);
          for (uint i = 0; i < xOUT.species_pp_vLDAU.size(); i++) {
            for (int j = 0; j < 4; j++) {
              data.vLDAU[j].push_back(xOUT.species_pp_vLDAU[i][j]);
            }
          }
          //ME20190124 END
          if(AFLOWLIB_VERBOSE && data.ldau_TLUJ.size()) cout << MESSAGE << " LDAU_string=" << data.ldau_TLUJ << endl;
        }
      }
    }

    if(str_relax.species.size()==0) { cerr << XPID << MESSAGE << " ERROR - OUTCAR/POTCAR not FOUND in " << directory_LIB << endl;exit(0); }

    for(uint j=0;j<str_relax.species.size();j++) {
      str_relax.species_volume.push_back(GetAtomVolume(str_relax.species.at(j)));
      str_relax.species_mass.push_back(GetAtomMass(str_relax.species.at(j)));
    }

    if(LDEBUG) cerr << soliloquy << ": " << data.nspecies << " " << str_relax.species.size() << " " << str_relax.species_pp.size() << " " << str_relax.species_pp_type.size() << " " << str_relax.species_pp_version.size() << " " << str_relax.species_pp_ZVAL.size() << " " << str_relax.species_pp_vLDAU.size() << " " << str_relax.species_volume.size() << " " << str_relax.species_mass.size() << endl;
    if(data.nspecies!=str_relax.species.size()) { cerr << XPID << MESSAGE << " ERROR [1] - data.nspecies[" << data.nspecies << "]!=str_relax.species.size()[" << str_relax.species.size() << "]" << endl << str_relax << endl;;exit(0); }
    if(data.nspecies!=str_relax.species_pp.size()) { cerr << XPID << MESSAGE << " ERROR [2] - data.nspecies[" << data.nspecies << "]!=str_relax.species_pp.size()[" << str_relax.species_pp.size() << "]" << endl;exit(0); }
    if(data.nspecies!=str_relax.species_pp_type.size()) { cerr << XPID << MESSAGE << " ERROR [3] - data.nspecies[" << data.nspecies << "]!=str_relax.species_pp_type.size()[" << str_relax.species_pp_type.size() << "]" << endl;exit(0); }
    if(data.nspecies!=str_relax.species_pp_version.size()) { cerr << XPID << MESSAGE << " ERROR [4] - data.nspecies[" << data.nspecies << "]!=str_relax.species_pp_version.size()[" << str_relax.species_pp_version.size() << "]" << endl;exit(0); }
    if(data.nspecies!=str_relax.species_pp_ZVAL.size()) { cerr << XPID << MESSAGE << " ERROR [5] - data.nspecies[" << data.nspecies << "]!=str_relax.species_pp_ZVAL.size()[" << str_relax.species_pp_ZVAL.size() << "]" << endl;exit(0); }
    if(data.nspecies!=data.vspecies_pp_AUID.size()) { cerr << XPID << MESSAGE << " ERROR [5] - data.nspecies[" << data.nspecies << "]!=data.vspecies_pp_AUID.size()[" << data.vspecies_pp_AUID.size() << "]" << endl;exit(0); }
    if(data.nspecies!=str_relax.species_volume.size()) { cerr << XPID << MESSAGE << " ERROR [6] - data.nspecies[" << data.nspecies << "]!=str_relax.species_volume.size()[" << str_relax.species_volume.size() << "]" << endl;exit(0); }
    if(data.nspecies!=str_relax.species_mass.size()) { cerr << XPID << MESSAGE << " ERROR [7] - data.nspecies[" << data.nspecies << "]!=str_relax.species_mass.size()[" << str_relax.species_mass.size() << "]" << endl;exit(0); }

    data.compound="";
    data.composition="";
    data.vcomposition.clear();
    data.density=0.0;
    data.stoichiometry="";
    //data.stoich=""; //CO20171026 - we do this here now, mostly obsolete, keep for legacy
    data.species="";
    data.species_pp="";
    data.species_pp_version="";
    data.species_pp_ZVAL="";
    data.species_pp_AUID="";
    data.valence_cell_iupac=0.0;
    data.valence_cell_std=0.0;
    for(uint i=0;i<(uint) data.nspecies;i++) {
      data.compound+=str_relax.species.at(i)+aurostd::utype2string(str_relax.num_each_type.at(i)); // if(i<data.nspecies-1) data.compound+=",";
      data.composition+=aurostd::utype2string(str_relax.num_each_type.at(i)); if(i<data.nspecies-1) data.composition+=",";
      data.vcomposition.push_back(str_relax.num_each_type.at(i));
      data.density+=(double) str_relax.num_each_type.at(i)*GetAtomMass(str_relax.species.at(i));
      data.stoichiometry+=aurostd::utype2string(data_vcomposition.at(i),9); if(i<data.nspecies-1) data.stoichiometry+=",";
      //data.stoich+=aurostd::utype2string(data_vcomposition.at(i),4); if(i<data.nspecies-1) data.stoichiometry+="   "; //mimic old BS format, 4 digits of accuracy and 3 spaces between, stoich=0.5000   0.1667   0.3333
      data.species+=str_relax.species.at(i);if(i<data.nspecies-1) data.species+=",";
      data.species_pp+=str_relax.species_pp.at(i);if(i<data.nspecies-1) data.species_pp+=",";
      // [UNUSED]    data.species_pp_type+=str_relax.species_pp_type.at(i);if(i<data.nspecies-1) data.species_pp_type+=",";
      data.species_pp_version+=str_relax.species_pp_version.at(i);if(i<data.nspecies-1) data.species_pp_version+=",";
      data.species_pp_ZVAL+=aurostd::utype2string(str_relax.species_pp_ZVAL.at(i));if(i<data.nspecies-1) data.species_pp_ZVAL+=",";
      data.species_pp_AUID+=data.vspecies_pp_AUID.at(i);if(i<data.nspecies-1) data.species_pp_AUID+=",";
      //  cerr << XPID << "SIMPLE=" << str_relax.species.at(i) << endl;
      data.valence_cell_iupac+=str_relax.num_each_type.at(i)*GetAtomValenceIupac(str_relax.species.at(i));
      data.valence_cell_std+=str_relax.num_each_type.at(i)*GetAtomValenceStd(str_relax.species.at(i));
    }

    // density
    data.density/=data.volume_cell;
    data.density*=1000.0; // grams instead of kilos
    data.density*=1e8*1e8*1e8; // cm^3 instead of A^3
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DENSITY (grams/cm^3) = " << data.density << endl;

    //DX20190124 - original density info - START
    data.density_orig=0.0;
    for(uint i=0;i<(uint) data.nspecies;i++) {
      data.density_orig+=(double) str_orig.num_each_type.at(i)*GetAtomMass(str_orig.species.at(i));
    }
    // density
    data.density_orig/=data.volume_cell_orig;
    data.density_orig*=1000.0; // grams instead of kilos
    data.density_orig*=1e8*1e8*1e8; // cm^3 instead of A^3
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DENSITY_ORIG (grams/cm^3) = " << data.density_orig << endl;
    //DX20190124 - original density info - START

    // scintillation_attenuation_length
    data.scintillation_attenuation_length=0.0;
    data.scintillation_attenuation_length=GetCompoundAttenuationLength(str_relax.species,str_relax.num_each_type,(double) data.density);
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SCINTILLATION_ATTENUATION_LENGTH (cm) = " << data.scintillation_attenuation_length << endl;

    // [UNUSED] if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PSEUDOPOTENTIAL species_pp_type = " << data.species_pp_type << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PSEUDOPOTENTIAL dft_type=" << data.dft_type << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PSEUDOPOTENTIAL species_pp_version = " << data.species_pp_version << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PSEUDOPOTENTIAL species_pp_ZVAL = " << data.species_pp_ZVAL << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PSEUDOPOTENTIAL species_pp_AUID = " << data.species_pp_AUID << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PSEUDOPOTENTIAL METAGGA = [" << data.METAGGA << "]" << endl;
    // reference
    bool FORMATION_CALC=TRUE;    
    bool isLDAUcalc=FALSE;
    if(str_relax.species_pp_vLDAU.at(0).size()>0) isLDAUcalc=TRUE;
    if(LDEBUG){cerr << soliloquy << " isLDAUcalc=" << isLDAUcalc << endl;}
    
    //CO20200624 START - CCE
    string functional="";
    if(str_relax.species_pp_version.size()>0){
      //CO20200624 - don't use else if's, we might prefer that which comes later
      if(str_relax.species_pp_version[0].find("PBE")!=string::npos){functional="PBE";}
      if(str_relax.species_pp_version[0].find("LDA")!=string::npos){functional="LDA";}
      if(str_relax.species_pp_version[0].find("SCAN")!=string::npos){functional="SCAN";}
      //double check all of the pp match
      for(uint i=0;i<str_relax.species_pp_version.size()&&!functional.empty();i++){
        if(functional=="PBE" && str_relax.species_pp_version[i].find("PBE")==string::npos){
          pflow::logger(_AFLOW_FILE_NAME_,soliloquy,"Mismatch in species_pp_version found ("+functional+"), not calculating CCE correction",_LOGGER_WARNING_);
          functional="";
        }
        if(functional=="LDA" && str_relax.species_pp_version[i].find("LDA")!=string::npos){
          pflow::logger(_AFLOW_FILE_NAME_,soliloquy,"Mismatch in species_pp_version found ("+functional+"), not calculating CCE correction",_LOGGER_WARNING_);
          functional="";
        }
        if(functional=="SCAN" && str_relax.species_pp_version[i].find("SCAN")!=string::npos){
          pflow::logger(_AFLOW_FILE_NAME_,soliloquy,"Mismatch in species_pp_version found ("+functional+"), not calculating CCE correction",_LOGGER_WARNING_);
          functional="";
        }
      }
      if(data.catalog=="ICSD" && isLDAUcalc==TRUE && functional=="PBE"){functional="PBE+U:ICSD";}
    }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " FUNCTIONAL=" << functional << endl;
    //CO20200624 STOP - CCE

    // LDAU ?
    if(FORMATION_CALC) { //CO20200624 - we can correct PBE now with CCE
      //[CO20200624 - we can correct PBE now with CCE]if(isLDAUcalc) FORMATION_CALC=FALSE; // no formation for LDAU
      if(functional.empty()) FORMATION_CALC=FALSE;  //CO20200624
    }

    // PP AVAILABLE ?
    if(FORMATION_CALC) {
      for(uint i=0;i<(uint) data.nspecies;i++) {
        if(!xPOTCAR_EnthalpyReference_AUID(data.vspecies_pp_AUID.at(i),data.METAGGA)) { 
          FORMATION_CALC=FALSE; // WORKS WITH SCANN TOO
          if(AFLOWLIB_VERBOSE) cout << MESSAGE << " REFERENCE NOT AVAILABLE species_pp=" << str_relax.species_pp_version.at(i) << "  species_pp_AUID=" << data.vspecies_pp_AUID.at(i)
            << (data.METAGGA.empty()?string(""):string("  METAGGA=["+data.METAGGA+"]")) << "   Href=nothing" << endl;
        }

      }
    }

    // OPERATE FORMATION
    if(FORMATION_CALC) { // no LDAU yet
      if(LDEBUG) cerr << soliloquy << " [FCALC=1]" << endl;
      vector<double> venthalpy_atom_ref;
      for(uint i=0;i<(uint) data.nspecies;i++) {
        //      string pseudopotential,string type,vector<double> LDAU
        double enthalpy_atom_ref=data.enthalpy_atom; // if there is 1 then there is only one
        string aus_gs_structure;
        double aus_gs_atom,aus_volume_atom,aus_spin_atom;
        // NEW STYLE
        // bool found=
        xPOTCAR_EnthalpyReference_AUID(data.vspecies_pp_AUID.at(i),data.METAGGA,aus_gs_structure,aus_gs_atom,aus_volume_atom,aus_spin_atom);
        enthalpy_atom_ref=aus_gs_atom;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " REFERENCE species=" << str_relax.species.at(i) << " species_pp_AUID=" << data.vspecies_pp_AUID.at(i) << (data.METAGGA.empty()?string(""):string("  METAGGA=["+data.METAGGA+"]")) << "   Href=" << enthalpy_atom_ref << endl;
        venthalpy_atom_ref.push_back(enthalpy_atom_ref);
      }

      if(LDEBUG) cerr << soliloquy << " [FCALC=2]" << endl;
      if(LDEBUG) cerr << soliloquy << " [FCALC=2] data.nspecies=" << data.nspecies << endl;
      if(LDEBUG) cerr << soliloquy << " [FCALC=2] venthalpy_atom_ref.size()=" << venthalpy_atom_ref.size() << endl;

      // calculation of REF
      //     cerr << XPID << data.enthalpy << endl;
      data.enthalpy_formation_cell=data.enthalpy_cell;
      data.enthalpy_formation_atom=data.enthalpy_atom;

      for(uint i=0;i<(uint) data.nspecies;i++) data.enthalpy_formation_cell=data.enthalpy_formation_cell-(double(venthalpy_atom_ref.at(i)*str_relax.num_each_type.at(i)));
      //   for(uint i=0;i<(uint) data.nspecies;i++) data.enthalpy_formation_atom=data.enthalpy_formation_atom-venthalpy_atom_ref.at(i)*double(str_relax.num_each_type.at(i))/double(str_relax.atoms.size());
      data.enthalpy_formation_atom=data.enthalpy_formation_cell/double(str_relax.atoms.size());
      
      //CO20200624 START - adding cce variants
      if(!functional.empty()){
        bool found_correctable=false;
        if(aurostd::WithinList(data.vspecies,"O")) found_correctable=true;
        //[CO20200624 - not yet]if(aurostd::WithinList(data.vspecies,"N")) found_correctable=true;
        if(found_correctable){
          data.enthalpy_formation_cce_300K_cell=data.enthalpy_formation_cce_0K_cell=data.enthalpy_formation_cell;
          vector<double> enthalpy_formation_cell_corrections_cce=cce::calculate_corrections(str_relax,functional);
          if(enthalpy_formation_cell_corrections_cce.size()==2){  //the first is at 300K, the second at 0K
            data.enthalpy_formation_cce_300K_cell-=enthalpy_formation_cell_corrections_cce[0];
            data.enthalpy_formation_cce_0K_cell-=enthalpy_formation_cell_corrections_cce[1];
            data.enthalpy_formation_cce_300K_atom=data.enthalpy_formation_cce_300K_cell/double(str_relax.atoms.size());
            data.enthalpy_formation_cce_0K_atom=data.enthalpy_formation_cce_0K_cell/double(str_relax.atoms.size());
          }
        }
      }
      if(LDEBUG){
        cerr << soliloquy << " data.enthalpy_formation_cce_300K_cell=" << data.enthalpy_formation_cce_300K_cell << endl;
        cerr << soliloquy << " data.enthalpy_formation_cce_0K_cell=" << data.enthalpy_formation_cce_0K_cell << endl;
        cerr << soliloquy << " data.enthalpy_formation_cce_300K_atom=" << data.enthalpy_formation_cce_300K_atom << endl;
        cerr << soliloquy << " data.enthalpy_formation_cce_0K_atom=" << data.enthalpy_formation_cce_0K_atom << endl;
      }
      //CO20200624 STOP - adding cce variants

      if(LDEBUG) cerr << soliloquy << " [FCALC=3]" << endl;

      data.entropic_temperature=0;
      if(data_vcomposition.size()>1) {
        for(uint i=0;i<(uint) data_vcomposition.size();i++)
          if(data_vcomposition.at(i)>_EPSILON_COMPOSITION_ && data_vcomposition.at(i)<1-_EPSILON_COMPOSITION_)
            data.entropic_temperature+=data_vcomposition.at(i)*logl(data_vcomposition.at(i));
        data.entropic_temperature=data.enthalpy_formation_atom/(data.entropic_temperature*KBOLTZEV);
      }
      // cerr << XPID << data.enthalpy_formation_cell << endl << data.enthalpy_formation_cell/str_relax.atoms.size() << endl << data.enthalpy_formation_atom << endl;
      if(LDEBUG) cerr << soliloquy << " [FCALC=4]" << endl;
    }

    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " VALENCE_IUPAC = " << data.valence_cell_iupac << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " VALENCE_STD = " << data.valence_cell_std << endl;

    //   aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "energyd=" << data_dE;
    //   aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "energyd_atom=" << data_dEN;

    if(FORMATION_CALC==TRUE) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY FORMATION total E0 (eV) = " << data.enthalpy_formation_cell << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY FORMATION per atom E0/N (eV) = " << data.enthalpy_formation_atom << "   " << directory_LIB << endl;
      //CO20200624 START - CCE
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY FORMATION CCE total E(300K) (eV) = " << data.enthalpy_formation_cce_300K_cell << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY FORMATION CCE per atom E(300K)/N (eV) = " << data.enthalpy_formation_cce_300K_atom << "   " << directory_LIB << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY FORMATION CCE total E(0K) (eV) = " << data.enthalpy_formation_cce_0K_cell << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTHALPY FORMATION CCE per atom E(0K)/N (eV) = " << data.enthalpy_formation_cce_0K_atom << "   " << directory_LIB << endl;
      //CO20200624 END - CCE
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTROPIC_TEMPERATURE (eV) = " << data.entropic_temperature*KBOLTZEV << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ENTROPIC_TEMPERATURE (K) = " << data.entropic_temperature << "   " << directory_LIB << endl;
    }
    // [FIX]  if(flag_ENERGY1) aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy1_cell=" << data1_energy_cell;
    // [FIX] if(flag_ENERGY1) aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "energy1_atom=" << data1_energy_atom;
    // DONE WITH THERMO

    // do the TIMING
    if(LDEBUG) cerr << soliloquy << " [9]" << endl;
    data.calculation_cores=1;
    data.calculation_time=0.0;
    data.calculation_memory=0.0;
    if(flag_TIMING) {  // OUTCAR.relax.EXT
      for(uint i=1;i<=relax_max+2;i++) {
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          fileE_LIB=directory_LIB+"/OUTCAR.relax"+aurostd::utype2string(i)+XHOST.vext.at(iext);
          if(i==relax_max+1) fileE_LIB=directory_LIB+"/OUTCAR.static"+XHOST.vext.at(iext);  // do the static
          if(i==relax_max+2) fileE_LIB=directory_LIB+"/OUTCAR.bands"+XHOST.vext.at(iext);  // do the bands
          // cerr << XPID << fileE_LIB << endl;
          if(aurostd::FileExist(fileE_LIB)) {
            xOUTCAR outcar_tmp; // cerr << XPID << fileE_LIB << endl;
            outcar_tmp.GetPropertiesFile(fileE_LIB); // OK
            data.calculation_cores=aurostd::max((int) data.calculation_cores,(int) outcar_tmp.calculation_cores);
            //	data.calculation_time+=double(data.calculation_cores)*outcar_tmp.calculation_time;
            data.calculation_time+=outcar_tmp.calculation_time;  // will multiply after
            data.calculation_memory=aurostd::max(data.calculation_memory,outcar_tmp.calculation_memory);
            xOUTCAR outcar;
            if(i==relax_max+2) { outcar.GetPropertiesFile(directory_LIB+"/OUTCAR.bands"+XHOST.vext.at(iext)); } // cerr << XPID << "xOUTCAR.Efermi=" << outcar.Efermi << endl;  //CO20200106 - patching for auto-indenting
            xEIGENVAL eigenval;
            if(i==relax_max+2) { eigenval.GetPropertiesFile(directory_LIB+"/EIGENVAL.bands"+XHOST.vext.at(iext)); }
            xDOSCAR doscar;
            if(i==relax_max+2) { doscar.GetPropertiesFile(directory_LIB+"/DOSCAR.bands"+XHOST.vext.at(iext)); } // cerr << XPID << "xDOSCAR.Efermi=" << doscar.Efermi << " " << outcar.Efermi-doscar.Efermi << endl; //CO20200106 - patching for auto-indenting
          }
        }
      }

      if(_APENNSY_STYLE_OLD_) aurostd::string2file(aurostd::utype2string((long(100*data.calculation_time))/100)+"\n",directory_RAW+"/"+DEFAULT_FILE_TIME_OUT);
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CALCULATION time (sec) = " << data.calculation_time << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CALCULATION mem (MB) = " << data.calculation_memory << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CALCULATION cores = " << data.calculation_cores << endl;
      //	exit(0);
    }
    // do the SG
    if(LDEBUG) cerr << soliloquy << " [10]" << endl;
    if(flag_SG1 || flag_SG2) {  // POSCAR.orig CONTCAR.relax1 CONTCAR.relax
#ifdef USE_PLATON_SG
      string space_group_calculator_sg1_cmd=" | aflow --platonSG=3.0,0.5,0.5,0.5";
      string space_group_calculator_sg2_cmd=" | aflow --platonSG=1.5,0.25,0.25,0.25";
      //DX20190319 - moved to end of external call lines [OBSOLETE] #endif

      //DX20190319 - moved further down [OBSOLETE] #ifdef USE_AFLOW_SG
      //DX20190319 - moved further down [OBSOLETE]      //DX+CO START
      //DX20190319 - moved further down [OBSOLETE]      //DX [OBSOLETE] string space_group_calculator_sg1_cmd=" | aflow --aflowSG=0.001";
      //DX20190319 - moved further down [OBSOLETE]      //DX [OBSOLETE] string space_group_calculator_sg2_cmd=" | aflow --aflowSG=0.00075";
      //DX20190319 - moved further down [OBSOLETE]      string space_group_calculator_sg1_cmd=" | aflow --aflowSG=loose";
      //DX20190319 - moved further down [OBSOLETE]      string space_group_calculator_sg2_cmd=" | aflow --aflowSG=tight";
      //DX20190319 - moved further down [OBSOLETE]      //DX+CO END
      //DX20190319 - moved further down [OBSOLETE] #endif
      string ilattice_cmd=" | aflow --ilattice 2.0";
      if(flag_SG1) {
        stringstream ssfile;
        cout << MESSAGE << " Space Group analyzer: " << space_group_calculator_sg1_cmd << endl;
        ssfile << "Space Group analyzer: RELAX: " << space_group_calculator_sg1_cmd << endl;
        ssfile << directory_RAW << endl;
        // I run outside so a segfault+core does not block the aflow
        // INSIDE no more inside to avoid BOMBING
        // str_orig.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,3.0,0.5,0.5,0.5);data_sg1_pre=str_orig.spacegroup;     //  aflow --platonSG 3.0 0.5 0.5 0.5
        // str_relax1.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,3.0,0.5,0.5,0.5);data_sg1_mid=str_relax1.spacegroup; //  aflow --platonSG 3.0 0.5 0.5 0.5
        // str_relax.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,3.0,0.5,0.5,0.5);data_sg1_post=str_relax.spacegroup;  //  aflow --platonSG 3.0 0.5 0.5 0.5
        // OUTSIDE
        // sg1_pre
        data_sg1_pre=aurostd::execute2string("cat \""+directory_RAW+"/POSCAR.orig"+"\""+space_group_calculator_sg1_cmd);
        aurostd::string2tokens(data_sg1_pre,tokens,"#");if(tokens.size()!=2) { data_sg1_pre=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_pre=NOSG; } else { aurostd::StringSubst(data_sg1_pre,"\n",""); }}
        if(data_sg1_pre==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg1_pre" << endl;data_sg1_pre=aurostd::execute2string("cat \""+directory_RAW+"/POSCAR.orig"+"\""+ilattice_cmd+space_group_calculator_sg1_cmd); }
        aurostd::string2tokens(data_sg1_pre,tokens,"#");if(tokens.size()!=2) { data_sg1_pre=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_pre=NOSG; } else { aurostd::StringSubst(data_sg1_pre,"\n",""); }}
        if(aurostd::substring2bool(data_sg1_pre,"SymbolnotKnown")) data_sg1_pre=NOSG;  // give up
        // sg1_mid
        data_sg1_mid=aurostd::execute2string("cat \""+directory_RAW+"/CONTCAR.relax1"+"\""+space_group_calculator_sg1_cmd);
        aurostd::string2tokens(data_sg1_mid,tokens,"#");if(tokens.size()!=2) { data_sg1_mid=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_mid=NOSG; } else { aurostd::StringSubst(data_sg1_mid,"\n",""); }}
        if(data_sg1_mid==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg1_mid" << endl;data_sg1_mid=aurostd::execute2string("cat \""+directory_RAW+"/CONTCAR.relax1"+"\""+ilattice_cmd+space_group_calculator_sg1_cmd); }
        aurostd::string2tokens(data_sg1_mid,tokens,"#");if(tokens.size()!=2) { data_sg1_mid=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_mid=NOSG; } else { aurostd::StringSubst(data_sg1_mid,"\n",""); }}
        if(aurostd::substring2bool(data_sg1_mid,"SymbolnotKnown")) data_sg1_mid=NOSG;  // give up
        // sg1_mid
        data_sg1_post=aurostd::execute2string("cat \""+directory_RAW+"/CONTCAR.relax"+"\""+space_group_calculator_sg1_cmd);
        aurostd::string2tokens(data_sg1_post,tokens,"#");if(tokens.size()!=2) { data_sg1_post=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_post=NOSG; } else { aurostd::StringSubst(data_sg1_post,"\n",""); }}
        if(data_sg1_post==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg1_post" << endl;data_sg1_post=aurostd::execute2string("cat \""+directory_RAW+"/CONTCAR.relax"+"\""+ilattice_cmd+space_group_calculator_sg1_cmd); }
        aurostd::string2tokens(data_sg1_post,tokens,"#");if(tokens.size()!=2) { data_sg1_post=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg1_post=NOSG; } else { aurostd::StringSubst(data_sg1_post,"\n",""); }}
        if(aurostd::substring2bool(data_sg1_post,"SymbolnotKnown")) data_sg1_post=NOSG;  // give up
        // DONE
        ssfile << "PRE  " << data_sg1_pre << endl;
        ssfile << "MID  " << data_sg1_mid << endl;
        ssfile << "POST " << data_sg1_post << endl;
        //      aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_SPACEGROUP1_OUT);
        data.sg=data_sg1_pre+","+data_sg1_mid+","+data_sg1_post;
        data.vsg.clear();data.vsg.push_back(data_sg1_pre);data.vsg.push_back(data_sg1_mid);data.vsg.push_back(data_sg1_post); //CO20171202
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP1 = " << data.sg << endl;
      }
      if(flag_SG2) {
        stringstream ssfile;
        cout << MESSAGE << " Space Group analyzer: " << space_group_calculator_sg2_cmd << endl;
        ssfile << "Space Group analyzer: RELAX: " << space_group_calculator_sg2_cmd << endl;
        ssfile << directory_RAW << endl;
        // I run outside so a segfault+core does not block the aflow
        // INSIDE
        // str_orig.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,1.5,0.25,0.25,0.25);data_sg2_pre=str_orig.spacegroup;     //  aflow --platonSG 1.5 0.25 0.25 0.25
        // str_relax1.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,1.5,0.25,0.25,0.25);data_sg2_mid=str_relax1.spacegroup; //  aflow --platonSG 1.5 0.25 0.25 0.25
        // str_relax.platon2sg(DEFAULT_PLATON_P_EQUAL,DEFAULT_PLATON_P_EXACT,1.5,0.25,0.25,0.25);data_sg2_post=str_relax.spacegroup;  //  aflow --platonSG 1.5 0.25 0.25 0.25
        // OUTSIDE
        // sg2_pre
        data_sg2_pre=aurostd::execute2string("cat \""+directory_RAW+"/POSCAR.orig"+"\""+space_group_calculator_sg2_cmd);
        aurostd::string2tokens(data_sg2_pre,tokens,"#");if(tokens.size()!=2) { data_sg2_pre=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_pre=NOSG; } else { aurostd::StringSubst(data_sg2_pre,"\n",""); }}
        if(data_sg2_pre==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg2_pre" << endl;data_sg2_pre=aurostd::execute2string("cat \""+directory_RAW+"/POSCAR.orig"+"\""+ilattice_cmd+space_group_calculator_sg2_cmd); }
        aurostd::string2tokens(data_sg2_pre,tokens,"#");if(tokens.size()!=2) { data_sg2_pre=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_pre=NOSG; } else { aurostd::StringSubst(data_sg2_pre,"\n",""); }}
        if(aurostd::substring2bool(data_sg2_pre,"SymbolnotKnown")) data_sg2_pre=NOSG;  // give up
        // sg2_mid
        data_sg2_mid=aurostd::execute2string("cat \""+directory_RAW+"/CONTCAR.relax1"+"\""+space_group_calculator_sg2_cmd);
        aurostd::string2tokens(data_sg2_mid,tokens,"#");if(tokens.size()!=2) { data_sg2_mid=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_mid=NOSG; } else { aurostd::StringSubst(data_sg2_mid,"\n",""); }}
        if(data_sg2_mid==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg2_mid" << endl;data_sg2_mid=aurostd::execute2string("cat \""+directory_RAW+"/CONTCAR.relax1"+"\""+ilattice_cmd+space_group_calculator_sg2_cmd); }
        aurostd::string2tokens(data_sg2_mid,tokens,"#");if(tokens.size()!=2) { data_sg2_mid=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_mid=NOSG; } else { aurostd::StringSubst(data_sg2_mid,"\n",""); }}
        if(aurostd::substring2bool(data_sg2_mid,"SymbolnotKnown")) data_sg2_mid=NOSG;  // give up
        // sg2_mid
        data_sg2_post=aurostd::execute2string("cat \""+directory_RAW+"/CONTCAR.relax"+"\""+space_group_calculator_sg2_cmd);
        aurostd::string2tokens(data_sg2_post,tokens,"#");if(tokens.size()!=2) { data_sg2_post=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_post=NOSG; } else { aurostd::StringSubst(data_sg2_post,"\n",""); }}
        if(data_sg2_post==NOSG) { if(AFLOWLIB_VERBOSE) cout << MESSAGE << " CORRECTING sg2_post" << endl;data_sg2_post=aurostd::execute2string("cat \""+directory_RAW+"/CONTCAR.relax"+"\""+ilattice_cmd+space_group_calculator_sg2_cmd); }
        aurostd::string2tokens(data_sg2_post,tokens,"#");if(tokens.size()!=2) { data_sg2_post=NOSG; } else { if(aurostd::string2utype<uint>(tokens.at(1))==0) { data_sg2_post=NOSG; } else { aurostd::StringSubst(data_sg2_post,"\n",""); }}
        if(aurostd::substring2bool(data_sg2_post,"SymbolnotKnown")) data_sg2_post=NOSG;  // give up
        // DONE
        ssfile << "PRE  " << data_sg2_pre << endl;
        ssfile << "MID  " << data_sg2_mid << endl;
        ssfile << "POST " << data_sg2_post << endl;
        // aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_SPACEGROUP2_OUT);
        data.sg2=data_sg2_pre+","+data_sg2_mid+","+data_sg2_post;
        data.vsg2.clear();data.vsg2.push_back(data_sg2_pre);data.vsg2.push_back(data_sg2_mid);data.vsg2.push_back(data_sg2_post); //CO20171202
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP2 = " << data.sg2 << endl;
      }
#endif //DX20190319 - moved endif to end of external functions

      //DX20190319 - now call aflowSG internally - START
#ifdef USE_AFLOW_SG
      //DX+CO START
      //DX [OBSOLETE] string space_group_calculator_sg1_cmd=" | aflow --aflowSG=0.001";
      //DX [OBSOLETE] string space_group_calculator_sg2_cmd=" | aflow --aflowSG=0.00075";
      //DX20190319 [OBSOLETE] string space_group_calculator_sg1_cmd=" | aflow --aflowSG=loose";
      //DX20190319 [OBSOLETE] string space_group_calculator_sg2_cmd=" | aflow --aflowSG=tight";

      string space_group_calculator_sg1_function = "aflowSG(loose)";
      string space_group_calculator_sg2_function = "aflowSG(tight)";
      bool no_scan = false;

      // pre
      xstructure str_orig(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO);
      double pre_default_tolerance=SYM::defaultTolerance(str_orig);
      // mid
      xstructure str_mid(directory_RAW+"/CONTCAR.relax1",IOAFLOW_AUTO);
      double mid_default_tolerance=SYM::defaultTolerance(str_mid);
      // post
      xstructure str_post(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);
      double post_default_tolerance=SYM::defaultTolerance(str_post);

      if(flag_SG1) {
        stringstream ssfile;
        cout << MESSAGE << " Space Group analyzer: " << space_group_calculator_sg1_function << endl;
        ssfile << "Space Group analyzer: RELAX: " << space_group_calculator_sg1_function << endl;
        ssfile << directory_RAW << endl;

        // sg1_pre
        xstructure str_sg1_pre=str_orig;
        double pre_loose_tolerance = pre_default_tolerance*10.0;
        uint sgroup_sg1_pre = str_sg1_pre.SpaceGroup_ITC(pre_loose_tolerance,no_scan);
        if(sgroup_sg1_pre<1 || sgroup_sg1_pre>230) {data_sg1_pre=NOSG;}
        else { data_sg1_pre=GetSpaceGroupName(sgroup_sg1_pre,str_sg1_pre.directory)+" #"+aurostd::utype2string(sgroup_sg1_pre); }

        // sg1_mid
        xstructure str_sg1_mid=str_mid;
        double mid_loose_tolerance = mid_default_tolerance*10.0;
        uint sgroup_sg1_mid = str_sg1_mid.SpaceGroup_ITC(mid_loose_tolerance,no_scan);
        if(sgroup_sg1_mid<1 || sgroup_sg1_mid>230) {data_sg1_mid=NOSG;}
        else { data_sg1_mid=GetSpaceGroupName(sgroup_sg1_mid,str_sg1_mid.directory)+" #"+aurostd::utype2string(sgroup_sg1_mid); }

        // sg1_post
        xstructure str_sg1_post=str_post;
        double post_loose_tolerance = post_default_tolerance*10.0;
        uint sgroup_sg1_post = str_sg1_post.SpaceGroup_ITC(post_loose_tolerance,no_scan);
        if(sgroup_sg1_post<1 || sgroup_sg1_post>230) {data_sg1_post=NOSG;}
        else { data_sg1_post=GetSpaceGroupName(sgroup_sg1_post,str_sg1_post.directory)+" #"+aurostd::utype2string(sgroup_sg1_post); }

        // DONE
        ssfile << "PRE  " << data_sg1_pre << endl;
        ssfile << "MID  " << data_sg1_mid << endl;
        ssfile << "POST " << data_sg1_post << endl;
        //      aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_SPACEGROUP1_OUT);
        data.sg=data_sg1_pre+","+data_sg1_mid+","+data_sg1_post;
        data.vsg.clear();data.vsg.push_back(data_sg1_pre);data.vsg.push_back(data_sg1_mid);data.vsg.push_back(data_sg1_post); //CO20171202
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP1 = " << data.sg << endl;
      }
      if(flag_SG2) {
        stringstream ssfile;
        cout << MESSAGE << " Space Group analyzer: " << space_group_calculator_sg2_function << endl;
        ssfile << "Space Group analyzer: RELAX: " << space_group_calculator_sg2_function << endl;
        ssfile << directory_RAW << endl;

        // sg2_pre
        xstructure str_sg2_pre=str_orig;
        double pre_tight_tolerance = pre_default_tolerance;
        uint sgroup_sg2_pre = str_sg2_pre.SpaceGroup_ITC(pre_tight_tolerance,no_scan);
        if(sgroup_sg2_pre<1 || sgroup_sg2_pre>230) {data_sg2_pre=NOSG;}
        else { data_sg2_pre=GetSpaceGroupName(sgroup_sg2_pre,str_sg2_pre.directory)+" #"+aurostd::utype2string(sgroup_sg2_pre); }

        // sg2_mid
        xstructure str_sg2_mid=str_mid;
        double mid_tight_tolerance = mid_default_tolerance;
        uint sgroup_sg2_mid = str_sg2_mid.SpaceGroup_ITC(mid_tight_tolerance,no_scan);
        if(sgroup_sg2_mid<1 || sgroup_sg2_mid>230) {data_sg2_mid=NOSG;}
        else { data_sg2_mid=GetSpaceGroupName(sgroup_sg2_mid,str_sg2_mid.directory)+" #"+aurostd::utype2string(sgroup_sg2_mid); }

        // sg2_post
        xstructure str_sg2_post=str_post;
        double post_tight_tolerance = post_default_tolerance;
        uint sgroup_sg2_post = str_sg2_post.SpaceGroup_ITC(post_tight_tolerance,no_scan);
        if(sgroup_sg2_post<1 || sgroup_sg2_post>230) {data_sg2_post=NOSG;}
        else { data_sg2_post=GetSpaceGroupName(sgroup_sg2_post,str_sg2_post.directory)+" #"+aurostd::utype2string(sgroup_sg2_post); }

        // DONE
        ssfile << "PRE  " << data_sg2_pre << endl;
        ssfile << "MID  " << data_sg2_mid << endl;
        ssfile << "POST " << data_sg2_post << endl;
        // aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_SPACEGROUP2_OUT);
        data.sg2=data_sg2_pre+","+data_sg2_mid+","+data_sg2_post;
        data.vsg2.clear();data.vsg2.push_back(data_sg2_pre);data.vsg2.push_back(data_sg2_mid);data.vsg2.push_back(data_sg2_post); //CO20171202
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP2 = " << data.sg2 << endl;
      }

      //DX20190131 [OBSOLETE] - we should use the self-consistent space group from the edata run; otherwise, Wyckoffs and SG may not match - START
      //DX20190131 [OBSOLETE] aurostd::string2tokens(data_sg2_pre,tokens,"#");tokens.at(tokens.size()-1);
      //DX20190131 [OBSOLETE] data.spacegroup_orig="0"; if(tokens.size()>1) { data.spacegroup_orig=tokens.at(1); }
      //DX20190131 [OBSOLETE] if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP_ORIG = " << data.spacegroup_orig << endl;

      //DX20190131 [OBSOLETE] aurostd::string2tokens(data_sg2_post,tokens,"#");tokens.at(tokens.size()-1);
      //DX20190131 [OBSOLETE] data.spacegroup_relax="0"; if(tokens.size()>1) { data.spacegroup_relax=tokens.at(1); }
      //DX20190131 [OBSOLETE] if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP_RELAX = " << data.spacegroup_relax << endl;
      //DX20190131 [OBSOLETE] - we should use the self-consistent space group from the edata run; otherwise, Wyckoffs and SG may not match - END
#endif
    }
    //DX20190319 - now call aflowSG internally - END
    if(LDEBUG) cerr << soliloquy << " [11]" << endl;
    // VOLDISTParams
    if(flag_VOLDISTPARAMS) {  // CONTCAR.relax
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [VOLDISTParams]" << endl;
      data.vnbondxx=GetNBONDXX(str_relax);  //CO20171024
      data.nbondxx=aurostd::joinWDelimiter(aurostd::vecDouble2vecString(data.vnbondxx,5),",");
      //[CO20200624 - OBSOLETE]if(data.vnbondxx.size())
      //[CO20200624 - OBSOLETE]  for(uint i=0;i<(uint) data.vnbondxx.size();i++)
      //[CO20200624 - OBSOLETE]    data.nbondxx+=aurostd::utype2string<double>(data.vnbondxx.at(i),5)+(i<data.vnbondxx.size()-1?",":"");
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NBONDXX = " << data.nbondxx << endl;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " MIN_DIST = " << SYM::minimumDistance(str_relax) << endl; //CO20171025

      if(aurostd::abs((double) data.vnbondxx.size()-data.nspecies*(data.nspecies+1.0)/2.0)>0.1){
        message << "incompatible data.vnbondxx.size()[" << data.vnbondxx.size() << "]!=data.nspecies*(data.nspecies+1)/2)[" << (data.nspecies*(data.nspecies+1.0)/2.0) << "]";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_RUNTIME_ERROR_);
      }
      //[CO20200624 - OBSOLETE]{ cerr << XPID << MESSAGE << " incompatible data.vnbondxx.size()[" << data.vnbondxx.size() << "]!=_nspecies*(_nspecies+1)/2)[" << (_nspecies*(_nspecies+1.0)/2.0) << "]" << endl;exit(0); }

      //CO20200624 START - mimic stoich from PrintData1(): aflow_pflow_print.cpp, this is really obsolete
      stringstream stoich_ss;stoich_ss.precision(4);
      for(uint it=0;it<str_relax.num_each_type.size();it++) {
        stoich_ss << setw(8) << (double)str_relax.num_each_type[it]/(double)str_relax.atoms.size() << " ";
      }
      data.stoich=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(stoich_ss.str());
      //CO20200624 END - mimic stoich from PrintData1(): aflow_pflow_print.cpp, this is really obsolete

      if(0){  //CO20200624 - OBSOLETE, we don't write this file out anymore
        stringstream ssfile;
        ssfile << "Volume Distance Bonds Analyzer: RELAX" << endl;
        ssfile << directory_RAW << endl;
        string data1=pflow::PrintData1(str_relax,-1.0);
        vector<string> vdata1;aurostd::string2vectorstring(data1,vdata1);
        for(uint i=0;i<vdata1.size();i++) {
          if(aurostd::substring2bool(vdata1.at(i),"Stoich Str1")) {
            data.stoich=aurostd::substring2string(vdata1.at(i),"Stoich Str1: ",FALSE);
            ssfile << "STOICH " << data.stoich << endl;
          }
          if(aurostd::substring2bool(vdata1.at(i),"Vol Str1")) {
            data_v_atom=aurostd::substring2string(vdata1.at(i),"Vol Str1: ",FALSE);
            ssfile << "V_ATOM " << data_v_atom << endl;
          }
          if(aurostd::substring2bool(vdata1.at(i),"Cell Str1")) {
            data_ucelld=aurostd::substring2string(vdata1.at(i),"Cell Str1: ",FALSE);
            ssfile << "UCELLD " << data_ucelld << endl;
          }
        }
        // aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "v_atom=" << data_v_atom;
        // aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "ucelld=" << data_ucelld;
        vector<double> vnbondxx_OLD;
        uint _nspecies=data.nspecies;
        for(uint isp1=0;isp1<_nspecies;isp1++)
          for(uint isp2=isp1;isp2<_nspecies;isp2++)
            for(uint i=0;i<vdata1.size();i++) {
              string string2find,string2search="Pairs between types: "+aurostd::utype2string(isp1)+" "+aurostd::utype2string(isp2);
              tokens.clear();tokens.push_back("");
              if(aurostd::substring2bool(vdata1.at(i),string2search)) {
                string2find=aurostd::substring2string(vdata1.at(i),string2search,FALSE);
                aurostd::string2tokens(string2find,tokens);
                vnbondxx_OLD.push_back(aurostd::string2utype<double>(tokens.at(0))); // will do better but for now it is OK
                //[CO20171024 OBSOLETE]data.vnbondxx.push_back(aurostd::string2utype<double>(tokens.at(0))); // will do better but for now it is OK
              }
            }
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " NBONDXX_OLD = " << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(vnbondxx_OLD,6),',') << endl; //CO20171025

        for(uint isp1=0,inbondxx=0;isp1<_nspecies;isp1++)
          for(uint isp2=isp1;isp2<_nspecies;isp2++)
            ssfile << "BOND_" << char('A'+isp1) << char('A'+isp2) << " " << data.vnbondxx.at(inbondxx++) << " [norm V_ATOM^0.33]" << endl;
        if(flag_LIB2==TRUE && _nspecies==1) { _nspecies++;data.vnbondxx.push_back(0.0);data.vnbondxx.push_back(0.0); } // FIX for purity control and alien generation
        // aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_VOLDISTPARAMS_OUT);
      }
    }
    // PRINT FORCES/POSITIONS_CARTESIAN

    if(LDEBUG) cerr << soliloquy << " [12]" << endl;
    // VOLDISTEvolution
    if(flag_VOLDISTEVOLUTION) {  // CONTCAR.relax
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [VOLDISTEvolution]" << endl;
      stringstream ssfile;
      xstructure str;
      bool VOLDISTEvolution_flag=FALSE;//TRUE;
      for(int istep=1;istep<=3;istep++) {
        if(LDEBUG) cerr << soliloquy << " [AUID=] istep=" << istep << endl;
        if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [start istep] " << istep << endl;
        if(istep==1) str=str_orig;
        if(istep==2) str=str_relax1;
        if(istep==3) str=str_relax;
        if(istep==1) {
          if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [Volume Distance Bonds Analyzer: ORIGINAL]" << endl;
          ssfile << "Volume Distance Bonds Analyzer: ORIGINAL" << endl;
          stringstream sss;vector<string> vline;
          pflow::PrintData(xstructure(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO),sss,"DATA"); // DATA
          aurostd::string2vectorstring(sss.str(),vline);
          for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"real space volume")) ssfile << vline.at(iline) << endl;
          for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"Real space a b c alpha beta gamma")) ssfile << vline.at(iline) << endl;
        }
        if(istep==2) {
          if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [Volume Distance Bonds Analyzer: RELAX1]" << endl;
          ssfile << "Volume Distance Bonds Analyzer: RELAX1" << endl;
          stringstream sss;vector<string> vline;
          pflow::PrintData(xstructure(directory_RAW+"/CONTCAR.relax1",IOAFLOW_AUTO),sss,"DATA"); // DATA
          aurostd::string2vectorstring(sss.str(),vline);
          for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"real space volume")) ssfile << vline.at(iline) << endl;
          for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"Real space a b c alpha beta gamma")) ssfile << vline.at(iline) << endl;
        }
        if(istep==3) {
          if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [Volume Distance Bonds Analyzer: RELAX]" << endl;
          ssfile << "Volume Distance Bonds Analyzer: RELAX" << endl;
          stringstream sss;vector<string> vline;
          pflow::PrintData(xstructure(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO),sss,"DATA"); // DATA
          aurostd::string2vectorstring(sss.str(),vline);
          for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"real space volume")) ssfile << vline.at(iline) << endl;
          for(uint iline=0;iline<vline.size();iline++) if(aurostd::substring2bool(vline.at(iline),"Real space a b c alpha beta gamma")) ssfile << vline.at(iline) << endl;
        }
        if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [stop istep] " << istep << endl;
        if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << str << endl;
        if(VOLDISTEvolution_flag) if(AFLOWLIB_VERBOSE) cout << ssfile.str() << endl;

        vector<string> aus_data_nbondxx;
        string data1=pflow::PrintData1(str,-1.0);
        vector<string> vdata1;aurostd::string2vectorstring(data1,vdata1);
        uint _nspecies=data.nspecies;
        if(VOLDISTEvolution_flag) cerr << XPID << _nspecies << " " << data1.size() << " " << vdata1.size() << endl << data1 << endl;
        for(uint isp1=0;isp1<_nspecies;isp1++)
          for(uint isp2=isp1;isp2<_nspecies;isp2++)
            for(uint i=0;i<vdata1.size();i++) {
              if(VOLDISTEvolution_flag) cerr << XPID << isp1 << " " << isp2 << " " << i << endl;
              string string2find,string2search="Pairs between types: "+aurostd::utype2string(isp1)+" "+aurostd::utype2string(isp2);
              tokens.clear();tokens.push_back("");
              if(aurostd::substring2bool(vdata1.at(i),string2search)) {
                string2find=aurostd::substring2string(vdata1.at(i),string2search,FALSE);
                aurostd::string2tokens(string2find,tokens);
                aus_data_nbondxx.push_back(tokens.at(0));
              }
            }
        if(aurostd::abs((double) aus_data_nbondxx.size()-_nspecies*(_nspecies+1.0)/2.0)>0.1)
        { cerr << XPID << MESSAGE << " incompatible aus_data_nbondxx.size()[" << aus_data_nbondxx.size() << "]!=_nspecies*(_nspecies+1)/2)[" << (_nspecies*(_nspecies+1.0)/2.0) << "]" << endl;exit(0); }
        if(flag_LIB2==TRUE && _nspecies==1) { _nspecies++;aus_data_nbondxx.push_back("");aus_data_nbondxx.push_back(""); } // FIX for purity control and alien generation
        for(uint isp1=0,inbondxx=0;isp1<_nspecies;isp1++)
          for(uint isp2=isp1;isp2<_nspecies;isp2++)
            ssfile << "BOND_" << char('A'+isp1) << char('A'+isp2) << " " << aus_data_nbondxx.at(inbondxx++) << " [norm V_ATOM^0.33]" << endl;
      }
      //   aurostd::stringstream2file(ssfile,directory_RAW+"/"+DEFAULT_FILE_VOLDISTEVOLUTION_OUT);
    } // END VOLDISTEvolution
    if(LDEBUG) cerr << soliloquy << " [13]" << endl;
    // check for ERRORS
    if(flag_ENERGY1) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [flag_ENERGY1]" << endl;
      if(aurostd::abs(data.energy_atom-data1_energy_atom)>ENERGY_ATOM_ERROR_meV/1000.0)
        if(data.energy_atom<0.0 && data1_energy_atom>0.0) flag_ERROR=TRUE;
    }
    if(LDEBUG) cerr << soliloquy << " [14]" << endl;
    // done now WRITE aflowlib.out
    if(flag_ERROR) data.error_status="ERROR";

    if(LDEBUG) cerr << soliloquy << " [15]" << endl;

    // PERFORM EDATA DATA AND CIF STEPS -------------------------------------------------------------------------
    xstructure str,str_sp,str_sc;
    vector<xstructure> vcif;

    if(LDEBUG) cerr << soliloquy << " [16]" << endl;
    // PERFORM EDATA STEP
    if(flag_EDATA_ORIG_ || flag_EDATA_RELAX_) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA start: " << directory_RAW << endl;
    if(flag_EDATA_ORIG_) { // ORIG
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.orig")) {
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing orig (POSCAR.orig) text format: " << directory_RAW << endl;
        // [OBSOLETE] aurostd::execute("cd \""+directory_RAW+"\" && cat POSCAR.orig | aflow --edata > "+DEFAULT_FILE_EDATA_ORIG_OUT);
        str=xstructure(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
        //DX START
        //DX20180526 [OBSOLETE] str.directory = str_sp.directory = str_sc.directory = aflags.Directory;
        //DX END
        stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym=str; //CO20171027 //DX20180226 - set equal str
        aurostd::xoption vpflow_edata_orig; //DX20180823 - added xoption
        pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"EDATA",vpflow_edata_orig,"txt",false); // 1=EDATA //CO20171025 //CO20171027
        aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_OUT);
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing orig (POSCAR.orig) json format: " << directory_RAW << endl;
        stringstream jjj; //CO20171025
        pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"EDATA",vpflow_edata_orig,"json",true); // 1=EDATA, already_calculated!  //CO20171025 //CO20171027
        aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_EDATA_ORIG_JSON); //CO20171025
        vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
        //CO+DX START 20170713 - adding symmetry output to RAW
        string new_sym_file;
        //txt variants
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,new_sym_file);
        } //DX20171207 - added pgroupk_xtal
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_PATTERSON_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT,new_sym_file);
        } //DX20200520 - added pgroupk_Patterson
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,new_sym_file);
        } //CO20171025
        //json variants
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,new_sym_file);
        } //DX20171207 - added pgroupk_xtal
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_PATTERSON_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON,new_sym_file);
        } //DX20200520 - added pgroupk_Patterson
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,"orig",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,new_sym_file);
        } //CO20171025
        //cout << message;
        //CO+DX STOP 20170713 - adding symmetry output to RAW
        // now extract info
        //DX20180823 - extract info from xoption - START
        if(vpflow_edata_orig.flag("EDATA::CALCULATED")) {
          if(data.Bravais_lattice_orig.empty()) { // Bravais_Lattice_orig
            data.Bravais_lattice_orig=vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_TYPE"); 
          }
          if(data.lattice_variation_orig.empty()) { // Bravais_Lattice_Variation_orig
            data.lattice_variation_orig=vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_VARIATION_TYPE"); 
          }
          if(data.lattice_system_orig.empty()) { // Lattice_System_orig
            data.lattice_system_orig=vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_SYSTEM"); 
          }
          if(data.Pearson_symbol_orig.empty()) { // Pearson_orig
            data.Pearson_symbol_orig=vpflow_edata_orig.getattachedscheme("EDATA::PEARSON_SYMBOL"); 
          }
          //DX20190124 - extract additional info from xoption - START
          if(data.crystal_system_orig.empty()) { // crystal_system_orig
            data.crystal_system_orig=vpflow_edata_orig.getattachedscheme("EDATA::CRYSTAL_SYSTEM"); 
          }
          if(data.crystal_family_orig.empty()) { // crystal_family_orig
            data.crystal_family_orig=vpflow_edata_orig.getattachedscheme("EDATA::CRYSTAL_FAMILY"); 
          }
          if(data.point_group_Hermann_Mauguin_orig.empty()) { // point_group_Hermann_Mauguin_orig
            data.point_group_Hermann_Mauguin_orig=vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_HERMANN_MAUGUIN"); 
          }
          if(data.crystal_class_orig.empty()) { // crystal_class_orig
            data.crystal_class_orig=vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_CRYSTAL_CLASS"); 
          }
          if(data.point_group_Schoenflies_orig.empty()) { // point_group_Schoenflies_orig
            data.point_group_Schoenflies_orig=vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_SCHOENFLIES"); 
          }
          if(data.point_group_orbifold_orig.empty()) { // point_group_orbifold_orig
            data.point_group_orbifold_orig=vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_ORBIFOLD"); 
          }
          if(data.point_group_type_orig.empty()) { // point_group_type_orig
            data.point_group_type_orig=vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_TYPE"); 
          }
          if(data.point_group_order_orig==AUROSTD_NAN) { // point_group_order_orig
            data.point_group_order_orig=vpflow_edata_orig.getattachedutype<uint>("EDATA::POINT_GROUP_ORDER"); 
          }
          if(data.point_group_structure_orig.empty()) { // point_group_structure_orig
            data.point_group_structure_orig=vpflow_edata_orig.getattachedscheme("EDATA::POINT_GROUP_STRUCTURE"); 
          }
          if(data.Bravais_lattice_lattice_type_orig.empty()) { // Bravais_lattice_lattice_type_orig
            data.Bravais_lattice_lattice_type_orig=vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE"); 
          }
          if(data.Bravais_lattice_lattice_variation_type_orig.empty()) { // Bravais_lattice_lattice_variation_type_orig
            data.Bravais_lattice_lattice_variation_type_orig=vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE"); 
          }
          if(data.Bravais_lattice_lattice_system_orig.empty()) { // Bravais_lattice_lattice_system_orig
            data.Bravais_lattice_lattice_system_orig=vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM"); 
          }
          if(data.Bravais_superlattice_lattice_type_orig.empty()) { // Bravais_superlattice_lattice_type_orig 
            data.Bravais_superlattice_lattice_type_orig=vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_TYPE"); 
          }
          if(data.Bravais_superlattice_lattice_variation_type_orig.empty()) { // Bravais_superlattice_lattice_variation_type_orig
            data.Bravais_superlattice_lattice_variation_type_orig=vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE"); 
          }
          if(data.Bravais_superlattice_lattice_system_orig.empty()) { // Bravais_superlattice_lattice_system_orig
            data.Bravais_superlattice_lattice_system_orig=vpflow_edata_orig.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM"); 
          }
          if(data.Pearson_symbol_superlattice_orig.empty()) { // Pearson_symbol_superlattice_orig
            data.Pearson_symbol_superlattice_orig=vpflow_edata_orig.getattachedscheme("EDATA::PEARSON_SYMBOL_SUPERLATTICE"); 
          }
          if(data.reciprocal_geometry_orig.empty()) { // reciprocal_geometry_orig
            data.reciprocal_geometry_orig=vpflow_edata_orig.getattachedscheme("EDATA::RECIPROCAL_LATTICE_PARAMETERS"); 
            vector<string> ktokens; 
            aurostd::string2tokens(data.reciprocal_geometry_orig,ktokens,",");
            for(uint t=0;t<ktokens.size();t++) { data.vreciprocal_geometry_orig.push_back(aurostd::string2utype<double>(ktokens[t])); }
          }
          if(data.reciprocal_volume_cell_orig==AUROSTD_NAN) { // reciprocal_volume_cell_orig
            data.reciprocal_volume_cell_orig=vpflow_edata_orig.getattachedutype<double>("EDATA::RECIPROCAL_SPACE_VOLUME"); 
          }
          if(data.reciprocal_lattice_type_orig.empty()) { // reciprocal_lattice_type_orig
            data.reciprocal_lattice_type_orig=vpflow_edata_orig.getattachedscheme("EDATA::RECIPROCAL_LATTICE_TYPE"); 
          }
          if(data.reciprocal_lattice_variation_type_orig.empty()) { // reciprocal_lattice_variation_type_orig
            data.reciprocal_lattice_variation_type_orig=vpflow_edata_orig.getattachedscheme("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE"); 
          }
          //DX20190131 - use self-consistent space group orig - START
          if(data.spacegroup_orig.empty()) {
            data.spacegroup_orig=vpflow_edata_orig.getattachedscheme("SGDATA::SPACE_GROUP_NUMBER");
            if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP_ORIG = " << data.spacegroup_orig << endl;
          } 
          //DX20190131 - use self-consistent space group orig - END
          if(data.Wyckoff_letters_orig.empty()) { // Wyckoff_letters_orig
            data.Wyckoff_letters_orig=vpflow_edata_orig.getattachedscheme("SGDATA::WYCKOFF_LETTERS"); 
          }
          if(data.Wyckoff_multiplicities_orig.empty()) { // Wyckoff_multiplicities_orig
            data.Wyckoff_multiplicities_orig=vpflow_edata_orig.getattachedscheme("SGDATA::WYCKOFF_MULTIPLICITIES"); 
          }
          if(data.Wyckoff_site_symmetries_orig.empty()) { // Wyckoff_site_symmetries_orig
            data.Wyckoff_site_symmetries_orig=vpflow_edata_orig.getattachedscheme("SGDATA::WYCKOFF_SITE_SYMMETRIES"); 
          }
        }
        //DX20180823 - extract info from xoption - END
        else { //DX20180823 - if no xoption, read from file (safety)
          vector<string> vline_edata;
          aurostd::string2vectorstring(sss.str(),vline_edata);
          for(uint iline=0;iline<vline_edata.size();iline++) {
            if(data.Bravais_lattice_orig.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Bravais Lattice Primitive")) { // Bravais_Lattice
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.Bravais_lattice_orig=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
            if(data.lattice_variation_orig.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice Variation")) { // Bravais_Lattice_Variation
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.lattice_variation_orig=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
            if(data.lattice_system_orig.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice System")) { // Lattice_System
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.lattice_system_orig=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
            if(data.Pearson_symbol_orig.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Pearson Symbol") && !aurostd::substring2bool(vline_edata.at(iline),"Superlattice")) { // Pearson
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.Pearson_symbol_orig=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
          }
        }
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Bravais Lattice Primitive = " << data.Bravais_lattice_orig << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice Variation = " << data.lattice_variation_orig << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice System = " << data.lattice_system_orig << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Pearson Symbol = " << data.Pearson_symbol_orig << endl;
        //DX20190208 - add ANRL label/parameters/parameter values - START
        double anrl_symmetry_tolerance = str_sym.sym_eps;
        xstructure str_anrl = str;
        uint setting=SG_SETTING_ANRL;
        anrl::structure2anrl(str_anrl, anrl_symmetry_tolerance, setting); 
        if(data.anrl_label_orig.empty()) { // anrl label 
          data.anrl_label_orig = str_anrl.prototype; 
        }
        if(data.anrl_parameter_list_orig.empty()) { // anrl parameter list
          data.anrl_parameter_list_orig = aurostd::joinWDelimiter(str_anrl.prototype_parameter_list,","); 
        }
        if(data.anrl_parameter_values_orig.empty()) { // anrl parameter values
          data.anrl_parameter_values_orig = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(str_anrl.prototype_parameter_values,6),","); 
        }
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] ANRL Label = " << data.anrl_label_orig << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] ANRL parameter list = " << data.anrl_parameter_list_orig << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.ORIG] ANRL parameter values = " << data.anrl_parameter_values_orig << endl;
        //DX20190208 - add ANRL label/parameters/parameter values - END
      }
    }

    // ONLY TO GET REFERENCE ENERGIES AND INTERCEPT PROTOTYPES
    if(1) { // ONLY TO GENERATE AUIDS - SC ONLY DONT TOUCH
      if(data.nspecies==1) {
        string s=str_relax.species.at(0);
        string TXT1,TXT2;
        TXT1= "XXX if(AUID==\""+data.vspecies_pp_AUID.at(0);
        if(!aurostd::substring2bool(str_relax.species_pp_version.at(0),"SCAN")) { TXT1+="\" && nKIN) {found=TRUE;groundstate_structure=\""; } else { TXT1+="\" && SCAN) {found=TRUE;groundstate_structure=\""; }
        TXT2="\";groundstate_energy="+aurostd::utype2string<double>(data.enthalpy_atom,7)+";volume_atom="+aurostd::utype2string<double>(data.volume_atom,7);
        if(aurostd::abs(data.spin_atom)<0.1) { TXT2+=";spin_atom=0.0;} // "; } else { TXT2+=";spin_atom="+aurostd::utype2string<double>(data.spin_atom,7)+";} // ";}
        TXT2+=str_relax.species_pp_version.at(0);
        cout << data.anrl_label_orig << endl;
        // A1
        if(data.anrl_label_orig=="A_cF4_225_a" && (s=="Ac" || s=="Ag" || s=="Al" || s=="Au" || s=="Ca" || s=="Cu" || s=="Ir" || s=="La" || s=="Ni" || s=="Pb" || s=="Pd" || s=="Pt" || s=="Rh" || s=="Sr" || s=="Yb" || s=="Ar" || s=="Ne" || s=="Xe" || s=="Kr")) { // A1
          cout << TXT1 << "A1" << TXT2 << endl;}
        // A2
        if(data.anrl_label_orig=="A_cI2_229_a" && (s=="Ba" || s=="Cr" || s=="Fe" || s=="K" || s=="Li" || s=="Mo" || s=="Na" || s=="Nb" || s=="Ta" || s=="V" || s=="W" || s=="Cs" || s=="Eu")) { // A2
          cout << TXT1 << "A2" << TXT2 << endl;}
        // A3
        if(data.anrl_label_orig=="A_hP2_194_c" && (s=="Be" || s=="Cd" || s=="Co" || s=="Dy" || s=="Hf" || s=="Hg" || s=="Ho" || s=="Mg" || s=="Os" || s=="Re" || s=="Ru" || s=="Sc" || s=="Tc" || s=="Ti" || s=="Tl" || s=="Y" || s=="Zn" || s=="Zr" || s=="He")) { // A3
          cout << TXT1 << "A3" << TXT2 << endl;}
        // A4
        if(data.anrl_label_orig=="A_cF8_227_a" && (s=="Ge" || s=="Si")) { // A4 
          cout << TXT1 << "A4" << TXT2 << endl;}
        // A5
        if(data.anrl_label_orig=="A_tI4_141_a" && (s=="Sn")) { // A5 
          cout << TXT1 << "A5" << TXT2 << endl;}
        // A6
        if(data.anrl_label_orig=="A_tI2_139_a" && (s=="In")) { // A6
          cout << TXT1 << "A6" << TXT2 << endl;}
        // A7
        if(data.anrl_label_orig=="A_hR2_166_c" && (s=="As" || s=="Bi" || s=="Sb" || s=="P")) { // A7 P??
          cout << TXT1 << "A7" << TXT2 << endl;}
        // A8
        if(data.anrl_label_orig=="A_hP3_152_a" && (s=="Se" || s=="Te")) { // A8
          cout << TXT1 << "A8" << TXT2 << endl;}
        // A9
        if(data.anrl_label_orig=="A_hP4_194_bc" && (s=="C")) { // A9
          cout << TXT1 << "A9" << TXT2 << endl;}
        // A10
        if(data.anrl_label_orig=="A_hR1_166_a" && (s=="Hg")) { // A10
          cout << TXT1 << "A10" << TXT2 << endl;}
        // A11
        if(data.anrl_label_orig=="A_oC8_64_f" && (s=="Ga" || s=="Br")) { // A11 
          cout << TXT1 << "A11" << TXT2 << endl;}
        // A12
        if(data.anrl_label_orig=="A_cI58_217_ac2g" && (s=="Mn")) { // A12
          cout << TXT1 << "A12" << TXT2 << endl;} 
        // diatom (A_tP2_123_g) just RICO choice
        if(data.anrl_label_orig=="A_tP2_123_g" && (s=="O" || s=="N" || s=="F" || s=="H" || s=="Cl")) { // diatom //  cat /tmp/xscrubber_ppAUID.LIB1 | grep diatom | grep '/O'
          cout << TXT1 << "diatom" << TXT2 << endl;} 
        // A14
        if(data.anrl_label_orig=="A_oC8_64_f" && (s=="I")) { // A14
          cout << TXT1 << "A14" << TXT2 << endl;}
        // A16
        if(data.anrl_label_orig=="A_oF128_70_4h" && (s=="S")) { // A16
          cout << TXT1 << "A16" << TXT2 << endl;}
        if(data.anrl_label_orig=="A_hR12_166_2h" && (s=="B")) { // ICSD_56992
          cout << TXT1 << "ICSD_56992" << TXT2 << endl;}
        if(data.anrl_label_orig=="A_hR3_166_ac" && (s=="Sm")) { // C19
          cout << TXT1 << "C19" << TXT2 << endl;}
      }
    }


    if(flag_EDATA_RELAX_) { // RELAX
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT) && aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing relax (CONTCAR.relax) text format: " << directory_RAW << endl;
        // [OBSOLETE] aurostd::execute("cd \""+directory_RAW+"\" && cat CONTCAR.relax | aflow --edata > "+DEFAULT_FILE_EDATA_RELAX_OUT);
        str=xstructure(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
        stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym=str; //CO20171027 //DX20180226 - set equal str
        aurostd::xoption vpflow_edata_relax; //DX20180823 - added xoption
        pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"EDATA",vpflow_edata_relax,"txt",false); // EDATA //CO20171025 //CO20171027
        aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_OUT);
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing relax (CONTCAR.relax) json format: " << directory_RAW << endl;
        stringstream jjj; //CO20171025
        pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"EDATA",vpflow_edata_relax,"json",true); // EDATA already_calculated! //CO20171025 //CO20171027
        aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_EDATA_RELAX_JSON); //CO20171025
        vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
        //CO+DX START 20170713 - adding symmetry output to RAW
        string new_sym_file;
        //txt variants
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,new_sym_file);
        } //DX20171207 - added pgroupk_xtal
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_PATTERSON_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT,new_sym_file);
        } //DX20200520 - added pgroupk_Patterson
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,new_sym_file);
        } //CO20171025
        //json variants
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,new_sym_file);
        } //DX20171207 - added pgroupk_xtal
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_PATTERSON_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON,new_sym_file);
        } //DX20200520 - added pgroupk_Patterson
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,"relax",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,new_sym_file);
        } //CO20171025
        //cout << message;
        //CO+DX STOP 20170713 - adding symmetry output to RAW
        // now extract info
        //DX20180823 - extract info from xoption - START
        if(vpflow_edata_relax.flag("EDATA::CALCULATED")) {
          if(data.Bravais_lattice_relax.empty()) { // Bravais_Lattice
            data.Bravais_lattice_relax=vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_TYPE"); 
          }
          if(data.lattice_variation_relax.empty()) { // Bravais_Lattice_Variation
            data.lattice_variation_relax=vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_VARIATION_TYPE"); 
          }
          if(data.lattice_system_relax.empty()) { // Lattice_System
            data.lattice_system_relax=vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_SYSTEM"); 
          }
          if(data.Pearson_symbol_relax.empty()) { // Pearson
            data.Pearson_symbol_relax=vpflow_edata_relax.getattachedscheme("EDATA::PEARSON_SYMBOL"); 
          }
          if(data.crystal_system.empty()) { // crystal_system
            data.crystal_system=vpflow_edata_relax.getattachedscheme("EDATA::CRYSTAL_SYSTEM"); 
          }
          if(data.crystal_family.empty()) { // crystal_family
            data.crystal_family=vpflow_edata_relax.getattachedscheme("EDATA::CRYSTAL_FAMILY"); 
          }
          //DX20190115 - typo in Hermann-Mauguin and crystal class entries - START
          if(data.point_group_Hermann_Mauguin.empty()) { // point_group_Hermann_Mauguin
            data.point_group_Hermann_Mauguin=vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_HERMANN_MAUGUIN"); 
          }
          if(data.crystal_class.empty()) { // crystal_class
            data.crystal_class=vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_CRYSTAL_CLASS"); 
          }
          //DX20190115 - typo in Hermann-Mauguin and crystal class entries - END
          if(data.point_group_Schoenflies.empty()) { // point_group_Schoenflies
            data.point_group_Schoenflies=vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_SCHOENFLIES"); 
          }
          if(data.point_group_orbifold.empty()) { // point_group_orbifold
            data.point_group_orbifold=vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_ORBIFOLD"); 
          }
          if(data.point_group_type.empty()) { // point_group_type 
            data.point_group_type=vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_TYPE"); 
          }
          if(data.point_group_order==AUROSTD_NAN) { // point_group_order 
            data.point_group_order=vpflow_edata_relax.getattachedutype<uint>("EDATA::POINT_GROUP_ORDER"); 
          }
          if(data.point_group_structure.empty()) { // point_group_structure 
            data.point_group_structure=vpflow_edata_relax.getattachedscheme("EDATA::POINT_GROUP_STRUCTURE"); 
          }
          if(data.Bravais_lattice_lattice_type.empty()) { // Bravais_lattice_lattice_type 
            data.Bravais_lattice_lattice_type=vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE"); 
          }
          //DX20190115 - typo, missing "variation" - START
          if(data.Bravais_lattice_lattice_variation_type.empty()) { // Bravais_lattice_lattice_variation_type 
            data.Bravais_lattice_lattice_variation_type=vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE"); 
          }
          //DX20190115 - typo, missing "variation" - END
          if(data.Bravais_lattice_lattice_system.empty()) { // Bravais_lattice_lattice_system 
            data.Bravais_lattice_lattice_system=vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM"); 
          }
          if(data.Bravais_superlattice_lattice_type.empty()) { // Bravais_superlattice_lattice_type 
            data.Bravais_superlattice_lattice_type=vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_TYPE"); 
          }
          if(data.Bravais_superlattice_lattice_variation_type.empty()) { // Bravais_superlattice_lattice_variation_type 
            data.Bravais_superlattice_lattice_variation_type=vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE"); 
          }
          if(data.Bravais_superlattice_lattice_system.empty()) { // Bravais_superlattice_lattice_system 
            data.Bravais_superlattice_lattice_system=vpflow_edata_relax.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM"); 
          }
          if(data.Pearson_symbol_superlattice.empty()) { // Pearson_symbol_superlattice 
            data.Pearson_symbol_superlattice=vpflow_edata_relax.getattachedscheme("EDATA::PEARSON_SYMBOL_SUPERLATTICE"); 
          }
          if(data.reciprocal_geometry.empty()) { // reciprocal_geometry
            data.reciprocal_geometry=vpflow_edata_relax.getattachedscheme("EDATA::RECIPROCAL_LATTICE_PARAMETERS"); 
            vector<string> ktokens; 
            aurostd::string2tokens(data.reciprocal_geometry,ktokens,",");
            for(uint t=0;t<ktokens.size();t++) { data.vreciprocal_geometry.push_back(aurostd::string2utype<double>(ktokens[t])); }
          }
          if(data.reciprocal_volume_cell==AUROSTD_NAN) { // reciprocal_volume_cell
            data.reciprocal_volume_cell=vpflow_edata_relax.getattachedutype<double>("EDATA::RECIPROCAL_SPACE_VOLUME"); 
          }
          if(data.reciprocal_lattice_type.empty()) { // reciprocal_lattice_type 
            data.reciprocal_lattice_type=vpflow_edata_relax.getattachedscheme("EDATA::RECIPROCAL_LATTICE_TYPE"); 
          }
          if(data.reciprocal_lattice_variation_type.empty()) { // reciprocal_lattice_variation_type 
            data.reciprocal_lattice_variation_type=vpflow_edata_relax.getattachedscheme("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE"); 
          }
          //DX20190131 - use self-consistent space group relax - START
          if(data.spacegroup_relax.empty()) {
            data.spacegroup_relax=vpflow_edata_relax.getattachedscheme("SGDATA::SPACE_GROUP_NUMBER");
            if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP_RELAX = " << data.spacegroup_relax << endl;
          } 
          //DX20190131 - use self-consistent space group orig - END
          if(data.Wyckoff_letters.empty()) { // Wyckoff_letters 
            data.Wyckoff_letters=vpflow_edata_relax.getattachedscheme("SGDATA::WYCKOFF_LETTERS"); 
          }
          if(data.Wyckoff_multiplicities.empty()) { // Wyckoff_multiplicities 
            data.Wyckoff_multiplicities=vpflow_edata_relax.getattachedscheme("SGDATA::WYCKOFF_MULTIPLICITIES"); 
          }
          if(data.Wyckoff_site_symmetries.empty()) { // Wyckoff_site_symmetries
            data.Wyckoff_site_symmetries=vpflow_edata_relax.getattachedscheme("SGDATA::WYCKOFF_SITE_SYMMETRIES"); 
          }
        }
        //DX20180823 - extract info from xoption - END
        else { //DX20180823 - if no xoption, read from file (safety)
          vector<string> vline_edata;
          aurostd::string2vectorstring(sss.str(),vline_edata);
          for(uint iline=0;iline<vline_edata.size();iline++) {
            if(data.Bravais_lattice_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Bravais Lattice Primitive")) { // Bravais_Lattice
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.Bravais_lattice_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
            if(data.lattice_variation_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice Variation")) { // Bravais_Lattice_Variation
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.lattice_variation_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
            if(data.lattice_system_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice System")) { // Lattice_System
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.lattice_system_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
            if(data.Pearson_symbol_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Pearson Symbol") && !aurostd::substring2bool(vline_edata.at(iline),"Superlattice")) { // Pearson
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.Pearson_symbol_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
          }
        } 
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Bravais Lattice Primitive = " << data.Bravais_lattice_relax << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice Variation = " << data.lattice_variation_relax << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice System = " << data.lattice_system_relax << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Pearson Symbol = " << data.Pearson_symbol_relax << endl;
        //DX20190208 - add ANRL label/parameters/parameter values - START
        double anrl_symmetry_tolerance = str_sym.sym_eps;
        xstructure str_anrl = str;
        uint setting=SG_SETTING_ANRL;
        anrl::structure2anrl(str_anrl, anrl_symmetry_tolerance, setting); 
        if(data.anrl_label_relax.empty()) { // anrl label 
          data.anrl_label_relax = str_anrl.prototype; 
        }
        if(data.anrl_parameter_list_relax.empty()) { // anrl parameter list
          data.anrl_parameter_list_relax = aurostd::joinWDelimiter(str_anrl.prototype_parameter_list,","); 
        }
        if(data.anrl_parameter_values_relax.empty()) { // anrl parameter values
          data.anrl_parameter_values_relax = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(str_anrl.prototype_parameter_values,6),","); 
        }
        //DX20190208 - add ANRL label/parameters/parameter values - END
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] ANRL Label = " << data.anrl_label_relax << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] ANRL parameter list = " << data.anrl_parameter_list_relax << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.RELAX] ANRL parameter values = " << data.anrl_parameter_values_relax << endl;
      }
    }
    if(LDEBUG) cerr << soliloquy << " [17]" << endl;
    if(flag_EDATA_BANDS_) { // BANDS
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_EDATA_BANDS_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.bands")) {
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing bands (POSCAR.bands) text format: " << directory_RAW << endl;
        // [OBSOLETE] aurostd::execute("cd \""+directory_RAW+"\" && cat POSCAR.bands | aflow --edata > "+DEFAULT_FILE_EDATA_BANDS_OUT);
        str=xstructure(directory_RAW+"/POSCAR.bands",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
        stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym=str; //CO20171027 //DX20180226 - set equal str
        aurostd::xoption vpflow_edata_bands; //DX20180823 - added xoption
        pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"EDATA",vpflow_edata_bands,"txt",false); // EDATA //CO20171025 //CO20171027
        aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_EDATA_BANDS_OUT);
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " EDATA doing bands (POSCAR.bands) json format: " << directory_RAW << endl;
        stringstream jjj; //CO20171025
        pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"EDATA",vpflow_edata_bands,"json",true); // EDATA already_calculated! //CO20171025 //CO20171027
        aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_EDATA_BANDS_JSON); //CO20171025
        vcif.clear();vcif.push_back(str);vcif.push_back(str_sp);vcif.push_back(str_sc);
        //CO+DX START 20170713 - adding symmetry output to RAW
        string new_sym_file;
        //txt variants
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_OUT,new_sym_file);
        } //DX20171207 - added pgroupk_xtal
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_PATTERSON_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_OUT,new_sym_file);
        } //DX20200520 - added pgroupk_Patterson
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_OUT,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"txt");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_OUT,new_sym_file);
        } //CO20171025
        //json variants
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_FGROUP_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_FGROUP_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUP_XTAL_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUP_XTAL_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_XTAL_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_XTAL_JSON,new_sym_file);
        } //DX20171207 - added pgroupk_xtal
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_PGROUPK_PATTERSON_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_PGROUPK_PATTERSON_JSON,new_sym_file);
        } //DX20200520 - added pgroupk_Patterson
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_IATOMS_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_IATOMS_JSON,new_sym_file);
        } //CO20171025
        KBIN_SymmetryWrite(FileMESSAGE,str_sp,aflags,_AGROUP_,false,message,"json");
        if(aurostd::FileExist(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON)) {
          AddFileNameBeforeExtension(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,"bands",new_sym_file);
          aurostd::file2file(directory_RAW+"/"+DEFAULT_AFLOW_AGROUP_JSON,new_sym_file);
        } //CO20171025
        //cout << message;
        //CO+DX STOP 20170713 - adding symmetry output to RAW
        // now extract info
        //DX20180823 - extract info from xoption - START
        if(vpflow_edata_bands.flag("EDATA::CALCULATED")) {
          if(data.Bravais_lattice_relax.empty()) { // Bravais_Lattice
            data.Bravais_lattice_relax=vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_TYPE"); 
          }
          if(data.lattice_variation_relax.empty()) { // Bravais_Lattice_Variation
            data.lattice_variation_relax=vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_VARIATION_TYPE"); 
          }
          if(data.lattice_system_relax.empty()) { // Lattice_System
            data.lattice_system_relax=vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_SYSTEM"); 
          }
          if(data.Pearson_symbol_relax.empty()) { // Pearson
            data.Pearson_symbol_relax=vpflow_edata_bands.getattachedscheme("EDATA::PEARSON_SYMBOL"); 
          }
          if(data.crystal_system.empty()) { // crystal_system
            data.crystal_system=vpflow_edata_bands.getattachedscheme("EDATA::CRYSTAL_SYSTEM"); 
          }
          if(data.crystal_family.empty()) { // crystal_family
            data.crystal_family=vpflow_edata_bands.getattachedscheme("EDATA::CRYSTAL_FAMILY"); 
          }
          //DX20190115 - typo in Hermann-Mauguin and crystal class entries - START
          if(data.point_group_Hermann_Mauguin.empty()) { // point_group_Hermann_Mauguin
            data.point_group_Hermann_Mauguin=vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_HERMANN_MAUGUIN"); 
          }
          if(data.crystal_class.empty()) { // crystal_class
            data.crystal_class=vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_CRYSTAL_CLASS"); 
          }
          //DX20190115 - typo in Hermann-Mauguin and crystal class entries - END
          if(data.point_group_Schoenflies.empty()) { // point_group_Schoenflies
            data.point_group_Schoenflies=vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_SCHOENFLIES"); 
          }
          if(data.point_group_orbifold.empty()) { // point_group_orbifold
            data.point_group_orbifold=vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_ORBIFOLD"); 
          }
          if(data.point_group_type.empty()) { // point_group_type 
            data.point_group_type=vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_TYPE"); 
          }
          if(data.point_group_order==AUROSTD_NAN) { // point_group_order 
            data.point_group_order=vpflow_edata_bands.getattachedutype<uint>("EDATA::POINT_GROUP_ORDER"); 
          }
          if(data.point_group_structure.empty()) { // point_group_structure 
            data.point_group_structure=vpflow_edata_bands.getattachedscheme("EDATA::POINT_GROUP_STRUCTURE"); 
          }
          if(data.Bravais_lattice_lattice_type.empty()) { // Bravais_lattice_lattice_type 
            data.Bravais_lattice_lattice_type=vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_TYPE"); 
          }
          //DX20190115 - typo, missing "variation" - START
          if(data.Bravais_lattice_lattice_variation_type.empty()) { // Bravais_lattice_lattice_variation_type 
            data.Bravais_lattice_lattice_variation_type=vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_VARIATION_TYPE"); 
          }
          //DX20190115 - typo, missing "variation" - END
          if(data.Bravais_lattice_lattice_system.empty()) { // Bravais_lattice_lattice_system 
            data.Bravais_lattice_lattice_system=vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_LATTICE_LATTICE_SYSTEM"); 
          }
          if(data.Bravais_superlattice_lattice_type.empty()) { // Bravais_superlattice_lattice_type 
            data.Bravais_superlattice_lattice_type=vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_TYPE"); 
          }
          if(data.Bravais_superlattice_lattice_variation_type.empty()) { // Bravais_superlattice_lattice_variation_type 
            data.Bravais_superlattice_lattice_variation_type=vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_VARIATION_TYPE"); 
          }
          if(data.Bravais_superlattice_lattice_system.empty()) { // Bravais_superlattice_lattice_system 
            data.Bravais_superlattice_lattice_system=vpflow_edata_bands.getattachedscheme("EDATA::BRAVAIS_SUPERLATTICE_SYSTEM"); 
          }
          if(data.Pearson_symbol_superlattice.empty()) { // Pearson_symbol_superlattice 
            data.Pearson_symbol_superlattice=vpflow_edata_bands.getattachedscheme("EDATA::PEARSON_SYMBOL_SUPERLATTICE"); 
          }
          if(data.reciprocal_geometry.empty()) { // reciprocal_geometry
            data.reciprocal_geometry=vpflow_edata_bands.getattachedscheme("EDATA::RECIPROCAL_LATTICE_PARAMETERS"); 
            vector<string> ktokens;
            aurostd::string2tokens(data.reciprocal_geometry,ktokens,",");
            for(uint t=0;t<ktokens.size();t++) { data.vreciprocal_geometry.push_back(aurostd::string2utype<double>(ktokens[t])); }
          }
          if(data.reciprocal_volume_cell==AUROSTD_NAN) { // reciprocal_volume_cell
            data.reciprocal_volume_cell=vpflow_edata_bands.getattachedutype<double>("EDATA::RECIPROCAL_SPACE_VOLUME"); 
          }
          if(data.reciprocal_lattice_type.empty()) { // reciprocal_lattice_type 
            data.reciprocal_lattice_type=vpflow_edata_bands.getattachedscheme("EDATA::RECIPROCAL_LATTICE_TYPE"); 
          }
          if(data.reciprocal_lattice_variation_type.empty()) { // reciprocal_lattice_variation_type 
            data.reciprocal_lattice_variation_type=vpflow_edata_bands.getattachedscheme("EDATA::RECIPROCAL_LATTICE_VARIATION_TYPE"); 
          }
          //DX20190131 - use self-consistent space group relax - START
          if(data.spacegroup_relax.empty()) {
            data.spacegroup_relax=vpflow_edata_bands.getattachedscheme("SGDATA::SPACE_GROUP_NUMBER");
            if(AFLOWLIB_VERBOSE) cout << MESSAGE << " SPACEGROUP_RELAX = " << data.spacegroup_relax << endl;
          } 
          //DX20190131 - use self-consistent space group orig - END
          if(data.Wyckoff_letters.empty()) { // Wyckoff_letters 
            data.Wyckoff_letters=vpflow_edata_bands.getattachedscheme("SGDATA::WYCKOFF_LETTERS"); 
          }
          if(data.Wyckoff_multiplicities.empty()) { // Wyckoff_multiplicities 
            data.Wyckoff_multiplicities=vpflow_edata_bands.getattachedscheme("SGDATA::WYCKOFF_MULTIPLICITIES"); 
          }
          if(data.Wyckoff_site_symmetries.empty()) { // Wyckoff_site_symmetries
            data.Wyckoff_site_symmetries=vpflow_edata_bands.getattachedscheme("SGDATA::WYCKOFF_SITE_SYMMETRIES"); 
          }
        }
        //DX20180823 - extract info from xoption - END
        else { //DX20180823 - if no xoption, read from file (safety)
          vector<string> vline_edata;
          aurostd::string2vectorstring(sss.str(),vline_edata);
          for(uint iline=0;iline<vline_edata.size();iline++) {
            if(data.Bravais_lattice_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Bravais Lattice Primitive")) { // Bravais_Lattice
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.Bravais_lattice_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
            if(data.lattice_variation_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice Variation")) { // Bravais_Lattice_Variation
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.lattice_variation_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
            if(data.lattice_system_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Lattice System")) { // Lattice_System
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.lattice_system_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
            if(data.Pearson_symbol_relax.empty() && aurostd::substring2bool(vline_edata.at(iline),"Real space: Pearson Symbol") && !aurostd::substring2bool(vline_edata.at(iline),"Superlattice")) { // Pearson
              aurostd::string2tokens(vline_edata.at(iline),tokens,"=");
              data.Pearson_symbol_relax=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)); }
          }
        }
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Bravais Lattice Primitive = " << data.Bravais_lattice_relax << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice Variation = " << data.lattice_variation_relax << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Lattice System = " << data.lattice_system_relax << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] BRAVAIS LATTICE OF THE CRYSTAL (pgroup_xtal) - Real space: Pearson Symbol = " << data.Pearson_symbol_relax << endl;
        //DX20190208 - add ANRL label/parameters/parameter values - START
        double anrl_symmetry_tolerance = str_sym.sym_eps;
        xstructure str_anrl = str;
        uint setting=SG_SETTING_ANRL;
        anrl::structure2anrl(str_anrl, anrl_symmetry_tolerance, setting); 
        if(data.anrl_label_relax.empty()) { // anrl label 
          data.anrl_label_relax = str_anrl.prototype; 
        }
        if(data.anrl_parameter_list_relax.empty()) { // anrl parameter list
          data.anrl_parameter_list_relax = aurostd::joinWDelimiter(str_anrl.prototype_parameter_list,","); 
        }
        if(data.anrl_parameter_values_relax.empty()) { // anrl parameter values
          data.anrl_parameter_values_relax = aurostd::joinWDelimiter(aurostd::vecDouble2vecString(str_anrl.prototype_parameter_values,6),","); 
        }
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] ANRL Label = " << data.anrl_label_relax << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] ANRL parameter list = " << data.anrl_parameter_list_relax << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " [EDATA.BANDS] ANRL parameter values = " << data.anrl_parameter_values_relax << endl;
        //DX20190208 - add ANRL label/parameters/parameter values - END
      }
    }
    if(LDEBUG) cerr << soliloquy << " [18]" << endl;
    // PERFORM DATA STEP
    if(flag_DATA_ORIG_ || flag_DATA_RELAX_) if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA start: " << directory_RAW << endl;
    if(flag_DATA_ORIG_) { // ORIG
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_ORIG_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.orig")) {
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing orig (POSCAR.orig) text format: " << directory_RAW << endl;
        str=xstructure(directory_RAW+"/POSCAR.orig",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
        stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym=str; //CO20171027 //DX20180226 - set equal str
        pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"DATA","txt",false); // DATA //CO20171025 //CO20171027
        aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_ORIG_OUT);
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing orig (POSCAR.orig) json format: " << directory_RAW << endl;
        stringstream jjj; //CO20171025
        pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"DATA","json",true); // DATA already_calculated! //CO20171025 //CO20171027
        aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_DATA_ORIG_JSON); //CO20171025
      }
    }
    if(LDEBUG) cerr << soliloquy << " [19]" << endl;
    if(flag_DATA_RELAX_) { // RELAX
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_RELAX_OUT) && aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing relax (CONTCAR.relax) text format: " << directory_RAW << endl;
        str=xstructure(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
        stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym=str; //CO20171027 //DX20180226 - set equal str
        pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"DATA","txt",false); // DATA //CO20171025 //CO20171027
        aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_RELAX_OUT);
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing relax (CONTCAR.relax) json format: " << directory_RAW << endl;
        stringstream jjj; //CO20171025
        pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"DATA","json",true); // DATA already_calculated! //CO20171025 //CO20171027
        aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_DATA_RELAX_JSON); //CO20171025
      }
    }
    if(LDEBUG) cerr << soliloquy << " [20]" << endl;
    if(flag_DATA_BANDS_) { // BANDS
      if(!aurostd::FileExist(directory_RAW+"/"+DEFAULT_FILE_DATA_BANDS_OUT) && aurostd::FileExist(directory_RAW+"/POSCAR.bands")) {
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing relax (POSCAR.bands) text format: " << directory_RAW << endl;
        str=xstructure(directory_RAW+"/POSCAR.bands",IOAFLOW_AUTO);str_sp.clear();str_sc.clear(); //DX20191220 - uppercase to lowercase clear
        stringstream sss; sss << aflow::Banner("BANNER_TINY") << endl;
        xstructure str_sym=str; //CO20171027 //DX20180226 - set equal str
        pflow::PrintData(str,str_sym,str_sp,str_sc,sss,"DATA","txt",false); // DATA //CO20171025 //CO20171027
        aurostd::stringstream2file(sss,directory_RAW+"/"+DEFAULT_FILE_DATA_BANDS_OUT);
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " DATA doing relax (POSCAR.bands) json format: " << directory_RAW << endl;
        stringstream jjj; //CO20171025
        pflow::PrintData(str_sym,str_sym,str_sp,str_sc,jjj,"DATA","json",true); // DATA already_calculated! //CO20171025 //CO20171027
        aurostd::stringstream2file(jjj,directory_RAW+"/"+DEFAULT_FILE_DATA_BANDS_JSON); //CO20171025
      }
    }

    if(LDEBUG) cerr << soliloquy << " [21]" << endl;
    // CREATING CIFS
    //      cout << MESSAGE << " CIF creation: " << directory_LIB << endl;
    if(vcif.size()==0 && aurostd::FileExist(directory_RAW+"/POSCAR.bands")) vcif.push_back(xstructure(directory_RAW+"/POSCAR.bands",IOVASP_AUTO));
    xvector<double> nvec(3);nvec(1)=1;nvec(2)=1;nvec(3)=1;
    double angle=45;

    for (uint j=0;j<vcif.size();j++) {
      vcif.at(j)=GetLTFVCell(nvec,angle,vcif.at(j));
      if(j==0) cout << MESSAGE << " CIF creation: data.spacegroup_relax=" << data.spacegroup_relax << endl;
      if(j==0) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing normal" << endl;
      if(j==1) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sprim" << endl;
      if(j==2) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sconv" << endl;
      stringstream oss;
      //    for(uint i=0;i<data.vspecies.size();i++) cerr << XPID << data.vspecies.at(i) << endl;
      vcif.at(j).species.clear();for(uint i=0;i<data.vspecies.size();i++) vcif.at(j).species.push_back(data.vspecies.at(i));
      vcif.at(j).species_pp.clear();for(uint i=0;i<data.vspecies.size();i++) vcif.at(j).species_pp.push_back(data.vspecies.at(i));
      vcif.at(j).species_pp_type.clear();for(uint i=0;i<data.vspecies.size();i++) vcif.at(j).species_pp_type.push_back("");
      vcif.at(j).species_pp_version.clear();for(uint i=0;i<data.vspecies.size();i++) vcif.at(j).species_pp_version.push_back("");
      vcif.at(j).species_pp_ZVAL.clear();for(uint i=0;i<data.vspecies.size();i++) vcif.at(j).species_pp_ZVAL.push_back(0.0);
      vcif.at(j).species_pp_vLDAU.clear();for(uint i=0;i<data.vspecies.size();i++) vcif.at(j).species_pp_vLDAU.push_back(deque<double>());
      for(uint i=0;i<vcif.at(j).atoms.size();i++) vcif.at(j).atoms.at(i).name=vcif.at(j).species.at(vcif.at(j).atoms.at(i).type);
      for(uint i=0;i<vcif.at(j).atoms.size();i++) vcif.at(j).atoms.at(i).cleanname=vcif.at(j).species.at(vcif.at(j).atoms.at(i).type);
      pflow::PrintCIF(oss,vcif.at(j),1);//aurostd::string2utype<int>(data.spacegroup_relax));
      if(j==0) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+".cif");
      if(j==1) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sprim.cif");
      if(j==2) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sconv.cif");
      vcif.at(j).AddCorners();
      oss.clear();oss.str("");
      pflow::PrintCIF(oss,vcif.at(j),1);//aurostd::string2utype<int>(data.spacegroup_relax));
      if(j==0) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing normal_corner" << endl;
      if(j==1) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sprim_corner" << endl;
      if(j==2) cout << MESSAGE << " CIF creation: " << directory_LIB << " doing sconv_corner" << endl;
      if(j==0) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_corner.cif");
      if(j==1) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sprim_corner.cif");
      if(j==2) aurostd::stringstream2file(oss,directory_RAW+"/"+KBIN::ExtractSystemName(directory_LIB)+"_sconv_corner.cif");
    }

    // CREATING GEOMETRY
    if(aurostd::FileExist(directory_RAW+"/CONTCAR.relax")) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " LOAD SPECIES=" << data.species << endl;
      xstructure xstr(directory_RAW+"/CONTCAR.relax",IOAFLOW_AUTO);
      //[CO20200624 - OBSOLETE]deque<string> vspecies;
      //[CO20200624 - OBSOLETE]aurostd::string2tokens(data.species,vspecies,",");
      xstr.SetSpecies(deq_species);  //CO20200624
      stringstream sss_vasp;sss_vasp << xstr;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " WRITE relax positions for VASP" << endl;
      aurostd::stringstream2file(sss_vasp,directory_RAW+"/CONTCAR.relax.vasp");
      //	cerr << XPID << xstr << endl;
      xstr.xstructure2qe();
      stringstream sss_qe;sss_qe << xstr;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " WRITE relax positions for QE" << endl;
      aurostd::stringstream2file(sss_qe,directory_RAW+"/CONTCAR.relax.qe");
      //	cerr << XPID << xstr << endl;
      xstr.xstructure2abinit();
      stringstream sss_abinit;sss_abinit << xstr;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " WRITE relax positions for ABINIT" << endl;
      aurostd::stringstream2file(sss_abinit,directory_RAW+"/CONTCAR.relax.abinit");
      //	cerr << XPID << xstr << endl;
      xstr.xstructure2aims();
      stringstream sss_aims;sss_aims << xstr;
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " WRITE relax positions for AIMS" << endl;
      aurostd::stringstream2file(sss_aims,directory_RAW+"/CONTCAR.relax.aims");
      //	cerr << XPID << xstr << endl;
    }

    if(flag_LIB2==TRUE) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " PATCHING for LIB2: " << directory_RAW << endl;
      if(_APENNSY_STYLE_OLD_) {
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          fileA_LIB=directory_LIB+"/aflow.agroup.out"+XHOST.vext.at(iext);
          fileA_RAW=directory_RAW+"/aflow.agroup.out"+XHOST.vext.at(iext);
          if(aurostd::FileExist(fileA_LIB)) {
            aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.agroup.out"+XHOST.vext.at(iext),directory_RAW,"aflow.agroup.out"+XHOST.vext.at(iext),vfile,MESSAGE);
          }
        }
      }
      if(_APENNSY_STYLE_OLD_) {
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          fileA_LIB=directory_LIB+"/aflow.fgroup.out"+XHOST.vext.at(iext);
          fileA_RAW=directory_RAW+"/aflow.fgroup.out"+XHOST.vext.at(iext);
          if(aurostd::FileExist(fileA_LIB)) {
            aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.fgroup.out"+XHOST.vext.at(iext),directory_RAW,"aflow.fgroup.out"+XHOST.vext.at(iext),vfile,MESSAGE);
          }
        }
      }
      if(_APENNSY_STYLE_OLD_) {
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          fileA_LIB=directory_LIB+"/aflow.pgroup.out"+XHOST.vext.at(iext);
          fileA_RAW=directory_RAW+"/aflow.pgroup.out"+XHOST.vext.at(iext);
          if(aurostd::FileExist(fileA_LIB)) {
            aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.pgroup.out"+XHOST.vext.at(iext),directory_RAW,"aflow.pgroup.out"+XHOST.vext.at(iext),vfile,MESSAGE);
          }
        }
      }
      if(_APENNSY_STYLE_OLD_) {
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          fileA_LIB=directory_LIB+"/aflow.iatoms.out"+XHOST.vext.at(iext);
          fileA_RAW=directory_RAW+"/aflow.iatoms.out"+XHOST.vext.at(iext);
          if(aurostd::FileExist(fileA_LIB)) {
            aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.iatoms.out"+XHOST.vext.at(iext),directory_RAW,"aflow.iatoms.out"+XHOST.vext.at(iext),vfile,MESSAGE);
          }
        }
      }
      for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
        fileA_LIB=directory_LIB+"/aflow.qmvasp.out"+XHOST.vext.at(iext);
        fileA_RAW=directory_RAW+"/aflow.qmvasp.out"+XHOST.vext.at(iext);
        if(aurostd::FileExist(fileA_LIB)) {
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.qmvasp.out"+XHOST.vext.at(iext),directory_RAW,"aflow.qmvasp.out"+XHOST.vext.at(iext),vfile,MESSAGE);
        }
      }
      //     for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
      // no vasp.out.relax1  fileA_LIB=directory_LIB+"/vasp.out.relax1"+XHOST.vext.at(iext);fileA_RAW=directory_RAW+"/vasp.out.relax1"+XHOST.vext.at(iext);if(aurostd::FileExist(fileA_LIB)) { aflowlib::LIB2RAW_FileNeeded(directory_LIB,"vasp.out.relax1"+XHOST.vext.at(iext),directory_RAW,"vasp.out.relax1"+XHOST.vext.at(iext),vfile,MESSAGE); }
      // no vasp.out.relax2  fileA_LIB=directory_LIB+"/vasp.out.relax2"+XHOST.vext.at(iext);fileA_RAW=directory_RAW+"/vasp.out.relax2"+XHOST.vext.at(iext);if(aurostd::FileExist(fileA_LIB)) { aflowlib::LIB2RAW_FileNeeded(directory_LIB,"vasp.out.relax2"+XHOST.vext.at(iext),directory_RAW,"vasp.out.relax2"+XHOST.vext.at(iext),vfile,MESSAGE); }
      // }
      if(_APENNSY_STYLE_OLD_) {
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          fileX_LIB=directory_LIB+"/POSCAR.relax1"+XHOST.vext.at(iext);
          fileX_RAW=directory_RAW+"/POSCAR.relax1"+XHOST.vext.at(iext);
          if(aurostd::FileExist(fileX_LIB)) {
            aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.relax1"+XHOST.vext.at(iext),directory_RAW,"POSCAR.relax1"+XHOST.vext.at(iext),vfile,MESSAGE);
          }
        }
      }
      for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
        fileX_LIB=directory_LIB+"/POSCAR.relax2"+XHOST.vext.at(iext);
        fileX_RAW=directory_RAW+"/POSCAR.relax2"+XHOST.vext.at(iext);
        if(aurostd::FileExist(fileX_LIB)) {
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,"POSCAR.relax2"+XHOST.vext.at(iext),directory_RAW,"POSCAR.relax2"+XHOST.vext.at(iext),vfile,MESSAGE);
        }
      }
      aurostd::CopyFile(directory_RAW+"/CONTCAR.relax",directory_RAW+"/CONTCAR.relax2");vfile.push_back("CONTCAR.relax2");
      // aurostd::CopyFile(directory_LIB+"/OSZICAR.relax1"+XHOST.vext.at(iext),directory_RAW+"/OSZICAR.relax1"+XHOST.vext.at(iext));vfile.push_back("OSZICAR.relax1");
      // aurostd::CopyFile(directory_LIB+"/OSZICAR.relax2"+XHOST.vext.at(iext),directory_RAW+"/OSZICAR.relax2"+XHOST.vext.at(iext));vfile.push_back("OSZICAR.relax2");
      if(_APENNSY_STYLE_OLD_) {
        aurostd::CopyFile(directory_RAW+"/OSZICAR.relax",directory_RAW+"/OSZICAR.relax2");vfile.push_back("OSZICAR.relax2");
      }
      if(_APENNSY_STYLE_OLD_) {
        stringstream aus_exec;
        aus_exec << "cd \"" << directory_RAW << "\"" << endl;
        aus_exec << "bin2ascii ./OSZICAR*" << endl << "bin2ascii ./OSZICAR*" << endl;
        aus_exec << "subst \"\\n\\n\" \"\\n\" ./OSZICAR*" << endl << "subst \"\\n\\n\" \"\\n\" ./OSZICAR*" << endl;
        aus_exec << "rm -f ./*~" << endl;
        aurostd::execute(aus_exec);
      }
    }
    // COMPRESS
    // aurostd::execute(DEFAULT_KZIP_BIN+" -9f "+directory_RAW+"/OSZICAR.relax");

    // LINKING
    if(FileName_OUTCAR_relax!="") {
      string FROM=FileName_OUTCAR_relax;
      string TO=FileName_OUTCAR_relax;
      aurostd::StringSubst(TO,"LIB/","RAW/");
      aurostd::StringSubst(TO,"relax1","relax");
      aurostd::StringSubst(TO,"relax2","relax");
      aurostd::StringSubst(TO,"relax3","relax");
      aurostd::StringSubst(TO,"static","relax");
      aurostd::StringSubst(TO,"bands","relax");
      cout << MESSAGE << " linking " << FROM << "->" << TO << endl;
      aurostd::LinkFile(FROM,TO);
      //      exit(0);
    }

    // DONE
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << soliloquy << " - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Magnetic
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_Magnetic(const string& directory_LIB,const string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,const string& MESSAGE) {
    //CO20180130 - note that here we extract STATIC/BANDS properties only
    //some spin properties (spin/cell, spin/atom, spinD) can be extracted from relax2
    //we do this in the thermo loop and DO NOT attempt to redo here
    //spinF is nonzero for magnetic metals, so we need to determine Egap here as well (need BANDS)
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::LIB2RAW_Loop_Magnetic [1]" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Magnetic - begin " << directory_LIB << endl;
    // [OBSOLETE]  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=magnetic";
    data.vloop.push_back("magnetic");
    xOUTCAR outcar_bands,outcar_static;
    xDOSCAR doscar;

    double MAG_EPS=1e-6;

    //try to grab STATIC properties first, Efermi and spin/cell, spin/atom, spinD
    //NB, spin properties from relax2 are okay (lower kppra), we grab it in thermo loop, so we don't attempt to here
    double EFERMI=AUROSTD_NAN;
    outcar_static.clear();
    if(aurostd::FileExist(directory_LIB+"/"+"OUTCAR.static") || aurostd::EFileExist(directory_LIB+"/"+"OUTCAR.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.static",directory_RAW,"OUTCAR.static",vfile,MESSAGE);  // OUTCAR.static
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"OUTCAR.static") << endl;
      if(outcar_static.GetPropertiesFile(directory_RAW+"/"+"OUTCAR.static")) {
        EFERMI=outcar_static.Efermi;
        data.spin_cell=outcar_static.mag_cell;
        data.spin_atom=outcar_static.mag_atom;
        data.spinD="";
        data.vspinD.clear();
        if(outcar_static.vmag.size()) {
          for(uint i=0;i<(uint) outcar_static.vmag.size();i++) {
            data.spinD+=aurostd::utype2string<double>(outcar_static.vmag.at(i),5)+(i<outcar_static.vmag.size()-1?",":"");
            data.vspinD.push_back(outcar_static.vmag.at(i));
          }
        } else {
          for(uint i=0;i<outcar_static.natoms;i++) {  //use outcar_static.natoms as there can be a primitivization between relax and static
            data.spinD+=aurostd::utype2string<double>(0)+(i<outcar_static.natoms-1?",":"");
            data.vspinD.push_back(0.0);
          }
        }
      } else { cout << MESSAGE << " ERROR OUTCAR.static properties cannot be extracted" << endl; }  //CO20200404 - removed xOUTCAR.ERROR
    } else { cout << MESSAGE << " MISSING OUTCAR.static" << endl; }
    if(EFERMI==AUROSTD_NAN) { cout << MESSAGE << " unable to load OUTCAR.static, using Efermi from bands" << endl; }

    //ideally we grab both bands and static
    //we want static because we trust this Efermi (self-consistent)
    //however, we fall back on Efermi from bands

    // BANDS CALCULATION
    data.spinF=0.0;// DEFAULT
    outcar_bands.clear();
    doscar.clear();
    if(aurostd::FileExist(directory_LIB+"/"+"OUTCAR.bands") || aurostd::EFileExist(directory_LIB+"/"+"OUTCAR.bands")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.bands",directory_RAW,"OUTCAR.bands",vfile,MESSAGE);  // OUTCAR.bands
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"OUTCAR.bands") << endl;
      if(outcar_bands.GetPropertiesFile(directory_RAW+"/"+"OUTCAR.bands")) {
        if(outcar_bands.GetBandGap(EFERMI)) {
          data.Egap=outcar_bands.Egap_net;
          data.Egap_fit=outcar_bands.Egap_fit_net;
          data.Egap_type=outcar_bands.Egap_type_net;
          if(aurostd::substring2bool(data.Egap_type,"metal")) data.Egap=0.0;      //half-metal
          if(aurostd::substring2bool(data.Egap_type,"metal")) data.Egap_fit=0.0;  //half-metal

          // SPIN POLARIZATION AT FERMI LEVEL
          if(data.Egap<MAG_EPS && aurostd::abs(data.spin_cell)>MAG_EPS) { //must be metal and magnetic
            if(doscar.content=="" && (aurostd::FileExist(directory_LIB+"/"+"DOSCAR.static") || aurostd::EFileExist(directory_LIB+"/"+"DOSCAR.static"))) {
              aflowlib::LIB2RAW_FileNeeded(directory_LIB,"DOSCAR.static",directory_RAW,"DOSCAR.static",vfile,MESSAGE);  // DOSCAR.static
              if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"DOSCAR.static") << endl;
              doscar.GetPropertiesFile(directory_RAW+"/"+"DOSCAR.static");
            }
            if(doscar.content=="" && (aurostd::FileExist(directory_LIB+"/"+"DOSCAR.relax2") || aurostd::EFileExist(directory_LIB+"/"+"DOSCAR.relax2"))) {
              aflowlib::LIB2RAW_FileNeeded(directory_LIB,"DOSCAR.relax2",directory_RAW,"DOSCAR.relax",vfile,MESSAGE);  // DOSCAR.relax2
              if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"DOSCAR.relax") << endl;
              doscar.GetPropertiesFile(directory_RAW+"/"+"DOSCAR.relax");
            }
            if(!doscar.content.empty()) {
              data.spinF=doscar.spinF;
            } else {
              cout << MESSAGE << " MISSING DOSCAR.static and DOSCAR.relax[2]" << endl;
            }
          }
        } else { cout << MESSAGE << " ERROR OUTCAR.bands BandGap() cannot be extracted" << endl; }  //CO20181129  //CO20200404 - removed xOUTCAR.ERROR
      } else { cout << MESSAGE << " ERROR OUTCAR.bands properties cannot be extracted" << endl; }  //CO20181129 //CO20200404 - removed xOUTCAR.ERROR
    } else { cout << MESSAGE << " MISSING OUTCAR.bands" << endl; }  //CO20181129

    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " spin/cell = "          << data.spin_cell << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " spin/atom = "          << data.spin_atom << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " spinD     = "          << data.spinD << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " spinF     = "          << ((data.spinF   !=AUROSTD_NAN)?aurostd::utype2string(data.spinF,    5):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " Egap (eV)     = "      << ((data.Egap    !=AUROSTD_NAN)?aurostd::utype2string(data.Egap,     5):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " Egap_fit (eV) = " << ((data.Egap_fit!=AUROSTD_NAN)?aurostd::utype2string(data.Egap_fit, 5):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " Egap_type     = "      << ((data.Egap_type.size())?data.Egap_type:"unavailable") << endl;

    // done
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_Magnetic - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_Bader  //CO WORK HERE
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_Bader(const string& directory_LIB,const string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,const string& MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "aflowlib::LIB2RAW_Loop_Bader():";
    if(LDEBUG) cerr << soliloquy << " [1]" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << soliloquy << " - begin " << directory_LIB << endl;
    // [OBSOLETE]  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=bader";
    data.vloop.push_back("bader");

    // [OBSOLETE]   xOUTCAR outcar;
    // [OBSOLETE]aflowlib::LIB2RAW_FileNeeded(directory_LIB,"OUTCAR.static",directory_RAW,"OUTCAR.static",vfile,MESSAGE);  // OUTCAR.static
    // [OBSOLETE]if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"OUTCAR.static") << endl;
    // [OBSOLETE]outcar.GetPropertiesFile(directory_RAW+"/"+"OUTCAR.static");
    if(aurostd::FileExist(directory_LIB+"/"+"CHGCAR.static") || aurostd::EFileExist(directory_LIB+"/"+"CHGCAR.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"CHGCAR.static",directory_RAW,"CHGCAR.static",vfile,MESSAGE);  // CHGCAR.static
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"CHGCAR.static") << endl;
    } else {
      cout << MESSAGE << " MISSING CHGCAR.static" << endl;
      return FALSE;
    }
    if(aurostd::FileExist(directory_LIB+"/"+"AECCAR0.static") || aurostd::EFileExist(directory_LIB+"/"+"AECCAR0.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AECCAR0.static",directory_RAW,"AECCAR0.static",vfile,MESSAGE);  // AECCAR0.static
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"AECCAR0.static") << endl;
    } else {
      cout << MESSAGE << " MISSING AECCAR0.static" << endl;
      return FALSE;
    }
    if(aurostd::FileExist(directory_LIB+"/"+"AECCAR2.static") || aurostd::EFileExist(directory_LIB+"/"+"AECCAR2.static")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AECCAR2.static",directory_RAW,"AECCAR2.static",vfile,MESSAGE);  // AECCAR2.static
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"AECCAR2.static") << endl;
    } else {
      cout << MESSAGE << " MISSING AECCAR2.static" << endl;
      return FALSE;
    }

    if(LDEBUG) cerr << soliloquy << " [2]" << endl;

    vector<string> vspecies;aurostd::string2tokens(data.species,vspecies,",");
    for(uint i=0;i<vspecies.size();i++)
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << vspecies.at(i) << endl;

    if(LDEBUG) cerr << soliloquy << " [3]" << endl;

    vector<double> vspecies_pp_ZVAL;
    for(uint i=0;i<data.vspecies_pp_ZVAL.size();i++) { vspecies_pp_ZVAL.push_back(data.vspecies_pp_ZVAL.at(i)); }//;aurostd::string2tokens(data.species_pp_ZVAL,vspecies_pp_ZVAL,",");
    for(uint i=0;i<vspecies_pp_ZVAL.size();i++) { if(AFLOWLIB_VERBOSE) { cout << MESSAGE << " " << vspecies_pp_ZVAL.at(i) << endl; } } 

    if(LDEBUG) cerr << soliloquy << " [4]" << endl;

    deque<int> num_each_type;
    //aurostd::string2tokens<int>(data.composition,num_each_type,",");  //CO20180626 - BAD! relax2 can be different than static via AFLOW primitization
    xstructure xstr_relaxed = KBIN::GetMostRelaxedStructure(directory_LIB);
    for(uint i=0;i<xstr_relaxed.num_each_type.size();i++) { num_each_type.push_back(xstr_relaxed.num_each_type[i]);}  //vector vs. deque
    if(LDEBUG) {for(uint i=0;i<num_each_type.size();i++) { cerr << XPID << soliloquy << " num_each_type[" << i << "]=" << num_each_type[i] << endl;}}

    vector<double> vZVAL;
    for(uint i=0;i<vspecies_pp_ZVAL.size();i++) {
      for(uint j=0;j<(uint)num_each_type.at(i);j++) {
        vZVAL.push_back(vspecies_pp_ZVAL.at(i));
        //vZVAL.push_back(aurostd::string2utype<double>(vspecies_pp_ZVAL.at(i)));
      }
    }

    if(LDEBUG) cerr << soliloquy << " [5]" << endl;

    //flags
    aurostd::xoption bader_flags;
    //DX+CO START
    bader_flags.flag("BADER::AFLOWLIB_LIBRARY",TRUE); //net charge file
    bader_flags.flag("BADER::JVXL_ALL_SPECIES",TRUE); //for flag automation, no need to specify cutoffs,downsample_ratios here
    bader_flags.flag("BADER::JVXL_CYCLIC",TRUE);      //CYCLIC MODE
    //bader_flags.flag("BADER::JVXL_CYCLIC",FALSE);   //SETS MODE
    //DX+CO END

    vector<double> cutoffs;
    vector<int> downsample_ratios;
    cutoffs.push_back(0.10);
    cutoffs.push_back(0.20);
    cutoffs.push_back(0.25);
    cutoffs.push_back(0.30);
    cutoffs.push_back(0.333333333333333);
    cutoffs.push_back(0.40);
    cutoffs.push_back(0.50);
    cutoffs.push_back(0.60);
    cutoffs.push_back(0.666666666666666);
    cutoffs.push_back(0.75);
    cutoffs.push_back(0.80);
    cutoffs.push_back(0.90);
    //downsample_ratios.push_back(1);
    downsample_ratios.push_back(2);   //for all
    /*if(data.natoms<=7) {              //for 0.5 and 0.75, downsample_ratio depends on natoms
      downsample_ratios.push_back(1);
    //downsample_ratios.push_back(1);
    } else {                          //8 or more
    downsample_ratios.push_back(2);
    //downsample_ratios.push_back(2);
    }
    //downsample_ratios.push_back(1);*/

    if(LDEBUG) cerr << soliloquy << " [6]" << endl;

    string bader_options;
    ostringstream oss;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " beginning BADER calculation--please be patient" << endl;
    bader_functions::Flags2BaderCommands(bader_flags,bader_options,oss);
    //DX+CO START
    bader_functions::BaderCalc(bader_flags,bader_options,KBIN::ExtractSystemName(directory_LIB),vspecies,num_each_type,vZVAL,cutoffs,downsample_ratios,directory_RAW,oss);  //CO+ME20200601 - no data.prototype, use ExtractSystemName() //data.prototype
    //DX+CO END
    if(AFLOWLIB_VERBOSE) cout << oss.str();

    if(LDEBUG) cerr << soliloquy << " [7]" << endl;

    vector<string> vline,tokens;
    stringstream abader_ss;
    string abader_out=data.prototype+"_abader.out";

    abader_ss.clear();
    if(aurostd::FileExist(directory_RAW+"/"+abader_out)) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+abader_out) << endl;
      aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(directory_RAW+"/"+abader_out),abader_ss,"[BADER_RESULTS]START","[BADER_RESULTS]STOP");
      aurostd::stream2vectorstring(abader_ss,vline);
      for (uint i=0;i<vline.size();i++) {
        aurostd::StringSubst(vline.at(i),"="," ");
        aurostd::string2tokens(vline.at(i),tokens," ");
        if(tokens.size()>=2) {
          if(tokens.at(0)=="bader_net_charges") data.bader_net_charges=tokens.at(1);
          if(tokens.at(0)=="bader_atomic_volumes") data.bader_atomic_volumes=tokens.at(1);
        }
      }
      aurostd::string2tokens<double>(data.bader_net_charges,data.vbader_net_charges,",");         //be careful, sets precision to 20
      aurostd::string2tokens<double>(data.bader_atomic_volumes,data.vbader_atomic_volumes,",");   //be careful, sets precision to 20
      if(data.vbader_net_charges.size()==0) {
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " data.vbader_net_charges was not extracted correctly" << endl;
        return FALSE;
      }
      if(data.vbader_net_charges.size()!=data.vbader_atomic_volumes.size()) {
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " length of data.vbader_net_charges does not match length of data.vbader_atomic_volumes" << endl;
        if(AFLOWLIB_VERBOSE) cout << MESSAGE << " there was a problem with the bader data extraction" << endl;
        return FALSE;
      }
      stringstream num_prec;
      if(AFLOWLIB_VERBOSE) {
        for(uint i=0;i<data.vbader_net_charges.size();i++) {
          num_prec << std::fixed << setprecision(4) << data.vbader_net_charges.at(i);
          cout << MESSAGE << " bader net charge for atom " << i+1 << " = " << num_prec.str() << endl;
          num_prec.str("");
        }
        for(uint i=0;i<data.vbader_atomic_volumes.size();i++) {
          num_prec << std::fixed << setprecision(4) << data.vbader_atomic_volumes.at(i);
          cout << MESSAGE << " bader atomic volume for atom " << i+1 << " = " << num_prec.str() << endl;
          num_prec.str("");
        }
      }
    } else {
      return FALSE;
    }
    // done
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << soliloquy << " - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_AGL  // CORMAC
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_AGL(const string& directory_LIB,const string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,const string& MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::LIB2RAW_Loop_AGL [1]" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_AGL - begin " << directory_LIB << endl;
    // [OBSOLETE]  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=agl";
    data.vloop.push_back("agl");

    vector<string> vline,tokens;
    stringstream aflow_agl_out;

    // AFLOW AGL
    aflow_agl_out.clear();
    if(aurostd::FileExist(directory_LIB+"/"+"aflow.agl.out") || aurostd::EFileExist(directory_LIB+"/"+"aflow.agl.out")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.agl.out",directory_RAW,"aflow.agl.out",vfile,MESSAGE);  // aflow.agl.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL.out",directory_RAW,"AGL.out",vfile,MESSAGE);  // AGL.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_energies_temperature.out",directory_RAW,"AGL_energies_temperature.out",vfile,MESSAGE);  // AGL.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_thermal_properties_temperature.out",directory_RAW,"AGL_thermal_properties_temperature.out",vfile,MESSAGE);  // AGL.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_edos_gap_pressure.out",directory_RAW,"AGL_edos_gap_pressure.out",vfile,MESSAGE);  // AGL_edos_gap_pressure.out //CT20181212
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_edos_gap_pressure.json",directory_RAW,"AGL_edos_gap_pressure.json",vfile,MESSAGE);  // AGL_edos_gap_pressure.json //CT20181212
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_energy.json",directory_RAW,"AGL_energy.json",vfile,MESSAGE);  // AGL_energy.out //CT20181212
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_energy_structures.json",directory_RAW,"AGL_energy_structures.json",vfile,MESSAGE);  // AGL_energy_structures.json //CT20181212
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_energy_volume.out",directory_RAW,"AGL_energy_volume.out",vfile,MESSAGE);  // AGL_energy_volume.out //CT20181212
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_gibbs_energy_pT.out",directory_RAW,"AGL_gibbs_energy_pT.out",vfile,MESSAGE);  // AGL_gibbs_energy_pT.out //CT20181212
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AGL_Hugoniot.out",directory_RAW,"AGL_Hugoniot.out",vfile,MESSAGE);  // AGL_Hugoniot.out //CT20181212
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"aflow.agl.out") << endl;
      aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(directory_RAW+"/"+"aflow.agl.out"),aflow_agl_out,"[AGL_RESULTS]START","[AGL_RESULTS]STOP");
      aurostd::stream2vectorstring(aflow_agl_out,vline);
      for (uint i=0;i<vline.size();i++) {
        aurostd::StringSubst(vline.at(i),"="," ");
        aurostd::string2tokens(vline.at(i),tokens," ");
        if(tokens.size()>=2) {
          if(tokens.at(0)=="agl_thermal_conductivity_300K") data.agl_thermal_conductivity_300K=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="agl_debye") data.agl_debye=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="agl_acoustic_debye") data.agl_acoustic_debye=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="agl_gruneisen") data.agl_gruneisen=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="agl_heat_capacity_Cv_300K") data.agl_heat_capacity_Cv_300K=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="agl_heat_capacity_Cp_300K") data.agl_heat_capacity_Cp_300K=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="agl_thermal_expansion_300K") data.agl_thermal_expansion_300K=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="agl_bulk_modulus_static_300K") data.agl_bulk_modulus_static_300K=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="agl_bulk_modulus_isothermal_300K") data.agl_bulk_modulus_isothermal_300K=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="agl_poisson_ratio_source") data.agl_poisson_ratio_source=tokens.at(1); //CT20181212
          if(tokens.at(0)=="agl_vibrational_free_energy_300K_cell") data.agl_vibrational_free_energy_300K_cell=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
          if(tokens.at(0)=="agl_vibrational_free_energy_300K_atom") data.agl_vibrational_free_energy_300K_atom=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
          if(tokens.at(0)=="agl_vibrational_entropy_300K_cell") data.agl_vibrational_entropy_300K_cell=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
          if(tokens.at(0)=="agl_vibrational_entropy_300K_atom") data.agl_vibrational_entropy_300K_atom=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
        }
      }
    } else {
      return FALSE;
    }

    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_thermal_conductivity_300K (W/m*K) = " << ((data.agl_thermal_conductivity_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_thermal_conductivity_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_debye (K) = " << ((data.agl_debye!=AUROSTD_NAN)?aurostd::utype2string(data.agl_debye,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_acoustic_debye (K) = " << ((data.agl_acoustic_debye!=AUROSTD_NAN)?aurostd::utype2string(data.agl_acoustic_debye,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_gruneisen = " << ((data.agl_gruneisen!=AUROSTD_NAN)?aurostd::utype2string(data.agl_gruneisen,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_heat_capacity_Cv_300K (kB/cell) = " << ((data.agl_heat_capacity_Cv_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_heat_capacity_Cv_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_heat_capacity_Cp_300K (kB/cell) = " << ((data.agl_heat_capacity_Cp_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_heat_capacity_Cp_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_thermal_expansion_300K (1/K) = " << ((data.agl_thermal_expansion_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_thermal_expansion_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_bulk_modulus_static_300K (GPa) = " << ((data.agl_bulk_modulus_static_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_bulk_modulus_static_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_bulk_modulus_isothermal_300K (GPa) = " << ((data.agl_bulk_modulus_isothermal_300K!=AUROSTD_NAN)?aurostd::utype2string(data.agl_bulk_modulus_isothermal_300K,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_poisson_ratio_source = " << ((data.agl_poisson_ratio_source.size())?(data.agl_poisson_ratio_source):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_vibrational_free_energy_300K_cell (meV/cell) = " << ((data.agl_vibrational_free_energy_300K_cell!=AUROSTD_NAN)?aurostd::utype2string(data.agl_vibrational_free_energy_300K_cell,10):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_vibrational_free_energy_300K_atom (meV/atom) = " << ((data.agl_vibrational_free_energy_300K_atom!=AUROSTD_NAN)?aurostd::utype2string(data.agl_vibrational_free_energy_300K_atom,10):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_vibrational_entropy_300K_cell (meV/cell*K) = " << ((data.agl_vibrational_entropy_300K_cell!=AUROSTD_NAN)?aurostd::utype2string(data.agl_vibrational_entropy_300K_cell,10):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " agl_vibrational_entropy_300K_atom (meV/atom*K) = " << ((data.agl_vibrational_entropy_300K_atom!=AUROSTD_NAN)?aurostd::utype2string(data.agl_vibrational_entropy_300K_atom,10):"unavailable") << endl; //CT20181212
    // done
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_AGL - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_AEL  // CORMAC
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_AEL(const string& directory_LIB,const string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,const string& MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function = XPID + "aflowlib::LIB2RAW_loop_AEL()";  //ME20191119
    if(LDEBUG) cerr << "aflowlib::LIB2RAW_Loop_AEL [1]" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_AEL - begin " << directory_LIB << endl;
    // [OBSOLETE]  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=ael";
    data.vloop.push_back("ael");

    vector<string> vline,tokens;
    stringstream aflow_ael_out;

    // AFLOW AEL
    aflow_ael_out.clear();
    if(aurostd::FileExist(directory_LIB+"/"+"aflow.ael.out") || aurostd::EFileExist(directory_LIB+"/"+"aflow.ael.out")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.ael.out",directory_RAW,"aflow.ael.out",vfile,MESSAGE);  // aflow.ael.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AEL_Elastic_constants.out",directory_RAW,"AEL_Elastic_constants.out",vfile,MESSAGE);  // AEL_Elastic_constants.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AEL_Compliance_tensor.out",directory_RAW,"AEL_Compliance_tensor.out",vfile,MESSAGE);  // AEL_Compliance_tensor.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AEL_elastic_tensor.json",directory_RAW,"AEL_elastic_tensor.json",vfile,MESSAGE);  // AEL_elastic_tensor.json //CT20181212
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"AEL_energy_structures.json",directory_RAW,"AEL_energy_structures.json",vfile,MESSAGE);  // AEL_energy_structures.json //CT20181212
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " loading " << string(directory_RAW+"/"+"aflow.ael.out") << endl;
      //ME20191105
      string ael_out_str = aurostd::efile2string(directory_RAW + "/aflow.ael.out");
      aurostd::ExtractToStringstreamEXPLICIT(ael_out_str, aflow_ael_out, "[AEL_RESULTS]START", "[AEL_RESULTS]STOP");
      //aurostd::ExtractToStringstreamEXPLICIT(aurostd::efile2string(directory_RAW+"/"+"aflow.ael.out"),aflow_ael_out,"[AEL_RESULTS]START","[AEL_RESULTS]STOP");  OBSOLETE ME20191105
      aurostd::stream2vectorstring(aflow_ael_out,vline);
      for (uint i=0;i<vline.size();i++) {
        aurostd::StringSubst(vline.at(i),"="," ");
        aurostd::string2tokens(vline.at(i),tokens," ");
        if(tokens.size()>=2) {
          if(tokens.at(0)=="ael_poisson_ratio") data.ael_poisson_ratio=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="ael_bulk_modulus_voigt" || tokens.at(0)=="ael_bulk_modulus_voight") data.ael_bulk_modulus_voigt=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="ael_bulk_modulus_reuss") data.ael_bulk_modulus_reuss=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="ael_shear_modulus_voigt" || tokens.at(0)=="ael_shear_modulus_voigth") data.ael_shear_modulus_voigt=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="ael_shear_modulus_reuss") data.ael_shear_modulus_reuss=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="ael_bulk_modulus_vrh") data.ael_bulk_modulus_vrh=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="ael_shear_modulus_vrh") data.ael_shear_modulus_vrh=aurostd::string2utype<double>(tokens.at(1));
          if(tokens.at(0)=="ael_elastic_anisotropy") data.ael_elastic_anisotropy=aurostd::string2utype<double>(tokens.at(1)); //CO20181129
          if(tokens.at(0)=="ael_youngs_modulus_vrh") data.ael_youngs_modulus_vrh=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
          if(tokens.at(0)=="ael_speed_sound_transverse") data.ael_speed_sound_transverse=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
          if(tokens.at(0)=="ael_speed_sound_longitudinal") data.ael_speed_sound_longitudinal=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
          if(tokens.at(0)=="ael_speed_sound_average") data.ael_speed_sound_average=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
          if(tokens.at(0)=="ael_pughs_modulus_ratio") data.ael_pughs_modulus_ratio=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
          if(tokens.at(0)=="ael_debye_temperature") data.ael_debye_temperature=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
          if(tokens.at(0)=="ael_applied_pressure") data.ael_applied_pressure=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
          if(tokens.at(0)=="ael_average_external_pressure") data.ael_average_external_pressure=aurostd::string2utype<double>(tokens.at(1)); //CT20181212
        }
      }

      //ME20191105 BEGIN
      xmatrix<double> tensor(6, 6);
      vector<double> row;
      aurostd::ExtractToStringstreamEXPLICIT(ael_out_str, aflow_ael_out, "[AEL_STIFFNESS_TENSOR]START", "[AEL_STIFFNESS_TENSOR]STOP");
      aurostd::stream2vectorstring(aflow_ael_out, vline);
      if (vline.size() > 0) {
        vline.pop_back();  // Remove empty line at the end
        try {
          if (vline.size() != 6) {
            stringstream message;
            message << "Could not read stiffness tensor: wrong number of lines"
              << " (found " << vline.size() << ", need 6).";
            throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
          }
          for (int i = 0; i < 6; i++) {
            aurostd::string2tokens(vline[i], row);
            if (row.size() != 6) {
              stringstream message;
              message << "Could not read stiffness tensor."
                << " Wrong number of columns in line " << (i + 1)
                << " (found " << row.size() << ", need 6).";
              throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
            }
            for (int j = 0; j < 6; j++) tensor[i + 1][j + 1] = row[j];
          }
          data.ael_stiffness_tensor = tensor;
        } catch (aurostd::xerror& e) {
          std::cout << MESSAGE << " ERROR - " << e.error_message << std::endl;
        }
      } else {
        std::cout << MESSAGE << " WARNING - No stiffness tensor found in aflow.ael.out." << std::endl;
      }

      tensor.clear();
      aurostd::ExtractToStringstreamEXPLICIT(ael_out_str, aflow_ael_out, "[AEL_COMPLIANCE_TENSOR]START", "[AEL_COMPLIANCE_TENSOR]STOP");
      aurostd::stream2vectorstring(aflow_ael_out, vline);
      if (vline.size() > 0) {
        vline.pop_back();  // Remove empty line at the end
        try {
          if (vline.size() != 6) {
            string function = "aflowlib::LIB2RAW_loop_AEL()";
            stringstream message;
            message << "Could not read compliance tensor: wrong number of lines"
              << " (found " << vline.size() << ", need 6).";
            throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
          }
          for (int i = 0; i < 6; i++) {
            aurostd::string2tokens(vline[i], row);
            if (row.size() != 6) {
              string function = "aflowlib::LIB2RAW_loop_AEL()";
              stringstream message;
              message << "Could not read compliance tensor:"
                << " wrong number of columns in line " << (i + 1)
                << " (found " << row.size() << ", need 6).";
              throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
            }
            for (int j = 0; j < 6; j++) tensor[i + 1][j + 1] = row[j];
          }
          data.ael_compliance_tensor = tensor;
        } catch (aurostd::xerror& e) {
          std::cout << MESSAGE << " ERROR - " << e.error_message << std::endl;
        }
      } else {
        std::cout << MESSAGE << " WARNING - No compliance tensor found in aflow.ael.out." << std::endl;
      }
      //ME20191105 END
    } else {
      return FALSE;
    }

    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_poisson_ratio = " << ((data.ael_poisson_ratio!=AUROSTD_NAN)?aurostd::utype2string(data.ael_poisson_ratio,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_bulk_modulus_voigt (GPa) = " << ((data.ael_bulk_modulus_voigt!=AUROSTD_NAN)?aurostd::utype2string(data.ael_bulk_modulus_voigt,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_bulk_modulus_reuss (GPa) = " << ((data.ael_bulk_modulus_reuss!=AUROSTD_NAN)?aurostd::utype2string(data.ael_bulk_modulus_reuss,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_shear_modulus_voigt (GPa) = " << ((data.ael_shear_modulus_voigt!=AUROSTD_NAN)?aurostd::utype2string(data.ael_shear_modulus_voigt,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_shear_modulus_reuss (GPa) = " << ((data.ael_shear_modulus_reuss!=AUROSTD_NAN)?aurostd::utype2string(data.ael_shear_modulus_reuss,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_bulk_modulus_vrh (GPa) = " << ((data.ael_bulk_modulus_vrh!=AUROSTD_NAN)?aurostd::utype2string(data.ael_bulk_modulus_vrh,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_shear_modulus_vrh (GPa) = " << ((data.ael_shear_modulus_vrh!=AUROSTD_NAN)?aurostd::utype2string(data.ael_shear_modulus_vrh,10):"unavailable") << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_elastic_anisotropy = " << ((data.ael_elastic_anisotropy!=AUROSTD_NAN)?aurostd::utype2string(data.ael_elastic_anisotropy,10):"unavailable") << endl; //CO20181129
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_youngs_modulus_vrh (GPa) = " << ((data.ael_youngs_modulus_vrh!=AUROSTD_NAN)?aurostd::utype2string(data.ael_youngs_modulus_vrh,10):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_speed_sound_transverse (m/s) = " << ((data.ael_speed_sound_transverse!=AUROSTD_NAN)?aurostd::utype2string(data.ael_speed_sound_transverse,10):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_speed_sound_longitudinal (m/s) = " << ((data.ael_speed_sound_longitudinal!=AUROSTD_NAN)?aurostd::utype2string(data.ael_speed_sound_longitudinal,10):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_speed_sound_average (m/s) = " << ((data.ael_speed_sound_average!=AUROSTD_NAN)?aurostd::utype2string(data.ael_speed_sound_average,10):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_pughs_modulus_ratio = " << ((data.ael_pughs_modulus_ratio!=AUROSTD_NAN)?aurostd::utype2string(data.ael_pughs_modulus_ratio,10):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_debye_temperature (K) = " << ((data.ael_debye_temperature!=AUROSTD_NAN)?aurostd::utype2string(data.ael_debye_temperature,10):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_applied_pressure (GPa) = " << ((data.ael_applied_pressure!=AUROSTD_NAN)?aurostd::utype2string(data.ael_applied_pressure,10):"unavailable") << endl; //CT20181212
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " ael_average_external_pressure (GPa) = " << ((data.ael_average_external_pressure!=AUROSTD_NAN)?aurostd::utype2string(data.ael_average_external_pressure,10):"unavailable") << endl; //CT20181212
    // done
    if (AFLOWLIB_VERBOSE) {  //ME20191105
      std::cout << MESSAGE << " ael_stiffness_tensor = ";
      if ((data.ael_stiffness_tensor.rows != 6) || (data.ael_stiffness_tensor.cols != 6)) std::cout << "unavailable" << std::endl;
      else std::cout << std::endl << data.ael_stiffness_tensor << std::endl;

      std::cout << MESSAGE << " ael_compliance_tensor = ";
      if ((data.ael_compliance_tensor.rows != 6) || (data.ael_compliance_tensor.cols != 6)) std::cout << "unavailable" << std::endl;
      else std::cout << std::endl << data.ael_compliance_tensor << std::endl;
    }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_AEL - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_POCC - CO20200624
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_POCC(const string& directory_LIB,const string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,const string& MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"aflowlib::LIB2RAW_Loop_POCC():";
    if(LDEBUG) cerr << soliloquy << " [1]" << endl;
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << soliloquy << " - begin " << directory_LIB << endl;
    data.vloop.push_back("pocc");
    
    if(aurostd::EFileExist(directory_LIB+"/"+POCC_FILE_PREFIX+POCC_OUT_FILE)) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,POCC_FILE_PREFIX+POCC_OUT_FILE,directory_RAW,POCC_FILE_PREFIX+POCC_OUT_FILE,vfile,MESSAGE);  // aflow.pocc.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE,directory_RAW,POCC_FILE_PREFIX+POCC_UNIQUE_SUPERCELLS_FILE,vfile,MESSAGE);  // aflow.pocc.structures_unique.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,POCC_FILE_PREFIX+POCC_ALL_HNF_MATRICES_FILE,directory_RAW,POCC_FILE_PREFIX+POCC_ALL_HNF_MATRICES_FILE,vfile,MESSAGE);  // aflow.pocc.hnf_matrices.out
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,POCC_FILE_PREFIX+POCC_ALL_SITE_CONFIGURATIONS_FILE,directory_RAW,POCC_FILE_PREFIX+POCC_ALL_SITE_CONFIGURATIONS_FILE,vfile,MESSAGE);  // aflow.pocc.site_configurations.out
      //get DOSCAR.pocc + png's
      vector<string> vfiles;
      aurostd::DirectoryLS(directory_LIB,vfiles);
      for(uint i=0;i<vfiles.size();i++){
        if(vfiles[i].find(POCC_DOSCAR_FILE)!=string::npos){
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,vfiles[i],directory_RAW,vfiles[i],vfile,MESSAGE);  // DOSCAR.pocc_
        }
        if(vfiles[i].find("dos_orbitals_")!=string::npos && vfiles[i].find(".png")!=string::npos){
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,vfiles[i],directory_RAW,vfiles[i],vfile,MESSAGE);  // dos_orbitals_
        }
        if(vfiles[i].find("dos_species_")!=string::npos && vfiles[i].find(".png")!=string::npos){
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,vfiles[i],directory_RAW,vfiles[i],vfile,MESSAGE);  // dos_species_
        }
        if(vfiles[i].find("dos_atoms_")!=string::npos && vfiles[i].find(".png")!=string::npos){
          aflowlib::LIB2RAW_FileNeeded(directory_LIB,vfiles[i],directory_RAW,vfiles[i],vfile,MESSAGE);  // dos_atoms_
        }
      }
    }
    if(aurostd::EFileExist(directory_LIB+"/"+"aflow.pocc_agl.out")) {
      aflowlib::LIB2RAW_FileNeeded(directory_LIB,"aflow.pocc_agl.out",directory_RAW,"aflow.pocc_agl.out",vfile,MESSAGE);  // aflow.pocc_agl.out
    }

    if(data.vspecies.size()==0){
    int comp_prec=(int)ceil(log10(1.0/DEFAULT_POCC_STOICH_TOL));  //ceil ensures we round up above 1 //CO20181226
    }

    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << soliloquy << " - end " << directory_LIB << endl;
    return true;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_LOCK
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_LOCK(const string& directory_LIB,const string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,const string& MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"aflowlib::LIB2RAW_Loop_LOCK():";
    if(LDEBUG) cerr << soliloquy << " [1]" << endl;
    // Stefano Curtarolo 2009-2010-2011-2012
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << soliloquy << " - begin " << directory_LIB << endl;
    //  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=LOCK";
    data.vloop.push_back("lock");

    vector<string> vlock,tokens;
    string tmp="";
    string::size_type loc;
    aflowlib::LIB2RAW_FileNeeded(directory_LIB,_AFLOWLOCK_,directory_RAW,_AFLOWLOCK_,vfile,MESSAGE);  // OUTCAR.static
    aurostd::file2vectorstring(directory_RAW+"/"+_AFLOWLOCK_,vlock) ;
    _XHOST aus_XHOST;
    // ---------------------------------------------------------------
    data.vaflowlib_date.clear(); //clear here, 
    for(uint iline=0;iline<vlock.size();iline++) //CO20200624 - adding lock date to vaflowlib_date  //grab LAST date - data.vaflowlib_date.empty()
      if(aurostd::substring2bool(vlock[iline],"date=") && aurostd::substring2bool(vlock[iline],"[") && aurostd::substring2bool(vlock[iline],"]")) {
        loc=vlock[iline].find("date=");
        tmp=vlock[iline].substr(loc,string::npos);
        loc=tmp.find("]");
        tmp=tmp.substr(0,loc);
        tmp=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(tmp);
        aurostd::StringSubst(tmp,"[","");aurostd::StringSubst(tmp,"]","");  //just in case
        aurostd::StringSubst(tmp,"date=",""); //just in case
        tmp=aflow_convert_time_ctime2aurostd(tmp);
        if(!tmp.empty()){data.vaflowlib_date.clear();data.vaflowlib_date.push_back(tmp+"_GMT-5");}
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " lock_date = " << ((data.vaflowlib_date.size())?data.vaflowlib_date.front():"unavailable") << endl;
    // ---------------------------------------------------------------
    for(uint iline=0;iline<vlock.size()&&data.aflow_version.empty();iline++)
      if(aurostd::substring2bool(vlock[iline],"NFS") && aurostd::substring2bool(vlock[iline],"(") && aurostd::substring2bool(vlock[iline],")")) {
        aurostd::string2tokens(vlock[iline],tokens);
        data.aflow_version=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1));
        aurostd::StringSubst(data.aflow_version,"(","");aurostd::StringSubst(data.aflow_version,")","");
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflow_version = " << ((data.aflow_version.size())?data.aflow_version:"unavailable") << endl;
    // XHOST.CPU_Model ---------------------------------------------------------------
    for(uint iline=0;iline<vlock.size()&&aus_XHOST.CPU_Model.empty();iline++)
      if(aurostd::substring2bool(vlock[iline],"XHOST.CPU_Model") && aurostd::substring2bool(vlock[iline],":")) {
        aurostd::string2tokens(vlock[iline],tokens,":");
        aus_XHOST.CPU_Model=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1));
        data.node_CPU_Model=aus_XHOST.CPU_Model;
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " XHOST.CPU_Model = " << ((aus_XHOST.CPU_Model.size())?aus_XHOST.CPU_Model:"unavailable") << endl;
    // XHOST.CPU_Cores ---------------------------------------------------------------
    for(uint iline=0;iline<vlock.size()&&aus_XHOST.CPU_Cores==0;iline++)
      if(aurostd::substring2bool(vlock[iline],"XHOST.CPU_Cores") && aurostd::substring2bool(vlock[iline],":")) {
        aurostd::string2tokens(vlock[iline],tokens,":");
        aus_XHOST.CPU_Cores=aurostd::string2utype<int>(aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)));
        aus_XHOST.CPU_Cores=aurostd::string2utype<int>(aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1)));
        data.node_CPU_Cores=aus_XHOST.CPU_Cores;
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " XHOST.CPU_Cores = " << ((aus_XHOST.CPU_Cores)?aurostd::utype2string<int>(aus_XHOST.CPU_Cores):"unavailable") << endl;
    // XHOST.CPU_MHz ---------------------------------------------------------------
    for(uint iline=0;iline<vlock.size()&&aus_XHOST.CPU_MHz.empty();iline++)
      if(aurostd::substring2bool(vlock[iline],"XHOST.CPU_MHz") && aurostd::substring2bool(vlock[iline],":")) {
        aurostd::string2tokens(vlock[iline],tokens,":");
        aus_XHOST.CPU_MHz=aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1));
        data.node_CPU_MHz=ceil(aurostd::string2utype<double>(aus_XHOST.CPU_MHz));
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " XHOST.CPU_MHz = " << ((aus_XHOST.CPU_MHz.size())?aus_XHOST.CPU_MHz:"unavailable") << endl;
    // XHOST.RAM_GB ---------------------------------------------------------------
    for(uint iline=0;iline<vlock.size()&&aus_XHOST.RAM_GB<0.001;iline++)
      if(aurostd::substring2bool(vlock[iline],"XHOST.RAM_GB") && aurostd::substring2bool(vlock[iline],":")) {
        aurostd::string2tokens(vlock[iline],tokens,":");
        aus_XHOST.RAM_GB=ceil(aurostd::string2utype<double>(aurostd::RemoveWhiteSpaces(tokens.at(tokens.size()-1))));
        data.node_RAM_GB=aus_XHOST.RAM_GB;
      }
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " XHOST.RAM_GB = " << ((aus_XHOST.RAM_GB)?aurostd::utype2string<double>(aus_XHOST.RAM_GB):"unavailable") << endl;
    // REMOVING  ---------------------------------------------------------------
    vector<string> vremove,vremovej;
    aurostd::string2tokens("aflow.in~,WAVECAR,REPORT,POTCAR.relax1,POTCAR.relax2,POTCAR.relax3,POTCAR.static,POTCAR.bands,AECCAR0,AECCAR1,AECCAR2,AECCAR1.static,AECCAR0.bands,AECCAR1.bands,AECCAR2.bands,AECCAR0.relax1,AECCAR1.relax1,AECCAR2.relax1,AECCAR0.relax2,AECCAR1.relax2,AECCAR2.relax2",vremove,",");
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      for(uint iremove=0;iremove<vremove.size();iremove++) {
        vremovej.clear();
        aurostd::string2vectorstring(aurostd::execute2string("find \""+directory_LIB+"\" -name \""+vremove.at(iremove)+XHOST.vext.at(iext)+"\""),vremovej);
        for(uint j=0;j<vremovej.size();j++) {
          if(aurostd::FileExist(vremovej.at(j))) {
            aurostd::RemoveFile(vremovej.at(j));
            if(AFLOWLIB_VERBOSE) 
              cout << MESSAGE << " removing file: " << aurostd::CleanFileName(vremovej.at(j)) << endl;
          }
        }
      } // iremove
    } // iext

    vector<string> vremoveERROR;
    aurostd::string2tokens("aflow.error.rotmat,aflow.error.nbands,aflow.error.symprec,aflow.error.read_kpoints_rd_sym,aflow.error.ibzkpt,aflow.error.gamma_shift,aflow.error.mpich11,aflow.error.mpich139,aflow.error.nkxyz_ikptd,aflow.error.eddrmm,aflow.error.lreal,aflow.error.exccor,aflow.error.brmix,aflow.error.dav,aflow.error.edddav,aflow.error.efield_pead,aflow.error.zpotrf,aflow.error.zpotrf_potim,aflow.error.natoms,aflow.error.psmaxn,aflow.error.npar,aflow.error.nparc,aflow.error.nparn,aflow.error.npar_remove,aflow.error.csloshing,aflow.error.dentet",vremoveERROR,",");
    for(uint iext=0;iext<XHOST.vext.size();iext++) {
      for(uint iremove=0;iremove<vremoveERROR.size();iremove++) {
        if(aurostd::FileExist(directory_LIB+"/"+vremoveERROR.at(iremove)+XHOST.vext.at(iext))) {
          aurostd::RemoveFile(directory_LIB+"/"+vremoveERROR.at(iremove)+XHOST.vext.at(iext));
        }
      } // iremove
    } // iext

    // done   ---------------------------------------------------------------
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " " << soliloquy << " - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::LIB2RAW_Loop_PATCH
// ***************************************************************************
namespace aflowlib {
  bool LIB2RAW_Loop_PATCH(const string& directory_LIB,const string& directory_RAW,vector<string> &vfile,aflowlib::_aflowlib_entry& data,const string& MESSAGE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::LIB2RAW_Loop_PATCH [1]" << endl;
    // Stefano Curtarolo 2009-2010-2011-2012
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_PATCH - begin " << directory_LIB << endl;
    if(AFLOWLIB_VERBOSE&&0) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_PATCH - begin " << directory_RAW << endl;  // to avoid unused
    if(AFLOWLIB_VERBOSE&&0) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_PATCH - begin " << vfile.size() << endl;  // to avoid unused
    //  aflowlib_out << _AFLOWLIB_ENTRY_SEPARATOR_ << "loop=PATCH";

    // --------------------------------------------------------------------
    // CHECK existence of DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT
    string file="";
    for(uint iext=1;iext<XHOST.vext.size()&&file.empty();iext++) // SKIM uncompressed
      if(aurostd::FileExist(directory_LIB+"/"+DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT+XHOST.vext.at(iext)))
        file=directory_LIB+"/"+DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT+XHOST.vext.at(iext);

    if(file!="") {
      // if(AFLOWLIB_VERBOSE)
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " FOUND " << file << endl;
    } else {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " WRITING " << directory_LIB << "/" << DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT << endl;
      KBIN::VASP_Write_ppAUID_FILE(directory_LIB,data.vspecies_pp_AUID,data.vspecies);
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " COMPRESSING " << directory_LIB << "/" << DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT << endl;
      aurostd::CompressFile(directory_LIB+"/"+DEFAULT_AFLOW_PSEUDOPOTENTIAL_AUID_OUT,DEFAULT_KZIP_BIN);
    }
    // --------------------------------------------------------------------
    // CHECK existence of VASP_POTCAR_AUID in aflow.in
    vector<string> vVASP_POTCAR_AUID;
    aurostd::string2vectorstring(aurostd::execute2string("cat "+directory_LIB+"/"+_AFLOWIN_+" | grep VASP_POTCAR_AUID"),vVASP_POTCAR_AUID);
    if(LDEBUG) cout << XPID << "aflowlib::LIB2RAW_Loop_PATCH: vVASP_POTCAR_AUID.size()=" << vVASP_POTCAR_AUID.size() << endl;

    if(vVASP_POTCAR_AUID.size()) {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " FOUND [VASP_POTCAR_AUID] in " << _AFLOWIN_ << ": " << vVASP_POTCAR_AUID.at(0) << endl;
    } else {
      if(AFLOWLIB_VERBOSE) cout << MESSAGE << " FIXING " << directory_LIB << "/" << _AFLOWIN_ << endl;
      KBIN::VASP_Write_ppAUID_AFLOWIN(directory_LIB,data.vspecies_pp_AUID,data.vspecies);
    } 
    // done   ---------------------------------------------------------------
    if(AFLOWLIB_VERBOSE) cout << MESSAGE << " aflowlib::LIB2RAW_Loop_PATCH - end " << directory_LIB << endl;
    return TRUE;
  }
}

// ***************************************************************************
// aflowlib::XPLUG
// ***************************************************************************
namespace aflowlib {
  bool XPLUG_Directory_ok(const string &dirIN,string &directory_OUT,const uint &i,const uint &j,const int& ithread,stringstream &oss);
  void *_threaded_interface_XPLUG_Directory(void *ptr);

  typedef struct {
    int      itbusy;          // FOR XPLUG (ALL)
    bool     VERBOSE;         // FOR XPLUG (ALL)
    //  ostringstream  oss;   // printout
    deque<bool>   *vok;       // directory to check
    deque<string> *vdirsIN;   // FOR
    deque<string> *vdirsOUT;  // FOR
    deque<uint>   *vrun;      // index to run
    int      ITHREAD;         // FOR
    int      THREADS_MAX;     // FOR
  } _threaded_XPLUG_params;

  //_threaded_XPLUG_params params[MAX_ALLOCATABLE_PTHREADS];
  pthread_mutex_t mutex_LIBRARIES=PTHREAD_MUTEX_INITIALIZER;

  //#define _PTHREAD_FLUSH_TIME_ 1
#define _XPLUG_FLUSH_THREAD_SLEEP_ 1

  void *_threaded_interface_XPLUG_Directory(void *ptr) {
    bool CONCURRENT=TRUE;
    _threaded_XPLUG_params* pparams=(_threaded_XPLUG_params*) ptr;
    string directory; uint ith;
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=TRUE;
    AFLOW_PTHREADS::RUNNING++;
    if((pparams->VERBOSE)) { pthread_mutex_lock(&mutex_LIBRARIES);cout << "_threaded_COMMANDS " << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << endl;pthread_mutex_unlock(&mutex_LIBRARIES); }
    if(CONCURRENT==FALSE) { // SPLITS tasks in MOD(THREAD_MAX);
      for(ith=(pparams->ITHREAD);ith<(*pparams->vdirsIN).size();ith+=(pparams->THREADS_MAX)) {
        if((pparams->VERBOSE)) { pthread_mutex_lock(&mutex_LIBRARIES);cout << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << " " << ith << " " << (*pparams->vdirsIN).at(ith) << endl;pthread_mutex_unlock(&mutex_LIBRARIES); }
        stringstream oss;
        (*pparams->vok).at(ith)=XPLUG_Directory_ok((*pparams->vdirsIN).at(ith),(*pparams->vdirsOUT).at(ith),ith,(*pparams->vdirsIN).size(),ith,oss); // do this
        pthread_mutex_lock(&mutex_LIBRARIES);
        cerr << oss.str();cerr.flush();
        pthread_mutex_unlock(&mutex_LIBRARIES);
      }
    } else { // RUNS in a queue
      while((*pparams->vrun).size()>0) {
        bool FRONT=TRUE;
        pthread_mutex_lock(&mutex_LIBRARIES);
        if(FRONT)  { ith=(*pparams->vrun).at(0);(*pparams->vrun).pop_front(); } // from the front
        if(!FRONT) { ith=(*pparams->vrun).at((*pparams->vrun).size()-1);(*pparams->vrun).pop_back(); } // from the back
        pthread_mutex_unlock(&mutex_LIBRARIES);
        if((pparams->VERBOSE)) { pthread_mutex_lock(&mutex_LIBRARIES);cout << (pparams->ITHREAD) << "/" << (pparams->THREADS_MAX) << ": " << ith << endl;pthread_mutex_unlock(&mutex_LIBRARIES); }
        stringstream oss;
        (*pparams->vok).at(ith)=XPLUG_Directory_ok((*pparams->vdirsIN).at(ith),(*pparams->vdirsOUT).at(ith),ith,(*pparams->vdirsIN).size(),ith,oss); // do this
        pthread_mutex_lock(&mutex_LIBRARIES);
        cerr << oss.str();cerr.flush();
        pthread_mutex_unlock(&mutex_LIBRARIES);
      }
    }
    AFLOW_PTHREADS::vpthread_busy[pparams->itbusy]=FALSE;
    AFLOW_PTHREADS::RUNNING--;
    return NULL;
  }
}

// DO ONE
namespace aflowlib {
  bool XPLUG_Directory_ok(const string &dirIN,string &dir,const uint &i,const uint &j,const int& ithread,stringstream &oss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"aflowlib::XPLUG_Directory_ok():"; //CO20200404
    dir=dirIN; // start from here
    bool ok=FALSE,print=TRUE;
    int answer=0;
    stringstream obb;

    // clean
    if(LDEBUG) cerr << soliloquy << " [1]: analyzing = " << dir << endl;

    dir=aurostd::RemoveSubString(dir,"/"+_AFLOWIN_);
    dir=aurostd::RemoveSubString(dir,"/aflow.end.out");

    deque<string> vrelax; aurostd::string2tokens(".relax1,.relax2,.relax3,.static,.bands",vrelax,",");
    deque<string> vremdirname;aurostd::string2tokens("OUTCAR,OSZICAR,CONTCAR,EIGENVAL",vremdirname,",");
    deque<string> vbroken;  // no more list... this checks everything

    for(uint iremdirname=0;iremdirname<vremdirname.size();iremdirname++) {
      for(uint irelax=0;irelax<vrelax.size();irelax++) {
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          dir=aurostd::RemoveSubString(dir,"/"+vremdirname.at(iremdirname)+vrelax.at(irelax)+XHOST.vext.at(iext));
        }
        dir=aurostd::RemoveSubString(dir,"/"+vremdirname.at(iremdirname)+vrelax.at(irelax));
      }
      dir=aurostd::RemoveSubString(dir,"/"+vremdirname.at(iremdirname));
    }

    if(LDEBUG) cerr << soliloquy << " [2]: analyzing = " << dir << endl;

    if(!aurostd::FileExist(dir+"/"+_AFLOWLOCK_)) {
      if(LDEBUG) cerr << soliloquy << " [3]: not finding " << _AFLOWLOCK_ << " = " << dir << endl;
      return FALSE;
    }
    if(!aurostd::EFileExist(dir+"/OUTCAR.relax1") && 
        !aurostd::EFileExist(dir+"/OUTCAR.relax2") &&
        !aurostd::EFileExist(dir+"/OUTCAR.relax3") &&
        !aurostd::EFileExist(dir+"/OUTCAR.static") &&
        !aurostd::EFileExist(dir+"/OUTCAR.bands")) {
      if(LDEBUG) cerr << soliloquy << " [4]: not finding any OUTCAR.* = " << dir << endl;
      return FALSE;
    }

    obb << "[" << i+1 << "/" << j << "] ";

    if(ithread<0)  obb << "DIR=" << dir << " ";
    if(ithread>=0) obb << "DIR(" << ithread << ")=" << dir << " ";
    // CHECK ALL FILES
    print=TRUE;
    ok=TRUE;
    // TEST INCOMPLETE
    if(ok) { obb << ".";if(!aurostd::FileExist(dir+"/"+_AFLOWLOCK_)) { ok=FALSE;obb << " no=LOCK"; }}
    if(ok) { obb << ".";if(aurostd::FileExist(dir+"/EIGENVAL")) { ok=FALSE;obb << " yes=EIGENVAL"; }}
    if(ok) { obb << ".";if(aurostd::FileExist(dir+"/vasp.out")) { ok=FALSE;obb << " yes=vasp.out"; }}
    if(ok) { obb << ".";if(aurostd::FileExist(dir+"/AECCAR0")) { ok=FALSE;obb << " yes=AECCAR0"; }}
    if(ok) { obb << ".";if(aurostd::FileExist(dir+"/AECCAR1")) { ok=FALSE;obb << " yes=AECCAR1"; }}
    if(ok) { obb << ".";if(aurostd::FileExist(dir+"/AECCAR2")) { ok=FALSE;obb << " yes=AECCAR2"; }}
    if(ok) { obb << ".";if(aurostd::FileExist(dir+"/CHGCAR")) { ok=FALSE;obb << " yes=CHGCAR"; }}
    if(ok) { obb << ".";if(aurostd::FileExist(dir+"/DOSCAR")) { ok=FALSE;obb << " yes=DOSCAR"; }}
    if(ok) { obb << ".";if(aurostd::FileExist(dir+"/OUTCAR")) { ok=FALSE;obb << " yes=OUTCAR"; }}
    if(ok) { obb << ".";if(aurostd::FileExist(dir+"/OSZICAR")) { ok=FALSE;obb << " yes=OSZICAR"; }}
    if(ok) { obb << ".";if(aurostd::FileExist(dir+"/POTCAR")) { ok=FALSE;obb << " yes=POTCAR"; }}
    for(uint irelax=0;irelax<vrelax.size();irelax++) {
      if(ok) { obb << ".";if(aurostd::FileExist(dir+"/OUTCAR"+vrelax.at(irelax))) { ok=FALSE;obb << " yes=OUTCAR"+vrelax.at(irelax); }}
      if(ok) { obb << ".";if(aurostd::FileExist(dir+"/OSZICAR"+vrelax.at(irelax))) { ok=FALSE;obb << " yes=OSZICAR"+vrelax.at(irelax); }}
      if(ok) { obb << ".";if(aurostd::FileExist(dir+"/POTCAR"+vrelax.at(irelax))) { ok=FALSE;obb << " yes=POTCAR"+vrelax.at(irelax); }}
    }
    // TEST RELAX1 RELAX2 RELAX3
    for(uint irelax=0;irelax<vrelax.size();irelax++) {
      if(ok && aurostd::EFileExist(dir+"/OUTCAR"+vrelax.at(irelax))) { // relax1 relax2 relax3
        if(ok) { obb << ".";if(!aurostd::EFileExist(dir+"/OUTCAR"+vrelax.at(irelax))) { ok=FALSE;print=FALSE;obb << " no=OUTCAR"+vrelax.at(irelax)+".EXT"; }}
        if(ok) { obb << ".";if(!aurostd::EFileExist(dir+"/OSZICAR"+vrelax.at(irelax))) { ok=FALSE;obb << " no=OSZICAR"+vrelax.at(irelax)+".EXT"; }}
        if(ok) { obb << ".";if(!aurostd::EFileExist(dir+"/vasp.out"+vrelax.at(irelax))) { ok=FALSE;obb << " no=vasp.out"+vrelax.at(irelax)+".EXT"; }}
        if(ok) { obb << ".";if(!aurostd::EFileExist(dir+"/INCAR"+vrelax.at(irelax))) { ok=FALSE;obb << " no=INCAR"+vrelax.at(irelax)+".EXT"; }}
        //      if(ok) { obb << ".";if(!aurostd::EFileExist(dir+"/EIGENVAL"+vrelax.at(irelax))) { ok=FALSE;obb << " no=EIGENVAL"+vrelax.at(irelax)+".EXT"; }}
        if(vrelax.at(irelax)==".static")
          if(ok) { obb << "d"; if(!aurostd::EFileExist(dir+"/DOSCAR"+vrelax.at(irelax))) { ok=FALSE;obb << " no=DOSCAR"+vrelax.at(irelax)+".EXT"; }}
        // "OUTCAR"+vrelax.at(irelax)+".EXT"  //CO20200106 - patching for auto-indenting (quotes)
        obb << ".";
        if(ok) { obb << "o";  // check Answer 4 or 5 in OUTCAR.RELAX.EXT
          for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
            if(aurostd::FileExist(dir+"/OUTCAR"+vrelax.at(irelax)+XHOST.vext.at(iext))) {
              answer=aurostd::execute2utype<int>(XHOST.vcat.at(iext)+" \""+dir+"/OUTCAR"+vrelax.at(irelax)+XHOST.vext.at(iext)+"\" | grep -c \"(sec)\" ");
              if(answer==0){answer=aurostd::execute2utype<int>(XHOST.vcat.at(iext)+" \""+dir+"/OUTCAR"+vrelax.at(irelax)+XHOST.vext.at(iext)+"\" | grep -c \"(sec)\" ");} //CO20200501 - soft patch, this could STILL break if xcat does NOT load anything
              if(answer!=4 && answer!=5) { ok=FALSE;obb << ". error(" << answer << ")=OUTCAR"+vrelax.at(irelax)+".EXT"; }}
          }
        }
        obb << ".";
        if(ok) { obb << "v"; // TEST VASP.OUT
          for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
            if(aurostd::FileExist(dir+"/vasp.out"+vrelax.at(irelax)+XHOST.vext.at(iext))) {
              answer=aurostd::execute2utype<int>(XHOST.vcat.at(iext)+" \""+dir+"/vasp.out"+vrelax.at(irelax)+XHOST.vext.at(iext)+"\" | grep -c \"The distance between some ions is very small\" ");
              if(answer!=0) { ok=FALSE;obb << ". ions=vasp.out"+vrelax.at(irelax)+".EXT"; }}
          }
        }
        obb << ".";
      }
      if(ok) {// obb << "b"; // TEST "+vbroken.at(ibroken)+"
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          if(aurostd::FileExist(dir+"/OUTCAR"+vrelax.at(irelax)+XHOST.vext.at(iext))) {
            //	    cerr << soliloquy << " " << string("ls "+dir+"/*"+vrelax.at(irelax)+XHOST.vext.at(iext)) << endl;
            aurostd::string2dequestring(aurostd::execute2string("ls "+dir+"/*"+vrelax.at(irelax)+XHOST.vext.at(iext)),vbroken);
            for(uint ibroken=0;ibroken<vbroken.size();ibroken++) {
              //      cerr << soliloquy << " " << "[" << vbroken.at(ibroken) << "]" << endl;
              if(ok && aurostd::FileExist(vbroken.at(ibroken))) {
                obb << "b";
                //		answer=aurostd::execute2utype<int>(XHOST.vcat.at(iext)+" \""+dir+"/"+vbroken.at(ibroken)+"\" 2>&1 | grep -c \"Unexpected end of input\" ");
                answer=aurostd::execute2utype<int>(XHOST.vcat.at(iext)+" \""+vbroken.at(ibroken)+"\" 2>&1 | grep -c \"Unexpected end of input\" ");
                if(answer!=0) { ok=FALSE;obb << ". Broken_file="+vbroken.at(ibroken)+".EXT"; }}
            }
          }
        } 
      }
    }
    // DONE
    if(ok==TRUE) obb << " good";
    if(ok==FALSE) obb << " bad";
    if(aurostd::FileExist(dir+"/OUTCAR.relax2.gz")) obb << "      gz";
    if(aurostd::FileExist(dir+"/OUTCAR.relax2.bz2")) obb << "      bz2";
    if(aurostd::FileExist(dir+"/OUTCAR.relax2.xz")) obb << "      xz";
    obb << endl;
    if(print==TRUE) oss << obb.str();
    oss.flush();

    //DX+CO20191030 - check for identical subdirectories - START
    // check if directory has an identical subdirectory
    vector<string> dir_tokens;
    aurostd::string2tokens(dir,dir_tokens,"/");
    string run_directory = dir_tokens[dir_tokens.size()-1];
    string sub_directory = dir+"/"+run_directory;
    if(aurostd::FileExist(sub_directory)) {
      obb << "This directory contains an identical subdirectory (subdir=" << sub_directory << "). Removing subdirectory." << endl;
      oss << obb.str();
      oss.flush();
      aurostd::RemoveDirectory(sub_directory);
    }
    //DX+CO20191030 - check for identical subdirectories - END

    return ok;
  }
}

namespace aflowlib {
  bool XPLUG_CHECK_ONLY(const vector<string>& argv) {  //CO20200501
    deque<string> vdirsOUT,vzips,vcleans;
    return XPLUG_CHECK_ONLY(argv,vdirsOUT,vzips,vcleans);
  }
  bool XPLUG_CHECK_ONLY(const vector<string>& argv,deque<string>& vdirsOUT,deque<string>& vzips,deque<string>& vcleans) {  //CO20200501
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy=XPID+"aflowlib::XPLUG_CHECK_ONLY():";  //CO20200404
    //   cerr << soliloquy << " is deprecated" << endl; exit(0);
    int NUM_THREADS=XHOST.CPU_Cores; //ME20181226
    //ME20181109 - Handle NCPUS=MAX
    if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS") && !(XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX")))
      NUM_THREADS=aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS"));

    if(NUM_THREADS<1) NUM_THREADS=1; //[CO20200502 - --np>1 is NOT guaranteed to work, prefer to run --np=1]NUM_THREADS=XHOST.CPU_Cores/2;

    if(LDEBUG) cerr << soliloquy << " NUM_THREADS=" << NUM_THREADS << endl;
    //  exit(0);
    // if(LDEBUG) cerr << soliloquy << " = " << XHOST.hostname << endl;
    deque<string> vdirsIN; //CO20200501,vdirsOUT,vzips,vcleans;
    vdirsOUT.clear();vzips.clear();vcleans.clear(); //CO20200501
    deque<bool> vok;
    deque<uint> vrun;

    for(uint i=1;i<argv.size();i++) {
      if(aurostd::substring2bool(argv.at(i),"/")) {
        //	cerr << soliloquy << " adding" << argv.at(i) << endl;
        vdirsIN.push_back(argv.at(i));
        vdirsOUT.push_back(argv.at(i));
      }
    }

    std::sort(vdirsIN.begin(),vdirsIN.end());
    std::sort(vdirsOUT.begin(),vdirsOUT.end());
    if(LDEBUG) cerr << soliloquy << " vdirsIN.size()=" << vdirsIN.size() << endl;
    //  if(LDEBUG) cerr << soliloquy << " [1]" << endl;
    for(uint i=0;i<vdirsIN.size();i++) { vok.push_back(TRUE);vrun.push_back(i); }
    //  if(LDEBUG) cerr << soliloquy << " [2]" << endl;

    if((int) vdirsIN.size()<=NUM_THREADS) NUM_THREADS=(uint) vdirsIN.size();        // SAFETY
    //  if(LDEBUG) cerr << soliloquy << " [3]" << endl;

    if(NUM_THREADS<=1) {                                                        // run singular
      //  if(LDEBUG) cerr << soliloquy << " [3b]" << endl;
      for(uint i=0;i<vdirsIN.size();i++) {
        stringstream oss;
        vok.at(i)=XPLUG_Directory_ok(vdirsIN.at(i),vdirsOUT.at(i),i,vdirsIN.size(),-1,oss);
        cerr << oss.str();cerr.flush();
      }
    }
    //  if(LDEBUG) cerr << soliloquy << " [4]" << endl;
    if(NUM_THREADS>=2) {                                                       // multithread
      AFLOW_PTHREADS::FLAG=TRUE;AFLOW_PTHREADS::MAX_PTHREADS=NUM_THREADS;      // prepare
      if(AFLOW_PTHREADS::MAX_PTHREADS>MAX_ALLOCATABLE_PTHREADS) AFLOW_PTHREADS::MAX_PTHREADS=MAX_ALLOCATABLE_PTHREADS; // check max
      //  if(LDEBUG) cerr << soliloquy << " [5]" << endl;
      AFLOW_PTHREADS::Clean_Threads();                                         // clean threads
      //    if(LDEBUG) cerr << soliloquy << " [6]" << endl;
      _threaded_XPLUG_params params[MAX_ALLOCATABLE_PTHREADS];                 // prepare
      //    if(LDEBUG) cerr << soliloquy << " [7]" << endl;
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) {      // prepare loop
        params[ithread].ITHREAD=ithread;                                       // prepare params
        params[ithread].THREADS_MAX=AFLOW_PTHREADS::MAX_PTHREADS;              // prepare params
        params[ithread].vok=&vok;                                              // prepare params
        params[ithread].vdirsIN=&vdirsIN;                                      // prepare params
        params[ithread].vdirsOUT=&vdirsOUT;                                    // prepare params
        params[ithread].vrun=&vrun;                                            // prepare params
        params[ithread].itbusy=ithread;                                        // prepare params
        params[ithread].VERBOSE=FALSE;                                         // prepare params
      }
      //    if(LDEBUG) cerr << soliloquy << " [8]" << endl;
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) AFLOW_PTHREADS::viret[ithread]=pthread_create(&AFLOW_PTHREADS::vpthread[ithread],NULL,_threaded_interface_XPLUG_Directory,(void*)&params[ithread]); // run threads
      //    if(LDEBUG) cerr << soliloquy << " [9]" << endl;
      for(int ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++) pthread_join(AFLOW_PTHREADS::vpthread[ithread],NULL); // flush threads
      //  if(LDEBUG) cerr << soliloquy << " [10]" << endl;
    }

    // has OK. now do the counting
    for(uint i=0;i<vdirsIN.size();i++) {
      if(vok.at(i)==TRUE) vzips.push_back(vdirsOUT.at(i));
      if(vok.at(i)==FALSE) vcleans.push_back(vdirsOUT.at(i));
    }
    if(LDEBUG) cerr << soliloquy << " vdirsIN.size()=" << vdirsIN.size() << endl;
    if(LDEBUG) cerr << soliloquy << " vzips.size()=" << vzips.size() << endl;
    if(LDEBUG) cerr << soliloquy << " vcleans.size()=" << vcleans.size() << endl;

    return true;
  }
  bool XPLUG(const vector<string>& argv) {
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy=XPID+"aflowlib::XPLUG():";  //CO20200404
    //   cerr << soliloquy << " is deprecated" << endl; exit(0);
    int NUM_ZIP=aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_ZIP"));
    int NUM_SIZE=aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_SIZE"));
    bool FLAG_DO_CLEAN=XHOST.vflag_control.flag("XPLUG_DO_CLEAN");
    bool FLAG_DO_ADD=XHOST.vflag_control.flag("XPLUG_DO_ADD");
    string PREFIX=XHOST.vflag_control.getattachedscheme("XPLUG_PREFIX");

    if(NUM_ZIP<1) NUM_ZIP=1;
    if(NUM_SIZE<1) NUM_SIZE=128;

    if(LDEBUG) cerr << soliloquy << " NUM_ZIP=" << NUM_ZIP << endl;
    if(LDEBUG) cerr << soliloquy << " NUM_SIZE=" << NUM_SIZE << endl;
    if(LDEBUG) cerr << soliloquy << " FLAG_DO_CLEAN=" << FLAG_DO_CLEAN << endl;
    if(LDEBUG) cerr << soliloquy << " FLAG_DO_ADD=" << FLAG_DO_ADD << endl;
    if(LDEBUG) cerr << soliloquy << " PREFIX=" << PREFIX << endl;

    deque<string> vdirsOUT,vzips,vcleans;
    XPLUG_CHECK_ONLY(argv,vdirsOUT,vzips,vcleans); //CO20200501

    if(vzips.size()>0) {
      stringstream command;
      if(aurostd::substring2bool(XHOST.hostname,"m7int0") || aurostd::substring2bool(XHOST.hostname,"m6int0")) XHOST.hostname="marylou";
      XHOST.hostname=aurostd::RemoveSubString(XHOST.hostname,".egr.duke.edu");
      XHOST.hostname=aurostd::RemoveSubString(XHOST.hostname,".mems.duke.edu");
      XHOST.hostname=aurostd::RemoveSubString(XHOST.hostname,".pratt.duke.edu");
      XHOST.hostname=aurostd::RemoveSubString(XHOST.hostname,".duke.edu");
      for(uint i=0;i<vzips.size();i+=AFLOW_MAX_ARGV) {
        command << "aflow --multi=zip " << (FLAG_DO_ADD?"--add ":"") << "--np=" << NUM_ZIP;
        command << " --prefix=update_" << (PREFIX!=""?string(PREFIX+"_"):string(""));
        //[CO20200501 - OBSOLETE]command << aurostd::get_date() <<  "-" << aurostd::get_hour() << aurostd::get_min() << aurostd::get_sec();
        command << aurostd::get_datetime_formatted("",true,"-",""); //CO20200501 - no delim for date and time, include time, "-" between date and time
        command << "_" << i/AFLOW_MAX_ARGV+1 << "_" << XHOST.hostname;
        command << " --size=" << NUM_SIZE << " --DIRECTORY ";
        for(uint j=i;j<i+AFLOW_MAX_ARGV && j<vzips.size();j++) command << " " << vzips.at(j);
        command << " " << endl << endl;
      }
      cerr << command.str() << endl;
      aurostd::execute(command);
    }

    //if(vzips.size()>0) {
    //stringstream command;
    //if(aurostd::substring2bool(XHOST.hostname,"m7int0") || aurostd::substring2bool(XHOST.hostname,"m6int0")) XHOST.hostname="marylou";
    //XHOST.hostname=aurostd::RemoveSubString(XHOST.hostname,".egr.duke.edu");
    //XHOST.hostname=aurostd::RemoveSubString(XHOST.hostname,".mems.duke.edu");
    //XHOST.hostname=aurostd::RemoveSubString(XHOST.hostname,".pratt.duke.edu");
    //XHOST.hostname=aurostd::RemoveSubString(XHOST.hostname,".duke.edu");
    //command << "aflow --multi=zip " << (FLAG_DO_ADD?"--add ":"") << "--np=" << NUM_ZIP
    // << " --prefix=update_" << (PREFIX!=""?string(PREFIX+"_"):string("")) << aurostd::get_date()
    // << "-" << aurostd::get_hour() << aurostd::get_min() << aurostd::get_sec() << "_" << XHOST.hostname
    // << " --size=" << NUM_SIZE << " --DIRECTORY ";
    //for(uint j=0;j<vzips.size()&&j<AFLOW_MAX_ARGV;j++)
    //command << " " << vzips.at(j);
    //cerr << command.str() << endl;
    //aurostd::execute(command);
    //}

    if(FLAG_DO_CLEAN && vcleans.size()>0) {
      for(uint i=0;i<vcleans.size();i++) {
        // cerr << "Cleaning=" << vcleans.at(i) << endl;
        KBIN::Clean(vcleans.at(i));
      }
    }
    return FALSE;
  }
}

// ***************************************************************************
// void aflowlib::vaspfile2stringstream(string& str_dir, const string& FILE, stringstream& ss_vaspfile)
// ***************************************************************************
namespace aflowlib {
  bool VaspFileExist(const string& str_dir, const string& FILE) {
    bool RUN_FLAG=FALSE;
    if(!RUN_FLAG && aurostd::FileExist(str_dir+"/"+FILE)) RUN_FLAG=TRUE;
    if(!RUN_FLAG && aurostd::EFileExist(str_dir+"/"+FILE)) RUN_FLAG=TRUE;
    if(!RUN_FLAG && aurostd::FileExist(str_dir+"/"+FILE+".static")) RUN_FLAG=TRUE;
    if(!RUN_FLAG && aurostd::EFileExist(str_dir+"/"+FILE+".static")) RUN_FLAG=TRUE;
    if(!RUN_FLAG && aurostd::FileExist(str_dir+"/"+FILE+".bands")) RUN_FLAG=TRUE;
    if(!RUN_FLAG && aurostd::EFileExist(str_dir+"/"+FILE+".bands")) RUN_FLAG=TRUE;
    if(!RUN_FLAG) {
      cerr << FILE+" or "+FILE+".static/bands or "+FILE+".static/bands.EXT not found in the directory!" << endl;
    }
    return RUN_FLAG;
  }
}

// ***************************************************************************
// void aflowlib::vaspfile2stringstream(string& str_dir, const string& FILE, stringstream& ss_vaspfile)
// ***************************************************************************
namespace aflowlib {
  uint _OLD_vaspfile2stringstream(const string& str_dir, const string& FILE, stringstream& sss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    if(LDEBUG) cerr << XPID << "vaspfile2stringstream: BEGIN" << endl;
    //if XXXX.bands exists, you may also use this function by setting FILE=vaspfile.bands
    sss.clear(); sss.str(std::string(""));
    bool gfound=FALSE;
    if(!gfound && (FILE=="EIGENVAL" || FILE=="KPOINTS")) {
      gfound=TRUE;
      if(LDEBUG) cerr << XPID << "vaspfile2stringstream: 1=" << FILE << endl;
      for(uint iext=0;iext<XHOST.vext.size();iext++) {
        if(aurostd::FileExist(str_dir+"/"+FILE+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
        if(aurostd::FileExist(str_dir+"/"+FILE+".bands"+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+".bands"+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
      }
      cerr << FILE+" or "+FILE+".bands or "+FILE+".bands.EXT not found in the directory, aborting!" << endl;
      exit(1);
    }
    if(!gfound && (FILE=="DOSCAR")) {
      gfound=TRUE;
      if(LDEBUG) cerr << XPID << "vaspfile2stringstream: 2=" << FILE << endl;
      for(uint iext=0;iext<XHOST.vext.size();iext++) {
        if(aurostd::FileExist(str_dir+"/"+FILE+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
        if(aurostd::FileExist(str_dir+"/"+FILE+".static"+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+".static"+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
      }
      cerr << FILE+" or "+FILE+".static or "+FILE+".static.EXT not found in the directory, aborting!" << endl;
      exit(1);
    }
    if(!gfound && (FILE=="POSCAR")) {
      gfound=TRUE;
      if(LDEBUG) cerr << XPID << "vaspfile2stringstream: 1=" << FILE << endl;
      for(uint iext=0;iext<XHOST.vext.size();iext++) {
        if(aurostd::FileExist(str_dir+"/"+FILE+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
        if(aurostd::FileExist(str_dir+"/"+FILE+".bands"+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+".bands"+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
        if(aurostd::FileExist(str_dir+"/"+FILE+".static"+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+".static"+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
        if(aurostd::FileExist(str_dir+"/"+FILE+".relax1"+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+".relax1"+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
      }
      cerr << FILE+" or "+FILE+".bands/static/relax or "+FILE+".bands./static/relax.EXT not found in the directory, aborting!" << endl;
      exit(1);
    }
    if(!gfound) {
      if(LDEBUG) cerr << XPID << "vaspfile2stringstream: 3=" << FILE << endl;
      for(uint iext=0;iext<XHOST.vext.size();iext++) {
        if(aurostd::FileExist(str_dir+"/"+FILE+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
        if(aurostd::FileExist(str_dir+"/"+FILE+".static"+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+".static"+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
        if(aurostd::FileExist(str_dir+"/"+FILE+".relax1"+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+".relax1"+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
        if(aurostd::FileExist(str_dir+"/"+FILE+".bands"+XHOST.vext.at(iext))) { sss << aurostd::execute2string(XHOST.vcat.at(iext)+" \""+ str_dir+"/"+FILE+".bands"+XHOST.vext.at(iext)+"\"");return sss.str().length(); }
      }
      cerr << FILE+" or "+FILE+".static/relax1/bands or "+FILE+".static/relax1/bands.EXT not found in the directory, aborting!" << endl;
      exit(1);
    }
    if(LDEBUG) cerr << XPID << "vaspfile2stringstream: END" << endl;
    return sss.str().length();
  }
}


// ***************************************************************************
// void aflowlib::vaspfile2stringstream
// ***************************************************************************
namespace aflowlib {
  string vaspfile2stringstream(const string& str_dir, const string& FILE, stringstream& sss) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: BEGIN" << endl;
    //if XXXX.bands exists, you may also use this function by setting FILE=vaspfile.bands
    sss.clear(); sss.str(std::string(""));
    bool gfound=FALSE;
    if(!gfound && (FILE=="EIGENVAL" || FILE=="KPOINTS")) {
      gfound=TRUE;
      if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: FILE=" << FILE << endl;
      string extension=".bz2,.gz,.xz,.bands,.bands.bz2,.bands.gz,.bands.xz";
      deque<string> vtype;aurostd::string2tokens(extension,vtype,",");vtype.push_front(""); // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las \""+str_dir+"\"");
      for(uint i=0;i<vtype.size();i++) {
        if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: TESTING FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
        if(aurostd::FileExist(str_dir+"/"+FILE+vtype.at(i))) {
          aurostd::efile2stringstream(str_dir+"/"+FILE+vtype.at(i),sss);
          if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: FOUND FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
          return aurostd::CleanFileName(str_dir+"/"+FILE+vtype.at(i));
        }
      }
      //ME20190627 BEGIN
      // [OBSOLETE] cerr << FILE+" or "+FILE+".bands or "+FILE+".bands.EXT not found in the directory, aborting!" << endl;
      // [OBSOLETE] exit(1);
      string message = FILE+" or "+FILE+".bands or "+FILE+".bands.EXT not found in the directory, aborting!";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,"aflowlib::vaspfile2stringstream()", message, _FILE_NOT_FOUND_);
      //ME20190627 END
    }
    if(!gfound && (FILE=="DOSCAR")) {
      gfound=TRUE;
      if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: FILE=" << FILE << endl;
      string extension=".bz2,.gz,.xz,.static,.static.bz2,.static.gz,.static.xz";
      deque<string> vtype;aurostd::string2tokens(extension,vtype,",");vtype.push_front("");  // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las \""+str_dir+"\"");
      for(uint i=0;i<vtype.size();i++) {
        if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: TESTING FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
        if(aurostd::FileExist(str_dir+"/"+FILE+vtype.at(i))) {
          aurostd::efile2stringstream(str_dir+"/"+FILE+vtype.at(i),sss);
          if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: FOUND FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
          return aurostd::CleanFileName(str_dir+"/"+FILE+vtype.at(i));
        }
      }
      //ME20190627 BEGIN
      // [OBSOLETE] cerr << FILE+" or "+FILE+".static or "+FILE+".static.EXT not found in the directory, aborting!" << endl;
      // [OBSOLETE] exit(1);
      string message = FILE+" or "+FILE+".static or "+FILE+".static.EXT not found in the directory, aborting!";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,"aflowlib::vaspfile2stringstream()", message, _FILE_NOT_FOUND_);
      //ME20190627 END
    }
    if(!gfound && (FILE=="POSCAR" || FILE=="CONTCAR")) { //RF20200409
      gfound=TRUE;
      if(LDEBUG) cerr << "aflowlib::vaspfile2stringstream: FILE=" << FILE << endl;
      string extension=".bz2,.gz,.xz,.bands,.bands.bz2,.bands.gz,.bands.xz,.static,.static.bz2,.static.gz,.static.xz,.relax2,.relax2.bz2,.relax2.gz,.relax2.xz,.relax1,.relax1.bz2,.relax1.gz,.relax1.xz"; //RF20200409
      deque<string> vtype;aurostd::string2tokens(extension,vtype,",");vtype.push_front(""); // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las \""+str_dir+"\"");
      for(uint i=0;i<vtype.size();i++) {
        if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: TESTING FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
        if(aurostd::FileExist(str_dir+"/"+FILE+vtype.at(i))) {
          aurostd::efile2stringstream(str_dir+"/"+FILE+vtype.at(i),sss);
          if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: FOUND FILE=" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
          return aurostd::CleanFileName(str_dir+"/"+FILE+vtype.at(i));
        }
      }
      //ME20190627 BEGIN
      // [OBSOLETE] cerr << FILE+" or "+FILE+".bands/static/relax or "+FILE+".bands./static/relax.EXT not found in the directory, aborting!" << endl;
      // [OBSOLETE] exit(1);
      string message = FILE+" or "+FILE+".bands/static/relax or "+FILE+".bands./static/relax.EXT not found in the directory, aborting!";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,"aflowlib::vaspfile2stringstream()", message, _FILE_NOT_FOUND_);
      //ME20190627 END
    }
    if(!gfound) {
      if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: FILE=" << FILE << endl;
      string extension=".bz2,.gz,.xz,.static,.static.bz2,.static.gz,.static.xz,.relax1,.relax1.bz2,.relax1.gz,.relax1.xz,.bands,.bands.bz2,.bands.gz,.bands.xz";
      deque<string> vtype;aurostd::string2tokens(extension,vtype,",");vtype.push_front(""); // have to add emptyness to vtype at the beginning
      //     if(LDEBUG) aurostd::execute("ls -las \""+str_dir+"\"");
      for(uint i=0;i<vtype.size();i++) {
        if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: TESTING FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
        if(aurostd::FileExist(str_dir+"/"+FILE+vtype.at(i))) {
          aurostd::efile2stringstream(str_dir+"/"+FILE+vtype.at(i),sss);
          if(LDEBUG) cerr << XPID << "aflowlib::vaspfile2stringstream: FOUND FILE=[" << str_dir+"/"+FILE+vtype.at(i) << "]" << endl;
          return aurostd::CleanFileName(str_dir+"/"+FILE+vtype.at(i));
        }
      }
      //ME20190627 BEGIN
      // [OBSOLETE] cerr << FILE+" or "+FILE+".static/relax1/bands or "+FILE+".static/relax1/bands.EXT not found in the directory, aborting!" << endl;
      // [OBSOLETE] exit(1);
      string message = FILE+" or "+FILE+".static/relax1/bands or "+FILE+".static/relax1/bands.EXT not found in the directory, aborting!";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,"aflowlib::vaspfile2stringstream()", message, _FILE_NOT_FOUND_);
      //ME20190627 END
    }
    if(LDEBUG) cerr << XPID << "vaspfile2stringstream: END" << endl;
    return string("");
  }
}

namespace aflowlib {
  string vaspfile2stringstream(const string& str_dir, const string& FILE) {
    stringstream sss;
    return vaspfile2stringstream(str_dir,FILE,sss);
  }
}

// ***************************************************************************
// aflowlib::LIB2LIB
// ***************************************************************************
// Added by CT20181212
namespace aflowlib {
  bool LIB2LIB(const string& options,bool flag_FORCE,bool LOCAL) {
    bool LDEBUG=(FALSE || XHOST.DEBUG); //CO20200624
    string soliloquy=XPID+"aflowlib::LIB2LIB():";
    if(LDEBUG) cerr << soliloquy << " BEGIN" << endl;
    if(LDEBUG) cerr << soliloquy << " options=" << options << endl;
    vector<string> tokens;

    if(!aurostd::substring2bool(options,"POCC")) {
      aurostd::string2tokens(options,tokens,","); // QUICK FIX
    } else {
      tokens.clear();
      tokens.push_back(options);
    }

    if(LDEBUG) cerr << soliloquy << " tokens.size()=" << tokens.size() << endl;
    if(tokens.size()==0 || tokens.size()>2) {
      init::ErrorOption(options,soliloquy,aurostd::liststring2string("aflow --lib2lib=directory","aflow --lib2lib=all[,dir]"));
      exit(0);
    }

    if(tokens.size()>=1) {
      if(tokens.at(0)=="all") {
        XHOST.sensors_allowed=FALSE;
        XHOST.vflag_pflow.flag("MULTI=SH",FALSE);
        return LIB2RAW_ALL(options,flag_FORCE);
        XHOST.sensors_allowed=TRUE;
      }
    }

    string directory_LIB="";
    // OBSOLETE bool flag_files_LIB=FALSE;
    string PROJECT_LIBRARY="";
    if(LOCAL) {
      flag_FORCE=true;
      string directory=options;
      CleanDirectoryLIB(directory);
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"./","");
      //[CO20200624 - do inside CleanDirectoryLIB()]if(directory=="." || directory.empty()) { directory=aurostd::getPWD(); } //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")
      // OBSOLETE flag_files_LIB=FALSE;
      directory_LIB=directory;
      PROJECT_LIBRARY=directory_LIB;
    } else {    //normal run
      //  cout << soliloquy << " AFLOW (" << VERSION << ")" << endl;
      CheckMaterialServer(soliloquy); // must be in AFLOW_MATERIALS_SERVER
      string directory=options;
      CleanDirectoryLIB(directory);
      //[CO20200624 - do inside CleanDirectoryLIB()]string directory=aurostd::CleanFileName(options);
      //[CO20200624 - do inside CleanDirectoryLIB()]if(directory.at(directory.size()-1)=='/')  directory=directory.substr(0,directory.size()-1);
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"/"+_AFLOWIN_,"");
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"/"+_AFLOWLOCK_,"");
      //[CO20200624 - do inside CleanDirectoryLIB()]if(directory.size()>=5 && directory[directory.size()-5]=='/' && directory[directory.size()-4]=='c' && directory[directory.size()-3]=='o' && directory[directory.size()-2]=='r' && directory[directory.size()-1]=='e'){directory=directory.substr(0,directory.size()-5);} //aurostd::StringSubst(directory,"/core",""); //CO20200624 - prevent /home/corey -> /homey
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"common/SCINT","common/ICSD");       // [HISTORIC]
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"SCINT","ICSD");                     // [HISTORIC]
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"common/ELPASOLITES","common/AURO"); // [HISTORIC]
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"ELPASOLITES","AURO");               // [HISTORIC]
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"LIBRARYX/RAW","LIBRARYX/LIB");      // [HISTORIC]
      //[CO20200624 - do inside CleanDirectoryLIB()]aurostd::StringSubst(directory,"LIB2/RAW","LIB2/LIB");              // [HISTORIC]

      PROJECT_LIBRARY=aflowlib::LIB2RAW_CheckProjectFromDirectory(directory);
      // cout << soliloquy << " FOUND Project= " << XHOST.hostname << ": " << PROJECT_LIBRARY << endl;
      cout << soliloquy << " AFLOW V" << string(AFLOW_VERSION) << endl;
      cout << soliloquy << " directory=" << directory << endl;
      cout << soliloquy << " hostname=" << XHOST.hostname << endl;
      directory_LIB=directory;
      // wait for temperature (to avoid overheat)
      // WAIT FOR TEMPERATURE
      vector<string> vmessage;vmessage.push_back(XHOST.sPID);vmessage.push_back("       dir="+directory_LIB);
      init::WaitTEMP(AFLOW_CORE_TEMPERATURE_LIB2RAW,cout,TRUE,vmessage);  // fixed at define at the beginning
    }
    
    bool run_directory=false;
    string AflowInName = _AFLOWIN_;
    string FileLockName = _AFLOWLOCK_;

    //[CO20200624 - OBSOLETE]if(pocc::structuresGenerated(directory_LIB)){KBIN::VASP_RunPOCC(directory_LIB);}  //CO20200624
    //[CO20200624 - OBSOLETE]else if(aurostd::FileExist(directory_LIB+"/agl_aflow.in"))
    if(pocc::structuresGenerated(directory_LIB)){ //CO20200624
      run_directory=true;
      if(!aurostd::FileExist(directory_LIB+"/"+_AFLOWLOCK_+".pocc.preprocessing")){
        if(aurostd::FileExist(directory_LIB+"/"+_AFLOWLOCK_+".old")){aurostd::file2file(directory_LIB+"/"+_AFLOWLOCK_+".old",directory_LIB+"/"+_AFLOWLOCK_+".pocc.preprocessing");}
        else if(aurostd::FileExist(directory_LIB+"/"+_AFLOWLOCK_+".OLD")){aurostd::file2file(directory_LIB+"/"+_AFLOWLOCK_+".OLD",directory_LIB+"/"+_AFLOWLOCK_+".pocc.preprocessing");}
        else if(aurostd::FileExist(directory_LIB+"/"+_AFLOWLOCK_)){aurostd::file2file(directory_LIB+"/"+_AFLOWLOCK_,directory_LIB+"/"+_AFLOWLOCK_+".pocc.preprocessing");}
      }
    }
    else if(aurostd::FileExist(directory_LIB+"/agl_aflow.in")) {
      run_directory=true;
      for(uint iext=0;iext<XHOST.vext.size();iext++) {
        aurostd::RemoveFile(directory_LIB+"/aflow.agl.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL_edos_gap_pressure.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL_edos_gap_pressure.json"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL_energies_temperature.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL_energy.json"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL_energy_structures.json"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL_energy_volume.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL_gibbs_energy_pT.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL_Hugoniot.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL_thermal_properties_temperature.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AGL_THERMO"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/aflow.ael.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AEL_Compliance_tensor.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AEL_Elastic_constants.out"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AEL_elastic_tensor.json"+XHOST.vext.at(iext));
        aurostd::RemoveFile(directory_LIB+"/AEL_energy_structures.json"+XHOST.vext.at(iext));
      }
      // CORMAC FIX BELOW I NEED TO COMMENT TO RUN THE agl_aflow.in containing only statics.
      //   cerr << XPID << "CORMAC" << endl;exit(0);
      if(aurostd::FileExist(directory_LIB+"/agl.LOCK")) {
        // [OBSOLETE] aurostd::file2file(directory_LIB+"/agl.LOCK",directory_LIB+"/agl.LOCK.run");   // LINK
        // [OBSOLETE] aurostd::CopyFile(directory_LIB+"/agl.LOCK",directory_LIB+"/agl.LOCK.run");   // LINK
        // [OBSOLETE] aurostd::RemoveFile(directory_LIB+"/agl.LOCK");   //ME20191105 - otherwise aflow will not run
        FileLockName = "agl.LOCK";
      }
      if(aurostd::FileExist(directory_LIB+"/agl_aflow.in")) {
        AflowInName="agl_aflow.in";
      }
      //CT20200624 Calls functions to run AEL and AGL in postprocessing mode instead of executing an AFLOW run
      //CT20200624 This should help prevent VASP from running when performing postprocessing, since we go direct to AEL/AGL routines
      //[CO20200624 - OBSOLETE]KBIN::VASP_RunPhonons_AGL_postprocess(directory_LIB, AflowInName, FileLockName); //CT20200624 
      // cerr << XPID << "CORMAC" << endl;exit(0);
      // if(_AFLOWIN_=="aflow.in")
      if(aurostd::FileExist(directory_LIB+"/ael.LOCK")) {
        //[CO20200624 - do later]aurostd::CopyFile(directory_LIB+"/ael.LOCK",directory_LIB+"/ael.LOCK.run");   // LINK
        //[CO20200624 - do later]aurostd::RemoveFile(directory_LIB+"/ael.LOCK");   //ME20191105 - otherwise aflow will not run
        FileLockName = "ael.LOCK";
      }
      if(aurostd::FileExist(directory_LIB+"/ael_aflow.in")) {
        AflowInName="ael_aflow.in";
      }      
      //[CO20200624 - OBSOLETE]KBIN::VASP_RunPhonons_AEL_postprocess(directory_LIB, AflowInName, FileLockName); //CT20200624
      // [OBSOLETE] if (XHOST.QUIET) {      //CT20190903
      // [OBSOLETE]   aurostd::execute("aflow --use_aflow.in=agl_aflow.in --use_LOCK=agl.LOCK --force --run=0 --postprocess --quiet -D \""+directory_LIB+"\""); // do not mess up subdirectories
      // [OBSOLETE] } else {
      // [OBSOLETE]   aurostd::execute("aflow --use_aflow.in=agl_aflow.in --use_LOCK=agl.LOCK --force --run=0 --postprocess -D \""+directory_LIB+"\"");	
      // [OBSOLETE] }
    }

    //KBIN::RUN_Directory() should be used instead of individual post-processing functions
    //otherwise we will have to program every single instance here
    //this also ensures we compress everything at the end
    if(run_directory){
      _aflags aflags;aflags.Directory=directory_LIB;
      
      //save originals
      string _AFLOWIN_orig=_AFLOWIN_;
      string _AFLOWLOCK_orig=_AFLOWLOCK_;
      
      //set env for RUN_Directory()
      _AFLOWIN_=AflowInName;
      _AFLOWLOCK_=FileLockName;
      if(aurostd::FileExist(directory_LIB+"/"+_AFLOWLOCK_)){aurostd::file2file(directory_LIB+"/"+_AFLOWLOCK_,directory_LIB+"/"+_AFLOWLOCK_+".run");} //keep original LOCK
      KBIN::RUN_Directory(aflags);

      //return to original
      _AFLOWIN_=_AFLOWIN_orig;
      _AFLOWLOCK_=_AFLOWLOCK_orig;
    }

    if(LDEBUG) cerr << soliloquy << " END" << endl;
    return TRUE;
  }
}

#endif //  _AFLOWLIB_LIBRARIES_CPP_
// ***************************************************************************
