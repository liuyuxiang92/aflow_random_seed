// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo - 2007-2020
#ifndef _AFLOW_XPSEUDOPOTENTIAL_CPP
#define _AFLOW_XPSEUDOPOTENTIAL_CPP
#include "aflow.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// _XPSEUDOPOTENTIAL_PROTOTYPES_
/*
  # commands to cut/paste to generate stuff
  rm -f aflow_xpseudopotentials_data.cpp
  touch aflow_xpseudopotentials.cpp aflow_xpseudopotentials_data.cpp
  make -j
  rm -f aflow_xpseudopotentials_data.cpp
  ./aflow --scrub=POTCAR --FILE `find /common/VASP/pot_LDA/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
  ./aflow --scrub=POTCAR --FILE `find /common/VASP/pot_GGA/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
  ./aflow --scrub=POTCAR --FILE `find /common/VASP/potpaw_LDA/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
  ./aflow --scrub=POTCAR --FILE `find /common/VASP/potpaw_GGA/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
  ./aflow --scrub=POTCAR --FILE `find /common/VASP/potpaw_PBE/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
  ./aflow --scrub=POTCAR --FILE `find /common/VASP/potpaw_LDA.54/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
  ./aflow --scrub=POTCAR --FILE `find /common/VASP/potpaw_PBE.54/ -name POTCAR |sort ` >> aflow_xpseudopotentials_data.cpp
  touch aflow_xpseudopotentials.cpp aflow_xpseudopotentials_data.cpp
  # subst "/common/VASP/" "" aflow_xpseudopotentials_data.cpp
  cat aflow_xpseudopotentials_data.cpp | grep -c TITEL
  cat aflow_xpseudopotentials_data.cpp | grep structure | grep -v N/A | wc
  make -j

  ./aflow --scrub=POTCAR --FILE /common/VASP/potpaw_PBE/current/Mo_pv/POTCAR
  ./aflow --scrub=OUTCAR --FILE /common/LIB3/LIB/AgCdCo/TFCC001.ABC/OUTCAR.relax2.xz
  ./aflow --pseudopotentials_check=/tmp/POTCAR1
  ./aflow --pseudopotentials_check=/tmp/POTCAR2
  ./aflow --pseudopotentials_check=/tmp/OUTCAR.relax2
  ./aflow --pseudopotentials_check=/common/LIB3/LIB/AgCdCo/TFCC001.ABC/OUTCAR.relax2.xz

  ./aflow --use_aflow.in=aflow.in --beep --force --lib2raw="/common/LIB3/LIB/TeW_pvY_sv/TFCC001.ABC"

  #!/bin/sh
  #echo "$1"
  #echo "$2"
  STR1=`cat "/common/LIB1/LIB/$1/A1/aflow.in" | grep AUID | head -1 | sed "s/\[VASP_POTCAR_AUID\]/if(AUID==\""/g`
  STR2="\") {"$2"}   //   "$1
  echo $STR1$STR2

*/

std::vector<xPOTCAR> vxpseudopotential;        // store starting from ONE

#define PSEUDOPOTENTIAL_GENERATOR_pad 70

bool xPOTCAR_FixBoot(xPOTCAR& xPOT) {
  bool fix=FALSE;
  if(!xPOT.vENMAX.size()) { fix=TRUE; xPOT.vENMAX.push_back(NNN);}                                    // if not identified for BOOT STRAP
  if(!xPOT.vENMIN.size()) { fix=TRUE; xPOT.vENMIN.push_back(NNN);}                                    // if not identified for BOOT STRAP
  if(!xPOT.vPOMASS.size()) { fix=TRUE; xPOT.vPOMASS.push_back(NNN);}                                  // if not identified for BOOT STRAP
  if(!xPOT.vZVAL.size()) { fix=TRUE; xPOT.vZVAL.push_back(NNN);}                                      // if not identified for BOOT STRAP
  if(!xPOT.vEATOM.size()) { fix=TRUE; xPOT.vEATOM.push_back(NNN);}                                    // if not identified for BOOT STRAP
  if(!xPOT.vRCORE.size()) { fix=TRUE; xPOT.vRCORE.push_back(NNN);}                                    // if not identified for BOOT STRAP
  if(!xPOT.vRWIGS.size()) { fix=TRUE; xPOT.vRWIGS.push_back(NNN);}                                    // if not identified for BOOT STRAP
  if(!xPOT.vEAUG.size()) { fix=TRUE; xPOT.vEAUG.push_back(NNN);}                                      // if not identified for BOOT STRAP
  if(!xPOT.vRAUG.size()) { fix=TRUE; xPOT.vRAUG.push_back(NNN);}                                      // if not identified for BOOT STRAP
  if(!xPOT.vRMAX.size()) { fix=TRUE; xPOT.vRMAX.push_back(NNN);}                                      // if not identified for BOOT STRAP
  if(!xPOT.vTITEL.size()) { fix=TRUE; xPOT.vTITEL.push_back("N/A");}                                  // if not identified for BOOT STRAP
  if(!xPOT.vLEXCH.size()) { fix=TRUE; xPOT.vLEXCH.push_back("N/A");}                                  // if not identified for BOOT STRAP
  if(!xPOT.species.size()) { fix=TRUE; xPOT.species.push_back("N/A");}                                // if not identified for BOOT STRAP
  if(!xPOT.species_Z.size()) { fix=TRUE; xPOT.species_Z.push_back(0);}                                // if not identified for BOOT STRAP
  if(!xPOT.species_pp.size()) { fix=TRUE; xPOT.species_pp.push_back("N/A");}                          // if not identified for BOOT STRAP
  if(!xPOT.species_pp_type.size()) { fix=TRUE; xPOT.species_pp_type.push_back("N/A");}                // if not identified for BOOT STRAP
  if(!xPOT.species_pp_version.size()) { fix=TRUE; xPOT.species_pp_version.push_back("N/A");}          // if not identified for BOOT STRAP
  if(!xPOT.species_pp_AUID.size()) { fix=TRUE; xPOT.species_pp_AUID.push_back("N/A");}                // if not identified for BOOT STRAP
  if(!xPOT.species_pp_groundstate_energy.size()) { fix=TRUE; xPOT.species_pp_groundstate_energy.push_back(NNN);}          // if not identified for BOOT STRAP
  if(!xPOT.species_pp_groundstate_structure.size()) { fix=TRUE; xPOT.species_pp_groundstate_structure.push_back("N/A");}  // if not identified for BOOT STRAP
  return fix;
}

xPOTCAR xPOTCAR_Finder(vector<string>& species_pp_AUID,vector<string>& species_pp_AUID_collisions,const string& TITEL,const string& LEXCH,const double& EATOM,const double& RMAX,bool LVERBOSE) {
  xPOTCAR xPOT;
  bool found=FALSE;
  for(uint ipp=0;ipp<vxpseudopotential.size();ipp++) {
    bool test=
      (TITEL==vxpseudopotential.at(ipp).vTITEL.at(0)) &&
      (LEXCH==vxpseudopotential.at(ipp).vLEXCH.at(0)) &&
      (aurostd::abs(EATOM-vxpseudopotential.at(ipp).vEATOM.at(0))<0.00001) &&
      (aurostd::abs(RMAX-vxpseudopotential.at(ipp).vRMAX.at(0))<0.0001);
    if(test && found) {
      if(!aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Xe") &&
	 !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Kr") &&
	 !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Li_pv") &&
	 !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_LDA/potcar.Apr00/Li_pv")) {
	//if(LVERBOSE)
	cerr << "xPOTCAR::xPOTCAR_Finder: COLLISION: POTCAR " << vxpseudopotential.at(ipp).filename << " " << xPOT.filename << endl;
	species_pp_AUID_collisions.push_back(vxpseudopotential.at(ipp).AUID);
      }
    }
    if(test && !found) {
      found=TRUE;
      if(LVERBOSE) cerr << "xPOTCAR::xPOTCAR_Finder: FOUND: POTCAR=" << vxpseudopotential.at(ipp).filename << endl;
      species_pp_AUID.push_back(vxpseudopotential.at(ipp).AUID);
      xPOT=vxpseudopotential.at(ipp);
    }
  }
  if(!found) {
    if(vxpseudopotential.size()) {
      // if(LVERBOSE)
      cerr << "xPOTCAR::xPOTCAR_Finder: NOT FOUND: TITEL=" << TITEL << " LEXCH=" <<  LEXCH << " EATOM=" <<  EATOM << " RMAX =" <<  RMAX << endl;
      //    cout << "xPOTCAR::xPOTCAR_Finder: NOT FOUND: TITEL=" << TITEL << " LEXCH=" <<  LEXCH << " EATOM=" <<  EATOM << " RMAX =" <<  RMAX << endl;
      // exit(0);
    }
    species_pp_AUID.push_back("N/A");
  }

  xPOTCAR_FixBoot(xPOT);
  return xPOT;
}

xPOTCAR xPOTCAR_Finder(const string& AUID,bool LVERBOSE) {
  xPOTCAR xPOT;
  bool found=FALSE;
  for(uint ipp=0;ipp<vxpseudopotential.size();ipp++) {
    bool test= (AUID==vxpseudopotential.at(ipp).AUID);
    if(test && found) {
      if(!aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Xe") &&
	 !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Kr") &&
	 !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_GGA/potcar.Apr00/Li_pv") &&
	 !aurostd::substring2bool(vxpseudopotential.at(ipp).filename,"pot_LDA/potcar.Apr00/Li_pv")) {
	//if(LVERBOSE)
	cerr << "xPOTCAR::xPOTCAR_Finder: COLLISION: POTCAR " << vxpseudopotential.at(ipp).filename << " " << xPOT.filename << endl;
      }
    }
    if(test && !found) {
      found=TRUE;
      if(LVERBOSE) cerr << "xPOTCAR::xPOTCAR_Finder: FOUND: POTCAR=" << vxpseudopotential.at(ipp).filename << endl;
      xPOT=vxpseudopotential.at(ipp);
    }
  }
  if(!found) {
    if(vxpseudopotential.size()) {
      // if(LVERBOSE)
      cerr << "xPOTCAR::xPOTCAR_Finder: NOT FOUND: AUID=" << AUID << endl;
      //   cout << "xPOTCAR::xPOTCAR_Finder: NOT FOUND: AUID=" << AUID << endl;
      //    exit(0);
    }
  }

  xPOTCAR_FixBoot(xPOT);
  return xPOT;
}

bool xPOTCAR_PURE_Printer(xPOTCAR& xPOT,ostream& oss,bool LVERBOSE) {
  if(XHOST.PSEUDOPOTENTIAL_GENERATOR && xPOT.species.size()==1) {  //SC20200326
    string comment="";
    // objects/functions for references energies  //SC20200326
    //   vector<string> tokens;
    //   aurostd::string2tokens(xPOT.species_pp_version.at(0),tokens,":");
    // [OBSOLETE] if(tokens.size()>0) {vdate.clear();vdate.push_back(tokens.at(tokens.size()-1));};

    double groundstate_energy=NNN;
    string groundstate_structure="N/A_"+xPOT.AUID;

    double volume_atom;
    double spin_atom;

    bool found=xPOTCAR_EnthalpyReference_AUID(xPOT.AUID,"",groundstate_structure,groundstate_energy,volume_atom,spin_atom);

    if(found)  cerr << "xPOTCAR_PURE_Printer:     FOUND AUID=" << xPOT.AUID << endl;
    if(!found) cerr << "xPOTCAR_PURE_Printer: NOT_FOUND AUID=" << xPOT.AUID << endl;

    comment=xPOT.species_pp_version.at(0);
    xPOT.vTITEL.at(0)=aurostd::RemoveWhiteSpaces(xPOT.vTITEL.at(0));
    xPOT.vLEXCH.at(0)=aurostd::RemoveWhiteSpaces(xPOT.vLEXCH.at(0));
    if(LVERBOSE) {cerr << "xPOTCAR_PURE_Printer: vTITEL.at(0)=" << xPOT.vTITEL.at(0) << endl;}   //SC20200326
    if(LVERBOSE) {cerr << "xPOTCAR_PURE_Printer: species.at(0)=" << xPOT.species.at(0) << endl;}   //SC20200326
    if(LVERBOSE) {cerr << "xPOTCAR_PURE_Printer: species_Z.at(0)=" << xPOT.species_Z.at(0) << endl;}    //SC20200326
    oss << "  " << endl;
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    //   oss << "  // [AFLOW]START=" << comment << " " << endl;
    oss << "  // " << comment << " " << comment << " " << comment << " " << comment << " " << endl;
    oss << "  // " << xPOT.filename << endl;    //SC20200326
    oss << "  " << aurostd::PaddedPOST("{",PSEUDOPOTENTIAL_GENERATOR_pad) << "      // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("xPOTCAR x;",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.filename=\""+xPOT.filename+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.AUID=\""+xPOT.AUID+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vTITEL.push_back(\""+xPOT.vTITEL.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.pp_type=\""+xPOT.pp_type+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species.push_back(\""+xPOT.species.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_Z.push_back("+aurostd::utype2string<int>(xPOT.species_Z.at(0))+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp.push_back(\""+xPOT.species_pp.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_type.push_back(\""+xPOT.species_pp_type.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_version.push_back(\""+xPOT.species_pp_version.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_AUID.push_back(\""+xPOT.AUID+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    if(!found) oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_energy.push_back("+aurostd::utype2string<double>(groundstate_energy,10)+");//"+xPOT.AUID,PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    if(found)  oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_energy.push_back("+aurostd::utype2string<double>(groundstate_energy,10)+");//",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_structure.push_back(\""+groundstate_structure+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_PAW="+aurostd::bool2string(xPOT.POTCAR_PAW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_TYPE=\""+xPOT.POTCAR_TYPE+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_KINETIC="+aurostd::bool2string(xPOT.POTCAR_KINETIC)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_GW="+aurostd::bool2string(xPOT.POTCAR_GW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_AE="+aurostd::bool2string(xPOT.POTCAR_AE)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // EATOM RCORE RWIGS EAUG RAUG ENMAX ENMIN POMASS ZVAL RMAX LEXCH
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vENMAX.push_back("+aurostd::utype2string<double>(xPOT.vENMAX.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vENMIN.push_back("+aurostd::utype2string<double>(xPOT.vENMIN.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vPOMASS.push_back("+aurostd::utype2string<double>(xPOT.vPOMASS.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vZVAL.push_back("+aurostd::utype2string<double>(xPOT.vZVAL.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vEATOM.push_back("+aurostd::utype2string<double>(xPOT.vEATOM.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vRCORE.push_back("+aurostd::utype2string<double>(xPOT.vRCORE.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vRWIGS.push_back("+aurostd::utype2string<double>(xPOT.vRWIGS.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vEAUG.push_back("+aurostd::utype2string<double>(xPOT.vEAUG.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vRAUG.push_back("+aurostd::utype2string<double>(xPOT.vRAUG.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vRMAX.push_back("+aurostd::utype2string<double>(xPOT.vRMAX.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vLEXCH.push_back(\""+xPOT.vLEXCH.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vxpseudopotential.push_back(x);",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "  " << aurostd::PaddedPOST("}",PSEUDOPOTENTIAL_GENERATOR_pad) << "      // " << comment << endl;    //SC20200326
    oss << "  // " << xPOT.filename << endl;    //SC20200326
    //   oss << "  // [AFLOW]STOP=" << comment << " " << endl;
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    oss << endl;    //SC20200326
    return TRUE;
  }
  return FALSE;
}

ostream& operator<<(ostream& oss,const xPOTCAR& xPOT) {
  oss.setf(std::ios::fixed,std::ios::floatfield);
  oss.precision(10);
  for(uint i=0;i<xPOT.species.size();i++) {
    string comment=xPOT.species_pp_version.at(i);
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    oss << "  // [AFLOW]START=" << comment << " " << endl;
    oss << "  // " << comment << " " << comment << " " << comment << " " << comment << " " << endl;
    oss << "  // " << xPOT.filename << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("filename=\""+xPOT.filename+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("AUID=\""+xPOT.AUID+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vTITEL.push_back(\""+xPOT.vTITEL.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("pp_type=\""+xPOT.pp_type+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species.push_back(\""+xPOT.species.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_Z.push_back("+aurostd::utype2string<int>(xPOT.species_Z.at(i))+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp.push_back(\""+xPOT.species_pp.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_type.push_back(\""+xPOT.species_pp_type.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_version.push_back(\""+xPOT.species_pp_version.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_AUID.push_back(\""+xPOT.species_pp_AUID.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_groundstate_energy.push_back("+aurostd::utype2string<double>(xPOT.species_pp_groundstate_energy.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_groundstate_structure.push_back(\""+xPOT.species_pp_groundstate_structure.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_PAW="+aurostd::bool2string(xPOT.POTCAR_PAW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_TYPE=\""+xPOT.POTCAR_TYPE+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_KINETIC="+aurostd::bool2string(xPOT.POTCAR_KINETIC)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_GW="+aurostd::bool2string(xPOT.POTCAR_GW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_AE="+aurostd::bool2string(xPOT.POTCAR_AE)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vENMAX.push_back("+aurostd::utype2string<double>(xPOT.vENMAX.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vENMIN.push_back("+aurostd::utype2string<double>(xPOT.vENMIN.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vPOMASS.push_back("+aurostd::utype2string<double>(xPOT.vPOMASS.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vZVAL.push_back("+aurostd::utype2string<double>(xPOT.vZVAL.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vEATOM.push_back("+aurostd::utype2string<double>(xPOT.vEATOM.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vRCORE.push_back("+aurostd::utype2string<double>(xPOT.vRCORE.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vRWIGS.push_back("+aurostd::utype2string<double>(xPOT.vRWIGS.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vEAUG.push_back("+aurostd::utype2string<double>(xPOT.vEAUG.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;   //SC20200326
    oss << "    " << aurostd::PaddedPOST("vRAUG.push_back("+aurostd::utype2string<double>(xPOT.vRAUG.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;   //SC20200326
    oss << "    " << aurostd::PaddedPOST("vRMAX.push_back("+aurostd::utype2string<double>(xPOT.vRMAX.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "    " << aurostd::PaddedPOST("vLEXCH.push_back(\""+xPOT.vLEXCH.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    //SC20200326
    oss << "  // " << xPOT.filename << endl;    //SC20200326
    oss << "  // [AFLOW]STOP=" << comment << " " << endl;
    oss << "  // ******************************************************************************************************************************************************** " << endl;
  }
  return oss;
}

uint xPOTCAR_Initialize(void) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  uint COLLISIONS=0;
  if(LDEBUG) cerr << "xpseudopotentials_initialize [BEGIN]" << endl;
  long double start=aurostd::get_useconds();
  //  if(LDEBUG) cerr << "xpseudopotentials_initialize: " << start << endl;
  if(LDEBUG) cerr << "vxpseudopotential.size()=" << vxpseudopotential.size() << endl;
#include "aflow_xpseudopotentials_data.cpp"
  long double stop=aurostd::get_useconds();
  if(LDEBUG) cerr << "xpseudopotentials_initialize: " << (stop-start)/1000  << endl;
  if(LDEBUG) cerr << "vxpseudopotential.size()=" << vxpseudopotential.size() << endl;
  if(LDEBUG) cerr << "xpseudopotentials_initialize [END]" << endl;

  for(uint i=0;i<vxpseudopotential.size();i++) {
    for(uint j=i+1;j<vxpseudopotential.size();j++) {
      if(vxpseudopotential.at(i).AUID!=vxpseudopotential.at(j).AUID) {
	if(vxpseudopotential.at(i).vTITEL.at(0)==vxpseudopotential.at(j).vTITEL.at(0)) {
	  if(vxpseudopotential.at(i).vLEXCH.at(0)==vxpseudopotential.at(j).vLEXCH.at(0)) {
	    if(abs(vxpseudopotential.at(i).vEATOM.at(0)-vxpseudopotential.at(j).vEATOM.at(0))<0.00001) {
	      if(abs(vxpseudopotential.at(i).vRMAX.at(0)-vxpseudopotential.at(j).vRMAX.at(0))<0.0001) {
		cerr << "COLLISION" << " "
		     << "FILENAME " << vxpseudopotential.at(i).filename << " " << vxpseudopotential.at(j).filename << " "
		     << "TITEL " << vxpseudopotential.at(i).vTITEL.at(0) << " "
		     << "LEXCH " << vxpseudopotential.at(i).vLEXCH.at(0) << " "
		     << "EATOM " << vxpseudopotential.at(i).vEATOM.at(0) << " "
		     << "RMAX " << vxpseudopotential.at(i).vRMAX.at(0) << " "
		     << endl;
		COLLISIONS++;
	      }
	    }
	  }
	}
      }
    }
  }
  //  cerr << "COLLISIONS=" << COLLISIONS << endl;
  return COLLISIONS;
}

bool xPOTCAR_EnthalpyReference_AUID(string AUID,string METAGGA) {
  string groundstate_structure="";
  double groundstate_energy=0.0,volume_atom=0.0,spin_atom=0.0;
  return  xPOTCAR_EnthalpyReference_AUID(AUID,METAGGA,groundstate_structure,groundstate_energy,volume_atom,spin_atom);
}

bool xPOTCAR_EnthalpyReference_AUID(string AUID,string METAGGA,string& groundstate_structure,double& groundstate_energy,double& volume_atom,double& spin_atom) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  bool VERBOSE=0;//TRUE;
  if(LDEBUG) cerr << "xPOTCAR_EnthalpyReference_AUID: [BEGIN]" << endl;
  bool found=FALSE;
  if(LDEBUG) cout <<"ERROR (xPOTCAR_EnthalpyReference_AUID): AUID=[" << AUID << "]" << endl; // exit(0);
  if(LDEBUG) cout <<"ERROR (xPOTCAR_EnthalpyReference_AUID): METAGGA=[" << METAGGA << "]" << endl; // exit(0);

  bool nKIN=FALSE,SCAN=FALSE;
  if(METAGGA.empty() || METAGGA=="none" || METAGGA=="NONE") {nKIN=TRUE;SCAN=FALSE;}
  if(METAGGA=="SCAN" || METAGGA=="scan") {nKIN=FALSE;SCAN=TRUE;}

  if(!XHOST.PSEUDOPOTENTIAL_GENERATOR) if(VERBOSE) cout <<"xPOTCAR_EnthalpyReference_AUID: AUID=[" << AUID << "]  METAGGA=[" << METAGGA << "]" << endl;

  // Ac
  if(AUID=="335f20e70f06b78b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.09443;volume_atom=45.4098;spin_atom=0.0;} // Ac:PAW_PBE:06Sep2000
  if(AUID=="b5bf833f40cc220c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.04579;volume_atom=44.9982;spin_atom=0.0;} // Ac_s:PAW_GGA:11Apr2000
  if(AUID=="f04c8b0982efae26" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.06438;volume_atom=45.1122;spin_atom=0.0;} // Ac:PAW_GGA:11Apr2000
  if(AUID=="f4d16e5d594ed3ef" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.64874;volume_atom=40.4514;spin_atom=0.0;} // Ac:PAW_LDA:12Apr2000
  if(AUID=="863169f4af8ed89f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.63742;volume_atom=40.4895;spin_atom=0.0;} // Ac_s:PAW_LDA:12Apr2000
  if(AUID=="06e8eff4a6b826a0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.03196;volume_atom=44.7612;spin_atom=0.0;} // Ac_s:PAW_PBE:06Sep2000
  if(AUID=="b35f9d074448e332" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-90.4912;volume_atom=43.9831;spin_atom=0.0;} // Ac:PAW_PBE_KIN:SCAN:06Sep2000

  // Ag
  if(AUID=="dc0f72d3f1bd620d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.72686;volume_atom=17.7851;spin_atom=0.0;} // Ag:PAW_GGA:18Jul2000
  if(AUID=="eca65d7d992efb48" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.82769;volume_atom=17.9065;spin_atom=0.0;} // Ag:PAW_PBE:06Sep2000
  if(AUID=="dac6684dedca50c0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.72042;volume_atom=17.9853;spin_atom=0.0;} // Ag:GGA:01Apr2000
  if(AUID=="639bb9917452d998" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.74450;volume_atom=16.2004;spin_atom=0.0;} // Ag:LDA:01Apr2000
  if(AUID=="9144029176631616" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.74599;volume_atom=16.0612;spin_atom=0.0;} // Ag:PAW_LDA:17Apr2000
  if(AUID=="cc0097ef8847f6c4" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-31.3840;volume_atom=16.9292;spin_atom=0.0;} // Ag:PAW_PBE_KIN:SCAN:02Apr2005
  if(AUID=="f4452364bfcfa5d4" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-31.3515;volume_atom=16.8463;spin_atom=0.0;} // Ag_pv:PAW_PBE_KIN:SCAN:09Dec2005

  // Al
  if(AUID=="e7afb3f96db37614" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.69194;volume_atom=16.5273;spin_atom=0.0;} // Al:GGA:01Apr2000
  if(AUID=="a0a3e933f8f5dcc7" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.72092;volume_atom=16.3640;spin_atom=0.0;} // Al:GGA:01Apr2000
  if(AUID=="ec638285f61e73c3" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.23100;volume_atom=15.6073;spin_atom=0.0;} // Al:LDA:01Apr2000
  if(AUID=="b9cdcfa63ddf1ea9" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.21254;volume_atom=15.1208;spin_atom=0.0;} // Al_h:PAW_GGA:08Apr2002
  if(AUID=="de9492ed3ab234ca" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.80021;volume_atom=14.8772;spin_atom=0.0;} // Al_h:PAW_PBE:08Apr2002
  if(AUID=="f2de0641bce99445" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.69580;volume_atom=16.5600;spin_atom=0.0;} // Al:PAW_GGA:05Jan2001
  if(AUID=="2b253cb3349b4835" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.19184;volume_atom=15.8018;spin_atom=0.0;} // Al:PAW_LDA:17Apr2000
  if(AUID=="624cbf44af81780c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.19464;volume_atom=15.7463;spin_atom=0.0;} // Al:LDA:01Apr2000
  if(AUID=="423324aa45d0242f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.74356;volume_atom=16.4719;spin_atom=0.0;} // Al:PAW_PBE:04Jan2001
  if(AUID=="35865a9b2b516833" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.74846;volume_atom=16.1309;spin_atom=0.0;} // Al:PAW_PBE_KIN:SCAN:04Jan2001

  // Ar
  if(AUID=="5ea4f89868f4ce2e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.145017;volume_atom=29.6839;spin_atom=0.0;} // Ar:PAW_LDA:07Sep2000
  if(AUID=="5adcf9fd506a6f8f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.079570;volume_atom=50.7986;spin_atom=0.0;} // Ar:PAW_GGA:06Sep2000
  if(AUID=="76ed84ac7c922a74" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.065575;volume_atom=44.0826;spin_atom=0.0;} // Ar:PAW_PBE:07Sep2000
  if(AUID=="c50dbc48fec4befc" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.917290;volume_atom=36.2195;spin_atom=0.0;} // Ar:PAW_PBE_KIN:SCAN:07Sep2000

  // As
  if(AUID=="f1de4433a638eae7" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-4.65243;volume_atom=22.6706;spin_atom=0.0;} // As:PAW_PBE:06Sep2000
  if(AUID=="1170cc199ae6ad79" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-20.2444;volume_atom=21.5214;spin_atom=0.0;} // As_d:PAW_PBE_KIN:SCAN:11Apr2003
  if(AUID=="b7b512c9cb6a6957" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-20.1938;volume_atom=21.8855;spin_atom=0.0;} // As:PAW_PBE_KIN:SCAN:22Sep2009

  // Au
  if(AUID=="57972e71335a410c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.20315;volume_atom=18.0489;spin_atom=0.0;} // Au:PAW_GGA:18Jul2000
  if(AUID=="c42e6bdfe7e7024e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.38948;volume_atom=16.7676;spin_atom=0.0;} // Au:LDA:01Apr2000
  if(AUID=="3477c9645a808524" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.27212;volume_atom=18.0368;spin_atom=0.0;} // Au:PAW_PBE:06Sep2000
  if(AUID=="fdbc4d4856dea2ca" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.19304;volume_atom=18.1997;spin_atom=0.0;} // Au:GGA:01Apr2000
  if(AUID=="834bf85c5ba24d77" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.39641;volume_atom=16.6632;spin_atom=0.0;} // Au:PAW_LDA:04Feb1998
  if(AUID=="e2936a3049f7d935" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-69.2026;volume_atom=17.0078;spin_atom=0.0;} // Au:PAW_PBE_KIN:SCAN:04Oct2007

  /*
  // B
  aflow --aflow_proto=ICSD_56992.A:B --potential=pot_LDA --potential_complete --run_relax_static_bands --kppra=32768  --kppra_static=32768
  aflow --aflow_proto=ICSD_56992.A:B --potential=pot_GGA --potential_complete --run_relax_static_bands --kppra=32768  --kppra_static=32768
  aflow --aflow_proto=ICSD_56992.A:B --potential=potpaw_LDA --potential_complete --run_relax_static_bands --kppra=32768  --kppra_static=32768
  aflow --aflow_proto=ICSD_56992.A:B --potential=potpaw_GGA --potential_complete --run_relax_static_bands --kppra=32768  --kppra_static=32768
  aflow --aflow_proto=ICSD_56992.A:B --potential=potpaw_PBE --potential_complete --run_relax_static_bands --kppra=32768  --kppra_static=32768

  aflow --aflow_proto=ICSD_94429.A:B --potential=pot_LDA --potential_complete --run_relax_static_bands --kppra=32768  --kppra_static=32768
  aflow --aflow_proto=ICSD_94429.A:B --potential=pot_GGA --potential_complete --run_relax_static_bands --kppra=32768  --kppra_static=32768
  aflow --aflow_proto=ICSD_94429.A:B --potential=potpaw_LDA --potential_complete --run_relax_static_bands --kppra=32768  --kppra_static=32768
  aflow --aflow_proto=ICSD_94429.A:B --potential=potpaw_GGA --potential_complete --run_relax_static_bands --kppra=32768  --kppra_static=32768
  aflow --aflow_proto=ICSD_94429.A:B --potential=potpaw_PBE --potential_complete --run_relax_static_bands --kppra=32768  --kppra_static=32768

B_h:PAW_GGA:18Jul2000/ICSD_56992.A/ 1 F= -.80301674E+02 E0= -.80301674E+02 d E =0.000000E+00 mag= -0.0000
B_h:PAW_GGA:18Jul2000/ICSD_94429.A/ 1 F= -.80301674E+02 E0= -.80301674E+02 d E =0.000000E+00 mag= -0.0000

B_h:PAW_LDA:17Apr2000/ICSD_56992.A/ 1 F= -.89629659E+02 E0= -.89629659E+02 d E =0.000000E+00 mag= 0.0000
B_h:PAW_LDA:17Apr2000/ICSD_94429.A/ 1 F= -.89629658E+02 E0= -.89629658E+02 d E =0.000000E+00 mag= 0.0000

B_h:PAW_PBE:07Sep2000/ICSD_56992.A/ 1 F= -.80319089E+02 E0= -.80319089E+02 d E =0.000000E+00 mag= -0.0000
B_h:PAW_PBE:07Sep2000/ICSD_94429.A/ 1 F= -.80319089E+02 E0= -.80319089E+02 d E =0.000000E+00 mag= -0.0000
 */

  // B
  if(AUID=="76781ebe8489383f" && nKIN) {found=TRUE;groundstate_structure="ICSD_56992";groundstate_energy=-6.69326;volume_atom=7.24167;spin_atom=0.0;} // B_h:PAW_PBE:07Sep2000
  if(AUID=="70110ee6c6cbaf90" && nKIN) {found=TRUE;groundstate_structure="ICSD_56992";groundstate_energy=-7.46914;volume_atom=6.99838;spin_atom=0.0;} // B_h:PAW_LDA:17Apr2000
  if(AUID=="00b4dfcc28b5887b" && nKIN) {found=TRUE;groundstate_structure="ICSD_56992";groundstate_energy=-6.69181;volume_atom=7.24862;spin_atom=0.0;} // B_h:PAW_GGA:18Jul2000
  // A3  if(AUID=="97c23347b69407ee" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.64898;volume_atom=8.99490;spin_atom=0.0;} // B_h:PAW_PBE_KIN:SCAN:06Feb2004
    
  // Ba
  if(AUID=="062068333d4e7e80" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.92400;volume_atom=63.1759;spin_atom=0.0;} // Ba_sv:PAW_PBE:06Sep2000
  if(AUID=="0ea886f29a8f18e6" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.92401;volume_atom=62.1307;spin_atom=0.0;} // Ba_sv:PAW_GGA:14Apr2000
  if(AUID=="b5d1dafcf2f8d798" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.19536;volume_atom=55.5803;spin_atom=0.0;} // Ba:LDA:01Apr2000
  if(AUID=="45efe879c954c89e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.89600;volume_atom=62.4601;spin_atom=0.0;} // Ba:GGA:01Apr2000
  if(AUID=="7cbe610c465db194" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.23839;volume_atom=53.9574;spin_atom=0.0;} // Ba_sv:PAW_LDA:17Apr2000
  if(AUID=="63fd0fe091a069b0" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-39.5316;volume_atom=63.8147;spin_atom=0.0;} // Ba_sv:PAW_PBE_KIN:SCAN:06Sep2000

  // Be
  if(AUID=="abd9038ab359fe8e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.19956;volume_atom=7.60317;spin_atom=0.0;} // Be_sv:PAW_LDA:23Feb1998
  if(AUID=="80136f0733d30e11" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.74102;volume_atom=7.92019;spin_atom=0.0;} // Be_sv:PAW_PBE:06Sep2000
  if(AUID=="250ccb49102b1ab7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.75366;volume_atom=7.92077;spin_atom=0.0;} // Be:PAW_PBE:06Sep2000
  if(AUID=="362b70d9e55ede03" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.72800;volume_atom=7.91229;spin_atom=0.0;} // Be:PAW_GGA:11Feb1998
  if(AUID=="0134da61389a9f60" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-3.73743;volume_atom=7.84295;spin_atom=0.0;} // Be:GGA:01Apr2000
  if(AUID=="653b8d1063d6b4e0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.19996;volume_atom=7.57850;spin_atom=0.0;} // Be:LDA:01Apr2000
  if(AUID=="d89972c63683236a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.21025;volume_atom=7.57396;spin_atom=0.0;} // Be:PAW_LDA:02Feb1998
  if(AUID=="b72167bbb69596d1" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.36112;volume_atom=7.92170;spin_atom=0.0;} // Be:PAW_PBE_KIN:SCAN:06Sep2000

  // Bi
  if(AUID=="8b1158914bff5430" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-3.87274;volume_atom=36.8656;spin_atom=0.0;} // Bi:PAW_PBE:08Apr2002
  if(AUID=="bba01b714d812459" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-4.03716;volume_atom=36.2786;spin_atom=0.0;} // Bi_d:PAW_PBE:06Sep2000
  if(AUID=="2efd860034c70e39" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-78.5112;volume_atom=35.6425;spin_atom=0.0;} // Bi:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="46860e6586f9539f" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-78.5664;volume_atom=34.5043;spin_atom=0.0;} // Bi_d:PAW_PBE_KIN:SCAN:06Sep2000

  // Br
  if(AUID=="c53210bda32c7827" && nKIN) {found=TRUE;groundstate_structure="A11";groundstate_energy=-1.5898;volume_atom=45.5907;spin_atom=0.0;} // Br:PAW_PBE:06Sep2000

  // C
  if(AUID=="f3800ea7e3e20d86" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.19568;volume_atom=10.6949;spin_atom=0.0;} // C_h:PAW_PBE:20Dec2001
  if(AUID=="a720f5418b5ad14f" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.22034;volume_atom=10.4453;spin_atom=0.0;} // C:PAW_PBE:08Apr2002
  if(AUID=="8d30df5759a1e14d" && nKIN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-9.19508;volume_atom=10.5847;spin_atom=0.0;} // C_s:PAW_PBE:06Sep2000
  if(AUID=="8397f0cab3a0348d" && SCAN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-10.1261;volume_atom=8.99501;spin_atom=0.0;} // C_h:PAW_PBE_KIN:SCAN:06Feb2004
  if(AUID=="f9bce748fcbf37b3" && SCAN) {found=TRUE;groundstate_structure="A9";groundstate_energy=-10.0941;volume_atom=9.02062;spin_atom=0.0;} // C:PAW_PBE_KIN:SCAN:08Apr2002

  // Ca
  if(AUID=="3e3d86ebf8f9e84a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.93199;volume_atom=41.5802;spin_atom=0.0;} // Ca:GGA:01Apr2000
  if(AUID=="572c23d8cdb07c91" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.94755;volume_atom=41.3321;spin_atom=0.0;} // Ca:GGA:01Apr2000
  if(AUID=="7f030fc2723afab5" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.90750;volume_atom=42.5547;spin_atom=0.0;} // Ca:PAW_GGA:10Feb1998
  if(AUID=="d85b55fde567529b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.24285;volume_atom=37.3248;spin_atom=0.0;} // Ca:LDA:01Apr2000
  if(AUID=="d0c3b753999f51e0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.91163;volume_atom=41.5636;spin_atom=0.0;} // Ca_pv:PAW_GGA:05May1998
  if(AUID=="e65bce3da30c5d46" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.18236;volume_atom=38.2930;spin_atom=0.0;} // Ca:LDA:01Apr2000
  if(AUID=="26e7ac00e1d3cd41" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.21057;volume_atom=37.7193;spin_atom=0.0;} // Ca_pv:PAW_LDA:05May1998
  if(AUID=="b216760dca178acb" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.97638;volume_atom=41.7973;spin_atom=0.0;} // Ca_pv:PAW_PBE:06Sep2000
  if(AUID=="560b47e7f6be795e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.91723;volume_atom=41.8249;spin_atom=0.0;} // Ca_sv:PAW_GGA:04May1998
  if(AUID=="a5df1335919961ad" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.00107;volume_atom=42.1521;spin_atom=0.0;} // Ca_sv:PAW_PBE:06Sep2000
  if(AUID=="3bdd926d04dbca1e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.18326;volume_atom=37.8791;spin_atom=0.0;} // Ca_sv:PAW_LDA:17Apr2000
  if(AUID=="818406d11e64bddc" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-9.66266;volume_atom=42.9481;spin_atom=0.0;} // Ca_sv:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="529aebe19b9289ea" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-9.66372;volume_atom=43.0071;spin_atom=0.0;} // Ca_pv:PAW_PBE_KIN:SCAN:06Sep2000

  // Cd
  if(AUID=="0b0758f9de1eb816" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.761895;volume_atom=22.4448;spin_atom=0.0;} // Cd:PAW_GGA:04May1998
  if(AUID=="3388a058d6fd344d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.906233;volume_atom=22.4532;spin_atom=0.0;} // Cd:PAW_PBE:06Sep2000
  if(AUID=="c5c179cf112468e6" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.503460;volume_atom=19.7986;spin_atom=0.0;} // Cd:PAW_LDA:03Mar1998
  if(AUID=="07875e5a0789cd01" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.530730;volume_atom=20.0760;spin_atom=0.0;} // Cd:LDA:01Apr2000
  if(AUID=="f0b4c3233e3b19d8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.778952;volume_atom=22.8753;spin_atom=0.0;} // Cd:GGA:01Apr2000
  if(AUID=="fdf368e153b383d1" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-30.40360;volume_atom=21.1426;spin_atom=0.0;} // Cd:PAW_PBE_KIN:SCAN:06Sep2000

  // Cl
  if(AUID=="6bf17162620b7ce3" && nKIN) {found=TRUE;groundstate_structure="A11";groundstate_energy=-1.81560;volume_atom=37.3299;spin_atom=0.0;} // Cl:PAW_PBE:17Jan2003

  // Co
  if(AUID=="41f6635e51e063f9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.13760;volume_atom=10.0906;spin_atom=1.51284;} // Co:LDA:01Apr2000
  if(AUID=="4fb418aa1d0f1607" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.05079;volume_atom=10.9889;spin_atom=1.62675;} // Co:GGA:01Apr2000
  if(AUID=="f3c42ce518194ac4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.99076;volume_atom=10.8013;spin_atom=1.57160;} // Co:PAW_GGA:03Mar1998
  if(AUID=="9c4653e42c9f870b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.10876;volume_atom=10.8414;spin_atom=1.60252;} // Co:PAW_PBE:06Sep2000
  if(AUID=="23d97019de916d58" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.11138;volume_atom=9.98379;spin_atom=1.49449;} // Co:PAW_LDA:03Mar1998
  if(AUID=="4cc3ed3cb4d06e66" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-17.2118;volume_atom=10.3970;spin_atom=1.72628;} // Co:PAW_PBE_KIN:SCAN:02Aug2007
  if(AUID=="6d93f7d0ea628829" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-17.2429;volume_atom=10.4421;spin_atom=1.72940;} // Co_pv:PAW_PBE_KIN:SCAN:23Apr2009
  if(AUID=="12d6a4385e456ff2" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-17.3028;volume_atom=10.2840;spin_atom=1.72626;} // Co_sv:PAW_PBE_KIN:SCAN:23Jul2007

  // Cr
  if(AUID=="e56aa7dff1571851" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.6156;volume_atom=10.8356;spin_atom=0.0;} // Cr_pv:PAW_LDA:07Sep2000
  if(AUID=="411ed623d1b5054d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.51285;volume_atom=11.3963;spin_atom=0.0;} // Cr:PAW_PBE:06Sep2000
  if(AUID=="1365c4634d5c2ea1" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.46687;volume_atom=11.3926;spin_atom=0.0;} // Cr:PAW_GGA:03Mar1998
  if(AUID=="f020a26862c0a532" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.62940;volume_atom=11.5175;spin_atom=0.0;} // Cr_pv:PAW_PBE:07Sep2000
  if(AUID=="7ab8ddc788e0ece9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.6036;volume_atom=10.7989;spin_atom=0.0;} // Cr:LDA:01Apr2000
  if(AUID=="235696d86e852bef" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.43370;volume_atom=11.5367;spin_atom=0.0;} // Cr:GGA:01Apr2000
  if(AUID=="caa6b9695b37b7c7" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.46956;volume_atom=11.5096;spin_atom=0.0;} // Cr_pv:PAW_GGA:07Sep2000
  if(AUID=="46e25663256eec1e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.6345;volume_atom=10.7192;spin_atom=0.0;} // Cr:PAW_LDA:03Mar1998
  if(AUID=="46684a319ca854a2" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.2521;volume_atom=11.1136;spin_atom=0.0;} // Cr:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="c349e20347a8f61f" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.3863;volume_atom=11.0185;spin_atom=0.0;} // Cr_sv:PAW_PBE_KIN:SCAN:23Jul2007
  if(AUID=="d84af42941e8d51d" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.3160;volume_atom=11.2542;spin_atom=0.0;} // Cr_pv:PAW_PBE_KIN:SCAN:02Aug2007

  // Cs
  if(AUID=="e33c8a4092f5c27c" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-2.52656;volume_atom=114.954;spin_atom=0.0;} // Cs_sv_GW:PAW_PBE_KIN:23Mar2010

  // Cu
  if(AUID=="1794dd824ade7bfb" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.61717;volume_atom=10.7855;spin_atom=0.0;} // Cu_pv:PAW_LDA:19Apr2000
  if(AUID=="cd77b2d6e432ab5d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.72886;volume_atom=11.9621;spin_atom=0.0;} // Cu:PAW_GGA:05Jan2001
  if(AUID=="82741f3fd74a85b4" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.76087;volume_atom=12.0520;spin_atom=0.0;} // Cu:GGA:01Apr2000
  if(AUID=="62263cf2c885649c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.09701;volume_atom=11.8086;spin_atom=0.0;} // Cu_pv:PAW_PBE:06Sep2000
  if(AUID=="c7a1b91d209a648b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.71935;volume_atom=11.9654;spin_atom=0.0;} // Cu:PAW_PBE:05Jan2001
  if(AUID=="920126d2172870b5" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.70071;volume_atom=10.8985;spin_atom=0.0;} // Cu:PAW_LDA:03Mar1998
  if(AUID=="88d5d10238695c60" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.75059;volume_atom=10.9687;spin_atom=0.0;} // Cu:LDA:01Apr2000
  if(AUID=="a81d4798925d371d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.63758;volume_atom=11.7904;spin_atom=0.0;} // Cu_pv:PAW_GGA:19Apr2000
  if(AUID=="2ce0157f813571ea" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-15.1713;volume_atom=11.1286;spin_atom=0.0;} // Cu_pv:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="34032e61f91e67dc" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-15.1334;volume_atom=11.2002;spin_atom=0.0;} // Cu:PAW_PBE_KIN:SCAN:22Jun2005

  // Dy
  if(AUID=="0373199190a8abb7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.58813;volume_atom=31.8570;spin_atom=0.0;} // Dy_3:PAW_PBE:06Sep2000
  if(AUID=="57f7009ab6dbe961" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.52172;volume_atom=31.5074;spin_atom=0.0;} // Dy_3:PAW_GGA:10May2000
  if(AUID=="b7037854f4315aa3" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-47.8191;volume_atom=30.5976;spin_atom=0.0;} // Dy_3:PAW_PBE_KIN:SCAN:06Sep2000

  // Fe
  if(AUID=="53451df4122ca435" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.27793;volume_atom=10.5121;spin_atom=2.03940;} // Fe:LDA:01Apr2000
  if(AUID=="8a9cdad8b2fd617c" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.16576;volume_atom=11.2687;spin_atom=2.15865;} // Fe:PAW_GGA:03Mar1998
  if(AUID=="ca7f0085c55a7aa1" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.30490;volume_atom=11.6137;spin_atom=2.31462;} // Fe:GGA:01Apr2000
  if(AUID=="f5bec9ed423a0ff3" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.31138;volume_atom=11.3320;spin_atom=2.19316;} // Fe:PAW_PBE:06Sep2000
  if(AUID=="24dd99c81a522c7f" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.23151;volume_atom=10.3385;spin_atom=1.95457;} // Fe:PAW_LDA:03Mar1998
  if(AUID=="8c4dfe9b38461b3c" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.11742;volume_atom=11.2608;spin_atom=2.15340;} // Fe_pv:PAW_GGA:06May1998
  if(AUID=="15e65f5e50047695" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.45502;volume_atom=11.3277;spin_atom=2.19491;} // Fe_pv:PAW_PBE:06Sep2000
  if(AUID=="ea29acbd64d22bd0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.17964;volume_atom=10.3583;spin_atom=1.95509;} // Fe_pv:PAW_LDA:03Mar1998
  if(AUID=="46955a9bf291f348" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.27136;volume_atom=11.2061;spin_atom=2.14058;} // Fe_sv:PAW_GGA:14Sep2000
  if(AUID=="4d78633077b48d89" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.0562;volume_atom=11.5112;spin_atom=2.62602;} // Fe:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="c3ba95aef6c439d4" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.1013;volume_atom=11.6239;spin_atom=2.66996;} // Fe_pv:PAW_PBE_KIN:SCAN:02Aug2007
  if(AUID=="dbda82e42dec4a56" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-18.1874;volume_atom=11.3325;spin_atom=2.62763;} // Fe_sv:PAW_PBE_KIN:SCAN:23Jul2007

  // Ga
  if(AUID=="567d70c3d85e2c7f" && nKIN) {found=TRUE;groundstate_structure="A11";groundstate_energy=-2.88121;volume_atom=18.6097;spin_atom=0.0;} // Ga_h:PAW_PBE:09Apr2002

  // Ge
  if(AUID=="e58c5d01a1b37232" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.49261;volume_atom=24.1781;spin_atom=0.0;} // Ge:PAW_PBE:05Jan2001
  if(AUID=="e8adcacdfca69c78" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.17465;volume_atom=22.4994;spin_atom=0.0;} // Ge:PAW_LDA:03Mar1998
  if(AUID=="fd03d8c50dd891a8" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.52888;volume_atom=23.8391;spin_atom=0.0;} // Ge:GGA:01Apr2000
  if(AUID=="23df10d879b3b58b" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.50372;volume_atom=23.6953;spin_atom=0.0;} // Ge_h:PAW_PBE:09Apr2002
  if(AUID=="1a9537cf1728245b" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.62213;volume_atom=23.7855;spin_atom=0.0;} // Ge_d:PAW_PBE:06Sep2000
  if(AUID=="60aa6eae51da5ce4" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.19244;volume_atom=22.2190;spin_atom=0.0;} // Ge:LDA:01Apr2000
  if(AUID=="8f36ae9587237414" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.52998;volume_atom=23.6635;spin_atom=0.0;} // Ge_h:PAW_RPBE:09Apr2002
  if(AUID=="4538ff07d45d4093" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.13954;volume_atom=22.1566;spin_atom=0.0;} // Ge_d:PAW_LDA:03Mar1998
  if(AUID=="584c68479b3b04b1" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-4.48878;volume_atom=23.7431;spin_atom=0.0;} // Ge_d:PAW_GGA:03Mar1998
  if(AUID=="efd53a76a6f5032e" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.18648;volume_atom=22.0633;spin_atom=0.0;} // Ge_h:PAW_LDA:21Jan2003
  if(AUID=="18fa5cbdea403800" && SCAN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-19.0570;volume_atom=22.4931;spin_atom=0.0;} // Ge_d:PAW_PBE_KIN:SCAN:03Jul2007
  if(AUID=="571d96beea462000" && SCAN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-19.0597;volume_atom=22.5845;spin_atom=0.0;} // Ge_h:PAW_PBE_KIN:SCAN:09Apr2002
  if(AUID=="39a1162598519002" && SCAN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-19.0406;volume_atom=22.9647;spin_atom=0.0;} // Ge:PAW_PBE_KIN:SCAN:05Jan2001

  // He
  if(AUID=="c135275f4a0c77c9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=0.481074;volume_atom=5.48651;spin_atom=0.0;} // He:PAW_GGA:05Jan2001
  if(AUID=="60575574eefd1ae1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-0.04401;volume_atom=9.55685;spin_atom=0.0;} // He:PAW_LDA:07Sep2000
  if(AUID=="d3ff93e340e01503" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=0.236864;volume_atom=6.96149;spin_atom=0.0;} // He:PAW_PBE:05Jan2001

  // Hf
  if(AUID=="bb7552942568c756" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.88402;volume_atom=21.9665;spin_atom=0.0;} // Hf:GGA:01Apr2000
  if(AUID=="1af9281ed7be5121" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.7693;volume_atom=20.3825;spin_atom=0.0;} // Hf:LDA:01Apr2000
  if(AUID=="1703a3c3d5bdc240" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.87902;volume_atom=22.3346;spin_atom=0.0;} // Hf:PAW_GGA:6May2002
  if(AUID=="2a020b71ae180392" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.95711;volume_atom=22.2563;spin_atom=0.0;} // Hf:PAW_PBE:20Jan2003
  if(AUID=="5e77afe83e49db7f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.7507;volume_atom=20.7174;spin_atom=0.0;} // Hf:PAW_LDA:21Jan2003
  if(AUID=="d56380d31571155d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.83075;volume_atom=22.3351;spin_atom=0.0;} // Hf_pv:PAW_GGA:17Apr2000
  if(AUID=="95ae304e242d499a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.6877;volume_atom=20.7631;spin_atom=0.0;} // Hf_pv:PAW_LDA:17Apr2000
  if(AUID=="15832a2c336c16e0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.95294;volume_atom=22.4349;spin_atom=0.0;} // Hf_pv:PAW_PBE:06Sep2000
  if(AUID=="b5f0d256277acd2b" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-62.7979;volume_atom=21.3757;spin_atom=0.0;} // Hf:PAW_PBE_KIN:SCAN:20Jan2003
  if(AUID=="d78f974a4747b7ff" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-65.4899;volume_atom=21.3256;spin_atom=0.0;} // Hf_sv:PAW_PBE_KIN:SCAN:10Jan2008
  if(AUID=="8be591141a9fad1a" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-62.7249;volume_atom=21.4698;spin_atom=0.0;} // Hf_pv:PAW_PBE_KIN:SCAN:06Sep2000

  // Hg - requires rechecking
  if(AUID=="85cf414c5473e942" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.917393;volume_atom=22.3518;spin_atom=0.0;} // Hg:PAW_LDA:04Feb1998
  if(AUID=="bc158f9698c492b9" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.232554;volume_atom=29.8458;spin_atom=0.0;} // Hg:GGA:01Apr2000
  if(AUID=="96247dcf8cda9590" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.300828;volume_atom=29.3560;spin_atom=0.0;} // Hg:PAW_PBE:06Sep2000
  if(AUID=="df9134b9ac23884a" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.197681;volume_atom=28.7964;spin_atom=0.0;} // Hg:PAW_GGA:04May1998
  if(AUID=="03e3289690e1f860" && nKIN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-0.958689;volume_atom=22.7055;spin_atom=0.0;} // Hg:LDA:01Apr2000
  if(AUID=="56eba2b878545a76" && SCAN) {found=TRUE;groundstate_structure="A10(A3)";groundstate_energy=-68.11140;volume_atom=26.8710;spin_atom=0.0;} // Hg:PAW_PBE_KIN:SCAN:06Sep2000

  // Ho
  if(AUID=="ac7fc7278038dec4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.56866;volume_atom=31.3589;spin_atom=0.0;} // Ho_3:PAW_PBE:06Sep2000
  if(AUID=="20bea66adf91e4a2" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.49473;volume_atom=31.0514;spin_atom=0.0;} // Ho_3:PAW_GGA:10May2000
  if(AUID=="141af062a506aee7" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-49.0018;volume_atom=30.1155;spin_atom=0.0;} // Ho_3:PAW_PBE_KIN:SCAN:06Sep2000

  // K
  if(AUID=="1faff02cfb784b29" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.02995;volume_atom=73.8776;spin_atom=0.0;} // K:GGA:01Apr2000
  if(AUID=="bc9d6fb83de4a3e6" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.17040;volume_atom=64.3539;spin_atom=0.0;} // K:LDA:01Apr2000
  if(AUID=="e9c22b07bc0ec1ed" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.02890;volume_atom=73.4124;spin_atom=0.0;} // K:GGA:01Apr2000
  if(AUID=="90f35bb229eef165" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.02436;volume_atom=73.3337;spin_atom=0.0;} // K_pv:PAW_GGA:11Feb1998
  if(AUID=="828b9d8df4c3212b" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.16884;volume_atom=63.8721;spin_atom=0.0;} // K:LDA:01Apr2000
  if(AUID=="57bec8a3e025c3e4" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.16059;volume_atom=63.6556;spin_atom=0.0;} // K_pv:PAW_LDA:02Feb1998
  if(AUID=="941fbf66b21efd05" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.02692;volume_atom=72.4893;spin_atom=0.0;} // K_pv:PAW_PBE:17Jan2003
  if(AUID=="5c717bdf563a1964" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.03858;volume_atom=72.8681;spin_atom=0.0;} // K_sv:PAW_GGA:04May1998
  if(AUID=="cc26596520b9b525" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.09647;volume_atom=73.1151;spin_atom=0.0;} // K_sv:PAW_PBE:06Sep2000
  if(AUID=="85c370c6dd5f769e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-1.15723;volume_atom=63.4202;spin_atom=0.0;} // K_sv:PAW_LDA:24Mar1998
  if(AUID=="ef04627355867a05" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.29912;volume_atom=71.3848;spin_atom=0.0;} // K_pv:PAW_PBE_KIN:SCAN:17Jan2003
  if(AUID=="3ac2263c1999fd5c" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.33540;volume_atom=74.0621;spin_atom=0.0;} // K_sv:PAW_PBE_KIN:SCAN:06Sep2000

  // In
  if(AUID=="5ce1eee07a5df3a7" && nKIN) {found=TRUE;groundstate_structure="A6";groundstate_energy=-2.72115;volume_atom=27.1064;spin_atom=0.0;} // In_d:PAW_PBE:06Sep2000
  if(AUID=="b26165280c5a6d7c" && SCAN) {found=TRUE;groundstate_structure="A6";groundstate_energy=-33.3068;volume_atom=26.1419;spin_atom=0.0;} // In:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="ff28085cf0ed0c70" && SCAN) {found=TRUE;groundstate_structure="A6";groundstate_energy=-33.3094;volume_atom=25.5954;spin_atom=0.0;} // In_d:PAW_PBE_KIN:SCAN:06Sep2000

  // Ir
  if(AUID=="78171971d3315e28" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-10.3309;volume_atom=13.8873;spin_atom=0.0;} // Ir:PAW_LDA:10Feb1998
  if(AUID=="8714595dd67f5cc6" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.72200;volume_atom=14.6380;spin_atom=0.0;} // Ir:GGA:01Apr2000
  if(AUID=="0d6f43382108d13d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-10.2667;volume_atom=13.9078;spin_atom=0.0;} // Ir:LDA:01Apr2000
  if(AUID=="b33a03214533a06a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.79336;volume_atom=14.5807;spin_atom=0.0;} // Ir:PAW_GGA:04May1998
  if(AUID=="9c41e654d31db16a" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.85711;volume_atom=14.5260;spin_atom=0.0;} // Ir:PAW_PBE:06Sep2000
  if(AUID=="be8b84b10c2543f2" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-71.2587;volume_atom=13.5393;spin_atom=0.0;} // Ir:PAW_PBE_KIN:SCAN:06Sep2000

  // La
  if(AUID=="d0888fff3f31a114" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.56862;volume_atom=31.6652;spin_atom=0.0;} // La:PAW_LDA:17Apr2000
  if(AUID=="ec88cfb1842bf27d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.90610;volume_atom=36.2172;spin_atom=0.0;} // La:PAW_GGA:14Apr2000
  if(AUID=="8353806ee117d7fd" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.91782;volume_atom=36.5244;spin_atom=0.0;} // La:PAW_PBE:06Sep2000
  if(AUID=="d04dad60260af779" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.86223;volume_atom=37.0366;spin_atom=0.0;} // La_s:PAW_PBE:06Sep2000
  if(AUID=="b2abe18f66fcdafe" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.86314;volume_atom=36.6893;spin_atom=0.0;} // La_s:PAW_GGA:17Apr2000
  if(AUID=="0549ee555a0b18c8" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.53774;volume_atom=32.3324;spin_atom=0.0;} // La_s:PAW_LDA:17Apr2000
  if(AUID=="450ae928112e37cc" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-43.0405;volume_atom=37.4721;spin_atom=0.0;} // La_s:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="4098c299a18eef08" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-43.1348;volume_atom=36.4040;spin_atom=0.0;} // La:PAW_PBE_KIN:SCAN:06Sep2000

  // Li
  if(AUID=="d51fe7240453d7c4" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.897390;volume_atom=20.3189;spin_atom=0.0;} // Li:PAW_PBE:17Jan2003
  if(AUID=="67021bad74ddb4c5" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.897720;volume_atom=20.3318;spin_atom=0.0;} // Li_sv:PAW_GGA:23Jan2001
  if(AUID=="99a2f30dac0e41d6" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.905030;volume_atom=20.2771;spin_atom=0.0;} // Li_sv:PAW_PBE:23Jan2001
  if(AUID=="e3ad72d8dbbd284c" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-2.04188;volume_atom=19.0006;spin_atom=0.0;} // Li:PAW_LDA:21Jan2003
  if(AUID=="1c81ce7af52fff4e" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-2.04908;volume_atom=18.8562;spin_atom=0.0;} // Li:LDA:01Apr2000
  if(AUID=="6076711c64467c4f" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.90718;volume_atom=19.7515;spin_atom=0.0;} // Li:GGA:01Apr2000
  if(AUID=="c680d72aeb02882e" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-2.04290;volume_atom=19.0034;spin_atom=0.0;} // Li_sv:PAW_LDA:19Jan2001
  if(AUID=="1c81ce7af52fff4e" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-2.04908;volume_atom=18.8562;spin_atom=0.0;} // Li:LDA:01Apr2000
  if(AUID=="0666c5aca67cb28d" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.89268;volume_atom=20.3860;spin_atom=0.0;} // Li:PAW_GGA:21Jan2003
  if(AUID=="6076711c64467c4f" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-1.90718;volume_atom=19.7515;spin_atom=0.0;} // Li:GGA:01Apr2000
  if(AUID=="84026814dcd62d7a" && SCAN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-2.33938;volume_atom=20.8216;spin_atom=0.0;} // Li:PAW_PBE_KIN:SCAN:17Jan2003

  // Mg
  if(AUID=="285c78eec2b06a80" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.51618;volume_atom=22.8934;spin_atom=0.0;} // Mg:GGA:01Apr2000
  if(AUID=="d8bb30571ef22203" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.52082;volume_atom=22.7555;spin_atom=0.0;} // Mg:GGA:01Apr2000
  if(AUID=="78efd19d72d98cb8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.78695;volume_atom=21.4123;spin_atom=0.0;} // Mg:LDA:01Apr2000
  if(AUID=="61fec4224ab4cab8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.78165;volume_atom=21.5278;spin_atom=0.0;} // Mg:LDA:01Apr2000
  if(AUID=="de98b20ad5c678cf" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.52338;volume_atom=22.8721;spin_atom=0.0;} // Mg:PAW_GGA:05Jan2001
  if(AUID=="7d88d84366cb7141" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.54148;volume_atom=22.8463;spin_atom=0.0;} // Mg:PAW_PBE:05Jan2001
  if(AUID=="b3d0cc0187ea8112" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.78256;volume_atom=21.5714;spin_atom=0.0;} // Mg:PAW_LDA:02Mar1998
  if(AUID=="7a5a79e9587d1be6" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.48173;volume_atom=22.8793;spin_atom=0.0;} // Mg:GGA:01Apr2000
  if(AUID=="1308594682e50447" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.78486;volume_atom=21.5561;spin_atom=0.0;} // Mg:LDA:01Apr2000
  if(AUID=="3e32a29daafd48c5" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.47720;volume_atom=22.7969;spin_atom=0.0;} // Mg_pv:PAW_GGA:10Feb1998
  if(AUID=="c2b2c35b8213a287" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.75631;volume_atom=21.4740;spin_atom=0.0;} // Mg_pv:PAW_LDA:02Mar1998
  if(AUID=="393ee52ce57367b8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.59348;volume_atom=22.7983;spin_atom=0.0;} // Mg_pv:PAW_PBE:06Sep2000
  if(AUID=="2c5dfb19cf1554e5" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.96828;volume_atom=22.6560;spin_atom=0.0;} // Mg:PAW_PBE_KIN:SCAN:13Apr2007
  if(AUID=="26735ac174693885" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-4.98465;volume_atom=22.7620;spin_atom=0.0;} // Mg_pv:PAW_PBE_KIN:SCAN:13Apr2007

  // Mn
  if(AUID=="8f2e213b3592f62e" && nKIN) {found=TRUE;groundstate_structure="A12";groundstate_energy=-9.02786;volume_atom=10.7275;spin_atom=0.0;} // Mn:PAW_PBE:06Sep2000
  if(AUID=="99be850476e2dfb3" && nKIN) {found=TRUE;groundstate_structure="A12";groundstate_energy=-9.15350;volume_atom=10.6973;spin_atom=0.0;} // Mn_pv:PAW_PBE:07Sep2000

  // Mo
  if(AUID=="8f74857b7edbc55b" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.8344;volume_atom=15.6219;spin_atom=0.0;} // Mo:GGA:01Apr2000
  if(AUID=="4690ce5c08300f84" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.9104;volume_atom=15.6337;spin_atom=0.0;} // Mo:PAW_GGA:08Jan2002
  if(AUID=="c56388dcfe8555c0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.1515;volume_atom=14.8749;spin_atom=0.0;} // Mo:LDA:01Apr2000
  if(AUID=="c90275c91471b466" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.9465;volume_atom=15.5844;spin_atom=0.0;} // Mo:PAW_PBE:08Apr2002
  if(AUID=="3e6c7aec9ade11fe" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.5202;volume_atom=16.1255;spin_atom=0.0;} // Mo:GGA:01Apr2000
  if(AUID=="097174f08d663402" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.7663;volume_atom=15.3630;spin_atom=0.0;} // Mo:LDA:01Apr2000
  if(AUID=="173298095bee51b0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.8042;volume_atom=15.9061;spin_atom=0.0;} // Mo_pv:PAW_GGA:08Jan2002
  if(AUID=="517b9d282b85b101" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.0734;volume_atom=15.1997;spin_atom=0.0;} // Mo_pv:PAW_LDA:08Jan2002
  if(AUID=="444a5f723c4c0576" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.8439;volume_atom=15.8555;spin_atom=0.0;} // Mo_pv:PAW_PBE:08Apr2002
  if(AUID=="027c5a9d41c6523f" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.1895;volume_atom=15.0334;spin_atom=0.0;} // Mo_sv:PAW_LDA:15Nov2001
  if(AUID=="5001e3854b3663c0" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-34.8382;volume_atom=15.5390;spin_atom=0.0;} // Mo_pv:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="0808a3790788f802" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-34.7737;volume_atom=15.3662;spin_atom=0.0;} // Mo:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="eb2ad6b5ae2c88f9" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-34.8748;volume_atom=15.5342;spin_atom=0.0;} // Mo_sv:PAW_PBE_KIN:SCAN:02Feb2006

  // Na
  if(AUID=="bd75e1ab544c0ed0" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.31096;volume_atom=36.2541;spin_atom=0.0;} // Na_pv:PAW_PBE:05Jan2001
  if(AUID=="192d4f6f863806ad" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.31443;volume_atom=35.0106;spin_atom=0.0;} // Na_sv:PAW_GGA:28Sep2000
  if(AUID=="637322e42c89bd5e" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.31537;volume_atom=35.6060;spin_atom=0.0;} // Na_sv:PAW_PBE:28Sep2000
  if(AUID=="2215677fd2f98aea" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.45287;volume_atom=33.2473;spin_atom=0.0;} // Na:PAW_LDA:24Mar1998
  if(AUID=="e7d71dba59886669" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.30398;volume_atom=37.1851;spin_atom=0.0;} // Na:PAW_GGA:05Jan2001
  if(AUID=="79f12c5d16cf98b8" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.30640;volume_atom=36.7352;spin_atom=0.0;} // Na:PAW_PBE:08Apr2002
  if(AUID=="f44fe5066a3bc573" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.30382;volume_atom=36.4114;spin_atom=0.0;} // Na_pv:PAW_GGA:05Jan2001
  if(AUID=="189d870b4537263b" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.45445;volume_atom=31.4092;spin_atom=0.0;} // Na_sv:PAW_LDA:28Sep2000
  if(AUID=="b994b2842f3b80ab" && nKIN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-1.45431;volume_atom=32.6970;spin_atom=0.0;} // Na_pv:PAW_LDA:28Sep2000
  if(AUID=="86f3cf1cb9c67ee6" && SCAN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-4.20653;volume_atom=37.2429;spin_atom=0.0;} // Na_pv:PAW_PBE_KIN:SCAN:19Sep2006
  if(AUID=="63e83f9cf84ac969" && SCAN) {found=TRUE;groundstate_structure="A2/A7";groundstate_energy=-4.18698;volume_atom=36.5426;spin_atom=0.0;} // Na:PAW_PBE_KIN:SCAN:08Apr2002

  // Nb
  if(AUID=="31fbbfcce7b1cf49" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.3250;volume_atom=16.8059;spin_atom=0.0;} // Nb:LDA:01Apr2000
  if(AUID=="724f3dcc974f7f25" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.1730;volume_atom=17.8196;spin_atom=0.0;} // Nb:GGA:01Apr2000
  if(AUID=="e1ae346e0f041fdc" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.99099;volume_atom=18.3604;spin_atom=0.0;} // Nb:GGA:01Apr2000
  if(AUID=="31ac6525497e7fb5" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.0804;volume_atom=17.3804;spin_atom=0.0;} // Nb:LDA:01Apr2000
  if(AUID=="3df5731180246390" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.0550;volume_atom=18.2510;spin_atom=0.0;} // Nb_pv:PAW_GGA:09Jan2002
  if(AUID=="73b4ad5a036338fa" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.1463;volume_atom=17.3220;spin_atom=0.0;} // Nb_pv:PAW_LDA:09Jan2002
  if(AUID=="d9711dcad08e9fe9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.0830;volume_atom=18.2314;spin_atom=0.0;} // Nb_pv:PAW_PBE:08Apr2002
  if(AUID=="79b1b1f832ef1be0" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.1943;volume_atom=18.0754;spin_atom=0.0;} // Nb_sv:PAW_GGA:14Nov2001
  if(AUID=="141a7bea5919e0b3" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-10.2253;volume_atom=18.0525;spin_atom=0.0;} // Nb_sv:PAW_PBE:17Jan2003
  if(AUID=="9fdcd7998786bd34" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.2626;volume_atom=17.0887;spin_atom=0.0;} // Nb_sv:PAW_LDA:15Nov2001
  if(AUID=="2d2f42f890bd3354" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-33.0916;volume_atom=18.0829;spin_atom=0.0;} // Nb_pv:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="b383f09c3645c1a5" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-33.2389;volume_atom=17.8122;spin_atom=0.0;} // Nb_sv:PAW_PBE_KIN:SCAN:25May2007

  // Ne
  if(AUID=="bf15cfac49144396" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.048593;volume_atom=13.7731;spin_atom=0.0;} // Ne:PAW_LDA:07Sep2000
  if(AUID=="78bf25649465a255" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.103775;volume_atom=13.5480;spin_atom=0.0;} // Ne:LDA:01Apr2000
  if(AUID=="4bda2c3e931df0db" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.102230;volume_atom=13.5480;spin_atom=0.0;} // Ne:LDA:01Apr2000
  if(AUID=="445c00adf152283d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.066215;volume_atom=20.2857;spin_atom=0.0;} // Ne:PAW_GGA:05Jan2001
  if(AUID=="91b32db4af57705e" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-0.032434;volume_atom=20.1045;spin_atom=0.0;} // Ne:PAW_PBE:05Jan2001
  if(AUID=="326c9a1ee646c679" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-2.394710;volume_atom=15.4322;spin_atom=0.0;} // Ne:PAW_PBE_KIN:SCAN:05Jan2001

  // Ni
  if(AUID=="eafb70a88aeb4012" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.47655;volume_atom=10.9943;spin_atom=0.636833;} // Ni:GGA:01Apr2000
  if(AUID=="1dc196d6de4f6b37" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.57108;volume_atom=10.9016;spin_atom=0.627335;} // Ni:PAW_PBE:06Sep2000
  if(AUID=="5737f6f3ba793dbe" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.56769;volume_atom=10.0249;spin_atom=0.578133;} // Ni:PAW_LDA:03Mar1998
  if(AUID=="ead949d4e6fa93d1" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.77783;volume_atom=10.7823;spin_atom=0.633813;} // Ni_pv:PAW_PBE:06Sep2000
  if(AUID=="134ab84006e2f008" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.58610;volume_atom=10.1132;spin_atom=0.586344;} // Ni:LDA:01Apr2000
  if(AUID=="f316c1e20cc629f5" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.46133;volume_atom=10.8760;spin_atom=0.608677;} // Ni:PAW_GGA:03Mar1998
  if(AUID=="55fd1c17c8c96b19" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.39117;volume_atom=10.7667;spin_atom=0.614692;} // Ni_pv:PAW_GGA:19Apr2000
  if(AUID=="bee457b2befa4114" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.50164;volume_atom=9.95066;spin_atom=0.583619;} // Ni_pv:PAW_LDA:19Apr2000
  if(AUID=="1f50d23b4fbbe9d3" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-16.0643;volume_atom=10.3563;spin_atom=0.696317;} // Ni:PAW_PBE_KIN:SCAN:02Aug2007
  if(AUID=="1141cd9c7118591b" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-16.1373;volume_atom=10.2771;spin_atom=0.738419;} // Ni_pv:PAW_PBE_KIN:SCAN:06Sep2000

  // O
  if(AUID=="85ded82734544fa9" && nKIN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-4.50622;volume_atom=12.5315;spin_atom=0.0;} // O:PAW_PBE:08Apr2002

  // Os
  if(AUID=="4cdd79622cbed39c" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.1466;volume_atom=14.3178;spin_atom=0.0;} // Os:PAW_GGA:06Feb2003
  if(AUID=="458230748f9ee8c4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.1170;volume_atom=14.3175;spin_atom=0.0;} // Os:GGA:01Apr2000
  if(AUID=="8c0f432aa073c44a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.6831;volume_atom=13.6759;spin_atom=0.0;} // Os:LDA:01Apr2000
  if(AUID=="b1954be3210262f1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.2440;volume_atom=14.2137;spin_atom=0.0;} // Os:PAW_PBE:17Jan2003
  if(AUID=="d38b52b3db36dc96" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.7031;volume_atom=13.7273;spin_atom=0.0;} // Os:PAW_LDA:21Jan2003
  if(AUID=="216ecda62b9b32bf" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.2193;volume_atom=14.2979;spin_atom=0.0;} // Os_pv:PAW_PBE:20Jan2003
  if(AUID=="d2665c78be0bcda4" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.1086;volume_atom=14.3251;spin_atom=0.0;} // Os_pv:PAW_GGA:10Feb1998
  if(AUID=="6948d5b374665a4d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.6613;volume_atom=13.7168;spin_atom=0.0;} // Os_pv:PAW_LDA:10Feb1998
  if(AUID=="0b4b416ff63ebb16" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-71.8359;volume_atom=13.6516;spin_atom=0.0;} // Os:PAW_PBE_KIN:SCAN:17Jan2003
  if(AUID=="6fba18d27313ce5e" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-71.7960;volume_atom=13.7403;spin_atom=0.0;} // Os_pv:PAW_PBE_KIN:SCAN:20Jan2003

  // P
  if(AUID=="0df9932192e3546b" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-5.32414;volume_atom=16.0927;spin_atom=0.0;} // P:PAW_PBE:17Jan2003
  if(AUID=="d64d49a8028ca2bc" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-10.3005;volume_atom=16.0865;spin_atom=0.0;} // P_h:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="3cafdbeb02f9e951" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-10.3006;volume_atom=16.0927;spin_atom=0.0;} // P:PAW_PBE_KIN:SCAN:06Sep2000

  // Pb
  if(AUID=="1c02af03a0cb18b2" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.59467;volume_atom=31.8712;spin_atom=0.0;} // Pb_d:GGA:01Apr2000
  if(AUID=="58341640c6c638ed" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.26624;volume_atom=29.0840;spin_atom=0.0;} // Pb_d:LDA:01Apr2000
  if(AUID=="7d26482cab914656" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.54563;volume_atom=31.4146;spin_atom=0.0;} // Pb_d:PAW_GGA:04May1998
  if(AUID=="0e47cf73a95cbfc2" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.70470;volume_atom=31.5298;spin_atom=0.0;} // Pb_d:PAW_PBE:06Sep2000
  if(AUID=="1010cb02b8b71033" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.20301;volume_atom=28.7142;spin_atom=0.0;} // Pb_d:PAW_LDA:30Apr1998
  if(AUID=="2f64d0d7f1c21e2c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.57056;volume_atom=31.7950;spin_atom=0.0;} // Pb:PAW_PBE:08Apr2002
  if(AUID=="899670f2ebab4069" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.58825;volume_atom=31.8825;spin_atom=0.0;} // Pb:GGA:01Apr2000
  if(AUID=="9568f688bb71b2b0" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-3.60201;volume_atom=31.7526;spin_atom=0.0;} // Pb:PAW_GGA:25Jul2001
  if(AUID=="626de6bdde387007" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.25933;volume_atom=29.0453;spin_atom=0.0;} // Pb:LDA:01Apr2000
  if(AUID=="381b125c2959cfe6" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-4.28081;volume_atom=28.9939;spin_atom=0.0;} // Pb:PAW_LDA:02Sep2001
  if(AUID=="727d5248f9fb2e47" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-75.7142;volume_atom=30.6165;spin_atom=0.0;} // Pb:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="3143a75eeaa29ea5" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-75.7300;volume_atom=30.1590;spin_atom=0.0;} // Pb_d:PAW_PBE_KIN:SCAN:06Sep2000

  // Pd
  if(AUID=="854a0d2f9191bbb2" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.19801;volume_atom=15.4798;spin_atom=0.0;} // Pd:GGA:01Apr2000
  if(AUID=="7e8d64cc173f9f0c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.21249;volume_atom=15.3824;spin_atom=0.0;} // Pd:PAW_GGA:05Jan2001
  if(AUID=="552d7d192268351d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.42893;volume_atom=14.3040;spin_atom=0.0;}  // Pd:LDA:01Apr2000
  if(AUID=="3b6d075d41b9645f" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.44052;volume_atom=14.2369;spin_atom=0.0;} // Pd:PAW_LDA:09Oct1998
  if(AUID=="6ade515555bfea29" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.17847;volume_atom=15.3433;spin_atom=0.0;} // Pd:PAW_PBE:05Jan2001
  if(AUID=="46106d1a706186be" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.13856;volume_atom=15.4912;spin_atom=0.0;} // Pd_pv:PAW_GGA:04Mar1998
  if(AUID=="4af3f66c14a34cee" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-5.38197;volume_atom=15.3938;spin_atom=0.0;} // Pd_pv:PAW_PBE:06Sep2000
  if(AUID=="4158852ee0320965" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.35943;volume_atom=14.2548;spin_atom=0.0;} // Pd_pv:PAW_LDA:17Apr2000
  if(AUID=="2ef0a2bc50307dc4" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-32.8128;volume_atom=14.7536;spin_atom=0.0;} // Pd:PAW_PBE_KIN:SCAN:04Jan2005
  if(AUID=="35cb7e1d80725a17" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-32.7762;volume_atom=14.6457;spin_atom=0.0;} // Pd_pv:PAW_PBE_KIN:SCAN:28Jan2005

  // Pt
  if(AUID=="4934ee186867e17b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.01467;volume_atom=15.8490;spin_atom=0.0;} // Pt:GGA:01Apr2000
  if(AUID=="59081d68fabd053c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.05482;volume_atom=15.6527;spin_atom=0.0;} // Pt:PAW_PBE:05Jan2001
  if(AUID=="129e32c460592b3d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.42489;volume_atom=14.9080;spin_atom=0.0;} // Pt:LDA:01Apr2000
  if(AUID=="2f9c34be5eee2bc7" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.45079;volume_atom=14.8389;spin_atom=0.0;} // Pt:PAW_LDA:17Apr2000
  if(AUID=="05e3fa0ef96c01d7" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-6.04404;volume_atom=15.7566;spin_atom=0.0;} // Pt:PAW_GGA:05Jan2001
  if(AUID=="c67ecad4e9259793" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-70.2902;volume_atom=14.7627;spin_atom=0.0;} // Pt:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="35f9da8fbc175cc5" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-70.1284;volume_atom=14.8589;spin_atom=0.0;} // Pt_pv:PAW_PBE_KIN:SCAN:12Dec2005

  // Rb
  if(AUID=="4a3726ab01a44892" && nKIN) {found=TRUE;groundstate_structure="A2*";groundstate_energy=-0.962354;volume_atom=90.2922;spin_atom=0.0;} // Rb_sv:PAW_PBE:06Sep2000

  // Re
  if(AUID=="7c2099bc993d23a3" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.3054;volume_atom=14.9209;spin_atom=0.0;} // Re:GGA:01Apr2000
  if(AUID=="98510e8bfa14709e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.7774;volume_atom=14.2528;spin_atom=0.0;} // Re:LDA:01Apr2000
  if(AUID=="1e59564cdfe6dd0a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.3475;volume_atom=14.9608;spin_atom=0.0;} // Re:PAW_GGA:05Jan2001
  if(AUID=="e5efc178b4648996" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.4113;volume_atom=14.8653;spin_atom=0.0;} // Re:PAW_PBE:17Jan2003
  if(AUID=="7289c4501ec2c23a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.7534;volume_atom=14.3332;spin_atom=0.0;} // Re:PAW_LDA:21Jan2003
  if(AUID=="754454a117ea525b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.2208;volume_atom=15.0395;spin_atom=0.0;} // Re_pv:PAW_GGA:11Feb1998
  if(AUID=="598a3a038e3785b3" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.6537;volume_atom=14.4140;spin_atom=0.0;} // Re_pv:PAW_LDA:11Feb1998
  if(AUID=="d72276b27490b853" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-12.4325;volume_atom=14.9985;spin_atom=0.0;} // Re_pv:PAW_PBE:06Sep2000
  if(AUID=="6c3f2ec3842638b8" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-71.0938;volume_atom=14.2825;spin_atom=0.0;} // Re:PAW_PBE_KIN:SCAN:17Jan2003
  if(AUID=="3f9c2349a9115560" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-71.0339;volume_atom=14.3427;spin_atom=0.0;} // Re_pv:PAW_PBE_KIN:SCAN:06Sep2000

  // Rh
  if(AUID=="e9aa6eb55d4587f1" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.13631;volume_atom=14.1993;spin_atom=0.0;} // Rh_pv:PAW_GGA:17Apr2000
  if(AUID=="c6d7ff86b44b5d4c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.52654;volume_atom=13.3693;spin_atom=0.0;} // Rh_pv:PAW_LDA:17Apr2000
  if(AUID=="b1547f3221c0c660" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.34058;volume_atom=14.1854;spin_atom=0.0;} // Rh_pv:PAW_PBE:06Sep2000
  if(AUID=="883dfe9dc76b6f31" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.58863;volume_atom=13.3502;spin_atom=0.0;} // Rh:LDA:01Apr2000
  if(AUID=="b02edafce0da47bc" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.17866;volume_atom=14.2151;spin_atom=0.0;} // Rh:GGA:01Apr2000
  if(AUID=="e06760ebe9e0e518" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.27044;volume_atom=14.1233;spin_atom=0.0;} // Rh:PAW_PBE:06Sep2000
  if(AUID=="2297d0177b7db39d" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-7.22290;volume_atom=14.1495;spin_atom=0.0;} // Rh:PAW_GGA:04May1998
  if(AUID=="b03cdc726a115d0b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-8.62793;volume_atom=13.3185;spin_atom=0.0;} // Rh:PAW_LDA:03Mar1998
  if(AUID=="dc3fbc9801231714" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-33.8609;volume_atom=13.4868;spin_atom=0.0;} // Rh:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="98605c42de12449b" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-33.7489;volume_atom=13.4857;spin_atom=0.0;} // Rh_pv:PAW_PBE_KIN:SCAN:25Jan2005

  // Ru
  if(AUID=="e2e4aaaa46da8131" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.20414;volume_atom=13.8077;spin_atom=0.0;} // Ru:PAW_PBE:06Sep2000
  if(AUID=="1fa7a7266ef21653" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.09581;volume_atom=13.8865;spin_atom=0.0;} // Ru:GGA:01Apr2000
  if(AUID=="24117781e96ec111" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.08258;volume_atom=13.9181;spin_atom=0.0;} // Ru_pv:PAW_GGA:10Feb1998
  if(AUID=="d45ec3c5462bef36" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.5801;volume_atom=13.1589;spin_atom=0.0;} // Ru:LDA:01Apr2000
  if(AUID=="48750268bde7c56d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.16053;volume_atom=13.8456;spin_atom=0.0;} // Ru:PAW_GGA:03Mar1998
  if(AUID=="64d0c8f1c723d0dc" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.6345;volume_atom=13.1584;spin_atom=0.0;} // Ru:PAW_LDA:17Apr2000
  if(AUID=="04267bcea324606c" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.5430;volume_atom=13.2301;spin_atom=0.0;} // Ru_pv:PAW_LDA:03Mar1998
  if(AUID=="a3d039d8d88a4a86" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.27135;volume_atom=13.8821;spin_atom=0.0;} // Ru_pv:PAW_PBE:06Sep2000
  if(AUID=="c6fbedbf10339a1d" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.16320;volume_atom=13.8588;spin_atom=0.0;} // Ru_sv:PAW_GGA:02Oct2001
  if(AUID=="6a824616a7311021" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-34.8327;volume_atom=13.3205;spin_atom=0.0;} // Ru_pv:PAW_PBE_KIN:SCAN:28Jan2005
  if(AUID=="ee46e79f573bbba1" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-34.9124;volume_atom=13.0842;spin_atom=0.0;} // Ru:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="1550a87ae48844e0" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-34.8652;volume_atom=13.3100;spin_atom=0.0;} // Ru_sv:PAW_GGA_KIN:SCAN:28Jan2005

  // Sb
  if(AUID=="6ec70613ee3528b2" && nKIN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-3.88932;volume_atom=27.1685;spin_atom=0.0;} // Sb:PAW_PBE:06Sep2000
  if(AUID=="bbf81761207d2817" && SCAN) {found=TRUE;groundstate_structure="A7";groundstate_energy=-37.1553;volume_atom=28.6389;spin_atom=0.0;} // Sb:PAW_PBE_KIN:SCAN:06Sep2000

  // Sc
  if(AUID=="4e27e0eebb4ca9ea" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.69169;volume_atom=22.2639;spin_atom=0.0;} // Sc:LDA:01Apr2000
  if(AUID=="11f9d8ee4231562f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.13734;volume_atom=24.1795;spin_atom=0.0;} // Sc:GGA:01Apr2000
  if(AUID=="1d27da013e36ce9b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.18297;volume_atom=24.0170;spin_atom=0.0;} // Sc:PAW_GGA:08Aug2001
  if(AUID=="90cbeb0ff7fd69f7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.15503;volume_atom=24.8607;spin_atom=0.0;} // Sc:GGA:01Apr2000
  if(AUID=="5b645c3a2198a968" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.67171;volume_atom=22.9075;spin_atom=0.0;} // Sc:LDA:01Apr2000
  if(AUID=="a563cb0b31f21d93" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.22496;volume_atom=24.2228;spin_atom=0.0;} // Sc_sv:PAW_GGA:07Sep2000
  if(AUID=="4b48722b5ae8f9e9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.33212;volume_atom=24.4214;spin_atom=0.0;} // Sc_sv:PAW_PBE:07Sep2000
  if(AUID=="dd64f21bfc2e4cab" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.72459;volume_atom=22.3699;spin_atom=0.0;} // Sc_sv:PAW_LDA:07Sep2000
  if(AUID=="69bd66903db1e199" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-14.1593;volume_atom=24.2183;spin_atom=0.619075;} // Sc:PAW_PBE_KIN:SCAN:04Feb2005  MAGNETIC ??
  if(AUID=="83d4f13df3294a60" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-14.2457;volume_atom=24.7117;spin_atom=0.597575;} // Sc_sv:PAW_PBE_KIN:SCAN:07Sep2000  MAGNETIC ??

  // Se
  if(AUID=="40060a482797271e" && nKIN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-3.48266;volume_atom=29.6441;spin_atom=0.0;} // Se:PAW_PBE:06Sep2000
  if(AUID=="00143077dee7333f" && SCAN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-20.0853;volume_atom=28.3006;spin_atom=0.0;} // Se:PAW_PBE_KIN:SCAN:06Sep2000

  // Si
  if(AUID=="096bd93da38a71ba" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.46219;volume_atom=20.3474;spin_atom=0.0;} // Si_h:PAW_GGA:08Apr2002
  if(AUID=="b19b07eb7ec794be" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.43180;volume_atom=20.2984;spin_atom=0.0;} // Si:GGA:01Apr2000
  if(AUID=="5d164d220fd7ceee" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.44184;volume_atom=20.3640;spin_atom=0.0;} // Si_h:PAW_PBE:08Apr2002
  if(AUID=="69684f2b4007eeda" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.44713;volume_atom=20.2871;spin_atom=0.0;} // Si:GGA:01Apr2000
  if(AUID=="530c405087cb5e2b" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.99306;volume_atom=19.5385;spin_atom=0.0;} // Si:LDA:01Apr2000
  if(AUID=="940bc6e41f3ea1ba" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.96148;volume_atom=19.6993;spin_atom=0.0;} // Si_h:PAW_LDA:21Jan2003
  if(AUID=="6b54cf2aa9fd187c" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.97452;volume_atom=19.5596;spin_atom=0.0;} // Si:LDA:01Apr2000
  if(AUID=="2e636f3ea3dd411d" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.95931;volume_atom=19.7009;spin_atom=0.0;} // Si:PAW_LDA:02Apr1999
  if(AUID=="3ea77bd2a2af2ad4" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.42373;volume_atom=20.4310;spin_atom=0.0;} // Si:PAW_PBE:05Jan2001
  if(AUID=="d47b3f10236456be" && nKIN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-5.43030;volume_atom=20.4209;spin_atom=0.0;} // Si:PAW_GGA:05Jan2001
  if(AUID=="afc7580b135e98ce" && SCAN) {found=TRUE;groundstate_structure="A4";groundstate_energy=-10.0057;volume_atom=19.9704;spin_atom=0.0;} // Si:PAW_PBE_KIN:SCAN:05Jan2001

  // Sn
  if(AUID=="ec3f84d587dd2653" && nKIN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-3.96372;volume_atom=28.1306;spin_atom=0.0;} // Sn_d:PAW_PBE:06Sep2000
  if(AUID=="37c5d78c5699ee73" && nKIN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-3.79514;volume_atom=28.3410;spin_atom=0.0;} // Sn:PAW_PBE:08Apr2002
  if(AUID=="4ddcdf83c7589054" && SCAN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-35.8121;volume_atom=27.4595;spin_atom=0.0;} // Sn:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="6431a791389bb4f7" && SCAN) {found=TRUE;groundstate_structure="A5";groundstate_energy=-35.7984;volume_atom=27.2303;spin_atom=0.0;} // Sn_d:PAW_PBE_KIN:SCAN:06Sep2000

  // Sr
  if(AUID=="82ec420862ed9432" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.73980;volume_atom=40.5798;spin_atom=0.0;} // Sr:GGA:01Apr2000
  if(AUID=="6fb0d22509f8b232" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.92701;volume_atom=41.8890;spin_atom=0.0;} // Sr:LDA:01Apr2000
  if(AUID=="24554a4dcb7ec11b" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.61166;volume_atom=54.1148;spin_atom=0.0;} // Sr_pv:GGA:01Apr2000
  if(AUID=="b3f00f0671094a1c" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.90232;volume_atom=48.0657;spin_atom=0.0;} // Sr_pv:LDA:01Apr2000
  if(AUID=="bbfa6aec08364643" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.62474;volume_atom=53.7508;spin_atom=0.0;} // Sr_sv:PAW_GGA:10Feb1998
  if(AUID=="3951da71de7f3413" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.88167;volume_atom=48.2867;spin_atom=0.0;} // Sr_sv:PAW_LDA:10Feb1998
  if(AUID=="68b7a60c49c97739" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.68354;volume_atom=54.6082;spin_atom=0.0631431;} /// Sr_sv:PAW_PBE:07Sep2000
  if(AUID=="403b3e0afc9132cc" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-22.3021;volume_atom=54.0609;spin_atom=0.156893;} // Sr_sv:PAW_PBE_KIN:SCAN:07Sep2000

  // Ta
  if(AUID=="ecb6cb1d273c3dcf" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.8043;volume_atom=17.8706;spin_atom=0.0;} // Ta:GGA:01Apr2000
  if(AUID=="c88d4f9a05e45d41" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.7802;volume_atom=18.0957;spin_atom=0.0;} // Ta:PAW_GGA:06Feb2003
  if(AUID=="d1561ec02d179810" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.8929;volume_atom=17.1116;spin_atom=0.0;} // Ta:PAW_LDA:21Jan2003
  if(AUID=="d83777bdbc002be9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.8600;volume_atom=18.0798;spin_atom=0.0;} // Ta:PAW_PBE:17Jan2003
  if(AUID=="de2a3eb2f27e9633" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.7319;volume_atom=18.2562;spin_atom=0.0;} // Ta_pv:PAW_GGA:07Sep2000
  if(AUID=="926ce9924f7a5fd2" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-11.8489;volume_atom=18.2778;spin_atom=0.0;} // Ta_pv:PAW_PBE:07Sep2000
  if(AUID=="6d000fcf1ccde6da" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.8285;volume_atom=17.2465;spin_atom=0.0;} // Ta_pv:PAW_LDA:07Sep2000
  if(AUID=="096d94d0306813d6" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-66.6819;volume_atom=17.5029;spin_atom=0.0;} // Ta:PAW_PBE_KIN:SCAN:17Jan2003
  if(AUID=="2287798c5502b426" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-66.6142;volume_atom=17.5162;spin_atom=0.0;} // Ta_pv:PAW_PBE_KIN:SCAN:07Sep2000

  // Tc
  if(AUID=="fd329b53fcbb5dbd" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.2088;volume_atom=14.4754;spin_atom=0.0;} // Tc:GGA:01Apr2000
  if(AUID=="aabb0d26abfc5599" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.3047;volume_atom=14.4300;spin_atom=0.0;} // Tc:PAW_PBE:17Jan2003
  if(AUID=="34c9772cc010d2cb" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.6391;volume_atom=13.7601;spin_atom=0.0;} // Tc:LDA:01Apr2000
  if(AUID=="80004c8d11aa64f2" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.2801;volume_atom=14.4848;spin_atom=0.0;} // Tc:PAW_GGA:21Dec2000
  if(AUID=="0eb01e3427b59ddc" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.2038;volume_atom=14.5814;spin_atom=0.0;} // Tc_pv:PAW_GGA:20Feb1998
  if(AUID=="43f93be60f147933" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.6972;volume_atom=13.8130;spin_atom=0.0;} // Tc:PAW_LDA:21Jan2003
  if(AUID=="c33964edb4f74ebf" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-10.3597;volume_atom=14.5515;spin_atom=0.0;} // Tc_pv:PAW_PBE:06Sep2000
  if(AUID=="aaf18f3b4fa31411" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-11.6002;volume_atom=13.9073;spin_atom=0.0;} // Tc_pv:PAW_LDA:03Mar1998
  if(AUID=="5751220a43dddd9b" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-35.0200;volume_atom=13.9455;spin_atom=0.0;} // Tc:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="5c4bf0cceb477531" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-35.0765;volume_atom=14.0558;spin_atom=0.0;} // Tc_pv:PAW_PBE_KIN:SCAN:04Feb2005
  if(AUID=="3816a915d663bd59" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-35.0857;volume_atom=14.0351;spin_atom=0.0;} // Tc_sv:PAW_PBE_KIN:SCAN:23Mar2010

  // Te
  if(AUID=="96c99e8488c592ae" && nKIN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-3.14141;volume_atom=34.8509;spin_atom=0.0;} // Te:PAW_PBE:08Apr2002
  if(AUID=="8d2f6568d1d06421" && SCAN) {found=TRUE;groundstate_structure="A8";groundstate_energy=-37.2600;volume_atom=33.9467;spin_atom=0.0;} // Te:PAW_PBE_KIN:SCAN:08Apr2002

  // Ti
  if(AUID=="e97d34c6b4cd62df" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.53337;volume_atom=15.9367;spin_atom=0.0;} // Ti:LDA:01Apr2000
  if(AUID=="8498e164009c60f7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.73332;volume_atom=17.2006;spin_atom=0.0;} // Ti:GGA:01Apr2000
  if(AUID=="264e8b27aa49e864" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.53301;volume_atom=15.9106;spin_atom=0.0;} // Ti:PAW_LDA:03Oct2001
  if(AUID=="47d8726ba10b866f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.74356;volume_atom=17.0410;spin_atom=0.0;} // Ti:PAW_GGA:08Aug2001
  if(AUID=="e18697ab9588a353" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.76396;volume_atom=17.1060;spin_atom=0.0;} // Ti:PAW_PBE:08Apr2002
  if(AUID=="257a6b2b7e91023a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.75293;volume_atom=17.5203;spin_atom=0.0;} // Ti:GGA:01Apr2000
  if(AUID=="8e962327621a06f8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.51437;volume_atom=16.2957;spin_atom=0.0;} // Ti:LDA:01Apr2000
  if(AUID=="8fecbeeb2db133cb" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.77528;volume_atom=17.1860;spin_atom=0.0;} // Ti_pv:PAW_GGA:07Sep2000
  if(AUID=="c8c8567f19acfb26" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.80066;volume_atom=17.1803;spin_atom=0.0;} // Ti_sv:PAW_GGA:07Sep2000
  if(AUID=="567ef720cface091" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.54336;volume_atom=16.0264;spin_atom=0.0;} // Ti_pv:PAW_LDA:07Sep2000
  if(AUID=="3f1fd5d5748fffbc" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.89056;volume_atom=17.2662;spin_atom=0.0;} // Ti_pv:PAW_PBE:07Sep2000
  if(AUID=="9cdd346558c528c6" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-7.93863;volume_atom=17.2475;spin_atom=0.0;} // Ti_sv:PAW_PBE:07Sep2000
  if(AUID=="5b9d5db647f6fc85" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.53889;volume_atom=16.0094;spin_atom=0.0;} // Ti_sv:PAW_LDA:07Sep2000
  if(AUID=="e19a9be220c9f36d" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-15.9523;volume_atom=16.9114;spin_atom=0.0;} // Ti:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="571036b918ee987a" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-15.9827;volume_atom=17.0821;spin_atom=0.0;} // Ti_pv:PAW_PBE_KIN:SCAN:07Sep2000
  if(AUID=="9695cfc266c5809c" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-16.0198;volume_atom=17.2317;spin_atom=0.0;} // Ti_sv:PAW_PBE_KIN:SCAN:26Sep2005

  // Tl
  if(AUID=="9f2d6c60104e71a3" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.25921;volume_atom=30.8473;spin_atom=0.0;} // Tl_d:GGA:01Apr2000
  if(AUID=="06fbb1e6c3452edb" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.93920;volume_atom=27.1009;spin_atom=0.0;} // Tl_d:LDA:01Apr2000
  if(AUID=="958e09d75d34dcd0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.23399;volume_atom=30.4705;spin_atom=0.0;} // Tl_d:PAW_GGA:11Feb1998
  if(AUID=="1a0488ce9d93df98" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.90255;volume_atom=26.8803;spin_atom=0.0;} // Tl_d:PAW_LDA:11Feb1998
  if(AUID=="b88a3147de55f4a8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.24350;volume_atom=31.1190;spin_atom=0.0;} // Tl:PAW_PBE:08Apr2002
  if(AUID=="f473cbb4831e0bc0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.36274;volume_atom=30.7601;spin_atom=0.0;} // Tl_d:PAW_PBE:06Sep2000
  if(AUID=="1298edb6f5cd44f1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.28269;volume_atom=30.9588;spin_atom=0.0;} // Tl:PAW_GGA:25Jul2001
  if(AUID=="2ad3380a31e82af9" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-2.95325;volume_atom=27.3640;spin_atom=0.0;} // Tl:PAW_LDA:03Oct2001
  if(AUID=="0e2337aac65480f2" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-72.3276;volume_atom=29.3500;spin_atom=0.0;} // Tl:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="4f0fac6e26fd86a9" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-72.3189;volume_atom=29.0596;spin_atom=0.0;} // Tl_d:PAW_PBE_KIN:SCAN:06Sep2000

  // V
  if(AUID=="aee485cb636ab31d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.90719;volume_atom=13.3310;spin_atom=0.0;} // V:GGA:01Apr2000
  if(AUID=="40048762e7327671" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.95900;volume_atom=12.3899;spin_atom=0.0;} // V:LDA:01Apr2000
  if(AUID=="dce68ed9d844d970" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.92373;volume_atom=13.1803;spin_atom=0.0;} // V:PAW_GGA:07Aug2001
  if(AUID=="a741261c8d2045c9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.95883;volume_atom=12.3350;spin_atom=0.0;} // V:PAW_LDA:07Aug2001
  if(AUID=="183846be5b03c64d" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.94402;volume_atom=13.2022;spin_atom=0.0;} // V:PAW_PBE:08Apr2002
  if(AUID=="ea4e243d3a0747f9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.99426;volume_atom=12.6635;spin_atom=0.0;} // V:LDA:01Apr2000
  if(AUID=="1dd37802a9ba101a" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.92582;volume_atom=13.3629;spin_atom=0.0;} // V_pv:PAW_GGA:07Sep2000
  if(AUID=="d16570ec1e90da05" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.99489;volume_atom=13.5698;spin_atom=0.0;} // V:GGA:01Apr2000
  if(AUID=="0725afc42e47381e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.93464;volume_atom=12.5085;spin_atom=0.0;} // V_pv:PAW_LDA:07Sep2000
  if(AUID=="b7d3211b63a6261e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.07822;volume_atom=13.3934;spin_atom=0.0;} // V_pv:PAW_PBE:07Sep2000
  if(AUID=="6555390e553533e4" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-8.95881;volume_atom=13.3628;spin_atom=0.0;} // V_sv:PAW_GGA:14Sep2000
  if(AUID=="c3b5d34259b9b5f9" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.93559;volume_atom=12.5083;spin_atom=0.0;} // V_sv:PAW_LDA:07Sep2000
  if(AUID=="94da8dda3e65c762" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-9.11496;volume_atom=13.3971;spin_atom=0.0;} // V_sv:PAW_PBE:07Sep2000
  if(AUID=="62ba629379463606" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-17.3619;volume_atom=12.8510;spin_atom=0.0;} // V:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="46928c16c2f8a488" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-17.4590;volume_atom=13.0021;spin_atom=0.0;} // V_sv:PAW_PBE_KIN:SCAN:02Aug2007
  if(AUID=="ff7f68300ba3e357" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-17.4221;volume_atom=13.0707;spin_atom=0.0;} // V_pv:PAW_PBE_KIN:SCAN:07Sep2000

  // W
  if(AUID=="fc6e347dc7ddd7ab" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.8623;volume_atom=15.9358;spin_atom=0.0;} // W:GGA:01Apr2000
  if(AUID=="4a36ee2cfb9fae34" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.9177;volume_atom=15.9658;spin_atom=0.0;} // W:PAW_GGA:21Dec2000
  if(AUID=="6a5db51b820aa7a5" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-14.2099;volume_atom=15.2787;spin_atom=0.0;} // W:PAW_LDA:19Jan2001
  if(AUID=="c895658f3e6606f4" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-14.1883;volume_atom=15.1794;spin_atom=0.0;} // W:LDA:01Apr2000
  if(AUID=="22a83c4c0d26e6eb" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.7816;volume_atom=16.1779;spin_atom=0.0;} // W_pv:PAW_GGA:15Jul1998
  if(AUID=="e05803166f7d0e3e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-14.0476;volume_atom=15.4745;spin_atom=0.0;} // W_pv:PAW_LDA:22Jul1998
  if(AUID=="cad81c4b63d2d96e" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-12.9546;volume_atom=16.1649;spin_atom=0.0;} // W_pv:PAW_PBE:06Sep2000
  if(AUID=="1febde78746c14e4" && nKIN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-13.0124;volume_atom=15.9045;spin_atom=0.0;} // W:PAW_PBE:08Apr2002
  if(AUID=="b99d5a263927c2ea" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-69.9743;volume_atom=15.5134;spin_atom=0.0;} // W:PAW_PBE_KIN:SCAN:08Apr2002
  if(AUID=="1d8c010cef748b9c" && SCAN) {found=TRUE;groundstate_structure="A2";groundstate_energy=-69.6572;volume_atom=15.6161;spin_atom=0.0;} // W_sv:PAW_PBE_KIN:SCAN:04Sep2015

  // Y
  if(AUID=="377b8f487d291e0e" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.33584;volume_atom=32.5652;spin_atom=0.0;} // Y:GGA:01Apr2000
  if(AUID=="a0835c00b2263c8b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.40419;volume_atom=31.9468;spin_atom=0.0;} // Y:GGA:01Apr2000
  if(AUID=="1844bc40a30b31c7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.90740;volume_atom=29.8089;spin_atom=0.0;} // Y_sv:PAW_LDA:10Feb1998
  if(AUID=="20813d90f7aff62f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.38136;volume_atom=32.4851;spin_atom=0.0;} // Y_sv:PAW_GGA:10Feb1998
  if(AUID=="18ed87089235423b" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.98897;volume_atom=29.2790;spin_atom=0.0;} // Y:LDA:01Apr2000
  if(AUID=="5c8efcbfc77486c0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.89015;volume_atom=29.8469;spin_atom=0.0;} // Y:LDA:01Apr2000
  if(AUID=="b07bf0a656aa22bc" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-6.46317;volume_atom=32.7578;spin_atom=0.0;} // Y_sv:PAW_PBE:06Sep2000
  if(AUID=="2fb876577258940d" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-27.8500;volume_atom=33.3097;spin_atom=0.602783;} // Y_sv:PAW_PBE_KIN:SCAN:25May2007   MAGNETIC ??

  // Yb
  if(AUID=="4e35fa207bf4aaa9" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.51978;volume_atom=39.7380;spin_atom=0.0;} // Yb_2:PAW_PBE:06Sep2000
  if(AUID=="b1f8bc4d83e97833" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.44398;volume_atom=39.1468;spin_atom=0.0;} // Yb_2:PAW_GGA:10May2000
  if(AUID=="de1c22a384a8d945" && nKIN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-1.67034;volume_atom=36.9163;spin_atom=0.0;} // Yb:PAW_PBE:24Feb2003
  if(AUID=="a2a66bbb705fd114" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-50.5319;volume_atom=37.1321;spin_atom=0.0;} // Yb_2:PAW_PBE_KIN:SCAN:06Sep2000
  if(AUID=="d4d1484e34700879" && SCAN) {found=TRUE;groundstate_structure="A1";groundstate_energy=-53.5119;volume_atom=28.3498;spin_atom=0.0;} // Yb_3:PAW_PBE_KIN:SCAN:08Jul2013

  // Zn
  if(AUID=="c765a5877b95887a" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.12141;volume_atom=15.3617;spin_atom=0.0;} // Zn:GGA:01Apr2000
  if(AUID=="e75d47ef817f48d8" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.11328;volume_atom=15.2205;spin_atom=0.0;} // Zn:PAW_GGA:03Mar1998
  if(AUID=="4b8b99f638e1a173" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.86830;volume_atom=13.5768;spin_atom=0.0;} // Zn:PAW_LDA:03Mar1998
  if(AUID=="d321d33f5619ef83" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.89186;volume_atom=13.5675;spin_atom=0.0;} // Zn:LDA:01Apr2000
  if(AUID=="d4ad14ef9da00329" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-1.26581;volume_atom=15.1693;spin_atom=0.0;} // Zn:PAW_PBE:06Sep2000
  if(AUID=="7ed9d78dd30715fd" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-13.4842;volume_atom=14.1685;spin_atom=0.0;} // Zn:PAW_PBE_KIN:SCAN:06Sep2000

  // Zr
  if(AUID=="c87b485bdb0c4683" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.36116;volume_atom=21.3443;spin_atom=0.0;} // Zr:LDA:01Apr2000
  if(AUID=="014ba38906f94787" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.47756;volume_atom=23.4023;spin_atom=0.0;} // Zr:PAW_PBE:08Apr2002
  if(AUID=="05c7a473e8b0fed1" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.45764;volume_atom=23.3305;spin_atom=0.0;} // Zr:PAW_GGA:08Aug2001
  if(AUID=="a2c98479cda145dd" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.49915;volume_atom=22.8558;spin_atom=0.0;} // Zr:GGA:01Apr2000
  if(AUID=="844977898c624f61" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.22027;volume_atom=21.8813;spin_atom=0.0;} // Zr_sv:PAW_LDA:10Feb1998
  if(AUID=="ae40f9ef1519f7f7" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.40052;volume_atom=23.4147;spin_atom=0.0;} // Zr:GGA:01Apr2000
  if(AUID=="7471a45b48d5cbca" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.43218;volume_atom=23.3391;spin_atom=0.0;} // Zr_sv:PAW_GGA:10Feb1998
  if(AUID=="3f478742b0dd98d0" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-9.21967;volume_atom=21.9211;spin_atom=0.0;} // Zr:LDA:01Apr2000
  if(AUID=="8f81b69844f3963f" && nKIN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-8.54365;volume_atom=23.4268;spin_atom=0.0;} // Zr_sv:PAW_PBE:07Sep2000
  if(AUID=="a3ef4a907c7b4c96" && SCAN) {found=TRUE;groundstate_structure="A3";groundstate_energy=-30.6972;volume_atom=23.1863;spin_atom=0.0;} // Zr_sv:PAW_PBE_KIN:SCAN:04Jan2005

  /*
    ./xgo B:PAW_GGA:18Jul2000 "found=TRUE;groundstate_structure=\"ICSD_240995\";groundstate_energy=0;volume_atom=7.76493;spin_atom=0.0;"// FIX
    ./xgo B_s:PAW_GGA:21Dec2000 "found=TRUE;groundstate_structure=\"ICSD_240995\";groundstate_energy=-6.52874;volume_atom=7.79803;spin_atom=0.0;"// FIX
    ./xgo B_h:PAW_GGA:18Jul2000 "found=TRUE;groundstate_structure=\"ICSD_108026\";groundstate_energy=-6.177995603333333;volume_atom=6.717450;spin_atom=0.0;"// FIX
    // OLD ./xgo B:PAW_PBE:06Sep2000 "found=TRUE;groundstate_structure=\"A3\";groundstate_energy=-5.96676;volume_atom=9.11282;spin_atom=0.0;"
    // /xgo B_s:PAW_PBE:22Jan2003 "found=TRUE;groundstate_structure=\"ICSD_240995\";groundstate_energy=-6.52874;volume_atom=7.79803;spin_atom=0.0;"// FIX
    // OLD ./xgo B_h:PAW_PBE:07Sep2000 "found=TRUE;groundstate_structure=\"A3\";groundstate_energy=-5.98337;volume_atom=9.39294;spin_atom=0.0;"
    // OLD ./xgo B_s:PAW_PBE:22Jan2003 "found=TRUE;groundstate_structure=\"A3\";groundstate_energy=-5.9736;volume_atom=8.8771;spin_atom=0.0;"

    ./xgo Sm_3:PAW_GGA:11May2000 "found=TRUE;groundstate_structure=\"ICSD_246657\";groundstate_energy=-4.621400;volume_atom=33.447633;spin_atom=0.0;"// FIX
    ./xgo Sm_3:PAW_GGA:11May2000 && 0 "found=TRUE;groundstate_structure=\"ICSD_652637\";groundstate_energy=-4.64136;volume_atom=33.5075;spin_atom=0.0;"// IT HAS LDAU
    ./xgo Sm_3:PAW_PBE:07Sep2000 && 0 "found=TRUE;groundstate_structure=\"A1\";groundstate_energy=-4.7062;volume_atom=33.8339;spin_atom=0.0;"// IT HAS LDAU

    // ./xgo Ce "found=TRUE;groundstate_structure=\"A1\";groundstate_energy=-5.92998;volume_atom=26.0579;spin_atom=0.0;"
    // ./xgo Ce "found=TRUE;groundstate_structure=\"ICSD_2284-mS4\";groundstate_energy=-5.93013;volume_atom=26.0697;spin_atom=0.0;"
    // ./xgo Cl_h:PAW_PBE:08Apr2002 "found=TRUE;groundstate_structure=\"A11\";groundstate_energy=-1.8156;volume_atom=37.3299;spin_atom=0.0;"WAITING

    */
  if(!found) { volume_atom=999999,spin_atom=999999;} // some defaults
  //  if(!found) cerr <<"ERROR (xPOTCAR_EnthalpyReference_AUID): NOT FOUND: AUID=" << AUID << endl;// exit(0);
  if(LDEBUG && !found) cout <<"ERROR (xPOTCAR_EnthalpyReference_AUID): NOT FOUND: AUID=" << AUID << endl;// exit(0);
  if(LDEBUG &&  found) cout <<"ERROR (xPOTCAR_EnthalpyReference_AUID): FOUND: AUID=" << AUID << endl;// exit(0);
  if(LDEBUG) cerr << "xPOTCAR_EnthalpyReference_AUID: [END]" << endl;
  return found;
};



#endif // _AFLOW_XPSEUDOPOTENTIAL_CPP

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2020              *
// *                                                                        *
// **************************************************************************
