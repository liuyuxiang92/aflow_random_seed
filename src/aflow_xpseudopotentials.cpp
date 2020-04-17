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
  make -j
 
  ./aflow --scrub=POTCAR --FILE /common/VASP/potpaw_PBE/current/Mo_pv/POTCAR
  ./aflow --scrub=OUTCAR --FILE /common/LIB3/LIB/AgCdCo/TFCC001.ABC/OUTCAR.relax2.xz 
  ./aflow --pseudopotentials_check=/tmp/POTCAR1
  ./aflow --pseudopotentials_check=/tmp/POTCAR2
  ./aflow --pseudopotentials_check=/tmp/OUTCAR.relax2
  ./aflow --pseudopotentials_check=/common/LIB3/LIB/AgCdCo/TFCC001.ABC/OUTCAR.relax2.xz 

  ./aflow --use_aflow.in=aflow.in --beep --force --lib2raw="/common/LIB3/LIB/TeW_pvY_sv/TFCC001.ABC"
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
  if(XHOST.PSEUDOPOTENTIAL_GENERATOR && xPOT.species.size()==1) {  // SC20200326
    string comment="";
    // objects/functions for references energies  // SC20200326
    //   vector<string> tokens;
    //   aurostd::string2tokens(xPOT.species_pp_version.at(0),tokens,":");
    // [OBSOLETE] if(tokens.size()>0) {vdate.clear();vdate.push_back(tokens.at(tokens.size()-1));};
    comment=xPOT.species_pp_version.at(0);
    xPOT.vTITEL.at(0)=aurostd::RemoveWhiteSpaces(xPOT.vTITEL.at(0));
    xPOT.vLEXCH.at(0)=aurostd::RemoveWhiteSpaces(xPOT.vLEXCH.at(0));
    if(LVERBOSE) {cerr << "xPOTCAR::GetProperties: vTITEL.at(0)=" << xPOT.vTITEL.at(0) << endl;}   // SC20200326
    if(LVERBOSE) {cerr << "xPOTCAR::GetProperties: species.at(0)=" << xPOT.species.at(0) << endl;}   // SC20200326
    if(LVERBOSE) {cerr << "xPOTCAR::GetProperties: species_Z.at(0)=" << xPOT.species_Z.at(0) << endl;}    // SC20200326
    oss << "  " << endl;
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    //   oss << "  // [AFLOW]START=" << comment << " " << endl;
    oss << "  // " << comment << " " << comment << " " << comment << " " << comment << " " << endl;
    oss << "  // " << xPOT.filename << endl;    // SC20200326
    oss << "  " << aurostd::PaddedPOST("{",PSEUDOPOTENTIAL_GENERATOR_pad) << "      // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("xPOTCAR x;",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.filename=\""+xPOT.filename+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
     oss << "    " << aurostd::PaddedPOST("x.AUID=\""+xPOT.AUID+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vTITEL.push_back(\""+xPOT.vTITEL.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.pp_type=\""+xPOT.pp_type+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species.push_back(\""+xPOT.species.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_Z.push_back("+aurostd::utype2string<int>(xPOT.species_Z.at(0))+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp.push_back(\""+xPOT.species_pp.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_type.push_back(\""+xPOT.species_pp_type.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_version.push_back(\""+xPOT.species_pp_version.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_AUID.push_back(\""+xPOT.AUID+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_energy.push_back("+aurostd::utype2string<double>(NNN)+");//"+xPOT.AUID,PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.species_pp_groundstate_structure.push_back(\"N/A_"+xPOT.AUID+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_PAW="+aurostd::bool2string(xPOT.POTCAR_PAW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_TYPE=\""+xPOT.POTCAR_TYPE+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_KINETIC="+aurostd::bool2string(xPOT.POTCAR_KINETIC)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_GW="+aurostd::bool2string(xPOT.POTCAR_GW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.POTCAR_AE="+aurostd::bool2string(xPOT.POTCAR_AE)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // EATOM RCORE RWIGS EAUG RAUG ENMAX ENMIN POMASS ZVAL RMAX LEXCH
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vENMAX.push_back("+aurostd::utype2string<double>(xPOT.vENMAX.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vENMIN.push_back("+aurostd::utype2string<double>(xPOT.vENMIN.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vPOMASS.push_back("+aurostd::utype2string<double>(xPOT.vPOMASS.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vZVAL.push_back("+aurostd::utype2string<double>(xPOT.vZVAL.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vEATOM.push_back("+aurostd::utype2string<double>(xPOT.vEATOM.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vRCORE.push_back("+aurostd::utype2string<double>(xPOT.vRCORE.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vRWIGS.push_back("+aurostd::utype2string<double>(xPOT.vRWIGS.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vEAUG.push_back("+aurostd::utype2string<double>(xPOT.vEAUG.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    // [NON_NECESSARY] oss << "    " << aurostd::PaddedPOST("x.vRAUG.push_back("+aurostd::utype2string<double>(xPOT.vRAUG.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vRMAX.push_back("+aurostd::utype2string<double>(xPOT.vRMAX.at(0),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("x.vLEXCH.push_back(\""+xPOT.vLEXCH.at(0)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vxpseudopotential.push_back(x);",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "  " << aurostd::PaddedPOST("}",PSEUDOPOTENTIAL_GENERATOR_pad) << "      // " << comment << endl;    // SC20200326
    oss << "  // " << xPOT.filename << endl;    // SC20200326
    //   oss << "  // [AFLOW]STOP=" << comment << " " << endl;
    oss << "  // ******************************************************************************************************************************************************** " << endl;
    oss << endl;    // SC20200326
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
    oss << "  // " << xPOT.filename << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("filename=\""+xPOT.filename+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("AUID=\""+xPOT.AUID+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vTITEL.push_back(\""+xPOT.vTITEL.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("pp_type=\""+xPOT.pp_type+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("species.push_back(\""+xPOT.species.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_Z.push_back("+aurostd::utype2string<int>(xPOT.species_Z.at(i))+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp.push_back(\""+xPOT.species_pp.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_type.push_back(\""+xPOT.species_pp_type.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_version.push_back(\""+xPOT.species_pp_version.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_AUID.push_back(\""+xPOT.species_pp_AUID.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_groundstate_energy.push_back("+aurostd::utype2string<double>(xPOT.species_pp_groundstate_energy.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("species_pp_groundstate_structure.push_back(\""+xPOT.species_pp_groundstate_structure.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_PAW="+aurostd::bool2string(xPOT.POTCAR_PAW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_TYPE=\""+xPOT.POTCAR_TYPE+"\";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_KINETIC="+aurostd::bool2string(xPOT.POTCAR_KINETIC)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_GW="+aurostd::bool2string(xPOT.POTCAR_GW)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("POTCAR_AE="+aurostd::bool2string(xPOT.POTCAR_AE)+";",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vENMAX.push_back("+aurostd::utype2string<double>(xPOT.vENMAX.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vENMIN.push_back("+aurostd::utype2string<double>(xPOT.vENMIN.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vPOMASS.push_back("+aurostd::utype2string<double>(xPOT.vPOMASS.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vZVAL.push_back("+aurostd::utype2string<double>(xPOT.vZVAL.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vEATOM.push_back("+aurostd::utype2string<double>(xPOT.vEATOM.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vRCORE.push_back("+aurostd::utype2string<double>(xPOT.vRCORE.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vRWIGS.push_back("+aurostd::utype2string<double>(xPOT.vRWIGS.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vEAUG.push_back("+aurostd::utype2string<double>(xPOT.vEAUG.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;   // SC20200326
    oss << "    " << aurostd::PaddedPOST("vRAUG.push_back("+aurostd::utype2string<double>(xPOT.vRAUG.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;   // SC20200326
    oss << "    " << aurostd::PaddedPOST("vRMAX.push_back("+aurostd::utype2string<double>(xPOT.vRMAX.at(i),10)+");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "    " << aurostd::PaddedPOST("vLEXCH.push_back(\""+xPOT.vLEXCH.at(i)+"\");",PSEUDOPOTENTIAL_GENERATOR_pad) << "    // " << comment << endl;    // SC20200326
    oss << "  // " << xPOT.filename << endl;    // SC20200326
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

#endif // _AFLOW_XPSEUDOPOTENTIAL_CPP

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2020              *
// *                                                                        *
// **************************************************************************
