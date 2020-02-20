//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
//****************************************************************************

#include "aflow_apl.h"

using std::string;
using std::vector;

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            INPUT FILE CREATORS                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //createAflowInPhonons////////////////////////////////////////////////////////
  // ME 181022 - New method to create the aflow.in files. Uses the aflow.in
  // creator in aflow_avasp.cpp
  bool createAflowInPhonons(const _aflags& _aflowFlags, const _kflags& _kbinFlags, const _xflags& _xFlags, _xinput& xinp) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="apl::createAflowInPhonons():";
    bool write = false;
    if (xinp.AFLOW_MODE_VASP) {
      if(LDEBUG){
        std::cerr << soliloquy << " BEFORE xinp.xvasp.Directory=" << xinp.xvasp.Directory << std::endl;
        std::cerr << soliloquy << " BEFORE xinp.getDirectory()=" << xinp.getDirectory() << std::endl;
      }
      AVASP_populateXVASP(_aflowFlags, _kbinFlags, _xFlags.vflags, xinp.xvasp);
      if(LDEBUG){
        std::cerr << soliloquy << " AFTER xinp.xvasp.Directory=" << xinp.xvasp.Directory << std::endl;
        std::cerr << soliloquy << " AFTER xinp.getDirectory()=" << xinp.getDirectory() << std::endl;
      }
      xinp.setDirectory(xinp.xvasp.Directory);
      if (!filesExistPhonons(xinp)) {
        stringstream aflowin;
        write = true;
        AVASP_MakeSingleAFLOWIN(xinp.xvasp, aflowin, write);
      } else {
        write = false;
      }
    }
    return write;
  }

  // ME181022 - Old method to create aflow.in files for AIMS
  void createAflowInPhononsAIMS(_aflags& _aflowFlags, _kflags& _kbinFlags, _xflags& _xFlags, string& _AflowIn, _xinput& xinp, ofstream& fileMessage) {
    if (!xinp.AFLOW_MODE_AIMS) {
      string function = "apl::createAflowInPhononsAIMS()";
      string message = "This function only works with AIMS.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }
    xinp.xaims.CONTROL.str(std::string());
    KBIN::AIMS_Produce_CONTROL(xinp.xaims,_AflowIn,fileMessage,_aflowFlags,_kbinFlags,_xFlags.aimsflags);  //DEFAULT
    KBIN::AIMS_Modify_CONTROL(xinp.xaims,fileMessage,_aflowFlags,_kbinFlags,_xFlags.aimsflags);            //DEFAULT

    // Write aflow.in

    //copying from createAFLOWIN
    _vflags vflags(_xFlags.vflags);
    _xaims xaims(xinp.xaims);
    _aimsflags aimsflags(_xFlags.aimsflags);

    string directory=xinp.getDirectory();
    if(directory.empty()){
      string function = "apl::createAflowInPhonons()";
      string message =  "no output directory found";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }

    if(!aurostd::FileExist(directory)){aurostd::DirectoryMake(directory);}  // Create directory if it is not created
    aurostd::DirectoryChmod("777", directory);                              // CHMOD Directory 777

    stringstream outfile;

    // OK, fill it...
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    outfile << "[AFLOW] _ ___ _" << std::endl;
    outfile << "[AFLOW] / \\| o \\ |" << std::endl;
    outfile << "[AFLOW] | o | _/ |_" << std::endl;
    outfile << "[AFLOW] |_n_|_| |___| automatic generated file" << std::endl;
    outfile << "[AFLOW]" << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    outfile << "[AFLOW_MODE=AIMS]" << std::endl;
    if(!_kbinFlags.KZIP_BIN.empty()){outfile << "[AFLOW_MODE_ZIP=" << _kbinFlags.KZIP_BIN << "]" << std::endl;}  //CO

    //CO 180130 - START
    //corey - at some point, fix alien mode for aims, for now omit!
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    outfile << "[AIMS_CONTROL_MODE_EXPLICIT]START " << std::endl;
    outfile << xaims.CONTROL.str();
    outfile << "[AIMS_CONTROL_MODE_EXPLICIT]STOP " << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;
    outfile << "[AIMS_GEOM_MODE_EXPLICIT]START " << std::endl;
    outfile << xaims.str;
    outfile << "[AIMS_GEOM_MODE_EXPLICIT]STOP " << std::endl;
    outfile << AFLOWIN_SEPARATION_LINE << std::endl;

    //also write out
    if(1){
      KBIN::AIMS_Write_CONTROL(xaims,aimsflags);
      xaims.GEOM.clear(); xaims.GEOM.str("");
      xaims.GEOM << xaims.str;

      string geom_filename = xaims.Directory + "/" + AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM;
      aurostd::stringstream2file(xaims.GEOM, geom_filename);
      if(!aurostd::FileExist(geom_filename)){
        string function = "apl::createAflowInPhonons()";
        string message = "Cannot create [" + AFLOWRC_DEFAULT_AIMS_EXTERNAL_GEOM + "] file.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
      }
      aurostd::ChmodFile("a+rw", geom_filename);
    }

    //CO - START
    string filename = directory + string("/") + _AFLOWIN_;
    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)){
      string function = "apl::createAflowInPhonons()";
      string message = "Cannot create [" + _AFLOWIN_ + "] file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
    aurostd::ChmodFile("a+rw", filename); // CHMOD a+rw _AFLOWIN_
    //CO - END
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           OUTPUT FILE READERS                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //outfileFoundAnywherePhonons/////////////////////////////////////////////////
  bool outfileFoundAnywherePhonons(vector<_xinput>& xinps) {
    for (uint idxRun = 0; idxRun < xinps.size(); idxRun++) {
      string dir = xinps[idxRun].getDirectory();

      if(xinps[idxRun].AFLOW_MODE_VASP) {
        if(aurostd::EFileExist(dir + string("/vasprun.xml.static")) ||
            aurostd::EFileExist(dir + string("/vasprun.xml")) ||
            aurostd::EFileExist(dir + "/" + DEFAULT_AFLOW_QMVASP_OUT)) {  // ME190607
          return true;
        }
      }
      if(xinps[idxRun].AFLOW_MODE_AIMS) {
        if(aurostd::EFileExist(xinps[idxRun].getDirectory() + string("/aims.out"))) {
          return true;
        }
      }
    }
    return false;
  }

  //filesExistPhonons/////////////////////////////////////////////////////////
  bool filesExistPhonons(_xinput& xinp) {
    string dir = xinp.getDirectory() + string("/");
    if (aurostd::FileExist(dir + _AFLOWIN_)) {
      return true;  //do not OVERWRITE an aflow.in
    }
    if (xinp.AFLOW_MODE_VASP){
      if(aurostd::EFileExist(dir + string("vasprun.xml.static")) ||
          aurostd::EFileExist(dir + string("vasprun.xml")) ||
          aurostd::EFileExist(dir + DEFAULT_AFLOW_QMVASP_OUT)) {  // ME200203 - Added qmvasp file
        return true;
      }
    }
    if(xinp.AFLOW_MODE_AIMS){
      if(aurostd::EFileExist(dir + string("aims.out"))) {
        return true;
      }
    }
    return false;
  }

  //outfileFoundEverywherePhonons/////////////////////////////////////////////
  void outfileFoundEverywherePhonons(vector<_xinput>& xinps, ofstream& fileMessage, bool contains_born) {
    string function = "apl::outfileFoundEverywherePhonons()";
    stringstream _logger;
    _logger << "Reading force files";
    pflow::logger(_AFLOW_FILE_NAME_, "AAPL", _logger, fileMessage, std::cout);
    uint ninps = xinps.size();
    if (contains_born) ninps--;
    for (uint idxRun = 0; idxRun < ninps; idxRun++) {
      _logger << "Reading force file " << (idxRun + 1) << "/" << ninps << "."; //CO190116  // ME190607
      pflow::logger(_AFLOW_FILE_NAME_, "AAPL", _logger, fileMessage, std::cout);
      string directory = xinps[idxRun].getDirectory();
      readForcesFromDirectory(xinps[idxRun]);

      // Was it all right?
      if (!xinps[idxRun].getXStr().qm_calculated) {
        string message = "The force file in " + directory + " is wrong.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
      }
    }
    _logger << "No errors caught, all force files read successfully."; //CO190116  // ME190607
    pflow::logger(_AFLOW_FILE_NAME_, "AAPL", _logger, fileMessage, std::cout, 'C');
  }

  void readForcesFromDirectory(_xinput& xinp) {
    string function = "apl::readForcesFromDirectory()";
    uint natoms = xinp.getXStr().atoms.size();
    xinp.getXStr().qm_forces.clear();
    // Load data....
    if(xinp.AFLOW_MODE_VASP) {
      if (!aurostd::EFileExist(xinp.getDirectory() + "/" + DEFAULT_AFLOW_QMVASP_OUT)
        && !aurostd::EFileExist(xinp.getDirectory() + "/vasprun.xml.static")
        && !aurostd::EFileExist(xinp.getDirectory() + "/vasprun.xml")) {
        string message = "The force file in " + xinp.getDirectory() + " is missing.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
      }
      // ME 190607 - BEGIN
      // Read forces from aflow qmvasp file - much faster
      string file = xinp.getDirectory() + "/" + DEFAULT_AFLOW_QMVASP_OUT;
      if (aurostd::EFileExist(file)) {
        xQMVASP qmvasp(file);
        xinp.getXStr().qm_forces = qmvasp.vforces;
      }
      if (xinp.getXStr().qm_forces.size() != natoms) {
        xinp.getXStr().qm_forces.clear();
        file = xinp.getDirectory() + string("/vasprun.xml.static");
        if(!aurostd::EFileExist(file)) {
          file = xinp.getDirectory() + string("/vasprun.xml");
          if(!aurostd::EFileExist(file)) {
            string message = "The force file in " + xinp.getDirectory() + " directory is missing.";
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
          }
        }
        //xVASPRUNXML vasprunxml(file); OBSOLETE ME 190204 - far too slow
        xVASPRUNXML vasprunxml;
        vasprunxml.GetForcesFile(file);
        for (uint i = 0; i < vasprunxml.vforces.size(); i++) xinp.getXStr().qm_forces.push_back(vasprunxml.vforces[i]);
      }
      if (xinp.getXStr().qm_forces.size() == natoms) xinp.getXStr().qm_calculated = true;
    }
    if(xinp.AFLOW_MODE_AIMS){
      if(!aurostd::EFileExist(xinp.getDirectory() + string("/aims.out"))) {
        string message = "The aims.out file in " + xinp.getDirectory() + " is missing.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
      }
      xAIMSOUT xaimsout(xinp.getDirectory() + "/aims.out");
      for (uint i = 0; i < xaimsout.vforces.size(); i++) xinp.getXStr().qm_forces.push_back(xaimsout.vforces[i]);
      if (xinp.getXStr().qm_forces.size() == natoms) xinp.getXStr().qm_calculated = true;
    }
  }

  //subtractZeroStateForces///////////////////////////////////////////////////
  void subtractZeroStateForces(vector<_xinput>& xinps, bool contains_born) {
    uint ninps = xinps.size() - 1;
    if (contains_born) ninps--;
    uint natoms = xinps[ninps].getXStr().atoms.size();
    for (uint idxRun = 0; idxRun < ninps; idxRun++) {
      if (xinps[idxRun].getXStr().atoms.size() != natoms) {
        string function = "apl::subtractZeroStateForces()";
        string message = "Structure and ZEROSTATE structure do not have the same number of atoms.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_MISMATCH_);
      }
      for (uint k = 0; k < natoms; k++) {
        xinps[idxRun].getXStr().qm_forces[k](1) = xinps[idxRun].getXStr().qm_forces[k](1) - xinps[ninps].getXStr().qm_forces[k](1);
        xinps[idxRun].getXStr().qm_forces[k](2) = xinps[idxRun].getXStr().qm_forces[k](2) - xinps[ninps].getXStr().qm_forces[k](2);
        xinps[idxRun].getXStr().qm_forces[k](3) = xinps[idxRun].getXStr().qm_forces[k](3) - xinps[ninps].getXStr().qm_forces[k](3);
      }
    }
  }

  // ME190114
  // Cannot use const reference for zerostate because of getXStr()
  void subtractZeroStateForces(vector<_xinput>& xinps, _xinput& zerostate) {
    string function = "apl::subtractZeroStateForces()";
    stringstream _logger;
    uint natoms = zerostate.getXStr().atoms.size();
    if (!zerostate.getXStr().qm_calculated) {
      readForcesFromDirectory(zerostate);
      if (!zerostate.getXStr().qm_calculated) {
        string message = "The force file in " + zerostate.getDirectory() + " is wrong.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
      }
    }
    for (uint idxRun = 0; idxRun < xinps.size(); idxRun++) {
      if (xinps[idxRun].getXStr().atoms.size() != natoms) {
        string function = "apl::subtractZeroStateForces()";
        string message = "Structure and ZEROSTATE structure do not have the same number of atoms.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_MISMATCH_);
      }
      for (uint at = 0; at < natoms; at++) {
        xinps[idxRun].getXStr().qm_forces[at](1) = xinps[idxRun].getXStr().qm_forces[at](1) - zerostate.getXStr().qm_forces[at](1);
        xinps[idxRun].getXStr().qm_forces[at](2) = xinps[idxRun].getXStr().qm_forces[at](2) - zerostate.getXStr().qm_forces[at](2);
        xinps[idxRun].getXStr().qm_forces[at](3) = xinps[idxRun].getXStr().qm_forces[at](3) - zerostate.getXStr().qm_forces[at](3);
      }
    }
  }

}  // namespace apl

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
//****************************************************************************
