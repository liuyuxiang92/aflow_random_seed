#include "aflow_apl.h"

namespace apl {

// ///////////////////////////////////////////////////////////////////////////
LinearResponsePC::LinearResponsePC(Supercell& sc, vector<ClusterSet>& clst,
                                   _xinput& xinput, _aflags& aflags, _kflags& kflags,
                                   _xflags& xflags, //_vflags& vflags, 
                                   string& AflowIn, Logger& l)
    : PhononCalculator(sc, clst, xinput, aflags, kflags, xflags, AflowIn, l) {
  _isPolarMaterial = false;
}

// ///////////////////////////////////////////////////////////////////////////
LinearResponsePC::~LinearResponsePC() {
  clear();
}

// ///////////////////////////////////////////////////////////////////////////
void LinearResponsePC::clear() {
}

//////////////////////////////////////////////////////////////////////////////
// ME190412 - Added aapl_stagebreak
void LinearResponsePC::calculateForceFields(bool aapl_stagebreak) {
  // Check if supercell is already built
  if (!_supercell.isConstructed()) {
    throw APLRuntimeError("apl::LinearResponsePC::calculateForceFields(); The supercell structure has not been initialized yet.");
  }

  // Do calculation an additional calculations for polar materials
  bool halfStageBreak = false;
  if (_isPolarMaterial) {
    try {
      // if tarred and compressed directory exists...


      // COREY CHECK THIS
      deque<string> vext; aurostd::string2tokens(".bz2,.xz,.gz",vext,",");vext.push_front(""); // cheat for void string
      string tarfilename =_AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ + ".tar";
      for(uint iext=0;iext<vext.size();iext++) {
	if (aurostd::FileExist(tarfilename+vext.at(iext))) {
	  // We will extract only what we need OUTCAR and LOCK...
	  if(_kbinFlags.AFLOW_MODE_VASP)
	    aurostd::execute(string("tar -xf ") + tarfilename + vext.at(iext) + " --wildcards "+"/OUTCAR*"); 
	  if(_kbinFlags.AFLOW_MODE_AIMS)
	    aurostd::execute(string("tar -xf ") + tarfilename + vext.at(iext) + " --wildcards "+"/aims.out*");
	}
      }
      //                 " --wildcards "
      //                 " " +
      //                 _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ + "/OUTCAR*" +
      //                 " " + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ + "/" + _AFLOWLOCK_);  //CO
      
      _xinput xinpBE(_xInput);  // ME190113
      // Calc. Born effective charge tensors and dielectric constant matrix
      runVASPCalculationsBE(xinpBE);  // ME190113

      // Parse it from OUTCAR
      if (!aapl_stagebreak) {  // ME190412 - Skip if AAPL calculations still need to be run
      if(_kbinFlags.AFLOW_MODE_VASP){readBornEffectiveChargesFromOUTCAR(xinpBE);}  // ME190113
      if(_kbinFlags.AFLOW_MODE_AIMS){readBornEffectiveChargesFromAIMSOUT();}

      // Enforce ASR (Acoustic sum rules)
      symmetrizeBornEffectiveChargeTensors();

      // Parser epsilon from OUTCAR
      if(_kbinFlags.AFLOW_MODE_VASP){readDielectricTensorFromOUTCAR(xinpBE);}  // ME190113
      if(_kbinFlags.AFLOW_MODE_AIMS){readDielectricTensorFromAIMSOUT();}

      //
      _logger << "Dielectric tensor: ";
      for (int a = 1; a <= 3; a++)
        for (int b = 1; b <= 3; b++)
          _logger << sf("%5.3f") << _dielectricTensor(a, b) << " ";
      _logger << apl::endl;

      // precompute
      _inverseDielectricTensor = inverse(_dielectricTensor);
      _recsqrtDielectricTensorDeterminant = 1.0 / sqrt(determinant(_dielectricTensor));
      }

//  OBSOLETE - ME181024
//      // Pack/Remove the whole directory...
//      // COREY CHECK THIS
//      if (DOtar) {
//        if (!aurostd::FileExist(tarfilename)) {
//	  aurostd::execute(string("tar -cf ") + tarfilename + " " + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ + "/");
//	  aurostd::CompressFile(tarfilename,_kbinFlags.KZIP_BIN);
//	}
//        if (aurostd::FileExist(tarfilename)) aurostd::execute(string("rm -rf ") + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ + "/");
//      }
    } catch (APLLogicError& e) {
      _logger << apl::error << e.what() << apl::endl;
      _logger << warning << "Switching the dipole-dipole correction off." << apl::endl;
      _isPolarMaterial = false;
    } catch (APLStageBreak& e) {
      halfStageBreak = true;
    }
  }

  // if tarred and compressed directory exists...
  string tarfilename = string(_AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_) + ".tar.bz2";
  if (aurostd::FileExist(tarfilename)) {
    // We will extract only what we need DYNMAT and LOCK...
    aurostd::execute(string("tar -xf ") + tarfilename +
                     " --wildcards "
                     " " +
                     _AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_ + "/DYNMAT*" +
                     " " + _AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_ + "/" + _AFLOWLOCK_);  //CO
  }

  // Call VASP to calculate forces by LR
  _xinput xinpFF(_xInput);
  runVASPCalculationsFF(xinpFF);

  // ME190412 - do not read files when AAPL calculations still need to run
  // or when Born charges still need to be calculated.
  if (aapl_stagebreak || halfStageBreak) throw APLStageBreak();
  // Get forces from the DYNMAT file generated by VASP
  readForceFieldsFromDYNMAT(xinpFF);

//  OBSOLETE - ME181024
//  // Pack/Remove the whole directory...
//  if (DOtar)
//    if (!aurostd::FileExist(tarfilename)) aurostd::execute(string("tar -cjf ") + tarfilename + " " + _AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_ + "/");
//  if (DOtar)
//    if (aurostd::FileExist(tarfilename)) aurostd::execute(string("rm -rf ") + _AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_ + "/");

  // It can happen the forcefields are computed, but the born effective charges not...
  // if (halfStageBreak) throw APLStageBreak();  OBSOLETE ME190412 - moved before readForceFieldsFromDYNMAT

  // Print some information
  int dof = 0;
  for (uint i = 0; i < _uniqueDistortions.size(); i++)
    dof += _uniqueDistortions[i].size();
  _logger << "Found " << dof << " degree(s) of freedom." << apl::endl;
  for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
    int id = _supercell.getUniqueAtomID(i);
    for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
      _logger << "Atom [" << sf("%03d") << id << "] (" << sf("%f")
              << sw(2) << _supercell.getSupercellStructure().atoms[id].cleanname
              << ") has been distorted in direction ["
              << sf("%5.3f") << _uniqueDistortions[i][j](1) << ","
              << sf("%5.3f") << _uniqueDistortions[i][j](2) << ","
              << sf("%5.3f") << _uniqueDistortions[i][j](3) << "]." << apl::endl;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// Wrapper function to run the linear response (force fields) calculation.
void LinearResponsePC::runVASPCalculationsFF(_xinput& xinp) {  // ME190113
  runVASPCalculationsLRBE(xinp, false);

//[ME181226 OBSOLETE]  // Our run;
//[ME181226 OBSOLETE]  _xinput xInput(_xInput); //_xvasp vaspRun(_vaspRun);
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  // Create run id name
//[ME181226 OBSOLETE]  string runname = string(_AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_);
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  // Setup working directory
//[ME181226 OBSOLETE]  xInput.setDirectory( _xInput.getDirectory() + runname );
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  // Generate POSCAR
//[ME181226 OBSOLETE]  xInput.setXStr(_supercell.getSupercellStructureLight());
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  // If there is already a LOCK file, it means this directory was already generated
//[ME181226 OBSOLETE]  // and computed, hence do not touch and leave, but store this structure in the
//[ME181226 OBSOLETE]  // list, hence it will be used in next part of code.
//[ME181226 OBSOLETE]  if (aurostd::FileExist(xInput.getDirectory() + string("/") + _AFLOWLOCK_)) return;
//[ME181226 OBSOLETE]  // If not, continue in this way, prepare generation of _AFLOWLIN ...
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  _logger << "Creating " << xInput.getDirectory() << apl::endl;
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  if(xInput.AFLOW_MODE_VASP){
//[ME181226 OBSOLETE]  // Common KPOINTS settings with PRIORITIES
//[ME181226 OBSOLETE]    xInput.xvasp.AVASP_KSCHEME = _xFlags.vflags.KBIN_VASP_KPOINTS_KSCHEME.content_string;
//[ME181226 OBSOLETE]    if (_xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.isentry)
//[ME181226 OBSOLETE]      xInput.xvasp.AVASP_KSCHEME = _xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.content_string;
//[ME181226 OBSOLETE]    xInput.xvasp.AVASP_value_KPPRA = _xFlags.vflags.KBIN_VASP_KPOINTS_KPPRA.content_int;
//[ME181226 OBSOLETE]    if (_xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.isentry)
//[ME181226 OBSOLETE]      xInput.xvasp.AVASP_value_KPPRA = _xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.content_uint;
//[ME181226 OBSOLETE]    // [OBSOLETE]  xInput.xvasp.AVASP_KSCHEME = _xFlags.vflags.KBIN_VASP_KPOINTS_KSCHEME.content_string;
//[ME181226 OBSOLETE]    // [OBSOLETE]  xInput.xvasp.AVASP_value_KPPRA = _xFlags.vflags.KBIN_VASP_KPOINTS_KPPRA.content_int;
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  // Switch off autotune, because....
//[ME181226 OBSOLETE]  _kbinFlags.KBIN_MPI_AUTOTUNE = false;
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  // Clear old INCAR and set it as we want...
//[ME181226 OBSOLETE]    xInput.xvasp.INCAR.str(std::string());
//[ME181226 OBSOLETE]  string system;
//[ME181226 OBSOLETE]    for (uint i = 0; i < xInput.getXStr().species.size(); i++)
//[ME181226 OBSOLETE]      system = system + xInput.getXStr().species_pp.at(i);
//[ME181226 OBSOLETE]  system = system + "@" + runname;
//[ME181226 OBSOLETE]    xInput.xvasp.INCAR << "SYSTEM=" << system << std::endl;
//[ME181226 OBSOLETE]    xInput.xvasp.INCAR << std::endl;
//[ME181226 OBSOLETE]    xInput.xvasp.INCAR << "# Added by [AFLOW_APL] begin" << std::endl;
//[ME181226 OBSOLETE]    xInput.xvasp.INCAR << "NELMIN=4         # The forces have to be well converged" << std::endl;
//[ME181226 OBSOLETE]    xInput.xvasp.INCAR << "NELM = 120       # Many electronic steps (SC2013)" << std::endl;
//[ME181226 OBSOLETE]    xInput.xvasp.INCAR << "ADDGRID=.TRUE.   # For finer forces" << std::endl;
//[ME181226 OBSOLETE]    xInput.xvasp.INCAR << "IBRION=8         # Linear Response method (needs VASP 5.2 or higher)" << std::endl;
//[ME181226 OBSOLETE]    xInput.xvasp.INCAR << "# Added by [AFLOW_APL] end" << std::endl;
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  // For this we need VASP5
//[ME181226 OBSOLETE]    xInput.getXStr().is_vasp4_poscar_format = false;
//[ME181226 OBSOLETE]    xInput.getXStr().is_vasp5_poscar_format = true;
//[ME181226 OBSOLETE]  }
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  if(xInput.AFLOW_MODE_AIMS){
//[ME181226 OBSOLETE]    xInput.xaims.CONTROL.str(std::string());
//[ME181226 OBSOLETE]    KBIN::AIMS_Produce_CONTROL(xInput.xaims,_AflowIn,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);  //DEFAULT
//[ME181226 OBSOLETE]    KBIN::AIMS_Modify_CONTROL(xInput.xaims,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);            //DEFAULT
//[ME181226 OBSOLETE]  }
//[ME181226 OBSOLETE]
//[ME181226 OBSOLETE]  // Create _AFLOWIN_
//[ME181226 OBSOLETE]  writeOUTPUT(xInput);

}

//////////////////////////////////////////////////////////////////////////////
// We will use VASP5.2+ to calculate Born effective charge tensors and
// dielectric constant matrix in the primitive cell with very high precision
// Both values are needed by non-analytical term of dynamic matrix for
// correct TO-LO splitting of optical phonon branches of polar systems

//////////////////////////////////////////////////////////////////////////////
// Wrapper function to calculate the Born effective charge tensor and the
// dilectric tensor.
void PhononCalculator::runVASPCalculationsBE(_xinput& xinp, uint ncalcs) { // ME190113
  runVASPCalculationsLRBE(xinp, true, ncalcs);

//[ME18226 OBSOLETE]  // Our run;
//[ME18226 OBSOLETE]  _xinput xInput(_xInput); //_xvasp vaspRun(_vaspRun);
//[ME18226 OBSOLETE]
//[ME18226 OBSOLETE]  // Create run id name
//[ME18226 OBSOLETE]  string runname = string(_AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_);
//[ME18226 OBSOLETE]
//[ME18226 OBSOLETE]  // Setup working directory
//[ME18226 OBSOLETE]  xInput.setDirectory( _xInput.getDirectory() + runname );
//[ME18226 OBSOLETE]
//[ME18226 OBSOLETE]  // Generate POSCAR
//[ME18226 OBSOLETE]  xInput.setXStr(_supercell.getInputStructure());
//[ME18226 OBSOLETE]
//[ME18226 OBSOLETE]  // If there is already a LOCK file, it means this directory was already generated
//[ME18226 OBSOLETE]  // and computed, hence do not touch and leave, but store this structure in the
//[ME18226 OBSOLETE]  // list, hence it will be used in next part of code.
//[ME18226 OBSOLETE]  if (aurostd::FileExist(xInput.getDirectory() + string("/") + _AFLOWLOCK_)) return;  //CO
//[ME18226 OBSOLETE]                                                                                  // If not, continue in this way, prepare generation of _AFLOWLIN ...
//[ME18226 OBSOLETE]  _logger << "Creating " << xInput.getDirectory() << apl::endl;                       //CO
//[ME18226 OBSOLETE]
//[ME18226 OBSOLETE]  if(xInput.AFLOW_MODE_VASP){
//[ME18226 OBSOLETE]  // Common KPOINTS settings with PRIORITIES
//[ME18226 OBSOLETE]    xInput.xvasp.AVASP_KSCHEME = _xFlags.vflags.KBIN_VASP_KPOINTS_KSCHEME.content_string;
//[ME18226 OBSOLETE]    if (_xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.isentry) {
//[ME18226 OBSOLETE]      xInput.xvasp.AVASP_KSCHEME = _xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.content_string;
//[ME18226 OBSOLETE]  }
//[ME18226 OBSOLETE]    //xInput.xvasp.AVASP_value_KPPRA = _xFlags.vflags.KBIN_VASP_KPOINTS_KPPRA.content_int;
//[ME18226 OBSOLETE]    xInput.xvasp.AVASP_value_KPPRA = 10000;  // JJPR Modification: Dielectric tensor and Born charges require a very dense q-point mesh.
//[ME18226 OBSOLETE]    if (_xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.isentry) {
//[ME18226 OBSOLETE]      xInput.xvasp.AVASP_value_KPPRA = _xFlags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.content_uint;
//[ME18226 OBSOLETE]  }
//[ME18226 OBSOLETE]  // [OBSOLETE]  xInput.xvasp.AVASP_KSCHEME = _xFlags.vflags.KBIN_VASP_KPOINTS_KSCHEME.content_string;
//[ME18226 OBSOLETE]  // [OBSOLETE]  xInput.xvasp.AVASP_value_KPPRA = _xFlags.vflags.KBIN_VASP_KPOINTS_KPPRA.content_int;
//[ME18226 OBSOLETE]
//[ME18226 OBSOLETE]#ifdef AFLOW_APL_VASP_USE_LCALCEPS
//[ME18226 OBSOLETE]    xInput.xvasp.AVASP_KSCHEME = "GAMMA    // forced by  AFLOW_APL_VASP_USE_LCALCEPS";
//[ME18226 OBSOLETE]#endif
//[ME18226 OBSOLETE]  // Switch off autotune, because....
//[ME18226 OBSOLETE]  _kbinFlags.KBIN_MPI_AUTOTUNE = false;
//[ME18226 OBSOLETE]
//[ME18226 OBSOLETE]  // Clear old INCAR and set it as we want...
//[ME18226 OBSOLETE]    xInput.xvasp.INCAR.str(std::string());
//[ME18226 OBSOLETE]  string system;
//[ME18226 OBSOLETE]    for (uint i = 0; i < xInput.getXStr().species.size(); i++)
//[ME18226 OBSOLETE]      system = system + xInput.getXStr().species_pp.at(i);
//[ME18226 OBSOLETE]  system = system + "@" + runname;
//[ME18226 OBSOLETE]    xInput.xvasp.INCAR << "SYSTEM=" << system << std::endl;
//[ME18226 OBSOLETE]    //    xInput.xvasp.INCAR << std::endl;
//[ME18226 OBSOLETE]    xInput.xvasp.INCAR << "# Added by [AFLOW_APL] begin" << std::endl;
//[ME18226 OBSOLETE]    xInput.xvasp.INCAR << "ADDGRID=.TRUE.   # For finer forces" << std::endl;
//[ME18226 OBSOLETE]#ifdef AFLOW_APL_VASP_USE_LEPSILON
//[ME18226 OBSOLETE]    xInput.xvasp.INCAR << "IBRION=8         # Linear Response method (needs VASP 5.2 or higher)" << std::endl;
//[ME18226 OBSOLETE]    xInput.xvasp.INCAR << "LEPSILON=.TRUE.  # Calculate Born effective charges and dielectric constants" << std::endl;
//[ME18226 OBSOLETE]#endif
//[ME18226 OBSOLETE]#ifdef AFLOW_APL_VASP_USE_LCALCEPS
//[ME18226 OBSOLETE]    xInput.xvasp.INCAR << "LCALCEPS=.TRUE.  # Calculate Born effective charges and dielectric constants" << std::endl;
//[ME18226 OBSOLETE]    xInput.xvasp.INCAR << "### not necessary EFIELD_PEAD= 0.001 0.001 0.001   # ten times smaller than the default" << std::endl;
//[ME18226 OBSOLETE]#endif
//[ME18226 OBSOLETE]    xInput.xvasp.INCAR << "# Added by [AFLOW_APL] end" << std::endl;
//[ME18226 OBSOLETE]
//[ME18226 OBSOLETE]  // For this we need VASP5
//[ME18226 OBSOLETE]    xInput.getXStr().is_vasp4_poscar_format = false;
//[ME18226 OBSOLETE]    xInput.getXStr().is_vasp5_poscar_format = true;
//[ME18226 OBSOLETE]  }
//[ME18226 OBSOLETE]
//[ME18226 OBSOLETE]  if(xInput.AFLOW_MODE_AIMS){
//[ME18226 OBSOLETE]    xInput.xaims.CONTROL.str(std::string());
//[ME18226 OBSOLETE]    KBIN::AIMS_Produce_CONTROL(xInput.xaims,_AflowIn,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);  //DEFAULT
//[ME18226 OBSOLETE]    KBIN::AIMS_Modify_CONTROL(xInput.xaims,_logger.getOutputStream(),_aflowFlags,_kbinFlags,_xFlags.aimsflags);            //DEFAULT
//[ME18226 OBSOLETE]  }
//[ME18226 OBSOLETE]
//[ME18226 OBSOLETE]  // Create _AFLOWIN_
//[ME18226 OBSOLETE]  writeOUTPUT(xInput);

}

//////////////////////////////////////////////////////////////////////////////
// A unified function to handle Born effective charge and force fields
// calculations. Both methods require similar treatment, so it is cleaner to
// use only one function for both methods.
void PhononCalculator::runVASPCalculationsLRBE(_xinput& xInput, bool born, uint ncalcs) { // ME190112

  //_xinput xInput(_xInput); OBSOLETE ME190113
  if (born) {
  xInput.setXStr(_supercell.getInputStructure());
  } else {
  xInput.setXStr(_supercell.getSupercellStructureLight());
  }

  // ME 190108 - Added title
  xInput.getXStr().title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xInput.getXStr().title);
  if(xInput.getXStr().title.empty()){xInput.getXStr().buildGenericTitle(true,false);}
  if (born) {
    xInput.getXStr().title+=" dielectric tensor";
  } else {
    xInput.getXStr().title+=" linear response";
  }

  // For VASP, use the standardized aflow.in creator
  if(xInput.AFLOW_MODE_VASP){
    if (born) {
      // ME190112 - add calculation index
      if (ncalcs == 0) {  // Linear response method
        xInput.xvasp.AVASP_arun_runname = "2_" + _AFLOW_APL_BORN_EPSILON_RUNNAME_;
      } else {  // Direct method
        xInput.xvasp.AVASP_arun_runname = stringify(ncalcs+1) + "_" + _AFLOW_APL_BORN_EPSILON_RUNNAME_;
      }
      xInput.xvasp.aopts.flag("FLAG::AVASP_LR", false);
      xInput.xvasp.aopts.flag("FLAG::AVASP_BORN", true);
    } else {
      xInput.xvasp.AVASP_arun_runname = "1_" + _AFLOW_APL_FORCEFIELDS_RUNNAME_;
      xInput.xvasp.aopts.flag("FLAG::AVASP_BORN", false);
      xInput.xvasp.aopts.flag("FLAG::AVASP_LR", true);
  }
    // Switch off autotune
  _kbinFlags.KBIN_MPI_AUTOTUNE = false;

    // Set POSCAR to VASP5 format
    xInput.getXStr().is_vasp4_poscar_format = false;
    xInput.getXStr().is_vasp5_poscar_format = true;
    createAflowInPhonons(xInput);
  }
  // For AIMS, use the old method until we have AVASP_populateXAIMS
  if(xInput.AFLOW_MODE_AIMS){
    string runname;
    if (born) {
      runname = _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_;
    } else {
      runname = _AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_;
  }
    xInput.setDirectory( _xInput.getDirectory() + "/" + runname );
    if (!filesExistPhonons(xInput)) {
  _logger << "Creating " << xInput.getDirectory() << apl::endl;
      createAflowInPhonons(xInput, runname);
  }
  }
}

//////////////////////////////////////////////////////////////////////////////
void PhononCalculator::readBornEffectiveChargesFromAIMSOUT(void) {
  string directory = string("./") + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_;
  string infilename = directory + string("/aims.out");

  if (!aurostd::EFileExist(infilename, infilename)) {
    //_logger << apl::warning << "The aims.out file in " << directory << " directory is missing." << apl::endl;
    //throw APLLogicError("apl::LinearResponsePC::readBornEffectiveChargesFromAIMSOUT(); Missing data from one job.");
    throw APLStageBreak();
  }

  throw APLLogicError("apl::LinearResponsePC::readBornEffectiveChargesFromAIMSOUT(); This functionality has yet to be implemented.");
}
//////////////////////////////////////////////////////////////////////////////
void PhononCalculator::readBornEffectiveChargesFromOUTCAR(const _xinput& xinp) {  // ME190113
  //string directory = string("./") + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_; OBSOLETE ME190113
  string directory = xinp.xvasp.Directory;  // ME190113

  //if (!aurostd::FileExist(directory + string("/") + _AFLOWLOCK_))  //CO
  //  throw APLStageBreak();

  //CO - START
  string infilename = directory + string("/OUTCAR.static");

  if (!aurostd::EFileExist(infilename, infilename)) {
    infilename = directory + string("/OUTCAR");
    if (!aurostd::EFileExist(infilename, infilename)) {
      //_logger << apl::warning << "The OUTCAR file in " << directory << " directory is missing." << apl::endl;
      //throw APLLogicError("apl::LinearResponsePC::readBornEffectiveChargesFromOUTCAR(); Missing data from one job.");
      throw APLStageBreak();
    }
  }

  // Open our file
  //CO - START
  vector<string> vlines;
  aurostd::efile2vectorstring(infilename, vlines);
  if (!vlines.size()) {
    //CO - END
    throw apl::APLLogicError("apl::LinearResponsePC::readBornEffectiveChargesFromOUTCAR(); Cannot open input OUTCAR.static file.");
  }

  string line;
  uint line_count = 0;  //CO
  string KEY; //ME181226
  if (DEFAULT_APL_USE_LEPSILON) { //ME181226
    KEY = string("BORN EFFECTIVE CHARGES (in e, cummulative output)");//ME181226
  } else { //ME181226
    KEY = string("BORN EFFECTIVE CHARGES (including local field effects)"); //ME181226
  } //ME181226

  while (true) {
    // Get line
    //CO - START
    //getline(infile,line);
    //if( infile.eof() ) {
    if (line_count == vlines.size()) {
      //CO - END
      throw apl::APLLogicError("apl::LinearResponsePC::readBornEffectiveChargesFromOUTCAR(); No information about Born effective charges in OUTCAR.");
    }
    line = vlines[line_count++];  //CO

    // Check for our key line
    if (line.size() < KEY.size()) continue;
    if (line.find(KEY) != string::npos) break;
  }
  // Read in all ...
  xmatrix<double> m(3, 3);
  vector<string> tokens;
  //CO - START
  //getline(infile,line); // Skip line "----------------...."
  line = vlines[line_count++];
  for (uint i = 0; i < _supercell.getInputStructure().atoms.size(); i++) {
    // Get atom ID but not use it...
    //getline(infile,line);
    line = vlines[line_count++];
    //tokenize(line,tokens,string(" "));
    //auto int id = aurostd::string2utype<int>(tokens.at(1)) - 1;
    //if( id != _supercell.getInputStructure().iatoms[i][0] ) continue;
    //tokens.clear();

    // Get its charge tensor
    for (int j = 1; j <= 3; j++) {
      //getline(infile,line);
      line = vlines[line_count++];
      tokenize(line, tokens, string(" "));
      m(j, 1) = aurostd::string2utype<double>(tokens.at(1));
      m(j, 2) = aurostd::string2utype<double>(tokens.at(2));
      m(j, 3) = aurostd::string2utype<double>(tokens.at(3));
      tokens.clear();
    }

    // Store it
    _bornEffectiveChargeTensor.push_back(m);
  }

  // Clear used stuff
  //infile.clear();
  //infile.close();
  //CO - END
}

//////////////////////////////////////////////////////////////////////////////
void PhononCalculator::symmetrizeBornEffectiveChargeTensors(void) {
  //CO - START
  // Test of stupidity...
  if (_supercell.getEPS() == AUROSTD_NAN) {
    throw APLRuntimeError("apl::PhononCalculator::buildForceConstantMatrices(); Need to define symmetry tolerance.");
  }
  //CO - END
  // Show charges
  _logger << "Input born effective charge tensors (for primitive cell):" << apl::endl;
  for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++) {
    int id = i;
    _logger << "Atom [" << sf("%03d") << id << "] (" << sf("%f")
            << sw(2) << _supercell.getInputStructure().atoms[id].cleanname
            << ") Born effective charge = ";
    for (int a = 1; a <= 3; a++)
      for (int b = 1; b <= 3; b++)
        _logger << sf("%+5.3f") << _bornEffectiveChargeTensor[i](a, b) << " ";
    _logger << apl::endl;
  }
  _logger << sf("%f");

  //_logger << "Symmetrization of born effective charge tensors (1)." << apl::endl;

  // Step1
  for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
    // Get the number of this atom in the whole list
    int basedUniqueAtomID = _supercell.getUniqueAtomID(i);

    xmatrix<double> sum(3, 3);
    for (int j = 0; j < _supercell.getNumberOfEquivalentAtomsOfType(i); j++) { //CO190218
      try {  //CO
        const _sym_op& symOp = _supercell.getSymOpWhichMatchAtoms(_supercell.getUniqueAtomID(i, j), basedUniqueAtomID, _FGROUP_);
        sum += inverse(symOp.Uc) * _bornEffectiveChargeTensor[_supercell.sc2pcMap(_supercell.getUniqueAtomID(i, j))] * symOp.Uc;
      }
      //CO - START
      catch (APLLogicError& e) {
        _logger << error << "Mapping problem " << _supercell.getUniqueAtomID(i, j) << " <-> " << basedUniqueAtomID << "?" << apl::endl;
        throw APLLogicError("apl::PhononCalculator::symmetrizeBornEffectiveChargeTensors(); Mapping failed.");
      }
      //CO - END
    }

    sum = (1.0 / _supercell.getNumberOfEquivalentAtomsOfType(i)) * sum; //CO190218

    for (int j = 0; j < _supercell.getNumberOfEquivalentAtomsOfType(i); j++) //CO190218
      _bornEffectiveChargeTensor[_supercell.sc2pcMap(_supercell.getUniqueAtomID(i, j))] = sum;
  }

  //_logger << "Symmetrization of born effective charge tensors (2)." << apl::endl;

  // Step2
  vector<xmatrix<double> > newbe = _bornEffectiveChargeTensor;
  const vector<vector<_sym_op> >& agroup = _supercell.getAGROUP();  //CO
  for (int i = 0; i < _supercell.getNumberOfAtoms(); i++) {
    // Translate the center to this atom
    _supercell.center(i);

    xmatrix<double> sum(3, 3);
    for (uint symOpID = 0; symOpID < agroup[i].size(); symOpID++) {
      const _sym_op& symOp = agroup[i][symOpID];
      sum = sum + (inverse(symOp.Uc) * _bornEffectiveChargeTensor[_supercell.sc2pcMap(i)] * symOp.Uc);
    }
    newbe[_supercell.sc2pcMap(i)] = (1.0 / agroup[i].size()) * sum;
    // Translate the center back
    //_supercell.center_original(); //CO
  }
  // Translate the center back
  //_supercell.center(0);
  _supercell.center_original();  //CO

  _bornEffectiveChargeTensor.clear();
  _bornEffectiveChargeTensor = newbe;
  newbe.clear();

  // Step 3
  _logger << "Forcing the acoustic sum rule (ASR). Resulted born effective charges (for supercell):" << apl::endl;

  xmatrix<double> sum(3, 3);
  for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++)
    sum += _bornEffectiveChargeTensor[i];
  sum = (1.0 / _bornEffectiveChargeTensor.size()) * sum;
  for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++)
    _bornEffectiveChargeTensor[i] -= sum;

  // Make list only for unique atoms
  for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++)
    newbe.push_back(_bornEffectiveChargeTensor[_supercell.sc2pcMap(_supercell.getUniqueAtomID(i))]);
  _bornEffectiveChargeTensor.clear();
  _bornEffectiveChargeTensor = newbe;
  newbe.clear();

  // Show charges
  for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
    int id = _supercell.getUniqueAtomID(i);
    _logger << "Atom [" << sf("%03d") << id << "] (" << sf("%f")
            << sw(2) << _supercell.getSupercellStructure().atoms[id].cleanname
            << ") Born effective charge = ";
    for (int a = 1; a <= 3; a++)
      for (int b = 1; b <= 3; b++)
        _logger << sf("%+5.3f") << _bornEffectiveChargeTensor[i](a, b) << " ";
    _logger << apl::endl;
  }
  _logger << sf("%f");
}

//////////////////////////////////////////////////////////////////////////////
void PhononCalculator::readDielectricTensorFromAIMSOUT(void) {
  string directory = string("./") + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_;

  //CO - START
  string infilename = directory + string("/aims.out");
  if (!aurostd::EFileExist(infilename, infilename)) {
    //_logger << apl::warning << "The OUTCAR file in " << directory << " directory is missing." << apl::endl;
    //throw APLLogicError("apl::LinearResponsePC::readDielectricTensorFromOUTCAR(); Missing data from one job.");
    throw APLStageBreak();
  }
  throw APLLogicError("apl::LinearResponsePC::readDielectricTensorFromAIMSOUT(); This functionality has yet to be implemented.");
}
//////////////////////////////////////////////////////////////////////////////
void PhononCalculator::readDielectricTensorFromOUTCAR(const _xinput& xinp) {  // ME190113
  //string directory = string("./") + _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_; OBSOLETE ME190113
  string directory = xinp.xvasp.Directory;  // ME190113

  //if (!aurostd::FileExist(directory + string("/") + _AFLOWLOCK_))  //CO
  //  throw APLStageBreak();

  //CO - START
  string infilename = directory + string("/OUTCAR.static");
  if (!aurostd::EFileExist(infilename, infilename)) {
    infilename = directory + string("/OUTCAR");
    if (!aurostd::EFileExist(infilename, infilename)) {
      //_logger << apl::warning << "The OUTCAR file in " << directory << " directory is missing." << apl::endl;
      //throw APLLogicError("apl::LinearResponsePC::readDielectricTensorFromOUTCAR(); Missing data from one job.");
      throw APLStageBreak();
    }
  }

  // Open our file
  vector<string> vlines;
  uint line_count = 0;
  string line;
  aurostd::efile2vectorstring(infilename, vlines);
  if (!vlines.size()) {
    throw apl::APLLogicError("LinearResponsePC::readDielectricTensorFromOUTCAR(); Cannot open input OUTCAR file.");
  }
//CO - END

// Find
  string KEY; //ME181226
  if (DEFAULT_APL_USE_LEPSILON) { //ME181226
    KEY = string("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)"); //ME181226
  } else { //ME181226
    KEY = string("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects)"); //ME181226
  }//ME181226

  while (true) {
    // Get line
    //CO - START
    //getline(infile,line);
    //if( infile.eof() )
    if (line_count == vlines.size()) {
      throw apl::APLLogicError("LinearResponsePC::readDielectricTensorFromOUTCAR(); No information about dielectric tensor in OUTCAR.");
    }
    line = vlines[line_count++];
    //CO - END

    // Check for our key line
    if (line.size() < KEY.size()) continue;
    if (line.find(KEY) != string::npos) break;
  }

  // Read in all ...
  //CO - START
  //getline(infile,line); // Skip line "----------------...."
  line = vlines[line_count++];

  // Get it
  vector<string> tokens;
  for (int j = 1; j <= 3; j++) {
    //getline(infile,line);
    line = vlines[line_count++];
    tokenize(line, tokens, string(" "));
    _dielectricTensor(j, 1) = aurostd::string2utype<double>(tokens.at(0));
    _dielectricTensor(j, 2) = aurostd::string2utype<double>(tokens.at(1));
    _dielectricTensor(j, 3) = aurostd::string2utype<double>(tokens.at(2));
    tokens.clear();
  }

  //???MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION
  /*
      const string KEY2 = string("MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION");
      while( true ) {
      // Get line
      //getline(infile,line);
      line=vlines[line_count++];
      //if( infile.eof() )
      if( line_count==vlines.size() ) {
      {
      throw apl::APLLogicError("LinearResponsePC::readDielectricTensorFromOUTCAR(); No information about dielectric tensor in OUTCAR.");
      }

      // Check for our key line
      if( line.size() < KEY2.size() ) continue;
      if( line.find(KEY2) != string::npos ) break;
      }

      // Read in all ...
      //getline(infile,line); // Skip line "----------------...."
      line=vlines[line_count++];

      // Get it
      for(int j = 1; j <= 3; j++) {
      //getline(infile,line);
      line=vlines[line_count++];
      tokenize(line,tokens,string(" "));
      _dielectricTensor(j,1) -= aurostd::string2utype<double>(tokens.at(0));
      _dielectricTensor(j,2) -= aurostd::string2utype<double>(tokens.at(1));
      _dielectricTensor(j,3) -= aurostd::string2utype<double>(tokens.at(2));
      tokens.clear();
      }
    */
  // Clear used stuff
  //infile.clear();
  //infile.close();
  //CO - END

  // Symmetrize
  /*
      const vector<_sym_op>& pgroup = _supercell.getSupercellStructure().pgroup;
      xmatrix<double> sum(3,3);
      for(uint symOpID = 0; symOpID < pgroup.size(); symOpID++) {
      const _sym_op& symOp = pgroup[symOpID];
      sum = sum + ( inverse(symOp.Uc) * _dielectricTensor * symOp.Uc );
      }
      _dielectricTensor = ( 1.0 / pgroup.size() ) * sum;
    */
}

//////////////////////////////////////////////////////////////////////////////
void LinearResponsePC::readForceFieldsFromDYNMAT(const _xinput& xinp) { // ME190113
  // string directory = string("./") + _AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_; ME190113
  string directory = xinp.xvasp.Directory; // ME190113

  //if (!aurostd::FileExist(directory + string("/") + _AFLOWLOCK_))  //CO
  //  throw APLStageBreak();

  string infilename = directory + string("/DYNMAT");
  if (!aurostd::EFileExist(infilename, infilename)) {  //CO
    infilename = directory + string("/DYNMAT.static");
    if (!aurostd::EFileExist(infilename, infilename))  //CO
    {
      throw APLStageBreak(); //ME181226
//      _logger << apl::warning << "The DYNMAT file in " << directory << " directory is missing." << apl::endl;  OBSOLETE - ME181024
//      throw APLLogicError("apl::LinearResponsePC::readForceFieldsFromDYNMAT(); Missing data from one job.");  OBSOLETE - ME181024
    }
  }

  // Clear old lists
  clear();

  // Open our file
  //CO - START
  vector<string> vlines;
  aurostd::efile2vectorstring(infilename, vlines);
  //ifstream infile(infilename.c_str());
  //if( !infile.is_open() ) {
  if (!vlines.size()) {
    throw apl::APLRuntimeError("LinearResponsePC::readForceFieldsFromDYNMAT(); Cannot open input DYNMAT file.");
  }

  // Header
  vector<int> values;
  uint line_count = 0;
  aurostd::string2tokens<int>(vlines[line_count++], values);
  //int ntypes=values[0]; //not used
  //int ntypes;
  //infile >> ntypes;

  int natoms = values[1];
  int nfree = values[2];
  //int natoms;
  //infile >> natoms;
  //int nfree;
  //infile >> nfree;
  //for(int i = 0; i < ntypes; i++) {
  //  double dtemp;
  //  infile >> dtemp;
  //}
  line_count++;
  //CO - END

  // Get forces
  vector<xvector<double> > distortionsOfOneAtom;
  vector<vector<xvector<double> > > forceFieldsOfOneAtom;
  int old_atom_id = 0;
  aurostd::string2tokens<int>(vlines[line_count++], values);  //CO
  for (int i = 0; i < nfree; i++) {
    int atom_id = values[0];  //CO
    //int atom_id;  //CO
    //infile >> atom_id;  //CO
    if (i == 0) old_atom_id = atom_id;

    //if(infile.eof()) break; //CO
    if (line_count == vlines.size()) break;  //CO

    // Distortion ID
    //CO - START
    //int idir=values[1]; //not used
    //int idir;
    //infile >> idir;

    // Distortion vector
    xvector<double> distortion(3);
    distortion(1) = values[2];
    distortion(2) = values[3];
    distortion(3) = values[4];
    //infile >> distortion(1);
    //infile >> distortion(2);
    //infile >> distortion(3);
    //CO - END
    double distortionLength = aurostd::modulus(distortion);
    DISTORTION_MAGNITUDE = distortionLength;
    distortion(1) = distortion(1) / distortionLength;
    distortion(2) = distortion(2) / distortionLength;
    distortion(3) = distortion(3) / distortionLength;

    // Get forces and calculate drift
    vector<xvector<double> > qm_forces;
    xvector<double> force(3);
    xvector<double> drift(3);
    for (int k = 0; k < natoms; k++) {
      //CO - START
      aurostd::string2tokens(vlines[line_count++], values);
      force(1) = values[0];
      force(2) = values[1];
      force(3) = values[2];
      //infile >> force(1);
      //infile >> force(2);
      //infile >> force(3);
      //CO - END
      qm_forces.push_back(force);
      drift = drift + force;
    }
    drift(1) = drift(1) / natoms;
    drift(2) = drift(2) / natoms;
    drift(3) = drift(3) / natoms;
    for (int k = 0; k < natoms; k++)
      qm_forces[k] = qm_forces[k] - drift;

    // Finish this atoms list and stasrt new one
    if (atom_id != old_atom_id) {
      _uniqueForces.push_back(forceFieldsOfOneAtom);
      _uniqueDistortions.push_back(distortionsOfOneAtom);
      old_atom_id = atom_id;
      for (uint i = 0; i < distortionsOfOneAtom.size(); i++)
        distortionsOfOneAtom[i].clear();
      distortionsOfOneAtom.clear();
      for (uint i = 0; i < forceFieldsOfOneAtom.size(); i++)
        forceFieldsOfOneAtom[i].clear();
      forceFieldsOfOneAtom.clear();
    }

    // Store distortions and forces
    distortionsOfOneAtom.push_back(distortion);
    forceFieldsOfOneAtom.push_back(qm_forces);

    // Final store ....
    if (i == nfree - 1) {
      _uniqueForces.push_back(forceFieldsOfOneAtom);
      _uniqueDistortions.push_back(distortionsOfOneAtom);
      for (uint i = 0; i < distortionsOfOneAtom.size(); i++)
        distortionsOfOneAtom[i].clear();
      distortionsOfOneAtom.clear();
      for (uint i = 0; i < forceFieldsOfOneAtom.size(); i++)
        forceFieldsOfOneAtom[i].clear();
      forceFieldsOfOneAtom.clear();
    }
  }

  // Clear used stuff
  //infile.clear(); //CO
  //infile.close(); //CO
}

//////////////////////////////////////////////////////////////////////////////

}  // namespace apl
