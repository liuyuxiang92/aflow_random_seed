// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************

#include "aflow_apl.h"

static const string _APL_LRPC_MODULE_ = "APL";  // for the logger

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         CONSTRUCTORS/DESTRUCTORS                         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  LinearResponsePC::LinearResponsePC(ostream& oss) : ForceConstantCalculator(oss) {
    free();
  }

  LinearResponsePC::LinearResponsePC(Supercell& sc, ofstream& mf, ostream& oss)
    : ForceConstantCalculator(sc, mf, oss) {
      free();
    }

  LinearResponsePC::LinearResponsePC(const LinearResponsePC& that)
    : ForceConstantCalculator(*that._supercell, *that.getOFStream(), *that.getOSS()) {
    free();
    copy(that);
  }

  LinearResponsePC& LinearResponsePC::operator=(const LinearResponsePC& that) {
    if (this != &that) {
      free();
      copy(that);
    }
    return *this;
  }

  LinearResponsePC::~LinearResponsePC() {
    xStream::free();
    free();
  }

  void LinearResponsePC::clear(Supercell& sc) {
    free();
    _supercell = &sc;
    _sc_set = true;
  }

  void LinearResponsePC::copy(const LinearResponsePC& that) {
    xStream::copy(that);
    _bornEffectiveChargeTensor = that._bornEffectiveChargeTensor;
    _dielectricTensor = that._dielectricTensor;
    _directory = that._directory;
    _forceConstantMatrices = that._forceConstantMatrices;
    _isPolarMaterial = that._isPolarMaterial;
    _supercell = that._supercell;
    _sc_set = that._sc_set;
    xInputs = that.xInputs;
  }

  void LinearResponsePC::free() {
    xInputs.clear();
    _bornEffectiveChargeTensor.clear();
    _dielectricTensor.clear();
    _directory = "";
    _forceConstantMatrices.clear();
    _isPolarMaterial = false;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            VASP CALCULATIONS                             //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  bool LinearResponsePC::runVASPCalculations(_xinput& xInput, _aflags& _aflowFlags,
      _kflags& _kbinFlags, _xflags& _xFlags, string& _AflowIn) {
    if (!_sc_set) {
      string function = "apl::LinearResponsePC::runVASPCalculations():";
      string message = "Supercell pointer not set.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }
    bool stagebreak = false;

    // Call VASP to calculate force constants using DFPT
    xInput.xvasp.AVASP_arun_mode = "APL";
    xInputs.clear();
    xInputs.push_back(xInput);
    stagebreak = runVASPCalculationsDFPT(xInputs[0], _aflowFlags, _kbinFlags, _xFlags, _AflowIn);
    if (_isPolarMaterial) {
      xInputs.push_back(xInput);
      stagebreak = (runVASPCalculationsBE(xInputs[1], _aflowFlags, _kbinFlags, _xFlags, _AflowIn, 2) || stagebreak);
    }
    return stagebreak;
  }

  //////////////////////////////////////////////////////////////////////////////
  // We will use VASP5.2+ to calculate Born effective charge tensors and
  // dielectric constant matrix in the primitive cell with very high precision
  // Both values are needed by the non-analytical term of the dynamical matrix
  // to capture the TO-LO splitting of optical phonon branches of polar systems.
  bool LinearResponsePC::runVASPCalculationsDFPT(_xinput& xInput, _aflags& _aflowFlags,
      _kflags& _kbinFlags, _xflags& _xFlags, string& _AflowIn) {
    bool stagebreak = false;

    xInput.setXStr(_supercell->getSupercellStructureLight());
    xInput.getXStr().title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xInput.getXStr().title);
    if(xInput.getXStr().title.empty()){xInput.getXStr().buildGenericTitle(true,false);}
    xInput.getXStr().title+=" linear response";

    // For VASP, use the standardized aflow.in creator
    if(xInput.AFLOW_MODE_VASP) {
      xInput.xvasp.AVASP_arun_runname = "1_" + _AFLOW_APL_DFPT_RUNNAME_;  // ME20200213
      xInput.xvasp.aopts.flag("APL_FLAG::AVASP_BORN", false);
      xInput.xvasp.aopts.flag("APL_FLAG::AVASP_LR", true);

      // Set POSCAR to VASP5 format
      xInput.getXStr().is_vasp4_poscar_format = false;
      xInput.getXStr().is_vasp5_poscar_format = true;
      stagebreak = (createAflowInPhonons(_aflowFlags, _kbinFlags, _xFlags, xInput) || stagebreak);
    }
    // For AIMS, use the old method until we have AVASP_populateXAIMS
    if(xInput.AFLOW_MODE_AIMS) {
      string runname = _AFLOW_APL_DFPT_DIRECTORY_NAME_;  // ME20200213
      xInput.setDirectory(_directory + "/" + runname);
      if (!filesExistPhonons(xInput)) {
        string message = "Creating " + xInput.getDirectory();
        pflow::logger(_AFLOW_FILE_NAME_, _APL_LRPC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
        createAflowInPhononsAIMS(_aflowFlags, _kbinFlags, _xFlags, _AflowIn, xInput, *p_FileMESSAGE);
        stagebreak = true;
      }
    }
    return stagebreak;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                             FORCE CONSTANTS                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  bool LinearResponsePC::calculateForceConstants() {
    // Check if supercell is already built
    if (!_supercell->isConstructed()) {
      string function = "apl::LinearResponsePC::calculateForceFields()";
      string message = "The supercell structure has not been initialized yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }

    // First pass - check if any of the calculations ran (gives no error message,
    // similar to DirectMethodPC).
    if (!outfileFoundAnywherePhonons(xInputs)) return false;
    if (!readForceConstantsFromVasprun(xInputs[0])) return false;
    if (_isPolarMaterial) {
      if (!calculateBornChargesDielectricTensor(xInputs[1])) return false;
    }
    return true;
  }

  //////////////////////////////////////////////////////////////////////////////
  //ME20200211
  bool LinearResponsePC::readForceConstantsFromVasprun(_xinput& xinp) {
    stringstream message;
    message << "Reading force constants from vasprun.xml";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_LRPC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);

    // Read vasprun.xml
    string filename = aurostd::CleanFileName(xinp.getDirectory() + "/vasprun.xml.static");
    if (!aurostd::EFileExist(filename)) {
      filename = aurostd::CleanFileName(xinp.getDirectory() + "/vasprun.xml");
      if (aurostd::EFileExist(filename)) {
        message << "Could not find vasprun.xml file for linear response calculations.";
        pflow::logger(_AFLOW_FILE_NAME_, _APL_LRPC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
        return false;
      }
    }
    vector<string> vlines;
    aurostd::efile2vectorstring(filename, vlines);
    uint nlines = vlines.size();

    // Read Hessian
    vector<vector<double> > hessian;
    vector<double> row;
    vector<string> line;
    uint iline = 0;
    for (iline = 0; iline < nlines; iline++) {
      if (aurostd::substring2bool(vlines[iline], "hessian")) {
        while ((++iline != nlines) && !aurostd::substring2bool(vlines[iline], "</varray>")) {
          aurostd::string2tokens(vlines[iline], line);
          row.clear();
          for (uint i = 1; i < (line.size() - 1); i++) {
            row.push_back(aurostd::string2utype<double>(line[i]));
          }
          hessian.push_back(row);
        }
        if (iline < nlines) break;
      }
    }

    string function = "apl::LinearResponsePC::readForceConstantsFromVasprun()";
    // Check that the file was read successfully.
    if (iline == nlines) {
      message << "Hessian tag not found or incomplete.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
    }
    uint natoms = _supercell->getSupercellStructure().atoms.size();
    uint nhessian = hessian.size();
    if (nhessian != 3 * natoms) {
      message << "Hessian matrix does not have the correct number of rows (has "
        << nhessian << ", should have " << (3 * natoms) << ")." << std::endl;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
    }
    uint i = 0;
    for (i = 0; i < nhessian; i++) {
      if (hessian[i].size() != 3 * natoms) break;
    }
    if (i != nhessian) {
      message << "Row " << (i + 1) << " of the Hessian matrix does not have the correct number of columns"
        << " (has " << hessian[i].size() << ", should have " << (3 * natoms) << ")." << std::endl;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
    }

    _forceConstantMatrices.clear();
    _forceConstantMatrices.resize(natoms, vector<xmatrix<double> >(natoms, xmatrix<double>(3, 3)));
    double mass = 0.0;
    // Convert Hessian matrix into force constants
    for (uint i = 0; i < natoms; i++) {
      for (uint j = 0; j < natoms; j++) {
        // Hessian matrix is normalized by masses, so multiply to get IFCs
        mass = std::sqrt(_supercell->getAtomMass(i) * _supercell->getAtomMass(j));
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            _forceConstantMatrices[i][j][k+1][l+1] = -mass * hessian[3 * i + k][3 * j + l];
          }
        }
      }
    }
    return true;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               FILE OUTPUT                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  void LinearResponsePC::saveState(const string& filename) {
    string function = "apl::LinearResponsePC::saveState()";
    string message = "";
    if (!_sc_set) {
      message = "Supercell pointer not set.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }
    if (xInputs.size() == 0) return;  // Nothing to write
    message = "Saving state of the force constant calculator into " + aurostd::CleanFileName(filename) + ".";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_LRPC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
    stringstream out;
    string tag = "[APL_FC_CALCULATOR]";
    out << AFLOWIN_SEPARATION_LINE << std::endl;
    out << tag << "ENGINE=LR" << std::endl;
    if (xInputs[0].AFLOW_MODE_VASP) out << tag << "AFLOW_MODE=VASP" << std::endl;
    else if (xInputs[0].AFLOW_MODE_AIMS) out << tag << "AFLOW_MODE=AIMS" << std::endl;
    out << AFLOWIN_SEPARATION_LINE << std::endl;
    out << tag << "SUPERCELL=" << _supercell->scell << std::endl;
    out << tag << "INPUT_STRUCTURE=START" << std::endl;
    out << _supercell->getInputStructure();  // No endl necessary
    out << tag << "INPUT_STRUCTURE=STOP" << std::endl;
    out << AFLOWIN_SEPARATION_LINE << std::endl;
    out << tag << "POLAR=" << _isPolarMaterial << std::endl;
    out << AFLOWIN_SEPARATION_LINE << std::endl;
    aurostd::stringstream2file(out, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Could not save state into file " + filename + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
  }

  // ME20200501
  // The state reader is designed to read the directory structure and supercell
  // structures from a prior run for post-processing. It cannot create APL
  // aflow.in files and should only be used to read forces for force constant
  // calculations. It is still in development and has only been tested with VASP.
  void LinearResponsePC::readFromStateFile(const string& filename) {
    string function = "apl::LinearResponsePC::readFromStateFile()";
    string message = "";
    if (!_sc_set) {
      message = "Supercell pointer not set.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }
    message = "Reading state of the phonon calculator from " + filename + ".";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_LRPC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
    if (!aurostd::EFileExist(filename)) {
      message = "Could not find file " + filename + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
    }

    // Find if the calculations are VASP or AIMS calculations
    vector<string> vlines, tokens;
    aurostd::efile2vectorstring(filename, vlines);
    uint nlines = vlines.size();
    uint iline = 0;
    _xinput xInput;
    while (++iline < nlines) {
      if (aurostd::substring2bool(vlines[iline], "AFLOW_MODE") ){
        aurostd::string2tokens(vlines[iline], tokens, "=");
        if (tokens.size() != 2) {
          message = "Tag for AFLOW_MODE is broken.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
        }
        if (tokens[1] == "VASP") {
          xInput.AFLOW_MODE_VASP = true;
        } else if (tokens[1] == "AIMS") {
          xInput.AFLOW_MODE_AIMS = true;
        } else {
          message = "Unknown AFLOW_MODE " + tokens[1] + ".";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
        }
      }
    }
    if (iline >= nlines) {
      message = "AFLOW_MODE tag is missing.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
    }

    // Defaults
    xInput.xvasp.AVASP_arun_mode = "APL";
    _isPolarMaterial = DEFAULT_APL_POLAR;

    // Set xInput for the linear response calculation
    xInputs.push_back(xInput);
    xInputs[0].setXStr(_supercell->getSupercellStructureLight());
    xInputs[0].xvasp.AVASP_arun_runname = "1_" + _AFLOW_APL_DFPT_RUNNAME_;

    // Read
    xInputs.clear();
    iline = 0;
    while (++iline < nlines) {
      if (aurostd::substring2bool(vlines[iline], "POLAR=")) {
        aurostd::string2tokens(vlines[iline], tokens, "=");
        if (tokens.size() != 2) {
          string message = "Tag for POLAR is broken.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
        }
        _isPolarMaterial = aurostd::string2utype<bool>(tokens[1]);
        if (_isPolarMaterial) {
          xInputs.push_back(xInput);
          xInputs[1].setXStr(_supercell->getInputStructureLight());
          xInputs[1].xvasp.AVASP_arun_runname = "2_" + _AFLOW_APL_BORN_EPSILON_RUNNAME_;
        }
      }
    }

    // Set directories
    string dir = "";
    for (uint i = 0; i < xInputs.size(); i++) {
      const _xvasp& xvasp = xInputs[i].xvasp;
      dir = _directory + "/ARUN." + xvasp.AVASP_arun_mode + "_" + xvasp.AVASP_arun_runname + "/";
      xInputs[i].setDirectory(aurostd::CleanFileName(dir));
    }
  }

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
