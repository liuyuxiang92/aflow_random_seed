#include "aflow_apl.h"

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////
  LinearResponsePC::LinearResponsePC(Supercell& sc,
      _xinput& xinput, _aflags& aflags, _kflags& kflags,
      _xflags& xflags, string& AflowIn, ofstream& mf)
    : ForceConstantCalculator(sc, xinput, aflags, kflags, xflags, AflowIn, mf) {
      free();
    }

  LinearResponsePC::LinearResponsePC(const LinearResponsePC& that)
    : ForceConstantCalculator(that._supercell, that._xInput, that._aflowFlags, that._kbinFlags,
      that._xFlags, that._AflowIn, *that.messageFile) {
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
    free();
  }

  void LinearResponsePC::clear(Supercell& sc, _xinput& xinput,
      _aflags& aflags, _kflags& kflags, _xflags& xflags, string& AflowIn, ofstream& mf) {
    LinearResponsePC that(sc, xinput, aflags, kflags, xflags, AflowIn, mf);
    copy(that);
  }

  void LinearResponsePC::copy(const LinearResponsePC& that) {
    _aflowFlags = that._aflowFlags;
    _AflowIn = that._AflowIn;
    _bornEffectiveChargeTensor = that._bornEffectiveChargeTensor;
    _dielectricTensor = that._dielectricTensor;
    messageFile = that.messageFile;
    _forceConstantMatrices = that._forceConstantMatrices;
    _isPolarMaterial = that._isPolarMaterial;
    _kbinFlags = that._kbinFlags;
    _supercell = that._supercell;
    _xFlags = that._xFlags;
    _xInput = that._xInput;
    xInputs = that.xInputs;
  }

  void LinearResponsePC::free() {
    xInputs.clear();
    _bornEffectiveChargeTensor.clear();
    _dielectricTensor.clear();
    _forceConstantMatrices.clear();
    _isPolarMaterial = false;
  }

  //////////////////////////////////////////////////////////////////////////////
  bool LinearResponsePC::runVASPCalculations(bool zerostate_chgcar) {
    if (zerostate_chgcar) {
      string message = "ZEROSTATE_CHGCAR not implemented for linear response calculations.";
      pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout, 'W');
    }
    bool stagebreak = false;

    // Call VASP to calculate forces by LR
    _xInput.xvasp.AVASP_arun_mode = "APL";
    xInputs.clear();
    xInputs.push_back(_xInput);
    stagebreak = runVASPCalculationsDFPT(xInputs[0]);
    if (_isPolarMaterial) {
      xInputs.push_back(_xInput);
      stagebreak = (runVASPCalculationsBE(xInputs[1], 2) || stagebreak);
    }
    return stagebreak;
  }

  //////////////////////////////////////////////////////////////////////////////
  // We will use VASP5.2+ to calculate Born effective charge tensors and
  // dielectric constant matrix in the primitive cell with very high precision
  // Both values are needed by non-analytical term of dynamic matrix for
  // correct TO-LO splitting of optical phonon branches of polar systems
  bool LinearResponsePC::runVASPCalculationsDFPT(_xinput& xInput) {
    bool stagebreak = false;

    xInput.setXStr(_supercell.getSupercellStructureLight());
    xInput.getXStr().title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xInput.getXStr().title);
    if(xInput.getXStr().title.empty()){xInput.getXStr().buildGenericTitle(true,false);}
    xInput.getXStr().title+=" linear response";

    // For VASP, use the standardized aflow.in creator
    if(xInput.AFLOW_MODE_VASP) {
      xInput.xvasp.AVASP_arun_runname = "1_" + _AFLOW_APL_DFPT_RUNNAME_;  // ME200213
      xInput.xvasp.aopts.flag("APL_FLAG::AVASP_BORN", false);
      xInput.xvasp.aopts.flag("APL_FLAG::AVASP_LR", true);
      _kbinFlags.KBIN_MPI_AUTOTUNE = false;

      // Set POSCAR to VASP5 format
      xInput.getXStr().is_vasp4_poscar_format = false;
      xInput.getXStr().is_vasp5_poscar_format = true;
      stagebreak = (createAflowInPhonons(_aflowFlags, _kbinFlags, _xFlags, xInput) || stagebreak);
    }
    // For AIMS, use the old method until we have AVASP_populateXAIMS
    if(xInput.AFLOW_MODE_AIMS) {
      string runname = _AFLOW_APL_DFPT_DIRECTORY_NAME_;  // ME200213
      xInput.setDirectory( _xInput.getDirectory() + "/" + runname );
      if (!filesExistPhonons(xInput)) {
        string message = "Creating " + xInput.getDirectory();
        pflow::logger(_AFLOW_FILE_NAME_, xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);
        createAflowInPhononsAIMS(_aflowFlags, _kbinFlags, _xFlags, _AflowIn, xInput, *messageFile);
        stagebreak = true;
      }
    }
    return stagebreak;
  }

  //////////////////////////////////////////////////////////////////////////////
  bool LinearResponsePC::calculateForceConstants() {
    // Check if supercell is already built
    if (!_supercell.isConstructed()) {
      // ME191031 - use xerror
      //throw APLRuntimeError("apl::LinearResponsePC::calculateForceFields(); The supercell structure has not been initialized yet.");
      string function = "apl::LinearResponsePC::calculateForceFields()";
      string message = "The supercell structure has not been initialized yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }

    // First pass - check if any of the calculations ran (gives no error message,
    // similar to DirectMethodPC).
    if (!outfileFoundAnywherePhonons(xInputs)) return false;
    if (!readForceConstantsFromVasprun(xInputs[0])) return false;
    if (_isPolarMaterial) {
      if (!calculateDielectricTensor(xInputs[1])) return false;
    }
    return true;
  }

  //////////////////////////////////////////////////////////////////////////////
  // ME200211
  bool LinearResponsePC::readForceConstantsFromVasprun(_xinput& xinp) {
    stringstream message;
    message << "Reading force constants from vasprun.xml";
    pflow::logger(_AFLOW_FILE_NAME_, xinp.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);
    string function = "apl::LinearResponsePC::readForceConstantsFromVasprun()";

    // Read vasprun.xml
    string filename = aurostd::CleanFileName(xinp.getDirectory() + "/vasprun.xml.static");
    if (!aurostd::EFileExist(filename)) {
      filename = aurostd::CleanFileName(xinp.getDirectory() + "/vasprun.xml");
      if (aurostd::EFileExist(filename)) {
        message << "Could not find vasprun.xml file for linear response calculations.";
        pflow::logger(_AFLOW_FILE_NAME_, xinp.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);
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

    // Check that the file was read successfully.
    if (iline == nlines) {
      message << "Hessian tag not found or incomplete.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
    }
    uint natoms = _supercell.getSupercellStructure().atoms.size();
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
        // Hessian matrix is normalized by masses, so multiply to get FCs
        mass = std::sqrt(_supercell.getAtomMass(i) * _supercell.getAtomMass(j));
        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            _forceConstantMatrices[i][j][k+1][l+1] = -mass * hessian[3 * i + k][3 * j + l];
          }
        }
      }
    }
    return true;
  }

  //////////////////////////////////////////////////////////////////////////////

  void LinearResponsePC::hibernate(const string& filename) {
    stringstream outfile;

    writeHibernateHeader(outfile);
    writeForceConstants(outfile);
    if (_isPolarMaterial) writePolar(outfile);
    outfile << "</apl>" << std::endl;

    aurostd::stringstream2file(outfile, filename); //ME181226
    if (!aurostd::FileExist(filename)) { //ME181226
      string function = "ForceConstantCalculator::hibernate()";
      string message = "Cannot open output file " + filename + "."; //ME181226
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
    }
  }


}  // namespace apl
