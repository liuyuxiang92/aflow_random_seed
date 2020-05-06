// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************

#include "aflow_apl.h"
#define _DEBUG_APL_HARM_IFCS_ false

using std::vector;
using std::string;

static const string _APL_FCCALC_MODULE_ = "APL";  // for the logger

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         CONSTRUCTORS/DESTRUCTORS                         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  ForceConstantCalculator::ForceConstantCalculator(ostream& oss) {
    free();
    xStream::initialize(oss);
    _directory = "./";
  }

  ForceConstantCalculator::ForceConstantCalculator(Supercell& sc, ofstream& mf, ostream& oss) {
    free();
    _supercell = &sc;
    _sc_set = true;
    xStream::initialize(mf, oss);
    _directory = "./";
  }

  ForceConstantCalculator::ForceConstantCalculator(const ForceConstantCalculator& that) {
    free();
    copy(that);
  }

  ForceConstantCalculator& ForceConstantCalculator::operator=(const ForceConstantCalculator& that) {
    if (this != &that) {
      free();
      copy(that);
    }
    return *this;
  }

  void ForceConstantCalculator::clear(Supercell& sc) {
    free();
    _supercell = &sc;
  }

  void ForceConstantCalculator::copy(const ForceConstantCalculator& that) {
    xStream::copy(that);
    _bornEffectiveChargeTensor = that._bornEffectiveChargeTensor;
    _dielectricTensor = that._dielectricTensor;
    _directory = that._directory;
    _forceConstantMatrices = that._forceConstantMatrices;
    _isPolarMaterial = that._isPolarMaterial;
    _sc_set = that._sc_set;
    _supercell = that._supercell;
    xInputs = that.xInputs;
  }

  void ForceConstantCalculator::free() {
    xInputs.clear();
    _bornEffectiveChargeTensor.clear();
    _dielectricTensor.clear();
    _directory = "";
    _forceConstantMatrices.clear();
    _isPolarMaterial = false;
    _sc_set = false;
    _supercell = NULL;
  }

}  // namespace apl

namespace apl {

  const vector<vector<xmatrix<double> > >& ForceConstantCalculator::getForceConstants() const {
    return _forceConstantMatrices;
  }

  const vector<xmatrix<double> >& ForceConstantCalculator::getBornEffectiveChargeTensor() const {
    return _bornEffectiveChargeTensor;
  }

  const xmatrix<double>& ForceConstantCalculator::getDielectricTensor() const {
    return _dielectricTensor;
  }

  bool ForceConstantCalculator::isPolarMaterial() const {
    return _isPolarMaterial;
  }

  const string& ForceConstantCalculator::getDirectory() const {
    return _directory;
  }

  void ForceConstantCalculator::setDirectory(const string& dir) {
    _directory = dir;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                             FORCE CONSTANTS                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // Runs the force constant calculator (main post-processing engine)
  bool ForceConstantCalculator::run() {
    string function = "apl::ForceConstantCalculator::run():";
    string message = "";
    if (!_sc_set) {
      message = "Supercell pointer not set.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }
    // Check if supercell is already built
    if (!_supercell->isConstructed()) {
      message = "The supercell structure has not been initialized yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }

    if (!calculateForceConstants()) return false;

    // ME20191219 - atomGoesTo and atomComesFrom can now use basis_atoms_map.
    // Calculating the full basis ahead of time is much faster than calculating all
    // symmetry operations on-the-fly.
    if (!_supercell->fullBasisCalculatedAGROUP()) _supercell->getFullBasisAGROUP();

    // Symmetrization of the force-constant matrices
    symmetrizeForceConstantMatrices();

    // Force the force-constant matrices to obey the sum-rule conditions
    correctSumRules();
    return true;
  }

  // Symmetrizes the force constant matrices using site point group symmetry
  void ForceConstantCalculator::symmetrizeForceConstantMatrices() {
    bool LDEBUG=(FALSE || _DEBUG_APL_HARM_IFCS_ || XHOST.DEBUG);
    string soliloquy="apl::ForceConstantCalculator::symmetrizeForceConstantMatrices()"; //CO20190218
    // Test of stupidity...
    if (!_supercell->getSupercellStructure().agroup_calculated) {
      string message = "The site groups have not been calculated yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, message, _RUNTIME_INIT_);
    }
    //CO - START
    if (_supercell->getEPS() == AUROSTD_NAN) {
      string message = "Need to define symmetry tolerance.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, message, _VALUE_ERROR_);
    }
    //CO - END

    string message = "Symmetrizing the force constant matrices.";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);

    vector<xmatrix<double> > row;
    for (int i = 0; i < _supercell->getNumberOfAtoms(); i++) {
      const vector<_sym_op>& agroup = _supercell->getAGROUP(i);  //CO //CO20190218
      if (agroup.size() == 0) {
        string message = "Site point group operations are missing.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, message, _RUNTIME_INIT_);
      }

      for (int j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        if(LDEBUG){ //CO20190218
          cerr << soliloquy << " compare original m=" << std::endl;
          cerr << _forceConstantMatrices[i][j] << std::endl;
        }
        xmatrix<double> m(3, 3); //CO20190218
        for (uint symOpID = 0; symOpID < agroup.size(); symOpID++) {
          const _sym_op& symOp = agroup[symOpID];

          try {
            //_AFLOW_APL_REGISTER_ int l = _supercell.atomComesFrom(symOp, j, i, FALSE);  //CO NEW //CO20190218
            // ME20191219 - atomGoesTo now uses basis_atoms_map; keep translation option in case
            // the basis has not been calculated for some reason
            _AFLOW_APL_REGISTER_ int l = _supercell->atomGoesTo(symOp, j, i, true); //JAHNATEK ORIGINAL //CO20190218
            m = m + (inverse(symOp.Uc) * _forceConstantMatrices[i][l] * symOp.Uc);  //JAHNATEK ORIGINAL //CO20190218
            //m = m + (symOp.Uc * _forceConstantMatrices[i][l] * inverse(symOp.Uc));  //CO NEW //CO20190218
            if(LDEBUG){ //CO20190218
              std::cerr << soliloquy << " atom[" << l << "].cpos=" << _supercell->getSupercellStructure().atoms[l].cpos << std::endl;
              std::cerr << soliloquy << " atom[" << j << "].cpos=" << _supercell->getSupercellStructure().atoms[j].cpos << std::endl;
              std::cerr << soliloquy << " agroup(" << l << " -> " << j << ")=" << std::endl;
              std::cerr << symOp.Uc << std::endl;
              std::cerr << soliloquy << " forceConstantMatrices[i=" << i << "][l=" << l << "]=" << std::endl;
              std::cerr << _forceConstantMatrices[i][l] << std::endl;
              std::cerr << soliloquy << " with new m=" << std::endl;
              std::cerr << (symOp.Uc * _forceConstantMatrices[i][l] * inverse(symOp.Uc)) << std::endl;
            }
          } catch (aurostd::xerror& e) {
            string message = "Mapping problem " + aurostd::utype2string<int>(j);
            throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, message, _RUNTIME_ERROR_);
          }
        }
        m = ( 1.0 / agroup.size() ) * m; //CO20190218
        row.push_back(m);
      }
      _forceConstantMatrices[i] = row;
      row.clear();
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME20200504 - this function needs to be rewritten to be more clear
  // Enfoces acoustic sum rules
  void ForceConstantCalculator::correctSumRules() {
    // ME20200504
    // sum appears to be the self-interaction term (diagonal terms) of the
    // force constant matrices. See http://cmt.dur.ac.uk/sjc/thesis_prt/node83.html
    xmatrix<double> sum(3, 3);//, sum2(3, 3); OBSOLETE ME20200504 - not used

    for (int i = 0; i < _supercell->getNumberOfAtoms(); i++) {
      // ME20200504 - sums are not used or cleared before they are used
      //[OBSOLETE] // Get SUMs
      //[OBSOLETE] for (int j = 0; j < _supercell->getNumberOfAtoms(); j++) {
      //[OBSOLETE]   if (i != j) {
      //[OBSOLETE]     sum = sum + _forceConstantMatrices[i][j];
      //[OBSOLETE]     //sum2 = sum2 + trasp(_forceConstantMatrices[j][i]);
      //[OBSOLETE]   }
      //[OBSOLETE] }

      // Correct SUM2
      // ME20200504 - This appears to enforce the invariance of the force constants
      // upon permutations
      for (int j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        if (i == j) continue;
        _forceConstantMatrices[i][j] = 0.5 * (_forceConstantMatrices[i][j] + trasp(_forceConstantMatrices[j][i]));
        _forceConstantMatrices[j][i] = trasp(_forceConstantMatrices[i][j]);
      }

      // Get SUMs again
      sum.clear();
      //sum2.clear(); OBSOLETE ME20200504 - not used
      for (int j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        if (i != j) {
          sum = sum + _forceConstantMatrices[i][j];
          //sum2 = sum2 + trasp(_forceConstantMatrices[j][i]);  // OBSOLETE ME20200504 - not used
        }
      }

      // Correct SUM1 to satisfied
      // ME20200504 - Self-interaction term
      _forceConstantMatrices[i][i] = -sum;
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                          BORN/DIELECTRIC TENSOR                          //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {
 
  // Sets up the calculations that determine the Born effective charges and
  // the dielectric tensor
  bool ForceConstantCalculator::runVASPCalculationsBE(_xinput& xInput, _aflags& _aflowFlags,
      _kflags& _kbinFlags, _xflags& _xFlags, string& _AflowIn, uint ncalcs) {
    bool stagebreak = false;

    xInput.setXStr(_supercell->getInputStructure());
    xInput.getXStr().title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xInput.getXStr().title);
    if(xInput.getXStr().title.empty()){xInput.getXStr().buildGenericTitle(true,false);}
    xInput.getXStr().title+=" Born effective charges/dielectric tensor";

    // For VASP, use the standardized aflow.in creator
    if(xInput.AFLOW_MODE_VASP) {
      xInput.xvasp.AVASP_arun_runname = aurostd::utype2string<uint>(ncalcs) + "_" + _AFLOW_APL_BORN_EPSILON_RUNNAME_;
      xInput.xvasp.aopts.flag("APL_FLAG::AVASP_LR", false);
      xInput.xvasp.aopts.flag("APL_FLAG::AVASP_BORN", true);
      // Set POSCAR to VASP5 format
      xInput.getXStr().is_vasp4_poscar_format = false;
      xInput.getXStr().is_vasp5_poscar_format = true;
      stagebreak = (createAflowInPhonons(_aflowFlags, _kbinFlags, _xFlags, xInput) || stagebreak);
    } else if (xInput.AFLOW_MODE_AIMS) {
      string runname = _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_;
      xInput.setDirectory(_directory + "/" + runname );
      if (!filesExistPhonons(xInput)) {
        string message = "Creating " + xInput.getDirectory();
        pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
        createAflowInPhononsAIMS(_aflowFlags, _kbinFlags, _xFlags, _AflowIn, xInput, *p_FileMESSAGE);
        stagebreak = true;
      }
    }
    return stagebreak;
  }

  //////////////////////////////////////////////////////////////////////////////
  
  // Calculates the dielectric tensor and Born effective charges
  bool ForceConstantCalculator::calculateBornChargesDielectricTensor(const _xinput& xinpBE) {
    stringstream message;
    // Parse effective charges from OUTCAR
    if (xinpBE.AFLOW_MODE_VASP) {
      string directory = xinpBE.xvasp.Directory;
      string infilename = directory + "/OUTCAR.static";

      if (!aurostd::EFileExist(infilename, infilename)) {
        infilename = directory + string("/OUTCAR");
        if (!aurostd::EFileExist(infilename, infilename)) {
          message << "The OUTCAR file in " << directory << " is missing.";
          pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
          return false;
        }
      }
    } else if (xinpBE.AFLOW_MODE_AIMS) {
      string function = "ForceConstantCalculator::calculateBornChargesDielectricTensor()";
      message << "This functionality has not been implemented yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    } else {
      return false;
    }

    if (xinpBE.AFLOW_MODE_VASP) readBornEffectiveChargesFromOUTCAR(xinpBE);
    else if (xinpBE.AFLOW_MODE_AIMS) readBornEffectiveChargesFromAIMSOUT();

    // Enforce ASR (Acoustic sum rules)
    symmetrizeBornEffectiveChargeTensors();

    // Parse epsilon from OUTCAR
    if(xinpBE.AFLOW_MODE_VASP) readDielectricTensorFromOUTCAR(xinpBE);
    if(xinpBE.AFLOW_MODE_AIMS) readDielectricTensorFromAIMSOUT();

    message << "Dielectric tensor: ";
    for (int a = 1; a <= 3; a++)
      for (int b = 1; b <= 3; b++)
        message << std::fixed << std::setw(5) << std::setprecision(3) << _dielectricTensor(a, b) << " ";

    pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
    return true;
  }

  //////////////////////////////////////////////////////////////////////////////
  void ForceConstantCalculator::readBornEffectiveChargesFromAIMSOUT(void) {
    string function = "ForceConstantCalculator::readBornEffectiveChargesFromAIMSOUT()";
    string message = "This functionality has not been implemented yet.";
    throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
  }

  void ForceConstantCalculator::readBornEffectiveChargesFromOUTCAR(const _xinput& xinp) {  // ME20190113
    string directory = xinp.xvasp.Directory;  // ME20190113
    string function = "apl::ForceConstantCalculator::readBornEffectiveChargesFromOUTCAR():";

    //CO - START
    string infilename = directory + string("/OUTCAR.static");

    if (!aurostd::EFileExist(infilename, infilename)) {
      // We already know that one of the files exists, so
      // if OUTCAR.static was not found, it must be OUTCAR
      infilename = directory + string("/OUTCAR");
    }

    // Open our file
    //CO - START
    vector<string> vlines;
    aurostd::efile2vectorstring(infilename, vlines);
    if (!vlines.size()) {
      //CO - END
      string message = "Cannot open input file OUTCAR.static.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }

    string line = "";
    uint line_count = 0;  //CO
    string KEY = ""; //ME20181226
    if (DEFAULT_APL_USE_LEPSILON) { //ME20181226
      KEY = string("BORN EFFECTIVE CHARGES (in e, cummulative output)");//ME20181226
    } else { //ME20181226
      KEY = string("BORN EFFECTIVE CHARGES (including local field effects)"); //ME20181226
    } //ME20181226

    while (true) {
      // Get line
      //CO - START
      if (line_count == vlines.size()) {
        string message = "No information on Born effective charges in OUTCAR file.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
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
    line = vlines[line_count++]; // Skip line "----------------...."
    for (uint i = 0; i < _supercell->getInputStructure().atoms.size(); i++) {
      // Get atom ID but not use it...
      line = vlines[line_count++];

      // Get its charge tensor
      for (int j = 1; j <= 3; j++) {
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
    //CO - END
  }

  //////////////////////////////////////////////////////////////////////////////
  void ForceConstantCalculator::symmetrizeBornEffectiveChargeTensors(void) {
    //CO - START
    // Test of stupidity...
    stringstream message;
    if (_supercell->getEPS() == AUROSTD_NAN) {
      string function = "apl::ForceConstantCalculator::symmetrizeEffectiveChargeTensors()";
      message << "Symmetry tolerance not defined.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ERROR_);
    }
    //CO - END
    // Show charges
    message << "Input born effective charge tensors (for primitive cell):";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
    for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++) {
      int id = i;
      message << "Atom [" << aurostd::PaddedNumString(id, 3) << "] ("
        << std::setw(2) << _supercell->getInputStructure().atoms[id].cleanname
        << ") Born effective charge = ";
      pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
      for (int a = 1; a <= 3; a++)
        for (int b = 1; b <= 3; b++)
          message << std::fixed << std::setw(5) << std::setprecision(3) << _bornEffectiveChargeTensor[i](a, b) << " ";
      pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
    }

    // Step1
    for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
      int basedUniqueAtomID = _supercell->getUniqueAtomID(i);

      xmatrix<double> sum(3, 3);
      for (int j = 0; j < _supercell->getNumberOfEquivalentAtomsOfType(i); j++) { //CO20190218
        try {  //CO
          const _sym_op& symOp = _supercell->getSymOpWhichMatchAtoms(_supercell->getUniqueAtomID(i, j), basedUniqueAtomID, _FGROUP_);
          sum += inverse(symOp.Uc) * _bornEffectiveChargeTensor[_supercell->sc2pcMap(_supercell->getUniqueAtomID(i, j))] * symOp.Uc;
        }
        //CO - START
        catch (aurostd::xerror& e) {
          string function = "apl::ForceConstantCalculator::symmetrizeBornEffectiveChargeTensors()";
          stringstream message;
          message << "Mapping problem " << _supercell->getUniqueAtomID(i, j) << " <-> " << basedUniqueAtomID << "?";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
        }
        //CO - END
      }

      sum = (1.0 / _supercell->getNumberOfEquivalentAtomsOfType(i)) * sum; //CO20190218

      for (int j = 0; j < _supercell->getNumberOfEquivalentAtomsOfType(i); j++) //CO20190218
        _bornEffectiveChargeTensor[_supercell->sc2pcMap(_supercell->getUniqueAtomID(i, j))] = sum;
    }

    // Step2
    vector<xmatrix<double> > newbe = _bornEffectiveChargeTensor;
    const vector<vector<_sym_op> >& agroup = _supercell->getAGROUP();  //CO
    for (int i = 0; i < _supercell->getNumberOfAtoms(); i++) {
      // Translate the center to this atom
      _supercell->center(i);

      xmatrix<double> sum(3, 3);
      for (uint symOpID = 0; symOpID < agroup[i].size(); symOpID++) {
        const _sym_op& symOp = agroup[i][symOpID];
        sum = sum + (inverse(symOp.Uc) * _bornEffectiveChargeTensor[_supercell->sc2pcMap(i)] * symOp.Uc);
      }
      newbe[_supercell->sc2pcMap(i)] = (1.0 / agroup[i].size()) * sum;
    }
    // Translate the center back
    _supercell->center_original();  //CO

    _bornEffectiveChargeTensor.clear();
    _bornEffectiveChargeTensor = newbe;
    newbe.clear();

    // Step 3
    message << "Forcing the acoustic sum rule (ASR). Resulting born effective charges (for the supercell):";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);

    xmatrix<double> sum(3, 3);
    for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++)
      sum += _bornEffectiveChargeTensor[i];
    sum = (1.0 / _bornEffectiveChargeTensor.size()) * sum;
    for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++)
      _bornEffectiveChargeTensor[i] -= sum;

    // Make list only for unique atoms
    for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++)
      newbe.push_back(_bornEffectiveChargeTensor[_supercell->sc2pcMap(_supercell->getUniqueAtomID(i))]);
    _bornEffectiveChargeTensor.clear();
    _bornEffectiveChargeTensor = newbe;
    newbe.clear();

    // Show charges
    for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
      int id = _supercell->getUniqueAtomID(i);
      message << "Atom [" << aurostd::PaddedNumString(id, 3) << "] ("
        << std::setw(2) << _supercell->getSupercellStructure().atoms[id].cleanname
        << ") Born effective charge = ";
      for (int a = 1; a <= 3; a++)
        for (int b = 1; b <= 3; b++)
          message << std::fixed << std::setw(5) << std::setprecision(3) << _bornEffectiveChargeTensor[i](a, b) << " ";
      pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  void ForceConstantCalculator::readDielectricTensorFromAIMSOUT(void) {
    string function = "apl::ForceConstantCalculator::readDielectricTensorFromAIMSOUT()";
    string message = "This functionality has not been implemented yet.";
    throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
  }

  //////////////////////////////////////////////////////////////////////////////
  void ForceConstantCalculator::readDielectricTensorFromOUTCAR(const _xinput& xinp) {  // ME20190113
    string directory = xinp.xvasp.Directory;  // ME20190113

    //CO - START
    string infilename = directory + string("/OUTCAR.static");
    if (!aurostd::EFileExist(infilename, infilename)) {
      // We already checked outside if one of the files exists, so if
      // it is not OUTCAR.static, it must be OUTCAR
      infilename = directory + string("/OUTCAR");
    }

    // Open our file
    vector<string> vlines;
    uint line_count = 0;
    string line;
    aurostd::efile2vectorstring(infilename, vlines);
    if (!vlines.size()) {
      string function = "apl::ForceConstantCalculator::readDielectricTensorFromOUTCAR():";
      string message = "Cannot open input file OUTCAR.";
      aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
    //CO - END

    // Find
    string KEY = ""; //ME20181226
    if (DEFAULT_APL_USE_LEPSILON) { //ME20181226
      KEY = string("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)"); //ME20181226
    } else { //ME20181226
      KEY = string("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects)"); //ME20181226
    }//ME20181226

    while (true) {
      // Get line
      //CO - START
      if (line_count == vlines.size()) {
      string function = "apl::ForceConstantCalculator::readDielectricTensorFromOUTCAR():";
        string message = "No information on dielectric tensor in OUTCAR.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
      }
      line = vlines[line_count++];
      //CO - END

      // Check for our key line
      if (line.size() < KEY.size()) continue;
      if (line.find(KEY) != string::npos) break;
    }

    // Read in all ...
    //CO - START
    line = vlines[line_count++]; // Skip line "----------------...."

    // Get it
    vector<string> tokens;
    for (int j = 1; j <= 3; j++) {
      line = vlines[line_count++];
      tokenize(line, tokens, string(" "));
      _dielectricTensor(j, 1) = aurostd::string2utype<double>(tokens.at(0));
      _dielectricTensor(j, 2) = aurostd::string2utype<double>(tokens.at(1));
      _dielectricTensor(j, 3) = aurostd::string2utype<double>(tokens.at(2));
      tokens.clear();
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               FILE OUTPUT                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // Writes the results into xml files
  void ForceConstantCalculator::hibernate() {
    string function = "ForceConstantCalculator::hibernate():";
    string message = "";
    if (!_sc_set) {
      message = "Supercell pointer not set.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }
    string base = _directory + "/" + DEFAULT_APL_FILE_PREFIX;
    string filename = aurostd::CleanFileName(base + DEFAULT_APL_HARMIFC_FILE);
    message = "Writing harmonic IFCs into file " + filename + ".";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
    writeHarmonicIFCs(filename);
    if (_isPolarMaterial) {
      filename = aurostd::CleanFileName(base + DEFAULT_APL_POLAR_FILE);
      message = "Writing harmonic IFCs into file " + filename + ".";
      writeBornChargesDielectricTensor(filename);
      pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, _directory, *p_FileMESSAGE, *p_oss);
    }
  }

  // Writes the harmonic force constants into an xml file
  void ForceConstantCalculator::writeHarmonicIFCs(const string& filename) {
    stringstream outfile;
    string tab = " ";

    // Header
    outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
    outfile << "<apl>" << std::endl;
    outfile << tab << "<generator>" << std::endl;
    outfile << tab << tab << "<i name=\"aflow_version\" type=\"string\">" << AFLOW_VERSION << "</i>" << std::endl;
    string time = aflow_get_time_string();
    if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
    outfile << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
    // OBSOLETE ME20200427 - we do not compare checksums anymore
    //outfile << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_ << "\" type=\"" << APL_CHECKSUM_ALGO << "\">"
    //  << std::hex << aurostd::getFileCheckSum(_directory + "/" + _AFLOWIN_ + "", APL_CHECKSUM_ALGO) << "</i>" << std::endl;  // ME20190219
    //outfile.unsetf(std::ios::hex); //CO20190116 - undo hex immediately
    outfile << tab << "</generator>" << std::endl;

    // Force constants
    outfile << tab << "<fcms units=\"eV/Angstrom^2\" cs=\"cartesian\" rows=\""
      << _forceConstantMatrices.size() << "\" cols=\""
      << _forceConstantMatrices[0].size() << "\">" << std::endl;

    outfile << tab << tab << "<varray>" << std::endl;
    for (uint i = 0; i < _forceConstantMatrices.size(); i++) {
      outfile << tab << tab << tab << "<varray row=\"" << i << "\">" << std::endl;
      for (uint j = 0; j < _forceConstantMatrices[i].size(); j++) {
        outfile << tab << tab << tab << tab << "<matrix row=\"" << i
          << "\" col=\"" << j << "\">" << std::endl;
        for (int k = 1; k <= 3; k++) {
          outfile << tab << tab << tab << tab << tab << "<v>";
          for (int l = 1; l <= 3; l++) {
            outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
            outfile << setprecision(15);
            outfile << setw(24) << std::scientific << _forceConstantMatrices[i][j](k, l) << " ";
          }
          outfile << "</v>" << std::endl;
        }
        outfile << tab << tab << tab << tab << "</matrix>" << std::endl;
      }
      outfile << tab << tab << tab << "</varray>" << std::endl;
    }
    outfile << tab << tab << "</varray>" << std::endl;
    outfile << tab << "</fcms>" << std::endl;
    outfile << "</apl>" << std::endl;

    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)) {
      string function = "apl::ForceConstantCalculator::writeHarmonicIFCs():";
      string message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
    }
  }

  // Writes the Born effective charges and the dielectric tensor into an xml file
  void ForceConstantCalculator::writeBornChargesDielectricTensor(const string& filename) {
    stringstream outfile;
    string tab = " ";

    // Header
    outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
    outfile << "<apl>" << std::endl;
    outfile << tab << "<generator>" << std::endl;
    string time = aflow_get_time_string();
    if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
    outfile << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
    // OBSOLETE ME20200428 - Checksums are not used anymore
    //outfile << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_ << "\" type=\"" << APL_CHECKSUM_ALGO << "\">"
    //  << std::hex << aurostd::getFileCheckSum(_directory + "/" + _AFLOWIN_ + "", APL_CHECKSUM_ALGO) << "</i>" << std::endl;  // ME20190219
    //outfile.unsetf(std::ios::hex); //CO20190116 - undo hex immediately
    outfile << tab << "</generator>" << std::endl;

    // Born effective charge tensors
    outfile << tab << "<born units=\"a.u.\" cs=\"cartesian\">" << std::endl;
    outfile << tab << tab << "<varray>" << std::endl;
    for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++) {
      int id = _supercell->getUniqueAtomID(i);
      outfile << tab << tab << tab << "<matrix type=\"" << _supercell->getSupercellStructure().atoms[id].cleanname << "\">" << std::endl;
      for (int k = 1; k <= 3; k++) {
        outfile << tab << tab << tab << tab << "<v>";
        for (int l = 1; l <= 3; l++) {
          outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
          outfile << setprecision(8);
          // ME20181030 - fixed prevents hexadecimal output
          outfile << setw(15) << std::fixed << _bornEffectiveChargeTensor[i](k, l) << " ";
        }
        outfile << "</v>" << std::endl;
      }
      outfile << tab << tab << tab << "</matrix>" << std::endl;
    }
    outfile << tab << tab << "</varray>" << std::endl;
    outfile << tab << "</born>" << std::endl;

    // Dielectric tensor
    outfile << tab << "<epsilon units=\"a.u.\" cs=\"cartesian\">" << std::endl;
    outfile << tab << tab << "<matrix>" << std::endl;
    for (int k = 1; k <= 3; k++) {
      outfile << tab << tab << tab << "<v>";
      for (int l = 1; l <= 3; l++) {
        outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        outfile << setprecision(8);
        // ME20181030 - fixed prevents hexadecimal output
        outfile << setw(15) << std::fixed << _dielectricTensor(k, l) << " ";
      }
      outfile << "</v>" << std::endl;
    }
    outfile << tab << tab << "</matrix>" << std::endl;
    outfile << tab << "</epsilon>" << std::endl;
    outfile << "</apl>" << std::endl;

    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)) {
      string function = "ForceConstantCalculator::writeHarmonicIFCs()";
      string message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
    }
  }

  // ///////////////////////////////////////////////////////////////////////////
  // OBSOLETE ME20200504  - not used
  //[OBSOLETE] void ForceConstantCalculator::printForceConstantMatrices(ostream& os) {
  //[OBSOLETE]   int units = 1;
  //[OBSOLETE]   double conversionFactor = 1.0;

  //[OBSOLETE]   switch (units) {
  //[OBSOLETE]     case (1):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in eV/A^2:" << std::endl;
  //[OBSOLETE]       conversionFactor = 1.0;
  //[OBSOLETE]       break;
  //[OBSOLETE]     case (2):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in 10 Dyn/cm:" << std::endl;
  //[OBSOLETE]       conversionFactor = 1602.17733;
  //[OBSOLETE]       break;
  //[OBSOLETE]     case (3):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in N/m:" << std::endl;
  //[OBSOLETE]       conversionFactor = 16.0217733;
  //[OBSOLETE]       break;
  //[OBSOLETE]   }
  //[OBSOLETE]   os << std::endl;

  //[OBSOLETE]   for (int i = 0; i < _supercell->getNumberOfAtoms(); i++) {
  //[OBSOLETE]     for (int k = 0; k < _supercell->getNumberOfAtoms(); k++) {
  //[OBSOLETE]       os << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  //[OBSOLETE]       os << setprecision(4);
  //[OBSOLETE]       os << "- MATRIX: " << i + 1 << "/" << k + 1 << " " << k + 1 << "/" << i + 1 << std::endl;
  //[OBSOLETE]       for (int m = 1; m <= 3; m++) {
  //[OBSOLETE]         for (int n = 1; n <= 3; n++)
  //[OBSOLETE]           os << setw(10) << (conversionFactor * _forceConstantMatrices[k][i](m, n)) << " ";
  //[OBSOLETE]         os << " ";
  //[OBSOLETE]         for (int n = 1; n <= 3; n++)
  //[OBSOLETE]           os << setw(10) << (conversionFactor * _forceConstantMatrices[i][k](n, m)) << " ";
  //[OBSOLETE]         os << std::endl;
  //[OBSOLETE]       }
  //[OBSOLETE]       os << std::endl;
  //[OBSOLETE]     }
  //[OBSOLETE]   }
  //[OBSOLETE] }

  //[OBSOLETE] // ///////////////////////////////////////////////////////////////////////////

  //[OBSOLETE] void ForceConstantCalculator::printFCShellInfo(ostream& os) {
  //[OBSOLETE]   int units = 4;
  //[OBSOLETE]   double conversionFactor = 1.0;

  //[OBSOLETE]   switch (units) {
  //[OBSOLETE]     case (1):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in eV/A^2:" << std::endl;
  //[OBSOLETE]       conversionFactor = 1.0;
  //[OBSOLETE]       break;
  //[OBSOLETE]     case (2):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in 10 Dyn/cm:" << std::endl;
  //[OBSOLETE]       conversionFactor = 1602.17733;
  //[OBSOLETE]       break;
  //[OBSOLETE]     case (3):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in N/m:" << std::endl;
  //[OBSOLETE]       conversionFactor = 16.0217733;
  //[OBSOLETE]       break;
  //[OBSOLETE]     case (4):
  //[OBSOLETE]       os << "FORCE CONSTANT MATRICES in 10^3 Dyn/cm:" << std::endl;
  //[OBSOLETE]       conversionFactor = 16.0217733;
  //[OBSOLETE]       break;
  //[OBSOLETE]   }
  //[OBSOLETE]   os << std::endl;

  //[OBSOLETE]   int maxshell = _supercell->getMaxShellID();
  //[OBSOLETE]   if (maxshell == -1) maxshell = 25;
  //[OBSOLETE]   std::vector<ShellHandle> sh;
  //[OBSOLETE]   for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
  //[OBSOLETE]     ShellHandle s;
  //[OBSOLETE]     sh.push_back(s);
  //[OBSOLETE]     sh.back().init(_supercell->getInputStructure(),
  //[OBSOLETE]         _supercell->getInputStructure().iatoms[i][0],
  //[OBSOLETE]         maxshell);
  //[OBSOLETE]     sh[i].splitBySymmetry();
  //[OBSOLETE]     sh[i].mapStructure(_supercell->getSupercellStructure(), _supercell->getUniqueAtomID(i));
  //[OBSOLETE]   }

  //[OBSOLETE]   for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
  //[OBSOLETE]     sh[i].printReport(cout);
  //[OBSOLETE]     for (int ishell = 0; ishell <= sh[i].getLastOccupiedShell(); ishell++) {
  //[OBSOLETE]       for (int isubshell = 0; isubshell < sh[i].getNumberOfSubshells(ishell); isubshell++) {
  //[OBSOLETE]         const deque<_atom>& atomsAtSameShell = sh[i].getAtomsAtSameShell(ishell, isubshell);
  //[OBSOLETE]         cout << "SHELL " << ishell << " " << isubshell << std::endl;

  //[OBSOLETE]         for (uint ai = 0; ai < atomsAtSameShell.size(); ai++) {
  //[OBSOLETE]           int nb = atomsAtSameShell[ai].number;
  //[OBSOLETE]           cout << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
  //[OBSOLETE]           cout << setprecision(4);
  //[OBSOLETE]           cout << "- MATRIX: " << i << "/" << nb << " " << nb << "/" << i << std::endl;
  //[OBSOLETE]           for (int m = 1; m <= 3; m++) {
  //[OBSOLETE]             for (int n = 1; n <= 3; n++)
  //[OBSOLETE]               cout << setw(10) << (conversionFactor * _forceConstantMatrices[nb][i](m, n)) << " ";
  //[OBSOLETE]             cout << " ";
  //[OBSOLETE]             for (int n = 1; n <= 3; n++)
  //[OBSOLETE]               cout << setw(10) << (conversionFactor * _forceConstantMatrices[i][nb](n, m)) << " ";
  //[OBSOLETE]             cout << std::endl;
  //[OBSOLETE]           }
  //[OBSOLETE]           cout << std::endl;
  //[OBSOLETE]         }
  //[OBSOLETE]       }
  //[OBSOLETE]     }
  //[OBSOLETE]   }

  //[OBSOLETE]   // Clear
  //[OBSOLETE]   for (uint i = 0; i < sh.size(); i++)
  //[OBSOLETE]     sh[i].clear();
  //[OBSOLETE]   sh.clear();
  //[OBSOLETE] }

}  // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
