//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
//****************************************************************************

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

  ForceConstantCalculator::ForceConstantCalculator() {
    free();
  }

  ForceConstantCalculator::ForceConstantCalculator(Supercell& sc, _xinput& xinput,
      _aflags& aflags, _kflags& kflags, _xflags& xflags, string& AflowIn, ofstream& mf, ostream& os) {
    free();
    _supercell = &sc;
    _xInput = &xinput;
    _aflowFlags = &aflags;
    _kbinFlags = &kflags;
    _xFlags = &xflags;
    _AflowIn = &AflowIn;
    messageFile = &mf;
    oss = &os;
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

  void ForceConstantCalculator::clear(Supercell& sc, _xinput& xinput,
      _aflags& aflags, _kflags& kflags, _xflags& xflags, string& AflowIn, ofstream& mf, ostream& os) {
    free();
    _supercell = &sc;
    _xInput = &xinput;
    _aflowFlags =  &aflags;
    _kbinFlags = &kflags;
    _xFlags = &xflags;
    _AflowIn = &AflowIn;
    messageFile = &mf;
    oss = &os;
  }

  void ForceConstantCalculator::copy(const ForceConstantCalculator& that) {
    _aflowFlags = that._aflowFlags;
    _AflowIn = that._AflowIn;
    _bornEffectiveChargeTensor = that._bornEffectiveChargeTensor;
    _dielectricTensor = that._dielectricTensor;
    messageFile = that.messageFile;
    oss = that.oss;
    _forceConstantMatrices = that._forceConstantMatrices;
    _isPolarMaterial = that._isPolarMaterial;
    _kbinFlags = that._kbinFlags;
    _supercell = that._supercell;
    _xFlags = that._xFlags;
    _xInput = that._xInput;
    xInputs = that.xInputs;
  }

  void ForceConstantCalculator::free() {
    xInputs.clear();
    _bornEffectiveChargeTensor.clear();
    _dielectricTensor.clear();
    _forceConstantMatrices.clear();
    _isPolarMaterial = false;
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

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                             FORCE CONSTANTS                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  bool ForceConstantCalculator::run() {
    // Check if supercell is already built
    if (!_supercell->isConstructed()) {
      string function = "apl::ForceConstantCalculator::run()";
      string message = "The supercell structure has not been initialized yet.";
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

  void ForceConstantCalculator::symmetrizeForceConstantMatrices() {
    bool LDEBUG=(FALSE || _DEBUG_APL_HARM_IFCS_ || XHOST.DEBUG);
    string soliloquy="apl::ForceConstantCalculator::symmetrizeForceConstantMatrices()"; //CO20190218
    // Test of stupidity...
    if (!_supercell->getSupercellStructure().agroup_calculated) {
      string function = "apl::ForceConstantCalculator::symmetrizeForceConstantMatrices()";
      string message = "The site groups have not been calculated yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }
    //CO - START
    if (_supercell->getEPS() == AUROSTD_NAN) {
      string function = "apl::ForceConstantCalculator::symmetrizeForceConstantMatrices()";
      string message = "Need to define symmetry tolerance.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ERROR_);
    }
    //CO - END

    string message = "Symmetrizing the force constant matrices.";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);

    vector<xmatrix<double> > row;
    for (int i = 0; i < _supercell->getNumberOfAtoms(); i++) {
      const vector<_sym_op>& agroup = _supercell->getAGROUP(i);  //CO //CO20190218
      if (agroup.size() == 0) {
        string function = "apl::ForceConstantCalculator::symmetrizeForceConstantMatrices()";
        string message = "Site point group operations are missing.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
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
            string function = "apl::ForceConstantCalculator::symmetrizeForceConstantMatrices()";
            string message = "Mapping problem " + aurostd::utype2string<int>(j);
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
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

  void ForceConstantCalculator::correctSumRules() {
    xmatrix<double> sum(3, 3), sum2(3, 3);

    for (int i = 0; i < _supercell->getNumberOfAtoms(); i++) {
      // Get SUMs
      for (int j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        if (i != j) {
          sum = sum + _forceConstantMatrices[i][j];
          sum2 = sum2 + trasp(_forceConstantMatrices[j][i]);
        }
      }

      // Correct SUM2
      for (int j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        if (i == j) continue;
        _forceConstantMatrices[i][j] = 0.5 * (_forceConstantMatrices[i][j] + trasp(_forceConstantMatrices[j][i]));
        _forceConstantMatrices[j][i] = trasp(_forceConstantMatrices[i][j]);
      }

      // Get SUMs again
      sum.clear();
      sum2.clear();
      for (int j = 0; j < _supercell->getNumberOfAtoms(); j++) {
        if (i != j) {
          sum = sum + _forceConstantMatrices[i][j];
          sum2 = sum2 + trasp(_forceConstantMatrices[j][i]);
        }
      }

      // Correct SUM1 to satisfied
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
 
  bool ForceConstantCalculator::runVASPCalculationsBE(_xinput& xInput, uint ncalcs) {
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
      // Switch off autotune
      _kbinFlags->KBIN_MPI_AUTOTUNE = false;
      // Set POSCAR to VASP5 format
      xInput.getXStr().is_vasp4_poscar_format = false;
      xInput.getXStr().is_vasp5_poscar_format = true;
      stagebreak = (createAflowInPhonons(*_aflowFlags, *_kbinFlags, *_xFlags, xInput) || stagebreak);
    } else if (xInput.AFLOW_MODE_AIMS) {
      string runname = _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_;
      xInput.setDirectory( _xInput->getDirectory() + "/" + runname );
      if (!filesExistPhonons(xInput)) {
        string message = "Creating " + xInput.getDirectory();
        pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);
        createAflowInPhononsAIMS(*_aflowFlags, *_kbinFlags, *_xFlags, *_AflowIn, xInput, *messageFile);
        stagebreak = true;
      }
    }
    return stagebreak;
  }

  //////////////////////////////////////////////////////////////////////////////
  
  bool ForceConstantCalculator::calculateDielectricTensor(const _xinput& xinpBE) {
    stringstream message;
    // Parse effective charges from OUTCAR
    if (_kbinFlags->AFLOW_MODE_VASP) {
      string directory = xinpBE.xvasp.Directory;
      string infilename = directory + "/OUTCAR.static";

      if (!aurostd::EFileExist(infilename, infilename)) {
        infilename = directory + string("/OUTCAR");
        if (!aurostd::EFileExist(infilename, infilename)) {
          message << "The OUTCAR file in " << directory << " is missing.";
          pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);
          return false;
        }
      }
    } else if (_kbinFlags->AFLOW_MODE_AIMS) {
      string function = "ForceConstantCalculator::readBornEffectiveChargesFromAIMSOUT()";
      message << "This functionality has not been implemented yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    } else {
      return false;
    }

    if(_kbinFlags->AFLOW_MODE_VASP) readBornEffectiveChargesFromOUTCAR(xinpBE);
    else if(_kbinFlags->AFLOW_MODE_AIMS) readBornEffectiveChargesFromAIMSOUT();

    // Enforce ASR (Acoustic sum rules)
    symmetrizeBornEffectiveChargeTensors();

    // Parse epsilon from OUTCAR
    if(_kbinFlags->AFLOW_MODE_VASP) readDielectricTensorFromOUTCAR(xinpBE);
    if(_kbinFlags->AFLOW_MODE_AIMS) readDielectricTensorFromAIMSOUT();

    message << "Dielectric tensor: ";
    for (int a = 1; a <= 3; a++)
      for (int b = 1; b <= 3; b++)
        message << std::fixed << std::setw(5) << std::setprecision(3) << _dielectricTensor(a, b) << " ";

    pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);
    return true;
  }

  void ForceConstantCalculator::readBornEffectiveChargesFromAIMSOUT(void) {
    string function = "ForceConstantCalculator::readBornEffectiveChargesFromAIMSOUT()";
    string message = "This functionality has not been implemented yet.";
    throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
  }

  //////////////////////////////////////////////////////////////////////////////
  void ForceConstantCalculator::readBornEffectiveChargesFromOUTCAR(const _xinput& xinp) {  // ME20190113
    string directory = xinp.xvasp.Directory;  // ME20190113

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
      string function = "apl::LinearResponsePC::readBornEffectiveChargesFromOUTCAR()";
      string message = "Cannot open input file OUTCAR.static.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }

    string line;
    uint line_count = 0;  //CO
    string KEY; //ME20181226
    if (DEFAULT_APL_USE_LEPSILON) { //ME20181226
      KEY = string("BORN EFFECTIVE CHARGES (in e, cummulative output)");//ME20181226
    } else { //ME20181226
      KEY = string("BORN EFFECTIVE CHARGES (including local field effects)"); //ME20181226
    } //ME20181226

    while (true) {
      // Get line
      //CO - START
      if (line_count == vlines.size()) {
        string function = "apl::LinearResponsePC::readBornEffectiveChargesFromOUTCAR()";
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
    pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);
    for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++) {
      int id = i;
      message << "Atom [" << aurostd::PaddedNumString(id, 3) << "] ("
        << std::setw(2) << _supercell->getInputStructure().atoms[id].cleanname
        << ") Born effective charge = ";
      pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);
      for (int a = 1; a <= 3; a++)
        for (int b = 1; b <= 3; b++)
          message << std::fixed << std::setw(5) << std::setprecision(3) << _bornEffectiveChargeTensor[i](a, b) << " ";
      pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);
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
    pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);

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
      pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);
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
      string function = "apl::LinearResponsePC::readDielectricTensorFromOUTCAR()";
      string message = "Cannot open input file OUTCAR.";
      aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
    //CO - END

    // Find
    string KEY; //ME20181226
    if (DEFAULT_APL_USE_LEPSILON) { //ME20181226
      KEY = string("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)"); //ME20181226
    } else { //ME20181226
      KEY = string("MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects)"); //ME20181226
    }//ME20181226

    while (true) {
      // Get line
      //CO - START
      if (line_count == vlines.size()) {
        string function = "apl::LinearResponsePC::readDielectricTensorFromOUTCAR()";
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

  void ForceConstantCalculator::hibernate() {
    string base = _aflowFlags->Directory + "/" + DEFAULT_APL_FILE_PREFIX;
    string filename = aurostd::CleanFileName(base + DEFAULT_APL_HARMIFC_FILE);
    string message = "Writing harmonic IFCs into file " + filename + ".";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);
    writeHarmonicIFCs(filename);
    if (_isPolarMaterial) {
      filename = aurostd::CleanFileName(base + DEFAULT_APL_POLAR_FILE);
      message = "Writing harmonic IFCs into file " + filename + ".";
      writeBornChargesDielectricTensor(filename);
      pflow::logger(_AFLOW_FILE_NAME_, _APL_FCCALC_MODULE_, message, *_aflowFlags, *messageFile, *oss);
    }
  }

  void ForceConstantCalculator::writeHarmonicIFCs(const string& filename) {
    stringstream outfile;
    string tab = " ";

    // Header
    outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
    outfile << "<apl>" << std::endl;
    outfile << tab << "<generator>" << std::endl;
    string time = aflow_get_time_string();
    if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
    outfile << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
    outfile << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_ << "\" type=\"" << APL_CHECKSUM_ALGO << "\">"
      << std::hex << aurostd::getFileCheckSum(_aflowFlags->Directory + "/" + _AFLOWIN_, APL_CHECKSUM_ALGO) << "</i>" << std::endl;  // ME20190219
    outfile.unsetf(std::ios::hex); //CO20190116 - undo hex immediately
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
      string function = "ForceConstantCalculator::writeHarmonicIFCs()";
      string message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
    }
  }

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
    outfile << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_ << "\" type=\"" << APL_CHECKSUM_ALGO << "\">"
      << std::hex << aurostd::getFileCheckSum(_aflowFlags->Directory + "/" + _AFLOWIN_ + "", APL_CHECKSUM_ALGO) << "</i>" << std::endl;  // ME20190219
    outfile.unsetf(std::ios::hex); //CO20190116 - undo hex immediately
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

  void ForceConstantCalculator::printForceConstantMatrices(ostream& os) {
    int units = 1;
    double conversionFactor = 1.0;

    switch (units) {
      case (1):
        os << "FORCE CONSTANT MATRICES in eV/A^2:" << std::endl;
        conversionFactor = 1.0;
        break;
      case (2):
        os << "FORCE CONSTANT MATRICES in 10 Dyn/cm:" << std::endl;
        conversionFactor = 1602.17733;
        break;
      case (3):
        os << "FORCE CONSTANT MATRICES in N/m:" << std::endl;
        conversionFactor = 16.0217733;
        break;
    }
    os << std::endl;

    for (int i = 0; i < _supercell->getNumberOfAtoms(); i++) {
      for (int k = 0; k < _supercell->getNumberOfAtoms(); k++) {
        os << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        os << setprecision(4);
        os << "- MATRIX: " << i + 1 << "/" << k + 1 << " " << k + 1 << "/" << i + 1 << std::endl;
        for (int m = 1; m <= 3; m++) {
          for (int n = 1; n <= 3; n++)
            os << setw(10) << (conversionFactor * _forceConstantMatrices[k][i](m, n)) << " ";
          os << " ";
          for (int n = 1; n <= 3; n++)
            os << setw(10) << (conversionFactor * _forceConstantMatrices[i][k](n, m)) << " ";
          os << std::endl;
        }
        os << std::endl;
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void ForceConstantCalculator::printFCShellInfo(ostream& os) {
    int units = 4;
    double conversionFactor = 1.0;

    switch (units) {
      case (1):
        os << "FORCE CONSTANT MATRICES in eV/A^2:" << std::endl;
        conversionFactor = 1.0;
        break;
      case (2):
        os << "FORCE CONSTANT MATRICES in 10 Dyn/cm:" << std::endl;
        conversionFactor = 1602.17733;
        break;
      case (3):
        os << "FORCE CONSTANT MATRICES in N/m:" << std::endl;
        conversionFactor = 16.0217733;
        break;
      case (4):
        os << "FORCE CONSTANT MATRICES in 10^3 Dyn/cm:" << std::endl;
        conversionFactor = 16.0217733;
        break;
    }
    os << std::endl;

    int maxshell = _supercell->getMaxShellID();
    if (maxshell == -1) maxshell = 25;
    std::vector<ShellHandle> sh;
    for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
      ShellHandle s;
      sh.push_back(s);
      sh.back().init(_supercell->getInputStructure(),
          _supercell->getInputStructure().iatoms[i][0],
          maxshell);
      sh[i].splitBySymmetry();
      sh[i].mapStructure(_supercell->getSupercellStructure(), _supercell->getUniqueAtomID(i));
    }

    //
    for (int i = 0; i < _supercell->getNumberOfUniqueAtoms(); i++) {
      sh[i].printReport(cout);
      for (int ishell = 0; ishell <= sh[i].getLastOccupiedShell(); ishell++) {
        for (int isubshell = 0; isubshell < sh[i].getNumberOfSubshells(ishell); isubshell++) {
          const deque<_atom>& atomsAtSameShell = sh[i].getAtomsAtSameShell(ishell, isubshell);
          cout << "SHELL " << ishell << " " << isubshell << std::endl;

          for (uint ai = 0; ai < atomsAtSameShell.size(); ai++) {
            int nb = atomsAtSameShell[ai].number;
            cout << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
            cout << setprecision(4);
            cout << "- MATRIX: " << i << "/" << nb << " " << nb << "/" << i << std::endl;
            for (int m = 1; m <= 3; m++) {
              for (int n = 1; n <= 3; n++)
                cout << setw(10) << (conversionFactor * _forceConstantMatrices[nb][i](m, n)) << " ";
              cout << " ";
              for (int n = 1; n <= 3; n++)
                cout << setw(10) << (conversionFactor * _forceConstantMatrices[i][nb](n, m)) << " ";
              cout << std::endl;
            }
            cout << std::endl;
          }
        }
      }
    }

    // Clear
    for (uint i = 0; i < sh.size(); i++)
      sh[i].clear();
    sh.clear();
  }

}  // namespace apl

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
//****************************************************************************
