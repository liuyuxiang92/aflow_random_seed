#include "aflow_apl.h"

#define _DEBUG_APL_PHONCALC_ false  //CO190116
#define ERROR_VERBOSE false  // CO

using namespace std;

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  PhononCalculator::PhononCalculator(Supercell& sc,
      _xinput& xinput, _aflags& aflags, _kflags& kflags,
      _xflags& xflags, //_vflags& vflags, 
      string& AflowIn, Logger& l)
    : _xInput(xinput), _aflowFlags(aflags), _kbinFlags(kflags), _xFlags(xflags),
    _AflowIn(AflowIn), _logger(l), _supercell(sc) {
      DISTORTION_MAGNITUDE = 0.015;
      _isGammaEwaldPrecomputed = false;
      //    DOtar = false;  OBSOLETE - ME 181024
      // ME190614 - Add system for VASP-style output files
      if ((_xFlags.AFLOW_MODE_VASP) && (!_xFlags.vflags.AFLOW_SYSTEM.content_string.empty())) {
        _system = _xFlags.vflags.AFLOW_SYSTEM.content_string;
      } else {
        _system = _supercell.getInputStructure().title;
      }
      zerostate_dir = "";  // ME191030
    }

  // ME191228 - BEGIN
  // Copy constructors
  PhononCalculator& PhononCalculator::operator=(const PhononCalculator& that) {
    if (this != &that) copy(that);
    return *this;
  }

  void PhononCalculator::copy(const PhononCalculator& that) {
    DISTORTION_MAGNITUDE = that.DISTORTION_MAGNITUDE;
    DISTORTION_INEQUIVONLY = that.DISTORTION_INEQUIVONLY;
    TCOND = that.TCOND;
    anharmonic_IFC_options = that.anharmonic_IFC_options;
    _xInput = that._xInput;
    _aflowFlags = that._aflowFlags;
    _xFlags = that._xFlags;
    _AflowIn = that._AflowIn;
    _system = that._system;
    _logger = that._logger;
    _supercell = that._supercell;
    xInputs = that.xInputs;
    _uniqueDistortions = that._uniqueDistortions;
    _uniqueForces = that._uniqueForces;
    _forceConstantMatrices = that._forceConstantMatrices;
    _calculateZeroStateForces = that._calculateZeroStateForces;
    _bornEffectiveChargeTensor = that._bornEffectiveChargeTensor;
    _dielectricTensor = that._dielectricTensor;
    _inverseDielectricTensor = that._inverseDielectricTensor;
    _recsqrtDielectricTensorDeterminant = that._recsqrtDielectricTensorDeterminant;
    _isGammaEwaldPrecomputed = that._isGammaEwaldPrecomputed;
    _gammaEwaldCorr = that._gammaEwaldCorr;
    zerostate_dir = that.zerostate_dir;
  }
  // ME191228 - END

  // ///////////////////////////////////////////////////////////////////////////

  PhononCalculator::~PhononCalculator() {
    clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::clear() {
    for (uint i = 0; i < _uniqueDistortions.size(); i++)
      _uniqueDistortions[i].clear();
    _uniqueDistortions.clear();

    for (uint i = 0; i < _uniqueForces.size(); i++) {
      for (uint j = 0; j < _uniqueForces[i].size(); j++)
        _uniqueForces[i][j].clear();
      _uniqueForces[i].clear();
    }
    _uniqueForces.clear();

    for (uint i = 0; i < _forceConstantMatrices.size(); i++)
      _forceConstantMatrices[i].clear();
    _forceConstantMatrices.clear();

    _isGammaEwaldPrecomputed = false;
    _gammaEwaldCorr.clear();

    _bornEffectiveChargeTensor.clear();
    zerostate_dir = "";  // ME191030
  }

  // INTERFACE /////////////////////////////////////////////////////////////////

  const Supercell& PhononCalculator::getSupercell() { //CO 180409
    return _supercell;
  }

  const xstructure& PhononCalculator::getInputCellStructure() {
    //     cerr << _supercell.getInputStructure() << std::endl;
    return _supercell.getInputStructure();
  }

  const xstructure& PhononCalculator::getSuperCellStructure() {
    return _supercell.getSupercellStructure();
  }

  //CO - START
  double PhononCalculator::getEPS() {
    return _supercell.getEPS();
  }
  //CO - END

  uint PhononCalculator::getNumberOfBranches() {
    return (3 * _supercell.getInputStructure().atoms.size());
  }

  // ME190614
  string PhononCalculator::getSystemName() {
    return _system;
  }

  // ME200206
  bool PhononCalculator::isPolarMaterial() {
    return _isPolarMaterial;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::run() {
    // Check if supercell is already built
    if (!_supercell.isConstructed()) {
      // ME191031 - use xerror
      //throw APLRuntimeError("apl::PhononCalculator::run(); The supercell structure has not been initialized yet.");
      string function = "apl::PhononCalculator::run()";
      string message = "The supercell structure has not been initialized yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }

    calculateForceConstants();

    // ME191219 - atomGoesTo and atomComesFrom can now use basis_atoms_map.
    // Calculating the full basis ahead of time is much faster than calculating all
    // symmetry operations on-the-fly.
    if (!_supercell.fullBasisCalculatedAGROUP()) _supercell.getFullBasisAGROUP();

    //cout << "NON-SYMMETRIZED:" << std::endl;
    //printForceConstantMatrices(cout);

    // Symmetrization of the force-constant matrices
    symmetrizeForceConstantMatrices();

    //cout << "AFTER SYMMETRIZATION:" << std::endl;
    //printForceConstantMatrices(cout);

    // Force the force-constant matrices to obey the sum-rule conditions
    correctSumRules();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::symmetrizeForceConstantMatrices() {
    bool LDEBUG=(FALSE || _DEBUG_APL_PHONCALC_ || XHOST.DEBUG);
    string soliloquy="apl::PhononCalculator::symmetrizeForceConstantMatrices()"; //CO190218
    // Test of stupidity...
    if (!_supercell.getSupercellStructure().agroup_calculated) {
      // ME191031 - use xerror
      //throw APLRuntimeError("apl::PhononCalculator::symmetrizeForceConstantMatrices(); The site groups have not been calculated yet.");
      string function = "apl::PhononCalculator::symmetrizeForceConstantMatrices()";
      string message = "The site groups have not been calculated yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
    }
    //CO - START
    if (_supercell.getEPS() == AUROSTD_NAN) {
      // ME191031 - use xerror
      //throw APLRuntimeError("apl::PhononCalculator::symmetrizeForceConstantMatrices(); Need to define symmetry tolerance.");
      string function = "apl::PhononCalculator::symmetrizeForceConstantMatrices()";
      string message = "Need to define symmetry tolerance.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ERROR_);
    }
    //CO - END

    //
    _logger << "Symmetrizing the force constant matrices." << apl::endl;

    // Get site symmetry group
    //const vector< vector<_sym_op> >& agroup = _supercell.getSupercellStructure().agroup; //JAHNATEK ORIGINAL
    //[CO190218 - moved down]const vector<vector<_sym_op> >& agroup = _supercell.getAGROUP();  //CO

    //
    vector<xmatrix<double> > row;
    //[CO190218 - moved down]xmatrix<double> m(3, 3);
    //[CO190218 - moved down]m.clear();
    //[CO181226 - OBSOLETE]uint agroup_size;  //CO
    for (int i = 0; i < _supercell.getNumberOfAtoms(); i++) {
      const vector<_sym_op>& agroup = _supercell.getAGROUP(i);  //CO //CO190218
      if (agroup.size() == 0) {
        // ME191031 - use xerror
        //throw APLRuntimeError("apl::PhononCalculator::symmetrizeForceConstantMatrices(); Site point group operations are missing.");
        string function = "apl::PhononCalculator::symmetrizeForceConstantMatrices()";
        string message = "Site point group operations are missing.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_INIT_);
      }

      // Translate the center to this atom
      //_supercell.center(i);  // OBSOLETE ME191219

      //
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
        //[CO181226 - OBSOLETE]agroup_size = agroup.size();  //CO
        if(LDEBUG){ //CO190218
          cerr << soliloquy << " compare original m=" << std::endl;
          cerr << _forceConstantMatrices[i][j] << std::endl;
        }
        xmatrix<double> m(3, 3); //CO190218
        for (uint symOpID = 0; symOpID < agroup.size(); symOpID++) {
          const _sym_op& symOp = agroup[symOpID];

          try {
            //_AFLOW_APL_REGISTER_ int l = _supercell.atomComesFrom(symOp, j, i, FALSE);  //CO NEW //CO190218
            // ME191219 - atomGoesTo now uses basis_atoms_map; keep translation option in case
            // the basis has not been calculated for some reason
            _AFLOW_APL_REGISTER_ int l = _supercell.atomGoesTo(symOp, j, i, true); //JAHNATEK ORIGINAL //CO190218
            //cout << "Mapping " << j << " <-> " << l << std::endl;
            m = m + (inverse(symOp.Uc) * _forceConstantMatrices[i][l] * symOp.Uc);  //JAHNATEK ORIGINAL //CO190218
            //m = m + (symOp.Uc * _forceConstantMatrices[i][l] * inverse(symOp.Uc));  //CO NEW //CO190218
            if(LDEBUG){ //CO190218
              cerr << soliloquy << " atom[" << l << "].cpos=" << _supercell.getSupercellStructure().atoms[l].cpos << std::endl;
              cerr << soliloquy << " atom[" << j << "].cpos=" << _supercell.getSupercellStructure().atoms[j].cpos << std::endl;
              cerr << soliloquy << " agroup(" << l << " -> " << j << ")=" << std::endl;
              cerr << symOp.Uc << std::endl;
              cerr << soliloquy << " forceConstantMatrices[i=" << i << "][l=" << l << "]=" << std::endl;
              cerr << _forceConstantMatrices[i][l] << std::endl;
              cerr << soliloquy << " with new m=" << std::endl;
              cerr << (symOp.Uc * _forceConstantMatrices[i][l] * inverse(symOp.Uc)) << std::endl;
            }
            //CO - START
          }
          // ME191031 - use xerror
          //catch (APLLogicError& e)
          catch (aurostd::xerror& e)
          { //CO200106 - patching for auto-indenting
            //_logger << error << "Mapping problem " << j << " <-> ?. Skipping." << apl::endl;
            //derivative structures are expected to lose symmetry, don't bother exiting
            //CO181226 - forget about this junk
            //if it's a derivative structure, we recalculate the symmetry for the supercell, it's necessary
            //[CO181226 - OBSOLETE]if (!_supercell.isDerivativeStructure()) {
            _logger << error << "Mapping problem " << j << " <-> ?." << apl::endl;
            // ME191031 - use xerror
            //throw APLLogicError("apl::PhononCalculator::symmetrizeForceConstantMatrices(); Mapping failed.");
            string function = "apl::PhononCalculator::symmetrizeForceConstantMatrices()";
            string message = "Mapping failed.";
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
            //[CO181226 - OBSOLETE]}
            //[CO181226 - OBSOLETE]agroup_size -= 1;  //CO, reduce agroup size
            //CO - END
          }
        }
        m = ( 1.0 / agroup.size() ) * m; //CO190218
        //CO - START
        //CO181226 - forget about this junk
        //if it's a derivative structure, we recalculate the symmetry for the supercell, it's necessary
        //[CO181226 - OBSOLETE]if (agroup_size) {
        //[CO181226 - OBSOLETE]m = (1.0 / agroup_size) * m;  //CO
        //[CO181226 - OBSOLETE]}
        //CO - END
        row.push_back(m);
        //m.clear();  //JAHNATEK ORIGINAL //CO190218
      }
      _forceConstantMatrices[i] = row;
      row.clear();

      // Translate the center back
      //_supercell.center_original();  //CO  // OBSOLETE ME191219
    }

    // Translate the center back
    //_supercell.center(0); //JAHNATEK ORIGINAL
    //[CO181226 - not necessary]_supercell.center_original();  //CO
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::correctSumRules() {
    xmatrix<double> sum(3, 3), sum2(3, 3);

    for (int i = 0; i < _supercell.getNumberOfAtoms(); i++) {
      // Get SUMs
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
        if (i == j) continue;
        //if( !TRUNCATE_DYNAMICAL_MATRIX || ( rshell <= maxShellRadiusOfType[isc1] + _AFLOW_APL_EPS_ ) )
        {
          sum = sum + _forceConstantMatrices[i][j];
          sum2 = sum2 + trasp(_forceConstantMatrices[j][i]);
        }
      }

      //cout << "SUM RULE 1:" << std::endl;
      //printXMatrix(_forceConstantMatrices[i][i]);
      //printXMatrix(-1.0*sum);
      //cout << "SUM RULE 2:" << std::endl;
      //printXMatrix(sum);
      //printXMatrix(sum2);

      // Correct SUM2
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
        if (i == j) continue;
        _forceConstantMatrices[i][j] = 0.5 * (_forceConstantMatrices[i][j] + trasp(_forceConstantMatrices[j][i]));
        _forceConstantMatrices[j][i] = trasp(_forceConstantMatrices[i][j]);
      }

      // Get SUMs again
      sum.clear();
      sum2.clear();
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
        if (i == j) continue;
        //if( !TRUNCATE_DYNAMICAL_MATRIX || ( rshell <= maxShellRadiusOfType[isc1] + _AFLOW_APL_EPS_ ) )
        {
          sum = sum + _forceConstantMatrices[i][j];
          sum2 = sum2 + trasp(_forceConstantMatrices[j][i]);
        }
      }

      // Correct SUM1 to satisfied
      _forceConstantMatrices[i][i] = -sum;

      //cout << "cSUM RULE 1:" << std::endl;
      //printXMatrix(_forceConstantMatrices[i][i]);
      //printXMatrix(-1.0*sum);
      //cout << "cSUM RULE 2:" << std::endl;
      //printXMatrix(sum);
      //printXMatrix(sum2);
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::printForceConstantMatrices(ostream& os) {
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

    for (int i = 0; i < _supercell.getNumberOfAtoms(); i++)
      //for(int ii = 0; ii < _supercell.getNumberOfUniqueAtoms(); ii++)
    {
      //int i = _supercell.getUniqueAtomID(ii);
      for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
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

  void PhononCalculator::printFCShellInfo(ostream& os) {
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

    int maxshell = _supercell.getMaxShellID();
    if (maxshell == -1) maxshell = 25;
    std::vector<ShellHandle> sh;
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      ShellHandle s;
      sh.push_back(s);
      sh.back().init(_supercell.getInputStructure(),
          _supercell.getInputStructure().iatoms[i][0],
          maxshell);
      sh[i].splitBySymmetry();
      sh[i].mapStructure(_supercell.getSupercellStructure(), _supercell.getUniqueAtomID(i));
    }

    //
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
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
            //cout << "atom " << setw(3) << nb << ": "; printXVector(atomsAtSameShell[ai].cpos);
            //cout << "atom " << setw(3) << i << ": "; printXVector(_supercell.getSupercellStructure().atoms[_supercell.getUniqueAtomID(i)].cpos);
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

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::hibernate() {
    _logger << "Hibernating..." << apl::endl;

    //
    //CO - START
    stringstream outfile;
    //ofstream outfile("apl.xml", ios_base::out);
    //if (!outfile.is_open())
    //  throw apl::APLRuntimeError("PhononCalculator::hibernate(); Cannot open output file.");
    //CO - END

    // XML declaration
    outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;

    // Our structure
    string tab = " ";
    outfile << "<apl>" << std::endl;

    // Info about calculation run
    outfile << tab << "<generator>" << std::endl;
    string time = aflow_get_time_string();
    if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
    outfile << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
    outfile << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_ << "\" type=\"" << APL_CHECKSUM_ALGO << "\">"
      << std::hex << aurostd::getFileCheckSum(_aflowFlags.Directory + "/" + _AFLOWIN_ + "", APL_CHECKSUM_ALGO) << "</i>" << std::endl;  // ME190219
    outfile.unsetf(std::ios::hex); //CO190116 - undo hex immediately
    outfile << tab << "</generator>" << std::endl;

    if (_uniqueDistortions.size() > 0) {  // ME200211 - linear response has no distortions
      // Unique distortions
      outfile << tab << "<distortions units=\"Angstrom\" cs=\"cartesian\">" << std::endl;

      outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      outfile << setprecision(8);
      outfile << tab << tab << "<i name=\"magnitude\">" << setw(15) << DISTORTION_MAGNITUDE << "</i>" << std::endl;

      outfile << tab << tab << "<varray>" << std::endl;
      for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
        outfile << tab << tab << tab << "<varray atomID=\"" << _supercell.getUniqueAtomID(i) << "\">" << std::endl;
        for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
          outfile << tab << tab << tab << tab << "<v>";
          outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
          outfile << setprecision(8);
          outfile << setw(15) << _uniqueDistortions[i][j](1) << " ";
          outfile << setw(15) << _uniqueDistortions[i][j](2) << " ";
          outfile << setw(15) << _uniqueDistortions[i][j](3);
          outfile << "</v>" << std::endl;
        }
        outfile << tab << tab << tab << "</varray>" << std::endl;
      }
      outfile << tab << tab << "</varray>" << std::endl;
      outfile << tab << "</distortions>" << std::endl;

      // Forces
      outfile << tab << "<forcefields units=\"eV/Angstrom\" cs=\"cartesian\">" << std::endl;
      outfile << tab << tab << "<varray>" << std::endl;
      for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
        outfile << tab << tab << tab << "<varray atomID=\"" << _supercell.getUniqueAtomID(i) << "\">" << std::endl;
        for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
          outfile << tab << tab << tab << tab << "<varray distortion=\"" << j << "\">" << std::endl;
          for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
            outfile << tab << tab << tab << tab << tab << "<v>";
            outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
            outfile << setprecision(15);
            outfile << setw(24) << std::scientific << _uniqueForces[i][j][k](1) << " ";
            outfile << setw(24) << std::scientific << _uniqueForces[i][j][k](2) << " ";
            outfile << setw(24) << std::scientific << _uniqueForces[i][j][k](3);
            outfile << "</v>" << std::endl;
          }
          outfile << tab << tab << tab << tab << "</varray>" << std::endl;
        }
        outfile << tab << tab << tab << "</varray>" << std::endl;
      }
      outfile << tab << tab << "</varray>" << std::endl;
      outfile << tab << "</forcefields>" << std::endl;
    }

    // Force constant matrices
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

    //_logger << "APL-DEBUG  Only for polar materials" << apl::endl;
    // Only for polar materials
    if (_isPolarMaterial) {
      _logger << "APL-DEBUG Born effective charge tensors" << apl::endl;
      // Born effective charge tensors
      outfile << tab << "<born units=\"a.u.\" cs=\"cartesian\">" << std::endl;
      outfile << tab << tab << "<varray>" << std::endl;
      for (uint i = 0; i < _bornEffectiveChargeTensor.size(); i++) {
        int id = _supercell.getUniqueAtomID(i);
        outfile << tab << tab << tab << "<matrix type=\"" << _supercell.getSupercellStructure().atoms[id].cleanname << "\">" << std::endl;
        for (int k = 1; k <= 3; k++) {
          outfile << tab << tab << tab << tab << "<v>";
          for (int l = 1; l <= 3; l++) {
            outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
            outfile << setprecision(8);
            // ME 181030 - fixed prevents hexadecimal output
            outfile << setw(15) << std::fixed << _bornEffectiveChargeTensor[i](k, l) << " ";
          }
          outfile << "</v>" << std::endl;
        }
        outfile << tab << tab << tab << "</matrix>" << std::endl;
      }
      outfile << tab << tab << "</varray>" << std::endl;
      outfile << tab << "</born>" << std::endl;

      // Dielectric constant matrix
      outfile << tab << "<epsilon units=\"a.u.\" cs=\"cartesian\">" << std::endl;
      outfile << tab << tab << "<matrix>" << std::endl;
      for (int k = 1; k <= 3; k++) {
        outfile << tab << tab << tab << "<v>";
        for (int l = 1; l <= 3; l++) {
          outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
          outfile << setprecision(8);
          // ME 181030 - fixed prevents hexadecimal output
          outfile << setw(15) << std::fixed << _dielectricTensor(k, l) << " ";
        }
        outfile << "</v>" << std::endl;
      }
      outfile << tab << tab << "</matrix>" << std::endl;
      outfile << tab << "</epsilon>" << std::endl;
    }

    //
    outfile << "</apl>" << std::endl;

    //
    // outfile.clear();
    // outfile.close();

    // COREY, KBIN_ZIP will compress the whole directory, so just leave it alone
    // if (DOtar) {
    //   aurostd::stringstream2compressfile(kbinFlags.KZIP_BIN,outfile,"apl.xml");
    //   if (!aurostd::EFileExist("apl.xml"))
    //    throw apl::APLRuntimeError("PhononCalculator::hibernate(); Cannot open output apl.xml.");
    //} else {
    string filename = aurostd::CleanFileName(_aflowFlags.Directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_HARMIFC_FILE); //ME181226
    aurostd::stringstream2file(outfile, filename); //ME181226
    if (!aurostd::FileExist(filename)) { //ME181226
      string function = "PhononCalculator::hibernate()";
      string message = "Cannot open output file " + filename + "."; //ME181226
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
      //      throw apl::APLRuntimeError("PhononCalculator::hibernate(); Cannot open output apl.xml.");
    }
    //}

    // Compress
    //if (DOtar) aurostd::execute(string("COMPRESS -fz9 apl.xml"));
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::awake(string hibfile, bool compare_checksum) {
    _logger << "Awakening..." << apl::endl;

    //CO, we already checked that it exists before, just open

    hibfile = aurostd::CleanFileName(hibfile); //ME200220
    vector<string> vlines;                           //CO
    aurostd::efile2vectorstring(hibfile, vlines);  //CO //ME181226
    // Decompress
    //bool isXMLCompressed = aurostd::FileExist(string("apl.xml.EXT")); //CO
    //if (isXMLCompressed) //CO
    //  aurostd::execute(string("EXT -d apl.xml.EXT")); //CO

    //
    //CO - START
    //ifstream infile("apl.xml", ios_base::in);
    //if (!infile.is_open())
    if (!vlines.size()) {
      // ME191031 - use xerror
      string function = "PhononCalculator::awake()";
      string message = "Cannot open output file " + hibfile + "."; //ME181226
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
      //      throw apl::APLRuntimeError("apl::PhononCalculator::awake(); Cannot open input apl.xml.");
    }

    string line;
    uint line_count = 0;
    vector<string> tokens;

    // Test of xml...
    line = vlines[line_count++];
    //getline(infile, line);
    if (line.find("xml") == string::npos) {
      // ME191031 - use xerror
      //throw APLLogicError("apl::PhononCalculator::awake(); Wrong xml file.");
      string function = "apl::PhononCalculator::awake()";
      string message = "Not an xml file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_WRONG_FORMAT_);
    }
    //CO - END

    if (compare_checksum) {
      // Get _AFLOWIN_ checksum and compare it to current
      while (true) {
        //getline(infile, line); //CO
        //if (infile.eof()) //CO
        if (line_count == vlines.size())  //CO
          // ME191031 - use xerror
          //throw APLLogicError("apl::PhononCalculator::awake(); Can not find <i name=\"checksum\" ...> tag.");
          throw aurostd::xerror(_AFLOW_FILE_NAME_, "apl::PhononCalculator::awake()", "Can not find <i name=\"checksum\" ...> tag.", _FILE_CORRUPT_);
        line = vlines[line_count++];  //CO
        if (line.find("checksum") != string::npos)
          break;
      }
      int t = line.find_first_of(">") + 1;
      tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
      if (strtoul(tokens[0].c_str(), NULL, 16) != aurostd::getFileCheckSum(_aflowFlags.Directory + "/" + _AFLOWIN_ + "", APL_CHECKSUM_ALGO)) {  // ME190219
        // ME191031 - use xerror
        //throw APLLogicError("apl::PhononCalculator::awake(); The " + _AFLOWIN_ + " file has been changed from the hibernated state.");
        string function = "apl::PhononCalculator::awake()";
        string message = "The " + _AFLOWIN_ + " file has been changed from the hibernated state.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
      }
    }
    tokens.clear();

    // Get force constant matrices
    while (true) {
      //getline(infile, line); //CO
      //if (infile.eof()) //CO
      if (line_count == vlines.size()) { //CO
        // ME191031 - use xerror
        //throw APLLogicError("apl::PhononCalculator::awake(); Can not find <fcms> tag.");
        string function = "apl::PhononCalculator::awake()";
        string message = "Cannot find <fcms> tag.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];  //CO
      if (line.find("fcms") != string::npos)
        break;
    }
    //CO - START
    line = vlines[line_count++];
    //getline(infile, line);
    line = vlines[line_count++];
    //getline(infile, line);
    //CO - END
    vector<xmatrix<double> > row;
    xmatrix<double> m(3, 3);
    while (true) {
      //getline(infile, line); //CO
      //if (infile.eof()) //CO
      if (line_count == vlines.size()) { //CO
        // ME191031 - use xerror
        //throw APLLogicError("apl::PhononCalculator::awake(); Incomplete <fcms> tag.");
        string function = "apl::PhononCalculator::awake()";
        string message = "Incomplete <fcms> tag.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];  //CO
      if (line.find("</varray>") != string::npos) {
        _forceConstantMatrices.push_back(row);
        row.clear();
        //getline(infile, line); //CO
        line = vlines[line_count++];  //CO
        if (line.find("</varray>") != string::npos)
          break;
        //getline(infile, line); //CO
        line = vlines[line_count++];  //CO
      }

      for (int k = 1; k <= 3; k++) {
        //getline(infile, line); //CO
        line = vlines[line_count++];  //CO
        int t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        m(k, 1) = aurostd::string2utype<double>(tokens.at(0));
        m(k, 2) = aurostd::string2utype<double>(tokens.at(1));
        m(k, 3) = aurostd::string2utype<double>(tokens.at(2));
        tokens.clear();
      }
      row.push_back(m);
      //getline(infile, line); //CO
      line = vlines[line_count++];  //CO
    }

    // Try to read born effective charges and dielectric constant
    try {  //CO
      //cerr << "APL-DEBUG Get born effective charge tensors" << std::endl; //CO
      // Get born effective charge tensors
      while (true) {
        //getline(infile, line); //CO
        //if (infile.eof()) //CO
        if (line_count == vlines.size())  //CO
        {  //CO200106 - patching for auto-indenting
          _isPolarMaterial = false;
          // ME191031 - use xerror
          //throw APLLogicError("apl::PhononCalculator::awake(); Can not find <born> tag.");
          string function = "apl::PhononCalculator::awake()";
          string message = "Cannot find <born> tag.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
        }
        line = vlines[line_count++];  //CO
        if (line.find("born") != string::npos)
          break;
      }
      //getline(infile, line); //CO
      line = vlines[line_count++];  //CO
      while (true) {
        //getline(infile, line); //CO
        //if (infile.eof()) //CO
        if (line_count == vlines.size()) { //CO
          // ME191031 - use xerror
          //throw APLLogicError("apl::PhononCalculator::awake(); Incomplete <born> tag.");
          string function = "apl::PhononCalculator::awake()";
          string message = "Incomplete <born> tag.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
        }
        line = vlines[line_count++];  //CO
        if (line.find("</varray>") != string::npos)
          break;
        for (int k = 1; k <= 3; k++) {
          //getline(infile, line); //CO
          line = vlines[line_count++];  //CO
          int t = line.find_first_of(">") + 1;
          tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
          m(k, 1) = aurostd::string2utype<double>(tokens.at(0));
          m(k, 2) = aurostd::string2utype<double>(tokens.at(1));
          m(k, 3) = aurostd::string2utype<double>(tokens.at(2));
          tokens.clear();
        }
        _bornEffectiveChargeTensor.push_back(m);
        //getline(infile, line); //CO
        line = vlines[line_count++];  //CO
      }

      // Get dielectric constant tensor
      while (true) {
        //getline(infile, line); //CO
        //if (infile.eof()) //CO
        if (line_count == vlines.size())  //CO
        {  //CO200106 - patching for auto-indenting
          _isPolarMaterial = false;
          // ME191031 - use xerror
          //throw APLLogicError("apl::PhononCalculator::awake(); Can not find <epsilon> tag.");
          string function = "apl::PhononCalculator::awake()";
          string message = "Can not find <epsilon> tag.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
        }
        line = vlines[line_count++];  //CO
        if (line.find("epsilon") != string::npos)
          break;
      }
      //getline(infile, line); //CO
      line = vlines[line_count++];  //CO
      for (int k = 1; k <= 3; k++) {
        //getline(infile, line); //CO
        line = vlines[line_count++];  //CO
        int t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        _dielectricTensor(k, 1) = aurostd::string2utype<double>(tokens.at(0));
        _dielectricTensor(k, 2) = aurostd::string2utype<double>(tokens.at(1));
        _dielectricTensor(k, 3) = aurostd::string2utype<double>(tokens.at(2));
        tokens.clear();
      }
      _inverseDielectricTensor = inverse(_dielectricTensor);
      _recsqrtDielectricTensorDeterminant = 1.0 / sqrt(determinant(_dielectricTensor));
      // ME191031 - use xerror
    }  //CO200106 - patching for auto-indenting
    //catch (APLLogicError& e)  //CO
    catch (aurostd::xerror& e)  //CO
    {  //CO200106 - patching for auto-indenting
      //_logger << apl::warning << e.what() << apl::endl;
    }

    //
    //infile.close(); //CO
    //infile.clear(); //CO

    // Compress
    //CO - START
    //if (DOtar)
    //  if (isXMLCompressed)
    //    aurostd::execute(string("COMPRESS -fz9 apl.xml"));
    //CO - END
  }

  // ///////////////////////////////////////////////////////////////////////////

  //CO 180214 - START

  void PhononCalculator::get_NCPUS(string& ncpus) {
    ncpus="MAX";
    if(_kbinFlags.KBIN_MPI_NCPUS>0){ncpus=aurostd::utype2string(_kbinFlags.KBIN_MPI_NCPUS);}
    if(XHOST.vflag_control.flag("XPLUG_NUM_THREADS")){ncpus=XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS");}
  }

  void PhononCalculator::get_NCPUS(int& ncpus) {
    string ncpus_str;
    get_NCPUS(ncpus_str);
    if(ncpus_str=="MAX"){ncpus=MPI_NCPUS_MAX;return;}
    ncpus=aurostd::string2utype<int>(ncpus_str);
    if(ncpus<1){ncpus=1;}
  }

  //CO 180214 - STOP

  // ///////////////////////////////////////////////////////////////////////////

  // ME180827 - Overloaded to calculate derivative and eigenvectors for AAPL
  // ME200206 - Added variants for the case near the Gamma point where the
  // non-analytical correction also needs a direction.
  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const IPCFreqFlags& flags) {
    return getFrequency(kpoint, kpoint, flags);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const xvector<double>& kpoint_nac, const IPCFreqFlags& flags) {
    const xstructure& pc = _supercell.getInputStructureLight();  //CO
    uint nBranches = 3 * pc.atoms.size();
    xmatrix<xcomplex<double> > placeholder_eigen(nBranches, nBranches, 1, 1);
    return getFrequency(kpoint, kpoint_nac, flags, placeholder_eigen);
  }

  // ME190624 - get eigenvectors and frequencies
  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const IPCFreqFlags& flags,
      xmatrix<xcomplex<double> >& eigenvectors) {
    return getFrequency(kpoint, kpoint, flags, eigenvectors);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const xvector<double>& kpoint_nac,
      const IPCFreqFlags& flags, xmatrix<xcomplex<double> >& eigenvectors) {
    vector<xmatrix<xcomplex<double> > > placeholder_mat;
    return getFrequency(kpoint, kpoint_nac, flags, eigenvectors, placeholder_mat, false);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint,
      const IPCFreqFlags& flags, xmatrix<xcomplex<double> >& eigenvectors,
      vector<xmatrix<xcomplex<double> > >& dDynMat, bool calc_derivative) {
    return getFrequency(kpoint, kpoint, flags, eigenvectors, dDynMat, calc_derivative);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const xvector<double>& kpoint_nac,
      const IPCFreqFlags& flags, xmatrix<xcomplex<double> >& eigenvectors,
      vector<xmatrix<xcomplex<double> > >& dDynMat, bool calc_derivative) {
    // Compute frequency(omega) from eigenvalues [in eV/A/A/atomic_mass_unit]
    xvector<double> omega = getEigenvalues(kpoint, kpoint_nac, eigenvectors, dDynMat, calc_derivative);

    // Get value of conversion factor
    double conversionFactor = getFrequencyConversionFactor(apl::RAW | apl::OMEGA, flags);

    // Transform values to desired format
    for (_AFLOW_APL_REGISTER_ int i = omega.lrows; i <= omega.urows; i++) {
      if (omega(i) < 0) {
        if (flags & ALLOW_NEGATIVE)
          omega(i) = -sqrt(-omega(i));
        else
          omega(i) = 0.0;
      } else {
        omega(i) = sqrt(omega(i));
      }

      // Convert to desired units
      omega(i) *= conversionFactor;
    }

    // Return
    return (omega);
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME200108 - replaced with constants in xscalar
  double PhononCalculator::getFrequencyConversionFactor(IPCFreqFlags inFlags, IPCFreqFlags outFlags) {
    double conversionFactor = 1.0;

    // Conversion from eV/A/A/atomic_mass_unit -> something
    if (inFlags & apl::RAW) {
      if (outFlags & apl::RAW) {
        // Transform to eV/A/A/atomic_mass_unit
        conversionFactor = 1.0;
      } else if (outFlags & apl::HERTZ) {
        // Transform to s-1; sqrt(EV_TO_JOULE / (ANGSTROM_TO_METER*ANGSTROM_TO_METER) / AMU_TO_KG);
        conversionFactor = au2Hz;
      } else if (outFlags & apl::THZ) {
        // Transform to THz; (in Hertz) / 1E12;
        conversionFactor = au2Hz * Hz2THz;
      } else if (outFlags & apl::RECIPROCAL_CM) {
        // Transform to cm-1; 1/lambda(m) = freq.(s-1) / light_speed(m/s);
        conversionFactor = au2rcm;
      } else if (outFlags & apl::MEV) {
        // Transform to meV; E(eV) = h(eV.s) * freq(s-1); h[(from J.s) -> (eV.s)] = 4.1356673310E-15
        conversionFactor = 1000 * au2eV;
      } else {
        // ME191031 - use xerror
        //throw APLRuntimeError("apl::PhononCalculator:convertFrequencyUnit(); Not implemented conversion.");
        string function = "apl::PhononCalculator:convertFrequencyUnit()";
        string message = "Not implemented conversion.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
      }
    }

    // Conversion from THz -> something
    else if (inFlags & apl::THZ) {
      if (outFlags & apl::RAW) {
        // Transform to eV/A/A/atomic_mass_unit
        conversionFactor = 1.0 / (au2Hz * Hz2THz);
      } else if (outFlags & apl::THZ) {
        conversionFactor = 1.0;
      } else if (outFlags & apl::MEV) {
        conversionFactor = 1000 * PLANCKSCONSTANTEV_h * THz2Hz;
      } else {
        // ME191031 - use xerror
        //throw APLRuntimeError("apl::PhononCalculator:convertFrequencyUnit(); Not implemented conversion.");
        string function = "apl::PhononCalculator:convertFrequencyUnit()";
        string message = "Not implemented conversion.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
      }
    }

    // Nothing suits?
    else {
      // ME191031 - use xerror
      //throw APLRuntimeError("apl::PhononCalculator:convertFrequencyUnit(); Not implemented conversion.");
      string function = "apl::PhononCalculator:convertFrequencyUnit()";
      string message = "Not implemented conversion.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
    }

    //
    if ((outFlags & OMEGA) && !(inFlags & OMEGA))
      conversionFactor *= 2.0 * M_PI;
    if (!(outFlags & OMEGA) && (inFlags & OMEGA))
      conversionFactor /= 2.0 * M_PI;

    //
    return (conversionFactor);
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME180827 - Overloaded to calculate derivative and eigenvectors for AAPL
  // OBSOLETE ME200206 - not used anywhere and notuseful for debugging (use get Frequency)
  //[OBSOLETE] xvector<double> PhononCalculator::getEigenvalues(const xvector<double>& kpoint) {
  //[OBSOLETE]   const xstructure& pc = _supercell.getInputStructureLight();  //CO
  //[OBSOLETE]   uint nBranches = 3 * pc.atoms.size();
  //[OBSOLETE]   xmatrix<xcomplex<double> > placeholder_eigen(nBranches, nBranches, 1, 1);
  //[OBSOLETE]   vector<xmatrix<xcomplex<double> > > placeholder_mat;
  //[OBSOLETE]   return getEigenvalues(kpoint, kpoint, placeholder_eigen, placeholder_mat, false);
  //[OBSOLETE] }

  xvector<double> PhononCalculator::getEigenvalues(const xvector<double>& kpoint,
      const xvector<double>& kpoint_nac,
      xmatrix<xcomplex<double> >& eigenvectors,
      vector<xmatrix<xcomplex<double> > >& dDynMat,
      bool calc_derivative) {
    // Get dynamical matrix
    xmatrix<xcomplex<double> > dynamicalMatrix = getDynamicalMatrix(kpoint, kpoint_nac, dDynMat, calc_derivative);

    // Diagonalize
    xvector<double> eigenvalues(dynamicalMatrix.rows, 1);
    //    xmatrix<xcomplex<double> > unitaryMatrix;  OBSOLETE ME 180827

    // OBSOLETE ME190815 - moved to aurostd::xmatrix
    //#ifdef USE_MKL
    //    zheevMKL(dynamicalMatrix, eigenvalues, eigenvectors);
    //#else
    //    //tred2(dynamicalMatrix);
    //    zheevByJacobiRotation(dynamicalMatrix2, eigenvalues2, eigenvectors2);
    //    eigenvectors2 = trasp(eigenvectors2);
    //#endif

    // ME 180828; OBSOLETE ME190815 - use Jacobi algorithm in aurostd::xmatrix, which
    // is much, much faster than aplEigensystems for large systems
    //    apl::aplEigensystems e;
    //    e.eigen_calculation(dynamicalMatrix, eigenvalues, eigenvectors, APL_MV_EIGEN_SORT_VAL_ASC);

    eigenvalues = jacobiHermitian(dynamicalMatrix, eigenvectors);  // ME190815

    return eigenvalues;
  }

  //  // ///////////////////////////////////////////////////////////////////////////
  // ME180827 - Overloaded to calculate derivative for AAPL
  // ME200206 - Added variants for the case near the Gamma point where the
  // non-analytical correction also needs a direction. While dynamical matrices
  // are not used directly, these functions are helpful debugging tools.
  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint) {
    return getDynamicalMatrix(kpoint, kpoint);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint, const xvector<double>& kpoint_nac) {
    vector<xmatrix<xcomplex<double> > > placeholder;
    return getDynamicalMatrix(kpoint, kpoint_nac, placeholder, false);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint,
      const xvector<double>& kpoint_nac,
      vector<xmatrix<xcomplex<double> > >& dDynMat,
      bool calc_derivative) {
    uint scAtomsSize = _supercell.getSupercellStructure().atoms.size();
    uint pcAtomsSize = _supercell.getInputStructure().atoms.size();

    uint nBranches = 3 * pcAtomsSize;
    xmatrix<xcomplex<double> > dynamicalMatrix(nBranches, nBranches, 1, 1);
    xmatrix<xcomplex<double> > dynamicalMatrix0(nBranches, nBranches, 1, 1);

    xcomplex<double> phase;
    double value;
    // ME 180828 - Prepare derivative calculation
    xvector<xcomplex<double> > derivative(3);
    vector<xmatrix<xcomplex<double> > > dDynMat_NAC;
    if (calc_derivative) {  // reset dDynMat
      dDynMat.clear();
      xmatrix<xcomplex<double> > mat(nBranches, nBranches, 1, 1);
      dDynMat.assign(3, mat);
    }

    // Calculate nonanalytical contribution
    xmatrix<xcomplex<double> > dynamicalMatrixNA(nBranches, nBranches, 1, 1);
    if (_isPolarMaterial)
      dynamicalMatrixNA = getNonanalyticalTermWang(kpoint_nac, dDynMat_NAC, calc_derivative);

    // Loop over primitive cell
    xcomplex<double> nac;
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      uint isc1 = _supercell.pc2scMap(ipc1);

      for (uint isc2 = 0; isc2 < scAtomsSize; isc2++) {
        uint ipc2 = _supercell.sc2pcMap(isc2);
        int neq;  // Important for NAC derivative
        if (_supercell.calcShellPhaseFactor(isc2, isc1, kpoint, phase, neq, derivative, calc_derivative)) {  // ME180827
          for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
            for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
              value = 0.5 * (_forceConstantMatrices[isc1][isc2](ix, iy) + _forceConstantMatrices[isc2][isc1](iy, ix));
              dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += value * phase;
              if (_isPolarMaterial)
                dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += dynamicalMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy) * phase;
              dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy) += value;
              if (_isPolarMaterial)
                dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy) += dynamicalMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy);
              if (calc_derivative) {
                for (int d = 0; d < 3; d++) {
                  dDynMat[d](3 * ipc1 + ix, 3 * ipc2 + iy) += value * derivative[d+1];
                  if (_isPolarMaterial && (aurostd::modulus(kpoint) > _AFLOW_APL_EPS_)) {
                    nac = ((double) neq) * phase * dDynMat_NAC[d](3 * ipc1 + ix, 3 * ipc2 + iy);
                    dDynMat[d](3 * ipc1 + ix, 3 * ipc2 + iy) += nac;
                  }
                }
              }
            }
          }
        }
      }
    }
    //printXMatrix2(dynamicalMatrix);

    // Subtract the sum of all "forces" from the central atom, this is like an automatic sum rule...
    for (uint i = 0; i < pcAtomsSize; i++) {
      for (uint j = 0; j < pcAtomsSize; j++) {
        for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
          for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
            dynamicalMatrix(3 * i + ix, 3 * i + iy) = dynamicalMatrix(3 * i + ix, 3 * i + iy) - dynamicalMatrix0(3 * i + ix, 3 * j + iy);
          }
        }
      }
    }

    // Get correction for polar materials
    //if( _isPolarMaterial )
    // dynamicMatrix += getNonanalyticalTermGonze(kpoint);

    // Make it hermitian
    for (uint i = 0; i <= pcAtomsSize - 1; i++) {
      for (uint j = 0; j <= i; j++) {
        for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
          for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
            dynamicalMatrix(3 * i + ix, 3 * j + iy) += conj(dynamicalMatrix(3 * j + iy, 3 * i + ix));
            dynamicalMatrix(3 * i + ix, 3 * j + iy) *= 0.5;
            dynamicalMatrix(3 * j + iy, 3 * i + ix) = conj(dynamicalMatrix(3 * i + ix, 3 * j + iy));
            if (calc_derivative) {
              for (int d = 0; d < 3; d++) {
                dDynMat[d](3 * i + ix, 3 * j + iy) += conj(dDynMat[d](3 * j + iy, 3 * i + ix));
                dDynMat[d](3 * i + ix, 3 * j + iy) *= 0.5;
                dDynMat[d](3 * j + iy, 3 * i + ix) = conj(dDynMat[d](3 * i + ix, 3 * j + iy));
              }
            }
          }
        }
      }
    }

    // Divide by masses
    for (uint i = 0; i < pcAtomsSize; i++) {
      double mass_i = _supercell.getAtomMass(_supercell.pc2scMap(i));
      for (uint j = 0; j < pcAtomsSize; j++) {
        double mass_j = _supercell.getAtomMass(_supercell.pc2scMap(j));
        for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
          for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
            dynamicalMatrix(3 * i + ix, 3 * j + iy) *= 1.0 / sqrt(mass_i * mass_j);
            if (calc_derivative) {
              for (int d = 0; d < 3; d++) {
                dDynMat[d](3 * i + ix, 3 * j + iy) *= 1.0/sqrt(mass_i * mass_j);
              }
            }
          }
        }
      }
    }

    return dynamicalMatrix;
  }

 ///////////////////////////////////////////////////////////////////////////

  // Y. Wang et.al, J. Phys.:Condens. Matter 22, 202201 (2010)
  // DOI: 10.1088/0953-8984/22/20/202201

  // ME180827 - Overloaded to calculate derivative for AAPL
  // ME200207 - This function assummed that Born charges were stored for each type,
  // but it is actually stored for each iatom.
  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermWang(const xvector<double>& _q) {
    vector<xmatrix<xcomplex<double> > > placeholder;
    return getNonanalyticalTermWang(_q, placeholder, false);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermWang(const xvector<double>& _q,
      vector<xmatrix<xcomplex<double> > >& derivative,
      bool calc_derivative) {
    const xstructure& sc = _supercell.getSupercellStructureLight();           //CO
    const xstructure& pc = _supercell.getInputStructure();  //CO  // ME200207 - grab input structure (need iatoms)

    // to correct the q=\Gamma as a limit
    xvector<double> q(_q);
    if (aurostd::modulus(q) < _AFLOW_APL_EPS_) {
      q(1) = _AFLOW_APL_EPS_ * 1.001;
    }

    uint pcAtomsSize = pc.atoms.size();
    uint pcIAtomsSize = pc.iatoms.size();
    uint nBranches = 3 * pcAtomsSize;

    xmatrix<xcomplex<double> > dynamicalMatrix(nBranches, nBranches);

    if (calc_derivative) {  // reset derivative
      derivative.clear();
      xmatrix<xcomplex<double> > mat(nBranches, nBranches, 1, 1);
      derivative.assign(3, mat);
    }

    // Calculation
    double fac0 = hartree2eV * bohr2angst;  // from a.u. to eV/A  // ME200206 - replaced with xscalar constants
    double volume = det(pc.lattice);
    double fac1 = 4.0 * PI / volume;
    double nbCells = det(sc.lattice) / volume;


    if (aurostd::modulus(q) > _AFLOW_APL_EPS_) {
      // Precompute product of q-point with charge tensor
      vector<xvector<double> > qZ(pcIAtomsSize);
      for (uint at = 0; at < pcIAtomsSize; at++) qZ[at] = q * _bornEffectiveChargeTensor[at];

      double dotprod = scalar_product(q, _dielectricTensor * q);
      double prefactor = fac0 * fac1/(dotprod * nbCells);
      for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
        int iat1 = pc.atoms[ipc1].index_iatoms;
        for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
          int iat2 = pc.atoms[ipc2].index_iatoms;
          for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
            for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
              //int typei = pc.atoms[ipc1].type;
              //int typej = pc.atoms[ipc2].type;
              //double borni = (q * _bornEffectiveChargeTensor[typei])(ix);
              //double bornj = (q * _bornEffectiveChargeTensor[typej])(iy);
              double borni = qZ[iat1][ix];
              double bornj = qZ[iat2][iy];
              dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) = prefactor * borni * bornj;
              if (calc_derivative) {
                for (int d = 0; d < 3; d++) {
                  xcomplex<double> coeff(0, 0);
                  //coeff += borni * _bornEffectiveChargeTensor[ipc1](iy, d + 1);
                  //coeff += bornj * _bornEffectiveChargeTensor[ipc2](ix, d + 1);
                  coeff += borni * _bornEffectiveChargeTensor[iat1](iy, d + 1);
                  coeff += bornj * _bornEffectiveChargeTensor[iat2](ix, d + 1);
                  coeff -= 2 * borni * bornj * scalar_product(_dielectricTensor(d + 1), q)/dotprod;
                  derivative[d](3 * ipc1 + ix, 3 * ipc2 + iy) = prefactor * coeff;
                }
              }
            }
          }
        }
      }
    }

    //
    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // X. Gonze et al., Phys. Rev. B 50, 13035 (1994)
  // X. Gonze and Ch. Lee, Phys. Rev. B 55, 10355 (1997)

  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermGonze(const xvector<double> kpoint) {
    uint pcAtomsSize = _supercell.getInputStructure().atoms.size();
    //    uint nBranches = 3 * pcAtomsSize; // not needed

    if (!_isGammaEwaldPrecomputed) {
      xvector<double> zero(3);
      xmatrix<xcomplex<double> > dynamicalMatrix0(getEwaldSumDipolDipolContribution(zero, false));

      _gammaEwaldCorr.clear();
      for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
        xmatrix<xcomplex<double> > sum(3, 3);
        for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
          for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
            for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
              sum(ix, iy) += dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy);
        }
        _gammaEwaldCorr.push_back(sum);
      }

      _isGammaEwaldPrecomputed = true;
    }

    //
    xmatrix<xcomplex<double> > dynamicalMatrix(getEwaldSumDipolDipolContribution(kpoint));

    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
        for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
          dynamicalMatrix(3 * ipc1 + ix, 3 * ipc1 + iy) -= _gammaEwaldCorr[ipc1](ix, iy);
    }

    //
    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME200207 - This function assummed that Born charges were stored for each type,
  // but it is actually stored for each iatom.
  xmatrix<xcomplex<double> > PhononCalculator::getEwaldSumDipolDipolContribution(const xvector<double> qpoint, bool includeTerm1) {
    // Definitions
    const xstructure& sc = _supercell.getSupercellStructureLight();           //CO
    const xstructure& pc = _supercell.getInputStructure();  //CO  // ME200207 - grab input structure (need iatoms)

    uint pcAtomsSize = pc.atoms.size();
    uint nBranches = 3 * pcAtomsSize;

    xmatrix<xcomplex<double> > dynamicalMatrix(nBranches, nBranches);

    double gmax = 14.0;
    double lambda = 1.0;
    double lambda2 = lambda * lambda;
    double lambda3 = lambda2 * lambda;
    double geg = gmax * lambda2 * 4.0;

    // Reciprocal Space
    xmatrix<double> klattice = trasp(ReciprocalLattice(pc.lattice));

    // Grid
    int n1 = (int)(sqrt(geg) / aurostd::modulus(klattice(1))) + 1;
    int n2 = (int)(sqrt(geg) / aurostd::modulus(klattice(2))) + 1;
    int n3 = (int)(sqrt(geg) / aurostd::modulus(klattice(3))) + 1;

    // Calculation
    double fac0 = hartree2eV * bohr2angst;  // from a.u. to eV/A  // ME200207 - replaced with xscalar constants
    double SQRTPI = sqrt(PI);
    double volume = det(pc.lattice);
    double fac = 4.0 * PI / volume;
    xcomplex<double> iONE(0.0, 1.0);

    // Term 1 - Reciprocal space sum

    if (includeTerm1) {
      for (_AFLOW_APL_REGISTER_ int m1 = -n1; m1 <= n1; m1++) {
        for (_AFLOW_APL_REGISTER_ int m2 = -n2; m2 <= n2; m2++) {
          for (_AFLOW_APL_REGISTER_ int m3 = -n3; m3 <= n3; m3++) {
            xvector<double> g = m1 * klattice(1) + m2 * klattice(2) + m3 * klattice(3) + qpoint;

            geg = scalar_product(g, _dielectricTensor * g);

            if (aurostd::abs(geg) > _AFLOW_APL_EPS_ && geg / lambda2 / 4.0 < gmax) {
              double fac2 = fac * exp(-geg / lambda2 / 4.0) / geg;

              for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
                //xvector<double> zag = g * _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
                int iat1 = pc.atoms[ipc1].index_iatoms;
                xvector<double> zag = g * _bornEffectiveChargeTensor[iat1];

                for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
                  //xvector<double> zbg = g * _bornEffectiveChargeTensor[pc.atoms[ipc2].type];
                  int iat2 = pc.atoms[ipc2].index_iatoms;
                  xvector<double> zbg = g * _bornEffectiveChargeTensor[iat2];

                  //xcomplex<double> e;
                  //(void)_supercell.calcShellPhaseFactor(ipc2,ipc1,g,e);
                  //xcomplex<double> facg = fac2 * e;
                  xcomplex<double> facg = fac2 * exp(iONE * scalar_product(g, sc.atoms[ipc2].cpos - sc.atoms[ipc1].cpos));

                  for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
                    for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
                      dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += fac0 * facg * zag(ix) * zbg(iy);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // Term 2 - Real space sum
    //for(_AFLOW_APL_REGISTER_ int m1 = -n1; m1 <= n1; m1++)
    //  for(_AFLOW_APL_REGISTER_ int m2 = -n2; m2 <= n2; m2++)
    //    for(_AFLOW_APL_REGISTER_ int m3 = -n2; m3 <= n3; m3++) {
    //      xvector<double> rc = m1 * pc.lattice(1) + m2 * pc.lattice(2)
    //        + m3 * pc.lattice(3);

    //      //xvector<double> zero(3);
    //      //xvector<double> rf = _supercell.getFPositionItsNearestImage(rc,zero,pc.lattice);
    //      //rc = F2C(pc.lattice,rf);

    //      if( aurostd::modulus(rc) < _AFLOW_APL_EPS_ ) continue;

    //      //
    //      xvector<double> delta = _inverseDielectricTensor * rc;
    //      double D = sqrt( scalar_product(delta,rc) );

    //      //
    //      xmatrix<double> H(3,3);
    //      xvector<double> x = lambda * delta;
    //      double y = lambda * D;
    //      double y2 = y * y;
    //      double ym2 = 1.0 / y2;
    //      double emy2dpi = 2.0 * exp( -y2 ) / SQRTPI;
    //      double erfcdy = erfc(y) / y;
    //      double c1 = ym2 * ( 3.0 * erfcdy * ym2 + ( emy2dpi * ( 3.0 * ym2 + 2.0 ) ) );
    //      double c2 = ym2 * ( erfcdy + emy2dpi );
    //      for(_AFLOW_APL_REGISTER_ int a = 1; a <= 3; a++)
    //        for(_AFLOW_APL_REGISTER_ int b = 1; b <= 3; b++) {
    //          H(a,b) = x(a) * x(b) * c1 - _inverseDielectricTensor(a,b) * c2;
    //        }

    //      //
    //      xcomplex<double> e = exp( iONE * scalar_product(qpoint,rc) );
    //      xcomplex<double> fac = fac0 * lambda3 * _recsqrtDielectricTensorDeterminant * e;

    //      //
    //      for(uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
    //        xmatrix<double> zh = _bornEffectiveChargeTensor[pc.atoms[ipc1].type] * H;

    //        for(uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
    //          xmatrix<double> zhz = zh * _bornEffectiveChargeTensor[pc.atoms[ipc2].type];

    //          for(_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
    //            for(_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
    //              dynamicalMatrix(3*ipc1+ix,3*ipc2+iy) -= fac * zhz(ix,iy);
    //        }
    //      }
    //    }

    // Term 2
    uint scAtomsSize = sc.atoms.size();
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      uint isc1 = _supercell.pc2scMap(ipc1);

      for (uint isc2 = 0; isc2 < scAtomsSize; isc2++) {
        uint ipc2 = _supercell.sc2pcMap(isc2);

        xvector<double> rf = _supercell.getFPositionItsNearestImage(isc2, isc1);
        xvector<double> rc = F2C(sc.lattice, rf);

        if (aurostd::modulus(rc) < _AFLOW_APL_EPS_) continue;

        //
        xvector<double> delta = _inverseDielectricTensor * rc;
        double D = sqrt(scalar_product(delta, rc));

        //
        xmatrix<double> H(3, 3);
        xvector<double> x = lambda * delta;
        double y = lambda * D;
        double y2 = y * y;
        double ym2 = 1.0 / y2;
        double emy2dpi = 2.0 * exp(-y2) / SQRTPI;
        double erfcdy = erfc(y) / y;
        double c1 = ym2 * (3.0 * erfcdy * ym2 + (emy2dpi * (3.0 * ym2 + 2.0)));
        double c2 = ym2 * (erfcdy + emy2dpi);
        for (_AFLOW_APL_REGISTER_ int a = 1; a <= 3; a++) {
          for (_AFLOW_APL_REGISTER_ int b = 1; b <= 3; b++) {
            H(a, b) = x(a) * x(b) * c1 - _inverseDielectricTensor(a, b) * c2;
          }
        }

        //xmatrix<double> za = _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
        //xmatrix<double> zb = _bornEffectiveChargeTensor[pc.atoms[ipc2].type];
        int iat1 = pc.atoms[ipc1].index_iatoms;
        int iat2 = pc.atoms[ipc2].index_iatoms;
        xmatrix<double> za = _bornEffectiveChargeTensor[iat1];
        xmatrix<double> zb = _bornEffectiveChargeTensor[iat2];
        xmatrix<double> zhz = za * H * zb;

        //
        xcomplex<double> e;  // = exp( iONE * scalar_product(qpoint,rc) );
        (void)_supercell.calcShellPhaseFactor(isc2, isc1, qpoint, e);

        //
        xcomplex<double> fac = fac0 * lambda3 * _recsqrtDielectricTensorDeterminant * e;
        for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
          for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
            dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) -= fac * zhz(ix, iy);
      }
    }

    // Term 3 - Limiting contribution

    double facterm3 = fac0 * 4.0 * lambda3 * _recsqrtDielectricTensorDeterminant / (3.0 * SQRTPI);
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      //xmatrix<double> z = _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
      int iat1 = pc.atoms[ipc1].index_iatoms;
      xmatrix<double> z = _bornEffectiveChargeTensor[iat1];
      xmatrix<double> zez = z * _inverseDielectricTensor * z;

      for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
        for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
          dynamicalMatrix(3 * ipc1 + ix, 3 * ipc1 + iy) -= facterm3 * zez(ix, iy);
    }

    //
    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PhononCalculator::readAnharmonicIFCs(string filename, bool compare_checksum) {
    filename = aurostd::CleanFileName(filename);
    string message = "Reading anharmonic IFCs from " + filename + ".";
    pflow::logger(_AFLOW_FILE_NAME_, "AAPL", message, _logger.getOutputStream(), std::cout);
    string function = "apl::PhononCalculator::readAnharmonicIFCs()";
    if (!aurostd::EFileExist(filename)) {
      message = "Could not open file " + filename + ". File not found.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
    }

    vector<string> vlines;
    aurostd::efile2vectorstring(filename, vlines);
    uint nlines = vlines.size();
    if (nlines == 0) {
      message = "Cannot open file " + filename + ". File empty or corrupt.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
    }

    uint line_count = 0;
    string line = vlines[line_count++];

    // Check that this is a valid xml file
    if (line.find("xml") == string::npos) {
      message = "File is not a valid xml file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_WRONG_FORMAT_);
    }

    int t = 0;
    vector<string> tokens;

    // Checksum
    if (compare_checksum) {
      while (true) {
        if (line_count == nlines) {
          message = "Checksum not found in hibernate file.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
        }
        line = vlines[line_count++];
        if (line.find("checksum") != string::npos) {
          t = line.find_first_of(">") + 1;
          aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, " ");
          if (strtoul(tokens[0].c_str(), NULL, 16) != aurostd::getFileCheckSum(_aflowFlags.Directory + "/" + _AFLOWIN_, APL_CHECKSUM_ALGO)) {
            message = "The " + _AFLOWIN_ + " file has been changed from the hibernated state.";
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
          }
          break;
        }
      }
    }

    // Read order
    uint order = 0;
    while (true) {
      if (line_count == nlines) {
        message = "order tag not found.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("order") != string::npos) {
        t = line.find_first_of(">") + 1;
        order = aurostd::string2utype<uint>(line.substr(t, line.find_last_of("<") - t));
      }
    }

    // Read IFCs
    while (true) {
      if (line_count == nlines) {
        message = "force_constants tag not found.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("force_constants") != string::npos) break;
    }

    vector<vector<double> > ifcs;
    vector<vector<int> > clst;
    uint nanharm = anharmonicIFCs.size();
    for (uint i = nanharm; i < order - 2; i++) {
      anharmonicIFCs.push_back(ifcs);
      clusters.push_back(clst);
    }
    vector<double> fc;
    vector<int> cl;
    while (line.find("/force_constants") == string::npos) {
      if (line_count == nlines) {
        message = "force_constants tag incomplete.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
      }
      line = vlines[line_count++];
      if (line.find("atoms") != string::npos) {
        // New tensor and cluster
        fc.clear();
        cl.clear();
        t = line.find_first_of("\"") + 1;
        aurostd::string2tokens(line.substr(t, line.find_last_of("\"") - t), cl, " ");
        clst.push_back(cl);
      } else if (line.find("slice") != string::npos) {
        // Populate tensor
        for (int i = 0; i < 3; i++) {
          tokens.clear();
          line = vlines[line_count++];
          t = line.find_first_of(">") + 1;
          aurostd::string2tokens(line.substr(t, line.find_last_of("<") - t), tokens, " ");
          for (int j = 0; j < 3; j++) fc.push_back(aurostd::string2utype<double>(tokens[j]));
        }
        line_count++;  // Skip /varray
      } else if (line.find("</varray>") != string::npos) {
        ifcs.push_back(fc);
      }
    }

    if (ifcs.size() != clst.size()) {
      message = "Number of IFCs is different from the number of clusters.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
    }

    // Populate
    anharmonicIFCs[order - 2] = ifcs;
    clusters[order - 2] = clst;
  }

}  // namespace apl
