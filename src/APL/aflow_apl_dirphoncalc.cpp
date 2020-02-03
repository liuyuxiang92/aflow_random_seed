#include "aflow_apl.h"

#define _DEBUG_APL_DIRPHONCALC_ false  //CO190116

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  DirectMethodPC::DirectMethodPC(Supercell& sc, vector<ClusterSet>& clst,
      _xinput& xinput, _aflags& aflags, _kflags& kflags,
      _xflags& xflags, //_vflags& vflags, 
      string& AflowIn, //this is the file CONTENTS, not the file path
      Logger& l)
    : PhononCalculator(sc, clst, xinput, aflags, kflags, xflags, AflowIn, l) {
      //GENERATE_PLUS_MINUS = false; //JAHNATEK ORIGINAL
      AUTO_GENERATE_PLUS_MINUS = true;   //CO
      USER_GENERATE_PLUS_MINUS = false;  //CO
      GENERATE_ONLY_XYZ = false;
      DISTORTION_SYMMETRIZE = true;   //CO190116
      DISTORTION_INEQUIVONLY = true;   //CO190116
      _isPolarMaterial = false;
    }

  // ///////////////////////////////////////////////////////////////////////////

  DirectMethodPC::~DirectMethodPC() {
    clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void DirectMethodPC::clear() {
  }

  //////////////////////////////////////////////////////////////////////////////

  void DirectMethodPC::runVASPCalculations(bool zerostate_chgcar) {
    string soliloquy="apl::DirectMethodPC::runVASPCalculations():"; //CO190218

    // Check if supercell is already built
    if (!_supercell.isConstructed()) {
      // ME191029 - use xerror
      //throw APLRuntimeError("apl::DirectMethodPC::calculateForceFields(); The supercell structure has not been initialized yet.");
      string message = "The supercell structure has not been initialized yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, message, _RUNTIME_INIT_);
    }
    _xInput.xvasp.AVASP_arun_mode = "APL";

    // Determine the distortion vectors
    estimateUniqueDistortions(_supercell.getSupercellStructure(), _uniqueDistortions);

    //CO - START
    vvgenerate_plus_minus.clear();  //CO //CO181226  // ME191029
    bool generate_plus_minus;           //CO
    //bool         check_minus_needed = ( AUTO_GENERATE_PLUS_MINUS && !USER_GENERATE_PLUS_MINUS && !_supercell.isDerivativeStructure() );
    //bool check_minus_needed = (AUTO_GENERATE_PLUS_MINUS && !USER_GENERATE_PLUS_MINUS);  OBSOLETE ME 181028 - this overrides DPM=OFF
    //[CO181212]bool check_minus_needed = AUTO_GENERATE_PLUS_MINUS;  // ME 181028
    int ncalcs = 0;  // ME 190107 - total number of calculations for padding
    if (_calculateZeroStateForces) ncalcs++;  // ME190112
    if (_isPolarMaterial) ncalcs++;  // ME190112

    for (uint i = 0; i < _uniqueDistortions.size(); i++) {
      vvgenerate_plus_minus.push_back(vector<bool>(0)); //CO181226
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        //      vvgenerate_plus_minus.push_back(true);  //assume we need plus/minus OBSOLETE ME 181028 - this overrides DPM=OFF
        // ME 190107 - Calculate "need minus" here
        if (AUTO_GENERATE_PLUS_MINUS) {
          vvgenerate_plus_minus.back().push_back(needMinus(i, j, DISTORTION_INEQUIVONLY)); //CO190218
        } else {
          vvgenerate_plus_minus.back().push_back(USER_GENERATE_PLUS_MINUS);  // ME 181028
        }
        if (vvgenerate_plus_minus[i][j]) {
          ncalcs += 2;
        } else {
          ncalcs += 1;
        }
      }
    }
    //CO - END
    // ME 181022 - START
    // Generate calculation directories
    string chgcar_file = "";
    if (zerostate_chgcar) {  // ME191029 - for ZEROSTATE CHGCAR
      zerostate_dir = "ARUN.APL_";
      int index = ncalcs;
      if (_isPolarMaterial) index--;
      zerostate_dir += aurostd::PaddedNumString(index, aurostd::getZeroPadding(ncalcs)) + "_ZEROSTATE";
      chgcar_file = "../" + zerostate_dir + "/CHGCAR.static";
      if (_kbinFlags.KZIP_COMPRESS) chgcar_file += "." + _kbinFlags.KZIP_BIN;
    }

    for (uint i = 0; i < _uniqueDistortions.size(); i++) {
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        //CO - START
        //[CO181212]if (check_minus_needed)
        //[ME190107] if(AUTO_GENERATE_PLUS_MINUS)  //CHECK -
        //[ME190107] {  //CO200106 - patching for auto-indenting
        //[ME190107]  vvgenerate_plus_minus[i][j] = needMinus(i, j);
        //[ME190107]  if (!vvgenerate_plus_minus[i][j]) {_logger << "No negative distortion needed for distortion [atom=" << i << ",direction=" << j << "]." << apl::endl;}
        //[ME190107]}
        generate_plus_minus = vvgenerate_plus_minus[i][j];
        if (AUTO_GENERATE_PLUS_MINUS && !generate_plus_minus) {
          _logger << "No negative distortion needed for distortion [atom=" << i << ",direction=" << j << "]." << apl::endl;
        }
        for (uint k = 0; k < (generate_plus_minus ? 2 : 1); k++) {
          //CO - END
          // Copy settings from common case
          xInputs.push_back(_xInput);
          int idxRun = xInputs.size() - 1;
          int idAtom = (DISTORTION_INEQUIVONLY ? _supercell.getUniqueAtomID(i) : i); //CO190218

          // Create run ID
          // ME190107 - added padding
          string runname = aurostd::PaddedNumString(idxRun + 1, aurostd::getZeroPadding(ncalcs)) + "_";  // ME190112
          runname += "A" + stringify(idAtom) + "D" + stringify(j); //CO190218

          if (generate_plus_minus) {  //CO
            runname = runname + ((k == 0) ? "P" : "M");
          }
          xInputs[idxRun].xvasp.AVASP_arun_runname = runname;
          // Apply the unique distortion to one inequvalent atom
          // This distortion vector is stored in Cartesian form, hence use C2F before applying
          xInputs[idxRun].setXStr(_supercell.getSupercellStructureLight()); //CO faster, only what's necessary here
          xstructure& xstr = xInputs[idxRun].getXStr(); // ME190109 - Declare to make code more legible
          //CO190114 - it is very silly to try to add in fpos
          //add to cpos, then convert to fpos
          xstr.atoms[idAtom].cpos += ((k == 0) ? 1.0 : -1.0) * DISTORTION_MAGNITUDE * _uniqueDistortions[i][j];
          xstr.atoms[idAtom].fpos = C2F(xstr.lattice, xstr.atoms[idAtom].cpos);
          //[CO190114 - OBSOLETE]xInputs[idxRun].getXStr().atoms[idAtom].fpos = xInputs[idxRun].getXStr().atoms[idAtom].fpos + C2F(xInputs[idxRun].getXStr().lattice, ((k == 0) ? 1.0 : -1.0) * DISTORTION_MAGNITUDE * _uniqueDistortions[i][j]);
          //[CO190114 - OBSOLETE]xInputs[idxRun].getXStr().atoms[idAtom].cpos = F2C(xInputs[idxRun].getXStr().lattice,
          //[CO190114 - OBSOLETE]                                                 xInputs[idxRun].getXStr().atoms[idAtom].fpos);

          //clean title //CO181226
          //[CO190131 - moved up]xstructure& xstr = xInputs[idxRun].getXStr(); // ME190109 - Declare to make code more legible
          xstr.title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr.title); //CO181226, ME 190109
          if(xstr.title.empty()){xstr.buildGenericTitle(true,false);} //CO181226, ME 190109
          xstr.title += " APL supercell=" + aurostd::joinWDelimiter(_supercell.scell, 'x'); //ME190109
          xstr.title += " atom=" + stringify(idAtom); //ME190109
          //xstr.title += " distortion=[" + aurostd::RemoveWhiteSpacesFromTheFrontAndBack(stringify(DISTORTION_MAGNITUDE*_uniqueDistortions[i][j])) + "]"; //ME190109 - OBSOLETE ME190112
          std::stringstream distortion; // ME190112 - need stringstream for nicer formatting
          xvector<double> dist_cart = DISTORTION_MAGNITUDE * _uniqueDistortions[i][j];  // ME190112
          distortion << " distortion=["
            << std::setprecision(3) << dist_cart[1] << ","
            << std::setprecision(3) << dist_cart[2] << ","
            << std::setprecision(3) << dist_cart[3] << "]"; // ME190112
          xstr.title += distortion.str();

          // For VASP, use the standardized aflow.in creator
          if (_kbinFlags.AFLOW_MODE_VASP){
            // ME191029
            xInputs[idxRun].xvasp.aopts.flag("APL_FLAG::ZEROSTATE_CHGCAR", zerostate_chgcar);
            if (zerostate_chgcar) {
              xInputs[idxRun].xvasp.aopts.push_attached("APL_FLAG::CHGCAR_FILE", chgcar_file);
            }

            _kbinFlags.KBIN_MPI_AUTOTUNE = true;
            // Change format of POSCAR
            // ME 190228 - OBSOLETE for two reasons:
            // 1. This method is not robust
            // 2. This will be taken care of when the actual POSCAR is generated
            // [OBSOLETE - 190228] if ((!_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_BIN.find("46") != string::npos)) ||
            // [OBSOLETE - 190228]    (_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_MPI_BIN.find("46") != string::npos))) {
            // [OBSOLETE - 190228]  xInputs[idxRun].getXStr().is_vasp5_poscar_format = false;
            // [OBSOLETE - 190228] }
            _stagebreak = (createAflowInPhonons(xInputs[idxRun]) || _stagebreak);
          }
          // For AIMS, use the old method until we have AVASP_populateXAIMS
          if (_kbinFlags.AFLOW_MODE_AIMS) {
            string runname = ARUN_DIRECTORY_PREFIX + "APL_" + stringify(idxRun) + "A" + stringify(_supercell.getUniqueAtomID(i)) + "D" + stringify(j);
            xInputs[idxRun].setDirectory(_xInput.getDirectory() + "/" + runname);
            if (!filesExistPhonons(xInputs[idxRun])) {
              _logger << "Creating " << xInputs[idxRun].getDirectory() << apl::endl;
              createAflowInPhonons(xInputs[idxRun], runname);
            }
          }
        }
      }
    }

    // Add zero state if requested
    if (_calculateZeroStateForces) {
      // Copy settings from common case
      xInputs.push_back(_xInput);
      int idxRun = xInputs.size() - 1;
      // Create run ID //ME181226
      xInputs[idxRun].xvasp.AVASP_arun_runname = aurostd::PaddedNumString(idxRun+1, aurostd::getZeroPadding(ncalcs)) + "_ZEROSTATE"; //ME181226, ME190112

      // Get structure
      xInputs[idxRun].setXStr(_supercell.getSupercellStructureLight()); //CO

      // ME 190108 - Set title
      xInputs[idxRun].getXStr().title=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xInputs[idxRun].getXStr().title); //CO181226
      if(xInputs[idxRun].getXStr().title.empty()){xInputs[idxRun].getXStr().buildGenericTitle(true,false);} //CO181226
      xInputs[idxRun].getXStr().title += " APL supercell=" + aurostd::joinWDelimiter(_supercell.scell, 'x'); //ME190112
      xInputs[idxRun].getXStr().title += " undistorted";
      xInputs[idxRun].xvasp.aopts.flag("APL_FLAG::IS_ZEROSTATE", true);  // ME191029
      // For VASP, use the standardized aflow.in creator
      if(_kbinFlags.AFLOW_MODE_VASP){
        xInputs[idxRun].xvasp.aopts.flag("APL_FLAG::ZEROSTATE_CHGCAR", zerostate_chgcar);
        _stagebreak = (createAflowInPhonons(xInputs[idxRun]) || _stagebreak); //ME181226
      }
      // For AIMS, use the old method until we have AVASP_populateXAIMS //ME181226
      if(_kbinFlags.AFLOW_MODE_AIMS){
        string runname = ARUN_DIRECTORY_PREFIX + "APL_" + stringify(idxRun) + "ZEROSTATE";
        xInputs[idxRun].setDirectory(_xInput.getDirectory() + "/" + runname);
        if (!filesExistPhonons(xInputs[idxRun])) {
          _logger << "Creating " << xInputs[idxRun].getDirectory() << apl::endl;
          createAflowInPhonons(xInputs[idxRun], runname);
        }
      }
    }

    // BEGIN STEFANO
    // Do an additional calculation for polar materials
    if (_isPolarMaterial) {
      // Calc. Born effective charge tensors and dielectric constant matrix
      _xinput xinpBE(_xInput);  // ME190113
      runVASPCalculationsBE(xinpBE, xInputs.size());  // ME190113
      xInputs.push_back(xinpBE);
    }
    // END STEFANO
  }

  void DirectMethodPC::estimateUniqueDistortions(const xstructure& xstr,
      vector<vector<xvector<double> > >& uniqueDistortions) {
    //COREY NOTES ON THIS FUNCTION
    // - this function creates symmetrically unique distortion vectors for each iatom
    // - you can have at most 3 unique (orthogonal) distortions per atom, but probably fewer considering symmetry
    // - distortions are relative to the lattice vectors
    // - we first generate different types of basic distortions (along lattice vectors, diagonal, body diagonal) in fractional coordinates
    // - convert fractional to cartesian, now distortions truly represent those along lattice vectors, etc.
    // - despite that the vector order changes with every loop, these basic distortions are not changed
    // - for each iatom,
    //   - we build three lists:
    //     - allDistortionsOfAtom (horrible name) is simply a list of symmetrically equivalent (orthogonal) distortions
    //       to each basic distortion (by rotations of the crystal) - this is true because we clear allDistortionsOfAtom with every loop,
    //       so original distortionVector is unchanged
    //     - uniqueDistortions does not play a role in this loop right now, so ignore
    //     - testVectorDim is a count of symmetrically equivalent distortions for each basic distortion
    //   - we sort the basic distortions by testVectorDim (largest first), i.e. the basic distortions that generate the highest count of
    //     equivalent distortions go first (you could consider it like the most "natural" choices for distortions based on crystal symmetry)

    // Is there a list of inequivalent atoms?
    if (DISTORTION_INEQUIVONLY && !xstr.iatoms_calculated) { //CO190218
      // ME191031 - use xerror
      //throw APLRuntimeError("apl::DirectMethodPC::estimateUniqueDistortions(); The list of the inequivalent atoms is missing.");
      string function = "apl::DirectMethodPC::estimateUniqueDistortions()";
      string message = "The list of the inequivalent atoms is missing.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }

    // Clear old stuff
    for (uint i = 0; i < uniqueDistortions.size(); i++) {
      uniqueDistortions[i].clear();
    }
    uniqueDistortions.clear();

    // Determine the distortion vectors for this system
    if (xstr.agroup_calculated) {
      // Test directions for distortion (in cartesian coordinates)
      vector<xvector<double> > testDistortions;
      xvector<double> testVec(3);

      if (GENERATE_ONLY_XYZ) {
        // Test distortion (in cartesian coordinates) along x, y, and z
        // - they have to be orthonormal
        testVec(1) = 1.0;
        testVec(2) = 0.0;
        testVec(3) = 0.0;
        testDistortions.push_back(testVec);
        testVec(1) = 0.0;
        testVec(2) = 1.0;
        testVec(3) = 0.0;
        testDistortions.push_back(testVec);
        testVec(1) = 0.0;
        testVec(2) = 0.0;
        testVec(3) = 1.0;
        testDistortions.push_back(testVec);
      } else {
        // Add lattice vectors
        testVec(1) = 1.0;
        testVec(2) = 0.0;
        testVec(3) = 0.0;
        testDistortions.push_back(testVec);
        testVec(1) = 0.0;
        testVec(2) = 1.0;
        testVec(3) = 0.0;
        testDistortions.push_back(testVec);
        testVec(1) = 0.0;
        testVec(2) = 0.0;
        testVec(3) = 1.0;
        testDistortions.push_back(testVec);

        // Add diagonal vectors
        testVec(1) = 1.0;
        testVec(2) = 1.0;
        testVec(3) = 0.0;
        testDistortions.push_back(testVec);
        testVec(1) = 1.0;
        testVec(2) = 0.0;
        testVec(3) = 1.0;
        testDistortions.push_back(testVec);
        testVec(1) = 0.0;
        testVec(2) = 1.0;
        testVec(3) = 1.0;
        testDistortions.push_back(testVec);

        // Add body diagonal vectors
        testVec(1) = 1.0;
        testVec(2) = 1.0;
        testVec(3) = 1.0;
        testDistortions.push_back(testVec);

        // Convert to cartesian representation
        for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
          testDistortions[idTestVector] = F2C(xstr.lattice, testDistortions[idTestVector]);
        }

        // Norm to one (so it will be scaled by user magnitude of distortion)
        for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
          double len = aurostd::modulus(testDistortions[idTestVector]);
          testDistortions[idTestVector](1) /= len;
          testDistortions[idTestVector](2) /= len;
          testDistortions[idTestVector](3) /= len;
        }
      }

      // Loop over inequivalent atoms
      vector<xvector<double> > uniqueDistortionsOfAtom;
      vector<xvector<double> > allDistortionsOfAtom;
      for (uint i = 0; i < (DISTORTION_INEQUIVONLY ? xstr.iatoms.size() : xstr.atoms.size()); i++) { //CO190218
        int atomID = (DISTORTION_INEQUIVONLY ? xstr.iatoms[i][0] : i); //CO190218
        //cout << "atomID = " << atomID << std::endl; //CO190218
        if (xstr.agroup[atomID].size() == 0) { //CO190218
          // ME191031 - use xerror
          //throw APLRuntimeError("apl::DirectMethodPC::estimateUniqueDistortions(); Site point group operations are missing.");
          string function = "apl::DirectMethodPC::estimateUniqueDistortions()";
          string message = "Site point group operations are missing.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
        }

        // Loop over test directions vectors - we count the number of unique
        // which can be obtain from the each test vector
        vector<int> testVectorDim;
        for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
          // Get the list of all equivalent distortions
          testDistortion(testDistortions[idTestVector],
              xstr.agroup[atomID], allDistortionsOfAtom,
              uniqueDistortionsOfAtom,true); //CO190114 - maximize distortion impact first
          testVectorDim.push_back(allDistortionsOfAtom.size());
          allDistortionsOfAtom.clear();
          uniqueDistortionsOfAtom.clear();
        }

        // Now, we order all test vectors according to their dimensionality
        // Those with the highest score go first
        for (uint j = 0; j < testDistortions.size() - 1; j++) {
          for (uint k = j + 1; k < testDistortions.size(); k++) {
            if (testVectorDim[k] > testVectorDim[j]) {
              xvector<double> temp = testDistortions[k];
              testDistortions[k] = testDistortions[j];
              testDistortions[j] = temp;
              int itemp = testVectorDim[k];
              testVectorDim[k] = testVectorDim[j];
              testVectorDim[j] = itemp;
            }
          }
        }

        // Now we are going again, but slightly different, we count all
        // generated directions together until the total count is lower than three
        for (uint idTestVector = 0; idTestVector < testDistortions.size(); idTestVector++) {
          testDistortion(testDistortions[idTestVector],
              xstr.agroup[atomID], allDistortionsOfAtom,
              uniqueDistortionsOfAtom,DISTORTION_SYMMETRIZE);  //CO190114 - if DISTORTION_SYMMETRIZE==false, we want all 3 directions (don't integrate equivalent distortions into count)
          if (allDistortionsOfAtom.size() >= 3) break;
        }

        //cout << "XXXXX  Number of unique distortion vectors for atom ["
        //     << atomID << "] = " << uniqueDistortionsOfAtom.size() << std::endl; //CO190218
        uniqueDistortions.push_back(uniqueDistortionsOfAtom);
        // Free useless stuff
        allDistortionsOfAtom.clear();
        uniqueDistortionsOfAtom.clear();
      }
      // Free useless stuff
      testDistortions.clear();
    } else {
      // ME191031
      //throw APLRuntimeError("apl::DirectMethodPC::estimateUniqueDistortions(); The list of the site point group operations is missing.");
      string function = "apl::DirectMethodPC::estimateUniqueDistortions()";
      string message = "The list of the site point group operations is missing.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }

    // Print some information
    int dof = 0;
    for (uint i = 0; i < _uniqueDistortions.size(); i++){ //CO200106 - wrapping with guard
      dof += _uniqueDistortions[i].size();
    }
    _logger << "Found " << dof << " degree(s) of freedom." << apl::endl;
    for (int i = 0; i < (DISTORTION_INEQUIVONLY ? _supercell.getNumberOfUniqueAtoms() : _supercell.getNumberOfAtoms()); i++) { //CO190218
      int id = (DISTORTION_INEQUIVONLY ? _supercell.getUniqueAtomID(i) : i); //CO190218
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        _logger << "Atom [" << sf("%03d") << id << "] (" << sf("%f")
          << sw(2) << _supercell.getSupercellStructure().atoms[id].cleanname
          << ") will be distorted in direction ["
          << sf("%5.3f") << _uniqueDistortions[i][j](1) << ","
          << sf("%5.3f") << _uniqueDistortions[i][j](2) << ","
          << sf("%5.3f") << _uniqueDistortions[i][j](3) << "]." << apl::endl;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  void DirectMethodPC::testDistortion(const xvector<double>& distortionVector,
      const vector<_sym_op>& symPool,
      vector<xvector<double> >& allDistortionsOfAtom,
      vector<xvector<double> >& uniqueDistortionsOfAtom,
      bool integrate_equivalent_distortions) {  //CO190114
    // Test if it is unique direction
    // Use the Gram–Schmidt method for orthogonalizing, if the final vectors
    // is non zero length -> is unique
    xvector<double> testDistortion = distortionVector;
    for (uint k = 0; k < allDistortionsOfAtom.size(); k++) {
      testDistortion = testDistortion - getVectorProjection(testDistortion, allDistortionsOfAtom[k]);
    }
    if (aurostd::modulus(testDistortion) > _AFLOW_APL_EPS_) {
      // Normalize vector
      testDistortion = testDistortion / aurostd::modulus(testDistortion);
      // Store this vector (in Cartesian form!)
      allDistortionsOfAtom.push_back(testDistortion);
      uniqueDistortionsOfAtom.push_back(testDistortion);
      //cout << "new unique distortion vector: " << testDistortion << std::endl;
    }

    if(integrate_equivalent_distortions) { //CO181226
      // Apply all symmetry operations on test vector and generate all
      // unique directions and store them for future comparison
      for (uint iSymOp = 0; iSymOp < symPool.size(); iSymOp++) {
        testDistortion = symPool[iSymOp].Uc * distortionVector;

        for (uint k = 0; k < allDistortionsOfAtom.size(); k++) {
          testDistortion = testDistortion - getVectorProjection(testDistortion, allDistortionsOfAtom[k]);
        }
        if (aurostd::modulus(testDistortion) > _AFLOW_APL_EPS_) {
          testDistortion = testDistortion / aurostd::modulus(testDistortion);
          allDistortionsOfAtom.push_back(testDistortion);
          //cout << "new symmetry generated distortion vector: " << testDistortion << std::endl;
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  //CO - START
  bool DirectMethodPC::needMinus(uint atom_index, uint distortion_index, bool inequiv_only) { //CO190218
    //bool need_minus = true;
    const vector<_sym_op>& agroup = _supercell.getAGROUP( inequiv_only ? _supercell.getUniqueAtomID(atom_index) : atom_index);  //CO190116
    //[CO190131 OBSOLETE]uint atom_index = _supercell.getUniqueAtomID(ineq_atom_indx);
    //[CO190131 OBSOLETE]vector<_sym_op> agroup = _supercell.getAGROUP(atom_index);

    xvector<double> rotated_distortion, distortion_sum;
    //cerr << agroup.size() << std::endl;
    //cerr << "distortion  : " << _uniqueDistortions[atom_index][distortion_index] << std::endl;
    for (uint i = 0; i < agroup.size(); i++) {
      rotated_distortion = agroup[i].Uc * _uniqueDistortions[atom_index][distortion_index];
      //cerr << "rdistortion : " << rotated_distortion << std::endl;
      if (identical(_uniqueDistortions[atom_index][distortion_index], -rotated_distortion, _AFLOW_APL_EPS_)) {  //just mimicking Jahnatek tolerance here
        return FALSE;
        //need_minus = false;
        //break;
      } //else {
      //  continue;
      //}
    }
    return TRUE;
    //return need_minus;
  }
  //CO - END

  //////////////////////////////////////////////////////////////////////////////

  void DirectMethodPC::calculateForceFields() {
    bool LDEBUG=(FALSE || _DEBUG_APL_DIRPHONCALC_ || XHOST.DEBUG);
    string soliloquy="apl::DirectMethodPC::runVASPCalculations():"; //CO190218
    // Extract all forces ////////////////////////////////////////////////////

    //first pass, just find if outfile is found ANYWHERE
    if(!outfileFoundAnywherePhonons(xInputs)){throw APLStageBreak();}

    //second pass, make sure it's everywhere!
    outfileFoundEverywherePhonons(xInputs, _isPolarMaterial);

    // Remove zero state forces if necessary
    if (_calculateZeroStateForces) {
      subtractZeroStateForces(xInputs);
    }

    // Store forces //////////////////////////////////////////////////////////

    bool generate_plus_minus = false;  // ME190129

    int idxRun = 0;
    for (int i = 0; i < (DISTORTION_INEQUIVONLY ? _supercell.getNumberOfUniqueAtoms() : _supercell.getNumberOfAtoms()); i++) {
      vector<vector<xvector<double> > > forcesForOneAtomAndAllDistortions;
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        generate_plus_minus = vvgenerate_plus_minus[i][j];  //CO181226
        vector<xvector<double> > forcefield;
        xvector<double> drift(3);
        for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
          xvector<double> force(3);
          force(1) = xInputs[idxRun].getXStr().qm_forces[k](1);
          force(2) = xInputs[idxRun].getXStr().qm_forces[k](2);
          force(3) = xInputs[idxRun].getXStr().qm_forces[k](3);
          if(LDEBUG) { cerr << soliloquy << " force[idistortion=" << i << ",atom=" << k << ",+]=" << xInputs[idxRun].getXStr().qm_forces[k] << std::endl;} //CO190218
          if (generate_plus_minus) {  //CO
            force(1) = 0.5 * (force(1) - xInputs[idxRun + 1].getXStr().qm_forces[k](1));
            force(2) = 0.5 * (force(2) - xInputs[idxRun + 1].getXStr().qm_forces[k](2));
            force(3) = 0.5 * (force(3) - xInputs[idxRun + 1].getXStr().qm_forces[k](3));
            if(LDEBUG) { cerr << soliloquy << " force[idistortion=" << i << ",atom=" << k << ",-]=" << xInputs[idxRun + 1].getXStr().qm_forces[k] << std::endl;} //CO190218
          }
          if(LDEBUG) { cerr << soliloquy << " force[idistortion=" << i << ",atom=" << k << ",AVG]=" << force << std::endl;} //CO190218
          forcefield.push_back(force);
          drift = drift + force;
        }

        if (generate_plus_minus) {  //CO
          idxRun += 2;
        } else {
          idxRun++;
        }

        // Remove drift
        if(LDEBUG) { cerr << soliloquy << " drift[idistortion=" << i << ",total]=" << drift << std::endl;} //CO190218
        drift(1) = drift(1) / forcefield.size();
        drift(2) = drift(2) / forcefield.size();
        drift(3) = drift(3) / forcefield.size();
        if(LDEBUG) { cerr << soliloquy << " drift[idistortion=" << i << ",AVG]=" << drift << std::endl;} //CO190218
        for (_AFLOW_APL_REGISTER_ uint k = 0; k < forcefield.size(); k++) {
          forcefield[k] = forcefield[k] - drift;
          if(LDEBUG) { cerr << soliloquy << " force[idistortion=" << i << ",atom=" << k << ",-drift]=" << forcefield[k] << std::endl;} //CO190218
        }

        // Store
        forcesForOneAtomAndAllDistortions.push_back(forcefield);
        forcefield.clear();
      }
      _uniqueForces.push_back(forcesForOneAtomAndAllDistortions);
      forcesForOneAtomAndAllDistortions.clear();
    }

    if (_isPolarMaterial) calculateDielectricTensor(xInputs.back());

    //******** BEGIN JJPR: Clean later to use ZERO STATE FORCES with TENSOR  *******
    // Clear useless stuff
    if (!TCOND) {
      for (uint i = 0; i < xInputs.size(); i++) {
        xInputs[i].clear();
      }
      xInputs.clear();
    }
    //****** END JJPR *******
  }

  //////////////////////////////////////////////////////////////////////////////

}  // namespace apl
