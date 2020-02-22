#include "aflow_apl.h"

#define _DEBUG_APL_DIRPHONCALC_ false  //CO190116

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  DirectMethodPC::DirectMethodPC(Supercell& sc,
      _xinput& xinput, _aflags& aflags, _kflags& kflags,
      _xflags& xflags, string& AflowIn, ofstream& mf)
    : ForceConstantCalculator(sc, xinput, aflags, kflags, xflags, AflowIn, mf) {
      free();
    }

  DirectMethodPC::DirectMethodPC(const DirectMethodPC& that)
    : ForceConstantCalculator(that._supercell, that._xInput, that._aflowFlags, that._kbinFlags,
      that._xFlags, that._AflowIn, *that.messageFile) {
    free();
    copy(that);
  }

  DirectMethodPC& DirectMethodPC::operator=(const DirectMethodPC& that) {
    if (this != &that) {
      free();
      copy(that);
    }
    return *this;
  }

  DirectMethodPC::~DirectMethodPC() {
    free();
  }

  void DirectMethodPC::clear(Supercell& sc, _xinput& xinput,
      _aflags& aflags, _kflags& kflags, _xflags& xflags, string& AflowIn, ofstream& mf) {
    DirectMethodPC that(sc, xinput, aflags, kflags, xflags, AflowIn, mf);
    copy(that);
  }

  void DirectMethodPC::copy(const DirectMethodPC& that) {
    AUTO_GENERATE_PLUS_MINUS = that.AUTO_GENERATE_PLUS_MINUS;
    DISTORTION_MAGNITUDE = that.DISTORTION_MAGNITUDE;
    DISTORTION_INEQUIVONLY = that.DISTORTION_INEQUIVONLY;
    DISTORTION_SYMMETRIZE = that.DISTORTION_SYMMETRIZE;
    GENERATE_ONLY_XYZ = that.GENERATE_ONLY_XYZ;
    USER_GENERATE_PLUS_MINUS = that.USER_GENERATE_PLUS_MINUS;
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

  void DirectMethodPC::free() {
    AUTO_GENERATE_PLUS_MINUS = true;   //CO
    DISTORTION_MAGNITUDE = 0.0;
    DISTORTION_INEQUIVONLY = true;   //CO190116
    DISTORTION_SYMMETRIZE = true;   //CO190116
    GENERATE_ONLY_XYZ = false;
    USER_GENERATE_PLUS_MINUS = false;  //CO
    xInputs.clear();
    _bornEffectiveChargeTensor.clear();
    _dielectricTensor.clear();
    _forceConstantMatrices.clear();
    _isPolarMaterial = false;
  }

  //////////////////////////////////////////////////////////////////////////////

  bool DirectMethodPC::runVASPCalculations(bool zerostate_chgcar) {
    string soliloquy="apl::DirectMethodPC::runVASPCalculations():"; //CO190218
    bool stagebreak = false;
    stringstream message;

    // Check if supercell is already built
    if (!_supercell.isConstructed()) {
      // ME191029 - use xerror
      //throw APLRuntimeError("apl::DirectMethodPC::calculateForceFields(); The supercell structure has not been initialized yet.");
      message << "The supercell structure has not been initialized yet.";
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
          message << "No negative distortion needed for distortion [atom=" << i << ",direction=" << j << "].";
          pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);
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
            stagebreak = (createAflowInPhonons(_aflowFlags, _kbinFlags, _xFlags, xInputs[idxRun]) || stagebreak);
          }
          // For AIMS, use the old method until we have AVASP_populateXAIMS
          if (_kbinFlags.AFLOW_MODE_AIMS) {
            string runname = ARUN_DIRECTORY_PREFIX + "APL_" + stringify(idxRun) + "A" + stringify(_supercell.getUniqueAtomID(i)) + "D" + stringify(j);
            xInputs[idxRun].setDirectory(_xInput.getDirectory() + "/" + runname);
            if (!filesExistPhonons(xInputs[idxRun])) {
              message << "Creating " << xInputs[idxRun].getDirectory();
              pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);
              createAflowInPhononsAIMS(_aflowFlags, _kbinFlags, _xFlags, _AflowIn, xInputs[idxRun], *messageFile);
              stagebreak = true;
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
        stagebreak = (createAflowInPhonons(_aflowFlags, _kbinFlags, _xFlags, xInputs[idxRun]) || stagebreak); //ME181226
      }
      // For AIMS, use the old method until we have AVASP_populateXAIMS //ME181226
      if(_kbinFlags.AFLOW_MODE_AIMS){
        string runname = ARUN_DIRECTORY_PREFIX + "APL_" + stringify(idxRun) + "ZEROSTATE";
        xInputs[idxRun].setDirectory(_xInput.getDirectory() + "/" + runname);
        if (!filesExistPhonons(xInputs[idxRun])) {
          message << "Creating " << xInputs[idxRun].getDirectory();
          pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);
          createAflowInPhononsAIMS(_aflowFlags, _kbinFlags, _xFlags, _AflowIn, xInputs[idxRun], *messageFile);
          stagebreak = true;
        }
      }
    }

    // BEGIN STEFANO
    // Do an additional calculation for polar materials
    if (_isPolarMaterial) {
      // Calc. Born effective charge tensors and dielectric constant matrix
      _xinput xinpBE(_xInput);  // ME190113
      stagebreak = (runVASPCalculationsBE(xinpBE, xInputs.size()) || stagebreak);  // ME190113
      xInputs.push_back(xinpBE);
    }
    return stagebreak;
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

    stringstream message;
    string function = "apl::DirectMethodPC::estimateUniqueDistortions()";
    // Is there a list of inequivalent atoms?
    if (DISTORTION_INEQUIVONLY && !xstr.iatoms_calculated) { //CO190218
      message << "The list of the inequivalent atoms is missing.";
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
      uint natoms = DISTORTION_INEQUIVONLY ? xstr.iatoms.size() : xstr.atoms.size();
      for (uint i = 0; i < natoms; i++) {
        int atomID = (DISTORTION_INEQUIVONLY ? xstr.iatoms[i][0] : i); //CO190218
        //cout << "atomID = " << atomID << std::endl; //CO190218
        if (xstr.agroup[atomID].size() == 0) { //CO190218
          message << "Site point group operations are missing.";
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
      message << "The list of the site point group operations is missing.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }

    // Print some information
    int dof = 0;
    for (uint i = 0; i < _uniqueDistortions.size(); i++){ //CO200106 - wrapping with guard
      dof += _uniqueDistortions[i].size();
    }
    message << "Found " << dof << " degree(s) of freedom.";
    pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);
    uint natoms = DISTORTION_INEQUIVONLY ? _supercell.getNumberOfUniqueAtoms() : _supercell.getNumberOfAtoms();
    for (uint i = 0; i < natoms; i++) {  //CO200212 - int->uint
      uint id = (DISTORTION_INEQUIVONLY ? _supercell.getUniqueAtomID(i) : i); //CO190218
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        message << "Atom [" << aurostd::PaddedNumString(id, 3) << "] ("
          << std::setw(2) << _supercell.getSupercellStructure().atoms[id].cleanname
          << ") will be distorted in direction ["
          << std::setw(5) << std::setprecision(3) << _uniqueDistortions[i][j](1) << ","
          << std::setw(5) << std::setprecision(3) << _uniqueDistortions[i][j](2) << ","
          << std::setw(5) << std::setprecision(3) << _uniqueDistortions[i][j](3) << "].";
        pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);
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

  void DirectMethodPC::calculateForceConstants() {
    // Get all forces required for the construction of force-constant matrices
    calculateForceFields();

    // ME191219 - atomGoesTo and atomComesFrom can now use basis_atoms_map.
    // Calculating the full basis ahead of time is much faster than calculating all
    // symmetry operations on-the-fly.
    if (!_supercell.fullBasisCalculatedAGROUP()) _supercell.getFullBasisAGROUP();

    // For construction of the force-constant matrices we need three
    // independent distortions. Hence, calculate the remaining distortions and
    // forces by the symmetry (if necessary).
    completeForceFields();

    // Ensure that all distortion vectors are along the cartesian directions
    projectToCartesianDirections();

    // Construct the matrix of force-constant matrices for all atoms based
    // on the force fields for the inequivalent atoms
    buildForceConstantMatrices();

    // Store data into DYNMAT file format - vasp like
    writeDYNMAT();
  }

  void DirectMethodPC::calculateForceFields() {
    bool LDEBUG=(FALSE || _DEBUG_APL_DIRPHONCALC_ || XHOST.DEBUG);
    string soliloquy="apl::DirectMethodPC::runVASPCalculations():"; //CO190218
    // Extract all forces ////////////////////////////////////////////////////

    //first pass, just find if outfile is found ANYWHERE
    if(!outfileFoundAnywherePhonons(xInputs)){throw APLStageBreak();}

    //second pass, make sure it's everywhere!
    outfileFoundEverywherePhonons(xInputs, _aflowFlags.Directory, *messageFile, _isPolarMaterial);

    // Remove zero state forces if necessary
    if (_calculateZeroStateForces) {
      subtractZeroStateForces(xInputs, _isPolarMaterial);
    }

    // Store forces //////////////////////////////////////////////////////////

    bool generate_plus_minus = false;  // ME190129

    int idxRun = 0;
    uint natoms = DISTORTION_INEQUIVONLY ? _supercell.getNumberOfUniqueAtoms() : _supercell.getNumberOfAtoms();
    for (uint i = 0; i < natoms; i++) {  //CO200212 - int->uint
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
  }

  // ///////////////////////////////////////////////////////////////////////////

  void DirectMethodPC::completeForceFields() {
    stringstream message;
    string function = "apl::DirectMethodPC::completeForceFields()";
    //CO - START
    // Test of stupidity...
    if (_supercell.getEPS() == AUROSTD_NAN) {
      message << "Need to define symmetry tolerance.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ERROR_);
    }
    //CO - END
    // Show info
    message << "Calculating the missing force fields by symmetry.";
    pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);

    // Let's go
    for (int i = 0; i < (DISTORTION_INEQUIVONLY ? _supercell.getNumberOfUniqueAtoms() : _supercell.getNumberOfAtoms()); i++) { //CO190218
      // We need to have 3 linearly independent distortions
      if (_uniqueDistortions[i].size() != 3) {
        vector<xvector<double> > allDistortionsOfAtom;
        vector<xvector<double> > testForce;
        vector<vector<xvector<double> > > forcePool;
        xvector<double> testVec(3), testVec0(3);

        int atomID = (DISTORTION_INEQUIVONLY ? _supercell.getUniqueAtomID(i) : i); //CO190218
        const vector<_sym_op>& agroup = _supercell.getAGROUP(atomID); //CO190218

        // Generate next independent distortion by symmetry operations...
        uint currentSizeDistortions = _uniqueDistortions[i].size(); //CO190218
        for (uint idistor = 0; idistor < currentSizeDistortions; idistor++) { //CO190218
          // Apply all symmetry operations and check if it is independent
          for (uint symOpID = 0; symOpID < agroup.size(); symOpID++) {
            const _sym_op& symOp = agroup[symOpID];

            testForce.clear();
            for (_AFLOW_APL_REGISTER_ int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
              try {
                // ME191219 - atomGoesTo now uses basis_atoms_map; keep translation option in case
                // the basis has not been calculated for some reason
                _AFLOW_APL_REGISTER_ int l = _supercell.atomComesFrom(symOp, k, atomID, true); //CO190218
                testForce.push_back(symOp.Uc * _uniqueForces[i][idistor][l]);
              } catch (aurostd::xerror& e) {
                message << "Mapping problem ? <-> " << k << ".";
                throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
              }
            }

            // Get distortion vector (it is in the cartesian form) and apply symmetry rotation
            testVec = symOp.Uc * _uniqueDistortions[i][idistor];

            // Orthogonalize new rotated distortion vector on all accepted distortion vectors
            for (uint k = 0; k < allDistortionsOfAtom.size(); k++) {
              for (_AFLOW_APL_REGISTER_ int l = 0; l < _supercell.getNumberOfAtoms(); l++) {
                testForce[l] = testForce[l] - getModeratedVectorProjection(forcePool[k][l], testVec, allDistortionsOfAtom[k]);
              }
              testVec = testVec - getVectorProjection(testVec, allDistortionsOfAtom[k]);
            }

            // If the remaining vector is non-zero length, it is new independent direction, hence store it...
            if (aurostd::modulus(testVec) > _AFLOW_APL_EPS_) {
              // Normalize to unit length
              double testVectorLength = aurostd::modulus(testVec);
              for (int l = 0; l < _supercell.getNumberOfAtoms(); l++) {
                testForce[l] = testForce[l] / testVectorLength;
              }
              testVec = testVec / testVectorLength;
              allDistortionsOfAtom.push_back(testVec);

              // We suppose the symOpID == 0 is E (Identity) operation, hence the first
              // independent vector is already the calculated vector, hence new forces are not need to store
              if (allDistortionsOfAtom.size() > 1) {
                // Store new distortion
                _uniqueDistortions[i].push_back(testVec);
                // Store new force field
                _uniqueForces[i].push_back(testForce);
              }

              // Store for next ortogonalization procedure
              forcePool.push_back(testForce);
            }
            if (_uniqueDistortions[i].size() == 3) break;
          }
          if (_uniqueDistortions[i].size() == 3) break;
        }
        allDistortionsOfAtom.clear();
        for (uint ii = 0; ii < forcePool.size(); ii++) forcePool[ii].clear();
        forcePool.clear();
      }

      // I hope this will never happen...
      if (_uniqueDistortions[i].size() != 3) {
        string message = "Cannot complete force fields by symmetry.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void DirectMethodPC::projectToCartesianDirections() {
    bool LDEBUG=(FALSE || _DEBUG_APL_DIRPHONCALC_ || XHOST.DEBUG);
    string soliloquy="apl::DirectMethodPC::projectToCartesianDirections():"; //CO190218
    for (int i = 0; i < (DISTORTION_INEQUIVONLY ? _supercell.getNumberOfUniqueAtoms() : _supercell.getNumberOfAtoms()); i++) { //CO190218
      if(LDEBUG) {cerr << soliloquy << " looking at distorted atom[idistortion=" << i << "]" << std::endl;} //CO190218
      // Construct transformation matrix A
      xmatrix<double> A(3, 3), U(3, 3);
      for (uint j = 0; j < 3; j++) {
        // Ensure it is unit length
        _uniqueDistortions[i][j] = _uniqueDistortions[i][j] / aurostd::modulus(_uniqueDistortions[i][j]);
        if(LDEBUG){ //CO190218
          cerr << soliloquy << " checking if uniqueDistortion[" << i << "][" << j << "] is unit length: ";
          cerr << "modulus(" << _uniqueDistortions[i][j] << ")=" << aurostd::modulus(_uniqueDistortions[i][j]) << std::endl;
        }

        // Copy to rows of U^T matrix
        for (uint k = 1; k <= 3; k++) {
          U(j + 1, k) = _uniqueDistortions[i][j](k);
        }
      }
      A = inverse(U);
      //CO190116 - I believe U is an orthonormal matrix, as it defines a 3D axis
      //hence A = trasp(U) as well (faster)
      //keep for now

      if(LDEBUG){ //CO190218
        cerr << soliloquy << " distortion matrix U(distortion,direction):" << std::endl;
        cerr << U << std::endl;
        cerr << soliloquy << " inverse matrix A:" << std::endl;
        cerr << A << std::endl;
      }

      // Update unique distortion vectors
      // CO190116 - using trasp(A) instead of A because _uniqueDistortions[i][0] is a vector, not a matrix (as m is below)
      // we are really applying A * U == I,
      // so use A below (not trasp(A))
      _uniqueDistortions[i][0] = trasp(A) * _uniqueDistortions[i][0];
      _uniqueDistortions[i][1] = trasp(A) * _uniqueDistortions[i][1];
      _uniqueDistortions[i][2] = trasp(A) * _uniqueDistortions[i][2];

      if(LDEBUG){ //CO190218
        cerr << soliloquy << " new cartesian-direction-projected uniqueDistortions[" << i << "][0]=" << _uniqueDistortions[i][0] << std::endl;
        cerr << soliloquy << " new cartesian-direction-projected uniqueDistortions[" << i << "][1]=" << _uniqueDistortions[i][1] << std::endl;
        cerr << soliloquy << " new cartesian-direction-projected uniqueDistortions[" << i << "][2]=" << _uniqueDistortions[i][2] << std::endl;
        //CO190116 - cerr << soliloquy << " testing: trasp(A) * U should give same as above: trasp(A) * U = " << std::endl;  //U ~ m below
        //CO190116 - cerr << trasp(A) * U << std::endl;
        cerr << soliloquy << " testing: A * U should give same as above: A * U = " << std::endl;  //U ~ m below //DUH A = inverse(U), so A*U = I
        cerr << A * U << std::endl;
      }

      // Update forces
      xmatrix<double> m(3, 3);
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
        if(LDEBUG) {cerr << soliloquy << " looking at supercell atom[" << j << "]" << std::endl;} //CO190218
        for (_AFLOW_APL_REGISTER_ int k = 0; k < 3; k++)
          for (_AFLOW_APL_REGISTER_ int l = 1; l <= 3; l++)
            m(k + 1, l) = _uniqueForces[i][k][j](l);
        if(LDEBUG){ //CO190218
          cerr << soliloquy << " BEFORE m = " << std::endl;
          cerr << m << std::endl;
        }
        // m = A * m * U; ??? I am not sure...
        m = A * m;
        // m = trasp(A) * m;  //CO NEW, treat forces exactly as distortion //CO190116 - wrong, see above, trasp(A) is only for vectors
        if(LDEBUG){ //CO190218
          cerr << soliloquy << " AFTER m = " << std::endl;
          cerr << m << std::endl;
        }
        for (_AFLOW_APL_REGISTER_ int k = 0; k < 3; k++)
          for (_AFLOW_APL_REGISTER_ int l = 1; l <= 3; l++)
            _uniqueForces[i][k][j](l) = m(k + 1, l);
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  void DirectMethodPC::buildForceConstantMatrices() {
    bool LDEBUG=(FALSE || _DEBUG_APL_DIRPHONCALC_ || XHOST.DEBUG);
    string soliloquy="apl::DirectMethodPC::buildForceConstantMatrices():"; //CO190218
    stringstream message;
    // Test of stupidity...
    if (DISTORTION_INEQUIVONLY && !_supercell.getSupercellStructure().fgroup_calculated) { //CO190218
      string message = "The factor group has not been calculated yet.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, message, _RUNTIME_INIT_);
    }
    //CO - START
    if (DISTORTION_INEQUIVONLY && _supercell.getEPS() == AUROSTD_NAN) { //CO190218
      string message = "Need to define symmetry tolerance.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, message, _VALUE_ERROR_);
    }
    //CO - END

    // Clear old matrices
    for (uint i = 0; i < _forceConstantMatrices.size(); i++)
      _forceConstantMatrices[i].clear();
    _forceConstantMatrices.clear();

    //CO190116 - BIG BUG
    //do NOT push_back() with forceConstantMatrices
    //Jahnatek assumed that iatoms structure was [0 1 2 3] [4 5 6] (in order)
    //therefore, pushing back meant keeping forceConstantMatrices in order of supercell atoms
    //this is not necessarily true, as the mappings could be out of order
    //therefore, we create the vector of the necessary dimensions, and put the row in the right place
    //CO190131 UPDATE - this is NOT the only part of the code for which this dependency (iatoms sorted) exists

    for (_AFLOW_APL_REGISTER_ int i = 0; i < _supercell.getNumberOfAtoms(); i++) {
      _forceConstantMatrices.push_back(vector<xmatrix<double> >(0));
      for (_AFLOW_APL_REGISTER_ int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
        _forceConstantMatrices.back().push_back(xmatrix<double>(3,3,1,1));
      }
    }


    //
    message << "Calculating the force constant matrices.";
    pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);

    // We have a party. Let's fun with us...
    //vector<xmatrix<double> > row; //JAHNATEK ORIGINAL //CO190218
    for (int i = 0; i < (DISTORTION_INEQUIVONLY ? _supercell.getNumberOfUniqueAtoms() : _supercell.getNumberOfAtoms()); i++) { //CO190218
      // Get the number of this atom in the whole list
      int basedAtomID = (DISTORTION_INEQUIVONLY ? _supercell.getUniqueAtomID(i) : i); //CO190218

      // This is easy. We know everything. Just construct a set of matrices.
      xmatrix<double> m(3, 3, 1, 1);
      for (int j = 0; j < _supercell.getNumberOfAtoms(); j++) {
        for (uint k = 0; k < _uniqueDistortions[i].size(); k++) {
          double distortionLength = aurostd::modulus(DISTORTION_MAGNITUDE * _uniqueDistortions[i][k]);
          // FCM element = -F/d, but we will omit minus, because next force transformations are better
          // done without it, and in construction of dyn. matrix we will add it to the sum
          // ME200212 - we store _forceConstantMatrices in a human-readable file, so they should
          // represent the actual force constants, not an AFLOW-customized construct
          m(k + 1, 1) = _uniqueForces[i][k][j](1) / distortionLength;
          m(k + 1, 2) = _uniqueForces[i][k][j](2) / distortionLength;
          m(k + 1, 3) = _uniqueForces[i][k][j](3) / distortionLength;
          //cout << i << " " << k << " " << j << " "; printXVector(_uniqueForces[i][k][j]);
        }
        //printXMatrix(m);
        //row.push_back(m); //JAHNATEK ORIGINAL //CO190218
        _forceConstantMatrices[basedAtomID][j] = m;  //CO NEW //CO190218
        if(LDEBUG){ //CO190218
          cerr << soliloquy << " adding m to forceConstantMatrices[" << basedAtomID << "][" << j << "]=" << std::endl;
          cerr << m << std::endl;
        }
      }
      //_forceConstantMatrices.push_back(row);  //JAHNATEK ORIGINAL //CO190218
      //row.clear();  //JAHNATEK ORIGINAL //CO190218

      if(DISTORTION_INEQUIVONLY){ //CO190218
        _sym_op symOp;  //CO
        // Calculate rows for next equivalent atoms starting 1 (structure of iatoms)... //CO190218
        for (int j = 1; j < _supercell.getNumberOfEquivalentAtomsOfType(i); j++) { //CO190218
          try {
            //CO190116 - we want to map the forces of the inequivalent atoms (for which we ran vasp) onto the equivalent ones
            //hence, we need the FGROUP that takes us from the inequivalent atom to the equivalent
            //then, we need to find the atom which, upon application of that symop, becomes k (below)
            //symOp = _supercell.getSymOpWhichMatchAtoms(basedAtomID, _supercell.getUniqueAtomID(i, j), _FGROUP_);  //CO NEW //CO190218
            symOp = _supercell.getSymOpWhichMatchAtoms(_supercell.getUniqueAtomID(i, j), basedAtomID, _FGROUP_);  //JAHNATEK ORIGINAL //CO190218
            //const _sym_op& symOp = _supercell.getSymOpWhichMatchAtoms(_supercell.getUniqueAtomID(i,j),basedAtomID,_FGROUP_); //JAHNATEK ORIGINAL
            //cout << basedAtomID << " -> " << _supercell.getUniqueAtomID(i,j) << " " << symOp.str_type << " shift:"; printXVector(symOp.ftau);
            //printXVector(_supercell.getSupercellStructure().atoms[basedAtomID].fpos);
            //printXVector(_supercell.getSupercellStructure().atoms[_supercell.getUniqueAtomID(i,j)].fpos);
          } catch (aurostd::xerror& e) {
            message << "Mapping problem " << _supercell.getUniqueAtomID(i, j) << " <-> " << basedAtomID << "?"; //CO190218
            throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, message, _RUNTIME_ERROR_);
          }

          for (_AFLOW_APL_REGISTER_ int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
            try {
              //CO190116 - read atomComesFrom() as: applying symOp to l makes k
              //_AFLOW_APL_REGISTER_ int l = _supercell.atomComesFrom(symOp, k, _supercell.getUniqueAtomID(i, j));  //CO NEW //CO190218
              _AFLOW_APL_REGISTER_ int l = _supercell.atomGoesTo(symOp, k, _supercell.getUniqueAtomID(i, j)); //JAHNATEK ORIGINAL //CO190218
              //cout << "MAP " << k << " <-> " << l << std::endl;
              //row.push_back(inverse(symOp.Uc) * _forceConstantMatrices[basedAtomID][l] * symOp.Uc); //JAHNATEK ORIGINAL //CO190218
              //row.push_back(symOp.Uc * _forceConstantMatrices[basedAtomID][l] * inverse(symOp.Uc)); //CO NEW  //JAHNATEK ORIGINAL //CO190218
              //m = symOp.Uc * _forceConstantMatrices[basedAtomID][l] * inverse(symOp.Uc);  //CO NEW //CO190218
              m = inverse(symOp.Uc) * _forceConstantMatrices[basedAtomID][l] * symOp.Uc;  //JAHNATEK ORIGINAL //CO190218
              _forceConstantMatrices[_supercell.getUniqueAtomID(i, j)][k] = m;  //CO NEW //CO190218
              if(LDEBUG){ //CO190218
                cerr << soliloquy << " adding m to forceConstantMatrices[" << _supercell.getUniqueAtomID(i, j) << "][" << k << "]=" << std::endl;
                cerr << m << std::endl;
              }
            } catch (aurostd::xerror& e) {
              message << "Mapping problem " << k << " <-> ?.";
              throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, "Mapping failed.");
            }
          }
          //_forceConstantMatrices.push_back(row);  //JAHNATEK ORIGINAL //CO190218
          //row.clear();  //JAHNATEK ORIGINAL //CO190218
        }
      }
      //row.clear();  //JAHNATEK ORIGINAL //CO190218
    }

    // Test of correctness
    if ((int)_forceConstantMatrices.size() != _supercell.getNumberOfAtoms()) {
      message << "Some problem with the application of factor group operations.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, soliloquy, message, _RUNTIME_ERROR_);
    }

    // ME200211 - force constants are -F/d, not F/d
    for (uint i = 0; i < _forceConstantMatrices.size(); i++) {
      for (uint j = 0; j < _forceConstantMatrices.size(); j++) {
        _forceConstantMatrices[i][j] = -_forceConstantMatrices[i][j];
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  void DirectMethodPC::hibernate(const string& filename) {
    stringstream outfile;

    writeHibernateHeader(outfile);
    writeForceConstants(outfile);
    writeForceField(outfile);
    if (_isPolarMaterial) writePolar(outfile);
    outfile << "</apl>" << std::endl;

    aurostd::stringstream2file(outfile, filename); //ME181226
    if (!aurostd::FileExist(filename)) { //ME181226
      string function = "ForceConstantCalculator::hibernate()";
      string message = "Cannot open output file " + filename + "."; //ME181226
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
    }
  }

  void DirectMethodPC::writeForceField(stringstream& outfile) {
    string tab = " ";
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

  //////////////////////////////////////////////////////////////////////////////

  void DirectMethodPC::writeDYNMAT() {
    string filename = aurostd::CleanFileName(_aflowFlags.Directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_DYNMAT_FILE);  //ME181226
    string message = "Writing forces into file " + filename + ".";
    pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);

    stringstream outfile;

    // 1st line
    outfile << _supercell.getNumberOfUniqueAtoms() << " ";
    outfile << _supercell.getNumberOfAtoms() << " ";
    int dof = 0;
    for (uint i = 0; i < _uniqueDistortions.size(); i++)
      dof += _uniqueDistortions[i].size();
    outfile << dof << std::endl;

    // 2nd line
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(3);
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      if (i != 0) outfile << " ";
      outfile << _supercell.getUniqueAtomMass(i);
    }
    outfile << std::endl;

    // forces + 1 line info about distortion
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        // line info
        outfile << (_supercell.getUniqueAtomID(i) + 1) << " ";
        outfile << (j + 1) << " ";
        xvector<double> shift(3);
        shift = DISTORTION_MAGNITUDE * _uniqueDistortions[i][j];
        outfile << setprecision(3);
        outfile << shift(1) << " " << shift(2) << " " << shift(3) << std::endl;
        // forces
        outfile << setprecision(6);
        for (int k = 0; k < _supercell.getNumberOfAtoms(); k++)
          outfile << setw(15) << _uniqueForces[i][j][k](1)
            << setw(15) << _uniqueForces[i][j][k](2)
            << setw(15) << _uniqueForces[i][j][k](3) << std::endl;
      }
    }

    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)) {
      string function = "DirectMethodPC::writeDYNMAT()";
      message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
    }
  }

   //////////////////////////////////////////////////////////////////////////////

  // This is the interface to phonopy code

  void DirectMethodPC::writeFORCES() {
    string function = "apl::DirectMethodPC::writeFORCES()";
    string message = "Writing forces into file FORCES.";
    pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);

    xstructure ix;
    string filename = "SPOSCAR";
    if (!aurostd::FileEmpty(filename)) {
      message = "Reading " + filename;
      pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);
      stringstream SPOSCAR;
      aurostd::efile2stringstream(filename, SPOSCAR);
      SPOSCAR >> ix;
    } else {
      ix = _supercell.getSupercellStructure();
    }

    stringstream outfile;

    // 1st line
    int dof = 0;
    for (uint i = 0; i < _uniqueDistortions.size(); i++)
      dof += _uniqueDistortions[i].size();
    outfile << dof << std::endl;

    // forces + 1 line info about distortion
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        // line info
        outfile << (_supercell.getUniqueAtomID(i) + 1) << " ";
        xvector<double> shift(3);
        shift = C2F(_supercell.getSupercellStructure().lattice, DISTORTION_MAGNITUDE * _uniqueDistortions[i][j]);
        outfile << setprecision(6);
        outfile << shift(1) << " " << shift(2) << " " << shift(3) << std::endl;
        // forces
        outfile << setprecision(6);
        for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
          int l = 0;
          for (; l < _supercell.getNumberOfAtoms(); l++)
            if ((aurostd::abs(ix.atoms[k].cpos(1) - _supercell.getSupercellStructure().atoms[l].cpos(1)) < _AFLOW_APL_EPS_) &&
                (aurostd::abs(ix.atoms[k].cpos(2) - _supercell.getSupercellStructure().atoms[l].cpos(2)) < _AFLOW_APL_EPS_) &&
                (aurostd::abs(ix.atoms[k].cpos(3) - _supercell.getSupercellStructure().atoms[l].cpos(3)) < _AFLOW_APL_EPS_))
              break;
          //CO, not really mapping error, just mismatch between structure read in (ix) and current supercell structure (should be exact)
          if (l == _supercell.getNumberOfAtoms()) {
            cout << k << std::endl;
            cout << ix.atoms[k].fpos << std::endl;
            cout << ix.atoms[k].cpos << std::endl;
            message = "Mapping error.";
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
          }

          outfile << setw(15) << _uniqueForces[i][j][l](1) << " "
            << setw(15) << _uniqueForces[i][j][l](2) << " "
            << setw(15) << _uniqueForces[i][j][l](3) << std::endl;
        }
      }
    }

    filename = "FORCES";
    aurostd::stringstream2file(outfile, filename);
    if (!aurostd::FileExist(filename)) {
      message = "Cannot open output file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }

  }

  // ///////////////////////////////////////////////////////////////////////////

  void DirectMethodPC::writeXCrysDenForces() {
    string function = "apl::DirectMethodPC::writeXCrysDenForces()";
    string message = "Writing forces into file XCrysDenForces.";
    pflow::logger(_AFLOW_FILE_NAME_, _xInput.xvasp.AVASP_arun_mode, message, _aflowFlags, *messageFile, std::cout);
    _supercell.center_original();  //COREY

    stringstream outfile;  //CO
    // forces + 1 line info about distortion
    for (int i = 0; i < _supercell.getNumberOfUniqueAtoms(); i++) {
      for (uint j = 0; j < _uniqueDistortions[i].size(); j++) {
        //string s = "FORCES_A" + stringify(_supercell.getUniqueAtomID(i)) + "D" + stringify(j) + ".xsf"; //CO
        outfile.str("");  //CO
        //ofstream outfile(s.c_str(), ios_base::out); //CO

        outfile << "CRYSTAL" << std::endl;
        outfile << "PRIMVEC 1" << std::endl;
        outfile << _supercell.getSupercellStructure().lattice << std::endl;
        outfile << "CONVEC 1" << std::endl;
        outfile << _supercell.getSupercellStructure().lattice << std::endl;
        outfile << "PRIMCOORD 1" << std::endl;
        outfile << _supercell.getNumberOfAtoms() << " 1" << std::endl;

        xvector<double> shift(3);
        shift = C2F(_supercell.getSupercellStructure().lattice, DISTORTION_MAGNITUDE * _uniqueDistortions[i][j]);

        outfile << setprecision(6);
        for (int k = 0; k < _supercell.getNumberOfAtoms(); k++) {
          outfile << _supercell.getAtomNumber(k) << " ";
          outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
          outfile << setprecision(8);
          xvector<double> r = F2C(_supercell.getSupercellStructure().lattice,
              _supercell.getSupercellStructure().atoms[k].fpos);
          outfile << setw(15) << r(1) << setw(15) << r(2) << setw(15) << r(3) << " ";
          // this is strange...
          //outfile << setw(15) << _superCellStructure.atoms[k].cpos << " ";

          // Scale force, it is expected in Hartree/Angs.
          xvector<double> f = 27.212 * _uniqueForces[i][j][k];

          outfile << setw(15) << f(1)
            << setw(15) << f(2)
            << setw(15) << f(3) << std::endl;
        }

        string filename = "FORCES_A" + stringify(_supercell.getUniqueAtomID(i)) + "D" + stringify(j) + ".xsf";
        aurostd::stringstream2file(outfile, filename);
        if (!aurostd::FileExist(filename)) {
          string message = "Cannot create " + filename + " file.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
        }

      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////

}  // namespace apl
