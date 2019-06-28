//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2017                    *
// *                                                                         *
//****************************************************************************
// Written by Marco Esters (ME), 2018. Based on work by Jose J. Plata (AFLOW
// AAPL, DOI: 10.1038/s41524-017-0046-7). 
//
// Routines to set up AAPL calculations.

#include "aflow_apl.h"

#define _DEBUG_AAPL_SETUP_ false

using std::vector;
using std::string;
using aurostd::xerror;

static const string _AAPL_FORCES_ERR_PREFIX_ = "apl::PhononCalculator::";

namespace apl {

//setAnharmonicOptions////////////////////////////////////////////////////////
// Sets the calculation options for the calculations of the anharmonic IFCs.
void PhononCalculator::setAnharmonicOptions(int iter, double mix, double threshold) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_SETUP_);
  // ME190501 - START
  // replaced old struct with xoption
  aurostd::xoption options;
  options.push_attached("MAX_ITER", aurostd::utype2string<int>(iter));
  options.push_attached("MIXING_COEFFICIENT", aurostd::utype2string<double>(mix));
  options.push_attached("THRESHOLD", aurostd::utype2string<double>(threshold));
  // ME190501 - END
  anharmonic_IFC_options = options;
  if (LDEBUG) {
    string function = _AAPL_FORCES_ERR_PREFIX_ + "setAnharmonicOptions(): ";
    // ME190501 - START
    std::cerr << function << "mixing coefficient = " << anharmonic_IFC_options.getattachedscheme("MIXING_COEFFICIENT") << std::endl;
    std::cerr << function << "max. number of iterations = " << anharmonic_IFC_options.getattachedscheme("MAX_ITER") << std::endl;
    std::cerr << function << "sumrule threshold = " << anharmonic_IFC_options.getattachedscheme("THRESHOLD") << std::endl;
    // ME190501 - END
  }
}

//buildVaspAAPL///////////////////////////////////////////////////////////////
// Creates the folders for the VASP calculations.
bool PhononCalculator::buildVaspAAPL(const ClusterSet& clst) {
  bool stagebreak = false;
  _xInput.xvasp.AVASP_arun_mode = "AAPL";
  _logger << "Managing directories for ";
  if (clst.order == 3) {
    _logger << "3rd";
  } else {
    _logger << "4th";
  }
  _logger << " order IFCs." << apl::endl;

  // Check if supercell is built
  if (!_supercell.isConstructed()) {
    string function = _AAPL_FORCES_ERR_PREFIX_ + "buildVaspAAPL";
    string message = "The supercell has not been initialized yet.";
    throw xerror(function, message, _RUNTIME_INIT_);
  }

  // Determine the number of runs so the run ID in the folder name can be
  // padded with the appropriate number of zeros.
  int nruns = 0;
  for (uint ineq = 0; ineq < clst.ineq_distortions.size(); ineq++) {
    nruns += clst.ineq_distortions[ineq].distortions.size();
  }
  if (clst.order == 4) {
    for (uint ineq = 0; ineq < clst.higher_order_ineq_distortions.size(); ineq++) {
      nruns += clst.higher_order_ineq_distortions[ineq].distortions.size();
    }
  }
  
  vector<_xinput> x;
  xInputsAAPL.push_back(x);
  vector<_xinput>& xinp = xInputsAAPL.back();
  xinp.assign(nruns, _xInput);
  
  std::cout << "allocated" << std::endl;

  int idxRun = 0;
  for (uint ineq = 0; ineq < clst.ineq_distortions.size(); ineq++) {
    const vector<int>& atoms = clst.ineq_distortions[ineq].atoms;  // ME190108 - Declare to make code more legible
    const _ineq_distortions& idist = clst.ineq_distortions[ineq];  // ME190108 - Declare to make code more legible
    for (uint dist = 0; dist < idist.distortions.size(); dist++) {
      const vector<int>& distortions = idist.distortions[dist][0];  // ME190108 Declare to make code more legible

      // ME 190109 - add title
      xstructure& xstr = xinp[idxRun].getXStr();
      xstr.title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr.title);
      if (xstr.title.empty()) {
        xstr.buildGenericTitle(true, false);
      }
      xstr.title += " AAPL supercell=" + aurostd::joinWDelimiter(clst.sc_dim, 'x');

      // ME 190113 - make sure that POSCAR has the correct format
      if ((!_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_BIN.find("46") != string::npos)) ||
          (_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_MPI_BIN.find("46") != string::npos))) {
        xstr.is_vasp5_poscar_format = false;
      }

      // Set up runname and generate distorted structure
      xinp[idxRun].xvasp.AVASP_arun_runname = buildRunNameAAPL(distortions, atoms, clst.order, idxRun, nruns);
      applyDistortionsAAPL(xinp[idxRun], clst.distortion_vectors, distortions, atoms);

      // Create aflow.in files if they don't exist. Stagebreak is true as soon
      // as one aflow.in file was created.
      stagebreak = (createAflowInPhonons(xinp[idxRun]) || stagebreak);
      idxRun++;
    }
  }
  if (clst.order == 4) {
    for (uint ineq = 0; ineq < clst.higher_order_ineq_distortions.size(); ineq++) {
      const _ineq_distortions& idist = clst.higher_order_ineq_distortions[ineq];
      const vector<int>& atoms = idist.atoms;
      for (uint dist = 0; dist < idist.distortions.size(); dist++) {
        xinp[idxRun] = _xInput;
        const vector<int>& distortions = idist.distortions[dist][0];
        xstructure& xstr = xinp[idxRun].getXStr();
        xstr.title = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(xstr.title);

        if (xstr.title.empty()) {
          xstr.buildGenericTitle(true, false);
        }
        xstr.title += " AAPL supercell=" + aurostd::joinWDelimiter(clst.sc_dim, 'x');
  
        // ME 190113 - make sure that POSCAR has the correct format
        if ((!_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_BIN.find("46") != string::npos)) ||
            (_kbinFlags.KBIN_MPI && (_kbinFlags.KBIN_MPI_BIN.find("46") != string::npos))) {
          xstr.is_vasp5_poscar_format = false;
        }
  
        // Set up runname and generate distorted structure
        xinp[idxRun].xvasp.AVASP_arun_runname = buildRunNameAAPL(distortions, atoms, clst.order, idxRun, nruns);
        xinp[idxRun].xvasp.AVASP_arun_runname += "_3rd";
        applyDistortionsAAPL(xinp[idxRun], clst.distortion_vectors, distortions, atoms, 2.0);
  
        // Create aflow.in files if they don't exist. Stagebreak is true as soon
        // as one aflow.in file was created.
        stagebreak = (createAflowInPhonons(xinp[idxRun]) || stagebreak);
        idxRun++;
      }
    }
  }
  return stagebreak;
}

//buildRunNameAAPL////////////////////////////////////////////////////////////
// Creates the name of the folder for the VASP calculation.
string PhononCalculator::buildRunNameAAPL(const vector<int>& distortions,
                                             const vector<int>& atoms, const int& ord,
                                             const int& run, const int& nruns) {
  std::stringstream runname;
  // runname << "C" << ord << "_"; ME 190112 - removed redundancy
  runname << ord << "_"; // ME 190112

  // ME 190109 - Run ID with padding
  runname << std::setfill('0') << std::setw(aurostd::getZeroPadding(nruns)) << run + 1 << "_";  // ME190112

  // Atom and distortion IDs
  for (uint at = 0; at < atoms.size(); at++) {
    // ME 190112 - made more compact
    if (at > 0) runname << "-";
    runname << "A" << atoms[at] << "D" << distortions[at];
    //runname << "_A" << (at+1) << "_" << atoms[at];
    //runname << "_D" << (at+1) << "_" << distortions[at];
  }
  return runname.str();
}

//applyDistortionsAAPL////////////////////////////////////////////////////////
// Applies the inequivalent distortions to the supercell structures.
void PhononCalculator::applyDistortionsAAPL(_xinput& xinp,
                                            const vector<xvector<double> >& distortion_vectors,
                                            const vector<int>& distortions,
                                            const vector<int>& atoms, double scale) {
  xinp.setXStr(_supercell.getSupercellStructureLight());  // Light copy because we don't need symmetry, etc.
  xstructure& xstr = xinp.getXStr();  // ME 190109
  for (uint at = 0; at < atoms.size(); at++) {
    int atsc, dist_index;
    atsc = atoms[at];
    dist_index = distortions[at];
    xvector<double> dist_cart = distortion_vectors[dist_index];
    while (((at + 1) < atoms.size()) && (atoms[at] == atoms[at+1])) {
       at++;
       dist_index = distortions[at];
       dist_cart += distortion_vectors[dist_index];
    }
    // Normalize dist_cart coordinates to 1 so that distortions have the same magnitude
    for (int i = 1; i < 4; i++) {
      if (abs(dist_cart(i)) < _ZERO_TOL_) {
        dist_cart(i) = 0.0;
      } else {
        dist_cart(i) *= scale/std::abs(dist_cart(i));
      }
    }
    dist_cart *= DISTORTION_MAGNITUDE;
    // ME 190109 - Add to title
    xstr.title += " atom=" + stringify(atoms[at]);
    //xstr.title += " distortion=[" + aurostd::RemoveWhiteSpacesFromTheFrontAndBack(stringify(dist_cart)) + "]"; OBSOLETE ME190112
    std::stringstream distortion; // ME190112 - need stringstream for nicer formatting
    distortion << " distortion=["
               << std::setprecision(3) << dist_cart[1] << ","
               << std::setprecision(3) << dist_cart[2] << ","
               << std::setprecision(3) << dist_cart[3] << "]"; // ME190112
    xstr.title += distortion.str();
    xinp.getXStr().atoms[atsc].cpos += dist_cart;
    xinp.getXStr().atoms[atsc].fpos = C2F(xinp.getXStr().lattice, xinp.getXStr().atoms[atsc].cpos);
  }
}

//calculateAnharmonicIFCs/////////////////////////////////////////////////////
// Calculates the anharmonic IFCs.
void PhononCalculator::calculateAnharmonicIFCs(ClusterSet& clst) {
  _logger << "Checking file integrity for anharmonic IFCs." << apl::endl;
  int o = clst.order - 3;
  if (!outfileFoundAnywherePhonons(xInputsAAPL[o])) {
    throw APLStageBreak();
  }
  outfileFoundEverywherePhonons(xInputsAAPL[o]);
  if (_calculateZeroStateForces) {
    // xInputs is empty if APL awakes, so generate the ZEROSTATE input
    if (xInputs.size() == 0) {
      vector<string> directory;
      aurostd::DirectoryLS(".", directory);
      for (uint d = 0; d < directory.size(); d++) {
        if (aurostd::IsDirectory("./" + directory[d]) &&
            aurostd::substring2bool(directory[d], "ZEROSTATE")) {
          xInputs.push_back(_xInput);
          xInputs.back().setDirectory("./" + directory[d]);
        }
      }
      if (xInputs.size() == 0) {
        string function = _AAPL_FORCES_ERR_PREFIX_ + "calculateAnharmonicIFCs()";
        string message = "Could not find ZEROSTATE directory.";
        throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
      }
    }
    subtractZeroStateForcesAAPL(xInputsAAPL[o], xInputs.back());
  }
  AnharmonicIFCs ifcs(xInputsAAPL[o], clst, DISTORTION_MAGNITUDE,
                      anharmonic_IFC_options, _logger);
  _anharmonicIFCs.push_back(ifcs);
  // Clear runs from memory because they are no longer needed
  xInputsAAPL[o].clear();
}


void PhononCalculator::readAnharmonicIFCs(const string& filename, ClusterSet& clst) {
  _logger << "Reading anharmonic IFCs from file " << filename << "." << apl::endl;
  int o = clst.order - 3;
  AnharmonicIFCs ifcs(filename, clst, DISTORTION_MAGNITUDE,
                      anharmonic_IFC_options, _logger);
  _anharmonicIFCs.push_back(ifcs);
  // Clear runs from memory because they are no longer needed
  xInputsAAPL[o].clear();
}

// ME190114
// Cannot use const reference for zerostate because of getXStr()
void PhononCalculator::subtractZeroStateForcesAAPL(vector<_xinput>& xinps, _xinput& zerostate) {
  if (!zerostate.getXStr().qm_calculated) {
    if (_kbinFlags.AFLOW_MODE_VASP) {
      string vasprunxml_file = zerostate.getDirectory() + string("/vasprun.xml.static");
      if (!aurostd::EFileExist(vasprunxml_file)) {
        vasprunxml_file = zerostate.getDirectory() + string("/vasprun.xml");
        if (!aurostd::EFileExist(vasprunxml_file)) {
          _logger << apl::warning << "The vasprun.xml file in " << zerostate.getDirectory() << " directory is missing." << apl::endl;
          throw APLRuntimeError("apl::PhononCalculator::subtractZeroStateForcesAAL(); Missing data from the ZEROSTATE calculation.");
        }
      }
      xVASPRUNXML vasprunxml;
      vasprunxml.GetForcesFile(vasprunxml_file);
      zerostate.getXStr().qm_forces = vasprunxml.vforces;
    } else if (_kbinFlags.AFLOW_MODE_AIMS) {
      if (!aurostd::EFileExist(zerostate.getDirectory() + string("/aims.out"))) {
        _logger << apl::warning << "The aims.out file in " << zerostate.getDirectory() << " directory is missing." << apl::endl;
        throw APLRuntimeError("apl::PhononCalculator::subtractZeroStateForcesAAPL; Missing data from one job.");
      }
      xAIMSOUT xaimsout(zerostate.getDirectory() + "/aims.out");
      zerostate.getXStr().qm_forces = xaimsout.vforces;
    }
    if ((int) zerostate.getXStr().qm_forces.size() != _supercell.getNumberOfAtoms()) {
      _logger << apl::warning << "The output file in " << zerostate.getDirectory() << " is wrong." << apl::endl;
      throw APLRuntimeError("apl::PhononCalculator::subtractZeroStateForcesAAPL(); Missing data from one job.");
    }
  }
  for (uint idxRun = 0; idxRun < xinps.size(); idxRun++) {
    for (int at = 0; at < _supercell.getNumberOfAtoms(); at++) {
	xinps[idxRun].getXStr().qm_forces[at](1) = xinps[idxRun].getXStr().qm_forces[at](1) - zerostate.getXStr().qm_forces[at](1);
	xinps[idxRun].getXStr().qm_forces[at](2) = xinps[idxRun].getXStr().qm_forces[at](2) - zerostate.getXStr().qm_forces[at](2);
	xinps[idxRun].getXStr().qm_forces[at](3) = xinps[idxRun].getXStr().qm_forces[at](3) - zerostate.getXStr().qm_forces[at](3);
    }
  }
}

} // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                Aflow Marco Esters - Duke University 2018                *
// *                                                                         *
// ***************************************************************************
