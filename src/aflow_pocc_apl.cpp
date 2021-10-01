//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                  Marco Esters - Duke University 2021                    *
// *                                                                         *
//****************************************************************************
// Written by Marco Esters, 2021.
//
// This file provides a framework to calculate phonon properties for
// disordered materials modeled using the POCC algorithm.

#include "aflow.h"
#include "aflow_pocc.h"
#include "APL/aflow_apl.h"

#define _DEBUG_POCC_APL_ false

// Some parts are written within the C++0x support in GCC, especially the std::thread,
// which is implemented in gcc 4.4 and higher.... For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400  // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will not be included, since they need gcc 4.4 and higher (C++0x support)."
#endif

using std::string;
using std::vector;
using std::deque;

static const string _POCC_APL_ERR_PREFIX_ = "pocc::";
static const string _POCC_APL_MODULE_ = "POCC-APL";
static const double _FREQ_WARNING_THRESHOLD_ = 1e-3;

namespace pocc {

  void POccCalculator::calculatePhononPropertiesAPL(const vector<double>& v_temperatures) {
    string function = _POCC_APL_ERR_PREFIX_ + "POccCalculator::calculatePhononProperties()";
    stringstream message;
    // Make sure that everything is consistent
    uint nruns = m_ARUN_directories.size();
    if (nruns == 0) {
      message << "Number of ARUN directories is zero.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }
    if (l_supercell_sets.size() != nruns) {
      message << "Number of directories and number of supercells are different.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }

    vector<int> vexclude;
    string aruns2skip_backup = "";
    if (m_kflags.KBIN_POCC_EXCLUDE_UNSTABLE) {
      // vflag_control (command line) has priority
      if (XHOST.vflag_control.flag("ARUNS2SKIP")) {
        aruns2skip_backup = XHOST.vflag_control.getattachedscheme("ARUNS2SKIP");
      } else if (!m_kflags.KBIN_POCC_ARUNS2SKIP_STRING.empty()) {
        aruns2skip_backup = m_kflags.KBIN_POCC_ARUNS2SKIP_STRING;
      }
    }

    // Parse aflow.in options
    const vector<aurostd::xoption>& vxopts = m_kflags.KBIN_MODULE_OPTIONS.aplflags;
    aurostd::xoption aplopts;
    for (uint i = 0; i < vxopts.size(); i++) {
      const string& key = vxopts[i].keyword;
      aplopts.push_attached(key, vxopts[i].xscheme);
      // Special case: boolean keyword
      if (key == "DOS_PROJECT") aplopts.flag(key, vxopts[i].option);
    }
    apl::validateParametersDosAPL(aplopts, m_aflags, *p_FileMESSAGE);

    vector<apl::PhononCalculator> vphcalc;
    initializePhononCalculators(vphcalc);

    // Get phonon DOS from each run
    vector<xDOSCAR> vxdos = getPhononDoscars(vphcalc, aplopts, vexclude);
    if (vxdos.size() == 0) {
      if (m_kflags.KBIN_POCC_EXCLUDE_UNSTABLE) {
        message << "No structures left after excluding dynamically unstable representatives.";
        pflow::logger(_AFLOW_FILE_NAME_, _POCC_APL_MODULE_, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        return;
      } else {
        message << "No phonon DOSCARs calculated.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
      }
    }
    vphcalc.clear();

    // Compress new files in subdirectories
    for (uint i = 0; i < nruns; i++) {
      if (m_kflags.KZIP_COMPRESS) KBIN::CompressDirectory(m_ARUN_directories[i], m_kflags);
      aurostd::DirectoryChmod("777", m_ARUN_directories[i]);
    }

    // Store PHPOSCAR for plotting
    stringstream phposcar;
    xstructure xstr_phposcar = xstr_pocc;
    xstr_phposcar.is_vasp4_poscar_format = false;
    xstr_phposcar.is_vasp5_poscar_format = true;
    phposcar << xstr_phposcar;
    string phposcar_file = aurostd::CleanFileName(m_aflags.Directory + "/" + DEFAULT_APL_PHPOSCAR_FILE);
    aurostd::stringstream2file(phposcar, phposcar_file);
    if (!aurostd::FileExist(phposcar_file)) {
      message << "Could not write file " << phposcar_file << ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }

    // Calculate vibrational properties
    double T = 0;
    xDOSCAR xdos_T;
    stringstream ossmain;
    apl::ThermalPropertiesCalculator tpc(*p_FileMESSAGE, *p_oss);
    tpc._directory = m_aflags.Directory;
    tpc.natoms = (uint) aurostd::sum(xstr_pocc.comp_each_type);  // Set to number of atoms in parent structure
    double tpt_start = aurostd::string2utype<double>(aplopts.getattachedscheme("TSTART"));
    double tpt_end = aurostd::string2utype<double>(aplopts.getattachedscheme("TEND"));
    double tpt_step = aurostd::string2utype<double>(aplopts.getattachedscheme("TSTEP"));
    for (uint t = 0; t < v_temperatures.size(); t++) {
      T = v_temperatures[t];
      xdos_T = getAveragePhononDos(T, vxdos);
      uint i = 0;
      for ( ; i < xdos_T.number_energies; i++) {
        if (xdos_T.venergy[i] > -_ZERO_TOL_) break;
      }
      if (i > 0) {
        double rel_idos = xdos_T.viDOS[0][i-1]/xdos_T.viDOS[0].back();
        if (rel_idos > _FREQ_WARNING_THRESHOLD_) {
          double percent = 100 * xdos_T.viDOS[0][i-1]/xdos_T.viDOS[0].back();
          message << "There are imaginary frequencies in the DOS at " << T << " K, covering "
            << std::dec << setprecision(3) << percent <<  "\% of the integrated DOS. These frequencies were"
            << " omitted in the calculation of thermodynamic properties.";
          pflow::logger(_AFLOW_FILE_NAME_, _POCC_APL_MODULE_, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        }
      }

      string tstring = getTemperatureString(T);
      string filename = aurostd::CleanFileName(m_aflags.Directory + "/" + POCC_PHDOSCAR_FILE + "_T" + tstring + "K");
      message << "Writing out " << filename << ".";
      pflow::logger(_AFLOW_FILE_NAME_, _POCC_APL_MODULE_, message, *p_FileMESSAGE, *p_oss);
      stringstream phdoscar;
      phdoscar << xdos_T;
      aurostd::stringstream2file(phdoscar, filename);
      if (!aurostd::FileExist(filename)) {
        message << "Could not write file " << filename << ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
      }

      tpc.clear();
      tpc.calculateThermalProperties(tpt_start, tpt_end, tpt_step);
      ossmain << "[POCC_APL_RESULTS]START_TEMPERATURE=" << tstring << "_K" << endl;
      tpc.addToAPLOut(ossmain);
      ossmain << "[POCC_APL_RESULTS]STOP_TEMPERATURE=" << tstring << "_K" << endl;
      tpc.writePropertiesToFile(m_aflags.Directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_THERMO_FILE);
    }

    // Restore ARUNS2SKIP
    if (vexclude.size() > 0) {
      if (XHOST.vflag_control.flag("ARUNS2SKIP")) {
        XHOST.vflag_control.push_attached("ARUNS2SKIP", aruns2skip_backup);
      } else if (!m_kflags.KBIN_POCC_ARUNS2SKIP_STRING.empty()) {
        m_kflags.KBIN_POCC_ARUNS2SKIP_STRING = aruns2skip_backup;
      } else {
        m_kflags.KBIN_POCC_ARUNS2SKIP_STRING = "";
      }
      loadDataIntoCalculator();
      setDFTEnergies();
    }
  }

  void POccCalculator::initializePhononCalculators(vector<apl::PhononCalculator>& vphcalc) {
    string function = "pocc::POccCalculator::initializePhononCalculators()";
    string message = "Initializing phonon calculators.";
    pflow::logger(_AFLOW_FILE_NAME_, _POCC_APL_MODULE_, message, m_aflags, *p_FileMESSAGE, *p_oss);
    unsigned long long int isupercell = 0;
    int imax = 0;
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
      isupercell = std::distance(l_supercell_sets.begin(), it);
      imax = std::max((int) (*it).getHNFIndex(), imax);
      string directory = aurostd::CleanFileName(m_aflags.Directory + "/" + m_ARUN_directories[isupercell]);
      message = "Descending into directory " + directory + ".";
      pflow::logger(_AFLOW_FILE_NAME_, _POCC_APL_MODULE_, message, m_aflags, *p_FileMESSAGE, *p_oss);

      string statefile = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_STATE_FILE);
      if (!aurostd::EFileExist(statefile)) {
        message = "Cannot find state file in directory " + directory + ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
      }

      // Initialize phonon calculator
      apl::PhononCalculator phcalc(*p_FileMESSAGE, *p_oss);
      phcalc.initialize_supercell(statefile);
      phcalc.setDirectory(directory);
      phcalc.setNCPUs(m_kflags);
      phcalc._system = m_vflags.AFLOW_SYSTEM.content_string + "." + m_ARUN_directories[isupercell];
      phcalc.setPolarMaterial(aurostd::EFileExist(directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_POLAR_FILE));

      // Get force constants
      string hibfile = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_HARMIFC_FILE);
      bool awakeHarmIFCs = aurostd::EFileExist(hibfile);
      if (awakeHarmIFCs) {
        try {
          message = "Reading force constant from hibernate file.";
          pflow::logger(_AFLOW_FILE_NAME_, _POCC_APL_MODULE_, message, directory, *p_FileMESSAGE, *p_oss);
          phcalc.awake();
        } catch (aurostd::xerror& e) {
          message = e.error_message + " Recalculating force constants.";
          pflow::logger(_AFLOW_FILE_NAME_, _POCC_APL_MODULE_, message, directory, *p_FileMESSAGE, *p_oss);
          awakeHarmIFCs = false;
        }
      }
      if (!awakeHarmIFCs) {
        // Calculate force constants
        _xvasp xvasp;
        xvasp.Directory = directory;
        _aflags _aflowFlags = m_aflags;
        _aflowFlags.Directory = directory;
        _xinput _xInput(xvasp);
        _xflags _xFlags(m_vflags);
        string AflowIn = "";  // Dummy for AIMS, which is not supported by POCC
        apl::ForceConstantCalculator fccalc(phcalc.getSupercell(), *p_FileMESSAGE, *p_oss);

        // Read parameters and calculate IFCs
        fccalc.readFromStateFile(statefile);
        if (fccalc.run()) {
          fccalc.hibernate();
          phcalc.setHarmonicForceConstants(fccalc);
        } else {
          message = "Could not calculate harmonic force constants for directory " + directory + ".";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
        }
      }

      // Convert to primitive cell for projections
      phcalc.getSupercell().projectToPrimitive();
      vphcalc[isupercell] = phcalc;
    }
  }

  vector<xDOSCAR> POccCalculator::getPhononDoscars(vector<apl::PhononCalculator>& vphcalc, xoption& aplopts, vector<int>& vexclude) {
    string function = "pocc::POccCalculator::getPhononDoscars()";
    string message = "Calculating phonon densities of states.";
    pflow::logger(_AFLOW_FILE_NAME_, _POCC_APL_MODULE_, message, m_aflags, *p_FileMESSAGE, *p_oss);

    // q-point mesh
    vector<int> qpt_mesh;
    aurostd::string2tokens(aplopts.getattachedscheme("DOSMESH"), qpt_mesh, " xX");

    // Calculate phonon DOS
    uint nruns = vphcalc.size();
    vector<xDOSCAR> vxdos(nruns);
    vector<apl::DOSCalculator> vphdos(nruns);
    double minfreq = 0.0, maxfreq = 0.0;
    unsigned long long int isupercell = 0;
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
      isupercell = std::distance(l_supercell_sets.begin(), it);

      string directory = aurostd::CleanFileName(m_aflags.Directory + "/" + m_ARUN_directories[isupercell]);
      vphcalc[isupercell].initialize_qmesh(qpt_mesh, true, true);
      if (!aplopts.flag("DOS_PROJECT")) vphcalc[isupercell].getQMesh().makeIrreducible();

      apl::DOSCalculator dosc(vphcalc[isupercell], aplopts);
      if (m_kflags.KBIN_POCC_EXCLUDE_UNSTABLE && dosc.hasImaginaryFrequencies()) {
        vexclude.push_back(isupercell);
      } else {
        if (dosc.getMinFreq() < minfreq) minfreq = dosc.getMinFreq();
        if (dosc.getMaxFreq() > maxfreq) maxfreq = dosc.getMaxFreq();
        vphdos[isupercell] = dosc;
      }
    }

    uint nexclude = vexclude.size();
    if (nexclude > 0) {
      string aruns2skip = "";
      for (uint i = 0; i < nexclude; i++) {
        aruns2skip += m_ARUN_directories[vexclude[i]];
        if (i < nexclude - 1) aruns2skip += ",";
      }
      if (XHOST.vflag_control.flag("ARUNS2SKIP")) {
        aruns2skip = XHOST.vflag_control.getattachedscheme("ARUNS2SKIP") + "," + aruns2skip;
        XHOST.vflag_control.push_attached("ARUNS2SKIP", aruns2skip);
      } else if (!m_kflags.KBIN_POCC_ARUNS2SKIP_STRING.empty()) {
        m_kflags.KBIN_POCC_ARUNS2SKIP_STRING += "," + aruns2skip;
      } else { // Nothing set in aflow.in or command line, so add to kflags
        m_kflags.KBIN_POCC_ARUNS2SKIP_STRING = aruns2skip;
      }
    }
    aplopts.push_attached("MINFREQ", aurostd::utype2string<double>(minfreq));
    aplopts.push_attached("MAXFREQ", aurostd::utype2string<double>(maxfreq));

    if (aplopts.getattachedscheme("DOSMETHOD") == "LT") message = "Calculating phonon DOS using the linear tetrahedron method.";
    else message = "Calculating phonon DOS using the root sampling method.";
    pflow::logger(_AFLOW_FILE_NAME_, _POCC_APL_MODULE_, message, m_aflags, *p_FileMESSAGE, *p_oss);

#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    int ncpus = KBIN::get_NCPUS(m_kflags);
    vector<std::thread*> threads;
    if (ncpus > (int) nruns) ncpus = (int) nruns;
    vector<vector<int> > thread_dist = getThreadDistribution((int) nruns, ncpus);
    for (int i = 0; i < ncpus; i++) {
      int startIndex = thread_dist[i][0];
      int endIndex = thread_dist[i][1];
      threads.push_back(new std::thread(&POccCalculator::calculatePhononDOSThread, this, startIndex, endIndex,
        std::ref(aplopts), std::ref(vphdos), std::ref(vxdos)));
    }

    for (uint i = 0; i < threads.size(); i++) {
      threads[i]->join();
      delete threads[i];
    }
#else
    calculatePhononDOSThread(0, (int) nruns, aplopts, vphdos, vxdos);
#endif

    return vxdos;
  }

  void POccCalculator::calculatePhononDOSThread(int startIndex, int endIndex,
      const aurostd::xoption& aplopts, vector<apl::DOSCalculator>& vphdos, vector<xDOSCAR>& vxdos) {
    string function = "pocc::POccCalculator::getPhononDoscars()";

    // Normalize to number of branches in the parent structure
    double pocc_sum = aurostd::sum(xstr_pocc.comp_each_type);  // Will be needed for projections
    double nbranches = 3 * pocc_sum;

    // DOS options
    int dos_npoints = aurostd::string2utype<int>(aplopts.getattachedscheme("DOSPOINTS"));
    double dos_smear = aurostd::string2utype<double>(aplopts.getattachedscheme("DOSSMEAR"));
    double minfreq = aurostd::string2utype<double>(aplopts.getattachedscheme("MINFREQ"));
    double maxfreq = aurostd::string2utype<double>(aplopts.getattachedscheme("MAXFREQ"));

    for (int i = startIndex; i < endIndex; i++) {
      vphdos[i].calc(dos_npoints, dos_smear, minfreq, maxfreq, false);
      xDOSCAR phdos = vphdos[i].createDOSCAR();

      stringstream dosfile_raw;
      dosfile_raw << phdos;
      aurostd::stringstream2file(dosfile_raw, m_ARUN_directories[i] + "/PHDOSCAR.pocc.raw");

      double norm_factor = nbranches/((double) vphdos[i].getNumberOfBranches());
      for (uint e = 0; e < phdos.number_energies; e++) phdos.viDOS[0][e] *= norm_factor;

      const xstructure& xstr_phcalc = vphdos[i].getInputStructure();
      uint natoms_phcalc = xstr_phcalc.atoms.size();
      uint nproj = phdos.vDOS[0].size();
      bool projected = (nproj > 1);
      for (uint at = 0; at < (projected?(natoms_phcalc + 1):1); at++) {  // +1 for totals
        for (uint p = 0; p < nproj; p++) {
          for (uint e = 0; e < phdos.number_energies; e++) {
            phdos.vDOS[at][p][0][e] *= norm_factor;
          }
        }
      }

      uint natoms_pocc = xstr_pocc.atoms.size();
      phdos.number_atoms = natoms_pocc;

      // If there are projected DOS, add the DOS contributions that belong to
      // each site in the parent structure. Since there have been relaxation
      // calculations and since AFLOW uses a large initial volume, the lattice
      // vectors are considerably different from the original POCC structure,
      // i.e. FPOSMatch will not help. First order approximation: scale the
      // volume of the initial POCC structure and find the site that is the
      // is closest to the atom.
      if (projected) {
        xstructure xstr_scaled = xstr_pocc;
        double V_new = GetVolume(xstr_phcalc);
        V_new *= pocc_sum/((double) natoms_phcalc);
        xstr_scaled.SetVolume(V_new);

        deque<deque<deque<deque<double> > > > vDOS_pocc(natoms_pocc + 1, deque<deque<deque<double> > >(nproj, deque<deque<double> >(1, deque<double>(phdos.number_energies, 0.0))));
        vDOS_pocc[0] = phdos.vDOS[0];  // Totals stay the same

        xvector<double> cpos(3);
        for (uint at = 0; at < natoms_phcalc; at++) {
          cpos = xstr_scaled.f2c * BringInCell(xstr_scaled.c2f * xstr_phcalc.atoms[at].cpos);
          int at_mapped = -1;
          double min_dist = AUROSTD_MAX_DOUBLE, dist = AUROSTD_MAX_DOUBLE;
          for (uint at_pocc = 0; at_pocc < natoms_pocc; at_pocc++) {
            if (xstr_phcalc.atoms[at].cleanname == xstr_scaled.atoms[at_pocc].cleanname) {
              dist = aurostd::modulus(SYM::minimizeDistanceCartesianMethod(xstr_scaled.atoms[at_pocc].cpos, cpos, xstr_scaled.lattice));
              if (dist < min_dist) {
                min_dist = dist;
                at_mapped = (int) at_pocc;
              }
            }
          }
          if (at_mapped < 0) {
            string message = "Could not map atom " + aurostd::utype2string<uint>(at) + " to parent structure.";
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
          } else {
            for (uint p = 0; p < nproj; p++) {
              for (uint e = 0; e < phdos.number_energies; e++) {
                vDOS_pocc[at_mapped + 1][p][0][e] += phdos.vDOS[at + 1][p][0][e];
              }
            }
          }
        }
        phdos.vDOS = vDOS_pocc;
      }
      vxdos[i] = phdos;
    }
  }

  xDOSCAR POccCalculator::getAveragePhononDos(double T, const vector<xDOSCAR>& vxdos) {
    string function = "pocc::POccCalculator::getAveragePhononDos()";
    stringstream message;
    xDOSCAR xdos;

    setPOccStructureProbabilities(T);
    unsigned long long int isupercell = 0;
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); it != l_supercell_sets.end(); ++it) {
      isupercell = std::distance(l_supercell_sets.begin(), it);
      if (isupercell == 0) { // reset xDOSCAR
        xdos = vxdos[isupercell];
        std::fill(xdos.viDOS[0].begin(), xdos.viDOS[0].end(), 0.0);
        for (uint at = 0; at < xdos.vDOS.size(); at++) {
          for (uint p = 0; p < xdos.vDOS[at].size(); p++) {
            std::fill(xdos.vDOS[at][p][0].begin(), xdos.vDOS[at][p][0].end(), 0.0);
          }
        }
        xdos.temperature = T;
        xdos.Efermi = 0.0;
        xdos.spin = 0;
        xdos.spinF = 0.0;
        if (m_vflags.AFLOW_SYSTEM.content_string.empty()) {
          xdos.title = "Unknown system";
        } else {
          xdos.title = m_vflags.AFLOW_SYSTEM.content_string;
        }
      } else { // check that all dimensions are the same
        bool mismatch = false;
        // Check frequencies
        if (!mismatch) {
          if ((xdos.energy_max != vxdos[isupercell].energy_max)
            || (xdos.energy_min != vxdos[isupercell].energy_min)
            || (xdos.number_energies != vxdos[isupercell].venergy.size())) {
            message << "Frequencies do not match for phonon DOS in " << m_ARUN_directories[isupercell] << ".";
            mismatch = true;
          }
        }
        if (!mismatch) {
          if (xdos.vDOS.size() != vxdos[isupercell].vDOS.size()) {
            message << "Phonon DOS in " << m_ARUN_directories[isupercell] << " has different number of atoms.";
            mismatch = true;
          }
        }
        if (!mismatch) {
          for (uint i = 0; i < xdos.vDOS.size() && !mismatch; i++) {
            if (xdos.vDOS[i].size() != vxdos[isupercell].vDOS[i].size()) {
              message << "Phonon DOS in " << m_ARUN_directories[isupercell] << " has different number of projections.";
              mismatch = true;
            }
          }
        }
        if (mismatch) throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_MISMATCH_);
      }
      uint nat = xdos.vDOS.size();
      uint nproj = xdos.vDOS[0].size();
      for (uint e = 0; e < xdos.number_energies; e++) {
        xdos.viDOS[0][e] += (*it).m_probability * vxdos[isupercell].viDOS[0][e];
        for (uint at = 0; at < nat; at++) {
          for (uint p = 0; p < nproj; p++) {
            xdos.vDOS[at][p][0][e] += (*it).m_probability * vxdos[isupercell].vDOS[at][p][0][e];
          }
        }
      }
    }
    return xdos;
  }

}  // namespace pocc

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                  Marco Esters - Duke University 2020                    *
// *                                                                         *
//****************************************************************************
