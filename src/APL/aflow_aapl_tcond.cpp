//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2018                    *
// *                                                                         *
//****************************************************************************
// Written by Marco Esters, 2018. Based on work by Jose J. Plata (AFLOW AAPL,
// DOI: 10.1038/s41524-017-0046-7) and Jesus Carrete (ShengBTE, 
// DOI: 10.1016/j.cpc.2014.02.015).
//
// This class calculates the thermal conductivity of a material using the
// Boltzmann Transport Equation (BTE). To determine the thermal conductivity,
// this class does the following steps:
//
// 1. It generates a q-point mesh according to the user input.
// 2. The frequencies, group velocities, and eigenvectors are calculated along
//    that mesh.
// 3. It calculates the intrinsic scattering rates, i.e. the parts that do not
//    depend on the temperature, of the anharmonic contributions.
// 4. It solves the BTE and calculates the thermal conductivity.
//
// See aflow_apl.h for descriptions of the classes and their members, and for
// the structs _kcell and _qpoint.

#include "aflow_apl.h"

// Some parts are written within the C++0x support in GCC, especially std::thread,
// which is implemented in gcc 4.4 and higher. For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400
#define AFLOW_APL_MULTITHREADS_ENABLE
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

#define _DEBUG_AAPL_TCOND_ false

// String constants for file output and exception handling
static const string _AAPL_TCOND_ERR_PREFIX_ = "apl::TCONDCalculator::";

static const int max_iter = 250;  // Maximum number of iterations for the iterative BTE solution - ME190408 updated

// Define constants and conversion factors. See AUROSTD/aurostd_xscalar.h for more.
static const double au2THz = 9.648553873170e+02;  // eV/(A amu) -> nm * THz^2
static const double hbar = PLANCKSCONSTANTEV_hbar;  // hbar in eVs
static const double hbar_J = E_ELECTRON * 1e12 * hbar;  // hbar in J/THz;
static const double BEfactor = hbar*1e12/KBOLTZEV;  // hbar/kB in K/THz
static const aurostd::xcomplex<double> iONE(0.0, 1.0);  // imaginary number

using aurostd::xcombos;
using aurostd::xerror;
using std::vector;

/************************************ MPI ***********************************/

#ifdef AFLOW_APL_MULTITHREADS_ENABLE

namespace apl {

//setupMPI////////////////////////////////////////////////////////////////////
// Sets up an MPI calculation.
vector<vector<int> > setupMPI(string message, Logger& log,
                              int nproc, int& ncpus) {
  if (ncpus < 1) ncpus = 1;

  if (ncpus > 1) {
    message += " (" + stringify(ncpus) + " threads)";
  }
  log.initProgressBar(message);

  return getThreadDistribution(nproc, ncpus);
}

//finishMPI///////////////////////////////////////////////////////////////////
// Finishes the MPI progress bar and deletes the threads.
void finishMPI(vector<std::thread*>& threads, Logger& log) {
  for (uint t = 0; t < threads.size(); t++) {
    threads[t]->join();
    delete threads[t];
  }
  log.finishProgressBar();
}

}  // namespace apl
#endif

/************************** CONSTRUCTOR/DESTRUCTOR **************************/

namespace apl {

//Constructor/////////////////////////////////////////////////////////////////
TCONDCalculator::TCONDCalculator(PhononCalculator& pc, QMesh& qm, 
                                 Logger& l) : _pc(pc), _qm(qm), _logger(l) {
  free();
  nBranches = _pc.getNumberOfBranches();
  nQPs = _qm.getnQPs();
  nIQPs = _qm.getnIQPs();
  pcell = _pc.getInputCellStructure();
}

TCONDCalculator::~TCONDCalculator() {
  free();
}

void TCONDCalculator::clear() {
  free();
}

void TCONDCalculator::free() {
  calc_options.clear();
  eigenvectors.clear();
  freq.clear();
  gvel.clear();
  intr_trans_probs.clear();
  irred_qpts_symops.clear();
  nBranches = 0;
  nIQPs = 0;
  nQPs = 0;
  processes.clear();
  temperatures.clear();
  thermal_conductivity.clear();
}

/*********************************** SETUP ***********************************/

//setCalculationOptions///////////////////////////////////////////////////////
// Sets all options for the thermal conductivity calculations, which includes
// the corrections to the scattering rates and the temperature steps.
void TCONDCalculator::setCalculationOptions(string USER_BTE, bool isotope,
                                            bool cumulative, bool fourth_order,
                                            bool boundary, double grain_size,
                                            double temp_start, double temp_end,
                                            double temp_step) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _DEBUG_AAPL_TCOND_);
  string function = _AAPL_TCOND_ERR_PREFIX_ + "setCalculationOptions()";
  if (USER_BTE == "RTA") {
    calc_options.flag("RTA", true);
  } else if (USER_BTE == "FULL") {
    calc_options.flag("RTA", false);
  } else {
    string message = "Illegal value for the flag BTE. Use RTA or FULL.";
    throw xerror(function, message, _INPUT_ILLEGAL_);
  }

  calc_options.flag("ISOTOPE", isotope);
  calc_options.flag("BOUNDARY", boundary);
  calc_options.flag("CUMULATIVE", cumulative);
  calc_options.flag("FOURTH_ORDER", fourth_order);
  calc_options.push_attached("GRAIN_SIZE", aurostd::utype2string<double>(grain_size));
  calc_options.push_attached("TSTART", aurostd::utype2string<double>(temp_start));
  calc_options.push_attached("TEND", aurostd::utype2string<double>(temp_end));
  calc_options.push_attached("TSTEP", aurostd::utype2string<double>(temp_step));

  if (LDEBUG) {
    function += ": ";
    std::cerr << function << "rta_only = " << calc_options.flag("RTA") << std::endl;
    std::cerr << function << "calc_isotopes = " << calc_options.flag("ISOTOPES") << std::endl;
    std::cerr << function << "calc_boundary = " << calc_options.flag("BOUNDARY") << std::endl;
    std::cerr << function << "calc_cumulative = " << calc_options.flag("CUMULATIVE") << std::endl;
    std::cerr << function << "fourth_order = " << calc_options.flag("FOURTH_ORDER") << std::endl;
    std::cerr << function << "grain_size = " << calc_options.getattachedutype<double>("GRAIN_SIZE") << std::endl;
    std::cerr << function << "temp_start = " << calc_options.getattachedutype<double>("TSTART") << std::endl;
    std::cerr << function << "temp_end = " << calc_options.getattachedutype<double>("TEND") << std::endl;
    std::cerr << function << "temp_step = " << calc_options.getattachedutype<double>("TSTEP") << std::endl;
  }
}

}  // namespace apl

/******************************* MAIN FUNCTION ******************************/

namespace apl {

void TCONDCalculator::calculateThermalConductivity() {
  calculateFrequenciesGroupVelocities();
}

} // namespace apl

/*********************** FREQUENCIES/GROUP VELOCITIES ***********************/

namespace apl {

//calculateFrequenciesGroupVelocities/////////////////////////////////////////
// Calculates the frequencies and group velocities for each q-point.
// This function is mostly overhead, the actual calculation happens in
// calculateFreqGvel.
void TCONDCalculator::calculateFrequenciesGroupVelocities() {
  // MPI variables
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  int ncpus, startIndex, endIndex;
  _pc.get_NCPUS(ncpus);
  string messageMPI;
  vector<vector<int> > thread_dist;
  vector<std::thread*> threads;
#endif

  _logger << "Calculating frequencies and group velocities." << apl::endl;
  // Prepare storage
  xmatrix<xcomplex<double> > eigen(nBranches, nBranches, 1, 1);
  eigenvectors.assign(nQPs, eigen);
  freq.assign(nQPs, vector<double>(nBranches));
  xvector<double> g(3);
  gvel.assign(nQPs, vector<xvector<double> >(nBranches, g));

  // Calculate frequencies and group velocities
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  messageMPI = "Frequencies and group velocities";
  thread_dist = setupMPI(messageMPI, _logger, nQPs, ncpus);
  threads.clear();
  for (int icpu = 0; icpu < ncpus; icpu++) {
    startIndex = thread_dist[icpu][0];
    endIndex = thread_dist[icpu][1];
    threads.push_back(new std::thread(&TCONDCalculator::calculateFreqGvel,
                                      this, startIndex, endIndex));
  }
  finishMPI(threads, _logger);
#else
  calculateFreqGvel(0, nQPs);
#endif
  writeFrequencies();
  writeGroupVelocities();
}

//calculateFreqGvel///////////////////////////////////////////////////////////
// Calculates the frequencies and group velocities using the eigenvalue
// solver implemented in apl::PhononCalculator.
void TCONDCalculator::calculateFreqGvel(int startIndex, int endIndex) {
  xmatrix<xcomplex<double> > eigen(nBranches, nBranches, 1, 1);
  vector<xmatrix<xcomplex<double> > > dDynMat(3, eigen);
  for (int q = startIndex; q < endIndex; q++) {
    // Frequency
    xvector<double> f = _pc.getFrequency(_qm.getQPoint(q).cpos, apl::THZ | apl::OMEGA,
                                         eigenvectors[q], dDynMat);
    freq[q] = aurostd::xvector2vector(f);  // Convert to vector to have same indexing as gvel
    // Group velocity
    for (int br = 0; br < nBranches; br++) {
      if (freq[q][br] > _AFLOW_APL_EPS_) {
        xvector<xcomplex<double> > eigenvec = eigenvectors[q].getcol(br+1);
        xvector<xcomplex<double> > eigenvec_conj = conj(eigenvec);
        for (int i = 1; i < 4; i++) {
          xcomplex<double> integral = eigenvec_conj * (dDynMat[i-1] * eigenvec);
          gvel[q][br][i] = au2THz * integral.re/(2.0 * freq[q][br]);
        }
      } else {
        for (int i = 1; i < 4; i++) {
          gvel[q][br][i] = 0.0;
        }
      }
    }
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    _logger.updateProgressBar(1.0/(nQPs - 1));
#endif
  }
}

}  // namespace apl

/******************************** FILE OUTPUT *******************************/

namespace apl {

//writeFrequencies////////////////////////////////////////////////////////////
// Writes the frequencies into a file. Each row belongs to a q-point, and
// each column belongs to a phonon branch.
void TCONDCalculator::writeFrequencies() {
  stringstream output;
  string filename = DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_FREQ_FILE;

  // Header
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  output << "[AAPL_FREQUENCY]START" << std::endl;
  output << std::setiosflags(std::ios::fixed | std::ios::right);
  output << std::setw(10) << "# Q-point";
  output << std::setw(20) << " ";
  output << "Frequencies (THz)" << std::endl;

  // Body
  for (int q = 0; q < nQPs; q++) {
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << std::setw(10) << q;
    for (int br = 0; br < nBranches; br++) {
      output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      output << std::setw(20) << std::setprecision(10) << std::scientific << freq[q][br];
    }
    output << std::endl;
  }

  output << "[AAPL_FREQUENCY]STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;

  // Write to file
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeFrequencies";
    string message = "Could not write frequencies to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

//writeGroupVelocities////////////////////////////////////////////////////////
// Writes the group velocities into a file. Each row belongs to a q-point,
// and each column triplet belongs to a phonon branch.
void TCONDCalculator::writeGroupVelocities() {
  stringstream output;
  string filename = DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_GVEL_FILE;

  // Header
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  output << "[AAPL_GROUP_VELOCITY]START" << std::endl;
  output << std::setiosflags(std::ios::fixed | std::ios::right);
  output << std::setw(10) << "# Q-point";
  output << std::setw(20) << " ";
  output << "Group Velocity (km/s)" << std::endl;

  // Body
  for (int q = 0; q < nQPs; q++) {
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << std::setw(10) << q;
    for (int br = 0; br < nBranches; br++) {
      for (int i = 1; i < 4; i++) {
        output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        output << std::setw(20) << std::setprecision(10) << std::scientific << gvel[q][br][i];
      }
      output << std::setw(5) << " ";
    }
    output << std::endl;
  }

  output << "[AAPL_GROUP_VELOCITY]STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;

  // Write to file
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeGroupVelocities";
    string message = "Could not write group velocities to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

}  // namespace apl
