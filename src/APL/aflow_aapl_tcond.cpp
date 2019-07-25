//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
//****************************************************************************
// Written by Marco Esters, 2019.
//
// This class calculates the thermal conductivity of a material using the
// Boltzmann Transport Equation (BTE).

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
//#define ALAMODE

// String constants for file output and exception handling
static const string _AAPL_TCOND_ERR_PREFIX_ = "apl::TCONDCalculator::";

static const int max_iter = 250;  // Maximum number of iterations for the iterative BTE solution

// Define constants and conversion factors. See AUROSTD/aurostd_xscalar.h for more.
static const double au2THz = 9.648553873170e+02;  // eV/(A amu) -> nm * THz^2
static const double hbar = PLANCKSCONSTANTEV_hbar;  // hbar in eVs
static const double hbar_J = E_ELECTRON * 1e12 * hbar;  // hbar in J/THz;
static const double hbar_amu = hbar * au2THz * 1e13;  // hbar in amu A^2 THz
static const double BEfactor = hbar*1e12/KBOLTZEV;  // hbar/kB in K/THz
static const aurostd::xcomplex<double> iONE(0.0, 1.0);  // imaginary number
static const double scatt_multi[2][3] = {{1.0, 0.5, 0},
                                          {0.5, 0.5, 1.0/6.0}};
static const double ps_prefactors[2] = {2.0/3.0, 6.0/7.0};

using aurostd::xcombos;
using aurostd::xcomplex;
using aurostd::xmatrix;
using aurostd::xvector;
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
  tmpdir = "";
}

}  // namespace apl

/******************************* MAIN FUNCTION ******************************/

namespace apl {

void TCONDCalculator::calculateThermalConductivity() {
  // Setup temporary directory
  //tmpdir = aurostd::TmpDirectoryCreate("TCOND");
  tmpdir = "./test";
  if (aurostd::IsDirectory(tmpdir)) aurostd::RemoveDirectory(tmpdir);
  aurostd::DirectoryMake(tmpdir);
  //
  // Setup temperatures
  double tstart, tend, tstep;
  tstart = aurostd::string2utype<double>(calc_options.getattachedscheme("TSTART"));
  tend = aurostd::string2utype<double>(calc_options.getattachedscheme("TEND"));
  tstep = aurostd::string2utype<double>(calc_options.getattachedscheme("TSTEP"));
  tstart = 300.0;
  tend = 300.0;
  for (double t = tstart; t <= tend; t += tstep) temperatures.push_back(t);

  calculateFrequenciesGroupVelocities();
  writeTempIndepOutput(DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_FREQ_FILE, "Frequency", "THz", freq);
  writeGroupVelocities();
//  vector<vector<vector<xcomplex<double> > > > phases = calculatePhases();
  vector<vector<vector<vector<xcomplex<double> > > > > phases = calculatePhases();
  vector<vector<int> > invariant_sym_ops;
  if (calc_options.flag("RTA")) invariant_sym_ops = calculateInvariantSymOps();
  for (int q = 0; q < nQPs; q++) {
    if (_qm.getIbzqpt(q) == 0) {
      int symop = _qm.getQPoint(q).symop;
      std::cout << q << std::endl;
      std::cout << _qm.getQPoint(q).cpos << std::endl;
      std::cout << "Uf" << std::endl;
      std::cout << _qm.getReciprocalCell().pgroup[symop].Uf << std::endl << std::endl;
      std::cout << "Uc" << std::endl;
      std::cout << _qm.getReciprocalCell().pgroup[symop].Uc << std::endl << std::endl;
      std::cout << eigenvectors[q] << std::endl << std::endl;
    }
  }

  for (int q = 0; q < nQPs; q++) {
    std::cout << q << " " << _qm.getIbzqpt(q) << std::endl;
  }
  calculateTransitionProbabilities(phases);

  xmatrix<double> xmtrx(3, 3);
  thermal_conductivity.assign(temperatures.size(), xmtrx);
  for (uint t = 0; t < temperatures.size(); t++) {
    thermal_conductivity[t] = calculateThermalConductivityTensor(temperatures[t], invariant_sym_ops);
  }
  //aurostd::RemoveDirectory(tmpdir);
}

} // namespace apl

namespace apl {

//vector<vector<vector<xcomplex<double> > > > TCONDCalculator::calculatePhases() {
vector<vector<vector<vector<xcomplex<double> > > > > TCONDCalculator::calculatePhases() {
  const xstructure& scell = _pc.getSuperCellStructure();
  const xstructure& pcell = _pc.getInputCellStructure();
  const vector<int>& sc2pcMap = _pc.getSupercell()._sc2pcMap;
  const vector<int>& pc2scMap = _pc.getSupercell()._pc2scMap;
  uint niatoms = pcell.atoms.size();
  uint natoms = scell.atoms.size();
//  vector<vector<vector<xcomplex<double> > > > phases(niatoms, vector<vector<xcomplex<double> > >(natoms, vector<xcomplex<double> >(nQPs)));
  vector<vector<vector<vector<xcomplex<double> > > > > phases(niatoms, vector<vector<vector<xcomplex<double> > > >(natoms, vector<vector<xcomplex<double> > >(nQPs, vector<xcomplex<double> >(2))));

  int at_eq, at_eq_sc, iat_sc;
  xvector<double> min_vec(3);
  for (uint iat = 0; iat < niatoms; iat++) {
    iat_sc = pc2scMap[iat];
    const xvector<double>& iat_cpos = scell.atoms[iat_sc].cpos;
    for (uint at = 0; at < natoms; at++) {
      min_vec = SYM::minimizeDistanceCartesianMethod(scell.atoms[at].cpos, iat_cpos, scell.lattice);
      at_eq = sc2pcMap[at];
      at_eq_sc = pc2scMap[at_eq];
      min_vec += iat_cpos;
      min_vec -= scell.atoms[at_eq_sc].cpos;
      for (int q = 0; q < nQPs; q++) {
//        phases[iat][at][q] = exp(iONE * scalar_product(_qm.getQPoint(q).cpos, min_vec));
        phases[iat][at][q][0] = exp(iONE * scalar_product(_qm.getQPoint(q).cpos, min_vec));
        phases[iat][at][q][1] = conj(phases[iat][at][q][0]);
      }
    }
  }
  return phases;
}

vector<vector<int> > TCONDCalculator::calculateInvariantSymOps() {
  vector<vector<int> > invar_symops(nIQPs, vector<int>(1, 0));  // Identity is always invariant
  vector<xvector<double> > fpos = _qm.getIrredQPointsFPOS();
  const _kcell& kcell = _qm.getReciprocalCell();
  const vector<_sym_op>& symops = kcell.pgroup;
  uint nsym = symops.size();
  double tol = _AFLOW_APL_EPS_;
  xvector<double> fpos_trans(3);
  for (int iq = 0; iq < nIQPs; iq++) {
    for (uint symop = 1; symop < nsym; symop++) {
      fpos_trans = symops[symop].Uf * fpos[iq];
      if (SYM::FPOSMatch(fpos_trans, fpos[iq], kcell.lattice, kcell.f2c, kcell.skewed, tol)) {
        invar_symops[iq].push_back(symop);
      }
    }
  }
  return invar_symops;
}

} // namespace apl

/*********************** FREQUENCIES/GROUP VELOCITIES ***********************/

namespace apl {

//calculateFrequenciesGroupVelocities/////////////////////////////////////////
// Calculates the frequencies and group velocities for each q-point.
// This function is mostly overhead, the actual calculation happens in
// calculateFreqGvel.
void TCONDCalculator::calculateFrequenciesGroupVelocities() {
  _logger << "Calculating frequencies and group velocities." << apl::endl;
  string message = "Frequencies and group velocities";

  // MPI variables
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  int ncpus;
  _pc.get_NCPUS(ncpus);
  vector<vector<int> > thread_dist;
  vector<std::thread*> threads;
#endif

  // Prepare storage
  xmatrix<xcomplex<double> > eigen(nBranches, nBranches, 1, 1);
  eigenvectors.assign(nQPs, eigen);
  freq.assign(nQPs, vector<double>(nBranches));
  xvector<double> g(3);
  gvel.assign(nQPs, vector<xvector<double> >(nBranches, g));

  // Calculate frequencies and group velocities
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  thread_dist = setupMPI(message, _logger, nQPs, ncpus);
  threads.clear();
  for (int icpu = 0; icpu < ncpus; icpu++) {
    threads.push_back(new std::thread(&TCONDCalculator::calculateFreqGvel, this,
                                      thread_dist[icpu][0], thread_dist[icpu][1]));
  }
  finishMPI(threads, _logger);
#else
  _logger.initProgressBar(message);
  calculateFreqGvel(0, nQPs);
  _logger.finishProgressBar();
#endif
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
    _logger.updateProgressBar(1.0/(nQPs - 1));
  }
}

}  // namespace apl

/*************************** SCATTERING PROCESSES ***************************/

namespace apl {

void TCONDCalculator::calculateScatteringProcesses(int order) {
  // First get all the scattering processes allowed by momemtnum conservation
  _logger << "Determining scattering processes." << apl::endl;
  string message;
  message = "Scattering processes";
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  vector<vector<int> > thread_dist;
  vector<std::thread*> threads;
  int ncpus;
  _pc.get_NCPUS(ncpus);
  thread_dist = setupMPI(message, _logger, nIQPs, ncpus);
  for (int icpu = 0; icpu < ncpus; icpu++) {
    threads.push_back(new std::thread(&TCONDCalculator::getLastQPoint, this,
                      order, thread_dist[icpu][0], thread_dist[icpu][1]));
  }
  finishMPI(threads, _logger);
#else
  _logger.initProgressBar(message);
  getLastQPoint(order, 0, nIQPs);
  _logger.finishProgressBar();
#endif
}

void TCONDCalculator::getLastQPoint(int order, int startIndex, int endIndex) {
  vector<int> qpts(order, 0);
  vector<int> signs(order - 2, 0);
  xvector<double> qsum(3);
  int q;
  string path = tmpdir + "/qpt." + aurostd::utype2string<int>(order) + ".";
  string filename;

  aurostd::xcombos qpt_combos(nQPs, order - 2, 'E', true);
  aurostd::xcombos signs_combos(2, order - 2, 'C', true);
  for (int i = startIndex; i < endIndex; i++) {
    filename = path + aurostd::utype2string<int>(i);
    ofstream out;
    openTmpFile(out, filename);
    qpts[0] = _qm.getIbzqpts()[i];
    signs_combos.reset();
    while (signs_combos.increment()) {
      signs = signs_combos.getCombo();    
      qpt_combos.reset();
      while (qpt_combos.increment()) {
        qsum = _qm.getQPoint(qpts[0]).fpos;
        for (int i = 0; i < order - 2; i++) {
          if (signs[i] == 0) qsum += _qm.getQPoint(qpt_combos.getCombo()[i]).fpos;
          else qsum -= _qm.getQPoint(qpt_combos.getCombo()[i]).fpos;
        }
        q = _qm.getQPointIndex(qsum);
        out.write((char*)(&q), sizeof(int));
      }
    }
    closeTmpFile(out, filename);
    _logger.updateProgressBar(1.0/(nIQPs - 1));
  }
}

} // namespace apl

/********************************** WEIGHTS *********************************/

namespace apl {

double TCONDCalculator::calculateIntegrationWeights(int order, const LTMethod& _lt) {
  _logger << "Calculating integration weights" << apl::endl;
  string message = "Calculating weights";
  double phase_space = 0.0;
  vector<double> ps(nIQPs);
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  vector<std::thread*> threads;
  int ncpus;
  _pc.get_NCPUS(ncpus);
  vector<vector<int> > thread_dist = setupMPI(message, _logger, nIQPs, ncpus);
  for (int icpu = 0; icpu < ncpus; icpu++) {
    threads.push_back(new std::thread(&TCONDCalculator::calculateWeightsLT, this,
                      std::ref(_lt), order, thread_dist[icpu][0], thread_dist[icpu][1], std::ref(ps)));
  }
  finishMPI(threads, _logger);
#else
  _logger.initProgressBar(message);
  calculateWeightsLT(_lt, order, 0, nIQPs, ps);
  _logger.finishProgressBar();
#endif
  for (int i = 0; i < nIQPs; i++) phase_space += ps[i];
  return phase_space;
}

void TCONDCalculator::calculateWeightsLT(const LTMethod& _lt, int order,
                                         int startIndex, int endIndex, vector<double>& ps) {
  vector<vector<int> > branches, qptstest;
  vector<vector<bool> > signs;
  vector<double> frequencies(nQPs), weights(nQPs);
  vector<int> lastq(nQPs);
  xvector<double> qsum(3);
  double freq_ref;
  bool sign;
  vector<bool> zerofreq(nQPs);
  uint b, s;
  int iq, qpt, q, o;
  const int o2 = order - 2;
  const int o3 = order - 3;
  const double prefactor = ps_prefactors[o3]/(std::pow(nBranches, order) * std::pow(nQPs, o2));

  string pathp, pathq, pathw, filenamep, filenameq, filenamew;
  pathp = tmpdir + "/proc." + aurostd::utype2string<int>(order) + ".";
  pathq = tmpdir + "/qpt." + aurostd::utype2string<int>(order) + ".";
  pathw = tmpdir + "/wght." + aurostd::utype2string<int>(order) + ".";

  aurostd::xcombos qpt_combos(nQPs, o3, 'E', true);
  while (qpt_combos.increment()) qptstest.push_back(qpt_combos.getCombo());
  aurostd::xcombos sign_combos(2, o2, 'C', true);
  while (sign_combos.increment()) {
    vector<bool> sgn(o2);
    for (o = 0; o < o2; o++) {
      if (sign_combos.getCombo()[o] == 0) sgn[o] = false;
      else sgn[o] = true;
    }
    signs.push_back(sgn);
  }
  uint nsigns = signs.size();
  aurostd::xcombos branch_combos(nBranches, order, 'E', true);
  while (branch_combos.increment()) branches.push_back(branch_combos.getCombo());
  uint nbr = branches.size();

  int nprocs = nIQPs;
  for (int i = 3; i < order; i++) nprocs *= nQPs;
  for (int i = startIndex; i < endIndex; i++) {
    int c = 0;
    iq = _qm.getIbzqpts()[i];
    qpt_combos.reset();
    bool increment = ((order == 3) || qpt_combos.increment());
    while (increment) {
      const vector<int>& qpoints = qpt_combos.getCombo();
      ofstream outp, outw, outq;
      filenamew = pathw + aurostd::utype2string<int>(i);
      filenamep = pathp + aurostd::utype2string<int>(i);
      filenameq = pathq + aurostd::utype2string<int>(i);
      if (order > 3) {
        string suffix = "." + aurostd::joinWDelimiter(qpoints, "_");
        filenamep += suffix;
        filenameq += suffix;
        filenamew += suffix;
      }
      openTmpFile(outp, filenamep);
      openTmpFile(outq, filenameq);
      openTmpFile(outw, filenamew);
      for (s = 0; s < nsigns; s++) {
        qsum = _qm.getQPoint(iq).fpos;
        for (o = 0; o < o3; o++) {
          if (!signs[s][o]) qsum += _qm.getQPoint(qpoints[o]).fpos;
          else qsum -= _qm.getQPoint(qpoints[o]).fpos;
        }
        for (q = 0; q < nQPs; q++) {
          if (!signs[s].back()) qpt = _qm.getQPointIndex(qsum + _qm.getQPoint(q).fpos);
          else qpt = _qm.getQPointIndex(qsum - _qm.getQPoint(q).fpos);
          lastq[q] = qpt;
          outq.write((char*)(&qpt), sizeof(int));
        }
        for (b = 0; b < nbr; b++) {
          freq_ref = freq[iq][branches[b][0]];
          if (freq_ref > _AFLOW_APL_EPS_) {
            for (o = 0; o < o3; o++) {
              if (freq[qpoints[o]][branches[b][o+1]] < _AFLOW_APL_EPS_) break;
              if (!signs[s][o]) freq_ref += freq[qpoints[o]][branches[b][o + 1]];
              else freq_ref -= freq[qpoints[o]][branches[b][o + 1]];
            }
            if (o == o3) {
              for (q = 0; q < nQPs; q++) {
                zerofreq[q] = ((freq[q][branches[b][o2]] < _AFLOW_APL_EPS_) || (freq[lastq[q]][branches[b].back()] < _AFLOW_APL_EPS_));
                if (!signs[s].back()) frequencies[q] = -freq[q][branches[b][o2]];
                else frequencies[q] = freq[q][branches[b][o2]];
                frequencies[q] += freq[lastq[q]][branches[b].back()];
              }
              weights = getWeightsLT(_lt, freq_ref, frequencies);
              for (q = 0; q < nQPs; q++) {
                if (!zerofreq[q] && (weights[q] > _ZERO_TOL_)) {
                  outw.write((char*)(&weights[q]), sizeof(double));
                  qpt = c + q;
                  for (o = 0; o < o2; o++) {
                    sign = signs[s][o];
                    outp.write((char*)(&sign), sizeof(bool));
                  }
                  outp.write((char*)(&qpt), sizeof(int));
                  outp.write((char*)(&b), sizeof(int));
                  ps[i] += scatt_multi[o3][s] * weights[q];
                }
              }
            }
          }
        }
      }
      closeTmpFile(outp, filenamep);
      closeTmpFile(outq, filenameq);
      closeTmpFile(outw, filenamew);
      increment = ((order > 3) && qpt_combos.increment());
      c += nQPs;
      _logger.updateProgressBar(1.0/nprocs);
    }
    ps[i] *= _qm.getWeights()[i] * prefactor;
  }
}

#ifdef ALAMODE

double fij(double ei, double ej, double e) {
  return (e - ej)/(ei - ej);
}

void insertion_sort(double *a, int *ind, int n) {
  int i;
  for (i = 0; i < n; ++i) ind[i] = i;
  for (i = 1; i < n; ++i) {
    double tmp = a[i];
    int j = i;
    while (j > 0 && tmp < a[j - 1]) {
      a[j] = a[j - 1];
      ind[j] = ind[j - 1];
      --j;
    }
    a[j] = tmp;
    ind[j] = i;
  }
}

vector<double> TCONDCalculator::getWeightsLT(const LTMethod& _lt, double e_ref,
                                             const vector<double>& energy) {
  vector<double> weight(nQPs);
  const vector<vector<int> >& corners = _lt.getTetrahedra();
  double vol = _lt.getVolumePerTetrahedron();

  double g;

  double e_tmp[4];
  int sort_arg[4], kindex[4];

  for (int i = 0; i < _lt.getnTetrahedra(); ++i) {
    for (int j = 0; j < 4; ++j) {
      e_tmp[j] = energy[corners[i][j]];
      kindex[j] = corners[i][j];
    }
    insertion_sort(e_tmp, sort_arg, 4);
    double e1 = e_tmp[0];
    double e2 = e_tmp[1];
    double e3 = e_tmp[2];
    double e4 = e_tmp[3];

    int k1 = kindex[sort_arg[0]];
    int k2 = kindex[sort_arg[1]];
    int k3 = kindex[sort_arg[2]];
    int k4 = kindex[sort_arg[3]];

    double I1 = 0.0;
    double I2 = 0.0;
    double I3 = 0.0;
    double I4 = 0.0;

    if (e3 <= e_ref && e_ref < e4) {
        g = std::pow(e4 - e_ref, 2) / ((e4 - e1) * (e4 - e2) * (e4 - e3));

        I1 = g * fij(e1, e4, e_ref);
        I2 = g * fij(e2, e4, e_ref);
        I3 = g * fij(e3, e4, e_ref);
        I4 = g * (fij(e4, e1, e_ref) + fij(e4, e2, e_ref) + fij(e4, e3, e_ref));

    } else if (e2 <= e_ref && e_ref < e3) {
        g = (e2 - e1 + 2.0 * (e_ref - e2) - (e4 + e3 - e2 - e1)
            * std::pow(e_ref - e2, 2) / ((e3 - e2) * (e4 - e2))) / ((e3 - e1) * (e4 - e1));

        I1 = g * fij(e1, e4, e_ref) + fij(e1, e3, e_ref) * fij(e3, e1, e_ref) * fij(e2, e3, e_ref) / (e4 - e1);
        I2 = g * fij(e2, e3, e_ref) + std::pow(fij(e2, e4, e_ref), 2) * fij(e3, e2, e_ref) / (e4 - e1);
        I3 = g * fij(e3, e2, e_ref) + std::pow(fij(e3, e1, e_ref), 2) * fij(e2, e3, e_ref) / (e4 - e1);
        I4 = g * fij(e4, e1, e_ref) + fij(e4, e2, e_ref) * fij(e2, e4, e_ref) * fij(e3, e2, e_ref) / (e4 - e1);

    } else if (e1 <= e_ref && e_ref < e2) {
        g = std::pow(e_ref - e1, 2) / ((e2 - e1) * (e3 - e1) * (e4 - e1));

        I1 = g * (fij(e1, e2, e_ref) + fij(e1, e3, e_ref) + fij(e1, e4, e_ref));
        I2 = g * fij(e2, e1, e_ref);
        I3 = g * fij(e3, e1, e_ref);
        I4 = g * fij(e4, e1, e_ref);

    }
    weight[k1] += vol * I1;
    weight[k2] += vol * I2;
    weight[k3] += vol * I3;
    weight[k4] += vol * I4;
  }
  return weight;
}

double TCONDCalculator::integrate(const LTMethod& _lt, double e_ref,
                                  const vector<double>& f, const vector<double>& energy) {
  double ret = 0.0;
  double I1, I2, I3, I4;

  double frac3 = 1.0 / 3.0;
  double g;

  const vector<vector<int> >& corners = _lt.getTetrahedra();
  double vol = _lt.getVolumePerTetrahedron();

  vector<std::pair<double, double> > tetra_data(4);

  for (int i = 0; i < _lt.getnTetrahedra(); ++i) {

    for (uint j = 0; j < 4; j++) {
      int knum = corners[i][j];
      tetra_data[j].first = energy[knum];
      tetra_data[j].second = f[knum];
    }

    std::sort(tetra_data.begin(), tetra_data.end());

    double e1 = tetra_data[0].first;
    double e2 = tetra_data[1].first;
    double e3 = tetra_data[2].first;
    double e4 = tetra_data[3].first;

    double f1 = tetra_data[0].second;
    double f2 = tetra_data[1].second;
    double f3 = tetra_data[2].second;
    double f4 = tetra_data[3].second;

    if (e3 <= e_ref && e_ref < e4) {
        g = 3.0 * std::pow(e4 - e_ref, 2) / ((e4 - e1) * (e4 - e2) * (e4 - e3));

        I1 = frac3 * fij(e1, e4, e_ref);
        I2 = frac3 * fij(e2, e4, e_ref);
        I3 = frac3 * fij(e3, e4, e_ref);
        I4 = frac3 * (fij(e4, e1, e_ref) + fij(e4, e2, e_ref) + fij(e4, e3, e_ref));

        ret += vol * g * (I1 * f1 + I2 * f2 + I3 * f3 + I4 * f4);

    } else if (e2 <= e_ref && e_ref < e3) {
        g = 3.0 * (e2 - e1 + 2.0 * (e_ref - e2) - (e4 + e3 - e2 - e1)
            * std::pow(e_ref - e2, 2) / ((e3 - e2) * (e4 - e2))) / ((e3 - e1) * (e4 - e1));

        I1 = frac3 * fij(e1, e4, e_ref) * g + fij(e1, e3, e_ref) * fij(e3, e1, e_ref) * fij(e2, e3, e_ref) / (e4 -
            e1);
        I2 = frac3 * fij(e2, e3, e_ref) * g + std::pow(fij(e2, e4, e_ref), 2) * fij(e3, e2, e_ref) / (e4 - e1);
        I3 = frac3 * fij(e3, e2, e_ref) * g + std::pow(fij(e3, e1, e_ref), 2) * fij(e2, e3, e_ref) / (e4 - e1);
        I4 = frac3 * fij(e4, e1, e_ref) * g + fij(e4, e2, e_ref) * fij(e2, e4, e_ref) * fij(e3, e2, e_ref) / (e4 -
            e1);

        ret += vol * (I1 * f1 + I2 * f2 + I3 * f3 + I4 * f4);

    } else if (e1 <= e_ref && e_ref < e2) {
        g = 3.0 * std::pow(e_ref - e1, 2) / ((e2 - e1) * (e3 - e1) * (e4 - e1));

        I1 = frac3 * (fij(e1, e2, e_ref) + fij(e1, e3, e_ref) + fij(e1, e4, e_ref));
        I2 = frac3 * fij(e2, e1, e_ref);
        I3 = frac3 * fij(e3, e1, e_ref);
        I4 = frac3 * fij(e4, e1, e_ref);

        ret += vol * g * (I1 * f1 + I2 * f2 + I3 * f3 + I4 * f4);

    }
    
  }
  return ret;
}

#else

vector<double> TCONDCalculator::getWeightsLT(const LTMethod& _lt, double freq_ref,
                                             const vector<double>& _freqs) {
  vector<double> weights(nQPs);
  double f[4], w[4];
  int k[4];
  const vector<vector<int> >& corners = _lt.getTetrahedra();
  double vol = _lt.getVolumePerTetrahedron();

  int i, j;
  double t, t2;

  // Shorthands for weight calculations. Will only populate necessary values.
  // For notation and formulas, see DOI 10.1007/978-3-642-25864-0_9.
  double E[4];
  double eps[4][4];

  for (int itet = 0; itet < _lt.getnTetrahedra(); itet++) {
    for (int icorner = 0; icorner < 4; icorner++) {
      f[icorner] = _freqs[corners[itet][icorner]];
      k[icorner] = icorner;
    }
    // Perform an explicit insertion sort to keep track of swaps
    for (i = 1; i < 4; i++) {
      t = f[i];
      j = i;
      while ((j > 0) && (t < f[j - 1])) {
        f[j] = f[j - 1];
        k[j] = k[j - 1];
        j--;
      }
      f[j] = t;
      k[j] = i;
    }
    if ((f[2] <= freq_ref) && (freq_ref < f[3])) {
      E[3] = freq_ref - f[3];
      eps[3][0] = f[3] - f[0];
      eps[3][1] = f[3] - f[1];
      eps[3][2] = f[3] - f[2];
      t = std::pow(E[3], 2)/(eps[3][0] * eps[3][1] * eps[3][2]);
      for (j = 0; j < 3; j++) w[j] = -t * E[3]/eps[3][j];
      w[3] = 3.0;
      for (i = 0; i < 3; i++) w[3] += E[3]/eps[3][i];
      w[3] *= t;
    } else if ((f[1] <= freq_ref) && (freq_ref < f[2])) {
      for (i = 0; i < 4; i++) E[i] = freq_ref - f[i];
      for (i = 1; i < 4; i++) {
        for (j = 0; j < i; j++) {
          eps[i][j] = f[i] - f[j];
        }
      }
      t = E[0] * E[2]/(eps[2][0] * eps[2][1] * eps[3][0]);
      t2 = E[1] * E[3]/(eps[2][1] * eps[3][0] * eps[3][1]);
      w[0] = t * (E[2]/eps[2][0] + E[3]/eps[3][0]) + t2 * E[3]/eps[3][0];
      w[1] = t * E[2]/eps[2][1] + t2 * (E[2]/eps[2][1] + E[3]/eps[3][1]);
      w[2] = -t * (E[0]/eps[2][0] + E[1]/eps[2][1]) - t2 * E[1]/eps[2][1];
      w[3] = -t * E[0]/eps[3][0] - t2 * (E[0]/eps[3][0] + E[1]/eps[3][1]);
    } else if ((f[0] <= freq_ref) && (freq_ref < f[1])) {
      E[0] = freq_ref - f[0];
      eps[1][0] = f[1] - f[0];
      eps[2][0] = f[2] - f[0];
      eps[3][0] = f[3] - f[0];
      t = std::pow(E[0], 2)/(eps[1][0] * eps[2][0] * eps[3][0]);
      w[0] = 3.0;
      for (i = 1; i < 4; i++) w[0] -= E[0]/eps[i][0];
      w[0] *= t;
      for (j = 1; j < 4; j++) w[j] = t * E[0]/eps[j][0];
    } else {
      for (int i = 0; i < 4; i++) w[i] = 0.0;
    }

    for (i = 0; i < 4; i++) weights[corners[itet][k[i]]] += vol * w[i];
  }
  return weights;
}

double TCONDCalculator::integrate(const LTMethod& _lt, double freq_ref,
                                  const vector<double>& _func, const vector<double>& _freqs) {
  double integ = 0.0;
  double c[4], f[4], w[4];
  int k[4];
  const vector<vector<int> >& corners = _lt.getTetrahedra();
  double vol = _lt.getVolumePerTetrahedron();

  int i, j;
  double t, t2;

  // Shorthands for weight calculations. Will only populate necessary values.
  // For notation and formulas, see DOI 10.1007/978-3-642-25864-0_9.
  double E[4];
  double eps[4][4];

  for (int itet = 0; itet < _lt.getnTetrahedra(); itet++) {
    for (int icorner = 0; icorner < 4; icorner++) {
      f[icorner] = _freqs[corners[itet][icorner]];
      c[icorner] = _func[corners[itet][icorner]];
      k[icorner] = icorner;
    }
    // Perform an explicit insertion sort to keep track of swaps
    for (i = 1; i < 4; i++) {
      t = f[i];
      j = i;
      while ((j > 0) && (t < f[j - 1])) {
        f[j] = f[j - 1];
        c[j] = c[j - 1];
        k[j] = k[j - 1];
        j--;
      }
      f[j] = t;
      c[j] = c[j - 1];
      k[j] = i;
    }
    if ((f[0] <= freq_ref) && (freq_ref < f[1])) {
      E[0] = freq_ref - f[0];
      eps[1][0] = f[1] - f[0];
      eps[2][0] = f[2] - f[0];
      eps[3][0] = f[3] - f[0];
      t = std::pow(E[0], 2)/(eps[1][0] * eps[2][0] * eps[3][0]);
      w[0] = 3.0;
      for (i = 1; i < 4; i++) w[0] -= E[0]/eps[i][0];
      w[0] *= t;
      for (j = 1; j < 4; j++) w[j] = t * E[0]/eps[j][0];
    } else if ((f[1] <= freq_ref) && (freq_ref < f[2])) {
      for (i = 0; i < 4; i++) E[i] = freq_ref - f[i];
      for (i = 1; i < 4; i++) {
        for (j = 0; j < i; j++) {
          eps[i][j] = f[i] - f[j];
        }
      }
      t = E[0] * E[2]/(eps[2][0] * eps[2][1] * eps[3][0]);
      t2 = E[1] * E[3]/(eps[2][1] * eps[3][0] * eps[3][1]);
      w[0] = t * (E[2]/eps[2][0] + E[3]/eps[3][0]) + t2 * E[3]/eps[3][0];
      w[1] = t * E[2]/eps[2][1] + t2 * (E[2]/eps[2][1] + E[3]/eps[3][1]);
      w[2] = -t * (E[0]/eps[2][0] + E[1]/eps[2][1]) - t2 * E[1]/eps[2][1];
      w[3] = -t * E[0]/eps[3][0] - t2 * (E[0]/eps[3][0] + E[1]/eps[3][1]);
    } else if ((f[2] <= freq_ref) && (freq_ref < f[3])) {
      E[3] = freq_ref - f[3];
      eps[3][0] = f[3] - f[0];
      eps[3][1] = f[3] - f[1];
      eps[3][2] = f[3] - f[2];
      t = std::pow(E[3], 2)/(eps[3][0] * eps[3][1] * eps[3][2]);
      for (j = 0; j < 3; j++) w[j] = -t * E[3]/eps[3][j];
      w[3] = 3.0;
      for (i = 0; i < 3; i++) w[3] += E[3]/eps[3][i];
      w[3] *= t;
    } else {
      for (i = 0; i < 4; i++) w[i] = 0;
    }
    for (i = 0; i < 4; i++) integ += vol * c[i] * w[i];
  }
  return integ;
}

#endif

}  // namespace apl

/************************* TRANSITION PROBABILITIES *************************/

namespace apl {

//void TCONDCalculator::calculateTransitionProbabilities(const vector<vector<vector<xcomplex<double> > > >& phases) {
void TCONDCalculator::calculateTransitionProbabilities(const vector<vector<vector<vector<xcomplex<double> > > > >& phases) {
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  int ncpus;
  _pc.get_NCPUS(ncpus);
  vector<std::thread*> threads;
  vector<vector<int> > thread_dist;
#endif
  string message;
  LTMethod _lt(_qm, _logger);

  int max_order;
  if (calc_options.flag("FOURTH_ORDER")) max_order = 4;
  else max_order = 3;
  vector<double> phase_space(max_order - 2);
  for (int order = 3; order <= max_order; order++) {
    _logger << "Calculating transition probabilities for " << order << "-phonon processes." << apl::endl;
    phase_space[order - 3] = calculateIntegrationWeights(order, _lt);
    _logger << "phase space = " << phase_space[order - 3] << apl::endl;

    message = "Transition Probabilities";
    _logger << "Calculating transition probabilities." << apl::endl;
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    thread_dist = setupMPI(message, _logger, nIQPs, ncpus);
    threads.clear();
    for (int icpu = 0; icpu < ncpus; icpu++) {
      threads.push_back(new std::thread(&TCONDCalculator::calculateTransitionProbabilitiesPhonon, this,
                                        order, thread_dist[icpu][0], thread_dist[icpu][1], std::ref(phases)));
    }
    finishMPI(threads, _logger);
#else
    _logger.initProgressBar(message);
    calculateTransitionProbabilitiesPhonon(order, 0, nIQPs, phases);
    _logger.finishProgressBar();
#endif
    _logger << "Finished" << apl::endl;
    for (int i = 0; i < 1; i++) {
      std::cout << i << std::endl;
      ifstream transfl, procfl, qptfl;
      double t;
      bool sign;
      int s, q, q2, br;
      vector<double> trans;
      vector<vector<int> > processes;
      vector<int> proc(4);
      string suffix = ".3." + aurostd::utype2string<int>(i);
      openTmpFile(transfl, tmpdir + "/trans" + suffix);
      openTmpFile(procfl, tmpdir + "/proc" + suffix);
      openTmpFile(qptfl, tmpdir + "/qpt" + suffix);
      while (transfl.peek() != EOF) {
        transfl.read((char*)(&t), sizeof(double));
        trans.push_back(t);
        procfl.read((char*)(&sign), sizeof(bool));
        procfl.read((char*)(&q), sizeof(int));
        procfl.read((char*)(&br), sizeof(int));
        s = (int) sign;
        qptfl.seekg(sizeof(int) * (s + q));
        qptfl.read((char*)(&q2), sizeof(int));
        proc[0] = s;
        proc[1] = q;
        proc[2] = q2;
        proc[3] = br;
        processes.push_back(proc);
        std::cout << aurostd::joinWDelimiter(processes.back(), " ") << " " << trans.back() << std::endl;
      }
      for (uint i = 0; i < trans.size(); i++) {
        std::cout << aurostd::joinWDelimiter(processes[i], " ") << " " << trans[i] << std::endl;
        if (trans[i] > _ZERO_TOL_) {
          for (uint j = i + 1; j < trans.size(); j++) {
            if (std::abs(trans[i] - trans[j]) < _ZERO_TOL_) {
              std::cout << aurostd::joinWDelimiter(processes[j], " ") << std::endl;
            }
          }
        }
      }
      closeTmpFile(transfl, tmpdir + "/trans" + suffix);
      closeTmpFile(procfl, tmpdir + "/proc" + suffix);
      closeTmpFile(procfl, tmpdir + "/qpt" + suffix);
    }
  }

  if (0 && calc_options.flag("ISOTOPE")) {
    vector<vector<double> > isotope_rates(nIQPs, vector<double>(nBranches, 0.0));
    _logger << "Calculating isotope transition probabilities." << apl::endl;
    message = "Isotope Transition Probabilities";
#ifdef AFLOW_APL_MULTITHREADS_ENABLE    
    thread_dist = setupMPI(message, _logger, nIQPs, ncpus);
    threads.clear();
    for (int icpu = 0; icpu < ncpus; icpu++) {
      threads.push_back(new std::thread(&TCONDCalculator::calculateTransitionProbabilitiesIsotope, this,
                                        std::ref(_lt), thread_dist[icpu][0], thread_dist[icpu][1], std::ref(isotope_rates)));
    }
    finishMPI(threads, _logger);
#else
    _logger.initProgressBar(message);
    calculateTransitionProbabilitiesIsotope(_lt, 0, nIQPs, isotope_rates);
    _logger.finishProgressBar();
#endif
    writeRatesToTmpFile("rates.iso", isotope_rates);
  }

  if (calc_options.flag("BOUNADRY")) {
    calculateTransitionProbabilitiesBoundary();
  }
}

void TCONDCalculator::calculateTransitionProbabilitiesPhonon(int order, int startIndex, int endIndex,
//                                                             const vector<vector<vector<xcomplex<double> > > >& phases) {
                                                             const vector<vector<vector<vector<xcomplex<double> > > > >& phases) {
  const int o3 = order - 3;
  const int o2 = order - 2;
  const int o1 = order - 1;

  const Supercell& scell = _pc.getSupercell();
  const vector<_cluster>& clusters = _pc._clusters[o3].clusters;
  uint nclusters = clusters.size();
  vector<double> invmasses(nclusters);
  for (uint c = 0; c < nclusters; c++) {
    double mass = 1.0;
    for (int o = 0; o < order; o++) mass *= scell.getAtomMass(clusters[c].atoms[o]);
    invmasses[c] = 1/sqrt(mass);
  }

  const vector<vector<double> >& ifcs = _pc._anharmonicIFCs[o3].force_constants;
  vector<vector<int> > cart_indices;
  aurostd::xcombos cart(3, order, 'E', true);
  while (cart.increment()) cart_indices.push_back(cart.getCombo());
  uint ncart = cart_indices.size();

  // Units are chosen so that probabilities are in THz (1/ps)
  //double probability_prefactor = std::pow(hbar_amu, o2) * PI/(std::pow(nQPs, o2) * std::pow(2.0, o1));
  double probability_prefactor = std::pow(hbar_amu, o2) * PI/(std::pow(nQPs, o3) * std::pow(2.0, o1));
  probability_prefactor *= std::pow(au2THz * 10.0, 2);

  int natoms = (int) _pc.getInputCellStructure().atoms.size();
  vector<int> atpowers(order, 1);
  vector<vector<int> > at_eigen;
  aurostd::xcombos at_combos(natoms, order, 'E' , true);
  while (at_combos.increment()) at_eigen.push_back(at_combos.getCombo());
  uint nateigen = at_eigen.size();
  vector<vector<xcomplex<double> > > eigenprods(nateigen, vector<xcomplex<double> >(ncart));
  for (int o = o1; o > 0; o--) atpowers[o1 - o] = (int) std::pow(natoms, o);

  xcomplex<double> matrix, prefactor, eigen;
  vector<int> qpts(order), branches(order);
  vector<bool> signs(o2);
  int iat, o, e;
  uint c, crt;
  double transprob, weight;

  string transfile, procfile, weightfile, qptfile;

  int nprocs = nIQPs;
  for (int i = 3; i < order; i++) nprocs *= nQPs;
  
  for (int i = startIndex; i < endIndex; i++) {
    ifstream procfl, weightfl, qptfl;
    ofstream transfl;
    qpts[0] = _qm.getIbzqpts()[i];

    transfile = tmpdir + "/trans." + aurostd::utype2string<int>(order) + "." + aurostd::utype2string<int>(i);
    procfile = tmpdir + "/proc." + aurostd::utype2string<int>(order) + "." + aurostd::utype2string<int>(i);
    qptfile = tmpdir + "/qpt." + aurostd::utype2string<int>(order) + "." + aurostd::utype2string<int>(i);
    weightfile = tmpdir + "/wght." + aurostd::utype2string<int>(order) + "." + aurostd::utype2string<int>(i);
    openTmpFile(qptfl, qptfile);
    openTmpFile(procfl, procfile);
    openTmpFile(transfl, transfile);
    openTmpFile(weightfl, weightfile);

    while (weightfl.peek() != EOF) {
      weightfl.read((char*)(&weight), sizeof(double));
      getProcess(order, procfl, qptfl, signs, qpts, branches);

      // Precompute eigenvalue products
      for (c = 0; c < nateigen; c++) {
        for (crt = 0; crt < ncart; crt++) {
          e = at_eigen[c][0] * 3 + cart_indices[crt][0] + 1;
          eigen = eigenvectors[qpts[0]][e][branches[0] + 1];
          for (o = 1; o < o1 ; o++) {
            e = at_eigen[c][o] * 3 + cart_indices[crt][o] + 1;
            if (signs[o - 1]) eigen *= conj(eigenvectors[qpts[o]][e][branches[o] + 1]);
            else eigen *= eigenvectors[qpts[o]][e][branches[o] + 1];
          }
          e = at_eigen[c][o1] * 3 + cart_indices[crt][o1] + 1;
          eigen *= conj(eigenvectors[qpts[o1]][e][branches[o1] + 1]);
          eigenprods[c][crt] = eigen;
        }
      }

      matrix.re = 0.0;
      matrix.im = 0.0;
      for (c = 0; c < nclusters; c++) {
        const vector<int>& atoms = clusters[c].atoms;
        iat = scell.sc2pcMap(atoms[0]);
        prefactor = invmasses[c] * phases[iat][atoms[o1]][qpts[o1]][1];
        for (o = 1; o < o1; o++) prefactor *= phases[iat][atoms[o]][qpts[o]][(int) signs[o - 1]];

        e = 0;
        for (o = 0; o < order; o++) e += scell.sc2pcMap(atoms[o]) * atpowers[o];
        for (crt = 0; crt < ncart; crt++) {
          if (std::abs(ifcs[c][crt]) > _ZERO_TOL_) {
            // Perform multiplication expliclty in place instead of using xcomplex.
            // This three times as fast because constructors and destructors are not called.
            matrix.re += ifcs[c][crt] * (prefactor.re * eigenprods[e][crt].re - prefactor.im * eigenprods[e][crt].im);
            matrix.im += ifcs[c][crt] * (prefactor.re * eigenprods[e][crt].im + prefactor.im * eigenprods[e][crt].re);
          }
        }
      }
      transprob = probability_prefactor * magsqr(matrix) * weight;
      for (o = 0; o < order; o++) transprob /= freq[qpts[o]][branches[o]];
      transfl.write((char*)(&transprob), sizeof(double));
    }
    aurostd::RemoveFile(weightfile);
    closeTmpFile(procfl, procfile);
    closeTmpFile(qptfl, qptfile);
    closeTmpFile(transfl, transfile);
    _logger.updateProgressBar(1.0/nprocs);
  }
}

void TCONDCalculator::getProcess(ifstream& procfl, vector<int>& qpts, vector<int>& branches) {
  int r;
  procfl.read((char*)(&r), sizeof(int));
  qpts[1] = r;
  procfl.read((char*)(&r), sizeof(int));
  branches[0] = r/nBranches;
  branches[1] = r % nBranches;
}

void TCONDCalculator::getProcess(int order, ifstream& procfl, ifstream& qptfl,
                                 vector<bool>& signs, vector<int>& qpts, vector<int>& branches) {
  const int o1 = order - 1;
  const int o2 = order - 2;
  const int o3 = order - 3;
  int s = 0, q, r;
  bool sign;
  for (q = 0; q < o2; q++) {
    procfl.read((char*)(&sign), sizeof(bool));
    signs[q] = sign;
    if (sign) s += nQPs;
  }
  procfl.read((char*)(&r), sizeof(int));
  qptfl.seekg(sizeof(int) * (s + r));
  qptfl.read((char*)(&qpts[o1]), sizeof(int));
  for (q = 0; q < o3; q++) {
    qpts[q + 1] = r/(int) std::pow(nQPs, o3 - q);
    r = r % (int) std::pow(nQPs, o3 - q);
  }
  qpts[o2] = r;
  procfl.read((char*)(&r), sizeof(int));
  for (q = 0; q < o1; q++) {
    branches[q] = r/(int) std::pow(nBranches, o1 - q);
    r = r % (int) std::pow(nBranches, o1 - q);
  }
  branches[o1] = r;
}

void TCONDCalculator::calculateTransitionProbabilitiesIsotope(const LTMethod& _lt, int startIndex,
                                                              int endIndex, vector<vector<double> >& rates) {
  const xstructure& pcell = _pc.getInputCellStructure();
  uint natoms = pcell.atoms.size();
  vector<double> pearson(natoms);
  uint at;
  for (at = 0; at < natoms; at++) pearson[at] = GetPearsonCoefficient(pcell.atoms[at].atomic_number);

  vector<double> freqs(nQPs), func(nQPs), weights(nQPs);
  double prefactor, rate;
  int q1, q2, br1, br2, b, e;
  xcomplex<double> eig;
  xvector<xcomplex<double> > eigen1(3), eigen2(3);

  string transfile, procfile;
  ofstream transfl, procfl;
  for (int iq = startIndex; iq < endIndex; iq++) {
    transfile = tmpdir + "/trans.iso." + aurostd::utype2string<int>(iq);
    procfile = tmpdir + "/proc.iso." + aurostd::utype2string<int>(iq);
    openTmpFile(transfl, transfile);
    openTmpFile(procfl, procfile);
    q1 = _qm.getIbzqpts()[iq];
    for (br1 = 0; br1 < nBranches; br1++) {
      //std::cout << "br1 = " << br1 << std::endl;
      //std::cout << freq[q1][br1] << std::endl;
      prefactor = freq[q1][br1] * freq[q1][br1] * PI/2.0;
      for (br2 = 0; br2 < nBranches; br2++) {
        for (q2 = 0; q2 < nQPs; q2++) {
          freqs[q2] = freq[q2][br2];
          func[q2] = 0.0;
          for (at = 0; at < natoms; at++) {
            if (pearson[at] != 0) {
              eig.re = 0.0;
              eig.im = 0.0;
              e = 3 * at;
              for (int i = 1; i < 4; i++) {
                eig.re += eigenvectors[q1][e + i][br1 + 1].re * eigenvectors[q2][e + i][br2 + 1].re;
                eig.re += eigenvectors[q1][e + i][br1 + 1].im * eigenvectors[q2][e + i][br2 + 1].im;
                eig.im += eigenvectors[q1][e + i][br1 + 1].im * eigenvectors[q2][e + i][br2 + 1].re;
                eig.im -= eigenvectors[q1][e + i][br1 + 1].re * eigenvectors[q2][e + i][br2 + 1].im;
              }
              func[q2] += pearson[at] * magsqr(eig);
            }
          }
        }
        weights = getWeightsLT(_lt, freq[q1][br1], freqs);
        for (q2 = 0; q2 < nQPs; q2++) {
          rate = prefactor * weights[q2] * func[q2];
          if (rate > _ZERO_TOL_) {
            rates[iq][br1] += rate;
            //if(br1 == 4) std::cout << q2 << " " << br1 << " " << br2 << " " << weights[q2] << " " << prefactor << " " << func[q2] << " " << rate << std::endl;
            transfl.write((char*)(&rate), sizeof(double));
            b = nBranches * br1 + br2;
            procfl.write((char*)(&q2), sizeof(int));
            procfl.write((char*)(&b), sizeof(int));
          }
        }
      }
//      for (int i = 0; i < nBranches; i++) std::cout << rates[iq][i] << " ";
//      std::cout << std::endl;
    }
    closeTmpFile(transfl, transfile);
    closeTmpFile(procfl, procfile);
    _logger.updateProgressBar(1.0/nIQPs);
  }
}

void TCONDCalculator::calculateTransitionProbabilitiesBoundary() {
  string transfile = tmpdir + "/trans.boundary";
  ofstream transfl;
  openTmpFile(transfl, transfile);
  int q;
  double rate;
  vector<vector<double> > rates(nIQPs, vector<double>(nBranches));
  double grain_size = aurostd::string2utype<double>(calc_options.getattachedscheme("GRAIN_SIZE"));
  _logger << "Calculating grain boundary transition probabilities with a grain size of " << grain_size << " nm." << apl::endl;
  for (int iq = 0; iq < nIQPs; iq++) {
    for (int br = 0; br < nBranches; br++) {
      q = _qm.getIbzqpts()[iq];
      rate = aurostd::modulus(gvel[q][br])/grain_size;
      rates[iq][br] = rate;
      transfl.write((char*)(&rate), sizeof(double));
    }
  }
  closeTmpFile(transfl, transfile);
  writeRatesToTmpFile("rates.boundary", rates);
}

}  // namespace apl

/************************************ BTE ***********************************/

namespace apl {

xmatrix<double> TCONDCalculator::calculateThermalConductivityTensor(double T,
                                                                    const vector<vector<int> >& invar_symops) {
  _logger << "Calculating thermal conductivity for " << T << " K." << apl::endl;
  vector<vector<double> > occ = getOccupationNumbers(T);
  for (int i = 0; i < nQPs; i++) {
    std::cout << i;
    for (int j = 0; j < nBranches; j++) {
      std::cout << " " << occ[i][j];
    }
    std::cout << std::endl;
  }
  _logger << "Calculating scattering rates." << apl::endl;
  int max_order;
  if (calc_options.flag("FOURTH_ORDER")) max_order = 4;
  else max_order = 3;
  for (int o = 3; o <= max_order; o++) calculateAnharmonicRates(o, T, occ);

  _logger << "Calculating RTA" << apl::endl;
  vector<vector<double> > rates = getRates(T);
  std::cout << "rates" << std::endl;
  vector<vector<xvector<double> > > mfd = getMeanFreeDispRTA(rates);
  std::cout << "mfd" << std::endl;
  xmatrix<double> tcond = calcTCOND(T, occ, mfd); // RTA solution
  std::cout << tcond << std::endl;

  if (0 && !calc_options.flag("RTA")) {
    xmatrix<double> tcond_prev(3, 3), diff(3, 3);
    int num_iter = 1;
    double norm;
    _logger << "Begin SCF for the Boltzmann transport equation." << apl::endl;
    std::cout << std::setiosflags(std::ios::fixed | std::ios::right);
    std::cout << std::setw(15) << "Iteration";
    std::cout << std::setiosflags(std::ios::fixed | std::ios::right);
    std::cout << std::setw(25) << "Norm" << std::endl;
    do {
      tcond_prev = tcond;
      getMeanFreeDispFull(rates, occ, invar_symops, mfd);
      tcond = calcTCOND(T, occ, mfd);
      norm = frobenius_norm(tcond_prev - tcond);
      std::cout << std::setiosflags(std::ios::fixed | std::ios::right);
      std::cout << std::setw(15) << num_iter;
      std::cout << std::setiosflags(std::ios::fixed | std::ios::right);
      std::cout << std::setw(25) << std::dec << (norm) << std::endl;
    } while((std::abs(norm) < 1e-5) && (num_iter <= max_iter));
    if (num_iter > max_iter) {
      string function = _AAPL_TCOND_ERR_PREFIX_ + "calculateThermalConductivityTensor()";
      stringstream message;
      message << "Thermal conductivity did not converge within " << max_iter << " iterations ";
      message << "at " << T << " K.";
      throw xerror(function, message, _RUNTIME_ERROR_);
    }
  }
  return tcond;
}

vector<vector<double> > TCONDCalculator::getOccupationNumbers(double T) {
  vector<vector<double> > occ(nQPs, vector<double>(nBranches, 0.0));
  for (int q = 0; q < nQPs; q++) {
    for (int br = 0; br < nBranches; br++) {
      occ[q][br] = 1.0/(exp(BEfactor * freq[q][br]/T) - 1.0);
    }
  }
  return occ;
}

void TCONDCalculator::calculateAnharmonicRates(int order, double T,
                                               const vector<vector<double> >& occ) {
  vector<vector<double> > rates(nIQPs, vector<double>(nBranches, 0.0));
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  int ncpus;
  _pc.get_NCPUS(ncpus);
  vector<std::thread*> threads;
  vector<vector<int> > thread_dist = getThreadDistribution(nIQPs, ncpus);
  threads.clear();
  for (int icpu = 0; icpu < ncpus; icpu++) {
    threads.push_back(new std::thread(&TCONDCalculator::calcAnharmRates, this,
                                      thread_dist[icpu][0], thread_dist[icpu][1],
                                      order, std::ref(occ), std::ref(rates)));
  }
  for (int icpu = 0; icpu < ncpus; icpu++) {
    threads[icpu]->join();
    delete threads[icpu];
  }
#else
  calcAnharmRates(0, nIQPs, order, occ, rates);
#endif
  for (int i = 0; i < nIQPs; i++) {
    std::cout << i;
    for (int j = 0; j < nBranches; j++) {
      std::cout << " " << rates[i][j];
    }
    std::cout << std::endl;
  }
  string ratefile = "rates." + aurostd::utype2string<int>(order) + "." + aurostd::utype2string<double>(T);;
  writeRatesToTmpFile(ratefile, rates);
}

void TCONDCalculator::calcAnharmRates(int startIndex, int endIndex, int order,
                                      const vector<vector<double> >& occ, vector<vector<double> >& rates) {
  vector<bool> signs(order - 2);
  vector<int> qpts(order), branches(order);
  string procfile, transfile, qptfile;
  double trans;
  for (int i = startIndex; i < endIndex; i++) {
    ifstream procfl, transfl, qptfl;
    qpts[0] = _qm.getIbzqpts()[i];
    procfile = tmpdir + "/" + "proc." + aurostd::utype2string<int>(order) + "." + aurostd::utype2string<int>(i);
    transfile = tmpdir + "/" + "trans." + aurostd::utype2string<int>(order) + "." + aurostd::utype2string<int>(i);
    qptfile = tmpdir + "/" + "qpt." + aurostd::utype2string<int>(order) + "." + aurostd::utype2string<int>(i);
    openTmpFile(procfl, procfile);
    openTmpFile(transfl, transfile);
    openTmpFile(qptfl, qptfile);
    int j = 0;
    while (transfl.peek() != EOF) {
      transfl.read((char*)(&trans), sizeof(double));
      getProcess(order, procfl, qptfl, signs, qpts, branches);
      j++;
      rates[i][branches[0]] += trans * getOccupationTerm(order, occ, signs, qpts, branches);
    }
    closeTmpFile(procfl, procfile);
    closeTmpFile(qptfl, qptfile);
    closeTmpFile(transfl, transfile);
  }
}

double TCONDCalculator::getOccupationTerm(int order,
                                          const vector<vector<double> >& occ,
                                          const vector<bool>& signs,
                                          const vector<int>& qpts,
                                          const vector<int>& branches) {
  double term = 0;
  if (order == 3) {
    if (signs[0]) {
      term = (1.0 + occ[qpts[1]][branches[1]] + occ[qpts[2]][branches[2]])/2.0;
    } else {
      term = occ[qpts[1]][branches[1]] - occ[qpts[2]][branches[2]];
    }
  } else if (order == 4) {
    if (!signs[0] && !signs[1]) {  // ++
      term = occ[qpts[1]][branches[1]] * occ[qpts[2]][branches[2]] * occ[qpts[3]][branches[3]]/6.0;
    } else if (signs[0] && signs[1]) {  // --
      term = (1 + occ[qpts[1]][branches[1]]) * (1 + occ[qpts[2]][branches[2]]) * occ[qpts[3]][branches[3]]/2.0;
    } else {  // +-
      term = (1 + occ[qpts[1]][branches[1]]) * occ[qpts[2]][branches[2]] * occ[qpts[3]][branches[3]]/2.0;
    }
    term /= occ[qpts[0]][branches[0]];
  }
  return term;
}

vector<vector<double> > TCONDCalculator::getRates(double T) {
  vector<vector<double> > rates(nIQPs, vector<double>(nBranches, 0.0));

  vector<vector<double> > rtmp = readRatesFromTmpFile("rates.3." + aurostd::utype2string<double>(T));
  for (int iq = 0; iq < nIQPs; iq++) {
    for (int br = 0; br < nBranches; br++) {
      rates[iq][br] += rtmp[iq][br];
    }
  }

  if (calc_options.flag("FOURTH_ORDER")) {
    rtmp = readRatesFromTmpFile("rates.4." + aurostd::utype2string<double>(T));
    for (int iq = 0; iq < nIQPs; iq++) {
      for (int br = 0; br < nBranches; br++) {
        rates[iq][br] += rtmp[iq][br];
      }
    }
  }

  if (0 && calc_options.flag("ISOTOPE")) {
    rtmp = readRatesFromTmpFile("rates.iso");
    for (int iq = 0; iq < nIQPs; iq++) {
      for (int br = 0; br < nBranches; br++) {
        rates[iq][br] += rtmp[iq][br];
      }
    }
  }

  if (calc_options.flag("BOUNDARY")) {
    rtmp = readRatesFromTmpFile("rates.boundary");
    for (int iq = 0; iq < nIQPs; iq++) {
      for (int br = 0; br < nBranches; br++) {
        rates[iq][br] += rtmp[iq][br];
      }
    }
  }

  return rates;
}

vector<vector<xvector<double> > > TCONDCalculator::getMeanFreeDispRTA(const vector<vector<double> >& rates) {
  xvector<double> xvec(3);
  vector<vector<xvector<double> > > mfd(nQPs, vector<xvector<double> >(nBranches, xvec));
  int iq;
  for (int q = 0; q < nQPs; q++) {
    iq = _qm.getIbzqpt(q);
    for (int br = 0; br < nBranches; br++) {
      if (rates[iq][br] > 0.0) {
        mfd[q][br] = gvel[q][br] * freq[q][br]/rates[iq][br];
      }
    }
  }
  return mfd;
}

xmatrix<double> TCONDCalculator::calcTCOND(double T, const vector<vector<double> >& occ,
                                           const vector<vector<xvector<double> > >& mfd) {
  xmatrix<double> tcond(3, 3);
  bool cumulative = calc_options.flag("CUMULATIVEK");
  double grain_size = aurostd::string2utype<double>(calc_options.getattachedscheme("GRAIN_SIZE"));
  double prefactor = 1E24 * std::pow(hbar_J, 2)/(KBOLTZ * std::pow(T, 2) * nQPs * Volume(_pc.getInputCellStructure()));
  for (int q = 0; q < nQPs; q++) {
    for (int br = 0; br < nBranches; br++) {
      bool include = true;
      if (freq[q][br] < _AFLOW_APL_EPS_) {
        include = false;
      } else if (cumulative) {
        double mfpath = scalar_product(mfd[q][br], gvel[q][br])/aurostd::modulus(gvel[q][br]);
        include = !(mfpath > grain_size);
      }
      if (include) {
        double x = occ[q][br] * (occ[q][br] + 1) * freq[q][br];
        for (int i = 1; i < 4; i++) {
          for (int j = 1; j < 4; j++) {
            tcond[i][j] += x * gvel[q][br][i] * mfd[q][br][j];
          }
        }
      }
    }
  }
  tcond = prefactor * tcond;
  return tcond;
}

void TCONDCalculator::getMeanFreeDispFull(const vector<vector<double> >& rates,
                                          const vector<vector<double> >& occ,
                                          const vector<vector<int> >& invar_symops,
                                          vector<vector<xvector<double> > >& mfd) {
  // MPI variables
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  int ncpus;
  _pc.get_NCPUS(ncpus);
  vector<vector<int> > thread_dist;
  vector<std::thread*> threads;
#endif

  vector<vector<xvector<double> > > delta(nIQPs, vector<xvector<double> >(nBranches, xvector<double>(3)));
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  thread_dist = getThreadDistribution(nIQPs, ncpus);
  threads.clear();
  for (int icpu = 0; icpu < ncpus; icpu++) {
    threads.push_back(new std::thread(&TCONDCalculator::calculateDelta, this,
                                      thread_dist[icpu][0], thread_dist[icpu][1], 
                                      std::ref(invar_symops), std::ref(occ), std::ref(mfd), std::ref(delta)));
  }
  for (int icpu = 0; icpu < ncpus; icpu++) {
    threads[icpu]->join();
    delete threads[icpu];
  }
#else
  calculateDelta(0, nIQPs, invar_symops, occ, mfd, delta);
#endif

  correctMFD(rates, delta, mfd);
}

void TCONDCalculator::calculateDelta(int startIndex, int endIndex, 
                                     const vector<vector<int> >& invar_symops,
                                     const vector<vector<double> >& occ,
                                     const vector<vector<xvector<double> > >& mfd,
                                     vector<vector<xvector<double> > >& delta) {
  xvector<double> correction(3);
  double trans;
  string procfile, transfile, qptfile;
  int max_order;
  if (calc_options.flag("FOURTH_ORDER")) max_order = 4;
  else max_order = 3;

  for (int i = startIndex; i < endIndex; i++) {
    for (int o = 3; o <= max_order; o++) {
      const int o2 = o - 2;
      vector<int> qpts(o), branches(o);
      vector<bool> signs(o2);
      qpts[0] = _qm.getIbzqpts()[i];
      ifstream procfl, transfl, qptfl;
      procfile = tmpdir + "/proc." + aurostd::utype2string<int>(o) + "." + aurostd::utype2string<int>(i);
      transfile = tmpdir + "/trans." + aurostd::utype2string<int>(o) + "." + aurostd::utype2string<int>(i);
      qptfile = tmpdir + "/qpt." + aurostd::utype2string<int>(o) + "." + aurostd::utype2string<int>(i);
      openTmpFile(procfl, procfile);
      openTmpFile(transfl, transfile);
      openTmpFile(qptfl, qptfile);
      while (transfl.peek() != EOF) {
        transfl.read((char*)(&trans), sizeof(double));
        getProcess(o, procfl, qptfl, signs, qpts, branches);
        correction = mfd[qpts.back()][branches.back()];
        for (int j = 0; j < o2; j++) {
          if (signs[j]) correction -= mfd[qpts[j + 1]][branches[j + 1]];
          else correction += mfd[qpts[j + 1]][branches[j + 1]];
        }
        delta[qpts[0]][branches[0]] += trans * correction * getOccupationTerm(o, occ, signs, qpts, branches);
      }
      closeTmpFile(procfl, procfile);
      closeTmpFile(transfl, transfile);
      closeTmpFile(qptfl, qptfile);
    }

    if (calc_options.flag("ISOTOPE")) {
      ifstream procfl_iso, transfl_iso;
      procfile = tmpdir + "/proc.iso." + aurostd::utype2string<int>(i);
      transfile = tmpdir + "/trans.iso." + aurostd::utype2string<int>(i);

      vector<int> qpts_iso(2), branches_iso(2);
      qpts_iso[i] = _qm.getIbzqpts()[i];
      openTmpFile(procfl_iso, procfile);
      openTmpFile(transfl_iso, transfile);
      while (transfl_iso.peek() != EOF) {
        transfl_iso.read((char*)(&trans), sizeof(double));
        getProcess(procfl_iso, qpts_iso, branches_iso);
        delta[qpts_iso[0]][branches_iso[0]] += trans * mfd[qpts_iso[1]][branches_iso[1]];
      }
      closeTmpFile(procfl_iso, procfile);
      closeTmpFile(transfl_iso, transfile);
    }

    // Symmetrize
    int symop;
    const vector<_sym_op>& pgroup = _qm.getReciprocalCell().pgroup;
    xmatrix<double> Uc(3, 3);
    uint nsym = invar_symops[i].size();
    for (uint isym = 0; isym < nsym; isym++) {
      symop = invar_symops[i][isym];
      Uc += pgroup[symop].Uc;
    }
    Uc = 1.0/nsym * Uc;
    for (int br = 0; br < nBranches; br++) {
      delta[i][br] = Uc * delta[i][br];
    }
  }
}

void TCONDCalculator::correctMFD(const vector<vector<double> >& rates,
                                 const vector<vector<xvector<double> > >& delta,
                                 vector<vector<xvector<double> > >& mfd) {
  const vector<_sym_op>& pgroup = _qm.getReciprocalCell().pgroup;
  xvector<double> xvec(3);
  int iq;
  for (int q = 0; q < nQPs; q++) {
    const xmatrix<double>& Uc = pgroup[_qm.getQPoint(q).symop].Uc;
    iq = _qm.getIbzqpt(q);
    for (int br = 0; br < nBranches; br++) {
      if (rates[iq][br] > 0.0) mfd[q][br] = (freq[q][br] * gvel[q][br] + Uc * delta[iq][br])/rates[iq][br];
      else mfd[q][br] = xvec;
    }
  }
}

}  // namespace apl

/********************************* FILE I/O *********************************/

namespace apl {

void TCONDCalculator::openTmpFile(ofstream& file, const string& filename) {
  file.open(filename, std::ios::out | std::ios::binary);
}

void TCONDCalculator::openTmpFile(ifstream& file, const string& filename) {
  aurostd::UncompressFile(filename + "." + calc_options.getattachedscheme("KZIP_BIN"));
  file.open(filename, std::ios::in | std::ios::binary);
}

void TCONDCalculator::closeTmpFile(ofstream& file, const string& filename) {
  file.close();
  aurostd::CompressFile(filename, calc_options.getattachedscheme("KZIP_BIN"));
}

void TCONDCalculator::closeTmpFile(ifstream& file, const string& filename) {
  file.close();
  aurostd::CompressFile(filename, calc_options.getattachedscheme("KZIP_BIN"));
}

void TCONDCalculator::writeRatesToTmpFile(const string& ratefile, 
                                          const vector<vector<double> >& rates) {
  ofstream ratefl;
  string filename = tmpdir + "/" + ratefile;
  openTmpFile(ratefl, filename);
  for (int q = 0; q < nIQPs; q++) {
    for (int br = 0; br < nBranches; br++) {
      ratefl.write((char*)(&rates[q][br]), sizeof(double));
    }
  }
  closeTmpFile(ratefl, filename);
}

vector<vector<double> > TCONDCalculator::readRatesFromTmpFile(const string& ratefile) {
  ifstream ratefl;
  string filename = tmpdir + "/" + ratefile;
  vector<vector<double> > rates(nIQPs, vector<double>(nBranches, 0.0));
  openTmpFile(ratefl, filename);
  for (int q = 0; q < nIQPs; q++) {
    for (int br = 0; br < nBranches; br++) {
      ratefl.read((char*)(&rates[q][br]), sizeof(double));
    }
  }
  closeTmpFile(ratefl, filename);
  return rates;
}

void TCONDCalculator::writeTempIndepOutput(const string& filename, string keyword,
                                           const string& unit, const vector<vector<double> >& data) {
  stringstream output;
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  string key = "[AAPL_" + aurostd::toupper(aurostd::StringSubst(keyword, " ", "_")) + "]";
  if (!_pc.getSystemName().empty()) {
    output << key << "SYSTEM=" << _pc.getSystemName() << std::endl;
  }
  output << key << "START" << std::endl;
  output << std::setiosflags(std::ios::fixed | std::ios::right);
  output << std::setw(10) << "# Q-point"
         << std::setw(20) << " "
         << keyword << " (" << unit << ")" << std::endl;
  writeDataBlock(output, data);
  output << key << "STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeTempIndepOutput()";
    stringstream message;
    message << "Could not write file " << filename << ".";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

void TCONDCalculator::writeTempDepOutput(const string& filename, string keyword, const string& unit,
                                         const vector<double>& temps, const vector<vector<vector<double> > >& data) {
  stringstream output;
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  string key = "[AAPL_" + aurostd::toupper(aurostd::StringSubst(keyword, " ", "_")) + "]";
  if (!_pc.getSystemName().empty()) {
    output << key << "SYSTEM=" << _pc.getSystemName() << std::endl;
  }
  output << key << "START" << std::endl;
  for (uint t = 0; t < temps.size(); t++) {
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << key << "T = " << temps[t] << " K" << std::endl;
    output << std::setw(10) << "# Q-point"
           << std::setw(20) << " "
           << keyword << " (" << unit << ")" << std::endl;
    writeDataBlock(output, data[t]);
  }
  output << key << "STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeTempDepOutput()";
    stringstream message;
    message << "Could not write file " << filename << ".";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

void TCONDCalculator::writeDataBlock(stringstream& output,
                                     const vector<vector<double> >& data) {
  for (uint q = 0; q < data.size(); q++) {
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << std::setw(10) << q;
    for (uint br = 0; br < data[q].size(); br++) {
      output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      output << std::setw(20) << std::setprecision(10) << std::scientific << data[q][br];
    }
    output << std::endl;
  }
}

//writeFrequencies////////////////////////////////////////////////////////////
// Writes the frequencies into a file. Each row belongs to a q-point, and
// each column belongs to a phonon branch.
void TCONDCalculator::writeFrequencies() {
  stringstream output;
  string filename = DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_FREQ_FILE;

  // Header
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  if (!_pc.getSystemName().empty()) {
    output << "[AAPL_FREQUENCY]SYSTEM=" << _pc.getSystemName() << std::endl;
  }
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
  if (!_pc.getSystemName().empty()) {
    output << "[AAPL_GROUP_VELOCITY]SYSTEM=" << _pc.getSystemName() << std::endl;
  }
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
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeGroupVelocities()";
    string message = "Could not write group velocities to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

void TCONDCalculator::writeThermalConductivity() {
  stringstream output;
  string filename = DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_TCOND_FILE;

  // Header
  output << AFLOWIN_SEPARATION_LINE << std::endl;
  if (!_pc.getSystemName().empty()) {
    output << "[AAPL_THERMAL_CONDUCTIVITY]SYSTEM=" << _pc.getSystemName() << std::endl;
  }
  output << "[AAPL_THERMAL_CONDUCTIVITY]START" << std::endl;
  output << std::setw(8) << "# T (K)"
         << std::setw(75) << " "
         << "Thermal Conductivity (W/m*K)" << std::endl;
  string xyz[3] = {"x", "y", "z"};
  output << "#" << std::setw(8) << " ";
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      output << std::setw(7) << " ";
      // Write columns first to make compatible with gnuplot
      string k = "k(" + xyz[j] + "," + xyz[i] + ")";
      output << k;
      output << std::setw(7) << " ";
    }
  }
  output << std::endl;

  // Body
  for (uint t = 0; t < temperatures.size(); t++) {
    output << std::setw(8) << std::fixed << std::setprecision(2) << temperatures[t];
    for (int i = 1; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        output << std::setw(2) << " ";
        // Write columns first to make compatible with gnuplot
        output << std::setprecision(10) << std::scientific << thermal_conductivity[t][j][i];
        output << std::setw(2) << " ";
      }
    }
    output << std::endl;
  }
  output << "[AAPL_THERMAL_CONDUCTIVITY]STOP" << std::endl;
  output << AFLOWIN_SEPARATION_LINE << std::endl;

  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_TCOND_ERR_PREFIX_ + "writeThermalConductivity";
    string message = "Could not write thermal conductivities to file.";
    throw xerror(function, message, _FILE_ERROR_);
  }
}

}  // namespace apl
