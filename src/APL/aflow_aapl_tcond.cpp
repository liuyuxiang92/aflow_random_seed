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
#define ALAMODE

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
static const double scatt_multi[2] = {1.0, 0.5};

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
  intr_trans_probs_iso.clear();
  nBranches = 0;
  nIQPs = 0;
  nQPs = 0;
  processes.clear();
  processes_iso.clear();
  rates_boundary.clear();
  rates_isotope.clear();
  temperatures.clear();
  thermal_conductivity.clear();
}

}  // namespace apl

/******************************* MAIN FUNCTION ******************************/

namespace apl {

void TCONDCalculator::calculateThermalConductivity() {
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

  calculateTransitionProbabilities();

  vector<vector<int> > small_groups = calculateSmallGroups();
  thermal_conductivity.assign(temperatures.size(), xmatrix<double>(3, 3));
  for (uint t = 0; t < temperatures.size(); t++) {
    thermal_conductivity[t] = calculateThermalConductivityTensor(temperatures[t], small_groups);
  }
}

vector<vector<int> > TCONDCalculator::calculateSmallGroups() {
  vector<vector<int> > small_groups(nIQPs, vector<int>(1, 0));  // Identity is always invariant
  const vector<_sym_op>& symops = _qm.getReciprocalCell().pgroup;
  const vector<int>& ibzqpts = _qm.getIbzqpts();
  int q;
  for (int iq = 0; iq < nIQPs; iq++) {
    q = ibzqpts[iq];
    const xvector<double>& fpos = _qm.getQPoint(q).fpos;
    for (uint isym = 1; isym < symops.size(); isym++) {
      if (_qm.getQPointIndex(symops[isym].Uf * fpos) == q) small_groups[iq].push_back(isym);
    }
  }
  return small_groups;
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

  // Prepare storage
  xmatrix<xcomplex<double> > eigen(nBranches, nBranches, 1, 1);
  eigenvectors.assign(nQPs, eigen);
  freq.assign(nQPs, vector<double>(nBranches));
  xvector<double> g(3);
  gvel.assign(nQPs, vector<xvector<double> >(nBranches, g));

  // Calculate frequencies and group velocities
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  int ncpus;
  _pc.get_NCPUS(ncpus);
  vector<vector<int> > thread_dist = setupMPI(message, _logger, nQPs, ncpus);
  vector<std::thread*> threads;
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

/********************************** WEIGHTS *********************************/

namespace apl {

double fij(double ei, double ej, double e) {
  return (e - ej)/(ei - ej);
}

void TCONDCalculator::getWeightsLT(const LTMethod& _lt, double e_ref,
                                   const vector<double>& energy, vector<double>& weights) {
  for (int q = 0; q < nQPs; q++) weights[q] = 0;
  const vector<vector<int> >& corners = _lt.getTetrahedra();
  double vol = _lt.getVolumePerTetrahedron();
  int ntet = _lt.getnTetrahedra();

  double g, tmp;
  int i, j, ii, jj;

  double e_tmp[4];
  int sort_arg[4], kindex[4];

  for (i = 0; i < ntet; i++) {
    for (j = 0; j < 4; j++) {
      e_tmp[j] = energy[corners[i][j]];
      kindex[j] = corners[i][j];
    }

    for (ii = 0; ii < 4; ii++) sort_arg[ii] = ii;

    for (ii = 0; ii < 4; ii++) {
      tmp = e_tmp[ii];
      jj = ii;
      while (jj > 0 && tmp < e_tmp[jj - 1]) {
        e_tmp[jj] = e_tmp[jj - 1];
        sort_arg[jj] = sort_arg[jj - 1];
        jj--;
      }
      e_tmp[jj] = tmp;
      sort_arg[jj] = ii;
    }

    double e1 = e_tmp[0];
    double e2 = e_tmp[1];
    double e3 = e_tmp[2];
    double e4 = e_tmp[3];

    if (aurostd::isequal(e1, e4)) continue;

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
    weights[k1] += vol * I1;
    weights[k2] += vol * I2;
    weights[k3] += vol * I3;
    weights[k4] += vol * I4;
  }
}

}  // namespace apl

/************************* TRANSITION PROBABILITIES *************************/

namespace apl {

void TCONDCalculator::calculateTransitionProbabilities() {
  _logger << "Calculating transition probabilities." << apl::endl;
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  int ncpus;
  _pc.get_NCPUS(ncpus);
  vector<std::thread*> threads;
  vector<vector<int> > thread_dist;
#endif
  string message;
  LTMethod _lt(_qm, _logger);
  vector<vector<vector<vector<xcomplex<double> > > > > phases = calculatePhases();

  message = "Transition Probabilities";
  _logger << "Calculating transition probabilities for 3-phonon scattering processes." << apl::endl;
  processes.resize(nIQPs);
  intr_trans_probs.resize(nIQPs);
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
  thread_dist = setupMPI(message, _logger, nIQPs, ncpus);
  threads.clear();
  for (int icpu = 0; icpu < ncpus; icpu++) {
    threads.push_back(new std::thread(&TCONDCalculator::calculateTransitionProbabilitiesPhonon, this,
                                      thread_dist[icpu][0], thread_dist[icpu][1], std::ref(_lt), std::ref(phases)));
  }
  finishMPI(threads, _logger);
#else
  _logger.initProgressBar(message);
  calculateTransitionProbabilitiesPhonon(0, nIQPs, phases);
  _logger.finishProgressBar();
#endif

  if (0 && calc_options.flag("ISOTOPE")) {
    _logger << "Calculating isotope transition probabilities." << apl::endl;
    message = "Isotope Transition Probabilities";
    processes_iso.resize(nIQPs);
    intr_trans_probs_iso.resize(nIQPs);

#ifdef AFLOW_APL_MULTITHREADS_ENABLE    
    thread_dist = setupMPI(message, _logger, nIQPs, ncpus);
    threads.clear();
    for (int icpu = 0; icpu < ncpus; icpu++) {
      threads.push_back(new std::thread(&TCONDCalculator::calculateTransitionProbabilitiesIsotope, this,
                                        thread_dist[icpu][0], thread_dist[icpu][1], std::ref(_lt)));
    }
    finishMPI(threads, _logger);
#else
    _logger.initProgressBar(message);
    calculateTransitionProbabilitiesIsotope(0, nIQPs, _lt);
    _logger.finishProgressBar();
#endif
  }

  if (calc_options.flag("BOUNADRY")) {
    rates_boundary = calculateTransitionProbabilitiesBoundary();
  }
}

vector<vector<vector<vector<xcomplex<double> > > > > TCONDCalculator::calculatePhases() {
  const xstructure& scell = _pc.getSuperCellStructure();
  const xstructure& pcell = _pc.getInputCellStructure();
  const vector<int>& sc2pcMap = _pc.getSupercell()._sc2pcMap;
  const vector<int>& pc2scMap = _pc.getSupercell()._pc2scMap;
  uint niatoms = pcell.atoms.size();
  uint natoms = scell.atoms.size();
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
        phases[iat][at][q][0] = exp(iONE * scalar_product(_qm.getQPoint(q).cpos, min_vec));
        phases[iat][at][q][1] = conj(phases[iat][at][q][0]);
      }
    }
  }
  return phases;
}

void TCONDCalculator::calculateTransitionProbabilitiesPhonon(int startIndex, int endIndex, const LTMethod& _lt,
                                                             const vector<vector<vector<vector<xcomplex<double> > > > >& phases) {
  const Supercell& scell = _pc.getSupercell();
  const vector<_cluster>& clusters = _pc._clusters[0].clusters;
  uint nclusters = clusters.size();
  vector<double> invmasses(nclusters);
  for (uint c = 0; c < nclusters; c++) {
    double mass = 1.0;
    for (int o = 0; o < 3; o++) mass *= scell.getAtomMass(clusters[c].atoms[o]);
    invmasses[c] = 1/sqrt(mass);
  }

  const vector<vector<double> >& ifcs = _pc._anharmonicIFCs[0].force_constants;
  vector<vector<int> > cart_indices;
  aurostd::xcombos cart(3, 3, 'E', true);
  while (cart.increment()) cart_indices.push_back(cart.getCombo());
  uint ncart = cart_indices.size();

  // Units are chosen so that probabilities are in THz (1/ps)
  const double probability_prefactor = std::pow(au2THz * 10.0, 2) * hbar_amu * PI/4.0;
  const double ps_prefactor = 2.0/(3.0 * std::pow(nBranches, 3) * nQPs);

  int natoms = (int) _pc.getInputCellStructure().atoms.size();
  vector<int> atpowers(3, 1);
  vector<vector<int> > at_eigen;
  aurostd::xcombos at_combos(natoms, 3, 'E' , true);
  while (at_combos.increment()) at_eigen.push_back(at_combos.getCombo());
  uint nateigen = at_eigen.size();
  vector<vector<xcomplex<double> > > eigenprods(nateigen, vector<xcomplex<double> >(ncart));
  for (int o = 0; o < 2; o++) atpowers[o] = (int) std::pow(natoms, 2 - o);

  vector<vector<int> > branches;
  aurostd::xcombos branch_combos(nBranches, 3, 'E', true);
  while (branch_combos.increment()) branches.push_back(branch_combos.getCombo());
  uint nbr = branches.size();

  xcomplex<double> matrix, prefactor, eigen;
  vector<double> weights(nQPs), frequencies(nQPs);
  vector<int> qpts(3), proc(3), lastq(nQPs);
  int iat, o, e, q, s;
  uint c, crt, br;
  double transprob, freq_ref;
  bool calc;

  for (int i = startIndex; i < endIndex; i++) {
    qpts[0] = _qm.getIbzqpts()[i];

    for (s = 0; s < 2; s++) {
      for (q = 0; q < nQPs; q++) {
        if (s) lastq[q] = _qm.getQPointIndex(_qm.getQPoint(qpts[0]).fpos - _qm.getQPoint(q).fpos);
        else lastq[q] = _qm.getQPointIndex(_qm.getQPoint(qpts[0]).fpos + _qm.getQPoint(q).fpos);
      }
      for (br = 0; br < nbr; br++) {
        freq_ref = freq[qpts[0]][branches[br][0]];
        for (q = 0; q < nQPs; q++) {
          if (s) frequencies[q] = freq[q][branches[br][1]] + freq[lastq[q]][branches[br][2]];
          else frequencies[q] = -freq[q][branches[br][1]] + freq[lastq[q]][branches[br][2]];
        }
        
        getWeightsLT(_lt, freq_ref, frequencies, weights);

        for (q = 0; q < nQPs; q++) {
          calc = (weights[q] > _ZERO_TOL_);
          if (calc) {
            qpts[1] = q;
            qpts[2] = lastq[q];
            for (o = 0; o < 3; o++) {
              if (freq[qpts[o]][branches[br][o]] < _ZERO_TOL_) break;
            }
            calc = (o == 3);
          }

          if (calc) {
            // Precompute eigenvalue products
            for (c = 0; c < nateigen; c++) {
              for (crt = 0; crt < ncart; crt++) {
                e = at_eigen[c][0] * 3 + cart_indices[crt][0] + 1;
                eigen = eigenvectors[qpts[0]][e][branches[br][0] + 1];
                e = at_eigen[c][o] * 3 + cart_indices[crt][o] + 1;
                if (s) eigen *= conj(eigenvectors[qpts[o]][e][branches[br][o] + 1]);
                else eigen *= eigenvectors[qpts[o]][e][branches[br][o] + 1];
                e = at_eigen[c][2] * 3 + cart_indices[crt][2] + 1;
                eigen *= conj(eigenvectors[qpts[2]][e][branches[br][2] + 1]);
                eigenprods[c][crt] = eigen;
              }
            }
      
            matrix.re = 0.0;
            matrix.im = 0.0;
            for (c = 0; c < nclusters; c++) {
              const vector<int>& atoms = clusters[c].atoms;
              iat = scell.sc2pcMap(atoms[0]);
              prefactor = invmasses[c] * phases[iat][atoms[2]][qpts[2]][1];
              for (o = 1; o < 2; o++) prefactor *= phases[iat][atoms[o]][qpts[o]][s];
              e = 0;
              for (o = 0; o < 3; o++) e += scell.sc2pcMap(atoms[o]) * atpowers[o];
              for (crt = 0; crt < ncart; crt++) {
                // Perform multiplication expliclty in place instead of using xcomplex.
                // This is three times as fast because constructors and destructors are not called.
                matrix.re += ifcs[c][crt] * (prefactor.re * eigenprods[e][crt].re - prefactor.im * eigenprods[e][crt].im);
                matrix.im += ifcs[c][crt] * (prefactor.re * eigenprods[e][crt].im + prefactor.im * eigenprods[e][crt].re);
              }
            }
            transprob = probability_prefactor * magsqr(matrix) * weights[q];
            for (o = 0; o < 3; o++) transprob /= freq[qpts[o]][branches[br][o]];
            if (transprob > _ZERO_TOL_) {
              proc[0] = s;
              proc[1] = q;
              proc[2] = br;
              processes[i].push_back(proc);
              intr_trans_probs[i].push_back(transprob);
            }
          }
        }
      }
    }
    _logger.updateProgressBar(1.0/nIQPs);
  }
}

void TCONDCalculator::calculateTransitionProbabilitiesIsotope(int startIndex, int endIndex, const LTMethod& _lt) {
  const xstructure& pcell = _pc.getInputCellStructure();
  uint natoms = pcell.atoms.size();
  vector<double> pearson(natoms);
  uint at;
  for (at = 0; at < natoms; at++) pearson[at] = GetPearsonCoefficient(pcell.atoms[at].atomic_number);

  vector<double> frequencies(nQPs), func(nQPs), weights(nQPs);
  vector<int> proc(2);
  double prefactor, rate;
  int q1, q2, br1, br2, b, e;
  xcomplex<double> eig;
  xvector<xcomplex<double> > eigen1(3), eigen2(3);

  for (int i = startIndex; i < endIndex; i++) {
    q1 = _qm.getIbzqpts()[i];
    for (br1 = 0; br1 < nBranches; br1++) {
      prefactor = freq[q1][br1] * freq[q1][br1] * PI/2.0;
      for (br2 = 0; br2 < nBranches; br2++) {
        for (q2 = 0; q2 < nQPs; q2++) {
          frequencies[q2] = freq[q2][br2];
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
        getWeightsLT(_lt, freq[q1][br1], frequencies, weights);
        for (q2 = 0; q2 < nQPs; q2++) {
          rate = prefactor * weights[q2] * func[q2];
          if (rate > _ZERO_TOL_) {
            b = nBranches * br1 + br2;
            proc[0] = q2;
            proc[1] = b;
            intr_trans_probs_iso[i].push_back(rate);
            processes_iso[i].push_back(proc);
            rates_isotope[i][br1] += rate;
          }
        }
      }
    }
    _logger.updateProgressBar(1.0/nIQPs);
  }
}

vector<vector<double> > TCONDCalculator::calculateTransitionProbabilitiesBoundary() {
  int br, iq, q;
  vector<vector<double> > rates(nIQPs, vector<double>(nBranches));
  double grain_size = aurostd::string2utype<double>(calc_options.getattachedscheme("GRAIN_SIZE"));

  _logger << "Calculating grain boundary transition probabilities with a grain size of " << grain_size << " nm." << apl::endl;
  for (iq = 0; iq < nIQPs; iq++) {
    for (br = 0; br < nBranches; br++) {
      q = _qm.getIbzqpts()[iq];
      rates[iq][br] = aurostd::modulus(gvel[q][br])/grain_size;
    }
  }
  return rates;
}

}  // namespace apl

/************************************ BTE ***********************************/

namespace apl {

void TCONDCalculator::getProcess(const vector<int>& process,
                                 vector<int>& qpts, vector<int>& branches) {
  xvector<double> qsum = _qm.getQPoint(qpts[0]).fpos;
  qpts[1] = process[1];
  if (process[0]) qpts[2] = _qm.getQPointIndex(qsum - _qm.getQPoint(qpts[1]).fpos);
  else qpts[2] = _qm.getQPointIndex(qsum + _qm.getQPoint(qpts[1]).fpos);

  int pw, br;
  br = process[2];
  for (int i = 0; i < 2; i++) {
    pw = (int) std::pow(nBranches, 2 - i);
    branches[i] = br/pw;
    br = br % pw;
  }
  branches[2] = br;
}

xmatrix<double> TCONDCalculator::calculateThermalConductivityTensor(double T,
                                                                    const vector<vector<int> >& small_groups) {
  _logger << "Calculating thermal conductivity for " << T << " K." << apl::endl;
  vector<vector<double> > occ = getOccupationNumbers(T);
  _logger << "Calculating scattering rates." << apl::endl;
  vector<vector<double> > rates = getRates(occ);

  _logger << "Calculating RTA" << apl::endl;
  vector<vector<xvector<double> > > mfd = getMeanFreeDispRTA(rates);
  xmatrix<double> tcond = calcTCOND(T, occ, mfd); // RTA solution
  std::cout << "tcond: " << std::endl;
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
      getMeanFreeDispFull(rates, occ, small_groups, mfd);
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

vector<vector<double> > TCONDCalculator::calculateAnharmonicRates(const vector<vector<double> >& occ) {
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
                                      std::ref(occ), std::ref(rates)));
  }
  for (int icpu = 0; icpu < ncpus; icpu++) {
    threads[icpu]->join();
    delete threads[icpu];
  }
#else
  calcAnharmRates(0, nIQPs, occ, rates);
#endif
  for (int i = 0; i < nIQPs; i++) {
    std::cout << i;
    for (int j = 0; j < nBranches; j++) {
      std::cout << " " << rates[i][j];
    }
    std::cout << std::endl;
  }
  return rates;
}

void TCONDCalculator::calcAnharmRates(int startIndex, int endIndex,
                                      const vector<vector<double> >& occ,
                                      vector<vector<double> >& rates) {
  vector<int> qpts(3), branches(3);

  for (int i = startIndex; i < endIndex; i++) {
    qpts[0] = _qm.getIbzqpts()[i];
    for (uint p = 0, nprocs = processes[i].size(); p < nprocs; p++) {
      getProcess(processes[i][p], qpts, branches);
      rates[i][branches[0]] += intr_trans_probs[i][p] * getOccupationTerm(occ, processes[i][p][0], qpts, branches);
    }
  }
}

double TCONDCalculator::getOccupationTerm(const vector<vector<double> >& occ, int sign,
                                          const vector<int>& qpts, const vector<int>& branches) {
  if (sign) {  // -
      return (1.0 + occ[qpts[1]][branches[1]] + occ[qpts[2]][branches[2]])/2.0;
  } else {  // +
      return occ[qpts[1]][branches[1]] - occ[qpts[2]][branches[2]];
  }
}

vector<vector<double> > TCONDCalculator::getRates(const vector<vector<double> >& occ) {
  vector<vector<double> > rates = calculateAnharmonicRates(occ);

  if (0 && calc_options.flag("ISOTOPE")) {
    for (int iq = 0; iq < nIQPs; iq++) {
      for (int br = 0; br < nBranches; br++) {
        rates[iq][br] += rates_isotope[iq][br];
      }
    }
  }

  if (calc_options.flag("BOUNDARY")) {
    for (int iq = 0; iq < nIQPs; iq++) {
      for (int br = 0; br < nBranches; br++) {
        rates[iq][br] += rates_boundary[iq][br];
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
                                          const vector<vector<int> >& small_groups,
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
                                      std::ref(small_groups), std::ref(occ), std::ref(mfd), std::ref(delta)));
  }
  for (int icpu = 0; icpu < ncpus; icpu++) {
    threads[icpu]->join();
    delete threads[icpu];
  }
#else
  calculateDelta(0, nIQPs, small_groups, occ, mfd, delta);
#endif

  correctMFD(rates, delta, mfd);
}

void TCONDCalculator::calculateDelta(int startIndex, int endIndex, 
                                     const vector<vector<int> >& small_groups,
                                     const vector<vector<double> >& occ,
                                     const vector<vector<xvector<double> > >& mfd,
                                     vector<vector<xvector<double> > >& delta) {
  xvector<double> correction(3);
  vector<int> qpts(3), branches(3);

  for (int i = startIndex; i < endIndex; i++) {
    qpts[0] = _qm.getIbzqpts()[i];
    for (uint p = 0, nprocs = processes[i].size(); p < nprocs; p++) {
      getProcess(processes[i][p], qpts, branches);
      correction = mfd[qpts[2]][branches[2]];
      if (processes[i][p][0]) correction -= mfd[qpts[1]][branches[1]];
      else correction += mfd[qpts[1]][branches[1]];
      correction *= getOccupationTerm(occ, processes[i][p][0], qpts, branches);
      delta[qpts[0]][branches[0]] += intr_trans_probs[i][p] * correction;
    }

    if (calc_options.flag("ISOTOPE")) {
      int q1, q2, br1, br2;
      q1 = _qm.getIbzqpts()[i];
      for (uint p = 0, nprocs = processes_iso[i].size(); p < nprocs; p++) {
        q2 = processes_iso[i][p][0];
        br1 = processes_iso[i][p][1]/nBranches;
        br2 = processes_iso[i][p][1] % nBranches;
        delta[q1][br1] += intr_trans_probs_iso[i][p] * mfd[q2][br2];
      }
    }

    // Symmetrize
    int symop;
    const vector<_sym_op>& pgroup = _qm.getReciprocalCell().pgroup;
    xmatrix<double> Uc(3, 3);
    uint nsym = small_groups[i].size();
    for (uint isym = 0; isym < nsym; isym++) {
      symop = small_groups[i][isym];
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
