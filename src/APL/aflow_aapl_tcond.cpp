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
  //tmpdir = aurostd::TmpDirectoryCreate("TCOND") + "/";
  tmpdir = "./test/";
  if (aurostd::IsDirectory(tmpdir)) aurostd::RemoveDirectory(tmpdir);
  aurostd::DirectoryMake(tmpdir);

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

  //int max_order = 3;
  //if (calc_options.flag("FOURTH_ORDER")) max_order = 4;
  //for (int order = 3; order <= max_order; order++) calculateScatteringProcesses(order, trans_map);

  calculateTransitionProbabilities();

  vector<vector<int> > small_groups = calculateSmallGroups();
  thermal_conductivity.assign(temperatures.size(), xmatrix<double>(3, 3));
  for (uint t = 0; t < temperatures.size(); t++) {
    thermal_conductivity[t] = calculateThermalConductivityTensor(temperatures[t], small_groups);
  }
  //aurostd::RemoveDirectory(tmpdir);
}

} // namespace apl

namespace apl {

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

/*
//calculateFrequenciesGroupVelocities/////////////////////////////////////////
// Calculates the frequencies and group velocities for each q-point.
void TCONDCalculator::calculateFrequenciesGroupVelocities() {
  _logger << "Calculating frequencies and group velocities." << apl::endl;
  eigenvectors.resize(nQPs, xmatrix<xcomplex<double> >(nBranches, nBranches));
  freq.resize(nQPs, vector<double>(nBranches));
  gvel.resize(nQPs, vector<xvector<double> >(nBranches, xvector<double>(3)));

  vector<vector<int> > irred_qpts;
  vector<vector<int> > symops;
  reduceByFGroup(irred_qpts, symops);
  uint nirred = irred_qpts.size();

  const vector<_sym_op>& fgroup = _pc.getInputCellStructure().fgroup;
  xvector<double> f;
  vector<xmatrix<xcomplex<double> > > dDynMat(3, xmatrix<xcomplex<double> >(nBranches, nBranches));
  for (uint i = 0; i < nirred; i++) {
    int iq = irred_qpts[i][0];
    // Frequencies
    f = _pc.getFrequency(_qm.getQPoint(iq).cpos, apl::THZ | apl::OMEGA, eigenvectors[iq], dDynMat);
    freq[iq] = aurostd::xvector2vector(f);  // Convert to vector to have same indexing as gvel
    // Group velocities
    vector<xvector<double> > gv(nBranches, xvector<double>(3));
    xmatrix<xcomplex<double> > eigen_conj = trasp(eigenvectors[iq]);
    for (int br = 0; br < nBranches; br++) {
      if (freq[iq][br] > _AFLOW_APL_EPS_) {
        for (int j = 1; j < 4; j++) {
          gv[br][j] = real(eigen_conj(br + 1) * (dDynMat[j - 1] * eigenvectors[iq].getcol(br + 1)));
        }
      }
      gvel[iq][br] = gv[br];
    }

    // Transform
    uint nq = irred_qpts[i].size();
    for (uint j = 1; j < nq; j++) {
      int q = irred_qpts[i][j];
      xmatrix<xcomplex<double> > Gamma = calculateGamma(q, symops[i][j]);
      // Frequencies
      freq[q] = freq[iq];
      // Eigenvectors
      eigenvectors[q] = Gamma * eigenvectors[iq];
      // Group velocities
      const xmatrix<double>& Uc = fgroup[symops[i][j]].Uc;
      for (int br = 0; br < nBranches; br++) {
        gvel[q][br] = Uc * gv[br];
      }
    }
  }
}
*/

//reduceByFGroup//////////////////////////////////////////////////////////////
// Reduce the q-point grid by the factor group. This is necessary to get the
// proper symmetry operations for the Gamma matrices (see calculateGamma).
void TCONDCalculator::reduceByFGroup(vector<vector<int> >& irred_qpts,
                                     vector<vector<int> >& symops) {
  irred_qpts.clear();
  symops.clear();
  const vector<_sym_op>& fgroup = _pc.getInputCellStructure().fgroup;
  uint nsym = fgroup.size();

  vector<xmatrix<double> > Uf_rec(nsym);
  for (uint isym = 0; isym < nsym; isym++) Uf_rec[isym] = trasp(inverse(fgroup[isym].Uf));

  vector<int> trans(nsym, -1);
  vector<vector<int> > irred_trans;
  for (int q = 0; q < nQPs; q++) {
    bool append = true;
    for (uint isym = 0; isym < nsym; isym++) {
      for (uint iq = 0; iq < irred_qpts.size(); iq++) {
        if (irred_trans[iq][isym] == q) {
          append = false;
          irred_qpts[iq].push_back(q);
          symops[iq].push_back(isym);
          isym = nsym;
          iq = irred_qpts.size();
        }
      }
    }
    if (append) {
      vector<int> ir(1, q), so(1, 0);
      irred_qpts.push_back(ir);
      symops.push_back(so);
      for (uint isym = 0; isym < nsym; isym++) {
        trans[isym] = _qm.getQPointIndex(Uf_rec[isym] * _qm.getQPoint(q).fpos);
      }
      irred_trans.push_back(trans);
    }
  }
}

//calculateGamma//////////////////////////////////////////////////////////////
// Calculates the Gamma matrix (see Eq. 2.37 in DOI: 10.1103/RevModPhys.40.1).
xmatrix<xcomplex<double> > TCONDCalculator::calculateGamma(int q, int isym) {
  const xstructure& pcell = _pc.getInputCellStructure();
  uint natoms = pcell.atoms.size();
  const _sym_op& fgroup = pcell.fgroup[isym];

  const xvector<double>& qpt = _qm.getQPoint(q).cpos;

  xmatrix<xcomplex<double> > Gamma(nBranches, nBranches);
  xvector<double> dx(3);
  xcomplex<double> phase;
  for (uint at = 0; at < natoms; at++) {
    int at_map = pcell.fgroup[isym].basis_atoms_map[at];
    dx = inverse(fgroup.Uc) * pcell.atoms[at].cpos - pcell.atoms[at_map].cpos;
    phase = exp(iONE * scalar_product(qpt, dx));
    for (int i = 1; i < 4; i++) {
      for (int j = 1; j < 4; j++) {
        Gamma[3 * at_map + i][3 * at + j] = fgroup.Uc[i][j] * phase;
      }
    }
  }
  return Gamma;
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

/*
void TCONDCalculator::getInequivalentQPointSet(int order) {
  const int o1 = order - 1;
  const int o2 = order - 2;

  const xstructure& pcell = _pc.getInputCellStructure();
  uint natoms = pcell.atoms.size();
  const vector<_sym_op>& fgroup = pcell.fgroup;

  int o;
  vector<vector<bool> > signs = getSigns(order);
  uint nsigns = signs.size();

  aurostd::xcombos qpts_combo(nQPs, o2, 'E', true);

  vector<int> set(o1), set_sorted(o1), qpts(o2);
  string path = "irred." + aurostd::utype2string<int>(order) + ".";
  xvector<double> qsum(3);
  int nprocs = nsigns * nIQPs;
  for (int i = startIndex; i < endIndex; i++) {
    int iq = _qm.getIbzqpts()[i];
    string irredfile = path + aurostd::utype2string<int>(i);
    ofstream irredfl;
    openTmpFile(irredfl, irredfile);

    vector<int> small_group;
    uint nsym = small_group.size();

    vector<vector<int> > trans_map(nsym, vector<int>(nQPs));
    for (uint isym = 0; isym < nsym; isym++) {
      for (int q = 0; q < nQPs; q++) {
        trans_map[isym][q]= _qm.getQPointIndex(fgroup[small_group[isym]].Uf * _qm.getQPoint(q).fpos);
      }
    }

    vector<xmatrix<double> > trans_eigen(nsym, xmatrix<double>(nBranches, nBranches));
    for (uint isym = 0; isym < nsym; isym++) {
      for (uint at = 0; at < natoms; at++) {
        int at_map = fgroup[isym].basis_atoms_map[at];
        for (int r = 1; r < 4; r++) {
          for (int c = 1; c < 4; c++) {
            trans_eigen[isym][3 * at_map + r][3 * at + c] = fgroup[small_group[isym]].Uc[r][c];
          }
        }
      }
    }

    for (uint s = 0; s < nsigns; s++) {
      vector<int> ineq_sets_index, weights;
      vector<vector<int> > ineq_sets;
      qpts_combo.reset();
      int q = 0;
      int nsets = 0;
      int iset;
      uint isym;
      while (qpts_combo.increment()) {
        qpts = qpts_combo.getCombo();
        qsum = _qm.getQPoint(iq).fpos;
        for (o = 0; o < o2; o++) {
          if (signs[s][o]) qsum -= _qm.getQPoint(qpts[o]).fpos;
          else qsum += _qm.getQPoint(qpts[o]).fpos;
          set[o] = qpts[o];
        }
        set[o2] = _qm.getQPointIndex(qsum);
        for (isym = 0; isym < nsym; isym++) {
          for (o = 0; o < o2; o++) {
            if (signs[s][o]) set_sorted[o] = -trans_map[small_groups[isym]][set[o]];
            else set_sorted[o] = trans_map[small_groups[i][isym]][set[o]];
          }
          set_sorted[o2] = -trans_map[small_groups[isym]][set[o2]];
          std::sort(set_sorted.begin(), set_sorted.end());
          for (iset = 0; iset < nsets; iset++) {
            for (o = 0; o < o1; o++) {
              if (set_sorted[o] != ineq_sets[iset][o]) break;
            }
            if (o == o1) break;
          }
          if (iset != nsets) {
            if (isym == 0) break;  // Identity means it is a permutation
            else {
              for (o = 0; o < o1; o++) {
                if (eigenvectors[std::abs(set_sorted[o])] != trans_eigen[isym] * eigenvectors[std::abs(ineq_sets[iset][o])]) break;
              }
              if (o != o1) break;
            }
          }
        }
        if (isym == nsym) {
          ineq_sets_index.push_back(q);
          for (o = 0; o < o2; o++) {
            if (signs[s][o]) set[o] = -set[o];
          }
          set[o2] = -set[o2];
          std::sort(set.begin(), set.end());
          ineq_sets.push_back(set);
          nsets++;
          weights.push_back(1);
        } else {
          weights[iset]++;
        }
        q++;
      }
      std::cout << "nsets = " << nsets << std::endl;
      exit(0);
      irredfl.write((char*)(&nsets), sizeof(int));
      for (uint q = 0; q < ineq_sets_index.size(); q++) {
        char w = (char) weights[q];
        irredfl.write((char*)(&ineq_sets_index[q]), sizeof(int));
        irredfl.write((char*)(&w), sizeof(char));
      }
      _logger.updateProgressBar(1.0/nprocs);
    }
    closeTmpFile(irredfl, irredfile);
  }
}
*/

vector<vector<bool> > TCONDCalculator::getSigns(int order) {
  const int o2 = order - 2;
  vector<vector<bool> > signs;
  aurostd::xcombos sign_combos(2, o2, 'C', true);
  while (sign_combos.increment()) {
    vector<bool> sgn(o2);
    for (int o = 0; o < o2; o++) {
      if (sign_combos.getCombo()[o] == 0) sgn[o] = false;
      else sgn[o] = true;
    }
    signs.push_back(sgn);
  }
  return signs;
}

}

/********************************** WEIGHTS *********************************/

namespace apl {

double TCONDCalculator::calculateIntegrationWeights(int order, const LTMethod& _lt) {
  _logger << "Calculating integration weights" << apl::endl;
  string message = "Calculating weights";
  double phase_space = 0.0;
  vector<double> ps(nIQPs);
  _logger.initProgressBar(message);
  calculateWeightsLT(_lt, order, 0, nIQPs, ps);
  _logger.finishProgressBar();
/*
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
*/
  for (int i = 0; i < nIQPs; i++) phase_space += ps[i];
  return phase_space;
}

void TCONDCalculator::calculateWeightsLT(const LTMethod& _lt, int order,
                                         int startIndex, int endIndex, vector<double>& ps) {
  vector<vector<int> > branches;
  vector<double> frequencies(nQPs), weights(nQPs);
  xvector<double> qsum;
  vector<int> lastq(nQPs);
  double freq_ref;
  vector<bool> zerofreq(nQPs);
  uint b, s;
  int iq, qpt, q, o;
  bool sign;

  const int o2 = order - 2;
  const int o3 = order - 3;
  const double prefactor = ps_prefactors[o3]/(std::pow(nBranches, order) * std::pow(nQPs, o2));

  string pathp, pathw, filenamep, filenamew;
  pathp = "proc.unfiltered." + aurostd::utype2string<int>(order) + ".";
  pathw = "wght." + aurostd::utype2string<int>(order) + ".";

  aurostd::xcombos qpt_combos(nQPs, o3, 'E', true);
  vector<vector<bool> > signs = getSigns(order);
  uint nsigns = signs.size();

  aurostd::xcombos branch_combos(nBranches, order, 'E', true);
  while (branch_combos.increment()) branches.push_back(branch_combos.getCombo());
  uint nbr = branches.size();

  int nprocs = nsigns * nIQPs;

  for (int i = startIndex; i < endIndex; i++) {
    ofstream outp, outw;
    iq = _qm.getIbzqpts()[i];
    filenamew = pathw + aurostd::utype2string<int>(i);
    filenamep = pathp + aurostd::utype2string<int>(i);
    openTmpFile(outp, filenamep);
    openTmpFile(outw, filenamew);
    for (s = 0; s < nsigns; s++) {
      int c = 0;
      qpt_combos.reset();
      bool increment = ((order == 3) || qpt_combos.increment());
      while (increment) {
        const vector<int>& qpoints = qpt_combos.getCombo();
        qsum = _qm.getQPoint(iq).fpos;
        for (o = 0; o < o3; o++) {
          if (signs[s][o]) qsum -= _qm.getQPoint(qpoints[o]).fpos;
          else qsum += _qm.getQPoint(qpoints[o]).fpos;
        }
        for (q = 0; q < nQPs; q++) {
          if (signs[s][o3]) lastq[q] = _qm.getQPointIndex(qsum - _qm.getQPoint(q).fpos);
          else lastq[q] = _qm.getQPointIndex(qsum + _qm.getQPoint(q).fpos);
        }
        for (b = 0; b < nbr; b++) {
          freq_ref = freq[iq][branches[b][0]];
          if (1 || freq_ref > _AFLOW_APL_EPS_) {
            for (o = 0; o < o3; o++) {
              if (freq[qpoints[o]][branches[b][o+1]] < _AFLOW_APL_EPS_) break;
              if (signs[s][o]) freq_ref -= freq[qpoints[o]][branches[b][o + 1]];
              else freq_ref += freq[qpoints[o]][branches[b][o + 1]];
            }
            if (o == o3) {
              for (q = 0; q < nQPs; q++) {
                zerofreq[q] = ((freq[q][branches[b][o2]] < _AFLOW_APL_EPS_) || (freq[lastq[q]][branches[b].back()] < _AFLOW_APL_EPS_));
                if (signs[s][o3]) frequencies[q] = freq[q][branches[b][o2]];
                else frequencies[q] = -freq[q][branches[b][o2]];
                frequencies[q] += freq[lastq[q]][branches[b].back()];
              }
              weights = getWeightsLT(_lt, freq_ref, frequencies);
              for (q = 0; q < nQPs; q++) {
                if (!zerofreq[q] && (weights[q] > _ZERO_TOL_)) {
                  outw.write((char*)(&weights[q]), sizeof(double));
                  for (o = 0; o < o2; o++) {
                    sign = signs[s][o];
                    outp.write((char*)(&sign), sizeof(bool));
                  }
                  qpt = c + q;
                  outp.write((char*)(&qpt), sizeof(int));
                  outp.write((char*)(&b), sizeof(int));
                  ps[i] += scatt_multi[o3][s] * weights[q];
                }
              }
            }
          }
        }
        increment = qpt_combos.increment();
        c += nQPs;
      }
      _logger.updateProgressBar(1.0/nprocs);
    }
    ps[i] *= _qm.getWeights()[i] * prefactor;
    closeTmpFile(outp, filenamep);
    closeTmpFile(outw, filenamew);
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
      for (i = 0; i < 4; i++) w[i] = 0.0;
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

  int max_order;
  if (calc_options.flag("FOURTH_ORDER")) max_order = 4;
  else max_order = 3;
  //max_order = 4;
  vector<double> phase_space(max_order - 2);
  for (int order = 3; order <= max_order; order++) {
    /*
    message = "Scattering Proccesses";
    _logger << "Calculating scattering processes" << apl::endl;
    getInequivalentQPointSet(order, 0, nIQPs, small_groups, trans_map);
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    thread_dist = setupMPI(message, _logger, nIQPs, ncpus);
    threads.clear();
    for (int icpu = 0; icpu < ncpus; icpu++) {
      threads.push_back(new std::thread(&TCONDCalculator::getInequivalentQPointSet, this,
                                        order, thread_dist[icpu][0], thread_dist[icpu][1],
                                        std::ref(small_groups), std::ref(trans_map)));
    }
    finishMPI(threads, _logger);
#else
    _logger.initProgressBar(message);
    getInequivalentQPointSet(order, 0, nIQPs, small_groups, trans_map);
    _logger.finishProgressBar();
#endif
    for (int i = 0; i < nIQPs; i++) {
      ifstream irredfl;
      string irredfile = "irred.3." + aurostd::utype2string<int>(i);
      openTmpFile(irredfl, irredfile);
      irredfl.seekg(0, irredfl.end);
      int n = irredfl.tellg()/sizeof(int);
      std::cout << i << " " << n << std::endl;
      closeTmpFile(irredfl, irredfile);
    }
    exit(0);
    */
    phase_space[order - 3] = calculateIntegrationWeights(order, _lt);
    _logger << "phase space = " << phase_space[order - 3] << apl::endl;
  }
  exit(0);

  for (int order = 3; order <= max_order; order++) {
    message = "Transition Probabilities";
    _logger << "Calculating transition probabilities for " << order << "-phonon processes." << apl::endl;
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
    ifstream transfl, procfl, weightfl;
    string transfile = "trans.3.0";
    string procfile = "proc.unfiltered.3.0";
    string weightfile = "wght.3.0";
    openTmpFile(transfl, transfile);
    openTmpFile(procfl, procfile);
    openTmpFile(weightfl, weightfile);
    vector<double> trans, weights;
    vector<vector<int> > qpts, branches;
    vector<int> qp(3), brnch(3);
    qp[0] = _qm.getIbzqpts()[0];
    vector<bool> sg(1);
    vector<vector<bool> > signs;
    double t;
    while (transfl.peek() != EOF) {
      transfl.read((char*)(&t), sizeof(double));
      trans.push_back(t);
      weightfl.read((char*)(&t), sizeof(double));
      weights.push_back(t);
      getProcess(3, procfl, sg, qp, brnch);
      signs.push_back(sg);
      qpts.push_back(qp);
      branches.push_back(brnch);
    }
    closeTmpFile(transfl, transfile);
    closeTmpFile(procfl, procfile);
    closeTmpFile(weightfl, weightfile);
    for (uint i = 0; i < trans.size(); i++) {
      std::cout << signs[i][0] << ", " << aurostd::joinWDelimiter(qpts[i], " ") << ", " << aurostd::joinWDelimiter(branches[i], " ") << ", " << weights[i] << " " << (trans[i]/weights[i]) << " " << trans[i] << std::endl;
      for (int q = 1; q < 3; q++) {
        std::cout << _qm.getQPoint(qpts[i][q]).fpos << std::endl;
      }
      if (trans[i] > _ZERO_TOL_ ) {
      for (uint j = 0; j < trans.size(); j++) {
//        if ((std::abs(trans[i] - trans[j]) < _ZERO_TOL_) && (std::abs((trans[i]/weights[i]) - (trans[j]/weights[j])) < _ZERO_TOL_)) {
        if ((i != j) && (std::abs(trans[i] - trans[j]) < _ZERO_TOL_)) {
//          for (int q = 1; q < 3; q++) {
//            std::cout << _qm.getQPoint(qpts[i][q]).cpos << std::endl;
//            std::cout << _qm.getQPoint(qpts[i][q]).fpos << std::endl;
//            std::cout << eigenvectors[qpts[i][q]] << std::endl << std::endl;
//          }
          std::cout << signs[j][0] << ", " << aurostd::joinWDelimiter(qpts[j], " ") << ", " << aurostd::joinWDelimiter(branches[j], " ") << ", " << weights[j] << " " << (trans[j]/weights[j]) << " " << trans[j] << std::endl;
          for (int q = 1; q < 3; q++) {
//            std::cout << _qm.getQPoint(qpts[j][q]).cpos << std::endl;
            std::cout << _qm.getQPoint(qpts[j][q]).fpos << std::endl;
//            std::cout << eigenvectors[qpts[j][q]] << std::endl << std::endl;
          }
        }
      }
      }
      std::cout << "--------------------------------------------------" << std::endl << std::endl;
    }
    exit(0);
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
                                                             const vector<vector<vector<vector<xcomplex<double> > > > >& phases) {
  const int o3 = order - 3;
  const int o2 = order - 2;
  const int o1 = order - 1;
  bool sign;

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
  int iat, o, e, q, br;
  uint c, crt;
  double transprob, weight;

  string transfile, procfile, weightfile, procunffile;
  string pathpunf = "proc.unfiltered." + aurostd::utype2string<int>(order) + ".";
  string pathp = "proc." + aurostd::utype2string<int>(order) + ".";
  string patht = "trans." + aurostd::utype2string<int>(order) + ".";
  string pathw = "wght." + aurostd::utype2string<int>(order) + ".";

  for (int i = startIndex; i < endIndex; i++) {
    qpts[0] = _qm.getIbzqpts()[i];
    ifstream procunffl, weightfl;
    ofstream transfl, procfl;
    procunffile = pathpunf + aurostd::utype2string<int>(i);
    procfile = pathp + aurostd::utype2string<int>(i);
    transfile = patht + aurostd::utype2string<int>(i);
    weightfile = pathw + aurostd::utype2string<int>(i);
    openTmpFile(procunffl, procunffile);
//    openTmpFile(procfl, procfile);
    openTmpFile(transfl, transfile);
    openTmpFile(weightfl, weightfile);

    while (!weightfl.eof()) {
      weightfl.read((char*)(&weight), sizeof(double));
      if (weightfl.eof()) break;
      getProcess(order, procunffl, signs, qpts, branches);
  
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
      vector<xcomplex<double> > m(nclusters, xcomplex<double>(0.0, 0.0));
      for (c = 0; c < nclusters; c++) {
        const vector<int>& atoms = clusters[c].atoms;
        iat = scell.sc2pcMap(atoms[0]);
        prefactor = invmasses[c] * phases[iat][atoms[o1]][qpts[o1]][1];
        for (o = 1; o < o1; o++) prefactor *= phases[iat][atoms[o]][qpts[o]][(int) signs[o - 1]];
        e = 0;
        for (o = 0; o < order; o++) e += scell.sc2pcMap(atoms[o]) * atpowers[o];
        for (crt = 0; crt < ncart; crt++) {
          // Perform multiplication expliclty in place instead of using xcomplex.
          // This is three times as fast because constructors and destructors are not called.
          matrix.re += ifcs[c][crt] * (prefactor.re * eigenprods[e][crt].re - prefactor.im * eigenprods[e][crt].im);
          matrix.im += ifcs[c][crt] * (prefactor.re * eigenprods[e][crt].im + prefactor.im * eigenprods[e][crt].re);
        }
      }
//      transprob = probability_prefactor * magsqr(matrix) * weight;
//      for (o = 0; o < order; o++) transprob /= freq[qpts[o]][branches[o]];
      transprob = magsqr(matrix);
      /*
      if (1 || transprob > _ZERO_TOL_) {
        q = 0;
        for (o = 0; o < o2; o++) q += qpts[o2 - o] * (int) std::pow(nQPs, o);
        br = 0;
        for (o = 0; o < order; o++) br += branches[o1 - o] * (int) std::pow(nBranches, o);
        for (o = 0; o < o2; o++) {
          sign = signs[o];
          procfl.write((char*)(&sign), sizeof(bool));
        }
        procfl.write((char*)(&q), sizeof(int));
        procfl.write((char*)(&br), sizeof(int));
        */
        transfl.write((char*)(&transprob), sizeof(double));
      //}
    }
//    weightfl.close();
//    aurostd::RemoveFile(tmpdir + weightfile);
//    aurostd::RemoveFile(tmpdir + procunffile);
    closeTmpFile(weightfl, weightfile);
    closeTmpFile(procfl, procfile);
    closeTmpFile(procunffl, procfile);
    closeTmpFile(transfl, transfile);
    _logger.updateProgressBar(1.0/nIQPs);
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

void TCONDCalculator::getProcess(int order, ifstream& procfl, vector<bool>& signs,
                                 vector<int>& qpts, vector<int>& branches) {
  const int o1 = order - 1;
  const int o2 = order - 2;
  const int o3 = order - 3;
  int o, r, pw;
  bool sign;

  for (o = 0; o < o2; o++) {
    procfl.read((char*)(&sign), sizeof(bool));
    signs[o] = sign;
  }

  procfl.read((char*)(&r), sizeof(int));
  xvector<double> qsum = _qm.getQPoint(qpts[0]).fpos;
  for (o = 0; o < o3; o++) {
    pw = (int) std::pow(nQPs, o3 - o);
    qpts[o + 1] = r/pw;
    r = r % pw;
    if (!signs[o]) qsum += _qm.getQPoint(qpts[o + 1]).fpos;
    else qsum -= _qm.getQPoint(qpts[o + 1]).fpos;
  }
  qpts[o2] = r;
  if (signs[o3]) qpts[o1] = _qm.getQPointIndex(qsum - _qm.getQPoint(qpts[o2]).fpos);
  else qpts[o1] = _qm.getQPointIndex(qsum + _qm.getQPoint(qpts[o2]).fpos);
  
  procfl.read((char*)(&r), sizeof(int));
  for (o = 0; o < o1; o++) {
    pw = (int) std::pow(nBranches, o1 - o);
    branches[o] = r/pw;
    r = r % pw;
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
    transfile = "trans.iso." + aurostd::utype2string<int>(iq);
    procfile = "proc.iso." + aurostd::utype2string<int>(iq);
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
  string transfile = "trans.boundary";
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
                                                                    const vector<vector<int> >& small_groups) {
  _logger << "Calculating thermal conductivity for " << T << " K." << apl::endl;
  vector<vector<double> > occ = getOccupationNumbers(T);
  _logger << "Calculating scattering rates." << apl::endl;
  int max_order;
  if (calc_options.flag("FOURTH_ORDER")) max_order = 4;
  else max_order = 3;
  for (int o = 3; o <= max_order; o++) calculateAnharmonicRates(o, T, occ);

  _logger << "Calculating RTA" << apl::endl;
  vector<vector<double> > rates = getRates(T);
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
  vector<int> qpts(order), branches(order);
  vector<bool> signs(order - 2);
  string procfile, transfile, qptfile;
  double trans;
  string pathp = "proc." + aurostd::utype2string<int>(order) + ".";
  string patht = "trans." + aurostd::utype2string<int>(order) + ".";

  for (int i = startIndex; i < endIndex; i++) {
    ifstream procfl, transfl;
    procfile = pathp + aurostd::utype2string<int>(i);
    transfile = patht + aurostd::utype2string<int>(i);
    openTmpFile(procfl, procfile);
    openTmpFile(transfl, transfile);
    qpts[0] = _qm.getIbzqpts()[i];
    while (transfl.peek() != EOF) {
      transfl.read((char*)(&trans), sizeof(double));
      getProcess(order, procfl, signs, qpts, branches);
      rates[i][branches[0]] += trans * getOccupationTerm(order, occ, signs, qpts, branches);
    }
    closeTmpFile(procfl, procfile);
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
      ifstream procfl, transfl;
      procfile = "proc." + aurostd::utype2string<int>(o) + "." + aurostd::utype2string<int>(i);
      transfile = "trans." + aurostd::utype2string<int>(o) + "." + aurostd::utype2string<int>(i);
      openTmpFile(procfl, procfile);
      openTmpFile(transfl, transfile);
      while (transfl.peek() != EOF) {
        transfl.read((char*)(&trans), sizeof(double));
        getProcess(o, procfl, signs, qpts, branches);
        correction = mfd[qpts.back()][branches.back()];
        for (int j = 0; j < o2; j++) {
          if (signs[j]) correction -= mfd[qpts[j + 1]][branches[j + 1]];
          else correction += mfd[qpts[j + 1]][branches[j + 1]];
        }
        delta[qpts[0]][branches[0]] += trans * correction * getOccupationTerm(o, occ, signs, qpts, branches);
      }
      closeTmpFile(procfl, procfile);
      closeTmpFile(transfl, transfile);
    }

    if (calc_options.flag("ISOTOPE")) {
      ifstream procfl_iso, transfl_iso;
      procfile = "proc.iso." + aurostd::utype2string<int>(i);
      transfile = "trans.iso." + aurostd::utype2string<int>(i);

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

void TCONDCalculator::openTmpFile(ofstream& file, const string& filename) {
  file.open((tmpdir + filename).c_str(), std::ios::out|std::ios::binary);
}

void TCONDCalculator::openTmpFile(ifstream& file, const string& filename) {
  aurostd::UncompressFile(tmpdir + filename + "." + calc_options.getattachedscheme("KZIP_BIN"));
  file.open((tmpdir + filename).c_str(), std::ios::in|std::ios::binary);
}

void TCONDCalculator::closeTmpFile(ofstream& file, const string& filename) {
  file.close();
  aurostd::CompressFile(tmpdir + filename, calc_options.getattachedscheme("KZIP_BIN"));
}

void TCONDCalculator::closeTmpFile(ifstream& file, const string& filename) {
  file.close();
  aurostd::CompressFile(tmpdir + filename, calc_options.getattachedscheme("KZIP_BIN"));
}

void TCONDCalculator::writeRatesToTmpFile(const string& ratefile, 
                                          const vector<vector<double> >& rates) {
  ofstream ratefl;
  openTmpFile(ratefl, ratefile);
  for (int q = 0; q < nIQPs; q++) {
    for (int br = 0; br < nBranches; br++) {
      ratefl.write((char*)(&rates[q][br]), sizeof(double));
    }
  }
  closeTmpFile(ratefl, ratefile);
}

vector<vector<double> > TCONDCalculator::readRatesFromTmpFile(const string& ratefile) {
  ifstream ratefl;
  vector<vector<double> > rates(nIQPs, vector<double>(nBranches, 0.0));
  openTmpFile(ratefl, ratefile);
  for (int q = 0; q < nIQPs; q++) {
    for (int br = 0; br < nBranches; br++) {
      ratefl.read((char*)(&rates[q][br]), sizeof(double));
    }
  }
  closeTmpFile(ratefl, ratefile);
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
