//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                  Simon Divilov - Duke University 2022                   *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov, 2022
// simon.divilov@duke.edu
// 
#ifndef _AFLOW_QCA_CPP_
#define _AFLOW_QCA_CPP_
#include "aflow.h"
#include "aflow_qca.h"
#include "aflow_chull.h"
#include "aflow_pocc.h"

// ###############################################################################
//            AFLOW Quasi-Chemical Approximation (QCA) (2022-)
// ###############################################################################

// **************************************************************************
// Class _qca data
// **************************************************************************
// Constructor
_qca_data::_qca_data() {
  // Input data
  num_threads = 0;
  min_sleep = 0;
  format_data = "";
  format_image = "";
  screen_only = false;
  image_only = false;
  calc_binodal = false;
  calc_spinodal = false;
  workdirpath = "";
  rootdirpath = "";
  plattice = "";
  elements.clear();
  aflow_max_num_atoms = 0;
  max_num_atoms = 0;
  cv_cut = 0.0;
  conc_npts = 0;
  conc_curve = false;
  conc_curve_range.clear();
  conc_macro.clear();
  temp_npts = 0;
  temp_range.clear();
  temp.clear();
  // Derived data
  alloyname = "";
  rundirpath = "";
  vstr_aflow.clear();
  lat_atat = "";
  vstr_atat.clear();
  mapstr.clear();
  // Cluster data
  cv_cluster = 0.0;
  num_atom_cluster.clear();
  degeneracy_cluster.clear();
  conc_cluster.clear();
  excess_energy_cluster.clear();
  // Thermo data
  prob_ideal_cluster.clear();
  prob_cluster.clear();
  rel_s_ec = 0.0;
  temp_ec = 0.0;
  rel_s.clear();
  binodal_boundary.clear();
}

// Destructor
_qca_data::~_qca_data() {
  free();
}
void _qca_data::free() {
}

// Copy constructor
const _qca_data& _qca_data::operator=(const _qca_data &b) {
  if (this != &b) {
    // Input data
    num_threads = b.num_threads;
    min_sleep = b.min_sleep;
    format_data = b.format_data;
    format_image = b.format_image;
    screen_only = b.screen_only;
    image_only = b.image_only;
    calc_binodal = b.calc_binodal;
    calc_spinodal = b.calc_spinodal;
    workdirpath = b.workdirpath;
    rootdirpath = b.rootdirpath;
    plattice = b.plattice;
    elements = b.elements;
    max_num_atoms = b.max_num_atoms;
    cv_cut = b.cv_cut;
    conc_npts = b.conc_npts;
    conc_curve = b.conc_curve;
    conc_curve_range = b.conc_curve_range;
    conc_macro = b.conc_macro;
    temp_npts = b.temp_npts;
    temp_range = b.temp_range;
    temp = b.temp;
    // Derived data
    alloyname = b.alloyname;
    rundirpath = b.rundirpath;
    vstr_aflow = b.vstr_aflow;
    lat_atat = b.lat_atat;
    vstr_atat = b.vstr_atat;
    mapstr = b.mapstr;
    // Cluster data
    cv_cluster = b.cv_cluster;
    num_atom_cluster = b.num_atom_cluster;
    degeneracy_cluster = b.degeneracy_cluster;
    conc_cluster = b.conc_cluster;
    excess_energy_cluster = b.excess_energy_cluster;
    // Thermo data
    prob_ideal_cluster = b.prob_ideal_cluster;
    prob_cluster = b.prob_cluster;
    rel_s_ec = b.rel_s_ec;
    temp_ec = b.temp_ec;
    rel_s = b.rel_s;
    binodal_boundary = b.binodal_boundary;
  }
  return *this;
}

// ***************************************************************************
// qca::quasiChemicalApprox
// ***************************************************************************
namespace qca {
 void quasiChemicalApprox(const aurostd::xoption& vpflow) {
    if (vpflow.flag("QCA::USAGE")) {
      if (!vpflow.flag("QCA::SCREEN_ONLY")) {displayUsage();}
      return;
    }
    _qca_data qca_data;
    qca_data.min_sleep = DEFAULT_QCA_MIN_SLEEP_SECONDS;
    qca_data.aflow_max_num_atoms = DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS;
    if (!vpflow.getattachedscheme("QCA::DIRECTORY").empty()) {
      qca_data.rootdirpath = vpflow.getattachedscheme("QCA::DIRECTORY");
    }
    else {
      qca_data.rootdirpath = aurostd::getPWD();
    }
    qca_data.plattice = vpflow.getattachedscheme("QCA::PLATTICE");
    if (!vpflow.getattachedscheme("QCA::ELEMENTS").empty()) {
      aurostd::string2tokens(vpflow.getattachedscheme("QCA::ELEMENTS"), qca_data.elements, ",");
    }
    if (!vpflow.getattachedscheme("QCA::MAX_NUM_ATOMS").empty()) {
      qca_data.max_num_atoms = aurostd::string2utype<int>(vpflow.getattachedscheme("QCA::MAX_NUM_ATOMS"));
    }
    else {
      qca_data.max_num_atoms = DEFAULT_QCA_MAX_NUM_ATOMS;
    }
    if (!vpflow.getattachedscheme("QCA::CV_CUTOFF").empty()) {
      qca_data.cv_cut = aurostd::string2utype<double>(vpflow.getattachedscheme("QCA::CV_CUTOFF"));
    }
    else {
      qca_data.cv_cut = DEFAULT_QCA_CV_CUTOFF;
    }
    if (!vpflow.getattachedscheme("QCA::CONC_CURVE_RANGE").empty()) {
      qca_data.conc_curve = true;
      aurostd::string2tokens(vpflow.getattachedscheme("QCA::CONC_CURVE_RANGE"), qca_data.conc_curve_range, ",");
    }
    if (!vpflow.getattachedscheme("QCA::CONC_NPTS").empty()) {
      qca_data.conc_npts = aurostd::string2utype<int>(vpflow.getattachedscheme("QCA::CONC_NPTS"));
    }
    else {
      qca_data.conc_npts = DEFAULT_QCA_CONC_NPTS;
    }
    if (!vpflow.getattachedscheme("QCA::TEMP_RANGE").empty()) {
      aurostd::string2tokens(vpflow.getattachedscheme("QCA::TEMP_RANGE"), qca_data.temp_range, ",");
    }
    else {
      qca_data.temp_range = {DEFAULT_QCA_TEMP_MIN, DEFAULT_QCA_TEMP_MAX};
    }
    if (!vpflow.getattachedscheme("QCA::TEMP_NPTS").empty()) {
      qca_data.temp_npts = aurostd::string2utype<int>(vpflow.getattachedscheme("QCA::TEMP_NPTS"));
    }
    else {
      qca_data.temp_npts = DEFAULT_QCA_TEMP_NPTS;
    }
    if (!vpflow.getattachedscheme("QCA::FORMAT_DATA").empty()) {
      qca_data.format_data = vpflow.getattachedscheme("QCA::FORMAT_DATA");
    }
    else {
      qca_data.format_data = DEFAULT_QCA_FORMAT_DATA;
    }
    if (!vpflow.flag("QCA::NO_PLOT")) {
      if (!vpflow.getattachedscheme("QCA::FORMAT_IMAGE").empty()) {
        qca_data.format_image = vpflow.getattachedscheme("QCA::FORMAT_IMAGE");
      }
      else {
        qca_data.format_image = DEFAULT_QCA_FORMAT_PLOT;
      }
    }
    if (vpflow.flag("QCA::SCREEN_ONLY")) {
      qca_data.screen_only = true;
      qca_data.format_data = DEFAULT_QCA_FORMAT_DATA;
    }
    if (vpflow.flag("QCA::IMAGE_ONLY")) {qca_data.image_only = true;}
    if (vpflow.flag("QCA::BINODAL")) {qca_data.calc_binodal = true;}
    if (vpflow.flag("QCA::SPINODAL")) {qca_data.calc_spinodal = true;}
    quasiChemicalApprox(qca_data);
    return;
 }

  void quasiChemicalApprox(_qca_data& qca_data) {
    // Clean-up input data and check for errors
    if (XHOST.vflag_control.flag("XPLUG_NUM_THREADS") && !(XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX"))) {
      qca_data.num_threads = aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS"));
    }
    qca_data.workdirpath = aurostd::getPWD();
    qca_data.rootdirpath = aurostd::CleanFileName(qca_data.rootdirpath);
    qca_data.plattice = aurostd::tolower(qca_data.plattice);
    qca_data.format_data = aurostd::tolower(qca_data.format_data);
    qca_data.format_image = aurostd::tolower(qca_data.format_image);
    aurostd::sort_remove_duplicates(qca_data.elements);
    errorChecks(qca_data);
    qca_data.rundirpath += qca_data.rootdirpath + "/" + pflow::arity_string(qca_data.elements.size(), false, false) + "/" + qca_data.plattice + "/" + qca_data.alloyname;
    aurostd::DirectoryMake(qca_data.rundirpath);
    // Only plot data from JSON file. This is useful when the plotting routine cannot be ran on the machine
    // running the calculation.
    if (qca_data.image_only) {
      readData(qca_data);
      plotData(qca_data);
      return;
    }
    // Calculations
    calcBinodalData(qca_data);
    // Write and plot results
    if (qca_data.calc_binodal) {
      writeData(qca_data);
      plotData(qca_data);
    }
    return;
  }
}

// ***************************************************************************
// qca::errorChecks
// ***************************************************************************
namespace qca {
  void errorChecks(_qca_data& qca_data) {
    // Check if number of threads is valid
    if (qca_data.num_threads < 1) {qca_data.num_threads = 1;}
    // Check if min sleep is at least 1 sec
    if (qca_data.min_sleep < 1) {qca_data.min_sleep = 1;}
    // Check if max number of atoms in cluster is valid
    if (qca_data.aflow_max_num_atoms < 1 || qca_data.max_num_atoms < 1) {
      string message = "Maximum number of atoms per cluster must be at least 1";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // Check if directory is writable
    if (!aurostd::DirectoryMake(qca_data.rootdirpath)) {
      string message = "Cannot create directory";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
    // Check if parent lattice is valid
    if (qca_data.plattice != "fcc" && qca_data.plattice != "bcc" && qca_data.plattice != "hcp") {
      string message = "Parent lattice \"" + qca_data.plattice + "\" is invalid";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    // Check if alloy is at least binary
    if (qca_data.elements.size() < 2) {
      string message = "Alloy must be at least binary";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // Check if HCP parent lattice is used for greater than binary
    if (qca_data.plattice == "hcp" && qca_data.elements.size() != 2) {
      string message = "HCP parent lattice of alloys greater than binary not supported";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    // Check if elements are valid, construct alloy name
    for (uint i = 0; i < qca_data.elements.size(); i++) {
      if (!xelement::xelement::isElement(qca_data.elements[i])) {
        string message = "Element \"" + qca_data.elements[i] + "\" is invalid";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      qca_data.alloyname += qca_data.elements[i];
    }
    // Check if concentration range has enough points
    if (qca_data.conc_npts < 2) {
      string message = "Number of points for the concentration range must be at least 2";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // Check if concentration curve format is valid
    if (qca_data.conc_curve) {
      if (qca_data.conc_curve_range.size() != 2 * qca_data.elements.size()) {
        string message = "Concentration curve must have format [X1_start, X1_end, X2_start, X2_end,...X(K)_start, X(K)_end]";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
      }
      // Check if concentration curve is within [0,1] and sums to 1
      double totconc_init = 0.0, totconc_final = 0.0;
      vector<double> vconc_init, vconc_final;
      vector<int> vzeros_init, vzeros_final;
      for (uint i = 0; i < qca_data.elements.size() && totconc_init <= 1.0 && totconc_final <= 1.0; i++) {
        if (qca_data.conc_curve_range[2 * i] < 0 || qca_data.conc_curve_range[2 * i] > 1 ||
            qca_data.conc_curve_range[2 * i + 1] < 0 || qca_data.conc_curve_range[2 * i + 1] > 1) {
          string message = "Concentration curve must be defined on [0,1]";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
        }
        totconc_init += qca_data.conc_curve_range[2 * i]; totconc_final += qca_data.conc_curve_range[2 * i + 1];
        vconc_init.push_back(qca_data.conc_curve_range[2 * i]); vconc_final.push_back(qca_data.conc_curve_range[2 * i + 1]);
      }
      if (!aurostd::isequal(totconc_init, 1.0) || !aurostd::isequal(totconc_final, 1.0)) {
        string message = "Total concentration must sum to 1";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
      }
      // Redefine concentration curve from [0,1] to (0,1)
      aurostd::WithinList(vconc_init, 0.0, vzeros_init); aurostd::WithinList(vconc_final, 0.0, vzeros_final);
      double cdelta_init = CONC_SHIFT * (double)vzeros_init.size() / ((double)vconc_init.size() - (double)vzeros_init.size());
      double cdelta_final = CONC_SHIFT * (double)vzeros_final.size() / ((double)vconc_final.size() - (double)vzeros_final.size());
      for (uint i = 0; i < qca_data.conc_curve_range.size(); i++) {
        if (aurostd::isequal(qca_data.conc_curve_range[i], 0.0)) {
          qca_data.conc_curve_range[i] += CONC_SHIFT;
        }
        else {
          qca_data.conc_curve_range[i] -= (i % 2) ? cdelta_final : cdelta_init;
        }
      }
    }
    // Check if temperature range format is valid
    if (qca_data.temp_range.size() != 2) {
      string message = "Temperature range must have format [T_start T_end]";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    // Check if temperature range has enough points
    if (qca_data.temp_npts < 2) {
      string message = "Number of points for the temperature range must be at least 2";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // Check if temperature values are valid
    if (qca_data.temp_range[0] < 0 || qca_data.temp_range[1] < 0) {
      string message = "Temperature cannot be below 0K";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    // Check if output format is valid
    if (qca_data.format_data != "txt" && qca_data.format_data != "json") {
      string message = "Format \"" + qca_data.format_data + "\" is invalid";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }
    if (qca_data.format_image != "" && qca_data.format_image != "pdf" && qca_data.format_image != "eps" && qca_data.format_image != "png") {
      string message = "Format \"" + qca_data.format_image + "\" is invalid";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }
  }
}

// ***************************************************************************
// qca::calcSpinodalData
// ***************************************************************************
namespace qca {
  void calcSpinodalData(_qca_data& qca_data) {
    cerr << qca_data.rundirpath << endl;
  }
}

// ***************************************************************************
// qca::calcBinodalData
// ***************************************************************************
// Binodal construction based on the method developed in Y. Lederer et al., Acta Materialia, 159 (2018)
namespace qca {
  void calcBinodalData(_qca_data& qca_data) {
    if (aurostd::FileExist(qca_data.rundirpath + "/fit.out") && aurostd::FileExist(qca_data.rundirpath + "/predstr.out")) { // read ATAT data
      qca_data.vstr_atat = getATATXstructures(qca_data.lat_atat, (uint)qca_data.max_num_atoms, qca_data.rundirpath);
    }
    else { // run ATAT
      qca_data.lat_atat = createLatForATAT(qca_data.plattice, qca_data.elements);
      qca_data.vstr_atat = getATATXstructures(qca_data.lat_atat, (uint)qca_data.max_num_atoms);
      qca_data.vstr_aflow = getAFLOWXstructures(qca_data.plattice, qca_data.elements, qca_data.num_threads);
      qca_data.mapstr = getMapForXstructures(getATATXstructures(qca_data.lat_atat, qca_data.aflow_max_num_atoms), qca_data.vstr_aflow, qca_data.num_threads); // map ATAT xstrs to AFLOW xstrs because ATAT cannot identify AFLOW xstrs
      generateFilesForATAT(qca_data.rundirpath, qca_data.lat_atat, qca_data.vstr_aflow, qca_data.vstr_atat, qca_data.mapstr);
      runATAT(qca_data.workdirpath, qca_data.rundirpath, qca_data.min_sleep);
    }
    qca_data.cv_cluster = getCVCluster(qca_data.rundirpath, qca_data.cv_cut);
    qca_data.num_atom_cluster = getNumAtomCluster(qca_data.vstr_atat);
    qca_data.conc_cluster = getConcentrationCluster(qca_data.rundirpath, qca_data.vstr_atat.size(), qca_data.elements.size());
    qca_data.excess_energy_cluster = getExcessEnergyCluster(qca_data.rundirpath, qca_data.conc_cluster, qca_data.max_num_atoms);
    setCongruentClusters(qca_data);
    qca_data.degeneracy_cluster = getDegeneracyCluster(qca_data.plattice, qca_data.vstr_atat, qca_data.elements, qca_data.max_num_atoms, true, qca_data.rundirpath);
    qca_data.conc_macro = getConcentrationMacro(qca_data.conc_curve_range, qca_data.conc_npts, qca_data.elements.size());
    qca_data.temp = getTemperature(qca_data.temp_range, qca_data.temp_npts);
    vector<double> data_ec = getRelativeEntropyEC(qca_data.conc_cluster, qca_data.degeneracy_cluster, qca_data.excess_energy_cluster, qca_data.temp, qca_data.max_num_atoms);
    qca_data.rel_s_ec = data_ec[0]; qca_data.temp_ec = data_ec[1];
    if (qca_data.calc_binodal) {
      qca_data.prob_ideal_cluster = getProbabilityIdealCluster(qca_data.conc_macro, qca_data.conc_cluster, qca_data.degeneracy_cluster, qca_data.max_num_atoms);
      checkProbability(qca_data.conc_macro, qca_data.conc_cluster, qca_data.prob_ideal_cluster);
      setProbabilityCluster(qca_data.conc_macro, qca_data.conc_cluster, qca_data.excess_energy_cluster, qca_data.prob_ideal_cluster, qca_data.temp, qca_data.max_num_atoms, qca_data.prob_cluster);
      try {
        checkProbability(qca_data.conc_macro, qca_data.conc_cluster, qca_data.prob_ideal_cluster, qca_data.prob_cluster, qca_data.temp);
      }
      catch (aurostd::xerror& err) {
        return;
      }
      qca_data.rel_s = getRelativeEntropy(qca_data.prob_cluster, qca_data.prob_ideal_cluster);
      qca_data.binodal_boundary = getBinodalBoundary(qca_data.rel_s, qca_data.rel_s_ec, qca_data.temp);
    }
    return;
  }
}

// ***************************************************************************
// qca::getBinodalBoundary
// ***************************************************************************
namespace qca {
  xvector<double> getBinodalBoundary(const xmatrix<double>& rel_s, const double rel_s_ec, const xvector<double>& temp) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    int nx = rel_s.rows, nt = rel_s.cols, n_fit = 8;
    xvector<double> binodal_boundary(nx), p, rr(n_fit), ri(n_fit);
    xvector<double> wts = aurostd::ones_xv<double>(nt);
    double temp_mean = aurostd::mean(temp), temp_std = aurostd::stddev(temp);
    xvector<double> temp_scaled = (temp - temp_mean) / temp_std; // scale for numerical stability
    for (int i = 1; i <= nx; i++) {
      p = aurostd::polynomialCurveFit(temp_scaled, rel_s(i) - rel_s_ec, n_fit, wts);
      aurostd::polynomialFindRoots(p, rr, ri);
      if (LDEBUG) {
        cerr << " i=" << i << " | p=" << p << endl;
        cerr << "   Real roots=" << rr << endl;
        cerr << "   Imag roots=" << ri << endl;
      }
      rr = temp_std * rr + temp_mean;
      for (int j = 1; j <= n_fit; j++) {
        if (rr(j) >= temp(1) && rr(j) <= temp(temp.rows) && aurostd::isequal(ri(j), 0.0) && binodal_boundary(i) < rr(j)) { // largest solution must be real and within temp range
          binodal_boundary(i) = rr(j);
        }
      }
    }
    return binodal_boundary;
  }
}

// ***************************************************************************
// qca::getRelativeEntropy
// ***************************************************************************
namespace qca {
  xmatrix<double> getRelativeEntropy(const vector<xmatrix<double>>& prob_cluster, const xmatrix<double>& prob_cluster_ideal) {
    int nx = prob_cluster_ideal.rows, ncl = prob_cluster_ideal.cols, nt = prob_cluster.size();
    xmatrix<double> rel_s(nx, nt);
    for (int i = 1; i <= nx; i++) {
      for (int j = 1; j <= nt; j++) {
        for (int k = 1; k <= ncl; k++) {
          rel_s(i, j) += prob_cluster[j - 1](i, k) * aurostd::log(prob_cluster[j - 1](i, k) / prob_cluster_ideal(i, k));
        }
      }
    }
    return rel_s;
  }
}

// ***************************************************************************
// qca::getRelativeEntropyEC
// ***************************************************************************
namespace qca {
  vector<double> getRelativeEntropyEC(const xmatrix<double>& conc_cluster, const xvector<int>& degeneracy_cluster, const xvector<double>& excess_energy_cluster, const xvector<double>& _temp, const int max_num_atoms, bool interp) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    double rel_s_ec = 0.0;
    int nelem = conc_cluster.cols, n_fit = 8;
    xvector<double> temp = _temp;
    xmatrix<double> conc_macro_ec(1, nelem);
    conc_macro_ec.set(1.0 / (double)nelem);
    xmatrix<double> prob_ideal_ec = getProbabilityIdealCluster(conc_macro_ec, conc_cluster, degeneracy_cluster, max_num_atoms);
    vector<xmatrix<double>> prob_ec;
    // Check that the probability is physical, if not, shift the temperature range upward
    bool shift_temp = false;
    double dtemp = temp(2) - temp(1);
    while (!setProbabilityCluster(conc_macro_ec, conc_cluster, excess_energy_cluster, prob_ideal_ec, temp, max_num_atoms, prob_ec) && 
           aurostd::min(temp) < 1e4) {
      shift_temp = true;
      temp += dtemp;
    }
    if (shift_temp && aurostd::min(temp) > DEFAULT_QCA_TEMP_MIN_LIMIT) {
      string message = "Could not find a temperature at equi-concentration that leads to a physical solution";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    xvector<double> order_param(prob_ec.size());
    xmatrix<double> m1, m2, m3 = prob_ideal_ec * aurostd::trasp(prob_ideal_ec);
    for (uint i = 0; i < prob_ec.size(); i++) {
      m1 = prob_ec[i] * aurostd::trasp(prob_ideal_ec);
      m2 = prob_ec[i] * aurostd::trasp(prob_ec[i]);
      order_param(i + 1) = m1(1, 1) / (aurostd::sqrt(m2(1, 1)) * aurostd::sqrt(m3(1, 1)) + _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);
    }
    // Check for divergence of the order parameter at low temperature
    bool low_t_div = false;
    if (order_param.rows > 3) {
      double d1 = (order_param(3) - order_param(1)) / (temp(3) - temp(1)), d2 = (order_param(4) - order_param(2)) / (temp(4) - temp(2));
      if (std::abs(d1) > std::abs(d2)) {low_t_div = true;}
    }
    double temp_mean = aurostd::mean(temp), temp_std = aurostd::stddev(temp);
    xvector<double> temp_scaled = (temp - temp_mean) / temp_std; // scale for numerical stability
    xvector<double> wts = aurostd::ones_xv<double>(order_param.rows);
    xvector<double> p = aurostd::polynomialCurveFit(temp_scaled, order_param, n_fit, wts);
    xvector<double> p1 = aurostd::evalPolynomialCoeff(p, 1), p2 = aurostd::evalPolynomialCoeff(p, 2);
    // If the order parameter diverges, then we can interpolate
    if (interp && low_t_div) {
      xvector<double> dadt = aurostd::evalPolynomial_xv(temp_scaled, p1), temp_interp(temp.rows + 1), dadt_interp(temp.rows + 1);
      for (int i = 1; i <= temp.rows; i++) {temp_interp(i + 1) = temp(i);} // add T = 0K
      for (int i = 1; i <= dadt.rows; i++) {dadt_interp(i + 1) = dadt(i);} // add Da/DT(0K) = 0
      temp_mean = aurostd::mean(temp_interp), temp_std = aurostd::stddev(temp_interp);
      temp_scaled = (temp_interp - temp_mean) / temp_std;
      wts = aurostd::exp(temp_scaled); // scale weights by temperature, since high-T is more accurate
      wts(1) = wts(wts.rows); // force convergence at T = 0K
      p1 = aurostd::polynomialCurveFit(temp_scaled, dadt_interp, n_fit - 1, wts);
      p2 = aurostd::evalPolynomialCoeff(p1, 1);
    }
    if (LDEBUG) {
      cerr << "shifted=" << shift_temp << endl;
      cerr << "interpolated=" << (interp && low_t_div) << endl;
      cerr << "alpha_orig=" << order_param << endl;
      cerr << "D[alpha, 0]=" << aurostd::evalPolynomial_xv(temp_scaled, p) << endl;
      cerr << "D[alpha, 1]=" << aurostd::evalPolynomial_xv(temp_scaled, p1) << endl;
      cerr << "D[alpha, 2]=" << aurostd::evalPolynomial_xv(temp_scaled, p2) << endl;
    }
    xvector<double> rr(n_fit - 2), ri(n_fit - 2);
    aurostd::polynomialFindRoots(p2, rr, ri);
    if (LDEBUG) {
      cerr << " p2=" << p2 << endl;
      cerr << "   Real roots=" << rr << endl;
      cerr << "   Imag roots=" << ri << endl;
    }
    vector<double> temp_ec;
    for (int j = 1; j <= n_fit; j++) {
      if (rr(j) >= temp_scaled(1) && rr(j) <= temp_scaled(temp.rows) && aurostd::isequal(ri(j), 0.0)) { // solution must be real and within temp range
        if (temp_ec.empty()) {
          temp_ec.push_back(rr(j));
        }
        else if (aurostd::evalPolynomial(temp_ec[0], p1) < aurostd::evalPolynomial(rr(j), p1)) { // keep largest gradient
          temp_ec[0] = rr(j);
        }
      }
    }
    temp_ec[0] = temp_std * temp_ec[0] + temp_mean;
    if (LDEBUG) {cerr << "T_ec=" << temp_ec[0] << "K" << endl;}
    setProbabilityCluster(conc_macro_ec, conc_cluster, excess_energy_cluster, prob_ideal_ec, aurostd::vector2xvector(temp_ec), max_num_atoms, prob_ec);
    for (int i = 1; i <= prob_ideal_ec.cols; i++) {
      rel_s_ec += prob_ec[0](1, i) * aurostd::log(prob_ec[0](1, i) / prob_ideal_ec(1, i));
    }
    return {rel_s_ec, temp_ec[0]};
  }
}

// ***************************************************************************
// qca::setProbabilityCluster
// ***************************************************************************
namespace qca {
  bool setProbabilityCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<double>& excess_energy_cluster, const xmatrix<double>& prob_ideal_cluster, const xvector<double>& temp, const int max_num_atoms, vector<xmatrix<double>>& prob_cluster) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    prob_cluster.clear();
    int nx = prob_ideal_cluster.rows, ncl = prob_ideal_cluster.cols, nt = temp.rows, neqs = conc_cluster.cols - 1;
    xmatrix<double> zeros(nx, ncl);
    for (uint it = 0; it < (uint)nt; it++) {prob_cluster.push_back(zeros);} // initialize
    xvector<double> beta = aurostd::pow(KBOLTZEV * temp, -1.0), p(max_num_atoms + 1), rr(max_num_atoms), ri(max_num_atoms), soln(neqs);
    xmatrix<int> natom_cluster = aurostd::xmatrixdouble2utype<int>((double)max_num_atoms * conc_cluster);
    bool soln_found = false;
    if (neqs == 1) {
      for (int it = 1; it <= nt; it++) {
        for (int i = 1; i <= nx; i++) {
          p.reset();
          soln.reset();
          for (int j = 1; j <= ncl; j++) {
            p(natom_cluster(j, 1) + 1) += prob_ideal_cluster(i, j) * std::exp(-beta(it) * excess_energy_cluster(j)) * (conc_cluster(j, 1) - conc_macro(i, 1));
          }
          aurostd::polynomialFindRoots(p, rr, ri);
          if (LDEBUG) {
            cerr << "it=" << it << " i=" << i << " | p=" << p << endl;
            cerr << "   Real roots=" << rr << endl;
            cerr << "   Imag roots=" << ri << endl;
          }
          for (int k = 1; k <= max_num_atoms; k++) {
            if (rr(k) > soln(1) && aurostd::isequal(ri(k), 0.0) && ncl * std::pow(rr(k), max_num_atoms) != INFINITY) { // solution must be positive, real and finite
              soln(1) = rr(k);
              soln_found = true;
            }
          }
          if (!soln_found) { // physical solution does not exist
            if (LDEBUG) {cerr << "Physical equilibrium probability does not exist for T=" << temp(it) << "K, X=[" << conc_macro(i) << " ]";}
            return false;
          }
          for (int j = 1; j <= ncl; j++) {
            prob_cluster[it - 1](i, j) = prob_ideal_cluster(i, j) * std::exp(-beta(it) * excess_energy_cluster(j)) * std::pow(soln(1), natom_cluster(j, 1));
          }
          prob_cluster[it - 1].setmat(prob_cluster[it - 1].getmat(i, i, 1, ncl) / aurostd::sum(prob_cluster[it - 1].getmat(i, i, 1, ncl)), i, 1); // normalize sum to 1
          if (LDEBUG) {cerr << "it=" << it << " i=" << i << " | SUM[P_cluster]=" << aurostd::sum(prob_cluster[it - 1].getmat(i, i, 1, ncl)) << endl;}
        }
      }
    }
    else { // homotopy continuation
    }
    return true;
  }
}

// ***************************************************************************
// qca::getProbabilityIdealCluster
// ***************************************************************************
// P_j(X) = g_j*(X1^N1_j)*(X2^N2_j)*...(X(K)^N(K)_j)
namespace qca {
  xmatrix<double> getProbabilityIdealCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<int>& degeneracy_cluster, const int max_num_atoms) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    int ncl = conc_cluster.rows, nx = conc_macro.rows, nelem = conc_macro.cols;
    xmatrix<double> prob_ideal_cluster(nx, ncl);
    for (int i = 1; i <= nx; i++) {
      for (int j = 1; j <= ncl; j++) {
        prob_ideal_cluster(i, j) = (double)degeneracy_cluster(j);
        for (int k = 1; k <= nelem; k++) {
          prob_ideal_cluster(i, j) *= std::pow(conc_macro(i, k), conc_cluster(j, k) * max_num_atoms);
        }
      }
      prob_ideal_cluster.setmat(prob_ideal_cluster.getmat(i, i, 1, ncl) / aurostd::sum(prob_ideal_cluster.getmat(i, i, 1, ncl)), i, 1); // normalize sum to 1
      if (LDEBUG) {cerr << "i=" << i << " | SUM[P_cluster]=" << aurostd::sum(prob_ideal_cluster.getmat(i, i, 1, ncl)) << endl;}
    }
    return prob_ideal_cluster;
  }
}

// ***************************************************************************
// qca::checkProbability
// ***************************************************************************
namespace qca {
  void checkProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob) {
    int nx = prob.rows;
    for (int i = 1; i <= nx; i++) {
      if (!aurostd::isequal(aurostd::sum(prob(i)), 1.0)) { // unnormalized
        stringstream message;
        message << "Ideal solution (high-T) probability is unnormalized for i=" << i;
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
      }
      else if (!aurostd::isequal(prob(i)*conc_cluster, conc_macro(i))) { // does not satisfy concentration constraints
        stringstream message;
        message << "Ideal solution (high-T) probability does not satisfy concentration contraint for X=[" << conc_macro(i) << " ]";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
      }
    }
  }

  void checkProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob_ideal_cluster, const vector<xmatrix<double>>& prob, const xvector<double>& temp) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    int nx = prob[0].rows;
    double diff = AUROSTD_MAX_DOUBLE, diff_old;
    for (int it = 1; it <= (int)prob.size(); it++) {
      for (int i = 1; i <= nx; i++) {
        if (!aurostd::isequal(aurostd::sum(prob[it - 1](i)), 1.0)) { // unnormalized
          stringstream message;
          message << "Equilibrium probability is unnormalized for T=" << temp(it) << "K, i=" << i;
          throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
        }
        else if (!aurostd::isequal(prob[it - 1](i)*conc_cluster, conc_macro(i))) { // does not satisfy concentration constraints
          stringstream message;
          message << "Equilibrium probability does not satisfy concentration contraint for T=" << temp(it) << "K, X=[" << conc_macro(i) << " ]";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
        }
      }
      diff_old = diff;
      diff = aurostd::sum(aurostd::abs(prob[it - 1] - prob_ideal_cluster));
      if (LDEBUG) {cerr << "|P - P0| = " << diff << endl;}
      if (diff > diff_old) { // P(T_inf) does not equal P0(T_inf)
        stringstream message;
        message << "Equilibrium probability does not converge to the high-T limit for T=" << temp(it) << "K";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }
  }
}

// ***************************************************************************
// qca::getConcentrationMacro
// ***************************************************************************
namespace qca {
  xmatrix<double> getConcentrationMacro(const vector<double>& conc_curve_range, const int conc_npts, const uint nelem) {
    xmatrix<double> conc_macro;
    if (!conc_curve_range.empty()) { // curve in concentration space
      xmatrix<double> _conc_macro(conc_npts, nelem);
      for (uint i = 0; i < nelem; i++) {
        _conc_macro.setcol(aurostd::linspace(conc_curve_range[2 * i], conc_curve_range[2 * i + 1], conc_npts), i + 1);
      }
      conc_macro = _conc_macro;
    }
    else {
      if (nelem == 2) { // trivial
        xmatrix<double> _conc_macro(conc_npts, nelem);
        xvector<double> x = aurostd::linspace(CONC_SHIFT, 1.0 - CONC_SHIFT, conc_npts);
        _conc_macro.setcol(x, 1);
        _conc_macro.setcol(1.0 - x, 2);
        conc_macro = _conc_macro;
      }
      else if (nelem == 3) { // barycentric coordinates
        xmatrix<double> _conc_macro((conc_npts * conc_npts - conc_npts) / 2, nelem); // always even
        xvector<double> x = aurostd::linspace(CONC_SHIFT, 1.0 - CONC_SHIFT, conc_npts);
        xvector<double> p1(2), p2(2), p3(2), p4(2);
        p1(1) = 0.0; p1(2) = 0.0;
        p2(1) = 1.0; p2(2) = 0.0;
        p3(1) = 1.0; p3(2) = 1.0;
        int ii = 1;
        for (int i = 1; i <= x.rows; i++) {
          for (int j = i + 1; j <= x.rows; j++) {
            p4(1) = x(j); p4(2) = x(i);
            _conc_macro(ii, 1) = (p1(1) - p3(1)) * (p4(2) - p1(2)) - (p1(1) - p4(1)) * (p3(2) - p1(2));
            _conc_macro(ii, 2) = (p1(1) - p4(1)) * (p2(2) - p1(2)) - (p1(1) - p2(1)) * (p4(2) - p1(2));
            _conc_macro(ii, 3) = (p4(1) - p3(1)) * (p2(2) - p4(2)) - (p4(1) - p2(1)) * (p3(2) - p4(2));
            ii++;
          }
        }
        conc_macro = _conc_macro;
      }
      else { // for greater than ternary only evaluate at equi-concentration
        xmatrix<double> _conc_macro(1, nelem);
        _conc_macro.set(1.0 / (double)nelem);
        conc_macro = _conc_macro;
      }
    }
    return conc_macro;
  }
}

// ***************************************************************************
// qca::getTemperature
// ***************************************************************************
namespace qca {
  xvector<double> getTemperature(const vector<double>& temp_range, const int temp_npts) {
    return aurostd::linspace(temp_range[0], temp_range[1], temp_npts);
  }
}

// ***************************************************************************
// qca::setCongruentClusters
// ***************************************************************************
namespace qca {
  void setCongruentClusters(_qca_data& qca_data) {
    vector<int> indx_cluster;
    for (int i = 1; i <= qca_data.num_atom_cluster.rows; i++) {
      if (!(qca_data.max_num_atoms % qca_data.num_atom_cluster(i))) {indx_cluster.push_back(i);}
    }
    int ncl = indx_cluster.size(), nelem = qca_data.elements.size();
    vector<xstructure> _vstr_atat(ncl);
    xvector<int> v1(ncl);
    xvector<double> v2(ncl);
    xmatrix<double> m1(ncl, nelem);
    for (int i = 0; i < ncl; i++) {
      _vstr_atat[i] = qca_data.vstr_atat[indx_cluster[i] - 1];
      m1.setmat(qca_data.conc_cluster.getmat(indx_cluster[i], indx_cluster[i], 1, nelem), i + 1, 1);
      v1(i + 1) = qca_data.num_atom_cluster(indx_cluster[i]);
      v2(i + 1) = qca_data.excess_energy_cluster(indx_cluster[i]);
    }
    qca_data.vstr_atat = _vstr_atat;
    qca_data.conc_cluster = m1;
    qca_data.num_atom_cluster = v1;
    qca_data.excess_energy_cluster = v2;
  }
}

// ***************************************************************************
// qca::getExcessEnergyCluster
// ***************************************************************************
namespace qca {
  xvector<double> getExcessEnergyCluster(const string& rundirpath, const xmatrix<double>& conc_cluster, const int max_num_atoms) {
    int ind, nstr = conc_cluster.rows, nelem = conc_cluster.cols;
    vector<string> vinput, tokens;
    aurostd::file2vectorstring(rundirpath + "/ref_energy.out", vinput);
    xvector<double> nrg_ref = aurostd::vector2xvector(aurostd::vectorstring2vectordouble(vinput));
    xvector<double> nrg(nstr);
    aurostd::file2vectorstring(rundirpath + "/fit.out", vinput);
    // Read energies per atom
    for (uint line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      ind = aurostd::string2utype<int>(tokens[tokens.size() - 1]) + 1;
      nrg(ind) = aurostd::string2utype<double>(tokens[nelem]);
    }
    aurostd::file2vectorstring(rundirpath + "/predstr.out", vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      ind = aurostd::string2utype<int>(tokens[tokens.size() - 2]) + 1;
      nrg(ind) = aurostd::string2utype<double>(tokens[nelem]);
    }
    return nrg * max_num_atoms; // per cell
  }
}

// ***************************************************************************
// qca::getConcentrationCluster
// ***************************************************************************
namespace qca {
  xmatrix<double> getConcentration(const vector<string>& elements, const vector<xstructure>& vstr) {
    uint nstr = vstr.size(), nelem = elements.size();
    xmatrix<double> conc_cluster(nstr, nelem);
    int ie = -1;
    xvector<double> stoich(nelem);
    vector<string> str_elements;
    for (uint i = 0; i < vstr.size(); i++) {
      if (nelem != vstr[i].stoich_each_type.size()) {
        str_elements = vstr[i].GetElements(true, true);
        stoich.reset();
        for (uint j = 0; j < nelem; j++) {
          if (aurostd::WithinList(str_elements, elements[j], ie)) {stoich(j + 1) = vstr[i].stoich_each_type[ie];}
        }
      }
      else {
        stoich = aurostd::vector2xvector(aurostd::deque2vector(vstr[i].stoich_each_type));
      }
      for (uint j = 0; j < nelem; j++) {conc_cluster(i + 1, j + 1) = stoich(j + 1);}
    }
    return conc_cluster;
  }

  xmatrix<double> getConcentrationCluster(const string& rundirpath, const int nstr, const int nelem) {
    xmatrix<double> conc_cluster(nstr, nelem);
    vector<string> vinput, tokens;
    aurostd::file2vectorstring(rundirpath + "/fit.out", vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      for (int i = 1; i <= nelem; i++) {
        conc_cluster(aurostd::string2utype<int>(tokens[tokens.size() - 1]) + 1, i) = aurostd::string2utype<double>(tokens[i - 1]);
      }
    }
    aurostd::file2vectorstring(rundirpath + "/predstr.out", vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      for (int i = 1; i <= nelem; i++) {
        conc_cluster(aurostd::string2utype<int>(tokens[tokens.size() - 2]) + 1, i) = aurostd::string2utype<double>(tokens[i - 1]);
      }
    }
    return conc_cluster;
  }
}

// ***************************************************************************
// qca::getNumAtomCluster
// ***************************************************************************
namespace qca {
  xvector<int> getNumAtomCluster(const vector<xstructure>& vstr) {
    int nstr = vstr.size(), natom;
    xvector<int> num_atom_cluster(nstr);
    for (int i = 1; i <= nstr; i++) {
      natom = 0;
      for (uint j = 0; j < vstr[i - 1].num_each_type.size(); j++) {natom += vstr[i - 1].num_each_type[j];}
      num_atom_cluster(i) = natom;
    }
    return num_atom_cluster;
  }
}

// ***************************************************************************
// qca::getDegeneracyCluster
// ***************************************************************************
namespace qca {
  xvector<int> getDegeneracyCluster(const string& plattice, const vector<xstructure>& _vstr, const vector<string>& elements, const int max_num_atoms, const bool shuffle, const string& rundirpath) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string filepath = "";
    if (!rundirpath.empty()) {
      vector<string> tokens;
      aurostd::string2tokens(rundirpath, tokens, "/");
      filepath = "/" + QCA_FILE_PREFIX + "degen_" + aurostd::utype2string<int>(max_num_atoms) + "_atom.txt";
      for (int i = tokens.size() - 2; i >= 0; i--) {filepath = "/" + tokens[i] + filepath;}
      if (aurostd::file2vectorstring(filepath, tokens)) {return aurostd::vector2xvector(aurostd::vectorstring2vectorint(tokens));}
    }
    vector<xstructure> vstr = _vstr;
    vector<uint> index;
    for (uint i = 0; i < vstr.size(); i++) {index.push_back(i);}
    // Shuffling the xstructures is important because it can avoid edge case scenarios where the reference centroid
    // in Xtalfinder does not pick up the proper neighbors
    if (shuffle) { // introduce randomness into grouping
      aurostd::random_shuffle(index);
      for (uint i = 0; i < index.size(); i++) {vstr[i] = _vstr[index[i]];}
    }
    xvector<int> degeneracy_cluster(vstr.size());
    xstructure str_plat, str;
    deque<_atom> atoms;
    _atom atom;
    xmatrix<double> lattice(3,3);
    // Define primitive parent lattice and atom atoms
    if (plattice == "fcc") {
      lattice(1, 1) = 0.0; lattice(1, 2) = 0.5; lattice(1, 3) = 0.5;
      lattice(2, 1) = 0.5; lattice(2, 2) = 0.0; lattice(2, 3) = 0.5;
      lattice(3, 1) = 0.5; lattice(3, 2) = 0.5; lattice(3, 3) = 0.0;
      str_plat.lattice = lattice;
      atom.clear();
      atom.fpos(1) = 0.0; atom.fpos(2) = 0.0; atom.fpos(3) = 0.0;
      atom.cpos = str_plat.f2c * atom.fpos;
      atoms.push_back(atom);
    }
    else if (plattice == "bcc") {
      lattice(1, 1) = -0.5; lattice(1, 2) = 0.5; lattice(1, 3) = 0.5;
      lattice(2, 1) = 0.5; lattice(2, 2) = -0.5; lattice(2, 3) = 0.5;
      lattice(3, 1) = 0.5; lattice(3, 2) = 0.5; lattice(3, 3) = -0.5;
      str_plat.lattice = lattice;
      atom.clear();
      atom.fpos(1) = 0.0; atom.fpos(2) = 0.0; atom.fpos(3) = 0.0;
      atom.cpos = str_plat.f2c * atom.fpos;
      atoms.push_back(atom);
    }
    else if (plattice == "hcp") {
      double sqrt3 = std::sqrt(3.0), sqrt8 = std::sqrt(8.0);
      lattice(1, 1) = 1.0; lattice(1, 2) = 0.0; lattice(1, 3) = 0.0;
      lattice(2, 1) = -0.5; lattice(2, 2) = 0.5 * sqrt3; lattice(2, 3) = 0.0;
      lattice(3, 1) = 0.0; lattice(3, 2) = 0.0; lattice(3, 3) = sqrt8 / sqrt3;
      str_plat.lattice = lattice;
      atom.clear();
      atom.fpos(1) = 0.0; atom.fpos(2) = 0.0; atom.fpos(3) = 0.0;
      atom.cpos = str_plat.f2c * atom.fpos;
      atoms.push_back(atom);
      atom.clear();
      atom.fpos(1) = 2.0 / 3.0; atom.fpos(2) = 1.0 / 3.0; atom.fpos(3) = 0.5;
      atom.cpos = str_plat.f2c * atom.fpos;
      atoms.push_back(atom);
    }
    uint natoms = atoms.size();
    for (uint i = 0; i < natoms; i++) {str_plat.AddAtom(atoms[i]);}
    str_plat.DecorateWithElements();
    str_plat.iomode = IOVASP_POSCAR;
    str_plat.partial_occupation_HNF = max_num_atoms;
    str_plat.partial_occupation_flag = true;
    str_plat.neg_scale_second = true;
    // Enumerate unique superlattices
    bool quiet = XHOST.QUIET;
    XHOST.QUIET = true;
    pocc::POccCalculator pcalc(str_plat);
    pcalc.calculateHNF();
    pcalc.getTotalPermutationsCount();
    pcalc.calculate();
    vector<xstructure> vstr_sup = pcalc.getUniqueDerivativeStructures();
    if (LDEBUG) {cerr << "Number of unique superlattices (HNF) = " << vstr_sup.size() << endl;}
    // Enumerate all indicies
    vector<vector<int>> all_indicies;
    aurostd::xcombos xc(elements.size(), max_num_atoms, 'E', true);
    while (xc.increment()) {all_indicies.push_back(xc.getCombo());}
    // Enumerate all derivative structures
    vector<xstructure> vstr_ds;
    uint itype;
    for (uint i = 0; i < vstr_sup.size(); i++) {
      lattice = vstr_sup[i].lattice;
      for (uint j = 0; j < all_indicies.size(); j++) {
        str.clear();
        atoms.clear();
        for (uint k = 0; k < vstr_sup[i].atoms.size(); k++) {
          atom.clear();
          atom.name = atom.cleanname = elements[all_indicies[j][k]];
          atom.cpos = vstr_sup[i].atoms[k].cpos;
          atom.fpos = vstr_sup[i].c2f * atom.cpos;
          atom.name_is_given = TRUE;
          atoms.push_back(atom);
        }
        std::stable_sort(atoms.begin(), atoms.end(), sortAtomsNames);
        itype = 0;
        atoms[0].type = itype;
        for (uint k = 1; k < atoms.size(); k++) {
          if (atoms[k].name != atoms[k - 1].name) {itype++;}
          atoms[k].type = itype;
        }
        str.scale = 1.0;
        str.neg_scale = false;
        str.lattice = lattice;
        str.AddAtom(atoms);
        vstr_ds.push_back(str);
      }
    }
    if (LDEBUG) {cerr << "Number of total derivative structures = " << vstr_ds.size() << endl;}
    // Find degenerate structures
    vstr_ds.insert(vstr_ds.begin(), vstr.begin(), vstr.end()); // concatenate xstructures, subtract by 1 in the end
    XtalFinderCalculator xtal_calc;
    vector<vector<uint>> dsg = xtal_calc.groupSimilarXstructures(vstr_ds); // costly part of the function
    for (uint i = 0; i < dsg.size(); i++) {
      std::sort(dsg[i].begin(), dsg[i].end()); // first index is the cluster index
      if (dsg[i][0] < (uint)degeneracy_cluster.rows) {
        degeneracy_cluster(dsg[i][0] + 1) = dsg[i].size() - 1;
        if (LDEBUG) {
          cerr << "i=" << i << " | DG=" << dsg[i].size() - 1 << " | ";
          for (uint j = 0; j < dsg[i].size(); j++) {cerr << dsg[i][j] << " ";}
          cerr << endl;
        }
      }
    }
    XHOST.QUIET = quiet;
    int sum_calculated = aurostd::sum(degeneracy_cluster), sum_accepted = (int)vstr_sup.size() * (int)std::pow(elements.size(), max_num_atoms);
    if (!aurostd::isequal(sum_calculated, sum_accepted)) {
      stringstream message;
      message << "Degeneracies do not satisfy the sum rule, the sum is " << sum_calculated << " but should be " << sum_accepted;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    if (shuffle) { // place back in correct order
      xvector<int> _degeneracy_cluster = degeneracy_cluster;
      for (uint i = 0; i < index.size(); i++) {degeneracy_cluster(index[i] + 1) = _degeneracy_cluster(i + 1);}
    }
    if (!filepath.empty()) {
      stringstream ss;
      for (int i = 1; i <= degeneracy_cluster.rows; i++) {ss << degeneracy_cluster(i) << endl;}
      aurostd::stringstream2file(ss, filepath);
    }
    return degeneracy_cluster;
  }
}

// ***************************************************************************
// qca::getCVCluster
// ***************************************************************************
namespace qca {
  double getCVCluster(const string& rundirpath, const double cv_cut) {
    vector<string> vinput, tokens;
    aurostd::file2vectorstring(rundirpath + "/maps.log", vinput);
    aurostd::string2tokens(vinput[vinput.size() - 1], tokens, " ");
    if (tokens[0] != "Crossvalidation") {
      string message = "ATAT run did not complete";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    double cv_cluster = aurostd::string2utype<double>(tokens[2]);
    if (cv_cluster > cv_cut) {
      stringstream message;
      message << "Cluster cross-validation score is above the cut-off, CV_cluster=" << cv_cluster << ", CV_cutoff=" << cv_cut;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    return cv_cluster;
  }
}

// ***************************************************************************
// qca::runATAT
// ***************************************************************************
namespace qca {
  void runATAT(const string& workdirpath, const string& rundirpath, const uint min_sleep) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    uint iter = 0;
    if (aurostd::substring2bool(aurostd::execute2string("mmaps", stdouterr_fsio), "command not found")) {
      string message = "Missing mmaps program";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    aurostd::execute("touch " + rundirpath + "/stop"); // pre-kill mmaps gracefully
    aurostd::RemoveFile(rundirpath + "/maps_is_running");
    aurostd::RemoveFile(rundirpath + "/maps.log");
    chdir(rundirpath.c_str());
    string tmpfile = aurostd::TmpStrCreate(), logstring = "";
    aurostd::execute("mmaps -d > " + tmpfile + " 2>&1 &");
    while (!aurostd::substring2bool(logstring, "true and predicted ground states agree") && !aurostd::substring2bool(logstring, "No other ground states")) {
      iter++;
      if (LDEBUG) {cerr << "Sleeping, iter=" << iter << endl;}
      if (iter > 30) { // wait 30 times the minimum sleep (60 seconds)
        string message = "mmaps is taking too long to predict structures, dir=" + rundirpath;
        aurostd::RemoveFile(tmpfile);
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      aurostd::Sleep(min_sleep);
      aurostd::file2string(rundirpath + "/maps.log", logstring);
    }
    aurostd::RemoveFile(tmpfile);
    chdir(workdirpath.c_str());
  }
}

// ***************************************************************************
// qca::generateFilesForATAT
// ***************************************************************************
namespace qca {
  void generateFilesForATAT(const string& rundirpath, const string& lat_atat, const vector<xstructure>& vstr_aflow, const vector<xstructure>& vstr_atat, const vector<int>& mapstr) {
    stringstream oss;
    vector<xstructure> vstr;
    // Generate lat.in file
    aurostd::string2file(lat_atat, rundirpath + "/lat.in");
    // Generate str.out and energy files
    uint msize = mapstr.size();
    for (uint i = 0; i < vstr_atat.size(); i++) {
      aurostd::DirectoryMake(rundirpath + "/" + aurostd::utype2string<uint>(i));
      oss << vstr_atat[i];
      aurostd::string2file(oss.str(), rundirpath + "/" + aurostd::utype2string<uint>(i) + "/str.out");
      aurostd::StringstreamClean(oss);
      if (i < msize && mapstr[i] != -1) {
        aurostd::string2file(aurostd::utype2string<double>(vstr_aflow[mapstr[i]].qm_E_cell) + "\n", rundirpath + "/" + aurostd::utype2string<uint>(i) + "/energy");
      }
    }
  }
}

// ***************************************************************************
// qca::getAFLOWXstructures
// ***************************************************************************
namespace qca {
  vector<xstructure> getAFLOWXstructures(const string& plattice, const vector<string>& elements, const int num_threads, bool use_sg) {
    vector<xstructure> vstr;
    vector<string> vstrlabel;
    string aflowlib, aflowurl;
    string alloyname = "";
    aflowlib::_aflowlib_entry entry;
    uint nelem = elements.size();
    stringstream oss;
    for (uint i = 0; i < nelem; i++) {alloyname += AVASP_Get_PseudoPotential_PAW_PBE(elements[i]);}
    aflowlib = "/common/LIB" + aurostd::utype2string<uint>(nelem) + "/RAW/" + alloyname;
    aflowurl = "aflowlib.duke.edu:AFLOWDATA/LIB" + aurostd::utype2string<uint>(nelem) + "_RAW/" + alloyname;
    if (plattice == "fcc") {
      if (nelem >= 2) {
        for (uint i = 1; i < 30; i++) {vstrlabel.push_back(aurostd::utype2string<uint>(i));}
      }
      if (nelem >= 3) {
        vstrlabel.push_back("TFCC001.ABC"); vstrlabel.push_back("TFCC002.ABC"); vstrlabel.push_back("TFCC003.ABC");
        for (uint i = 4; i < 17; i++) {
          vstrlabel.push_back("TFCC00" + aurostd::utype2string<uint>(i) + ".ABC");
          vstrlabel.push_back("TFCC00" + aurostd::utype2string<uint>(i) + ".BCA");
          vstrlabel.push_back("TFCC00" + aurostd::utype2string<uint>(i) + ".CAB");
        }
      }
    }
    else if (plattice == "bcc") {
      if (nelem >= 2) {
        for (uint i = 58; i < 87; i++) {vstrlabel.push_back(aurostd::utype2string<uint>(i));}
      }
      if (nelem >= 3) {
        vstrlabel.push_back("TBCC001.ABC"); vstrlabel.push_back("TBCC002.ABC"); vstrlabel.push_back("TBCC003.ABC");
        for (uint i = 4; i < 17; i++) {
          vstrlabel.push_back("TBCC00" + aurostd::utype2string<uint>(i) + ".ABC");
          vstrlabel.push_back("TBCC00" + aurostd::utype2string<uint>(i) + ".BCA");
          vstrlabel.push_back("TBCC00" + aurostd::utype2string<uint>(i) + ".CAB");
        }
      }
    }
    else if (plattice == "hcp") {
      if (nelem >= 2) {
        for (uint i = 115; i < 178; i++) {vstrlabel.push_back(aurostd::utype2string<uint>(i));}
      }
    }
    for (uint i = 0; i < vstrlabel.size(); i++) {
        entry.Load(aurostd::file2string(aflowlib + "/" + vstrlabel[i] + "/aflowlib.out"), oss);
        entry.aurl = aflowurl + "/" + vstrlabel[i];
        if (pflow::loadXstructures(entry, oss, false)) { // initial = unrelaxed; final = relaxed
          entry.vstr[0].qm_E_cell = entry.enthalpy_cell; // ATAT needs energy per cell
          if (use_sg && (entry.spacegroup_orig == entry.spacegroup_relax)) {
            vstr.push_back(entry.vstr[0]);
          }
          else if (compare::structuresMatch(entry.vstr[0], entry.vstr[entry.vstr.size() - 1], true, num_threads)) {
            vstr.push_back(entry.vstr[0]);
          }
        }
        entry.clear();
    }
    return vstr;
  }
}

// ***************************************************************************
// qca::createLatForATAT
// ***************************************************************************
namespace qca {
  string createLatForATAT(const string& plattice, const vector<string>& elements) {
    stringstream oss;
    uint nelem = elements.size();
    xmatrix<double> lattice(3,3);
    xvector<double> angles = 90.0 * aurostd::ones_xv<double>(3);
    xvector<double> coorsys = aurostd::ones_xv<double>(3);
    oss.precision(_DOUBLE_WRITE_PRECISION_MAX_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
    if (plattice == "fcc") {
      lattice(1, 1) = 0.0; lattice(2, 1) = 0.5; lattice(3, 1) = 0.5;
      lattice(1, 2) = 0.5; lattice(2, 2) = 0.0; lattice(3, 2) = 0.5;
      lattice(1, 3) = 0.5; lattice(2, 3) = 0.5; lattice(3, 3) = 0.0;
    }
    else if (plattice == "bcc") {
      lattice(1, 1) = -0.5; lattice(2, 1) = 0.5; lattice(3, 1) = 0.5;
      lattice(1, 2) = 0.5; lattice(2, 2) = -0.5; lattice(3, 2) = 0.5;
      lattice(1, 3) = 0.5; lattice(2, 3) = 0.5; lattice(3, 3) = -0.5;
    }
    else if (plattice == "hcp") {
      lattice = aurostd::eye<double>(3, 3);
      coorsys(3) = std::sqrt(8.0 / 3.0);
      angles(3) = 120.0;
    }
    for (uint i = 1; i <= 3; i++) {oss << DEFAULT_ATAT_ALAT * coorsys(i) << " ";}
    for (uint i = 1; i <= 3; i++) {oss << angles(i) << " ";}
    oss << endl;
    for (uint i = 1; i <= 3; i++) {
      for (uint j = 1; j <= 3; j++) {
        oss << lattice(i, j) << " ";
      }
      oss << endl;
    }
    oss << 0.0 << " " << 0.0 << " " << 0.0 << " ";
    for (uint i = 0; i < nelem; i++) {
      oss << elements[i] << ",";
    }
    if (plattice == "hcp") {
      oss << endl << 2.0 / 3.0 << " " << 1.0 / 3.0 << " " << 0.5 << " ";
      for (uint i = 0; i < nelem; i++) {
        oss << elements[i] << ",";
      }
    }
    oss << endl;
    return oss.str();
  }
}

// ***************************************************************************
// qca::getATATXstructures
// ***************************************************************************
namespace qca {
  vector<xstructure> getATATXstructures(const string& lat, const uint max_num_atoms, const string& rundirpath) {
    vector<xstructure> vstr;
    stringstream oss;
    if (!rundirpath.empty()) {
      vector<string> files;
      vector<xstructure> vstr_tmp;
      vector<uint> index;
      aurostd::DirectoryLS(rundirpath, files);
      for (uint i = 0; i < files.size(); i++) {
        if (aurostd::FileExist(rundirpath + "/" + files[i] + "/str.out")) {
          aurostd::file2stringstream(rundirpath + "/" + files[i] + "/str.out", oss);
          vstr_tmp.push_back(xstructure(oss, IOATAT_STR));
          index.push_back(aurostd::string2utype<uint>(files[i]));
          aurostd::StringstreamClean(oss);
        }
      }
      vstr.resize(index.size());
      for (uint i = 0; i < index.size(); i++) {vstr[index[i]] = vstr_tmp[i];}
      return vstr;
    }
    vector<string> vinput, tokens;
    string tmpfile = aurostd::TmpStrCreate();
    aurostd::string2file(lat, tmpfile);
    string sstr = aurostd::execute2string("genstr -n " + aurostd::utype2string<uint>(max_num_atoms) + " -l " + tmpfile, stdouterr_fsio);
    aurostd::RemoveFile(tmpfile);
    if (sstr.size() == 0 || aurostd::substring2bool(sstr, "Unable to open lattice file")) {
      string message = "Invalid lat.in file";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    else if (aurostd::substring2bool(sstr, "command not found")) {
      string message = "Missing genstr program";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    aurostd::string2vectorstring(sstr, vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      if (aurostd::substring2bool(vinput[line], "end")) {
        vstr.push_back(xstructure(oss, IOATAT_STR));
        aurostd::StringstreamClean(oss);
      }
      else {
        oss << vinput[line] << endl;
      }
    }
    return vstr;
  }
}

// ***************************************************************************
// qca::getMapForXstructures
// ***************************************************************************
namespace qca {
  vector<int> getMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, const int num_threads) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    vector<int> mapstr;
    bool match;
    for (uint i = 0; i < vstr1.size(); i++) {
      match = false;
      for (uint j = 0; j < vstr2.size() && !match; j++) {
        if (compare::structuresMatch(vstr1[i], vstr2[j], true, num_threads)) {
          if (LDEBUG) {cerr << "str1 i=" << i << " matched str2 j=" << j << endl << "str1=" << vstr1[i] << endl << "str2=" << vstr2[j] << endl;}
          mapstr.push_back(j); 
          match = true;
        }
      }
      if (!match) {mapstr.push_back(-1);}
    }
    return mapstr;
  }
}

// ***************************************************************************
// qca::displayUsage
// ***************************************************************************
namespace qca {
  void displayUsage(void) {
    vector<string> usage_options;
    usage_options.push_back("aflow --quasi_chem_approx|--qca --plattice=|--plat=fcc --elements=|--elem=Au,Pt[,Zn] [qca_options] [--directory=[DIRECTORY]]");
    usage_options.push_back(" ");
    usage_options.push_back("qca_options:");
    usage_options.push_back(" ");
    usage_options.push_back("GENERAL OPTIONS:");
    usage_options.push_back("--binodal");
    usage_options.push_back("--usage");
    usage_options.push_back("--screen_only");
    usage_options.push_back("--image_only|--image");
    usage_options.push_back("--no_plot|--noplot");
    usage_options.push_back(" ");
    usage_options.push_back("BINODAL OPTIONS:");
    usage_options.push_back("--max_num_atoms=|--mna=8");
    usage_options.push_back("--cv_cutoff=|--cv_cut=0.05");
    usage_options.push_back("--conc_curve_range=|--conc_curve=0,1,1,0");
    usage_options.push_back("--conc_npts=20");
    usage_options.push_back("--temp_range=|--temp=300,5000");
    usage_options.push_back("--temp_npts=150");
    usage_options.push_back(" ");
    usage_options.push_back("SPINODAL OPTIONS:");
    usage_options.push_back("--spinodal");
    usage_options.push_back(" ");
    usage_options.push_back("FORMAT OPTIONS:");
    usage_options.push_back("--format_data=|--data_format=txt|json");
    usage_options.push_back("--format_image=|--image_format=pdf|eps|png");
    usage_options.push_back(" ");
    init::MessageOption("--usage", "QCA()", usage_options);
    return;
  }
}

// ***************************************************************************
// qca::writeData
// ***************************************************************************
namespace qca {
  void writeData(const _qca_data& qca_data) {
    string filepath = qca_data.rundirpath + "/" + QCA_FILE_PREFIX + "output." + qca_data.format_data;
    stringstream output;
    if (qca_data.format_data == "txt") {
      string info_prefix = "";
      output << "QCA DATA" << endl;
      info_prefix = "Input data ";
      output << " " << info_prefix << "Alloy name                     = " << qca_data.alloyname << endl;
      output << " " << info_prefix << "Parent lattice                 = " << qca_data.plattice << endl;
      output << " " << info_prefix << "Macroscopic concentration      = " << endl << qca_data.conc_macro << endl;
      output << " " << info_prefix << "Temperature range (K)          = " << endl << trasp(qca_data.temp) << endl;
      output << " " << info_prefix << "Max atoms per cell             = " << qca_data.max_num_atoms << endl;
      info_prefix = "Cluster data ";
      output << " " << info_prefix << "CV score (eV)                  = " << qca_data.cv_cluster << endl;
      output << " " << info_prefix << "Number of atoms                = " << endl << trasp(qca_data.num_atom_cluster) << endl;
      output << " " << info_prefix << "Degeneracy                     = " << endl << trasp(qca_data.degeneracy_cluster) << endl;
      output << " " << info_prefix << "Concentration                  = " << endl << qca_data.conc_cluster << endl;
      output << " " << info_prefix << "Excess energy (eV)             = " << endl << trasp(qca_data.excess_energy_cluster) << endl;
      info_prefix = "Thermo data ";
      output << " " << info_prefix << "EC transition temperature (K)  = " << qca_data.temp_ec << endl;
      output << " " << info_prefix << "Binodal boundary (K)           = " << endl << trasp(qca_data.binodal_boundary) << endl;
      if (!qca_data.screen_only) {
        aurostd::stringstream2file(output, filepath);
        return;
      }
    }
    else if (qca_data.format_data == "json") {
      aurostd::JSONwriter json;
      json.addString("Alloy name", qca_data.alloyname);
      json.addVector("Elements", qca_data.elements);
      json.addString("Parent lattice", qca_data.plattice);
      json.addNumber("Concentration curve", qca_data.conc_curve);
      json.addMatrix("Macroscopic concentration", qca_data.conc_macro);
      json.addVector("Temperature range (K)", qca_data.temp);
      json.addNumber("Max atoms per cell", qca_data.max_num_atoms);
      json.addNumber("Cluster CV score (eV)", qca_data.cv_cluster);
      json.addVector("Cluster number of atoms", aurostd::xvectorutype2double(qca_data.num_atom_cluster));
      json.addVector("Cluster degeneracy", aurostd::xvectorutype2double(qca_data.degeneracy_cluster));
      json.addMatrix("Cluster concentration", qca_data.conc_cluster);
      json.addVector("Cluster excess energy (eV)", qca_data.excess_energy_cluster);
      json.addNumber("EC transition temperature (K)", qca_data.temp_ec);
      json.addVector("Binodal boundary (K)", qca_data.binodal_boundary);
      aurostd::string2file(json.toString(), filepath);
      return;
    }
    cout << output.str() << endl;
    return;
  }
}

// ***************************************************************************
// qca::readData
// ***************************************************************************
namespace qca {
  void readData(_qca_data& qca_data) {
    string filepath = qca_data.rundirpath + "/" + QCA_FILE_PREFIX + "output.json";
    if (!aurostd::FileExist(filepath)) {
      string message = "The JSON file does not exist, filepath=" + filepath; 
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _FILE_ERROR_);
    }
    string jsonfile = aurostd::file2string(filepath);
    qca_data.alloyname = aurostd::extractJsonValueAflow(jsonfile, "Alloy name");
    qca_data.conc_curve = aurostd::string2utype<bool>(aurostd::extractJsonValueAflow(jsonfile, "Concentration curve"));
    qca_data.conc_macro = aurostd::vectorvector2xmatrix<double>(aurostd::extractJsonMatrixAflow(jsonfile, "Macroscopic concentration"));
    qca_data.temp = aurostd::vector2xvector<double>(aurostd::extractJsonVectorAflow(jsonfile, "Temperature range (K)"));
    qca_data.binodal_boundary = aurostd::vector2xvector<double>(aurostd::extractJsonVectorAflow(jsonfile, "Binodal boundary (K)"));
  }
}

// ***************************************************************************
// qca::plotData
// ***************************************************************************
namespace qca {
  void plotData(const _qca_data& qca_data) {
    if (qca_data.format_image.empty()) {return;}
    stringstream gpfile;
    aurostd::xoption plotoptions;
    vector<vector<double>> data;
    // BEGIN - set standard
    plotoptions.push_attached("BACKGROUND_COLOR", "#FFFFFF");
    plotoptions.push_attached("GRID_COLOR", "#808080");
    plotoptions.push_attached("GRID_WIDTH", "1");
    plotoptions.push_attached("GRID_LINE_TYPE", "0");
    plotoptions.push_attached("PLOT_SIZE", "5.333, 3");
    // END - set standard
    // START - set default
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    plotoptions.push_attached("FILE_NAME", qca_data.rundirpath + "/" + QCA_FILE_PREFIX + "plot");
    plotoptions.push_attached("FILE_NAME_LATEX", qca_data.rundirpath + "/" + QCA_FILE_PREFIX + "plot");
    plotoptions.push_attached("IMAGE_FORMAT", qca_data.format_image);
    plotoptions.push_attached("PLOT_TITLE", qca_data.alloyname);
    plotoptions.push_attached("XLABEL", "Concentration");
    plotoptions.push_attached("XMIN", "0.0");
    plotoptions.push_attached("XMAX", "1.0");
    plotoptions.push_attached("XUNIT", "");
    plotoptions.push_attached("YLABEL", "Temperature");
    plotoptions.push_attached("YMIN", "1.0");
    plotoptions.push_attached("YMAX", aurostd::utype2string(1.2 * aurostd::max(qca_data.binodal_boundary)));
    plotoptions.push_attached("YUNIT", "K");
    plotoptions.push_attached("TITLES", "Binodal,Spinodal");
    plotoptions.push_attached("COLORS", "#000000");
    plotoptions.push_attached("LINETYPES", "-1,7");
    // END - set default
    plotter::generateHeader(gpfile, plotoptions, false);
    if (qca_data.conc_curve) { // generic x-axis for the concentration curve
      xvector<double> a = aurostd::linspace(0.0, 1.0, qca_data.conc_macro.rows);
      for (int i = 1; i <= qca_data.conc_macro.rows; i++) {data.push_back({a(i), qca_data.binodal_boundary(i)});}
        // START - custom
        gpfile << "# Custom" << endl;
        gpfile << "set label '\\shortstack{c}{$" << qca_data.elements[0] << "$=" << qca_data.conc_macro(1, 1);
        for (uint i = 1; i < qca_data.elements.size(); i++) {
          gpfile << "\\\\$" << qca_data.elements[i] << "$=" << qca_data.conc_macro(1, i + 1);
        }
        gpfile << "}' at 0.0,0.0 center offset 0.0,-2.5" << endl;
        gpfile << "set label '\\shortstack{c}{$" << qca_data.elements[0] << "$=" << qca_data.conc_macro(qca_data.conc_macro.rows, 1);
        for (uint i = 1; i < qca_data.elements.size(); i++) {
          gpfile << "\\\\$" << qca_data.elements[i] << "$=" << qca_data.conc_macro(qca_data.conc_macro.rows, i + 1);
        }
        gpfile << "}' at 0.0,0.0 center offset 0.0,-2.5" << endl;
        // END - custom
    }
    else if (qca_data.elements.size() == 2) {
      xvector<double> a = qca_data.conc_macro.getcol(1);
      for (int i = 1; i <= qca_data.conc_macro.rows; i++) {data.push_back({a(i), qca_data.binodal_boundary(i)});}
      // START - custom
      gpfile << "# Custom" << endl;
      gpfile << "set label '$" << qca_data.elements[0] << "$' at 0.0,0.0 center offset 0.0,-2.5" << endl;
      gpfile << "set label '$" << qca_data.elements[1] << "$' at 1.0,0.0 center offset 0.0,-2.5" << endl;
      // END - custom
    }
    else {
      return; // no plotting support for alloys higher than binary
    }
    plotter::generatePlotGNUPLOT(gpfile, plotoptions, data);
    if (qca_data.screen_only) {
      cout << gpfile.str() << endl;
    }
    else {
      plotter::savePlotGNUPLOT(plotoptions, gpfile);
    }
    return;
  }
}

#endif
