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

#include "aflow_qca.h"

// ###############################################################################
//            AFLOW Quasi-Chemical Approximation (QCA) (2022-)
// ###############################################################################

// ***************************************************************************
// qca::quasiChemicalApprox
// ***************************************************************************
namespace qca {
  /// @brief modifies the QCA variables based on the input flags and starts the module
  ///
  /// @param vpflow command line options
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void quasiChemicalApprox(const aurostd::xoption& vpflow) {
    if (vpflow.flag("QCA::USAGE")) {
      if (!vpflow.flag("QCA::SCREEN_ONLY")) {displayUsage();}
      return;
    }
    stringstream message;
    _aflags aflags;
    ofstream FileMESSAGE;
    ostream& p_oss = cout;
    _qca_data qca_data;
    initQCA(qca_data);
    qca_data.cdirpath = aurostd::getPWD();
    qca_data.rootdirpath = aurostd::CleanFileName(qca_data.rootdirpath);
    qca_data.plattice = aurostd::tolower(qca_data.plattice);
    qca_data.print = aurostd::tolower(qca_data.print);
    aurostd::sort_remove_duplicates(qca_data.elements);
    qca_data.min_sleep = DEFAULT_QCA_MIN_SLEEP_SECONDS;
    if (!vpflow.getattachedscheme("QCA::DIRECTORY").empty()) {
      qca_data.rootdirpath = vpflow.getattachedscheme("QCA::DIRECTORY");
    }
    else {
      qca_data.rootdirpath = aurostd::getPWD();
    }
    message << "Root directory = " << qca_data.rootdirpath;
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (!vpflow.getattachedscheme("QCA::AFLOWLIB_DIRECTORY").empty()) {qca_data.aflowlibpath = vpflow.getattachedscheme("QCA::AFLOWLIB_DIRECTORY");}
    message << "AFLOWLIB directory = " << qca_data.aflowlibpath;
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (!vpflow.getattachedscheme("QCA::AFLOW_MAX_NUM_ATOMS").empty()) {
      qca_data.aflow_max_num_atoms = aurostd::string2utype<int>(vpflow.getattachedscheme("QCA::AFLOW_MAX_NUM_ATOMS"));
    }
    else {
      qca_data.aflow_max_num_atoms = DEFAULT_QCA_AFLOW_MAX_NUM_ATOMS;
    }
    message << "Maximum number of atoms in AFLOW = " << qca_data.aflow_max_num_atoms;
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    qca_data.plattice = vpflow.getattachedscheme("QCA::PLATTICE");
    message << "Parent lattice = " << qca_data.plattice;
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (!vpflow.getattachedscheme("QCA::ELEMENTS").empty()) {
      aurostd::string2tokens(vpflow.getattachedscheme("QCA::ELEMENTS"), qca_data.elements, ",");
    }
    message << "Alloy elements = " << aurostd::joinWDelimiter(qca_data.elements, ",");
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (!vpflow.getattachedscheme("QCA::MAX_NUM_ATOMS").empty()) {
      qca_data.max_num_atoms = aurostd::string2utype<int>(vpflow.getattachedscheme("QCA::MAX_NUM_ATOMS"));
    }
    else {
      qca_data.max_num_atoms = DEFAULT_QCA_MAX_NUM_ATOMS;
    }
    message << "Maximum number of atoms in cluster expansion = " << qca_data.max_num_atoms;
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (!vpflow.getattachedscheme("QCA::CV_CUTOFF").empty()) {
      qca_data.cv_cut = aurostd::string2utype<double>(vpflow.getattachedscheme("QCA::CV_CUTOFF"));
    }
    else {
      qca_data.cv_cut = DEFAULT_QCA_CV_CUTOFF;
    }
    message << "Cross validation cutoff = " << qca_data.cv_cut;
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (!vpflow.getattachedscheme("QCA::CONC_CURVE_RANGE").empty()) {
      qca_data.conc_curve = true;
      aurostd::string2tokens(vpflow.getattachedscheme("QCA::CONC_CURVE_RANGE"), qca_data.conc_curve_range, ",");
    }
    message << "Concentration curve range = " << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < qca_data.conc_curve_range.size(); i++) {message << qca_data.conc_curve_range[i] << " ";}
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (!vpflow.getattachedscheme("QCA::CONC_NPTS").empty()) {
      qca_data.conc_npts = aurostd::string2utype<int>(vpflow.getattachedscheme("QCA::CONC_NPTS"));
    }
    else {
      qca_data.conc_npts = DEFAULT_QCA_CONC_NPTS;
    }
    message << "Number of points (concentration) = " << qca_data.conc_npts;
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (!vpflow.getattachedscheme("QCA::TEMP_RANGE").empty()) {
      aurostd::string2tokens(vpflow.getattachedscheme("QCA::TEMP_RANGE"), qca_data.temp_range, ",");
    }
    else {
      qca_data.temp_range = {DEFAULT_QCA_TEMP_MIN, DEFAULT_QCA_TEMP_MAX};
    }
    message << "Temperature range (K) = " <<  std::fixed << std::setprecision(2) 
            << qca_data.temp_range[0] << " " << qca_data.temp_range[1];
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (!vpflow.getattachedscheme("QCA::TEMP_NPTS").empty()) {
      qca_data.temp_npts = aurostd::string2utype<int>(vpflow.getattachedscheme("QCA::TEMP_NPTS"));
    }
    else {
      qca_data.temp_npts = DEFAULT_QCA_TEMP_NPTS;
    }
    message << "Number of points (temperature) = " << qca_data.temp_npts;
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (!vpflow.getattachedscheme("QCA::PRINT").empty()) {
      qca_data.print = vpflow.getattachedscheme("QCA::PRINT");
    }
    else {
      qca_data.print = DEFAULT_QCA_PRINT;
    }
    message << "Output format = " << qca_data.print;
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (vpflow.flag("QCA::SCREEN_ONLY")) {qca_data.screen_only = true;}
    if (vpflow.flag("QCA::IMAGE_ONLY")) {qca_data.image_only = true;}
    if (vpflow.flag("QCA::BINODAL")) {qca_data.calc_binodal = true;}
    if (vpflow.flag("QCA::USE_SG")) {qca_data.use_sg = true;}
    errorFix(qca_data);
    runQCA(qca_data);
   }
}

// ***************************************************************************
// qca::initQCA
// ***************************************************************************
namespace qca {
  /// @brief initializes the QCA variables
  ///
  /// @param qca_data object that contains all the QCA data needed to perform the calculation
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void initQCA(_qca_data& qca_data) {
    // Input data
    qca_data.min_sleep = 0;
    qca_data.print = "";
    qca_data.screen_only = false;
    qca_data.image_only = false;
    qca_data.calc_binodal = false;
    qca_data.use_sg = false;
    qca_data.cdirpath = "";
    qca_data.rootdirpath = "";
    qca_data.aflowlibpath = "";
    qca_data.plattice = "";
    qca_data.elements.clear();
    qca_data.aflow_max_num_atoms = 0;
    qca_data.max_num_atoms = 0;
    qca_data.cv_cut = 0.0;
    qca_data.conc_npts = 0;
    qca_data.conc_curve = false;
    qca_data.conc_curve_range.clear();
    qca_data.conc_macro.clear();
    qca_data.temp_npts = 0;
    qca_data.temp_range.clear();
    qca_data.temp.clear();
    // Derived data
    qca_data.alloyname = "";
    qca_data.rundirpath = "";
    qca_data.vstr_aflow.clear();
    qca_data.lat_atat = "";
    qca_data.vstr_atat.clear();
    qca_data.mapstr.clear();
    // Cluster data
    qca_data.cv_cluster = 0.0;
    qca_data.num_atom_cluster.clear();
    qca_data.degeneracy_cluster.clear();
    qca_data.conc_cluster.clear();
    qca_data.excess_energy_cluster.clear();
    // Thermo data
    qca_data.prob_ideal_cluster.clear();
    qca_data.prob_cluster.clear();
    qca_data.param_ec.first = 0.0; qca_data.param_ec.second = 0.0;
    qca_data.rel_s.clear();
    qca_data.binodal_boundary.clear();
  }
}

// ***************************************************************************
// qca::runQCA
// ***************************************************************************
namespace qca {
  /// @brief runs the QCA module
  ///
  /// @param qca_data object that contains all the QCA data needed to perform the calculation
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void runQCA(_qca_data& qca_data) {
    qca_data.rundirpath += qca_data.rootdirpath + "/" + pflow::arity_string(qca_data.elements.size(), false, false) + "/" + qca_data.plattice + "/" + qca_data.alloyname;
    aurostd::DirectoryMake(qca_data.rundirpath);
    // Only plot data from JSON file. This is useful when the plotting routine cannot be run on the machine
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
  }
}

// ***************************************************************************
// qca::errorFix
// ***************************************************************************
namespace qca {
  /// @brief checks and fixes errors in the input variables
  ///
  /// @param qca_data object that contains all the QCA data needed to perform the calculation
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void errorFix(_qca_data& qca_data) {
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
      vector<size_t> vzeros_init, vzeros_final;
      for (size_t i = 0; i < qca_data.elements.size() && totconc_init <= 1.0 && totconc_final <= 1.0; i++) {
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
      for (size_t i = 0; i < qca_data.conc_curve_range.size(); i++) {
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
    if (qca_data.print != "txt" && qca_data.print != "json" &&
        qca_data.print != "pdf" && qca_data.print != "eps" && qca_data.print != "png") {
      string message = "Format \"" + qca_data.print + "\" is invalid";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }
  }
}

// ***************************************************************************
// qca::calcBinodalData
// ***************************************************************************
namespace qca {
  /// @brief calculates the binodal curve and the quantities needed for the calculation
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  ///
  /// @see
  /// @doi{10.1016/j.actamat.2018.07.042}
  void calcBinodalData(_qca_data& qca_data) {
    if (aurostd::FileExist(qca_data.rundirpath + "/fit.out") && aurostd::FileExist(qca_data.rundirpath + "/predstr.out")) { // read ATAT data
      qca_data.vstr_atat = getATATXstructures(qca_data.lat_atat, qca_data.plattice, qca_data.elements, (uint)qca_data.max_num_atoms, qca_data.rundirpath);
    }
    else { // run ATAT
      qca_data.lat_atat = createLatForATAT(qca_data.plattice, qca_data.elements);
      qca_data.vstr_atat = getATATXstructures(qca_data.lat_atat, qca_data.plattice, qca_data.elements, (uint)qca_data.max_num_atoms);
      qca_data.vstr_aflow = getAFLOWXstructures(qca_data.plattice, qca_data.elements, qca_data.use_sg);
      if (!qca_data.aflowlibpath.empty()) { // add custom xstrs
        vector<xstructure> vstr_add = getAFLOWXstructures(qca_data.aflowlibpath, qca_data.use_sg);
        qca_data.vstr_aflow.insert(qca_data.vstr_aflow.end(), vstr_add.begin(), vstr_add.end());
      }
      qca_data.mapstr = calcMapForXstructures(getATATXstructures(qca_data.lat_atat, qca_data.plattice, qca_data.elements, qca_data.aflow_max_num_atoms), qca_data.vstr_aflow); // map ATAT xstrs to AFLOW xstrs because ATAT cannot identify AFLOW xstrs
      generateFilesForATAT(qca_data.rundirpath, createLatForATAT(qca_data.plattice, qca_data.elements, true), qca_data.vstr_aflow, qca_data.vstr_atat, qca_data.mapstr);
      runATAT(qca_data.cdirpath, qca_data.rundirpath, qca_data.min_sleep);
    }
    qca_data.cv_cluster = getCVCluster(qca_data.rundirpath, qca_data.cv_cut);
    qca_data.num_atom_cluster = getNumAtomCluster(qca_data.vstr_atat);
    qca_data.conc_cluster = getConcentrationCluster(qca_data.rundirpath, qca_data.vstr_atat.size(), qca_data.elements.size());
    qca_data.excess_energy_cluster = getExcessEnergyCluster(qca_data.rundirpath, qca_data.conc_cluster, qca_data.max_num_atoms);
    setCongruentClusters(qca_data);
    qca_data.degeneracy_cluster = calcDegeneracyCluster(qca_data.plattice, qca_data.vstr_atat, qca_data.elements, qca_data.max_num_atoms, qca_data.rundirpath);
    qca_data.conc_macro = getConcentrationMacro(qca_data.conc_curve_range, qca_data.conc_npts, qca_data.elements.size());
    qca_data.temp = getTemperatureRange(qca_data.temp_range, qca_data.temp_npts);
    qca_data.param_ec = calcRelativeEntropyEC(qca_data.conc_cluster, qca_data.degeneracy_cluster, qca_data.excess_energy_cluster, qca_data.temp, qca_data.max_num_atoms);
    if (qca_data.calc_binodal) {
      qca_data.prob_ideal_cluster = calcProbabilityIdealCluster(qca_data.conc_macro, qca_data.conc_cluster, qca_data.degeneracy_cluster, qca_data.max_num_atoms);
      checkProbability(qca_data.conc_macro, qca_data.conc_cluster, qca_data.prob_ideal_cluster);
      calcProbabilityCluster(qca_data.conc_macro, qca_data.conc_cluster, qca_data.excess_energy_cluster, qca_data.prob_ideal_cluster, qca_data.temp, qca_data.max_num_atoms, qca_data.prob_cluster);
      try {
        checkProbability(qca_data.conc_macro, qca_data.conc_cluster, qca_data.prob_ideal_cluster, qca_data.prob_cluster, qca_data.temp);
      }
      catch (aurostd::xerror& err) {
        return;
      }
      qca_data.rel_s = calcRelativeEntropy(qca_data.prob_cluster, qca_data.prob_ideal_cluster);
      qca_data.binodal_boundary = calcBinodalBoundary(qca_data.rel_s, qca_data.param_ec.first, qca_data.temp);
    }
    return;
  }
}

// ***************************************************************************
// qca::calcBinodalBoundary
// ***************************************************************************
namespace qca {
  /// @brief gets the binodal boundary
  ///
  /// @param rel_s relative entropy as a function of concentration and temperature
  /// @param rel_s_ec relative entropy at the equi-concentration and the transition temperature
  /// @param temp temperature range
  ///
  /// @return binodal boundary as a function of the concentration
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  xvector<double> calcBinodalBoundary(const xmatrix<double>& rel_s, const double rel_s_ec, const xvector<double>& temp) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    _aflags aflags;
    ofstream FileMESSAGE;
    ostream& p_oss = cout;
    int n_fit = 8;
    xvector<double> binodal_boundary(rel_s.rows), p, rr(n_fit), ri(n_fit);
    xvector<double> wts = aurostd::ones_xv<double>(rel_s.cols);
    double temp_mean = aurostd::mean(temp), temp_std = aurostd::stddev(temp);
    xvector<double> temp_scaled = (temp - temp_mean) / temp_std; // scale for numerical stability
    message << "Calculating binodal curve";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    for (int i = rel_s.lrows; i <= rel_s.urows; i++) {
      p = aurostd::polynomialCurveFit(temp_scaled, rel_s(i) - rel_s_ec, n_fit, wts);
      aurostd::polynomialFindRoots(p, rr, ri);
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " i=" << i << " | p=" << p << endl;
        cerr << __AFLOW_FUNC__ << "   Real roots=" << rr << endl;
        cerr << __AFLOW_FUNC__ << "   Imag roots=" << ri << endl;
      }
      rr = temp_std * rr + temp_mean;
      for (int j = 1; j <= n_fit; j++) {
        if (rr(j) >= temp(1) && rr(j) <= temp(temp.rows) && aurostd::isequal(ri(j), 0.0) && binodal_boundary(i) < rr(j)) { // largest solution must be real and within temp range
          binodal_boundary(i) = rr(j);
        }
      }
      pflow::updateProgressBar(i - rel_s.lrows + 1, rel_s.rows, p_oss);
    }
    return binodal_boundary;
  }
}

// ***************************************************************************
// qca::calcRelativeEntropy
// ***************************************************************************
namespace qca {
  /// @brief calculates the relative entropy.
  ///
  /// @param prob_cluster equilibrium probability of the clusters as a function of concentration and temperature.
  /// @param prob_cluster_ideal ideal (high-T) probability of the clusters as a function of concentration and temperature.
  ///
  /// @return relative entropy as a function of concentration and temperature.
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  xmatrix<double> calcRelativeEntropy(const vector<xmatrix<double>>& prob_cluster, const xmatrix<double>& prob_cluster_ideal) {
    xmatrix<double> rel_s(prob_cluster_ideal.rows, prob_cluster.size());
    for (int i = prob_cluster_ideal.lrows; i <= prob_cluster_ideal.urows; i++) {
      for (int j = 1; j <= (int)prob_cluster.size(); j++) {
        for (int k = prob_cluster_ideal.lcols; k <= prob_cluster_ideal.ucols; k++) {
          rel_s(i, j) += prob_cluster[j - 1](i, k) * aurostd::log((prob_cluster[j - 1](i, k) + _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_) / (prob_cluster_ideal(i, k) + _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_));
        }
      }
    }
    return rel_s;
  }
}

// ***************************************************************************
// qca::calcRelativeEntropyEC
// ***************************************************************************
namespace qca {
  /// @brief calculates the relative entropy evaluated at equi-concentration and the transition temperature
  ///
  /// @param conc_cluster concentration of the clusters
  /// @param degeneracy_cluster degeneracy of the clusters
  /// @param excess_energy_cluster excess energy of the clusters
  /// @param temp temperature range
  /// @param max_num_atoms maximum number of atoms in the cluster expansion
  /// @param interp interpolate order parameter to 0K
  ///
  /// @return relative entropy at the equi-concentration and the transition temperature
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  std::pair<double, double> calcRelativeEntropyEC(const xmatrix<double>& conc_cluster, const xvector<long int>& degeneracy_cluster, const xvector<double>& excess_energy_cluster, const xvector<double>& _temp, const int max_num_atoms, bool interp) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    _aflags aflags;
    ofstream FileMESSAGE;
    ostream& p_oss = cout;
    std::pair<double, double> param_ec(0.0, 0.0);
    int n_fit = 8;
    xvector<double> temp = _temp;
    xmatrix<double> conc_macro_ec(1, conc_cluster.cols);
    conc_macro_ec.set(1.0 / (double)conc_cluster.cols);
    xmatrix<double> prob_ideal_ec = calcProbabilityIdealCluster(conc_macro_ec, conc_cluster, degeneracy_cluster, max_num_atoms);
    vector<xmatrix<double>> prob_ec;
    // Check that the probability is physical, if not, shift the temperature range upward
    bool shift_temp = false;
    double dtemp = temp(2) - temp(1);
    message << "Checking for physical probabilites at equi-composition";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    while (!calcProbabilityCluster(conc_macro_ec, conc_cluster, excess_energy_cluster, prob_ideal_ec, temp, max_num_atoms, prob_ec) && 
           aurostd::min(temp) < DEFAULT_QCA_TEMP_MIN_LIMIT) {
      shift_temp = true;
      temp += dtemp;
    }
    message << "Equi-composition temperature range = [" << temp(1) << "K, " << temp(temp.urows) << "K]";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (shift_temp && aurostd::min(temp) > DEFAULT_QCA_TEMP_MIN_LIMIT) {
      message << "Could not find a temperature at equi-concentration that leads to a physical solution";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    xvector<double> order_param(prob_ec.size());
    xmatrix<double> m1, m2, m3 = prob_ideal_ec * aurostd::trasp(prob_ideal_ec);
    for (size_t i = 0; i < prob_ec.size(); i++) {
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
      message << "Order parameter diverges at low temperature, using interpolation";
      pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_WARNING_);
      xvector<double> dadt = aurostd::evalPolynomial_xv(temp_scaled, p1), temp_interp(temp.rows + 1), dadt_interp(temp.rows + 1);
      for (int i = temp.lrows; i <= temp.urows; i++) {temp_interp(i + 1) = temp(i);} // add T = 0K
      for (int i = dadt.lrows; i <= dadt.urows; i++) {dadt_interp(i + 1) = dadt(i);} // add Da/DT(0K) = 0
      temp_mean = aurostd::mean(temp_interp), temp_std = aurostd::stddev(temp_interp);
      temp_scaled = (temp_interp - temp_mean) / temp_std;
      wts = aurostd::exp(temp_scaled); // scale weights by temperature, since high-T is more accurate
      wts(1) = wts(wts.rows); // force convergence at T = 0K
      p1 = aurostd::polynomialCurveFit(temp_scaled, dadt_interp, n_fit - 1, wts);
      p2 = aurostd::evalPolynomialCoeff(p1, 1);
    }
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " alpha_orig=" << order_param << endl;
      cerr << __AFLOW_FUNC__ << " D[alpha, 0]=" << aurostd::evalPolynomial_xv(temp_scaled, p) << endl;
      cerr << __AFLOW_FUNC__ << " D[alpha, 1]=" << aurostd::evalPolynomial_xv(temp_scaled, p1) << endl;
      cerr << __AFLOW_FUNC__ << " D[alpha, 2]=" << aurostd::evalPolynomial_xv(temp_scaled, p2) << endl;
    }
    xvector<double> rr(n_fit - 2), ri(n_fit - 2);
    aurostd::polynomialFindRoots(p2, rr, ri);
    if (LDEBUG) {
      cerr << __AFLOW_FUNC__ << " p2=" << p2 << endl;
      cerr << __AFLOW_FUNC__ << "   Real roots=" << rr << endl;
      cerr << __AFLOW_FUNC__ << "   Imag roots=" << ri << endl;
    }
    bool unset = true;
    for (int j = 1; j <= n_fit; j++) {
      if (rr(j) >= temp_scaled(1) && rr(j) <= temp_scaled(temp.rows) && aurostd::isequal(ri(j), 0.0)) { // solution must be real and within temp range
        if (unset) {
          param_ec.second = rr(j);
          unset = false;
        }
        else if (aurostd::evalPolynomial(param_ec.second, p1) < aurostd::evalPolynomial(rr(j), p1)) { // keep largest gradient
          param_ec.second = rr(j);
        }
      }
    }
    param_ec.second = temp_std * param_ec.second + temp_mean;
    message << "Equi-composition transition temperature = " << param_ec.second << "K";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_COMPLETE_);
    calcProbabilityCluster(conc_macro_ec, conc_cluster, excess_energy_cluster, prob_ideal_ec, aurostd::vector2xvector(vector<double> {param_ec.second}), max_num_atoms, prob_ec);
    for (int i = prob_ideal_ec.lcols; i <= prob_ideal_ec.ucols; i++) {
      param_ec.first += prob_ec[0](1, i) * aurostd::log((prob_ec[0](1, i) + _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_) / (prob_ideal_ec(1, i) + _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_));
    }
    return param_ec;
  }
}

// ***************************************************************************
// qca::calcProbabilityCluster
// ***************************************************************************
namespace qca {
  /// @brief calculates the equilibrium cluster probability
  ///
  /// @param conc_macro macroscopic concentration of the alloy
  /// @param conc_cluster concentration of the clusters
  /// @param excess_energy_cluster excess energy of the clusters
  /// @param prob_cluster_ideal ideal (high-T) probability of the clusters as a function of concentration and temperature
  /// @param temp temperature range
  /// @param max_num_atoms maximum number of atoms in the cluster expansion
  /// @param prob_cluster equilibrium probability of the clusters as a function of concentration and temperature
  ///
  /// @return whether the calculation finished successfully
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  bool calcProbabilityCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<double>& excess_energy_cluster, const xmatrix<double>& prob_ideal_cluster, const xvector<double>& temp, const int max_num_atoms, vector<xmatrix<double>>& prob_cluster) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    _aflags aflags;
    ofstream FileMESSAGE;
    ostream& p_oss = cout;
    prob_cluster.clear();
    int neq = conc_cluster.cols - 1;
    xmatrix<double> zeros(prob_ideal_cluster.rows, prob_ideal_cluster.cols), natom_cluster = (double)max_num_atoms * conc_cluster;
    prob_cluster.assign(temp.rows, zeros); // initialize
    xvector<double> beta = aurostd::pow(KBOLTZEV * temp, -1.0), soln(neq);
    bool soln_found = false;
    message << "Calculating cluster probabilities";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    if (neq == 1) {
      xvector<double> p(max_num_atoms + 1), rr(max_num_atoms), ri(max_num_atoms);
      for (int it = temp.lrows; it <= temp.urows; it++) {
        for (int ix = prob_ideal_cluster.lrows; ix <= prob_ideal_cluster.urows; ix++) {
          p.reset();
          soln.reset();
          for (int j = prob_ideal_cluster.lcols; j <= prob_ideal_cluster.ucols; j++) {
            p((int)natom_cluster(j, 1) + 1) += prob_ideal_cluster(ix, j) * std::exp(-beta(it) * excess_energy_cluster(j)) * (conc_cluster(j, 1) - conc_macro(ix, 1));
          }
          if (LDEBUG) {cerr << __AFLOW_FUNC__ << " it=" << it << " ix=" << ix << " | p=" << p << endl;}
          aurostd::polynomialFindRoots(p, rr, ri);
          if (LDEBUG) {
            cerr << __AFLOW_FUNC__ << "   real roots=" << rr << endl;
            cerr << __AFLOW_FUNC__ << "   imag roots=" << ri << endl;
          }
          for (int k = 1; k <= max_num_atoms; k++) {
            if (rr(k) > soln(1) && 
                aurostd::isequal(ri(k), 0.0) && 
                prob_ideal_cluster.cols * std::pow(rr(k), max_num_atoms) != INFINITY) { // solution must be positive, real, and finite
              soln(1) = rr(k);
              soln_found = true;
            }
          }
          if (!soln_found) { // physical solution does not exist
            message << " Physical equilibrium probability does not exist for T=" << temp(it) << "K, X=[" << conc_macro(ix) << " ]" << endl;
            pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_WARNING_);
            return false;
          }
          for (int j = prob_ideal_cluster.lcols; j <= prob_ideal_cluster.ucols; j++) {
            prob_cluster[it - 1](ix, j) = prob_ideal_cluster(ix, j) * std::exp(-beta(it) * excess_energy_cluster(j)) * std::pow(soln(1), natom_cluster(j, 1));
          }
          prob_cluster[it - 1].setmat(prob_cluster[it - 1].getxmat(ix, ix, prob_ideal_cluster.lcols, prob_ideal_cluster.ucols) / 
                                      aurostd::sum(prob_cluster[it - 1].getxmat(ix, ix, prob_ideal_cluster.lcols, prob_ideal_cluster.ucols)), ix, 1); // normalize sum to 1
          pflow::updateProgressBar((it - temp.lrows) * prob_ideal_cluster.rows + (ix - prob_ideal_cluster.lrows) + 1, temp.rows * prob_ideal_cluster.rows, p_oss);
        }
      }
    }
    else {
      xmatrix<double> soln0 = aurostd::ones_xm<double>(neq, prob_ideal_cluster.rows), msoln;
      std::function<double(int it, int ix, int ieq, int ideq, xvector<double>)> 
        poly = [conc_macro, conc_cluster, excess_energy_cluster, prob_ideal_cluster, beta, natom_cluster]
               (int it, int ix, int ieq, int ideq, xvector<double> xvar) {
                 return calcProbabilityConstraint(conc_macro, conc_cluster, excess_energy_cluster, prob_ideal_cluster, beta, natom_cluster, it, ix, ieq, ideq, xvar);
               };
      vector<std::function<double(xvector<double>)>> vpoly, vdpoly;
      vector<vector<std::function<double(xvector<double>)>>> jac;
      for (int it = temp.urows; it >= temp.lrows; it--) { // go backwards
        for (int ix = prob_ideal_cluster.lrows; ix <= prob_ideal_cluster.urows; ix++) {
          vpoly.clear();
          vdpoly.clear();
          jac.clear();
          for (int ieq = 1; ieq <= neq; ieq++) {
            vpoly.push_back([poly, it, ix, ieq](xvector<double> xvar) {return poly(it, ix, ieq, 0, xvar);});
            for (int ideq = 1; ideq <= neq; ideq++) {
              vdpoly.push_back([poly, it, ix, ieq, ideq](xvector<double> xvar) {return poly(it, ix, ieq, ideq, xvar);});
            }
            jac.push_back(vdpoly);
            vdpoly.clear();
          }
          if (LDEBUG) {cerr << __AFLOW_FUNC__ << " it=" << it << " ix=" << ix << endl;}
          soln_found = aurostd::findZeroDeflation(soln0.getcol(ix), vpoly, jac, msoln);
          if (LDEBUG) {cerr << __AFLOW_FUNC__ << "   exist=" << soln_found << " | real roots=" << endl << msoln << endl;}
          if (soln_found) { // if solution is found, it must also be positive, real, and finite
            for (int isol = 1; isol <= msoln.cols; isol++) {
              soln_found = true;
              for (int ieq = 1; ieq <= neq && soln_found; ieq++) {
                soln_found = !(msoln.getcol(isol)(ieq) <= 0.0);
              }
              if (soln_found) {
                soln = msoln.getcol(isol);
                break;
              }
            }
          }
          if (!soln_found) {
            message << " Physical equilibrium probability does not exist for T=" << temp(it) << "K, X=[" << conc_macro(ix) << " ]" << endl;
            pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_WARNING_);
            return false;
          }
          soln0.setcol(soln, ix); // use the higher temperature solution as an initial guess for lower temperature solution
          for (int j = prob_ideal_cluster.lcols; j <= prob_ideal_cluster.ucols; j++) {
            prob_cluster[it - 1](ix, j) = prob_ideal_cluster(ix, j) * std::exp(-beta(it) * excess_energy_cluster(j)) * 
                                          aurostd::elements_product(aurostd::pow(soln, natom_cluster.getxvec(j, j, 1, neq)));
          }
          prob_cluster[it - 1].setmat(prob_cluster[it - 1].getxmat(ix, ix, prob_ideal_cluster.lcols, prob_ideal_cluster.ucols) / 
                                      aurostd::sum(prob_cluster[it - 1].getxmat(ix, ix, prob_ideal_cluster.lcols, prob_ideal_cluster.ucols)), ix, 1); // normalize sum to 1
          pflow::updateProgressBar((it - temp.lrows) * prob_ideal_cluster.rows + (ix - prob_ideal_cluster.lrows) + 1, temp.rows * prob_ideal_cluster.rows, p_oss);
        }
      }
    }
    return true;
  }
}

// ***************************************************************************
// qca::calcProbabilityConstraint
// ***************************************************************************
  /// @brief calculates the probability constraint for a particular set of indices
  ///
  /// @param conc_macro macroscopic concentration of the alloy
  /// @param conc_cluster concentration of the clusters
  /// @param excess_energy_cluster excess energy of the clusters
  /// @param prob_cluster_ideal ideal (high-T) probability of the clusters as a function of concentration and temperature
  /// @param beta inverse temperature range
  /// @param natom_cluster number of atoms in the clusters
  /// @param it temperature index
  /// @param ix macroscopic concentration index
  /// @param ik element index
  /// @param ideq index of the equation that is differentiated
  ///
  /// @return value of the constraint for a particular set of indices
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
namespace qca {
  double calcProbabilityConstraint(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<double>& excess_energy_cluster, const xmatrix<double>& prob_ideal_cluster, const xvector<double>& beta, const xmatrix<double>& natom_cluster, const int it, const int ix, const int ik, const int ideq, const xvector<double>& xvar) {
    double totsum = 0.0, prodx;
    int neq = conc_cluster.cols - 1;
    for (int j = prob_ideal_cluster.lcols; j <= prob_ideal_cluster.ucols; j++) {
      prodx = 1.0;
      for (int k = 1; k <= neq; k++) {
        if (k == ideq) {
          prodx *= natom_cluster(j, k) * std::pow(xvar(k), natom_cluster(j, k) - 1);
        }
        else {
          prodx *= std::pow(xvar(k), natom_cluster(j, k));
        }
      }
      totsum += prob_ideal_cluster(ix, j) * std::exp(-beta(it) * excess_energy_cluster(j)) * (conc_cluster(j, ik) - conc_macro(ix, ik)) * prodx;
    }
    return totsum;
  }
}

// ***************************************************************************
// qca::calcProbabilityIdealCluster
// ***************************************************************************
namespace qca {
  /// @brief calculates the ideal (high-T) cluster probability
  ///
  /// @param conc_macro macroscopic concentration of the alloy
  /// @param conc_cluster concentration of the clusters
  /// @param degeneracy_cluster degeneracy of the clusters
  /// @param max_num_atoms maximum number of atoms in the cluster expansion
  ///
  /// @return ideal (high-T) probability of the clusters as a function of concentration and temperature.
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  xmatrix<double> calcProbabilityIdealCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<long int>& degeneracy_cluster, const int max_num_atoms) {
    xmatrix<double> prob_ideal_cluster(conc_macro.rows, conc_cluster.rows);
    for (int i = conc_macro.lrows; i <= conc_macro.urows; i++) {
      for (int j = conc_cluster.lrows; j <= conc_cluster.urows; j++) {
        prob_ideal_cluster(i, j) = (double)degeneracy_cluster(j);
        for (int k = conc_macro.lcols; k <= conc_macro.ucols; k++) {
          prob_ideal_cluster(i, j) *= std::pow(conc_macro(i, k), conc_cluster(j, k) * max_num_atoms);
        }
      }
      prob_ideal_cluster.setmat(prob_ideal_cluster.getxmat(i, i, conc_cluster.lrows, conc_cluster.urows) /
                                aurostd::sum(prob_ideal_cluster.getxmat(i, i, conc_cluster.lrows, conc_cluster.urows)), i, 1); // normalize sum to 1
    }
    return prob_ideal_cluster;
  }
}

// ***************************************************************************
// qca::checkProbability
// ***************************************************************************
namespace qca {
  /// @brief checks the ideal (high-T) probability for errors
  ///
  /// @param conc_macro macroscopic concentration of the alloy
  /// @param conc_cluster concentration of the clusters
  /// @param prob_cluster_ideal ideal (high-T) probability of the clusters as a function of concentration and temperature
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void checkProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob_cluster_ideal) {
    for (int i = prob_cluster_ideal.lrows; i <= prob_cluster_ideal.urows; i++) {
      if (!aurostd::isequal(aurostd::sum(prob_cluster_ideal(i)), 1.0)) { // unnormalized
        stringstream message;
        message << "Ideal solution (high-T) probability is unnormalized for i=" << i << " | SUM[P_cluster]=" << aurostd::sum(prob_cluster_ideal(i));
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
      }
      else if (!aurostd::isequal(prob_cluster_ideal(i) * conc_cluster, conc_macro(i))) { // does not satisfy concentration constraints
        stringstream message;
        message << "Ideal solution (high-T) probability does not satisfy concentration contraint X=[" << conc_macro(i) << " ], X_calc=[" << prob_cluster_ideal(i) * conc_cluster << "]";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
      }
    }
  }

  /// @brief checks the equilibrium probability for errors
  ///
  /// @param conc_macro macroscopic concentration of the alloy
  /// @param conc_cluster concentration of the clusters
  /// @param prob_cluster_ideal ideal (high-T) probability of the clusters as a function of concentration and temperature
  /// @param prob_cluster equilibrium probability of the clusters as a function of concentration and temperature
  /// @param temp temperature range
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void checkProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob_ideal_cluster, const vector<xmatrix<double>>& prob_cluster, const xvector<double>& temp) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if (prob_cluster.empty()) {return;}
    double diff = AUROSTD_MAX_DOUBLE, diff_old;
    for (int it = 1; it <= (int)prob_cluster.size(); it++) {
      for (int i = prob_cluster[0].lrows; i <= prob_cluster[0].urows; i++) {
        if (!aurostd::isequal(aurostd::sum(prob_cluster[it - 1](i)), 1.0)) { // unnormalized
          stringstream message;
          message << "Equilibrium probability is unnormalized for T=" << temp(it) << "K, i=" << i << " | SUM[P_cluster]=" << aurostd::sum(prob_cluster[it - 1](i));
          throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
        }
        else if (!aurostd::isequal(prob_cluster[it - 1](i) * conc_cluster, conc_macro(i))) { // does not satisfy concentration constraints
          stringstream message;
          message << "Equilibrium probability does not satisfy concentration contraint for T=" << temp(it) << "K, X=[" << conc_macro(i) << " ], X_calc=[" << prob_cluster[it - 1](i) * conc_cluster << "]";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
        }
      }
      diff_old = diff;
      diff = aurostd::sum(aurostd::abs(prob_cluster[it - 1] - prob_ideal_cluster));
      if (LDEBUG) {cerr << __AFLOW_FUNC__ << " |P - P0| = " << diff << endl;}
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
  /// @brief gets the macroscopic concentration
  ///
  /// @param conc_curve_range concentration endpoints used to evaluate the macroscopic concentration
  /// @param conc_npts number of points used to evaluate the macroscopic concentration
  /// @param nelem number of elements in the alloy
  ///
  /// @return macroscopic concentration of the alloy
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
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
        for (int i = x.lrows; i <= x.rows; i++) {
          for (int j = x.lrows + 1; j <= x.rows; j++) {
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
// qca::getTemperatureRange
// ***************************************************************************
namespace qca {
  /// @brief gets the temperature profile
  ///
  /// @param temp_range temperature endpoints used to evaluate the temperature range
  /// @param temp_npts number of points used to evaluate the temperature range
  ///
  /// @return temperature range
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  xvector<double> getTemperatureRange(const vector<double>& temp_range, const int temp_npts) {
    return aurostd::linspace(temp_range[0], temp_range[1], temp_npts);
  }
}

// ***************************************************************************
// qca::setCongruentClusters
// ***************************************************************************
namespace qca {
  /// @brief prune the clusters to only include clusters congruent with the maximum number of atoms in the cluster expansion
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void setCongruentClusters(_qca_data& qca_data) {
    vector<int> index_cluster;
    for (int i = 1; i <= qca_data.num_atom_cluster.rows; i++) {
      if (!(qca_data.max_num_atoms % qca_data.num_atom_cluster(i))) {index_cluster.push_back(i);}
    }
    int ncl = index_cluster.size(), nelem = qca_data.elements.size();
    vector<xstructure> _vstr_atat(ncl);
    xvector<int> v1(ncl);
    xvector<double> v2(ncl);
    xmatrix<double> m1(ncl, nelem);
    for (size_t i = 0; i < index_cluster.size(); i++) {
      _vstr_atat[i] = qca_data.vstr_atat[index_cluster[i] - 1];
      m1.setmat(qca_data.conc_cluster.getxmat(index_cluster[i], index_cluster[i], 1, nelem), i + 1, 1);
      v1(i + 1) = qca_data.num_atom_cluster(index_cluster[i]);
      v2(i + 1) = qca_data.excess_energy_cluster(index_cluster[i]);
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
  /// @brief gets the cluster excess energy
  ///
  /// @param rundirpath path to the directory where AFLOW is running
  /// @param conc_cluster concentration of the clusters
  /// @param max_num_atoms maximum number of atoms in the cluster expansion
  ///
  /// @return excess energy of the clusters
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  xvector<double> getExcessEnergyCluster(const string& rundirpath, const xmatrix<double>& conc_cluster, const int max_num_atoms) {
    int ind, nstr = conc_cluster.rows, nelem = conc_cluster.cols;
    vector<string> vinput, tokens;
    aurostd::file2vectorstring(rundirpath + "/ref_energy.out", vinput);
    xvector<double> nrg_ref = aurostd::vector2xvector(aurostd::vectorstring2vectorutype<double>(vinput));
    xvector<double> nrg(nstr);
    // Read energies per atom
    aurostd::file2vectorstring(rundirpath + "/fit.out", vinput);
    aurostd::string2tokens(vinput[0], tokens, " ");
    if (nelem > (int)tokens.size()) {
      string message = "fit.out has the wrong output format";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    for (size_t line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      ind = aurostd::string2utype<int>(tokens[tokens.size() - 1]) + 1;
      nrg(ind) = aurostd::string2utype<double>(tokens[nelem]);
    }
    aurostd::file2vectorstring(rundirpath + "/predstr.out", vinput);
    aurostd::string2tokens(vinput[0], tokens, " ");
    if (nelem > (int)tokens.size()) {
      string message = "predstr.out has the wrong output format";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    for (size_t line = 0; line < vinput.size(); line++) {
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
  /// @brief gets the cluster concentration
  ///
  /// @param elements elements in the alloy
  /// @param vstr xstructures of the clusters
  ///
  /// @return concentration of the clusters
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  xmatrix<double> getConcentrationCluster(const vector<string>& elements, const vector<xstructure>& vstr) {
    uint nstr = vstr.size(), nelem = elements.size();
    xmatrix<double> conc_cluster(nstr, nelem);
    size_t ie;
    xvector<double> stoich(nelem);
    vector<string> str_elements;
    for (size_t i = 0; i < vstr.size(); i++) {
      if (nelem != vstr[i].stoich_each_type.size()) {
        str_elements = vstr[i].GetElements(true, true);
        stoich.reset();
        for (size_t j = 0; j < nelem; j++) {
          if (aurostd::WithinList(str_elements, elements[j], ie)) {stoich(j + 1) = vstr[i].stoich_each_type[ie];}
        }
      }
      else {
        stoich = aurostd::vector2xvector(aurostd::deque2vector(vstr[i].stoich_each_type));
      }
      for (size_t j = 0; j < nelem; j++) {conc_cluster(i + 1, j + 1) = stoich(j + 1);}
    }
    return conc_cluster;
  }

  /// @brief gets the cluster concentration
  ///
  /// @param rundirpath path to the directory where AFLOW is running
  /// @param nstr number of clusters
  /// @param nelem number of elements in the alloy
  ///
  /// @return concentration of the clusters
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
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
  /// @brief gets the cluster number of atoms
  ///
  /// @param vstr xstructures of the clusters
  ///
  /// @return number of atoms of the clusters
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  xvector<int> getNumAtomCluster(const vector<xstructure>& vstr) {
    int natom = 0;
    xvector<int> num_atom_cluster(vstr.size());
    for (size_t i = 0; i < vstr.size(); i++) {
      natom = 0;
      for (size_t j = 0; j < vstr[i].num_each_type.size(); j++) {natom += vstr[i].num_each_type[j];}
      num_atom_cluster(i + 1) = natom;
    }
    return num_atom_cluster;
  }
}

// ***************************************************************************
// qca::calcDegeneracyCluster
// ***************************************************************************
namespace qca {
  /// @brief calculates the cluster degeneracy
  ///
  /// @param plattice parent lattice of the alloy
  /// @param vstr xstructures of the clusters
  /// @param elements elements in the alloy
  /// @param max_num_atoms maximum number of atoms in the cluster expansion
  /// @param rundirpath path to the directory where AFLOW is running
  ///
  /// @return degeneracy of the clusters
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  xvector<long int> calcDegeneracyCluster(const string& plattice, const vector<xstructure>& _vstr, const vector<string>& elements, const int max_num_atoms, const string& rundirpath) {
    stringstream message;
    _aflags aflags;
    ofstream FileMESSAGE;
    ostream& p_oss = cout;
    string filepath = "";
    if (!rundirpath.empty()) {
      vector<string> tokens;
      aurostd::string2tokens(rundirpath, tokens, "/");
      filepath = "/" + QCA_FILE_PREFIX + "degen_" + aurostd::utype2string<int>(max_num_atoms) + "_atom.txt";
      for (int i = tokens.size() - 2; i >= 0; i--) {filepath = "/" + tokens[i] + filepath;} // well defined path due to how we define rundirpath
      if (aurostd::file2vectorstring(filepath, tokens)) {return aurostd::vector2xvector(aurostd::vectorstring2vectorutype<long int>(tokens));}
    }
    message << "Generating derivative structures";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    vector<xstructure> vstr = _vstr;
    xvector<long int> degeneracy_cluster(vstr.size());
    xstructure str_pocc, str;
    deque<_atom> patoms, atoms;
    _atom atom;
    uint itype;
    xmatrix<double> lattice(3,3);
    // Define primitive parent lattice and pocc atoms
    if (plattice == "fcc") {
      lattice(1, 1) = 0.0; lattice(1, 2) = 0.5; lattice(1, 3) = 0.5;
      lattice(2, 1) = 0.5; lattice(2, 2) = 0.0; lattice(2, 3) = 0.5;
      lattice(3, 1) = 0.5; lattice(3, 2) = 0.5; lattice(3, 3) = 0.0;
      atom.clear();
      atom.fpos(1) = 0.0; atom.fpos(2) = 0.0; atom.fpos(3) = 0.0;
      patoms.push_back(atom);
    }
    else if (plattice == "bcc") {
      lattice(1, 1) = -0.5; lattice(1, 2) = 0.5; lattice(1, 3) = 0.5;
      lattice(2, 1) = 0.5; lattice(2, 2) = -0.5; lattice(2, 3) = 0.5;
      lattice(3, 1) = 0.5; lattice(3, 2) = 0.5; lattice(3, 3) = -0.5;
      atom.clear();
      atom.fpos(1) = 0.0; atom.fpos(2) = 0.0; atom.fpos(3) = 0.0;
      patoms.push_back(atom);
    }
    else if (plattice == "hcp") {
      double sqrt3 = std::sqrt(3.0), sqrt8 = std::sqrt(8.0);
      lattice(1, 1) = 1.0; lattice(1, 2) = 0.0; lattice(1, 3) = 0.0;
      lattice(2, 1) = -0.5; lattice(2, 2) = 0.5 * sqrt3; lattice(2, 3) = 0.0;
      lattice(3, 1) = 0.0; lattice(3, 2) = 0.0; lattice(3, 3) = sqrt8 / sqrt3;
      atom.clear();
      atom.fpos(1) = 0.0; atom.fpos(2) = 0.0; atom.fpos(3) = 0.0;
      patoms.push_back(atom);
      atom.clear();
      atom.fpos(1) = 2.0 / 3.0; atom.fpos(2) = 1.0 / 3.0; atom.fpos(3) = 0.5;
      patoms.push_back(atom);
    }
    aurostd::xcombos xcpocc(max_num_atoms + 1, elements.size(), 'E', false);
    vector<xvector<double>> vpocc;
    while (xcpocc.increment()) {
      if (aurostd::sum(xcpocc.getCombo()) == max_num_atoms) {
        vpocc.push_back(aurostd::xvectorutype2xvectorvtype<int, double>(aurostd::vector2xvector(xcpocc.getCombo())) / (double)max_num_atoms);
      }
    }
    vector<xstructure> vstr_ds, vstr_u;
    long int sum_calculated = 0, sum_accepted = 0;
    for (size_t i = 0; i < vpocc.size(); i++) {
      str_pocc.clear();
      str_pocc.lattice = lattice;
      str_pocc.partial_occupation_HNF = max_num_atoms;
      str_pocc.partial_occupation_flag = true;
      str_pocc.neg_scale_second = true;
      atoms.clear();
      for (size_t j = 0; j < elements.size(); j++) {
        for (size_t k = 0; k < patoms.size(); k++) {
          atom.clear();
          atom.name = atom.cleanname = elements[j];
          atom.name_is_given = true;
          atom.partial_occupation_flag = true;
          atom.partial_occupation_value = vpocc[i](j + 1);
          atom.fpos = patoms[k].fpos;
          atom.cpos = str_pocc.f2c * patoms[k].fpos;
          atoms.push_back(atom);
        }
      }
      std::stable_sort(atoms.begin(), atoms.end(), sortAtomsNames);
      itype = 0;
      atoms[0].type = itype;
      for (size_t j = 1; j < atoms.size(); j++) {
        if (atoms[j].name != atoms[j - 1].name) {itype++;}
        atoms[j].type = itype;
      }
      str_pocc.AddAtom(atoms);
      // Enumerate unique superlattices
      pocc::POccCalculator pcalc(str_pocc);
      pcalc.calculateHNF();
      pcalc.getTotalPermutationsCount();
      pcalc.calculate();
      vstr_u = pcalc.getUniqueDerivativeStructures();
      for (vector<xstructure>::iterator it = vstr_u.begin(); it != vstr_u.end(); it++) {vstr_ds.push_back(*it);}
      sum_accepted += pcalc.total_permutations_count;
    }
    message << "Generated all unique derivative structures, count = " << vstr_ds.size();
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_COMPLETE_);
    // Find degenerate structures
    message << "Finding the degeneracy of clusters";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    // Shuffling the xstructures is important because it can avoid edge case scenarios where the reference centroid
    // in Xtalfinder does not pick up the proper neighbors
    vector<uint> index;
    for (size_t i = 0; i < vstr.size(); i++) {index.push_back(i);}
    aurostd::random_shuffle(index);
    for (size_t i = 0; i < index.size(); i++) {vstr[i] = _vstr[index[i]];}
    vstr_ds.insert(vstr_ds.begin(), vstr.begin(), vstr.end()); // concatenate xstructures
    XtalFinderCalculator xtal_calc;
    vector<vector<uint>> vindex = xtal_calc.groupSimilarXstructures(vstr_ds); // costly part of the function
    for (size_t i = 0; i < vindex.size(); i++) {
      std::sort(vindex[i].begin(), vindex[i].end()); // first index is the cluster index
      if (vindex[i][0] < vstr.size()) {
        for (size_t j = 1; j < vindex[i].size(); j++) {
          degeneracy_cluster(vindex[i][0] + 1) += pocc::getDGFromXStructureTitle(vstr_ds[vindex[i][j]].title);
        }
      }
    }
    xvector<long int> _degeneracy_cluster = degeneracy_cluster;
    for (size_t i = 0; i < index.size(); i++) {degeneracy_cluster(index[i] + 1) = _degeneracy_cluster(i + 1);} // reorder back to the original
    sum_calculated = aurostd::sum(degeneracy_cluster);
    if (!aurostd::isequal(sum_calculated, sum_accepted)) {
      stringstream message;
      message << "Degeneracies do not satisfy the sum rule, the sum is " << sum_calculated << " but should be " << sum_accepted;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
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
  /// @brief gets the coefficient of variation
  ///
  /// @param rundirpath path to the directory where AFLOW is running
  /// @param cv_cut coefficient of variation cut-off
  ///
  /// @return coefficient of variation of the cluster expansion
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  double getCVCluster(const string& rundirpath, const double cv_cut) {
    vector<string> vinput, tokens;
    aurostd::file2vectorstring(rundirpath + "/maps.log", vinput);
    aurostd::string2tokens(vinput.back(), tokens, " ");
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
  /// @brief runs the ATAT program
  ///
  /// @param cdirpath path to the directory where the original command was run
  /// @param rundirpath path to the directory where AFLOW is running
  /// @param min_sleep minimum number of seconds to sleep
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void runATAT(const string& cdirpath, const string& rundirpath, const uint min_sleep) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    _aflags aflags;
    ofstream FileMESSAGE;
    ostream& p_oss = cout;
    uint iter = 0;
    if (aurostd::substring2bool(aurostd::execute2string("mmaps", stdouterr_fsio), "command not found")) {
      message << "Missing mmaps program";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    aurostd::execute("touch " + rundirpath + "/stop"); // pre-kill mmaps gracefully
    aurostd::RemoveFile(rundirpath + "/maps_is_running");
    aurostd::RemoveFile(rundirpath + "/maps.log");
    chdir(rundirpath.c_str());
    string tmpfile = aurostd::TmpStrCreate(), logstring = "";
    aurostd::execute("mmaps -d > " + tmpfile + " 2>&1 &");
    message << "Running ATAT cluster expansion";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    while (!aurostd::substring2bool(logstring, "true and predicted ground states agree") && !aurostd::substring2bool(logstring, "No other ground states")) {
      iter++;
      if (LDEBUG) {cerr << __AFLOW_FUNC__ << " Sleeping, iter=" << iter << endl;}
      if (iter > 60) { // wait 60 times the minimum sleep
        message << "mmaps is taking too long to predict structures, directory = " << rundirpath;
        aurostd::RemoveFile(tmpfile);
        throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      aurostd::Sleep(min_sleep);
      aurostd::file2string(rundirpath + "/maps.log", logstring);
    }
    aurostd::RemoveFile(tmpfile);
    chdir(cdirpath.c_str());
  }
}

// ***************************************************************************
// qca::generateFilesForATAT
// ***************************************************************************
namespace qca {
  /// @brief generates the files for the ATAT program
  ///
  /// @param rundirpath path to the directory where AFLOW is running
  /// @param vstr_aflow xstructures from AFLOW runs
  /// @param vstr_atat xstructures from ATAT runs
  /// @param mapstr xstructure map between AFLOW and ATAT
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void generateFilesForATAT(const string& rundirpath, const string& lat_atat, const vector<xstructure>& vstr_aflow, const vector<xstructure>& vstr_atat, const vector<int>& mapstr) {
    stringstream oss;
    vector<xstructure> vstr;
    // Generate lat.in file
    aurostd::string2file(lat_atat, rundirpath + "/lat.in");
    // Generate str.out and energy files
    uint msize = mapstr.size();
    for (size_t i = 0; i < vstr_atat.size(); i++) {
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
  /// @brief gets the AFLOW xstructures from the AFLOW database
  ///
  /// @param plattice parent lattice of the alloy
  /// @param elements elements in the alloy
  /// @param use_sg compare initial and final xstructures only using their space groups
  ///
  /// @return xstructures from AFLOW runs
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  vector<xstructure> getAFLOWXstructures(const string& plattice, const vector<string>& elements, bool use_sg) {
    stringstream message;
    _aflags aflags;
    ofstream FileMESSAGE;
    ostream& p_oss = cout;
    vector<xstructure> vstr;
    vector<string> vstrlabel;
    string aflowlib, aflowurl;
    string alloyname = "";
    aflowlib::_aflowlib_entry entry;
    uint nelem = elements.size();
    stringstream oss;
    for (size_t i = 0; i < nelem; i++) {alloyname += AVASP_Get_PseudoPotential_PAW_PBE(elements[i]);}
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
    message << "Reading AFLOW xstructures from AFLOW database";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    for (size_t i = 0; i < vstrlabel.size(); i++) {
        entry.Load(aurostd::file2string(aflowlib + "/" + vstrlabel[i] + "/aflowlib.out"), oss);
        if (pflow::loadXstructures(entry, oss, false)) { // initial = unrelaxed; final = relaxed
          entry.aurl = aflowurl + "/" + vstrlabel[i];
          entry.vstr[0].qm_E_cell = entry.enthalpy_cell; // ATAT needs energy per cell
          if (use_sg && (entry.spacegroup_orig == entry.spacegroup_relax)) {
            vstr.push_back(entry.vstr[0]);
          }
          else if (compare::structuresMatch(entry.vstr[0], entry.vstr.back(), true, true, false)) {
            vstr.push_back(entry.vstr[0]);
          }
        }
        entry.clear();
        pflow::updateProgressBar(i, vstrlabel.size() - 1, p_oss);
    }
    return vstr;
  }

  /// @brief gets the AFLOW xstructures from a given directory
  ///
  /// @param aflowlibpath path to the parent directory where the aflowlib output files are stored
  /// @param use_sg compare initial and final xstructures only using their space groups
  ///
  /// @return xstructures from AFLOW runs
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  vector<xstructure> getAFLOWXstructures(const string& aflowlibpath, bool use_sg) {
    stringstream message;
    _aflags aflags;
    ofstream FileMESSAGE;
    ostream& p_oss = cout;
    vector<xstructure> vstr;
    vector<string> vdir;
    aurostd::SubDirectoryLS(aflowlibpath, vdir);
    string data;
    stringstream oss;
    aflowlib::_aflowlib_entry entry;
    message << "Reading AFLOW xstructures from directory = " << aflowlibpath;
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    for (size_t i = 0; i < vdir.size(); i++) {
      if (!aurostd::FileExist(vdir[i] + "/aflowlib.out") ||
          !aurostd::FileExist(vdir[i] + "/POSCAR.orig") ||
          !aurostd::FileExist(vdir[i] + "/POSCAR.relax2")) {
        continue;
      }
      entry.Load(aurostd::file2string(vdir[i] + "/aflowlib.out"), oss);
      oss.str(aurostd::file2string(vdir[i] + "/POSCAR.orig"));
      entry.vstr.push_back(xstructure(oss));
      oss.str(aurostd::file2string(vdir[i] + "/POSCAR.relax2"));
      entry.vstr.push_back(xstructure(oss));
      entry.vstr[0].qm_E_cell = entry.enthalpy_cell; // ATAT needs energy per cell
      if (use_sg && (entry.spacegroup_orig == entry.spacegroup_relax)) {
        vstr.push_back(entry.vstr[0]);
      }
      else if (compare::structuresMatch(entry.vstr[0], entry.vstr.back(), true, true, false)) {
        vstr.push_back(entry.vstr[0]);
      }
      entry.clear();
      pflow::updateProgressBar(i, vdir.size() - 1, p_oss);
    }
    return vstr;
  }
}

// ***************************************************************************
// qca::createLatForATAT
// ***************************************************************************
namespace qca {
  /// @brief gets the lattice file for ATAT
  ///
  /// @param plattice parent lattice of the alloy
  /// @param elements elements in the alloy
  ///
  /// @return lattice file for ATAT
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  string createLatForATAT(const string& plattice, const vector<string>& elements, bool scale) {
    stringstream oss;
    uint nelem = elements.size();
    double alat = 0.0;
    for (size_t i = 0; i < nelem; i++) {alat += GetAtomRadiusCovalent(elements[i]);}
    alat /= nelem;
    xmatrix<double> lattice(3,3);
    xvector<double> angles = 90.0 * aurostd::ones_xv<double>(3);
    xvector<double> coorsys = aurostd::ones_xv<double>(3);
    oss.precision(_DOUBLE_WRITE_PRECISION_MAX_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
    if (plattice == "fcc") {
      alat *= 2.0 * std::sqrt(2.0);
      lattice(1, 1) = 0.0; lattice(2, 1) = 0.5; lattice(3, 1) = 0.5;
      lattice(1, 2) = 0.5; lattice(2, 2) = 0.0; lattice(3, 2) = 0.5;
      lattice(1, 3) = 0.5; lattice(2, 3) = 0.5; lattice(3, 3) = 0.0;
    }
    else if (plattice == "bcc") {
      alat *= 4.0 / std::sqrt(3.0);
      lattice(1, 1) = -0.5; lattice(2, 1) = 0.5; lattice(3, 1) = 0.5;
      lattice(1, 2) = 0.5; lattice(2, 2) = -0.5; lattice(3, 2) = 0.5;
      lattice(1, 3) = 0.5; lattice(2, 3) = 0.5; lattice(3, 3) = -0.5;
    }
    else if (plattice == "hcp") {
      alat *= 2.0;
      lattice = aurostd::eye<double>(3, 3);
      coorsys(3) = std::sqrt(8.0 / 3.0);
      angles(3) = 120.0;
    }
    alat = (scale) ? aurostd::round(alat, 3) : 1.0;
    for (uint i = 1; i <= 3; i++) {oss << alat * coorsys(i) << " ";}
    for (uint i = 1; i <= 3; i++) {oss << angles(i) << " ";}
    oss << endl;
    for (uint i = 1; i <= 3; i++) {
      for (uint j = 1; j <= 3; j++) {
        oss << lattice(i, j) << " ";
      }
      oss << endl;
    }
    oss << 0.0 << " " << 0.0 << " " << 0.0 << " ";
    for (size_t i = 0; i < nelem; i++) {
      oss << elements[i] << ",";
    }
    if (plattice == "hcp") {
      oss << endl << 2.0 / 3.0 << " " << 1.0 / 3.0 << " " << 0.5 << " ";
      for (size_t i = 0; i < nelem; i++) {
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
  /// @brief gets the ATAT xstructures
  ///
  /// @param lat lattice file for ATAT
  /// @param max_num_atoms maximum number of atoms in the cluster expansion
  /// @param elements elements in the alloy
  /// @param rundirpath path to the directory where AFLOW is running
  ///
  /// @return xstructures from ATAT runs
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  vector<xstructure> getATATXstructures(const string& lat, const string& plattice, const vector<string>& elements, const uint max_num_atoms, const string& rundirpath) {
    stringstream message;
    _aflags aflags;
    ofstream FileMESSAGE;
    ostream& p_oss = cout;
    vector<xstructure> vstr;
    stringstream oss;
    if (!rundirpath.empty()) {
      message << "Reading ATAT xstructures from directory = " << rundirpath;
      pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
      vector<string> files;
      vector<xstructure> vstr_tmp;
      vector<uint> index;
      aurostd::DirectoryLS(rundirpath, files);
      for (size_t i = 0; i < files.size(); i++) {
        if (aurostd::FileExist(rundirpath + "/" + files[i] + "/str.out")) {
          aurostd::file2stringstream(rundirpath + "/" + files[i] + "/str.out", oss);
          vstr_tmp.push_back(xstructure(oss, IOATAT_STR));
          index.push_back(aurostd::string2utype<uint>(files[i]));
          aurostd::StringstreamClean(oss);
        }
        pflow::updateProgressBar(i, files.size() - 1, p_oss);
      }
      vstr.resize(index.size());
      for (size_t i = 0; i < index.size(); i++) {vstr[index[i]] = vstr_tmp[i];}
      return vstr;
    }
    uint nelem = elements.size();
    double scale = 0.0;
    for (size_t i = 0; i < nelem; i++) {scale += GetAtomRadiusCovalent(elements[i]);}
    scale /= nelem;
    if (plattice == "fcc") {
      scale *= 2.0 * std::sqrt(2.0);
    }
    else if (plattice == "bcc") {
      scale *= 4.0 / std::sqrt(3.0);
    }
    else if (plattice == "hcp") {
      scale *= 2.0;
    }
    scale = aurostd::round(scale, 3);
    vector<string> vinput, tokens;
    string tmpfile = aurostd::TmpStrCreate();
    aurostd::string2file(lat, tmpfile);
    string sstr = aurostd::execute2string("genstr -n " + aurostd::utype2string<uint>(max_num_atoms) + " -l " + tmpfile, stdouterr_fsio);
    aurostd::RemoveFile(tmpfile);
    if (sstr.size() == 0 || aurostd::substring2bool(sstr, "Unable to open lattice file")) {
      message << "Invalid lat.in file";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _FILE_CORRUPT_);
    }
    else if (aurostd::substring2bool(sstr, "command not found")) {
      message << "Missing genstr program";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
    aurostd::string2vectorstring(sstr, vinput);
    message << "Reading ATAT xstructures from genstr";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    for (size_t line = 0; line < vinput.size(); line++) {
      if (aurostd::substring2bool(vinput[line], "end")) {
        vstr.push_back(xstructure(oss, IOATAT_STR));
        vstr.back().scale = scale;
        vstr.back().ReScale(1.0);
        aurostd::StringstreamClean(oss);
      }
      else {
        oss << vinput[line] << endl;
      }
      pflow::updateProgressBar(line, vinput.size() - 1, p_oss);
    }
    return vstr;
  }
}

// ***************************************************************************
// qca::calcMapForXstructures
// ***************************************************************************
namespace qca {
  /// @brief calculates the map between two groups of xstructures
  ///
  /// @param vstr1 first group of xstructures
  /// @param vstr2 second group of xstructures
  ///
  /// @return xstructure map between the two groups
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  vector<int> calcMapForXstructures(const vector<xstructure>& _vstr1, const vector<xstructure>& vstr2) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    stringstream message;
    _aflags aflags;
    ofstream FileMESSAGE;
    ostream& p_oss = cout;
    vector<int> mapstr;
    vector<xstructure> vstr1 = _vstr1;
    vector<uint> index;
    message << "Mapping xstructures between AFLOW and ATAT";
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, FileMESSAGE, p_oss, _LOGGER_MESSAGE_);
    for (uint i = 0; i < vstr1.size(); i++) {index.push_back(i);}
    // Shuffling the xstructures is important because it can avoid edge case scenarios where the reference centroid
    // in Xtalfinder does not pick up the proper neighbors
    aurostd::random_shuffle(index); 
    for (size_t i = 0; i < index.size(); i++) {vstr1[i] = _vstr1[index[i]];}
    vector<xstructure> vstr3 = vstr1;
    vstr3.insert(vstr3.end(), vstr2.begin(), vstr2.end());
    XtalFinderCalculator xtal_calc;
    vector<vector<uint>> gsx = xtal_calc.groupSimilarXstructures(vstr3);
    for (size_t i = 0; i < gsx.size(); i++) {
      std::sort(gsx[i].begin(), gsx[i].end());
      if (gsx[i][0] > vstr1.size() -  1) { // index is not part of vstr1
        continue;
      }
      else if (gsx[i].size() == 1) { // no match for xstr in vstr1
        mapstr.push_back(-1);
      }
      else { // take first match of xstr in vstr2 to xstr in vstr1
        for (size_t j = 1; j < gsx[i].size(); j++) {
          if (gsx[i][j] > vstr1.size() - 1) {
            mapstr.push_back((int)gsx[i][j] - (int)vstr1.size());
            break;
          }
        }
      }
      if (LDEBUG) {
        cerr << __AFLOW_FUNC__ << " i=" << i << " | ";
        for (size_t j = 0; j < gsx[i].size(); j++) {cerr << gsx[i][j] << " ";}
        cerr << endl;
      }
    }
    if (mapstr.size() != index.size()) {
      message << "Something went wrong in the mapping of the xstructures; index.size()=" << index.size() << " | mapstr.size()=" << mapstr.size() << endl;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, _VALUE_ERROR_);
    }
    vector<int> _mapstr = mapstr;
    for (size_t i = 0; i < index.size(); i++) {mapstr[index[i]] = _mapstr[i];} // reorder back to the original
    return mapstr;
  }
}

// ***************************************************************************
// qca::displayUsage
// ***************************************************************************
namespace qca {
  /// @brief displays the usage commands and options
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void displayUsage(void) {
    vector<string> usage_options;
    usage_options.push_back("aflow --quasi_chem_approx|--qca --plattice=|--plat=fcc --elements=|--elem=Au,Pt[,Zn] [qca_options] [--directory=[DIRECTORY]]");
    usage_options.push_back(" ");
    usage_options.push_back("qca_options:");
    usage_options.push_back(" ");
    usage_options.push_back("GENERAL OPTIONS:");
    usage_options.push_back("--usage");
    usage_options.push_back("--screen_only");
    usage_options.push_back("--image_only|--image");
    usage_options.push_back("--aflowlib_directory=|--aflowlib_dir=...");
    usage_options.push_back("--print=|--p=|--output=|--o=txt");
    usage_options.push_back(" ");
    usage_options.push_back("BINODAL OPTIONS:");
    usage_options.push_back("--binodal");
    usage_options.push_back("--use_sg");
    usage_options.push_back("--aflow_max_num_atoms=4");
    usage_options.push_back("--max_num_atoms=|--mna=8");
    usage_options.push_back("--cv_cutoff=|--cv_cut=0.05");
    usage_options.push_back("--conc_curve_range=|--conc_curve=0,1,1,0");
    usage_options.push_back("--conc_npts=20");
    usage_options.push_back("--temp_range=|--temp=300,5000");
    usage_options.push_back("--temp_npts=150");
    usage_options.push_back(" ");
    init::MessageOption("--usage", "QCA()", usage_options);
  }
}

// ***************************************************************************
// qca::writeData
// ***************************************************************************
namespace qca {
  /// @brief write the QCA data to screen or to file
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void writeData(const _qca_data& qca_data) {
    string filepath = qca_data.rundirpath + "/" + QCA_FILE_PREFIX + "output." + qca_data.print;
    stringstream output;
    if (qca_data.print == "txt") {
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
      output << " " << info_prefix << "EC transition temperature (K)  = " << qca_data.param_ec.second << endl;
      output << " " << info_prefix << "Binodal boundary (K)           = " << endl << trasp(qca_data.binodal_boundary) << endl;
      if (!qca_data.screen_only) {
        aurostd::stringstream2file(output, filepath);
        return;
      }
    }
    else if (qca_data.print == "json") {
      aurostd::JSONwriter json;
      json.addString("Alloy name", qca_data.alloyname);
      json.addVector("Elements", qca_data.elements);
      json.addString("Parent lattice", qca_data.plattice);
      json.addNumber("Concentration curve", qca_data.conc_curve);
      json.addMatrix("Macroscopic concentration", qca_data.conc_macro);
      json.addVector("Temperature range (K)", qca_data.temp);
      json.addNumber("Max atoms per cell", qca_data.max_num_atoms);
      json.addNumber("Cluster CV score (eV)", qca_data.cv_cluster);
      json.addVector("Cluster number of atoms", aurostd::xvectorutype2xvectorvtype<int, double>(qca_data.num_atom_cluster));
      json.addVector("Cluster degeneracy", aurostd::xvectorutype2xvectorvtype<long int, double>(qca_data.degeneracy_cluster));
      json.addMatrix("Cluster concentration", qca_data.conc_cluster);
      json.addVector("Cluster excess energy (eV)", qca_data.excess_energy_cluster);
      json.addNumber("EC transition temperature (K)", qca_data.param_ec.second);
      json.addVector("Binodal boundary (K)", qca_data.binodal_boundary);
      aurostd::string2file(json.toString(), filepath);
      return;
    }
    cout << output.str() << endl;
  }
}

// ***************************************************************************
// qca::readData
// ***************************************************************************
namespace qca {
  /// @brief reads the QCA data
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
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
  /// @brief plots the QCA data
  ///
  /// @authors
  /// @mod{SD,20220718,created function}
  void plotData(const _qca_data& qca_data) {
    if (qca_data.print == "txt" || qca_data.print == "json") {return;}
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
    plotoptions.push_attached("IMAGE_FORMAT", qca_data.print);
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
      for (int i = qca_data.conc_macro.lrows; i <= qca_data.conc_macro.urows; i++) {data.push_back({a(i), qca_data.binodal_boundary(i)});}
        // START - custom
        gpfile << "# Custom" << endl;
        gpfile << "set label '\\shortstack{c}{$" << qca_data.elements[0] << "$=" << qca_data.conc_macro(1, 1);
        for (size_t i = 1; i < qca_data.elements.size(); i++) {
          gpfile << "\\\\$" << qca_data.elements[i] << "$=" << qca_data.conc_macro(1, i + 1);
        }
        gpfile << "}' at 0.0,0.0 center offset 0.0,-2.5" << endl;
        gpfile << "set label '\\shortstack{c}{$" << qca_data.elements[0] << "$=" << qca_data.conc_macro(qca_data.conc_macro.rows, 1);
        for (size_t i = 1; i < qca_data.elements.size(); i++) {
          gpfile << "\\\\$" << qca_data.elements[i] << "$=" << qca_data.conc_macro(qca_data.conc_macro.rows, i + 1);
        }
        gpfile << "}' at 0.0,0.0 center offset 0.0,-2.5" << endl;
        // END - custom
    }
    else if (qca_data.elements.size() == 2) {
      xvector<double> a = qca_data.conc_macro.getcol(1);
      for (int i = qca_data.conc_macro.lrows; i <= qca_data.conc_macro.urows; i++) {data.push_back({a(i), qca_data.binodal_boundary(i)});}
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
  }
}

#endif
