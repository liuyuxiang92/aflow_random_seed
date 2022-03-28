//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                  Simon Divilov - Duke University 2022                   *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov, 2022
// simon.divilov@duke.edu
// 
#ifndef _AFLOW_APDC_CPP_
#define _AFLOW_APDC_CPP_
#include "aflow.h"
#include "aflow_apdc.h"

// ###############################################################################
//            AFLOW Alloy Phase Diagram Constructor (APDC) (2022-)
// ###############################################################################

// **************************************************************************
// Class _apdc data
// **************************************************************************
// Constructor
_apdc_data::_apdc_data() {
  // Input data
  num_threads = 0;
  min_sleep = 0;
  workdirpath = "";
  rootdirpath = "";
  plattice = "";
  elements.clear();
  aflow_max_num_atoms = 0;
  max_num_atoms = 0;
  conc_npts = 0;
  conc_range.clear();
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
  mult_cluster.clear();
  natom_cluster.clear();
  conc_cluster.clear();
  excess_energy_cluster.clear();
  // Thermo data
  prob_ideal_cluster.clear();
  prob_cluster.clear();
}

// Destructor
_apdc_data::~_apdc_data() {
  free();
}
void _apdc_data::free() {
}

// Copy constructor
const _apdc_data& _apdc_data::operator=(const _apdc_data &b) {
  if (this != &b) {
    // Input data
    num_threads = b.num_threads;
    min_sleep = b.min_sleep;
    workdirpath = b.workdirpath;
    rootdirpath = b.rootdirpath;
    plattice = b.plattice;
    elements = b.elements;
    max_num_atoms = b.max_num_atoms;
    conc_npts = b.conc_npts;
    conc_range = b.conc_range;
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
    mult_cluster = b.mult_cluster;
    natom_cluster = b.natom_cluster;
    conc_cluster = b.conc_cluster;
    excess_energy_cluster = b.excess_energy_cluster;
    // Thermo data
    prob_ideal_cluster = b.prob_ideal_cluster;
    prob_cluster = b.prob_cluster;
  }
  return *this;
}

// ***************************************************************************
// apdc::GetPhaseDiagram
// ***************************************************************************
namespace apdc {
  void GetPhaseDiagram(_apdc_data& apdc_data) {
    // Clean-up input data and check for errors
    if (XHOST.vflag_control.flag("XPLUG_NUM_THREADS") && !(XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX"))) {
      apdc_data.num_threads = aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS"));
    }
    apdc_data.workdirpath = aurostd::getPWD();
    apdc_data.rootdirpath = aurostd::CleanFileName(apdc_data.rootdirpath);
    apdc_data.plattice = aurostd::tolower(apdc_data.plattice);
    aurostd::sort_remove_duplicates(apdc_data.elements);
    ErrorChecks(apdc_data);
    // Binodal
    apdc_data.rundirpath += apdc_data.rootdirpath + "/" + pflow::arity_string(apdc_data.elements.size(), false, false) + "/" + apdc_data.plattice + "/" + apdc_data.alloyname;
    aurostd::DirectoryMake(apdc_data.rundirpath);
    GetBinodal(apdc_data);
  }

 void GetPhaseDiagram(const string& aflowin, bool command_line_call) {
    _apdc_data apdc_data;
    string function_name = XPID + "GetPhaseDiagram():";
    if (command_line_call) {
      // FORMAT: <plattice>:<element1>,<element2>,...<element(K)>:<conc1_i>,<conc1_f>,<conc2_i><conc2_f>,..<conc(K),i><conc(K)_f>
      vector<string> tokens;
      aurostd::string2tokens(aflowin, tokens, ":");
      if (tokens.empty()) {
        string message = "Missing input";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _INPUT_ILLEGAL_);
      }
      if (tokens.size() != 3) {
        string message = "Invalid input";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _INPUT_ILLEGAL_);
      }
      apdc_data.min_sleep = DEFAULT_APDC_MIN_SLEEP_SECONDS;
      apdc_data.aflow_max_num_atoms = DEFAULT_APDC_AFLOW_MAX_NUM_ATOMS;
      apdc_data.max_num_atoms = DEFAULT_APDC_MAX_NUM_ATOMS;
      apdc_data.rootdirpath = aurostd::getPWD();
      apdc_data.plattice = tokens[0];
      aurostd::string2tokens(tokens[1], apdc_data.elements, ",");
      apdc_data.conc_npts = DEFAULT_APDC_CONC_NPTS;
      vector<double> vc;
      aurostd::string2tokens(tokens[2], vc, ",");
      apdc_data.conc_range = aurostd::vector2xvector(vc);
      apdc_data.temp_npts = DEFAULT_APDC_TEMP_NPTS;
      vector<double> vt = {DEFAULT_APDC_TEMP_MIN, DEFAULT_APDC_TEMP_MAX};
      apdc_data.temp_range = aurostd::vector2xvector(vt);
    }
    else {
      cerr<<aflowin<<endl;
    }
    GetPhaseDiagram(apdc_data);
 }

 void GetPhaseDiagram(istream& infile) {string aflowin; std::getline(infile, aflowin); GetPhaseDiagram(aflowin, false);}
}

// ***************************************************************************
// apdc::ErrorChecks
// ***************************************************************************
namespace apdc {
  void ErrorChecks(_apdc_data& apdc_data) {
    string function_name = XPID + "ErrorChecks():";
    // Check if number of threads is valid
    if (apdc_data.num_threads < 1) {apdc_data.num_threads = 1;}
    // Check if min sleep is at least 1 sec
    if (apdc_data.min_sleep < 1) {apdc_data.min_sleep = 1;}
    // Check if directory is writable
    if (!aurostd::DirectoryMake(apdc_data.rootdirpath)) {
      string message = "Cannot create directory";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _FILE_ERROR_);
    }
    // Check if parent lattice is valid
    if (apdc_data.plattice != "fcc" && apdc_data.plattice != "bcc" && apdc_data.plattice != "hcp") {
      string message = "Invalid parent lattice";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _FILE_ERROR_);
    }
    // Check if alloy is at least binary
    if (apdc_data.elements.size() < 2) {
      string message = "Alloy must be at least binary";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _FILE_ERROR_);
    }
    // Check if elements are valid, construct alloy name
    for (uint i = 0; i < apdc_data.elements.size(); i++) {
      if (!xelement::xelement::isElement(apdc_data.elements[i])) {
        string message = "Element \"" + apdc_data.elements[i] + "\" is invalid";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _INPUT_ILLEGAL_);
      }
      apdc_data.alloyname += apdc_data.elements[i];
    }
    // Check if concentration range format is valid
    if (apdc_data.conc_range.rows != 2 * (int)apdc_data.elements.size()) {
      string message = "Concentration range must have format [X1_start, X1_end, X2_start, X2_end,...X(K)_start, X(K)_end]";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _INPUT_ILLEGAL_);
    }
    // Check if concentration is within [0,1] and sums to 1
    double totconc_init = 0.0, totconc_final = 0.0;
    vector<double> vconc_init, vconc_final;
    vector<int> vzeros_init, vzeros_final;
    for (int i = 1; i <= (int)apdc_data.elements.size() && totconc_init <= 1.0 && totconc_final <= 1.0; i++) {
      if (apdc_data.conc_range(2 * i - 1) < 0 || apdc_data.conc_range(2 * i - 1) > 1 ||
          apdc_data.conc_range(2 * i) < 0 || apdc_data.conc_range(2 * i) > 1) {
        string message = "Concentration range must be defined on [0,1]";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _VALUE_ERROR_);
      }
      totconc_init += apdc_data.conc_range(2 * i - 1); totconc_final += apdc_data.conc_range(2 * i);
      vconc_init.push_back(apdc_data.conc_range(2 * i - 1)); vconc_final.push_back(apdc_data.conc_range(2 * i));
    }
    if (!aurostd::isequal(totconc_init, 1.0) || !aurostd::isequal(totconc_final, 1.0)) {
      string message = "Total concentration must sum to 1";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _VALUE_ERROR_);
    }
    // Redefine concentration range from [0,1] to (0,1)
    aurostd::WithinList(vconc_init, 0.0, vzeros_init); aurostd::WithinList(vconc_final, 0.0, vzeros_final);
    double cdelta_init = CONC_SHIFT * (double)vzeros_init.size() / ((double)vconc_init.size() - (double)vzeros_init.size());
    double cdelta_final = CONC_SHIFT * (double)vzeros_final.size() / ((double)vconc_final.size() - (double)vzeros_final.size());
    for (int i = 1; i <= apdc_data.conc_range.rows; i++) {
      if (aurostd::isequal(apdc_data.conc_range(i), 0.0)) {
        apdc_data.conc_range(i) += CONC_SHIFT;
      }
      else {
        apdc_data.conc_range(i) -= (i % 2) ? cdelta_init : cdelta_final;
      }
    }
    // Check if temperature range format is valid
    if (apdc_data.temp_range.rows != 2) {
      string message = "Temperature range must have format [T_start T_end]";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _INPUT_ILLEGAL_);
    }
    // Check if temperature values are valid
    if (apdc_data.temp_range(1) < 0 || apdc_data.temp_range(2) < 0) {
      string message = "Temperature cannot be below 0K";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _VALUE_ERROR_);
    }
  }
}

// ***************************************************************************
// apdc::GetBinodal
// ***************************************************************************
// Binodal construction based on method developed in Y. Lederer et al., Acta Materialia, 159 (2018)
namespace apdc {
  void GetBinodal(_apdc_data& apdc_data) {
    apdc_data.lat_atat = CreateLatForATAT(apdc_data.plattice, apdc_data.elements);
    apdc_data.vstr_atat = GetATATXstructures(apdc_data.lat_atat, (uint)apdc_data.max_num_atoms);
    apdc_data.vstr_aflow = GetAFLOWXstructures(apdc_data.plattice, apdc_data.elements, apdc_data.num_threads, false);
    //apdc_data.vstr_aflow = GetAFLOWXstructures(apdc_data.plattice, apdc_data.elements, apdc_data.num_threads);
    // map ATAT xstrs to AFLOW xstrs because ATAT cannot identify AFLOW xstrs
    apdc_data.mapstr = GetMapForXstructures(GetATATXstructures(apdc_data.lat_atat, apdc_data.aflow_max_num_atoms), apdc_data.vstr_aflow, apdc_data.num_threads);
    GenerateFilesForATAT(apdc_data.rundirpath, apdc_data.lat_atat, apdc_data.vstr_aflow, apdc_data.vstr_atat, apdc_data.mapstr);
    RunATAT(apdc_data.workdirpath, apdc_data.rundirpath, apdc_data.min_sleep);
    vector<xvector<int>> mult_cluster = GetMultiplicityCluster(apdc_data.vstr_atat);
    apdc_data.natom_cluster = mult_cluster[0];
    apdc_data.mult_cluster = mult_cluster[1];
    apdc_data.conc_cluster = GetConcentrationCluster(apdc_data.rundirpath, apdc_data.vstr_atat.size(), apdc_data.elements.size());
    apdc_data.excess_energy_cluster = GetExcessEnergyCluster(apdc_data.rundirpath, apdc_data.conc_cluster, apdc_data.natom_cluster);
    apdc_data.conc_macro = GetConcentrationMacro(apdc_data.conc_range, apdc_data.conc_npts, apdc_data.elements.size());
    SetCongruentClusters(apdc_data);
    apdc_data.temp = GetTemperature(apdc_data.temp_range, apdc_data.temp_npts);
    apdc_data.prob_ideal_cluster = GetProbabilityIdealCluster(apdc_data.conc_macro, apdc_data.conc_cluster, apdc_data.mult_cluster, apdc_data.max_num_atoms);
    CheckProbability(apdc_data.conc_macro, apdc_data.conc_cluster, apdc_data.prob_ideal_cluster);
    return;
    apdc_data.prob_cluster = GetProbabilityCluster(apdc_data.conc_macro, apdc_data.conc_cluster, apdc_data.excess_energy_cluster, apdc_data.prob_ideal_cluster, apdc_data.temp, apdc_data.max_num_atoms);
  }
}

// ***************************************************************************
// apdc::GetSpinodal
// ***************************************************************************
namespace apdc {
  void GetSpinodal(_apdc_data& apdc_data) {
    cerr << apdc_data.rundirpath << endl;
  }
}

// ***************************************************************************
// apdc::GetProbabilityCluster
// ***************************************************************************
namespace apdc {
  xtensor<double> GetProbabilityCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, xvector<double>& excess_energy_cluster, const xmatrix<double>& prob_ideal_cluster, const xvector<double>& temp, const int max_num_atoms) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    LDEBUG=TRUE;
    int nx = prob_ideal_cluster.rows, ncl = prob_ideal_cluster.cols, nt = temp.rows, neqs = conc_cluster.cols - 1;
    xtensor<double> prob_cluster({nx, ncl, nt});
    xvector<double> beta = aurostd::pow(CONSTANT_BOLTZMANN * temp, -1.0);
    xmatrix<int> natom_cluster = aurostd::xmatrixdouble2utype<int>((double)max_num_atoms * conc_cluster);
    xvector<double> rr, ri;
    if (neqs == 1) {
      for (int i = 1; i <= nx; i++) {
        for (int j = 1; j <= nt; j++) {
          xvector<double> coeff = 0.0 * aurostd::ones_xv<double>(max_num_atoms + 1);
          for (int k = 1; k <= ncl; k++) {
            coeff(natom_cluster(k, 1) + 1) += prob_ideal_cluster(i, k) * std::exp(-beta(j) * excess_energy_cluster(k)) * (conc_cluster(k, 1) - conc_macro(i, 1));
          }
          aurostd::polynomialFindRoots(coeff, rr, ri);
          if (LDEBUG) {
            cerr << "p=" << coeff << endl;
            cerr << "i=" << i << " j=" << j << " | Real roots=" << rr << endl;
            cerr << "i=" << i << " j=" << j << " | Imag roots=" << ri << endl;
          }
        }
      }
    }
    else { // homotopy continuation
    }
    return prob_cluster;
  }
}

// ***************************************************************************
// apdc::GetProbabilityIdealCluster
// ***************************************************************************
// P_j(X) = g_j*(X1^N1_j)*(X2^N2_j)*...(X(K)^N(K)_j)
namespace apdc {
  xmatrix<double> GetProbabilityIdealCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<int>& mult_cluster, const int max_num_atoms) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    int ncl = conc_cluster.rows, nx = conc_macro.rows, nelem = conc_macro.cols;
    xmatrix<double> prob_ideal_cluster(nx, ncl);
    for (int i = 1; i <= nx; i++) {
      for (int j = 1; j <= ncl; j++) {
        prob_ideal_cluster(i, j) = mult_cluster(j);
        for (int k = 1; k <= nelem; k++) {prob_ideal_cluster(i, j) *= std::pow(conc_macro(i, k), conc_cluster(j, k) * max_num_atoms);}
      }
      prob_ideal_cluster.setmat(prob_ideal_cluster.getmat(i, i, 1, ncl) / aurostd::sum(prob_ideal_cluster.getmat(i, i, 1, ncl)), i, 1); // normalize sum to 1
      if (LDEBUG) {cerr << "i=" << i << " | SUM[P_j]=" << aurostd::sum(prob_ideal_cluster.getmat(i, i, 1, ncl)) << endl;}
    }
    return prob_ideal_cluster;
  }
}

// ***************************************************************************
// apdc::CheckProbability
// ***************************************************************************
namespace apdc {
  bool CheckProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob) {
    int nx = prob.rows;
    for (int i = 1; i <= nx; i++) {
      if (!aurostd::isequal(aurostd::sum(prob(i)), 1.0)) {return false;}
      cerr<<prob(i)*conc_cluster<<endl;
      cerr<<conc_macro(i)<<endl;
      cerr<<"++++"<<endl;
    }
    return true;
  }
}

// ***************************************************************************
// apdc::GetConcentrationMacro
// ***************************************************************************
namespace apdc {
  xmatrix<double> GetConcentrationMacro(const xvector<double>& conc_range, const int conc_npts, const int nelem) {
    string function_name = XPID + "GetConcentrationMacro():";
    xmatrix<double> conc_macro(conc_npts, nelem);
    for (int i = 1; i <= nelem; i++) {
      conc_macro.setcol(aurostd::linspace(conc_range(2 * i - 1), conc_range(2 * i), conc_npts), i);
    }
    return conc_macro;
  }
}

// ***************************************************************************
// apdc::GetTemperature
// ***************************************************************************
namespace apdc {
  xvector<double> GetTemperature(const xvector<double>& temp_range, const int temp_npts) {
    return aurostd::linspace(temp_range(1), temp_range(2), temp_npts);
  }
}

// ***************************************************************************
// apdc::SetCongruentClusters
// ***************************************************************************
namespace apdc {
  void SetCongruentClusters(_apdc_data& apdc_data) {
    vector<int> indx_cluster;
    for (int i = 1; i <= apdc_data.natom_cluster.rows; i++) {
      if (!(apdc_data.max_num_atoms % apdc_data.natom_cluster(i))) {indx_cluster.push_back(i);}
    }
    int ncl = indx_cluster.size(), nelem = apdc_data.elements.size();
    xvector<int> v1(ncl), v2(ncl);
    xvector<double> v3(ncl);
    xmatrix<double> m1(ncl, nelem);
    for (int i = 0; i < ncl; i++) {
      m1.setmat(apdc_data.conc_cluster.getmat(indx_cluster[i], indx_cluster[i], 1, nelem), i + 1, 1);
      v1(i + 1) = apdc_data.mult_cluster(indx_cluster[i]);
      v2(i + 1) = apdc_data.natom_cluster(indx_cluster[i]);
      v3(i + 1) = apdc_data.excess_energy_cluster(indx_cluster[i]);
    }
    apdc_data.conc_cluster = m1;
    apdc_data.mult_cluster = v1;
    apdc_data.natom_cluster = v2;
    apdc_data.excess_energy_cluster = v3;
  }
}

// ***************************************************************************
// apdc::GetMultiplicityCluster
// ***************************************************************************
// N_j = N1_j+N2_j+...N(K)_j
// g_j = N_j!/(N1_j!*N2_j!*...N(K)_j)
// return [N,g]
namespace apdc {
  vector<xvector<int> > GetMultiplicityCluster(const vector<xstructure>& vstr) {
    int nstr = vstr.size();
    vector<xvector<int> > mult_cluster;
    xvector<int> g(nstr), n(nstr);
    uint natom, fact_prod;
    for (int i = 1; i <= nstr; i++) {
      natom = 0;
      fact_prod = 1;
      for (uint j = 0; j < vstr[i - 1].num_each_type.size(); j++) {
        natom += vstr[i - 1].num_each_type[j];
        fact_prod *= aurostd::factorial(vstr[i - 1].num_each_type[j]);
      }
      n(i) = natom;
      g(i) = aurostd::factorial(natom) / fact_prod;
    }
    mult_cluster.push_back(n); mult_cluster.push_back(g);
    return mult_cluster;
  }
}

// ***************************************************************************
// apdc::GetConcentrationCluster
// ***************************************************************************
namespace apdc {
  xmatrix<double> GetConcentration(const vector<string>& elements, const vector<xstructure>& vstr) {
    uint nstr = vstr.size(), nelem = elements.size();
    xmatrix<double> conc_cluster(nstr, nelem);
    int ie = -1;
    xvector<double> stoich;
    vector<string> str_elements;
    for (uint i = 0; i < vstr.size(); i++) {
      if (nelem != vstr[i].stoich_each_type.size()) {
        str_elements = vstr[i].GetElements(true, true);
        stoich = 0.0 * aurostd::ones_xv<double>(nelem);
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

  xmatrix<double> GetConcentrationCluster(const string& rundirpath, const int nstr, const int nelem) {
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
// apdc::GetExcessEnergyCluster
// ***************************************************************************
namespace apdc {
  xvector<double> GetExcessEnergyCluster(const string& rundirpath, const xmatrix<double>& conc_cluster, const xvector<int>& natom) {
    int nstr = conc_cluster.rows, nelem = conc_cluster.cols;
    vector<string> vinput, tokens;
    aurostd::file2vectorstring(rundirpath + "/ref_energy.out", vinput);
    xvector<double> nrg_ref = aurostd::vector2xvector(aurostd::vectorstring2vectordouble(vinput));
    xvector<double> nrg(nstr);
    aurostd::file2vectorstring(rundirpath + "/fit.out", vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      nrg(aurostd::string2utype<int>(tokens[tokens.size() - 1]) + 1) = aurostd::string2utype<double>(tokens[nelem]) * natom(aurostd::string2utype<int>(tokens[tokens.size() - 1]) + 1);
    }
    aurostd::file2vectorstring(rundirpath + "/predstr.out", vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      nrg(aurostd::string2utype<int>(tokens[tokens.size() - 2]) + 1) = aurostd::string2utype<double>(tokens[nelem]) * natom(aurostd::string2utype<int>(tokens[tokens.size() - 2]) + 1);
    }
    return nrg;
  }
}

// ***************************************************************************
// apdc::RunATAT
// ***************************************************************************
namespace apdc {
  void RunATAT(const string& workdirpath, const string& rundirpath, const uint min_sleep) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "RunATAT():";
    uint iter = 0;
    if (aurostd::substring2bool(aurostd::execute2string("mmaps", stdouterr_fsio), "command not found")) {
      string message = "Missing mmaps program";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _RUNTIME_ERROR_);
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
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _RUNTIME_ERROR_);
      }
      aurostd::Sleep(min_sleep);
      aurostd::file2string(rundirpath + "/maps.log", logstring);
    }
    aurostd::RemoveFile(tmpfile);
    chdir(workdirpath.c_str());
  }
}

// ***************************************************************************
// apdc::GenerateFilesForATAT
// ***************************************************************************
namespace apdc {
  void GenerateFilesForATAT(const string& rundirpath, const string& lat_atat, const vector<xstructure>& vstr_aflow, const vector<xstructure>& vstr_atat, const vector<int>& mapstr) {
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
// apdc::GetAFLOWXstructures
// ***************************************************************************
namespace apdc {
  vector<xstructure> GetAFLOWXstructures(const string& plattice, const vector<string>& elements, const int num_proc, bool use_sg) {
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
        vstrlabel.push_back("TFCC001.ABC");vstrlabel.push_back("TFCC002.ABC");vstrlabel.push_back("TFCC003.ABC");
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
        vstrlabel.push_back("TBCC001.ABC");vstrlabel.push_back("TBCC002.ABC");vstrlabel.push_back("TBCC003.ABC");
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
          else if (compare::structuresMatch(entry.vstr[0], entry.vstr[entry.vstr.size() - 1], true, num_proc)) {
            vstr.push_back(entry.vstr[0]);
          }
        }
        entry.clear();
    }
    return vstr;
  }
}

// ***************************************************************************
// apdc::CreateLatForATAT
// ***************************************************************************
namespace apdc {
  string CreateLatForATAT(const string& plattice, const vector<string>& elements) {
    stringstream oss;
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
    for (uint i = 0; i < elements.size(); i++) {
      oss << elements[i] << ",";
    }
    if (plattice == "hcp") {
      oss << endl << 2.0 / 3.0 << " " << 1.0 / 3.0 << " " << 0.5 << " ";
      for (uint i = 0; i < elements.size(); i++) {
        oss << elements[i] << ",";
      }
    }
    oss << endl;
    return oss.str();
  }
}

// ***************************************************************************
// apdc::GetATATXstructures
// ***************************************************************************
namespace apdc {
  vector<xstructure> GetATATXstructures(const string& lat, const uint max_num_atoms) {
    string function_name = XPID + "GetATATXstructures():";
    vector<xstructure> vstr;
    stringstream oss;
    vector<string> vinput, tokens;
    string tmpfile = aurostd::TmpStrCreate();
    aurostd::string2file(lat, tmpfile);
    string sstr = aurostd::execute2string("genstr -n " + aurostd::utype2string<uint>(max_num_atoms) + " -l " + tmpfile, stdouterr_fsio);
    aurostd::RemoveFile(tmpfile);
    if (sstr.size() == 0 || aurostd::substring2bool(sstr, "Unable to open lattice file")) {
      string message = "Invalid lat.in file";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _FILE_CORRUPT_);
    }
    else if (aurostd::substring2bool(sstr, "command not found")) {
      string message = "Missing genstr program";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, message, _RUNTIME_ERROR_);
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
// apdc::GetMapForXstructures
// ***************************************************************************
namespace apdc {
  vector<int> GetMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, const int num_proc) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    vector<int> mapstr;
    bool match;
    for (uint i = 0; i < vstr1.size(); i++) {
      match = false;
      for (uint j = 0; j < vstr2.size() && !match; j++) {
        if (compare::structuresMatch(vstr1[i], vstr2[j], true, num_proc)) {
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






#endif
