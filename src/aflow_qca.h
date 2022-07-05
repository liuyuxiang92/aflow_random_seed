//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                  Simon Divilov - Duke University 2022                   *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov, 2022
// simon.divilov@duke.edu
// 
#ifndef _AFLOW_QCA_H_
#define _AFLOW_QCA_H_

#define CONC_SHIFT 0.01 // concentration shift away from 0
#define QCA_FILE_PREFIX string("aflow_qca_")

// Struct _qca data
struct _qca_data {
  // Input data
  int num_threads;
  uint min_sleep;
  string format_data;
  string format_image;
  bool screen_only;
  bool image_only;
  bool calc_binodal;
  bool use_sg;
  string cdirpath;
  string rootdirpath;
  string aflowlibpath;
  string plattice;
  vector<string> elements;
  int aflow_max_num_atoms;
  int max_num_atoms;
  int conc_npts;
  bool conc_curve;
  vector<double> conc_curve_range; // DIM: 2*Nk
  xmatrix<double> conc_macro; // UNIT: unitless | DIM: Nc, Nk
  int temp_npts;
  vector<double> temp_range; // UNIT: K | DIM: 2
  xvector<double> temp; // UNIT: K | DIM: Nt
  double cv_cut; // UNIT: eV

  // Derived data
  string alloyname;
  string rundirpath;
  vector<xstructure> vstr_aflow;
  string lat_atat;
  vector<xstructure> vstr_atat; // DIM: Nj
  vector<int> mapstr;

  // Cluster data
  double cv_cluster; // UNIT: eV
  xvector<int> num_atom_cluster; // DIM: Nj
  xvector<int> degeneracy_cluster; // DIM: Nj
  xmatrix<double> conc_cluster; // UNIT: unitless | DIM: Nj, Nk
  xvector<double> excess_energy_cluster; // UNIT: eV | DIM: Nj

  // Thermo data
  xmatrix<double> prob_ideal_cluster; // DIM: Nc, Nj
  vector<xmatrix<double>> prob_cluster; // DIM: Nc, Nj, Nt
  std::pair<double, double> param_ec; // UNIT: unitless, K
  xmatrix<double> rel_s; // UNIT: unitless | DIM: Nc, Nt
  xvector<double> binodal_boundary; // UNIT: K | DIM: Nc
};

// Namespace for functions used by QCA
namespace qca {
  void quasiChemicalApprox(const aurostd::xoption& vpflow);
  void initQCA(_qca_data& qca_data);
  void runQCA(_qca_data& qca_data);
  void errorFix(_qca_data& qca_data);
  void calcSpinodalData(_qca_data& qca_data);
  void calcBinodalData(_qca_data& qca_data);
  xvector<double> calcBinodalBoundary(const xmatrix<double>& rel_s, const double rel_s_ec, const xvector<double>& temp);
  xmatrix<double> calcRelativeEntropy(const vector<xmatrix<double>>& prob_cluster, const xmatrix<double>& prob_cluster_ideal);
  std::pair<double, double> calcRelativeEntropyEC(const xmatrix<double>& conc_cluster, const xvector<int>& degeneracy_cluster, const xvector<double>& excess_energy_cluster, const xvector<double>& temp, const int max_num_atoms, bool interp=true);
  bool calcProbabilityCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<double>& excess_energy_cluster, const xmatrix<double>& prob_ideal_cluster, const xvector<double>& temp, const int max_num_atoms, vector<xmatrix<double>>& prob_cluster);
  double calcProbabilityConstraint(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<double>& excess_energy_cluster, const xmatrix<double>& prob_ideal_cluster, const xvector<double>& beta, const xmatrix<double>& natom_cluster, const int it, const int ix, const int ik, const int ideq, const xvector<double>& xvar);
  xmatrix<double> calcProbabilityIdealCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<int>& degeneracy_cluster, const int max_num_atoms);
  void checkProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob_ideal_cluster);
  void checkProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob_ideal_cluster, const vector<xmatrix<double>>& prob_cluster, const xvector<double>& temp);
  xmatrix<double> getConcentrationMacro(const vector<double>& conc_curve_range, const int conc_npts, const uint nelem);
  xvector<double> getTemperatureRange(const vector<double>& temp_range, const int temp_npts);
  void setCongruentClusters(_qca_data& qca_data);
  xvector<double> getExcessEnergyCluster(const string& rundirpath, const xmatrix<double>& conc_cluster, const int max_num_atoms);
  xmatrix<double> getConcentrationCluster(const vector<xstructure>& vstr, const vector<string>& elements);
  xmatrix<double> getConcentrationCluster(const string& rundirpath, const int nstr, const int nelem);
  xvector<int> getNumAtomCluster(const vector<xstructure>& vstr);
  xvector<int> calcDegeneracyCluster(const string& plattice, const vector<xstructure>& vstr, const vector<string>& elements, const int max_num_atoms, const int num_threads, const string& rundirpath="", const string& algo="SLOW");
  double getCVCluster(const string& rundirpath, const double cv_cut);
  void runATAT(const string& cdirpath, const string& rundirpath, const uint min_sleep);
  void generateFilesForATAT(const string& rundirpath, const string& lat_atat, const vector<xstructure>& vstr_aflow, const vector<xstructure>& vstr_atat, const vector<int>& mapstr);
  vector<xstructure> getAFLOWXstructures(const string& plattice, const vector<string>& elements, const int num_threads, bool use_sg);
  vector<xstructure> getAFLOWXstructures(const string& aflowlibpath, const int num_threads, bool use_sg);
  string createLatForATAT(const string& plattice, const vector<string>& elements, bool scale=false);
  vector<xstructure> getATATXstructures(const string& lat, const string& plattice, const vector<string>& elements, const uint max_num_atoms, const string& rundirpath="");
  vector<int> calcMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, const int num_threads, const string& algo="SLOW");
  void displayUsage(void);
  void writeData(const _qca_data& qca_data);
  void readData(_qca_data& qca_data);
  void plotData(const _qca_data& qca_data);
}

#endif
