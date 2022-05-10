//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                  Simon Divilov - Duke University 2022                   *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov, 2022
// simon.divilov@duke.edu
// 
#ifndef _AFLOW_APEC_H_
#define _AFLOW_APEC_H_

#define DEFAULT_ATAT_ALAT 4.0 // Ang
#define CONC_SHIFT 0.01 // concentration shift away from 0
#define APEC_FILE_PREFIX string("apec_")

// Class _apec data
class _apec_data {
  public:
    _apec_data();
    ~_apec_data();
    const _apec_data& operator=(const _apec_data &b);
    
    // Input data
    int num_threads;
    uint min_sleep;
    string format_data;
    string format_image;
    bool screen_only;
    bool image_only;
    string workdirpath;
    string rootdirpath;
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
    double rel_s_ec; // UNIT: unitless
    double temp_ec; // UNIT: K
    xmatrix<double> rel_s; // UNIT: unitless | DIM: Nc, Nt
    xvector<double> binodal_boundary; // UNIT: K | DIM: Nc
    
  private:
  void free();
};

// Namespace for functions used by APEC
namespace apec {
  void getPhaseDiagram(const aurostd::xoption& vpflow);
  void getPhaseDiagram(_apec_data& apec_data);
  void errorChecks(_apec_data& apec_data);
  void getSpinodalData(_apec_data& apec_data);
  void getBinodalData(_apec_data& apec_data);
  xvector<double> getBinodalBoundary(const xmatrix<double>& rel_s, const double rel_s_ec, const xvector<double>& temp);
  xmatrix<double> getRelativeEntropy(const vector<xmatrix<double>>& prob_cluster, const xmatrix<double>& prob_cluster_ideal);
  vector<double> getRelativeEntropyEC(const xmatrix<double>& conc_cluster, const xvector<int>& degeneracy_cluster, const xvector<double>& excess_energy_cluster, const xvector<double>& temp, const int max_num_atoms);
  vector<xmatrix<double>> getProbabilityCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<double>& excess_energy_cluster, const xmatrix<double>& prob_ideal_cluster, const xvector<double>& temp, const int max_num_atoms);
  xmatrix<double> getProbabilityIdealCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<int>& degeneracy_cluster, const int max_num_atoms);
  void checkProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob);
  void checkProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob0, const vector<xmatrix<double>>& prob);
  xmatrix<double> getConcentrationMacro(const vector<double>& conc_curve_range, const int conc_npts, const uint nelem);
  xvector<double> getTemperature(const vector<double>& temp_range, const int temp_npts);
  void setCongruentClusters(_apec_data& apec_data);
  xvector<double> getExcessEnergyCluster(const string& rundirpath, const xmatrix<double>& conc_cluster, const int max_num_atoms);
  xmatrix<double> getConcentrationCluster(const vector<string>& elements, const vector<xstructure>& vstr);
  xmatrix<double> getConcentrationCluster(const string& rundirpath, const int nstr, const int nelem);
  xvector<int> getNumAtomCluster(const vector<xstructure>& vstr);
  xvector<int> getDegeneracyCluster(const string& plattice, const vector<xstructure>& vstr, const vector<string>& elements, const int max_num_atoms, const bool shuffle, const string& rundirpath="");
  double getCVCluster(const string& rundirpath, const double cv_cut);
  void runATAT(const string& workdirpath, const string& rundirpath, const uint min_sleep);
  void generateFilesForATAT(const string& rundirpath, const string& lat_atat, const vector<xstructure>& vstr_aflow, const vector<xstructure>& vstr_atat, const vector<int>& mapstr);
  vector<xstructure> getAFLOWXstructures(const string& plattice, const vector<string>& elements, const int num_threads, bool use_sg=false);
  string createLatForATAT(const string& plattice, const vector<string>& elements);
  vector<xstructure> getATATXstructures(const string& lat, const uint max_num_atoms, const string& rundirpath="");
  vector<int> getMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, const int num_threads);
  void displayUsage(void);
  void writeData(const _apec_data& apec_data);
  void readData(_apec_data& apec_data);
  void plotData(const _apec_data& apec_data);
}

#endif
