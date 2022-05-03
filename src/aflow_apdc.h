//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                  Simon Divilov - Duke University 2022                   *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov, 2022
// simon.divilov@duke.edu
// 
#ifndef _AFLOW_APDC_H_
#define _AFLOW_APDC_H_

#define DEFAULT_ATAT_ALAT 4.0 // Ang
#define CONC_SHIFT 0.01 // concentration shift away from 0
#define APDC_FILE_PREFIX string("apdc_")

// Class _apdc data
class _apdc_data {
  public:
    _apdc_data();
    ~_apdc_data();
    const _apdc_data& operator=(const _apdc_data &b);
    
    // Input data
    int num_threads;
    uint min_sleep;
    string format_data;
    string format_plot;
    bool image_only;
    string workdirpath;
    string rootdirpath;
    string plattice;
    vector<string> elements;
    int aflow_max_num_atoms;
    int max_num_atoms;
    int conc_npts;
    vector<double> conc_curve; // DIM: 2*Nk
    xmatrix<double> conc_macro; // UNIT: unitless | DIM: Nc, Nk
    int temp_npts;
    vector<double> temp_range; // UNIT: K | DIM: 2
    xvector<double> temp; // UNIT: K | DIM: Nt

    // Derived data
    string alloyname;
    string rundirpath;
    vector<xstructure> vstr_aflow;
    string lat_atat;
    vector<xstructure> vstr_atat; // DIM: Nj
    vector<int> mapstr;

    // Cluster data
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

// Namespace for functions used by APDC
namespace apdc {
  void GetPhaseDiagram(const aurostd::xoption& vpflow);
  void GetPhaseDiagram(_apdc_data& apdc_data);
  void ErrorChecks(_apdc_data& apdc_data);
  void GetSpinodalData(_apdc_data& apdc_data);
  void GetBinodalData(_apdc_data& apdc_data);
  xvector<double> GetBinodalBoundary(const xmatrix<double>& rel_s, const double& rel_s_ec, const xvector<double>& temp);
  xmatrix<double> GetRelativeEntropy(const vector<xmatrix<double>>& prob_cluster, const xmatrix<double>& prob_cluster_ideal);
  vector<double> GetRelativeEntropyEC(const xmatrix<double>& conc_cluster, const xvector<int>& degeneracy_cluster, const xvector<double>& excess_energy_cluster, const xvector<double>& temp, const int max_num_atoms);
  vector<xmatrix<double>> GetProbabilityCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<double>& excess_energy_cluster, const xmatrix<double>& prob_ideal_cluster, const xvector<double>& temp, const int max_num_atoms);
  xmatrix<double> GetProbabilityIdealCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<int>& degeneracy_cluster, const int max_num_atoms);
  void CheckProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob);
  void CheckProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob0, const vector<xmatrix<double>>& prob);
  xmatrix<double> GetConcentrationMacro(const vector<double>& conc_curve, const int conc_npts, const uint nelem);
  xvector<double> GetTemperature(const vector<double>& temp_range, const int temp_npts);
  void SetCongruentClusters(_apdc_data& apdc_data);
  xvector<int> GetNumAtomCluster(const vector<xstructure>& vstr);
  xvector<int> GetDegeneracyCluster(const string& plattice, const vector<xstructure>& vstr, const vector<string>& elements, const int max_num_atoms, const bool shuffle, const string& rundirpath="");
  xmatrix<double> GetConcentrationCluster(const vector<string>& elements, const vector<xstructure>& vstr);
  xmatrix<double> GetConcentrationCluster(const string& rundirpath, const int nstr, const int nelem);
  xvector<double> GetExcessEnergyCluster(const string& rundirpath, const xmatrix<double>& conc_cluster, const xvector<int>& natom);
  void RunATAT(const string& workdirpath, const string& rundirpath, const uint min_sleep);
  void GenerateFilesForATAT(const string& rundirpath, const string& lat_atat, const vector<xstructure>& vstr_aflow, const vector<xstructure>& vstr_atat, const vector<int>& mapstr);
  vector<xstructure> GetAFLOWXstructures(const string& plattice, const vector<string>& elements, const int num_threads, bool use_sg=true);
  string CreateLatForATAT(const string& plattice, const vector<string>& elements);
  vector<xstructure> GetATATXstructures(const string& lat, const uint max_num_atoms, const string& rundirpath="");
  vector<int> GetMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, const int num_threads);
  void DisplayUsage(void);
  void PrintData(const _apdc_data& apdc_data);
  void PlotData(const _apdc_data& apdc_data);
  void PlotData(const string& rundirpath);
}

#endif
