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
#define CONSTANT_BOLTZMANN 8.6173324e-5 // eV/K
#define CONC_DELTA 0.01
#define _APDC_STR_OPT_ string("[AFLOW_APDC]")

// Class _apdc data
class _apdc_data {
  public:
    _apdc_data();
    ~_apdc_data();
    const _apdc_data& operator=(const _apdc_data &b);
    
    // Input data
    int num_threads;
    uint min_sleep;
    string workdirpath;
    string rootdirpath;
    string plattice;
    vector<string> elements;
    int aflow_max_num_atoms;
    int max_num_atoms;
    int conc_npts;
    xvector<double> conc_range; // DIM: 2*K
    xmatrix<double> conc_macro; // UNIT: unitless | DIM: Nc, K
    int temp_npts;
    xvector<double> temp_range; // DIM: 2
    xvector<double> temp; // UNIT: K | DIM: Nt

    // Derived data
    string alloyname;
    string rundirpath;
    vector<xstructure> vstr_aflow;
    string lat_atat;
    vector<xstructure> vstr_atat; // DIM: J
    vector<int> mapstr;

    // Cluster data
    xvector<int> multiplicity; // DIM: J
    xmatrix<double> conc_cluster; // UNIT: unitless | DIM: J, K
    xvector<double> excess_energy_cluster; // UNIT: eV | DIM: J, K

    // Thermo data
    xmatrix<double> prob_ideal_cluster; // DIM: Nc, J
    xtensor<double> prob_cluster; // DIM: Nc, Nt, J
    
  private:
  void free();
};

// Namespace for functions used by APDC
namespace apdc {
  void GetPhaseDiagram(_apdc_data& apdc_data);
  void GetPhaseDiagram(const string& aflowin, bool command_line_call);
  void GetPhaseDiagram(istream& infile);
  void ErrorChecks(_apdc_data& apdc_data);
  void GetBinodal(_apdc_data& apdc_data);
  void GetSpinodal(_apdc_data& apdc_data);
  xtensor<double> GetProbabilityCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, xvector<double>& excess_energy_cluster, const xmatrix<double>& prob_ideal_cluster, const xvector<double>& temp, const int max_num_atoms);
  xmatrix<double> GetProbabilityIdealCluster(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xvector<int>& multiplicity, const int max_num_atoms);
  bool CheckProbability(const xmatrix<double>& conc_macro, const xmatrix<double>& conc_cluster, const xmatrix<double>& prob);
  xmatrix<double> GetConcentrationMacro(const xvector<double>& conc_range, const int conc_npts, const int nelem);
  xvector<double> GetTemperature(const xvector<double>& temp_range, const int temp_npts);
  vector<xvector<int> > GetMultiplicity(const vector<xstructure>& vstr);
  xmatrix<double> GetConcentrationCluster(const vector<string>& elements, const vector<xstructure>& vstr);
  xmatrix<double> GetConcentrationCluster(const string& rundirpath, const int nstr, const int nelem);
  xvector<double> GetExcessEnergyCluster(const string& rundirpath, const xmatrix<double>& conc_cluster, const xvector<int>& natom);
  void RunATAT(const string& workdirpath, const string& rundirpath, const uint min_sleep);
  void GenerateFilesForATAT(const string& rundirpath, const string& lat_atat, const vector<xstructure>& vstr_aflow, const vector<xstructure>& vstr_atat, const vector<int>& mapstr);
  vector<xstructure> GetAFLOWXstructures(const string& plattice, const vector<string>& elements, const int num_threads, bool keep_all=true);
  string CreateLatForATAT(const string& plattice, const vector<string>& elements);
  vector<xstructure> GetATATXstructures(const string& lat, const uint max_num_atoms);
  vector<int> GetMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, const int num_threads);
}

#endif
