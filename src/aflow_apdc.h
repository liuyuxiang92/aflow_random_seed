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

#define _APDC_NUM_PROC_ 8
#define _APDC_MIN_SLEEP_ 2

// Class _apdc data
class _apdc_data {
  public:
    _apdc_data();
    ~_apdc_data();
    const _apdc_data& operator=(const _apdc_data &b);
    
    // Input data
    int num_threads;
    string workdirpath;
    string rootdirpath;
    string plattice;
    vector<string> elements;
    int aflow_max_num_atoms;
    int max_num_atoms;
    vector<int> conc_macro_npts;
    xmatrix<double> conc_macro; // unitless
    int temp_npts;
    vector<double> temp;

    // Derived data
    string alloyname;
    string rundirpath;
    vector<xstructure> vstr_aflow;
    string lat_atat;
    vector<xstructure> vstr_atat;
    vector<int> mapstr;

    // Structure data
    xvector<int> multiplicity;
    xmatrix<double> conc; // unitless
    xvector<double> excess_energies; // eV

    // Thermo data
    xvector<double> prob_rand;
    xvector<double> prob;
    
  private:
  void free();
};

// Namespace for functions used by APDC
namespace apdc {
  void GetPhaseDiagram(_apdc_data& apdc_data);
  void GetPhaseDiagram(const string& aflowin, bool elements_only);
  void GetPhaseDiagram(istream& infile);
  void ErrorChecks(_apdc_data& apdc_data);
  void GetBinodal(_apdc_data& apdc_data);
  void GetSpinodal(_apdc_data& apdc_data);
  void RunATAT(const string& workdirpath, const string& rundirpath);
  void GenerateFilesForATAT(const string& rundirpath, const string& lat_atat, const vector<xstructure>& vstr_aflow, const vector<xstructure>& vstr_atat, const vector<int>& mapstr);
  vector<xvector<int> > GetMultiplicity(const vector<xstructure>& vstr);
  xvector<double> GetProbabilityRandom(const xmatrix<double>& conc, const xvector<int>& multiplicity, const int max_num_atoms);
  xmatrix<double> GetConcentration(const vector<string>& elements, const vector<xstructure>& vstr);
  xmatrix<double> GetConcentration(const string& rundirpath, const int nstr, const int nelem);
  xvector<double> GetExcessEnergy(const string& rundirpath, const xmatrix<double>& conc, const xvector<int>& natom);
  vector<xstructure> GetAFLOWXstructures(const string& plattice, const vector<string>& elements, const int num_threads, bool keep_all=true);
  string CreateLatForATAT(const string& plattice, const vector<string>& elements);
  vector<xstructure> GetATATXstructures(const string& rundirpath, const uint max_num_atoms);
  vector<int> GetMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, const int num_threads);
}

#endif
