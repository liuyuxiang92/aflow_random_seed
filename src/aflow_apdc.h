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

#define _AFLOW_APDC_NUM_PROC 8

// Class _apdc data
class _apdc_data {
  public:
    _apdc_data();
    ~_apdc_data();
    const _apdc_data& operator=(const _apdc_data &b);
    
    // Input data
    string workdirpath;
    string rootdirpath;
    string plattice;
    vector<string> elements;
    int max_num_atoms;

    // Derived data
    string alloyname;
    string rundirpath;
    vector<xstructure> vstr_aflow;
    string lat_atat;
    vector<xstructure> vstr_atat;
    vector<int> mapstr;

    // Structure data
    xvector<int> multiplicity;
    xmatrix<double> conc;
    xvector<double> excess_energies; // eV
    xvector<double> prob_rand;

    // Thermo data
    xvector<double> prob;
    
  

  private:
  void free();
};

// Namespace for functions used by APDC
namespace apdc {
  void GetPhaseDiagram(_apdc_data& apdc_data);
  void GetBinodal(_apdc_data& apdc_data);
  void GetSpinodal(_apdc_data& apdc_data);
  void RunATAT(const string& workdirpath, const string& rundirpath);
  void GenerateFilesForATAT(const string& rundirpath, const string& lat_atat, const vector<xstructure>& vstr_aflow, const vector<xstructure>& vstr_atat, const vector<int>& mapstr);
  vector<xvector<int> > GetMultiplicity(const vector<xstructure>& vstr);
  xvector<double> GetProbabilityRandom(const xmatrix<double>& conc, const xvector<int>& multiplicity, const int max_num_atoms);
  xmatrix<double> GetConcentration(const vector<string>& elements, const vector<xstructure>& vstr);
  xmatrix<double> GetConcentration(const string& rundirpath, const int nstr, const int nelem);
  xvector<double> GetExcessEnergy(const string& rundirpath, const xmatrix<double>& conc, const xvector<int>& natom);
  vector<xstructure> GetAFLOWXstructures(const string& plattice, const vector<string>& elements, bool keep_all=true, uint num_proc=_AFLOW_APDC_NUM_PROC);
  string CreateLatForATAT(const string& plattice, const vector<string>& elements);
  vector<xstructure> GetATATXstructures(const string& rundirpath, const uint max_num_atoms);
  vector<int> GetMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, uint num_proc=_AFLOW_APDC_NUM_PROC);
}

#endif
