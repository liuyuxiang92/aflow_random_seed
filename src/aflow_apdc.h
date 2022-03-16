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

#define _MAX_NUM_ATOMS 4
#define _NUM_PROC 8

// Class _apdc data
class _apdc_data {
  public:
    _apdc_data();
    ~_apdc_data();
    const _apdc_data& operator=(const _apdc_data &b);
    
    // Input data
    string rootdirpath;
    string plattice;
    vector<string> elements;

    // Derived data
    string alloyname;
    string rundirpath;
    vector<xstructure> vstr;

    // Xstructure data
    vector<uint> multiplicity;
    vector<xvector<double> > composition;
  

  private:
    void free();
};

// Namespace for functions used by APDC
namespace apdc {
  void GetPhaseDiagram(_apdc_data& apdc_data);
  void GetBinodal(_apdc_data& apdc_data);
  void GetSpinodal(_apdc_data& apdc_data);
  void RunATAT(const string& rundirpath);
  void GenerateFilesForATAT(const string& rundirpath, const string& plattice, const vector<string>& elements, const vector<xstructure>& _vstr, bool use_atat_xstr=true);
  vector<uint> GetMultiplicity(const vector<xstructure>& vstr);
  vector<xvector<double> > GetComposition(const vector<string>& elements, const vector<xstructure>& vstr);
  vector<xstructure> GetAFLOWXstructures(const string& plattice, const vector<string>& elements, uint num_proc=_NUM_PROC);
  vector<xstructure> GetATATXstructures(const string& rundirpath, uint max_num_atoms=_MAX_NUM_ATOMS);
  vector<uint> GetDictionaryForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, uint num_proc=_NUM_PROC);
}

#endif
