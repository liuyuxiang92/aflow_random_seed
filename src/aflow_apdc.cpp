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
//            AFLOW Automatic Phase Diagram Constructor (APDC) (2022-)
// ###############################################################################

// **************************************************************************
// Class _apdc data
// **************************************************************************
// Constructor
_apdc_data::_apdc_data() {
  // Input data
  rootdirpath = "";
  plattice = "";
  elements.clear();

  // Calculated data
  alloyname = "";
  rundirpath = "";
  vstr.clear();
  
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
    rootdirpath = b.rootdirpath;
    plattice = b.plattice;
    elements = b.elements;
    alloyname = b.alloyname;
    rundirpath = b.rundirpath;
    vstr = b.vstr;
  }
  return *this;
}

// ***************************************************************************
// apdc::GetPhaseDiagram
// ***************************************************************************
namespace apdc {
  void GetPhaseDiagram(_apdc_data& apdc_data) {
    // Clean-up input data
    apdc_data.rundirpath = aurostd::CleanFileName(apdc_data.rundirpath);
    apdc_data.plattice = aurostd::tolower(apdc_data.plattice);
    aurostd::sort_remove_duplicates(apdc_data.elements);
    // Binodal
    for (uint i = 0; i < apdc_data.elements.size(); i++) {apdc_data.alloyname += apdc_data.elements[i];}
    apdc_data.rundirpath = apdc_data.rundirpath + "/" + pflow::arity_string(apdc_data.elements.size(), false, false) + "/" + apdc_data.plattice + "/" + apdc_data.alloyname;
    GetBinodal(apdc_data);
  }
}

// ***************************************************************************
// apdc::GetBinodal
// ***************************************************************************
namespace apdc {
  void GetBinodal(_apdc_data& apdc_data) {
    apdc_data.vstr = GetXstructuresForATAT(apdc_data.plattice, apdc_data.elements);
    GenerateFilesForATAT(apdc_data.rundirpath, apdc_data.plattice, apdc_data.vstr);
  }
}

// ***************************************************************************
// apdc::GetSpinodal
// ***************************************************************************
namespace apdc {
  void GetSpinodal(_apdc_data& apdc_data) {
  }
}

// ***************************************************************************
// apdc::GetXstructuresForATAT
// ***************************************************************************
namespace apdc {
  vector<xstructure> GetXstructuresForATAT(const string& plattice, const vector<string>& elements) {
    vector<xstructure> vstr;
    aflowlib::_aflowlib_entry entry;
    uint istart, iend;
    stringstream oss;
    int nary = elements.size();
    string alloyname = "";
    for (uint i = 0; i < elements.size(); i++) {alloyname += AVASP_Get_PseudoPotential_PAW_PBE(elements[i]);}
    string aflowlib = "/common/LIB" + aurostd::utype2string<int>(nary) + "/RAW/" + alloyname;
    string aflowurl = "aflowlib.duke.edu:AFLOWDATA/LIB" + aurostd::utype2string<int>(nary) + "_RAW/" + alloyname;
    if (plattice == "fcc") {
      if (nary == 2) {
        istart = 1;
        iend = 29;
      }
    }
    else if (plattice == "bcc") {
      if (nary == 2) {
        istart = 58;
        iend = 86;
      }
    }
    else if (plattice == "hcp") {
      if (nary == 2) {
        istart = 115;
        iend = 177;
      }
    }
    for (uint i = istart; i <= iend; i++) {
        entry.Load(aflowlib + "/" + aurostd::utype2string<uint>(i), oss);
        entry.aurl = aflowurl + "/" + aurostd::utype2string<uint>(i);
        pflow::loadXstructures(entry, oss, false); // initial = unrelaxed; final = relaxed
        entry.vstr[0].iomode = IOATAT_STR;
        entry.vstr[0].qm_E_cell = entry.enthalpy_cell; // ATAT needs energy per cell
        if (entry.spacegroup_orig == entry.spacegroup_relax) {vstr.push_back(entry.vstr[0]);}
        entry.clear();
    }
    return vstr;
  }
}

// ***************************************************************************
// apdc::GenerateFilesForATAT
// ***************************************************************************
namespace apdc {
  void GenerateFilesForATAT(const string& rundirpath, const string& plattice, const vector<xstructure>& vstr) {
  }
}









#endif
