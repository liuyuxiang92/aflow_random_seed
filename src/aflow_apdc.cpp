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

#define _AFLOW_APDC_ALAT 4.0

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
    // Clean-up input data and check for errors
    apdc_data.rundirpath = aurostd::CleanFileName(apdc_data.rundirpath);
    apdc_data.plattice = aurostd::tolower(apdc_data.plattice);
    aurostd::sort_remove_duplicates(apdc_data.elements);
    if (apdc_data.plattice != "fcc" && apdc_data.plattice != "bcc" && apdc_data.plattice != "hcp") {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "GetPhaseDiagram():", "Invalid parent lattice", _INPUT_ILLEGAL_);
    }
    if (apdc_data.elements.size() < 2) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "GetPhaseDiagram():", "Alloy must be at least binary", _VALUE_ERROR_);
    }
    // Binodal
    for (uint i = 0; i < apdc_data.elements.size(); i++) {apdc_data.alloyname += apdc_data.elements[i];}
    apdc_data.rundirpath += apdc_data.rootdirpath + "/" + pflow::arity_string(apdc_data.elements.size(), false, false) + "/" + apdc_data.plattice + "/" + apdc_data.alloyname;
    GetBinodal(apdc_data);
  }
}

// ***************************************************************************
// apdc::GetBinodal
// ***************************************************************************
namespace apdc {
  void GetBinodal(_apdc_data& apdc_data) {
    apdc_data.vstr = GetXstructuresForATAT(apdc_data.plattice, apdc_data.elements);
    GenerateFilesForATAT(apdc_data.rundirpath, apdc_data.plattice, apdc_data.elements, apdc_data.vstr);
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
// apdc::GenerateFilesForATAT
// ***************************************************************************
namespace apdc {
  void GenerateFilesForATAT(const string& rundirpath, const string& plattice, const vector<string>& elements, const vector<xstructure>& vstr) {
    stringstream oss;
    xmatrix<double> lattice(3,3);
    if (!aurostd::DirectoryMake(rundirpath)) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "GenerateFilesForATAT():", "Cannot create directory", _INPUT_ILLEGAL_); 
    }
    // Generate lat.in file
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
      lattice(1, 1) = 0.5; lattice(2, 1) = -std::sqrt(3.0) / 2.0; lattice(3, 1) = 0.0;
      lattice(1, 2) = 0.5; lattice(2, 2) = std::sqrt(3.0) / 2.0; lattice(3, 2) = 0.0;
      lattice(1, 3) = 0.0; lattice(2, 3) = 0.0; lattice(3, 3) = std::sqrt(8.0 / 3.0);
    }
    lattice *= _AFLOW_APDC_ALAT;
    oss << aurostd::eye<double>(3, 3) << endl << lattice << endl << 0.0 * aurostd::ones_xv<double>(3) << " ";
    for (uint i = 0; i < elements.size(); i++) {
      oss << elements[i] << ",";
    }
    if (plattice == "hcp") {
      oss << endl << 0.5 * aurostd::ones_xv<double>(3) << " ";
      for (uint i = 0; i < elements.size(); i++) {
        oss << elements[i] << ",";
      }
    }
    oss << endl;
    aurostd::string2file(oss.str(), rundirpath + "/lat.in");
    aurostd::StringstreamClean(oss);
    // Generate str.out and energy
    for (uint i = 0; i < vstr.size(); i++) {
      aurostd::DirectoryMake(rundirpath + "/" + aurostd::utype2string<uint>(i));
      oss << vstr[i] << endl;
      aurostd::string2file(oss.str(), rundirpath + "/" + aurostd::utype2string<uint>(i) + "/str.out");
      aurostd::string2file(aurostd::utype2string<double>(vstr[i].qm_E_cell) + "\n", rundirpath + "/" + aurostd::utype2string<uint>(i) + "/energy");
      aurostd::StringstreamClean(oss);
    }
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
        cerr << entry.spacegroup_orig << " " << entry.spacegroup_relax << " " << entry.enthalpy_cell << endl;
        entry.aurl = aflowurl + "/" + aurostd::utype2string<uint>(i);
        if (pflow::loadXstructures(entry, oss, false)) { //initial = unrelaxed; final = relaxed
          entry.vstr[0].iomode = IOATAT_STR;
          entry.vstr[0].qm_E_cell = entry.enthalpy_cell; // ATAT needs energy per cell
          if (entry.spacegroup_orig == entry.spacegroup_relax) {vstr.push_back(entry.vstr[0]);}
        }
        entry.clear();
    }
    return vstr;
  }
}










#endif
