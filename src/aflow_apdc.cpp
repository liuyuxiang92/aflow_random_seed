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
  rundirpath = "";
  plattice = "";
  elements.clear();

  // Calculated data
  alloyname = "";
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
    rundirpath = b.rundirpath;
    plattice = b.plattice;
    elements = b.elements;
    alloyname = b.alloyname;
    vstr = b.vstr;
  }
  return *this;
}

// ***************************************************************************
// apdc::GetPhaseDiagram
// ***************************************************************************
namespace apdc {
  void GetPhaseDiagram(_apdc_data& apdc_data) {
    GetBinodal(apdc_data);
  }
}

// ***************************************************************************
// apdc::GetBinodal
// ***************************************************************************
namespace apdc {
  void GetBinodal(_apdc_data& apdc_data) {
    GenerateFilesForATAT(apdc_data.rundirpath, apdc_data.plattice, apdc_data.elements);
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
// apdc::GenerateFilesForATAT
// ***************************************************************************
namespace apdc {
  void GenerateFilesForATAT(const string& rundirpath, const string& plattice, const vector<string>& elements) {
  }
}









#endif
