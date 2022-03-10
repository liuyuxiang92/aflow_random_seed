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

  // Derived data
  alloyname = "";
  rundirpath = "";
  vstr.clear();

  // Xstrucutre data
  multiplicity.clear();
  composition.clear();
  
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
    apdc_data.rootdirpath = aurostd::CleanFileName(apdc_data.rootdirpath);
    apdc_data.plattice = aurostd::tolower(apdc_data.plattice);
    aurostd::sort_remove_duplicates(apdc_data.elements);
    if (!aurostd::DirectoryMake(apdc_data.rootdirpath)) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "GetPhaseDiagram():", "Cannot create directory", _FILE_ERROR_);
    }
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
    apdc_data.vstr = GetAFLOWXstructuresForATAT(apdc_data.plattice, apdc_data.elements);
    apdc_data.multiplicity = GetMultiplicity(apdc_data.vstr);
    apdc_data.composition = GetComposition(apdc_data.elements, apdc_data.vstr);
    //for (uint i = 0; i < apdc_data.composition.size(); i++){cerr << apdc_data.composition[i] << endl;}
    GenerateFilesForATAT(apdc_data.rundirpath, apdc_data.plattice, apdc_data.elements, apdc_data.vstr);
    RunATAT(apdc_data.rundirpath);
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
// apdc::RunATAT
// ***************************************************************************
namespace apdc {
  void RunATAT(const string& rundirpath) {
    cerr << rundirpath << endl;
  }
}

// ***************************************************************************
// apdc::GenerateFilesForATAT
// ***************************************************************************
namespace apdc {
  void GenerateFilesForATAT(const string& rundirpath, const string& plattice, const vector<string>& elements, const vector<xstructure>& vstr) {
    stringstream oss;
    xmatrix<double> lattice(3,3);
    xvector<double> angles = 90.0 * aurostd::ones_xv<double>(3);
    xvector<double> coorsys = aurostd::ones_xv<double>(3);
    oss.precision(_DOUBLE_WRITE_PRECISION_MAX_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
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
      lattice = aurostd::eye<double>(3, 3);
      coorsys(3) = std::sqrt(8.0 / 3.0);
      angles(3) = 120.0;
    }
    for (uint i = 1; i <= 3; i++) {oss << _AFLOW_APDC_ALAT * coorsys(i) << " ";}
    for (uint i = 1; i <= 3; i++) {oss << angles(i) << " ";}
    oss << endl;
    for (uint i = 1; i <= 3; i++) {
      for (uint j = 1; j <= 3; j++) {
        oss << lattice(i, j) << " ";
      }
      oss << endl;
    }
    oss << 0.0 << " " << 0.0 << " " << 0.0 << " ";
    for (uint i = 0; i < elements.size(); i++) {
      oss << elements[i] << ",";
    }
    if (plattice == "hcp") {
      oss << endl << 2.0 / 3.0 << " " << 1.0 / 3.0 << " " << 0.5 << " ";
      for (uint i = 0; i < elements.size(); i++) {
        oss << elements[i] << ",";
      }
    }
    oss << endl;
    aurostd::string2file(oss.str(), rundirpath + "/lat.in");
    aurostd::StringstreamClean(oss);
    // Generate str.out and energy files
    for (uint i = 0; i < vstr.size(); i++) {
      aurostd::DirectoryMake(rundirpath + "/" + aurostd::utype2string<uint>(i));
      oss << vstr[i];
      aurostd::string2file(oss.str(), rundirpath + "/" + aurostd::utype2string<uint>(i) + "/str.out");
      aurostd::string2file(aurostd::utype2string<double>(vstr[i].qm_E_cell) + "\n", rundirpath + "/" + aurostd::utype2string<uint>(i) + "/energy");
      aurostd::StringstreamClean(oss);
    }
  }
}

// ***************************************************************************
// apdc::GetMultiplicity
// ***************************************************************************
namespace apdc {
  vector<uint> GetMultiplicity(const vector<xstructure>& vstr) {
    vector<uint> multiplicity;
    uint natom, fact_prod;
    for (uint i = 0; i < vstr.size(); i++) {
      natom = 0;
      fact_prod = 1;
      for (uint j = 0; j < vstr[i].num_each_type.size(); j++) {
        natom += vstr[i].num_each_type[j];
        fact_prod *= aurostd::factorial(vstr[i].num_each_type[j]);
      }
      multiplicity.push_back(aurostd::factorial(natom) / fact_prod);
    }
    return multiplicity;
  }
}

// ***************************************************************************
// apdc::GetComposition
// ***************************************************************************
namespace apdc {
  vector<xvector<double> > GetComposition(const vector<string>& elements, const vector<xstructure>& vstr) {
    vector<xvector<double> > composition;
    uint nary = elements.size();
    int ie = -1;
    xvector<double> stoich;
    vector<string> str_elements;
    xstructure tmp = vstr[0];
    for (uint i = 0; i < vstr.size(); i++) {
      if (nary != vstr[i].stoich_each_type.size()) {
        str_elements = vstr[i].GetElements(true, true);
        stoich = 0.0 * aurostd::ones_xv<double>(nary);
        for (uint j = 0; j < nary; j++) {
          if (aurostd::WithinList(str_elements, elements[j], ie)) {stoich(j + 1) = vstr[i].stoich_each_type[ie];}
        }
      }
      else {
        stoich = aurostd::vector2xvector(aurostd::deque2vector(vstr[i].stoich_each_type));
      }
      composition.push_back(stoich);
    }
    return composition;
  }
}

// ***************************************************************************
// apdc::GetAFLOWXstructuresForATAT
// ***************************************************************************
namespace apdc {
  vector<xstructure> GetAFLOWXstructuresForATAT(const string& plattice, const vector<string>& elements, bool use_xstr_atat) {
    vector<xstructure> vstr, vstr_atat;
    string aflowlib, aflowurl;
    string alloyname = "";
    aflowlib::_aflowlib_entry entry;
    uint istart, iend, nary = elements.size();
    stringstream oss;
    for (uint i = 0; i < nary; i++) {alloyname += AVASP_Get_PseudoPotential_PAW_PBE(elements[i]);}
    aflowlib = "/common/LIB" + aurostd::utype2string<uint>(nary) + "/RAW/" + alloyname;
    aflowurl = "aflowlib.duke.edu:AFLOWDATA/LIB" + aurostd::utype2string<uint>(nary) + "_RAW/" + alloyname;
    if (use_xstr_atat) {
      vstr_atat = GetATATXstructuresForAFLOW(plattice, elements);} // read the ATAT geometries
      vector<uint> dict = DictAFLOW2ATAT(plattice, nary); // dictionary ATAT to AFLOW
    if (nary == 2) {
      if (plattice == "fcc") {
        istart = 1;
        iend = 29;
      }
      else if (plattice == "bcc") {
        istart = 58;
        iend = 86;
      }
      else if (plattice == "hcp") {
        istart = 115;
        iend = 177;
      }
    }
    else if (nary == 3) {
    }
    for (uint i = istart; i <= iend; i++) {
        entry.Load(aurostd::file2string(aflowlib + "/" + aurostd::utype2string<uint>(i) + "/aflowlib.out"), oss);
        entry.aurl = aflowurl + "/" + aurostd::utype2string<uint>(i);
        if (!use_xstr_atat && pflow::loadXstructures(entry, oss, false)) { // initial = unrelaxed; final = relaxed
          entry.vstr[0].FixLattices();
          entry.vstr[0].ReScale(1.0);
          entry.vstr[0].BringInCell();
          entry.vstr[0].iomode = IOATAT_STR;
          entry.vstr[0].qm_E_cell = entry.enthalpy_cell; // ATAT needs energy per cell
          if (entry.spacegroup_orig == entry.spacegroup_relax) {vstr.push_back(entry.vstr[0]);}
        }
        else {
          vstr_atat[dict[i]].iomode = IOATAT_STR;
          vstr_atat[dict[i]].qm_E_cell = entry.enthalpy_cell; // ATAT needs energy per cell
          if (entry.spacegroup_orig == entry.spacegroup_relax) {vstr.push_back(vstr_atat[dict[i]]);}
        }
        entry.clear();
    }
    return vstr;
  }
}

// ***************************************************************************
// apdc::GetATATXstructuresForAFLOW
// ***************************************************************************
namespace apdc {
  vector<xstructure> GetATATXstructuresForAFLOW(const string& plattice, const vector<string>& elements) {
    vector<xstructure> vstr;
    xstructure str;
    _atom atom;
    uint nary = elements.size();
    xmatrix<double> lat = aurostd::eye<double>(3, 3);
    xvector<double> coord = aurostd::ones_xv<double>(3);
    if (nary == 2) {
      if (plattice == "fcc") {
        // #1
        lat(1,1) =  0.50; lat(1,2) =  0.00; lat(1,3) =  0.50;
        lat(2,1) =  0.50; lat(2,2) =  0.50; lat(2,3) =  0.00;
        lat(3,1) =  0.00; lat(3,2) =  0.50; lat(3,3) =  0.50;
        str.lattice = _AFLOW_APDC_ALAT * lat; str.ReScale(1.0);
        coord(1) = 0.00; coord(2) = 0.00; coord(3) = 0.00;
        atom.name = atom.cleanname = elements[0]; atom.name_is_given = TRUE; atom.fpos = coord; atom.cpos = str.f2c * atom.fpos;
        str.AddAtom(atom); atom.clear();
        str.SpeciesPutAlphabetic(); str.MakeBasis(); str.MakeTypes(); str.buildGenericTitle(); vstr.push_back(str); str.clear();
        // #2
        lat(1,1) =  0.50; lat(1,2) =  0.00; lat(1,3) =  0.50;
        lat(2,1) =  0.50; lat(2,2) =  0.50; lat(2,3) =  0.00;
        lat(3,1) =  0.00; lat(3,2) =  0.50; lat(3,3) =  0.50;
        str.lattice = _AFLOW_APDC_ALAT * lat; str.ReScale(1.0);
        coord(1) = 0.00; coord(2) = 0.00; coord(3) = 0.00;
        atom.name = atom.cleanname = elements[1]; atom.name_is_given = TRUE; atom.fpos = coord; atom.cpos = str.f2c * atom.fpos;
        str.AddAtom(atom); atom.clear();
        str.SpeciesPutAlphabetic(); str.MakeBasis(); str.MakeTypes(); str.buildGenericTitle(); vstr.push_back(str); str.clear();
        // #3
        lat(1,1) = -0.50; lat(1,2) =  0.50; lat(1,3) =  1.00;
        lat(2,1) = -0.50; lat(2,2) =  1.00; lat(2,3) =  0.50;
        lat(3,1) = -1.00; lat(3,2) =  0.50; lat(3,3) =  0.50;
        str.lattice = _AFLOW_APDC_ALAT * lat; str.ReScale(1.0);
        coord(1) = 0.50; coord(2) = 0.50; coord(3) = 0.50;
        atom.name = atom.cleanname = elements[0]; atom.name_is_given = TRUE; atom.fpos = coord; atom.cpos = str.f2c * atom.fpos;
        str.AddAtom(atom); atom.clear();
        coord(1) = 0.00; coord(2) = 0.00; coord(3) = 0.00;
        atom.name = atom.cleanname = elements[1]; atom.name_is_given = TRUE; atom.fpos = coord; atom.cpos = str.f2c * atom.fpos;
        str.AddAtom(atom); atom.clear();
        str.SpeciesPutAlphabetic(); str.MakeBasis(); str.MakeTypes(); str.buildGenericTitle(); vstr.push_back(str); str.clear();
        // #4
        lat(1,1) =  0.00; lat(1,2) = -0.50; lat(1,3) = -0.50;
        lat(2,1) =  0.00; lat(2,2) = -0.50; lat(2,3) =  0.50;
        lat(3,1) =  1.00; lat(3,2) =  0.00; lat(3,3) =  0.00;
        str.lattice = _AFLOW_APDC_ALAT * lat; str.ReScale(1.0);
        coord(1) = 0.00; coord(2) = 0.00; coord(3) = 0.00;
        atom.name = atom.cleanname = elements[0]; atom.name_is_given = TRUE; atom.fpos = coord; atom.cpos = str.f2c * atom.fpos;
        str.AddAtom(atom); atom.clear();
        coord(1) = 0.50; coord(2) = 0.50; coord(3) = 0.50;
        atom.name = atom.cleanname = elements[1]; atom.name_is_given = TRUE; atom.fpos = coord; atom.cpos = str.f2c * atom.fpos;
        str.AddAtom(atom); atom.clear();
        str.SpeciesPutAlphabetic(); str.MakeBasis(); str.MakeTypes(); str.buildGenericTitle(); vstr.push_back(str); str.clear();
        // #5
        lat(1,1) = -0.50; lat(1,2) =  0.00; lat(1,3) = -0.50;
        lat(2,1) =  0.00; lat(2,2) = -0.50; lat(2,3) =  0.50;
        lat(3,1) = -1.00; lat(3,2) =  1.00; lat(3,3) =  1.00;
        str.lattice = _AFLOW_APDC_ALAT * lat; str.ReScale(1.0);
        coord(1) = 2.0/3.0; coord(2) = 1.0/3.0; coord(3) = 2.0/3.0;
        atom.name = atom.cleanname = elements[0]; atom.name_is_given = TRUE; atom.fpos = coord; atom.cpos = str.f2c * atom.fpos;
        str.AddAtom(atom); atom.clear();
        coord(1) = 1.0/3.0; coord(2) = 2.0/3.0; coord(3) = 1.0/3.0;
        atom.name = atom.cleanname = elements[0]; atom.name_is_given = TRUE; atom.fpos = coord; atom.cpos = str.f2c * atom.fpos;
        str.AddAtom(atom); atom.clear();
        coord(1) = 0.00; coord(2) = 0.00; coord(3) = 0.00;
        atom.name = atom.cleanname = elements[1]; atom.name_is_given = TRUE; atom.fpos = coord; atom.cpos = str.f2c * atom.fpos;
        str.AddAtom(atom); atom.clear();
        str.SpeciesPutAlphabetic(); str.MakeBasis(); str.MakeTypes(); str.buildGenericTitle(); vstr.push_back(str); str.clear();
      }
      else if (plattice == "bcc") {
      }
      else if (plattice == "hcp") {
      }
    }
    return vstr;
  }
}

// ***************************************************************************
// apdc::DictAFLOW2ATAT
// ***************************************************************************
namespace apdc {
  vector<uint> DictAFLOW2ATAT(const string& plattice, const uint nary) {
    vector<uint> dict;
    if (nary == 2) {
      if (plattice == "fcc") {
        //dict = {1, 0, 3, 2, 9, 8};
        dict = {0, 1, 2, 3, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
      }
      else if (plattice == "bcc") {
      }
      else if (plattice == "hcp") {
      }
    }
    else if (nary == 3) {
    }
    return dict;
  }
}







#endif
