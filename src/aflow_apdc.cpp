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
#define _AFLOW_APDC_MAX_NUM_ATOMS 4
#define _AFLOW_APDC_MIN_SLEEP 2

// ###############################################################################
//            AFLOW Automatic Phase Diagram Constructor (APDC) (2022-)
// ###############################################################################

// **************************************************************************
// Class _apdc data
// **************************************************************************
// Constructor
_apdc_data::_apdc_data() {
  // Input data
  workdirpath = "";
  rootdirpath = "";
  plattice = "";
  elements.clear();
  max_num_atoms = 0;

  // Derived data
  alloyname = "";
  rundirpath = "";
  vstr_aflow.clear();
  lat_atat = "";
  vstr_atat.clear();
  mapstr.clear();

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
    // Input data
    rootdirpath = b.rootdirpath;
    plattice = b.plattice;
    elements = b.elements;
    max_num_atoms = b.max_num_atoms;
    // Derived data
    alloyname = b.alloyname;
    rundirpath = b.rundirpath;
    vstr_aflow = b.vstr_aflow;
    lat_atat = b.lat_atat;
    vstr_atat = b.vstr_atat;
    mapstr = b.mapstr;
    // Xstructure data
    multiplicity = b.multiplicity;
    composition = b.composition;
  }
  return *this;
}

// ***************************************************************************
// apdc::GetPhaseDiagram
// ***************************************************************************
namespace apdc {
  void GetPhaseDiagram(_apdc_data& apdc_data) {
    // Clean-up input data and check for errors
    apdc_data.workdirpath = aurostd::getPWD();
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
    aurostd::DirectoryMake(apdc_data.rundirpath);
    GetBinodal(apdc_data);
  }
}

// ***************************************************************************
// apdc::GetBinodal
// ***************************************************************************
namespace apdc {
  void GetBinodal(_apdc_data& apdc_data) {
    apdc_data.vstr_aflow = GetAFLOWXstructures(apdc_data.plattice, apdc_data.elements);
    apdc_data.lat_atat = CreateLatForATAT(apdc_data.plattice, apdc_data.elements);
    apdc_data.vstr_atat = GetATATXstructures(apdc_data.lat_atat, apdc_data.max_num_atoms);
    // map ATAT xstrs to AFLOW xstrs because ATAT cannot identify AFLOW xstrs
    apdc_data.mapstr = GetMapForXstructures(GetATATXstructures(apdc_data.lat_atat, _AFLOW_APDC_MAX_NUM_ATOMS), apdc_data.vstr_aflow);
    GenerateFilesForATAT(apdc_data.rundirpath, apdc_data.lat_atat, apdc_data.vstr_aflow, apdc_data.vstr_atat, apdc_data.mapstr);
    RunATAT(apdc_data.workdirpath, apdc_data.rundirpath);
    //apdc_data.multiplicity = GetMultiplicity(apdc_data.vstr_atat);
    //apdc_data.composition = GetComposition(apdc_data.elements, apdc_data.vstr_atat);
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
  void RunATAT(const string& workdirpath, const string& rundirpath) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    uint iter = 0;
    if (aurostd::substring2bool(aurostd::execute2string("mmaps", stdouterr_fsio), "command not found")) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "RunATAT():", "Missing mmaps program", _RUNTIME_ERROR_);
    }
    aurostd::execute("touch " + rundirpath + "/stop"); // pre-kill mmaps gracefully
    aurostd::RemoveFile(rundirpath + "/maps_is_running");
    aurostd::RemoveFile(rundirpath + "/predstr.out");
    chdir(rundirpath.c_str());
    string tmpfile = aurostd::TmpStrCreate();
    aurostd::execute("mmaps -d > " + tmpfile + " 2>&1 &");
    while (!aurostd::FileExist(rundirpath + "/predstr.out")) {
      iter++;
      if (LDEBUG) {cerr << "Sleeping, iter=" << iter << endl;}
      if (iter > 30) { // wait 30 times the minimum sleep (60 seconds)
        aurostd::RemoveFile(tmpfile);
        throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "RunATAT():", "mmaps is taking too long to predict structures", _RUNTIME_ERROR_);
      }
      aurostd::Sleep(_AFLOW_APDC_MIN_SLEEP);
    }
    aurostd::RemoveFile(tmpfile);
    chdir(workdirpath.c_str());
  }
}

// ***************************************************************************
// apdc::GenerateFilesForATAT
// ***************************************************************************
namespace apdc {
  void GenerateFilesForATAT(const string& rundirpath, const string& lat_atat, const vector<xstructure>& vstr_aflow, const vector<xstructure>& vstr_atat, const vector<int>& mapstr) {
    stringstream oss;
    vector<xstructure> vstr;
    // Generate lat.in file
    aurostd::string2file(lat_atat, rundirpath + "/lat.in");
    // Generate str.out and energy files
    uint msize = mapstr.size();
    for (uint i = 0; i < vstr_atat.size(); i++) {
      aurostd::DirectoryMake(rundirpath + "/" + aurostd::utype2string<uint>(i));
      oss << vstr_atat[i];
      aurostd::string2file(oss.str(), rundirpath + "/" + aurostd::utype2string<uint>(i) + "/str.out");
      aurostd::StringstreamClean(oss);
      if (i < msize && mapstr[i] != -1) {
        aurostd::string2file(aurostd::utype2string<double>(vstr_aflow[mapstr[i]].qm_E_cell) + "\n", rundirpath + "/" + aurostd::utype2string<uint>(i) + "/energy");
      }
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
// apdc::GetAFLOWXstructures
// ***************************************************************************
namespace apdc {
  vector<xstructure> GetAFLOWXstructures(const string& plattice, const vector<string>& elements, bool keep_all, uint num_proc) {
    vector<xstructure> vstr;
    vector<string> vstrlabel;
    string aflowlib, aflowurl;
    string alloyname = "";
    aflowlib::_aflowlib_entry entry;
    uint nary = elements.size();
    stringstream oss;
    for (uint i = 0; i < nary; i++) {alloyname += AVASP_Get_PseudoPotential_PAW_PBE(elements[i]);}
    aflowlib = "/common/LIB" + aurostd::utype2string<uint>(nary) + "/RAW/" + alloyname;
    aflowurl = "aflowlib.duke.edu:AFLOWDATA/LIB" + aurostd::utype2string<uint>(nary) + "_RAW/" + alloyname;
    if (plattice == "fcc") {
      if (nary >= 2) {
        for (uint i = 1; i < 30; i++) {vstrlabel.push_back(aurostd::utype2string<uint>(i));}
      }
      if (nary >= 3) {
        vstrlabel.push_back("TFCC001.ABC");vstrlabel.push_back("TFCC002.ABC");vstrlabel.push_back("TFCC003.ABC");
        for (uint i = 4; i < 17; i++) {
          vstrlabel.push_back("TFCC00" + aurostd::utype2string<uint>(i) + ".ABC");
          vstrlabel.push_back("TFCC00" + aurostd::utype2string<uint>(i) + ".BCA");
          vstrlabel.push_back("TFCC00" + aurostd::utype2string<uint>(i) + ".CAB");
        }
      }
    }
    else if (plattice == "bcc") {
      if (nary >= 2) {
        for (uint i = 58; i < 87; i++) {vstrlabel.push_back(aurostd::utype2string<uint>(i));}
      }
      if (nary >= 3) {
        vstrlabel.push_back("TBCC001.ABC");vstrlabel.push_back("TBCC002.ABC");vstrlabel.push_back("TBCC003.ABC");
        for (uint i = 4; i < 17; i++) {
          vstrlabel.push_back("TBCC00" + aurostd::utype2string<uint>(i) + ".ABC");
          vstrlabel.push_back("TBCC00" + aurostd::utype2string<uint>(i) + ".BCA");
          vstrlabel.push_back("TBCC00" + aurostd::utype2string<uint>(i) + ".CAB");
        }
      }
    }
    else if (plattice == "hcp") {
      if (nary >= 2) {
        for (uint i = 115; i < 178; i++) {vstrlabel.push_back(aurostd::utype2string<uint>(i));}
      }
    }
    for (uint i = 0; i < vstrlabel.size(); i++) {
        entry.Load(aurostd::file2string(aflowlib + "/" + vstrlabel[i] + "/aflowlib.out"), oss);
        entry.aurl = aflowurl + "/" + vstrlabel[i];
        if (pflow::loadXstructures(entry, oss, false)) { // initial = unrelaxed; final = relaxed
          entry.vstr[0].qm_E_cell = entry.enthalpy_cell; // ATAT needs energy per cell
          if (keep_all || compare::structuresMatch(entry.vstr[0], entry.vstr[entry.vstr.size() - 1], true, num_proc)) {vstr.push_back(entry.vstr[0]);}
        }
        entry.clear();
    }
    return vstr;
  }
}

// ***************************************************************************
// apdc::CreateLatForATAT
// ***************************************************************************
namespace apdc {
  string CreateLatForATAT(const string& plattice, const vector<string>& elements) {
    stringstream oss;
    xmatrix<double> lattice(3,3);
    xvector<double> angles = 90.0 * aurostd::ones_xv<double>(3);
    xvector<double> coorsys = aurostd::ones_xv<double>(3);
    oss.precision(_DOUBLE_WRITE_PRECISION_MAX_);
    oss.setf(std::ios::fixed,std::ios::floatfield);
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
    return oss.str();
  }
}

// ***************************************************************************
// apdc::GetATATXstructures
// ***************************************************************************
namespace apdc {
  vector<xstructure> GetATATXstructures(const string& lat, const uint max_num_atoms) {
    vector<xstructure> vstr;
    stringstream oss;
    vector<string> vinput, tokens;
    string tmpfile = aurostd::TmpStrCreate();
    aurostd::string2file(lat, tmpfile);
    string sstr = aurostd::execute2string("genstr -n " + aurostd::utype2string<uint>(max_num_atoms) + " -l " + tmpfile, stdouterr_fsio);
    aurostd::RemoveFile(tmpfile);
    if (sstr == "" || aurostd::substring2bool(sstr, "Unable to open lattice file")) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "GetATATXstructures():", "Invalid lat.in file", _FILE_CORRUPT_);
    }
    else if (aurostd::substring2bool(sstr, "command not found")) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "GetATATXstructures():", "Missing genstr program", _RUNTIME_ERROR_);
    }
    aurostd::string2vectorstring(sstr, vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      if (aurostd::substring2bool(vinput[line], "end")) {
        vstr.push_back(xstructure(oss, IOATAT_STR));
        aurostd::StringstreamClean(oss);
      }
      else {
        oss << vinput[line] << endl;
      }
    }
    return vstr;
  }
}

// ***************************************************************************
// apdc::GetMapForXstructures
// ***************************************************************************
namespace apdc {
  vector<int> GetMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, uint num_proc) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    vector<int> mapstr;
    bool match;
    for (uint i = 0; i < vstr1.size(); i++) {
      match = false;
      for (uint j = 0; j < vstr2.size() && !match; j++) {
        if (compare::structuresMatch(vstr1[i], vstr2[j], true, num_proc)) {
          if (LDEBUG) {cerr << "str1 i=" << i << " matched str2 j=" << j << endl << "str1=" << vstr1[i] << endl << "str2=" << vstr2[j] << endl;}
          mapstr.push_back(j); 
          match = true;
        }
      }
      if (!match) {mapstr.push_back(-1);}
    }
    return mapstr;
  }
}






#endif
