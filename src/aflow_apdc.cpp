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
// Strings for I/O
#define _APDC_STR_OPT_ string("[AFLOW_APDC]")

// ###############################################################################
//            AFLOW Automatic Phase Diagram Constructor (APDC) (2022-)
// ###############################################################################

// **************************************************************************
// Class _apdc data
// **************************************************************************
// Constructor
_apdc_data::_apdc_data() {
  // Input data
  num_threads = 0;
  workdirpath = "";
  rootdirpath = "";
  plattice = "";
  elements.clear();
  aflow_max_num_atoms = 0;
  max_num_atoms = 0;
  conc_macro_npts.clear();
  conc_macro.clear();
  temp_npts = 0;
  temp.clear();
  // Derived data
  alloyname = "";
  rundirpath = "";
  vstr_aflow.clear();
  lat_atat = "";
  vstr_atat.clear();
  mapstr.clear();
  // Structure data
  multiplicity.clear();
  conc.clear();
  excess_energies.clear();
  // Thermo data
  prob_rand.clear();
  prob.clear();
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
    num_threads = b.num_threads;
    workdirpath = b.workdirpath;
    rootdirpath = b.rootdirpath;
    plattice = b.plattice;
    elements = b.elements;
    max_num_atoms = b.max_num_atoms;
    conc_macro_npts = b.conc_macro_npts;
    conc_macro = b.conc_macro;
    temp_npts = b.temp_npts;
    temp = b.temp;
    // Derived data
    alloyname = b.alloyname;
    rundirpath = b.rundirpath;
    vstr_aflow = b.vstr_aflow;
    lat_atat = b.lat_atat;
    vstr_atat = b.vstr_atat;
    mapstr = b.mapstr;
    // Structure data
    multiplicity = b.multiplicity;
    conc = b.conc;
    excess_energies = b.excess_energies;
    // Thermo data
    prob_rand = b.prob_rand;
    prob = b.prob;
  }
  return *this;
}

// ***************************************************************************
// apdc::GetPhaseDiagram
// ***************************************************************************
namespace apdc {
  void GetPhaseDiagram(_apdc_data& apdc_data) {
    // Clean-up input data and check for errors
    if (XHOST.vflag_control.flag("XPLUG_NUM_THREADS") && !(XHOST.vflag_control.flag("XPLUG_NUM_THREADS_MAX"))) {
      apdc_data.num_threads = aurostd::string2utype<int>(XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS"));
    }
    if (apdc_data.num_threads < 1) {apdc_data.num_threads = 1;}
    apdc_data.workdirpath = aurostd::getPWD();
    apdc_data.rootdirpath = aurostd::CleanFileName(apdc_data.rootdirpath);
    apdc_data.plattice = aurostd::tolower(apdc_data.plattice);
    aurostd::sort_remove_duplicates(apdc_data.elements);
    ErrorChecks(apdc_data);
    // Binodal
    apdc_data.rundirpath += apdc_data.rootdirpath + "/" + pflow::arity_string(apdc_data.elements.size(), false, false) + "/" + apdc_data.plattice + "/" + apdc_data.alloyname;
    aurostd::DirectoryMake(apdc_data.rundirpath);
    GetBinodal(apdc_data);
  }

 void GetPhaseDiagram(const string& aflowin, bool elements_only) {
    _apdc_data apdc_data;
    if (elements_only) { // command line call
      aurostd::string2tokens(aflowin, apdc_data.elements, ",");
      apdc_data.rootdirpath = aurostd::getPWD();
      apdc_data.plattice = DEFAULT_APDC_PLATTICE;
      apdc_data.aflow_max_num_atoms = DEFAULT_APDC_AFLOW_MAX_NUM_ATOMS;
      apdc_data.max_num_atoms = DEFAULT_APDC_MAX_NUM_ATOMS;
    }
    else {
      cerr<<aflowin<<endl;
    }
    GetPhaseDiagram(apdc_data);
 }

 void GetPhaseDiagram(istream& infile) {string aflowin; std::getline(infile, aflowin); GetPhaseDiagram(aflowin, false);}
}

// ***************************************************************************
// apdc::ErrorChecks
// ***************************************************************************
namespace apdc {
  void ErrorChecks(_apdc_data& apdc_data) {
    if (!aurostd::DirectoryMake(apdc_data.rootdirpath)) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "ErrorChecks():", "Cannot create directory", _FILE_ERROR_);
    }
    if (apdc_data.plattice != "fcc" && apdc_data.plattice != "bcc" && apdc_data.plattice != "hcp") {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "ErrorChecks():", "Invalid parent lattice", _INPUT_ILLEGAL_);
    }
    if (apdc_data.elements.size() < 2) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "ErrorChecks():", "Alloy must be at least binary", _VALUE_ERROR_);
    }
    for (uint i = 0; i < apdc_data.elements.size(); i++) {
      if (!xelement::xelement::isElement(apdc_data.elements[i])) {
        throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "ErrorChecks():", "Element \"" + apdc_data.elements[i] + "\" is invalid", _INPUT_ILLEGAL_);
      }
      apdc_data.alloyname += apdc_data.elements[i];
    }
  }
}

// ***************************************************************************
// apdc::GetBinodal
// ***************************************************************************
namespace apdc {
  void GetBinodal(_apdc_data& apdc_data) {
    apdc_data.lat_atat = CreateLatForATAT(apdc_data.plattice, apdc_data.elements);
    apdc_data.vstr_atat = GetATATXstructures(apdc_data.lat_atat, (uint)apdc_data.max_num_atoms);
    apdc_data.vstr_aflow = GetAFLOWXstructures(apdc_data.plattice, apdc_data.elements, apdc_data.num_threads);
 //   // map ATAT xstrs to AFLOW xstrs because ATAT cannot identify AFLOW xstrs
    apdc_data.mapstr = GetMapForXstructures(GetATATXstructures(apdc_data.lat_atat, apdc_data.aflow_max_num_atoms), apdc_data.vstr_aflow, apdc_data.num_threads);
    GenerateFilesForATAT(apdc_data.rundirpath, apdc_data.lat_atat, apdc_data.vstr_aflow, apdc_data.vstr_atat, apdc_data.mapstr);
    RunATAT(apdc_data.workdirpath, apdc_data.rundirpath);
    apdc_data.conc = GetConcentration(apdc_data.rundirpath, apdc_data.vstr_atat.size(), apdc_data.elements.size());
    vector<xvector<int> > multiplicity = GetMultiplicity(apdc_data.vstr_atat);
    apdc_data.multiplicity = multiplicity[1];
    apdc_data.excess_energies = GetExcessEnergy(apdc_data.rundirpath, apdc_data.conc, multiplicity[0]);
    apdc_data.prob_rand = GetProbabilityRandom(apdc_data.conc, apdc_data.multiplicity, apdc_data.max_num_atoms);
    
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
// apdc::GetProbabilityRandom
// ***************************************************************************
// P(X) = g*(Xa^Na)*(Xb*Nb)*(Xc*Nc)*...
namespace apdc {
  xvector<double> GetProbabilityRandom(const xmatrix<double>& conc, const xvector<int>& multiplicity, const int max_num_atoms) {
    int nstr = conc.rows, nelem = conc.cols;
    xvector<double> prob_rand(nstr);
    for (int i = 1; i <= nstr; i++) {
      prob_rand(i) = multiplicity(i);
      for (int j = 1; j <= nelem; j++) {prob_rand(i) *= std::pow(conc(i, j), max_num_atoms);}
    }
    return prob_rand / aurostd::sum(prob_rand);
  }
}

// ***************************************************************************
// apdc::GetMultiplicity
// ***************************************************************************
// N = Na+Nb+Nc+...
// g = N!/(Na!*Nb!*Nc!*...)
// return [N,g]
namespace apdc {
  vector<xvector<int> > GetMultiplicity(const vector<xstructure>& vstr) {
    int nstr = vstr.size();
    vector<xvector<int> > multiplicity;
    xvector<int> g(nstr), n(nstr);
    uint natom, fact_prod;
    for (int i = 1; i <= nstr; i++) {
      natom = 0;
      fact_prod = 1;
      for (uint j = 0; j < vstr[i - 1].num_each_type.size(); j++) {
        natom += vstr[i - 1].num_each_type[j];
        fact_prod *= aurostd::factorial(vstr[i - 1].num_each_type[j]);
      }
      n(i) = natom;
      g(i) = aurostd::factorial(natom) / fact_prod;
    }
    multiplicity.push_back(n); multiplicity.push_back(g);
    return multiplicity;
  }
}

// ***************************************************************************
// apdc::GetConcentration
// ***************************************************************************
namespace apdc {
  xmatrix<double> GetConcentration(const vector<string>& elements, const vector<xstructure>& vstr) {
    uint nstr = vstr.size(), nelem = elements.size();
    xmatrix<double> conc(nstr, nelem);
    int ie = -1;
    xvector<double> stoich;
    vector<string> str_elements;
    for (uint i = 0; i < vstr.size(); i++) {
      if (nelem != vstr[i].stoich_each_type.size()) {
        str_elements = vstr[i].GetElements(true, true);
        stoich = 0.0 * aurostd::ones_xv<double>(nelem);
        for (uint j = 0; j < nelem; j++) {
          if (aurostd::WithinList(str_elements, elements[j], ie)) {stoich(j + 1) = vstr[i].stoich_each_type[ie];}
        }
      }
      else {
        stoich = aurostd::vector2xvector(aurostd::deque2vector(vstr[i].stoich_each_type));
      }
      for (uint j = 0; j < nelem; j++) {conc(i + 1, j + 1) = stoich(j + 1);}
    }
    return conc;
  }

  xmatrix<double> GetConcentration(const string& rundirpath, const int nstr, const int nelem) {
    xmatrix<double> conc(nstr, nelem);
    vector<string> vinput, tokens;
    aurostd::file2vectorstring(rundirpath + "/fit.out", vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      for (int i = 1; i <= nelem; i++) {
        conc(aurostd::string2utype<int>(tokens[tokens.size() - 1]) + 1, i) = aurostd::string2utype<double>(tokens[i - 1]);
      }
    }
    aurostd::file2vectorstring(rundirpath + "/predstr.out", vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      for (int i = 1; i <= nelem; i++) {
        conc(aurostd::string2utype<int>(tokens[tokens.size() - 2]) + 1, i) = aurostd::string2utype<double>(tokens[i - 1]);
      }
    }
    return conc;
  }
}

// ***************************************************************************
// apdc::GetExcessEnergy
// ***************************************************************************
namespace apdc {
  xvector<double> GetExcessEnergy(const string& rundirpath, const xmatrix<double>& conc, const xvector<int>& natom) {
    int nstr = conc.rows, nelem = conc.cols;
    vector<string> vinput, tokens;
    aurostd::file2vectorstring(rundirpath + "/ref_energy.out", vinput);
    xvector<double> nrg_ref = aurostd::vector2xvector(aurostd::vectorstring2vectordouble(vinput));
    xvector<double> nrg(nstr);
    aurostd::file2vectorstring(rundirpath + "/fit.out", vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      nrg(aurostd::string2utype<int>(tokens[tokens.size() - 1]) + 1) = aurostd::string2utype<double>(tokens[nelem]);
    }
    aurostd::file2vectorstring(rundirpath + "/predstr.out", vinput);
    for (uint line = 0; line < vinput.size(); line++) {
      aurostd::string2tokens(vinput[line], tokens, " ");
      nrg(aurostd::string2utype<int>(tokens[tokens.size() - 2]) + 1) = aurostd::string2utype<double>(tokens[nelem]);
    }
    // Convert D_j = E_j/n_j - \sum(X_k*E0_k) ==> D_j = E_j - \sum(X_k*E0_k)
    double nrg_macro;
    for (int i = 1; i <= nstr; i++) {
      nrg_macro = aurostd::sum(conc.getmat(i, i, 1, nelem, 1, 1) * nrg_ref);
      nrg(i) += nrg_macro;
      nrg(i) *= natom(i);
      nrg(i) -= nrg_macro;
    }
    return nrg;
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
    aurostd::RemoveFile(rundirpath + "/maps.log");
    chdir(rundirpath.c_str());
    string tmpfile = aurostd::TmpStrCreate(), logstring = "";
    aurostd::execute("mmaps -d > " + tmpfile + " 2>&1 &");
    while (!aurostd::substring2bool(logstring, "true and predicted ground states agree") && !aurostd::substring2bool(logstring, "No other ground states")) {
      iter++;
      if (LDEBUG) {cerr << "Sleeping, iter=" << iter << endl;}
      if (iter > 30) { // wait 30 times the minimum sleep (60 seconds)
        aurostd::RemoveFile(tmpfile);
        throw aurostd::xerror(_AFLOW_FILE_NAME_, XPID + "RunATAT():", "mmaps is taking too long to predict structures, dir=" + rundirpath, _RUNTIME_ERROR_);
      }
      aurostd::Sleep(_APDC_MIN_SLEEP_);
      aurostd::file2string(rundirpath + "/maps.log", logstring);
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
// apdc::GetAFLOWXstructures
// ***************************************************************************
namespace apdc {
  vector<xstructure> GetAFLOWXstructures(const string& plattice, const vector<string>& elements, const int num_proc, bool keep_all) {
    vector<xstructure> vstr;
    vector<string> vstrlabel;
    string aflowlib, aflowurl;
    string alloyname = "";
    aflowlib::_aflowlib_entry entry;
    uint nelem = elements.size();
    stringstream oss;
    for (uint i = 0; i < nelem; i++) {alloyname += AVASP_Get_PseudoPotential_PAW_PBE(elements[i]);}
    aflowlib = "/common/LIB" + aurostd::utype2string<uint>(nelem) + "/RAW/" + alloyname;
    aflowurl = "aflowlib.duke.edu:AFLOWDATA/LIB" + aurostd::utype2string<uint>(nelem) + "_RAW/" + alloyname;
    if (plattice == "fcc") {
      if (nelem >= 2) {
        for (uint i = 1; i < 30; i++) {vstrlabel.push_back(aurostd::utype2string<uint>(i));}
      }
      if (nelem >= 3) {
        vstrlabel.push_back("TFCC001.ABC");vstrlabel.push_back("TFCC002.ABC");vstrlabel.push_back("TFCC003.ABC");
        for (uint i = 4; i < 17; i++) {
          vstrlabel.push_back("TFCC00" + aurostd::utype2string<uint>(i) + ".ABC");
          vstrlabel.push_back("TFCC00" + aurostd::utype2string<uint>(i) + ".BCA");
          vstrlabel.push_back("TFCC00" + aurostd::utype2string<uint>(i) + ".CAB");
        }
      }
    }
    else if (plattice == "bcc") {
      if (nelem >= 2) {
        for (uint i = 58; i < 87; i++) {vstrlabel.push_back(aurostd::utype2string<uint>(i));}
      }
      if (nelem >= 3) {
        vstrlabel.push_back("TBCC001.ABC");vstrlabel.push_back("TBCC002.ABC");vstrlabel.push_back("TBCC003.ABC");
        for (uint i = 4; i < 17; i++) {
          vstrlabel.push_back("TBCC00" + aurostd::utype2string<uint>(i) + ".ABC");
          vstrlabel.push_back("TBCC00" + aurostd::utype2string<uint>(i) + ".BCA");
          vstrlabel.push_back("TBCC00" + aurostd::utype2string<uint>(i) + ".CAB");
        }
      }
    }
    else if (plattice == "hcp") {
      if (nelem >= 2) {
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
    if (sstr.size() == 0 || aurostd::substring2bool(sstr, "Unable to open lattice file")) {
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
  vector<int> GetMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2, const int num_proc) {
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
