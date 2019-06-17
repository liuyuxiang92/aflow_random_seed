// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2018                    *
// *                                                                         *
// ***************************************************************************

// This file provides serves as an AflowIn reader for AFLOW modules (e.g. APL).
// It loads all parameters from the aflow.in file (or the defaults if no file
// is present).

#include "aflow.h"

#define _ASTROPT_APL_OLD_ string("[AFLOW_PHONONS]")
#define _ASTROPT_APL_ string("[AFLOW_APL]")
#define _ASTROPT_QHA_ string("[AFLOW_QHA]")
#define _ASTROPT_AAPL_ string("[AFLOW_AAPL]")

#define FLAG_PRECISION 7  //CO181226

#define DEBUG_MODULES false

using aurostd::utype2string;

namespace KBIN {

//setModules//////////////////////////////////////////////////////////////////
// Set all module flags to their default values. This is used when no
// aflow.in file is availabe.
void setModules(_xvasp& xvasp) {
  _xinput xinput(xvasp);
  setModules(xinput);
  xvasp = xinput.xvasp;
}
  
// General case
void setModules(_xinput& xinput) {
  _moduleOptions module_opts;
  module_opts.aplflags = loadDefaultsAPL();
  module_opts.aaplflags = loadDefaultsAAPL();
  // The readParameters functions are necessary to set xvasp
  string placeholder = "";  // acts as pseudo-aflow.in
  readParametersAPL(placeholder, module_opts, xinput);
  readParametersAAPL(placeholder, module_opts, xinput);
}

//readModulesfromAflowIn//////////////////////////////////////////////////////
// Reads all module flags from an AflowIn file.
void readModulesFromAflowIn(const string& AflowIn,
                            _kflags& kflags, _xvasp& xvasp) {
  _xinput xinput(xvasp);
  readModulesFromAflowIn(AflowIn, kflags, xinput);
  xvasp = xinput.xvasp;
}

// General case
void readModulesFromAflowIn(const string& AflowIn,
                            _kflags& kflags, _xinput& xinput) {
  _moduleOptions module_opts;
  module_opts.aplflags = loadDefaultsAPL();
  module_opts.aaplflags = loadDefaultsAAPL();
  readParametersAPL(AflowIn, module_opts, xinput);
  readParametersAAPL(AflowIn, module_opts, xinput);
  kflags.KBIN_MODULE_OPTIONS = module_opts;
}                 

// APL-related functions -----------------------------------------------------

//loadDefaultsAPL/////////////////////////////////////////////////////////////
// Sets all APL flags to their default values.
vector<aurostd::xoption> loadDefaultsAPL() {
  bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
  string soliloquy="loadDefaultsAPL():";
  vector<aurostd::xoption> aplflags;
  aurostd::xoption opt;
  opt.keyword="RELAX"; opt.option = DEFAULT_APL_RELAX; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
  opt.keyword="HIBERNATE"; opt.option = DEFAULT_APL_HIBERNATE; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
  opt.keyword="ENGINE"; opt.xscheme = DEFAULT_APL_ENGINE; aplflags.push_back(opt); opt.clear();
  opt.keyword="SUPERCELL"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
  opt.keyword="MINATOMS"; opt.xscheme = utype2string<int>(DEFAULT_APL_MINATOMS); aplflags.push_back(opt); opt.clear();
  opt.keyword="MINATOMS_RESTRICTED"; opt.xscheme = utype2string<int>(DEFAULT_APL_MINATOMS); aplflags.push_back(opt); opt.clear();
  opt.keyword="MINSHELL"; opt.xscheme = utype2string<int>(DEFAULT_APL_MINSHELL); aplflags.push_back(opt); opt.clear();
  opt.keyword="POLAR"; opt.option = DEFAULT_APL_POLAR; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
  opt.keyword="DMAG"; opt.xscheme = utype2string<double>(DEFAULT_APL_DMAG, FLAG_PRECISION); aplflags.push_back(opt); opt.clear();
  opt.keyword="DXYZONLY"; opt.option = DEFAULT_APL_DXYZONLY; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
  opt.keyword="DSYMMETRIZE"; opt.option = DEFAULT_APL_DSYMMETRIZE; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
  opt.keyword="DINEQUIV_ONLY"; opt.option = DEFAULT_APL_DINEQUIV_ONLY; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear(); //CO190131
  //[ME181226 - now a default in .aflow.rc]// Special case: DPM can be true, false, or empty
  opt.keyword="DPM"; opt.xscheme = DEFAULT_APL_DPM; opt.option = (opt.xscheme=="ON"?true:false); aplflags.push_back(opt); opt.clear();  //CO181226
  //[ME181226 - now a default in .aflow.rc]// Special case: k-points options can be empty
  opt.keyword="KPPRA"; opt.xscheme = utype2string<int>(DEFAULT_PHONONS_KPPRA); aplflags.push_back(opt); opt.clear(); //CO181226 // ME 190112
  opt.keyword="KSCHEME"; opt.xscheme = DEFAULT_PHONONS_KSCHEME; aplflags.push_back(opt); opt.clear();  // ME 190109 - KPPRA can be taken from STATIC, but KSCHEME should default to G
  opt.keyword="KPOINTS"; aplflags.push_back(opt); opt.clear();
  opt.keyword="PREC"; opt.xscheme = DEFAULT_APL_PREC; aplflags.push_back(opt); opt.clear();
  opt.keyword="ZEROSTATE"; opt.option = DEFAULT_APL_ZEROSTATE; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
  opt.keyword="FREQFORMAT"; opt.xscheme = DEFAULT_APL_FREQFORMAT; aplflags.push_back(opt); opt.clear();
  opt.keyword="DC"; opt.option = DEFAULT_APL_DC; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
  opt.keyword="DCPOINTS"; opt.xscheme = utype2string<int>(DEFAULT_APL_DCPOINTS); aplflags.push_back(opt); opt.clear();
  opt.keyword="DCPATH"; opt.xscheme = DEFAULT_APL_DCPATH; aplflags.push_back(opt); opt.clear();
  opt.keyword="DCINITCOORDSFRAC"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
  opt.keyword="DCINITCOORDSCART"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
  opt.keyword="DCINITCOORDSLABELS"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
  opt.keyword="DCUSERPATH"; opt.xscheme = ""; aplflags.push_back(opt); opt.clear();
  opt.keyword="DOS"; opt.option = DEFAULT_APL_DOS; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
  opt.keyword="DOSMETHOD"; opt.xscheme = DEFAULT_APL_DOSMETHOD; aplflags.push_back(opt); opt.clear();
  opt.keyword="DOSMESH"; opt.xscheme = DEFAULT_APL_DOSMESH; aplflags.push_back(opt); opt.clear();
  opt.keyword="DOSSMEAR"; opt.xscheme = utype2string<double>(DEFAULT_APL_DOSSMEAR, FLAG_PRECISION); aplflags.push_back(opt); opt.clear();
  opt.keyword="DOSPOINTS"; opt.xscheme = utype2string<int>(DEFAULT_APL_DOSPOINTS); aplflags.push_back(opt); opt.clear();
  opt.keyword="TP"; opt.option = DEFAULT_APL_TP; opt.xscheme = (opt.option?"ON":"OFF"); aplflags.push_back(opt); opt.clear();
  opt.keyword="TPT"; opt.xscheme = DEFAULT_APL_TPT; aplflags.push_back(opt); opt.clear();
  if (LDEBUG) {
    for (uint i = 0; i < aplflags.size(); i++) {
      std::cerr << soliloquy << " key: " << aplflags[i].keyword << ", xscheme: " << aplflags[i].xscheme << ", option: " << aplflags[i].option << std::endl;
    }
  }
  return aplflags;
}

//writeFlagAPL///////////////////////////////////////////////////////////////
// Determines whether flag should be written to aflow.in
// CO181226
bool writeFlagAPL(const string& key,const xoption& xopt){  // ME190113
  // return true;  OBSOLETE ME190113
  if (xopt.isentry) {return true;}  // ME190116 - Do not remove user entries
  if(key=="RELAX"){return true;}
  if(key=="HIBERNATE"){if((xopt.option == AFLOWRC_DEFAULT_APL_HIBERNATE) && (xopt.option == DEFAULT_APL_HIBERNATE)) {return false;}}  // ME190113
  if(key=="ENGINE"){return true;}
  if(key=="SUPERCELL"){return true;}
  if(key=="MINATOMS"){return true;}
  if(key=="MINSHELL"){return true;}
  if(key=="POLAR"){return true;}
  if(key=="DMAG"){return true;}
  if(key=="DXYZONLY"){if((xopt.option == AFLOWRC_DEFAULT_APL_DXYZONLY) && (xopt.option == DEFAULT_APL_DXYZONLY)) {return false;}}  // ME190113
  if(key=="DSYMMETRIZE"){if((xopt.option == AFLOWRC_DEFAULT_APL_DSYMMETRIZE) && (xopt.option == DEFAULT_APL_DSYMMETRIZE)) {return false;}}  // ME190113
  if(key=="DINEQUIV_ONLY"){if((xopt.option == AFLOWRC_DEFAULT_APL_DINEQUIV_ONLY) && (xopt.option == DEFAULT_APL_DINEQUIV_ONLY)) {return false;}}  // ME190113
  if(key=="DPM"){if(AFLOWRC_DEFAULT_APL_DPM==xopt.xscheme && DEFAULT_APL_DPM==xopt.xscheme){return false;}}  // ME190113
  if(key=="PREC"){if(AFLOWRC_DEFAULT_APL_PREC==xopt.xscheme && DEFAULT_APL_PREC==xopt.xscheme){return false;}}  // ME190113
  if(key=="ZEROSTATE"){return true;}
  if(key=="FREQFORMAT"){if(AFLOWRC_DEFAULT_APL_FREQFORMAT==xopt.xscheme && DEFAULT_APL_FREQFORMAT==xopt.xscheme){return false;}}  // ME190113
  if(key=="DC"){return true;}
  if(key=="DCPOINTS"){if(utype2string<int>(AFLOWRC_DEFAULT_APL_DCPOINTS)==xopt.xscheme && utype2string<int>(DEFAULT_APL_DCPOINTS)==xopt.xscheme){return false;}}  // ME190113
  if(key=="DCPATH"){if(AFLOWRC_DEFAULT_APL_DCPATH==xopt.xscheme && DEFAULT_APL_DCPATH==xopt.xscheme){return false;}}  // ME190113
  if(key=="DOS"){return true;}
  if(key=="DOSMETHOD"){if(AFLOWRC_DEFAULT_APL_DOSMETHOD==xopt.xscheme && DEFAULT_APL_DOSMETHOD==xopt.xscheme){return false;}}  // ME190113
  if(key=="DOSMESH"){return true;}  // ME190113 - should always write
  if(key=="DOSSMEAR"){if(utype2string<double>(AFLOWRC_DEFAULT_APL_DOSSMEAR, FLAG_PRECISION)==xopt.xscheme && utype2string<double>(DEFAULT_APL_DOSSMEAR, FLAG_PRECISION)==xopt.xscheme){return false;}}  // ME190113
  if(key=="DOSPOINTS"){if(utype2string<int>(AFLOWRC_DEFAULT_APL_DOSPOINTS)==xopt.xscheme && utype2string<int>(DEFAULT_APL_DOSPOINTS)==xopt.xscheme){return false;}}  // ME190113
  if(key=="TP"){return true;}
  if(key=="TPT"){return true;}
  return true;
}

//readParametersAPL///////////////////////////////////////////////////////////
// Reads APL flags from an aflow.in file.
void readParametersAPL(const string& AflowIn,
                       _moduleOptions& module_opts, _xinput& xinput) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
  string soliloquy="readParametersAPL():";
  string key, entry, xvaspflag;
  module_opts.supercell_method.assign(4, false);
  for (uint i = 0; i < module_opts.aplflags.size(); i++) {
    key = module_opts.aplflags[i].keyword;
    entry = _ASTROPT_APL_ + key + "=|" + _ASTROPT_APL_OLD_ + key + "="; //CO181226
    entry += "|" + _ASTROPT_AAPL_ + key + "=|" + _ASTROPT_QHA_ + key + "="; //CO181226
    if (key == "DMAG") { // for backwards compatibility
      entry += "|" + _ASTROPT_APL_ + "DISMAG"+ "=|" + _ASTROPT_APL_OLD_ + "DISMAG" + "="; //CO181226
      entry += "|" + _ASTROPT_AAPL_ + "DISMAG" + "=|" + _ASTROPT_QHA_ + "DISMAG" + "="; //CO181226
      entry += "|" + _ASTROPT_APL_ + "TDMAG"+ "=|" + _ASTROPT_APL_OLD_ + "TDMAG" + "="; //CO181226
      entry += "|" + _ASTROPT_AAPL_ + "TDMAG" + "="; //CO181226
      entry += "|" + _ASTROPT_APL_ + "TDISMAG"+ "=|" + _ASTROPT_APL_OLD_ + "TDISMAG" + "="; //CO181226
      entry += "|" + _ASTROPT_AAPL_ + "TDISMAG" + "="; //CO181226
    }
    if (key == "DSYMMETRIZE") {  // for backwards compatibility
      entry += "|" + _ASTROPT_APL_ + "SYMMETRIZE=";  //CO181226
      entry += "|" + _ASTROPT_APL_ + "SYM=";  //CO181226
    }
    if (key == "DINEQUIV_ONLY") {  // for backwards compatibility
      entry += "|" + _ASTROPT_APL_ + "INEQUIVONLY=";  //CO181226
    }
    module_opts.aplflags[i].options2entry(AflowIn, entry, module_opts.aplflags[i].option, module_opts.aplflags[i].xscheme);

    // options2entry sets the keyword to _ASTROPT_ + key + "=", so reset
    module_opts.aplflags[i].keyword = key;

    //[ME181226 - now a default in .aflow.rc]// Special case: set DPM to AUTO if not an entry
    //[ME181226 - now a default in .aflow.rc]if ((key == "DPM") && (!module_opts.aplflags[i].isentry)) module_opts.aplflags[i].xscheme = "AUTO";

    // Write xvasp
    if(xinput.AFLOW_MODE_VASP) {
      xvaspflag = "AFLOWIN_FLAG::APL_" + key;
      //[CO181226 - need to revise]if(writeFlagAPL(key,module_opts.aplflags[i].xscheme)){xinput.xvasp.aplopts.flag(xvaspflag, TRUE);}  //CO181226
      xinput.xvasp.aplopts.flag(xvaspflag, TRUE);  //CO181226
      xinput.xvasp.aplopts.push_attached(xvaspflag, module_opts.aplflags[i].xscheme); //this should become pop/push or changeattachedscheme (eventually)
    }

    // Supercell options
    if ((key == "SUPERCELL") && (module_opts.aplflags[i].isentry)) {
      module_opts.supercell_method[0] = true;
      continue;
    }
    if ((key == "MINATOMS") && (module_opts.aplflags[i].isentry)) {
      module_opts.supercell_method[1] = true;
      continue;
    }
    if ((key == "MINATOMS_RESTRICTED") && (module_opts.aplflags[i].isentry)) {
      module_opts.supercell_method[1] = true;
      module_opts.minatoms_restricted = true;
      continue;
    }
    // MAXSHELL will be skipped because it's not documented at all
    if ((key == "MINSHELL") && (module_opts.aplflags[i].isentry)) {
      module_opts.supercell_method[3] = true;
      continue;
    }
  }

  // Was a supercell entry found? If not, switch to MINATOMS
  bool supercell_found = false;
  for (uint i = 0; i < module_opts.supercell_method.size(); i++) {
    if (module_opts.supercell_method[i]) {
      supercell_found = true;
      break;
    }
  }
  if (!supercell_found) module_opts.supercell_method[1] = true;

  // Unset contradicting/unnecessary xvasp flags
  if (xinput.AFLOW_MODE_VASP) {
    // Engine-related parameters
    if (xinput.xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_ENGINE") == "LR") {
      // Unset DM parameters - do not unset DMAG because AAPL may need it
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DPM", false);
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DSYMMETRIZE", false); //CO190131
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DINEQUIV_ONLY", false); //CO190131
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DXYZONLY", false);
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_ZEROSTATE", false);
    } else {
      // DPM
      //[ME181226 - now a default in .aflow.rc]if (xinput.xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_DPM") == "AUTO") {
      //[ME181226 - now a default in .aflow.rc]  xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DPM", false);
      //[ME181226 - now a default in .aflow.rc]}
    }

    // Supercell
    xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_SUPERCELL", module_opts.supercell_method[0]);
    xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_MINATOMS", (module_opts.supercell_method[1] && !module_opts.minatoms_restricted));
    xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_MINATOMS_RESTRICTED", (module_opts.supercell_method[1] && module_opts.minatoms_restricted));
    xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_MINSHELL", module_opts.supercell_method[3]);

    // Phonon dispersion
    if (xinput.xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_DCPATH") == "LATTICE") {
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DCINITCOORDSFRAC", false);
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DCINITCOORDSCART", false);
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DCINITCOORDSLABELS", false);
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DCUSERPATH", false);
    }

    // DOS
    if (xinput.xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_DOSMETHOD") == "LT") {
      xinput.xvasp.aplopts.flag("AFLOWIN_FLAG::APL_DOSSMEAR", false);
    }
  }
  if (LDEBUG) {
    for (uint i = 0; i < module_opts.aplflags.size(); i++) {
      std::cerr << soliloquy << "  " << module_opts.aplflags[i].keyword << " = " << module_opts.aplflags[i].xscheme << std::endl;
    }
  }
}

// AAPL-related functions ----------------------------------------------------

//loadDefaultsAAPL////////////////////////////////////////////////////////////
// Sets all AAPL flags to their default values.
vector<aurostd::xoption> loadDefaultsAAPL() {
  bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
  string soliloquy="loadDefaultsAAPL():";
  vector<aurostd::xoption> aaplflags;
  aurostd::xoption opt;
  opt.keyword="BTE"; opt.xscheme = DEFAULT_AAPL_BTE; aaplflags.push_back(opt); opt.clear();
  opt.keyword="FOURTH_ORDER"; opt.option = DEFAULT_AAPL_FOURTH_ORDER; opt.xscheme = (opt.option?"ON":"OFF"); aaplflags.push_back(opt); opt.clear();
  opt.keyword="CUT_RAD"; opt.xscheme = DEFAULT_AAPL_CUT_RAD; aaplflags.push_back(opt); opt.clear();
  opt.keyword="CUT_SHELL"; opt.xscheme = DEFAULT_AAPL_CUT_SHELL; aaplflags.push_back(opt); opt.clear();
  opt.keyword="THERMALGRID"; opt.xscheme = DEFAULT_AAPL_THERMALGRID; aaplflags.push_back(opt); opt.clear();
  opt.keyword="TCT"; opt.xscheme = DEFAULT_AAPL_TCT; aaplflags.push_back(opt); opt.clear();
  opt.keyword="SUMRULE"; opt.xscheme = utype2string<double>(DEFAULT_AAPL_SUMRULE, FLAG_PRECISION); aaplflags.push_back(opt); opt.clear();
  opt.keyword="SUMRULE_MAX_ITER"; opt.xscheme = utype2string<int>(DEFAULT_AAPL_SUMRULE_MAX_ITER); aaplflags.push_back(opt); opt.clear();
  opt.keyword="MIXING_COEFFICIENT"; opt.xscheme = utype2string<double>(DEFAULT_AAPL_MIXING_COEFFICIENT, FLAG_PRECISION); aaplflags.push_back(opt); opt.clear();
  opt.keyword="ISOTOPE"; opt.option = DEFAULT_AAPL_ISOTOPE; opt.xscheme = (opt.option?"ON":"OFF"); aaplflags.push_back(opt); opt.clear();
  opt.keyword="BOUNDARY"; opt.option = DEFAULT_AAPL_BOUNDARY; opt.xscheme = (opt.option?"ON":"OFF"); aaplflags.push_back(opt); opt.clear();
  opt.keyword="CUMULATIVEK"; opt.option = DEFAULT_AAPL_CUMULATIVEK; opt.xscheme = (opt.option?"ON":"OFF"); aaplflags.push_back(opt); opt.clear();
  opt.keyword="NANO_SIZE"; opt.xscheme = utype2string<double>(DEFAULT_AAPL_NANO_SIZE, FLAG_PRECISION); aaplflags.push_back(opt); opt.clear();
  if (LDEBUG) {
    for (uint i = 0; i < aaplflags.size(); i++) {
      std::cerr << soliloquy << " key: " << aaplflags[i].keyword << ", xscheme: " << aaplflags[i].xscheme << ", option: " << aaplflags[i].option << std::endl;
    }
  }
  return aaplflags;
}

//writeFlagAAPL///////////////////////////////////////////////////////////////
// Determines whether flag should be written to aflow.in
// CO181226
bool writeFlagAAPL(const string& key,const xoption& xopt){
  // return true; OBSOLETE ME190113
  if (xopt.isentry) {return true;}  // ME190116 - Do not remove user entries
  if(key=="BTE"){return true;}
  if(key=="FOURTH_ORDER"){return true;}  // ME190113 - should always write to be explicit
  if(key=="CUT_RAD"){return true;}
  if(key=="CUT_SHELL"){return true;}
  if(key=="THERMALGRID"){return true;}
  if(key=="TCT"){return true;}
  if(key=="SUMRULE"){if(utype2string<double>(AFLOWRC_DEFAULT_AAPL_SUMRULE, FLAG_PRECISION)==xopt.xscheme && utype2string<double>(DEFAULT_AAPL_SUMRULE, FLAG_PRECISION)==xopt.xscheme){return false;}}  // ME190113
  if(key=="SUMRULE_MAX_ITER"){if(utype2string<int>(AFLOWRC_DEFAULT_AAPL_SUMRULE_MAX_ITER)==xopt.xscheme && utype2string<int>(DEFAULT_AAPL_SUMRULE_MAX_ITER)==xopt.xscheme){return false;}}  // ME190113
  if(key=="MIXING_COEFFICIENT"){if(utype2string<double>(AFLOWRC_DEFAULT_AAPL_MIXING_COEFFICIENT, FLAG_PRECISION)==xopt.xscheme && utype2string<double>(DEFAULT_AAPL_MIXING_COEFFICIENT, FLAG_PRECISION)==xopt.xscheme){return false;}}  // ME190113
  if(key=="ISOTOPE"){return true;}
  if(key=="BOUNDARY"){return true;}
  if(key=="CUMULATIVEK"){if((xopt.option == AFLOWRC_DEFAULT_AAPL_CUMULATIVEK) && (xopt.option == DEFAULT_AAPL_CUMULATIVEK)) {return false;}}  // ME190113
  if(key=="NANO_SIZE"){return true;}
  return true;
}

//readParametersAAPL//////////////////////////////////////////////////////////
// Reads AAPL flags from an aflow.in file.
void readParametersAAPL(const string& AflowIn,
                         _moduleOptions& module_opts, _xinput& xinput) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || DEBUG_MODULES);
  string soliloquy="readParametersAAPL():";
  string key, entry, xvaspflag;
  module_opts.cut_rad_shell.assign(2, false);
  for (uint i = 0; i < module_opts.aaplflags.size(); i++) {
    key = module_opts.aaplflags[i].keyword;
    entry = _ASTROPT_AAPL_ + key + "=|" + _ASTROPT_APL_OLD_ + key + "=";  //CO181226
    module_opts.aaplflags[i].options2entry(AflowIn, entry, module_opts.aaplflags[i].option, module_opts.aaplflags[i].xscheme);
    // options2entry sets the keyword to _ASTROPT_ + key + "=", so reset
    module_opts.aaplflags[i].keyword = key;
    if (xinput.AFLOW_MODE_VASP) {
      xvaspflag = "AFLOWIN_FLAG::AAPL_" + key;
      //[CO181226 - need to revise]if(writeFlagAAPL(key,module_opts.aaplflags[i].xscheme)){xinput.xvasp.aaplopts.flag(xvaspflag, TRUE);}  //CO181226
      xinput.xvasp.aaplopts.flag(xvaspflag, TRUE);  //CO181226
      xinput.xvasp.aaplopts.push_attached(xvaspflag, module_opts.aaplflags[i].xscheme); //this should become pop/push or changeattachedscheme (eventually)
    }
    // Special rules for certain keywords
    if (key == "CUT_RAD") {
      module_opts.cut_rad_shell[0] = module_opts.aaplflags[i].isentry;
      continue;
    }
    if (key == "CUT_SHELL") {
      module_opts.cut_rad_shell[1] = module_opts.aaplflags[i].isentry;
      continue;
    }
  }
  if (module_opts.cut_rad_shell[0] != module_opts.cut_rad_shell[1]) {
    if (xinput.AFLOW_MODE_VASP) {
      xinput.xvasp.aaplopts.flag("AFLOWIN_FLAG::AAPL_CUT_SHELL", module_opts.cut_rad_shell[0]);
      xinput.xvasp.aaplopts.flag("AFLOWIN_FLAG::AAPL_CUT_RAD", module_opts.cut_rad_shell[1]);
    }
  }
  if (LDEBUG) {
    for (uint i = 0; i < module_opts.aaplflags.size(); i++) {
      std::cerr << soliloquy << "  " << module_opts.aaplflags[i].keyword << " = " << module_opts.aaplflags[i].xscheme << std::endl;
    }
  }
}

}  // namespace KBIN

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2018                    *
// *                                                                         *
//****************************************************************************
