// ***************************************************************************
// *                                                                         *
// *             STEFANO CURTAROLO - Duke University 2003-2019               *
// *                                                                         *
// ***************************************************************************

#include "../aflow.h"
#include "aflow_apl.h"

#define _ASTROPT_APL_OLD_ string("[AFLOW_PHONONS]") //CO 170601, ensure backwards compatibility (we ALWAYS support LEGACY)
#define _ASTROPT_APL_ string("[AFLOW_APL]") //CO 170601
#define _ASTROPT_QHA_ string("[AFLOW_QHA]") //CO 170601
#define _ASTROPT_AAPL_ string("[AFLOW_AAPL]") //CO 170601
//temporary directory for storing QHA files
#define _TMPDIR_ string("ARUN.APL.QH.TMP")  //[PINKU]
#define _STROPT_ string("[VASP_FORCE_OPTION]") //ME181226

 //CO181226
#define COMPILE_SLIM

bool _WITHIN_DUKE_ = false;

//CO fixing cpp version issues with auto_ptr (depreciated)
#if __cplusplus >= 201103L
template <typename T>
using auto_ptr = std::unique_ptr<T>;
#else
using std::auto_ptr;
#endif

static const string _ANHARMONIC_IFCS_FILE_[2] = {"anharmonicIFCs_3rd.xml", "anharmonicIFCs_4th.xml"};
static const string _CLUSTER_SET_FILE_[2] = {"clusterSet_3rd.xml", "clusterSet_4th.xml"};
static const int _NUM_RELAX_ = 2; //ME181226
static const string _APL_RELAX_PREFIX_ = "relax_apl"; //ME181226  // ME190125

namespace apl {
  // ME190119
  // Old apl::tokenize is essentially just aurostd::string2tokens
  void tokenize(const string& strin, vector<string>& tokens, const string& del) {
    aurostd::string2tokens(strin, tokens, del);
  }
}

namespace KBIN {
// ME 181107 - Relax structure with PREC=PHONONS before running APL
bool relaxStructureAPL_VASP(int start_relax,
                            const string& AflowIn,
                             _xvasp& xvasp,
                             _aflags& aflags,
                             _kflags& kflags,
                             _vflags& vflags,
                             ofstream& fileMessage) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  ostringstream aus;
  // Store original settings
  string prec = vflags.KBIN_VASP_FORCE_OPTION_PREC.xscheme;
  bool prec_entry = vflags.KBIN_VASP_FORCE_OPTION_PREC.isentry;
  string kppra = vflags.KBIN_VASP_KPOINTS_KPPRA.xscheme;
  bool kppra_entry = vflags.KBIN_VASP_KPOINTS_KPPRA.isentry;
  string kscheme = vflags.KBIN_VASP_KPOINTS_KSCHEME.xscheme;
  bool kscheme_entry = vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry;
  string kpoints = vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.xscheme;
  bool kpoints_entry = vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.isentry;

  // Change to phonon settings
  string prec_phonons = xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_PREC");
  vflags.KBIN_VASP_FORCE_OPTION_PREC.pop(prec);
  vflags.KBIN_VASP_FORCE_OPTION_PREC.push(prec_phonons);
  vflags.KBIN_VASP_FORCE_OPTION_PREC.isentry = true;

  string kppra_phonons = xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_KPPRA");
  if (!kppra_phonons.empty()) {
    vflags.KBIN_VASP_KPOINTS_KPPRA.pop(kppra);
    vflags.KBIN_VASP_KPOINTS_KPPRA.push(kppra_phonons);
    vflags.KBIN_VASP_KPOINTS_KPPRA.isentry = true;
  }

  string kscheme_phonons = xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_KSCHEME");
  if (!kscheme_phonons.empty()) {
    vflags.KBIN_VASP_KPOINTS_KSCHEME.pop(kscheme);
    vflags.KBIN_VASP_KPOINTS_KSCHEME.push(kscheme_phonons);
    vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry = true;
  }

  string kpoints_phonons = xvasp.aplopts.getattachedscheme("AFLOWIN_FLAG::APL_KPOINTS");
  if (!kpoints_phonons.empty()) {
    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.pop(kpoints);
    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.push(kpoints_phonons);
    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.isentry = true;
  }
  
  bool Krun = VASP_Produce_and_Modify_INPUT(xvasp, AflowIn, fileMessage, aflags, kflags, vflags);
  Krun = (Krun && VASP_Write_INPUT(xvasp, vflags));

  if (Krun) {
    int i;
    for (i = start_relax; Krun && i <= _NUM_RELAX_; i++) {
      aus << 11111*i << " RELAXATION APL - " << xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
      aurostd::PrintMessageStream(fileMessage, aus, XHOST.QUIET);
      if (i < _NUM_RELAX_) {
        Krun = VASP_Run(xvasp, aflags, kflags, vflags, _APL_RELAX_PREFIX_ + aurostd::utype2string<int>(i), _APL_RELAX_PREFIX_ + aurostd::utype2string<int>(i), true, fileMessage);
        XVASP_INCAR_SPIN_REMOVE_RELAX(xvasp, aflags, vflags, i, fileMessage);
        XVASP_KPOINTS_IBZKPT_UPDATE(xvasp, aflags, vflags, i, fileMessage);
      } else { 
        Krun = VASP_Run(xvasp, aflags, kflags, vflags, _APL_RELAX_PREFIX_ + aurostd::utype2string<int>(i), true, fileMessage);
      }
      }
      if (Krun && (i == _NUM_RELAX_)) {
        aus << 11111*i << " RELAXATION APL END - " << xvasp.Directory << " - K=[" << xvasp.str.kpoints_k1 << " " << xvasp.str.kpoints_k2 << " " << xvasp.str.kpoints_k3 << "]" << " - " << kflags.KBIN_BIN << " - " << Message("user,host,time") << endl;
      }
  }

  // Update structure - do not set xvasp.str = str_fin or all other
  // information (species, PPs) will be deleted!
  stringstream xstr;
  string filename = "CONTCAR." + _APL_RELAX_PREFIX_ + aurostd::utype2string<int>(_NUM_RELAX_);  // ME190308
  if (!aurostd::FileExist(filename)) return false;
  aurostd::file2stringstream(filename, xstr);
  xstructure str_fin(xstr, IOVASP_AUTO);
  // Restore species and PP information
  str_fin.species = xvasp.str.species; // ME190109
  str_fin.species_pp = xvasp.str.species_pp; // ME190109
  str_fin.species_pp_type = xvasp.str.species_pp_type; // ME190109
  str_fin.species_pp_version = xvasp.str.species_pp_version; // ME190109
  str_fin.species_pp_ZVAL = xvasp.str.species_pp_ZVAL; // ME190109
  str_fin.species_pp_vLDAU = xvasp.str.species_pp_vLDAU; // ME190109

//[ME190109] Can't just set the updated positions because the symmetry may have changed
//[ME190109]  xvasp.str.lattice = str_fin.lattice;
//[ME190109]  for (uint i = 0; i < xvasp.str.atoms.size(); i++) {
//[ME190109]    xvasp.str.atoms[i].cpos = str_fin.atoms[i].cpos;
//[ME190109]    xvasp.str.atoms[i].fpos = str_fin.atoms[i].fpos;
//[ME190109]  }

  xvasp.str = str_fin; // ME190109
  if(LDEBUG){std::cout << xvasp.str << std::endl;} // ME190308
  pflow::fixEmptyAtomNames(xvasp.str,true);  // ME190308

  // Safeguard for when CONVERT is set in the aflow.in file
  if (!vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.flag("PRESERVE") &&
      !vflags.KBIN_VASP_FORCE_OPTION_CONVERT_UNIT_CELL.xscheme.empty()) { // ME190109
  xvasp.str.Standard_Lattice_primitive = false;
  xvasp.str.Standard_Lattice_conventional = false;
  VASP_Convert_Unit_Cell(xvasp, vflags, aflags, fileMessage, aus);
  }

  // Restore original settings
  vflags.KBIN_VASP_FORCE_OPTION_PREC.pop(prec_phonons);
  vflags.KBIN_VASP_FORCE_OPTION_PREC.push(prec);
  vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.isentry = prec_entry;

  if (!kppra_phonons.empty()) {
    vflags.KBIN_VASP_KPOINTS_KPPRA.pop(kppra_phonons);
    vflags.KBIN_VASP_KPOINTS_KPPRA.push(kppra);
    vflags.KBIN_VASP_KPOINTS_KPPRA.isentry = kppra_entry;
  }

  if (!kscheme_phonons.empty()) {
    vflags.KBIN_VASP_KPOINTS_KSCHEME.pop(kscheme_phonons);
    vflags.KBIN_VASP_KPOINTS_KSCHEME.push(kscheme);
    vflags.KBIN_VASP_KPOINTS_KSCHEME.isentry = kscheme_entry;
  }

  if (!kpoints_phonons.empty()) {
    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.pop(kpoints_phonons);
    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.push(kpoints);
    vflags.KBIN_VASP_FORCE_OPTION_KPOINTS.isentry = kpoints_entry;
  }
  return Krun;
}

void VASP_RunPhonons_APL(_xvasp& xvasp,
                         string AflowIn,
                         _aflags& aflags,
                         _kflags& kflags,
                         _vflags& vflags, ofstream& messageFile) {
  _xinput xinput(xvasp);
  _xflags xflags(vflags);
  return RunPhonons_APL(xinput,AflowIn,aflags,kflags,xflags,messageFile);
}

void RunPhonons_APL(_xinput& xinput,
                    string AflowIn,
                    _aflags& aflags,
                    _kflags& kflags,
                    _xflags& xflags, 
                    ofstream& messageFile) {
  return RunPhonons_APL_181216(xinput,AflowIn,aflags,kflags,xflags,messageFile);
}
void RunPhonons_APL_181216(_xinput& xinput,
                    string AflowIn,
                    _aflags& aflags,
                    _kflags& kflags,
                    _xflags& xflags, 
                    ofstream& messageFile) {

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                               INITIALIZE                                //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  xinput.xvasp.AVASP_arun = true;
  // Test
    //if (!(kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_AAPL)) return; //PN180705
    if (!(kflags.KBIN_PHONONS_CALCULATION_APL ||
          kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_QHA_A || kflags.KBIN_PHONONS_CALCULATION_QHA_B || kflags.KBIN_PHONONS_CALCULATION_QHA_C ||
          kflags.KBIN_PHONONS_CALCULATION_SCQHA || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C ||
          kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C || //PN180717
          kflags.KBIN_PHONONS_CALCULATION_AAPL)) return; //PN180705

  //we make certain automatic fixes if we're within our domain, otherwise issue warning/error
  // ME 190109 - OBSOLETE because this doesn't make any sense. The automatic fix that is applied is
  // the VASP version, which can be done anywhere. Moreover, the hostname becomes the name of the node
  //  whenever a job is sent to the queue.
//  ::_WITHIN_DUKE_ = (aurostd::substring2bool(XHOST.hostname, "nietzsche") || aurostd::substring2bool(XHOST.hostname, "aflowlib") || aurostd::substring2bool(XHOST.hostname, "qrats") || aurostd::substring2bool(XHOST.hostname, "habana") || aurostd::substring2bool(XHOST.hostname, "quser"));

  //determine if we have a consistent mode defined between input and flags
  if(xinput.AFLOW_MODE_VASP){
    if(!xflags.AFLOW_MODE_VASP){
      cerr << "ERROR: KBIN::RunPhonons_APL: mismatch types between xinput(VASP) and xflags(!VASP)" << endl;
      return;
    }
  } else if(xinput.AFLOW_MODE_AIMS){
    if(!xflags.AFLOW_MODE_AIMS){
      cerr << "ERROR: KBIN::RunPhonons_APL: mismatch types between xinput(AIMS) and xflags(!AIMS)" << endl;
      return;
    }
  //} else if(xinput.AFLOW_MODE_ALIEN){  //alien doesn't have xstr, so we ignore
  } else {
    cerr << "ERROR: KBIN::RunPhonons_APL: unknown input type" << endl;
    return;
  }

  //corey
  //fix names if necessary, but do not remove pp info, we need mass!
  if(xinput.AFLOW_MODE_VASP){pflow::fixEmptyAtomNames(xinput.xvasp.str,true);}

  // Setup our logger
  apl::Logger logger(messageFile, aflags);
  //logger.setModuleName("PHONONS");  //will rename later
  
  string _ASTROPT_; //CO 170601
  if (kflags.KBIN_PHONONS_CALCULATION_AAPL) {
    logger.setModuleName("AAPL");  //CO 170601
    _ASTROPT_ = _ASTROPT_AAPL_;    //CO 170601
    } else if (kflags.KBIN_PHONONS_CALCULATION_QHA || 
               kflags.KBIN_PHONONS_CALCULATION_QHA_A || 
               kflags.KBIN_PHONONS_CALCULATION_QHA_B || 
               kflags.KBIN_PHONONS_CALCULATION_QHA_C || 
               kflags.KBIN_PHONONS_CALCULATION_SCQHA ||
               kflags.KBIN_PHONONS_CALCULATION_SCQHA_A ||
               kflags.KBIN_PHONONS_CALCULATION_SCQHA_B ||
               kflags.KBIN_PHONONS_CALCULATION_SCQHA_C ||
               kflags.KBIN_PHONONS_CALCULATION_QHA3P ||
               kflags.KBIN_PHONONS_CALCULATION_QHA3P_A ||
               kflags.KBIN_PHONONS_CALCULATION_QHA3P_B ||
               kflags.KBIN_PHONONS_CALCULATION_QHA3P_C) { //PN180705
    logger.setModuleName("QHA");  //CO 170601
    _ASTROPT_ = _ASTROPT_QHA_;    //CO 170601
  } else {
    logger.setModuleName("APL");  //CO 170601
    _ASTROPT_ = _ASTROPT_APL_;    //CO 170601
  }

  logger << "RUNNING..." << apl::endl;

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                      READ PARAMETERS FROM AFLOWIN                       //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  // ME 181019 - Overwrite defaults for CHGCAR, WAVECAR, etc. Only write
  // these files if the user explicitly sets these flags. Otherwise, APL will
  // use too much disk space.
  if (xinput.AFLOW_MODE_VASP) {
    xflags.vflags.KBIN_VASP_FORCE_OPTION_BADER.options2entry(AflowIn,_STROPT_+"BADER=",false);
    xinput.xvasp.aopts.flag("FLAG::AVASP_BADER",xflags.vflags.KBIN_VASP_FORCE_OPTION_BADER.option);
    if (!xflags.vflags.KBIN_VASP_FORCE_OPTION_BADER.isentry && DEFAULT_VASP_FORCE_OPTION_BADER) {
      logger << "Switching OFF BADER for APL calculations (default: OFF)." << apl::endl; //CO181226 - default OFF
      }

    xflags.vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.options2entry(AflowIn,_STROPT_+"CHGCAR=",false);
    xinput.xvasp.aopts.flag("FLAG::AVASP_CHGCAR",xflags.vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.option);
    if (!xflags.vflags.KBIN_VASP_FORCE_OPTION_CHGCAR.isentry && DEFAULT_VASP_FORCE_OPTION_CHGCAR) {
      logger << "Switching OFF CHGCAR for APL calculations (default: OFF)." << apl::endl;  //CO181226 - default OFF
      }

    xflags.vflags.KBIN_VASP_FORCE_OPTION_ELF.options2entry(AflowIn,_STROPT_+"ELF=",false);
    xinput.xvasp.aopts.flag("FLAG::AVASP_ELF",xflags.vflags.KBIN_VASP_FORCE_OPTION_ELF.option);
    if (!xflags.vflags.KBIN_VASP_FORCE_OPTION_ELF.isentry && DEFAULT_VASP_FORCE_OPTION_ELF) {
      logger << "Switching OFF ELF for APL calculations (default: OFF)." << apl::endl; //CO181226 - default OFF
      }

    xflags.vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.options2entry(AflowIn,_STROPT_+"WAVECAR=",false);
    xinput.xvasp.aopts.flag("FLAG::AVASP_WAVECAR",xflags.vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.option);
    if (!xflags.vflags.KBIN_VASP_FORCE_OPTION_WAVECAR.isentry && DEFAULT_VASP_FORCE_OPTION_WAVECAR) {
      logger << "Switching OFF WAVECAR for APL calculations (default: OFF)." << apl::endl;  //CO181226 - default OFF
        }
    }

  // APL ----------------------------------------------------------------------

  // ME 181026 START

  /***************************** READ PARAMETERS *****************************/

  string USER_ENGINE="", USER_FREQFORMAT="", USER_SUPERCELL="", DOS_MESH_SCHEME="", USER_DOS_METHOD="", USER_TPT="", USER_DC_METHOD=""; //CO190114 - initialize everything
  string USER_DC_INITLATTICE="", USER_DC_INITCOORDS_FRAC="", USER_DC_INITCOORDS_CART="", USER_DC_INITCOORDS_LABELS="", USER_DC_USERPATH=""; //CO190114 - initialize everything
  bool USER_DPM=false, USER_AUTO_DISTORTIONS=false, USER_DISTORTIONS_XYZ_ONLY=false, USER_DISTORTIONS_SYMMETRIZE=false, USER_DISTORTIONS_INEQUIVONLY=false, USER_RELAX=false, USER_ZEROSTATE=false; //CO190114 - initialize everything
  bool USER_HIBERNATE=false, USER_POLAR=false, USER_DC=false, USER_DOS=false, USER_TP=false;  //CO190114 - initialize everything
  double USER_DISTORTION_MAGNITUDE=false, USER_DOS_SMEAR=false, USER_TP_TSTART=false, USER_TP_TEND=false, USER_TP_TSTEP=false;  //CO190114 - initialize everything  
  int USER_MAXSHELL = 0, USER_MINSHELL = 0, USER_MINATOMS = 0, USER_MINATOMS_RESTRICTED = 0, USER_DC_NPOINTS = 0, USER_DOS_NPOINTS = 0, START_RELAX = 0;  //CO190114 - initialize everything
  vector<int> USER_DOS_MESH(3);
  for (uint i = 0; i < kflags.KBIN_MODULE_OPTIONS.aplflags.size(); i++) {
    const string& key = kflags.KBIN_MODULE_OPTIONS.aplflags[i].keyword;
    logger << (kflags.KBIN_MODULE_OPTIONS.aplflags[i].isentry? "Setting" : "DEFAULT") << " " << _ASTROPT_ << key << "=" << kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme << apl::endl;
    if (key == "RELAX") {USER_RELAX = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option; continue;}
    if (key == "HIBERNATE") {USER_HIBERNATE = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option; continue;}
    if (key == "ENGINE") {USER_ENGINE = kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme; continue;}
    if (key == "SUPERCELL") {USER_SUPERCELL = kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme; continue;}
    if (key == "MINATOMS") {USER_MINATOMS = kflags.KBIN_MODULE_OPTIONS.aplflags[i].content_int; continue;}
    if (key == "MINATOMS_RESTRICTED") {USER_MINATOMS_RESTRICTED = kflags.KBIN_MODULE_OPTIONS.aplflags[i].content_int; continue;}
    if (key == "MINSHELL") {USER_MINSHELL = kflags.KBIN_MODULE_OPTIONS.aplflags[i].content_int; continue;}
    if (key == "POLAR") {USER_POLAR = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option; continue;}
    if (key == "DMAG") {USER_DISTORTION_MAGNITUDE = kflags.KBIN_MODULE_OPTIONS.aplflags[i].content_double; continue;}
    if (key == "DXYZONLY") {USER_DISTORTIONS_XYZ_ONLY = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option; continue;}
    if (key == "DSYMMETRIZE") {USER_DISTORTIONS_SYMMETRIZE = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option; continue;}
    if (key == "DINEQUIV_ONLY") {USER_DISTORTIONS_INEQUIVONLY = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option; continue;} //CO190131
    // Special case: DPM can be true, false, or empty
    if (key == "DPM") {
      USER_DPM = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option;
      USER_AUTO_DISTORTIONS=(!kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme.empty() && (kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme[0] == 'A' || kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme[0] == 'a')); //CO190131
      continue;
      }
    if (key == "ZEROSTATE") {USER_ZEROSTATE = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option; continue;}
    if (key == "FREQFORMAT") {USER_FREQFORMAT = aurostd::toupper(kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme); continue;}
    if (key == "DC") {USER_DC = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option; continue;}
    if (key == "DCPATH") {USER_DC_METHOD = kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme; continue;}
    if (key == "DCPOINTS") {USER_DC_NPOINTS = kflags.KBIN_MODULE_OPTIONS.aplflags[i].content_int; continue;}
    if (key == "DCINITCOORDSFRAC") {USER_DC_INITCOORDS_FRAC = kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme; continue;}
    if (key == "DCINITCOORDSCART") {USER_DC_INITCOORDS_CART = kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme; continue;}
    if (key == "DCINITCOORDSLABELS") {USER_DC_INITCOORDS_LABELS = kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme; continue;}
    if (key == "DCUSERPATH") {USER_DC_USERPATH = kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme; continue;}
    if (key == "DOS") {USER_DOS = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option; continue;}
    if (key == "DOSMESH") {DOS_MESH_SCHEME = kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme; continue;}
    if (key == "DOSMETHOD") {USER_DOS_METHOD = kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme; continue;}
    if (key == "DOSSMEAR") {USER_DOS_SMEAR = kflags.KBIN_MODULE_OPTIONS.aplflags[i].content_double; continue;}
    if (key == "DOSPOINTS") {USER_DOS_NPOINTS = kflags.KBIN_MODULE_OPTIONS.aplflags[i].content_int; continue;}
    if (key == "TP") {USER_TP = kflags.KBIN_MODULE_OPTIONS.aplflags[i].option; continue;}
    if (key == "TPT") {USER_TPT = kflags.KBIN_MODULE_OPTIONS.aplflags[i].xscheme; continue;}
      }

  /***************************** CHECK PARAMETERS *****************************/

  try {
    string message;

    // ME 190313 - START
    // Do not relax with --generate_aflowin_only option
    if (XHOST.GENERATE_AFLOWIN_ONLY && USER_RELAX) {
      USER_RELAX = false;
      logger << apl::warning << "RELAX will be switched OFF for generate_aflowin_only." << apl::endl;
    }
    // ME 190313 - END
    // Relax structure
    if (USER_RELAX) {
      START_RELAX = 1;
      if(xinput.AFLOW_MODE_VASP){
        // Check if the structure has already been relaxed
        for (int i = 1; i <= _NUM_RELAX_; i++) {
          string contcar = "CONTCAR." + _APL_RELAX_PREFIX_ + aurostd::utype2string<int>(i);
          if (aurostd::EFileExist(contcar)) {
            START_RELAX++;
          } else {
            break;
      }
        }
        if (START_RELAX == _NUM_RELAX_ + 1) {
          logger << "Structure has already been relaxed. Relaxation will be skipped." << apl::endl;
          USER_RELAX = false;
        } else if (START_RELAX != 1) {
          logger << "APL has already performed " << (START_RELAX + 1) << " relaxations.";
          logger << "Number of relaxations remaining: " << (_NUM_RELAX_ - START_RELAX + 1) << "." << apl::endl;
      }
      } else {
        USER_RELAX = false;
        logger << apl::warning << "RELAX option only supported for VASP. Relaxations will be skipped." << apl::endl;
        }
          }

    // Correct user engine
    if (USER_ENGINE == "GSA") {
      logger << "The Generalized Supercell Approach (GSA) is deprecated - replaced with the Direct Method (DM)." << apl::endl;
      USER_ENGINE = "DM";
            }
    if ((USER_ENGINE != "DM") && (USER_ENGINE != "LR")) {
      message = "Wrong setting in " + _ASTROPT_ + "ENGINE. Use either DM or LR. ";
      message += "See README_AFLOW_APL.TXT for more information.";
      throw apl::APLRuntimeError(message);
    }

    // Correct MINATOMS if restricted
    if (kflags.KBIN_MODULE_OPTIONS.minatoms_restricted) USER_MINATOMS = USER_MINATOMS_RESTRICTED;

    vector<string> tokens;
    // Check supercell settings
    if (!USER_SUPERCELL.empty()) {
        tokens.clear();
      apl::tokenize(USER_SUPERCELL, tokens, string(" xX"));
      if (tokens.size() != 3) {
        message = "Wrong setting in " + _ASTROPT_ + "SUPERCELL. ";
        message += "See README_AFLOW_APL.TXT for the correct format.";
        throw apl::APLRuntimeError(message);
      }
    }

    if (USER_DC) {
      if (USER_DC_METHOD == "LATTICE") {
        USER_DC_INITLATTICE = xinput.getXStr().bravais_lattice_type;
      } else if (USER_DC_METHOD == "MANUAL") {
        // Make sure that the number of coordinates and labels agree
      tokens.clear();
        if (!USER_DC_INITCOORDS_FRAC.empty()) {
          apl::tokenize(USER_DC_INITCOORDS_FRAC, tokens, string(" ;"));
        } else {
          apl::tokenize(USER_DC_INITCOORDS_CART, tokens, string(" ;"));
    }
        uint ncoords = tokens.size();
        tokens.clear();
        apl::tokenize(USER_DC_INITCOORDS_LABELS, tokens, string(" ,;"));  // ME190427 - also break along semicolon
        if (tokens.size() != ncoords) {
          message = "Mismatch between the number of points and the number of labels for the phonon dispersions. ";
          message += "Check the parameters DCINITCOORDS" + string(USER_DC_INITCOORDS_FRAC.empty()?"CART":"FRAC") + " and DCINITCOORDSLABELS.";
          message += "See README_AFLOW_APL.TXT for more information.";
          throw apl::APLRuntimeError(message);
      }
      } else {
        message = "Wrong setting in " + _ASTROPT_ + "DCPATH. Use either LATTICE or MANUAL. ";
        message += "See README_AFLOW_APL.TXT for more information.";
        throw apl::APLRuntimeError(message);
    }
    }

    // DOS
    if (USER_DOS || USER_TP) {  // ME190423
      tokens.clear();
      if (USER_DOS_METHOD != "LT" && USER_DOS_METHOD != "RS") {
        message = "Wrong setting in " + _ASTROPT_ + "DOSMETHOD. Use either LT or RS. ";
        message += "See README_AFLOW_APL.TXT for more information.";
        throw apl::APLRuntimeError(message);
      }
      if ((USER_DOS_METHOD == "RS") && (USER_DOS_SMEAR < _ZERO_TOL_)) {
        logger << apl::warning << "Smearing value for DOS not set or set to zero. ";
        logger << "APL will overwrite the smearing value to 0.05 eV." << apl::endl;
        USER_DOS_SMEAR = 0.05;
      }
      apl::tokenize(DOS_MESH_SCHEME, tokens, string(" xX"));
      if (tokens.size() != 3) {
        message = "Wrong setting in " + _ASTROPT_ + "DOSMESH. ";
        message += "See README_AFLOW_APL.TXT for the correct format.";
        throw apl::APLRuntimeError(message);
      } else {
        USER_DOS_MESH[0] = aurostd::string2utype<int>(tokens[0]);
        USER_DOS_MESH[1] = aurostd::string2utype<int>(tokens[1]);
        USER_DOS_MESH[2] = aurostd::string2utype<int>(tokens[2]);
      }
    }

    // TPT
    if (USER_TP) {
      tokens.clear();
      apl::tokenize(USER_TPT, tokens, string(" :"));
      if (tokens.size() != 3) {
        message = "Wrong setting in " + _ASTROPT_ + "TPT. ";
        message += "See README_AFLOW_APL.TXT for the correct format.";
        throw apl::APLRuntimeError(message);
      }
      USER_TP_TSTART = aurostd::string2utype<double>(tokens[0]);
      USER_TP_TEND = aurostd::string2utype<double>(tokens[1]);
      USER_TP_TSTEP = aurostd::string2utype<double>(tokens[2]);
    }
  } catch (std::exception& e) {
    logger << apl::error << e.what() << apl::endl;
    return;
  }

  /****************************** OUTPUT SUMMARY ******************************/

  logger << "Parameters for the Automatic Phonon Library successfully read." << apl::endl;
  logger << "The structure will " << (USER_RELAX?"":"NOT ") << "be relaxed before running APL." << apl::endl;
  logger << "The hibernate feature is switched " << (USER_HIBERNATE?"ON":"OFF") << "." << apl::endl;
  logger << "Phonons will be calculated using the " << (USER_ENGINE=="DM"?"Direct":"Linear Response") << " Method." << apl::endl;

  if (USER_ENGINE == "DM") {
    logger << "The distortion magnitude will be " << USER_DISTORTION_MAGNITUDE << " Angstrom." << apl::endl;
    if (USER_DISTORTIONS_XYZ_ONLY) {
      logger << "Only distortions along the lattice vectors will be used." << apl::endl;
    } else {
      logger << "Atoms will be distorted along the lattice vectors, face diagonals, and body diagonals." << apl::endl;
    }
    if (USER_DISTORTIONS_SYMMETRIZE) {
      logger << "Non-symmetric distortion directions will be determined for each site." << apl::endl;
    } else {
      logger << "Distortions will be generated in three independent directions." << apl::endl;
    }
    if (USER_DISTORTIONS_INEQUIVONLY) { //CO190131
      logger << "Distortion directions will be determined for inequivalent site only." << apl::endl;
    } else {
      logger << "Distortion directions will be determined for ALL site." << apl::endl;
    }
    if (USER_AUTO_DISTORTIONS) {
      logger << "Positive/negative distortion directions will be determined for each site." << apl::endl; 
    } else if (USER_DPM) {
      logger << "Distortions will be generated in both the positive and negative direction." << apl::endl;
    } else {
      logger << apl::warning << "Distortions will only be generated in the positive direction - this is NOT recommended." << apl::endl;
    }
    logger << "Forces from the undistored state will " << (USER_ZEROSTATE?"":"NOT ") << "be used." << apl::endl;
  }

  logger << "Polar corrections will " << (USER_POLAR?"":"NOT ") << "be employed." << apl::endl;
  logger << "Frequencies will be returned in this format: " << USER_FREQFORMAT << "." << apl::endl;

  logger << "The supercell will be built using ";
  if (kflags.KBIN_MODULE_OPTIONS.supercell_method[0]) {
    logger << "the dimensions " << USER_SUPERCELL << "." << apl::endl;
  } else if (kflags.KBIN_MODULE_OPTIONS.supercell_method[1]) {
    logger << "at least " << USER_MINATOMS << " atoms." << apl::endl;
  } else if (kflags.KBIN_MODULE_OPTIONS.supercell_method[2]) {
    logger << "at most " << USER_MAXSHELL << " shells." << apl::endl;
  } else {
    logger << "at least " << USER_MINSHELL << " shells." << apl::endl;
  }

  if (USER_DC) {
    logger << "Phonon dispersion curves will be calculated ";
    if (USER_DC_METHOD == "LATTICE") {
        logger << "using the default path of the lattice (see DOI 10.1016/j.commatsci.2010.05.010). ";
    } else {
      logger << "along the " << (USER_DC_INITCOORDS_FRAC.empty()?"Cartesian":"fractional") << " coordinates ";
      logger << "[" << (USER_DC_INITCOORDS_FRAC.empty()?USER_DC_INITCOORDS_CART:USER_DC_INITCOORDS_FRAC) << "].";
    }
    logger << " Each subpath will be divided into " << USER_DC_NPOINTS << " points." << apl::endl;
  } else {
    logger << "Phonon dispersion curves will NOT be calculated." << apl::endl;
  }

  if (USER_DOS || USER_TP) {  // ME190423
    logger << "Phonon DOS will be calculated using the ";
    logger << (USER_DOS_METHOD == "LT"?"Linear Tetrahedron":"Root Sampling") << " method ";
    logger << "along a " << USER_DOS_MESH[0] << "x" << USER_DOS_MESH[1] << "x" << USER_DOS_MESH[2];
    logger << " mesh with " << USER_DOS_NPOINTS << " bins.";
    if (USER_DOS_METHOD == "RS")
      logger << " A smearing value of " << USER_DOS_SMEAR << " eV will be used.";
    logger << apl::endl;
  } else {
    logger << "Phonon DOS will NOT be calculated." << apl::endl;
    }

  if (USER_TP) {
    logger << "Thermodynamic properties will be calculated between ";
    logger << USER_TP_TSTART << " K and " << USER_TP_TEND << " K ";
    logger << "in " << USER_TP_TSTEP << " K steps." << apl::endl;
  } else {
    logger << "Thermodynamic properties will NOT be calculated." << apl::endl;
    }
  // ME 181026 END

  // AAPL ----------------------------------------------------------------------

  // ME 181027 START

  /***************************** READ PARAMETERS *****************************/

  string USER_BTE;
  bool USER_TCOND, USER_ISOTOPE, USER_BOUNDARY, USER_CUMULATIVEK, USER_AAPL_FOURTH_ORDER;
  double USER_NANO_SIZE, USER_EPS_SUM, USER_AAPL_MIX, USER_TCT_TSTART, USER_TCT_TEND, USER_TCT_TSTEP;
  int USER_AAPL_MAX_ITER, USER_KPPRA_AAPL = -1;  // ME190408 - Added KPPRA_AAPL
  vector<double> USER_CUTOFF_DISTANCE(2);
  vector<int> USER_CUTOFF_SHELL(2), USER_THERMALGRID(3);
  if (kflags.KBIN_PHONONS_CALCULATION_AAPL) {
    USER_TCOND = true;
    string CUTOFF_SCHEME, SHELL_SCHEME, THERMALGRID_SCHEME, TCT_SCHEME;
    for (uint i = 0; i < kflags.KBIN_MODULE_OPTIONS.aaplflags.size(); i++) {
      const string& key = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].keyword;
      logger << (kflags.KBIN_MODULE_OPTIONS.aaplflags[i].isentry? "Setting" : "DEFAULT") << " " << _ASTROPT_ << key << "=" << kflags.KBIN_MODULE_OPTIONS.aaplflags[i].xscheme << apl::endl;
      if (key == "BTE") {USER_BTE = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].xscheme; continue;}
      if (key == "FOURTH_ORDER") {USER_AAPL_FOURTH_ORDER = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].option; continue;}
      if (key == "CUT_RAD") {CUTOFF_SCHEME = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].xscheme; continue;}
      if (key == "CUT_SHELL") {SHELL_SCHEME = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].xscheme; continue;}
      if (key == "THERMALGRID") {THERMALGRID_SCHEME = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].xscheme; continue;}
      if (key == "TCT") {TCT_SCHEME = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].xscheme; continue;}
      if (key == "SUMRULE") {USER_EPS_SUM = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].content_double; continue;}
      if (key == "SUMRULE_MAX_ITER") {USER_AAPL_MAX_ITER = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].content_int; continue;}
      if (key == "MIXING_COEFFICIENT") {USER_AAPL_MIX = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].content_double; continue;}
      if (key == "ISOTOPE") {USER_ISOTOPE = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].option; continue;}
      if (key == "BOUNDARY") {USER_BOUNDARY = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].option; continue;}
      if (key == "CUMULATIVEK") {USER_CUMULATIVEK = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].option; continue;}
      if (key == "NANO_SIZE") {USER_NANO_SIZE = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].content_double; continue;}
      if (key == "KPPRA_AAPL") {USER_KPPRA_AAPL = kflags.KBIN_MODULE_OPTIONS.aaplflags[i].content_int; continue;}  // ME190408
    }

  /***************************** CHECK PARAMETERS *****************************/

      try {
      string message;
      // Correct BTE
      if (USER_BTE != "RTA" && USER_BTE != "FULL") {
        message = "Wrong setting in " + _ASTROPT_ + "BTE. Use either RTA or FULL.";
        message += "See README_AFLOW_APL.TXT for more information.";
        throw apl::APLRuntimeError(message);
      }

      vector<string> tokens;
      // CUT_SHELL and CUT_RAD
      bool defaults = (!kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[0] && !kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[1]);
      if (defaults || kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[0]) {
        tokens.clear();
        apl::tokenize(CUTOFF_SCHEME, tokens, string(" ,"));
        if (tokens.size() < 1) {
          message = "Not enought entries in " + _ASTROPT_ + "CUT_RAD. ";
          message += "See README_AFLOW_APL.TXT for more information.";
          throw apl::APLRuntimeError(message);
        } else if (tokens.size() > 2) {
          logger << apl::warning << "Too many entries for " << _ASTROPT_ << "CUT_RAD. ";
          logger << "Excess entries will be ignored." << apl::endl;
        }
        USER_CUTOFF_DISTANCE[0] = aurostd::string2utype<double>(tokens[0]);
        if (USER_AAPL_FOURTH_ORDER) {
          if (tokens.size() == 1) {
            logger << apl::warning << "Only one entry found for the cutoff radius. ";
            logger << "3rd and 4th order anharmonic IFCs will use the same value." << apl::endl;
            USER_CUTOFF_DISTANCE[1] = USER_CUTOFF_DISTANCE[0];
            } else {
            USER_CUTOFF_DISTANCE[1] = aurostd::string2utype<double>(tokens[1]);
            }
        }
      }

      if (defaults || kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[1]) {
        tokens.clear();
        apl::tokenize(SHELL_SCHEME, tokens, string(" ,"));
        if (tokens.size() < 1) {
          message = "Not enought entries in " + _ASTROPT_ + "CUT_SHELL. ";
          message += "See README_AFLOW_APL.TXT for more information.";
          throw apl::APLRuntimeError(message);
        } else if (tokens.size() > 2) {
          logger << apl::warning << "Too many entries for " << _ASTROPT_ << "CUT_SHELL. ";
          logger << "Excess entries will be ignored." << apl::endl;
        }
        USER_CUTOFF_SHELL[0] = aurostd::string2utype<int>(tokens[0]);
        if (USER_AAPL_FOURTH_ORDER) {
          if (tokens.size() == 1) {
            logger << apl::warning << "Only one entry found for the number of coordination shells. ";
            logger << "3rd and 4th order anharmonic IFCs will use the same value." << apl::endl;
            USER_CUTOFF_SHELL[1] = USER_CUTOFF_SHELL[0];
        } else {
            USER_CUTOFF_SHELL[1] = aurostd::string2utype<int>(tokens[1]);
        }
      }
    }
      // If only one parameter is set in the aflow.in file, unset the other  
      if (kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[0] != kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[1]) {
        if (!kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[0]) USER_CUTOFF_DISTANCE.assign(2, 0.0);
        if (!kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[1]) USER_CUTOFF_SHELL.assign(2, 0);
    }

      // THERMALGRID
      tokens.clear();
      apl::tokenize(THERMALGRID_SCHEME, tokens, string(" xX"));
      if (tokens.size() == 3) {
        USER_THERMALGRID[0] = aurostd::string2utype<int>(tokens[0]);
        USER_THERMALGRID[1] = aurostd::string2utype<int>(tokens[1]);
        USER_THERMALGRID[2] = aurostd::string2utype<int>(tokens[2]);
      } else {
        message = "Wrong setting in " + _ASTROPT_ + "THERMALGRID. ";
        message += "See README_AFLOW_APL.TXT for the correct format.";
        throw apl::APLRuntimeError(message);
    }

      // TCT
      tokens.clear();
      apl::tokenize(TCT_SCHEME, tokens, string(" :"));
      if (tokens.size() == 3) {
        USER_TCT_TSTART = aurostd::string2utype<double>(tokens[0]);
        USER_TCT_TEND = aurostd::string2utype<double>(tokens[1]);
        USER_TCT_TSTEP = aurostd::string2utype<double>(tokens[2]);
        if (USER_TCT_TSTART == 0) {
          logger << apl::warning << "Thermal conductivity is infinite at 0 K and will be skipped." << apl::endl;
          USER_TCT_TSTART += USER_TCT_TSTEP;
        }
      } else {
        message = "Wrong setting in " + _ASTROPT_ + "TCT. ";
        message += "See README_AFLOW_APL.TXT for the correct format.";
        throw apl::APLRuntimeError(message);
    }

      // BOUNDARY and CUMULATIVEK
      if (USER_BOUNDARY && USER_CUMULATIVEK) {
        USER_CUMULATIVEK = false;
        logger << apl::warning << "Both boundary effects and cumulative thermal conductivity cannot be ";
        logger << "set at the same time. Cumulative thermal conductivity has been switched off." << apl::endl;
      }
    } catch (std::exception& e) {
      logger << apl::error << e.what() << apl::endl;
      return;
    }

  /****************************** OUTPUT SUMMARY ******************************/

    logger << "Parameters for the Automatic Anharmonic Phonon Library successfully read." << apl::endl;
    logger << "Four-phonon processes will " << (USER_AAPL_FOURTH_ORDER?"":"NOT ") << "be included in the calculations." << apl::endl;

    bool defaults = (!kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[0] && !kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[1]);
    if (defaults || kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[0]) {
      logger << "The cutoff to compute the 3rd order anharmonic IFCs will be ";
      logger << USER_CUTOFF_DISTANCE[0] << " Angstrom." << apl::endl;
      if (USER_AAPL_FOURTH_ORDER) {
        logger << "The cutoff to compute the 4th order anharmonic IFCs will be ";
        logger << USER_CUTOFF_DISTANCE[1] << " Angstrom." << apl::endl;
      }
    }
    if (defaults || kflags.KBIN_MODULE_OPTIONS.cut_rad_shell[1]) {
      logger << "The calculation of 3rd order anharmonic IFCs will consider up to ";
      logger << USER_CUTOFF_SHELL[0] << " coordination shells." << apl::endl;
      if (USER_AAPL_FOURTH_ORDER) {
        logger << "The calculation of 4th order anharmonic IFCs will consider up to ";
        logger << USER_CUTOFF_SHELL[1] << " coordination shells." << apl::endl;
      }
    }

    if (USER_KPPRA_AAPL > 0) logger << "AAPL will use a KPPRA of " << USER_KPPRA_AAPL << " for static calculations." << apl::endl;  // ME190408

    logger << "Anharmonic IFCs will be calculated with a convergence criterion of " << USER_EPS_SUM << "." << apl::endl;
    logger << "A mixing coefficient of " << USER_AAPL_MIX << " will be used." << apl::endl;
    logger << "Anharmonic IFCs need to be converged within " << USER_AAPL_MAX_ITER << " iterations." << apl::endl;

    logger << "Thermal conductivity will be calculated between ";
    logger << USER_TCT_TSTART << " K and " << USER_TCT_TEND << "K ";
    logger << "in " << USER_TCT_TSTEP << " K steps." << apl::endl;

    logger << "The Boltzmann Transport Equation will be solved using ";
    if (USER_BTE == "RTA") {
      logger << "the Relaxation Time Approximation Approximation (RTA)." << apl::endl;
      } else {
      logger << "an iterative scheme." << apl::endl;
      }
    logger << "The equation will be solved using the adaptive broadening method along a ";
    logger << USER_THERMALGRID[0] << "x" << USER_THERMALGRID[1] << "x" << USER_THERMALGRID[2] << " q-point mesh." << apl::endl;
    logger << "Isotope effects will " << (USER_ISOTOPE?"":"NOT ") << "be included." << apl::endl;
    if (USER_BOUNDARY || USER_CUMULATIVEK) {
      logger << "Boundary effects will be included via ";
      if (USER_BOUNDARY) {
        logger << "boundary scattering with a grain size of ";
      } else {
        logger << "cumulative thermal conductivity and a mean free path of at most ";
      }
      logger << USER_NANO_SIZE << " nm." << apl::endl;
    } else {
      logger << "Boundary effects will NOT be included." << apl::endl;
    }
  } else {
    USER_TCOND = false;
    logger << "Anharmonic force constants and thermal conductivity will NOT be calculated." << apl::endl;
    }
  // ME 181027 START

  vector<string> tokens;

  // QHA ----------------------------------------------------------------------

    //  //PINKU QUASI-HARMONIC START
    aurostd::xoption CALCULATE_GROUPVELOCITY_OPTION; CALCULATE_GROUPVELOCITY_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_ATOMIC_DISPLACEMENT_OPTION; CALCULATE_ATOMIC_DISPLACEMENT_OPTION.option = false; //PN180705
  aurostd::xoption CALCULATE_GRUNEISEN_OPTION; CALCULATE_GRUNEISEN_OPTION.option = false;
  aurostd::xoption CALCULATE_DISPLACEMENTS_OPTION; CALCULATE_DISPLACEMENTS_OPTION.option = false;
  aurostd::xoption CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION; CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option = false;
  aurostd::xoption CALCULATE_EOS_OPTION; CALCULATE_EOS_OPTION.option = false;
  aurostd::xoption CALCULATE_EOS_SUBDIRECTORIES_OPTION; CALCULATE_EOS_SUBDIRECTORIES_OPTION.option = false;
    aurostd::xoption EDOS_ACURATE_OPTION; EDOS_ACURATE_OPTION.option = false; //PN180705
    aurostd::xoption INCLUDE_ELE_OPTION;  INCLUDE_ELE_OPTION.option = false; //PN180705
    //Anisotropic Gruneisen and EOS //PN180705
    //in the a-direction //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_A_OPTION; CALCULATE_GRUNEISEN_A_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION; CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option = false; //PN180705
    aurostd::xoption GP_DISTORTION_OPTION; GP_DISTORTION_OPTION.xscheme = "0.03"; double GP_DISTORTION = 0.03; //PN180705
      
    //in the b-direction //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_B_OPTION; CALCULATE_GRUNEISEN_B_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION; CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option = false; //PN180705

    //in the c-direction //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_C_OPTION; CALCULATE_GRUNEISEN_C_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION; CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option = false; //PN180705

    //SC-QHA //PN180705
    aurostd::xoption CALCULATE_SCQHA_OPTION; CALCULATE_SCQHA_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_SCQHA_SUBDIRECTORIES_OPTION; CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option = false; //PN180705
    //QHA3P option
    aurostd::xoption CALCULATE_QHA3P_OPTION; CALCULATE_QHA3P_OPTION.option = false; //PN180705

    //Anisotropic SCQHA EOS //PN180705
    //in the a-direction //PN180705
    aurostd::xoption CALCULATE_SCQHA_A_OPTION; CALCULATE_SCQHA_A_OPTION.option = false; //PN180705
    aurostd::xoption SCQHA_DISTORTION_OPTION; SCQHA_DISTORTION_OPTION.xscheme = "3.0"; double SCQHA_DISTORTION = 3.0; //PN180705
    aurostd::xoption CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION; CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option = false; //PN180705

    //in the b-direction //PN180705
    aurostd::xoption CALCULATE_SCQHA_B_OPTION; CALCULATE_SCQHA_B_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION; CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option = false; //PN180705

    //in the c-direction //PN180705
    aurostd::xoption CALCULATE_SCQHA_C_OPTION; CALCULATE_SCQHA_C_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION; CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.option = false; //PN180705

    //QHA3P a direction
    aurostd::xoption CALCULATE_QHA3P_A_OPTION; CALCULATE_QHA3P_A_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_QHA3P_A_SUBDIRECTORIES_OPTION; CALCULATE_QHA3P_A_SUBDIRECTORIES_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_QHA3P_B_OPTION; CALCULATE_QHA3P_B_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_QHA3P_B_SUBDIRECTORIES_OPTION; CALCULATE_QHA3P_B_SUBDIRECTORIES_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_QHA3P_C_OPTION; CALCULATE_QHA3P_C_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_QHA3P_C_SUBDIRECTORIES_OPTION; CALCULATE_QHA3P_C_SUBDIRECTORIES_OPTION.option = false; //PN180705


    //QHA, QHA3P and SCQHA options initializing from previous options
    CALCULATE_GRUNEISEN_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA;
    CALCULATE_GRUNEISEN_A_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA_A;
    CALCULATE_GRUNEISEN_B_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA_B;
    CALCULATE_GRUNEISEN_C_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA_C;

    CALCULATE_SCQHA_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_SCQHA;
    CALCULATE_SCQHA_A_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_SCQHA_A;
    CALCULATE_SCQHA_B_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_SCQHA_B;
    CALCULATE_SCQHA_C_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_SCQHA_C;

    CALCULATE_QHA3P_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA3P;
    CALCULATE_QHA3P_A_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA3P_A;
    CALCULATE_QHA3P_B_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA3P_B;
    CALCULATE_QHA3P_C_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA3P_C;

  // different types of fitting options for EOS calculations
  // (1) BM1 => Murnaghan EOS
  // (2) BM2 => Birch-Murnaghan 3rd-order EOS
  // (3) BM3 => Birch-Murnaghan 4th-order EOS
  //[OBSOLETE PN180705]aurostd::xoption GP_VOL_DISTORTION_OPTION; GP_VOL_DISTORTION_OPTION.xscheme = "0.03"; double GP_VOL_DISTORTION = 0.03;
  aurostd::xoption USER_PROJECTION_DIR_OPTION; USER_PROJECTION_DIR_OPTION.xscheme = "1:1:1"; vector<double> directions(3, 0); directions[0] = 1; directions[1] = 1; directions[2] = 1;   // 3 Miller indices
  aurostd::xoption CUTOFF_FREQ_OPTION; CUTOFF_FREQ_OPTION.xscheme="1e-5"; double CUTOFF_FREQ = 1e-5;  //in amu
    aurostd::xoption EOS_DISTORTION_RANGE_OPTION; EOS_DISTORTION_RANGE_OPTION.xscheme = "-3:6:1"; double EOS_DISTORTION_START = -3; double EOS_DISTORTION_END = 6; double EOS_DISTORTION_DISTORTION_INC=1; //PN180705
    aurostd::xoption EOS_STATIC_KPPRA_OPTION; EOS_STATIC_KPPRA_OPTION.xscheme = "10000"; int EOS_STATIC_KPPRA = 10000; //PN180705
    aurostd::xoption NEDOS_OPTION; NEDOS_OPTION.xscheme = "5000"; int NEDOS = 5000; //PN180705
  aurostd::xoption FITTING_TYPE_OPTION; FITTING_TYPE_OPTION.xscheme = "BM1"; string FITTING_TYPE = "BM1";
    aurostd::xoption SCQHA_PDIS_T_OPTION; SCQHA_PDIS_T_OPTION.xscheme = "100,400,600"; std::vector<double> scqha_pdis_T; //PN180705
  //PINKU QUASI-HARMONIC END
  
//PINKU PHONON START
    //GROUPVELOCITY=y/n
    CALCULATE_GROUPVELOCITY_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GROUP_VELOCITY=" + "|" + _ASTROPT_APL_OLD_ + "GROUP_VELOCITY="), CALCULATE_GROUPVELOCITY_OPTION.option,  CALCULATE_GROUPVELOCITY_OPTION.xscheme);
    logger << (CALCULATE_GROUPVELOCITY_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "GROUP_VELOCITY=" << (CALCULATE_GROUPVELOCITY_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

    //ATOMIC_DISPLACEMENT=y/n
    CALCULATE_ATOMIC_DISPLACEMENT_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "ATOMIC_DISPLACEMENT=" + "|" + _ASTROPT_APL_OLD_ + "ATOMIC_DISPLACEMENT="), CALCULATE_ATOMIC_DISPLACEMENT_OPTION.option,  CALCULATE_ATOMIC_DISPLACEMENT_OPTION.xscheme);
    logger << (CALCULATE_ATOMIC_DISPLACEMENT_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "ATOMIC_DISPLACEMENT=" << (CALCULATE_ATOMIC_DISPLACEMENT_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    //PINKU PHONON END

    //PINKU QUASI-HARMONIC START
    if(!USER_TCOND){
      if(kflags.KBIN_PHONONS_CALCULATION_QHA){
        CALCULATE_GRUNEISEN_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_QHA; //recycle what we parsed earlier
        logger << _ASTROPT_ << "CALC is" << ( CALCULATE_GRUNEISEN_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;
        logger << "The Gruneisen parameter will" << ( CALCULATE_GRUNEISEN_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-QHA." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA_A){
        CALCULATE_GRUNEISEN_A_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_QHA_A;
        logger << _ASTROPT_ << "CALC is" << ( CALCULATE_GRUNEISEN_A_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;
        logger << "The Gruneisen_A parameter will" << ( CALCULATE_GRUNEISEN_A_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-QHA." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA_B){
        CALCULATE_GRUNEISEN_B_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_QHA_B;
        logger << _ASTROPT_ << "CALC is" << ( CALCULATE_GRUNEISEN_B_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;
        logger << "The Gruneisen_B parameter will" << ( CALCULATE_GRUNEISEN_B_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-QHA." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA_C){
        CALCULATE_GRUNEISEN_C_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_QHA_C;
        logger << _ASTROPT_ << "CALC is" << ( CALCULATE_GRUNEISEN_C_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;
        logger << "The Gruneisen_C parameter will" << ( CALCULATE_GRUNEISEN_C_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-QHA." << apl::endl;
      }

      //QHA, QHA3P and SCQHA INCLUDE ELECTRONIC OPTION
      if(kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_QHA_A || kflags.KBIN_PHONONS_CALCULATION_QHA_B || kflags.KBIN_PHONONS_CALCULATION_QHA_C ||
          kflags.KBIN_PHONONS_CALCULATION_SCQHA|| kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C){
        INCLUDE_ELE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "INCLUDE_ELE=" + "|" + _ASTROPT_ + "INCLUDE_ELE="), INCLUDE_ELE_OPTION.option, INCLUDE_ELE_OPTION.xscheme);
        logger << (INCLUDE_ELE_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "INCLUDE_ELE=" << (INCLUDE_ELE_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      }

      //QHA3P and SCQHA temperature dependent phonon dispersion option
      if(kflags.KBIN_PHONONS_CALCULATION_SCQHA|| kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C){
        SCQHA_PDIS_T_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "SCQHA_PDIS_T=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_PDIS_T="), SCQHA_PDIS_T_OPTION.option, SCQHA_PDIS_T_OPTION.xscheme);
        tokens.clear(); scqha_pdis_T.clear();
        apl::tokenize(SCQHA_PDIS_T_OPTION.content_string, tokens, string(" ,"));
        if (tokens.size() == 0) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"SCQHA_PDIS_T. Specify as SCQHA_PDIS_T=-100.0, 300.0, 600.0");}
        if(tokens.size()!=0){
          for (uint i=0; i<tokens.size(); i++){
            scqha_pdis_T.push_back(aurostd::string2utype<double>(tokens.at(i)));
          }
        }
      }
      //rectifying possible user QHA-input errors
      if(CALCULATE_GRUNEISEN_OPTION.option){
        CALCULATE_SCQHA_A_OPTION.option=false;
        CALCULATE_SCQHA_B_OPTION.option=false;
        CALCULATE_SCQHA_C_OPTION.option=false;
      } else if(CALCULATE_GRUNEISEN_A_OPTION.option){
        CALCULATE_SCQHA_OPTION.option=false;
        CALCULATE_SCQHA_B_OPTION.option=false;
        CALCULATE_SCQHA_C_OPTION.option=false;
      } else if(CALCULATE_GRUNEISEN_B_OPTION.option){
        CALCULATE_SCQHA_OPTION.option=false;
        CALCULATE_SCQHA_A_OPTION.option=false;
        CALCULATE_SCQHA_C_OPTION.option=false;
      } else if(CALCULATE_GRUNEISEN_C_OPTION.option){
        CALCULATE_SCQHA_OPTION.option=false;
        CALCULATE_SCQHA_A_OPTION.option=false;
        CALCULATE_SCQHA_B_OPTION.option=false;
      }
      
      //Writing to log
      if(kflags.KBIN_PHONONS_CALCULATION_SCQHA){
        logger << (CALCULATE_SCQHA_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA=" << (CALCULATE_SCQHA_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_SCQHA_A){
        logger << (CALCULATE_SCQHA_A_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_A=" << (CALCULATE_SCQHA_A_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_SCQHA_B){
        logger << (CALCULATE_SCQHA_B_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_B=" << (CALCULATE_SCQHA_B_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_SCQHA_C){
        logger << (CALCULATE_SCQHA_C_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_C=" << (CALCULATE_SCQHA_C_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      }

      //Writing to log
      if(kflags.KBIN_PHONONS_CALCULATION_QHA3P){
        logger << (CALCULATE_QHA3P_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "QHA3P=" << (CALCULATE_QHA3P_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA3P_A){
        logger << (CALCULATE_QHA3P_A_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "QHA3P_A=" << (CALCULATE_QHA3P_A_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA3P_B){
        logger << (CALCULATE_QHA3P_B_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "QHA3P_B=" << (CALCULATE_QHA3P_B_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA3P_C){
        logger << (CALCULATE_QHA3P_C_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "QHA3P_C=" << (CALCULATE_QHA3P_C_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    }

      if(CALCULATE_SCQHA_OPTION.option || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C){
        SCQHA_DISTORTION_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "SCQHA_DISTORTION=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_DISTORTION="), SCQHA_DISTORTION_OPTION.option, SCQHA_DISTORTION_OPTION.xscheme);
        SCQHA_DISTORTION=SCQHA_DISTORTION_OPTION.content_double;
        if (SCQHA_DISTORTION_OPTION.isentry) {
          tokens.clear();
          apl::tokenize(SCQHA_DISTORTION_OPTION.content_string, tokens, string(" "));
          if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"SCQHA_DISTORTION. Specify as SCQHA_DISTORTION_OPTION=3.0.");}
        }
        logger << (SCQHA_DISTORTION_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "SCQHA_DISTORTION=" << SCQHA_DISTORTION << "." << apl::endl;
      }

      if(CALCULATE_GRUNEISEN_OPTION.option || kflags.KBIN_PHONONS_CALCULATION_QHA_A || kflags.KBIN_PHONONS_CALCULATION_QHA_B || kflags.KBIN_PHONONS_CALCULATION_QHA_C || kflags.KBIN_PHONONS_CALCULATION_SCQHA ||
          kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C || kflags.KBIN_PHONONS_CALCULATION_QHA3P ||
          kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C){

        CUTOFF_FREQ_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "CUTOFF_FREQ=" + "|" + _ASTROPT_APL_OLD_ + "CUTOFF_FREQ="), CUTOFF_FREQ_OPTION.option, CUTOFF_FREQ_OPTION.xscheme); //CO 170601
        CUTOFF_FREQ = CUTOFF_FREQ_OPTION.content_double;
        if (CUTOFF_FREQ_OPTION.isentry) {
      tokens.clear();
          apl::tokenize(CUTOFF_FREQ_OPTION.content_string, tokens, string(" "));
          if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"CUTOFF_FREQ. Specify as CUTOFF_FREQ=0.01.");}
        }
        logger << (CUTOFF_FREQ_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "CUTOFF_FREQ=" << CUTOFF_FREQ_OPTION.content_string << "." << apl::endl;


        CALCULATE_EOS_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "EOS=" + "|" + _ASTROPT_APL_OLD_ + "EOS="), CALCULATE_EOS_OPTION.option, CALCULATE_EOS_OPTION.xscheme); //CO 170601
        logger << (CALCULATE_EOS_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "EOS=" << (CALCULATE_EOS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

        if ((CALCULATE_EOS_OPTION.option) && (CALCULATE_GRUNEISEN_OPTION.option || kflags.KBIN_PHONONS_CALCULATION_QHA_A || kflags.KBIN_PHONONS_CALCULATION_QHA_B ||
              kflags.KBIN_PHONONS_CALCULATION_QHA_C || kflags.KBIN_PHONONS_CALCULATION_SCQHA || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A ||
              kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C || kflags.KBIN_PHONONS_CALCULATION_QHA3P ||
              kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C))

        {
          EOS_DISTORTION_RANGE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "EOS_DISTORTION_RANGE=" + "|" + _ASTROPT_APL_OLD_ + "EOS_DISTORTION_RANGE="), EOS_DISTORTION_RANGE_OPTION.option, EOS_DISTORTION_RANGE_OPTION.xscheme);
          tokens.clear();
          apl::tokenize(EOS_DISTORTION_RANGE_OPTION.content_string, tokens, string(" :"));
          if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"EOS_DISTORTION_RANGE. Specify as EOS_DISTORTION_RANGE=-3:6:1.");}
          EOS_DISTORTION_START = aurostd::string2utype<double>(tokens.at(0));
          EOS_DISTORTION_END = aurostd::string2utype<double>(tokens.at(1));
          EOS_DISTORTION_DISTORTION_INC = aurostd::string2utype<double>(tokens.at(2));
          logger << (EOS_DISTORTION_RANGE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "EOS_DISTORTION_RANGE=" << EOS_DISTORTION_RANGE_OPTION.content_string << "." << apl::endl;
          logger << "The EOS properties will be calculated in distortion range <" << EOS_DISTORTION_START << "," << EOS_DISTORTION_END << "," << EOS_DISTORTION_DISTORTION_INC << "." << apl::endl;

          FITTING_TYPE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "FITTING_TYPE=" + "|" + _ASTROPT_APL_OLD_ + "FITTING_TYPE="), FITTING_TYPE_OPTION.option, FITTING_TYPE_OPTION.xscheme); //CO 170601
          FITTING_TYPE = FITTING_TYPE_OPTION.content_string;
          if (FITTING_TYPE_OPTION.isentry) {
            tokens.clear();
            apl::tokenize(FITTING_TYPE_OPTION.content_string, tokens, string(" "));
            if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"FITTING_TYPE. Specify as FITTING_TYPE=BM2.");}
    }
          logger << (FITTING_TYPE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "FITTING_TYPE=" << FITTING_TYPE_OPTION.content_string << "." << apl::endl;
          logger << "EOS fitting type found = " << FITTING_TYPE << "." << apl::endl;

          EOS_STATIC_KPPRA_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "EOS_STATIC_KPPRA=" + "|" + _ASTROPT_APL_OLD_ + "EOS_STATIC_KPPRA="), EOS_STATIC_KPPRA_OPTION.option, EOS_STATIC_KPPRA_OPTION.xscheme);
          EOS_STATIC_KPPRA = EOS_STATIC_KPPRA_OPTION.content_int;
          if(CALCULATE_GRUNEISEN_OPTION.option){
            if (EOS_STATIC_KPPRA_OPTION.isentry) {
              tokens.clear();
              apl::tokenize(EOS_STATIC_KPPRA_OPTION.content_string, tokens, string(" "));
              if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"EOS_STATIC_KPPRA. Specify as EOS_STATIC_KPPRA=10000.");}
            }
            logger << (EOS_STATIC_KPPRA_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "EOS_STATIC_KPPRA=" << EOS_STATIC_KPPRA_OPTION.content_string << "." << apl::endl;
          }
          NEDOS_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "NEDOS=" + "|" + _ASTROPT_APL_OLD_ + "NEDOS="), NEDOS_OPTION.option, NEDOS_OPTION.xscheme);
          NEDOS = NEDOS_OPTION.content_int;
          if (NEDOS_OPTION.isentry) {
            tokens.clear();
            apl::tokenize(NEDOS_OPTION.content_string, tokens, string(" "));
            if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"NEDOS. Specify as NEDOS=5000.");}
          }
          logger << (NEDOS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "NEDOS=" << NEDOS_OPTION.content_string << "." << apl::endl;
        }
      }
    }
    //GP SUBDIRECTORY OPTIONS. These are automic options and not controlled by users
    CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GRUNEISEN_SD=" + "|" + _ASTROPT_APL_OLD_ + "GRUNEISEN_SD="), CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option, CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.xscheme);
    if(!CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option){
      CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GRUNEISEN_A_SD=" + "|" + _ASTROPT_APL_OLD_ + "GRUNEISEN_A_SD="), CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option, CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.xscheme);
      if(!CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option){
        CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GRUNEISEN_B_SD=" + "|" + _ASTROPT_APL_OLD_ + "GRUNEISEN_B_SD="), CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option, CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.xscheme);
        if(!CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option){
          CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GRUNEISEN_C_SD=" + "|" + _ASTROPT_APL_OLD_ + "GRUNEISEN_C_SD="), CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option, CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.xscheme);
        }
      }
    }

    //Writing to log
    if(CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GRUNEISEN_SD=" << (CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GRUNEISEN_A_SD=" << (CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GRUNEISEN_B_SD=" << (CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GRUNEISEN_C_SD=" << (CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    }

    if(CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option|| CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option || CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option || CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option){
      GP_DISTORTION_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "GP_DISTORTION=" + "|" + _ASTROPT_APL_OLD_ + "GP_DISTORTION="), GP_DISTORTION_OPTION.option, GP_DISTORTION_OPTION.xscheme);
      GP_DISTORTION=GP_DISTORTION_OPTION.content_double;
      if (GP_DISTORTION_OPTION.isentry) {
      tokens.clear();
        apl::tokenize(GP_DISTORTION_OPTION.content_string, tokens, string(" "));
        if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"GP_DISTORTION. Specify as GP_DISTORTION=0.03.");}
    }
      logger << (GP_DISTORTION_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GP_DISTORTION=" << GP_DISTORTION << "." << apl::endl;
    }
    //EOS SUBDIRECTORY OPTIONS. These are automic options and not controlled by users
    CALCULATE_EOS_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "EOS_SD=" + "|" + _ASTROPT_APL_OLD_ + "EOS_SD="), CALCULATE_EOS_SUBDIRECTORIES_OPTION.option, CALCULATE_EOS_SUBDIRECTORIES_OPTION.xscheme); //CO 170601
    if(CALCULATE_EOS_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_EOS_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "EOS_SD=" << (CALCULATE_EOS_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      EOS_DISTORTION_RANGE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "EOS_DISTORTION_RANGE=" + "|" + _ASTROPT_APL_OLD_ + "EOS_DISTORTION_RANGE="), EOS_DISTORTION_RANGE_OPTION.option, EOS_DISTORTION_RANGE_OPTION.xscheme); //CO 170601
      tokens.clear();
      apl::tokenize(EOS_DISTORTION_RANGE_OPTION.content_string, tokens, string(" :"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"EOS_DISTORTION_RANGE. Specify as EOS_DISTORTION_RANGE=-3:6:1");}
      EOS_DISTORTION_START = aurostd::string2utype<double>(tokens.at(0));
      EOS_DISTORTION_END = aurostd::string2utype<double>(tokens.at(1));
      EOS_DISTORTION_DISTORTION_INC = aurostd::string2utype<double>(tokens.at(2));
      logger << (EOS_DISTORTION_RANGE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "EOS_DISTORTION_RANGE=" << EOS_DISTORTION_RANGE_OPTION.content_string << "." << apl::endl;
      logger << "The EOS properties will be calculated in distortion range <" << EOS_DISTORTION_START << "," << EOS_DISTORTION_END << "," << EOS_DISTORTION_DISTORTION_INC << "." << apl::endl;
    }
    //}

    //SCQHA SUBDIRECTORY OPTIONS. These are automic options and not controlled by users
    CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "SCQHA_SD=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_SD="), CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option,  CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.xscheme);
    if(!CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option){
      CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "SCQHA_A_SD=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_A_SD="), CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option,  CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.xscheme);
      if(!CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option){
        CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "SCQHA_B_SD=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_B_SD="), CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option,  CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.xscheme);
        if(!CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option){
          CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "SCQHA_C_SD=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_C_SD="), CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.option,  CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.xscheme);
        }
      }
  }

    //Writing to log
    if(CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_SD=" << (CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_A_SD=" << (CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_B_SD=" << (CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_C_SD=" << (CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    }

    if(CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option || CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option || CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option || CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.isentry){
      SCQHA_DISTORTION_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "SCQHA_DISTORTION=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_DISTORTION="), SCQHA_DISTORTION_OPTION.option, SCQHA_DISTORTION_OPTION.xscheme);
      SCQHA_DISTORTION=SCQHA_DISTORTION_OPTION.content_double;
      if (SCQHA_DISTORTION_OPTION.isentry) {
        tokens.clear();
        apl::tokenize(SCQHA_DISTORTION_OPTION.content_string, tokens, string(" "));
        if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"SCQHA_DISTORTION. Specify as SCQHA_DISTORTION=3.0.");}
      }
      logger << (SCQHA_DISTORTION_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "SCQHA_DISTORTION=" << SCQHA_DISTORTION << "." << apl::endl;
    }
    // PINKU QUASI-HARMONIC END


//  ME181026 - OBSOLETE until it's properly documented
//    // Get the users maximum shell which will be included into calculation
//    // CO, not sure how maxshell works here (F option), need to investigate further and add to README
//    // also seems USER_WANTS_FULL_SHELL applies for both MAX and MIN shell settings, should one take precedence? should they be separate flags?
//    if(!found_supercell){
//      USER_MAXSHELL_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "MAXSHELL=" + "|" + _ASTROPT_QHA_ + "MAXSHELL=" + "|" + _ASTROPT_AAPL_ + "MAXSHELL=" + "|" + _ASTROPT_APL_OLD_ + "MAXSHELL="), USER_MAXSHELL_OPTION.option, USER_MAXSHELL_OPTION.xscheme); //CO 170601
//      test = USER_MAXSHELL_OPTION.content_string;
//      if (test[test.size() - 1] == 'f' || test[test.size() - 1] == 'F') {
//        USER_MAXSHELL = aurostd::string2utype<int>(test.substr(0, test.size() - 1));
//        USER_WANTS_FULL_SHELL = true;
//      } else {
//        USER_MAXSHELL = USER_MAXSHELL_OPTION.content_int;
//        USER_WANTS_FULL_SHELL = false;
//      }
//      if(USER_MAXSHELL_OPTION.isentry){
//        logger << (USER_MAXSHELL_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "MAXSHELL=" << USER_MAXSHELL << (USER_WANTS_FULL_SHELL ? " (FULL)" : "") << "." << apl::endl;
//        logger << "Supercell will be built with at most " << USER_MAXSHELL << " shells." << apl::endl;
//        found_supercell = true;
//      }
//    }
//
//    ME181026 - F option OBSOLETE until it's properly documented
//    // Get the users minimum shell which will be included into calculation
//    // CO, not sure how minshell works here (F option), need to investigate further and add to README
//    if(!found_supercell){
//      USER_MINSHELL_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "MINSHELL=" + "|" + _ASTROPT_QHA_ + "MINSHELL=" + "|" + _ASTROPT_AAPL_ + "MINSHELL=" + "|" + _ASTROPT_APL_OLD_ + "MINSHELL="), USER_MINSHELL_OPTION.option, USER_MINSHELL_OPTION.xscheme); //CO 170601
//      test = USER_MINSHELL_OPTION.content_string;
//      if (test[test.size() - 1] == 'f' || test[test.size() - 1] == 'F') {
//        USER_MINSHELL = aurostd::string2utype<int>(test.substr(0, test.size() - 1));
//        USER_WANTS_FULL_SHELL = true;
//      } else {
//        USER_MINSHELL = USER_MINSHELL_OPTION.content_int;
//        USER_WANTS_FULL_SHELL = false;
//      }
//      if(USER_MINSHELL_OPTION.isentry){
//        logger << (USER_MINSHELL_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "MINSHELL=" << USER_MINSHELL << (USER_WANTS_FULL_SHELL ? " (FULL)" : "") << "." << apl::endl;
//        logger << "Supercell will be built with at least " << USER_MINSHELL << " shells." << apl::endl;
//        found_supercell = true;
//      }
//    }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                          PREPARE CALCULATIONS                           //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

    //fix vasp bin for LR or DM+POLAR
    if (USER_ENGINE == string("LR") || (USER_ENGINE == string("DM") && USER_POLAR)) {
      if(xflags.AFLOW_MODE_VASP && !XHOST.GENERATE_AFLOWIN_ONLY){  // ME190313 - Do not check the VASP binary for generate_aflowin_only
  try {
        // Check the version of VASP binary
        logger << "Checking VASP version for linear response calculations.";
        string vaspVersion;
        vaspVersion = getVASPVersionString( (kflags.KBIN_MPI ? kflags.KBIN_MPI_BIN : kflags.KBIN_BIN ) );
        if (!vaspVersion.empty()) {
          logger << "[" << vaspVersion << "]";
          if ((vaspVersion[0] - '0') < 5) { //cool way of getting ascii value:  https://stackoverflow.com/questions/36310181/char-subtraction-in-c
            logger << apl::warning << "." << apl::endl;
            // if(_WITHIN_DUKE_){ OBSOLETE - ME 190108
              // ME 190107 - fix both serial and MPI binaries
              kflags.KBIN_SERIAL_BIN = DEFAULT_VASP5_BIN;
              kflags.KBIN_MPI_BIN = DEFAULT_VASP5_MPI_BIN;
              if (kflags.KBIN_MPI) {
                kflags.KBIN_BIN = kflags.KBIN_MPI_BIN;
            } else {
                kflags.KBIN_BIN = kflags.KBIN_SERIAL_BIN;
            }
              logger << apl::warning << "Modifying VASP bin to " << kflags.KBIN_BIN << " (AUTO modification)." << apl::endl;  // ME 190109
//[ME190109]            } else {
//[ME190109]              throw apl::APLRuntimeError("The LR engine needs VASP5 or higher.");
//[ME190109]            }
          } else {logger << " OK." << apl::endl;}
        } else {
          logger << "Failed." << apl::warning << apl::endl; 
          throw apl::APLLogicError("Unexpected binary format.");
        }
      } catch (apl::APLLogicError& e) {
        logger << apl::warning << "Failed to identify the version of VASP binary." << apl::endl;
        logger << apl::warning << e.what() << apl::endl;
      }
    }
    }

    // ME 190107 - Relax after fixing the vasp bin to make version consistent.
    // Run relaxations if necessary
    if (USER_RELAX) {
      bool Krun=true;
      string function;
      if (xinput.AFLOW_MODE_VASP) {
        function = "KBIN::relaxStructureAPL_VASP";
        Krun = relaxStructureAPL_VASP(START_RELAX, AflowIn, xinput.xvasp, aflags,
                                      kflags, xflags.vflags, messageFile);
      }
      if (!Krun) {
        string message = "Relaxation calculations did not run successfully.";
        throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
      }
    }

    // SUPERCELL ---------------------------------------------------------

    // Construct the working supercell ////////////////////////////////////

    //apl::Supercell supercell(xvasp.str,logger),supercell_test(xvasp.str,logger);    //corey, slow
  try {

    apl::Supercell supercell(xinput.getXStr(),aflags,logger); //xvasp.str, logger);  //CO  //CO181226
    apl::Supercell supercell_test = supercell;    //CO

    //   pflow::PrintDist(xinput.getXStr(),20.0,cerr);
    // ME 181026 - Added default case: if there are no supercell entries,
    // use MINATOMS (default value defined in aflow.rc)
    if ((USER_SUPERCELL.empty() && kflags.KBIN_MODULE_OPTIONS.supercell_method[1])) {
      stringstream aus;
      if(kflags.KBIN_MODULE_OPTIONS.minatoms_restricted){
        for (int Ni=1; USER_SUPERCELL == ""; Ni++) {
          aus.str("");
          aus << "Ni=" << Ni
              << " "
              << "supercell=" << Ni << "x" << Ni << "x" << Ni << "  natoms=" << Ni * Ni * Ni * xinput.getXStr().atoms.size(); //xvasp.str.atoms.size();
          //	  logger << aus.getXStr() << apl::endl;
          if (Ni * Ni * Ni * ((int)xinput.getXStr().atoms.size()) > (int)USER_MINATOMS) { // xvasp.str.atoms.size()) > (int)USER_MINATOMS) {
            USER_MINATOMS = 0;
            USER_SUPERCELL = aurostd::utype2string<uint>(Ni) + "X" + aurostd::utype2string<uint>(Ni) + "X" + aurostd::utype2string<uint>(Ni);
            logger << aus.str() << apl::endl;
          }
        }
      } else {
      for (double radius = 0.01; USER_SUPERCELL == ""; radius += 0.01) {
        xvector<int> dims(3);
        dims = LatticeDimensionSphere(xinput.getXStr().lattice,radius); //xvasp.str.lattice, radius);
          aus.str("");
        aus << "Radius=" << aurostd::PaddedPOST(aurostd::utype2string<double>(radius, 3), 4)
            << " "
            << " supercell=" << dims(1) << "x" << dims(2) << "x" << dims(3) << "  natoms=" << dims(1) * dims(2) * dims(3) * xinput.getXStr().atoms.size(); //xvasp.str.atoms.size();
        //	  logger << aus.getXStr() << apl::endl;
        if (dims(1) * dims(2) * dims(3) * ((int)xinput.getXStr().atoms.size()) > (int)USER_MINATOMS) { // xvasp.str.atoms.size()) > (int)USER_MINATOMS) {
          USER_MINATOMS = 0;
          USER_SUPERCELL = aurostd::utype2string<uint>(dims(1)) + "X" + aurostd::utype2string<uint>(dims(2)) + "X" + aurostd::utype2string<uint>(dims(3));
          logger << aus.str() << apl::endl;
          }
        }
      }
    }

    //      cerr << "USER_WANTS_FULL_SHELL=" << USER_WANTS_FULL_SHELL << endl;
    //    for(int i=2;i<20;i+=2) cerr << "try " << i << ": " << supercell_test.buildSuitableForShell(i,USER_WANTS_FULL_SHELL,FALSE) << endl;

    if (USER_SUPERCELL.empty() && kflags.KBIN_MODULE_OPTIONS.supercell_method[2]) {  // Not documented - ME181026
      logger << "a Searching for suitable cell to handle " << USER_MAXSHELL << " shells..." << apl::endl;
      bool USER_WANTS_FULL_SHELL = false;  // ME 181026
      supercell.buildSuitableForShell(USER_MAXSHELL, USER_WANTS_FULL_SHELL, TRUE);
      supercell.setupShellRestrictions(USER_MAXSHELL);
    } else if (USER_SUPERCELL.empty() && kflags.KBIN_MODULE_OPTIONS.supercell_method[3]) {
      bool USER_WANTS_FULL_SHELL = false;  // ME 181026
      logger << "b Searching for suitable cell to handle " << USER_MINSHELL << " shells..." << apl::endl;
      supercell.buildSuitableForShell(USER_MINSHELL, USER_WANTS_FULL_SHELL, TRUE);
    } else if (USER_SUPERCELL.find_first_of("xX") != string::npos) {
      // OK, user wants their own supercell...
      tokens.clear();
      apl::tokenize(USER_SUPERCELL, tokens, string(" xX"));
      supercell.build(aurostd::string2utype<int>(tokens.at(0)),
                      aurostd::string2utype<int>(tokens.at(1)),
                      aurostd::string2utype<int>(tokens.at(2)));
      // Did he specify also regular restriction for max shell included
      // in calculation?
      if (USER_MAXSHELL > 0)
        supercell.setupShellRestrictions(USER_MAXSHELL);
//    } else {  OBSOLETE ME181026
//      throw apl::APLRuntimeError("The settings for supercell construction are confusing.");
    }

    // CLUSTERS ---------------------------------------------------------

    // ME 180925
    // Calculate the clusters for thermal conductivity calculations
    vector<apl::ClusterSet> clusters;  // ME, default, only allocates to be passed into functions
    if (USER_TCOND) {
      int max_order;
      if (USER_AAPL_FOURTH_ORDER) {
        max_order = 4;
      } else {
        max_order = 3;
    }

      for (int o = 3; o <= max_order; o++) {
        apl::ClusterSet clst(logger);
        bool awakeClusterSet;
        string clust_hib_file = DEFAULT_AAPL_FILE_PREFIX + _CLUSTER_SET_FILE_[o-3];
        if (USER_HIBERNATE) {
          awakeClusterSet = (aurostd::EFileExist(clust_hib_file) ||
                             aurostd::FileExist(clust_hib_file));
        } else {
          awakeClusterSet = false;
        }

        if (awakeClusterSet) {
          try {
            clst = apl::ClusterSet(clust_hib_file, supercell, USER_CUTOFF_SHELL[o-3],
                                   USER_CUTOFF_DISTANCE[o-3], o, logger);
          } catch (aurostd::xerror excpt) {
            logger << apl::warning << excpt.where() << " " << excpt.error_message << std::endl;
            logger << apl::warning << "Skipping awakening of anharmonic IFCs." << apl::endl;
            awakeClusterSet = false;
          }
        }

        if (!awakeClusterSet) {
          clst = apl::ClusterSet(supercell, USER_CUTOFF_SHELL[o-3],
                                 USER_CUTOFF_DISTANCE[o-3], logger);
          clst.build(o);
          clst.buildDistortions();
          if (USER_HIBERNATE) {
            clst.writeClusterSetToFile(clust_hib_file);
          }
        }
        clusters.push_back(clst);
      }
    }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                           CALCULATE PHONONS                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

    auto_ptr<apl::PhononCalculator> phcalc;
    if (USER_ENGINE == string("DM")) {
      apl::DirectMethodPC* phcalcdm = new apl::DirectMethodPC(supercell, clusters, xinput, aflags,
                                                              kflags, xflags, AflowIn, logger);
      phcalcdm->isPolarMaterial(USER_POLAR);                                                       // TRY POLAR [STEFANO]
      phcalcdm->setTCOND(USER_TCOND);
      //phcalcdm->setGeneratePlusMinus(USER_DISTORTIONS_PLUS_MINUS_OPTION.option); //CO auto
      phcalcdm->setGeneratePlusMinus(USER_AUTO_DISTORTIONS, USER_DPM);  //CO auto
      phcalcdm->setGenerateOnlyXYZ(USER_DISTORTIONS_XYZ_ONLY);
      phcalcdm->setDistortionSYMMETRIZE(USER_DISTORTIONS_SYMMETRIZE);
      phcalcdm->setDistortionINEQUIVONLY(USER_DISTORTIONS_INEQUIVONLY); //CO190131
      phcalcdm->setDistortionMagnitude(USER_DISTORTION_MAGNITUDE);
      phcalcdm->setCalculateZeroStateForces(USER_ZEROSTATE);
      phcalcdm->get_special_inputs(AflowIn);  //PINKU, to include PSTRESS and LDAU_PARAMETERS in the SUPERCELL files
      phcalc.reset(phcalcdm);
    //CO generally redirects to DM, the distinction between DM and GSA is obsolete
    //} else if (USER_ENGINE == string("GSA")) {
    //  apl::GeneralizedSupercellApproach* gsa = new apl::GeneralizedSupercellApproach(supercell, strPair, xinput, aflags, kflags, xflags, logger);//xvasp, aflags, kflags, vflags, logger);  //Modified JJPR
    //  //gsa->setGeneratePlusMinus(USER_DISTORTIONS_PLUS_MINUS_OPTION.option); //CO auto
    //  gsa->setGeneratePlusMinus(AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option, USER_DISTORTIONS_PLUS_MINUS_OPTION.option);  //CO auto
    //  gsa->setGenerateOnlyXYZ(USER_DISTORTIONS_XYZ_ONLY_OPTION.option);
    //  gsa->setDistortionMagnitude(USER_DISTORTION_MAGNITUDE);
    //  gsa->setTensor(CALCULATE_TCOND_OPTION.option);  // TCOND JJPR
    //  gsa->setSumRule(USER_EPS_SUM);           // TCOND JJPR
    //  //phcalcdm->setCalculateZeroStateForces(USER_ZEROSTATE_OPTION.option);
    //  phcalc.reset(gsa);
    } else {
      phcalc.reset(new apl::LinearResponsePC(supercell, clusters, xinput, aflags,
                                             kflags, xflags, AflowIn, logger));
      phcalc->setTCOND(USER_TCOND);
      phcalc->isPolarMaterial(USER_POLAR);
      //phcalcdm->setCalculateZeroStateForces(USER_ZEROSTATE_OPTION.option);
    }

      //QHA/SCQHA/QHA3P  START //PN180705
      // Create directories for QHA/SCQHA/QHA3P
      // The pointer pheos should be called before creation of apl.xml
      auto_ptr<apl::QHA_AFLOWIN_CREATOR> pheos;
      if(CALCULATE_GRUNEISEN_OPTION.option || CALCULATE_GRUNEISEN_A_OPTION.option || CALCULATE_GRUNEISEN_B_OPTION.option || CALCULATE_GRUNEISEN_C_OPTION.option ||
         CALCULATE_SCQHA_OPTION.option || CALCULATE_SCQHA_A_OPTION.option || CALCULATE_SCQHA_B_OPTION.option || CALCULATE_SCQHA_C_OPTION.option ||
         CALCULATE_QHA3P_OPTION.option || CALCULATE_QHA3P_A_OPTION.option || CALCULATE_QHA3P_B_OPTION.option || CALCULATE_QHA3P_C_OPTION.option)
	{

        pheos.reset(new apl::QHA_AFLOWIN_CREATOR(supercell, clusters, xinput, aflags,
                                                 kflags, xflags, AflowIn, logger));

        pheos->setGP(CALCULATE_GRUNEISEN_OPTION.option, CALCULATE_GRUNEISEN_A_OPTION.option, CALCULATE_GRUNEISEN_B_OPTION.option, CALCULATE_GRUNEISEN_C_OPTION.option);
          if( CALCULATE_SCQHA_OPTION.option || CALCULATE_SCQHA_A_OPTION.option || CALCULATE_SCQHA_B_OPTION.option || CALCULATE_SCQHA_C_OPTION.option )
	    {
        pheos->setSCGP(CALCULATE_SCQHA_OPTION.option, CALCULATE_SCQHA_A_OPTION.option, CALCULATE_SCQHA_B_OPTION.option, CALCULATE_SCQHA_C_OPTION.option);
	    }
          if( CALCULATE_QHA3P_OPTION.option || CALCULATE_QHA3P_A_OPTION.option || CALCULATE_QHA3P_B_OPTION.option || CALCULATE_QHA3P_C_OPTION.option )
	    {
	      pheos->setSCGP(CALCULATE_QHA3P_OPTION.option, CALCULATE_QHA3P_A_OPTION.option, CALCULATE_QHA3P_B_OPTION.option, CALCULATE_QHA3P_C_OPTION.option);
	    }

        pheos->setGP_VOL_DISTORTION(GP_DISTORTION);
	  if(CALCULATE_SCQHA_OPTION.option || CALCULATE_SCQHA_A_OPTION.option || CALCULATE_SCQHA_B_OPTION.option || CALCULATE_SCQHA_C_OPTION.option ||
	     CALCULATE_QHA3P_OPTION.option || CALCULATE_QHA3P_A_OPTION.option || CALCULATE_QHA3P_B_OPTION.option|| CALCULATE_QHA3P_C_OPTION.option){
          pheos->setSCGP_VOL_DISTORTION(SCQHA_DISTORTION);
        }
        if(CALCULATE_EOS_OPTION.option){
      pheos->setEOS(CALCULATE_EOS_OPTION.option);
          pheos->setEOS_distortion_range(EOS_DISTORTION_START, EOS_DISTORTION_END, EOS_DISTORTION_DISTORTION_INC);
        pheos->setEOS_STATIC_KPPRA(EOS_STATIC_KPPRA);
          pheos->setEOS_NEDOS(NEDOS);
          pheos->set_edos_accurate(EDOS_ACURATE_OPTION.option);
      }
        pheos->run_qha();
        pheos->close_log();
    }
      //QHA/SCQHA/QHA3P END

    // ME180820 - set up VASP calculations for thermal conductivity calculations
    bool aapl_stagebreak;
    if (USER_TCOND) {
      aapl_stagebreak = phcalc->buildVaspAAPL(phcalc->_clusters[0]);
      if (USER_AAPL_FOURTH_ORDER) {
        aapl_stagebreak = (phcalc->buildVaspAAPL(phcalc->_clusters[1]) || aapl_stagebreak);
      }
    } else {
      aapl_stagebreak = false;
    }

    // Run or awake
    bool isHibFileAvailable = aurostd::EFileExist(DEFAULT_APL_FILE_PREFIX+DEFAULT_APL_HARMIFC_FILE);  //|| //CO
    //aurostd::FileExist(string("apl.xml")); //CO

    if (USER_HIBERNATE && isHibFileAvailable) {
      if (aapl_stagebreak) {
        throw apl::APLStageBreak();  // ME 180830
      }
      try {
        phcalc->awake();
      } catch (apl::APLLogicError& e) {
        logger << apl::warning << e.what() << apl::endl;
        logger << apl::warning << "Skipping awakening..." << apl::endl;
        isHibFileAvailable = false;
      }
    }

    if (!isHibFileAvailable) {
      phcalc->run(aapl_stagebreak);  // ME180830 -- added stagebreak bool
      if (USER_HIBERNATE)
        phcalc->hibernate();
    }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                            QUASI-HARMONIC                               //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

      //QHA/SCQHA/QHA3P START //PN180705
      //Store synamical matrics and PDOS from different distorted directores
      if(CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option)
        { 
          apl::QHAsubdirectoryData store(*phcalc, logger);
          store.setdir_prefix(_TMPDIR_);
      string dirname = store.getdir_name(aflags.Directory);
          store.set_gp_vol_distortion(GP_DISTORTION);
        //create uniform q-mesh
        store.createMPmesh(USER_DOS_MESH[0], USER_DOS_MESH[1], USER_DOS_MESH[2],
                           phcalc->getInputCellStructure());

          //check the distorted directort contains Gruneisen ON //PN180705
          if(store.check_GP()){ //PN180705
            //store dynamical matrices //PN180705
            store.create_dm(); //PN180705
        apl::PhononDispersionCalculator pdisc(*phcalc, logger);
        
        // Init path according to the aflow's definition for elec. struc.
        // ME 181029 - Restructured
        if (USER_DC_METHOD == "LATTICE") {
          pdisc.initPathLattice(USER_DC_INITLATTICE,USER_DC_NPOINTS);
        } else {
          if (!USER_DC_INITCOORDS_LABELS.empty() && !USER_DC_INITCOORDS_FRAC.empty()) {
            pdisc.initPathCoords(USER_DC_INITCOORDS_FRAC, USER_DC_INITCOORDS_LABELS, USER_DC_NPOINTS, false);
          } else if (!USER_DC_INITCOORDS_LABELS.empty() && !USER_DC_INITCOORDS_CART.empty()) {
            pdisc.initPathCoords(USER_DC_INITCOORDS_CART, USER_DC_INITCOORDS_LABELS, USER_DC_NPOINTS, true);
          }
        }
        // ME190501 Allow user to override path
          if(!USER_DC_USERPATH.empty()){  // Set path
            pdisc.setPath(USER_DC_USERPATH);
          }

        std::vector<xvector<double> > qpoints = pdisc.get_qpoints();
            //store dynamical matrices along path //PN180705
        store.create_pdispath(qpoints);
        qpoints.clear();
        pdisc.clear();
      }
          store.clear(); //PN180705
          phcalc->clear(); //PN180705
          return; //PN180705
        } //PN180705
      //store PDOS from different distorted directories
      if(CALCULATE_EOS_SUBDIRECTORIES_OPTION.option)
        {
          apl::QHAsubdirectoryData store(*phcalc, logger);
          store.setdir_prefix(_TMPDIR_);
          string dirname=store.getdir_name(aflags.Directory);
          // ME190428 - START
          // MonkhorstPackMesh replaced by qmesh
          //apl::MonkhorstPackMesh qmesh(USER_DOS_MESH[0], USER_DOS_MESH[1], USER_DOS_MESH[2],
          //  phcalc->getInputCellStructure(), logger);
         apl::QMesh qmesh(USER_DOS_MESH, phcalc->getInputCellStructure(), logger);
         qmesh.makeIrreducible();

        // OBSOLETE - DOSCalculator is not an auto_ptr anymore
        //auto_ptr<apl::DOSCalculator> dosc;

        //if (USER_DOS_METHOD == string("LT")) {
        //  apl::LTMethod _lt(qmesh, logger);
          //dosc.reset(new apl::LinearTetrahedronMethod(*phcalc, qmesh, logger));
        //  dosc.reset(new apl::DOSLinearTetrahedronMethod(_lt, *phcalc, qmesh, logger));
        //} else {
          //dosc.reset(new apl::RootSamplingMethod(*phcalc, qmesh, logger));
        //  dosc.reset(new apl::DOSRootSamplingMethod(*phcalc, qmesh, logger));
        //}
        
        apl::DOSCalculator dosc(*phcalc, qmesh, logger, USER_DOS_METHOD);

        // ME190428 - END
        // Calculate DOS
        dosc.calc(USER_DOS_NPOINTS, USER_DOS_SMEAR);
        if (USER_DOS) dosc.writePDOS(_TMPDIR_, dirname);
        dosc.clear();
          qmesh.clear();
          store.clear();
          phcalc->clear();
          return;
        }
      //SCQHA and QHA3P save dynamical matrices and PDOS from different distorted directories
      if(CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.option)
        {
          {
            apl::QHAsubdirectoryData store(*phcalc, logger);
            store.setdir_prefix(_TMPDIR_);
            string dirname=store.getdir_name(aflags.Directory);
            store.set_sc_vol_distortion(SCQHA_DISTORTION);

            //create uniform q-mesh
            store.createMPmesh(USER_DOS_MESH[0], USER_DOS_MESH[1], USER_DOS_MESH[2],
                               phcalc->getInputCellStructure());

            if(store.check_SCQHA())
              {
                store.create_dm();
                apl::PhononDispersionCalculator pdisc(*phcalc,logger);

		// Init path according to the aflow's definition for elec. struc.
                // ME 181029 - Restructured
                if (USER_DC_METHOD == "LATTICE") {
                  pdisc.initPathLattice(USER_DC_INITLATTICE,USER_DC_NPOINTS);
                } else {
                  if (!USER_DC_INITCOORDS_LABELS.empty() && !USER_DC_INITCOORDS_FRAC.empty()) {
                    pdisc.initPathCoords(USER_DC_INITCOORDS_FRAC, USER_DC_INITCOORDS_LABELS, USER_DC_NPOINTS, false);
                  } else if (!USER_DC_INITCOORDS_LABELS.empty() && !USER_DC_INITCOORDS_CART.empty()) {
                    pdisc.initPathCoords(USER_DC_INITCOORDS_CART, USER_DC_INITCOORDS_LABELS, USER_DC_NPOINTS, true);
                  }
                }
                // ME190501 Allow user to override path
                if(!USER_DC_USERPATH.empty()) {  // Set path
                  pdisc.setPath(USER_DC_USERPATH);
                }

                std::vector< xvector<double> > qpoints=pdisc.get_qpoints();
                store.create_pdispath(qpoints);
                qpoints.clear();
                pdisc.clear();
              }
            {
              // ME190428 - START
              // MonkhorstPackMesh replaced by qmesh
              //apl::MonkhorstPackMesh qmesh(USER_DOS_MESH[0], USER_DOS_MESH[1], USER_DOS_MESH[2],
              //                             phcalc->getInputCellStructure(),logger);
              apl::QMesh qmesh(USER_DOS_MESH, phcalc->getInputCellStructure(), logger);
              qmesh.makeIrreducible();
              
              // OBSOLETE - DOSCalculator is not an auto_ptr anymore
              //auto_ptr<apl::DOSCalculator> dosc;
              //if( USER_DOS_METHOD == string("LT") ) {
              //  apl::LTMethod _lt(qmesh, logger);
                //dosc.reset( new apl::LinearTetrahedronMethod(*phcalc,qmesh,logger) );
              //  dosc.reset( new apl::DOSLinearTetrahedronMethod(_lt, *phcalc,qmesh,logger) );
              //} else {
                //dosc.reset( new apl::RootSamplingMethod(*phcalc,qmesh,logger) );
              //  dosc.reset( new apl::DOSRootSamplingMethod(*phcalc,qmesh,logger) );
              //}
              
              apl::DOSCalculator dosc(*phcalc, qmesh, logger, USER_DOS_METHOD);
              // ME190428 - END
              // Calculate DOS
              dosc.calc(USER_DOS_NPOINTS,USER_DOS_SMEAR);
              if(USER_DOS)dosc.writePDOS(_TMPDIR_, dirname);
              dosc.clear();
              qmesh.clear();
      }
      store.clear();
      phcalc->clear();
          }
      return;
    }
      //PINKU QHA/SCQHA/QHA3P  END


  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                           PHONON PROPERTIES                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

    // Get the format of frequency desired by user ///////////////////////

    apl::IPCFreqFlags frequencyFormat = apl::NONE;

    if (!USER_FREQFORMAT.empty()) {
     try {
      // Convert format to machine representation
      tokens.clear();
      apl::tokenize(USER_FREQFORMAT, tokens, string(" |:;,"));
      for (uint i = 0; i < tokens.size(); i++) {
        if (tokens.at(i) == string("OMEGA")) {
          frequencyFormat |= apl::OMEGA;
          continue;
        }
        if (tokens.at(i) == string("HERTZ")) {
          frequencyFormat |= apl::HERTZ;
          continue;
        } else if (tokens.at(i) == string("THZ")) {
          frequencyFormat |= apl::THZ;
          continue;
        } else if (tokens.at(i) == string("CM-1") || tokens.at(i) == string("RECIPROCAL_CM")) {
          frequencyFormat |= apl::RECIPROCAL_CM;
          continue;
        } else if (tokens.at(i) == string("MEV")) {
          frequencyFormat |= apl::MEV;
          continue;
        }
        if (tokens.at(i) == string("ALLOW_NEGATIVE")) {
          frequencyFormat |= apl::ALLOW_NEGATIVE;
          continue;
        }
      }
      // Check if there was specified unit keyword...
      if (((frequencyFormat & ~apl::OMEGA) & ~apl::ALLOW_NEGATIVE) == apl::NONE)
        throw apl::APLLogicError("The mishmash frequency format.");
     } catch (std::exception& e) {
       logger << apl::error << e.what() << apl::endl;
     }
    } else {
      frequencyFormat = apl::THZ | apl::ALLOW_NEGATIVE;
    }

      //high-symmery qpoint auto pointer [PINKU] //PN180705
      auto_ptr<apl::PhononHSQpoints> ptr_hsq;
      bool is_negative_freq=false;
      bool scqha_is_vol_err=false;
      //high-symmery qpoint auto pointer END [PINKU]

    // PHONON DISPERSIONS ---------------------------------------------------------

    if (USER_DC) {
      apl::PhononDispersionCalculator pdisc(*phcalc, logger);

      // Init path according to the aflow's definition for elec. struc.
      // ME 181029 - Restructured
      if (USER_DC_METHOD == "LATTICE") {
        pdisc.initPathLattice(USER_DC_INITLATTICE,USER_DC_NPOINTS);
      } else {
        if (!USER_DC_INITCOORDS_LABELS.empty() && !USER_DC_INITCOORDS_FRAC.empty()) {
          pdisc.initPathCoords(USER_DC_INITCOORDS_FRAC, USER_DC_INITCOORDS_LABELS, USER_DC_NPOINTS, false);
        } else if (!USER_DC_INITCOORDS_LABELS.empty() && !USER_DC_INITCOORDS_CART.empty()) {
          pdisc.initPathCoords(USER_DC_INITCOORDS_CART, USER_DC_INITCOORDS_LABELS, USER_DC_NPOINTS, true);
        }
      }
      // ME190501 Allow user to override path
      if(!USER_DC_USERPATH.empty()){  // Set path
        pdisc.setPath(USER_DC_USERPATH);
      }

      // Calculate frequencies on path
      pdisc.calc(frequencyFormat);

      // Write results into PDIS file
      pdisc.writePDIS();
      pdisc.writePHEIGENVAL();  // ME190614
	//QHA/SCQHA/QHA3P  START //PN180705
	//////////////////////////////////////////////////////////////////////
        ptr_hsq.reset(new apl::PhononHSQpoints(logger));
        ptr_hsq->read_qpointfile();
      //compute Gruneisen dispersion curve
        if (CALCULATE_GRUNEISEN_OPTION.option ||
            CALCULATE_GRUNEISEN_A_OPTION.option ||
            CALCULATE_GRUNEISEN_B_OPTION.option ||
            CALCULATE_GRUNEISEN_C_OPTION.option)
          {
            apl::QHA qha(*phcalc, *pheos, logger);
            qha.get_tmp_dir_name(_TMPDIR_);
            qha.set_cutoff_freq(CUTOFF_FREQ);
            if(qha.set_imported_variables())
              {
                if(qha.calculation_gruneisen(ptr_hsq->get_qpoints()))
                  {
                    qha.write_gruneisen_parameter_path(ptr_hsq->get_path(), ptr_hsq->get_path_segment());
                    is_negative_freq=qha.get_is_negative_freq();
        }
        //[OBSOLETE PN180705]//clear used variables
        //[OBSOLETE PN180705]path.clear();
        //[OBSOLETE PN180705]path_segment.clear();
        //[OBSOLETE PN180705]qpoints.clear();
        //[OBSOLETE PN180705]qh->clear();
      }
            qha.clear();
          }
	//QHA/SCQHA/QHA3P  END
    }

    // PHONON DOS AND THERMODYNAMIC PROPERTIES ----------------------------------------

    if (USER_DOS || USER_TP) {
      // Generate mesh for calculation of DOS...
      // ME190428 - START
      // MonkhorstPackMesh replaced by qmesh
      //apl::MonkhorstPackMesh qmesh(USER_DOS_MESH[0], USER_DOS_MESH[1], USER_DOS_MESH[2],
      //                             phcalc->getInputCellStructure(), logger);

      apl::QMesh qmesh(USER_DOS_MESH, phcalc->getInputCellStructure(), logger);
      qmesh.makeIrreducible();
      // Setup the DOS engine which is used also for thermodynamic properties
      // OBSOLETE - DOSCalculator is not an auto_ptr anymore
      //auto_ptr<apl::DOSCalculator> dosc;
      //if (USER_DOS_METHOD == string("LT")) {
      //  apl::LTMethod _lt(qmesh, logger);
        //dosc.reset(new apl::LinearTetrahedronMethod(*phcalc, qmesh, logger));
      //  dosc.reset(new apl::DOSLinearTetrahedronMethod(_lt, *phcalc, qmesh, logger));
      //} else {
        //dosc.reset(new apl::RootSamplingMethod(*phcalc, qmesh, logger));
      //  dosc.reset(new apl::DOSRootSamplingMethod(*phcalc, qmesh, logger));
      //}

      // Calculate DOS
      apl::DOSCalculator dosc(*phcalc, qmesh, logger, USER_DOS_METHOD);
      // ME190428 - END
      dosc.calc(USER_DOS_NPOINTS, USER_DOS_SMEAR);
      if (USER_DOS) {
        dosc.writePDOS();
        dosc.writePHDOSCAR();  // ME190614
      }

      // Calculate thermal properties
      if (USER_TP) {
        if (!dosc.hasNegativeFrequencies()) {  // ME190423
          apl::ThermalPropertiesCalculator tpc(dosc, logger);  // ME190423
          tpc.writeTHERMO(USER_TP_TSTART, USER_TP_TEND, USER_TP_TSTEP);
	    //QHA/SCQHA/QHA3P START //PN180705
          //calculate Gruneisen
          // ME 190428 - START
            // Uniform Mesh is obsolete - replaced by QMesh
            //vector<int> sc_size(3,0);
            //for (int i = 0; i < 3; i++) sc_size[i] = USER_DOS_MESH[i];
            //apl::UniformMesh umesh(logger);
            //calculate group velocities
            //umesh.create_uniform_mesh(sc_size[0],sc_size[1],sc_size[2],phcalc->getInputCellStructure());
            if(CALCULATE_GROUPVELOCITY_OPTION.option){
              //apl::GroupVelocity vg(*phcalc, umesh, logger);
              apl::GroupVelocity vg(*phcalc, qmesh, logger);
              if(vg.check_negative_frequencies()){
                vg.write();
                vg.clear();
              }}
                //atomic displacement calculations
            if(CALCULATE_ATOMIC_DISPLACEMENT_OPTION.option){
              // Uniform Mesh is obsolete - replaced by QMesh
              //apl::AtomicDisplacements ad(*phcalc, umesh, logger);
              apl::AtomicDisplacements ad(*phcalc, qmesh, logger);
              ad.set_frequency_cutoff(CUTOFF_FREQ);
              ad.populate_variables(phcalc->getInputCellStructure());
              if(ad.eigen_solver()){
                ad.thermal_displacements(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP);
                ad.write_normal_mode_direction(ptr_hsq->get_hs_kpoints());
              }}
          // ME 190428 - END
            //QHA calculate Gruneisen
            std::vector< std::vector< double> > scqha_tv;
            if (CALCULATE_GRUNEISEN_OPTION.option   ||
                CALCULATE_GRUNEISEN_A_OPTION.option ||
                CALCULATE_GRUNEISEN_B_OPTION.option ||
                CALCULATE_GRUNEISEN_C_OPTION.option) {
              if(!is_negative_freq){
                // ME190428 - UniformMesh is obsolete - replaced by QMesh
                //if(umesh.get_kpoints().size()==0){
                //  umesh.create_uniform_mesh(sc_size[0],sc_size[1],sc_size[2],phcalc->getInputCellStructure());
                //}
                apl::QHA qha(*phcalc, *pheos, logger);
                qha.get_tmp_dir_name(_TMPDIR_);
                qha.set_cutoff_freq(CUTOFF_FREQ);
                if(qha.set_imported_variables())
                  {
                    //QHA Grunneisen parameter calculations     
                    //if(qha.calculation_gruneisen(&umesh))  OBSOLETE ME190428
                    if(qha.calculation_gruneisen(&qmesh))  // ME190428
                      {
                        qha.write_gruneisen_parameter_mesh();
                        qha.Writeaverage_gp(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP);
                        is_negative_freq=qha.get_is_negative_freq();
                  }
                  //[OBSOLETE PN180705]eos.clear();
                }
		if(CALCULATE_EOS_OPTION.option)
                  {
                    //QHA EOS calculations
                    //[CO181202 - NOT USED]apl::QH_ENERGIES eos_ens(*phcalc, *pheos, logger);
                    apl::QH_ENERGIES eos_ens(*pheos, logger);
                    eos_ens.get_tmp_dir_name(_TMPDIR_);
                    eos_ens.get_xtracture(phcalc->getInputCellStructure());
                    if(eos_ens.get_qha_energies())
                      {     
                        apl::QHAEOS qheos(qha, eos_ens, logger);
                        qheos.set_fitting_type(FITTING_TYPE);
                        if(qheos.setvariables())
                          { 
                            qheos.set_include_ele(INCLUDE_ELE_OPTION.option);
                            qheos.cal_qheos(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP, tpc);
              }
                        qheos.clear();
                        qha.clear();
                        //QHA3P and SCQHA calculations
                        if( (CALCULATE_SCQHA_OPTION.option || 
                             CALCULATE_SCQHA_A_OPTION.option || 
                             CALCULATE_SCQHA_B_OPTION.option || 
                             CALCULATE_SCQHA_C_OPTION.option ||
                             CALCULATE_QHA3P_OPTION.option || 
                             CALCULATE_QHA3P_A_OPTION.option || 
                             CALCULATE_QHA3P_B_OPTION.option || 
                             CALCULATE_QHA3P_C_OPTION.option) && (!scqha_is_vol_err) )
                          {         
                            if(!is_negative_freq){
                              apl::SCQHA_QHA3P scqha(*phcalc, *pheos, logger);
                              scqha.get_tmp_dir_name(_TMPDIR_);
                              scqha.set_cutoff_freq(CUTOFF_FREQ);
                              if(scqha.set_imported_variables())
                                {    
                                  //QHA3P Gruneisen parameter calculations 
                                  //if(scqha.calculation_gruneisen(&umesh))  OBSOLETE ME190428
                                  if(scqha.calculation_gruneisen(&qmesh))  // ME190428
                                    {
				      if(kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C)
					{ //PN 180719
                                      scqha.write_gruneisen_parameter_mesh();
                                      scqha.Writeaverage_gp(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP);
                                      is_negative_freq=scqha.get_is_negative_freq();
            }
				    }
            //[OBSOLETE PN180705]pheos->clear();
          }
                              if(!is_negative_freq) //PN180705
                                {   
                                  //SCQHA EOS
				  if(kflags.KBIN_PHONONS_CALCULATION_SCQHA || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A)
				    {
                                  apl::SCQHAEOS scqhaeos(scqha, eos_ens, logger);
                                  if(scqhaeos.import_variables())
                                    { 
                                      scqhaeos.set_input_temperature(scqha_pdis_T);
                                      scqhaeos.sccycle(USER_TP_TSTART,USER_TP_TEND, 0.1);
                                      scqha_tv=scqhaeos.get_TV_data();
                                    } 
                                  scqhaeos.clear();
				    }
                                  //QHA3P EOS
				  if(kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C)
				    {
                                  apl::QHA3POINTS qha3points(scqha, eos_ens, logger);
                                  if(qha3points.import_variables())
                                    {
                                      qha3points.set_include_ele(INCLUDE_ELE_OPTION.option);
                                      qha3points.qha3pts_temperature_loop(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP, tpc);
                                    }
                                  qha3points.clear();
                                }
                                }
                              eos_ens.clear();
                              scqha.clear();
                            }
                          }
                      }
                  }
	      } //QHA3P and SCQHA calculations    
	    } else if(kflags.KBIN_PHONONS_CALCULATION_SCQHA || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A ||
	             kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C)
	      {
              //[CO181202 - NOT USED]apl::QH_ENERGIES eos_ens(*phcalc, *pheos, logger);
              apl::QH_ENERGIES eos_ens(*pheos, logger);
              eos_ens.get_xtracture(phcalc->getInputCellStructure());
              eos_ens.get_tmp_dir_name(_TMPDIR_);
              if(eos_ens.get_scqha_energies())
                {
                  if(!is_negative_freq){
                    // ME190428 - UniformMesh is obsolete - replaced by QMesh
                    //if(umesh.get_kpoints().size()==0){
                    //  umesh.create_uniform_mesh(sc_size[0],sc_size[1],sc_size[2],phcalc->getInputCellStructure());
                    //}
                    apl::SCQHA_QHA3P scqha(*phcalc, *pheos, logger);
                    scqha.get_tmp_dir_name(_TMPDIR_);
                    scqha.set_cutoff_freq(CUTOFF_FREQ);
                    if(scqha.set_imported_variables())
                      {
                        //QHA3P Gruneisen parameter calculation
                        //if(scqha.calculation_gruneisen(&umesh))  OBSOLETE M190428
                        if(scqha.calculation_gruneisen(&qmesh))  // ME190428
                          {
			      if(kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C)
				{ //PN 180719
                            scqha.write_gruneisen_parameter_mesh();
                            scqha.Writeaverage_gp(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP);
                            is_negative_freq=scqha.get_is_negative_freq();
				}
                          }
                      }
                    if(!is_negative_freq)
                      {
                        //SCQHA calculations
			  if(kflags.KBIN_PHONONS_CALCULATION_SCQHA || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C)
			    {
                        apl::SCQHAEOS scqhaeos(scqha, eos_ens, logger);
                        if(scqhaeos.import_variables())
                          {
                            scqhaeos.sccycle(USER_TP_TSTART,USER_TP_TEND, 0.1);
                          }
                        scqhaeos.clear();
			    }
                        //QHA3P calculations
			  if(kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C)
			    {
                        apl::QHA3POINTS qha3points(scqha, eos_ens, logger);
                        if(qha3points.import_variables())
                          {
                            qha3points.set_include_ele(INCLUDE_ELE_OPTION.option);
                            qha3points.qha3pts_temperature_loop(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP, tpc);
                          }
                        qha3points.clear();
                      }
			}
                    scqha.clear();
                  }
                }
              eos_ens.clear();
	    }
	    //compute SCQHA temperature dependent dispersion curve
	    if((scqha_tv.size()!=0) && (SCQHA_PDIS_T_OPTION.option)){
	      apl::T_spectra_SCQHA_QHA3P scqha_T(*phcalc, *pheos, logger);
	      scqha_T.get_tmp_dir_name(_TMPDIR_);
	      scqha_T.set_cutoff_freq(CUTOFF_FREQ);
	      scqha_T.get_input_data(scqha_tv);
	      if(scqha_T.set_imported_variables())
		{
		  if(scqha_T.calculation_freqs(ptr_hsq->get_qpoints()))
		    {
		      scqha_T.calculate_pdis_T(ptr_hsq->get_path(), ptr_hsq->get_path_segment());
		    }
		}
	      scqha_T.clear();
	    }
            if(ptr_hsq.get()) ptr_hsq->clear();  // ME190423
            // umesh.clear();  OBSOLETE ME190428
	    //QHA/SCQHA/QHA3P END
          tpc.clear();
        } else {
          logger << apl::warning << "There are negative frequencies in DOS. The calculation of thermal properties has been skipped." << apl::endl;
        }
      }

      // Clear old stuff
      dosc.clear();  // ME190423
      //delete dosc; //auto_ptr will do
      qmesh.clear();
    }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                         THERMAL CONDUCTIVITY                            //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

    if (USER_TCOND) {
      // Get anharmonic force constants
      phcalc->setAnharmonicOptions(USER_AAPL_MAX_ITER, USER_AAPL_MIX, USER_EPS_SUM);
      bool awakeAnharmIFCs;
      for (uint i = 0; i < phcalc->_clusters.size(); i++) {
        string ifcs_hib_file = DEFAULT_AAPL_FILE_PREFIX + _ANHARMONIC_IFCS_FILE_[i];
        if (USER_HIBERNATE) {
          awakeAnharmIFCs = (aurostd::EFileExist(ifcs_hib_file) ||
                             aurostd::FileExist(ifcs_hib_file));
        } else {
          awakeAnharmIFCs = false;
        }

        if (awakeAnharmIFCs) {
          try {
            phcalc->readAnharmonicIFCs(ifcs_hib_file, phcalc->_clusters[i]);
          } catch (aurostd::xerror excpt) {
            logger << apl::warning << excpt.where() << " " << excpt.error_message << std::endl;
            logger << apl::warning << "Skipping awakening of anharmonic IFCs." << apl::endl;
            awakeAnharmIFCs = false;
          }
        }

        if (!awakeAnharmIFCs) {
          phcalc->calculateAnharmonicIFCs(phcalc->_clusters[i]);
          if (USER_HIBERNATE) {
            phcalc->_anharmonicIFCs[i].writeIFCsToFile(ifcs_hib_file);
          }
        }
      }

      // Do the thermal conductivity calculation
      logger << "Starting thermal conductivity calculations." << apl::endl;
      apl::TCONDCalculator tcond(*phcalc, supercell, logger);

      tcond.setCalculationOptions(USER_BTE, USER_ISOTOPE,
                                  USER_CUMULATIVEK, USER_AAPL_FOURTH_ORDER,
                                  USER_BOUNDARY, USER_NANO_SIZE, USER_TCT_TSTART,
                                  USER_TCT_TEND, USER_TCT_TSTEP);

      // Get q-points
      logger << "Preparing a q-mesh of " << USER_THERMALGRID[0] << "x" << USER_THERMALGRID[1] << "x" << USER_THERMALGRID[2] << "." << apl::endl;
      tcond.buildQpoints(aurostd::vector2xvector(USER_THERMALGRID));

      // Calculate lattice thermal conductivity
      tcond.calculateFrequenciesGroupVelocities();
      tcond.calculateTransitionProbabilities(3);
      if (USER_AAPL_FOURTH_ORDER) {
        tcond.calculateTransitionProbabilities(4);
    }
      tcond.calculateThermalConductivity();
 
      tcond.clear();
    }
    /*************** End Thermal Conductivity Calculations ****************/

    phcalc->clear();
    // delete phcalc; //auto_ptr will do
    supercell.clear();
  } catch (apl::APLStageBreak& e) {
    logger << apl::notice << "Stopped. Waiting for required calculations..." << apl::endl;  //CO181226
  }
}
}

#ifndef COMPILE_SLIM
void RunPhonons_APL_180101(_xinput& xinput,
                    string AflowIn,
                    _aflags& aflags,
                    _kflags& kflags,
                    _xflags& xflags, 
                    ofstream& messageFile) {
  // Test
    //if (!(kflags.KBIN_PHONONS_CALCULATION_APL || kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_AAPL)) return; //PN180705
    if (!(kflags.KBIN_PHONONS_CALCULATION_APL ||
          kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_QHA_A || kflags.KBIN_PHONONS_CALCULATION_QHA_B || kflags.KBIN_PHONONS_CALCULATION_QHA_C ||
          kflags.KBIN_PHONONS_CALCULATION_SCQHA || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C ||
          kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C || //PN180717
          kflags.KBIN_PHONONS_CALCULATION_AAPL)) return; //PN180705

  //we make certain automatic fixes if we're within our domain, otherwise issue warning/error
  ::_WITHIN_DUKE_ = (aurostd::substring2bool(XHOST.hostname, "nietzsche") || aurostd::substring2bool(XHOST.hostname, "aflowlib") || aurostd::substring2bool(XHOST.hostname, "qrats") || aurostd::substring2bool(XHOST.hostname, "habana") || aurostd::substring2bool(XHOST.hostname, "quser"));

  //determine if we have a consistent mode defined between input and flags
  if(xinput.AFLOW_MODE_VASP){
    if(!xflags.AFLOW_MODE_VASP){
      cerr << "ERROR: KBIN::RunPhonons_APL: mismatch types between xinput(VASP) and xflags(!VASP)" << endl;
      return;
    }
  } else if(xinput.AFLOW_MODE_AIMS){
    if(!xflags.AFLOW_MODE_AIMS){
      cerr << "ERROR: KBIN::RunPhonons_APL: mismatch types between xinput(AIMS) and xflags(!AIMS)" << endl;
      return;
    }
  //} else if(xinput.AFLOW_MODE_ALIEN){  //alien doesn't have xstr, so we ignore
  } else {
    cerr << "ERROR: KBIN::RunPhonons_APL: unknown input type" << endl;
    return;
  }

  //corey
  //fix names regardless of POSCAR style, we need mass!
  if(xinput.AFLOW_MODE_VASP){pflow::fixEmptyAtomNames(xinput.xvasp.str, true);}

  // Setup our logger
  apl::Logger logger(messageFile, aflags);
  //logger.setModuleName("PHONONS");  //will rename later
  
  string _ASTROPT_; //CO 170601
  if (kflags.KBIN_PHONONS_CALCULATION_AAPL) {
    logger.setModuleName("AAPL");  //CO 170601
    _ASTROPT_ = _ASTROPT_AAPL_;    //CO 170601
    } else if (kflags.KBIN_PHONONS_CALCULATION_QHA || 
               kflags.KBIN_PHONONS_CALCULATION_QHA_A || 
               kflags.KBIN_PHONONS_CALCULATION_QHA_B || 
               kflags.KBIN_PHONONS_CALCULATION_QHA_C || 
               kflags.KBIN_PHONONS_CALCULATION_SCQHA ||
               kflags.KBIN_PHONONS_CALCULATION_SCQHA_A ||
               kflags.KBIN_PHONONS_CALCULATION_SCQHA_B ||
               kflags.KBIN_PHONONS_CALCULATION_SCQHA_C ||
               kflags.KBIN_PHONONS_CALCULATION_QHA3P ||
               kflags.KBIN_PHONONS_CALCULATION_QHA3P_A ||
               kflags.KBIN_PHONONS_CALCULATION_QHA3P_B ||
               kflags.KBIN_PHONONS_CALCULATION_QHA3P_C) { //PN180705
    logger.setModuleName("QHA");  //CO 170601
    _ASTROPT_ = _ASTROPT_QHA_;    //CO 170601
  } else {
    logger.setModuleName("APL");  //CO 170601
    _ASTROPT_ = _ASTROPT_APL_;    //CO 170601
  }

  logger << "RUNNING..." << apl::endl;

  // CONTROL PARAMETERS FROM _AFLOWIN_ --------------------------------------

  // General switches what to calculate
  aurostd::xoption CALCULATE_PHONON_DISPERSIONS_OPTION; CALCULATE_PHONON_DISPERSIONS_OPTION.option = DEFAULT_APL_DC;
  aurostd::xoption CALCULATE_PHONON_DOS_OPTION; CALCULATE_PHONON_DOS_OPTION.option = DEFAULT_APL_DOS;
  aurostd::xoption CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION; CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option = DEFAULT_APL_TP;
  aurostd::xoption CALCULATE_POLAR_CORRECTIONS_OPTION; CALCULATE_POLAR_CORRECTIONS_OPTION.option = DEFAULT_APL_POLAR;
  // BEGIN ME - AAPL
  aurostd::xoption CALCULATE_APL_OPTION; CALCULATE_APL_OPTION.option = false;
  aurostd::xoption CALCULATE_TCOND_OPTION; CALCULATE_TCOND_OPTION.option = false;
  aurostd::xoption USER_BTE_OPTION; USER_BTE_OPTION.xscheme = DEFAULT_AAPL_BTE; string USER_BTE = DEFAULT_AAPL_BTE;
  aurostd::xoption CALCULATE_ISOTOPE_OPTION; CALCULATE_ISOTOPE_OPTION.option = DEFAULT_AAPL_ISOTOPE;
  aurostd::xoption CALCULATE_CUMULATIVEK_OPTION; CALCULATE_CUMULATIVEK_OPTION.option = DEFAULT_AAPL_CUMULATIVEK;
  aurostd::xoption CALCULATE_BOUNDARY_OPTION; CALCULATE_BOUNDARY_OPTION.option = DEFAULT_AAPL_CUMULATIVEK;
//  aurostd::xoption USER_TDISTORTION_MAGNITUDE_OPTION; USER_TDISTORTION_MAGNITUDE_OPTION.xscheme = "0.015"; double USER_TDISTORTION_MAGNITUDE = 0.015;  OBSOLETE ME 181018
  aurostd::xoption USER_THERMALGRID_OPTION; USER_THERMALGRID_OPTION.xscheme = DEFAULT_AAPL_THERMALGRID; string USER_THERMALGRID = DEFAULT_AAPL_THERMALGRID;
  aurostd::xoption USER_CUTOFF_DISTANCE_OPTION; USER_CUTOFF_DISTANCE_OPTION.xscheme = DEFAULT_AAPL_CUT_RAD; vector<double> USER_CUTOFF_DISTANCE; //CO 180409
  vector<string> default_tokens;
  default_tokens.clear();
  apl::tokenize(DEFAULT_AAPL_CUT_RAD, default_tokens, string(" ,"));
  for (uint i = 0; i < default_tokens.size(); i++) USER_CUTOFF_DISTANCE.push_back(aurostd::string2utype<double>(default_tokens[i]));
  stringstream default_stream;
  default_stream << DEFAULT_AAPL_NANO_SIZE;
  aurostd::xoption USER_NANO_SIZE_OPTION; USER_NANO_SIZE_OPTION.xscheme = default_stream.str(); double USER_NANO_SIZE = DEFAULT_AAPL_NANO_SIZE;
  default_stream.str("");
  default_stream << DEFAULT_AAPL_SUMRULE;
  aurostd::xoption USER_EPS_SUM_OPTION; USER_EPS_SUM_OPTION.xscheme = default_stream.str(); double USER_EPS_SUM = DEFAULT_AAPL_SUMRULE;
  default_stream.str("");
  default_stream << DEFAULT_AAPL_SUMRULE_MAX_ITER;
  aurostd::xoption USER_AAPL_MAX_ITER_OPTION; USER_AAPL_MAX_ITER_OPTION.xscheme = default_stream.str(); int USER_AAPL_MAX_ITER = DEFAULT_AAPL_SUMRULE_MAX_ITER;
  default_stream.str("");
  default_stream << DEFAULT_AAPL_MIXING_COEFFICIENT;
  aurostd::xoption USER_AAPL_MIX_OPTION; USER_AAPL_MIX_OPTION.xscheme = default_stream.str(); double USER_AAPL_MIX = DEFAULT_AAPL_MIXING_COEFFICIENT;
  aurostd::xoption USER_AAPL_FOURTH_ORDER_OPTION; USER_AAPL_FOURTH_ORDER_OPTION.option = DEFAULT_AAPL_FOURTH_ORDER;
  aurostd::xoption USER_CUTOFF_SHELL_OPTION; USER_CUTOFF_SHELL_OPTION.xscheme = DEFAULT_AAPL_CUT_SHELL; vector<int> USER_CUTOFF_SHELL; //CO 180409
  default_tokens.clear();
  apl::tokenize(DEFAULT_AAPL_CUT_SHELL, default_tokens, string(" ,"));
  for (uint i = 0; i < default_tokens.size(); i++) USER_CUTOFF_SHELL.push_back(aurostd::string2utype<int>(default_tokens[i]));
  aurostd::xoption USER_TCT_OPTION; USER_TCT_OPTION.xscheme = DEFAULT_AAPL_TCT;
  default_tokens.clear(); apl::tokenize(DEFAULT_AAPL_TCT, default_tokens, string(" :"));
  double USER_TCT_TSTART = aurostd::string2utype<double>(default_tokens[0]);
  double USER_TCT_TEND = aurostd::string2utype<double>(default_tokens[1]);
  double USER_TCT_TSTEP = aurostd::string2utype<double>(default_tokens[2]);
  // END ME
  
  //PINKU QUASI-HARMONIC START
    aurostd::xoption CALCULATE_GROUPVELOCITY_OPTION; CALCULATE_GROUPVELOCITY_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_ATOMIC_DISPLACEMENT_OPTION; CALCULATE_ATOMIC_DISPLACEMENT_OPTION.option = false; //PN180705
  aurostd::xoption CALCULATE_GRUNEISEN_OPTION; CALCULATE_GRUNEISEN_OPTION.option = false;
  aurostd::xoption CALCULATE_DISPLACEMENTS_OPTION; CALCULATE_DISPLACEMENTS_OPTION.option = false;
  aurostd::xoption CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION; CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option = false;
  aurostd::xoption CALCULATE_EOS_OPTION; CALCULATE_EOS_OPTION.option = false;
  aurostd::xoption CALCULATE_EOS_SUBDIRECTORIES_OPTION; CALCULATE_EOS_SUBDIRECTORIES_OPTION.option = false;
    aurostd::xoption EDOS_ACURATE_OPTION; EDOS_ACURATE_OPTION.option = false; //PN180705
    aurostd::xoption INCLUDE_ELE_OPTION;  INCLUDE_ELE_OPTION.option = false; //PN180705
    //Anisotropic Gruneisen and EOS //PN180705
    //in the a-direction //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_A_OPTION; CALCULATE_GRUNEISEN_A_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION; CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option = false; //PN180705
    aurostd::xoption GP_DISTORTION_OPTION; GP_DISTORTION_OPTION.xscheme = "0.03"; double GP_DISTORTION = 0.03; //PN180705

    //in the b-direction //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_B_OPTION; CALCULATE_GRUNEISEN_B_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION; CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option = false; //PN180705

    //in the c-direction //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_C_OPTION; CALCULATE_GRUNEISEN_C_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION; CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option = false; //PN180705

    //SC-QHA //PN180705
    aurostd::xoption CALCULATE_SCQHA_OPTION; CALCULATE_SCQHA_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_SCQHA_SUBDIRECTORIES_OPTION; CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option = false; //PN180705
    //QHA3P option
    aurostd::xoption CALCULATE_QHA3P_OPTION; CALCULATE_QHA3P_OPTION.option = false; //PN180705

    //Anisotropic SCQHA EOS //PN180705
    //in the a-direction //PN180705
    aurostd::xoption CALCULATE_SCQHA_A_OPTION; CALCULATE_SCQHA_A_OPTION.option = false; //PN180705
    aurostd::xoption SCQHA_DISTORTION_OPTION; SCQHA_DISTORTION_OPTION.xscheme = "3.0"; double SCQHA_DISTORTION = 3.0; //PN180705
    aurostd::xoption CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION; CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option = false; //PN180705

    //in the b-direction //PN180705
    aurostd::xoption CALCULATE_SCQHA_B_OPTION; CALCULATE_SCQHA_B_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION; CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option = false; //PN180705

    //in the c-direction //PN180705
    aurostd::xoption CALCULATE_SCQHA_C_OPTION; CALCULATE_SCQHA_C_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION; CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.option = false; //PN180705

    //QHA3P a direction
    aurostd::xoption CALCULATE_QHA3P_A_OPTION; CALCULATE_QHA3P_A_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_QHA3P_A_SUBDIRECTORIES_OPTION; CALCULATE_QHA3P_A_SUBDIRECTORIES_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_QHA3P_B_OPTION; CALCULATE_QHA3P_B_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_QHA3P_B_SUBDIRECTORIES_OPTION; CALCULATE_QHA3P_B_SUBDIRECTORIES_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_QHA3P_C_OPTION; CALCULATE_QHA3P_C_OPTION.option = false; //PN180705
    aurostd::xoption CALCULATE_QHA3P_C_SUBDIRECTORIES_OPTION; CALCULATE_QHA3P_C_SUBDIRECTORIES_OPTION.option = false; //PN180705


    //QHA, QHA3P and SCQHA options initializing from previous options
    CALCULATE_GRUNEISEN_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA;
    CALCULATE_GRUNEISEN_A_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA_A;
    CALCULATE_GRUNEISEN_B_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA_B;
    CALCULATE_GRUNEISEN_C_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA_C;

    CALCULATE_SCQHA_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_SCQHA;
    CALCULATE_SCQHA_A_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_SCQHA_A;
    CALCULATE_SCQHA_B_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_SCQHA_B;
    CALCULATE_SCQHA_C_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_SCQHA_C;

    CALCULATE_QHA3P_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA3P;
    CALCULATE_QHA3P_A_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA3P_A;
    CALCULATE_QHA3P_B_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA3P_B;
    CALCULATE_QHA3P_C_OPTION.option=kflags.KBIN_PHONONS_CALCULATION_QHA3P_C;

  // different types of fitting options for EOS calculations
  // (1) BM1 => Murnaghan EOS
  // (2) BM2 => Birch-Murnaghan 3rd-order EOS
  // (3) BM3 => Birch-Murnaghan 4th-order EOS
  //[OBSOLETE PN180705]aurostd::xoption GP_VOL_DISTORTION_OPTION; GP_VOL_DISTORTION_OPTION.xscheme = "0.03"; double GP_VOL_DISTORTION = 0.03;
  aurostd::xoption USER_PROJECTION_DIR_OPTION; USER_PROJECTION_DIR_OPTION.xscheme = "1:1:1"; vector<double> directions(3, 0); directions[0] = 1; directions[1] = 1; directions[2] = 1;   // 3 Miller indices
  aurostd::xoption CUTOFF_FREQ_OPTION; CUTOFF_FREQ_OPTION.xscheme="1e-5"; double CUTOFF_FREQ = 1e-5;  //in amu
    aurostd::xoption EOS_DISTORTION_RANGE_OPTION; EOS_DISTORTION_RANGE_OPTION.xscheme = "-3:6:1"; double EOS_DISTORTION_START = -3; double EOS_DISTORTION_END = 6; double EOS_DISTORTION_DISTORTION_INC=1; //PN180705
    aurostd::xoption EOS_STATIC_KPPRA_OPTION; EOS_STATIC_KPPRA_OPTION.xscheme = "10000"; int EOS_STATIC_KPPRA = 10000; //PN180705
    aurostd::xoption NEDOS_OPTION; NEDOS_OPTION.xscheme = "5000"; int NEDOS = 5000; //PN180705
  aurostd::xoption FITTING_TYPE_OPTION; FITTING_TYPE_OPTION.xscheme = "BM1"; string FITTING_TYPE = "BM1";
    aurostd::xoption SCQHA_PDIS_T_OPTION; SCQHA_PDIS_T_OPTION.xscheme = "100,400,600"; std::vector<double> scqha_pdis_T; //PN180705
  //PINKU QUASI-HARMONIC END
  
  // User's control about general phonon engine
  aurostd::xoption USER_ENGINE_OPTION; USER_ENGINE_OPTION.xscheme = DEFAULT_APL_ENGINE; string USER_ENGINE = DEFAULT_APL_ENGINE;
  aurostd::xoption AUTO_DISTORTIONS_PLUS_MINUS_OPTION; AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option = true;  //CO
  aurostd::xoption USER_DISTORTIONS_PLUS_MINUS_OPTION; USER_DISTORTIONS_PLUS_MINUS_OPTION.option = DEFAULT_APL_DPM;  //CO
  aurostd::xoption USER_DISTORTIONS_XYZ_ONLY_OPTION; USER_DISTORTIONS_XYZ_ONLY_OPTION.option = DEFAULT_APL_DXYZONLY;
  default_stream.str("");
  default_stream << DEFAULT_APL_DMAG;
  aurostd::xoption USER_DISTORTION_MAGNITUDE_OPTION; USER_DISTORTION_MAGNITUDE_OPTION.xscheme = default_stream.str(); double USER_DISTORTION_MAGNITUDE = DEFAULT_APL_DMAG;
  aurostd::xoption USER_ZEROSTATE_OPTION; USER_ZEROSTATE_OPTION.option = DEFAULT_APL_ZEROSTATE;
  aurostd::xoption USER_HIBERNATE_OPTION; USER_HIBERNATE_OPTION.option = DEFAULT_APL_HIBERNATE;
  
  aurostd::xoption USER_FREQFORMAT_OPTION; USER_FREQFORMAT_OPTION.xscheme=DEFAULT_APL_FREQFORMAT; string USER_FREQFORMAT = DEFAULT_APL_FREQFORMAT;
  // User's control of supercell used for calculation
  aurostd::xoption USER_WANTS_RELAX_OPTION; USER_WANTS_RELAX_OPTION.option = false;
  bool USER_WANTS_FULL_SHELL = false;
  //  bool   USER_WANTS_FULL_ATOMS            = false;
  aurostd::xoption USER_MAXSHELL_OPTION; USER_MAXSHELL_OPTION.xscheme = "-1"; int USER_MAXSHELL = -1;  //CO set by OPTION
  default_stream.str("");
  default_stream << DEFAULT_APL_MINSHELL;
  aurostd::xoption USER_MINSHELL_OPTION; USER_MINSHELL_OPTION.xscheme = default_stream.str(); int USER_MINSHELL = DEFAULT_APL_MINSHELL;  //CO set by OPTION
  default_stream.str("");
  default_stream << DEFAULT_APL_MINATOMS;
  aurostd::xoption USER_MINATOMS_OPTION; USER_MINATOMS_OPTION.xscheme = default_stream.str(); int USER_MINATOMS = DEFAULT_APL_MINATOMS;  //CO set by OPTION
  bool USER_MINATOMS_RESTICTED_OPTION=0;  //CO 180404 - adding restricted variant to minatoms
  aurostd::xoption USER_SUPERCELL_OPTION; USER_SUPERCELL_OPTION.xscheme = ""; string USER_SUPERCELL = "";

  // User's control about calculation of dispersion curves
  // there are several ways a user can change the default path
  // a USERPATH approach is either a subset or a swapping of the default path
  // OWNPATH is a true user defined path, either by providing frac or cartesian coordinates
  // OWNPATH can be used in combination with LABELS
  // useful for matching with experimental results
  aurostd::xoption USER_DC_INITLATTICE_OPTION; USER_DC_INITLATTICE_OPTION.xscheme = ""; string USER_DC_INITLATTICE = "";  //CO set by OPTION
  aurostd::xoption USER_DC_INITSG_OPTION; USER_DC_INITSG_OPTION.xscheme = ""; int USER_DC_DCINITSG = 1;
  aurostd::xoption USER_DC_INITCOORDS_FRAC_OPTION; USER_DC_INITCOORDS_FRAC_OPTION.xscheme = ""; string USER_DC_INITCOORDS_FRAC = ""; //CO set by OPTION
  aurostd::xoption USER_DC_INITCOORDS_CART_OPTION; USER_DC_INITCOORDS_CART_OPTION.xscheme = ""; string USER_DC_INITCOORDS_CART = ""; //CO set by OPTION
  aurostd::xoption USER_DC_INITCOORDS_LABELS_OPTION; USER_DC_INITCOORDS_LABELS_OPTION.xscheme = ""; string USER_DC_INITCOORDS_LABELS = ""; //CO set by OPTION
  aurostd::xoption USER_DC_USERPATH_OPTION; USER_DC_USERPATH_OPTION.xscheme = ""; string USER_DC_USERPATH = ""; //CO set by OPTION
  default_stream.str("");
  default_stream << DEFAULT_APL_DCPOINTS;
  aurostd::xoption USER_DC_NPOINTS_OPTION; USER_DC_NPOINTS_OPTION.xscheme = default_stream.str(); int USER_DC_NPOINTS = DEFAULT_APL_DCPOINTS;

  // User's control about phonon DOS
  aurostd::xoption USER_DOS_MESH_OPTION; USER_DOS_MESH_OPTION.xscheme = DEFAULT_APL_DOSMESH; string USER_DOS_MESH = DEFAULT_APL_DOSMESH;
  default_stream.str("");
  default_stream << DEFAULT_APL_DOSPOINTS;
  aurostd::xoption USER_DOS_NPOINTS_OPTION; USER_DOS_NPOINTS_OPTION.xscheme = default_stream.str(); int USER_DOS_NPOINTS = DEFAULT_APL_DOSPOINTS;
  aurostd::xoption USER_DOS_METHOD_OPTION; USER_DOS_METHOD_OPTION.xscheme = DEFAULT_APL_DOSMETHOD; string USER_DOS_METHOD = DEFAULT_APL_DOSMETHOD;
  aurostd::xoption USER_DOS_SMEAR_OPTION; USER_DOS_SMEAR_OPTION.xscheme = "0.0"; double USER_DOS_SMEAR = 0.0;

  // User's control of thermodynamic properties calculation
  aurostd::xoption USER_TPT_OPTION; USER_TPT_OPTION.xscheme=DEFAULT_APL_TPT;
  default_tokens.clear(); apl::tokenize(DEFAULT_APL_TPT, default_tokens, string(" :"));
  double USER_TP_TSTART = aurostd::string2utype<double>(default_tokens[0]);
  double USER_TP_TEND = aurostd::string2utype<double>(default_tokens[1]);
  double USER_TP_TSTEP = aurostd::string2utype<double>(default_tokens[2]);

  //if current dir DCUSERPATH is same as other PHONON sub-directories, if not make changes accordingly
  //PINKU QUASI-HARMONIC START
    //[OBSOLETE PN180705]apl::check_consistency_aflow_apl check_consistency(logger);
    //[OBSOLETE PN180705]check_consistency.getdir_name(aflags.Directory);
    //[OBSOLETE PN180705]check_consistency.check_consistency_in_aflow(AflowIn);
  //PINKU QUASI-HARMONIC END

  // TAR
  //  bool DO_TAR                             = false;

  // Get user's parameters from _AFLOWIN_ ///////////////////////////////////

  //for when we need double or string variable extraction
  vector<string> tokens;
  string test;
  bool override_option;
  bool found_supercell;
  try {
    // BEGIN JJPR
    // Anharmonic and thermal conductivity options

    if(kflags.KBIN_PHONONS_CALCULATION_AAPL){
    // TCOND, e.g. TCOND = y
      CALCULATE_TCOND_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_AAPL; //recycle what we parsed earlier
      logger << _ASTROPT_ << "CALC is" << ( CALCULATE_TCOND_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;  //CO 170601
      logger << "Anharmonic force constants will" << ( CALCULATE_TCOND_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-AAPL." << apl::endl;   //CO 170601

      // Get the user's magnitute of the threshold for the sumrule
      USER_EPS_SUM_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "SUMRULE=" + "|" + _ASTROPT_APL_OLD_ + "SUMRULE="), USER_EPS_SUM_OPTION.option, USER_EPS_SUM_OPTION.xscheme); //CO 170601
      USER_EPS_SUM = USER_EPS_SUM_OPTION.content_double; //CO 170601
      logger << (USER_EPS_SUM_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "SUMRULE=" << USER_EPS_SUM << "." << apl::endl;
      logger << "Convergence criterion for the sumrules is set to " <<  abs(USER_EPS_SUM) << " eV/Angs.^3." << apl::endl;

      // ME 180821 - Get the mixing coefficient for the SCF procedure in AAPL
      USER_AAPL_MIX_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "MIXING_COEFFICIENT=" + "|" + _ASTROPT_APL_OLD_ + "MIXING_COEFFICIENT="), USER_AAPL_MIX_OPTION.option, USER_AAPL_MIX_OPTION.xscheme);
      USER_AAPL_MIX = USER_AAPL_MIX_OPTION.content_double;
      logger << (USER_AAPL_MIX_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "MIXING_COEFFICIENT=" << USER_AAPL_MIX << "." << apl::endl;
      logger << "The SCF for AAPL will use a mixing coefficient of " << USER_AAPL_MIX << "." << apl::endl;

      // ME 180622 - Get the user's numnber of iterations for the sumrule
      USER_AAPL_MAX_ITER_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "SUMRULE_MAX_ITER=" + "|" + _ASTROPT_APL_OLD_ + "SUMRULE_MAX_ITER="), USER_AAPL_MAX_ITER_OPTION.option, USER_AAPL_MAX_ITER_OPTION.xscheme);
      USER_AAPL_MAX_ITER = USER_AAPL_MAX_ITER_OPTION.content_double;
      logger << (USER_AAPL_MAX_ITER_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "SUMRULE_MAX_ITER=" << USER_AAPL_MAX_ITER << "." << apl::endl;
      logger << "Anharmonic IFCs need to be converged within " << abs(USER_AAPL_MAX_ITER) << " iterations." << apl::endl;

      // ME 180913 - Calculate fourth order correction for thermal conductivity
      USER_AAPL_FOURTH_ORDER_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "FOURTH_ORDER=" + "|" +_ASTROPT_APL_OLD_ + "FOURTH_ORDER="), USER_AAPL_FOURTH_ORDER_OPTION.option, USER_AAPL_FOURTH_ORDER_OPTION.xscheme);
      logger << (USER_AAPL_FOURTH_ORDER_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "FOURTH_ORDER=" << (USER_AAPL_FOURTH_ORDER_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      logger << "Thermal conductivity will be calculated " << (USER_AAPL_FOURTH_ORDER_OPTION.option ? "with" : "without") << " four-phonon processes." << apl::endl;

/* OBSOLETE ME 181018
      // Get the users magnitude of the distortion vector (in Angs. and real space)
      USER_TDISTORTION_MAGNITUDE_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "TDMAG=" + "|" +_ASTROPT_APL_OLD_ + "TDMAG=" + "|" + _ASTROPT_ + "TDISMAG=" + "|" +_ASTROPT_APL_OLD_ + "TDISMAG="), USER_TDISTORTION_MAGNITUDE_OPTION.option, USER_TDISTORTION_MAGNITUDE_OPTION.xscheme);  //CO 170621 - TDISMAG legacy
      USER_TDISTORTION_MAGNITUDE = USER_TDISTORTION_MAGNITUDE_OPTION.content_double;
      logger << (USER_TDISTORTION_MAGNITUDE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "TDMAG=" << USER_TDISTORTION_MAGNITUDE << "." << apl::endl;
      logger << "The distortion magnitude for anharmonic IFCs will be " << USER_TDISTORTION_MAGNITUDE << " Angs." << apl::endl;
*/

      // ISOTOPE, e.g. ISOTOPE = y
      CALCULATE_ISOTOPE_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "ISOTOPE=" + "|" +_ASTROPT_APL_OLD_ + "ISOTOPE="), CALCULATE_ISOTOPE_OPTION.option, CALCULATE_ISOTOPE_OPTION.xscheme); //CO 170601
      logger << (CALCULATE_ISOTOPE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "ISOTOPE=" << (CALCULATE_ISOTOPE_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      logger << "Isotope effects " << (CALCULATE_ISOTOPE_OPTION.option ? "will" : "will NOT") << " be considered in the calculation." << apl::endl;

      //scattering at boundaries is handled either with BOUNDARY_OPTION + NANO_SIZE or CUMULATIVEK, if handled AT ALL
      // BOUNDARY, e.g. BOUNDARY = y
      CALCULATE_BOUNDARY_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "BOUNDARY=" + "|" +_ASTROPT_APL_OLD_ + "BOUNDARY="), CALCULATE_BOUNDARY_OPTION.option, CALCULATE_BOUNDARY_OPTION.xscheme); //CO 170601
      logger << (CALCULATE_BOUNDARY_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "BOUNDARY=" << (CALCULATE_BOUNDARY_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      logger << "The boundary term will " << (CALCULATE_BOUNDARY_OPTION.option ? "": "NOT ") << "be considered in the scattering time." << apl::endl;
      // CUMULATIVEK, e.g. CUMULATIVEK = y
      CALCULATE_CUMULATIVEK_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "CUMULATIVEK=" + "|" +_ASTROPT_APL_OLD_ + "CUMULATIVEK="), CALCULATE_CUMULATIVEK_OPTION.option, CALCULATE_CUMULATIVEK_OPTION.xscheme); //CO 170601
      if(!(CALCULATE_BOUNDARY_OPTION.option || CALCULATE_CUMULATIVEK_OPTION.option) ){
        logger << (CALCULATE_CUMULATIVEK_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "CUMULATIVEK=" << (CALCULATE_CUMULATIVEK_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
        logger << "Boundary effects will NOT be considered in calculating the lattice thermal conductivity." << apl::endl;
      } else {
        // BOUNDARY takes precedence
        override_option = false;
        if (CALCULATE_BOUNDARY_OPTION.option && CALCULATE_CUMULATIVEK_OPTION.option) {
          CALCULATE_CUMULATIVEK_OPTION.option=false;
          override_option = true;
          logger << apl::warning << "Both boundary effects and cumulative thermal conductivity cannot be ";
          logger << "set at the same time. Cumulative thermal conductivity has been switched off." << apl::endl;
        }
        logger << (CALCULATE_CUMULATIVEK_OPTION.isentry || override_option ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "CUMULATIVEK=" << (CALCULATE_CUMULATIVEK_OPTION.option ? "ON" : "OFF") << (override_option ? " (overridden by BOUNDARY)": "") << "." << apl::endl;
        logger << "The cumulative lattice thermal conductivity will " << (CALCULATE_CUMULATIVEK_OPTION.option ? "" : "NOT") << " be calculated." << apl::endl;
      }
      // Grain size
      if(CALCULATE_BOUNDARY_OPTION.option){
        USER_NANO_SIZE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "NANO_SIZE=" + "|" + _ASTROPT_APL_OLD_ + "NANO_SIZE="), USER_NANO_SIZE_OPTION.option, USER_NANO_SIZE_OPTION.xscheme); //CO 170601
        USER_NANO_SIZE = USER_NANO_SIZE_OPTION.content_double; //CO 170601
        logger << (USER_NANO_SIZE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "NANO_SIZE=" << USER_NANO_SIZE << "." << apl::endl;
        logger << "Boundary scattering will be computed for grain sizes of " << abs(USER_NANO_SIZE) << " nm." << apl::endl;
      }

      // BTE, e.g. BTE = RTA
      USER_BTE_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "BTE=" + "|" +_ASTROPT_APL_OLD_ + "BTE="), USER_BTE_OPTION.option, USER_BTE_OPTION.xscheme); //CO 170601
      USER_BTE = USER_BTE_OPTION.content_string;
      transform(USER_BTE.begin(), USER_BTE.end(), USER_BTE.begin(), toupper);
      logger << (USER_BTE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "BTE=" << USER_BTE << "." << apl::endl;
      if (USER_BTE == string("RTA")){logger << "The Boltzmann Transport Equation will be solved using the Relaxation Time Approximation Approximation (RTA)." << apl::endl;}
      else if (USER_BTE == string("FULL")) {logger << "The Boltzmann Transport Equation will be solved using an iterative scheme (FULL)." << apl::endl;}
      else {
        string function = "apl::RunPhonons_APL";
        stringstream message;
        message << "Wrong setting in "+_ASTROPT_+"BTE. Specify as BTE=RTA or FULL.";
        throw aurostd::xerror(function, message, _INPUT_ILLEGAL_);
      }

      // THERMALGRID, e.g., THERMALGRID = 2x2x2
      USER_THERMALGRID_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "THERMALGRID=" + "|" + _ASTROPT_APL_OLD_ + "THERMALGRID="), USER_THERMALGRID_OPTION.option, USER_THERMALGRID_OPTION.xscheme); //CO 170601
      USER_THERMALGRID = USER_THERMALGRID_OPTION.content_string; //CO 170601
      tokens.clear();
      apl::tokenize(USER_THERMALGRID, tokens, string(" xX"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in "+_ASTROPT_+"THERMALGRID. Specify as THERMALGRID=10x10x10.");}
      logger << (USER_THERMALGRID_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "THERMALGRID=" << USER_THERMALGRID << "." << apl::endl;

      // TCT, e.g., TCT = 1000:2000:10 -> temperature from 1000 to 2000 K by step 10 K
      // no default here
      USER_TCT_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "TCT=" + "|" + _ASTROPT_APL_OLD_ + "TCT="), USER_TCT_OPTION.option, USER_TCT_OPTION.xscheme); //CO 170601
      tokens.clear();
      apl::tokenize(USER_TCT_OPTION.content_string, tokens, string(" :"));
      if (tokens.size() != 3) {
        string function = "apl::RunPhonons_APL";
        stringstream message;
        message << "Wrong setting in "+_ASTROPT_+"TCT. Specify as TCT=1000:2000:10.";
        throw aurostd::xerror(function, message, _INPUT_NUMBER_);
      }
      USER_TCT_TSTART = aurostd::string2utype<double>(tokens.at(0));
      USER_TCT_TEND = aurostd::string2utype<double>(tokens.at(1));
      USER_TCT_TSTEP = aurostd::string2utype<double>(tokens.at(2));
      override_option = false;
      if (USER_TCT_TSTART == 0) {
        logger << apl::warning << "Thermal conductivity is infinite at 0 K and will be skipped." << apl::endl;
        USER_TCT_TSTART += USER_TCT_TSTEP;
        override_option = true;
      }
      logger << (USER_TCT_OPTION.isentry || override_option ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "TCT=" << USER_TCT_OPTION.content_string << "." << apl::endl;
      logger << "Thermal conductivity will be calculated in the temperature range <" << USER_TCT_TSTART << "," << USER_TCT_TEND << "> with step size " << USER_TCT_TSTEP << " K." << apl::endl;

      // Cutoff Radius (in Angstrom)
      USER_CUTOFF_DISTANCE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "CUT_RAD=" + "|" + _ASTROPT_APL_OLD_ + "CUT_RAD="), USER_CUTOFF_DISTANCE_OPTION.option, USER_CUTOFF_DISTANCE_OPTION.xscheme); //CO 170601
      if (USER_CUTOFF_SHELL_OPTION.isentry) {
        if (!USER_AAPL_FOURTH_ORDER_OPTION.option) {
          USER_CUTOFF_DISTANCE.resize(1);
        }
        tokens.clear();
        apl::tokenize(USER_CUTOFF_DISTANCE_OPTION.content_string, tokens, string(" ,"));
        if (tokens.size() < 1) {
          string function = "apl::RunPhonons_APL";
          stringstream message;
          message << "Not enough entries in "+_ASTROPT_+"CUT_RAD.";
          throw aurostd::xerror(function, message, _INPUT_NUMBER_);
        } else if (tokens.size() > USER_CUTOFF_DISTANCE.size()) {
          logger << apl::warning << "Too many entries for " << _ASTROPT_ << "CUT_RAD. ";
          logger << "Excess entries will be ignored." << apl::endl;
        }
        USER_CUTOFF_DISTANCE[0] = aurostd::string2utype<double>(tokens.at(0));
        logger << "The cutoff to compute the 3rd order anharmonic IFCs will be ";
        logger << abs(USER_CUTOFF_DISTANCE[0]) << " Angstrom." << apl::endl;
        if (USER_AAPL_FOURTH_ORDER_OPTION.option) {
          if (tokens.size() == 1) {
            logger << apl::warning << "Only one entry found for the cutoff radius. ";
            logger << "3rd and 4th order anharmonic IFCs will use the same value." << apl::endl;
            USER_CUTOFF_DISTANCE[1] = USER_CUTOFF_DISTANCE[0];
          } else {
            USER_CUTOFF_DISTANCE[1] = aurostd::string2utype<double>(tokens.at(1));
            logger << "The cutoff to compute the 4th order anharmonic IFCs will be ";
            logger << abs(USER_CUTOFF_DISTANCE[1]) << " Angstrom." << apl::endl;
          }
        }
    }

      // Cutoff radius (in coordination shells)
      USER_CUTOFF_SHELL_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "CUT_SHELL=" + "|" + _ASTROPT_APL_OLD_ + "CUT_SHELL="), USER_CUTOFF_SHELL_OPTION.option, USER_CUTOFF_SHELL_OPTION.xscheme); //CO 170601
      if (USER_CUTOFF_SHELL_OPTION.isentry) {
        if (!USER_AAPL_FOURTH_ORDER_OPTION.option) {
          USER_CUTOFF_SHELL.resize(1);
        }
        tokens.clear();
        apl::tokenize(USER_CUTOFF_SHELL_OPTION.content_string, tokens, string(" ,"));
        if (tokens.size() < 1) {
          string function = "apl::RunPhonons_APL";
          stringstream message;
          message << "Not enough entries in "+_ASTROPT_+"CUT_SHELL.";
          throw aurostd::xerror(function, message, _INPUT_NUMBER_);
        } else if (tokens.size() > USER_CUTOFF_SHELL.size()) {
          logger << apl::warning << "Too many entries for " << _ASTROPT_ << "CUT_SHELL. ";
          logger << "Excess entries will be ignored." << apl::endl;
        }
        USER_CUTOFF_SHELL[0] = aurostd::string2utype<int>(tokens.at(0));
        logger << "The calculation of 3rd order anharmonic IFCs will consider up to ";
        logger << USER_CUTOFF_SHELL[0] << " coordination shells." << apl::endl;
        if (USER_AAPL_FOURTH_ORDER_OPTION.option) {
          if (tokens.size() == 1) {
            logger << apl::warning << "Only one entry found for the number of coordination shells. ";
            logger << "3rd and 4th order anharmonic IFCs will use the same value." << apl::endl;
            USER_CUTOFF_SHELL[1] = USER_CUTOFF_SHELL[0];
          } else {
            USER_CUTOFF_SHELL[1] = aurostd::string2utype<int>(tokens.at(1));
            logger << "The calculation of 4th order anharmonic IFCs will consider up to ";
            logger << USER_CUTOFF_SHELL[1] << " coordination shells." << apl::endl;
          }
        }
      }

      // ME 180501 - If the user only specifies CUT_SHELL or CUT_RAD, unset the default values
      if (USER_CUTOFF_SHELL_OPTION.isentry && ! USER_CUTOFF_DISTANCE_OPTION.isentry){
        if (USER_AAPL_FOURTH_ORDER_OPTION.option) {
          USER_CUTOFF_DISTANCE.assign(2, 0.0);
        } else {
          USER_CUTOFF_DISTANCE.assign(1, 0.0);
        }
      }
      if (! USER_CUTOFF_SHELL_OPTION.isentry && USER_CUTOFF_DISTANCE_OPTION.isentry) {
        if (USER_AAPL_FOURTH_ORDER_OPTION.option) {
          USER_CUTOFF_SHELL.assign(2, 0);
        } else {
          USER_CUTOFF_SHELL.assign(1, 0);
        }
      }
    }
    //Anharmonic forces and thermal conductivity
    //END JJPR ANHARMONIC

    //PINKU PHONON START
    //GROUPVELOCITY=y/n
    CALCULATE_GROUPVELOCITY_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GROUP_VELOCITY=" + "|" + _ASTROPT_APL_OLD_ + "GROUP_VELOCITY="), CALCULATE_GROUPVELOCITY_OPTION.option,  CALCULATE_GROUPVELOCITY_OPTION.xscheme);
    logger << (CALCULATE_GROUPVELOCITY_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "GROUP_VELOCITY=" << (CALCULATE_GROUPVELOCITY_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

    //ATOMIC_DISPLACEMENT=y/n
    CALCULATE_ATOMIC_DISPLACEMENT_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "ATOMIC_DISPLACEMENT=" + "|" + _ASTROPT_APL_OLD_ + "ATOMIC_DISPLACEMENT="), CALCULATE_ATOMIC_DISPLACEMENT_OPTION.option,  CALCULATE_ATOMIC_DISPLACEMENT_OPTION.xscheme);
    logger << (CALCULATE_ATOMIC_DISPLACEMENT_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "ATOMIC_DISPLACEMENT=" << (CALCULATE_ATOMIC_DISPLACEMENT_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    //PINKU PHONON END

    //PINKU QUASI-HARMONIC START
    if(!CALCULATE_TCOND_OPTION.option){
      if(kflags.KBIN_PHONONS_CALCULATION_QHA){
        CALCULATE_GRUNEISEN_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_QHA; //recycle what we parsed earlier
        logger << _ASTROPT_ << "CALC is" << ( CALCULATE_GRUNEISEN_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;
        logger << "The Gruneisen parameter will" << ( CALCULATE_GRUNEISEN_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-QHA." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA_A){
        CALCULATE_GRUNEISEN_A_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_QHA_A;
        logger << _ASTROPT_ << "CALC is" << ( CALCULATE_GRUNEISEN_A_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;
        logger << "The Gruneisen_A parameter will" << ( CALCULATE_GRUNEISEN_A_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-QHA." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA_B){
        CALCULATE_GRUNEISEN_B_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_QHA_B;
        logger << _ASTROPT_ << "CALC is" << ( CALCULATE_GRUNEISEN_B_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;
        logger << "The Gruneisen_B parameter will" << ( CALCULATE_GRUNEISEN_B_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-QHA." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA_C){
        CALCULATE_GRUNEISEN_C_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_QHA_C;
        logger << _ASTROPT_ << "CALC is" << ( CALCULATE_GRUNEISEN_C_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;
        logger << "The Gruneisen_C parameter will" << ( CALCULATE_GRUNEISEN_C_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-QHA." << apl::endl;
      }

      //QHA, QHA3P and SCQHA INCLUDE ELECTRONIC OPTION
      if(kflags.KBIN_PHONONS_CALCULATION_QHA || kflags.KBIN_PHONONS_CALCULATION_QHA_A || kflags.KBIN_PHONONS_CALCULATION_QHA_B || kflags.KBIN_PHONONS_CALCULATION_QHA_C ||
          kflags.KBIN_PHONONS_CALCULATION_SCQHA|| kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C){
        INCLUDE_ELE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "INCLUDE_ELE=" + "|" + _ASTROPT_ + "INCLUDE_ELE="), INCLUDE_ELE_OPTION.option, INCLUDE_ELE_OPTION.xscheme);
        logger << (INCLUDE_ELE_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "INCLUDE_ELE=" << (INCLUDE_ELE_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      }

      //QHA3P and SCQHA temperature dependent phonon dispersion option
      if(kflags.KBIN_PHONONS_CALCULATION_SCQHA|| kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C){
        SCQHA_PDIS_T_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "SCQHA_PDIS_T=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_PDIS_T="), SCQHA_PDIS_T_OPTION.option, SCQHA_PDIS_T_OPTION.xscheme);
        tokens.clear(); scqha_pdis_T.clear();
        apl::tokenize(SCQHA_PDIS_T_OPTION.content_string, tokens, string(" ,"));
        if (tokens.size() == 0) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"SCQHA_PDIS_T. Specify as SCQHA_PDIS_T=-100.0, 300.0, 600.0");}
        if(tokens.size()!=0){
          for (uint i=0; i<tokens.size(); i++){
            scqha_pdis_T.push_back(aurostd::string2utype<double>(tokens.at(i)));
          }
        }
      }
      //rectifying possible user QHA-input errors
      if(CALCULATE_GRUNEISEN_OPTION.option){
        CALCULATE_SCQHA_A_OPTION.option=false;
        CALCULATE_SCQHA_B_OPTION.option=false;
        CALCULATE_SCQHA_C_OPTION.option=false;
      } else if(CALCULATE_GRUNEISEN_A_OPTION.option){
        CALCULATE_SCQHA_OPTION.option=false;
        CALCULATE_SCQHA_B_OPTION.option=false;
        CALCULATE_SCQHA_C_OPTION.option=false;
      } else if(CALCULATE_GRUNEISEN_B_OPTION.option){
        CALCULATE_SCQHA_OPTION.option=false;
        CALCULATE_SCQHA_A_OPTION.option=false;
        CALCULATE_SCQHA_C_OPTION.option=false;
      } else if(CALCULATE_GRUNEISEN_C_OPTION.option){
        CALCULATE_SCQHA_OPTION.option=false;
        CALCULATE_SCQHA_A_OPTION.option=false;
        CALCULATE_SCQHA_B_OPTION.option=false;
      }

      //Writing to log
      if(kflags.KBIN_PHONONS_CALCULATION_SCQHA){
        logger << (CALCULATE_SCQHA_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA=" << (CALCULATE_SCQHA_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_SCQHA_A){
        logger << (CALCULATE_SCQHA_A_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_A=" << (CALCULATE_SCQHA_A_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_SCQHA_B){
        logger << (CALCULATE_SCQHA_B_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_B=" << (CALCULATE_SCQHA_B_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_SCQHA_C){
        logger << (CALCULATE_SCQHA_C_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_C=" << (CALCULATE_SCQHA_C_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      }

      //Writing to log
      if(kflags.KBIN_PHONONS_CALCULATION_QHA3P){
        logger << (CALCULATE_QHA3P_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "QHA3P=" << (CALCULATE_QHA3P_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA3P_A){
        logger << (CALCULATE_QHA3P_A_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "QHA3P_A=" << (CALCULATE_QHA3P_A_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA3P_B){
        logger << (CALCULATE_QHA3P_B_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "QHA3P_B=" << (CALCULATE_QHA3P_B_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      } else if(kflags.KBIN_PHONONS_CALCULATION_QHA3P_C){
        logger << (CALCULATE_QHA3P_C_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "QHA3P_C=" << (CALCULATE_QHA3P_C_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      }

      if(CALCULATE_SCQHA_OPTION.option || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C){
        SCQHA_DISTORTION_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "SCQHA_DISTORTION=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_DISTORTION="), SCQHA_DISTORTION_OPTION.option, SCQHA_DISTORTION_OPTION.xscheme);
        SCQHA_DISTORTION=SCQHA_DISTORTION_OPTION.content_double;
        if (SCQHA_DISTORTION_OPTION.isentry) {
          tokens.clear();
          apl::tokenize(SCQHA_DISTORTION_OPTION.content_string, tokens, string(" "));
          if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"SCQHA_DISTORTION. Specify as SCQHA_DISTORTION_OPTION=3.0.");}
        }
        logger << (SCQHA_DISTORTION_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "SCQHA_DISTORTION=" << SCQHA_DISTORTION << "." << apl::endl;
      }

      if(CALCULATE_GRUNEISEN_OPTION.option || kflags.KBIN_PHONONS_CALCULATION_QHA_A || kflags.KBIN_PHONONS_CALCULATION_QHA_B || kflags.KBIN_PHONONS_CALCULATION_QHA_C || kflags.KBIN_PHONONS_CALCULATION_SCQHA ||
          kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C || kflags.KBIN_PHONONS_CALCULATION_QHA3P ||
          kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C){

        CUTOFF_FREQ_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "CUTOFF_FREQ=" + "|" + _ASTROPT_APL_OLD_ + "CUTOFF_FREQ="), CUTOFF_FREQ_OPTION.option, CUTOFF_FREQ_OPTION.xscheme); //CO 170601
        CUTOFF_FREQ = CUTOFF_FREQ_OPTION.content_double;
        if (CUTOFF_FREQ_OPTION.isentry) {
          tokens.clear();
          apl::tokenize(CUTOFF_FREQ_OPTION.content_string, tokens, string(" "));
          if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"CUTOFF_FREQ. Specify as CUTOFF_FREQ=0.01.");}
        }
        logger << (CUTOFF_FREQ_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "CUTOFF_FREQ=" << CUTOFF_FREQ_OPTION.content_string << "." << apl::endl;


        CALCULATE_EOS_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "EOS=" + "|" + _ASTROPT_APL_OLD_ + "EOS="), CALCULATE_EOS_OPTION.option, CALCULATE_EOS_OPTION.xscheme); //CO 170601
        logger << (CALCULATE_EOS_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "EOS=" << (CALCULATE_EOS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

        if ((CALCULATE_EOS_OPTION.option) && (CALCULATE_GRUNEISEN_OPTION.option || kflags.KBIN_PHONONS_CALCULATION_QHA_A || kflags.KBIN_PHONONS_CALCULATION_QHA_B ||
              kflags.KBIN_PHONONS_CALCULATION_QHA_C || kflags.KBIN_PHONONS_CALCULATION_SCQHA || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A ||
              kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C || kflags.KBIN_PHONONS_CALCULATION_QHA3P ||
              kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C))

        {
          EOS_DISTORTION_RANGE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "EOS_DISTORTION_RANGE=" + "|" + _ASTROPT_APL_OLD_ + "EOS_DISTORTION_RANGE="), EOS_DISTORTION_RANGE_OPTION.option, EOS_DISTORTION_RANGE_OPTION.xscheme);
          tokens.clear();
          apl::tokenize(EOS_DISTORTION_RANGE_OPTION.content_string, tokens, string(" :"));
          if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"EOS_DISTORTION_RANGE. Specify as EOS_DISTORTION_RANGE=-3:6:1.");}
          EOS_DISTORTION_START = aurostd::string2utype<double>(tokens.at(0));
          EOS_DISTORTION_END = aurostd::string2utype<double>(tokens.at(1));
          EOS_DISTORTION_DISTORTION_INC = aurostd::string2utype<double>(tokens.at(2));
          logger << (EOS_DISTORTION_RANGE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "EOS_DISTORTION_RANGE=" << EOS_DISTORTION_RANGE_OPTION.content_string << "." << apl::endl;
          logger << "The EOS properties will be calculated in distortion range <" << EOS_DISTORTION_START << "," << EOS_DISTORTION_END << "," << EOS_DISTORTION_DISTORTION_INC << "." << apl::endl;

          FITTING_TYPE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "FITTING_TYPE=" + "|" + _ASTROPT_APL_OLD_ + "FITTING_TYPE="), FITTING_TYPE_OPTION.option, FITTING_TYPE_OPTION.xscheme); //CO 170601
          FITTING_TYPE = FITTING_TYPE_OPTION.content_string;
          if (FITTING_TYPE_OPTION.isentry) {
            tokens.clear();
            apl::tokenize(FITTING_TYPE_OPTION.content_string, tokens, string(" "));
            if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"FITTING_TYPE. Specify as FITTING_TYPE=BM2.");}
          }
          logger << (FITTING_TYPE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "FITTING_TYPE=" << FITTING_TYPE_OPTION.content_string << "." << apl::endl;
          logger << "EOS fitting type found = " << FITTING_TYPE << "." << apl::endl;

          EOS_STATIC_KPPRA_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "EOS_STATIC_KPPRA=" + "|" + _ASTROPT_APL_OLD_ + "EOS_STATIC_KPPRA="), EOS_STATIC_KPPRA_OPTION.option, EOS_STATIC_KPPRA_OPTION.xscheme);
          EOS_STATIC_KPPRA = EOS_STATIC_KPPRA_OPTION.content_int;
          if(CALCULATE_GRUNEISEN_OPTION.option){
            if (EOS_STATIC_KPPRA_OPTION.isentry) {
              tokens.clear();
              apl::tokenize(EOS_STATIC_KPPRA_OPTION.content_string, tokens, string(" "));
              if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"EOS_STATIC_KPPRA. Specify as EOS_STATIC_KPPRA=10000.");}
            }
            logger << (EOS_STATIC_KPPRA_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "EOS_STATIC_KPPRA=" << EOS_STATIC_KPPRA_OPTION.content_string << "." << apl::endl;
          }
          NEDOS_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "NEDOS=" + "|" + _ASTROPT_APL_OLD_ + "NEDOS="), NEDOS_OPTION.option, NEDOS_OPTION.xscheme);
          NEDOS = NEDOS_OPTION.content_int;
          if (NEDOS_OPTION.isentry) {
            tokens.clear();
            apl::tokenize(NEDOS_OPTION.content_string, tokens, string(" "));
            if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"NEDOS. Specify as NEDOS=5000.");}
          }
          logger << (NEDOS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "NEDOS=" << NEDOS_OPTION.content_string << "." << apl::endl;
        }
      }
    }
    //GP SUBDIRECTORY OPTIONS. These are automic options and not controlled by users
    CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GRUNEISEN_SD=" + "|" + _ASTROPT_APL_OLD_ + "GRUNEISEN_SD="), CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option, CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.xscheme);
    if(!CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option){
      CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GRUNEISEN_A_SD=" + "|" + _ASTROPT_APL_OLD_ + "GRUNEISEN_A_SD="), CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option, CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.xscheme);
      if(!CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option){
        CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GRUNEISEN_B_SD=" + "|" + _ASTROPT_APL_OLD_ + "GRUNEISEN_B_SD="), CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option, CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.xscheme);
        if(!CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option){
          CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "GRUNEISEN_C_SD=" + "|" + _ASTROPT_APL_OLD_ + "GRUNEISEN_C_SD="), CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option, CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.xscheme);
        }
      }
    }

    //Writing to log
    if(CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GRUNEISEN_SD=" << (CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GRUNEISEN_A_SD=" << (CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GRUNEISEN_B_SD=" << (CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GRUNEISEN_C_SD=" << (CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    }

    if(CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option|| CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option || CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option || CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option){
      GP_DISTORTION_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "GP_DISTORTION=" + "|" + _ASTROPT_APL_OLD_ + "GP_DISTORTION="), GP_DISTORTION_OPTION.option, GP_DISTORTION_OPTION.xscheme);
      GP_DISTORTION=GP_DISTORTION_OPTION.content_double;
      if (GP_DISTORTION_OPTION.isentry) {
        tokens.clear();
        apl::tokenize(GP_DISTORTION_OPTION.content_string, tokens, string(" "));
        if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"GP_DISTORTION. Specify as GP_DISTORTION=0.03.");}
      }
      logger << (GP_DISTORTION_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "GP_DISTORTION=" << GP_DISTORTION << "." << apl::endl;
    }
    //EOS SUBDIRECTORY OPTIONS. These are automic options and not controlled by users
    CALCULATE_EOS_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "EOS_SD=" + "|" + _ASTROPT_APL_OLD_ + "EOS_SD="), CALCULATE_EOS_SUBDIRECTORIES_OPTION.option, CALCULATE_EOS_SUBDIRECTORIES_OPTION.xscheme); //CO 170601
    if(CALCULATE_EOS_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_EOS_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "EOS_SD=" << (CALCULATE_EOS_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      EOS_DISTORTION_RANGE_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "EOS_DISTORTION_RANGE=" + "|" + _ASTROPT_APL_OLD_ + "EOS_DISTORTION_RANGE="), EOS_DISTORTION_RANGE_OPTION.option, EOS_DISTORTION_RANGE_OPTION.xscheme); //CO 170601
      tokens.clear();
      apl::tokenize(EOS_DISTORTION_RANGE_OPTION.content_string, tokens, string(" :"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"EOS_DISTORTION_RANGE. Specify as EOS_DISTORTION_RANGE=-3:6:1");}
      EOS_DISTORTION_START = aurostd::string2utype<double>(tokens.at(0));
      EOS_DISTORTION_END = aurostd::string2utype<double>(tokens.at(1));
      EOS_DISTORTION_DISTORTION_INC = aurostd::string2utype<double>(tokens.at(2));
      logger << (EOS_DISTORTION_RANGE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "EOS_DISTORTION_RANGE=" << EOS_DISTORTION_RANGE_OPTION.content_string << "." << apl::endl;
      logger << "The EOS properties will be calculated in distortion range <" << EOS_DISTORTION_START << "," << EOS_DISTORTION_END << "," << EOS_DISTORTION_DISTORTION_INC << "." << apl::endl;
    }
    //}

    //SCQHA SUBDIRECTORY OPTIONS. These are automic options and not controlled by users
    CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "SCQHA_SD=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_SD="), CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option,  CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.xscheme);
    if(!CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option){
      CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "SCQHA_A_SD=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_A_SD="), CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option,  CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.xscheme);
      if(!CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option){
        CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "SCQHA_B_SD=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_B_SD="), CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option,  CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.xscheme);
        if(!CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option){
          CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.options2entry(AflowIn, string(_ASTROPT_ + "SCQHA_C_SD=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_C_SD="), CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.option,  CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.xscheme);
        }
      }
    }

    //Writing to log
    if(CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_SD=" << (CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_A_SD=" << (CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_B_SD=" << (CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    } else if(CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.option){
      logger << (CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.isentry ? "Setting " : "DEFAULT ") << _ASTROPT_ << "SCQHA_C_SD=" << (CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    }

    if(CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option || CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option || CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option || CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.isentry){
      SCQHA_DISTORTION_OPTION.options2entry(AflowIn, string( _ASTROPT_ + "SCQHA_DISTORTION=" + "|" + _ASTROPT_APL_OLD_ + "SCQHA_DISTORTION="), SCQHA_DISTORTION_OPTION.option, SCQHA_DISTORTION_OPTION.xscheme);
      SCQHA_DISTORTION=SCQHA_DISTORTION_OPTION.content_double;
      if (SCQHA_DISTORTION_OPTION.isentry) {
        tokens.clear();
        apl::tokenize(SCQHA_DISTORTION_OPTION.content_string, tokens, string(" "));
        if (tokens.size() != 1){throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"SCQHA_DISTORTION. Specify as SCQHA_DISTORTION=3.0.");}
      }
      logger << (SCQHA_DISTORTION_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "SCQHA_DISTORTION=" << SCQHA_DISTORTION << "." << apl::endl;
    }
    // PINKU QUASI-HARMONIC END

    if( !(CALCULATE_TCOND_OPTION.option||CALCULATE_GRUNEISEN_OPTION.option) ){
      CALCULATE_APL_OPTION.option = kflags.KBIN_PHONONS_CALCULATION_APL;
      logger << _ASTROPT_ << "CALC is" << ( CALCULATE_APL_OPTION.option ? "" : " NOT" ) << " set." << apl::endl;  //CO 170601
      logger << "Harmonic force constants will" << ( CALCULATE_APL_OPTION.option ? "" : " NOT" ) << " be computed via AFLOW-APL." << apl::endl;   //CO 170601
    }  //CO 170601

    //BELOW HERE, all tags all shared, so we need to look for all combinations

    // HIBERNATE, e.g. HIBERNATE = y
    USER_HIBERNATE_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "HIBERNATE=" + "|" + _ASTROPT_QHA_ + "HIBERNATE=" + "|" + _ASTROPT_AAPL_ + "HIBERNATE=" + "|" + _ASTROPT_APL_OLD_ + "HIBERNATE="), USER_HIBERNATE_OPTION.option, USER_HIBERNATE_OPTION.xscheme); //CO 170601
    logger << (USER_HIBERNATE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "HIBERNATE=" << (USER_HIBERNATE_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    logger << "The hibernate feature is switched " << (USER_HIBERNATE_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

    // FREQFORMAT, e.g. FREQFORMAT = "THz | allow_negative"
    //COREY, would help if you had some sort of stupidity test here
    USER_FREQFORMAT_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "FREQFORMAT=" + "|" + _ASTROPT_QHA_ + "FREQFORMAT=" + "|" + _ASTROPT_AAPL_ + "FREQFORMAT=" + "|" + _ASTROPT_APL_OLD_ + "FREQFORMAT="), USER_FREQFORMAT_OPTION.option, USER_FREQFORMAT_OPTION.xscheme); //CO 170601
    USER_FREQFORMAT = USER_FREQFORMAT_OPTION.content_string;  //CO 170601
    transform(USER_FREQFORMAT.begin(), USER_FREQFORMAT.end(), USER_FREQFORMAT.begin(), toupper);
    logger << (USER_FREQFORMAT_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "FREQFORMAT=" << USER_FREQFORMAT << "." << apl::endl;
    logger << "The frequency will be returned in this format: " << USER_FREQFORMAT << "." << apl::endl;

    // ENGINE, e.g., ENGINE = DM or LR
    USER_ENGINE_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "ENGINE=" + "|" + _ASTROPT_QHA_ + "ENGINE=" + "|" + _ASTROPT_AAPL_ + "ENGINE=" + "|" + _ASTROPT_APL_OLD_ + "ENGINE="), USER_ENGINE_OPTION.option, USER_ENGINE_OPTION.xscheme); //CO 170601
    USER_ENGINE = USER_ENGINE_OPTION.content_string;
    transform(USER_ENGINE.begin(), USER_ENGINE.end(), USER_ENGINE.begin(), toupper);
    logger << (USER_ENGINE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "ENGINE=" << USER_ENGINE << "." << apl::endl;
    if (USER_ENGINE == string("DM")) {logger << "The phonon calculator engine: Direct Method (DM)." << apl::endl;}
    else if (USER_ENGINE == string("GSA")) {
      //CO generally redirects to DM, the distinction between DM and GSA is obsolete
      logger << "The Generalized Supercell Approach (GSA) phonon calculator now directs to another engine: Direct Method (DM)." << apl::endl;
      USER_ENGINE = "DM";
    } else if (USER_ENGINE == string("LR")) {logger << "The phonon calculator engine: Linear Response (LR)." << apl::endl;}
    else {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"ENGINE. Specify as ENGINE=DM or LR.");}

    // Get the users magnitute of the distortion vector (in Angs. and real space)
    if (USER_ENGINE == string("DM")) {
      USER_DISTORTION_MAGNITUDE_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DMAG=" + "|" + _ASTROPT_QHA_ + "DMAG=" + "|" + _ASTROPT_AAPL_ + "DMAG=" + "|" + _ASTROPT_APL_OLD_ + "DMAG=" + "|" + _ASTROPT_APL_ + "DISMAG=" + "|" + _ASTROPT_QHA_ + "DISMAG=" + "|" + _ASTROPT_AAPL_ + "DISMAG=" + "|" + _ASTROPT_APL_OLD_ + "DISMAG="), USER_DISTORTION_MAGNITUDE_OPTION.option, USER_DISTORTION_MAGNITUDE_OPTION.xscheme); //CO 170601, 170621 DISMAG legacy
      USER_DISTORTION_MAGNITUDE = USER_DISTORTION_MAGNITUDE_OPTION.content_double;
      logger << (USER_DISTORTION_MAGNITUDE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DMAG=" << USER_DISTORTION_MAGNITUDE << "." << apl::endl;
      logger << "The distortion magnitude will be " << USER_DISTORTION_MAGNITUDE << " Angs." << apl::endl;

      // Get flag for generation of displacements with positive/negative magnitude
      USER_DISTORTIONS_PLUS_MINUS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DPM=" + "|" + _ASTROPT_QHA_ + "DPM=" + "|" + _ASTROPT_AAPL_ + "DPM=" + "|" + _ASTROPT_APL_OLD_ + "DPM="), USER_DISTORTIONS_PLUS_MINUS_OPTION.option, USER_DISTORTIONS_PLUS_MINUS_OPTION.xscheme); //CO 170601
      AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option = !USER_DISTORTIONS_PLUS_MINUS_OPTION.option;
      if (USER_DISTORTIONS_PLUS_MINUS_OPTION.isentry){logger << (USER_DISTORTIONS_PLUS_MINUS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DPM=" << (USER_DISTORTIONS_PLUS_MINUS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;}
      if (USER_DISTORTIONS_PLUS_MINUS_OPTION.option){logger << "Each distortion will be generated with both positive and negative magnitudes." << apl::endl;}
      if (USER_DISTORTIONS_PLUS_MINUS_OPTION.isentry && !USER_DISTORTIONS_PLUS_MINUS_OPTION.option){logger << apl::warning << "Distortions will only be considered in one direction (no \"negative\" distortions) - this is NOT recommended." << apl::endl;}
      if (AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option){logger << "DEFAULT " << _ASTROPT_ << "DPM=AUTO which considers negative distortions on a per-site basis." << apl::endl;} //corey AUTO

      // Get flag for generation of displacements only along the x, y, and z axis
      USER_DISTORTIONS_XYZ_ONLY_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DXYZONLY=" + "|" + _ASTROPT_QHA_ + "DXYZONLY=" + "|" + _ASTROPT_AAPL_ + "DXYZONLY=" + "|" + _ASTROPT_APL_OLD_ + "DXYZONLY="), USER_DISTORTIONS_XYZ_ONLY_OPTION.option, USER_DISTORTIONS_XYZ_ONLY_OPTION.xscheme); //CO 170601
      logger << (USER_DISTORTIONS_XYZ_ONLY_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DXYZONLY=" << (USER_DISTORTIONS_XYZ_ONLY_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      if (USER_DISTORTIONS_XYZ_ONLY_OPTION.option){logger << "Only distortions along the x, y, and z axis will be used." << apl::endl;}
      else {logger << "Distortions will be created along lattice vectors including faces and body diagonals." << apl::endl;}

      // ZEROSTATE; One next calculation will be done with no distortion, a such obtained
      // forces will be subtracted from the all forces obtained with distortions
      USER_ZEROSTATE_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "ZEROSTATE=" + "|" + _ASTROPT_QHA_ + "ZEROSTATE=" + "|" + _ASTROPT_AAPL_ + "ZEROSTATE=" + "|" + _ASTROPT_APL_OLD_ + "ZEROSTATE="), USER_ZEROSTATE_OPTION.option, USER_ZEROSTATE_OPTION.xscheme); //CO 170601
      logger << (USER_ZEROSTATE_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "ZEROSTATE=" << (USER_ZEROSTATE_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      if (USER_ZEROSTATE_OPTION.option){logger << "The zero state forces will be also calculated and subtracted from all forces." << apl::endl;}
    }

    // Do polar correction
    if (USER_ENGINE == string("LR") || USER_ENGINE == string("DM")) {
      CALCULATE_POLAR_CORRECTIONS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "POLAR=" + "|" + _ASTROPT_QHA_ + "POLAR=" + "|" + _ASTROPT_AAPL_ + "POLAR=" + "|" + _ASTROPT_APL_OLD_ + "POLAR="), CALCULATE_POLAR_CORRECTIONS_OPTION.option, CALCULATE_POLAR_CORRECTIONS_OPTION.xscheme); //CO 170601
      logger << (CALCULATE_POLAR_CORRECTIONS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "POLAR=" << (CALCULATE_POLAR_CORRECTIONS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      logger << "The calculation of POLAR MATERIALS corrections is switched " << (CALCULATE_POLAR_CORRECTIONS_OPTION.option? "ON": "OFF") << "." << apl::endl;
    }

    //fix vasp bin for LR or DM+POLAR
    if (USER_ENGINE == string("LR") || (USER_ENGINE == string("DM") && CALCULATE_POLAR_CORRECTIONS_OPTION.option)) {
      if(xflags.AFLOW_MODE_VASP){
      try {
        // Check the version of VASP binary
        logger << "Checking VASP version ... ";
        string vaspVersion;
        vaspVersion = apl::getVASPVersionString( (kflags.KBIN_MPI ? kflags.KBIN_MPI_BIN : kflags.KBIN_BIN ) );
        if (!vaspVersion.empty()) {
          logger << "[" << vaspVersion << "]";
          if ((vaspVersion[0] - '0') < 5) { //cool way of getting ascii value:  https://stackoverflow.com/questions/36310181/char-subtraction-in-c
            logger << apl::warning << "." << apl::endl;
            if(_WITHIN_DUKE_){
              kflags.KBIN_BIN = DEFAULT_VASP5_BIN;
              kflags.KBIN_MPI_BIN = DEFAULT_VASP5_MPI_BIN;
              logger << apl::warning << "Modifying VASP bin to " << kflags.KBIN_BIN << " (Duke machine AUTO modification)." << apl::endl;
            } else {
              throw apl::APLRuntimeError("The LR engine needs VASP5 or higher version.");
            }
          } else {logger << " OK." << apl::endl;}
        } else {
          logger << "Failed." << apl::warning << apl::endl; 
          throw apl::APLLogicError("Unexpected binary format.");
        }
      } catch (apl::APLLogicError& e) {
        logger << apl::warning << "Failed to identify the version of VASP binary." << apl::endl;
        logger << apl::warning << e.what() << apl::endl;
      }
    }
    }

    // SUPERCELL ---------------------------------------------------------

    // RELAX, Should be the primitive structure relax before any supercell is build?
    //CO looks like this has yet to be implemented, need to figure this out
    USER_WANTS_RELAX_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "RELAX=" + "|" + _ASTROPT_QHA_ + "RELAX=" + "|" + _ASTROPT_AAPL_ + "RELAX=" + "|" + _ASTROPT_APL_OLD_ + "RELAX="), USER_WANTS_RELAX_OPTION.option, USER_WANTS_RELAX_OPTION.xscheme); //CO 170601
    if(USER_WANTS_RELAX_OPTION.isentry){
      logger << "Setting " << _ASTROPT_ << "RELAX=" << (USER_WANTS_RELAX_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
      if (USER_WANTS_RELAX_OPTION.option){logger << "The primitive cell is going to relax before the supercell build." << apl::endl;}
    }

    found_supercell = false;
    // SUPERCELL, e.g., SUPERCELL = 2x2x2
    USER_SUPERCELL_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "SUPERCELL=" + "|" + _ASTROPT_QHA_ + "SUPERCELL=" + "|" + _ASTROPT_AAPL_ + "SUPERCELL=" + "|" + _ASTROPT_APL_OLD_ + "SUPERCELL="), USER_SUPERCELL_OPTION.option, USER_SUPERCELL_OPTION.xscheme); //CO 170601
    USER_SUPERCELL = USER_SUPERCELL_OPTION.content_string;
    if(USER_SUPERCELL_OPTION.isentry){
      tokens.clear();
      apl::tokenize(USER_SUPERCELL, tokens, string(" xX"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"SUPERCELL. Specify as SUPERCELL=2x2x2.");}
      // If the MAXSHELL was not specify, clear the default setting for shell restriction
      //if (!aurostd::substring2bool(AflowIn, _ASTROPT_ + "MINSHELL=", TRUE)){USER_MINSHELL = -1;}
      //CO, I think this should be MAX, not MIN
      //if(!USER_MAXSHELL_OPTION.isentry){USER_MAXSHELL = -1;}  //not really important
      logger << (USER_SUPERCELL_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "SUPERCELL=" << USER_SUPERCELL << "." << apl::endl;
      logger << "Supercell will be built with dimensions " << USER_SUPERCELL << "." << apl::endl;
      found_supercell = true;
    }

    // Get the users minimum atoms which will be included into calculation
    // takes back seat to SUPERCELL
    if(!found_supercell){
      USER_MINATOMS_OPTION.options2entry(AflowIn, string( 
            _ASTROPT_APL_ + "MINATOMS=" + "|" + _ASTROPT_QHA_ + "MINATOMS=" + "|" + _ASTROPT_AAPL_ + "MINATOMS=" + "|" + _ASTROPT_APL_OLD_ + "MINATOMS=" + "|" +
            _ASTROPT_APL_ + "MINATOMS_RESTRICTED=" + "|" + _ASTROPT_QHA_ + "MINATOMS_RESTRICTED=" + "|" + _ASTROPT_AAPL_ + "MINATOMS_RESTRICTED=" + "|" + _ASTROPT_APL_OLD_ + "MINATOMS_RESTRICTED="  //CO 180418 - restricted means all dims are equal
            ), USER_MINATOMS_OPTION.option, USER_MINATOMS_OPTION.xscheme); //CO 170601
      USER_MINATOMS_RESTICTED_OPTION=aurostd::substring2bool(USER_MINATOMS_OPTION.keyword,"_RESTRICTED");
      USER_MINATOMS = USER_MINATOMS_OPTION.content_int;
      if(USER_MINATOMS_OPTION.isentry){
        logger << (USER_MINATOMS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "MINATOMS" << (USER_MINATOMS_RESTICTED_OPTION?string("_RESTRICTED"):string("")) << "=" << USER_MINATOMS << "." << apl::endl;
        logger << "Supercell will be built with at least " << USER_MINATOMS << " atoms." << apl::endl;
        found_supercell = true;
      }
    }

    // Get the users maximum shell which will be included into calculation
    // CO, not sure how maxshell works here (F option), need to investigate further and add to README
    // also seems USER_WANTS_FULL_SHELL applies for both MAX and MIN shell settings, should one take precedence? should they be separate flags?
    if(!found_supercell){
      USER_MAXSHELL_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "MAXSHELL=" + "|" + _ASTROPT_QHA_ + "MAXSHELL=" + "|" + _ASTROPT_AAPL_ + "MAXSHELL=" + "|" + _ASTROPT_APL_OLD_ + "MAXSHELL="), USER_MAXSHELL_OPTION.option, USER_MAXSHELL_OPTION.xscheme); //CO 170601
      test = USER_MAXSHELL_OPTION.content_string;
      if (test[test.size() - 1] == 'f' || test[test.size() - 1] == 'F') {
        USER_MAXSHELL = aurostd::string2utype<int>(test.substr(0, test.size() - 1));
        USER_WANTS_FULL_SHELL = true;
      } else {
        USER_MAXSHELL = USER_MAXSHELL_OPTION.content_int;
        USER_WANTS_FULL_SHELL = false;
      }
      if(USER_MAXSHELL_OPTION.isentry){
        logger << (USER_MAXSHELL_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "MAXSHELL=" << USER_MAXSHELL << (USER_WANTS_FULL_SHELL ? " (FULL)" : "") << "." << apl::endl;
        logger << "Supercell will be built with at most " << USER_MAXSHELL << " shells." << apl::endl;
        found_supercell = true;
      }
    }

    // Get the users minimum shell which will be included into calculation
    // CO, not sure how minshell works here (F option), need to investigate further and add to README
    if(!found_supercell){
      USER_MINSHELL_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "MINSHELL=" + "|" + _ASTROPT_QHA_ + "MINSHELL=" + "|" + _ASTROPT_AAPL_ + "MINSHELL=" + "|" + _ASTROPT_APL_OLD_ + "MINSHELL="), USER_MINSHELL_OPTION.option, USER_MINSHELL_OPTION.xscheme); //CO 170601
      test = USER_MINSHELL_OPTION.content_string;
      if (test[test.size() - 1] == 'f' || test[test.size() - 1] == 'F') {
        USER_MINSHELL = aurostd::string2utype<int>(test.substr(0, test.size() - 1));
        USER_WANTS_FULL_SHELL = true;
      } else {
        USER_MINSHELL = USER_MINSHELL_OPTION.content_int;
        USER_WANTS_FULL_SHELL = false;
      }
      if(USER_MINSHELL_OPTION.isentry){
        logger << (USER_MINSHELL_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "MINSHELL=" << USER_MINSHELL << (USER_WANTS_FULL_SHELL ? " (FULL)" : "") << "." << apl::endl;
        logger << "Supercell will be built with at least " << USER_MINSHELL << " shells." << apl::endl;
        found_supercell = true;
      }
    }

    // Get the users KPPRA which will be included into calculation, get from AFLOW machinery
    if(xflags.AFLOW_MODE_VASP){
      if (xflags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.isentry) {logger << "Overriding with " << _ASTROPT_ << "KPPRA=" << xflags.vflags.KBIN_VASP_KPOINTS_PHONONS_KPPRA.content_uint << "." << apl::endl;}
    // Get the users KSCHEME which will be included into calculation, get from AFLOW machinery
      if (xflags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.isentry) {logger << "Overriding with " << _ASTROPT_ << "KSCHEME=" << xflags.vflags.KBIN_VASP_KPOINTS_PHONONS_KSCHEME.content_string << "." << apl::endl;}
      if (xflags.vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.flag("EVEN")) {logger << "Overriding with " << _ASTROPT_ << "KPOINTS=EVEN" << "." << apl::endl;}
      if (xflags.vflags.KBIN_VASP_FORCE_OPTION_KPOINTS_PHONONS_PARITY.flag("ODD")) {logger << "Overriding with " << _ASTROPT_ << "KPOINTS=ODD" << "." << apl::endl;}
    }

    // ADDITIONAL OPTIONS ------------------------------------------------

    // DC. e.g., DC = yes
    CALCULATE_PHONON_DISPERSIONS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DC=" + "|" + _ASTROPT_QHA_ + "DC=" + "|" + _ASTROPT_AAPL_ + "DC=" + "|" + _ASTROPT_APL_OLD_ + "DC="), CALCULATE_PHONON_DISPERSIONS_OPTION.option, CALCULATE_PHONON_DISPERSIONS_OPTION.xscheme); //CO 170601
    logger << (CALCULATE_PHONON_DISPERSIONS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DC=" << (CALCULATE_PHONON_DISPERSIONS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    logger << "The calculation of PHONON DISPERSION curves is switched " << (CALCULATE_PHONON_DISPERSIONS_OPTION.option? "ON" : "OFF" ) << "." << apl::endl;

    if (CALCULATE_PHONON_DISPERSIONS_OPTION.option) {
      //CO 180406 - only choose one way to define path, INITLATTICE/INITSG vs. INITCOORDS_FRAC/INITCOORDS_CART + INITCOORDS_LABELS
      bool found_user_path=false;
      
      // DCINITLATTICE, e.g.,       DCINITLATTICE = RHL
      // DCINITSSG,  e.g., DCINITSG = 166
      // DCINITCOORDSFRAC, e.g.,    DCINITCOORDSFRAC = 0,0,0;0.5,0.5,0.5
      // DCINITCOORDSCART, e.g.,    DCINITCOORDSCART = 0,0,0;0.5,0.5,0.5
      // DCINITCOORDSLABELS, e.g.,  DCINITCOORDSLABELS = G,L,M,X
      // defaults are empty strings, so don't set unless isentry
      USER_DC_INITLATTICE_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCINITLATTICE=" + "|" + _ASTROPT_QHA_ + "DCINITLATTICE=" + "|" + _ASTROPT_AAPL_ + "DCINITLATTICE=" + "|" + _ASTROPT_APL_OLD_ + "DCINITLATTICE="), USER_DC_INITLATTICE_OPTION.option, USER_DC_INITLATTICE_OPTION.xscheme); //CO 170601
      USER_DC_INITSG_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCINITSG=" + "|" + _ASTROPT_QHA_ + "DCINITSG=" + "|" + _ASTROPT_AAPL_ + "DCINITSG=" + "|" + _ASTROPT_APL_OLD_ + "DCINITSG="), USER_DC_INITSG_OPTION.option, USER_DC_INITSG_OPTION.xscheme); //CO 170601
      USER_DC_INITCOORDS_FRAC_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCINITCOORDSFRAC=" + "|" + _ASTROPT_QHA_ + "DCINITCOORDSFRAC=" + "|" + _ASTROPT_AAPL_ + "DCINITCOORDSFRAC=" + "|" + _ASTROPT_APL_OLD_ + "DCINITCOORDSFRAC="), USER_DC_INITCOORDS_FRAC_OPTION.option, USER_DC_INITCOORDS_FRAC_OPTION.xscheme); //CO 170601
      USER_DC_INITCOORDS_CART_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCINITCOORDSCART=" + "|" + _ASTROPT_QHA_ + "DCINITCOORDSCART=" + "|" + _ASTROPT_AAPL_ + "DCINITCOORDSCART=" + "|" + _ASTROPT_APL_OLD_ + "DCINITCOORDSCART="), USER_DC_INITCOORDS_CART_OPTION.option, USER_DC_INITCOORDS_CART_OPTION.xscheme); //CO 170601
      USER_DC_INITCOORDS_LABELS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCINITCOORDSLABELS=" + "|" + _ASTROPT_QHA_ + "DCINITCOORDSLABELS=" + "|" + _ASTROPT_AAPL_ + "DCINITCOORDSLABELS=" + "|" + _ASTROPT_APL_OLD_ + "DCINITCOORDSLABELS="), USER_DC_INITCOORDS_LABELS_OPTION.option, USER_DC_INITCOORDS_LABELS_OPTION.xscheme); //CO 170601
      if(!found_user_path && USER_DC_INITLATTICE_OPTION.isentry) {
        USER_DC_INITLATTICE.clear(); USER_DC_INITCOORDS_FRAC.clear(); USER_DC_INITCOORDS_CART.clear();
        USER_DC_INITLATTICE = USER_DC_INITLATTICE_OPTION.content_string;
        logger << "Setting " << _ASTROPT_ << "DCINITLATTICE=" << USER_DC_INITLATTICE << " (via DCINITLATTICE)." << apl::endl;
        found_user_path=true;
      }  //corey
      if(!found_user_path && USER_DC_INITSG_OPTION.isentry){
        USER_DC_INITLATTICE.clear(); USER_DC_INITCOORDS_FRAC.clear(); USER_DC_INITCOORDS_CART.clear();
        USER_DC_DCINITSG = USER_DC_INITSG_OPTION.content_int;
        USER_DC_INITLATTICE = LATTICE::SpaceGroup2LatticeVariation(USER_DC_DCINITSG, xinput.getXStr()); //xvasp.str);
        logger << "Setting " << _ASTROPT_ << "DCINITLATTICE=" << USER_DC_INITLATTICE << " (via DCINITSG)." << apl::endl;
        found_user_path=true;
      }  //corey
      if(!found_user_path && USER_DC_INITCOORDS_FRAC_OPTION.isentry){
        USER_DC_INITLATTICE.clear(); USER_DC_INITCOORDS_FRAC.clear(); USER_DC_INITCOORDS_CART.clear();
        USER_DC_INITCOORDS_FRAC = USER_DC_INITCOORDS_FRAC_OPTION.content_string;
        logger << "Setting " << _ASTROPT_ << "DCINITCOORDSFRAC=" << USER_DC_INITCOORDS_FRAC << "." << apl::endl;
        logger << "User's q-point path will be calculated along fractional coordinates [" << USER_DC_INITCOORDS_FRAC << "]." << apl::endl;
        found_user_path=true;
      }  //corey
      if(!found_user_path && USER_DC_INITCOORDS_CART_OPTION.isentry){
        USER_DC_INITLATTICE.clear(); USER_DC_INITCOORDS_FRAC.clear(); USER_DC_INITCOORDS_CART.clear();
        USER_DC_INITCOORDS_CART = USER_DC_INITCOORDS_CART_OPTION.content_string;
        logger << "Setting " << _ASTROPT_ << "DCINITCOORDSCART=" << USER_DC_INITCOORDS_CART << "." << apl::endl;
        logger << "User's q-point path will be calculated along cartesian coordinates [" << USER_DC_INITCOORDS_CART << "]." << apl::endl;
        found_user_path=true;
      }  //corey

      if(!found_user_path){
        USER_DC_INITLATTICE.clear();
        USER_DC_INITCOORDS_FRAC.clear();
        USER_DC_INITCOORDS_CART.clear();
        USER_DC_INITCOORDS_LABELS.clear();
      }

      if(found_user_path && ( USER_DC_INITCOORDS_FRAC_OPTION.isentry || USER_DC_INITCOORDS_CART_OPTION.isentry )){
        USER_DC_INITCOORDS_LABELS = USER_DC_INITCOORDS_LABELS_OPTION.content_string;
        if(USER_DC_INITCOORDS_LABELS.empty()){throw apl::APLRuntimeError("INITCOORDS set but no corresponding labels (DCINITCOORDSLABELS) found.");}
        logger << "Setting " << _ASTROPT_ << "DCUSERPATHLABELS=" << USER_DC_INITCOORDS_LABELS << "." << apl::endl;
        logger << "User's q-point path will be labeled [" << USER_DC_INITCOORDS_LABELS << "]." << apl::endl;
      }

      // DCUSERPATH, e.g., DCUSERPATH = G-X|X-U|K-G|G-L
      // default is empty string, so don't set unless isentry
      USER_DC_USERPATH_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCUSERPATH=" + "|" + _ASTROPT_QHA_ + "DCUSERPATH=" + "|" + _ASTROPT_AAPL_ + "DCUSERPATH=" + "|" + _ASTROPT_APL_OLD_ + "DCUSERPATH="), USER_DC_USERPATH_OPTION.option, USER_DC_USERPATH_OPTION.xscheme); //CO 170601
      if(USER_DC_USERPATH_OPTION.isentry){
        USER_DC_USERPATH = USER_DC_USERPATH_OPTION.content_string;
        logger << "Setting " << _ASTROPT_ << "DCUSERPATH=" << USER_DC_USERPATH << "." << apl::endl;
        logger << "User's q-point path will be calculated along [" << USER_DC_USERPATH << "]." << apl::endl;
      }
      
      //ASSUME SINGLE POINT OTHERWISE
      //if we initcoords, we also need to specify the path
      //if(found_user_path && ( USER_DC_INITCOORDS_FRAC_OPTION.isentry || USER_DC_INITCOORDS_CART_OPTION.isentry )){
      //  if(USER_DC_USERPATH.empty()){throw apl::APLRuntimeError("INITCOORDS set but no corresponding path (DCUSERPATH) found.");}
      //}

      // DCPOINTS, e.g., DCPOINTS = 100
      USER_DC_NPOINTS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DCPOINTS=" + "|" + _ASTROPT_QHA_ + "DCPOINTS=" + "|" + _ASTROPT_AAPL_ + "DCPOINTS=" + "|" + _ASTROPT_APL_OLD_ + "DCPOINTS="), USER_DC_NPOINTS_OPTION.option, USER_DC_NPOINTS_OPTION.xscheme); //CO 170601
      USER_DC_NPOINTS = USER_DC_NPOINTS_OPTION.content_int;
      logger << (USER_DC_NPOINTS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DCPOINTS=" << USER_DC_NPOINTS << "." << apl::endl;
      logger << "Each subpath will be divided into " << USER_DC_NPOINTS << " points." << apl::endl;
    }

    // DOS, e.g., DOS = yes
    CALCULATE_PHONON_DOS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DOS=" + "|" + _ASTROPT_QHA_ + "DOS=" + "|" + _ASTROPT_AAPL_ + "DOS=" + "|" + _ASTROPT_APL_OLD_ + "DOS="), CALCULATE_PHONON_DOS_OPTION.option, CALCULATE_PHONON_DOS_OPTION.xscheme); //CO 170601
    logger << (CALCULATE_PHONON_DOS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DOS=" << (CALCULATE_PHONON_DOS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    logger << "The calculation of PHONON DENSITY of states switched " << (CALCULATE_PHONON_DOS_OPTION.option ? "ON" : "OFF") << "." << apl::endl;

    if (CALCULATE_PHONON_DOS_OPTION.option){
      // DOSMESH, e.g., DOSMESH = 20x20x20
      USER_DOS_MESH_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DOSMESH=" + "|" + _ASTROPT_QHA_ + "DOSMESH=" + "|" + _ASTROPT_AAPL_ + "DOSMESH=" + "|" + _ASTROPT_APL_OLD_ + "DOSMESH="), USER_DOS_MESH_OPTION.option, USER_DOS_MESH_OPTION.xscheme); //CO 170601
      USER_DOS_MESH = USER_DOS_MESH_OPTION.content_string;
      tokens.clear();
      apl::tokenize(USER_DOS_MESH, tokens, string(" xX"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"DOSMESH. Specify as DOSMESH=20x20x20");}
      logger << (USER_DOS_MESH_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DOSMESH=" << USER_DOS_MESH << "." << apl::endl;

      // DOSPOINTS, e.g., DOSPOINTS = 2000
      USER_DOS_NPOINTS_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DOSPOINTS=" + "|" + _ASTROPT_QHA_ + "DOSPOINTS=" + "|" + _ASTROPT_AAPL_ + "DOSPOINTS=" + "|" + _ASTROPT_APL_OLD_ + "DOSPOINTS="), USER_DOS_NPOINTS_OPTION.option, USER_DOS_NPOINTS_OPTION.xscheme); //CO 170601
      USER_DOS_NPOINTS = USER_DOS_NPOINTS_OPTION.content_int;
      logger << (USER_DOS_NPOINTS_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "USER_DOS_NPOINTS=" << USER_DOS_NPOINTS << "." << apl::endl;
      logger << "The phonon density of states will be calculated for " << USER_DOS_NPOINTS << " bins." << apl::endl;

      // DOSMETHOD, e.g., DOSMETHOD = LT
      USER_DOS_METHOD_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DOSMETHOD=" + "|" + _ASTROPT_QHA_ + "DOSMETHOD=" + "|" + _ASTROPT_AAPL_ + "DOSMETHOD=" + "|" + _ASTROPT_APL_OLD_ + "DOSMETHOD="), USER_DOS_METHOD_OPTION.option, USER_DOS_METHOD_OPTION.xscheme); //CO 170601
      USER_DOS_METHOD = USER_DOS_METHOD_OPTION.content_string;
      transform(USER_DOS_METHOD.begin(), USER_DOS_METHOD.end(), USER_DOS_METHOD.begin(), toupper);
      logger << (USER_DOS_METHOD_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DOSMETHOD=" << USER_DOS_METHOD << "." << apl::endl;
      if (USER_DOS_METHOD == string("LT")) {logger << "The phonon density of states will be calculated by Linear Tetrahedron Method (LT)." << apl::endl;}
      else if (USER_DOS_METHOD == string("RS")){logger << "The phonon density of states will be calculated by Root Sampling Method (RS)." << apl::endl;}
      else {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"DOSMETHOD. Specify as DOSMETHOD=LT.");}

      // DOSSMEAR, e.g., DOSSMEAR = 0.05
      USER_DOS_SMEAR_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "DOSSMEAR=" + "|" + _ASTROPT_QHA_ + "DOSSMEAR=" + "|" + _ASTROPT_AAPL_ + "DOSSMEAR=" + "|" + _ASTROPT_APL_OLD_ + "DOSSMEAR="), USER_DOS_SMEAR_OPTION.option, USER_DOS_SMEAR_OPTION.xscheme); //CO 170601
      USER_DOS_SMEAR = USER_DOS_SMEAR_OPTION.content_double;
      override_option=false;
      if (USER_DOS_METHOD == string("RS")) {USER_DOS_SMEAR = DEFAULT_APL_DOSSMEAR; override_option=true;} // Default value, it is better with this...
      logger << ((USER_DOS_SMEAR_OPTION.isentry || override_option) ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "DOSSMEAR=" << USER_DOS_SMEAR << " (overridden by RS DOSMETHOD)." << apl::endl;
      if (USER_DOS_SMEAR > 1E-6){logger << "The phonon density of states will be smooth by gaussians with sigma = " << USER_DOS_SMEAR << "." << apl::endl;}
    }

    // TP, e.g., TP = yes
    CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "TP=" + "|" + _ASTROPT_QHA_ + "TP=" + "|" + _ASTROPT_AAPL_ + "TP=" + "|" + _ASTROPT_APL_OLD_ + "TP="), CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option, CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.xscheme); //CO 170601
    logger << (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "TP=" << (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    logger << "The calculation of THERMODYNAMIC properties switched " << (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option ? "ON" : "OFF") << "." << apl::endl;
    if (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option && !CALCULATE_PHONON_DOS_OPTION.option){logger << apl::warning << "The thermodynamic properties may be calculated by default settings for DOS calculation." << apl::endl;}

    // TPT, e.g., TPT = 1000:2000:10 -> temperature from 1000 to 2000 K by step 10 K
    if (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option){
      USER_TPT_OPTION.options2entry(AflowIn, string( _ASTROPT_APL_ + "TPT=" + "|" + _ASTROPT_QHA_ + "TPT=" + "|" + _ASTROPT_AAPL_ + "TPT=" + "|" + _ASTROPT_APL_OLD_ + "TPT="), USER_TPT_OPTION.option, USER_TPT_OPTION.xscheme); //CO 170601
      tokens.clear();
      apl::tokenize(USER_TPT_OPTION.content_string, tokens, string(" :"));
      if (tokens.size() != 3) {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"TPT. Specify as TPT=1000:2000:10.");}
      USER_TP_TSTART = aurostd::string2utype<double>(tokens.at(0));
      USER_TP_TEND = aurostd::string2utype<double>(tokens.at(1));
      USER_TP_TSTEP = aurostd::string2utype<double>(tokens.at(2));
      logger << (USER_TPT_OPTION.isentry ? "Setting" : "DEFAULT") << " " << _ASTROPT_ << "TPT=" << USER_TPT_OPTION.content_string << "." << apl::endl;
      logger << "The thermodynamic properties will be calculated in temperature range <" << USER_TP_TSTART << "," << USER_TP_TEND << "> with step " << USER_TP_TSTEP << " K." << apl::endl;
    }

  } catch (std::exception& e) {
    logger << apl::error << e.what() << apl::endl;
    return;
  }

  // ///////////////////////////////////////////////////////////////////////

  try {
    // Contruct the working supercell ////////////////////////////////////

    //apl::Supercell supercell(xvasp.str,logger),supercell_test(xvasp.str,logger);    //corey, slow
    apl::Supercell supercell(xinput.getXStr(),aflags,logger); //xvasp.str, logger);  //CO  //CO181226
    apl::Supercell supercell_test = supercell;    //CO

    //   pflow::PrintDist(xinput.getXStr(),20.0,cerr);
    if (USER_SUPERCELL.empty() && USER_MINATOMS > 0) {
      stringstream aus;
      if(USER_MINATOMS_RESTICTED_OPTION){
        for (int Ni=1; USER_SUPERCELL == ""; Ni++) {
          aus.str("");
          aus << "Ni=" << Ni
              << " "
              << "supercell=" << Ni << "x" << Ni << "x" << Ni << "  natoms=" << Ni * Ni * Ni * xinput.getXStr().atoms.size(); //xvasp.str.atoms.size();
          //	  logger << aus.getXStr() << apl::endl;
          if (Ni * Ni * Ni * ((int)xinput.getXStr().atoms.size()) > (int)USER_MINATOMS) { // xvasp.str.atoms.size()) > (int)USER_MINATOMS) {
            USER_MINATOMS = 0;
            USER_SUPERCELL = aurostd::utype2string<uint>(Ni) + "X" + aurostd::utype2string<uint>(Ni) + "X" + aurostd::utype2string<uint>(Ni);
            logger << aus.str() << apl::endl;
          }
        }
      } else {
      for (double radius = 0.01; USER_SUPERCELL == ""; radius += 0.01) {
        xvector<int> dims(3);
        dims = LatticeDimensionSphere(xinput.getXStr().lattice,radius); //xvasp.str.lattice, radius);
          aus.str("");
        aus << "Radius=" << aurostd::PaddedPOST(aurostd::utype2string<double>(radius, 3), 4)
            << " "
            << " supercell=" << dims(1) << "x" << dims(2) << "x" << dims(3) << "  natoms=" << dims(1) * dims(2) * dims(3) * xinput.getXStr().atoms.size(); //xvasp.str.atoms.size();
        //	  logger << aus.getXStr() << apl::endl;
        if (dims(1) * dims(2) * dims(3) * ((int)xinput.getXStr().atoms.size()) > (int)USER_MINATOMS) { // xvasp.str.atoms.size()) > (int)USER_MINATOMS) {
          USER_MINATOMS = 0;
          USER_SUPERCELL = aurostd::utype2string<uint>(dims(1)) + "X" + aurostd::utype2string<uint>(dims(2)) + "X" + aurostd::utype2string<uint>(dims(3));
          logger << aus.str() << apl::endl;
          }
        }
      }
    }

    //      cerr << "USER_WANTS_FULL_SHELL=" << USER_WANTS_FULL_SHELL << endl;
    //    for(int i=2;i<20;i+=2) cerr << "try " << i << ": " << supercell_test.buildSuitableForShell(i,USER_WANTS_FULL_SHELL,FALSE) << endl;

    if (USER_SUPERCELL.empty() && USER_MAXSHELL > 0) {
      logger << "a Searching for suitable cell to handle " << USER_MINSHELL << " shells..." << apl::endl;
      supercell.buildSuitableForShell(USER_MAXSHELL, USER_WANTS_FULL_SHELL, TRUE);
      supercell.setupShellRestrictions(USER_MAXSHELL);
    } else if (USER_SUPERCELL.empty() && USER_MINSHELL > 0) {
      logger << "b Searching for suitable cell to handle " << USER_MINSHELL << " shells..." << apl::endl;
      supercell.buildSuitableForShell(USER_MINSHELL, USER_WANTS_FULL_SHELL, TRUE);
    } else if (USER_SUPERCELL.find_first_of("xX") != string::npos) {
      // OK, user wants their own supercell...
      tokens.clear();
      apl::tokenize(USER_SUPERCELL, tokens, string(" xX"));
      supercell.build(aurostd::string2utype<int>(tokens.at(0)),
                      aurostd::string2utype<int>(tokens.at(1)),
                      aurostd::string2utype<int>(tokens.at(2)));
      // Did he specify also regular restriction for max shell included
      // in calculation?
      if (USER_MAXSHELL > 0)
        supercell.setupShellRestrictions(USER_MAXSHELL);
    } else {
      throw apl::APLRuntimeError("The settings for supercell construction are confusing.");
    }

    // ME 180925
    // Calculate the clusters for thermal conductivity calculations
    vector<apl::ClusterSet> clusters;  // ME, default, only allocates to be passed into functions
    if (CALCULATE_TCOND_OPTION.option) {
      int max_order;
      if (USER_AAPL_FOURTH_ORDER_OPTION.option) {
        max_order = 4;
      } else {
        max_order = 3;
    }

      for (int o = 3; o <= max_order; o++) {
        apl::ClusterSet clst(logger);
        bool awakeClusterSet;
        string clust_hib_file = DEFAULT_AAPL_FILE_PREFIX + _CLUSTER_SET_FILE_[o-3];
        if (USER_HIBERNATE_OPTION.option) {
          awakeClusterSet = (aurostd::EFileExist(clust_hib_file) ||
                             aurostd::FileExist(clust_hib_file));
        } else {
          awakeClusterSet = false;
        }

        if (awakeClusterSet) {
          try {
            clst = apl::ClusterSet(clust_hib_file, supercell, USER_CUTOFF_SHELL[o-3],
                                   USER_CUTOFF_DISTANCE[o-3], o, logger);
          } catch (aurostd::xerror excpt) {
            logger << apl::warning << excpt.where() << " " << excpt.error_message << std::endl;
            logger << apl::warning << "Skipping awakening of anharmonic IFCs." << apl::endl;
            awakeClusterSet = false;
          }
        }

        if (!awakeClusterSet) {
          clst = apl::ClusterSet(supercell, USER_CUTOFF_SHELL[o-3],
                                 USER_CUTOFF_DISTANCE[o-3], logger);
          clst.build(o);
          clst.buildDistortions();
          if (USER_HIBERNATE_OPTION.option) {
            clst.writeClusterSetToFile(clust_hib_file);
          }
        }
        clusters.push_back(clst);
      }
    }

    // Calculate phonons /////////////////////////////////////////////////

    auto_ptr<apl::PhononCalculator> phcalc;
    if (USER_ENGINE == string("DM")) {
      apl::DirectMethodPC* phcalcdm = new apl::DirectMethodPC(supercell, clusters, xinput, aflags,
                                                              kflags, xflags, AflowIn, logger);
      phcalcdm->isPolarMaterial(CALCULATE_POLAR_CORRECTIONS_OPTION.option);                                                       // TRY POLAR [STEFANO]
      phcalcdm->setTCOND(CALCULATE_TCOND_OPTION.option);  // TCOND JJPR
      //phcalcdm->setGeneratePlusMinus(USER_DISTORTIONS_PLUS_MINUS_OPTION.option); //CO auto
      phcalcdm->setGeneratePlusMinus(AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option, USER_DISTORTIONS_PLUS_MINUS_OPTION.option);  //CO auto
      phcalcdm->setGenerateOnlyXYZ(USER_DISTORTIONS_XYZ_ONLY_OPTION.option);
      phcalcdm->setDistortionMagnitude(USER_DISTORTION_MAGNITUDE);
      phcalcdm->setCalculateZeroStateForces(USER_ZEROSTATE_OPTION.option);
      phcalcdm->get_special_inputs(AflowIn);  //PINKU, to include PSTRESS and LDAU_PARAMETERS in the SUPERCELL files
      phcalc.reset(phcalcdm);
    //CO generally redirects to DM, the distinction between DM and GSA is obsolete
    //} else if (USER_ENGINE == string("GSA")) {
    //  apl::GeneralizedSupercellApproach* gsa = new apl::GeneralizedSupercellApproach(supercell, strPair, xinput, aflags, kflags, xflags, logger);//xvasp, aflags, kflags, vflags, logger);  //Modified JJPR
    //  //gsa->setGeneratePlusMinus(USER_DISTORTIONS_PLUS_MINUS_OPTION.option); //CO auto
    //  gsa->setGeneratePlusMinus(AUTO_DISTORTIONS_PLUS_MINUS_OPTION.option, USER_DISTORTIONS_PLUS_MINUS_OPTION.option);  //CO auto
    //  gsa->setGenerateOnlyXYZ(USER_DISTORTIONS_XYZ_ONLY_OPTION.option);
    //  gsa->setDistortionMagnitude(USER_DISTORTION_MAGNITUDE);
    //  gsa->setTensor(CALCULATE_TCOND_OPTION.option);  // TCOND JJPR
    //  gsa->setSumRule(USER_EPS_SUM);           // TCOND JJPR
    //  //phcalcdm->setCalculateZeroStateForces(USER_ZEROSTATE_OPTION.option);
    //  phcalc.reset(gsa);
    } else if (USER_ENGINE == string("LR")) {
      phcalc.reset(new apl::LinearResponsePC(supercell, clusters, xinput, aflags,
                                             kflags, xflags, AflowIn, logger));
      phcalc->setTCOND(CALCULATE_TCOND_OPTION.option);  // TCOND JJPR
      phcalc->isPolarMaterial(CALCULATE_POLAR_CORRECTIONS_OPTION.option);
      //phcalcdm->setCalculateZeroStateForces(USER_ZEROSTATE_OPTION.option);
    } else {throw apl::APLRuntimeError("Wrong setting in the "+_ASTROPT_+"ENGINE. Set DM or LR only.");}

      //QHA/SCQHA/QHA3P  START //PN180705
      // Create directories for QHA/SCQHA/QHA3P
      // The pointer pheos should be called before creation of apl.xml
      auto_ptr<apl::QHA_AFLOWIN_CREATOR> pheos;
      if(CALCULATE_GRUNEISEN_OPTION.option || CALCULATE_GRUNEISEN_A_OPTION.option || CALCULATE_GRUNEISEN_B_OPTION.option || CALCULATE_GRUNEISEN_C_OPTION.option ||
         CALCULATE_SCQHA_OPTION.option || CALCULATE_SCQHA_A_OPTION.option || CALCULATE_SCQHA_B_OPTION.option || CALCULATE_SCQHA_C_OPTION.option ||
         CALCULATE_QHA3P_OPTION.option || CALCULATE_QHA3P_A_OPTION.option || CALCULATE_QHA3P_B_OPTION.option || CALCULATE_QHA3P_C_OPTION.option)
	{

        pheos.reset(new apl::QHA_AFLOWIN_CREATOR(supercell, clusters, xinput, aflags,
                                                 kflags, xflags, AflowIn, logger));

        pheos->setGP(CALCULATE_GRUNEISEN_OPTION.option, CALCULATE_GRUNEISEN_A_OPTION.option, CALCULATE_GRUNEISEN_B_OPTION.option, CALCULATE_GRUNEISEN_C_OPTION.option);
          if( CALCULATE_SCQHA_OPTION.option || CALCULATE_SCQHA_A_OPTION.option || CALCULATE_SCQHA_B_OPTION.option || CALCULATE_SCQHA_C_OPTION.option )
	    {
        pheos->setSCGP(CALCULATE_SCQHA_OPTION.option, CALCULATE_SCQHA_A_OPTION.option, CALCULATE_SCQHA_B_OPTION.option, CALCULATE_SCQHA_C_OPTION.option);
	    }
          if( CALCULATE_QHA3P_OPTION.option || CALCULATE_QHA3P_A_OPTION.option || CALCULATE_QHA3P_B_OPTION.option || CALCULATE_QHA3P_C_OPTION.option )
	    {
	      pheos->setSCGP(CALCULATE_QHA3P_OPTION.option, CALCULATE_QHA3P_A_OPTION.option, CALCULATE_QHA3P_B_OPTION.option, CALCULATE_QHA3P_C_OPTION.option);
	    }

        pheos->setGP_VOL_DISTORTION(GP_DISTORTION);
	  if(CALCULATE_SCQHA_OPTION.option || CALCULATE_SCQHA_A_OPTION.option || CALCULATE_SCQHA_B_OPTION.option || CALCULATE_SCQHA_C_OPTION.option ||
	     CALCULATE_QHA3P_OPTION.option || CALCULATE_QHA3P_A_OPTION.option || CALCULATE_QHA3P_B_OPTION.option|| CALCULATE_QHA3P_C_OPTION.option){
          pheos->setSCGP_VOL_DISTORTION(SCQHA_DISTORTION);
        }
        if(CALCULATE_EOS_OPTION.option){
      pheos->setEOS(CALCULATE_EOS_OPTION.option);
          pheos->setEOS_distortion_range(EOS_DISTORTION_START, EOS_DISTORTION_END, EOS_DISTORTION_DISTORTION_INC);
        pheos->setEOS_STATIC_KPPRA(EOS_STATIC_KPPRA);
          pheos->setEOS_NEDOS(NEDOS);
          pheos->set_edos_accurate(EDOS_ACURATE_OPTION.option);
      }
        pheos->run_qha();
        pheos->close_log();
    }
      //QHA/SCQHA/QHA3P END

    // ME180820 - set up VASP calculations for thermal conductivity calculations
    bool aapl_stagebreak;
    if (CALCULATE_TCOND_OPTION.option) {
      aapl_stagebreak = phcalc->buildVaspAAPL(phcalc->_clusters[0]);
      if (USER_AAPL_FOURTH_ORDER_OPTION.option) {
        aapl_stagebreak = (phcalc->buildVaspAAPL(phcalc->_clusters[1]) || aapl_stagebreak);
      }
    } else {
      aapl_stagebreak = false;
    }

    // Run or awake
    bool isHibFileAvailable = aurostd::EFileExist(DEFAULT_APL_HARMIFC_FILE);  //|| //CO
    //aurostd::FileExist(string("apl.xml")); //CO

    if (USER_HIBERNATE_OPTION.option && isHibFileAvailable) {
      if (aapl_stagebreak) {
        throw apl::APLStageBreak();  // ME 180830
      }
      try {
        phcalc->awake();
      } catch (apl::APLLogicError& e) {
        logger << apl::warning << e.what() << apl::endl;
        logger << apl::warning << "Skipping awakening..." << apl::endl;
        isHibFileAvailable = false;
      }
    }

    if (!isHibFileAvailable) {
      phcalc->run(aapl_stagebreak);  // ME180830 -- added stagebreak bool
      if (USER_HIBERNATE_OPTION.option)
        phcalc->hibernate();
    }
      //QHA/SCQHA/QHA3P START //PN180705
      //Store synamical matrics and PDOS from different distorted directores
      if(CALCULATE_GRUNEISEN_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_GRUNEISEN_A_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_GRUNEISEN_B_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_GRUNEISEN_C_SUBDIRECTORIES_OPTION.option)
        { 
          apl::QHAsubdirectoryData store(*phcalc, logger);
          store.setdir_prefix(_TMPDIR_);
      string dirname = store.getdir_name(aflags.Directory);
          store.set_gp_vol_distortion(GP_DISTORTION);
          vector<string> tokens;
        apl::tokenize(USER_DOS_MESH, tokens, string(" xX"));
        //create uniform q-mesh
        store.createMPmesh(aurostd::string2utype<int>(tokens.at(0)),
                           aurostd::string2utype<int>(tokens.at(1)),
                           aurostd::string2utype<int>(tokens.at(2)),
                           phcalc->getInputCellStructure());
          tokens.clear(); //PN180705

          //check the distorted directort contains Gruneisen ON //PN180705
          if(store.check_GP()){ //PN180705
            //store dynamical matrices //PN180705
            store.create_dm(); //PN180705
        apl::PhononDispersionCalculator pdisc(*phcalc, logger);
        
        // Init path according to the aflow's definition for elec. struc.
        if((!USER_DC_INITCOORDS_LABELS.empty())&&(!USER_DC_INITCOORDS_FRAC.empty())){pdisc.initPathCoords(USER_DC_INITCOORDS_FRAC,USER_DC_INITCOORDS_LABELS,USER_DC_NPOINTS,false);}
        else if((!USER_DC_INITCOORDS_LABELS.empty())&&(!USER_DC_INITCOORDS_CART.empty())){pdisc.initPathCoords(USER_DC_INITCOORDS_CART,USER_DC_INITCOORDS_LABELS,USER_DC_NPOINTS,true);}
        else {pdisc.initPathLattice(USER_DC_INITLATTICE,USER_DC_NPOINTS);} //default!
        if(!USER_DC_USERPATH.empty()){pdisc.setPath(USER_DC_USERPATH);}   //Does user want his own path?
        std::vector<xvector<double> > qpoints = pdisc.get_qpoints();
            //store dynamical matrices along path //PN180705
        store.create_pdispath(qpoints);
        qpoints.clear();
        pdisc.clear();
      }
          store.clear(); //PN180705
          phcalc->clear(); //PN180705
          return; //PN180705
        } //PN180705
      //store PDOS from different distorted directories
      if(CALCULATE_EOS_SUBDIRECTORIES_OPTION.option)
        {
          apl::QHAsubdirectoryData store(*phcalc, logger);
          store.setdir_prefix(_TMPDIR_);
          string dirname=store.getdir_name(aflags.Directory);
          vector<string> tokens;
        apl::tokenize(USER_DOS_MESH, tokens, string(" xX"));
          apl::MonkhorstPackMesh qmesh(aurostd::string2utype<int>(tokens.at(0)),
            aurostd::string2utype<int>(tokens.at(1)),
            aurostd::string2utype<int>(tokens.at(2)),
            phcalc->getInputCellStructure(), logger);
          tokens.clear();
        auto_ptr<apl::DOSCalculator> dosc;

        if (USER_DOS_METHOD == string("LT"))
          dosc.reset(new apl::LinearTetrahedronMethod(*phcalc, qmesh, logger));
        else if (USER_DOS_METHOD == string("RS"))
          dosc.reset(new apl::RootSamplingMethod(*phcalc, qmesh, logger));
        else
          throw apl::APLRuntimeError("Unknown DOS method. Check "+_ASTROPT_+"DOSMETHOD command.");

        // Calculate DOS
        dosc->calc(USER_DOS_NPOINTS, USER_DOS_SMEAR);
        if (CALCULATE_PHONON_DOS_OPTION.option) dosc->writePDOS(_TMPDIR_, dirname);
        dosc->clear();
          qmesh.clear();
          store.clear();
          phcalc->clear();
          return;
        }
      //SCQHA and QHA3P save dynamical matrices and PDOS from different distorted directories
      if(CALCULATE_SCQHA_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_SCQHA_A_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_SCQHA_B_SUBDIRECTORIES_OPTION.option ||
         CALCULATE_SCQHA_C_SUBDIRECTORIES_OPTION.option)
        {
          {
            apl::QHAsubdirectoryData store(*phcalc, logger);
            store.setdir_prefix(_TMPDIR_);
            string dirname=store.getdir_name(aflags.Directory);
            store.set_sc_vol_distortion(SCQHA_DISTORTION);

            vector<string> tokens;
            apl::tokenize(USER_DOS_MESH,tokens,string(" xX"));
            //create uniform q-mesh
            store.createMPmesh(aurostd::string2utype<int>(tokens.at(0)),
                               aurostd::string2utype<int>(tokens.at(1)),
                               aurostd::string2utype<int>(tokens.at(2)),
                               phcalc->getInputCellStructure());
            tokens.clear();

            if(store.check_SCQHA())
              {
                store.create_dm();
                apl::PhononDispersionCalculator pdisc(*phcalc,logger);

		// Init path according to the aflow's definition for elec. struc.
		if((!USER_DC_INITCOORDS_LABELS.empty())&&(!USER_DC_INITCOORDS_FRAC.empty())){pdisc.initPathCoords(USER_DC_INITCOORDS_FRAC,USER_DC_INITCOORDS_LABELS,USER_DC_NPOINTS,false);}
		else if((!USER_DC_INITCOORDS_LABELS.empty())&&(!USER_DC_INITCOORDS_CART.empty())){pdisc.initPathCoords(USER_DC_INITCOORDS_CART,USER_DC_INITCOORDS_LABELS,USER_DC_NPOINTS,true);}
		else {pdisc.initPathLattice(USER_DC_INITLATTICE,USER_DC_NPOINTS);} //default!
		if(!USER_DC_USERPATH.empty()){pdisc.setPath(USER_DC_USERPATH);}   //Does user want his own path?
                std::vector< xvector<double> > qpoints=pdisc.get_qpoints();
                store.create_pdispath(qpoints);
                qpoints.clear();
                pdisc.clear();
              }
            {
              vector<string> tokens;
              apl::tokenize(USER_DOS_MESH,tokens,string(" xX"));
              apl::MonkhorstPackMesh qmesh(aurostd::string2utype<int>(tokens.at(0)),
                                           aurostd::string2utype<int>(tokens.at(1)),
                                           aurostd::string2utype<int>(tokens.at(2)),
                                           phcalc->getInputCellStructure(),logger);
              tokens.clear();
              auto_ptr<apl::DOSCalculator> dosc;
              if( USER_DOS_METHOD == string("LT") )
                dosc.reset( new apl::LinearTetrahedronMethod(*phcalc,qmesh,logger) );
              else if( USER_DOS_METHOD == string("RS") )
                dosc.reset( new apl::RootSamplingMethod(*phcalc,qmesh,logger) );
              else
                throw apl::APLRuntimeError("Unknown DOS method. Check [AFLOW_PHONONS]DOSMETHOD command.");
              // Calculate DOS
              dosc->calc(USER_DOS_NPOINTS,USER_DOS_SMEAR);
              if(CALCULATE_PHONON_DOS_OPTION.option)dosc->writePDOS(_TMPDIR_, dirname);
              dosc->clear();
              qmesh.clear();
      }
      store.clear();
      phcalc->clear();
          }
      return;
    }
      //PINKU QHA/SCQHA/QHA3P  END

    // Get the format of frequency desired by user ///////////////////////

    apl::IPCFreqFlags frequencyFormat = apl::NONE;

    if (!USER_FREQFORMAT.empty()) {
      // Convert format to machine representation
      tokens.clear();
      apl::tokenize(USER_FREQFORMAT, tokens, string(" |:;,"));
      for (uint i = 0; i < tokens.size(); i++) {
        if (tokens.at(i) == string("OMEGA")) {
          frequencyFormat |= apl::OMEGA;
          continue;
        }
        if (tokens.at(i) == string("HERTZ")) {
          frequencyFormat |= apl::HERTZ;
          continue;
        } else if (tokens.at(i) == string("THZ")) {
          frequencyFormat |= apl::THZ;
          continue;
        } else if (tokens.at(i) == string("CM-1") || tokens.at(i) == string("RECIPROCAL_CM")) {
          frequencyFormat |= apl::RECIPROCAL_CM;
          continue;
        } else if (tokens.at(i) == string("MEV")) {
          frequencyFormat |= apl::MEV;
          continue;
        }
        if (tokens.at(i) == string("ALLOW_NEGATIVE")) {
          frequencyFormat |= apl::ALLOW_NEGATIVE;
          continue;
        }
      }
      // Check if there was specified unit keyword...
      if (((frequencyFormat & ~apl::OMEGA) & ~apl::ALLOW_NEGATIVE) == apl::NONE)
        throw apl::APLLogicError("The mishmash frequency format.");
    } else {
      frequencyFormat = apl::THZ | apl::ALLOW_NEGATIVE;
    }

    //////////////////////////////////////////////////////////////////////
      //high-symmery qpoint auto pointer [PINKU] //PN180705
      auto_ptr<apl::PhononHSQpoints> ptr_hsq;
      bool is_negative_freq=false;
      bool scqha_is_vol_err=false;
      //high-symmery qpoint auto pointer END [PINKU]

    if (CALCULATE_PHONON_DISPERSIONS_OPTION.option) {
      apl::PhononDispersionCalculator pdisc(*phcalc, logger);

      // Init path according to the aflow's definition for elec. struc.
      if((!USER_DC_INITCOORDS_LABELS.empty())&&(!USER_DC_INITCOORDS_FRAC.empty())){pdisc.initPathCoords(USER_DC_INITCOORDS_FRAC,USER_DC_INITCOORDS_LABELS,USER_DC_NPOINTS,false);}
      else if((!USER_DC_INITCOORDS_LABELS.empty())&&(!USER_DC_INITCOORDS_CART.empty())){pdisc.initPathCoords(USER_DC_INITCOORDS_CART,USER_DC_INITCOORDS_LABELS,USER_DC_NPOINTS,true);}
      else {pdisc.initPathLattice(USER_DC_INITLATTICE,USER_DC_NPOINTS);} //default!
      if(!USER_DC_USERPATH.empty()){pdisc.setPath(USER_DC_USERPATH);}   //Does user want his own path?

      // Calculate frequencies on path
      pdisc.calc(frequencyFormat);

      // Write results into PDIS file
      pdisc.writePDIS();
	//QHA/SCQHA/QHA3P  START //PN180705
	//////////////////////////////////////////////////////////////////////
        ptr_hsq.reset(new apl::PhononHSQpoints(logger));
        ptr_hsq->read_qpointfile();
      //compute Gruneisen dispersion curve
        if (CALCULATE_GRUNEISEN_OPTION.option ||
            CALCULATE_GRUNEISEN_A_OPTION.option ||
            CALCULATE_GRUNEISEN_B_OPTION.option ||
            CALCULATE_GRUNEISEN_C_OPTION.option)
          {
            apl::QHA qha(*phcalc, *pheos, logger);
            qha.get_tmp_dir_name(_TMPDIR_);
            qha.set_cutoff_freq(CUTOFF_FREQ);
            if(qha.set_imported_variables())
              {
                if(qha.calculation_gruneisen(ptr_hsq->get_qpoints()))
                  {
                    qha.write_gruneisen_parameter_path(ptr_hsq->get_path(), ptr_hsq->get_path_segment());
                    is_negative_freq=qha.get_is_negative_freq();
        }
        //[OBSOLETE PN180705]//clear used variables
        //[OBSOLETE PN180705]path.clear();
        //[OBSOLETE PN180705]path_segment.clear();
        //[OBSOLETE PN180705]qpoints.clear();
        //[OBSOLETE PN180705]qh->clear();
      }
            qha.clear();
          }
	//QHA/SCQHA/QHA3P  END
    }

    //////////////////////////////////////////////////////////////////////

    if (CALCULATE_PHONON_DOS_OPTION.option || CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option) {
      // Generate mesh for calculation of DOS...
      tokens.clear();
      apl::tokenize(USER_DOS_MESH, tokens, string(" xX"));
      apl::MonkhorstPackMesh qmesh(aurostd::string2utype<int>(tokens.at(0)),
                                   aurostd::string2utype<int>(tokens.at(1)),
                                   aurostd::string2utype<int>(tokens.at(2)),
                                   phcalc->getInputCellStructure(), logger);

      // Setup the DOS engine which is used also for thermodynamic properties
      auto_ptr<apl::DOSCalculator> dosc;
      if (USER_DOS_METHOD == string("LT"))
        dosc.reset(new apl::LinearTetrahedronMethod(*phcalc, qmesh, logger));
      else if (USER_DOS_METHOD == string("RS"))
        dosc.reset(new apl::RootSamplingMethod(*phcalc, qmesh, logger));
      else
        throw apl::APLRuntimeError("Unknown DOS method. Check "+_ASTROPT_+"DOSMETHOD command.");

      // Calculate DOS
      dosc->calc(USER_DOS_NPOINTS, USER_DOS_SMEAR);
      if (CALCULATE_PHONON_DOS_OPTION.option) {
        dosc->writePDOS();
      }

      // Calculate thermal properties
      if (CALCULATE_THERMODYNAMIC_PROPERTIES_OPTION.option) {
        if (!dosc->hasNegativeFrequencies()) {
          apl::ThermalPropertiesCalculator tpc(*dosc, logger);
          tpc.writeTHERMO(USER_TP_TSTART, USER_TP_TEND, USER_TP_TSTEP);
	    //QHA/SCQHA/QHA3P START //PN180705
          //calculate Gruneisen
            vector<string> tokens;
            apl::tokenize(USER_DOS_MESH,tokens,string(" xX"));
            vector<int> sc_size(3,0);
            sc_size[0]=aurostd::string2utype<int>(tokens.at(0));
            sc_size[1]=aurostd::string2utype<int>(tokens.at(1));
            sc_size[2]=aurostd::string2utype<int>(tokens.at(2));
              tokens.clear();
            apl::UniformMesh umesh(logger);
            //calculate group velocities
            umesh.create_uniform_mesh(sc_size[0],sc_size[1],sc_size[2],phcalc->getInputCellStructure());
            if(CALCULATE_GROUPVELOCITY_OPTION.option){
              apl::GroupVelocity vg(*phcalc, umesh, logger);
              if(vg.check_negative_frequencies()){
                vg.write();
                vg.clear();
              }}
                //atomic displacement calculations
            if(CALCULATE_ATOMIC_DISPLACEMENT_OPTION.option){
              apl::AtomicDisplacements ad(*phcalc, umesh, logger);
              ad.set_frequency_cutoff(CUTOFF_FREQ);
              ad.populate_variables(phcalc->getInputCellStructure());
              if(ad.eigen_solver()){
                ad.thermal_displacements(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP);
                ad.write_normal_mode_direction(ptr_hsq->get_hs_kpoints());
              }}
            //QHA calculate Gruneisen
            std::vector< std::vector< double> > scqha_tv;
            if (CALCULATE_GRUNEISEN_OPTION.option   ||
                CALCULATE_GRUNEISEN_A_OPTION.option ||
                CALCULATE_GRUNEISEN_B_OPTION.option ||
                CALCULATE_GRUNEISEN_C_OPTION.option) {
              if(!is_negative_freq){
                if(umesh.get_kpoints().size()==0){
                  umesh.create_uniform_mesh(sc_size[0],sc_size[1],sc_size[2],phcalc->getInputCellStructure());
                }
                apl::QHA qha(*phcalc, *pheos, logger);
                qha.get_tmp_dir_name(_TMPDIR_);
                qha.set_cutoff_freq(CUTOFF_FREQ);
                if(qha.set_imported_variables())
                  {
                    //QHA Grunneisen parameter calculations     
                    if(qha.calculation_gruneisen(&umesh))
                      {
                        qha.write_gruneisen_parameter_mesh();
                        qha.Writeaverage_gp(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP);
                        is_negative_freq=qha.get_is_negative_freq();
                  }
                  //[OBSOLETE PN180705]eos.clear();
                }
		if(CALCULATE_EOS_OPTION.option)
                  {
                    //QHA EOS calculations
                    //[CO181202 - NOT USED]apl::QH_ENERGIES eos_ens(*phcalc, *pheos, logger);
                    apl::QH_ENERGIES eos_ens(*pheos, logger);
                    eos_ens.get_tmp_dir_name(_TMPDIR_);
                    eos_ens.get_xtracture(phcalc->getInputCellStructure());
                    if(eos_ens.get_qha_energies())
                      {     
                        apl::QHAEOS qheos(qha, eos_ens, logger);
                        qheos.set_fitting_type(FITTING_TYPE);
                        if(qheos.setvariables())
                          { 
                            qheos.set_include_ele(INCLUDE_ELE_OPTION.option);
                            qheos.cal_qheos(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP, tpc);
              }
                        qheos.clear();
                        qha.clear();
                        //QHA3P and SCQHA calculations
                        if( (CALCULATE_SCQHA_OPTION.option || 
                             CALCULATE_SCQHA_A_OPTION.option || 
                             CALCULATE_SCQHA_B_OPTION.option || 
                             CALCULATE_SCQHA_C_OPTION.option ||
                             CALCULATE_QHA3P_OPTION.option || 
                             CALCULATE_QHA3P_A_OPTION.option || 
                             CALCULATE_QHA3P_B_OPTION.option || 
                             CALCULATE_QHA3P_C_OPTION.option) && (!scqha_is_vol_err) )
                          {         
                            if(!is_negative_freq){
                              apl::SCQHA_QHA3P scqha(*phcalc, *pheos, logger);
                              scqha.get_tmp_dir_name(_TMPDIR_);
                              scqha.set_cutoff_freq(CUTOFF_FREQ);
                              if(scqha.set_imported_variables())
                                {    
                                  //QHA3P Gruneisen parameter calculations 
                                  if(scqha.calculation_gruneisen(&umesh))
                                    {
				      if(kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C)
					{ //PN 180719
                                      scqha.write_gruneisen_parameter_mesh();
                                      scqha.Writeaverage_gp(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP);
                                      is_negative_freq=scqha.get_is_negative_freq();
            }
				    }
            //[OBSOLETE PN180705]pheos->clear();
          }
                              if(!is_negative_freq) //PN180705
                                {   
                                  //SCQHA EOS
				  if(kflags.KBIN_PHONONS_CALCULATION_SCQHA || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A)
				    {
                                  apl::SCQHAEOS scqhaeos(scqha, eos_ens, logger);
                                  if(scqhaeos.import_variables())
                                    { 
                                      scqhaeos.set_input_temperature(scqha_pdis_T);
                                      scqhaeos.sccycle(USER_TP_TSTART,USER_TP_TEND, 0.1);
                                      scqha_tv=scqhaeos.get_TV_data();
                                    } 
                                  scqhaeos.clear();
				    }
                                  //QHA3P EOS
				  if(kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C)
				    {
                                  apl::QHA3POINTS qha3points(scqha, eos_ens, logger);
                                  if(qha3points.import_variables())
                                    {
                                      qha3points.set_include_ele(INCLUDE_ELE_OPTION.option);
                                      qha3points.qha3pts_temperature_loop(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP, tpc);
                                    }
                                  qha3points.clear();
                                }
                                }
                              eos_ens.clear();
                              scqha.clear();
                            }
                          }
                      }
                  }
	      } //QHA3P and SCQHA calculations    
	    } else if(kflags.KBIN_PHONONS_CALCULATION_SCQHA || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A ||
	             kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C)
	      {
              //[CO181202 - NOT USED]apl::QH_ENERGIES eos_ens(*phcalc, *pheos, logger);
              apl::QH_ENERGIES eos_ens(*pheos, logger);
              eos_ens.get_xtracture(phcalc->getInputCellStructure());
              eos_ens.get_tmp_dir_name(_TMPDIR_);
              if(eos_ens.get_scqha_energies())
                {
                  if(!is_negative_freq){
                    if(umesh.get_kpoints().size()==0){
                      umesh.create_uniform_mesh(sc_size[0],sc_size[1],sc_size[2],phcalc->getInputCellStructure());
                    }
                    apl::SCQHA_QHA3P scqha(*phcalc, *pheos, logger);
                    scqha.get_tmp_dir_name(_TMPDIR_);
                    scqha.set_cutoff_freq(CUTOFF_FREQ);
                    if(scqha.set_imported_variables())
                      {
                        //QHA3P Gruneisen parameter calculation
                        if(scqha.calculation_gruneisen(&umesh))
                          {
			      if(kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C)
				{ //PN 180719
                            scqha.write_gruneisen_parameter_mesh();
                            scqha.Writeaverage_gp(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP);
                            is_negative_freq=scqha.get_is_negative_freq();
				}
                          }
                      }
                    if(!is_negative_freq)
                      {
                        //SCQHA calculations
			  if(kflags.KBIN_PHONONS_CALCULATION_SCQHA || kflags.KBIN_PHONONS_CALCULATION_SCQHA_A || kflags.KBIN_PHONONS_CALCULATION_SCQHA_B || kflags.KBIN_PHONONS_CALCULATION_SCQHA_C)
			    {
                        apl::SCQHAEOS scqhaeos(scqha, eos_ens, logger);
                        if(scqhaeos.import_variables())
                          {
                            scqhaeos.sccycle(USER_TP_TSTART,USER_TP_TEND, 0.1);
                          }
                        scqhaeos.clear();
			    }
                        //QHA3P calculations
			  if(kflags.KBIN_PHONONS_CALCULATION_QHA3P || kflags.KBIN_PHONONS_CALCULATION_QHA3P_A || kflags.KBIN_PHONONS_CALCULATION_QHA3P_B || kflags.KBIN_PHONONS_CALCULATION_QHA3P_C)
			    {
                        apl::QHA3POINTS qha3points(scqha, eos_ens, logger);
                        if(qha3points.import_variables())
                          {
                            qha3points.set_include_ele(INCLUDE_ELE_OPTION.option);
                            qha3points.qha3pts_temperature_loop(USER_TP_TSTART,USER_TP_TEND,USER_TP_TSTEP, tpc);
                          }
                        qha3points.clear();
                      }
			}
                    scqha.clear();
                  }
                }
              eos_ens.clear();
	    }
	    //compute SCQHA temperature dependent dispersion curve
	    if((scqha_tv.size()!=0) && (SCQHA_PDIS_T_OPTION.option)){
	      apl::T_spectra_SCQHA_QHA3P scqha_T(*phcalc, *pheos, logger);
	      scqha_T.get_tmp_dir_name(_TMPDIR_);
	      scqha_T.set_cutoff_freq(CUTOFF_FREQ);
	      scqha_T.get_input_data(scqha_tv);
	      if(scqha_T.set_imported_variables())
		{
		  if(scqha_T.calculation_freqs(ptr_hsq->get_qpoints()))
		    {
		      scqha_T.calculate_pdis_T(ptr_hsq->get_path(), ptr_hsq->get_path_segment());
		    }
		}
	      scqha_T.clear();
	    }
            ptr_hsq->clear();
            umesh.clear();
	    //QHA/SCQHA/QHA3P END
          tpc.clear();
        } else {
          logger << apl::warning << "There are negative frequencies in DOS. The calculation of thermal properties has been skipped." << apl::endl;
        }
      }

      // Clear old stuff
      dosc->clear();
      //delete dosc; //auto_ptr will do
      qmesh.clear();
    }

    //////////////////////////////////////////////////////////////////////

    /***************** Thermal Conductivity Calculations *****************/
    if (CALCULATE_TCOND_OPTION.option) {
      // Get anharmonic force constants
      phcalc->setAnharmonicOptions(USER_AAPL_MAX_ITER, USER_AAPL_MIX, USER_EPS_SUM);
      bool awakeAnharmIFCs;
      for (uint i = 0; i < phcalc->_clusters.size(); i++) {
        string ifcs_hib_file = DEFAULT_AAPL_FILE_PREFIX + _ANHARMONIC_IFCS_FILE_[i];
        if (USER_HIBERNATE_OPTION.option) {
          awakeAnharmIFCs = (aurostd::EFileExist(ifcs_hib_file) ||
                             aurostd::FileExist(ifcs_hib_file));
        } else {
          awakeAnharmIFCs = false;
        }

        if (awakeAnharmIFCs) {
          try {
            phcalc->readAnharmonicIFCs(ifcs_hib_file, phcalc->_clusters[i]);
          } catch (aurostd::xerror excpt) {
            logger << apl::warning << excpt.where() << " " << excpt.error_message << std::endl;
            logger << apl::warning << "Skipping awakening of anharmonic IFCs." << apl::endl;
            awakeAnharmIFCs = false;
          }
        }

        if (!awakeAnharmIFCs) {
          phcalc->calculateAnharmonicIFCs(phcalc->_clusters[i]);
          if (USER_HIBERNATE_OPTION.option) {
            phcalc->_anharmonicIFCs[i].writeIFCsToFile(ifcs_hib_file);
          }
        }
      }

      // Do the thermal conductivity calculation
      logger << "Starting thermal conductivity calculations." << apl::endl;
      apl::TCONDCalculator tcond(*phcalc, supercell, logger);

      tcond.setCalculationOptions(USER_BTE, CALCULATE_ISOTOPE_OPTION.option,
                                  CALCULATE_CUMULATIVEK_OPTION.option,
                                  USER_AAPL_FOURTH_ORDER_OPTION.option,
                                  CALCULATE_BOUNDARY_OPTION.option,
                                  USER_NANO_SIZE, USER_TCT_TSTART,
                                  USER_TCT_TEND, USER_TCT_TSTEP);

      // Get q-points
      if (USER_THERMALGRID.find_first_of("xX") != string::npos) {
        logger << "Preparing a q-mesh of " << USER_THERMALGRID << "." << apl::endl;
        tokens.clear();
        apl::tokenize(USER_THERMALGRID, tokens, string(" xX"));
        xvector<int> grid(3);
        for (int i = 0; i < 3; i++) {
          grid[i+1] = aurostd::string2utype<int>(tokens.at(i));
        }
        tcond.buildQpoints(grid);
      } else {
        string function = "apl::RunPhonons_APL";
        string message = "Incorrect q-point settings in the THERMALGRID flag.";
        throw aurostd::xerror(function, message, _INPUT_ILLEGAL_);
      }

      // Calculate lattice thermal conductivity
      tcond.calculateFrequenciesGroupVelocities();
      tcond.calculateTransitionProbabilities(3);
      if (USER_AAPL_FOURTH_ORDER_OPTION.option) {
        tcond.calculateTransitionProbabilities(4);
    }
      tcond.calculateThermalConductivity();
 
      tcond.clear();
    }
    /*************** End Thermal Conductivity Calculations ****************/

    phcalc->clear();
    // delete phcalc; //auto_ptr will do
    supercell.clear();
  } catch (apl::APLStageBreak& e) {
    logger << apl::notice << "Stopped. Waiting for required calculations..." << apl::endl;  //CO181226
  } catch (std::exception& e) {
    logger << apl::error << e.what() << apl::endl;
  }
}
}
#endif

//////////////////////////////////////////////////////////////////////////////

bool PHON_RunPhonons(const xstructure& _str,
                     _aflags& aflags,
                     const double& _radius,
                     const bool& osswrite, ostream& oss) {
  bool Krun = FALSE;
  oss << "not implemented" << endl;
  if (0) {  // to avoid warnings
    oss << _str << endl
        << aflags.QUIET << _radius << osswrite << endl;
  }
  return Krun;
}

// ***************************************************************************
// *                                                                         *
// *             STEFANO CURTAROLO - Duke University 2003-2019              *
// *                                                                         *
// ***************************************************************************
