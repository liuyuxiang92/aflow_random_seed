// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************

#define __XLIBS_LINK
#include "aflow.h"
#define cdebug cerr
#include <algorithm>

//bool nocase_compare(char c1,char c2) {return toupper(c1)==toupper(c2);}

//#define MaxAflowInSize 65535
//string AflowIn; //[MaxAflowInSize];

#define VRUNS_MAX_CUTOFF 32768
#define DUKE_BETANEW_DEFAULT_KILL_MEM_CUTOFF 1.50
#define DUKE_QRATS_DEFAULT_KILL_MEM_CUTOFF 1.50
#define DUKE_QFLOW_DEFAULT_KILL_MEM_CUTOFF 1.50
#define MPCDF_EOS_DEFAULT_KILL_MEM_CUTOFF 1.50
#define MPCDF_DRACO_DEFAULT_KILL_MEM_CUTOFF 1.50
#define MPCDF_COBRA_DEFAULT_KILL_MEM_CUTOFF 1.50
#define MPCDF_HYDRA_DEFAULT_KILL_MEM_CUTOFF 1.50
#define MACHINE001_DEFAULT_KILL_MEM_CUTOFF 1.50  //DX20190509 - MACHINE001
#define MACHINE002_DEFAULT_KILL_MEM_CUTOFF 1.50  //DX20190509 - MACHINE002
#define CMU_EULER_DEFAULT_KILL_MEM_CUTOFF 1.50   //DX20190107 - CMU EULER

namespace aurostd {
  // ***************************************************************************
  // Function DirectoryAlreadyInDabatase
  // ***************************************************************************
  bool DirectoryAlreadyInDatabase(string directory,bool FORCE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(FORCE) return FALSE; // neglect already in the database

    // already scanned
    if(aurostd::FileExist(directory+"/ALREADY_IN_DATABASE") || aurostd::EFileExist(directory+"/ALREADY_IN_DATABASE")) return TRUE;

    // no, then scan
    if(LDEBUG) cerr << "SEARCHING" << endl;
    bool already_in_database=FALSE;
    string tmp_directory=directory;
    aurostd::StringSubst(tmp_directory,_AFLOWIN_," ");
    aurostd::StringSubst(tmp_directory,"./"," ");
    aurostd::StringSubst(tmp_directory,"/"," ");
    vector<string> tokens;
    aurostd::string2tokens(tmp_directory,tokens);
    if(tokens.size()>=2) tmp_directory=tokens.at(tokens.size()-2)+"/"+tokens.at(tokens.size()-1);
    if(tokens.size()==1) tmp_directory=tokens.at(tokens.size()-1);

    uint library=LIBRARY_NOTHING;
    // XHOST_LIBRARY_LIB0
    if(aurostd::substring2bool(directory,"LIB0"))     library=XHOST_LIBRARY_LIB0;
    // XHOST_LIBRARY_LIB1
    if(aurostd::substring2bool(directory,"LIB1"))     library=XHOST_LIBRARY_LIB1;
    // XHOST_LIBRARY_AURO
    if(aurostd::substring2bool(directory,"AURO"))     library=XHOST_LIBRARY_LIB1;    // [HISTORIC]
    // XHOST_LIBRARY_LIB2
    if(aurostd::substring2bool(directory,"LIBRARYU")) library=XHOST_LIBRARY_LIB2;
    if(aurostd::substring2bool(directory,"LIBRARYX")) library=XHOST_LIBRARY_LIB2;    // [HISTORIC]
    if(aurostd::substring2bool(directory,"LIB2"))     library=XHOST_LIBRARY_LIB2;
    // XHOST_LIBRARY_ICDS
    if(aurostd::substring2bool(directory,"ICSD"))     library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"SCINT"))    library=XHOST_LIBRARY_ICSD;    // [HISTORIC]
    if(aurostd::substring2bool(directory,"BCC"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"BCT"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"CUB"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"FCC"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"HEX"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"MCL"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"MCLC"))     library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"ORC"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"ORCC"))     library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"ORCF"))     library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"ORCI"))     library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"RHL"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"TET"))      library=XHOST_LIBRARY_ICSD;
    if(aurostd::substring2bool(directory,"TRI"))      library=XHOST_LIBRARY_ICSD;
    // XHOST_LIBRARY_ICSD
    if(aurostd::substring2bool(directory,"MAGNETIC")) library=XHOST_LIBRARY_LIB3;    // [HISTORIC]
    if(aurostd::substring2bool(directory,"LIB3"))     library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0001"))    library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0002"))    library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0001"))    library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0004"))    library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0005"))    library=XHOST_LIBRARY_LIB3;
    if(aurostd::substring2bool(directory,"T0006"))    library=XHOST_LIBRARY_LIB3;
    // XHOST_LIBRARY_LIB4
    if(aurostd::substring2bool(directory,"LIB4"))     library=XHOST_LIBRARY_LIB4;
    if(aurostd::substring2bool(directory,"Q0001"))    library=XHOST_LIBRARY_LIB4;
    // XHOST_LIBRARY_LIB5
    if(aurostd::substring2bool(directory,"LIB5"))     library=XHOST_LIBRARY_LIB5;
    if(aurostd::substring2bool(directory,"P0001"))    library=XHOST_LIBRARY_LIB5;
    // XHOST_LIBRARY_LIB6
    if(aurostd::substring2bool(directory,"LIB6"))     library=XHOST_LIBRARY_LIB6;
    if(aurostd::substring2bool(directory,"H0001"))    library=XHOST_LIBRARY_LIB6;
    // XHOST_LIBRARY_LIB7
    if(aurostd::substring2bool(directory,"LIB7"))     library=XHOST_LIBRARY_LIB7;
    // XHOST_LIBRARY_LIB8
    if(aurostd::substring2bool(directory,"LIB8"))     library=XHOST_LIBRARY_LIB8;
    // XHOST_LIBRARY_LIB9
    if(aurostd::substring2bool(directory,"LIB9"))     library=XHOST_LIBRARY_LIB9;

    // found something
    if(library!=LIBRARY_NOTHING) {
      init::InitLoadString("vLIBS");
      tokens.clear();
      string tmp;
      vector<string> vLibrary;
      if(library==XHOST_LIBRARY_LIB0) {
        if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB0" << endl;
        tokens=XHOST_Library_CALCULATED_LIB0_RAW;
      }
      if(library==XHOST_LIBRARY_LIB1) {
        if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB1" << endl;
        tokens=XHOST_Library_CALCULATED_LIB1_RAW;
      }
      if(library==XHOST_LIBRARY_LIB2) {
        if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB2" << endl;
        tokens=XHOST_Library_CALCULATED_LIB2_RAW;
      }
      if(library==XHOST_LIBRARY_LIB3) {
        if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB3" << endl;
        tokens=XHOST_Library_CALCULATED_LIB3_RAW;
      }
      if(library==XHOST_LIBRARY_LIB4) {
        if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB4" << endl;
        tokens=XHOST_Library_CALCULATED_LIB4_RAW;
      }
      if(library==XHOST_LIBRARY_LIB5) {
        if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB5" << endl;
        tokens=XHOST_Library_CALCULATED_LIB5_RAW;
      }
      if(library==XHOST_LIBRARY_LIB6) {
        if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB6" << endl;
        tokens=XHOST_Library_CALCULATED_LIB6_RAW;
      }
      if(library==XHOST_LIBRARY_LIB7) {
        if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB7" << endl;
        tokens=XHOST_Library_CALCULATED_LIB7_RAW;
      }
      if(library==XHOST_LIBRARY_LIB8) {
        if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB8" << endl;
        tokens=XHOST_Library_CALCULATED_LIB8_RAW;
      }
      if(library==XHOST_LIBRARY_LIB9) {
        if(LDEBUG) cerr << "library==XHOST_LIBRARY_LIB9" << endl;
        tokens=XHOST_Library_CALCULATED_LIB9_RAW;
      }
      for(uint i=0;i<tokens.size();i++) {
        if(aurostd::substring2bool(tokens.at(i),"/")) {
          tmp=tokens.at(i);
          aurostd::StringSubst(tmp," ","");
          vLibrary.push_back(tmp);
        }
      }
      for(int i=vLibrary.size()-1;i>=0;i--) {
        if(tmp_directory==vLibrary.at(i)) {
          already_in_database=TRUE;
          if(LDEBUG) cerr << vLibrary.at(i) << " FOUND .." << endl; // NEW
        }
      }
    }

    if(LDEBUG) cerr << "DONE..." << endl;
    if(LDEBUG) cerr << tmp_directory << endl;

    if(already_in_database) {
      //    cerr << directory << " already in database" << endl;
      aurostd::execute("cp -f "+directory+"/"+_AFLOWIN_+" "+directory+"/ALREADY_IN_DATABASE");
      aurostd::execute(DEFAULT_KZIP_BIN+" -f "+directory+"/"+_AFLOWIN_);
    }

    return already_in_database;
  }
}

using aurostd::DirectorySkipped;
using aurostd::DirectoryAlreadyInDatabase;
using aurostd::DirectoryUnwritable;

// int KBIN_MODE;

//#define KBIN_VOID_MODE 0             // just a shift
//#define KBIN_VASP_MODE 2             // for ab-initio VASP mode
//#define KBIN_XXXX_MODE 3             // for XXXX program
//#define KBIN_GRND_MODE 4             // for classical monte carlo

// GND MODE
#define KBIN_VASP_N_VPARS 32
#define _KBIN_VASP_SLEEP_ 2
#define _KBIN_LOOP_SLEEP_ 300
// PRIORITY
#define PRIORITY_PROBABILITY 0.2000
//#define PRIORITY_GREP_STRING string("grep -vi xxxxx ")
#define PRIORITY_GREP_STRING string("grep -vi PRIORITY ")

// ***************************************************************************
// KBIN::Legitimate_aflowin
// ***************************************************************************
namespace KBIN {
  bool Legitimate_aflowin(string _aflowindir,const bool& osswrite,ostringstream& oss) {
    string aflowindir=_aflowindir;
    aurostd::StringSubst(aflowindir,"//","/");

    if(aurostd::FileExist(aflowindir)) {  // file must exist
      if(!aurostd::FileEmpty(aflowindir)) { // must not be empty
        if(aurostd::substring2bool(aflowindir,_AFLOWIN_)) {  // there must be an _AFLOWIN_
          aurostd::StringSubst(aflowindir,_AFLOWIN_,"");
          if(!aurostd::FileExist(aflowindir+_AFLOWLOCK_)) { // it should be UNLOCKED
            //  if(osswrite) {oss << "MMMMM  Loading Valid File Entry = " << aflowindir << MessageTime(aflags);aurostd::PrintMessageStream(oss,XHOST.QUIET);};
            return TRUE;	
          } else { // must be unlocked
            if(osswrite) {oss << "MMMMM  Directory locked = " << aflowindir << Message(_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(oss,XHOST.QUIET);};
            return FALSE;
          }
        } else { // must contain _AFLOWIN_
          if(osswrite) {oss << "MMMMM  Not loading file without " << _AFLOWIN_ << " = " << aflowindir << Message(_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(oss,XHOST.QUIET);};
          return FALSE;
        }
      } else { // empty
        if(osswrite) {oss << "MMMMM  Not loading empty file = " << aflowindir << Message(_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(oss,XHOST.QUIET);};
        return FALSE;
      }
    } else { // unexisting
      // if(osswrite) {oss << "MMMMM  Not loading unexisting file = " << aflowindir << Message(_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(oss,XHOST.QUIET);};
      return FALSE;
    }
    return FALSE;
  }
}

namespace KBIN {
  bool Legitimate_aflowin(string aflowindir) { ostringstream aus; return KBIN::Legitimate_aflowin(aflowindir,FALSE,aus);};
}

namespace KBIN {
  void getAflowInFromAFlags(const _aflags& aflags,string& AflowIn_file,string& AflowIn,ostream& oss) {ofstream FileMESSAGE;return getAflowInFromAFlags(aflags,AflowIn_file,AflowIn,FileMESSAGE,oss);}  //CO20191110
  void getAflowInFromAFlags(const _aflags& aflags,string& AflowIn_file,string& AflowIn,ofstream& FileMESSAGE,ostream& oss) { //CO20191110
    return getAflowInFromDirectory(aflags.Directory,AflowIn_file,AflowIn,FileMESSAGE,oss);
  }
  void getAflowInFromDirectory(const string& directory,string& AflowIn_file,string& AflowIn,ostream& oss) {ofstream FileMESSAGE;return getAflowInFromDirectory(directory,AflowIn_file,AflowIn,FileMESSAGE,oss);}  //CO20191110
  void getAflowInFromDirectory(const string& directory,string& AflowIn_file,string& AflowIn,ofstream& FileMESSAGE,ostream& oss) { //CO20191110
    string soliloquy = XPID + "KBIN::getAflowInFromDirectory():";
    AflowIn_file=aurostd::CleanFileName(directory+"/"+_AFLOWIN_); //CO20200624
    if(!aurostd::FileExist(AflowIn_file)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Input file does not exist: "+AflowIn_file,_INPUT_ERROR_);}
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,"Using input file: "+AflowIn_file,FileMESSAGE,oss,_LOGGER_MESSAGE_);
    aurostd::file2string(AflowIn_file,AflowIn);
    AflowIn=aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,"AflowIn.size()="+aurostd::utype2string(AflowIn.size()),FileMESSAGE,oss,_LOGGER_MESSAGE_); //CO20200624 - check size!=0
  }
}

// ***************************************************************************
// KBIN::Main
// ***************************************************************************
namespace KBIN {
  int KBIN_Main(vector<string> argv) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "KBIN::KBIN_Main():";
    string GENERIC;
    //  string Directory;
    int i;
    ostringstream aus;
    ifstream FileAUS;
    _aflags aflags;
    aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");  //CO20190629 - a good default until we get down into vDirectory
    aurostd::StringstreamClean(aus);
    bool _VERBOSE_=FALSE;
    int XHOST_AFLOW_RUNXnumber_multiplier=3;

    std::deque<_aflags> qaflags;

    AFLOW_PTHREADS::Clean_Threads();                                    // clean threads
    // _aflags taflags[MAX_ALLOCATABLE_PTHREADS];
    // _threaded_KBIN_params params[MAX_ALLOCATABLE_PTHREADS];

    // cerr << "GMODE" << endl;
    // check BlackList **************************************************
    if(AFLOW_BlackList(XHOST.hostname)) {
      aus << "MMMMM  HOSTNAME BLACKLISTED = " << XHOST.hostname << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
      return 0;
    }
    // get KBIN **************************************************
    //  aus << "MMMMM  AFLOW: running " << _AFLOWIN_ << ": " << " "  << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;   // too much verbosity is annoying
    //  aurostd::PrintMessageStream(aus,XHOST.QUIET);                                                  // too much verbosity is annoying

    // do some running


    //   cerr << soliloquy << " XHOST.AFLOW_RUNXflag=" << XHOST.AFLOW_RUNXflag << endl; //

    if(XHOST.AFLOW_RUNXflag) {
      vector<string> tokens;
      if(aurostd::args2attachedflag(argv,"--run=")) {
        aurostd::string2tokens(aurostd::args2attachedstring(argv,"--run=","1"),tokens,"=");
        XHOST.AFLOW_RUNXnumber=aurostd::string2utype<uint>(tokens.at(tokens.size()-1));
      }
      if(aurostd::args2attachedflag(argv,"-run=")) {
        aurostd::string2tokens(aurostd::args2attachedstring(argv,"-run=","1"),tokens,"=");
        XHOST.AFLOW_RUNXnumber=aurostd::string2utype<uint>(tokens.at(tokens.size()-1));
      }
      if(XHOST.AFLOW_RUNXnumber<1) XHOST.AFLOW_RUNXnumber=1;
    }

    if(aurostd::args2flag(argv,"-runone|--runone|-run_one|--run_one|-run1|--run1|--run=1|-run=1")) { // RUNONE COMPATIBILITY
      XHOST.AFLOW_RUNXflag=TRUE;XHOST.AFLOW_RUNXnumber=1;
    }

    if(XHOST.AFLOW_RUNXflag)  {
      XHOST.AFLOW_RUNDIRflag=FALSE;
      XHOST.AFLOW_MULTIflag=FALSE;
    }
    if(XHOST.AFLOW_MULTIflag)  {
      XHOST.AFLOW_RUNDIRflag=FALSE;
      XHOST.AFLOW_RUNXflag=FALSE;
    }

    //    cerr << soliloquy << " XHOST.AFLOW_RUNXflag=" << XHOST.AFLOW_RUNXflag << endl; //

    if(XHOST.vflag_aflow.flag("LOOP")) {aus << "MMMMM  KBIN option XHOST.vflag_aflow.flag(\"LOOP\") - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(XHOST.AFLOW_RUNDIRflag) {aus << "MMMMM  KBIN option [--run] (XHOST.AFLOW_RUNDIRflag=TRUE) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(XHOST.AFLOW_MULTIflag) {aus << "MMMMM  KBIN option [--run=multi] (XHOST.AFLOW_MULTIflag=TRUE) - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(XHOST.AFLOW_RUNXflag && XHOST.AFLOW_RUNXnumber==1) {aus << "MMMMM  KBIN option [--run=1] (XHOST.AFLOW_RUNXflag=TRUE, XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(XHOST.AFLOW_RUNXflag && XHOST.AFLOW_RUNXnumber>1) {aus << "MMMMM  KBIN option [--run=N] (XHOST.AFLOW_RUNXflag=TRUE, XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}

    aflags.KBIN_RUN_AFLOWIN=TRUE;
    aflags.KBIN_GEN_VASP_FROM_AFLOWIN=FALSE;
    aflags.KBIN_GEN_AFLOWIN_FROM_VASP=FALSE;
    //DX
    aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN=FALSE;
    //DX
    aflags.KBIN_DELETE_AFLOWIN=FALSE;

    //--generate could be for VASP/AIMS/etc., need to read aflow.in later
    aflags.KBIN_GEN_GENERAL                  = aurostd::args2flag(argv,"--generate");                             //CO20180402 - we will use this to modify other GEN flags later, when reading AFLOWIN
    aflags.KBIN_GEN_VASP_FROM_AFLOWIN        = aurostd::args2flag(argv,"--generate_vasp_from_aflowin|--generate");//CO20180402 - --generate assumes vasp generation FOR NOW
    if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
      aflags.KBIN_RUN_AFLOWIN=FALSE;
    }
    aflags.KBIN_GEN_AIMS_FROM_AFLOWIN        = aurostd::args2flag(argv,"--generate_aims_from_aflowin");//CO20180402
    aflags.KBIN_GEN_AFLOWIN_FROM_VASP        = aurostd::args2flag(argv,"--generate_aflowin_from_vasp");
    aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN      = aurostd::args2flag(argv,"--generate_symmetry|--generate_sym"); //DX

    //DX
    //DX if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN || aflags.KBIN_GEN_AFLOWIN_FROM_VASP)
    if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN || aflags.KBIN_GEN_AFLOWIN_FROM_VASP || aflags.KBIN_GEN_SYMMETRY_OF_AFLOWIN || aflags.KBIN_GEN_AIMS_FROM_AFLOWIN) //CO20180409
    { //CO20200106 - patching for auto-indenting
      //DX
      XHOST.AFLOW_RUNXflag=TRUE;XHOST.AFLOW_RUNXnumber=1;
      XHOST.AFLOW_RUNDIRflag=TRUE;
      if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP) aflags.KBIN_RUN_AFLOWIN=FALSE;
    }

    aflags.KBIN_DELETE_AFLOWIN=aurostd::args2flag(argv,"--delete_aflowin");

    aflags.AFLOW_MODE_QSUB_MODE1 = aurostd::args2flag(argv,"--qsub1|-qsub1");
    aflags.AFLOW_MODE_QSUB_MODE2 = aurostd::args2flag(argv,"--qsub2|-qsub2");
    aflags.AFLOW_MODE_QSUB_MODE3 = aurostd::args2flag(argv,"--qsub3|-qsub3");

    // [OBSOLETE]  MPI=aurostd::args2flag(argv,"--MPI|--mpi");  // ABSOLUTELY otherwise the multithreads kicks
    aflags.AFLOW_FORCE_MPI = aurostd::args2flag(argv,"--MPI|--mpi");
    aflags.AFLOW_FORCE_SERIAL = aurostd::args2flag(argv,"--nompi|-nompi|--serial|-serial");
    // [OBSOLETE]  aflags.AFLOW_GLOBAL_NCPUS=aurostd::args2utype(argv,"--np",(int) 0);
    aflags.AFLOW_GLOBAL_NCPUS=aurostd::args2attachedutype<int>(argv,"--np=",(int) 0);
    // if(aflags.AFLOW_GLOBAL_NCPUS && !MPI && !aflags.AFLOW_FORCE_MPI) {}
    if(XHOST.MPI || aflags.AFLOW_FORCE_MPI) AFLOW_PTHREADS::No_Threads();

    aflags.AFLOW_MACHINE_GLOBAL.clear();
    // "MACHINE::DUKE_BETA_MPICH"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_MPICH",
        aurostd::args2flag(argv,"--machine=beta|--machine=duke_beta|--beta|--duke_beta|--machine=beta_mpich|--machine=duke_beta_mpich"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_MPICH")) XHOST.maxmem=DUKE_BETANEW_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::DUKE_BETA_OPENMPI"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_OPENMPI",aurostd::args2flag(argv,"--machine=beta_openmpi|--machine=duke_beta_openmpi"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_BETA_OPENMPI")) XHOST.maxmem=DUKE_BETANEW_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::DUKE_QRATS_MPICH"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QRATS_MPICH",aurostd::args2flag(argv,"--machine=qrats|--machine=duke_qrats|--machine=qrats_mpich|--machine=duke_qrats_mpich"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QRATS_MPICH")) XHOST.maxmem=DUKE_QRATS_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::DUKE_QFLOW_OPENMPI"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QFLOW_OPENMPI",aurostd::args2flag(argv,"--machine=qflow|--machine=duke_qflow|--machine=qflow_openmpi|--machine=duke_qflow_openmpi|--machine=quser|--machine=duke_quser|--machine=quser_openmpi|--machine=duke_quser_openmpi")); //backwards compatible //CO20180409
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_QFLOW_OPENMPI")) XHOST.maxmem=DUKE_QFLOW_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::MPCDF_EOS"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_EOS",aurostd::args2flag(argv,"--machine=eos|--machine=mpcdf_eos|--machine=eos_mpiifort|--machine=mpcdf_eos_mpiifort"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_EOS")) XHOST.maxmem=MPCDF_EOS_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::MPCDF_DRACO"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_DRACO",aurostd::args2flag(argv,"--machine=draco|--machine=mpcdf_draco|--machine=draco_mpiifort|--machine=mpcdf_draco_mpiifort"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_DRACO")) XHOST.maxmem=MPCDF_DRACO_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::MPCDF_COBRA"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_COBRA",aurostd::args2flag(argv,"--machine=cobra|--machine=mpcdf_cobra|--machine=cobra_mpiifort|--machine=mpcdf_cobra_mpiifort"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_COBRA")) XHOST.maxmem=MPCDF_COBRA_DEFAULT_KILL_MEM_CUTOFF;
    // "MACHINE::MPCDF_HYDRA"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_HYDRA",aurostd::args2flag(argv,"--machine=hydra|--machine=mpcdf_hydra|--machine=hydra_mpiifort|--machine=mpcdf_hydra_mpiifort"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MPCDF_HYDRA")) XHOST.maxmem=MPCDF_HYDRA_DEFAULT_KILL_MEM_CUTOFF;
    //DX20190509 - MACHINE001 - START
    // "MACHINE::MACHINE001"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE001",aurostd::args2flag(argv,"--machine=machine001"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE001")) XHOST.maxmem=MACHINE001_DEFAULT_KILL_MEM_CUTOFF;
    //DX20190509 - MACHINE001 - END
    //DX20190509 - MACHINE002 - START
    // "MACHINE::MACHINE002"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE002",aurostd::args2flag(argv,"--machine=machine002"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::MACHINE002")) XHOST.maxmem=MACHINE002_DEFAULT_KILL_MEM_CUTOFF;
    //DX20190509 - MACHINE002 - END
    // DUKE_MATERIALS
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_MATERIALS",aurostd::args2flag(argv,"--machine=materials|--machine=duke_materials"));
    // DUKE_AFLOWLIB
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_AFLOWLIB",aurostd::args2flag(argv,"--machine=aflowlib|--machine=duke_aflowlib"));
    // DUKE_HABANA
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::DUKE_HABANA",aurostd::args2flag(argv,"--machine=habana|--machine=duke_habana"));
    // FULTON_MARYLOU
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::FULTON_MARYLOU",aurostd::args2flag(argv,"--machine=marylou|--machine=fulton_marylou"));
    // MACHINE2
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::OHAD",aurostd::args2flag(argv,"--machine=ohad|--machine=machine2")); //CO20181113
    // MACHINE1
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::HOST1",aurostd::args2flag(argv,"--machine=host1|--machine=machine1")); //CO20181113
    //DX20190107 - CMU EULER - START
    // "MACHINE::CMU_EULER"
    aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::CMU_EULER",aurostd::args2flag(argv,"--machine=euler|--machine=cmu_euler"));
    if(aflags.AFLOW_MACHINE_GLOBAL.flag("MACHINE::CMU_EULER")) XHOST.maxmem=CMU_EULER_DEFAULT_KILL_MEM_CUTOFF;
    //DX20190107 - CMU EULER - END


    // turn off multi pthreads on specific machines
    if(aflags.AFLOW_MACHINE_GLOBAL.flag()) {
      AFLOW_PTHREADS::No_Threads();
      AFLOW_PTHREADS::FLAG=FALSE; // safety...
      AFLOW_PTHREADS::MAX_PTHREADS=1; // safety...
      //    kflags.KBIN_MPI=TRUE; // overrides the MPI for machines
      XHOST.MPI=TRUE;
    }

    vector<string> vruns;

    aflags.AFLOW_PERFORM_CLEAN=XHOST.vflag_aflow.flag("CLEAN");// || XHOST.vflag_aflow.flag("XCLEAN"));
    // [OBSOLETE] aflags.AFLOW_PERFORM_DIRECTORY=aurostd::args2flag(argv,"--DIRECTORY|--D|--d|./");
    aflags.AFLOW_PERFORM_DIRECTORY=XHOST.vflag_control.flag("VDIR");
    // [OBSOLETE] aflags.AFLOW_PERFORM_FILE=aurostd::args2flag(argv,"--FILE|--F|--f");
    aflags.AFLOW_PERFORM_FILE=XHOST.vflag_control.flag("FILE");
    aflags.AFLOW_PERFORM_ORDER_SORT=aurostd::args2flag(argv,"--sort|-sort");                    // Sorts the _AFLOWIN_ in the list
    aflags.AFLOW_PERFORM_ORDER_REVERSE=aurostd::args2flag(argv,"--reverse|--rsort|-reverse|-rsort"); // Reverse the _AFLOWIN_ in the list
    aflags.AFLOW_PERFORM_ORDER_RANDOM=aurostd::args2flag(argv,"--random|--rnd|-random|-rnd"); // Randomize the _AFLOWIN_ in the list
    aflags.AFLOW_FORCE_RUN = aurostd::args2flag(argv,"--force|-force");

    if(aflags.AFLOW_PERFORM_DIRECTORY) {aus << "MMMMM  KBIN option PERFORM_DIRECTORY - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(aflags.AFLOW_PERFORM_FILE) {aus << "MMMMM  KBIN option PERFORM_FILE - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(aflags.AFLOW_PERFORM_ORDER_SORT) {aus << "MMMMM  KBIN option ORDER_SORT - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(aflags.AFLOW_PERFORM_ORDER_REVERSE) {aus << "MMMMM  KBIN option ORDER_REVERSE - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(aflags.AFLOW_PERFORM_ORDER_RANDOM) {aus << "MMMMM  KBIN option ORDER_RANDOM - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
    if(aflags.AFLOW_FORCE_RUN) {aus << "MMMMM  KBIN option FORCE_RUN - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}

    uint maxcheck=VRUNS_MAX_CUTOFF;
    if(XHOST.AFLOW_MULTIflag) maxcheck=VRUNS_MAX_CUTOFF;
    if(XHOST.AFLOW_RUNXflag) maxcheck=XHOST.AFLOW_RUNXnumber;
    if(maxcheck==0) maxcheck=1; // safety

    // simple commands
    if(XHOST.vflag_aflow.flag("XCLEAN")) {
      KBIN::XClean(XHOST.vflag_aflow.getattachedscheme("XCLEAN"));
      return 1;
    }

    // for directory mode load them all
    if(aflags.AFLOW_PERFORM_DIRECTORY) {
      // [OBSOLETE] vruns=aurostd::args2vectorstring(argv,"--DIRECTORY|--D|--d","./");
      aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("VDIR"),vruns,",");
      if(LDEBUG) { for(uint i=0;i<vruns.size();i++) cerr << XPID << "KBIN::Main: vruns.at(i)=" << vruns.at(i) << endl;}
    } else {
      if(!aflags.AFLOW_PERFORM_FILE)
        vruns.push_back(aurostd::getPWD()); //CO20191112
      //[CO20191112 - OBSOLETE]vruns.push_back(aurostd::execute2string(XHOST.command("pwd")));
      // vruns.push_back(aurostd::execute2string(XHOST.command("pwd"))+" ./");
    }
    // if file found
    if(aflags.AFLOW_PERFORM_FILE) {
      // [OBSOLETE]   string file_name=aurostd::args2string(argv,"--FILE|--F|--f","xxxx");
      string file_name=XHOST.vflag_control.getattachedscheme("FILE");
      if(!aurostd::FileExist(file_name)) {
        aus << "EEEEE  FILE_NOT_FOUND = " << file_name  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);return 1;  //CO20200624 - previously exit
      }
      if(aurostd::FileEmpty(file_name)) {
        aus << "EEEEE  FILE_EMPTY = " << file_name  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);return 1;  //CO20200624 - previously exit
      }
      aus << "MMMMM  Loading File = " << file_name << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      vector<string> vlines;
      vlines.clear();
      aurostd::file2vectorstring(file_name,vlines);

      aus << "MMMMM  " <<  aurostd::PaddedPOST("Legitimate VLINES = "+aurostd::utype2string(vlines.size()),40) << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      aus << "MMMMM  " <<  aurostd::PaddedPOST("         maxcheck = "+aurostd::utype2string(maxcheck),40) << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);

      if(aflags.AFLOW_PERFORM_ORDER_SORT) {   // SORT do something
        aus << "MMMMM  Requested SORT [aflags.AFLOW_PERFORM_ORDER_SORT=1] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
        aurostd::sort(vlines);}
      if(aflags.AFLOW_PERFORM_ORDER_REVERSE) { // REVERSE do something
        aus << "MMMMM  Requested REVERSE_SORT [aflags.AFLOW_PERFORM_ORDER_REVERSE=1] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
        aurostd::rsort(vlines);}
      if(aflags.AFLOW_PERFORM_ORDER_RANDOM) { // RANDOM do something
        aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] start (XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(aus,XHOST.QUIET);
        aurostd::random_shuffle(vlines);  // uses the std library but the seed is initialized in xrandom too
        aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] stop (XHOST.AFLOW_RUNXnumber=" << XHOST.AFLOW_RUNXnumber << ") - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(aus,XHOST.QUIET);
      }

      for(uint i=0;(i<vlines.size() && vruns.size()<XHOST_AFLOW_RUNXnumber_multiplier*maxcheck && vruns.size()<VRUNS_MAX_CUTOFF);i++) {  // XHOST_AFLOW_RUNXnumber_multiplier times more... for safety
        if(KBIN::Legitimate_aflowin(vlines.at(i),FALSE,aus)) vruns.push_back(vlines.at(i));             // TRUE puts too much verbosity
      }

      aus << "MMMMM  " <<  aurostd::PaddedPOST("Legitimate VRUNS = "+aurostd::utype2string(vruns.size()),40) << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      // aus << "MMMMM  Legitimate VRUNS = " << vruns.size() << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      // cerr << vlines.size() << endl;
    }

    if(aflags.AFLOW_PERFORM_CLEAN && !XHOST.AFLOW_RUNDIRflag && !XHOST.AFLOW_MULTIflag && !XHOST.AFLOW_RUNXflag) {
      XHOST.AFLOW_RUNDIRflag=TRUE; // give something to clean
    }

    if(aflags.AFLOW_PERFORM_CLEAN && !aflags.AFLOW_PERFORM_DIRECTORY) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"to use --clean, you must specify one or more directories",_INPUT_MISSING_);  //CO20200624
    }
    if(aflags.AFLOW_PERFORM_CLEAN && aflags.AFLOW_PERFORM_DIRECTORY) {
      //    cout << "DEBUG CLEAN = " << vruns.size() << endl;
    }

    bool STOP_DEBUG=aurostd::args2flag(argv,"--STOP|--stop");

    // cdebug << "aflags.KBIN_RUN_AFLOWIN=" << aflags.KBIN_RUN_AFLOWIN << endl;
    // cdebug << "aflags.KBIN_GEN_VASP_FROM_AFLOWIN=" << aflags.KBIN_GEN_VASP_FROM_AFLOWIN << endl;
    // cdebug << "aflags.KBIN_GEN_AFLOWIN_FROM_VASP=" << aflags.KBIN_GEN_AFLOWIN_FROM_VASP << endl;

    // ------------------------------------------------------------------------------------------------------------------------------------
    // nothing to run
    if(!XHOST.AFLOW_RUNDIRflag && !XHOST.AFLOW_MULTIflag && !XHOST.AFLOW_RUNXflag && !aflags.AFLOW_PERFORM_CLEAN && !aflags.AFLOW_PERFORM_DIRECTORY) {
      aus << "MMMMM  KBIN option nothing to run [--run , --run=multi, --run=1, --run=N] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //    return 0;
    }

    // ------------------------------------------------------------------------------------------------------------------------------------
    // nothing specified : CHECK IF DAEMON
    if(XHOST.AFLOW_RUNDIRflag) { // check if daemaon aflowd
      string progname=argv.at(0);
      if(aurostd::substring2bool(progname,"aflowd")) {
        XHOST.AFLOW_MULTIflag=TRUE;
        XHOST.AFLOW_RUNXflag=FALSE;
        XHOST.vflag_aflow.flag("LOOP",TRUE); // add automatically
        aus << "MMMMM  AFLOW: running as DAEMON (aflow --kmode --multi --loop): " << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(aus,XHOST.QUIET);
      }
    }

    // ------------------------------------------------------------------------------------------------------------------------------------
    // run specified directory: XHOST.AFLOW_RUNDIRflag
    if(XHOST.AFLOW_RUNDIRflag) {
      //    bool krun=TRUE;
      if(LDEBUG) cerr << soliloquy << " STEP0b" << endl;
      // [OBSOLETE]   vector<string> vDirectory(aurostd::args2vectorstring(argv,"--DIRECTORY|--D|--d","./"));
      vector<string> vDirectory=vruns;
      // fix the RUNS
      for(uint ii=0;ii<vDirectory.size();ii++)
        if(aurostd::substring2bool(vDirectory.at(ii),_AFLOWIN_))
          aurostd::StringSubst(vDirectory.at(ii),_AFLOWIN_,"");

      // if(aflags.AFLOW_PERFORM_ORDER_SORT) {   // SORT do something
      //   aus << "MMMMM  Requested SORT [aflags.AFLOW_PERFORM_ORDER_SORT=1] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //   aurostd::sort(vDirectory);}
      // if(aflags.AFLOW_PERFORM_ORDER_REVERSE) { // REVERSE do something
      //   aus << "MMMMM  Requested REVERSE_SORT [aflags.AFLOW_PERFORM_ORDER_REVERSE=1] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //   aurostd::rsort(vDirectory);}
      // if(aflags.AFLOW_PERFORM_ORDER_RANDOM) { // RANDOM do something
      //   aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
      //   aurostd::random_shuffle(vDirectory);}  // uses the std library but the seed is initialized in xrandom too

      // aurostd::random_shuffle(vDirectory);
      // std::random_shuffle(vDirectory.begin(),vDirectory.end());

      for(uint idir=0;idir<vDirectory.size();idir++) {
        bool krun=TRUE;
        aflags.Directory=vDirectory.at(idir);
        aus << "MMMMM  AFLOW: running " << _AFLOWIN_ << ", directory" << "=\""  << aflags.Directory << "\" - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(aus,XHOST.QUIET);
        // If necessary PERFORM CLEAN
        if(krun && aflags.AFLOW_PERFORM_CLEAN) {
          aflags.Directory=vDirectory.at(idir);
          aurostd::StringSubst(aflags.Directory,"/"+_AFLOWLOCK_,"");
          for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
            aurostd::StringSubst(aflags.Directory,"/OUTCAR.relax1"+XHOST.vext.at(iext),"");    // CLEAN UP A LITTLE
            aurostd::StringSubst(aflags.Directory,"/OUTCAR.relax2"+XHOST.vext.at(iext),"");    // CLEAN UP A LITTLE
            aurostd::StringSubst(aflags.Directory,"/OUTCAR.static"+XHOST.vext.at(iext),"");    // CLEAN UP A LITTLE
            aurostd::StringSubst(aflags.Directory,"/OUTCAR.bands"+XHOST.vext.at(iext),"");     // CLEAN UP A LITTLE
          }
          for(uint iext=1;iext<XHOST.vext.size();iext++) {  // SKIP uncompressed
            aurostd::StringSubst(aflags.Directory,"/EIGENVAL.relax1"+XHOST.vext.at(iext),"");  // CLEAN UP A LITTLE
            aurostd::StringSubst(aflags.Directory,"/EIGENVAL.relax2"+XHOST.vext.at(iext),"");  // CLEAN UP A LITTLE
            aurostd::StringSubst(aflags.Directory,"/EIGENVAL.static"+XHOST.vext.at(iext),"");  // CLEAN UP A LITTLE
            aurostd::StringSubst(aflags.Directory,"/EIGENVAL.bands"+XHOST.vext.at(iext),"");   // CLEAN UP A LITTLE
          }
          aurostd::StringSubst(aflags.Directory,"/OUTCAR","");  // so it is easier to search
          aurostd::StringSubst(aflags.Directory,"/"+_AFLOWIN_,"");  // so it is easier to search
          //  cerr << aflags.Directory << endl;
          KBIN::Clean(aflags);
          krun=FALSE;
        }
        // RUN
        // cerr << "STEP0b" << endl;
        if(krun) {
          if(LDEBUG) cerr << "STEP1b" << endl;  //CO20170622 - should be debug
          ifstream FileCHECK;string FileNameCHECK;
          // check for directory
          FileNameCHECK=aflags.Directory;
          FileCHECK.open(FileNameCHECK.c_str(),std::ios::in);
          FileCHECK.clear();FileCHECK.close();
          if(!FileCHECK) {                                                                        // ******* Directory is non existent
            aus << "EEEEE  DIRECTORY_NOT_FOUND = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
            aurostd::PrintMessageStream(aus,XHOST.QUIET);
          } else {                                                                                // ******* Directory EXISTS
            if(LDEBUG) cerr << soliloquy << " STEP1c" << endl;
            if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_)) {                                               // ******* Directory is locked
              aus << "LLLLL  DIRECTORY_LOCKED ...bzzzz... !MULTI = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
              aurostd::PrintMessageStream(aus,XHOST.QUIET);
            } else {
              if(LDEBUG) cerr << soliloquy << " STEP1d" << endl;
              if(DirectorySkipped(aflags.Directory)) {                                            // ******* Directory is skipped
                aus << "LLLLL  DIRECTORY_SKIPPED ...bzzzz... !MULTI = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                aurostd::PrintMessageStream(aus,XHOST.QUIET);
              } else {        
                if(LDEBUG) cerr << soliloquy << " STEP1e" << endl;
                if(DirectoryAlreadyInDatabase(aflags.Directory,aflags.AFLOW_FORCE_RUN)) {         // ******* Directory is already in the database
                  aus << "LLLLL  DIRECTORY_ALREADY_IN_DATABASE ...bzzzz... !MULTI = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                  aus << "LLLLL  DIRECTORY_ALREADY_IN_DATABASE ... use \"aflow --multi\" to force the calculation of this entry "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                  aurostd::PrintMessageStream(aus,XHOST.QUIET);
                } else {        
                  if(LDEBUG) cerr << soliloquy << " STEP1f" << endl;
                  if(DirectoryUnwritable(aflags.Directory)) {                                     // ******* Directory is unwritable
                    aus << "LLLLL  DIRECTORY_UNWRITABLE ...bzzzz... !MULTI = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                    aurostd::PrintMessageStream(aus,XHOST.QUIET);
                  } else {                                                                        // ******* Directory is ok
                    if(LDEBUG) cerr << soliloquy << " STEP1g" << endl;
                    if(_VERBOSE_) aus << "LLLLL  GOOD ...bzzz... !MULTI = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                    if(_VERBOSE_) aurostd::PrintMessageStream(aus,XHOST.QUIET);
                    if(aflags.KBIN_RUN_AFLOWIN)               FileNameCHECK=aflags.Directory+"/"+_AFLOWIN_;
                    if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)     FileNameCHECK=aflags.Directory+"/"+_AFLOWIN_;
                    if(aflags.KBIN_GEN_AIMS_FROM_AFLOWIN)     FileNameCHECK=aflags.Directory+"/"+_AFLOWIN_; //CO20180409
                    if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP)     FileNameCHECK=aflags.Directory+"/INCAR";
                    FileCHECK.open(FileNameCHECK.c_str(),std::ios::in);
                    FileCHECK.clear();FileCHECK.close();
                    if(!FileCHECK) {                                                                    // ******* _AFLOWIN_ does not exist
                      if(LDEBUG) cerr << soliloquy << " STEP1h" << endl;
                      if(aflags.KBIN_RUN_AFLOWIN)               aus << "EEEEE  " << _AFLOWIN_ << " not found  = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                      if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)     aus << "EEEEE  " << _AFLOWIN_ << " not found  = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                      if(aflags.KBIN_GEN_AIMS_FROM_AFLOWIN)     aus << "EEEEE  " << _AFLOWIN_ << " not found  = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl; //CO20180409
                      if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP)     aus << "EEEEE  INCAR not found     = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                      aurostd::PrintMessageStream(aus,XHOST.QUIET);
                    } else {                                                                            // ******* _AFLOWIN_ exists RUN
                      if(LDEBUG) cerr << soliloquy << " STEP1i" << endl;
                      if(aflags.KBIN_RUN_AFLOWIN)  {
                        KBIN::RUN_Directory(aflags);
                      }
                      // if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)    KBIN::RUN_Directory(aflags);
                      if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP)     KBIN::GenerateAflowinFromVASPDirectory(aflags);
                      aus << "MMMMM  AFLOW: Run Done " << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                      aurostd::PrintMessageStream(aus,XHOST.QUIET);
                    }
                  } // DIRECTORY WRITABLE
                } // DIRECTORY NOT IN THE DATABASE
              } // DIRECTORY UN SKIPPED
            } // DIRECTORY UN LOCKED
          } // DIRECTORY FOUND
          aus << "MMMMM  AFLOW: Done " << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(aus,XHOST.QUIET);
        }
      } // idir
      // aus << "MMMMM  AFLOW: Done " << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      // aurostd::PrintMessageStream(aus,XHOST.QUIET);
      if(LDEBUG) cerr << soliloquy << " STEP2" << endl;
    }

    // ERRORS ------------------------------------------------------------------------------------------------

    // ------------------------------------------------------------------------------------------------------------------------------------
    // run MULTI and XHOST.AFLOW_RUNXflag (in XHOST.AFLOW_RUNXflag, runs only XHOST.AFLOW_RUNXnumber and then dies)
    // MULTI with SINGLE AND MULTI THREAD VERSION -----------------------------------------------------------------------------------------
    if(XHOST.AFLOW_MULTIflag || XHOST.AFLOW_RUNXflag) {
      uint RUN_times=0;
      bool MULTI_DEBUG=FALSE;
      if(MULTI_DEBUG) {aus << "MMMMM  AFLOW_PTHREADS::FLAG=" << AFLOW_PTHREADS::FLAG << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
      if(MULTI_DEBUG) {aus << "MMMMM  AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
      bool EMPTY=FALSE,FOUND=FALSE;
      bool free_thread;
      int ithread=0;
      if(AFLOW_PTHREADS::FLAG && AFLOW_PTHREADS::MAX_PTHREADS<=1) {
        aus << "EEEEE  ERROR PTHREADS" << endl;
        aus << "MMMMM  AFLOW_PTHREADS::FLAG=" << AFLOW_PTHREADS::FLAG << endl;
        aus << "MMMMM  AFLOW_PTHREADS::MAX_PTHREADS=" << AFLOW_PTHREADS::MAX_PTHREADS << endl;
        aurostd::PrintMessageStream(aus,XHOST.QUIET);
        return 1; //CO20200624 - previously exit
      }
      if(AFLOW_PTHREADS::MAX_PTHREADS>1) {
        aus << "MMMMM  AFLOW: MULTI THREAD START: phread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(aus,XHOST.QUIET);
      }
      aus << "MMMMM  AFLOW: searching subdirectories [d1] " << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
      while(!EMPTY) {
        vector<string> vaflowin;
        //   cerr << "aflags.AFLOW_PERFORM_DIRECTORY=" << aflags.AFLOW_PERFORM_DIRECTORY << endl;
        // cerr << "aflags.AFLOW_PERFORM_FILE=" << aflags.AFLOW_PERFORM_FILE << endl;

        if(aflags.AFLOW_PERFORM_FILE==FALSE) { // NO FILE SPECIFIED = standard  RUN DIRECTORY
          aus << "MMMMM  AFLOW: aflags.AFLOW_PERFORM_FILE==FALSE" << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);

          string FileNameSUBDIR=aurostd::TmpFileCreate("RUN");
          stringstream strstream;
          aus << "rm -f " <<  FileNameSUBDIR << endl;
          aurostd::execute(aus);  // RESET  // RESET

          bool isPRIORITY=FALSE;
          // NEW, the sorting is done internally (speed and reliability)
          for(int ifind=0;ifind<(int)vruns.size();ifind++) {
            isPRIORITY=(aurostd::substring2bool(vruns.at(ifind),"PRIORITY") || aurostd::substring2bool(vruns.at(ifind),"priority"));
            if((isPRIORITY && aurostd::uniform(1.0)<=PRIORITY_PROBABILITY) || !isPRIORITY) {
              aus << "find " << vruns.at(ifind) << " " << XHOST.Find_Parameters;
              if(aflags.KBIN_RUN_AFLOWIN) aus << " -name \"" << _AFLOWIN_ << "\" ";
              // if(aflags.KBIN_RUN_AFLOWIN) aus << " -name \"" << _AFLOWIN_ << "\" ";
              // if(aflags.KBIN_RUN_AFLOWIN) aus << " -name \"" << _AFLOWIN_ << "\" ";
              // if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)  aus << " -name \"" << _AFLOWIN_ << "\" ";   // IS THIS CORRECT ?? CHECK !!!
              if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP)     aus << " -name \"INCAR\" ";
              aus << " | " << PRIORITY_GREP_STRING << " >> ";
              aus << FileNameSUBDIR << endl;
            }
          }
          // perform the command
          aurostd::execute(aus);  // RESET  // RESET

          if(_VERBOSE_) aus << "MMMMM  AFLOW: searching subdirectories [d2] " << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          if(_VERBOSE_) aurostd::PrintMessageStream(aus,XHOST.QUIET);
          aurostd::string2tokens(aurostd::file2string(FileNameSUBDIR),vaflowin,"\n");
          for(i=0;i<(int) vaflowin.size();i++) aurostd::StringSubst(vaflowin.at(i),_AFLOWIN_,"");
          for(i=0;i<(int) vaflowin.size();i++) aurostd::StringSubst(vaflowin.at(i),"INCAR","");
          if(STOP_DEBUG) for(i=0;i<(int) vaflowin.size();i++) cout << vaflowin.at(i) << endl;
          // RANDOMIZING priority
          if(vaflowin.size()>1) { // only if I can poll
            // if(aurostd::substring2bool(vaflowin.at(0),"PRIORITY") || aurostd::substring2bool(vaflowin.at(0),"priority")) aurostd::random_shuffle(vaflowin);
          }
          // loaded up
          aus << "rm -f " << FileNameSUBDIR << endl;
          aurostd::execute(aus);  // RESET  // RESET
        }

        // FILE SPECIFIED
        if(aflags.AFLOW_PERFORM_FILE) {
          aus << "MMMMM  AFLOW: aflags.AFLOW_PERFORM_FILE==TRUE" << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
          vaflowin.clear();
          for(uint i=0;(i<vruns.size() && vaflowin.size()<maxcheck);i++) {
            if(KBIN::Legitimate_aflowin(vruns.at(i),FALSE,aus)) vaflowin.push_back(vruns.at(i)); // TRUE puts too much verbosity
            // vaflowin.push_back(vruns.at(i)); // just load them up... they were checked before //OLD MUST RECHECH THEM as things change on the fly
          }
        }
        // NOW TIME OF SORTING/RANDOMIZING
        if(aflags.AFLOW_PERFORM_ORDER_SORT) {   // SORT do something
          aus << "MMMMM  Requested SORT [aflags.AFLOW_PERFORM_ORDER_SORT=1] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
          aurostd::sort(vaflowin);}
        if(aflags.AFLOW_PERFORM_ORDER_REVERSE) { // REVERSE do something
          aus << "MMMMM  Requested REVERSE_SORT [aflags.AFLOW_PERFORM_ORDER_REVERSE=1] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
          aurostd::sort(vaflowin); //PN modification 
          aurostd::rsort(vaflowin);}
        if(aflags.AFLOW_PERFORM_ORDER_RANDOM) { // RANDOM do something
          aus << "MMMMM  Requested RANDOM [aflags.AFLOW_PERFORM_ORDER_RANDOM=1] - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
          aurostd::random_shuffle(vaflowin);}  // uses the std library but the seed is initialized in xrandom too
        aus << "MMMMM  " <<  aurostd::PaddedPOST("Legitimate VAFLOWIN = "+aurostd::utype2string(vaflowin.size()),40) << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);
        //      aus << "MMMMM  Legitimate VAFLOWIN = " << vaflowin.size() << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);

        //	for(uint i=0;i<vaflowin.size();i++) cerr << i << " " << vaflowin.at(i) << endl;

        // clean AFLOWIN // SAFETY
        for(uint i=0;i<vaflowin.size();i++) aurostd::StringSubst(vaflowin.at(i),_AFLOWIN_,"");  

        if(MULTI_DEBUG) {aus << "MMMMM  SIZE vaflowin=" << vaflowin.size() << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}
        if(MULTI_DEBUG) {for(uint i=0;i<vaflowin.size();i++) {aus << "MMMMM  vaflowin.at(i)=" << vaflowin.at(i) << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;aurostd::PrintMessageStream(aus,XHOST.QUIET);}}
        // cerr << vaflowin.at(10) << " " << vaflowin.at(20) << " " << vaflowin.at(30) << endl;
        // if(aurostd::substring2bool(vaflowin.at(0),"PRIORITY") || aurostd::substring2bool(vaflowin.at(0),"priority")) aurostd::random_shuffle(vaflowin);

        FOUND=FALSE;
        for(uint i=0;i<vaflowin.size()&& !FOUND;i++) {
          aflags.Directory="NULL";
          if(aurostd::DirectoryLocked(vaflowin.at(i),_AFLOWLOCK_)) {
            if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  LOCKED ...bzzz... MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
            if(_VERBOSE_ || STOP_DEBUG) aurostd::PrintMessageStream(aus,XHOST.QUIET);
            FOUND=FALSE;
          } else {
            if(DirectorySkipped(vaflowin.at(i))) {
              if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  SKIPPED ...bzzz... MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
              if(_VERBOSE_ || STOP_DEBUG) aurostd::PrintMessageStream(aus,XHOST.QUIET);
              FOUND=FALSE;
            } else {
              if(DirectoryAlreadyInDatabase(vaflowin.at(i),aflags.AFLOW_FORCE_RUN)) {
                if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  DIRECTORY_ALREADY_IN_DATABASE ...bzzz... MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
                if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  DIRECTORY_ALREADY_IN_DATABASE ... use \"aflow --multi\" to force the calculation of this entry...  MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
                if(_VERBOSE_ || STOP_DEBUG) aurostd::PrintMessageStream(aus,XHOST.QUIET);
                FOUND=FALSE;
              } else {
                if(DirectoryUnwritable(vaflowin.at(i))) {
                  if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  UNWRITABLE ...bzzz... MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
                  if(_VERBOSE_ || STOP_DEBUG) aurostd::PrintMessageStream(aus,XHOST.QUIET);
                  FOUND=FALSE;
                } else {
                  if(_VERBOSE_ || STOP_DEBUG) aus << "LLLLL  GOOD ...bzzz... MULTI "  << vaflowin.at(i) << " " << XHOST.hostname << " " << aflow_get_time_string() << endl;
                  if(_VERBOSE_ || STOP_DEBUG) aurostd::PrintMessageStream(aus,XHOST.QUIET);
                  aflags.Directory=vaflowin.at(i);
                  FOUND=TRUE;
                } // DIRECTORY WRITABLE
              } // DIRECTORY NOT IN THE DATABASE
            } // DIRECTORY UN SKIPPED
          } // DIRECTORY UN LOCKED
        } // DIRECTORY FOUND
        // exiting if STOP_DEBUG
        if(STOP_DEBUG) {cout << "aflow_kbin.cpp: STOP_DEBUG" << endl; return 1;}  //CO20200624 - previously exit

        // FOUND SOMETHING
        if(FOUND==FALSE) {
          EMPTY=TRUE;
          if(AFLOW_PTHREADS::FLAG) {
            aus << "MMMMM  AFLOW: MULTI-THREADED: FLUSHING PTHREADS - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
            aurostd::PrintMessageStream(aus,XHOST.QUIET);
            for(ithread=0;ithread<AFLOW_PTHREADS::MAX_PTHREADS;ithread++)
              if(AFLOW_PTHREADS::vpthread_busy[ithread]) {
                aus << "MMMMM  AFLOW: MULTI-THREADED: Flushing   pthread=" << ithread << "   pthread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                aurostd::PrintMessageStream(aus,XHOST.QUIET);
                pthread_join(AFLOW_PTHREADS::vpthread[ithread],NULL);
              }
          }
        }
        // again another check for LOCK, because NFS (network file system might be slow in concurrent seaches
        //     if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_) && FOUND) {cerr << "AFLOW EXCEPTION on concurrent LOCK: " << aflags.Directory << endl; FOUND=FALSE;}
        if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_) && FOUND) {
          aus << "AFLOW EXCEPTION on concurrent LOCK: " << aflags.Directory << endl;
          aurostd::PrintMessageStream(aus,XHOST.QUIET);
          FOUND=FALSE;
        }

        //  ---------------------------------------------------------------------------- START RUNNING
        if(FOUND) {
          // again another check for LOCK, because NFS (network file system might be slow in concurrent seaches
          EMPTY=FALSE;
          //  ---------------------------------------------------------------------------- KIN_RUN_AFLOWIN
          if(aflags.KBIN_RUN_AFLOWIN) {
            //	  bool PHONONS;
            //  -------------------------------------------------------------------------- KIN_RUN_AFLOWIN multithreaded
            // 	  bool found;
            // 	  for(uint ii=0;ii<qaflags.size()&&!found;ii++)
            // 	    found=(qaflags[ii].Directory==aflags.Directory);       // look in all the list of operations
            // 	  if(found==FALSE) {                                 // new operation, generate and save it
            // 	    qaflags.push_back(aflags);
            // 	  }
            //	  cerr << qaflags.size() << endl;
            if(AFLOW_PTHREADS::FLAG) {
              // there is something to run in aflags.
              // wait and put in ithread there is the number of the thread
              free_thread=AFLOW_PTHREADS::Wait_Available_Free_Threads(ithread,_VERBOSE_);        // WAIT A WHILE !!
              if(free_thread) {
                aus << "MMMMM  AFLOW: Found subdirectory to run " << aflags.Directory<< " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                aus << "MMMMM  AFLOW: MULTI-THREADED: Starting    pthread_free=" << ithread << "   pthread_max=" << AFLOW_PTHREADS::MAX_PTHREADS << " - " << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                aurostd::PrintMessageStream(aus,XHOST.QUIET);
                aflags.AFLOW_PTHREADS_NUMBER=ithread;
                KBIN::RUN_Directory_PTHREADS(aflags);
                RUN_times++;
                if(XHOST.AFLOW_RUNXflag) {
                  aus << "MMMMM  AFLOW: RUNFX finished running " << RUN_times <<  "/" << XHOST.AFLOW_RUNXnumber << " " << aflags.Directory<< " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                  aurostd::PrintMessageStream(aus,XHOST.QUIET);}
                if(XHOST.AFLOW_RUNXflag && RUN_times==XHOST.AFLOW_RUNXnumber) EMPTY=TRUE; // force to end if RUXN reached
              }
            }
            //  -------------------------------------------------------------------------- KIN_RUN_AFLOWIN normal
            if(!AFLOW_PTHREADS::FLAG) {
              KBIN::RUN_Directory(aflags);
              RUN_times++;
              if(XHOST.AFLOW_RUNXflag) {
                aus << "MMMMM  AFLOW: RUNFX finished running " << RUN_times <<  "/" << XHOST.AFLOW_RUNXnumber << " " << aflags.Directory<< " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
                aurostd::PrintMessageStream(aus,XHOST.QUIET);}
              if(XHOST.AFLOW_RUNXflag && RUN_times==XHOST.AFLOW_RUNXnumber) EMPTY=TRUE; // force to end if RUXN reached
            }
          }
          //  ---------------------------------------------------------------------------- KBIN_GEN_VASP_FROM_AFLOWIN normal
          //	if(aflags.KBIN_GEN_VASP_FROM_AFLOWIN)     KBIN::RUN_Directory(argv,aflags);  // IS THIS CORRECT ?? CHECK !!!
          //  ---------------------------------------------------------------------------- KBIN_GEN_AFLOWIN_FROM_VASP normal
          if(aflags.KBIN_GEN_AFLOWIN_FROM_VASP) {
            KBIN::GenerateAflowinFromVASPDirectory(aflags);
          }
        }
        if(XHOST.vflag_aflow.flag("LOOP") && EMPTY && XHOST.AFLOW_RUNXflag==FALSE) {
          EMPTY=FALSE;
          aus << "MMMMM  AFLOW: waiting for new subdirectories: " << (int) _KBIN_LOOP_SLEEP_/60 << "mins ";
          aus << " - " << XHOST.hostname << " - " << aflow_get_time_string() << endl;//endl;
          aurostd::PrintMessageStream(aus,XHOST.QUIET);
          aurostd::Sleep(_KBIN_LOOP_SLEEP_);
        }
      }
      aus << "MMMMM  AFLOW: no more subdirectories to run " << " - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
    }

    // ------------------------------------------------------------------------------------------------------------------------------------

    return 1;
  }
} // namespace end of MAIN


// ***************************************************************************
// KBIN::MPI_Extract
// ***************************************************************************
// This function extracts from _AFLOWIN_ the parameters for MPI run
namespace KBIN {
  void MPI_Extract(string AflowIn,ofstream &FileMESSAGE,_aflags &aflags,_kflags &kflags) {
    ostringstream aus;
    bool Kmpi=TRUE;
    kflags.KBIN_MPI_NCPUS=0;
    aus << "00000  [AFLOW_MODE_MPI] found in " << _AFLOWIN_ << " " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
    aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);
    // get (integer) kflags.KBIN_MPI_NCPUS

    if(aflags.AFLOW_GLOBAL_NCPUS<1) {
      if(Kmpi && !aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=",TRUE) && !aflags.AFLOW_MACHINE_LOCAL.flag()) {                      // DEFAULT NO CPU SPECIFIED
        kflags.KBIN_MPI_NCPUS=MPI_NCPUS_DEFAULT;
        aus << "00000  MESSAGE MPI: NCPUS=NNNN is missing, taking NCPUS=" << kflags.KBIN_MPI_NCPUS << "  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
        Kmpi=FALSE;
      }
      if(Kmpi && (aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=MAX",TRUE) || aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=AUTO",TRUE) || aflags.AFLOW_MACHINE_LOCAL.flag())) { // DEFAULT NCPUS=MAX
        kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;
        if (aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=MAX",TRUE)) kflags.KBIN_MPI_NCPUS_STRING = "MAX";
        if (aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=AUTO",TRUE)) kflags.KBIN_MPI_NCPUS_STRING = "AUTO";

        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_MPICH")) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;	//CO
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_BETA_OPENMPI")) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;        // with DUKE_BETA force NCPUS from QUEUE
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_QRATS_MPICH")) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN; //CO
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_MATERIALS")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;   // with DUKE_MATERIALS force NCPUS
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_AFLOWLIB")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;   // with DUKE_AFLOWLIB force NCPUS
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::DUKE_HABANA")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;   // with DUKE_HABANA force NCPUS
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_EOS")) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_NTASKS; // [OBSOLETE] XHOST.SLURM_CPUS_ON_NODE; no CPUS because it gets fooled by HT
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_DRACO")) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_NTASKS; // [OBSOLETE] XHOST.SLURM_CPUS_ON_NODE; no CPUS because it gets fooled by HT
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_COBRA")) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_NTASKS; // [OBSOLETE] XHOST.SLURM_CPUS_ON_NODE; no CPUS because it gets fooled by HT
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MPCDF_HYDRA")) kflags.KBIN_MPI_NCPUS=XHOST.SLURM_NTASKS; // [OBSOLETE] XHOST.SLURM_CPUS_ON_NODE; no CPUS because it gets fooled by HT
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE001")) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;   // with MACHINE001; DX added 20190509
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::MACHINE002")) kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;   // with MACHINE002; DX added 20190509
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::CMU_EULER"))  kflags.KBIN_MPI_NCPUS=XHOST.PBS_NUM_PPN;;  //DX20190107 - CMU EULER // with CMU_EULER force NCPUS //DX20181113
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::OHAD")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;           // MACHINE2 has only NCPUS //CO20181113
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::HOST1")) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores;          // MACHINE1 has only NCPUS //CO20181113
        if(aflags.AFLOW_MACHINE_LOCAL.flag("MACHINE::FULTON_MARYLOU"))  kflags.KBIN_MPI_NCPUS=XHOST.SLURM_NTASKS;        // with FULTON_MARYLOU force NCPUS
        if(kflags.KBIN_MPI_NCPUS<1) kflags.KBIN_MPI_NCPUS=XHOST.CPU_Cores; // SAFE

        aus << "00000  MESSAGE MPI: found NCPUS=MAX  NCPUS="<<kflags.KBIN_MPI_NCPUS<<" " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
        Kmpi=FALSE;
      }
      if(Kmpi && aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=",TRUE)) {                     // DEFAULT NCPUS=XXX
        kflags.KBIN_MPI_NCPUS=aurostd::substring2utype<int>(AflowIn,"[AFLOW_MODE_MPI_MODE]NCPUS=",TRUE);
        if(kflags.KBIN_MPI_NCPUS>0) {
          kflags.KBIN_MPI_NCPUS_STRING = aurostd::utype2string<int>(kflags.KBIN_MPI_NCPUS); //ME20181113
          aus << "00000  MESSAGE MPI: found NCPUS="<<kflags.KBIN_MPI_NCPUS<<" " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
        }
        Kmpi=FALSE;
      }
    } else {
      kflags.KBIN_MPI_NCPUS=aflags.AFLOW_GLOBAL_NCPUS;
      aus << "00000  MESSAGE MPI: NCPUS is overriden, taking NCPUS=" << kflags.KBIN_MPI_NCPUS << "  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);  
    }

    if(kflags.KBIN_MPI_NCPUS<1) kflags.KBIN_MPI_NCPUS=1;                                              // DEFAULT NCPUS=troubles

    if(kflags.KBIN_MPI_NCPUS==1) {
      kflags.KBIN_MPI=FALSE;
      aus << "00000  MESSAGE MPI: found NCPUS=1 " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aus << "00000  MESSAGE MPI: going back to SERIAL execution " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
    }
    if(kflags.KBIN_MPI) {
      // get (string) kflags.KBIN_MPI_START
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]START=",TRUE)) {
        kflags.KBIN_MPI_START=MPI_START_DEFAULT;
        aus << "00000  MESSAGE MPI: START string is missing, taking START=\"" << kflags.KBIN_MPI_START << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
        kflags.KBIN_MPI_START=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_MPI_MODE]START=",TRUE),'"');
        aus << "00000  MESSAGE MPI: found START=\"" << kflags.KBIN_MPI_START << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      // get (string) kflags.KBIN_MPI_STOP
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]STOP=",TRUE)) {
        kflags.KBIN_MPI_STOP=MPI_STOP_DEFAULT;
        aus << "00000  MESSAGE MPI: STOP string is missing, taking STOP=\"" << kflags.KBIN_MPI_STOP << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
        kflags.KBIN_MPI_STOP=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_MPI_MODE]STOP=",TRUE),'"');
        aus << "00000  MESSAGE MPI: found STOP=\"" << kflags.KBIN_MPI_STOP << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      // get (string) kflags.KBIN_MPI_COMMAND
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]COMMAND=",TRUE)) {
        kflags.KBIN_MPI_COMMAND=MPI_COMMAND_DEFAULT;
        aus << "00000  MESSAGE MPI: COMMAND string is missing, taking COMMAND=\"" << kflags.KBIN_MPI_COMMAND << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
        kflags.KBIN_MPI_COMMAND=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_MPI_MODE]COMMAND=",TRUE),'"');
        aus << "00000  MESSAGE MPI: found COMMAND=\"" << kflags.KBIN_MPI_COMMAND << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      kflags.KBIN_MPI_AUTOTUNE=aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]AUTOTUNE",TRUE);
      if(kflags.KBIN_MPI_AUTOTUNE) {
        aus << "00000  MESSAGE MPI: found AUTOTUNE option " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "00000  MESSAGE MPI: input files WILL be auto-tuned for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
        aus << "00000  MESSAGE MPI: AUTOTUNE option NOT found " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aus << "00000  MESSAGE MPI: input files MUST be appropriate for PARALLEL execution with " << kflags.KBIN_MPI_NCPUS << " CPUs "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      // get (string) kflags.KBIN_MPI_BIN
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]BINARY=",TRUE)) {
        kflags.KBIN_MPI_BIN=DEFAULT_VASP_MPI_BIN;
        aus << "00000  MESSAGE MPI: BINARY string is missing, taking BIN=\"" << kflags.KBIN_MPI_BIN << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
        kflags.KBIN_MPI_BIN=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_MPI_MODE]BINARY=",TRUE),'"');
        aus << "00000  MESSAGE MPI: found BINARY=\"" << kflags.KBIN_MPI_BIN << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
      //ME20190107 - Grab the serial binary to propagate into child aflow.in files
      if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY]")) {
        kflags.KBIN_SERIAL_BIN = aurostd::substring2string(AflowIn, "[AFLOW_MODE_BINARY]");
      } else if (aurostd::substring2bool(AflowIn, "[AFLOW_MODE_BINARY=")) {
        kflags.KBIN_SERIAL_BIN = aurostd::RemoveCharacter(aurostd::substring2string(AflowIn, "[AFLOW_MODE_BINARY="), ']');
      }
      aus << "00000  MESSAGE MPI: Overriding BINARY=\"" << kflags.KBIN_BIN << "\" to BINARY =\"" << kflags.KBIN_MPI_BIN << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      kflags.KBIN_BIN=kflags.KBIN_MPI_BIN;

      // get (string) kflags.KBIN_MPI_OPTIONS
      if(!aurostd::substring2bool(AflowIn,"[AFLOW_MODE_MPI_MODE]OPTIONS=",TRUE)) {
        kflags.KBIN_MPI_OPTIONS=VASP_OPTIONS_MPI_DEFAULT;
        aus << "00000  MESSAGE MPI: OPTIONS string is missing, taking OPTIONS=\"" << kflags.KBIN_MPI_OPTIONS << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      } else {
        kflags.KBIN_MPI_OPTIONS=aurostd::RemoveCharacter(aurostd::substring2string(AflowIn,"[AFLOW_MODE_MPI_MODE]OPTIONS=",TRUE),'"');
        aus << "00000  MESSAGE MPI: found OPTIONS=\"" << kflags.KBIN_MPI_OPTIONS << "\"  " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
        aurostd::PrintMessageStream(FileMESSAGE,aus,XHOST.QUIET);    
      }
    }
  }
} // namespace KBIN


// ***************************************************************************
// KBIN::StartStopCheck
// ***************************************************************************
namespace KBIN {
  void StartStopCheck(const string &AflowIn,string str1,string str2,bool &flag,bool &flagS) {
    flag =
      aurostd::substring2bool(AflowIn,str1) || aurostd::substring2bool(AflowIn,str2) ;
    flagS= (aurostd::substring2bool(AflowIn,str1+"START") && aurostd::substring2bool(AflowIn,str1+"STOP")) ||
      (aurostd::substring2bool(AflowIn,str2+"_START") && aurostd::substring2bool(AflowIn,str2+"_STOP"));
    if(flagS) flag=FALSE;
  }
}

namespace KBIN {
  void StartStopCheck(const string &AflowIn,string str1,bool &flag,bool &flagS) {
    flag = aurostd::substring2bool(AflowIn,str1);
    flagS= aurostd::substring2bool(AflowIn,str1+"START") && aurostd::substring2bool(AflowIn,str1+"STOP");
    if(flagS) flag=FALSE;
  }
}

// ***************************************************************************
// KBIN::RUN_Directory
// ***************************************************************************
namespace KBIN {
  void RUN_Directory(_aflags& aflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"KBIN::RUN_Directory():";
    ostringstream aus;
    ifstream FileSUBDIR;string FileNameSUBDIR;
    FileNameSUBDIR=aflags.Directory;
    FileSUBDIR.open(FileNameSUBDIR.c_str(),std::ios::in);
    FileSUBDIR.clear();FileSUBDIR.close();
    // string::size_type sub_size1,sub_size2;
    string AflowIn,AflowInMode,subS,subS1,subS2;
    string::iterator pos;
    bool Krun=TRUE;
    //  int i;

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    if(aflags.Directory.at(0)!='/' && aflags.Directory.at(0)!='.' && aflags.Directory.at(0)!=' ') aflags.Directory="./"+aflags.Directory;

    if(!FileSUBDIR) {                                                                                           // ******* Directory is non existent
      aus << "EEEEE  DIRECTORY_NOT_FOUND = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
      aurostd::PrintMessageStream(aus,XHOST.QUIET);
    } else {                                                                                                    // ******* Directory EXISTS
      // ***************************************************************************
      // Check LOCK again
      if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_) || DirectorySkipped(aflags.Directory) || DirectoryAlreadyInDatabase(aflags.Directory,aflags.AFLOW_FORCE_RUN) || DirectoryUnwritable(aflags.Directory)) {
        // ******* Directory is locked/skipped/unwritable
        // LOCK/SKIP/UNWRITABLE exist, then RUN already RUN
        if(aurostd::DirectoryLocked(aflags.Directory,_AFLOWLOCK_)) {
          aus << "LLLLL  LOCKED ... bzzz ... KBIN::RUN_Directory "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "LLLLL  LOCKED ... Probably other aflows are concurring with this. KBIN::RUN_Directory " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(aus,XHOST.QUIET);
        }
        if(DirectorySkipped(aflags.Directory)) {
          aus << "LLLLL  SKIPPED ... bzzz ... KBIN::RUN_Directory "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "LLLLL  SKIPPED ... Probably other aflows are concurring with this. KBIN::RUN_Directory " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(aus,XHOST.QUIET);
        }
        if(DirectoryAlreadyInDatabase(aflags.Directory,aflags.AFLOW_FORCE_RUN)) {
          aus << "LLLLL  ALREADY_IN_DATABASE ... bzzz ... KBIN::RUN_Directory "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "LLLLL  ALREADY_IN_DATABASE ... use \"aflow --multi\" to force the calculation of this entry "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "LLLLL  ALREADY_IN_DATABASE ... Probably other aflows are concurring with this. KBIN::RUN_Directory " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(aus,XHOST.QUIET);
        }
        if(DirectoryUnwritable(aflags.Directory)) {
          aus << "LLLLL  UNWRITABLE ... bzzz ... KBIN::RUN_Directory "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus << "LLLLL  UNWRITABLE ... Probably other aflows are concurring with this. KBIN::RUN_Directory " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(aus,XHOST.QUIET);
        }
      } else {                                                                                                  // ******* Directory is fine
        // make a dumb lock as soon as possible -------------------------------------
        aus.clear();aus.str(std::string());
        aus << "echo \"NNNNN  KBIN LOCK ASAP for NFS concurrent jobs (aflow" << string(AFLOW_VERSION) << ")\" >> " << aflags.Directory+"/"+_AFLOWLOCK_ << endl;
        // aus << "/home/auro/bin/aflow -machine >> " << aflags.Directory+"/"+_AFLOWLOCK_ << endl;
        // aus << XHOST.command("sensors") << " >> " << aflags.Directory+"/"+_AFLOWLOCK_ << endl;
        aurostd::execute(aus);
        // now change its permission
        aurostd::ChmodFile("664",string(aflags.Directory+"/"+_AFLOWLOCK_));
        // now the lock should be done ----------------------------------------------
        ifstream FileAFLOWIN;string FileNameAFLOWIN;
        FileNameAFLOWIN=aflags.Directory+"/"+_AFLOWIN_;
        FileAFLOWIN.open(FileNameAFLOWIN.c_str(),std::ios::in);
        if(!FileAFLOWIN) {                                                                                      // ******* _AFLOWIN_ does not exist
          aus << "EEEEE  " << _AFLOWIN_ << " ABSENT   = "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(aus,XHOST.QUIET);
        } else {                                                                                                // ******* _AFLOWIN_ exists RUN    
          // ***************************************************************************
          // RESET LOCK
          ofstream FileLOCK;
          string FileNameLOCK=aflags.Directory+"/"+_AFLOWLOCK_;
          //	FileLOCK.open(FileNameLOCK.c_str(),std::ios::out);
          FileLOCK.open(FileNameLOCK.c_str(),std::ios::app);
          // ***************************************************************************
          // WRITE LOCK
          if(0) {
            aus <<    "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-Flow - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
            aus <<    "MMMMM  (C) "<<XHOST.Copyright_Years<<", Stefano Curtarolo - Duke University   - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
            aus <<    "MMMMM  High-Throughput ab-initio Computing - " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
            aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
            // ***************************************************************************
            // WRITE AFLOW VERSION
            // aus << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-Flow " << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
            // aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
          }
          // ***************************************************************************
          // START DIRECTORY
          aus      << "XXXXX  KBIN DIRECTORY BEGIN (aflow" << string(AFLOW_VERSION) << ")  "  << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          //	aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
          aus      << "XXXXX  KBIN XHOST.CPU_Model : "<<  XHOST.CPU_Model << "" << endl;// << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus      << "XXXXX  KBIN XHOST.CPU_Cores : "<<  XHOST.CPU_Cores << "" << endl;// << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus      << "XXXXX  KBIN XHOST.CPU_MHz   : "<<  XHOST.CPU_MHz << "" << endl;// << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aus      << "XXXXX  KBIN XHOST.RAM_GB    : "<<  XHOST.RAM_GB << "" << endl;// << Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_) << endl;
          aurostd::PrintMessageStream(FileLOCK,aus,XHOST.QUIET);
          // ***************************************************************************
          // FLUSH & REOPEN to avoid double writing
          FileLOCK.flush();FileLOCK.close();FileLOCK.open(FileNameLOCK.c_str(),std::ios::app);
          // ***************************************************************************
          // NOW Digest AFLOWIN
          FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
          //DX20190125 [OBSOLETE] - need to remove null bytes : AflowIn.clear();char c; while (FileAFLOWIN.get(c)) AflowIn+=c;               // READ _AFLOWIN_ and put into AflowIn
          AflowIn.clear();char c; while (FileAFLOWIN.get(c)) if(c!='\0'){AflowIn+=c;}               // READ _AFLOWIN_ and put into AflowIn //DX20190125 - remove null bytes from AflowIn
          FileAFLOWIN.clear();FileAFLOWIN.seekg(0);
          AflowIn=aurostd::RemoveComments(AflowIn); // NOW Clean AFLOWIN
          vector<string> vAflowIn;aurostd::string2vectorstring(AflowIn,vAflowIn); //CO20181226

          _kflags kflags=KBIN::VASP_Get_Kflags_from_AflowIN(AflowIn,FileLOCK,aflags); //CO20200624 - made separate function for getting kflags

          if(kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP) {  // [AFLOW_MODE_PRESCRIPT] construction
            aurostd::stringstream2file(kflags.AFLOW_MODE_PRESCRIPT,string(aflags.Directory+"/"+DEFAULT_AFLOW_PRESCRIPT_COMMAND));
            aurostd::ChmodFile("755",string(aflags.Directory+"/"+DEFAULT_AFLOW_PRESCRIPT_COMMAND));
          }
          if(kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP) {  // [AFLOW_MODE_POSTSCRIPT] construction
            aurostd::stringstream2file(kflags.AFLOW_MODE_POSTSCRIPT,string(aflags.Directory+"/"+DEFAULT_AFLOW_POSTSCRIPT_COMMAND));
            aurostd::ChmodFile("755",string(aflags.Directory+"/"+DEFAULT_AFLOW_POSTSCRIPT_COMMAND));
          }
          // ************************************************************************************************************************************
          // ALIEN MODE
          if(Krun && kflags.AFLOW_MODE_ALIEN) {
            // ***************************************************************************
            // ALIEN MODE  // must contain EMAIL perform
            if(Krun) {
              Krun=(Krun && ALIEN::Run_Directory(FileLOCK,aflags,kflags));
            }
            // ***************************************************************************
            // COMPRESS
            if(Krun && kflags.KZIP_COMPRESS) {
              Krun=(Krun && KBIN::CompressDirectory(aflags,kflags));
            }
          }
          // ************************************************************************************************************************************
          // MATLAB MODE
          if(Krun && kflags.AFLOW_MODE_MATLAB) {
            if(Krun) {
              aurostd::CommandRequired(DEFAULT_KBIN_MATLAB_BIN); // MATLAB MUST BE AVAILABLE
              Krun=(Krun && KBIN_MATLAB_Directory(FileLOCK,aflags,kflags));
            }
            // ***************************************************************************
            // COMPRESS
            if(Krun && kflags.KZIP_COMPRESS)
              Krun=(Krun && KBIN::CompressDirectory(aflags,kflags));
            Krun=FALSE;
          }
          // ************************************************************************************************************************************
          // AIMS MODE
          if(Krun && kflags.AFLOW_MODE_AIMS) {
            // ***************************************************************************
            // AIMS MODE  // must contain EMAIL perform
            if(Krun) {
              Krun=(Krun && KBIN::AIMS_Directory(FileLOCK,aflags,kflags));
            }
            // ***************************************************************************
            // COMPRESS
            if(Krun && kflags.KZIP_COMPRESS)
              Krun=(Krun && KBIN::CompressDirectory(aflags,kflags));
          }
          // ************************************************************************************************************************************
          // ************************************************************************************************************************************
          // VASP MODE
          if(Krun && kflags.AFLOW_MODE_VASP) {
            // ***************************************************************************
            // VASP MODE  // must contain EMAIL perform
            if(Krun) {
              Krun=(Krun && KBIN::VASP_Directory(FileLOCK,aflags,kflags));
            }
            // ***************************************************************************
            // COMPRESS	    
            if(Krun && kflags.KZIP_COMPRESS) {
              // cerr << aurostd::execute2string("ls -las "+aflags.Directory) << endl;
              Krun=(Krun && KBIN::CompressDirectory(aflags,kflags));
              // cerr << aurostd::execute2string("ls -las "+aflags.Directory) << endl;
            }
          }
          // ************************************************************************************************************************************
          // MATLAB MODE
          if(Krun && kflags.KBIN_PHONONS_CALCULATION_FROZSL && !kflags.AFLOW_MODE_VASP) {
            // PRESCRIPT
            if(Krun && kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT || kflags.AFLOW_MODE_PRESCRIPT_EXPLICIT_START_STOP)
              KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_PRESCRIPT_COMMAND,DEFAULT_AFLOW_PRESCRIPT_OUT);
            // POSTSCRIPT
            if(Krun && kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT || kflags.AFLOW_MODE_POSTSCRIPT_EXPLICIT_START_STOP)
              KBIN::RUN_DirectoryScript(aflags,DEFAULT_AFLOW_POSTSCRIPT_COMMAND,DEFAULT_AFLOW_POSTSCRIPT_OUT);
          }
          // ***************************************************************************
          // FINALIZE AFLOWIN
          // DONE, turn OFF the flag
          kflags.AFLOW_MODE_VASP=FALSE;    
          // ***************************************************************************
          // NFS cache-cleaning HACK, opendir() and closedir() to invalidate cache
          // https://stackoverflow.com/questions/8311710/nfs-cache-cleaning-command
          //CO20171106
          vector<string> vfiles_dummyls;
          aurostd::DirectoryLS(aflags.Directory,vfiles_dummyls);
          // ***************************************************************************
          // WRITE END
          aurostd::string2file(string(Message(aflags,_AFLOW_MESSAGE_DEFAULTS_,_AFLOW_FILE_NAME_)+"\n"),string(aflags.Directory+"/"+DEFAULT_AFLOW_END_OUT));
          // ***************************************************************************
          // MAKE READEABLE
          aurostd::ChmodFile("664",string(aflags.Directory+"/"+_AFLOWLOCK_));
          aus.clear();aus.str(std::string());
          aus << "cp " << aflags.Directory << "/" << _AFLOWLOCK_ << " " << aflags.Directory << "/" << _AFLOWLOCK_ << "." << string(AFLOW_VERSION) << endl;  
          aurostd::execute(aus);  
          aurostd::ChmodFile("664",string(aflags.Directory+"/*"));
          aurostd::ChmodFile("777",string(aflags.Directory+"/*"));
        }
        FileAFLOWIN.clear();FileAFLOWIN.close();
      }
    }
  };
} // namespace

// *******************************************************************************************
namespace KBIN {
  void AFLOW_RUN_Directory(const _aflags& aflags) {
    aurostd::execute(XHOST.command("aflow")+" --run=1 --DIRECTORY="+aflags.Directory);  // run it OUTSIDE
    // this is cool as if the particula program explodes, there is still aflow running
  }
} // namespace

// *******************************************************************************************
namespace KBIN {
  void RUN_DirectoryScript(const _aflags& aflags,const string& script,const string& output) {        // AFLOW_FUNCTION_IMPLEMENTATION
    ostringstream aus;
    aurostd::StringstreamClean(aus);
    aus << "cd " << aflags.Directory << endl;
    aus <<  "chmod 755 " << script << endl;
    aus << "./" << script << " >> " << output << endl;
    aurostd::execute(aus);
  }
} // namespace

// *******************************************************************************************
namespace KBIN {
  bool CompressDirectory(const _aflags& aflags,const _kflags& kflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy = XPID + "KBIN::CompressDirectory():";
    //DX+CO START
    // [OBSOLETE] aus << "cd " << aflags.Directory << " && ";
    // [OBSOLETE]  aus << "ls | grep -v .EXT | ";                                //CO, skip anything with bzip extension
    // [OBSOLETE]  aus << "grep -v LOCK | grep -v " << _AFLOWLOCK_ << " | ";     //CO never zip LOCK or .LOCK (agl.LOCK) newly defined LOCK
    // [OBSOLETE] aus << "grep -v SKIP | ";
    // [OBSOLETE] aus << "grep -v " << KBIN_SUBDIRECTORIES << " | ";
    // [OBSOLETE]  aus << "grep -v aflow.in | grep -v " << _AFLOWIN_ << " ";;    //CO, never zip aflow.in or _aflow.in (agl_aflow.in) or newly defined aflow.in
    vector<string> _vfiles,vfiles;
    string compressed_variant;
    aurostd::DirectoryLS(aflags.Directory,_vfiles);
    string file_path;
    for(uint i=0;i<_vfiles.size();i++){
      if(aurostd::IsCompressed(_vfiles[i])){continue;}  //doesn't need full path, just substring2bool for zip variants, e.g., .EXT
      if(_vfiles[i].size() && _vfiles[i][0]=='.'){continue;}  //do not try to compress any hidden files, .nfs stuff is particularly problematic to compress
      if(aurostd::substring2bool(_vfiles[i],KBIN_SUBDIRECTORIES)){continue;}
      if(aurostd::substring2bool(_vfiles[i],"LOCK")){continue;}
      if(aurostd::substring2bool(_vfiles[i],_AFLOWLOCK_)){continue;}
      if(aurostd::substring2bool(_vfiles[i],"SKIP")){continue;}
      if(aurostd::substring2bool(_vfiles[i],"aflow.in")){continue;}
      if(aurostd::substring2bool(_vfiles[i],DEFAULT_AFLOW_END_OUT) || aurostd::substring2bool(_vfiles[i],"aflow.end.out")){continue;}  //CO20170613, file is special because it gets written after compression
      if(aurostd::substring2bool(_vfiles[i],_AFLOWIN_)){continue;}
      file_path=aflags.Directory + "/" + _vfiles[i];
      if(LDEBUG) {cerr << soliloquy << " file_path=" << file_path << endl;}
      if(aurostd::IsDirectory(file_path)){continue;}  //compress files only
      // [OBSOLETE]  if(aurostd::EFileExist(file_path,compressed_variant)){ //SC20200408
      // [OBSOLETE]  //need full path here, also, notice the placement here, actual compressed variant would have been skipped, this is for the uncompressed variant //SC20200408
      // [OBSOLETE]  //both compressed and uncompressed variants exists //SC20200408
      // [OBSOLETE]  //assume compressed is from a previous run, hence obsolete //SC20200408
      // [OBSOLETE]  //delete it before compressing uncompressed variant //SC20200408
      // [OBSOLETE]  aurostd::RemoveFile(compressed_variant);  //  //SC20200408 = no needed because it will overwrite it -f  //SC20200408
      // [OBSOLETE] }
      vfiles.push_back(_vfiles[i]);
    }
    // [OBSOLETE] aurostd::string2vectorstring(aurostd::execute2string(aus),vfiles);
    // [OBSOLETE] cerr << vfiles.size() << endl;
    // [OBSOLETE] aurostd::StringstreamClean(aus);
    // [OBSOLETE] cerr << "CO " << aurostd::joinWDelimiter(vfiles," ") << endl;
    if(vfiles.size()){
      ostringstream aus;
      //aurostd::StringstreamClean(aus);
      aus << "cd " << aflags.Directory << " && " << endl;
      for(uint i=0;i<vfiles.size();i++){ //better than doing it all in one shot
        aus << kflags.KZIP_BIN << " -9f " << vfiles[i] << "; " << endl;  // semi-colon is important, keeps going if it stalls on one
      }
      // aus << kflags.KZIP_BIN << " " << aurostd::joinWDelimiter(vfiles," ") << endl; //AVOID, because if one fails, the whole command stops
      aurostd::execute(aus);
      // cerr << aus.str() << endl;
    }
    // [OBSOLETE] aus << kflags.KZIP_BIN << " `find . " << XHOST.Find_Parameters << " -name \"*\" | grep -v LOCK | grep -v SKIP | grep \"./\"` " << endl;
    // [OBSOLETE] aus << kflags.KZIP_BIN << " `ls | grep -v LOCK | grep -v SKIP | grep -v " << KBIN_SUBDIRECTORIES << "| grep -v " << _AFLOWIN_ << " | grep -v apl.xml ` " << endl;
    // [OBSOLETE] apl.xml can now be zipped
    // [OBSOLETE] set up is slightly redundant (LOCK vs. agl.LOCK), but very safe
    // [OBSOLETE] aus << kflags.KZIP_BIN << " `ls | grep -v .EXT | ";                 //CO, skip anything with bzip extension
    // [OBSOLETE] aus << "grep -v LOCK | grep -v " << _AFLOWLOCK_ << " | ";           //CO never zip LOCK or .LOCK (agl.LOCK) newly defined LOCK
    // [OBSOLETE] aus << "grep -v SKIP | ";
    // [OBSOLETE] aus << "grep -v " << KBIN_SUBDIRECTORIES << " | ";
    // [OBSOLETE] aus << "grep -v aflow.in | grep -v " << _AFLOWIN_ << " ` " << endl; //CO, never zip aflow.in or _aflow.in (agl_aflow.in) or newly defined aflow.in
    // [OBSOLETE] cerr << aus.str() << endl;
    // [OBSOLETE] aurostd::execute(aus);
    //DX+CO END
    return TRUE;
  }
}

// *******************************************************************************************
namespace KBIN {
  bool CompressDirectory(const _aflags& aflags) {        // AFLOW_FUNCTION_IMPLEMENTATION
    _kflags kflags;
    kflags.KZIP_BIN=DEFAULT_KZIP_BIN+" -9q";
    return KBIN::CompressDirectory(aflags,kflags);
  }
}

// *******************************************************************************************
// KBIN::Clean
// *******************************************************************************************
namespace KBIN {
  void Clean(const _aflags& aflags) {          // AFLOW_FUNCTION_IMPLEMENTATION
    //    cerr << XPID << "KBIN::Clean: aflags.Directory=" << aflags.Directory << endl;
    KBIN::Clean(aflags.Directory);
  }
}

namespace KBIN {
  void Clean(const string _directory) {        // AFLOW_FUNCTION_IMPLEMENTATION
    string directory=_directory;
    //    cerr << XPID << "KBIN::Clean: directory=" << aflags.Directory << endl;

    aurostd::StringSubst(directory,"/"+_AFLOWIN_,"");  // so it is easier to search

    for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
      aurostd::StringSubst(directory,"/"+XHOST.vext.at(iext),"");  // so it is easier to search
      aurostd::StringSubst(directory,"/"+XHOST.vext.at(iext),"");  // so it is easier to search    
      aurostd::StringSubst(directory,"/"+XHOST.vext.at(iext),"");  // so it is easier to search
    }   

    for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
      aurostd::StringSubst(directory,"/agl_aflow.in"+XHOST.vext.at(iext),"");  // so it is easier to search
      aurostd::StringSubst(directory,"/ael_aflow.in"+XHOST.vext.at(iext),"");  // so it is easier to search    
      aurostd::StringSubst(directory,"/aflow.in"+XHOST.vext.at(iext),"");      // so it is easier to search
    }   
    aurostd::StringSubst(directory,"/agl_aflow.in","");  // so it is easier to search
    aurostd::StringSubst(directory,"/ael_aflow.in","");  // so it is easier to search    
    aurostd::StringSubst(directory,"/aflow.in","");      // so it is easier to search

    if(!aurostd::FileExist(string(directory+"/"+"NOCLEAN")) &&
        !aurostd::FileExist(string(directory+"/"+"NOCLEAR")) &&
        !aurostd::FileExist(string(directory+"/"+"noclean")) &&
        !aurostd::FileExist(string(directory+"/"+"noclear")) &&
        !aurostd::FileExist(string(directory+"/"+"Makefile")) &&
        !aurostd::FileExist(string(directory+"/"+"aflow.h")) &&
        !aurostd::FileExist(string(directory+"/"+"paper.tex")) &&
        !aurostd::FileExist(string(directory+"/"+"manuscript.tex")) &&
        !aurostd::FileExist(string(directory+"/"+"supplementary_information.tex")) &&
        !aurostd::FileExist(string(directory+"/"+"supplementary_materials.tex")) &&
        !aurostd::FileExist(string(directory+"/"+"review.tex"))) {
      if(aurostd::FileExist(string(directory+"/"+_AFLOWIN_)) ||    // normal aflow.in or specified it
          aurostd::FileExist(string(directory+"/aflow.in")) ||      // normal aflow.in
          aurostd::FileExist(string(directory+"/agl_aflow.in")) ||  // normal agl_aflow.in
          aurostd::FileExist(string(directory+"/ael_aflow.in")) ) { // normal ael_aflow.in

        // CLEAN directory
        //DX+CO START
        vector<string> vfiles;  //not only files, includes EVERYTHING
        string file_path;
        aurostd::DirectoryLS(directory,vfiles);
        for(uint i=0;i<vfiles.size();i++) {
          file_path=directory + "/" + vfiles.at(i);
          if(aurostd::substring2bool(vfiles.at(i),_AFLOWIN_)){continue;}
          if(aurostd::substring2bool(vfiles.at(i),"aflow.in")){continue;}
          if(aurostd::substring2bool(vfiles.at(i),"agl_aflow.in")){continue;}
          if(aurostd::substring2bool(vfiles.at(i),"ael_aflow.in")){continue;}
          if(aurostd::substring2bool(vfiles.at(i),DEFAULT_AFLOW_FROZSL_INPUT_OUT)){continue;}
          if(aurostd::IsDirectory(file_path)){                 
            if(aurostd::substring2bool(vfiles.at(i),KBIN_SUBDIRECTORIES)){ // only directories we don't ignore
              aurostd::RemoveDirectory(file_path);
            }
            continue;                                                   // ignore all other directories
          }
          if(aurostd::IsFile(file_path)){                 
            if(vfiles.at(i).size() && vfiles.at(i)[0]=='.'){
              if(aurostd::substring2bool(vfiles.at(i),".nfs")){            // only hidden files we don't ignore
                aurostd::RemoveFile(file_path);
              }
              continue;                                                 // ignore all other hidden files
            }
            aurostd::RemoveFile(file_path);
          }
        }
        aurostd::execute("rm -f "+directory+"/.pam*");
        aurostd::execute("rm -f "+directory+"/*~");

        //DX+CO END
        // now DECOMPRESS _AFLOWIN_.EXT if a mistake was made
        for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
          ifstream FileCHECK;string FileNameCHECK;
          FileNameCHECK=directory+"/" + _AFLOWIN_ + XHOST.vext.at(iext);                    // _AFLOWIN_.EXT
          FileCHECK.open(FileNameCHECK.c_str(),std::ios::in);                         // _AFLOWIN_.EXT
          FileCHECK.clear();FileCHECK.close();                                        // _AFLOWIN_.EXT
          if(FileCHECK) {                                                             // _AFLOWIN_.EXT
            aurostd::execute(XHOST.vzip.at(iext)+" -dqf "+_AFLOWIN_+XHOST.vext.at(iext)); // _AFLOWIN_.EXT
          }
        } // _AFLOWIN_.EXT
      }
    }
  }
}

// *******************************************************************************************
// KBIN::XClean
// *******************************************************************************************
namespace KBIN {
  void XClean(string options) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy=XPID+"KBIN::XClean():";
    if(LDEBUG) cerr << soliloquy << " BEGIN" << endl;  
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=0) {
      init::ErrorOption(options,soliloquy,"aflow --xclean");
    }

    vector<string> vcheck1;aurostd::string2tokens(string("OUTCAR.static,OUTCAR.relax2,OUTCAR.relax1"),vcheck1,",");
    vector<string> vcheck2;aurostd::string2tokens(string("OUTCAR,OUTCAR,OUTCAR"),vcheck2,",");

    for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
      vcheck1.push_back("OUTCAR.relax1"+XHOST.vext.at(iext));
    }

    for(uint iext=1;iext<XHOST.vext.size();iext++) { // SKIP uncompressed
      vcheck2.push_back("OUTCAR.relax2"+XHOST.vext.at(iext));
    }

    vector<string> vfile;
    bool test=false;

    cout << soliloquy << " checking missing " << "OUTCAR*" << " with " << _AFLOWLOCK_ << endl;  // check OUTCAR.static
    aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find")+" ./ -name "+_AFLOWLOCK_),vfile);
    for(uint j=0;j<vfile.size();j++) {
      aurostd::StringSubst(vfile.at(j),_AFLOWLOCK_,"");
      if(!aurostd::FileExist(vfile.at(j)+"OUTCAR") && !aurostd::EFileExist(vfile.at(j)+"OUTCAR.relax1")) {
        cout << soliloquy << " cleaning=" << vfile.at(j) << endl;
        if(!test) KBIN::Clean(vfile.at(j));
      }
    }
    for(uint i=0;i<vcheck1.size();i++) {
      cout << soliloquy << " checking missing " << vcheck2.at(i) << " with " << vcheck1.at(i) << endl;  // check OUTCAR.static
      aurostd::string2vectorstring(aurostd::execute2string(XHOST.command("find")+" ./ -name "+vcheck1.at(i)),vfile);
      for(uint j=0;j<vfile.size();j++) {
        aurostd::StringSubst(vfile.at(j),vcheck1.at(i),"");
        if(!aurostd::FileExist(vfile.at(j)+vcheck2.at(i))) {
          cout << soliloquy << " cleaning=" << vfile.at(j) << endl;
          if(!test) KBIN::Clean(vfile.at(j));
        }
      }
    }    
    if(LDEBUG) cerr << soliloquy << " END" << endl;  
  }
} // namespace KBIN

// *******************************************************************************************
// KBIN::XClean
// *******************************************************************************************

//ME20200219 - based on CO's code in old PhononCalculator
namespace KBIN {
  int get_NCPUS() {
    _kflags kflags;
    return get_NCPUS(kflags);
  }

  int get_NCPUS(const _kflags& kflags) {
    string ncpus_str = "MAX";
    int ncpus = 1;
    if (kflags.KBIN_MPI_NCPUS > 0) ncpus_str = aurostd::utype2string<int>(kflags.KBIN_MPI_NCPUS);
    if (XHOST.vflag_control.isdefined("XPLUG_NUM_THREADS")) ncpus_str = XHOST.vflag_control.getattachedscheme("XPLUG_NUM_THREADS");
    if (ncpus_str == "MAX") ncpus = MPI_NCPUS_MAX;
    else ncpus = aurostd::string2utype<int>(ncpus_str);
    if (ncpus < 1) ncpus = 1;
    return ncpus;
  }
}  // namespace KBIN

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
