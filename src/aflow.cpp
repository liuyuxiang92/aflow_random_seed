// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *             Stefano Curtarolo - Duke University - 2003-2021             *
// *                                                                         *
// ***************************************************************************
//
//  Copyright 2003-2021 - Stefano Curtarolo - AFLOW.ORG consortium
//
//  This file is part of AFLOW software.
//
//  AFLOW is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ***************************************************************************

#include "aflow.h"
#include "aflow_pflow.h"
#include "aflowlib_entry_loader.h"
#include "aflow_pocc.h"  //CO20200624
#include "aflow_anrl.h"  //DX20201104

#include <chrono> // benchmarking HE20220222

//#define  __XOPTIMIZE
//#include "aflow_array.h"

// ***************************************************************************

// bool AFLOW_PTHREADS::FLAG;
// int  AFLOW_PTHREADS::MAX_PTHREADS;
// int  AFLOW_PTHREADS::RUNNING;
// pthread_t thread[MAX_ALLOCATABLE_PTHREADS];
// int iret[MAX_ALLOCATABLE_PTHREADS];
// bool thread_busy[MAX_ALLOCATABLE_PTHREADS];

#include "aflow_test.cpp"
#include <sys/ioctl.h>
// #include <linux/kd.h>
//   0x00004B2F   KIOCSOUND     int

//CO20180729 - OBSOLETE - use xerror
//[OBSOLETE]//CO20180419 - global exception handling - START
//[OBSOLETE]AFLOWRuntimeError::AFLOWRuntimeError(const std::string& function,const std::string& message) : std::runtime_error(message),f_name(function) {}  // I/O or computer type errors (no entries loaded)
//[OBSOLETE]AFLOWRuntimeError::AFLOWRuntimeError(const std::string& function,std::stringstream& message) : std::runtime_error(message.str()),f_name(function) {message.str("");}  // I/O or computer type errors (no entries loaded)
//[OBSOLETE]string AFLOWRuntimeError::where(){return f_name;}
//[OBSOLETE]AFLOWLogicError::AFLOWLogicError(const std::string& function,const std::string& message) : std::logic_error(message),f_name(function) {}    //errors in logic, unintended (and insurmountable) use of functionality
//[OBSOLETE]AFLOWLogicError::AFLOWLogicError(const std::string& function,std::stringstream& message) : std::logic_error(message.str()),f_name(function) {message.str("");}    //errors in logic, unintended (and insurmountable) use of functionality
//[OBSOLETE]string AFLOWLogicError::where(){return f_name;}
//[OBSOLETE]//CO20180419 - global exception handling - STOP





//
//  string sinput="",soutput="";
//  vector<aflowlib::_aflowlib_entry> entries;
//  aurostd::xoption elflags;
//  aflowlib::EntryLoader el;
//
////  if(0){
////    sinput="species(Mn,Pd)";
////    el.initialize(sinput,elflags,FileMESSAGE,oss);
////    el.retrieveOutput(soutput);
////    if(soutput.empty()){
////      message << "string output empty";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
////      return false;
////    }
////    message << "string output = \""+soutput+"\"";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
////    el.retrieveOutput(entries);
////  }
//
//
//  //sinput="MnPdPt";
//  //el.initialize(sinput,elflags,FileMESSAGE,oss);
//  sinput="MnPdPt";  //CMoTaTiW  //CMoTaTiWZr
//  el.initialize(sinput,elflags,FileMESSAGE,oss);
//  //el.retrieveOutput(soutput);
//  //if(soutput.empty()){
//  //  message << "string output empty";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
//  //  return false;
//  //}
//  //message << "string output = \""+soutput+"\"";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
//
//  el.retrieveOutput(entries);
//  cout << entries.size() << std::endl;
//
//  return true;

bool SchemaTest(ostream& oss){ofstream FileMESSAGE;return SchemaTest(FileMESSAGE,oss);}
bool SchemaTest(ofstream& FileMESSAGE,ostream& oss) {
  string function = XPID+"SchemaTest()";
  _aflags aflags; aflags.Directory = ".";
  stringstream message;
  bool all_passed = true, check_passed = true;

  message << "Performing schema test.";
  pflow::logger(_AFLOW_FILE_NAME_, function, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);

  message << "Checking for internal consistency of XHOST.vschema.";
  pflow::logger(_AFLOW_FILE_NAME_, function, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
  vector<string> vschema_keys, vschema_types;
  string schema_types = "UNIT,TYPE";
  aurostd::string2tokens(schema_types, vschema_types, ",");
  string key = "";
  check_passed = true;
  for (uint i = 0; i < XHOST.vschema.vxsghost.size(); i+= 2) {
    if(XHOST.vschema.vxsghost[i].find("::NAME:") != string::npos) {
      key=aurostd::RemoveSubString(XHOST.vschema.vxsghost[i], "SCHEMA::NAME:");
      vschema_keys.push_back(XHOST.vschema.getattachedscheme("SCHEMA::NAME:" + key));
      for (uint j = 0; j < vschema_types.size(); j++) {
        if (!XHOST.vschema.isdefined("SCHEMA::" + vschema_types[j] + ":" + key)) {
          check_passed = false;
          message << "SCHEMA::" << vschema_types[j] << ":" << key << " not found.";
          pflow::logger(_AFLOW_FILE_NAME_, function, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
        }
      }
    }
  }
  message << "Schema internal consistency check " << (check_passed?"passed":"failed") << ".";
  pflow::logger(_AFLOW_FILE_NAME_, function, message, aflags, FileMESSAGE, oss, (check_passed?_LOGGER_COMPLETE_:_LOGGER_ERROR_));
  all_passed = (all_passed && check_passed);

  message << "Checking for consistency between _aflowlib_entry json and schema.";
  pflow::logger(_AFLOW_FILE_NAME_, function, message, aflags, FileMESSAGE, oss, _LOGGER_MESSAGE_);
  aflowlib::_aflowlib_entry aentry;
  string aflowlib_json = aentry.aflowlib2string("JSON", true);
  vector<string> json_keys = aurostd::extractJsonKeysAflow(aflowlib_json);

  string keys_ignore = "data_language,error_status,natoms_orig,density_orig,volume_cell_orig,volume_atom_orig,spinD_magmom_orig";
  vector<string> vkeys_ignore;
  aurostd::string2tokens(keys_ignore, vkeys_ignore, ",");
  check_passed = true;
  for (uint i = 0; i < json_keys.size(); i++) {
    if (!aurostd::WithinList(vkeys_ignore, json_keys[i]) && !aurostd::WithinList(vschema_keys, json_keys[i])) {
      check_passed = false;
      message << json_keys[i] << " not found in schema.";
      pflow::logger(_AFLOW_FILE_NAME_, function, message, aflags, FileMESSAGE, oss, _LOGGER_ERROR_);
    }
  }
  message << "Consistency check between schema and aflowlib.json " << (check_passed?"passed":"failed") << ".";
  pflow::logger(_AFLOW_FILE_NAME_, function, message, aflags, FileMESSAGE, oss, (check_passed?_LOGGER_COMPLETE_:_LOGGER_ERROR_));
  all_passed = (all_passed && check_passed);

  return all_passed;
}

bool CeramGenTest(ostream& oss){ofstream FileMESSAGE;return CeramGenTest(FileMESSAGE,oss);}  //CO20190520
bool CeramGenTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy=XPID+"CeramGenTest():";
  //bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  message << "Performing ceramics generation test";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  //./aflow --generate_ceramics --nm=N,C --m=Co,Mo,Fe,Ru,Ni,Rh,Pt,Cu,Cr,V --N=5
  vector<string> vnonmetals,vmetals;
  aurostd::string2tokens("N,C",vnonmetals,",");aurostd::string2tokens("Co,Mo,Fe,Ru,Ni,Rh,Pt,Cu,Cr,V",vmetals,",");

  vector<string> commands=pflow::GENERATE_CERAMICS(vnonmetals,vmetals,5);
  if(commands.size()!=6){
    message << "commands.size()!=6";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N
  //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V
  //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
  //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
  //C:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
  //C:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V

  if(!aurostd::WithinList(commands,"C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N")){
    message << "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  if(!aurostd::WithinList(commands,"C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V")){
    message << "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  if(!aurostd::WithinList(commands,"C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V")){
    message << "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  if(!aurostd::WithinList(commands,"C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V")){
    message << "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  if(!aurostd::WithinList(commands,"C:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V")){
    message << "C:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  if(!aurostd::WithinList(commands,"C:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V")){
    message << "C:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }

  message << "Ceramics generation test successful";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
  return true;
}
bool EgapTest(ostream& oss){ofstream FileMESSAGE;return EgapTest(FileMESSAGE,oss);}  //CO20190520
bool EgapTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy=XPID+"EgapTest():";
  bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  message << "Performing Egap test";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  string system="",path="",query="",file="",efile="",ext="",tfile="",Egap_type="";
  vector<string> files;
  xOUTCAR xout(FileMESSAGE,oss);
  double EFERMI=AUROSTD_MAX_DOUBLE,Egap=0.0;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  //FCC/Si1_ICSD_150530
  system="ICSD_WEB/FCC/Si1_ICSD_150530";

  path=AFLOWLIB_SERVER_DEFAULT+"/AFLOWDATA/"+system;
  query=path+"/?files";
  message << "Fetching: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  aurostd::url2tokens(query,files,",");
  if(files.size()==0){
    message << "Could not fetch query: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  //OUTCAR.static
  file="OUTCAR.static";
  if(!aurostd::EWithinList(files,file,efile)){
    message << "No " << file << " found within " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  ext=aurostd::GetCompressionExtension(efile);
  tfile=aurostd::TmpFileCreate("Egap_file1")+ext;
  query=path+"/"+efile;
  message << "Fetching: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  if(!(aurostd::url2file(query,tfile,LDEBUG) && aurostd::FileExist(tfile))){
    message << "Could not fetch query: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  message << "Loaded file to: " << tfile;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  if(!xout.GetPropertiesFile(tfile,!LDEBUG)){
    message << "xOUTCAR::GetProperties() failed";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
#ifndef _AFLOW_TEMP_PRESERVE_
  aurostd::RemoveFile(tfile);
#endif
  EFERMI=xout.Efermi;
  //OUTCAR.bands
  file="OUTCAR.bands";
  if(!aurostd::EWithinList(files,file,efile)){
    query=path+"/?files"; //reload query for error message
    message << "No " << file << " found within " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  ext=aurostd::GetCompressionExtension(efile);
  tfile=aurostd::TmpFileCreate("Egap_file1")+ext;
  query=path+"/"+efile;
  message << "Fetching: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  if(!(aurostd::url2file(query,tfile,LDEBUG) && aurostd::FileExist(tfile))){
    message << "Could not fetch query: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  message << "Loaded file to: " << tfile;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  if(!xout.GetPropertiesFile(tfile,!LDEBUG)){
    message << "xOUTCAR::GetProperties() failed";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
#ifndef _AFLOW_TEMP_PRESERVE_
  aurostd::RemoveFile(tfile);
#endif
  //GetBandGap
  message << "Running bandgap code";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  if(!xout.GetBandGap(EFERMI)){
    message << "xOUTCAR::GetBandGap() failed";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  Egap=+6.1000e-01;
  if(!aurostd::isequal(Egap,xout.Egap[0])){
    message << "xOUTCAR::GetBandGap() did not find Egap==" << Egap << ", found instead xOUTCAR.Egap[0]==" << xout.Egap[0];pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  Egap_type="insulator-indirect";
  if(xout.Egap_type[0]!="insulator-indirect"){
    message << "xOUTCAR::GetBandGap() did not find type==" << Egap_type << ", found instead xOUTCAR.Egap_type[0]==" << xout.Egap_type[0];pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////

  message << "Egap test successful";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
  return true;
}

bool gcdTest(ostream& oss){ofstream FileMESSAGE;return gcdTest(FileMESSAGE,oss);}  //CO20190520
bool gcdTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy = XPID + "gcdTest():";
  bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  message << "Performing gcd test";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  int a=0,b=0,x1=0,y1=0,gcd=0;

  a=25;b=15;
  aurostd::GCD(a,b,gcd,x1,y1);
  if(!(gcd==5 && x1==-1 && y1==2)){
    if(LDEBUG){
      cerr << soliloquy << " gcd(25,15) failed" << endl;
      cerr << soliloquy << " gcd=" << gcd << endl;
      cerr << soliloquy << " x=" << x1 << endl;
      cerr << soliloquy << " y=" << y1 << endl;
    }
    return FALSE;
  }

  a=25;b=0;
  aurostd::GCD(a,b,gcd,x1,y1);
  if(!(gcd==25 && x1==1 && y1==0)){
    if(LDEBUG){
      cerr << soliloquy << " gcd(25,0) failed" << endl;
      cerr << soliloquy << " gcd=" << gcd << endl;
      cerr << soliloquy << " x=" << x1 << endl;
      cerr << soliloquy << " y=" << y1 << endl;
    }
    return FALSE;
  }

  a=0;b=15;
  aurostd::GCD(a,b,gcd,x1,y1);
  if(!(gcd==15 && x1==0 && y1==1)){
    if(LDEBUG){
      cerr << soliloquy << " gcd(0,15) failed" << endl;
      cerr << soliloquy << " gcd=" << gcd << endl;
      cerr << soliloquy << " x=" << x1 << endl;
      cerr << soliloquy << " y=" << y1 << endl;
    }
    return FALSE;
  }

  a=-5100;b=30450;
  aurostd::GCD(a,b,gcd,x1,y1);
  if(!(gcd==150 && x1==-6 && y1==-1)){
    if(LDEBUG){
      cerr << soliloquy << " gcd(-5100,30450) failed" << endl;
      cerr << soliloquy << " gcd=" << gcd << endl;
      cerr << soliloquy << " x=" << x1 << endl;
      cerr << soliloquy << " y=" << y1 << endl;
    }
    return FALSE;
  }

  message << "gcd test successful";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
  return TRUE; //CO20180419
}

bool smithTest(ostream& oss){ofstream FileMESSAGE;return smithTest(FileMESSAGE,oss);}  //CO20190520
bool smithTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy = XPID + "smithTest():";
  bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  //test ehermite
  xmatrix<int> ehermite(2,2);
  aurostd::getEHermite(5,12,ehermite);
  if(!(
        ehermite[1][1]==5 &&
        ehermite[1][2]==-2 &&
        ehermite[2][1]==-12 &&
        ehermite[2][2]==5 &&
        TRUE
      )
    ){
    if(LDEBUG){cerr << soliloquy << " getEHermite(5,12) failed" << endl;}
    return FALSE;
  }

  xmatrix<int> A1(3,3),U1,V1,S1;
  A1[1][1]=3;A1[1][2]=2;A1[1][3]=1;
  A1[2][1]=5;A1[2][2]=3;A1[2][3]=1;
  A1[3][1]=6;A1[3][2]=8;A1[3][3]=9;

  aurostd::getSmithNormalForm(A1,U1,V1,S1);

  if(LDEBUG){
    cerr << soliloquy << " A=" << endl;cerr << A1 << endl;
    cerr << soliloquy << " U=" << endl;cerr << U1 << endl;
    cerr << soliloquy << " V=" << endl;cerr << V1 << endl;
    cerr << soliloquy << " S=" << endl;cerr << S1 << endl;
  }

  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]if(!(
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      U1[1][1]==24 && U1[1][2]==-13 && U1[1][3]==-1 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      U1[2][1]==13 && U1[2][2]==-7  && U1[2][3]==-1 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      U1[3][1]==2  && U1[3][2]==-1  && U1[3][3]==0  && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      TRUE
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]    )
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  ){
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  if(LDEBUG){cerr << soliloquy << " U1(1) failed of getSmithNormalForm()" << endl;}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  return FALSE;
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]if(!(
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      V1[1][1]==0  && V1[1][2]==1  && V1[1][3]==3  && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      V1[2][1]==-1 && V1[2][2]==-1 && V1[2][3]==-1 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      V1[3][1]==1  && V1[3][2]==0  && V1[3][3]==-1 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      TRUE
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]    )
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  ){
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  if(LDEBUG){cerr << soliloquy << " V1(1) failed of getSmithNormalForm()" << endl;}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  return FALSE;
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]if(!(
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      S1[1][1]==1 && S1[1][2]==0 && S1[1][3]==0 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      S1[2][1]==0 && S1[2][2]==1 && S1[2][3]==0 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      S1[3][1]==0 && S1[3][2]==0 && S1[3][3]==1 &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      TRUE
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]    )
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  ){
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  if(LDEBUG){cerr << soliloquy << " S1(1) failed of getSmithNormalForm()" << endl;}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  return FALSE;
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]}

  xmatrix<long long int> A2(5,5),U2,V2,S2;  //long long int is CRUCIAL, Matlab actually gets this wrong because it uses long int by default
  A2[1][1]=25;    A2[1][2]=-300;   A2[1][3]=1050;    A2[1][4]=-1400;   A2[1][5]=630;
  A2[2][1]=-300;  A2[2][2]=4800;   A2[2][3]=-18900;  A2[2][4]=26880;   A2[2][5]=-12600;
  A2[3][1]=1050;  A2[3][2]=-18900; A2[3][3]=79380;   A2[3][4]=-117600; A2[3][5]=56700;
  A2[4][1]=-1400; A2[4][2]=26880;  A2[4][3]=-117600; A2[4][4]=179200;  A2[4][5]=-88200;
  A2[5][1]=630;   A2[5][2]=-12600; A2[5][3]=56700;   A2[5][4]=-88200;  A2[5][5]=44100;

  aurostd::getSmithNormalForm(A2,U2,V2,S2);

  if(LDEBUG){ //COME BACK AND PATCH FOR ANSWERS
    cerr << soliloquy << " A=" << endl;cerr << A2 << endl;
    cerr << soliloquy << " U=" << endl;cerr << U2 << endl;
    cerr << soliloquy << " V=" << endl;cerr << V2 << endl;
    cerr << soliloquy << " S=" << endl;cerr << S2 << endl;
  }

  message << "smith test successful";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
  return TRUE; //CO20180419
}

bool coordinationTest(ostream& oss){ofstream FileMESSAGE;return coordinationTest(FileMESSAGE,oss);}  //CO20190520
bool coordinationTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy=XPID+"coordinationTest():";
  bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  xstructure str("aflowlib.duke.edu:AFLOWDATA/ICSD_WEB/FCC/Cl1Na1_ICSD_240599","CONTCAR.relax.vasp",IOAFLOW_AUTO);
  deque<deque<uint> > coordinations;
  str.GetCoordinations(coordinations);
  if(coordinations.size()<2){
    if(LDEBUG){cerr << soliloquy << " coordinations not found" << endl;}
    return FALSE;
  }
  if(coordinations[0].size()<2){
    if(LDEBUG){cerr << soliloquy << " coordinations[0] not found" << endl;}
    return FALSE;
  }
  if(coordinations[1].size()<2){
    if(LDEBUG){cerr << soliloquy << " coordinations[1] not found" << endl;}
    return FALSE;
  }
  //first iatom
  //first shell
  if(coordinations[0][0]!=6){
    if(LDEBUG){cerr << soliloquy << " coordinations[0][0]!=6 (==" << coordinations[0][0] << ")" << endl;}
    return FALSE;
  }
  //second shell
  if(coordinations[0][1]!=12){
    if(LDEBUG){cerr << soliloquy << " coordinations[0][1]!=12 (==" << coordinations[0][1] << ")" << endl;}
    return FALSE;
  }
  //second iatom
  //first shell
  if(coordinations[1][0]!=6){
    if(LDEBUG){cerr << soliloquy << " coordinations[1][0]!=6 (==" << coordinations[1][0] << ")" << endl;}
    return FALSE;
  }
  //second shell
  if(coordinations[1][1]!=12){
    if(LDEBUG){cerr << soliloquy << " coordinations[1][1]!=12 (==" << coordinations[1][1] << ")" << endl;}
    return FALSE;
  }

  message << "coordination test successful";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
  return TRUE; //CO20180419
}

bool PrototypeGeneratorTest(ostream& oss, bool check_symmetry, bool check_uniqueness){ofstream FileMESSAGE;return PrototypeGeneratorTest(FileMESSAGE,oss,check_symmetry,check_uniqueness);} //DX20200925
bool PrototypeGeneratorTest(ofstream& FileMESSAGE,ostream& oss,bool check_symmetry, bool check_uniqueness){  //DX20200925
  string function_name=XPID+"PrototypeGeneratorTest():";
  bool LDEBUG=FALSE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=aurostd::getPWD();

  message << "Testing generation of all AFLOW prototypes" << (check_symmetry?" AND checking symmetry of all generated AFLOW prototypes":check_uniqueness?" AND checking all AFLOW prototypes are unique":"");
  pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  vector<string> prototype_labels, compositions;
  vector<uint> space_group_numbers;
  vector<vector<vector<string> > > grouped_Wyckoff_letters;
  string library = "anrl";

  uint num_protos = aflowlib::GetAllPrototypeLabels(prototype_labels,
      compositions,
      space_group_numbers,
      grouped_Wyckoff_letters,
      library);

  message << "Number of prototype labels = " << num_protos << " (each may have multiple parameter sets)";
  pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  string catalog="anrl";
  for(uint i=0;i<num_protos;i++){
    // get parameters
    vector<string> parameter_sets = anrl::getANRLParameters(prototype_labels[i],"all");
    if(LDEBUG){ cerr << "Number of parameters for label=" << prototype_labels[i] << ": " << parameter_sets.size() << endl; }

    for(uint j=0;j<parameter_sets.size();j++){
      xstructure xstr;
      try{
        xstr = aflowlib::PrototypeLibraries(oss,prototype_labels[i],parameter_sets[j],1);
      }
      catch(aurostd::xerror& excpt){
        message << "Could not generate prototype=" << prototype_labels[i] << " given parameters=" << parameter_sets[j] << "; check inputs or the symbolic generator.";
        pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
        return false;
      }

      // check symmetry
      if(check_symmetry){
        if(LDEBUG){ cerr << "Check that the generated structure is consistent with the label=" << prototype_labels[i] << ": " << parameter_sets.size() << endl; }

        // symmetry tolerances
        // some prototype require special tolerance values
        stringstream label_input_ss; label_input_ss << prototype_labels[i] << "-" << std::setw(3) << std::setfill('0') << j+1;
        string label_input = label_input_ss.str();
        double tolerance_sym = anrl::specialCaseSymmetryTolerances(label_input);

        string updated_label_and_params = "";
        if(!anrl::structureAndLabelConsistent(xstr, prototype_labels[i], updated_label_and_params, tolerance_sym)){ //DX20201105 - added symmetry tolerance
          // if changes symmetry, give the appropriate label
          message << "The structure has a higher symmetry than indicated by the label ";
          message << "(orig: proto=" << prototype_labels[i] << " and " << parameter_sets[j] << "). ";
          message << "The correct label and parameters for this structure are:" << endl;
          message << updated_label_and_params << endl;
          message << "Please feed this label and set of parameters into the prototype generator.";
          pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
          return false;
        }
      }
      // check uniqueness
      if(check_uniqueness){
        aurostd::xoption vpflow;
        stringstream label_input_ss; label_input_ss << prototype_labels[i] << "-" << std::setw(3) << std::setfill('0') << j+1;
        string label_input = label_input_ss.str();
        // use special symmetry tolerance if necessary (otherwise, we won't check the prototypes with the correct symmetry)
        double sym_eps = anrl::specialCaseSymmetryTolerances(label_input);
        if(sym_eps!=AUROSTD_MAX_DOUBLE){;
          xstr.sym_eps = sym_eps;
          xstr.sym_eps_calculated = true;
        }
        // check if the prototype matches to more than one prototype
        // (i.e., a prototype should match with itself, but no others)
        vector<string> protos_matching = compare::getMatchingPrototypes(xstr, vpflow, catalog);
        // if it matches to more than one
        if(protos_matching.size()>1 && !anrl::isSpecialCaseEquivalentPrototypes(protos_matching)){
          message << "ERROR: " << prototype_labels[i] << " given parameters=" << parameter_sets[j] << " matches to more than one prototype (and not a documented special case): ";
          message << aurostd::joinWDelimiter(protos_matching,",") << ". ";
          message << "If the prototype was newly added, ONLY include it in the encyclopedia for a valid reason (e.g., historical, special designation, etc.)";
          message << " and document this in anrl::isSpecialCaseEquivalentPrototypes().";
          pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
          return false;
        }
        // if it doesn't match with ITSELF
        if(protos_matching.size()==0){
          message << "ERROR: " << prototype_labels[i] << " given parameters=" << parameter_sets[j] << " does NOT match to any prototypes ";
          message << "(either this system requires a special symmetry tolerance or there is a bug with XtalFinder)." << endl;
          pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
        }
      }
    }
  }
  message << "Successfully generated all prototypes!";
  pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);

  return true;
}

bool FoldAtomsInCellTest(ostream& oss){ofstream FileMESSAGE;return FoldAtomsInCellTest(FileMESSAGE,oss);} //DX20210129
bool FoldAtomsInCellTest(ofstream& FileMESSAGE,ostream& oss){ //DX20210129
  string function_name=XPID+"FoldAtomsInCellTest():";
  //bool LDEBUG=FALSE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=aurostd::getPWD();

  // ---------------------------------------------------------------------------
  // generate rocksalt structure
  string prototype_label = "AB_cF8_225_a_b"; 
  vector<string> parameter_sets = anrl::getANRLParameters(prototype_label,"all");

  if(parameter_sets.size() != 1){
    message << "Expected only one parameter set for the rocksalt structure (" << prototype_label << ") # different parameter sets=" << parameter_sets.size() << ".";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }

  xstructure xstr;
  try{
    xstr = aflowlib::PrototypeLibraries(oss,prototype_label,parameter_sets[0],1);
  }
  catch(aurostd::xerror& excpt){
    message << "Could not generate prototype=" << prototype_label << " given parameters=" << parameter_sets[0] << "; check inputs or the symbolic generator.";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }

  // set fold atoms in cell variables
  bool skew = false;
  double tol = 0.01;
  bool check_min_dists = false;

  // ---------------------------------------------------------------------------
  // test 1: expand cell
  // create 3x1x1 supercell expansion matrix
  xmatrix<double> supercell_matrix = aurostd::eye<double>(3,3);
  supercell_matrix(1,1)=3.0;
  xmatrix<double> lattice_new = supercell_matrix*xstr.lattice; // transform lattice

  xstructure xstr_supercell = xstr;
  xstr_supercell.foldAtomsInCell(lattice_new, skew, tol, check_min_dists);

  bool same_species = true;
  if(compare::structuresMatch(xstr,xstr_supercell,same_species)){
    message << "Successfully expanded rocksalt structure into a 3x1x1 supercell via foldAtomsInCell().";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  }
  else{
    message << "Expanded rocksalt structure is not equivalent to the original structure; bad expansion.";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);
    return false;
  }

  // ---------------------------------------------------------------------------
  // test 2: reduce cell
  // convert supercell back to original lattice
  xstructure xstr_reduced = xstr_supercell;
  xstr_reduced.foldAtomsInCell(xstr.lattice, skew, tol, check_min_dists);

  if(compare::structuresMatch(xstr,xstr_reduced,same_species)){
    message << "Successfully reduced 3x1x1 rocksalt structure into a primitive form via foldAtomsInCell().";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  }
  else{
    message << "Reduced 3x1x1 rocksalt structure is not equivalent to the original structure; bad reduction.";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);
    return false;
  }

  return true;
}

// Collection of generic check functions, to streamline testing.
// HE20210616
template <typename utype>
void check(const bool &passed, const utype &calculated, const utype &expected, const string &check_function,
    const string check_description, const uint check_num, uint &passed_checks, vector<string> &results){
  stringstream result;
  if (passed) {
    passed_checks++;
    if (check_function.empty()){
      result << std::setw(3) << check_num << " | pass | " << check_description;
    } else {
      result << std::setw(3) << check_num << " | pass | " << check_function << " | " << check_description;
    }
  }
  else {
    if (check_function.empty()) {
      result << std::setw(3) << check_num << " | FAIL | " << check_description
        << " (result: " << calculated << " | expected: " << expected << ")";
    } else {
      result << std::setw(3) << check_num << " | FAIL | " << check_function << " | " << check_description
        << " (result: " << calculated << " | expected: " << expected << ")";
    }
  }
  results.push_back(result.str());
}

template <typename utype>
void check_equal(const utype &calculated, const utype &expected, const string &check_function,
    const string check_description, const uint check_num, uint &passed_checks, vector<string> &results){
  bool passed = false;
  if (aurostd::isequal(calculated, expected)) passed = true;
  check(passed, calculated, expected, check_function, check_description, check_num, passed_checks, results);
}
void check_equal(const string &calculated, const string &expected, const string &check_function, 
    const string check_description, const uint check_num, uint &passed_checks, vector<string> &results){
  bool passed = false;
  if (calculated == expected) passed = true;
  check(passed, calculated.substr(0,250), expected, check_function, check_description, check_num, passed_checks, results); //HE20220124 Prevent a screen overflow if the received answer is too long
}
void check_equal(const bool &calculated, const bool &expected, const string &check_function,
    const string check_description, const uint check_num, uint &passed_checks, vector<string> &results){
  bool passed = false;
  if (calculated == expected) passed = true;
  check(passed, calculated, expected, check_function, check_description, check_num, passed_checks, results);
}

void check_similar(const double &calculated, const double &expected, const string &check_function,
    const string check_description, const uint check_num, uint &passed_checks, vector<string> &results, const double &relative=1E-10){
  bool passed = false;
  if (std::abs(expected - calculated) <= expected * relative) passed = true;
  check(passed, calculated, expected, check_function, check_description, check_num, passed_checks, results);
}

bool display_result(const uint passed_checks, const uint check_num, const string & task_description, const vector<string> & results,
    const string & function_name, ofstream & FileMESSAGE, ostream & oss){
  stringstream message;
  _aflags aflags;
  aflags.Directory=aurostd::getPWD();
  if (passed_checks == check_num) {
    message << "SUCCESS " << task_description << " (passing " << check_num << " checks)" << endl;
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
    message << "\t" << aurostd::joinWDelimiter(results, "\n\t");
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);
  }
  else {
    message << "FAIL " << task_description << " (" << check_num-passed_checks << " of " << check_num << " checks failed)" << endl;
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    message << "\t" << aurostd::joinWDelimiter(results, "\n\t");
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);
    return false;
  }
  return true;
}


bool EntryLoaderTest(ostream& oss){ofstream FileMESSAGE;return EntryLoaderTest(FileMESSAGE,oss);}  //CO20200520
bool EntryLoaderTest(ofstream& FileMESSAGE,ostream& oss) {  //CO20200520
  string function_name = XPID + __func__;

  // setup test environment
  string task_description = "Testing EntryLoader";
  vector<string> results;
  stringstream result;
  stringstream check_description;
  uint passed_checks = 0;
  string check_function = "";
  uint check_num = 0;

  std::string test_alloy = "MnPdPt";
  bool recursive = false;
  aflowlib::EntryLoader el;
  aflowlib::_aflowlib_entry test_entry;
  xstructure test_structure;


  size_t expected = 0;
  std::vector<std::string> test_AUIDs = {
      "aflow:2de63b1ebe0a1a83",
      "4d8cf7edb50d1901",
      "auid:6d47aa3f4f1286d0",
      "aflow:7dd846bc04c764e8",
      "9d84facf8161aa60",
      "broken"
  };
  std::string test_AUID = "aflow:0d16c1946df2435c";

  std::vector<std::string> test_AURLs = {
      "aflowlib.duke.edu:AFLOWDATA/LIB2_WEB/Ca_svCu_pv/84",
      "AFLOWDATA/LIB2_WEB/Ca_svCu_pv/546",
      "aurl:AFLOWDATA/LIB2_WEB/Ca_svCu_pv/724.BA",
      "LIB2_WEB/Ca_svCu_pv/253",
      "aflowlib.duke.edu:AFLOWDATA/LIB2_WEB/Ca_svCu_pv/230"
  };
  std::string test_AURL ="aflowlib.duke.edu:AFLOWDATA/LIB2_WEB/Ca_svCu_pv/539";

  std::map<std::string, aflowlib::EntryLoader::Source> test_sources = {
      {"AUTO SELECT", aflowlib::EntryLoader::Source::NONE},
      {"SQLITE", aflowlib::EntryLoader::Source::SQLITE},
      {"AFLUX", aflowlib::EntryLoader::Source::AFLUX},
      {"FILESYSTEM", aflowlib::EntryLoader::Source::FILESYSTEM},
      {"RESTAPI", aflowlib::EntryLoader::Source::RESTAPI},
      {"FILESYSTEM_RAW", aflowlib::EntryLoader::Source::FILESYSTEM_RAW},
      {"RESTAPI_RAW", aflowlib::EntryLoader::Source::RESTAPI_RAW}
      };

  std::map<std::string, aflowlib::EntryLoader::Source> short_test_sources = {
      {"SQLITE", aflowlib::EntryLoader::Source::SQLITE},
      {"AFLUX", aflowlib::EntryLoader::Source::AFLUX},
      {"FILESYSTEM", aflowlib::EntryLoader::Source::FILESYSTEM},
      {"RESTAPI", aflowlib::EntryLoader::Source::RESTAPI},
  };

  // ---------------------------------------------------------------------------
  // Check | load alloys

  for (auto source: test_sources) {
    aurostd::StringstreamClean(check_description);
    check_num++;
    check_function = "EntryLoader::loadAlloy()";
    if (source.first == "RESTAPI" || source.first == "RESTAPI_RAW") recursive = false;
    else recursive = true;
    check_description << source.first << " - " << test_alloy;
    if (recursive) {
      check_description << " - recursive";
      expected = 2500;
    } else expected = 90;
    el.clear();
    el.m_out_silent = true;
    {
      auto start = std::chrono::high_resolution_clock::now();
      if (el.setSource(source.second)) {
        el.loadAlloy(test_alloy, recursive);
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);
        check_description << " | speed " << el.m_entries_flat->size() / (duration.count() / 1000.0) << " entries/s; "
                          << el.m_entries_flat->size() << " entries";
        check((expected < el.m_entries_flat->size()), el.m_entries_flat->size(), expected, check_function,
              check_description.str(), check_num, passed_checks, results);
      } else {
        check_description << " | failed to load " << source.first;
        check(false, 0, 0, check_function, check_description.str(), check_num, passed_checks, results);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Check | load AUID + Xstructure

  for (auto source: short_test_sources) {
    aurostd::StringstreamClean(check_description);
    check_num++;
    check_function = "EntryLoader::loadAUID()";
    check_description << source.first << " + xstructure";
    expected = 6;
    el.clear();
    el.m_out_silent = true;
    el.m_xstructure_original = true;
    el.m_xstructure_relaxed = true;
    if (source.first == "RESTAPI" || source.first == "RESTAPI_RAW") el.m_filesystem_path="/fake/"; // force xstructure test to use REST API
    auto start = std::chrono::high_resolution_clock::now();
    if (el.setSource(source.second)) {
      el.loadAUID(test_AUID);
      if (source.first == "AFLUX") test_entry = *el.m_entries_flat->back();
      el.loadAUID(test_AUIDs);

      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::high_resolution_clock::now() - start);
      check_description << " | speed " << el.m_entries_flat->size() / (duration.count() / 1000.0) << " entries/s; "
                        << el.m_entries_flat->size() << " entries";
      check_equal(el.m_entries_flat->size(), expected, check_function,
            check_description.str(), check_num, passed_checks, results);
    } else {
      check_description << " | failed to load " << source.first;
      check(false, 0, 0, check_function, check_description.str(), check_num, passed_checks, results);
    }
  }

  // ---------------------------------------------------------------------------
  // Check | load xstructure from file
  aurostd::StringstreamClean(check_description);
  check_num++;
  check_function = "EntryLoader::loadXstructureFile()";
  el.loadXstructureFile(test_entry, test_structure);
  check_description << "load xstructure extern";
  check_equal(test_structure.atoms.size(), (size_t) 6, check_function, check_description.str(), check_num, passed_checks, results);


  // ---------------------------------------------------------------------------
  // Check | load AURL

  for (auto source: short_test_sources) {
    aurostd::StringstreamClean(check_description);
    check_num++;
    check_function = "EntryLoader::loadAURL()";
    check_description << source.first;
    expected = 6;
    el.clear();
    el.m_out_silent = true;
    {
      auto start = std::chrono::high_resolution_clock::now();
      if (el.setSource(source.second)) {
        el.loadAURL(test_AURLs);
        el.loadAURL(test_AURL);
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start);
        check_description << " | speed " << el.m_entries_flat->size() / (duration.count() / 1000.0) << " entries/s; "
                          << el.m_entries_flat->size() << " entries";
        check_equal(el.m_entries_flat->size(), expected, check_function,
                    check_description.str(), check_num, passed_checks, results);
      } else {
        check_description << " | failed to load " << source.first;
        check(false, 0, 0, check_function, check_description.str(), check_num, passed_checks, results);
      }
    }
  }



  // present overall result
  return display_result(passed_checks, check_num, task_description, results, function_name, FileMESSAGE, oss);
}

// This should become a collection of tests regarding aurostd.
// At the moment, just the functions aurostd::volume and aurostd::area are tested here.
bool aurostdTest(ostream& oss){ofstream FileMESSAGE; return aurostdTest(FileMESSAGE,oss);} //HE20210511
bool aurostdTest(ofstream& FileMESSAGE, ostream& oss) { //HE20210511

  string function_name = XPID + "aurostdTest():";
//  stringstream message;


  // setup test environment
  string task_description = "Testing aurostd";
  vector<string> results;
  stringstream result;
  uint passed_checks = 0;
  string check_function = "";
  string check_description = "";
  uint check_num = 0;

  double expected = 0.0;
  double calculated = 0.0;

  int expected_int = 0;
  string expected_error = "";
  vector<xvector<double> > points;
  vector<xvector<int> > ipoints;
  vector<vector<uint> > facets;
  vector<uint> facet;

  // variables to store examples as doubles (p#) and int (p#i) variants
  xvector<double> p0(3,1); xvector<int> p0i(3,1);
  xvector<double> p1(3,1); xvector<int> p1i(3,1);
  xvector<double> p2(3,1); xvector<int> p2i(3,1);
  xvector<double> p3(3,1); xvector<int> p3i(3,1);
  xvector<double> p4(3,1); xvector<int> p4i(3,1);
  xvector<double> p5(3,1); xvector<int> p5i(3,1);
  xvector<double> p6(3,1); xvector<int> p6i(3,1);
  xvector<double> p7(3,1); xvector<int> p7i(3,1);
  xvector<double> p8(3,1); xvector<int> p8i(3,1);
  xvector<double> p9(3,1); xvector<int> p9i(3,1);
  xvector<double> p10(3,1); xvector<int> p10i(3,1);
  xvector<double> p11(3,1); xvector<int> p11i(3,1);

  // define convex solid
  p0i(1) = p0(1) = 0.0; p0i(2) = p0(2) = 0.0; p0i(3) = p0(3) = 0.0;
  p1i(1) = p1(1) = 1.0; p1i(2) = p1(2) = 0.0; p1i(3) = p1(3) = 0.0;
  p2i(1) = p2(1) = 1.0; p2i(2) = p2(2) = 1.0; p2i(3) = p2(3) = 0.0;
  p3i(1) = p3(1) = 0.0; p3i(2) = p3(2) = 1.0; p3i(3) = p3(3) = 0.0;
  p4i(1) = p4(1) = 0.0; p4i(2) = p4(2) = 0.0; p4i(3) = p4(3) = 2.0;
  p5i(1) = p5(1) = 1.0; p5i(2) = p5(2) = 0.0; p5i(3) = p5(3) = 2.0;
  p6i(1) = p6(1) = 1.0; p6i(2) = p6(2) = 1.0; p6i(3) = p6(3) = 3.0;
  p7i(1) = p7(1) = 0.0; p7i(2) = p7(2) = 1.0; p7i(3) = p7(3) = 3.0;

  // transfer data into vectors
  points.clear(); ipoints.clear(); facets.clear();
  points.push_back(p0); points.push_back(p1); points.push_back(p2); points.push_back(p3); points.push_back(p4);
  points.push_back(p5); points.push_back(p6); points.push_back(p7);
  ipoints.push_back(p0i); ipoints.push_back(p1i); ipoints.push_back(p2i); ipoints.push_back(p3i); ipoints.push_back(p4i);
  ipoints.push_back(p5i); ipoints.push_back(p6i); ipoints.push_back(p7i);
  facet.clear(); facet.resize(4); facet[0]=0; facet[1]=1; facet[2]=2; facet[3]=3; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=4; facet[1]=5; facet[2]=6; facet[3]=7; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=1; facet[1]=2; facet[2]=6; facet[3]=5; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=0; facet[1]=3; facet[2]=7; facet[3]=7; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=0; facet[1]=1; facet[2]=5; facet[3]=4; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=3; facet[1]=2; facet[2]=6; facet[3]=7; facets.push_back(facet);

  // ---------------------------------------------------------------------------
  // Check | convex solid volume (double)
  check_num++;
  check_function = "aurostd::volume()";
  check_description = "convex solid, points as doubles";
  expected = 2.5;

  calculated = aurostd::volume(points, facets, true);
  check_equal(calculated, expected, check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | convex solid volume (int)
  check_num++;
  check_function = "aurostd::volume()";
  check_description = "convex solid, points as int";

  calculated = aurostd::volume(ipoints, facets, true);
  check_equal(calculated, expected, check_function, check_description, check_num, passed_checks, results);


  // define non convex solid
  p0i(1) = p0(1)   = 0.0; p0i(2) = p0(2)   = 0.0; p0i(3) = p0(3)   = 0.0;
  p1i(1) = p1(1)   = 0.0; p1i(2) = p1(2)   = 4.0; p1i(3) = p1(3)   = 0.0;
  p2i(1) = p2(1)   = 2.0; p2i(2) = p2(2)   = 4.0; p2i(3) = p2(3)   = 0.0;
  p3i(1) = p3(1)   = 1.0; p3i(2) = p3(2)   = 1.0; p3i(3) = p3(3)   = 0.0;
  p4i(1) = p4(1)   = 4.0; p4i(2) = p4(2)   = 2.0; p4i(3) = p4(3)   = 0.0;
  p5i(1) = p5(1)   = 4.0; p5i(2) = p5(2)   = 0.0; p5i(3) = p5(3)   = 0.0;
  p6i(1) = p6(1)   = 0.0; p6i(2) = p6(2)   = 0.0; p6i(3) = p6(3)   = 4.0;
  p7i(1) = p7(1)   = 0.0; p7i(2) = p7(2)   = 4.0; p7i(3) = p7(3)   = 4.0;
  p8i(1) = p8(1)   = 2.0; p8i(2) = p8(2)   = 4.0; p8i(3) = p8(3)   = 4.0;
  p9i(1) = p9(1)   = 1.0; p9i(2) = p9(2)   = 1.0; p9i(3) = p9(3)   = 4.0;
  p10i(1) = p10(1) = 4.0; p10i(2) = p10(2) = 2.0; p10i(3) = p10(3) = 4.0;
  p11i(1) = p11(1) = 4.0; p11i(2) = p11(2) = 0.0; p11i(3) = p11(3) = 4.0;

  // transfer data into vectors
  points.clear(); ipoints.clear(); facets.clear();
  points.push_back(p0); points.push_back(p1); points.push_back(p2); points.push_back(p3); points.push_back(p4);
  points.push_back(p5); points.push_back(p6); points.push_back(p7); points.push_back(p8); points.push_back(p9);
  points.push_back(p10); points.push_back(p11);
  ipoints.push_back(p0i); ipoints.push_back(p1i); ipoints.push_back(p2i); ipoints.push_back(p3i); ipoints.push_back(p4i);
  ipoints.push_back(p5i); ipoints.push_back(p6i); ipoints.push_back(p7i); ipoints.push_back(p8i); ipoints.push_back(p9i);
  ipoints.push_back(p10i); ipoints.push_back(p11i);

  facet.clear(); facet.resize(6); facet[0]=5; facet[1]=4; facet[2]=3; facet[3]=2; facet[4]=1; facet[5]=0; facets.push_back(facet);
  facet.clear(); facet.resize(6); facet[0]=6; facet[1]=7; facet[2]=8; facet[3]=9; facet[4]=10; facet[5]=11; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=0; facet[1]=6; facet[2]=11; facet[3]=5; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=4; facet[1]=5; facet[2]=11; facet[3]=10; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=3; facet[1]=4; facet[2]=10; facet[3]=9; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=3; facet[1]=9; facet[2]=8; facet[3]=2; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=1; facet[1]=2; facet[2]=8; facet[3]=7; facets.push_back(facet);
  facet.clear(); facet.resize(4); facet[0]=0; facet[1]=1; facet[2]=7; facet[3]=6; facets.push_back(facet);

  // ---------------------------------------------------------------------------
  // Check | non convex solid volume (double)
  check_num++;
  check_function = "aurostd::volume()";
  check_description = "non convex solid, points as doubles";
  expected = 40.0;

  calculated = aurostd::volume(points, facets);
  check_equal(calculated, expected, check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | error facet/normals mismatch
  check_num++;
  check_function = "aurostd::volume()";
  check_description = "error: facet/normals mismatch";
  vector<xvector<double> > normals;
  expected_error = "xerror code 30 (VALUE_ERROR)";
  expected_int = _VALUE_ERROR_;

  try {
    calculated = aurostd::volume(points, facets, normals);
    check(false, std::string("no error"), expected_error, check_function, check_description, check_num, passed_checks, results);
  }
  catch (aurostd::xerror e)
  {
    if (e.error_code == expected_int) check(true, "", "", check_function, check_description, check_num, passed_checks, results);
    else check(false, aurostd::utype2string(e.error_code), expected_error, check_function, check_description, check_num, passed_checks, results);
  }
  catch (...) {
    check(false, std::string("not an xerror"), expected_error, check_function, check_description, check_num, passed_checks, results);
  }

  // ---------------------------------------------------------------------------
  // Check | non convex solid volume (int)
  check_num++;
  check_function = "aurostd::volume()";
  check_description = "non convex solid, points as int";
  expected = 40.0;

  calculated = aurostd::volume(ipoints, facets);
  check_equal(calculated, expected, check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | error facet size
  check_num++;
  check_function = "aurostd::volume()";
  check_description = "error: wrong facet size";
  expected_error = "xerror code 30 (VALUE_ERROR)";
  expected_int = _VALUE_ERROR_;

  facet.clear(); facet.resize(2); facet[0]=1; facet[1]=2; facets.push_back(facet);
  try {
    calculated = aurostd::volume(points, facets);
    check(false, std::string("no error"), expected_error, check_function, check_description, check_num, passed_checks, results);
  }
  catch (aurostd::xerror e)
  {
    if (e.error_code == expected_int) check(true, "", "", check_function, check_description, check_num, passed_checks, results);
    else check(false, aurostd::utype2string(e.error_code), expected_error, check_function, check_description, check_num, passed_checks, results);
  }
  catch (...) {
    check(false, std::string("not an xerror"), expected_error, check_function, check_description, check_num, passed_checks, results);
  }

  // ---------------------------------------------------------------------------
  // Check | non convex area (double)
  check_num++;
  check_function = "aurostd::areaPointsOnPlane()";
  check_description = "non convex area; points as double";
  expected = 10.0;

  //fill vectors with data
  points.clear(); ipoints.clear(); facets.clear();
  points.push_back(p0); points.push_back(p1); points.push_back(p2); points.push_back(p3); points.push_back(p4);
  points.push_back(p5);
  ipoints.push_back(p0i); ipoints.push_back(p1i); ipoints.push_back(p2i); ipoints.push_back(p3i); ipoints.push_back(p4i);
  ipoints.push_back(p5i);

  calculated = aurostd::areaPointsOnPlane(points);
  check_equal(calculated, expected, check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | non convex area (int)
  check_num++;
  check_function = "aurostd::areaPointsOnPlane()";
  check_description = "non convex area; points as int";
  expected = 10.0;

  calculated = aurostd::areaPointsOnPlane(ipoints);
  check_equal(calculated, expected, check_function, check_description, check_num, passed_checks, results);


  // define triangle in 3D to better test int handling
  p0i(1) = p0(1) = 0.0; p0i(2) = p0(2) = 0.0; p0i(3) = p0(3) = 0.0;
  p1i(1) = p1(1) = 1.0; p1i(2) = p1(2) = 1.0; p1i(3) = p1(3) = 1.0;
  p2i(1) = p2(1) = 5.0; p2i(2) = p2(2) = 0.0; p2i(3) = p2(3) = 5.0;

  points.clear(); ipoints.clear(); facets.clear();
  points.push_back(p0); points.push_back(p1); points.push_back(p2);
  ipoints.push_back(p0i); ipoints.push_back(p1i); ipoints.push_back(p2i);

  // ---------------------------------------------------------------------------
  // Check | 3d triangle area (double)
  check_num++;
  check_function = "aurostd::areaPointsOnPlane()";
  check_description = "3d triangle; points as double";
  expected = 3.5355339059;

  calculated = aurostd::areaPointsOnPlane(points);
  check_similar(calculated, expected, check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | 3d triangle area (int)
  check_num++;
  check_function = "aurostd::areaPointsOnPlane()";
  check_description = "3d triangle; points as int";
  expected = 3.5355339059;

  calculated = aurostd::areaPointsOnPlane(ipoints);
  check_similar(calculated, expected, check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | double2fraction conversion //DX20210908
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  check_function = "aurostd::double2fraction()";
  check_description = "convert a double to a fraction.";

  double test_double = 1.625;
  int numerator=1, denominator=1;
  string answer = "13/8";
  aurostd::double2fraction(test_double,numerator,denominator);
  stringstream result_ss; result_ss << numerator << "/" << denominator;

  check_num++;
  check_equal(result_ss.str(), answer, check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check http get function
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  check_function = "aurostd::httpGet()";

  vector<std::string> urls;
  urls.push_back("http://aflowlib.duke.edu/test/?echo=Compicated1773");
  urls.push_back("http://aflowlib.duke.edu/test/?re=301");
  urls.push_back("http://aflowlib.duke.edu/test/?re=302");
  urls.push_back("http://aflowlib.duke.edu/test/?re=307");
  urls.push_back("http://aflowlib.duke.edu/test/?re=308");

  vector<std::string> http_results;
  http_results.push_back("Compicated1773");
  http_results.push_back("Redirected");
  http_results.push_back("Redirected");
  http_results.push_back("Redirected");
  http_results.push_back("Redirected");

  vector<std::string> http_description;
  http_description.push_back("Echo test");
  http_description.push_back("301 forward test");
  http_description.push_back("302 forward test");
  http_description.push_back("307 forward test");
  http_description.push_back("308 forward test");

  std::string response;
  for (uint i_task=0; i_task<urls.size(); i_task++){
    response.clear();
    response = aurostd::httpGet(urls[i_task]);
    check_num++;
    check_equal(response, http_results[i_task], check_function, http_description[i_task], check_num, passed_checks, results);
  }

  // Test different call pattern
  response.clear();
  std::string http_single_description;
  int status_code=-1;
  std::map<std::string, std::string> http_header;

  std::string http_host = "aflowlib.duke.edu";
  std::string http_path = "/test/";
  std::string http_query = "?echo=";
  std::string url = "http://" + http_host + http_path + http_query;

  http_single_description = "output=httpGet(url,status_code)";
  response = aurostd::httpGet(url+http_single_description, status_code);
  check_num++;
  check_equal(response, http_single_description, check_function, http_single_description, check_num, passed_checks, results);

  http_single_description = "output=httpGet(url,status_code,header)";
  response = aurostd::httpGet(url+http_single_description, status_code, http_header);
  check_num++;
  check_equal(response, http_single_description, check_function, http_single_description, check_num, passed_checks, results);

  http_single_description = "output=httpGet(host,path,query)";
  response = aurostd::httpGet(http_host, http_path, http_query+http_query+http_single_description);
  check_num++;
  check_equal(response, http_single_description, check_function, http_single_description, check_num, passed_checks, results);

  http_single_description = "output=httpGet(host,path,query,status_code)";
  response = aurostd::httpGet(http_host, http_path, http_query+http_single_description, status_code);
  check_num++;
  check_equal(response, http_single_description, check_function, http_single_description, check_num, passed_checks, results);

  http_single_description = "output=httpGet(host,path,query,status_code,header)";
  response = aurostd::httpGet(http_host, http_path, http_query+http_single_description, status_code, http_header);
  check_num++;
  check_equal(response, http_single_description, check_function, http_single_description, check_num, passed_checks, results);

  check_function = "aurostd::httpGetStatus()";
  http_single_description = "status_code=httpGetStatus(url,output)";
  status_code = aurostd::httpGetStatus(url+http_single_description, response);
  check_num++;
  check_equal(response, http_single_description, check_function, http_single_description, check_num, passed_checks, results);

  http_single_description = "status_code=httpGetStatus(url,output,header)";
  status_code = aurostd::httpGetStatus(url+http_single_description, response, http_header);
  check_num++;
  check_equal(response, http_single_description, check_function, http_single_description, check_num, passed_checks, results);

  http_single_description = "status_code=httpGetStatus(host,path,query,output)";
  status_code = aurostd::httpGetStatus(http_host, http_path, http_query+http_single_description, response);
  check_num++;
  check_equal(response, http_single_description, check_function, http_single_description, check_num, passed_checks, results);

  http_single_description = "status_code=httpGetStatus(host,path,query,output,header)";
  status_code = aurostd::httpGetStatus(http_host, http_path, http_query+http_single_description, response, http_header);
  check_num++;
  check_equal(response, http_single_description, check_function, http_single_description, check_num, passed_checks, results);

  check_function = "aurostd::httpPercentEncodingSelected()";
  response = aurostd::httpPercentEncodingSelected("Crazy!? _&@String", "@&");
  check_num++;
  check_equal(response, "Crazy!? _%26%40String", check_function, "just escape @ and &", check_num, passed_checks, results);

  check_function = "aurostd::httpPercentEncodingFull()";
  response = aurostd::httpPercentEncodingFull("[(Crazy!? _&@String}/|\\äßæ");
  check_num++;
  check_equal(response, "%5B%28Crazy%21%3F%20_%26%40String%7D%2F%7C%5C%C3%A4%C3%9F%C3%A6", check_function, "escape all chars", check_num, passed_checks, results);


  // present overall result
  return display_result(passed_checks, check_num, task_description, results, function_name, FileMESSAGE, oss);
}

bool AtomicEnvironmentTest(ostream& oss){ofstream FileMESSAGE;return AtomicEnvironmentTest(FileMESSAGE,oss);} //HE20210511
bool AtomicEnvironmentTest(ofstream& FileMESSAGE, ostream& oss){ //HE20210511

  string function_name = XPID + "AtomicEnvironmentTest():";

  // setup test environment
  string task_description = "Creating a atomic environment [aflow:d912e209c81aeb94]";
  vector<string> results;
  stringstream result;
  uint passed_checks = 0;
  string check_function = "";
  string check_description = "";
  uint check_num = 0;

  // ---------------------------------------------------------------------------
  // Test 1: create AE - mode 1
  // ---------------------------------------------------------------------------

  // load test system
  xstructure str;
  {
    aflowlib::EntryLoader el;
    el.m_xstructure_relaxed = true;
    el.m_out_silent = true;
    el.loadAUID("aflow:d912e209c81aeb94"); // aflowlib.duke.edu:AFLOWDATA/LIB2_RAW/Ca_svCu_pv/138
    str = el.m_entries_flat->operator[](0)->vstr.back();
  }

  vector<AtomEnvironment> AE = getAtomEnvironments(str, 1);

  // ---------------------------------------------------------------------------
  // Check | number of created AEs
  check_num++;
  check_description = "number of created AEs";
  check_equal(uint(AE.size()), uint(6), check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | point index mapping
  check_num++;
  check_description = "point index mapping";
  check_equal(AE[1].index2Point(10), AE[1].coordinates_neighbor[1][2], check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | coordinate matching
  check_num++;
  check_description = "coordinate matching";
  xvector<double> compare_point(3,1);
  compare_point(1) = -2.95227319118354;
  compare_point(2) = 1.02047725834487;
  compare_point(3) = 2.42721573754061;

  check_equal(AE[1].index2Point(2), compare_point, check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | center element
  check_num++;
  check_description = "center element";
  check_equal(std::string(AE[2].element_center), std::string("Ca"), check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | center element id
  check_num++;
  check_description = "center element id";
  check_equal(AE[4].type_center, uint(1), check_function, check_description, check_num, passed_checks, results);

  // present overall result
  bool isPassed = display_result(passed_checks, check_num, task_description, results, function_name, FileMESSAGE, oss);
  if (!isPassed) return isPassed;
  // ---------------------------------------------------------------------------
  // Test 2: create AE convex hull
  // ---------------------------------------------------------------------------

  // setup test environment
  results.clear();
  result.str("");
  result.clear();
  task_description = "Creating convex hull with constructAtomEnvironmentHull() [aflow:d912e209c81aeb94]";
  passed_checks = 0;
  check_num = 0;
  const uint test_AE = 4;

  // create hull
  AE[test_AE].constructAtomEnvironmentHull();

  // ---------------------------------------------------------------------------
  // Check | hull bit set
  check_num++;
  check_description = "hull bit set";
  check_equal(AE[test_AE].has_hull, true, check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | hull volume
  check_num++;
  check_description = "hull volume";
  check_similar(AE[test_AE].volume, 33.0352927028572, check_function, check_description, check_num, passed_checks, results);
  cout << std::setprecision(15) << AE[test_AE].volume << endl;
  // ---------------------------------------------------------------------------
  // Check | hull area
  check_num++;
  check_description = "hull area";
  check_similar(AE[test_AE].area, 62.4980204713889, check_function, check_description, check_num, passed_checks, results);
  cout << std::setprecision(15) << AE[test_AE].area << endl;
  // ---------------------------------------------------------------------------
  // Check | triangle count
  check_num++;
  check_description = "triangle count";
  check_equal(AE[test_AE].facet_order[0], uint(6), check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | tetragon count
  check_num++;
  check_description = "tetragon count";
  check_equal(AE[test_AE].facet_order[1], uint(2), check_function, check_description, check_num, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | pentagon count
  check_num++;
  check_description = "pentagon count";
  check_equal(AE[test_AE].facet_order[2], uint(0), check_function, check_description, check_num, passed_checks, results);

  // present overall result
  return display_result(passed_checks, check_num, task_description, results, function_name, FileMESSAGE, oss);
}

int main(int _argc,char **_argv) {
  string soliloquy = XPID + "main():"; //CO20180419
  ostream& oss=cout;  //CO20180419
  try{
    bool LDEBUG=FALSE; // TRUE;
    int return_code = 0;  //ME20200901
    if(LDEBUG) cerr << "AFLOW-MAIN [1]" << endl;
    std::vector<string> argv(aurostd::get_arguments_from_input(_argc,_argv));
    if(LDEBUG) cerr << "AFLOW-MAIN [2]" << endl;
    std::vector<string> cmds;

    // MACHINE
    //ME20200724
    int code = init::InitMachine(FALSE,argv,cmds,cerr);
    if (code >= 0) return code;
    if(LDEBUG || XHOST.DEBUG) cerr << "AFLOW-MAIN [3]" << endl;

    // aurostd::TmpDirectoryCreate("test");
    // cerr << args2flag(argv,"--aaa|--bbb |--ccc") << endl;
    // CHECK USERS MACHINES - DEBUG

    // initialize_templates_never_call_this_procedure(1);

    // INITIALIZE ***************************************************
    // INIT LOOK UP TABLES
    atoms_initialize();
    xelement::Initialize();
    xPOTCAR_Initialize();
    // spacegroup::SpaceGroupInitialize(); only if necessary
    // INFORMATION **************************************************
    AFLOW_PTHREADS::FLAG=AFLOW_PTHREADS::Check_Threads(argv,!XHOST.QUIET);

    bool Arun=FALSE;
    if(!Arun && aurostd::args2flag(argv,cmds,"--pocc_old2new|--pocc_o2n"))  {Arun=TRUE;pocc::poccOld2New();} //CO20200624
    if(!Arun && aurostd::args2flag(argv,cmds,"--prx|--prx="))  {Arun=TRUE;PERFORM_PRX(cout);}
    if(!Arun && aurostd::args2flag(argv,cmds,"--generate_makefile|--makefile"))  {Arun=TRUE;makefile::createMakefileAFLOW(".");}  //CO20200508 - if calling from command-line, you should be sitting inside aflow directory (will write out Makefile.aflow)
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_getpp")) {
      if(KBIN::VASP_PseudoPotential_CleanName_TEST()){return 0;}
      return 1;
    }
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=Init|--test=init")) {
      init::InitLoadString("vLIBS",1);
      return 1;
    }

    if (0)     { // works with pointers
      std::filebuf* fb_pre = new std::filebuf; fb_pre->open("friscosity_pre.txt",std::ios_base::out|std::ios_base::trunc);  // need to be generated before any assignments
      std::filebuf* fb_post = new std::filebuf; fb_post->open("friscosity_post.txt",std::ios_base::out|std::ios_base::trunc); // need to be generated before any assignments
      std::ostream* oss;
#define _oss (*oss)
      // PRE
      //   oss = fb_pre;
      //      _oss << "FRISCOSITY_PRE" << endl;
      //  _oss.flush(); 
      // COUT
      oss = &std::cout;
      _oss << "assigned COUT" << endl;
      _oss.flush();
      // CERR
      oss = &std::cerr;
      _oss << "assigned CERR" << endl;
      _oss.flush();
      return 1; //CO20200624 - debug mode
    }
    if(0) { // https://stdcxx.apache.org/doc/stdlibug/34-2.html
      std::ofstream oss;
      // default COUT
      // COUT DEFAULT
      oss.copyfmt(std::cout);                             
      oss.clear(std::cout.rdstate());                     
      oss.std::basic_ios<char>::rdbuf(std::cout.rdbuf());  //CO20210707 - patching for older TX machines
      oss << "COUT DEFAULT" << std::endl;
      if(1) { // somebody chose a file
        // FILE PRE
        std::ofstream ofs_pre("friscosity_pre.txt", std::ofstream::out);
        oss.copyfmt(ofs_pre);                             
        oss.clear(ofs_pre.rdstate());                     
        oss.std::basic_ios<char>::rdbuf(ofs_pre.rdbuf());  //CO20210707 - patching for older TX machines
      }
      oss << "FRISCOSITY_PRE" << std::endl;
      if(1) { // put it back on COUT
        // COUT
        oss.copyfmt(std::cout);                             
        oss.clear(std::cout.rdstate());                     
        oss.std::basic_ios<char>::rdbuf(std::cout.rdbuf());  //CO20210707 - patching for older TX machines
      }
      oss << "COUT" << std::endl;
      if(1) { // try CERR
        // CERR
        oss.copyfmt(std::cerr);                             
        oss.clear(std::cerr.rdstate());                     
        oss.std::basic_ios<char>::rdbuf(std::cerr.rdbuf());  //CO20210707 - patching for older TX machines
      }
      oss << "CERR" << std::endl;
      if(1) { // switch to file 
        // FILE POST
        std::ofstream ofs_post("friscosity_post.txt", std::ofstream::out);
        oss.copyfmt(ofs_post);                             
        oss.clear(ofs_post.rdstate());                     
        oss.std::basic_ios<char>::rdbuf(ofs_post.rdbuf());  //CO20210707 - patching for older TX machines
      }
      oss << "FRISCOSITY_POST" << std::endl;

      return 1; //CO20200624 - debug mode
    }


    if(!Arun && aurostd::args2flag(argv,cmds,"--test_xmatrix")) { //CO20190911
      string soliloquy = XPID + "test_xmatrix()::";
      bool LDEBUG=TRUE;// TRUE;
      xmatrix<double> mat;
      mat(1,1)=5;mat(1,2)=9;mat(1,3)=12;
      mat(2,1)=7;mat(2,2)=10;mat(2,3)=13;
      mat(3,1)=8;mat(3,2)=11;mat(3,3)=14;
      if(LDEBUG){cerr << soliloquy << " mat=" << endl;cerr << mat << endl;}
      //getmat()
      xmatrix<double> submat;
      mat.getmatInPlace(submat,2,3,2,3);
      if(LDEBUG){cerr << soliloquy << " submat=" << endl;cerr << submat << endl;}
      //setmat()
      mat.setmat(submat,1,1); //do nothing
      if(LDEBUG){
        cerr << soliloquy << " replacing with submat at 1,1" << endl;
        cerr << soliloquy << " mat=" << endl;cerr << mat << endl;
      }
      xvector<double> xv;
      xv(1)=2;xv(2)=3;xv(3)=4;
      if(LDEBUG){cerr << soliloquy << " xv=" << xv << endl;}
      mat.setmat(xv,1,false); //row
      if(LDEBUG){
        cerr << soliloquy << " replacing with xv at row=1" << endl;
        cerr << soliloquy << " mat=" << endl;cerr << mat << endl;
      }
      mat.setmat(xv,2,true); //col
      if(LDEBUG){
        cerr << soliloquy << " replacing with xv at col=2" << endl;
        cerr << soliloquy << " mat=" << endl;cerr << mat << endl;
      }
      //setrow()
      mat.setrow(xv,2);
      if(LDEBUG){
        cerr << soliloquy << " replacing with xv at row=2" << endl;
        cerr << soliloquy << " mat=" << endl;cerr << mat << endl;
      }
      //setcol()
      mat.setcol(xv,3);
      if(LDEBUG){
        cerr << soliloquy << " replacing with xv at col=3" << endl;
        cerr << soliloquy << " mat=" << endl;cerr << mat << endl;
      }
      return 1;
    }
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_stefano")) {
      uint y=2017,m=11;
      m+=1;
      for(uint i=0;i<200;i++) {
        if(m==0) {
          cout << "mv \"unknown.pdf\" stefano_" << y << m << ".pdf" << endl;
        } else {
          if(m<10) {
            cout << "mv \"unknown(" << i << ").pdf\" stefano_" << y << "0" << m << ".pdf" << endl;
          } else {
            cout << "mv \"unknown(" << i << ").pdf\" stefano_" << y << m << ".pdf" << endl;
          }
        }
        m--;
        if(m==0) {y--;m+=12;} 
      }
      return 0; //CO20180419
    }
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_EntryLoader|--EntryLoader_test")) {return (EntryLoaderTest()?0:1);}  //CO20190601
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_schema|--schema_test")) {return (SchemaTest()?0:1);}  //ME20210408
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_CeramGen|--CeramGen_test")) {return (CeramGenTest()?0:1);}  //CO20190601
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_Egap|--Egap_test")) {return (EgapTest()?0:1);}  //CO20190601
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_gcd|--gcd_test")) {return (gcdTest()?0:1);}  //CO20190601
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_smith|--smith_test")) {return (smithTest()?0:1);}  //CO20190601
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_coordination|--coordination_test")) {return (coordinationTest()?0:1);}  //CO20190601
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_PrototypeGenerator|--PrototypeGenerator_test")) {return (PrototypeGeneratorTest()?0:1);}  //DX20200928
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_PrototypeSymmetry|--PrototypeSymmetry_test")) {return (PrototypeGeneratorTest(cout,true)?0:1);}  //DX20201105
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_PrototypeUniqueness|--PrototypeUniqueness_test")) {return (PrototypeGeneratorTest(cout,false,true)?0:1);}  //DX20210429
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_FoldAtomsInCell|--FoldAtomsInCell_test")) {return (FoldAtomsInCellTest(cout)?0:1);}  //DX20210129
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_AtomicEnvironment|--AtomicEnvironment_test")) {return (AtomicEnvironmentTest(cout)?0:1);}  //HE20210511
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_aurostd|--aurostd_test")) {return (aurostdTest(cout)?0:1);} //HE20210512
    if(!Arun && aurostd::args2flag(argv,cmds,"--test")) {

      if(XHOST.vext.size()!=XHOST.vcat.size()) {throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"XHOST.vext.size()!=XHOST.vcat.size(), aborting.",_RUNTIME_ERROR_);}

      for(uint iext=0;iext<XHOST.vext.size();iext++) { 
        cout << "\"" << XHOST.vext.at(iext) << "\"" << " " << "\"" << XHOST.vcat.at(iext) << "\"" << endl;
      }

      int dim=7;
      cout << "dimm=" << dim << endl;

      xmatrix<double> m(dim,dim),mi(dim,dim);
      for(int i=1;i<=dim;i++) 
        for(int j=1;j<=dim;j++)
          m(i,j)=aurostd::ran0();
      cout << "m=" << endl << m << endl;
      mi=inverse(m);
      cout << "mi=" << endl << mi << endl;
      cout << "mi*m=" << endl << det(mi*m) << endl;

      //CO how to create 64bit string from binary file
      //string b64String;
      //aurostd::bin2base64("aflow_logo.pdf",b64String);
      //cout << b64String << endl;

      string test="2.730747137  -2.730747137-12.397646334";
      vector<string> _tokens;
      aurostd::string2tokens(test,_tokens,"-");
      for(uint i=0;i<_tokens.size();i++){
        cerr << _tokens[i] << endl;
      }
      return 0; //CO20180419
      //CO START 20170614 - some SQLITE tests
      //http://zetcode.com/db/sqlitec/ - more tests here
      //this will create test.db file
      sqlite3 *db;
      char *err_msg = 0;
      int rc = sqlite3_open("test.db", &db);
      if(rc != SQLITE_OK) {
        fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return 1;
      }
      char const *sql = "DROP TABLE IF EXISTS Cars;" 
        "CREATE TABLE Cars(Id INT, Name TEXT, Price INT);" 
        "INSERT INTO Cars VALUES(1, 'Audi', 52642);" 
        "INSERT INTO Cars VALUES(2, 'Mercedes', 57127);" 
        "INSERT INTO Cars VALUES(3, 'Skoda', 9000);" 
        "INSERT INTO Cars VALUES(4, 'Volvo', 29000);" 
        "INSERT INTO Cars VALUES(5, 'Bentley', 350000);" 
        "INSERT INTO Cars VALUES(6, 'Citroen', 21000);" 
        "INSERT INTO Cars VALUES(7, 'Hummer', 41400);" 
        "INSERT INTO Cars VALUES(8, 'Volkswagen', 21600);";
      rc = sqlite3_exec(db, sql, 0, 0, &err_msg);
      if(rc != SQLITE_OK ) {
        fprintf(stderr, "SQL error: %s\n", err_msg);
        sqlite3_free(err_msg);       
        sqlite3_close(db);
        return 1;
      } 
      sqlite3_close(db);
      //return 0;

      //MORE TESTS
      //printf("%s\n,sqlite3_libversion()");
      //    sqlite3 *db;
      //sqlite3_stmt *res;
      //int rc = sqlite3_open(":memory:", &db);
      //if(rc != SQLITE_OK) {
      //    fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
      //    sqlite3_close(db);
      //    return 1;
      //}
      //rc = sqlite3_prepare_v2(db, "SELECT SQLITE_VERSION()", -1, &res, 0);   
      //if(rc != SQLITE_OK) {
      //    fprintf(stderr, "Failed to fetch data: %s\n", sqlite3_errmsg(db));
      //    sqlite3_close(db);
      //    return 1;
      //}    
      //rc = sqlite3_step(res);
      //if(rc == SQLITE_ROW) {
      //    printf("%s\n", sqlite3_column_text(res, 0));
      //}
      //sqlite3_finalize(res);
      //sqlite3_close(db);
      //return 0;

      //quick easy test
      cerr << sqlite3_libversion() << endl;
      //CO END 20170614 - some SQLITE tests
      aurostd::xcomplex<double> x(123.0,456.0);
      cout << x.re << "," << x.im << " - " << x.real() << "," << x.imag() << " - " << x << endl;
      x.re=111;x.im=222;
      cout << x.re << "," << x.im << " - " << x.real() << "," << x.imag() << " - " << x << endl;
      cout << aurostd::PaddedPOST("EMIN= -30.0",10) << endl;;
      stringstream for_corey;
      for_corey << "scatter/use mapped color={draw=black,fill=mapped color,solid}";
      string corey=for_corey.str();
      cout << corey << endl;
      stringstream aus;
      aus << "************************   00000  MESSAGE KPOINTS KSHIFT=[" << 1 << " " << 2 << " " << 3 << "]" << " ************************ " << endl;
      cout << aus.str() << endl;
      return 0; //CO20180419
    }

    if(!Arun && aurostd::args2attachedflag(argv,"--bin2base64=")) {
      string b64String;
      aurostd::bin2base64(aurostd::args2attachedstring(argv,"--bin2base64=",""),b64String);
      cout << b64String << endl;
      return 0; //CO20180419
    }

    if(!Arun && aurostd::args2flag(argv,cmds,"--test=POTCAR|--test=POTCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=POTCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=POTCAR.static"+DEFAULT_KZIP_EXT+"|--test=POTCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xPOTCAR(aurostd::args2attachedstring(argv,"--test=",""));return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=DOSCAR|--test=DOSCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=DOSCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=DOSCAR.static"+DEFAULT_KZIP_EXT+"|--test=DOSCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xDOSCAR(aurostd::args2attachedstring(argv,"--test=",""));return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=EIGENVAL|--test=EIGENVAL.relax1"+DEFAULT_KZIP_EXT+"|--test=EIGENVAL.relax2"+DEFAULT_KZIP_EXT+"|--test=EIGENVAL.static"+DEFAULT_KZIP_EXT+"|--test=EIGENVAL.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xEIGENVAL(aurostd::args2attachedstring(argv,"--test=",""));return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=OUTCAR|--test=OUTCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=OUTCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=OUTCAR.static"+DEFAULT_KZIP_EXT+"|--test=OUTCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xOUTCAR(aurostd::args2attachedstring(argv,"--test=",""));return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=KPOINTS|--test=KPOINTS.relax1"+DEFAULT_KZIP_EXT+"|--test=KPOINTS.relax2"+DEFAULT_KZIP_EXT+"|--test=KPOINTS.static"+DEFAULT_KZIP_EXT+"|--test=KPOINTS.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xKPOINTS(aurostd::args2attachedstring(argv,"--test=",""));return 0;}  //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=vasprun|--test=vasprun.xml.relax1"+DEFAULT_KZIP_EXT+"|--test=vasprun.xml.relax2"+DEFAULT_KZIP_EXT+"|--test=vasprun.xml.static"+DEFAULT_KZIP_EXT+"|--test=vasprun.xml.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xVASPRUNXML(aurostd::args2attachedstring(argv,"--test=",""));return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=IBZKPT|--test=IBZKPT.relax1"+DEFAULT_KZIP_EXT+"|--test=IBZKPT.relax2"+DEFAULT_KZIP_EXT+"|--test=IBZKPT.static"+DEFAULT_KZIP_EXT+"|--test=IBZKPT.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xIBZKPT(aurostd::args2attachedstring(argv,"--test=",""));return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=CHGCAR|--test=CHGCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=CHGCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=CHGCAR.static"+DEFAULT_KZIP_EXT+"|--test=CHGCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xCHGCAR(aurostd::args2attachedstring(argv,"--test=",""));return 0;} //CO20180419

    // SCRUB things
    if(!Arun && aurostd::args2flag(argv,cmds,"--scrub=POTCAR")) { 
      XHOST.DEBUG=FALSE;
      XHOST.PSEUDOPOTENTIAL_GENERATOR=TRUE;
      vector<string> vfile(aurostd::args2vectorstring(argv,"--FILE|--file|--F|--f","./"));
      for(uint ifile=0;ifile<vfile.size();ifile++) {
        cerr << "PROCESSING = " << vfile.at(ifile) << endl;
        xPOTCAR xPOT(vfile.at(ifile));
      }
      return 0; //CO20180419
    }
    if(!Arun && aurostd::args2flag(argv,cmds,"--scrub=OUTCAR")) { 
      XHOST.DEBUG=FALSE;
      XHOST.PSEUDOPOTENTIAL_GENERATOR=FALSE;
      vector<string> vfile(aurostd::args2vectorstring(argv,"--FILE|--file|--F|--f","./"));
      for(uint ifile=0;ifile<vfile.size();ifile++) {
        cerr << "PROCESSING = " << vfile.at(ifile) << endl;
        xOUTCAR xOUT(vfile.at(ifile));
        // cout << xOUT << endl;
      }
      return 0; //CO20180419
    }

    if(!Arun && (aurostd::args2flag(argv,cmds,"--scrub") || aurostd::args2attachedflag(argv,"--scrub="))) {
      //  XHOST.DEBUG=TRUE;
      aflowlib::LIB2SCRUB(aurostd::args2attachedstring(argv,"--scrub=","ALL"),TRUE);
      return 0; //CO20180419
    }
    if(!Arun && (aurostd::args2attachedflag(argv,"--lib2auid="))) {
      //  XHOST.DEBUG=TRUE;
      aflowlib::LIB2AUID(aurostd::args2attachedstring(argv,"--lib2auid=","ALL"),FALSE,TRUE); // no test, act
      return 0; //CO20180419
    }
    if(!Arun && (aurostd::args2flag(argv,cmds,"--mosfet") || aurostd::args2attachedflag(argv,"--mosfet="))) {
      //  XHOST.DEBUG=TRUE;
      aflowlib::MOSFET(aurostd::args2attachedutype<int>(argv,"--mosfet=",0),TRUE);
      return 0; //CO20180419
    }
    if(!Arun && (aurostd::args2flag(argv,cmds,"--mail2scan") || aurostd::args2attachedflag(argv,"--mail2scan="))) {
      //  XHOST.DEBUG=TRUE;
      aflowlib::MAIL2SCAN(aurostd::args2attachedstring(argv,"--mail2scan=","/var/mail/auro"),TRUE);
      return 0; //CO20180419
    }

    if(!Arun && aurostd::args2flag(argv,cmds,"--test_proto1")) {
      vector<xstructure> vstr;
      vector<string> vlattice;aurostd::string2tokens("BCC,BCT,CUB,FCC,HEX,MCL,MCLC,ORC,ORCC,ORCF,ORCI,RHL,TET,TRI",vlattice,",");
      aflowlib::_aflowlib_entry data;
      vector<aflowlib::_aflowlib_entry> vdata;
      for(uint i=4;i<vlattice.size();i++) {
        data.clear();
        data.url2aflowlib("materials.duke.edu:AFLOWDATA/ICSD_WEB/"+vlattice.at(i)+"/?format=text",cout,FALSE);vdata.push_back(data);
        cout << "AFLOWLIB " << vlattice.at(i) << "=" << data.vaflowlib_entries.size() << endl;
        for(uint j=0;j<data.vaflowlib_entries.size();j++) {
          aflowlib::_aflowlib_entry dataj;
          dataj.url2aflowlib("materials.duke.edu:AFLOWDATA/ICSD_WEB/"+vlattice.at(i)+"/"+data.vaflowlib_entries.at(j),cout,TRUE);
          aurostd::StringSubst(dataj.aurl,"aflowlib","materials");
          if(dataj.aurl!="") {
            xstructure str(dataj.aurl,"CONTCAR.relax.vasp",IOAFLOW_AUTO);
            xEIGENVAL xEIGENVAL;xEIGENVAL.GetPropertiesUrlFile(dataj.aurl,"EIGENVAL.bands"+DEFAULT_KZIP_EXT+"",FALSE);
            xOUTCAR xOUTCAR;xOUTCAR.GetPropertiesUrlFile(dataj.aurl,"OUTCAR.static"+DEFAULT_KZIP_EXT+"",FALSE);
            xDOSCAR xDOSCAR;xDOSCAR.GetPropertiesUrlFile(dataj.aurl,"DOSCAR.static"+DEFAULT_KZIP_EXT+"",FALSE);
            // if(aurostd::args2flag(argv,cmds,"--vasp")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.vasp",stream);
            // if(aurostd::args2flag(argv,cmds,"--qe")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.qe",stream);
            // if(aurostd::args2flag(argv,cmds,"--abinit")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.abinit",stream);
            // if(aurostd::args2flag(argv,cmds,"--aims")) aurostd::url2stringstream(dataj.aurl+"/CONTCAR.relax.aims",stream);
            vstr.push_back(str);
            cerr << "vstr.size()=" << vstr.size() << "  "
              << "str.atoms.size()=" << str.atoms.size() << "  "
              << "OUTCAR.static"+DEFAULT_KZIP_EXT+".size()=" << xOUTCAR.vcontent.size() << "  "
              << "DOSCAR.static"+DEFAULT_KZIP_EXT+".size()=" << xDOSCAR.vcontent.size() << "  "
              << "EIGENVAL.static"+DEFAULT_KZIP_EXT+".size()=" << xEIGENVAL.vcontent.size() << "  "
              << endl;
            //	  cerr << str << endl;
          }
        }
      }
      return 0; //CO20180419
    }

    if(!Arun && aurostd::args2flag(argv,cmds,"--testJ")) {Arun=TRUE;PERFORM_TESTJ(cout);}
    // [OBSOLETE] if(!Arun && aurostd::args2flag(argv,cmds,"--test1")) {Arun=TRUE;PERFORM_TEST1(cout);}
    if(!Arun && aurostd::args2flag(argv,cmds,"--test3")) {Arun=TRUE;PERFORM_TEST3(cout);}
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_slab|--slab_test")) {return (slab::slabTest()?0:1);}  //CO20190601  //CO20190629 if TRUE(==1), return 0 (normal)
    if(!Arun && XHOST.vflag_control.flag("MACHINE")) {
      //ME20200724
      int code = init::InitMachine(FALSE,argv,cmds,cerr);
      if (code >= 0) return code;
    }

    // **************************************************************
    // INTERCEPT AFLOW
    if(!Arun && XHOST.vflag_control.flag("SWITCH_AFLOW")) {Arun=TRUE;AFLOW_main(argv);}
    //DX
    if(!Arun && XHOST.vflag_control.flag("AFLOWIN_SYM")) {Arun=TRUE;AFLOW_main(argv);} 
    //DX
    if(!Arun && (XHOST.vflag_aflow.flag("CLEAN") || XHOST.vflag_aflow.flag("XCLEAN") || XHOST.AFLOW_RUNDIRflag || XHOST.AFLOW_MULTIflag || XHOST.AFLOW_RUNXflag)) {
      Arun=TRUE;AFLOW_main(argv);
    }
    //  // **************************************************************
    // // INTERCEPT AFLOWLIB
    // MOVED INSIDE PFLOW
    // if(!Arun && XHOST.vflag_control.flag("SWITCH_AFLOWLIB") && !XHOST.vflag_pflow.flag("PROTOS") && !XHOST.vflag_pflow.flag("PROTOS_ICSD") && !XHOST.vflag_pflow.flag("PROTO")) {Arun=TRUE;aflowlib::WEB_Aflowlib_Entry_PHP(argv,cout);}
    // **************************************************************
    // INTERCEPT APENNSY
    if(!Arun && XHOST.vflag_control.flag("SWITCH_APENNSY1")) {Arun=TRUE;Apennsymain(argv,cmds);}
    if(!Arun && XHOST.vflag_control.flag("SWITCH_APENNSY2")) {Arun=TRUE;Apennsymain(argv,cmds);}

    // **************************************************************
    // INTERCEPT aconvasp/aqe/apennsy by title
    if(!Arun && aurostd::substring2bool(XHOST.progname,"aconvasp","convasp")) {Arun=TRUE;pflow::main(argv,cmds);}
    if(!Arun && aurostd::substring2bool(XHOST.progname,"aqe")) {Arun=TRUE;pflow::main(argv,cmds);}
    if(!Arun && aurostd::substring2bool(XHOST.progname,"apennsy")) {Arun=TRUE;Apennsymain(argv,cmds);}

    // **************************************************************
    // intercept commands
    if(!Arun && XHOST.vflag_control.flag("MULTI=SH")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_sh(argv);return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=BZIP2")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("bzip2",argv);return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=BUNZIP2")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("bunzip2",argv);return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=GZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("gzip",argv);return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=GUNZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("gunzip",argv);return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=XZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("xz",argv);return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=XUNZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("xunzip",argv);return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=BZ2XZ")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_bz2xz(argv);return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=GZ2XZ")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_gz2xz(argv);return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=ZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_zip(argv);return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MONITOR")) {Arun=TRUE;AFLOW_monitor(argv);return 0;} //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MONITOR_VASP")) {Arun=TRUE;AFLOW_monitor_VASP();return 0;} //CO20180419
    if(!Arun && XHOST.vflag_control.flag("GETTEMP")) {Arun=TRUE;AFLOW_getTEMP(argv);return 0;} //CO20180419

    // **************************************************************
    // INTERCEPT HELP
    // ME20200921 - Restructured to make web processing easier
    stringstream banner_message;
    if(XHOST.vflag_control.flag("AFLOW_HELP")) {
      banner_message << aflow::Banner("BANNER_BIG") << endl << aflow::Intro_HELP("aflow") << aflow::Banner("BANNER_BIG") << endl;
    } else if(XHOST.vflag_control.flag("AFLOW_EXCEPTIONS")) {
      banner_message << aflow::Banner("BANNER_BIG") << endl << aflow::Banner("EXCEPTIONS") << endl;
    } else if(XHOST.vflag_control.flag("README_AFLOW_LICENSE_GPL3")) {
      banner_message << aflow::License_Preamble_aflow() << endl;
      banner_message << " " << endl;
      banner_message << init::InitGlobalObject("README_AFLOW_LICENSE_GPL3_TXT") << endl;
      banner_message << " " << endl;
      banner_message << "*************************************************************************** " << endl;
    } else if(XHOST.vflag_control.flag("README_AFLOW"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_AFLOW_VERSIONS_HISTORY"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_VERSIONS_HISTORY_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_AFLOW_PFLOW"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_PFLOW_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_FROZSL"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_FROZSL_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_APL"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_APL_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_QHA"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_QHA_SCQHA_QHA3P_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_AAPL"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_APL_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_AGL"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_AGL_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_AEL"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_AEL_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_ANRL"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_ANRL_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_COMPARE"))  { //CO20190401
      banner_message << init::InitGlobalObject("README_AFLOW_COMPARE_TXT") << endl; //CO20190401
    } else if(XHOST.vflag_control.flag("README_GFA"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_GFA_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_SYMMETRY"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_SYM_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_CCE"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_CCE_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_CHULL"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_CHULL_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_PARTIAL_OCCUPATION")) {
      banner_message << init::InitGlobalObject("README_AFLOW_POCC_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_APENNSY"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_APENNSY_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_SCRIPTING"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_SCRIPTING_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_EXCEPTIONS"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_EXCEPTIONS_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_XAFLOW"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_XAFLOW_TXT") << endl;
    } else if(XHOST.vflag_control.flag("README_AFLOWRC"))  {
      banner_message << init::InitGlobalObject("README_AFLOW_AFLOWRC_TXT") << endl;
    }

    if (!banner_message.str().empty()) {
      std::cout << (XHOST.vflag_control.flag("WWW")?aurostd::text2html(banner_message.str()):banner_message.str()) << std::endl;
      return 0;
    }

    // **************************************************************
    // PHP-WEB AND CURRICULUM AND HIGH-THROUGHPUT STUFF
    if (ProcessPhpLatexCv()) return 0;
    //  ProcessSecurityOptions(argv,cmds); OLD STUFF AFLOW SECURITY

    // **************************************************************
    bool VVERSION=aurostd::args2flag(argv,cmds,"-v|--version");
    //ME20200921 - Added web mode
    if(!Arun && VVERSION)  {
      // look for version IMMEDIATELY //CO20180419
      Arun=TRUE;
      cout << (XHOST.vflag_control.flag("WWW")?aurostd::text2html(aflow::Banner("AFLOW_VERSION")):aflow::Banner("AFLOW_VERSION"));
      return 0;
    }
    if(!Arun && XHOST.TEST) { Arun=TRUE;cerr << "test" << endl;return 0;} //CO20180419

    // [OBSOLETE]  if(!Arun && (aurostd::substring2bool(XHOST.progname,"aflow1") || aurostd::substring2bool(XHOST.progname,"aflowd1"))) {
    // [OBSOLETE]  Arun=TRUE;AFLOW_main1(argv,cmds);}
    if(!Arun && XHOST.argv.size()==1 && (aurostd::substring2bool(XHOST.progname,"aflow")  || aurostd::substring2bool(XHOST.progname,"aflowd"))) {   
      //   Arun=TRUE;AFLOW_main(argv);
      Arun=TRUE;
      //    cout << "******************************************************************************************************" << endl;
      //  cout << aflow::Banner("BANNER_TINY") << endl;
      cout << aflow::Banner("BANNER_BIG") << endl;
      cout << aflow::Intro_aflow("aflow") << endl;
      cout << pflow::Intro_pflow("aflow") << endl;
      cout << aflow::Intro_sflow("aflow") << endl;
      cout << aflow::Intro_HELP("aflow") << endl;
      cout << aflow::Banner("BANNER_BIG") << endl;
      //    cout << "XHOST.argv.size()=" << XHOST.argv.size()<< endl;
      //    cout << "******************************************************************************************************" << endl;
    }

    // **************************************************************
    // LAST RESOURCE PFLOW
    if(!Arun) { 
      Arun=TRUE;
      return_code = pflow::main(argv,cmds);   //ME20200901 - use pflow::main return code for database handling
    }
    // **************************************************************
    // END
    return (Arun?return_code:1); //Arun==TRUE is 1, so flip because return 0 is normal  //CO20190629 - more explicit return 0//ME20200901 - use return_code
  }
  //CO20180729 - OBSOLETE - use xerror
  //[OBSOLETE]catch(AFLOWRuntimeError& re){
  //[OBSOLETE]  pflow::logger(_AFLOW_FILE_NAME_, soliloquy, "AFLOWRuntimeError detected. Report on the AFLOW Forum: aflow.org/forum.", oss, _LOGGER_ERROR_);
  //[OBSOLETE]  pflow::logger(_AFLOW_FILE_NAME_, re.where(), re.what(), oss, _LOGGER_ERROR_);
  //[OBSOLETE]  return 1;
  //[OBSOLETE]}
  //[OBSOLETE]catch(AFLOWLogicError& le){
  //[OBSOLETE]  pflow::logger(_AFLOW_FILE_NAME_, soliloquy, "AFLOWLogicError detected. Adjust your inputs accordingly.", oss, _LOGGER_ERROR_);
  //[OBSOLETE]  pflow::logger(_AFLOW_FILE_NAME_, le.where(), le.what(), oss, _LOGGER_ERROR_);
  //[OBSOLETE]  return 1;
  //[OBSOLETE]}
  catch (aurostd::xerror& excpt) {
    pflow::logger(excpt.whereFileName(), excpt.whereFunction(), excpt.error_message, oss, _LOGGER_ERROR_);
    return excpt.error_code;
  }
}


// ***************************************************************************
// AFLOW_main
// ***************************************************************************
int AFLOW_main(vector<string> &argv) {
  if(!XHOST.QUIET) cout << aflow::Banner("INTRODUCTION");// << endl;
  KBIN::KBIN_Main(argv);
  // if(!XHOST.QUIET) cout << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << "  " << endl;
  return 0; //1; //CO20180419 - return 0 is normal
}

// ***************************************************************************
// aflow::License_aflow
// ***************************************************************************
namespace aflow {
  string License_Preamble_aflow(void) {
    //(C) 2003-2021 Stefano Curtarolo, MIT-Duke University stefano@duke.edu
    stringstream strstream;
    strstream << endl;
    strstream << "***************************************************************************" << endl;
    strstream << "*                                                                         *" << endl;
    strstream << "*                    AFLOW - Duke University 2003-2021                    *" << endl; //CO20200502 - SC -> AFLOW consortium
    strstream << "*                                                                         *" << endl;
    strstream << "***************************************************************************" << endl;
    strstream << "Copyright 2003-2021 - AFLOW.ORG consortium" << endl;  //CO20200502 - SC -> AFLOW consortium
    strstream << endl;
    strstream << "This file is part of AFLOW software." << endl;
    strstream << endl;
    strstream << "AFLOW is free software: you can redistribute it and/or modify" << endl;
    strstream << "it under the terms of the GNU General Public License as published by" << endl;
    strstream << "the Free Software Foundation, either version 3 of the License, or" << endl;
    strstream << "(at your option) any later version." << endl;
    strstream << endl;
    strstream << "This program is distributed in the hope that it will be useful," << endl;
    strstream << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
    strstream << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
    strstream << "GNU General Public License for more details." << endl;
    strstream << endl;
    strstream << "You should have received a copy of the GNU General Public License" << endl;
    strstream << "along with this program.  If not, see <http://www.gnu.org/licenses/>." << endl;
    strstream << endl;
    strstream << "***************************************************************************" << endl;
    strstream << endl;
    return strstream.str();
  }
} // namespace aflow


// ***************************************************************************
// aflow::Intro_aflow
// ***************************************************************************
// patched by CO20200106 to avoid long string construction (patching for indents)
namespace aflow {
  string Intro_aflow(string x) {
    //(C) 2003-2021 Stefano Curtarolo, MIT-Duke University stefano@duke.edu
    stringstream strstream;
    string tab="  ";
    string xspaces="";for(uint i=0;i<x.size();i++){xspaces+=" ";} //spaces size of x
    strstream << endl;
    strstream << "******* BEGIN INFORMATION MODE *********************************************************************" << endl;
    strstream << tab << x << " -h|--help|--readme_aflow   CHECK" << endl;
    strstream << tab << x << " --version|-v                 CHECK" << endl;
    strstream << tab << x << " --machine                      CHECK" << endl;
    strstream << "******* END INFORMATION MODE ***********************************************************************" << endl;
    strstream << endl;
    strstream << "******* BEGIN RUNNING MODE *************************************************************************" << endl;
    strstream << tab << x << " --run|--run=multi|--run=N   COMMANDS for running" << endl;
    strstream << tab << x << " --clean|-c                    COMMAND for cleaning  " << endl;
    strstream << endl;
    strstream << " MODIFIERS" << endl;
    strstream << tab << " --DIRECTORY[=| ]dir|--D[=| ]dir|--d[=| ]dir" << endl;
    strstream << tab << " --FILE[=| ]file|--F[=| ]file|--f[=| ]file " << endl;
    strstream << tab << " --quiet|-quiet|-q " << endl;
    strstream << tab << " --loop " << endl;
    strstream << tab << " --sort|-sort" << endl;
    strstream << tab << " --reverse|-rsort" << endl;
    strstream << tab << " --random|-rnd" << endl;
    strstream << tab << " --force|-force" << endl;
    strstream << tab << " --mem=XX|--maxmem=XX" << endl;
    strstream << tab << " --readme= xaflow|aflow|aconvasp|aflowrc|scripting|apennsy|apl|agl|ael|anrl|compare|gfa|symmetry|chull|errors|exceptions|frozsl CHECK !!!!" << endl;
    strstream << tab << " --np=NUMBER|--npmax" << endl;
    strstream << tab << " --generate_aflowin_from_vasp" << endl;
    strstream << tab << " --generate_vasp_from_aflowin|--generate" << endl;
    strstream << tab << " --use_aflow.in=XXX" << endl;
    strstream << tab << " --use_LOCK=XXX" << endl;
    strstream << tab << " --use_tmpfs=XXXX " << endl;
    strstream << endl;
    strstream << " MODIFIERS MPI/SERIAL PARAMETERS" << endl;
    strstream << tab << " --mpi|-nompi|--serial " << endl;
    strstream << endl;
    strstream << " HOST ORIENTED OPTION" << endl;
    strstream << tab << " --machine=beta|beta_openmpi|qrats|qflow|x|conrad|eos|materials|habana|aflowlib|ranger|kraken" << endl;
    strstream << tab << "           marylou|parsons|jellium|ohad|host1" << endl;
    strstream << tab << "           raptor --np=N|diamond --np=N" << endl;
    strstream << "******* END RUNNING MODE ***************************************************************************" << endl;
    strstream << endl;
    // --readme=htresources
    return strstream.str();
  }
} // namespace aflow

// ***************************************************************************
// aflow::Intro_sflow
// ***************************************************************************
// patched by CO20200106 to avoid long string construction (patching for indents)
namespace aflow {
  string Intro_sflow(string x) {
    stringstream strstream;
    string tab="  ";
    //string xspaces="";for(uint i=0;i<x.size();i++){xspaces+=" ";} //spaces size of x
    strstream << "******* BEGIN SCRIPTING MODE ***********************************************************************" << endl;
    strstream << " AFLOW SCRIPTING COMMANDS" << endl;
    strstream << tab << x << " --justafter=string" << endl;
    strstream << tab << x << " --justbefore=string" << endl;
    strstream << tab << x << " --justbetween=string_from[,string_to]" << endl;
    strstream << tab << x << " --qsub=N,file" << endl;
    strstream << tab << x << " --qdel=aaa,nnn:mmm,aaa,bbb,ccc" << endl;
    strstream << tab << x << " --bsub=N,file" << endl;
    strstream << tab << x << " --bkill=aaa,nnn:mmm,aaa,bbb,ccc" << endl;
    strstream << tab << x << " --sbatch=N,file" << endl;
    strstream << tab << x << " --scancel=aaa,nnn:mmm,aaa,bbb,ccc" << endl;
    strstream << tab << x << " --kill=aaa,nnn:mmm,aaa,bbb,ccc" << endl;
    strstream << tab << x << " --multi=sh [--np=NUMBER|npmax|nothing] [--F[ILE]] file" << endl;
    strstream << tab << x << " --multi=zip [--prefix PREFIX] [--size SSSS] --F[ILE] file|--D[IRECTORY directory1 directory2 ...." << endl;
    strstream << tab << x << " --multi=bzip2 [--np=NUMBER|npmax|nothing] --F[ILE] file1 file2 file3 ...." << endl;
    strstream << tab << x << " --multi=bunzip2 [--np=NUMBER|npmax|nothing] --F[ILE] file1.bz2 file2.bz2 file3.bz2 ...." << endl;
    strstream << tab << x << " --multi=gzip [--np=NUMBER|npmax|nothing] --F[ILE] file1 file2 file3 ...." << endl;
    strstream << tab << x << " --multi=gunzip [--np=NUMBER|npmax|nothing] --F[ILE] file1.gz file2.gz file3.gz ...." << endl;
    strstream << tab << x << " --multi=xzip [--np=NUMBER|npmax|nothing] --F[ILE] file1 file2 file3 ...." << endl;
    strstream << tab << x << " --multi=xunzip [--np=NUMBER|npmax|nothing] --F[ILE] file1.xz file2.xz file3.xz ...." << endl;
    strstream << tab << x << " --multi=bz2xz [--np=NUMBER|npmax|nothing] --F[ILE] file1.bz2 file2.bz2 file3.bz2 ...." << endl;
    strstream << tab << x << " --multi=gz2xz [--np=NUMBER|npmax|nothing] --F[ILE] file1.gz file2.gz file3.gz ...." << endl;
    strstream << tab << x << " --getTEMP [--runstat|--runbar|--refresh=X|--warning_beep=T|--warning_halt=T|--mem=XX]" << endl;
    strstream << tab << x << " --monitor [--mem=XX]" << endl;
    strstream << "******* END SCRIPTING MODE *************************************************************************" << endl;
    strstream << endl;
    return strstream.str();
  }
} // namespace aflow

// ***************************************************************************
// aflow::Intro_HELP
// ***************************************************************************
// patched by CO20200106 to avoid long string construction (patching for indents)
namespace aflow {
  string Intro_HELP(string x) {
    stringstream strstream;
    string tab="  ";
    string xspaces="";for(uint i=0;i<x.size();i++){xspaces+=" ";} //xspacess size of x
    strstream << "******* BEGIN HELP MODE ****************************************************************************" << endl;
    strstream << " AFLOW HELP AVAILABLE HELPS" << endl;
    strstream << tab << x << " --license" << endl;
    strstream << tab << xspaces << " " << tab << "License information." << endl;
    strstream << tab << x << " --help" << endl;
    strstream << tab << xspaces << " " << tab << "This help." << endl;
    strstream << tab << x << " --readme" << endl;
    strstream << tab << xspaces << " " << tab << "The list of all the commands available." << endl;
    strstream << tab << x << " --readme=xaflow" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the installation of aflow." << endl;
    strstream << tab << x << " --readme=aflow|--readme=run|--readme_aflow" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"running machinery\"." << endl;
    strstream << tab << x << " --readme=aflowrc" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the installation of aflow." << endl;
    strstream << tab << x << " --readme=pflow|--readme=processor|--readme=aconvasp|--readme_aconvasp" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"processing machinery\"." << endl;
    strstream << tab << x << " --readme=scripting|--readme_scripting" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"scripting\" operations." << endl;
    strstream << tab << x << " --readme=apl|--readme_apl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-harmonic-phonon-library\"." << endl;
    strstream << tab << x << " --readme=qha|--readme_qha|--readme=qha3p|--readme_qha3p|--readme=scqha|--readme_scqha" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-quasi-harmonic-library\"." << endl;
    strstream << tab << x << " --readme=aapl|--readme_aapl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-anharmonic-phonon-library (AFLOW-AAPL)\"." << endl;
    strstream << tab << x << " --readme=agl|--readme_agl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-gibbs-library (AFLOW-AGL)\"." << endl;
    strstream << tab << x << " --readme=ael|--readme_ael" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-elastic-library (AFLOW-AEL)\"." << endl;
    strstream << tab << x << " --readme=prototypes|--readme_prototypes|--readme=anrl|--readme_anrl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow library of prototypes\"." << endl;
    strstream << tab << x << " --readme=xtalfinder|--readme_xtalfinder|--readme=compare|--readme_compare" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow-crystal-finder (AFLOW-XtalFinder) code\"." << endl;
    strstream << tab << x << " --readme=gfa|--readme_gfa" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"glass-forming-ability code\"." << endl;
    strstream << tab << x << " --readme=symmetry|--readme_symmetry" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"symmetry library (AFLOW-SYM)\"." << endl;
    strstream << tab << x << " --readme=chull|--readme_chull" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"convex hull library (AFLOW-hull)\"." << endl;
    strstream << tab << x << " --readme=errors|--readme=exceptions|--readme_errors|--readme_exceptions" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for exception handling in AFLOW." << endl;
    strstream << tab << x << " --readme=partial_occupation|--readme=pocc|--readme_pocc" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"partial occupation library\"." << endl;
    strstream << tab << x << " --readme=apennsy|--readme_apennsy" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"apennsy\" binary phase diagram generator." << endl;
    strstream << tab << x << " --readme=frozsl|--readme_frozsl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"frozsl\" add ons." << endl;
    strstream << "******* END HELP MODE ******************************************************************************" << endl;
    strstream << endl;
    return strstream.str();
  }
} // namespace aflow


// ***************************************************************************
// aflow::Banner
// ***************************************************************************
namespace aflow {
  string Banner(string type) {
    stringstream oss;
    if(type=="VERSION" || type=="AFLOW_VERSION") {
      oss << "" << endl;
      oss << "AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-FLOW [(C) "<<XHOST.Copyright_Years<<" aflow.org consortium]" << endl;
      oss << "New versions are available here: <http://" << XHOST.AFLOW_MATERIALS_SERVER << "/AFLOW/>" << endl;
      oss << "" << endl;
      oss << "AFLOW is free software: you can redistribute it and/or modify it under the terms of the" << endl;
      oss << "GNU General Public License as published by the Free Software Foundation, either version 3" << endl;
      oss << "of the License, or (at your option) any later version." << endl;
      oss << "" << endl;
      oss << "This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;" << endl;
      oss << "without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl;
      oss << "See the GNU General Public License for more details." << endl;
      oss << "" << endl;
      oss << "You should have received a copy of the GNU General Public License along with this program." << endl;
      oss << "If not, see <http://www.gnu.org/licenses/>." << endl;
      oss << "" << endl;
      oss << "AFLOW V" << string(AFLOW_VERSION) << " [" << XHOST.hostname << "] [" << XHOST.machine_type << "] ["<< XHOST.CPU_Cores << "] [" << XHOST.Find_Parameters << "]" << endl;
      //     oss << endl;
      return oss.str();
    }
    if(type=="INTRODUCTION") {
      stringstream oss;
      oss << "MMMMM  AFLOW VERSION " << string(AFLOW_VERSION) << " Automatic-FLOW for materials discovery - [" << TODAY << "] -" << aflow_get_time_string() << endl;
      //    oss << "MMMMM  AFLOW VERSION " << aurostd::PaddedPOST(string(AFLOW_VERSION),5) <<" - BUILT ["<<TODAY<<"] - (C) " << XHOST.Copyright_Years << "    " << endl;
      oss << "MMMMM  AFLOW.org consortium - High-Throughput ab-initio Computing Project - (C) " << XHOST.Copyright_Years << "    " << endl;
      oss << "MMMMM  ";// << endl;
      oss << "MMMMM  AFLOW is free software: you can redistribute it and/or modify it under the terms of the" << endl;
      oss << "MMMMM  GNU General Public License as published by the Free Software Foundation, either version 3" << endl;
      oss << "MMMMM  of the License, or (at your option) any later version." << endl;
      oss << "MMMMM  " << endl;
      oss << "MMMMM  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;" << endl;
      oss << "MMMMM  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl;
      oss << "MMMMM  See the GNU General Public License for more details." << endl;
      oss << "MMMMM  " << endl;
      oss << "MMMMM  You should have received a copy of the GNU General Public License along with this program." << endl;
      oss << "MMMMM  If not, see <http://www.gnu.org/licenses/>." << endl;
      oss << "MMMMM " << endl;
      return oss.str();
    }
    if(type=="BANNER_NORMAL") {
      oss << "****************************************************************************************************" << endl;
      oss << "MMMMM  AFLOW V" << string(AFLOW_VERSION) << " Automatic-FLOW [" << TODAY << "] - " << endl; // << aflow_get_time_string() << endl; //CO
      oss << "****************************************************************************************************" << endl;
      return oss.str();
    }
    if(type=="BANNER_BIG") {
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*                          aflow - Automatic-FLOW for materials discovery                          *" << endl;
      oss << "*                aflow.org consortium - High-Throughput ab-initio Computing Project                *" << endl;
      oss << "*" << aurostd::PaddedCENTER(string("version "+string(AFLOW_VERSION)+" - g++/gcc "+aurostd::utype2string(__GNUC__)+"."+aurostd::utype2string(__GNUC_MINOR__)+"."+aurostd::utype2string(__GNUC_PATCHLEVEL__)+" - built ["+string(TODAY)+"] - (C) " +XHOST.Copyright_Years),100) << "*" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*    AFLOW is free software: you can redistribute it and/or modify it under the terms of the       *" << endl;
      oss << "*    GNU General Public License as published by the Free Software Foundation, either version 3     *" << endl;
      oss << "*    of the License, or (at your option) any later version.                                        *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;     *" << endl;
      oss << "*    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     *" << endl;
      oss << "*    See the GNU General Public License for more details.                                          *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*    You should have received a copy of the GNU General Public License along with this program.    *" << endl;
      oss << "*    If not, see <http://www.gnu.org/licenses/>.                                                   *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*     Use of AFLOW software and repositories welcomes references to the following publications:    *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*  Friedrich et al. npj Comput. Mater. 5, 59 (2019)  10.1038/s41524-019-0192-1       (CCE)         *" << endl;
      oss << "*  Hicks et al.     Comp. Mat. Sci. 161, S1 (2019)   10.1016/j.commatsci.2018.10.043 (ANRL proto2) *" << endl;
      oss << "*  Oses et al.      J. Chem. Inf. Model. (2018)      10.1021/acs.jcim.8b00393        (AFLOW-CHULL) *" << endl;
      oss << "*  Gossett et al.   Comp. Mat. Sci. 152, 134 (2018)  10.1016/j.commatsci.2018.03.075 (AFLOW-ML)    *" << endl;
      oss << "*  Hicks et al.     Acta Cryst. A74, 184-203 (2018)  10.1107/S2053273318003066       (AFLOW-SYM)   *" << endl;
      oss << "*  MBNardelli et al Comp. Mat. Sci. 143, 462 (2018)  10.1016/j.commatsci.2017.11.034 (PAOFLOW)     *" << endl;
      oss << "*  Rose et al.      Comp. Mat. Sci. 137, 362 (2017)  10.1016/j.commatsci.2017.04.036 (AFLUX lang)  *" << endl;
      oss << "*  Supka et al.     Comp. Mat. Sci. 136, 76 (2017)   10.1016/j.commatsci.2017.03.055 (AFLOWpi)     *" << endl;
      oss << "*  Plata et al.     npj Comput. Mater. 3, 45 (2017)  10.1038/s41524-017-0046-7       (AAPL kappa)  *" << endl;
      oss << "*  Toher et al.     Phys. Rev.Mater.1, 015401 (2017) 10.1103/PhysRevMaterials.1.015401 (AEL elast) *" << endl;
      oss << "*  Mehl et al.      Comp. Mat. Sci. 136, S1 (2017)   10.1016/j.commatsci.2017.01.017 (ANRL proto1) *" << endl;
      oss << "*  Calderon et al.  Comp. Mat. Sci. 108A, 233 (2015) 10.1016/j.commatsci.2015.07.019 (standard)    *" << endl;
      oss << "*  Toher et al.     Phys. Rev. B 90, 174107 (2014)   10.1103/PhysRevB.90.174107      (AGL Gibbs)   *" << endl;
      oss << "*  Taylor et al.    Comp. Mat. Sci. 93, 178 (2014)   10.1016/j.commatsci.2014.05.014 (REST-API)    *" << endl;
      oss << "*  Curtarolo et al. Comp. Mat. Sci. 58, 227 (2012)   10.1016/j.commatsci.2012.02.002 (AFLOW.org)   *" << endl;
      oss << "*  Curtarolo et al. Comp. Mat. Sci. 58, 218 (2012)   10.1016/j.commatsci.2012.02.005 (AFLOW C++)   *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "*" << aurostd::PaddedCENTER(string("aflow/aflow.org - CONTRIBUTORS"),100) << "*" << endl;
      oss << "*  2000-2019 Stefano Curtarolo (aflow); 2002-2004 Dane Morgan (convasp); 2007-2011 Wahyu Setyawan  *" << endl;
      oss << "*  (--rsm --edos --kband --icsd*); 2008-2011 Roman Chepulskyy (--edos --kband  surfaces);          *" << endl;
      oss << "*  2008 Gus Hart (lattice reductions - prototypes); 2009-2011, Ohad Levy (prototypes);             *" << endl;
      oss << "*  2009-2010, Michal Jahnatek (APL); 2010-2013 Shidong Wang (cluster expansion); 2010-2013         *" << endl;
      oss << "*  Richard Taylor (surfaces, apennsy); 2010-2013 Junkai Xue (prototyper); 2010-2013 Kesong Yang    *" << endl;
      oss << "*  (findsym, frozsl, plotband/dos); 2013-2019 Cormac Toher (AGL Debye-Gruneisen, AEL elastic);     *" << endl;
      oss << "*  2013-2019 Frisco Rose (API, Aflux); 2013-2018 Pinku Nath (Quasi-harmonic approximation);        *" << endl;
      oss << "*  2013-2017 Jose J. Plata (AAPL, thermal cond.); 2014-2019 David Hicks (symmetry, structure       *" << endl;
      oss << "*  comparison, prototypes); 2014-2019 Corey Oses (Egap, bader, chull, APL, pocc); 2018-2019 Marco  *" << endl;
      oss << "*  Esters (AAPL, thermal cond.); 2016-2019 Denise Ford (GFA); 2018-2019 Rico Friedrich (CCE);      *" << endl;
      oss << "*                                                                                                  *" << endl;
      oss << "****************************************************************************************************" << endl;
      oss << "*" << aurostd::PaddedCENTER(string("version "+string(AFLOW_VERSION)+" - g++/gcc "+aurostd::utype2string(__GNUC__)+"."+aurostd::utype2string(__GNUC_MINOR__)+"."+aurostd::utype2string(__GNUC_PATCHLEVEL__)+" - built ["+string(TODAY)+"] - (C) " +XHOST.Copyright_Years),100) << "*" << endl;
      oss << "****************************************************************************************************" << endl;
      return oss.str();
    }
    if(type=="BANNER_TINY") {
      // called 63 times
      oss << "AFLOW VERSION "<<AFLOW_VERSION<<":  [aflow.org consortium - 2003-2021] ";
      return oss.str();
    }
    if(type == "EXCEPTIONS") {
      oss << "List of AFLOW exceptions with error codes. See README_AFLOW_EXCEPTIONS.TXT for more information" << endl;
      oss << endl;
      oss << "----------------------------------------------------------------------------------------------------" << endl;
      oss << "Error Code      Error Type             Error                              Name of Constant          " << endl;
      oss << "----------------------------------------------------------------------------------------------------" << endl;
      oss << "        1       N/A                    Generic error                      _GENERIC_ERROR_           " << endl;
      oss << "        2                              Illegal error code                 _ILLEGAL_CODE_            " << endl;
      oss << "       10       Input Error            generic                            _INPUT_ERROR_             " << endl;
      oss << "       11                              unknown flag                       _INPUT_UNKNOWN_           " << endl;
      oss << "       12                              missing flag                       _INPUT_MISSING_           " << endl;
      oss << "       13                              input ambiguous                    _INPUT_AMBIGUOUS_         " << endl;
      oss << "       14                              illegal parameter                  _INPUT_ILLEGAL_           " << endl;
      oss << "       15                              number of parameters               _INPUT_NUMBER_            " << endl;
      oss << "       20       File Error             generic                            _FILE_ERROR_              " << endl;
      oss << "       21                              file not found                     _FILE_NOT_FOUND_          " << endl;
      oss << "       22                              wrong format                       _FILE_WRONG_FORMAT_       " << endl;
      oss << "       23                              file corrupt                       _FILE_CORRUPT_            " << endl;
      oss << "       30       Value Error            generic                            _VALUE_ERROR_             " << endl;
      oss << "       31                              illegal value                      _VALUE_ILLEGAL_           " << endl;
      oss << "       32                              out of range                       _VALUE_RANGE_             " << endl;
      oss << "       40       Index Error            generic                            _INDEX_ERROR_             " << endl;
      oss << "       41                              illegal value                      _INDEX_ILLEGAL_           " << endl;
      oss << "       42                              out of bounds                      _INDEX_BOUNDS_            " << endl;
      oss << "       43                              mismatch                           _INDEX_MISMATCH_          " << endl;
      oss << "       50       Runtime Error          generic                            _RUNTIME_ERROR_           " << endl;
      oss << "       51                              not initialized                    _RUNTIME_INIT_            " << endl;
      oss << "       52                              SQL error                          _RUNTIME_SQL_             " << endl;
      oss << "       53                              busy                               _RUNTIME_BUSY_            " << endl;
      oss << "       54                              external command not found         _RUNTIME_EXTERNAL_MISS_   " << endl;  //CO20200531
      oss << "       55                              external command failed            _RUNTIME_EXTERNAL_FAIL_   " << endl;  //CO20200531
      oss << "       60       Allocation Error       generic                            _ALLOC_ERROR_             " << endl;
      oss << "       61                              could not allocate memory          _ALLOC_ALLOCATE_          " << endl;
      oss << "       62                              insufficient memory                _ALLOC_INSUFFICIENT_      " << endl;
      oss << "----------------------------------------------------------------------------------------------------" << endl;
      return oss.str();
    }
    cerr << XPID << "aflow::Banner type=" << type << " not found..." << endl;
    oss << "aflow::Banner type=" << type << " not found..." << endl;
    return 0; //CO20180419
    return oss.str();
  }
} // namespace aflow

// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *                                                                         *
// ***************************************************************************

// Update Mon Mar  7 14:05:40 EST 2011 (by WSETYAWAN):
// ICSD prototypes of compounds used in aflow are generated using the following
// command:
// xzcat /common/NIST/$1ary.icsd.xz | pflow --icsd_nobrokenbasis |
// pflow --icsd_nopartialocc | pflow --icsd2proto > README_LIBRARY_ICSD$1.TXT


//#include "../AFLOW3_AURO/aflow_auro.cpp"

// ***************************************************************************
#ifndef _AFLOW_AURO_CPP_
namespace aflowlib {
  uint MOSFET(int mode,bool VERBOSE) {
    if(VERBOSE) cerr << XPID << "aflowlib::MOSFET mode=" << mode << endl;
    return 0;
  }
}
namespace aflowlib {
  uint MAIL2SCAN(string library,bool VERBOSE) {
    if(VERBOSE) cerr << XPID << "aflowlib::MAIL2SCAN library=" << library << endl;
    return 0;
  }
}
#endif


// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *                                                                         *
// ***************************************************************************
