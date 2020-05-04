// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *             Stefano Curtarolo - Duke University - 2003-2020             *
// *                                                                         *
// ***************************************************************************
//
//  Copyright 2003-2019 - Stefano Curtarolo - AFLOW.ORG consortium
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
// [OBSOLETE] #include "aflow_aqe.h"

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


bool EgapTest(ostream& oss){ofstream FileMESSAGE;return EgapTest(FileMESSAGE,oss);}  //CO20190520
bool EgapTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy="EgapTest():";
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
  string soliloquy="gcdTest():";
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
  string soliloquy="smithTest():";
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

int main(int _argc,char **_argv) {
  string soliloquy="main():"; //CO20180419
  ostream& oss=cout;  //CO20180419
  try{
    bool LDEBUG=FALSE; // TRUE;
    if(LDEBUG) cerr << "AFLOW-MAIN [1]" << endl;
    std::vector<string> argv(aurostd::get_arguments_from_input(_argc,_argv));
    if(LDEBUG) cerr << "AFLOW-MAIN [2]" << endl;
    std::vector<string> cmds;  
    // MACHINE
    init::InitMachine(FALSE,argv,cmds,cerr);    
    if(LDEBUG || XHOST.DEBUG) cerr << "AFLOW-MAIN [3]" << endl;

    // aurostd::TmpDirectoryCreate("test");
    // cerr << args2flag(argv,"--aaa|--bbb |--ccc") << endl; exit(0);
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
    if(!Arun && aurostd::args2flag(argv,cmds,"--prx|--prx="))  {Arun=TRUE;PERFORM_PRX(cout);}
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_getpp")) {
      if(KBIN::VASP_PseudoPotential_CleanName_TEST()){return 0;}
      return 1;
    }
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=Init|--test=init")) {
      init::InitLoadString("vLIBS",1);
      return 1;
    }

    if(!Arun && aurostd::args2flag(argv,cmds,"--test_xmatrix")) { //CO20190911
      string soliloquy="test_xmatrix()::";
      bool LDEBUG=TRUE; // TRUE;
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
      //exit(0);
      return 0; //CO20180419
    }
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_Egap|--Egap_test")) {return (EgapTest()?0:1);}  //CO20190601
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_gcd|--gcd_test")) {return (gcdTest()?0:1);}  //CO20190601
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_smith|--smith_test")) {return (smithTest()?0:1);}  //CO20190601
    if(!Arun && aurostd::args2flag(argv,cmds,"--test")) {
      
      if(XHOST.vext.size()!=XHOST.vcat.size()) { cerr << "ERROR - aflow.cpp:main: XHOST.vext.size()!=XHOST.vcat.size(), aborting." << endl; exit(0); }

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
      //exit(0);
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
      //exit(0);
      return 0; //CO20180419
    }

    if(!Arun && aurostd::args2attachedflag(argv,"--bin2base64=")) {
      string b64String;
      aurostd::bin2base64(aurostd::args2attachedstring(argv,"--bin2base64=",""),b64String);
      cout << b64String << endl;
      exit(0);
    }
    
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=POTCAR|--test=POTCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=POTCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=POTCAR.static"+DEFAULT_KZIP_EXT+"|--test=POTCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xPOTCAR(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=DOSCAR|--test=DOSCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=DOSCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=DOSCAR.static"+DEFAULT_KZIP_EXT+"|--test=DOSCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xDOSCAR(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=EIGENVAL|--test=EIGENVAL.relax1"+DEFAULT_KZIP_EXT+"|--test=EIGENVAL.relax2"+DEFAULT_KZIP_EXT+"|--test=EIGENVAL.static"+DEFAULT_KZIP_EXT+"|--test=EIGENVAL.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xEIGENVAL(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=OUTCAR|--test=OUTCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=OUTCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=OUTCAR.static"+DEFAULT_KZIP_EXT+"|--test=OUTCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xOUTCAR(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=KPOINTS|--test=KPOINTS.relax1"+DEFAULT_KZIP_EXT+"|--test=KPOINTS.relax2"+DEFAULT_KZIP_EXT+"|--test=KPOINTS.static"+DEFAULT_KZIP_EXT+"|--test=KPOINTS.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xKPOINTS(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=vasprun|--test=vasprun.xml.relax1"+DEFAULT_KZIP_EXT+"|--test=vasprun.xml.relax2"+DEFAULT_KZIP_EXT+"|--test=vasprun.xml.static"+DEFAULT_KZIP_EXT+"|--test=vasprun.xml.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xVASPRUNXML(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=IBZKPT|--test=IBZKPT.relax1"+DEFAULT_KZIP_EXT+"|--test=IBZKPT.relax2"+DEFAULT_KZIP_EXT+"|--test=IBZKPT.static"+DEFAULT_KZIP_EXT+"|--test=IBZKPT.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xIBZKPT(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} //CO20180419
    if(!Arun && aurostd::args2flag(argv,cmds,"--test=CHGCAR|--test=CHGCAR.relax1"+DEFAULT_KZIP_EXT+"|--test=CHGCAR.relax2"+DEFAULT_KZIP_EXT+"|--test=CHGCAR.static"+DEFAULT_KZIP_EXT+"|--test=CHGCAR.bands"+DEFAULT_KZIP_EXT+"")) {
      XHOST.DEBUG=TRUE;xCHGCAR(aurostd::args2attachedstring(argv,"--test=",""));/*exit(0)*/return 0;} //CO20180419

    // SCRUB things
    if(!Arun && aurostd::args2flag(argv,cmds,"--scrub=POTCAR")) { 
      XHOST.DEBUG=FALSE;
      XHOST.PSEUDOPOTENTIAL_GENERATOR=TRUE;
      vector<string> vfile(aurostd::args2vectorstring(argv,"--FILE|--file|--F|--f","./"));
      for(uint ifile=0;ifile<vfile.size();ifile++) {
	cerr << "PROCESSING = " << vfile.at(ifile) << endl;
	xPOTCAR xPOT(vfile.at(ifile));
      }
      exit(0);
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
      exit(0);
    }

    if(!Arun && (aurostd::args2flag(argv,cmds,"--scrub") || aurostd::args2attachedflag(argv,"--scrub="))) {
      //  XHOST.DEBUG=TRUE;
      aflowlib::LIB2SCRUB(aurostd::args2attachedstring(argv,"--scrub=","ALL"),TRUE);
      exit(0);
    }
    if(!Arun && (aurostd::args2attachedflag(argv,"--lib2auid="))) {
      //  XHOST.DEBUG=TRUE;
      aflowlib::LIB2AUID(aurostd::args2attachedstring(argv,"--lib2auid=","ALL"),FALSE,TRUE); // no test, act
      exit(0);
    }
    if(!Arun && (aurostd::args2flag(argv,cmds,"--mosfet") || aurostd::args2attachedflag(argv,"--mosfet="))) {
      //  XHOST.DEBUG=TRUE;
      aflowlib::MOSFET(aurostd::args2attachedutype<int>(argv,"--mosfet=",0),TRUE);
      exit(0);
    }
    if(!Arun && (aurostd::args2flag(argv,cmds,"--mail2scan") || aurostd::args2attachedflag(argv,"--mail2scan="))) {
      //  XHOST.DEBUG=TRUE;
      aflowlib::MAIL2SCAN(aurostd::args2attachedstring(argv,"--mail2scan=","/var/mail/auro"),TRUE);
      exit(0);
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
            xEIGENVAL xEIGENVAL; xEIGENVAL.GetPropertiesUrlFile(dataj.aurl,"EIGENVAL.bands"+DEFAULT_KZIP_EXT+"",FALSE);
            xOUTCAR xOUTCAR; xOUTCAR.GetPropertiesUrlFile(dataj.aurl,"OUTCAR.static"+DEFAULT_KZIP_EXT+"",FALSE);
            xDOSCAR xDOSCAR; xDOSCAR.GetPropertiesUrlFile(dataj.aurl,"DOSCAR.static"+DEFAULT_KZIP_EXT+"",FALSE);
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
      //exit(0);
      return 0; //CO20180419
    }

    if(!Arun && aurostd::args2flag(argv,cmds,"--testJ")) {Arun=TRUE;PERFORM_TESTJ(cout);}
    // [OBSOLETE] if(!Arun && aurostd::args2flag(argv,cmds,"--test1")) {Arun=TRUE;PERFORM_TEST1(cout);}
    if(!Arun && aurostd::args2flag(argv,cmds,"--test3")) {Arun=TRUE;PERFORM_TEST3(cout);}
    if(!Arun && aurostd::args2flag(argv,cmds,"--test_slab|--slab_test")) {return (slab::slabTest()?0:1);}  //CO20190601  //CO20190629 if TRUE(==1), return 0 (normal)
    if(!Arun && XHOST.vflag_control.flag("MACHINE"))  {Arun=TRUE;init::InitMachine(TRUE,argv,cmds,cout);}

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
    if(!Arun && XHOST.vflag_control.flag("MULTI=SH")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_sh(argv);/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=BZIP2")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("bzip2",argv);/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=BUNZIP2")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("bunzip2",argv);/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=GZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("gzip",argv);/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=GUNZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("gunzip",argv);/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=XZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("xz",argv);/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=XUNZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_compress("xunzip",argv);/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=BZ2XZ")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_bz2xz(argv);/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=GZ2XZ")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_gz2xz(argv);/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MULTI=ZIP")) {Arun=TRUE;AFLOW_PTHREADS::MULTI_zip(argv);/*exit(0)*/return 0;}  //CO20180419
    if(!Arun && XHOST.vflag_control.flag("MONITOR")) {Arun=TRUE;AFLOW_monitor(argv);/*exit(0)*/return 0;} //CO20180419
    if(!Arun && XHOST.vflag_control.flag("GETTEMP")) {Arun=TRUE;AFLOW_getTEMP(argv);/*exit(0)*/return 0;} //CO20180419

    // **************************************************************
    // INTERCEPT HELP
    if(XHOST.vflag_control.flag("AFLOW_HELP")) {
      cout << aflow::Banner("BANNER_BIG") << endl << aflow::Intro_HELP("aflow") << aflow::Banner("BANNER_BIG") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("AFLOW_EXCEPTIONS")) {
      cout << aflow::Banner("BANNER_BIG") << endl << aflow::Banner("EXCEPTIONS") << endl;
      return 0;
    }  //ME20180531
    if(XHOST.vflag_control.flag("README_AFLOW_LICENSE_GPL3"))  {
      cout << aflow::License_Preamble_aflow() << endl;
      cout << " " << endl;
      cout << init::InitGlobalObject("README_AFLOW_LICENSE_GPL3_TXT") << endl;
      cout << " " << endl;
      cout << "*************************************************************************** " << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_AFLOW"))  {
      cout << init::InitGlobalObject("README_AFLOW_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_AFLOW_VERSIONS_HISTORY"))  {
      cout << init::InitGlobalObject("README_AFLOW_VERSIONS_HISTORY_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_AFLOW_PFLOW"))  {
      cout << init::InitGlobalObject("README_AFLOW_PFLOW_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_FROZSL"))  {
      cout << init::InitGlobalObject("README_AFLOW_FROZSL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_APL"))  {
      cout << init::InitGlobalObject("README_AFLOW_APL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_QHA"))  {
      cout << init::InitGlobalObject("README_AFLOW_QHA_SCQHA_QHA3P_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_AAPL"))  {
      cout << init::InitGlobalObject("README_AFLOW_APL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_AGL"))  {
      cout << init::InitGlobalObject("README_AFLOW_AGL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_AEL"))  {
      cout << init::InitGlobalObject("README_AFLOW_AEL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_ANRL"))  {
      cout << init::InitGlobalObject("README_AFLOW_ANRL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_COMPARE"))  { //CO20190401
      cout << init::InitGlobalObject("README_AFLOW_COMPARE_TXT") << endl; //CO20190401
      /*exit(1)*/return 0;} // << endl; //CO20190401
    if(XHOST.vflag_control.flag("README_GFA"))  {
      cout << init::InitGlobalObject("README_AFLOW_GFA_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_SYMMETRY"))  {
      cout << init::InitGlobalObject("README_AFLOW_SYM_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_CCE"))  {
      cout << init::InitGlobalObject("README_AFLOW_CCE_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_CHULL"))  {
      cout << init::InitGlobalObject("README_AFLOW_CHULL_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_PARTIAL_OCCUPATION")) {
      cout << init::InitGlobalObject("README_AFLOW_POCC_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_APENNSY"))  {
      cout << init::InitGlobalObject("README_AFLOW_APENNSY_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_SCRIPTING"))  {
      cout << init::InitGlobalObject("README_AFLOW_SCRIPTING_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_EXCEPTIONS"))  {
      cout << init::InitGlobalObject("README_AFLOW_EXCEPTIONS_TXT") << endl;
      return 0;}  //ME20180531
    if(XHOST.vflag_control.flag("README_XAFLOW"))  {
      cout << init::InitGlobalObject("README_AFLOW_XAFLOW_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419
    if(XHOST.vflag_control.flag("README_AFLOWRC"))  {
      cout << init::InitGlobalObject("README_AFLOW_AFLOWRC_TXT") << endl;
      /*exit(1)*/return 0;} // << endl; //CO20180419

    // **************************************************************
    // PHP-WEB AND CURRICULUM AND HIGH-THROUGHPUT STUFF
    ProcessPhpLatexCv();
    //  ProcessSecurityOptions(argv,cmds); OLD STUFF AFLOW SECURITY

    // **************************************************************
    bool VVERSION=aurostd::args2flag(argv,cmds,"-v|--version");
    if(!Arun && VVERSION)  {Arun=TRUE; cout << aflow::Banner("AFLOW_VERSION");/*exit(0)*/return 0;} // look for version IMMEDIATELY //CO20180419
    if(!Arun && XHOST.TEST) { Arun=TRUE;cerr << "test" << endl; /*exit(0)*/return 0;} //CO20180419

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
      Arun=TRUE;pflow::main(argv,cmds);
    }
    // **************************************************************
    // END
    return (Arun?0:1); //Arun==TRUE is 1, so flip because return 0 is normal  //CO20190629 - more explicit return 0
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
  return 0; //1;  //CO20180419 - return 0 is normal
}

// ***************************************************************************
// aflow::License_aflow
// ***************************************************************************
namespace aflow {
  string License_Preamble_aflow(void) {
    //( C) 2003-2019 Stefano Curtarolo, MIT-Duke University SC@duke.edu
    stringstream strstream;
    strstream << endl;
    strstream << "***************************************************************************" << endl;
    strstream << "*                                                                         *" << endl;
    strstream << "*           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *" << endl;
    strstream << "*                                                                         *" << endl;
    strstream << "***************************************************************************" << endl;
    strstream << "Copyright 2003-2019 - Stefano Curtarolo - AFLOW.ORG consortium" << endl;
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
    //( C) 2003-2019 Stefano Curtarolo, MIT-Duke University SC@duke.edu
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
    strstream << tab << " --machine=beta|beta_openmpi|qrats|qflow|conrad|eos|materials|habana|aflowlib|ranger|kraken" << endl;
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
    strstream << tab << x << " --readme=anrl|--readme_anrl" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"aflow library of prototypes\"." << endl;
    strstream << tab << x << " --readme=compare|--readme_compare" << endl;
    strstream << tab << xspaces << " " << tab << "Returns the HELP information for the \"structure geometry comparison code\"." << endl;
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
      oss << "AFLOW VERSION "<<AFLOW_VERSION<<":  [aflow.org consortium - 2003-2019] ";
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
      oss << "       60       Allocation Error       generic                            _ALLOC_ERROR_             " << endl;
      oss << "       61                              could not allocate memory          _ALLOC_ALLOCATE_          " << endl;
      oss << "       62                              insufficient memory                _ALLOC_INSUFFICIENT_      " << endl;
      oss << "----------------------------------------------------------------------------------------------------" << endl;
      return oss.str();
    }
    cerr << "aflow::Banner type=" << type << " not found..." << endl;
    oss << "aflow::Banner type=" << type << " not found..." << endl;
    //exit(0);
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
    if(VERBOSE) cerr << "aflowlib::MOSFET mode=" << mode << endl;
    return 0;
  }
}
namespace aflowlib {
  uint MAIL2SCAN(string library,bool VERBOSE) {
    if(VERBOSE) cerr << "aflowlib::MAIL2SCAN library=" << library << endl;
    return 0;
  }
}
#endif


// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *                                                                         *
// ***************************************************************************
