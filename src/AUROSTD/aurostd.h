// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************

#ifndef _AUROSTD_MAIN_H_
#define _AUROSTD_MAIN_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <fstream>
#include <functional>  //ME20220127 - for unit tests and multithreading
#include <grp.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <pthread.h>
#include <pwd.h>
#include <queue>
#include <sstream>
#include <stdarg.h>
#include <stdexcept>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#ifndef __CYGWIN__  //ME20190327 - Cygwin support
//#include <sys/sysctl.h>
#endif
// #include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <sys/wait.h>
#ifdef AFLOW_MULTITHREADS_ENABLE
#include <mutex>
#include <thread>
#endif
#include <time.h>
#include <typeinfo>
#include <unistd.h>
#include <signal.h>  //ME20191125 - needed for AflowDB
#include <vector>
#include <list>          //CO20170806 - need for POCC
#include <utility>       //HE2021069 - for pairs in chull (C++98 changes, already included in SYMBOLICCPLUSPLUS)
#include <netdb.h>       //CO20180321 - frisco needs for AFLUX + for EntryLoader
#include <fts.h>         //HE20220222 - for EntryLoader (effective filesystem tree walk)
#include <regex>         //HE20220222 - for EntryLoader (faster match of complex patterns like alloy matching)



#define GCC_VERSION (__GNUC__ * 10000  + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)  //CO20200502 - moved from aflow.h

//CO20200502 START - including gettid()
#ifdef __GLIBC__
#define GLIBC_VERSION (__GLIBC__ * 100 + __GLIBC_MINOR__)
#if (GLIBC_VERSION < 230) //CO20200502 - apparently they patched at 230 - https://stackoverflow.com/questions/30680550/c-gettid-was-not-declared-in-this-scope
//[CO20200502 - too many warnings]#warning "defining getid() with syscall(SYS_gettid)"
#include <sys/syscall.h>  //CO20200502 - need for gettid()
#define gettid() syscall(SYS_gettid)
#endif  //GLIBC_VERSION
#endif  //__GLIBC__
//CO20200502 END - including gettid()

#ifdef _USE_AFLOW_H_
//#include "aflow.h"
#endif

using std::abs;
using std::cerr;
using std::cin;
using std::cout;
using std::deque;
using std::endl;
using std::ends;
using std::ifstream;
using std::ios_base;
using std::istream;
using std::istringstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::setprecision;
using std::setw;
using std::string;
using std::stringstream;
using std::vector;

#ifndef SWAP
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#endif

#ifndef NNN
#define NNN -123456
#endif

#ifndef AUROSTD_NAN
#define AUROSTD_NAN 1E9
#endif

#ifndef AUROSTD_DEFAULT_PRECISION
#define AUROSTD_DEFAULT_PRECISION 20
#endif

//CO20171215 - more global control
//this is for PRINTING, not for algorithmic zeroing
#ifndef AUROSTD_ROUNDOFF_TOL
#define AUROSTD_ROUNDOFF_TOL 1e-6
#endif

//CO20171215 - more global control
//this is SMALLER than _ZERO_TOL_ in aflow.h 
//(aurostd default is less conservative to match read/write roundoff error)
#ifndef AUROSTD_IDENTITY_TOL
#define AUROSTD_IDENTITY_TOL 1e-6
#endif

//CO20171002 - USEFUL!
#ifndef AUROSTD_MAX_INT
#define AUROSTD_MAX_INT std::numeric_limits<int>::max()
#endif

//CO20171002 - USEFUL!
#ifndef AUROSTD_MAX_UINT
#define AUROSTD_MAX_UINT std::numeric_limits<uint>::max()
#endif

//CO20171002 - USEFUL!
#ifndef AUROSTD_MAX_ULLINT
#define AUROSTD_MAX_ULLINT std::numeric_limits<unsigned long long int>::max()
#endif

//CO20171002 - USEFUL!
#ifndef AUROSTD_MAX_DOUBLE
#define AUROSTD_MAX_DOUBLE std::numeric_limits<double>::max()
#endif

//CO20180101 - stream2stream modes
#ifndef DEFAULT_STREAM
#define DEFAULT_STREAM 'D'
#endif
#ifndef FIXED_STREAM
#define FIXED_STREAM 'F'
#endif
#ifndef SCIENTIFIC_STREAM
#define SCIENTIFIC_STREAM 'S'
#endif


#define AUROSTD_ZIP_BIN "xz"
#define AUROSTD_ZIP_EXT ".xz"

#define __GNU_CPP
#ifndef __xprototype
#ifdef GNU
#define __xprototype __attribute__((const))
#else
#define __xprototype
#endif
#endif

#include "aurostd_xscalar.h"
#include "aurostd_xcomplex.h"
#include "aurostd_xvector.h"
#include "aurostd_xmatrix.h"
#include "aurostd_xtensor.h"
#include "aurostd_xrandom.h"
#include "aurostd_xoption.h"
#include "aurostd_argv.h"
#include "aurostd_xparser.h" //CO20200624
#include "aurostd_xcombos.h"
#include "aurostd_xerror.h" //ME20180627
#include "aurostd_xfit.h" //AS20200824
#include "aurostd_xhttp.h" //HE20220121

using aurostd::min;
using aurostd::max;
using aurostd::mod;
using aurostd::_isodd;
using aurostd::_iseven;
using aurostd::_isfloat;
using aurostd::_iscomplex;
using aurostd::sign;
using aurostd::nint;
using aurostd::xcomplex;
using aurostd::xmatrix;
using aurostd::clear;
using aurostd::xvector;
using aurostd::xtensor; //ME20180627
//[ME20180627 OBSOLETE]using aurostd::xtensor3;
//[ME20180627 OBSOLETE]using aurostd::xtensor4;
//[ME20180627 OBSOLETE]using aurostd::xtensor5;
//[ME20180627 OBSOLETE]using aurostd::xtensor6;
//[ME20180627 OBSOLETE]using aurostd::xtensor7;
//[ME20180627 OBSOLETE]using aurostd::xtensor8;
using aurostd::xoption;

using aurostd::xvector;
using aurostd::xmatrix;
using aurostd::xcombos;

#ifndef uint
typedef unsigned uint;
#endif

// ----------------------------------------------------------------------------

#ifndef __STRICT_ANSI__
#define __STRICT_ANSI__
#endif

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

#ifndef OFF
#define OFF FALSE
#endif
#ifndef ON
#define ON TRUE
#endif


#define _AUROSTD_XLIBS_ERROR_    string("ERROR - AUROSTD_XLIBS++ [***]:")
#define _AUROSTD_XLIBS_WARNING_  string("WARNING - AUROSTD_XLIBS++ [***]:")

#define EMPTY_WORDING string("blablabla")

//extern bool QUIET,DEBUG;
//extern class _XHOST XHOST;
#include "../aflow.h"     //needed for XHOST //SD20220224 - also for _LOCK_LINK_SUFFIX_
//#include "../SQLITE/sqlite3.h"  // OBSOLETE ME20191228 - not used

//CO20200624 START - adding from Jahnatek
//http://ascii-table.com/ansi-escape-sequences.php
//http://ascii-table.com/ansi-escape-sequences-vt-100.php
#define cursor_moveyx(y, x, fstr) fprintf(fstr,"\033[%d;%dH", y, x) //Move cursor to position y,x (rows, columns) with (1,1) as origin    
#define cursor_moveup(y, fstr) fprintf(fstr,"\033[%dA", y)          //Move cursor up y    
#define cursor_movedown(y, fstr) fprintf(fstr,"\033[%dB", y)        //Move cursor down y    
#define cursor_moveright(x, fstr) fprintf(fstr,"\033[%dC", x)       //Move cursor right x    
#define cursor_moveleft(x, fstr) fprintf(fstr,"\033[%dD", x)        //Move cursor left x    
#define cursor_store(fstr) fprintf(fstr,"\033[s")                 //Store current cursor position and color    
#define cursor_restore(fstr) fprintf(fstr,"\033[u")               //Restore cursor position and color from cursor_store()    
#define cursor_clear(fstr) fprintf(fstr,"\033[2J")                //Clear screen and leave cursor where is    
#define cursor_clearline(fstr) fprintf(fstr,"\033[K")             //Clear to end of line and leave cursor where is    
#define cursor_fore_black(fstr) fprintf(fstr,"\033[30m")          //Change foreground color to black    
#define cursor_fore_red(fstr) fprintf(fstr,"\033[31m")            //Change foreground color to red    
#define cursor_fore_green(fstr) fprintf(fstr,"\033[32m")          //Change foreground color to green    
#define cursor_fore_orange(fstr) fprintf(fstr,"\033[33m")         //Change foreground color to orange    
#define cursor_fore_blue(fstr) fprintf(fstr,"\033[34m")           //Change foreground color to blue    
#define cursor_fore_magenta(fstr) fprintf(fstr,"\033[35m")        //Change foreground color to magenta    
#define cursor_fore_cyan(fstr) fprintf(fstr,"\033[36m")           //Change foreground color to cyan    
#define cursor_fore_yellow(fstr) fprintf(fstr,"\033[33m\033[1m")  //Change foreground color to yellow (add bold to help visibility) 
#define cursor_fore_white(fstr) fprintf(fstr,"\033[37m")          //Change foreground color to white    
#define cursor_back_black(fstr) fprintf(fstr,"\033[40m")          //Change background color to black    
#define cursor_back_red(fstr) fprintf(fstr,"\033[41m")            //Change background color to red    
#define cursor_back_green(fstr) fprintf(fstr,"\033[42m")          //Change background color to green    
#define cursor_back_orange(fstr) fprintf(fstr,"\033[43m")         //Change background color to orange    
#define cursor_back_blue(fstr) fprintf(fstr,"\033[44m")           //Change background color to blue    
#define cursor_back_magenta(fstr) fprintf(fstr,"\033[45m")        //Change background color to magenta    
#define cursor_back_cyan(fstr) fprintf(fstr,"\033[46m")           //Change background color to cyan    
#define cursor_back_white(fstr) fprintf(fstr,"\033[47m")          //Change background color to white    
#define cursor_attr_none(fstr) fprintf(fstr,"\033[0m")            //Turn off all cursor attributes    
#define cursor_attr_bold(fstr) fprintf(fstr,"\033[1m")            //Make test bold    
#define cursor_attr_underline(fstr) fprintf(fstr,"\033[4m")       //Underline text    
#define cursor_attr_blink(fstr) fprintf(fstr,"\033[5m")           //Supposed to make text blink, usually bolds it instead    
#define cursor_attr_reverse(fstr) fprintf(fstr,"\033[7m")         //Swap background and foreground colors
//CO20200624 END - adding from Jahnatek

template<class utype> std::ostream& operator<<(std::ostream&,const std::vector<utype>&);// __xprototype;
template<class utype> std::ostream& operator<<(std::ostream&,const std::deque<utype>&);// __xprototype;

// ----------------------------------------------------------------------------
// TIME stuff
namespace aurostd {
  int get_day(void);
  int get_day(const tm& tstruct); //CO20200624
  int get_month(void);
  int get_month(const tm& tstruct);  //CO20200624
  int get_year(void);
  int get_year(const tm& tsruct);  //CO20200624
  void get_offset_utc(int& offset_hours,int& offset_mins); //CO20210601
  void get_offset_utc(const tm& tstruct,int& offset_hours,int& offset_mins);  //CO20210601
  long int get_date(void);
  long int get_date(const tm& tsruct);  //CO20200624
  int get_hour(void);
  int get_hour(const tm& tsruct);  //CO20200624
  int get_min(void);
  int get_min(const tm& tsruct); //CO20200624
  int get_sec(void);
  int get_sec(const tm& tsruct); //CO20200624
  long double get_seconds(void);
  long double get_seconds(long double reference_seconds);
  long double get_delta_seconds(long double& seconds_begin);
  long double get_mseconds(void);
  long double get_mseconds(long double reference_useconds);
  long double get_delta_mseconds(long double& useconds_begin);
  long double get_useconds(void);
  long double get_useconds(long double reference_useconds);
  long double get_delta_useconds(long double& useconds_begin);
  string get_time(void);
  string get_time(const tm& tsruct); //CO20200624
  string get_datetime(bool include_utc_offset=false);
  string get_datetime(const tm& tsruct,bool include_utc_offset=false); //CO20200624
  string get_datetime_formatted(const string& date_delim="/",bool include_time=true,const string& date_time_sep=" ",const string& time_delim=":");  //CO20171215
  string get_datetime_formatted(const tm& tsruct,const string& date_delim="/",bool include_time=true,const string& date_time_sep=" ",const string& time_delim=":");  //CO20171215  //CO20200624
  bool beep(uint=2000,uint=100); // standard values
}
// ----------------------------------------------------------------------------
// threadID stuff
namespace aurostd {
  unsigned long long int getTID(void); //CO20200502 - threadID
}
// ----------------------------------------------------------------------------

//FSIO: file system IO
enum FSIO {  //CO20200624 - for execute2string
  stdout_fsio,
  stderr_fsio,
  stdouterr_fsio
};

namespace aurostd {
  void sizes(void) __xprototype;
  // aflow_aurostd.cpp
  //int Nint(double x);
  //int Sign(const double& x);
  //int SignNoZero(const double& x);
  template<class utype> void aswap(utype &a,utype &b);                     // some dumb algebra
  template<class utype> utype max(const std::vector<utype> vec);           // some dumb algebra
  template<class utype> utype max(const std::deque<utype> vec);            // some dumb algebra
  template<class utype> utype max(const std::vector<vector<utype> >mat);   // some dumb algebra
  template<class utype> utype min(const std::vector<utype> vec);           // some dumb algebra
  template<class utype> utype min(const std::deque<utype> vec);            // some dumb algebra
  template<class utype> utype min(const std::vector<vector<utype> >mat);   // some dumb algebra
  template<class utype> utype sum(const std::vector<utype> vec);           // some dumb algebra
  template<class utype> utype sum(const std::deque<utype> vec);            // some dumb algebra
  template<class utype> utype sum(const std::vector<vector<utype> >mat);   // some dumb algebra
  template<class utype> utype mean(const std::vector<utype> vec);          // some dumb algebra
  template<class utype> utype mean(const std::deque<utype> vec);           // some dumb algebra
  template<class utype> utype mean(const std::vector<vector<utype> >mat);  // some dumb algebra
  template<class utype> vector<utype> reset(std::vector<utype>& v);        // some dumb algebra
  template<class utype> deque<utype> reset(std::deque<utype>& v);          // some dumb algebra
  template<class utype> vector<vector<utype> > reset(std::vector<std::vector<utype> > m);  // some dumb algebra
  template<class utype> vector<utype> clear(std::vector<utype>& v);        // some dumb algebra
  template<class utype> deque<utype> clear(std::deque<utype>& v);          // some dumb algebra
  template<class utype> vector<vector<utype> > clear(std::vector<std::vector<utype> > m);  // some dumb algebra
  template<class utype> void random_shuffle(std::vector<utype>& v);        // some dumb algebra
  template<class utype> void random_shuffle(std::deque<utype>& v);         // some dumb algebra
  template<class utype> bool identical(std::vector<utype> v1,std::vector<utype> v2,utype epsilon);
  template<class utype> bool identical(std::deque<utype> v1,std::deque<utype> v2,utype epsilon);
  bool identical(std::vector<int> v1,std::vector<int> v2,int epsilon);
  bool identical(std::deque<int> v1,std::deque<int> v2,int epsilon);
  template<class utype> bool identical(const std::vector<vector<utype> >& m1,const std::vector<vector<utype> >& m2,utype epsilon);
  template<class utype> bool identical(const vector<utype>& vec, utype eps=(utype)AUROSTD_IDENTITY_TOL); //DX20210422 - checks if all values are the same
  template<class utype> bool identical(const deque<utype>& vec, utype eps=(utype)AUROSTD_IDENTITY_TOL); //DX20210422 - checks if all values are the same
  string toupper(const string& in)  __xprototype;
  string tolower(const string& in)  __xprototype;
  char toupper(const char& in)  __xprototype;
  char tolower(const char& in)  __xprototype;
  string getPWD();  //CO20191112
  int GetNumFields(const string& s);
  string GetNextVal(const string& s,int& id);
  string PaddedNumString(const int num,const int ndigits);
  int getZeroPadding(double num);  //CO20191217
  int getZeroPadding(int num);  //CO20191217
  int getZeroPadding(uint num); //CO20191217
  int getZeroPadding(long int num); //CO20191217
  int getZeroPadding(unsigned long int num);  //CO20191217
  int getZeroPadding(long long int num);  //CO20191217
  int getZeroPadding(unsigned long long int num);  //ME20190108
  template<class utype> string PaddedPRE(utype,int,string=" ");
  string PaddedPRE(string,int,string=" ");
  template<class utype> string PaddedPOST(utype,int,string=" ");
  string PaddedPOST(string,int,string=" ");
  template<class utype> string PaddedCENTER(utype,int,string=" ");
  //write progresses
  string PaddedCENTER(string,int,string=" ");
  uint ProgressBar(std::ostream& oss,string prelim,uint j,uint jmax,bool VERBOSE_PERCENTAGE,bool VERBOSE_ROLLER,bool VERBOSE_CURSOR);
  uint ProgressBar(std::ostream& oss,string prelim,uint j,uint jmax);
  uint ProgressBar(std::ostream& oss,string prelim,double j,bool VERBOSE_PERCENTAGE,bool VERBOSE_ROLLER,bool VERBOSE_CURSOR);
  uint ProgressBar(std::ostream& oss,string prelim,double j);  
  //about cleaning up strings
  bool RemoveControlCodeCharactersFromString(const string& in, string& out); //DX20190516  //CO20190620
  bool RemoveControlCodeCharactersFromStringstream(std::stringstream& ss_in, std::stringstream& ss_out); //DX20190516
  bool RemoveControlCodeCharactersFromFile(const string& directory,const string& filename, bool keep_orig_file=true); //DX20190516
  bool isNullByte(char c); //DX20190131
  string removeNullBytes(string in); //DX20190131
  bool RemoveBinaryCharactersFromFile(const string& directory,const string& filename); //DX20190211 //CO20210315
  string PercentEncodeASCII(const char c) __xprototype; //DX20210706
  string CleanStringASCII(const string& s) __xprototype;
  string CleanStringASCII_20190712(const string& s) __xprototype; //CO20190712
  string CleanStringASCII_20190101(const string& s) __xprototype; //CO20190712
  void CleanStringASCII_InPlace(string& s) __xprototype;  //CO20190712
  string RemoveTrailingCharacter(const string& s,char c); //CO+ME20200825
  void RemoveTrailingCharacter_InPlace(string& s,char c); //CO+ME20200825
  string CGI_StringClean(const string& stringIN) __xprototype;
  string RemoveWhiteSpaces(const string& s) __xprototype;
  string RemoveWhiteSpaces(const string& s, const char toogle) __xprototype;
  string RemoveWhiteSpacesFromTheBack(const string& s) __xprototype;
  string RemoveWhiteSpacesFromTheFront(const string& s) __xprototype;
  string RemoveWhiteSpacesFromTheFrontAndBack(const string& s) __xprototype;
  string RemoveSpaces(const string& s) __xprototype;
  string RemoveSpaces(const string& s, const char toogle) __xprototype;
  string RemoveSpacesFromTheBack(const string& s) __xprototype;
  string RemoveTabs(const string& s) __xprototype;
  string RemoveTabs(const string& s, const char toogle) __xprototype;
  string RemoveTabsFromTheBack(const string& s) __xprototype;
  string RemoveComments(const string& s) __xprototype;
  vector<string> RemoveComments(const vector<string>&) __xprototype;  //ME20190614
  deque<string> RemoveComments(const deque<string>&) __xprototype;  //ME20190614
  string RemoveCharacter(const string& s, const char character) __xprototype;
  void RemoveCharacterInPlace(string& s, const char character) __xprototype;  //CO20190712
  string RemoveCharacterFromTheBack(const string& s, const char character); __xprototype; //DX20190708
  string RemoveCharacterFromTheFront(const string& s, const char character); __xprototype; //DX20190708
  string RemoveCharacterFromTheFrontAndBack(const string& s, const char character); __xprototype; //DX20190708
  string RemoveNumbers(const string& s) __xprototype; //CO20190712
  string RemoveNumbers_20190712(const string& s) __xprototype;  //CO10712
  string RemoveNumbers_20190101(const string& s) __xprototype;  //CO20190712
  void RemoveNumbersInPlace(string& s) __xprototype;  //CO20190712
  string RemoveRounding(const string& s) __xprototype;
  //string RemoveCharacter(const string& s, const char character, const char toogle) __xprototype;
  string RemoveSubStringFirst(const string& str_orig, const string& str_rm) __xprototype;
  void RemoveSubStringFirstInPlace(string& str_orig, const string& str_rm) __xprototype;  //CO20190712
  string RemoveSubString(const string& str_orig, const string& str_rm) __xprototype;
  void RemoveSubStringInPlace(string& str_orig, const string& str_rm) __xprototype; //CO20190712
  double VersionString2Double(const string& version_str); //SD20220331
  vector<string> ProcessPIDs(const string& process,bool user_specific=true); //CO20210315
  vector<string> ProcessPIDs(const string& process,string& output_syscall,bool user_specific=true); //CO20210315
  vector<string> ProcessPIDs(const string& process,const string& pgid,string& output_syscall,bool user_specific=true); //SD20220329
  bool ProcessRunning(const string& process,bool user_specific=true); //CO20210315
  bool ProcessRunning(const string& process,const string& pgid,bool user_specific=true); //SD20220329
  void ProcessKill(const string& process,bool user_specific=true,bool sigkill=true); //CO20210315
  void ProcessKill(const string& process,const string& pgid,bool user_specific=true,bool sigkill=true); //SD20220329
  void ProcessRenice(const string& process,int nvalue,bool user_specific=true); //CO20210315
  // about directories and file existing or not
  bool DirectoryMake(string Directory);
  bool SSH_DirectoryMake(string user, string machine,string Directory);
  bool DirectoryChmod(string chmod_string,string Directory);
  bool SubDirectoryLS(const string& _Directory,vector<string>& vsubd);  //CO20200731
  bool DirectoryLS(const string& Directory,vector<string> &vfiles);
  bool DirectoryLS(const string& Directory,deque<string> &vfiles);
  string dirname(const string& _file);  //CO20210315
  string basename(const string& _file); //CO20210315
  bool DirectoryLocked(string directory,string="LOCK");
  bool DirectorySkipped(string directory);
  bool DirectoryWritable(string directory);
  bool DirectoryUnwritable(string directory);
  string TmpStrCreate(const string& _identifier="",const string& _tmpdir="",bool hidden=false,bool directory=false);  //CO20210624
  string TmpFileCreate(const string& identifier="",const string& tmpdir="",bool hidden=false);  //CO20210315 - empty tmpdir means use XHOST.tmpfs
  string TmpDirectoryCreate(const string& identifier="",const string& tmpdir="",bool hidden=false);  //CO20210315 - empty tmpdir means use XHOST.tmpfs
  string CleanFileName(const string& fileIN);
  string ProperFileName(const string& fileIN);
  bool CopyFile(const string& file_from,const string& file_to);
  bool LinkFile(const string& file_from,const string& file_to);
  bool LinkFileAtomic(const string& file_from,const string& file_to,bool soft=true); //SD20220208
  bool UnlinkFile(const string& file_link); //SD20220208
  //CO START
  bool MatchCompressed(const string& CompressedFileName,const string& FileNameOUT);
  // [OBSOLETE]  bool DecompressFile(const string& CompressedFileName);
  bool efile2tempfile(const string& _FileNameIN, string& FileNameOUT); //CO20180220
  bool efile2tempfile(const string& _FileNameIN, string& FileNameOUT,bool& tempfile_created); //CO20180220
  bool IsCompressed(const string& FileNameIn,string& FileNameOut);
  bool IsCompressed(const string& FileNameIn);
  string GetCompressionExtension(const string& CompressedFileName);
  string GetCatCommand(const string& CompressedFileName); //CO20210315
  //CO END
  bool UncompressFile(const string& FileName,const string& command);  bool UncompressFile(const string& FileName); // with guess  
  bool CompressFile(const string& FileName,const string& command=AUROSTD_ZIP_BIN); // with default
  bool ZIP2ZIP(string dir,string from,string to,bool=TRUE,const string& = "");
  bool BZ2XZ(string dir,bool=TRUE,const string& = "");
  bool GZ2XZ(string dir,bool=TRUE,const string& = "");
  // [OBSOLETE]  bool BunzipFile(const string& FileName); bool BzipFile(const string& FileName);
  // [OBSOLETE]  bool GunzipFile(const string& FileName); bool GzipFile(const string& FileName); 
  // [OBSOLETE]  bool XunzipFile(const string& FileName); bool XzipFile(const string& FileName); 
  // [OBSOLETE]  bool UnzipFile(const string& FileName);  bool ZipFile(const string& FileName);
  bool FileExist(const string& FileName);  bool FileExist(const string& FileName,string &FileNameOut);
  bool EFileExist(const string& FileName); bool EFileExist(const string& FileName,string &FileNameOut);
  unsigned long long int FileSize(const string& FileName);  //ME20191001
  bool GetMemoryUsagePercentage(double& usage_percentage_ram,double& usage_percentage_swap);  //CO20210601
  bool GetMemory(unsigned long long int& free_ram,unsigned long long int& total_ram,unsigned long long int& free_swap,unsigned long long int& total_swap); //CO20210315
  bool FileEmpty(const string& FileName);
  bool FileNotEmpty(const string& FileName);
  bool EFileEmpty(const string& FileName); //CO20190808
  bool EFileNotEmpty(const string& FileName); //CO20190808
  long int GetTimestampModified(const string&);  //ME20180712
  long int SecondsSinceFileModified(const string&);  //CO20210315
  unsigned int getFileCheckSum(const string&, const string&);  //ME20190219
  unsigned int getFletcher32(unsigned short*, size_t);  //ME20190219
  string FileToString(const string& FileName);
  void InFileExistCheck(const string& routine,const string& filename,ifstream& file_to_check);
  bool IsCommandAvailable(const string& command);
  bool IsCommandAvailable(const string& command, string& position);
  bool IsCommandAvailableModify(string& command);
  bool CommandRequired(const string& command);
  bool CommandRequired(const string& command, string& position);
  bool IsExecutableAvailable(const string& executable);
  bool IsExecutableAvailable(const string& executable, string& position);
  bool ExecutableRequired(const string& executable);
  bool ExecutableRequired(const string& executable, string& position);
  // about cleaning, higiene is important
  void StringstreamClean(ostringstream& stream);
  void StringstreamClean(stringstream& stream);
  int FindIfStringInStream(const string& key,std::istream& instream);
  // remove comments

  // about printing
  //[CO20200624 - OBSOLETE]void PrintMessageStream(ofstream& FileERROR,ostringstream& stream,bool quiet);
  //[CO20200624 - OBSOLETE]void PrintMessageStream(std::ostream& FileERROR,ostringstream& stream,bool quiet);
  void PrintANSIEscapeSequence(const aurostd::xoption& color,FILE* fstr);
  void PrintMessageStream(ostringstream& stream,bool quiet,std::ostream& oss=cout); //CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE,ostringstream& stream,bool quiet,std::ostream& oss=cout); //CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE,ostringstream& stream,bool quiet,bool osswrite,std::ostream& oss=cout);
  //[CO20200624 - OBSOLETE]void PrintMessageStream(ostringstream& stream,bool quiet);
  //[CO20200624 - OBSOLETE]void PrintErrorStream(ofstream& FileERROR,ostringstream& stream,bool quiet);
  //[CO20200624 - OBSOLETE]void PrintErrorStream(std::ostream& FileERROR,ostringstream& stream,bool quiet);
  void PrintErrorStream(ostringstream& stream,bool quiet); //CO20200624
  void PrintErrorStream(ofstream& FileERROR,ostringstream& stream,bool quiet); //CO20200624
  void PrintErrorStream(ofstream& FileERROR,ostringstream& stream,bool quiet,bool osswrite);
  //[CO20200624 - OBSOLETE]void PrintErrorStream(ostringstream& stream,bool quiet);
  //[CO20200624 - OBSOLETE]void PrintWarningStream(ofstream& FileERROR,ostringstream& stream,bool quiet);
  //[CO20200624 - OBSOLETE]void PrintWarningStream(std::ostream& FileERROR,ostringstream& stream,bool quiet);
  void PrintWarningStream(ostringstream& stream,bool quiet); //CO20200624
  void PrintWarningStream(ofstream& FileWARNING,ostringstream& stream,bool quiet); //CO20200624
  void PrintWarningStream(ofstream& FileWARNING,ostringstream& stream,bool quiet,bool osswrite);
  //[CO20200624 - OBSOLETE]void PrintWarningStream(ostringstream& stream,bool quiet);

  //[CO20200624 - OBSOLETE]void PrintMessageStream(ofstream& FileERROR,stringstream& stream,bool quiet);
  //[CO20200624 - OBSOLETE]void PrintMessageStream(std::ostream& FileERROR,stringstream& stream,bool quiet);
  void PrintMessageStream(stringstream& stream,bool quiet,std::ostream& oss=cout);  //CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE,stringstream& stream,bool quiet,std::ostream& oss=cout);  //CO20200624
  void PrintMessageStream(ofstream& FileMESSAGE,stringstream& stream,bool quiet,bool osswrite,std::ostream& oss=cout);
  //[CO20200624 - OBSOLETE]void PrintMessageStream(stringstream& stream,bool quiet);
  //[CO20200624 - OBSOLETE]void PrintErrorStream(ofstream& FileERROR,stringstream& stream,bool quiet);
  //[CO20200624 - OBSOLETE]void PrintErrorStream(std::ostream& FileERROR,stringstream& stream,bool quiet);
  void PrintErrorStream(stringstream& stream,bool quiet);  //CO20200624
  void PrintErrorStream(ofstream& FileERROR,stringstream& stream,bool quiet);  //CO20200624
  void PrintErrorStream(ofstream& FileERROR,stringstream& stream,bool quiet,bool osswrite);
  //[CO20200624 - OBSOLETE]void PrintErrorStream(stringstream& stream,bool quiet);
  //[CO20200624 - OBSOLETE]void PrintWarningStream(ofstream& FileERROR,stringstream& stream,bool quiet);
  //[CO20200624 - OBSOLETE]void PrintWarningStream(std::ostream& FileERROR,stringstream& stream,bool quiet);
  void PrintWarningStream(stringstream& stream,bool quiet);  //CO20200624
  void PrintWarningStream(ofstream& FileWARNING,stringstream& stream,bool quiet);  //CO20200624
  void PrintWarningStream(ofstream& FileWARNING,stringstream& stream,bool quiet,bool osswrite);
  //[CO20200624 - OBSOLETE]void PrintWarningStream(stringstream& stream,bool quiet);

  // about executing
  bool execute(ostringstream& command);
  bool execute(stringstream& command);
  bool execute(const string& command);
  bool execute(const vector<string>& vcommand);
  bool execute(const deque<string>& dcommand);
#ifdef _stringcharstar_
  bool execute(char* command);
#endif
  // Execute and report
  string execute2string(ostringstream& command,FSIO fsio=stdout_fsio);  //CO20200624 - added file system IO mode
  string execute2string(stringstream& command,FSIO fsio=stdout_fsio);  //CO20200624 - added file system IO mode
  string execute2string(const string& command,FSIO fsio=stdout_fsio);  //CO20200624 - added file system IO mode
  vector<string> execute2string(const vector<string>& vcommand,FSIO fsio=stdout_fsio);  //CO20200624 - added file system IO mode
  deque<string> execute2string(const deque<string>& dcommand,FSIO fsio=stdout_fsio);  //CO20200624 - added file system IO mode
#ifdef _stringcharstar_
  string execute2string(char* command,FSIO fsio=stdout_fsio);  //CO20200624 - added file system IO mode
#endif
  string CleanCommand4Execute(const string& _command); //CO20200624
  template<class utype> utype execute2utype(ostringstream& command);
  template<class utype> utype execute2utype(stringstream& command);
  template<class utype> utype execute2utype(string command);
  template<class utype> vector<utype> execute2utype(vector<utype> vcommand);
  template<class utype> deque<utype> execute2utype(deque<utype> dcommand);
#ifdef _stringcharstar_
  template<class utype> utype execute2utype(char* command);
#endif
  // about sleeping
  unsigned int Sleep(unsigned int seconds);
  // about extracting from to files
  vector<string> GrepFile(const string& filename,const string& keyword,bool RemoveWS=false,bool RemoveComments=true); //CO20210623
  bool ExtractToFileEXPLICIT(ifstream& FileIN,string FileNameOUTPUT,const string& Keyword);
  bool ExtractToFileEXPLICIT(const string& StringIN,string FileNameOUTPUT,const string& Keyword);
  bool ExtractToFileEXPLICIT(ifstream& FileIN,string FileNameOUTPUT,const string& Keyword_start,const string& Keyword_stop);
  bool ExtractToFileEXPLICIT(const string& StringIN,string FileNameOUTPUT,const string& Keyword_start,const string& Keyword_stop);
  bool ExtractToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,const string& Keyword);
  bool ExtractToStringEXPLICIT(const string& StringIN,string& StringOUTPUT,const string& Keyword);
  bool ExtractToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,const string& Keyword_start,const string& Keyword_stop);
  bool ExtractToStringEXPLICIT(const string& StringIN,string& StringOUTPUT,const string& Keyword_start,const string& Keyword_stop);
  bool ExtractToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,const string& Keyword);
  bool ExtractToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop);
  bool ExtractToStringstreamEXPLICIT(stringstream& StringStreamIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop);
  bool ExtractToStringstreamEXPLICIT(const string& StringIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop);
  bool ExtractToStringstreamEXPLICIT(const string& StringIN,stringstream& StringstreamOUTPUT,const string& Keyword);
  // take the nth //SD20220301
  uint ConvertNegativeIndex(const int index,const uint _size);
  bool ExtractNthToStringstreamEXPLICIT(ifstream &FileIN,stringstream& StringstreamOUTPUT,const string& Keyword,const int index);
  bool ExtractNthToStringstreamEXPLICIT(ifstream &FileIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop,const int index);
  bool ExtractNthToStringstreamEXPLICIT(stringstream& StringStreamIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop,int index);
  bool ExtractNthToStringstreamEXPLICIT(const string& StringIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop,const int index);
  // take the last
  bool ExtractLastToStringstreamEXPLICIT(ifstream &FileIN,stringstream& StringstreamOUTPUT,const string& Keyword);
  bool ExtractLastToStringstreamEXPLICIT(ifstream &FileIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop);
  bool ExtractLastToStringstreamEXPLICIT(stringstream& StringStreamIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop);
  bool ExtractLastToStringstreamEXPLICIT(const string& StringIN,stringstream& StringstreamOUTPUT,const string& Keyword_start,const string& Keyword_stop);
  // take just after
  bool ExtractJustAfterToFileEXPLICIT(ifstream& FileIN,string FileNameOUTPUT,const string& Keyword_start);
  bool ExtractJustAfterToStringEXPLICIT(ifstream& FileIN,string& StringOUTPUT,const string& Keyword_start);
  bool ExtractJustAfterToStringstreamEXPLICIT(ifstream& FileIN,stringstream& StringstreamOUTPUT,const string& Keyword_start);
  bool ExtractJustAfterToStringstreamEXPLICIT(stringstream& StringStreamIN,stringstream& StringstreamOUTPUT,const string& Keyword_start);
  bool ExtractJustAfterToStringstreamEXPLICIT(const string& StringIN,stringstream& StringstreamOUTPUT,const string& Keyword_start);
  bool ExtractJustAfterToStringEXPLICIT(const string& StringIN,string& StringOUTPUT,const string& Keyword_start);
  // about taking in istreams and stringstream and strings
  uint stream2vectorstring(std::istream& istreamIN,vector<string> &vstringout);
  uint stream2vectorstring(std::ifstream& ifstreamIN,vector<string> &vstringout);
  uint stream2vectorstring(stringstream& stringstreamIN,vector<string> &vstringout);
  uint string2vectorstring(const string& stringIN,vector<string> &vstringout,bool consecutive=false,bool trim_edges=false); //CO20170613, defaults to usual string2tokens() behavior
  vector<string> stream2vectorstring(std::istream& istreamIN);
  vector<string> stream2vectorstring(std::ifstream& ifstreamIN);
  vector<string> stream2vectorstring(stringstream& stringstreamIN);
  vector<string> string2vectorstring(const string& stringIN,bool consecutive=false,bool trim_edges=false);  //CO20170613, defaults to usual string2tokens() behavior
  string liststring2string(string="",string="",string="",string="",string="",string="",string="",string="",
      string="",string="",string="",string="",string="",string="",string="",string="",
      string="",string="",string="",string="",string="",string="",string="",string="",
      string="",string="",string="",string="",string="",string="",string="",string="",
      string="",string="",string="",string="",string="",string="",string="",string="",
      string="",string="",string="",string="",string="",string="",string="",string="");
  uint stream2dequestring(std::istream& istreamIN,deque<string> &vstringout);
  uint stream2dequestring(std::ifstream& ifstreamIN,deque<string> &vstringout);
  uint stream2dequestring(stringstream& stringstreamIN,deque<string> &vstringout);
  uint string2dequestring(const string& stringIN,deque<string> &vstringout);
  deque<string> stream2dequestring(std::istream& istreamIN);
  deque<string> stream2dequestring(std::ifstream& ifstreamIN);
  deque<string> stream2dequestring(stringstream& stringstreamIN);
  deque<string> string2dequestring(const string& stringIN);
  // about writing files
  bool string2file(const string& StringOUTPUT,const string& FileNameOUTPUT,const string& mode="");
  bool string2compressfile(const string& commamnd,const string& StringOUTPUT,const string& FileNameOUTPUT,const string& mode="");
  bool string2gzfile(const string& StringOUTPUT,const string& FileNameOUTPUT,const string& mode="");
  bool string2bz2file(const string& StringOUTPUT,const string& FileNameOUTPUT,const string& mode="");
  bool string2xzfile(const string& StringOUTPUT,const string& FileNameOUTPUT,const string& mode="");
  bool stringstream2file(const stringstream& StringstreamOUTPUT,const string& FileNameOUTPUT,const string& mode="");
  bool stringstream2compressfile(const string& commamnd,const stringstream& StringstreamOUTPUT,const string& FileNameOUTPUT,const string& mode="");
  bool stringstream2gzfile(const stringstream& StringstreamOUTPUT,const string& FileNameOUTPUT,const string& mode="");
  bool stringstream2bz2file(const stringstream& StringstreamOUTPUT,const string& FileNameOUTPUT,const string& mode="");
  bool stringstream2xzfile(const stringstream& StringstreamOUTPUT,const string& FileNameOUTPUT,const string& mode="");
  // file2string
  uint file2string(const string& FileNameIN,string& StringIN);  //CO20210624
  uint bz2file2string(const string& FileNameIN,string& StringIN); //CO20210624
  uint gzfile2string(const string& FileNameIN,string& StringIN);  //CO20210624
  uint xzfile2string(const string& FileNameIN,string& StringIN);  //CO20210624
  uint zipfile2string(const string& FileNameIN,string& StringIN); //CO  //CO20210624
  uint efile2string(const string& FileNameIN,string& StringIN);  //CO20191110 //CO20210624
  // file2vectorstring  
  uint file2vectorstring(const string& FileNameIN,vector<string>& vline,bool consecutive=false,bool trim_edges=false); //CO20170613, defaults to usual string2tokens() behavior //CO20210624
  uint bz2file2vectorstring(const string& FileNameIN,vector<string>& vline,bool consecutive=false,bool trim_edges=false);  //CO20170613, defaults to usual string2tokens() behavior //CO20210624
  uint gzfile2vectorstring(const string& FileNameIN,vector<string>& vline,bool consecutive=false,bool trim_edges=false); //CO20170613, defaults to usual string2tokens() behavior //CO20210624
  uint xzfile2vectorstring(const string& FileNameIN,vector<string>& vline,bool consecutive=false,bool trim_edges=false); //CO20170613, defaults to usual string2tokens() behavior //CO20210624
  uint efile2vectorstring(const string& FileNameIN,vector<string>& vline,bool consecutive=false,bool trim_edges=false);  //CO20170613, defaults to usual string2tokens() behavior //CO20210624
  bool vectorstring2file(const vector<string>& vline,string FileNameOUT);
  // file2dequestring  
  uint file2dequestring(const string& FileNameIN,deque<string>& vline);  //CO20210624
  uint bz2file2dequestring(const string& FileNameIN,deque<string>& vline); //CO20210624
  uint gzfile2dequestring(const string& FileNameIN,deque<string>& vline);  //CO20210624
  uint xzfile2dequestring(const string& FileNameIN,deque<string>& vline);  //CO20210624
  uint efile2dequestring(const string& FileNameIN,deque<string>& vline); //CO20210624
  bool dequestring2file(const deque<string>& vline,string FileNameOUT);
  // file2vectorstring overloading with deque
  uint file2vectorstring(const string& FileNameIN,deque<string>& vline);  //CO20210624
  uint bz2file2vectorstring(const string& FileNameIN,deque<string>& vline); //CO20210624
  uint gzfile2vectorstring(const string& FileNameIN,deque<string>& vline);  //CO20210624
  uint xzfile2vectorstring(const string& FileNameIN,deque<string>& vline);  //CO20210624
  uint efile2vectorstring(const string& FileNameIN,deque<string>& vline);   //CO20210624
  bool vectorstring2file(const deque<string>& vline,string FileNameOUT);
  // file2stringstream
  bool file2stringstream(const string& FileNameIN,stringstream& StringstreamIN);  //CO20210624
  bool bz2file2stringstream(const string& FileNameIN,stringstream& StringstreamIN); //CO20210624
  bool gzfile2stringstream(const string& FileNameIN,stringstream& StringstreamIN);  //CO20210624
  bool xzfile2stringstream(const string& FileNameIN,stringstream& StringstreamIN);  //CO20210624
  bool zipfile2stringstream(const string& _FileNameIN,stringstream& StringstreamIN); //CO //CO20210624
  bool efile2stringstream(const string& FileNameIN,stringstream& StringstreamIN); //CO20210624
  // return directly the string
  string file2string(const string& FileNameIN); //CO20210624
  string bz2file2string(const string& FileNameIN);  //CO20210624
  string gzfile2string(const string& FileNameIN); //CO20210624
  string xzfile2string(const string& FileNameIN); //CO20210624
  string zipfile2string(const string& FileNameIN); //CO20210624
  string efile2string(const string& FileNameIN);  //CO20210624
  // reading url to string/stringstream/tokens/vector/deque
  bool url2file(string url,string& fileIN,bool=FALSE)  __xprototype;   // bool = verbose
  bool eurl2string(const string& url,string& stringIN,bool verbose=FALSE);  //CO20200223
  bool url2string(const string& url,string& stringIN,bool=FALSE)  __xprototype;   // bool = verbose
  bool eurl2stringstream(const string& url,stringstream& stringstreamIN,bool=FALSE)  __xprototype;  // bool = verbose  //CO20200223
  bool url2stringstream(const string& url,stringstream& stringstreamIN,bool=FALSE)  __xprototype;  // bool = verbose
  bool eurl2vectorstring(const string& url,vector<string>& vlines,bool=FALSE)  __xprototype;  // bool = verbose  //CO20200223
  bool url2vectorstring(const string& url,vector<string>& vlines,bool=FALSE)  __xprototype;  // bool = verbose
  bool eurl2dequestring(const string& url,deque<string>& vlines,bool=FALSE)  __xprototype;  // bool = verbose  //CO20200223
  bool url2dequestring(const string& url,deque<string>& vlines,bool=FALSE)  __xprototype;  // bool = verbose
  template<typename utype> uint eurl2tokens(const string& url,vector<utype>& tokens,const string& delimiters = " ")  __xprototype; //CO20200223
  template<typename utype> uint url2tokens(const string& url,vector<utype>& tokens,const string& delimiters = " ")  __xprototype;
  template<typename utype> uint eurl2tokens(const string& url,deque<utype>& tokens,const string& delimiters = " ")  __xprototype;  //CO20200223
  template<typename utype> uint url2tokens(const string& url,deque<utype>& tokens,const string& delimiters = " ")  __xprototype;
  string eurl2string(const string& url)  __xprototype; //CO20200223
  string url2string(const string& url)  __xprototype;
  // about istream/ostream
  string ostream2string(std::ostream& oss);
  uint stream2string(std::istream& istreamIN,string &vstringout);
  uint stream2string(std::ifstream& ifstreamIN,string &vstringout);
  uint stream2string(stringstream& stringstreamIN,string &vstringout);
  // about environments
  string getenv2string(const string& str);
  int getenv2int(const string& str);
  uint getenv2uint(const string& str);
  double getenv2double(const string& str);
  // about operating on files
  bool ChmodFile(string chmod_string,string FileNameOUTPUT);
  //CO START
  bool file2directory(const string& _file,const string& destination);
  bool file2directory(const vector<string>& files,const string& destination);
  bool file2file(const string& _file,const string& destination); //CO20171025
  string file2md5sum(const string& _file); //SC20200326
  string file2auid(const string& _file); //SC20200326  
  bool IsDirectory(const string& path);
  bool IsFile(const string& path);
  //CO END
  bool RemoveFile(string FileNameOUTPUT);
  bool RemoveFile(const vector<string>& files); //CO
  bool RemoveDirectory(const string& dir); //CO
  // about getting info from strings
  uint string2tokens(const string& str,vector<string>& tokens,const string& delimiters = " ",bool consecutive=false) __xprototype;  //CO20170613, defaults to usual string2tokens() behavior
  uint string2tokens(const string& str,deque<string>& tokens,const string& delimiters = " ",bool consecutive=false) __xprototype; //CO20170613, defaults to usual string2tokens() behavior
  template<class utype> uint string2tokens(const string& str,std::vector<utype>& tokens,const string& delimiters = " ",bool consecutive=false) __xprototype;  //CO20170613, defaults to usual string2tokens() behavior
  template<class utype> uint string2tokens(const string& str,std::deque<utype>& tokens,const string& delimiters = " ",bool consecutive=false) __xprototype; //CO20170613, defaults to usual string2tokens() behavior
  uint string2tokensAdd(const string& str,vector<string>& tokens,const string& delimiters = " ") __xprototype;
  uint string2tokensAdd(const string& str,deque<string>& tokens,const string& delimiters = " ") __xprototype;
  template<class utype> uint string2tokensAdd(const string& str,std::vector<utype>& tokens,const string& delimiters = " ") __xprototype;
  template<class utype> uint string2tokensAdd(const string& str,std::deque<utype>& tokens,const string& delimiters = " ") __xprototype;

  //[CO20210315 - OBSOLETE use stream2stream()]template<typename typeTo, typename typeFrom> typeTo StringStreamConvert(const typeFrom& from);  //CO20210315 - cleaned up
  //[CO20210315 - OBSOLETE use stream2stream()]template<typename typeFrom> string StringConvert(const typeFrom& from);  //CO20210315 - cleaned up
  //[CO20210315 - not defined]template<typename typeTo, typename typeFrom> typeTo NumberStreamConvert(const typeFrom& from);  //CO20210315 - cleaned up

  // [OBSOLETE]  double string2double(const string& from) __xprototype;
  vector<double> vectorstring2vectordouble(const vector<string>& from); //CO20210315 - cleaned up
  // [OBSOLETE]  long double string2longdouble(const string& from) __xprototype;
  // [OBSOLETE]  int string2int(const string& from) __xprototype;
  string string2string(const string& from) __xprototype;
  template<typename utype> utype string2utype(const string& from, const uint base=10);  //CO20210315 - cleaned up //HE20220324 add base option
  vector<int> vectorstring2vectorint(const vector<string>& from); //CO20210315 - cleaned up
  // [OBSOLETE] uint string2uint(const string& from) __xprototype;
  vector<uint> vectorstring2vectoruint(const vector<string>& from); //CO20210315 - cleaned up
  // [OBSOLETE] long int string2longint(const string& from) __xprototype;
  // [OBSOLETE] float string2float(const string& from) __xprototype;

  vector<float> vectorstring2vectorfloat(const vector<string>& from);  //CO20210315 - cleaned up
  string vectorstring2string(const vector<string>& vstrings);
  string vectorstring2string(const deque<string>& vstrings);

  // [OBSOLETE] string double2string(double from) __xprototype;
  // [OBSOLETE] string double2string(double from,int precision) __xprototype;
  // [OBSOLETE] string longdouble2string(long double from) __xprototype;
  // [OBSOLETE] string int2string(int from) __xprototype;
  // [OBSOLETE] string uint2string(uint from) __xprototype;
  // [OBSOLETE] string longint2string(long int from) __xprototype;
  // [OBSOLETE] string float2string(float from) __xprototype;
  //  string utype2string(const std::basic_string<char, std::char_traits<char>, std::allocator<char> >& from) __xprototype;
  //  string utype2string(std::basic_string<char, std::char_traits<char>, std::allocator<char> > from) __xprototype;
  //string utype2string(string from) __xprototype;
  //  std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > utype2string(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > from) __xprototype;
  //DX20210128 [OBSOLETE - use default arguments] template<typename utype> string utype2string(const utype& from) __xprototype;
  //DX20210128 [OBSOLETE - use default arguments] template<typename utype> string utype2string(const utype& from,int precision) __xprototype;
  template<typename utype> string utype2string(const utype& from,int precision=AUROSTD_DEFAULT_PRECISION,char FORMAT=DEFAULT_STREAM) __xprototype; //DX20201028 - this declaration was missing //DX20210128 - add defaults
  string utype2string(double from,bool roff);
  string utype2string(double from,int precision,bool roff);
  string utype2string(double from,bool roff,double tol);
  string utype2string(double from,int precision,bool roff,double tol);
  string utype2string(double from,bool roff,char FORMAT);
  string utype2string(double from,int precision,char FORMAT,bool roff=false); //CO20200624
  string utype2string(double from,int precision,bool roff,char FORMAT);
  string utype2string(double from,bool roff,double tol,char FORMAT);
  string utype2string(double from,int precision,bool roff,double tol,char FORMAT);
  string bool2string(bool from);
  // string utype2string(const string& from) __xprototype;
  // [OBSOLETE]  template<class u1> string utype2string(u1) __xprototype;
  // [OBSOLETE]  template<class u1, class u2> string utype2string(u1,u2) __xprototype;
  // [OBSOLETE]  template<class u1, class u2, class u3> string utype2string(u1,u2,u3) __xprototype;
  // [OBSOLETE]  template<class u1, class u2, class u3, class u4> string utype2string(u1,u2,u3,u4) __xprototype;

  template<class utype> deque<utype> utypes2deque(utype u1) __xprototype;
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2) __xprototype;
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2,utype u3) __xprototype;
  template<class utype> deque<utype> utypes2deque(utype u1,utype u2,utype u3,utype u4) __xprototype;

  void StringCommasColumsVectorInt(string vstring,vector<int> &vint) __xprototype;
  void StringCommasColumsVectorUnsignedInt(string vstring,vector<uint> &vuint) __xprototype;
  void StringCommasColumsVectorFloat(string vstring,vector<float> &vfloat) __xprototype;
  void StringCommasColumsVectorDouble(string vstring,vector<double> &vdouble) __xprototype;
  int GetNLinesString(const string& strstream) __xprototype;
  int GetNLinesString(const stringstream& strstream) __xprototype;
  int GetNLinesFile(const string& file_name) __xprototype;
  string GetLineString(const string& strstream,int line);
  string GetLineString(const stringstream& strstream,int line);
  // substitute strings in strings and stringstreams
  bool StringsAlphabetic(const string& A,const string& B,bool allow_identical=true); //CO20180801
  bool StringsAlphabetic(const vector<string>& input,bool allow_identical=true);  //CO20180801
  bool StringsAlphabetic(const deque<string>& input,bool allow_identical=true);  //CO20180801
  string StringSubst(string &strstring, const string &strfind, const string &strreplace);
  string StringSubst(const string &strstring, const string &strfind, const string &strreplace); //HE20220321
  //  string StringSubst(string &strstring, const string &strfind0, const string &strfind1, const string &strfind2, const string &strfind3, const string &strreplace);
  string StringSubst(string &strstring, const char &charfind, const char &charreplace);
  string StringSubst(const string &strstring, const char &charfind, const char &charreplace);
  void StringStreamSubst(stringstream &strstring, const string &strfind, const string &strreplace);  //ME20190128 - fixed type declaration
  // about present substrings
  bool substring2bool(const string& strstream,const string& strsub1,bool RemoveWS=false,bool RemoveComments=true);  //CO20210315 - cleaned up
  bool substring2bool(const vector<string>& vstrstream,const string& strsub1,bool RemoveWS=false,bool RemoveComments=true); //CO20210315 - cleaned up
  bool substring2bool(const deque<string>& vstrstream,const string& strsub1,bool RemoveWS=false,bool RemoveComments=true);  //CO20210315 - cleaned up
  bool substring2bool(const stringstream& strstream,const string& strsub1,bool RemoveWS=false,bool RemoveComments=true);  //CO20210315 - cleaned up
  bool substring_present_file(const string& FileName,const string& strsub1,bool RemoveWS=false,bool RemoveComments=true); //CO20210315 - cleaned up
  bool substring_present_file_FAST(const string& FileName,const string& strsub1,bool RemoveWS=false,bool case_insensitive=false,bool expect_near_end=false,unsigned long long int size_max=AUROSTD_MAX_ULLINT);  //CO20210315 - cleaned up
  bool WithinList(const vector<string>& list,const string& input,bool sorted=false);  //CO20181010
  bool WithinList(const deque<string>& list,const string& input,bool sorted=false);  //CO20181010
  bool WithinList(const vector<int>& list,int input,bool sorted=false); //CO20181010
  bool WithinList(const vector<uint>& list,uint input,bool sorted=false); //CO20181010
  bool WithinList(const vector<string>&, const string&, int&,bool sorted=false);  //ME20190905
  bool WithinList(const vector<int>&, int, int&,bool sorted=false);  //ME20190905
  bool WithinList(const vector<uint>&, uint, int&,bool sorted=false);  //ME20190905
  bool EWithinList(const vector<string>& list,const string& input); //CO20200223
  bool EWithinList(const vector<string>& list, const string& input, string& output); //CO20200223
  // about present substrings and taking off the value
  string substring2string(const string& strstream, const string& strsub1, bool RemoveWS=false,bool RemoveComments=true);  //CO20210315 - cleaned up
  string substring2string(const stringstream& strstream, const string& strsub1, bool RemoveWS=false,bool RemoveComments=true);  //CO20210315 - cleaned up
  //[CO20210315 - not used, not sure the purpose of strsub2]string substring2string(const string& strstream, const string& strsub1, const string& strsub2, bool RemoveWS=false,bool RemoveComments=true); //CO20210315 - cleaned up
  template<typename utype> utype substring2utype(const string& strstream,const string& strsub1,bool RemoveWS=false,bool RemoveComments=true); //CO20210315 - cleaned up
  template<typename utype> utype substring2utype(const stringstream& strstream,const string& strsub1,bool RemoveWS=false,bool RemoveComments=true); //CO20210315 - cleaned up
  //[CO20210315 - not used, not sure the purpose of strsub2]template<typename utype> utype substring2utype(const string& strstream, const string& strsub1, const string& strsub2, bool RemoveWS=false,bool RemoveComments=true);  //CO20210315 - cleaned up

  bool kvpairfound(const string& strstream,const string& keyword,const string& delim,bool RemoveWS=false,bool RemoveComments=true);  //CO20210315
  bool kvpairfound(const stringstream& strstream,const string& keyword,const string& delim,bool RemoveWS=false,bool RemoveComments=true);  //CO20210315
  string kvpair2value(const string& strstream,const string& keyword,const string& delim,bool RemoveWS=false,bool RemoveComments=true);  //CO20210315
  string kvpair2value(const stringstream& strstream,const string& keyword,const string& delim,bool RemoveWS=false,bool RemoveComments=true);  //CO20210315
  template<typename utype> utype kvpair2utype(const string& strstream,const string& keyword,const string& delim,bool RemoveWS=false,bool RemoveComments=true); //CO20210315 - cleaned up
  template<typename utype> utype kvpair2utype(const stringstream& strstream,const string& keyword,const string& delim,bool RemoveWS=false,bool RemoveComments=true); //CO20210315 - cleaned up

  uint substring2strings(const string& strstream, vector<string> &vstringout, const string& strsub1, bool RemoveWS=false,bool RemoveComments=true); //CO20210315 - cleaned up
  //[CO20210315 - not used, not sure the purpose of strsub2]uint substring2strings(const string& strstream, vector<string> &vstringout, const string& strsub1, const string& strsub2, bool RemoveWS=false,bool RemoveComments=true);  //CO20210315 - cleaned up
  template<typename utype> uint substring2utypes(const string& strstream, vector<string> &vstringout, const string& strsub1, bool RemoveWS=false,bool RemoveComments=true); //CO20210315 - cleaned up
  template<typename utype> uint substring2utypes(const stringstream& strstream, vector<string> &vstringout, const string& strsub1, bool RemoveWS=false,bool RemoveComments=true); //CO20210315 - cleaned up
  //[CO20210315 - not used, not sure the purpose of strsub2]template<typename utype> uint substring2utypes(const string& strstream, vector<string> &vstringout, const string& strsub1, const string& strsub2, bool RemoveWS=false,bool RemoveComments=true);  //CO20210315 - cleaned up
}

// ***************************************************************************
namespace aurostd { // aurostd_crc64.cpp
  uint64_t crc64(uint64_t crc, const unsigned char *s, uint64_t l);
  uint64_t crc64(uint64_t crc, const string s);
  string crc2string(uint64_t crc);
  int crc64_main(void);
}

// ***************************************************************************

namespace aurostd {
  string text2html(const string& str) __xprototype;  //ME20200921
  string html2latex(const string& str) __xprototype;
  string html2txt(const string& str) __xprototype;
  string string2latex(const string& str) __xprototype;
  string latex2html(const string& str) __xprototype;
  string latex2txt(const string& str) __xprototype;
  string fixStringLatex(const string& input, bool double_back_slash=false,bool symmetry_string=false);  //CO20190419
}


// ***************************************************************************
// SORT WORLD
// double
// ----------------------------------------------------------------------------
// sort for vectors

namespace aurostd {
  template<class utype1> void sort(vector<utype1>& arr);
  template<class utype1> void sort_remove_duplicates(vector<utype1>& arr);
  template<class utype1,class utype2> void sort(vector<utype1>& arr, vector<utype2>& brr);
  template<class utype1,class utype2> void sort(deque<utype1>& arr, deque<utype2>& brr);  //CO20200915
  template<class utype1,class utype2,class utype3> void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr);
  template<class utype1,class utype2,class utype3,class utype4> void sort(vector<utype1>& arr, vector<utype2>& brr, vector<utype3>& crr, vector<utype4>& drr);
}

namespace aurostd {   // DOUBLE
  class _sort_double_value0 {                    // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]<v2[0]); } };
  class _isort_double_value0 {                   // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]>v2[0]); } };
  class _sort_double_value1 {                    // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[1]<v2[1]); } };
  class _isort_double_value1 {                   // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[1]>v2[1]); } };
  class _sort_double_value2 {                    // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[2]<v2[2]); } };
  class _isort_double_value2 {                   // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[2]>v2[2]); } };
  class _sort_double_value3 {                    // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[3]<v2[3]); } };
  class _isort_double_value3 {                   // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[3]>v2[3]); } };
  class _sort_double_value4 {                    // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[4]<v2[4]); } };
  class _isort_double_value4 {                   // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[4]>v2[4]); } };
  class _sort_double_value5 {                    // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[5]<v2[5]); } };
  class _isort_double_value5 {                   // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[5]>v2[5]); } };
  class _sort_double_value6 {                    // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[6]<v2[6]); } };
  class _isort_double_value6 {                   // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[6]>v2[6]); } };
  class _sort_double_value7 {                    // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[7]<v2[7]); } };
  class _isort_double_value7 {                   // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[7]>v2[7]); } };
  class _sort_double_value8 {                    // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[8]<v2[8]); } };
  class _isort_double_value8 {                   // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[8]>v2[8]); } };
  class _sort_double_value9 {                    // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[9]<v2[9]); } };
  class _isort_double_value9 {                   // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[9]>v2[9]); } };
  class _sort_double_value01 {                // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]<v2[0]+v2[1]); } };
  class _isort_double_value01 {               // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]>v2[0]+v2[1]); } };
  class _sort_double_value012 {                // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]<v2[0]+v2[1]+v2[2]); } };
  class _isort_double_value012 {               // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]>v2[0]+v2[1]+v2[2]); } };
  class _sort_double_value0123 {                // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]<v2[0]+v2[1]+v2[2]+v2[3]); } };
  class _isort_double_value0123 {               // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]>v2[0]+v2[1]+v2[2]+v2[3]); } };
  class _sort_double_value01234 {                // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]+v1[4]<v2[0]+v2[1]+v2[2]+v2[3]+v2[4]); } };
  class _isort_double_value01234 {               // sorting through reference
    public:
      bool operator()(const vector<double>& v1, const vector<double>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]+v1[4]>v2[0]+v2[1]+v2[2]+v2[3]+v2[4]); } };
}
// int
namespace aurostd {   // INT
  class _sort_int_value0 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]<v2[0]); } };
  class _isort_int_value0 {                   // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]>v2[0]); } };
  class _sort_int_value1 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[1]<v2[1]); } };
  class _isort_int_value1 {                   // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[1]>v2[1]); } };
  class _sort_int_value2 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[2]<v2[2]); } };
  class _isort_int_value2 {                   // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[2]>v2[2]); } };
  class _sort_int_value3 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[3]<v2[3]); } };
  class _isort_int_value3 {                   // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[3]>v2[3]); } };
  class _sort_int_value4 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[4]<v2[4]); } };
  class _isort_int_value4 {                   // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[4]>v2[4]); } };
  class _sort_int_value5 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[5]<v2[5]); } };
  class _isort_int_value5 {                   // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[5]>v2[5]); } };
  class _sort_int_value6 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[6]<v2[6]); } };
  class _isort_int_value6 {                   // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[6]>v2[6]); } };
  class _sort_int_value7 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[7]<v2[7]); } };
  class _isort_int_value7 {                   // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[7]>v2[7]); } };
  class _sort_int_value8 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[8]<v2[8]); } };
  class _isort_int_value8 {                   // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[8]>v2[8]); } };
  class _sort_int_value9 {                    // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[9]<v2[9]); } };
  class _isort_int_value9 {                   // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[9]>v2[9]); } };
  class _sort_int_value01 {                // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]<v2[0]+v2[1]); } };
  class _isort_int_value01 {               // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]>v2[0]+v2[1]); } };
  class _sort_int_value012 {                // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]<v2[0]+v2[1]+v2[2]); } };
  class _isort_int_value012 {               // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]>v2[0]+v2[1]+v2[2]); } };
  class _sort_int_value0123 {                // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]<v2[0]+v2[1]+v2[2]+v2[3]); } };
  class _isort_int_value0123 {               // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]>v2[0]+v2[1]+v2[2]+v2[3]); } };
  class _sort_int_value01234 {                // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]+v1[4]<v2[0]+v2[1]+v2[2]+v2[3]+v2[4]); } };
  class _isort_int_value01234 {               // sorting through reference
    public:
      bool operator()(const vector<int>& v1, const vector<int>& v2) const { return (bool) (v1[0]+v1[1]+v1[2]+v1[3]+v1[4]>v2[0]+v2[1]+v2[2]+v2[3]+v2[4]); } };
  // STRING
  void sort(vector<string>& arg);
  void sort(deque<string>& arg);
  void sort_remove_duplicates(vector<string>& arg);
  void sort_remove_duplicates(deque<string>& arg);
  class _sort_string_ {                    // sorting through reference
    public:
      bool operator()(const string& str1, const string& str2) const { return (bool) (str1<str2);}
  };
  void rsort(vector<string>& arg);
  void rsort(deque<string>& arg);
  void rsort_remove_duplicates(vector<string>& arg);
  void rsort_remove_duplicates(deque<string>& arg);


  // _STRING_INT_
  void sort(vector<string>& varg1,vector<int>& varg2);
  void sort(deque<string>& varg1,deque<int>& varg2);
  struct _string_int_ {
    string arg1;
    int arg2;
  };
  class _sort_string_int_ {                    // sorting through reference
    public:
      bool operator()(const _string_int_& x1, const _string_int_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_DOUBLE_
  void sort(vector<string>& varg1,vector<double>& varg2);
  void sort(deque<string>& varg1,deque<double>& varg2);
  struct _string_double_ {
    string arg1;
    double arg2;
  };
  class _sort_string_double_ {                    // sorting through reference
    public:
      bool operator()(const _string_double_& x1, const _string_double_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_STRING_
  void sort(vector<string>& varg1,vector<string>& varg2);
  void sort(deque<string>& varg1,deque<string>& varg2);
  struct _string_string_ {
    string arg1;
    string arg2;
  };
  class _sort_string_string_ {                    // sorting through reference
    public:
      bool operator()(const _string_string_& x1, const _string_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _DOUBLE_INT_
  void sort(vector<double>& varg1,vector<int>& varg2);
  void sort(deque<double>& varg1,deque<int>& varg2);
  struct _double_int_ {
    double arg1;
    int arg2;
  };
  class _sort_double_int_ {                    // sorting through reference
    public:
      bool operator()(const _double_int_& x1, const _double_int_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _DOUBLE_DOUBLE_
  void sort(vector<double>& varg1,vector<double>& varg2);
  void sort(deque<double>& varg1,deque<double>& varg2);
  struct _double_double_ {
    double arg1;
    double arg2;
  };
  class _sort_double_double_ {                    // sorting through reference
    public:
      bool operator()(const _double_double_& x1, const _double_double_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _DOUBLE_STRING_
  void sort(vector<double>& varg1,vector<string>& varg2);
  void sort(deque<double>& varg1,deque<string>& varg2);
  struct _double_string_ {
    double arg1;
    string arg2;
  };
  class _sort_double_string_ {                    // sorting through reference
    public:
      bool operator()(const _double_string_& x1, const _double_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_INT_STRING
  void sort(vector<string>& varg1,vector<int>& varg2,vector<string>& varg3);
  void sort(deque<string>& varg1,deque<int>& varg2,deque<string>& varg3);
  struct _string_int_string_ {
    string arg1;
    int arg2;
    string arg3;
  };
  class _sort_string_int_string_ {             // sorting through reference
    public:
      bool operator()(const _string_int_string_& x1, const _string_int_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_DOUBLE_STRING
  void sort(vector<string>& varg1,vector<double>& varg2,vector<string>& varg3);
  void sort(deque<string>& varg1,deque<double>& varg2,deque<string>& varg3);
  struct _string_double_string_ {
    string arg1;
    double arg2;
    string arg3;
  };
  class _sort_string_double_string_ {             // sorting through reference
    public:
      bool operator()(const _string_double_string_& x1, const _string_double_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_STRING_STRING
  void sort(vector<string>& varg1,vector<string>& varg2,vector<string>& varg3);
  void sort(deque<string>& varg1,deque<string>& varg2,deque<string>& varg3);
  struct _string_string_string_ {
    string arg1;
    string arg2;
    string arg3;
  };
  class _sort_string_string_string_ {             // sorting through reference
    public:
      bool operator()(const _string_string_string_& x1, const _string_string_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_STRING_DOUBLE_STRING
  void sort(vector<string>& varg1,vector<string>& varg2,vector<double>& varg3,vector<string>& varg4);
  void sort(deque<string>& varg1,deque<string>& varg2,deque<double>& varg3,deque<string>& varg4);
  struct _string_string_double_string_ {
    string arg1;
    string arg2;
    double arg3;
    string arg4;
  };
  class _sort_string_string_double_string_ {             // sorting through reference
    public:
      bool operator()(const _string_string_double_string_& x1, const _string_string_double_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
  // _STRING_STRING_DOUBLE_DOUBLE_STRING
  void sort(vector<string>& varg1,vector<string>& varg2,vector<double>& varg3,vector<double>& varg4,vector<string>& varg5);
  void sort(deque<string>& varg1,deque<string>& varg2,deque<double>& varg3,deque<double>& varg4,deque<string>& varg5);
  struct _string_string_double_double_string_ {
    string arg1;
    string arg2;
    double arg3;
    double arg4;
    string arg5;
  };
  class _sort_string_string_double_double_string_ {             // sorting through reference
    public:
      bool operator()(const _string_string_double_double_string_& x1, const _string_string_double_double_string_& x2) const { return (bool) (x1.arg1<x2.arg1);}
  };
}

// ***************************************************************************
// some statistical stuff
template<class utype> utype combinations(utype n,utype k) __xprototype; // http://en.wikipedia.org/wiki/Combination
template<class utype> utype Cnk(utype n,utype k) __xprototype; // http://en.wikipedia.org/wiki/Combination

// ***************************************************************************
namespace aurostd {
  //template <typename utype> utype sum(vector<utype>& a) {
  //utype result = 0;
  //for (unsigned int i=0; i<a.size();i++) result += a.at(i);
  //return result;
  //}
  vector<vector<double> > ShiftFirstColumn(const vector<vector<double> >& vva, const double& value);
  vector<vector<double> > ShrinkValuesExceptFirstColumn(const vector<vector<double> >& vva, const double& Fi);
  vector<vector<double> > NormalizeAndSum3DVector(const vector<vector<vector<double> > >& vvva, const vector<double>& vFi);
  vector<vector<double> > Sum3DVectorAndReduce2D(const vector<vector<vector<double> > >& vvva);
  vector<vector<double> > Sum2DVectorExceptFirstColumn(const vector<vector<double> >& vva, const vector<vector<double> >& vvb);
  string vector2string(const vector<vector<double> >& vva);
  template<typename utype> deque<utype> vector2deque(const vector<utype>& vin); //CO20181226
  template<typename utype> vector<utype> deque2vector(const deque<utype>& din); //CO20181226
  vector<vector<double> > ReduceVector(const vector<vector<double> >& vva, const int& n);
  double CalculateIntegrate(const vector<vector<double> >& vva, const int& n);
  double CalculateIntegrate(const vector<vector<double> >& vva, const int& n, const double& Emin, const double& Emax);
  double CalculateIntegrate(const vector<vector<double> >& vva);
  double CalculateIntegrate(const vector<vector<double> >& vva, const double& Emin, const double& Emax);
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva);
  double FindMaxIn2DvectorExcept1stColumn(const vector<vector<double> >& vva, const double& min, const double& max);
  double FindMaxInTDOS(const vector<vector<double> >& vva, const double& min, const double& max);
}

// ***************************************************************************


// ***************************************************************************
bool initialize_templates_never_call_this_procedure(bool flag);

////////////////////////////////////////////////////////////////////////////////
namespace aurostd {
  // joins int/string type of objects together by a delimiter
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const stringstream& delimiter);
  string joinWDelimiter(const xvector<int>& ientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter);
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter);
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter);
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<int>& ientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<uint>& uientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const stringstream& delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const vector<string>& _sentries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter);
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter);
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter);
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const deque<int>& ientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const deque<uint>& uientries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter,
      const char& _l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const char& _delimiter,
      const char& _m_delimiter, const char& _l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter,
      const string& _l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const string& _delimiter,
      const string& _m_delimiter, const string& _l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter,
      const stringstream& l_delimiter);
  string joinWDelimiter(const deque<string>& _sentries, const stringstream& delimiter,
      const stringstream& m_delimiter,
      const stringstream& l_delimiter);
}
////////////////////////////////////////////////////////////////////////////////

//DX20180118 - Add xcomplex to json
namespace aurostd {
  template<typename utype> string _xcomplex2json(xcomplex<utype>& number) __xprototype;
  string xcomplex2json(xcomplex<double>& number);
}

//DX20170803 - Add Matrix print out
namespace aurostd {
  // [OBSOLETE] string xmatDouble2String(const xmatrix<double>& xmat_in, bool roff=false);
  string xmatDouble2String(const xmatrix<double>& xmat_in, int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM);
}

//CO20171215 - more json functionality
namespace aurostd {
  string wrapString(const string& input,const string& wrapper);
  string wrapString(const string& input,const string& wrapper_start,const string& wrapper_end);
}

namespace aurostd {
  // [OBSOLETE] vector<string> vecDouble2vecString(const vector<double>& vin, bool roff=false);
  vector<string> vecDouble2vecString(const vector<double>& vin,int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM);
  // [OBSOLETE] vector<string> xvecDouble2vecString(const xvector<double>& vin, bool roff=false);
  vector<string> xvecDouble2vecString(const xvector<double>& vin,int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM);
  // [OBSOLETE] deque<string> deqDouble2deqString(const deque<double>& vin, bool roff=false);
  deque<string> vecDouble2vecString(const deque<double>& vin,int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM);
  // deque<string> deqDouble2deqString(const deque<double>& vin,int precision=AUROSTD_DEFAULT_PRECISION, bool roff=false, double tol=AUROSTD_ROUNDOFF_TOL, char FORMAT=DEFAULT_STREAM);
}

namespace aurostd {
  vector<string> wrapVecEntries(const vector<string>& vin,string wrap);
  vector<string> wrapVecEntries(const vector<string>& vin,string wrap_start,string wrap_end);
  deque<string> wrapVecEntries(const deque<string>& vin,string wrap);                          //SC20200329 nice overload to deal with ME
  deque<string> wrapVecEntries(const deque<string>& vin,string wrap_start,string wrap_end);    //SC20200329 nice overload to deal with ME
}

//base64 stuff
//CO START
namespace aurostd {
  //http://www.adp-gmbh.ch/cpp/common/base64.html for base64 encoding/decoding
  //http://stackoverflow.com/questions/535444/custom-manipulator-for-c-iostream for fancy stream manipulators
  //static inline bool isBase64(unsigned char c); //determines if char is base64
  inline bool isBase64(unsigned char c); //determines if char is base64
  std::string base64Encoder(unsigned char const* bytes_to_encode, unsigned int in_len); //encodes bytes to base64
  std::string base64Decoder(std::string const& encoded_string); //decodes base64 to bytes
  bool bin2base64(const std::string& b_file, std::string& b64String); //converts binary file to base64 string
  bool base642bin(const std::string& b64String, const std::string& b_file); //converts base64 string to binary file

  //structure that allows "cout << b64_encoder << ifstream"
  struct b64_encoder_proxy {
    explicit b64_encoder_proxy(std::ostream & os):os(os){}

    template<typename Rhs>
      friend std::ostream & operator<<(b64_encoder_proxy const& b, 
          Rhs const& rhs) {
        return b.os << rhs;
      }

    friend std::ostream & operator<<(b64_encoder_proxy const& b, 
        ifstream& file) {
      // Stop eating new lines in binary mode!!!
      file.unsetf(std::ios::skipws);

      // get its size:
      std::streampos fileSize;

      file.seekg(0, std::ios::end);
      fileSize = file.tellg();
      file.seekg(0, std::ios::beg);

      // read the data:
      vector<unsigned char> vec((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));

      return b.os << base64Encoder(reinterpret_cast<const unsigned char*>(&vec[0]), vec.size());
    }

    private:
    std::ostream & os;
  };

  //structure that allows "ofstream << b64_decoder << string" 
  //WARNING, SPITS BINARY CHARACTERS OUT IF COUT
  struct b64_decoder_proxy {
    explicit b64_decoder_proxy(std::ostream & os):os(os){}

    template<typename Rhs>
      friend std::ostream & operator<<(b64_decoder_proxy const& b, 
          Rhs const& rhs) {
        return b.os << rhs;
      }

    friend std::ostream & operator<<(b64_decoder_proxy const& b, 
        std::string const& rhs) {
      return b.os << base64Decoder(rhs);
    }

    private:
    std::ostream & os;
  };

  //structure that allows "cout << b64_encoder << ifstream"
  struct b64_encoder_creator { };
  extern b64_encoder_creator b64_encoder;
  b64_encoder_proxy operator<<(std::ostream & os, b64_encoder_creator);

  //structure that allows "ofstream << b64_decoder << string" 
  //WARNING, SPITS BINARY CHARACTERS OUT IF COUT
  struct b64_decoder_creator { };
  extern b64_decoder_creator b64_decoder;
  b64_decoder_proxy operator<<(std::ostream & os, b64_decoder_creator);
}

//binary to base64 conversion
const std::string base64_chars = 
"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
"abcdefghijklmnopqrstuvwxyz"
"0123456789+/";
//CO END

#endif

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************
