// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XSCALAR_CPP_
#define _AUROSTD_XSCALAR_CPP_

#ifndef XXEND
#define XXEND 1
#endif

#ifndef _AUROSTD_XSCALAR_H_
#include "aurostd_xscalar.h"
#endif

// ----------------------------------------------------------------------------
//namespace aurostd {
//template<class utype> utype                        
//nint(utype x) {
//if(x>=0) return (int)(x+0.5);
//else return (int)(x-0.5);
//}
//}

// namespace aurostd {
//   double atof(string str) { return atof(str.c_str());}
//   int atoi(string str) { return atoi(str.c_str());}
//   long atol(string str) { return atol(str.c_str());}
//   // long long atoll(string str) { return atoll(str.c_str());}
//   // long long atoq(string str) { return atoq(str.c_str());}
// }

// ----------------------------------------------------------------------------
// _isfloat  _isfloat  _isfloat  _isfloat  _isfloat
namespace aurostd {
  // namespace aurostd
  bool _isfloat(bool)          { return (bool) FALSE; }
  bool _isfloat(char)          { return (bool) FALSE; }
  bool _isfloat(uint)          { return (bool) FALSE; }
  bool _isfloat(int)           { return (bool) FALSE; }
  bool _isfloat(long int)      { return (bool) FALSE; }
  bool _isfloat(unsigned long int)      { return (bool) FALSE; }  //CO20191201
  bool _isfloat(long long int) { return (bool) FALSE; }
  bool _isfloat(unsigned long long int) { return (bool) FALSE; }
  bool _isfloat(float)         { return (bool) TRUE;  }
  bool _isfloat(double)        { return (bool) TRUE;  }
  bool _isfloat(long double)   { return (bool) TRUE;  }
#ifdef _AUROSTD_XCOMPLEX_
  bool _isfloat(xcomplex<float>) { return (bool) TRUE;}
  bool _isfloat(xcomplex<double>) { return (bool) TRUE;}
  bool _isfloat(xcomplex<long double>) { return (bool) TRUE;}
#endif
  bool isfloat(const string& in){ //CO20180729
    stringstream ss;
    double num;
    ss << in;
    return ((bool)(ss >> num));
  }
}

// ----------------------------------------------------------------------------
// _iscomplex  _iscomplex  _iscomplex  _iscomplex  _iscomplex
namespace aurostd {
  bool _iscomplex(bool)        { return (bool) FALSE; }
  bool _iscomplex(char)        { return (bool) FALSE; }
  bool _iscomplex(uint)        { return (bool) FALSE; }
  bool _iscomplex(int)         { return (bool) FALSE; }
  bool _iscomplex(long int)    { return (bool) FALSE; }
  bool _iscomplex(unsigned long int)    { return (bool) FALSE; }  //CO20191201
  bool _iscomplex(long long int) { return (bool) FALSE; }
  bool _iscomplex(unsigned long long int) { return (bool) FALSE; }
  bool _iscomplex(float)       { return (bool) FALSE; }
  bool _iscomplex(double)      { return (bool) FALSE; }
  bool _iscomplex(long double) { return (bool) FALSE; }
#ifdef _AUROSTD_XCOMPLEX_
  bool _iscomplex(xcomplex<float>) { return (bool) TRUE;}
  bool _iscomplex(xcomplex<double>) { return (bool) TRUE;}
  bool _iscomplex(xcomplex<long double>) { return (bool) TRUE;}
#endif
}

// ----------------------------------------------------------------------------
// _isreal  _isreal  _isreal  _isreal  _isreal
namespace aurostd {
  bool _isreal(bool)           { return (bool) TRUE; }
  bool _isreal(char)           { return (bool) TRUE; }
  bool _isreal(uint)           { return (bool) TRUE; }
  bool _isreal(int)            { return (bool) TRUE; }
  bool _isreal(long int)       { return (bool) TRUE; }
  bool _isreal(unsigned long int)       { return (bool) TRUE; } //CO20191201
  bool _isreal(long long int)  { return (bool) TRUE; }
  bool _isreal(unsigned long long int)  { return (bool) TRUE; }
  bool _isreal(float)          { return (bool) TRUE; }
  bool _isreal(double)         { return (bool) TRUE; }
  bool _isreal(long double)    { return (bool) TRUE; }
#ifdef _AUROSTD_XCOMPLEX_
  bool _isreal(xcomplex<float>) { return (bool) FALSE;}
  bool _isreal(xcomplex<double>) { return (bool) FALSE;}
  bool _isreal(xcomplex<long double>) { return (bool) FALSE;}
#endif
}

// ----------------------------------------------------------------------------
// _size  _size  _size  _size  _size
namespace aurostd {
  int _size(bool)              { return (int) sizeof(bool);}
  int _size(char)              { return (int) sizeof(char);}
  int _size(uint)              { return (int) sizeof(uint);}
  int _size(int)               { return (int) sizeof(int);}
  int _size(long int)          { return (int) sizeof(long int);}
  int _size(unsigned long int)          { return (int) sizeof(unsigned long int);} //CO20191201
  int _size(long long int)     { return (int) sizeof(long long int);}
  int _size(unsigned long long int)     { return (int) sizeof(unsigned long long int);}
  int _size(float)             { return (int) sizeof(float);}
  int _size(double)            { return (int) sizeof(double);}
  int _size(long double)       { return (int) sizeof(long double);}
#ifdef _AUROSTD_XCOMPLEX_
  int _size(xcomplex<float>)   { return (int) 2*sizeof(float);}
  int _size(xcomplex<double>)  { return (int) 2*sizeof(double);}
  int _size(xcomplex<long double>) { return (int) 2*sizeof(long double);}
#endif
}

// ----------------------------------------------------------------------------
// _real  _real  _real  _real  _real
namespace aurostd {
  bool _real(bool x) { return (bool) x;}
  char _real(char x) { return (char) x;}
  uint _real(uint x) { return (uint) x;}
  int _real(int x) { return (int) x;}
  long int _real(long int x) { return (long int) x;}
  unsigned long int _real(unsigned long int x) { return (unsigned long int) x;} //CO20191201
  long long int _real(long long int x) { return (long long int) x;}
  unsigned long long int _real(unsigned long long int x) { return (unsigned long long int) x;}
  float _real(float x) { return (float) x;}
  double _real(double x) { return (double) x;}
  long double _real(long double x) { return (long double) x;}
#ifdef _AUROSTD_XCOMPLEX_
  int _real(xcomplex<int> x) { return (int) real(x);}
  float _real(xcomplex<float> x) { return (float) real(x);}
  double _real(xcomplex<double> x) { return (double) real(x);}
  long double _real(xcomplex<long double> x) { return (long double) real(x);}
#endif
}



// ----------------------------------------------------------------------------
// _real _real _real _real _real _real _real _real
namespace aurostd {
  // namespace aurostd
  bool _real(bool) __xprototype;
  char _real(char) __xprototype;
  int _real(int) __xprototype;
  uint _real(uint) __xprototype;
  float _real(float) __xprototype;
  double _real(double) __xprototype;
  long int _real(long int) __xprototype;
  long long int _real(long long int) __xprototype;
  unsigned long long int _real(unsigned long long int) __xprototype;
  long double _real(long double) __xprototype;
#ifdef _AUROSTD_XCOMPLEX_
  float _real(xcomplex<float>)  __xprototype;
  double _real(xcomplex<double>) __xprototype;
  long double _real(xcomplex<long double>) __xprototype;
#endif
}

// ----------------------------------------------------------------------------
// abs  abs  abs  abs  abs
namespace aurostd {
  // namespace aurostd
  // ABS(X)
  // int abs(int x) { return (int) (x<0? -x:x);}
  char abs(char x) { return (char) std::abs(x);}
  int abs(int x) { return (int) std::abs(x);}
  uint abs(uint x) { return (uint) x;}
  float abs(float x) { return (float) fabsf(x);}
  double abs(double x) { return (double) std::fabs(x);}
  long int abs(long int x) { return (long int) std::labs(x);}
  long long int abs(long long int x) { return (long long int) std::llabs(x);}
  unsigned long int abs(unsigned long int x) { return (unsigned long int) x;} //std::llabs(x);  //CO, unsigned is already abs
  unsigned long long int abs(unsigned long long int x) { return (unsigned long long int) x;} //std::llabs(x);  //CO, unsigned is already abs
  long double abs(long double x) { return (long double) fabsl(x);}
#ifdef _AUROSTD_XCOMPLEX_
  float abs(xcomplex<float> x) { return (float) sqrtf(x.re*x.re+x.im*x.im);}
  double abs(xcomplex<double> x) { return (double) std::sqrt(x.re*x.re+x.im*x.im);}
  long double abs(xcomplex<long double> x) { return (long double) sqrtl(x.re*x.re+x.im*x.im);}
#endif
}

// ----------------------------------------------------------------------------
// sqrt  sqrt  sqrt  sqrt  sqrt
namespace aurostd {
  // namespace aurostd
  // SQRT(X)
  // int sqrt(int x) { return (int) (x<0? -x:x);}
  char sqrt(char x) { return (char) std::sqrt((double) x);}
  int sqrt(int x) { return (int) std::sqrt((double) x);}
  uint sqrt(uint x) { return (uint) std::sqrt((double) x);}
  float sqrt(float x) { return (float) sqrtf(x);}
  double sqrt(double x) { return (double) std::sqrt(x);}
  long int sqrt(long int x) { return (long int) std::sqrt((double) x);}
  long long int sqrt(long long int x) { return (long long int) std::sqrt((double) x);}
  unsigned long long int sqrt(unsigned long long int x) { return (unsigned long long int) std::sqrt((double) x);}
  long double sqrt(long double x) { return (long double) sqrtl(x);}
#ifdef _AUROSTD_XCOMPLEX_
  float sqrt(xcomplex<float> x) { return (float) sqrtf(x.re*x.re+x.im*x.im);}  // not defined Dec09
  double sqrt(xcomplex<double> x) { return (double) std::sqrt(x.re*x.re+x.im*x.im);}  // not defined Dec09
  long double sqrt(xcomplex<long double> x) { return (long double) sqrtl(x.re*x.re+x.im*x.im);}  // not defined Dec09
#endif
}

// ----------------------------------------------------------------------------
// round  floor ceil trunc
// [OBSOLETE] namespace aurostd {  // namespace aurostd
// [OBSOLETE]   // ROUND(X)
// [OBSOLETE]  double round(double x) { return (double) std::round(double(x));}
// [OBSOLETE]  float round(float x) { return (float) std::roundf(float(x));}
// [OBSOLETE]  long double round(long double x) { return (long double) std::roundl((long double) x);}
// [OBSOLETE]  int round(int x) { return (int) std::round(double(x));}
// [OBSOLETE]  long round(long x) { return (long) std::round(double(x));}
// [OBSOLETE]  // FLOOR(X)
// [OBSOLETE]  //  double floor(double x) { return (double) std::floor(double(x));}
// [OBSOLETE]  float floor(float x) { return (float) std::floorf(float(x));}
// [OBSOLETE]  long double floor(long double x) { return (long double) std::floorl((long double) x);}
// [OBSOLETE]  int floor(int x) { return (int) std::floor(double(x));}
// [OBSOLETE]  long floor(long x) { return (long) std::floor(double(x));}
// [OBSOLETE]  // CEIL(X)
// [OBSOLETE]  double ceil(double x) { return (double) std::ceil(double(x));}
// [OBSOLETE]  float ceil(float x) { return (float) std::ceilf(float(x));}
// [OBSOLETE]  long double ceil(long double x) { return (long double) std::ceill((long double) x);}
// [OBSOLETE]  int ceil(int x) { return (int) std::ceil(double(x));}
// [OBSOLETE]  long ceil(long x) { return (long) std::ceil(double(x));}
// [OBSOLETE]  // TRUNC(X)
// [OBSOLETE]  double trunc(double x) { return (double) std::trunc(double(x));}
// [OBSOLETE]  float trunc(float x) { return (float) std::truncf(float(x));}
// [OBSOLETE]  long double trunc(long double x) { return (long double) std::truncl((long double) x);}
// [OBSOLETE]  int trunc(int x) { return (int) std::trunc(double(x));}
// [OBSOLETE]  long trunc(long x) { return (long) std::trunc(double(x));}
// [OBSOLETE] }

namespace aurostd {
  double ln(double x) { return (double) std::log(x);};
  // float lnf(float x) { return (float) std::logf(x);};
  // long double lnl(long double x) { return (long double) std::logl(x);};
  // float ln(float x) { return (float) std::logf(x);};
  // long double ln(long double x) { return (long double) sstd::logl(x);};
  double log(double x) { return(double) std::log(x);};
  // float logf(float x) { return (float) std::logf(x);};
  // float log(float x) { return (float) std::logf(x);};
  //  long double logl(long double x) { return (long double) std::logl(x);};
  // long double log(long double x) { return (long double) std::logl(x);};
  double log10(double x) { return (double) std::log10(x);};
  // float log10f(float x) { return (float) std::log10f(x);};
  // float log10(float x) { return (float) std::log10f(x);};
  // long double log10l(long double x) { return (long double) std::log10l(x);};
  // long double log10(long double x) { return (long double) std::log10l(x);};
}

// ***************************************************************************
// Function sign
// ***************************************************************************
namespace aurostd {
  char sign(char x) {
    if(x>0) return (char) 1;
    if(x<0) return (char) -1;
    return (char) 0;
  }
  int sign(int x) {
    if(x>0) return (int) 1;
    if(x<0) return (int) -1;
    return (int) 0;
  }
  float sign(float x) {
    if(x>0) return (float) 1;
    if(x<0) return (float) -1;
    return (float) 0;
  }
  double sign(double x) {
    if(x>0) return (double) 1;
    if(x<0) return (double) -1;
    return (double) 0;
  }
  long int sign(long int x) {
    if(x>0) return (long int) 1;
    if(x<0) return (long int) -1;
    return (long int) 0;
  }
  long long int sign(long long int x) {
    if(x>0) return (long long int) 1;
    if(x<0) return (long long int) -1;
    return (long long int) 0;
  }
  long double sign(long double x) {
    if(x>0) return (long double) 1;
    if(x<0) return (long double) -1;
    return (long double) 0;
  }
}

// ***************************************************************************
// Function signnozero
// ***************************************************************************
namespace aurostd {
  char signnozero(char x) {
    if(x>=0) return (char) 1;
    return (char) -1;
  }
  int signnozero(int x) {
    if(x>=0) return (int) 1;
    return (int) -1;
  }
  float signnozero(float x) {
    if(x>=0) return (float) 1;
    return (float) -1;
  }
  double signnozero(double x) {
    if(x>=0) return (double) 1;
    return (double) -1;
  }
  long int signnozero(long int x) {
    if(x>=0) return (long int) 1;
    return (long int) -1;
  }
  long long int signnozero(long long int x) {
    if(x>=0) return (long long int) 1;
    return (long long int) -1;
  }
  long double signnozero(long double x) {
    if(x>=0) return (long double) 1;
    return (long double) -1;
  }
}

// ***************************************************************************
// Function nint
// ***************************************************************************
namespace aurostd {
  //int nint(double x) {        // AFLOW_FUNCTION_IMPLEMENTATION
  //  if(x>=0) return (int)(x+0.5);
  // else return (int)(x-0.5);
  //}
  // namespace aurostd
  template<class utype>
    utype nint(utype x) {
      if(x>=0) return (utype) std::floor((double)   x+0.5);
      else      return (utype) -std::floor((double) -x+0.5);
    }
  char nint(char x) {
    return (char) x;
  }
  int nint(int x) {
    if(x>=0) return (int) std::floor((double)   x+0.5);
    else      return (int) -std::floor((double) -x+0.5);
  }
  uint nint(uint x) {
    return (int) std::floor((double)   x+0.5);
    //if(x>=0) return (int) std::floor((double)   x+0.5);
    //else      return (int) -std::floor((double) -x+0.5);
  }
  float nint(float x) {
    if(x>=0) return (float) std::floor((double)   x+0.5);
    else      return (float) -std::floor((double) -x+0.5);
  }
  double nint(double x) {
    if(x>=0) return (double) std::floor((double)   x+0.5);
    else      return (double) -std::floor((double) -x+0.5);
  }
  long int nint(long int x) {
    if(x>=0) return (long int) std::floor((double)   x+0.5);
    else      return (long int) -std::floor((double) -x+0.5);
  }
  long long int nint(long long int x) {
    if(x>=0) return (long long int) std::floor((double)   x+0.5);
    else      return (long long int) -std::floor((double) -x+0.5);
  }
  long long unsigned int nint(long long unsigned int x) {
    return (long long unsigned int) std::floor((double)   x+0.5);
    //if(x>=0) return (long long int) std::floor((double)   x+0.5);//CO20180719
    //else      return (long long int) -std::floor((double) -x+0.5);//CO20180719
  }
  long double nint(long double x) {
    if(x>=0) return (long double) std::floor((double)   x+0.5);
    else      return (long double) -std::floor((double) -x+0.5);
  }
}

// ***************************************************************************
// Function GCD
// ***************************************************************************
namespace aurostd {
  //CO20191112 - extended GCD, get Bezout coefficients
  //algorithm inspired by python solution of https://brilliant.org/wiki/extended-euclidean-algorithm/
  //ax+by=gcd(a,b)
  //a=Z,b=0;gcd=Z;x=1;y=0;
  //a=0,b=Z;gcd=Z;x=0;y=1;
  //a=0,b=0;gcd=undefined;x=undefined;y=undefined  //all integers are common divisors of 0 and 0, so there is no greatest one.
  //note this implementation gives different results from matlab for gcd(1,1): either x or y can be 1, the other is zero
  template<class utype>
    void _GCD(utype a,utype b,utype& gcd,utype& x,utype& y){ //CO20180409
      if(!a && !b){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","gcd(0,0) is undefined",_INPUT_ILLEGAL_);} //only special case needed, all other cases work perfectly
      if(false && isequal(a,(utype)1) && isequal(b,(utype)1)){gcd=1;x=0;y=1;return;} //matlab implementation, this algorithm gives x=1,y=0 which IS valid
      utype a_orig=a,b_orig=b;
      x=(utype)0;y=(utype)1;
      utype u=(utype)1,v=(utype)0;
      utype q=(utype)0,r=(utype)0,m=(utype)0,n=(utype)0;
      while(a){ //is not 0
        q=(utype)std::floor(b/a); //beware of negative numbers, not just /
        r=b%a;
        m=x-u*q;
        n=y-v*q;
        b=a;a=r;
        x=u;y=v;
        u=m;v=n;
      }
      gcd=b;
      if(std::signbit(gcd)){gcd=-gcd;x=-x;y=-y;}  //GCD(a,-b)==GCD(a,b), so flip all the signs
      if(!isequal(a_orig*x+b_orig*y,gcd)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","Bezout's identity not satisfied",_RUNTIME_ERROR_);}
    }
  template<class utype>
    void _GCD(utype a,utype b,utype& gcd){ //CO20180409  //keep this one too, fewer operations than if you need x and y too
      // added for safety, will always give nonzero result, important for division!
      if(a==(utype)0 && b==(utype)0) {throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","gcd(0,0) is undefined",_INPUT_ILLEGAL_);}  //special case
      else if(a==(utype)0) {gcd=b;return;} //special case
      else if(b==(utype)0) {gcd=a;return;} //special case
      // borrowed from KY aflow_contrib_kesong_pocc_basic.cpp
      // calculate greatest common denominator of two integers
      if(a % b == (utype)0) {gcd=b;return;}
      else {return GCD(b, a % b, gcd);}
    }
  void GCD(int a,int b,int& gcd,int& x,int& y){return _GCD(a,b,gcd,x,y);} //CO20191201
  void GCD(int a,int b,int& gcd){return _GCD(a,b,gcd);} //CO20191201
  void GCD(uint a,uint b,uint& gcd,uint& x,uint& y){return _GCD(a,b,gcd,x,y);}  //CO20191201
  void GCD(uint a,uint b,uint& gcd){return _GCD(a,b,gcd);}  //CO20191201
  void GCD(long int a,long int b,long int& gcd,long int& x,long int& y){return _GCD(a,b,gcd,x,y);}  //CO20191201
  void GCD(long int a,long int b,long int& gcd){return _GCD(a,b,gcd);}  //CO20191201
  void GCD(unsigned long int a,unsigned long int b,unsigned long int& gcd,unsigned long int& x,unsigned long int& y){return _GCD(a,b,gcd,x,y);} //CO20191201
  void GCD(unsigned long int a,unsigned long int b,unsigned long int& gcd){return _GCD(a,b,gcd);} //CO20191201
  void GCD(long long int a,long long int b,long long int& gcd,long long int& x,long long int& y){return _GCD(a,b,gcd,x,y);} //CO20191201
  void GCD(long long int a,long long int b,long long int& gcd){return _GCD(a,b,gcd);} //CO20191201
  void GCD(unsigned long long int a,unsigned long long int b,unsigned long long int& gcd,unsigned long long int& x,unsigned long long int& y){return _GCD(a,b,gcd,x,y);}  //CO20191201
  void GCD(unsigned long long int a,unsigned long long int b,unsigned long long int& gcd){return _GCD(a,b,gcd);}  //CO20191201

  void GCD(float a,float b,float& gcd,float& x,float& y,float tolerance){  //CO20191201
    if(!isinteger(a,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","a[="+aurostd::utype2string(a)+"] is not an integer",_INPUT_ILLEGAL_);}
    if(!isinteger(b,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","b[="+aurostd::utype2string(b)+"] is not an integer",_INPUT_ILLEGAL_);}
    long long int igcd=0,ix=0,iy=0; //CO20191201 - long long int as SNF matrices can get big
    GCD((long long int)aurostd::nint(a),(long long int)aurostd::nint(b),igcd,ix,iy);  //CO20191201 - long long int as SNF matrices can get big
    gcd=(float)igcd;
    x=(float)ix;
    y=(float)iy;
  }
  void GCD(float a,float b,float& gcd,float tolerance){  //CO20191201
    if(!isinteger(a,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","a[="+aurostd::utype2string(a)+"] is not an integer",_INPUT_ILLEGAL_);}
    if(!isinteger(b,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","b[="+aurostd::utype2string(b)+"] is not an integer",_INPUT_ILLEGAL_);}
    long long int igcd=0; //CO20191201 - long long int as SNF matrices can get big
    GCD((long long int)aurostd::nint(a),(long long int)aurostd::nint(b),igcd);  //CO20191201 - long long int as SNF matrices can get big
    gcd=(float)igcd;
  }
  void GCD(double a,double b,double& gcd,double& x,double& y,double tolerance){  //CO20191201
    if(!isinteger(a,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","a[="+aurostd::utype2string(a)+"] is not an integer",_INPUT_ILLEGAL_);}
    if(!isinteger(b,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","b[="+aurostd::utype2string(b)+"] is not an integer",_INPUT_ILLEGAL_);}
    long long int igcd=0,ix=0,iy=0; //CO20191201 - long long int as SNF matrices can get big
    GCD((long long int)aurostd::nint(a),(long long int)aurostd::nint(b),igcd,ix,iy);  //CO20191201 - long long int as SNF matrices can get big
    gcd=(double)igcd;
    x=(double)ix;
    y=(double)iy;
  }
  void GCD(double a,double b,double& gcd,double tolerance){  //CO20191201
    if(!isinteger(a,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","a[="+aurostd::utype2string(a)+"] is not an integer",_INPUT_ILLEGAL_);}
    if(!isinteger(b,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","b[="+aurostd::utype2string(b)+"] is not an integer",_INPUT_ILLEGAL_);}
    long long int igcd=0; //CO20191201 - long long int as SNF matrices can get big
    GCD((long long int)aurostd::nint(a),(long long int)aurostd::nint(b),igcd);  //CO20191201 - long long int as SNF matrices can get big
    gcd=(double)igcd;
  }
  void GCD(long double a,long double b,long double& gcd,long double& x,long double& y,long double tolerance){  //CO20191201
    if(!isinteger(a,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","a[="+aurostd::utype2string(a)+"] is not an integer",_INPUT_ILLEGAL_);}
    if(!isinteger(b,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","b[="+aurostd::utype2string(b)+"] is not an integer",_INPUT_ILLEGAL_);}
    long long int igcd=0,ix=0,iy=0; //CO20191201 - long long int as SNF matrices can get big
    GCD((long long int)aurostd::nint(a),(long long int)aurostd::nint(b),igcd,ix,iy);  //CO20191201 - long long int as SNF matrices can get big
    gcd=(long double)igcd;
    x=(long double)ix;
    y=(long double)iy;
  }
  void GCD(long double a,long double b,long double& gcd,long double tolerance){  //CO20191201
    if(!isinteger(a,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","a[="+aurostd::utype2string(a)+"] is not an integer",_INPUT_ILLEGAL_);}
    if(!isinteger(b,tolerance)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::GCD():","b[="+aurostd::utype2string(b)+"] is not an integer",_INPUT_ILLEGAL_);}
    long long int igcd=0;
    GCD((long long int)aurostd::nint(a),(long long int)aurostd::nint(b),igcd);
    gcd=(long double)igcd;
  }

  //https://www.geeksforgeeks.org/program-to-find-lcm-of-two-numbers/
  int LCM(int a,int b){ //CO20190520
    if(!a || !b){return 0;} //special, trivial case: lcm must really be positive
    int gcd=0;  //CO20191201
    GCD(a,b,gcd); //CO20191201
    return a*b/gcd; //CO20191201
  } //CO20190520
} // namespace aurostd

// ***************************************************************************
// Function isinteger
// ***************************************************************************
// namespace aurostd {
//   // namespace aurostd
//   template<class utype>
//   bool isinteger(utype x) {
//     if(aurostd::abs(x-((utype)aurostd::nint(x)))<0.01) return TRUE;
//     return FALSE;
//   }

//   // bool _aurostd_initialize_isinteger(bool x) { return isinteger(x);}
//   // bool _aurostd_initialize_isinteger(char x) { return isinteger(x);}
//   bool _aurostd_initialize_isinteger(int x) { return isinteger(x);}
//   // bool _aurostd_initialize_isinteger(long x) { return isinteger(x);}
//   // bool _aurostd_initialize_isinteger(uint x) { return isinteger(x);}
//   bool _aurostd_initialize_isinteger(float x) { return isinteger(x);}
//   bool _aurostd_initialize_isinteger(double x) { return isinteger(x);}
//   bool _aurostd_initialize_isinteger(long int x) { return isinteger(x);}
//   bool _aurostd_initialize_isinteger(long long int x) { return isinteger(x);}
//   bool _aurostd_initialize_isinteger(unsigned long long int x) { return isinteger(x);}
//   bool _aurostd_initialize_isinteger(long double x) { return isinteger(x);}
// }

namespace aurostd {
  // namespace aurostd
  template<class utype>
    bool _isinteger(utype x,utype tolerance) {
      //DX20191125 [OBSOLETE] if(aurostd::abs(x-((utype)aurostd::nint(x)))<tolerance) return TRUE;
      if(aurostd::abs(x-((utype)aurostd::nint(x)))<=tolerance) return TRUE; //DX20191125 - added <= to account for int and uint and tolerance=0, e.g., 3-nint(3)=0 is not less than 0
      return FALSE;
    }

  //bool isinteger(bool x,bool tolerance){return _isinteger(x,tolerance);}
  //bool isinteger(char x,char tolerance){return _isinteger(x,tolerance);}
  bool isinteger(uint,uint){return true;}  //CO20191201 - obvious but define for boot
  bool isinteger(int,int){return true;}  //CO20191201 - obvious but define for boot
  bool isinteger(long int,long int){return true;}  //CO20191201 - obvious but define for boot
  bool isinteger(unsigned long int,long int){return true;}  //CO20191201 - obvious but define for boot
  bool isinteger(long long int,long long int){return true;}  //CO20191201 - obvious but define for boot
  bool isinteger(unsigned long long int,unsigned long long int){return true;}  //CO20191201 - obvious but define for boot
  bool isinteger(float x,float tolerance){return _isinteger(x,tolerance);}  //CO20191201
  bool isinteger(double x,double tolerance){return _isinteger(x,tolerance);}  //CO20191201
  bool isinteger(long double x,long double tolerance){return _isinteger(x,tolerance);}  //CO20191201

  // bool _aurostd_initialize_isinteger(bool x) { return _isinteger(x);}
  // bool _aurostd_initialize_isinteger(char x) { return _isinteger(x);}
  //bool _aurostd_initialize_isinteger(int x) { return _isinteger(x);}
  // bool _aurostd_initialize_isinteger(long x) { return _isinteger(x);}
  // bool _aurostd_initialize_isinteger(uint x) { return _isinteger(x);}
  //bool _aurostd_initialize_isinteger(float x) { return _isinteger(x);}
  //bool _aurostd_initialize_isinteger(double x) { return _isinteger(x);}
  //bool _aurostd_initialize_isinteger(long int x) { return _isinteger(x);}
  //bool _aurostd_initialize_isinteger(long long int x) { return _isinteger(x);}
  //bool _aurostd_initialize_isinteger(unsigned long long int x) { return _isinteger(x);}
  //bool _aurostd_initialize_isinteger(long double x) { return _isinteger(x);}

  // bool _aurostd_initialize_isinteger(bool x,bool tolerance) { return _isinteger(x,tolerance);}
  // bool _aurostd_initialize_isinteger(char x,char tolerance) { return _isinteger(x,tolerance);}
  //bool _aurostd_initialize_isinteger(int x,int tolerance) { return _isinteger(x,tolerance);}
  // bool _aurostd_initialize_isinteger(long x,long tolerance) { return _isinteger(x,tolerance);}
  // bool _aurostd_initialize_isinteger(uint x,uint tolerance) { return _isinteger(x,tolerance);}
  //bool _aurostd_initialize_isinteger(float x,float tolerance) { return _isinteger(x,tolerance);}
  //bool _aurostd_initialize_isinteger(double x,double tolerance) { return _isinteger(x,tolerance);}
  //bool _aurostd_initialize_isinteger(long int x,long int tolerance) { return _isinteger(x,tolerance);}
  //bool _aurostd_initialize_isinteger(long long int x,long long int tolerance) { return _isinteger(x,tolerance);}
  //bool _aurostd_initialize_isinteger(unsigned long long int x,unsigned long long int tolerance) { return _isinteger(x,tolerance);}
  //bool _aurostd_initialize_isinteger(long double x,long double tolerance) { return isinteger(x,_tolerance);}
}

// ***************************************************************************
// Function iszero
// ***************************************************************************
namespace aurostd {
  template<class utype> bool _iszero(utype a,utype tol) {return abs(a) <= tol;}  //CO20191201

  //bool iszero(bool x,bool tolerance){return _iszero(x,tolerance);}  //CO20191201
  //bool iszero(char x,char tolerance){return _iszero(x,tolerance);}  //CO20191201
  bool iszero(uint x,uint){return (bool) x==(uint)0;}  //CO20191201
  bool iszero(int x,int){return (bool) x==(int)0;}  //CO20191201
  bool iszero(long int x,long int){return (bool) x==(long int)0;}  //CO20191201
  bool iszero(unsigned long int x,long int){return (bool) x==(unsigned long int)0;}  //CO20191201
  bool iszero(long long int x,long long int){return (bool) x==(long long int)0;}  //CO20191201
  bool iszero(unsigned long long int x,unsigned long long int){return (bool) x==(unsigned long long int)0;}  //CO20191201
  bool iszero(float x,float tolerance){return _iszero(x,tolerance);}  //CO20191201
  bool iszero(double x,double tolerance){return _iszero(x,tolerance);}  //CO20191201
  bool iszero(long double x,long double tolerance){return _iszero(x,tolerance);}  //CO20191201
}

// ***************************************************************************
// Function factorial
// ***************************************************************************
namespace aurostd {
  // namespace aurostd
  bool factorial(bool x) {
    if(_isfloat(x)) {
      string function = XPID + "aurostd::factorial():";
      string message = "factorial(bool) implemented only for bool !";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
    }
    return TRUE;
  }

  template<class utype>
    utype factorial(utype x) {
      utype out=1;  
      for(int i=1;(utype) i<=x;i++) out=out*(utype) i;
      return out;
    }

  bool _aurostd_initialize_factorial(bool x) { return factorial(x);}
  char _aurostd_initialize_factorial(char x) { return factorial(x);}
  int _aurostd_initialize_factorial(int x) { return factorial(x);}
  uint _aurostd_initialize_factorial(uint x) { return factorial(x);}
  float _aurostd_initialize_factorial(float x) { return factorial(x);}
  double _aurostd_initialize_factorial(double x) { return factorial(x);}
  long int _aurostd_initialize_factorial(long int x) { return factorial(x);}
  long long int _aurostd_initialize_factorial(long long int x) { return factorial(x);}
  unsigned long long int _aurostd_initialize_factorial(unsigned long long int x) { return factorial(x);}
  long double _aurostd_initialize_factorial(long double x) { return factorial(x);}

  template<class utype> utype fact(utype x) { return factorial(x);}
  bool _aurostd_initialize_fact(bool x) { return factorial(x);}
  char _aurostd_initialize_fact(char x) { return factorial(x);}
  int _aurostd_initialize_fact(int x) { return factorial(x);}
  uint _aurostd_initialize_fact(uint x) { return factorial(x);}
  float _aurostd_initialize_fact(float x) { return factorial(x);}
  double _aurostd_initialize_fact(double x) { return factorial(x);}
  long int _aurostd_initialize_fact(long int x) { return fact(x);}
  long long int _aurostd_initialize_fact(long long int x) { return fact(x);}
  unsigned long long int _aurostd_initialize_fact(unsigned long long int x) { return fact(x);}
  long double _aurostd_initialize_fact(long double x) { return fact(x);}

}


// ***************************************************************************
// Function mod
// ***************************************************************************
namespace aurostd {
  // namespace aurostd
  template<class utype> utype mod(utype x,utype y) {
    utype xx;
    if(x== y) return (utype) 0.0;
    if(x== 0) return (utype) 0.0;
    if(x==-y) return (utype) 0.0;
    xx=x+y;
    // returns (xx mod y)
    //  if(xx>=0) return (utype) (xx-std::floor((double) xx/y)*y);
    // else return (utype) (y+xx-std::floor((double) xx/y)*y);
    if(xx>(utype) 0.0) {
      xx=(utype) (xx-std::floor((double) xx/y)*y);
      return xx;
    }
    else {
      xx=(utype) (y+xx-std::floor((double) xx/y)*y);  // FOR G++ 3.2
      return xx;    
    }
  }
}

// ***************************************************************************
// Function ishex
// ***************************************************************************
// ME20200707
// Quick check if the string is a hexadecimal string
namespace aurostd {
  bool _ishex(const string& hexstr) {
    uint istart = 0;
    uint str_len = hexstr.size();
    // Also process hexadecimal numbers with the '0x' prefix.
    if ((str_len > 2) && (hexstr[0] == '0') && hexstr[1] == 'x') istart = 2;
    for (uint i = istart; i < str_len; i++) {
      // Return false if chars aren't '0-9' or 'a-f' - https://www.asciitable.com/
      if ((hexstr[i] < 48) || ((hexstr[i] > 57) && (hexstr[i] < 97)) || (hexstr[i] > 102)) return false;
    }
    return true;
  }
}

// ***************************************************************************
// Function isodd
// ***************************************************************************
namespace aurostd {
  // namespace aurostd
  template<class utype>
    bool _isodd(utype x) {
      if(_isfloat(x)) {
        string function = XPID + "aurostd::factorial():";
        string message = "_isodd implemented only for integers !";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
      }
      if(!mod(x,(utype) 2)) return FALSE;
      else return TRUE;
    }
  bool _isodd(int x) {
    if(!mod(x,(int) 2)) return FALSE;
    else return TRUE;
  }
  bool _isodd(uint x) {
    if(!mod(x,(uint) 2)) return FALSE;
    else return TRUE;
  }
  bool _isodd(long int x) {
    if(!mod(x,(long int) 2)) return FALSE;
    else return TRUE;
  }
  bool _isodd(long long int x) {
    if(!mod(x,(long long int) 2)) return FALSE;
    else return TRUE;
  }
  bool _isodd(unsigned long long int x) {
    if(!mod(x,(unsigned long long int) 2)) return FALSE;
    else return TRUE;
  }
}

// ***************************************************************************
// Function iseven
// ***************************************************************************
namespace aurostd {
  // namespace aurostd
  template<class utype>
    bool _iseven(utype x) {
      if(_isfloat(x)) {
        string function = XPID + "aurostd::factorial():";
        string message = "_iseven implemented only for integers !";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
      }
      if(mod(x,(utype) 2)) return FALSE;
      else return TRUE;
    }
  bool _iseven(int x) {
    if(mod(x,(int) 2)) return FALSE;
    else return TRUE;
  }
  bool _iseven(uint x) {
    if(mod(x,(uint) 2)) return FALSE;
    else return TRUE;
  }
  bool _iseven(long int x) {
    if(mod(x,(long int) 2)) return FALSE;
    else return TRUE;
  }
  bool _iseven(long long int x) {
    if(mod(x,(long long int) 2)) return FALSE;
    else return TRUE;
  }
  bool _iseven(unsigned long long int x) {
    if(mod(x,(unsigned long long int) 2)) return FALSE;
    else return TRUE;
  }
}

// ----------------------------------------------------------------------------
// --------------------------------------------- simple operations modulus/angle
namespace aurostd {
  // namespace aurostd
  template<class utype>
    utype angle(utype x1,utype x2,utype x3,utype y1,utype y2,utype y3) {
      return (utype) std::acos((x1*y1+x2*y2+x3*y3)/(std::sqrt(x1*x1+x2*x2+x3*x3)*std::sqrt(y1*y1+y2*y2+y3*y3)))*rad2deg;
    }

  template<class utype>
    utype angle(utype x1,utype x2,utype y1,utype y2) {
      return (utype) std::acos((x1*y1+x2*y2)/(std::sqrt(x1*x1+x2*x2)*std::sqrt(y1*y1+y2*y2)))*rad2deg;
    }

  template<class utype>
    utype modulus(utype x1,utype x2,utype x3) {
      return (utype) std::sqrt(x1*x1+x2*x2+x3*x3);
    }

  template<class utype>
    utype modulus(utype x1,utype x2) {
      return (utype) std::sqrt(x1*x1+x2*x2);
    }
}

// ----------------------------------------------------------------------------
//--------------------------------------------------------------- extra_minmax min/max
namespace aurostd {
  // with const utype&
  template<class utype> bool                             // is scalar == scalar ?
    identical(const utype& a,const utype& b,const utype& _tol_) {
      if(abs(a-b)<=_tol_) return TRUE;
      return FALSE;
    }
  template<class utype> bool                             // is scalar == scalar ?
    identical(const utype& a,const utype& b) {
      return (bool) identical(a,b,(utype) _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);
    } 
  template<class utype> bool                             // is scalar != scalar ?
    isdifferent(const utype& a,const utype& b,const utype& _tol_) {
      return (bool) !identical(a,b,_tol_);
    }
  template<class utype> bool                             // is scalar != scalar ?
    isdifferent(const utype& a,const utype& b) {
      return (bool) !identical(a,b,(utype) _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);
    }
  template<class utype> bool                             // is scalar == scalar ?
    isequal(const utype& a,const utype& b,const utype& _tol_) {
      return (bool) identical(a,b,_tol_);
    }
  template<class utype> bool                             // is scalar == scalar ?
    isequal(const utype& a,const utype& b) {
      return (bool) identical(a,b,(utype) _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);
    }
  //// with utype
  //template<class utype> bool                             // is scalar == scalar ?
  //identical(utype a,utype b,utype _tol_) {
  //if(abs(a-b)<=_tol_) return TRUE;
  //return FALSE;
  //}
  //template<class utype> bool                             // is scalar == scalar ?
  //identical(utype a,utype b) {
  //return (bool) identical(a,b,(utype) _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);
  //} 
  //template<class utype> bool                             // is scalar != scalar ?
  //isdifferent(utype a,utype b,utype _tol_) {
  //return (bool) !identical(a,b,_tol_);
  //}
  //template<class utype> bool                             // is scalar != scalar ?
  //isdifferent(utype a,utype b) {
  //return (bool) !identical(a,b,(utype) _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);
  //}
  //template<class utype> bool                             // is scalar == scalar ?
  //isequal(utype a,utype b,utype _tol_) {
  //return (bool) identical(a,b,_tol_);
  //}
  //template<class utype> bool                             // is scalar == scalar ?
  //isequal(utype a,utype b) {
  //return (bool) identical(a,b,(utype) _AUROSTD_XSCALAR_TOLERANCE_IDENTITY_);
  //}
}

// ----------------------------------------------------------------------------
//--------------------------------------------------------------- extra_minmax min/max
#define _max(a,b) (a>b?a:b)
#define _min(a,b) (a<b?a:b)
// ----------------------------------------------------------------------------

namespace aurostd {
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2) {
    return _min(x1,x2);}
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2,utype x3) {
    return _min(min(x1,x2),x3);}
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2,utype x3,utype x4) {
    return _min(min(x1,x2,x3),x4);}
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2,utype x3,utype x4,utype x5) {
    return _min(min(x1,x2,x3,x4),x5);}
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6) {
    return _min(min(x1,x2,x3,x4,x5),x6);}
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7) {
    return _min(min(x1,x2,x3,x4,x5,x6),x7);}
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7,utype x8) {
    return _min(min(x1,x2,x3,x4,x5,x6,x7),x8);}
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7,utype x8,utype x9) {
    return _min(min(x1,x2,x3,x4,x5,x6,x7,x8),x9);}
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7,utype x8,utype x9,utype x10) {
    return _min(min(x1,x2,x3,x4,x5,x6,x7,x8,x9),x10);}
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7,utype x8,utype x9,utype x10,utype x11) {
    return _min(min(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10),x11);}
  // namespace aurostd
  template<class utype> utype min(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7,utype x8,utype x9,utype x10,utype x11,utype x12) {
    return _min(min(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11),x12);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2) {
    return _max(x1,x2);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2,utype x3) {
    return _max(max(x1,x2),x3);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2,utype x3,utype x4) {
    return _max(max(x1,x2,x3),x4);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2,utype x3,utype x4,utype x5) {
    return _max(max(x1,x2,x3,x4),x5);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6) {
    return _max(max(x1,x2,x3,x4,x5),x6);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7) {
    return _max(max(x1,x2,x3,x4,x5,x6),x7);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7,utype x8) {
    return _max(max(x1,x2,x3,x4,x5,x6,x7),x8);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7,utype x8,utype x9) {
    return _max(max(x1,x2,x3,x4,x5,x6,x7,x8),x9);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7,utype x8,utype x9,utype x10) {
    return _max(max(x1,x2,x3,x4,x5,x6,x7,x8,x9),x10);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7,utype x8,utype x9,utype x10,utype x11) {
    return _max(max(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10),x11);}
  // namespace aurostd
  template<class utype> utype max(utype x1,utype x2,utype x3,utype x4,utype x5,utype x6,utype x7,utype x8,utype x9,utype x10,utype x11,utype x12) {
    return _max(max(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11),x12);}
}



// ----------------------------------------------------------------------------

//ROUNDOFF for scalars
namespace aurostd {
  // namespace aurostd
  template<class utype>
    utype _roundoff(const utype& x,utype tolerance) {
      return ((abs(x)<(utype) tolerance) ? (utype) 0.0 : x);
    }

  int roundoff(int x,int tolerance){return _roundoff(x,tolerance);}
  long roundoff(long x,long tolerance){return _roundoff(x,tolerance);}
  uint roundoff(uint x,uint tolerance){return _roundoff(x,tolerance);}
  float roundoff(float x,float tolerance){return _roundoff(x,tolerance);}
  double roundoff(double x,double tolerance){return _roundoff(x,tolerance);}
  long long int roundoff(long long int x,long long int tolerance){return _roundoff(x,tolerance);}
  unsigned long int roundoff(unsigned long int x,unsigned long int tolerance){return _roundoff(x,tolerance);}
  unsigned long long int roundoff(unsigned long long int x,unsigned long long int tolerance){return _roundoff(x,tolerance);}
  long double roundoff(long double x,long double tolerance){return _roundoff(x,tolerance);}
}

namespace aurostd { //CO20190419
  int boundary_conditions_periodic(int lrows,int urows,int i){  //CO20190419 - taken from xvector BOUNDARY_CONDITIONS_PERIODIC
    int ii=i;
    if(ii==urows+1) ii=lrows;
    if(ii==lrows-1) ii=urows;
    if(ii>urows) ii=lrows+mod(i-lrows,urows-lrows+1);
    if(ii<lrows) ii=urows-mod(urows-i,urows-lrows+1);
    return ii;
  }
}

namespace aurostd { //CO20191201
  uint powint(uint x,uint exp){uint y=1;for(uint i=0;i<exp;i++){y*=x;}return y;}  //CO20191201
  int powint(int x,uint exp){int y=1;for(uint i=0;i<exp;i++){y*=x;}return y;} //CO20191201
  long int powint(long int x,uint exp){int y=1;for(uint i=0;i<exp;i++){y*=x;}return y;} //CO20191201
  unsigned long int powint(unsigned long int x,uint exp){int y=1;for(uint i=0;i<exp;i++){y*=x;}return y;} //CO20191201
  long long int powint(long long int x,uint exp){int y=1;for(uint i=0;i<exp;i++){y*=x;}return y;} //CO20191201
  unsigned long long int powint(unsigned long long int x,uint exp){int y=1;for(uint i=0;i<exp;i++){y*=x;}return y;} //CO20191201
}

//AS20200513 BEGIN
namespace aurostd{
  double FermiDirac(double E, double mu, double T){
    if (T<0) return 0;

    if (T>_ZERO_TOL_){
      return 1/(1+std::exp((E-mu)/(KBOLTZEV*T)));
    }
    //[//CO20200731 - what else is there?]else
    // At T=0 FD transforms to Heaviside step function
    if (E<mu) return 1;
    else if (E>mu) return 0;
    //CO20200731 - what else is there?else 
    return 0.5;
  }
}
//AS20200513 END

#endif // _AUROSTD_XSCALAR_CPP_



// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2020              *
// *                                                                        *
// **************************************************************************

