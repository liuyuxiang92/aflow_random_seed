// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo 1994-2011

#ifndef _AUROSTD_XVECTOR_CPP_
#define _AUROSTD_XVECTOR_CPP_

#ifndef XXEND
#define XXEND 1
#endif

#ifndef _AUROSTD_XCOMPLEX_H_
#include "aurostd_xcomplex.h"
#endif
#ifndef _AUROSTD_XSCALAR_H_
#include "aurostd_xscalar.h"
#endif
#ifndef _AUROSTD_XVECTOR_H_
#include "aurostd_xvector.h"
#endif
#ifndef _AUROSTD_XMATRIX_H_
#include "aurostd_xmatrix.h"
#endif
#ifndef _AUROSTD_XTENSOR_H_
#include "aurostd_xtensor.h"
#endif

// ----------------------------------------------------------------------------
// --------------------------------------------------------------- constructors
namespace aurostd {  // namespace aurostd
  template<class utype>                                    // default constructor
    xvector<utype>::xvector(int nh,int nl) : vsize(0) {        
      // allocate a xvector with subscript range [nl..nh]
      lrows=std::min(nl,nh);// if(!nh)lrows=0; this messes up convasp
      urows=std::max(nl,nh);// if(!nh)urows=0; this messes up convasp
      refresh(); //CO191110
#ifdef _AUROSTD_XVECTOR_DEBUG_CONSTRUCTORS
      cerr << "xxvector -> default constructor: lrows=" << lrows << ", urows=" << urows << ", rows=" << rows << endl;
#endif
      if(vsize>0) {
        corpus=new utype[rows+XXEND];
        if(!corpus) {throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::xvector<utype>::xvector():","allocation failure in default constructor",_ALLOC_ERROR_);}
        corpus+= -lrows+XXEND;
        reset(); //CO191110
      }
#ifdef _AUROSTD_XVECTOR_DEBUG_CONSTRUCTORS
      cerr << " isfloat=" << isfloat << ", iscomplex=" << iscomplex << ", sizeof=" << size << ", vsize=" << vsize << endl;
#endif
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                       // copy constructor
    xvector<utype>::xvector(const xvector<utype>& b) : vsize(0) {copy(b);} //CO191110
  template<class utype>                                       // copy constructor
    xvector<utype>::xvector(const xmatrix<utype>& b) : vsize(0) {copy(b);} //CO191110
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------- destructor

namespace aurostd {  // namespace aurostd
  template<class utype>                                     // default destructor
  xvector<utype>::~xvector() {
#ifdef _AUROSTD_XVECTOR_DEBUG_DESTRUCTORS
    cerr << "xxvector -> default destructor: lrows=" << lrows << ", urows=" << urows << ", rows=" << rows << endl;
#endif
    free(); //CO190808
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
  void xvector<utype>::free() { //CO190808
    if(vsize>0) delete [] (corpus+lrows-XXEND);
    lrows=urows=0;refresh();
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
    void xvector<utype>::copy(const xvector<utype>& b) { //CO190808
      //[CO190808 - this assumes you have corpus, perhaps it is not constructed yet]if(corpus!=b.corpus||rows!=b.rows||lrows!=b.lrows||urows!=b.urows) {     // check  for a=a */
      //CO 170803 - ODD CORNER CASE, same corpus and rows, but different lrows and urows
      //if(rows!=b.rows) {    // if dims(this)!=dims(a) => build a new xvector !!!
      if(lrows!=b.lrows||urows!=b.urows||vsize!=b.vsize) {    // if dims(this)!=dims(a) => build a new xvector !!!  //CO190808 - VERY IMPORTANT that we not only check lrows and urows, but vsize, as xvector could just have been initialized (vsize==0)
        free();
        lrows=b.lrows;
        urows=b.urows;
        //[simply copy instead]refresh();
        rows=b.rows;
        isfloat=b.isfloat;
        iscomplex=b.iscomplex;
        size=b.size;
        vsize=b.vsize;
#ifdef _AUROSTD_XVECTOR_DEBUG_CONSTRUCTORS
        printf("xxvector -> default constructor: lrows=%i, urows=%i,",lrows,urows);
#endif
        if(vsize>0) {
          corpus=new utype[rows+XXEND];
          if(!corpus) {throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::xvector<utype>::copy():","allocation failure in COPY",_ALLOC_ERROR_);}
          corpus+= -lrows+XXEND;
        }
#ifdef _AUROSTD_XVECTOR_DEBUG_CONSTRUCTORS
        printf(" isfloat=%i, iscomplex=%i, sizeof=%i, vsize=%i\n",isfloat,iscomplex,size,vsize);
#endif
        //[CO190808 - this assumes you have corpus, perhaps it is not constructed yet]}
        //[CO190808 OBSOLETE]for(int i=0;i<rows;i++)
        //[CO190808 OBSOLETE]  this->corpus[i+lrows] =
        //[CO190808 OBSOLETE]    (utype) b.corpus[i+b.lrows];
    }
      if(corpus!=b.corpus){
        for(int i=lrows;i<=urows;i++){this->corpus[i]=b.corpus[i];}
      } //CO190808 - we definitely have corpus now
    }
  template<class utype>
    void xvector<utype>::copy(const xmatrix<utype>& b) { //CO190808
      if(b.rows==1){return copy(b(b.lrows));}
      else if(b.cols==1)return copy(b.getcol(b.lcols));
      throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::xvector<utype>::copy():","xmatrix input cannot be converted to xvector",_VALUE_ILLEGAL_);
    }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
    void xvector<utype>::refresh(void) { //CO190808
      rows=urows-lrows+1;   // if(!nh) rows=0; this messes up convasp
      // [[BUG]    if(nh==0 && nl==0) {
      if(rows==0) {
        cerr << "XVECTOR constructor: creating EMPTY xvector<utype>" << endl;
        lrows=0;urows=0;rows=0;
      };
      // cerr << "XVECTOR constructor:  nh=" << nh << " nl=" << nl << " lrows=" << lrows << " urows=" << urows << " rows=" << rows << endl;
      isfloat=aurostd::_isfloat((utype) 0);      
      iscomplex=aurostd::_iscomplex((utype) 0);
      size=(char) (sizeof(utype));
      vsize=(long int) size*rows;
    }
}

// ----------------------------------------------------------------------------
// -------------------------------------------------------- assigment operators

namespace aurostd {  // namespace aurostd
  template<class utype>                                             // operator =
    xvector<utype>& xvector<utype>::operator=(const xvector<utype>& b) { //CO191110
      if(this!=&b) {copy(b);}
      return *this;
    }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------------ index operators

namespace aurostd {  // namespace aurostd
  template<class utype>                                            // operator []
  // removed inline
  utype& xvector<utype>::operator[](int i) const {
#ifndef __XOPTIMIZE
    if(i>urows) {
      cerr << "xvector[1]=" << corpus[1] << endl;
      cerr << "xvector[2]=" << corpus[2] << endl;
      cerr << "xvector[3]=" << corpus[3] << endl;
      cerr << _AUROSTD_XLIBS_ERROR_ << " xvector[] -> i=" << i << " > urows=" << urows << " lrows=" << lrows << " float=" << isfloat << endl;
      exit(0);
    }
    if(i<lrows) {
      cerr << "xvector[1]=" << corpus[1] << endl;
      cerr << "xvector[2]=" << corpus[2] << endl;
      cerr << "xvector[3]=" << corpus[3] << endl;
      cerr << _AUROSTD_XLIBS_ERROR_ << " xvector[] -> i=" << i << " < lrows=" << lrows << " urows=" << urows << " float=" << isfloat << endl;
      exit(0);
    }
#endif
    return corpus[i];
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                            // operator ()
  // removed inline
  utype& xvector<utype>::operator()(int i) const {
    //#ifndef BOUNDARY_CONDITIONS_PERIODIC
#ifndef __XVECTOR_IGNORE_BOUNDARIES
    if(i>urows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << " xvector() -> i=" << i << " > urows=" << urows << " lrows=" << lrows << " float=" << isfloat << endl;
      exit(0);
    }
    if(i<lrows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << " xvector() -> i=" << i << " < lrows=" << lrows << " urows=" << urows << " float=" << isfloat << endl;
      exit(0);
    }
#endif  // __XVECTOR_IGNORE_BOUNDARIES
    return corpus[i];
  }
}
/*
  #else
  #ifdef XVECTOR_WARNING
  #warning "BOUNDARY_CONDITIONS_PERIODIC"
  #endif
  int ii=i;
  // if(ii>urows) ii=lrows+mod(i-lrows,urows-lrows+1);
  // if(ii<lrows) ii=urows-mod(urows-i,urows-lrows+1);
  if(ii>urows) ii-=rows;
  if(ii<lrows) ii+=rows;
  return corpus[ii];
  #endif
  }
  }
*/

// ----------------------------------------------------------------------------
// ----------------------------------- index operators with boundary conditions

namespace aurostd {  // namespace aurostd
  template<class utype>                        // operator () boundary conditions
  inline utype& xvector<utype>::operator()(int i,bool bc) const {
    if(bc==BOUNDARY_CONDITIONS_NONE) {
#ifndef __XOPTIMIZE
      if(i>urows) {
	cerr << _AUROSTD_XLIBS_ERROR_ << " i > xvector<utype>.urows, BC=" << bc << endl;
	exit(0);
      }
      if(i<lrows) {
	cerr << _AUROSTD_XLIBS_ERROR_ << " i < xvector<utype>.lrows, BC=" << bc << endl;
	exit(0);
      }
#endif
      return corpus[i];
    }
    if(bc==BOUNDARY_CONDITIONS_PERIODIC) {
      int ii=boundary_conditions_periodic(lrows,urows,i); //CO190520
      //[CO190419 - moved to xscalar]int ii=i; // ,jj=j; //CO190520
      //[CO190419 - moved to xscalar]if(ii==urows+1) ii=lrows; //CO190520
      //[CO190419 - moved to xscalar]if(ii==lrows-1) ii=urows; //CO190520
      //[CO190419 - moved to xscalar]if(ii>urows) ii=lrows+mod(i-lrows,urows-lrows+1); //CO190520
      //[CO190419 - moved to xscalar]if(ii<lrows) ii=urows-mod(urows-i,urows-lrows+1); //CO190520
      return corpus[ii];
    }
    return corpus[i];  //CO190419 - needs to return something
  }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------- math unary operators

// --------------------------------------------------------- operator += xvector
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>&
  // removed inline
  xvector<utype>::operator +=(const xvector<utype>& r)
  {
#ifdef _AUROSTD_XVECTOR_DEBUG_OPERATORS
    printf("xxvector -> operator xvector+=xvector: ");
    printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
    printf("r.lrows=%i, r.urows=%i\n",r.lrows,r.urows);
#endif
    if(this->rows!=r.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xvector<utype>: failure in operator+=:" << " (this->rows!=r.rows)=FALSE" << endl;
      exit(0);
    }
    for(int i=0;i<rows;i++)
      corpus[i+lrows]+=r[i+r.lrows];
    return *this;
  }
}

// --------------------------------------------------------- operator += xvector
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>&
  // removed inline
  xvector<utype>::operator +=(utype r)
  {
    for(int i=0;i<rows;i++)
      corpus[i+lrows]+=r;
    return *this;
  }
}

// --------------------------------------------------------- operator -= xvector
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>&  
  // removed inline
  xvector<utype>::operator -=(const xvector<utype>& r)
  {
#ifdef _AUROSTD_XVECTOR_DEBUG_OPERATORS
    printf("xxvector -> operator xvector-=xvector: ");
    printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
    printf("r.lrows=%i, r.urows=%i\n",r.lrows,r.urows);
#endif
    if(this->rows!=r.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xvector<utype>: failure in operator-=:" << " (this->rows!=r.rows)=FALSE" << endl;
      exit(0);
    }
    for(int i=0;i<rows;i++)
      corpus[i+lrows]-=r[i+r.lrows];
    return *this;
  }
}

// --------------------------------------------------------- operator -= xvector
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>&  
  // removed inline
  xvector<utype>::operator -=(utype r)
  {
    for(int i=0;i<rows;i++)
      corpus[i+lrows]-=r;
    return *this;
  }
}

// --------------------------------------------------------- operator *= double
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>&
  // removed inline
  xvector<utype>::operator *=(utype r) //(const utype& r) //CO 171130 - v*=v[1] doesn't work
  {
#ifdef _AUROSTD_XVECTOR_DEBUG_OPERATORS
    printf("xxvector -> operator xvector*=xvector: ");
    printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
    printf("r.lrows=%i, r.urows=%i\n",r.lrows,r.urows);
#endif
    for(int i=0;i<rows;i++)
      corpus[i+lrows]*=r;
    return *this;
  }
}

// --------------------------------------------------------- operator /= double
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>&
  // removed inline
  xvector<utype>::operator /=(utype r) //(const utype& r) //CO 171130 - v/=v[1] doesn't work
  {
#ifdef _AUROSTD_XVECTOR_DEBUG_OPERATORS
    printf("xxvector -> operator xvector/=xvector: ");
    printf("this->lrows=%i, this->urows=%i, ",this->lrows,this->urows);
    printf("r.lrows=%i, r.urows=%i\n",r.lrows,r.urows);
#endif
    for(int i=0;i<rows;i++)
      corpus[i+lrows]/=r;
    return *this;
  }
}

// ---------------------------------------------------------- operator+xvector
namespace aurostd {  // namespace aurostd
  template<class utype>                          
  xvector<utype> operator+(const xvector<utype>& a) {
    return a;
  }
}

// ---------------------------------------------------------- operator-xvector
namespace aurostd {  // namespace aurostd
  template<class utype>                          
  xvector<utype> operator-(const xvector<utype>& a) {
    xvector<utype> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      c[i]=-a[i];
    return c;
  }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------ math binary operators

// SCALAR PRODUCT OPERATOR
namespace aurostd {  // namespace aurostd
  template<class utype> utype                         // operator xvector * xvector
  operator*(const xvector<utype>& a,const xvector<utype>& b) {
    if(a.rows!=b.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xvector<utype>: failure in operator*" << endl;
      exit(0);
    }
    utype out=(utype) 0.0;
    for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
      out += a[i]*b[ii];
    return (utype) out;
  }
}

// SCALAR PRODUCT FUNCTION
namespace aurostd {  // namespace aurostd
  template<class utype> utype                   // scalar_product xvector * xvector
  scalar_product(const xvector<utype>& a,const xvector<utype>& b) {
    if(a.rows!=b.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xvector<utype>: failure in operator*" << endl;
      exit(0);
    }
    utype out=(utype) 0.0;
    for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++)
      out += a[i]*b[ii];
    return (utype) out;
  }
}

// VECTOR PRODUCT OPERATOR
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                 // operator xvector % xvector
  operator %(const xvector<utype>& a,const xvector<utype>& b) {
    xvector<utype> c(3);
    if(a.rows!=3) {
      cerr << _AUROSTD_XLIBS_ERROR_ << " XVECTOR.CPP: xvector product (a%b) a.rows=" << a.rows << " !=3" << endl;
      exit(0);
    }
    if(b.rows!=3) {
      cerr << _AUROSTD_XLIBS_ERROR_ << " XVECTOR.CPP: xvector product (a%b) b.rows=" << b.rows << " !=3" << endl;
      exit(0);
    }
    c(1)=a(2)*b(3)-a(3)*b(2);
    c(2)=a(3)*b(1)-a(1)*b(3);
    c(3)=a(1)*b(2)-a(2)*b(1);
    return c;
  }
}

// VECTOR PRODUCT FUNCTION
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>          // vector_product xvector % xvector
  vector_product(const xvector<utype>& a,const xvector<utype>& b) {
    xvector<utype> c(3);
    if(a.rows!=3) {
      cerr << _AUROSTD_XLIBS_ERROR_ << " XVECTOR.CPP: xvector product (a%b) a.rows=" << a.rows << " !=3" << endl;
      exit(0);
    }
    if(b.rows!=3) {
      cerr << _AUROSTD_XLIBS_ERROR_ << " XVECTOR.CPP: xvector product (a%b) b.rows=" << b.rows << " !=3" << endl;
      exit(0);
    }

    c(1)=a(2)*b(3)-a(3)*b(2);
    c(2)=a(3)*b(1)-a(1)*b(3);
    c(3)=a(1)*b(2)-a(2)*b(1);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                 // operator xvector+xvector
  operator+(const xvector<utype>& a,const xvector<utype>& b) {
    //xvector<utype> c(std::max(a.rows,b.rows));  //CO 170910 - doesn't work with 0 starting index
    xvector<utype> c(std::max(a.urows,b.urows),std::min(a.lrows,b.lrows));  //CO 170910 - max/min
    for(int i=0;i<c.rows;i++)
      if(a.urows-i>=a.lrows && b.urows-i>=b.lrows)
	c[c.urows-i]=a[a.urows-i]+b[b.urows-i];
      else if(a.urows-i>=a.lrows)
	c[c.urows-i]=a[a.urows-i];
      else c[c.urows-i]=b[b.urows-i];
    //   xvector<utype> cc(c); return cc; VISUALAGE ???
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                 // operator xvector-xvector
  operator-(const xvector<utype>& a,const xvector<utype>& b) {
    //xvector<utype> c(std::max(a.rows,b.rows));  //CO 170910 - doesn't work with 0 starting index
    xvector<utype> c(std::max(a.urows,b.urows),std::min(a.lrows,b.lrows));  //CO 170910 - max/min
    for(int i=0;i<c.rows;i++)
      if(a.urows-i>=a.lrows && b.urows-i>=b.lrows)
	c[c.urows-i]=a[a.urows-i]-b[b.urows-i];
      else if(a.urows-i>=a.lrows)
	c[c.urows-i]=a[a.urows-i];
      else c[c.urows-i]=-b[b.urows-i];
    return c;
  }
}

// ------------------------------------------------------ operator scalar+xvector
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                 // operator scalar+xvector
  operator+(const utype s,const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) s+a[i];
    return c;
  }
  // namespace aurostd
  template<class utype> xvector<utype>                 // operator xvector+scalar
  operator+(const xvector<utype>& a,const utype s) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) a[i]+s;
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype,class stype> xvector<utype>     // operator scalar+xvector
  operator+(const stype s,const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) s+a[i];
    return c;
  }
  // namespace aurostd
  template<class utype,class stype> xvector<utype>     // operator xvector+scalar
  operator+(const xvector<utype>& a,const stype s) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) a[i]+s;
    return c;
  }
}

// ------------------------------------------------------ operator scalar-xvector
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                 // operator scalar-xvector
  operator-(const utype s,const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) s-a[i];
    return c;
  }
  // namespace aurostd
  template<class utype> xvector<utype>                 // operator xvector-scalar
  operator-(const xvector<utype>& a,const utype s) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) a[i]-s;
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype,class stype> xvector<utype>     // operator scalar-xvector
  operator-(const stype s,const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) s-a[i];
    return c;
  }
  // namespace aurostd
  template<class utype,class stype> xvector<utype>     // operator xvector-scalar
  operator-(const xvector<utype>& a,const stype s) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) a[i]-s;
    return c;
  }
}

// ------------------------------------------------------ operator scalar * xvector
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                 // operator scalar * xvector
  operator*(const utype s,const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) s*a[i];
    return c;
  }
  // namespace aurostd
  template<class utype> xvector<utype>                 // operator xvector * scalar
  operator*(const xvector<utype>& a,const utype s) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) a[i]*((utype)s);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype,class stype> xvector<utype>     // operator scalar * xvector
  operator*(const stype s,const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=((utype) s)*a[i];
    return c;
  }
  // namespace aurostd
  template<class utype,class stype> xvector<utype>     // operator scalar * xvector
  operator*(const xvector<utype>& a,const stype s) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) a[i]*((utype) s);
    return c;
  }
}

// ------------------------------------------------------ operator scalar / xvector
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                 // operator scalar / xvector
  operator/(const utype s,const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) s/a[i];
    return c;
  }
  // namespace aurostd
  template<class utype> xvector<utype>                 // operator xvector / scalar
  operator/(const xvector<utype>& a,const utype s) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=(utype) a[i]/s;
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype,class stype> xvector<utype>     // operator scalar / xvector
  operator/(const stype s,const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=((utype) s)/a[i];
    return c;
  }
  // namespace aurostd
  template<class utype,class stype> xvector<utype>     // operator xvector / scalar
  operator/(const xvector<utype>& a,const stype s) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) c[i]=a[i]/((utype) s);
    return c;
  }
}

// --------------------------------------------------------------------------------
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>               // operator xvector << xvector
  operator<<(const xvector<utype>& a,const xvector<utype>& b) {
    xvector<utype> c(a.rows+b.rows);
    for(int i=1;i<=c.rows;i++)
      if(i<=a.rows) c[i]=a[a.lrows-1+i];
      else c[i]=b[b.lrows-1+i-a.rows];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>               // operator xvector << scalar
  operator<<(const xvector<utype>& a,const utype s) {
    xvector<utype> c(a.rows+1);
    for(int i=1;i<=a.rows;i++)
      c[i]=a[a.lrows-1+i];
    c[a.rows+1]=(utype) s;
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>               // operator scalar << xvector
  operator<<(const utype s,const xvector<utype>& a) {
    xvector<utype> c(a.rows+1);
    c[1]=(utype) s;
    for(int i=1;i<=a.rows;i++)
      c[i+1]=a[a.lrows-1+i];
    return c;
  }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------------- conditionals !

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<char>                    // is xvector > scalar ?
  operator>(const xvector<utype>& a,const utype& b) {
    xvector<char> c(a.lrows,a.urows);
    for(int i=a.lrows;i<=a.urows;i++) {
      c[i]=FALSE;
      if(a[i]>b) c[i]=TRUE;
    }
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<char>                    // is xvector < scalar ?
  operator<(const xvector<utype>& a,const utype& b) {
    xvector<char> c(a.lrows,a.urows);
    for(int i=a.lrows;i<=a.urows;i++) {
      c[i]=FALSE;
      if(a[i]<b) c[i]=TRUE;
    }
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<char>                    // is xvector == scalar ?
  operator==(const xvector<utype>& a,const utype& b) {
    xvector<char> c(a.lrows,a.urows);
    for(int i=a.lrows;i<=a.urows;i++) {
      c[i]=FALSE;
      if(a[i]==b) c[i]=TRUE;
    }
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<char>                    // is xvector > xvector ?
  operator>(const xvector<utype>& a,const xvector<utype>& b) {
    if(a.rows!=b.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xvector<utype>: failure in operator> (xvector > xvector)" << endl;
      exit(0);
    }
    xvector<char> c(a.lrows,a.urows);
    for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++) {
      c[i]=FALSE;
      if(a[i]>b[ii]) c[i]=TRUE;
    }
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<char>                    // is xvector < xvector ?
  operator<(const xvector<utype>& a,const xvector<utype>& b) {
    if(a.rows!=b.rows)  {
      cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xvector<utype>: failure in operator< (xvector < xvector)" << endl;
      exit(0);
    }
    xvector<char> c(a.lrows,a.urows);
    for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++) {
      c[i]=FALSE;
      if(a[i]<b[ii]) c[i]=TRUE;
    }
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> bool                             // is xvector == xvector ?
  identical(const xvector<utype>& a,const xvector<utype>& b,const utype& _tol_) {
    if(a.rows!=b.rows)  {
      cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xvector<utype>: failure in function identical (xvector == xvector)[1]" << endl;
      exit(0);
    }
    bool output=TRUE;
    if(a.isfloat || a.iscomplex) {
      for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++) {
	if((abs(a[i]-b[ii])/(abs(a[i])+abs(b[ii])+_tol_))>_tol_)  output=FALSE; // SC 20180115
	// output=output*(((abs(a[i]-b[ii]))/(abs(a[i])+abs(b[ii])+_tol_))<=_tol_);// SC 20180115
	// output=output*(abs(a[i]-b[ii])<=_tol_); // SC pre 20180115
      }
      if(output==FALSE) return (bool) output;
    } else {
      for(int i=a.lrows,ii=b.lrows;i<=a.urows;i++,ii++) {
	if(a[i]!=b[ii]) output=FALSE; // SC 20180115
	// output=output*(a[i]==b[ii]); // SC 20180115
      }
      if(output==FALSE) return (bool) output;
    }    
    return (bool) output;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> bool                             // is xvector == xvector ?
  identical(const xvector<utype>& a,const xvector<utype>& b) {
    return (bool) identical(a,b,(utype) _AUROSTD_XVECTOR_TOLERANCE_IDENTITY_);
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> bool                             // is xvector == xvector ?
  operator==(const xvector<utype>& a,const xvector<utype>& b) {
    return (bool) identical(a,b,(utype) _AUROSTD_XVECTOR_TOLERANCE_IDENTITY_);
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> bool                             // is xvector!=xvector ?
  isdifferent(const xvector<utype>& a,const xvector<utype>& b,const utype& _tol_) {
    return (bool) !identical(a,b,_tol_);
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> bool                             // is xvector!=xvector ?
  isdifferent(const xvector<utype>& a,const xvector<utype>& b) {
    return (bool) !identical(a,b,(utype) _AUROSTD_XVECTOR_TOLERANCE_IDENTITY_);
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> bool                             // is xvector == xvector ?
  isequal(const xvector<utype>& a,const xvector<utype>& b,const utype& _tol_) {
    return (bool) identical(a,b,_tol_);
  }
  bool _aurostd_initialize_isequal(const xvector<int>& a,const xvector<int>& b,const int& _tol_) { return isequal(a,b,_tol_);}
  bool _aurostd_initialize_isequal(const xvector<uint>& a,const xvector<uint>& b,const uint& _tol_) { return isequal(a,b,_tol_);} //CO 180409
  bool _aurostd_initialize_isequal(const xvector<float>& a,const xvector<float>& b,const float& _tol_) { return isequal(a,b,_tol_);}
  bool _aurostd_initialize_isequal(const xvector<double>& a,const xvector<double>& b,const double& _tol_) { return isequal(a,b,_tol_);}
}

namespace aurostd {  // namespace aurostd
  template<class utype> bool                             // is xvector == xvector ?
  isequal(const xvector<utype>& a,const xvector<utype>& b) {
    return (bool) identical(a,b,(utype) _AUROSTD_XVECTOR_TOLERANCE_IDENTITY_);
  }
  bool _aurostd_initialize_isequal(const xvector<int>& a,const xvector<int>& b) { return isequal(a,b);}
  bool _aurostd_initialize_isequal(const xvector<uint>& a,const xvector<uint>& b) { return isequal(a,b);} //CO 180409
  bool _aurostd_initialize_isequal(const xvector<float>& a,const xvector<float>& b) { return isequal(a,b);}
  bool _aurostd_initialize_isequal(const xvector<double>& a,const xvector<double>& b) { return isequal(a,b);}
}

namespace aurostd {  // namespace aurostd
  template<class utype> bool                             // is xvector!=xvector ?
  operator!=(const xvector<utype>& a,const xvector<utype>& b) {
    return (bool) !identical(a,b,(utype) _AUROSTD_XVECTOR_TOLERANCE_IDENTITY_);
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> bool
  isinteger(const xvector<utype>& a,const utype& tol) { //CO 180409
    for(int i=a.lrows;i<=a.urows;i++)
      if(isinteger(a[i],tol)==FALSE) return FALSE; //CO 180409
    return TRUE;
  }
  bool _aurostd_initialize_isinteger(const xvector<int>& a) { return isinteger(a);}
  bool _aurostd_initialize_isinteger(const xvector<uint>& a) { return isinteger(a);} //CO 180409
  bool _aurostd_initialize_isinteger(const xvector<float>& a) { return isinteger(a);}
  bool _aurostd_initialize_isinteger(const xvector<double>& a) { return isinteger(a);}
}

// ME 180702 - Tests if vector is the zero vector
namespace aurostd {
  template<class utype> bool
  iszero(const xvector<utype>& a, const double& tol) {
    for (int i = a.lrows; i <= a.urows; i++) {
      if (abs(a(i)) > tol) {
        return false;
      }
    }
    return true; 
  }
  bool _aurostd_initialize_iszero(const xvector<int>& a) { return iszero(a);}
  bool _aurostd_initialize_iszero(const xvector<uint>& a) { return iszero(a);}
  bool _aurostd_initialize_iszero(const xvector<float>& a) { return iszero(a);}
  bool _aurostd_initialize_iszero(const xvector<double>& a) { return iszero(a);}
}

// ------------------------------------------------------------ std::cout operations

/*
  namespace aurostd {  // namespace aurostd
  template<class utype>
  std::ostream& operator<< (std::ostream& buf,const xvector<utype>& x) {
  char buf1[80],*iobuf;                                    
  utype xi;
  if(!x.isfloat) {
  if(x.size==sizeof(char))
  iobuf= " %4d";        
  if(x.size==sizeof(int))
  iobuf= " %4i";        
  if(x.size==sizeof(long int))
  iobuf= " %9i";        
  } else {
  if(!x.iscomplex && x.size==sizeof(float))
  iobuf= " %7.3le";
  if(!x.iscomplex && x.size==sizeof(double))
  iobuf= " %11.4le";
  if(!x.iscomplex && x.size==sizeof(long double))
  iobuf= " %13.5le";
  if(x.iscomplex && x.size==sizeof(complex<float>))
  iobuf= " (%7.3le,%7.3le)";
  if(x.iscomplex && x.size==sizeof(complex<double>))
  iobuf= " (%11.4le,%11.4le)";
  if(x.iscomplex && x.size==sizeof(complex<long double>))
  iobuf= " (%13.5le,%13.5le)";
  }
  for(int i=x.lrows;i<=x.urows;i++) {
  xi=x[i];
  if(!x.isfloat) {  
  sprintf(buf1,"[%2i]  ",i);
  buf << buf1;
  }
  //    sprintf(buf1,iobuf+((xi >=0) ? 0 :1 ),xi);
  if(!x.iscomplex)
  sprintf(buf1,iobuf,xi);
  else
  sprintf(buf1,iobuf,real(xi),imag(xi));
  buf << buf1;
  if(!x.isfloat)
  buf << endl;
  else
  if(i<x.urows) buf << "  ";
  }
  return buf;
  }
  }
*/

namespace aurostd {  // namespace aurostd
  template<class utype>                            // operator <<  xvector<>
  std::ostream& operator<< (std::ostream& buf,const xvector<utype>& x) {
    char buf1[80];
    string iobuf;                                    
    utype xi=0;
    if(!aurostd::_isfloat(xi)) {
      if(aurostd::_size(xi)==sizeof(char))
	iobuf= "%4d";        
      if(aurostd::_size(xi)==sizeof(int))
	//	iobuf= " %4i";        
	iobuf="%11i";        
      if(aurostd::_size(xi)==sizeof(long int))
	iobuf= "%9i";        
      if(aurostd::_size(xi)==sizeof(uint)) //CO 180409
        //	iobuf= " %4i";        
        iobuf="%11i";        
      if(aurostd::_size(xi)==sizeof(unsigned long int)) //CO 180409
	iobuf= "%9i";        
    } else {
      if(!aurostd::_iscomplex(xi) && aurostd::_size(xi)==sizeof(float))
	iobuf="%7.3le";
      if(!aurostd::_iscomplex(xi) && aurostd::_size(xi)==sizeof(double))
	iobuf="%11.4le";
      if(!aurostd::_iscomplex(xi) && aurostd::_size(xi)==sizeof(long double))
	iobuf="%13.5le";
      //    if(aurostd::_iscomplex(xi) && aurostd::_size(xi)==sizeof(complex<float>))
      //      iobuf= " (%7.3le,%7.3le)";
      //    if(aurostd::_iscomplex(xi) && aurostd::_size(xi)==sizeof(complex<double>))
      //      iobuf= " (%11.4le,%11.4le)";
      //    if(aurostd::_iscomplex(xi) && aurostd::_size(xi)==sizeof(complex<long double>))
      //      iobuf= " (%13.5le,%13.5le)";
    }
    if(x.rows>0) {
	for(int i=x.lrows;i<=x.urows;i++) {
	  xi=x[i];
	  // if(!aurostd::_isfloat(xi)){                              // for index [n]
	  //	sprintf(buf1,"[%2i]  ",i);                // for index [n]
	  //	buf << buf1;                              // for index [n]
	  // }                                                // for index [n]
	  //    sprintf(buf1,iobuf.c_str()+((xi >=0) ? 0 :1 ),xi);
	  if(!aurostd::_iscomplex(xi)) {
	    sprintf(buf1,iobuf.c_str(),aurostd::_real(xi));
	    buf << buf1;
	  } else {
	    //      sprintf(buf1,iobuf.c_str(),real(xi),imag(xi));
	    buf << xi << "";  // problem is in << of xcomplex
	  }
	  //      if(!aurostd::_isfloat(xi))                          // remove newline for integers
	  //	buf << endl;                              // remove newline for integers
	  //  else                                            // remove newline for integers
	  if(i<x.urows) buf << " ";
	}
    } else {
      buf << " xvector=empty ";
    }
    return buf;
  }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------ xvector constrtuction
// reshape from scalars
namespace aurostd {  
  template<class utype>
  xvector<utype> reshape(const utype& s1) {
    xvector<utype> v(1);
    v(1)=s1;
    return v;
  }

  template<class utype>
  xvector<utype> reshape(const utype& s1,const utype& s2) {
    xvector<utype> v(2);
    v(1)=s1;v(2)=s2;
    return v;
  }

  template<class utype>
  xvector<utype> reshape(const utype& s1,const utype& s2,const utype& s3) {
    xvector<utype> v(3);
    v(1)=s1;v(2)=s2;v(3)=s3;
    return v;
  }

  template<class utype>
  xvector<utype> reshape(const utype& s1,const utype& s2,const utype& s3,const utype& s4) {
    xvector<utype> v(4);
    v(1)=s1;v(2)=s2;v(3)=s3;v(4)=s4;
    return v;
  }

  template<class utype>
  xvector<utype> reshape(const utype& s1,const utype& s2,const utype& s3,const utype& s4,const utype& s5) {
    xvector<utype> v(5);
    v(1)=s1;v(2)=s2;v(3)=s3;v(4)=s4;v(5)=s5;
    return v;
  }

  template<class utype>
  xvector<utype> reshape(const utype& s1,const utype& s2,const utype& s3,const utype& s4,const utype& s5,const utype& s6) {
    xvector<utype> v(6);
    v(1)=s1;v(2)=s2;v(3)=s3;v(4)=s4;v(5)=s5;v(6)=s6;
    return v;
  }

}  

namespace aurostd {
  template<class utype> xvector<utype> ones_xv(int nh,int nl) __xprototype { //CO190419
    xvector<utype> a(nh,nl);
    for(int i=a.lrows;i<=a.urows;i++){a[i]=(utype)1;}
    return a;
  }
  template<class utype> xvector<utype> box_filter_xv(int window,int lrows) __xprototype {  //CO190419
    xvector<utype> filter=aurostd::ones_xv<utype>(window+(lrows-1),lrows);filter/=(utype)window;
    return filter;
  }
#define STDDEV_TRUNCATE_GAUSSIAN 4 //after 4 stddev's, gaussian is effectively 0
  template<class utype> xvector<utype> gaussian_filter_xv(utype sigma) __xprototype { //CO190419
    int half_width=(int) (sigma * STDDEV_TRUNCATE_GAUSSIAN);
    int window=2*half_width+1;
    if(window%2==0){window++;}
    return gaussian_filter_xv<utype>(sigma,window); //if you need lrows!=1, use shiftlrows()
  }
  template<class utype> xvector<utype> gaussian_filter_xv(utype sigma,int window,int lrows) __xprototype { //CO190419
    if(window%2==0){throw aurostd::xerror(_AFLOW_FILE_NAME_,"aurostd::gaussian_filter_xv():","window should NOT be even (window="+aurostd::utype2string(window)+")",_INPUT_ILLEGAL_);}
    xvector<utype> filter(window+(lrows-1),lrows);
    int ind=lrows;utype x=0.0;
    for(int val=-window/2;val<=window/2;val++){
      x=(utype)val; //will look like -3 -2 -1 0 1 2 3
      filter[ind++]=std::exp(-(x*x)/((utype)2*sigma*sigma));
    }
    filter/=(utype)aurostd::sum(filter);
    return filter;
  }
}

// ----------------------------------------------------------------------------
// -------------------------------------------------------------- xvector casts

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to long double
  xvector<long double> xlongdouble(const xvector<utype> &a) {
    xvector<long double> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      c[i]=(long double) a[i];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to double
  xvector<double> xdouble(const xvector<utype> &a) {
    xvector<double> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      c[i]=(double) a[i];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                          // function floor
  xvector<double> floor(const xvector<utype> &a) {
    xvector<double> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      c[i]=std::floor(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                           // function ceil
  xvector<double> ceil(const xvector<utype> &a) {
    xvector<double> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      c[i]=std::ceil(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                           // function round
  xvector<double> round(const xvector<utype> &a) {
    xvector<double> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      // c[i]=std::round(a[i]);
      c[i]=round(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                           // function trunc
  xvector<double> trunc(const xvector<utype> &a) {
    xvector<double> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      //   c[i]=std::trunc(a[i]);
       c[i]=trunc(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to float
  xvector<float> xfloat(const xvector<utype> &a) {
    xvector<float> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      c[i]=(float) a[i];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to long int
  xvector<long int> xlongint(const xvector<utype> &a) {
    xvector<long int> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      c[i]=(long int) a[i];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to int
  xvector<int> xint(const xvector<utype> &a) {
    xvector<int> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      c[i]=(int) a[i];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // conversion to char
  xvector<char> xchar(const xvector<utype> &a) {
    xvector<char> c(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++)
      c[i]=(char) a[i];
    return c;
  }
}

namespace aurostd {                   // conversion to vector<utype>
  template<class utype> vector<utype>
  xvector2vector(const xvector<utype> & xvec) {
    int isize=xvec.rows;
    vector<utype> vvector(isize);
    for(int i=0;i<isize;i++)
      vvector[i]=xvec(i+xvec.lrows);
    return vvector;
  }
}

namespace aurostd {                   // conversion to xvector<utype>
  template<class utype> xvector<utype>
  vector2xvector(const vector<utype> & vec,int lrows) { //CO 180409
    int isize=vec.size();
    xvector<utype> xvec(isize+lrows-1,lrows); //CO 180409
    for(int i=lrows;i<=isize+lrows-1;i++) //CO 180409
      xvec(i)=vec.at(i-lrows); //CO 180409
    return xvec;
  }
}

//CO190516
namespace aurostd {                   // conversion from xvector<int> to xvector<double>
  xvector<double> xvectorint2double(const xvector<int>& a){
    xvector<double> b(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++){b[i]=(double)a[i];}  //nint is for safety
    return b;
  }
}

//CO190516
namespace aurostd {                   // conversion from xvector<int> to xvector<double>
  xvector<int> xvectordouble2int(const xvector<double>& a){
    xvector<int> b(a.urows,a.lrows);
    for(int i=a.lrows;i<=a.urows;i++){b[i]=(int)nint(a[i]);}  //nint is for safety
    return b;
  }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------- set reset operations

namespace aurostd {  // namespace aurostd
  template<class utype>
  void xvector<utype>::reset(void) {
    for(int i=lrows;i<=urows;i++) corpus[i]=(utype) 0.0;
  }
  template<class utype>                              // function reset xvector<>
  void reset(xvector<utype>& a) {
    for(int i=a.lrows;i<=a.urows;i++) a[i]=(utype) 0.0;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
  void xvector<utype>::clear(void) {
    //[CO190808 - this is ideal behavior of clear, but to avoid seg faults with size changes, simply reset() instead]xvector<utype> a;copy(a);
    reset(); //CO191110
  }
  template<class utype>                              // function clear xvector<>
  void clear(xvector<utype>& b) {
    //[CO190808 - this is ideal behavior of clear, but to avoid seg faults with size changes, simply reset() instead]xvector<utype> a;b=a;
    b.reset(); //CO191110
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>
  void xvector<utype>::set(const utype& s) {
    for(int i=lrows;i<=urows;i++) corpus[i]=(utype) s;
  }
  template<class utype>                                // function set xvector<>
  void set(xvector<utype>& a, const utype& s) {
    for(int i=a.lrows;i<=a.urows;i++) a[i]=(utype) s;
  }
}

// ----------------------------------------------------------------------------
// --------------------------------------- abs/sign/modulus/sum/nint operations

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                // function vabs xvector<>
  vabs(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) aurostd::abs(a[i]);  // (if any abs is defined on utype)
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                // function abs xvector<>
  abs(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) aurostd::abs(a[i]);  // (if any abs is defined on utype)
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                // function sign xvector<>
  sign(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) aurostd::sign(a[i]);  // (if any abs is defined on utype)
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                // function nint xvector<>
  nint(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) aurostd::nint(a[i]);  // (if any nint is defined on utype)
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                            // function modulus xvector<>
  utype modulus(const xvector<utype>& a) {
    utype c=(utype) 0.0;
    for(int i=a.lrows;i<=a.urows;i++)
      //      c+=(utype) std::abs(a[i]*a[i]);   // ABS FIX
      c+=(utype) abs(a[i]*a[i]);
    c=(utype) std::sqrt((double) c);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                           // function modulussquare xvector<>
  utype modulussquare(const xvector<utype>& a) {
    utype c=(utype) 0.0;
    for(int i=a.lrows;i<=a.urows;i++)
      //     c+=(utype) std::abs(a[i]*a[i]); // ABS FIX
      c+=(utype) abs(a[i]*a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                           // function modulus2 xvector<>
  utype modulus2(const xvector<utype>& a) {
    return modulussquare(a);
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                // function sum xvector<>
  utype sum(const xvector<utype>& a) {
    utype c=a[a.lrows];
    for(int i=a.lrows+1;i<=a.urows;i++)
      c+=a[i];
    return c;
  }
}

// ----------------------------------------------------------------------------
// --------------------------------------------------------- min max operations

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function min xvector<>
  utype min(const xvector<utype>& a) {
    utype c=a[a.lrows];
    for(int i=a.lrows+1;i<=a.urows;i++)
      c=c < a[i] ? c:a[i];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function min xvector<>
  utype min(const xvector<utype>& a,int& index) {
    utype c=a[a.lrows];
    index=a.lrows;
    for(int i=a.lrows+1;i<=a.urows;i++)
      if(a[i]<c) {
	c=a[i];
	index=i;
      }
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                // function mini xvector<>
  int mini(const xvector<utype>& a) {
    utype c=a[a.lrows];
    int index=a.lrows;
    for(int i=a.lrows+1;i<=a.urows;i++)
      if(a[i]<c) {
	c=a[i];
	index=i;
      }
    return index;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function max xvector<>
  utype max(const xvector<utype>& a) {
    utype c=a[a.lrows];
    for(int i=a.lrows+1;i<=a.urows;i++)
      c=c > a[i] ? c:a[i];
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                 // function max xvector<>
  utype max(const xvector<utype>& a,int& index) {
    utype c=a[a.lrows];
    index=a.lrows;
    for(int i=a.lrows+1;i<=a.urows;i++)
      if(a[i]>c) {
	c=a[i];
	index=i;
      }
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype>                                // function maxi xvector<>
  int maxi(const xvector<utype>& a) {
    utype c=a[a.lrows];
    int index=a.lrows;
    for(int i=a.lrows+1;i<=a.urows;i++)
      if(a[i]>c) {
	c=a[i];
	index=i;
      }
    return index;
  }
}

// ----------------------------------------------------------------------------
// roundoff operations
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>  // function roundoff clear small elements
  roundoff(const xvector<utype>& a,utype _tol_) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) {
      c[i]=roundoff(a[i],_tol_); //CO 180409
      //if(abs(a[i])<(utype) _tol_) c[i]=a[i]=(utype) 0.0; else  c[i]=a[i];
      // c[i]=nint(a[i]/_tol_)*_tol_;
    }
    return c;
  }

  xvector<float>  _aurostd_initialize_roundoff(const xvector<float>& a) { return roundoff(a);} //CO 180409
  xvector<double> _aurostd_initialize_roundoff(const xvector<double>& a) { return roundoff(a);} //CO 180409

}

//[OBSOLETE CO 180409]namespace aurostd {  // namespace aurostd
//[OBSOLETE CO 180409]  template<class utype> xvector<utype>  // function roundoff clear small elements
//[OBSOLETE CO 180409]  roundoff(const xvector<utype>& a) {
//[OBSOLETE CO 180409]    return roundoff(a,(utype) _AUROSTD_XVECTOR_TOLERANCE_ROUNDOFF_);
//[OBSOLETE CO 180409]  }
//[OBSOLETE CO 180409]}

// ----------------------------------------------------------------------------
// GCD //CO 180409
namespace aurostd {
  int GCD(const xvector<int>& in_V){
    // find first nonzero entry
    int counter;
    bool set=false;
    for(int i=in_V.lrows;i<=in_V.urows&&!set;i++) {
      if(in_V[i]) {
        counter=i;
        set=true;
      }
    }
    if(!set) {return 1;}  // always works
    int gcd=in_V[counter];
    for(int i=counter+1;i<=in_V.urows;i++){if(in_V[i]){gcd=GCD(gcd,in_V[i]);}}// if we use chullpoint, there will be 0's!
    return gcd;
  }
  int LCM(const xvector<int>& in_V){
    // find first nonzero entry
    int counter;
    bool set=false;
    for(int i=in_V.lrows;i<=in_V.urows&&!set;i++) {
      if(in_V[i]) {
        counter=i;
        set=true;
      }
    }
    if(!set) {return 0;}  //special, trivial case: lcm must really be positive
    int lcm=in_V[counter];
    for(int i=counter+1;i<=in_V.urows;i++){if(in_V[i]){lcm=LCM(lcm,in_V[i]);}}// if we use chullpoint, there will be 0's!
    return lcm;
  }
}
 
namespace aurostd {
  template<class utype> void reduceByGCD(const xvector<utype>& in_V, xvector<utype>& out_V, utype tol){
  // DX 20191125 [OBSOLETE]   reduceByGCD(const xvector<utype>& in_V,const utype& tol){
    //DX 20191125 [OBSOLETE] xvector<utype> out_V=in_V;
    out_V=in_V;
    if(!isinteger(out_V,tol)){return;} //nothing to reduce //DX 20191125 - return type is void

    xvector<int> v1(in_V.lrows,in_V.urows); //cast to xvector of ints
    for(int i=in_V.lrows;i<=in_V.urows;i++){v1[i]=nint(in_V[i]);}
    int denom=GCD(v1);
    if(denom!=0){v1/=denom;}  //safety
    for(int i=v1.lrows;i<=v1.urows;i++){out_V[i]=(utype)v1[i];}  //cast back
    
    //DX 20191125 [OBSOLETE] return out_V;
  }
  //xvector<int> reduceByGCD(const xvector<int>& in_V){
  //  xvector<int> out_V=in_V;
  //  int denom=GCD(in_V);
  //  if(denom!=1){out_V/=denom;}
  //  return out_V;
  //}
}

// ----------------------------------------------------------------------------
// GCD //DX 20191122
// vector version (modeled after CO's xvector version)
namespace aurostd {
  int GCD(const vector<int>& in_V){
    // find first nonzero entry
    int counter;
    bool set=false;
    // REMOVE for(int i=in_V.lrows;i<=in_V.urows&&!set;i++) {
    for(uint i=0;i<in_V.size()&&!set;i++) {
      if(in_V[i]) {
        counter=i;
        set=true;
      }
    }
    if(!set) {return 1;}  // always works
    int gcd=in_V[counter];
    for(uint i=counter+1;i<in_V.size();i++){if(in_V[i]){gcd=GCD(gcd,in_V[i]);}}// if we use chullpoint, there will be 0's!
    return gcd;
  }
  int LCM(const vector<int>& in_V){
    // find first nonzero entry
    int counter;
    bool set=false;
    for(uint i=0;i<in_V.size()&&!set;i++) {
      if(in_V[i]) {
        counter=i;
        set=true;
      }
    }
    if(!set) {return 0;}  //special, trivial case: lcm must really be positive
    int lcm=in_V[counter];
    for(uint i=counter+1;i<in_V.size();i++){if(in_V[i]){lcm=LCM(lcm,in_V[i]);}}// if we use chullpoint, there will be 0's!
    return lcm;
  }
}

namespace aurostd {
  template<class utype> void reduceByGCD(const vector<utype>& in_V, vector<utype>& out_V, utype tol){
    //DX 20191125 [OBSOLETE] vector<utype> out_V=in_V;
    out_V=in_V;
    for(uint i=0;i<out_V.size();i++){ 
      if(!isinteger(out_V[i],tol)){return;}
    }

    vector<int> v1; v1.resize(in_V.size()); //cast to vector of ints
    for(uint i=0;i<in_V.size();i++){v1[i]=nint(in_V[i]);}
    int denom=GCD(v1);
    if(denom!=0){
      for(uint i=0;i<v1.size();i++){v1[i]/=denom;}
    }
    for(uint i=0;i<v1.size();i++){out_V[i]=(utype)v1[i];}  //cast back
    
    //DX 20191125 [OBSOLETE] return out_V;
  }
}

// ----------------------------------------------------------------------------
// GCD //DX 20191122
// deque version (modeled after CO's xvector version)
namespace aurostd {
  int GCD(const deque<int>& in_V){
    // find first nonzero entry
    int counter;
    bool set=false;
    // REMOVE for(int i=in_V.lrows;i<=in_V.urows&&!set;i++) {
    for(uint i=0;i<in_V.size()&&!set;i++) {
      if(in_V[i]) {
        counter=i;
        set=true;
      }
    }
    if(!set) {return 1;}  // always works
    int gcd=in_V[counter];
    for(uint i=counter+1;i<in_V.size();i++){if(in_V[i]){gcd=GCD(gcd,in_V[i]);}}// if we use chullpoint, there will be 0's!
    return gcd;
  }
  int LCM(const deque<int>& in_V){
    // find first nonzero entry
    int counter;
    bool set=false;
    for(uint i=0;i<in_V.size()&&!set;i++) {
      if(in_V[i]) {
        counter=i;
        set=true;
      }
    }
    if(!set) {return 0;}  //special, trivial case: lcm must really be positive
    int lcm=in_V[counter];
    for(uint i=counter+1;i<in_V.size();i++){if(in_V[i]){lcm=LCM(lcm,in_V[i]);}}// if we use chullpoint, there will be 0's!
    return lcm;
  }
}

namespace aurostd {
  template<class utype> void reduceByGCD(const deque<utype>& in_V, deque<utype>& out_V, utype tol){
    //DX 20191125 [OBSOLETE] deque<utype> out_V=in_V;
    out_V=in_V;
    for(uint i=0;i<out_V.size();i++){ 
      if(!isinteger(out_V[i],tol)){return;}
    }

    vector<int> v1; v1.resize(in_V.size()); //cast to vector of ints
    for(uint i=0;i<in_V.size();i++){v1[i]=nint(in_V[i]);}
    int denom=GCD(v1);
    if(denom!=0){
      for(uint i=0;i<v1.size();i++){v1[i]/=denom;}
    }
    for(uint i=0;i<v1.size();i++){out_V[i]=(utype)v1[i];}  //cast back
    
    //DX 20191125 [OBSOLETE] return out_V;
  }
}

// ----------------------------------------------------------------------------
// normalizeSumToOne //CO 180101
namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype> 
    normalizeSumToOne(const xvector<utype>& in_V,const utype& tol){
      utype s=sum(in_V);
      if(abs(s)<tol){return (utype)0.0;}
      return in_V/s;
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// swap operations
namespace aurostd {  // namespace aurostd
  template<class utype> void                                   // swap
  swap(xvector<utype>& a,const int& i,const int& j) {
    if(i<a.lrows || i>a.urows) return; // nothing to do, out of boundaries
    if(j<a.lrows || j>a.urows) return; // nothing to do, out of boundaries
    if(i==j) return; // nothing to do, no swap
    utype temp=a[i];a[i]=a[j];a[j]=temp;
  }
}

// ----------------------------------------------------------------------------
// shiftlrows operations  //CO 171128
namespace aurostd {  // namespace aurostd
  template<class utype> void  // function lrows shift lrows so first index is i
  shiftlrows(xvector<utype>& a,const int& i) {
    if(a.lrows==i){return;}
    xvector<utype> b(a.rows+i-1,i);
    int j=i;
    for(int ii=a.lrows;ii<=a.urows;ii++){b[j++]=a[ii];}
    a=b;
  }
}

// ----------------------------------------------------------------------------
// ---- Operations on complex vectors

namespace aurostd { // namespace aurostd
  template<class utype> xvector<utype>                 // function conj xvector<>
  conj (const xvector<utype>& a) {
    if (a.iscomplex) {
      xvector<utype> c(a.lrows, a.urows);
      for (int i = c.lrows; i <= c.urows; i++) {
        c[i] = conj(a[i]);
      }
      return c;
    } else {
      return a;
    }
  }
}

// ----------------------------------------------------------------------------
// ---- exponential operations on namespace template<class utype> xvector<utype>
// EXPONENTIAL OPERATIONS

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                    // function exp xvector<>
  exp(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      //      c[i]=(utype) std::exp(a[i]);
      c[i]=(utype) exp(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                    // function log xvector<>
  log(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      //      c[i]=(utype) std::log(a[i]);
      c[i]=(utype) log(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                   // function log10 xvector<>
  log10(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) std::log10(a[i]);
    return c;
  }
}

// ----------------------------------------------------------------------------
// ---- tridimensional operations on namespace template<class utype> xvector<utype>
// TRIDIMENSIONAL OPERATIONS
namespace aurostd {  // namespace aurostd
  template<class utype> utype
  distance(const xvector<utype>& v1,const xvector<utype>& v2) {
    return (utype) modulus(v1-v2);
  }
}

// ----------------------------------------------------------------------------
// ---- trigonometric operations on namespace template<class utype> xvector<utype>
// TRIGONOMETRIC OPERATIONS

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                    // function sin xvector<>
  sin(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++) {
      // if(aurostd::_isreal(c[i])) c[i]=(utype) std::sin(aurostd::_real(a[i]));
      // else c[i]=(utype) aurostd::sin(a[i]);
      c[i]=(utype) sin(a[i]);
    }
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                    // function cos xvector<>
  cos(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      //    c[i]=(utype) std::cos(a[i]);
      c[i]=(utype) cos(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                   // function asin xvector<>
  asin(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) std::asin(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                   // function acos xvector<>
  acos(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) std::acos(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                    // function tan xvector<>
  tan(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) std::tan(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                   // function atan xvector<>
  atan(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) std::atan(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                   // function sinh xvector<>
  sinh(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      //     c[i]=(utype) std::sinh(a[i]);
      c[i]=(utype) sinh(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                  // function asinh xvector<>
  asinh(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) asinh(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                   // function cosh xvector<>
  cosh(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      //     c[i]=(utype) std::cosh(a[i]);
      c[i]=(utype) cosh(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                  // function acosh xvector<>
  acosh(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) acosh(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                   // function tanh xvector<>
  tanh(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) std::tanh(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                  // function atanh xvector<>
  atanh(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) atanh(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                 // function cotanh xvector<>
  cotanh(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) cotanh(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                // function acotanh xvector<>
  acotanh(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) acotanh(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                    // function sec xvector<>
  sec(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) sec(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                  // function cosec xvector<>
  cosec(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) cosec(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                   // function sech xvector<>
  sech(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) sech(a[i]);
    return c;
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>                 // function cosech xvector<>
  cosech(const xvector<utype>& a) {
    xvector<utype> c(a.lrows,a.urows);
    for(int i=c.lrows;i<=c.urows;i++)
      c[i]=(utype) cosech(a[i]);
    return c;
  }
}

// ----------------------------------------------------------------------------
// ---------------------- trigonometric operations between two three xvectors<>
// TRIGONOMETRIC OPERATIONS BETWEEN TWO/THREE VECTORS

namespace aurostd {  // namespace aurostd
  template<class utype> double   // cos of angle between two vectors
  cos(const xvector<utype>& v1,const xvector<utype>& v2) {
    if(v1.rows!=v2.rows) {
      cerr << _AUROSTD_XLIBS_ERROR_ << " XVECTOR.CPP: cos(xvector,xvector) v1.rows,v2.rows=" << v1.rows << "," << v2.rows << endl;
      exit(0);
    }
    double out=0.0,n_v1=0.0,n_v2=0.0;
    int size=v1.rows,i;
    xvector<double> _v1(size),_v2(size);
    for(i=0;i<size;i++) {
      _v1(i+_v1.lrows)=(double) v1(i+v1.lrows);
      _v2(i+_v2.lrows)=(double) v2(i+v2.lrows);
    }
    out=0.0;
    for(i=1;i<=size;i++)
      out+=_v1[i]*_v2[i]; // scalar product
    n_v1=modulus(_v1);
    n_v2=modulus(_v2);
    // cerr << n_v1 << " " << n_v2 << endl;
    if(n_v1==0.0) {
      cerr << _AUROSTD_XLIBS_ERROR_ << " XVECTOR.CPP: cos(xvector,xvector)=modulus(v1)=0" << endl;
      exit(0);
    } 
    if(n_v2==0.0) {
      cerr << _AUROSTD_XLIBS_ERROR_ << " XVECTOR.CPP: cos(xvector,xvector)=modulus(v2)=0" << endl;
      exit(0);
    }
    // assert(n_v1>0 && n_v2>0);
    out/=(n_v1*n_v2);
    if(out >=1.0) out =1.0; //make sure numerical errors don't place the arguments outside of the [-1,1] interval
    if(out<=-1.0) out=-1.0; //make sure numerical errors don't place the arguments outside of the [-1,1] interval
    return out;  
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> double   // cos of angle between two vectors
  getcos(const xvector<utype>& v1,const xvector<utype>& v2) {            // for convasp
    return cos(v1,v2);
  }
}

namespace aurostd {
  template<class utype> double   // sin of angle between two vectors
  // namespace aurostd
  sin(const xvector<utype>& v1,const xvector<utype>& v2) {
    double _cos=0.0;
    _cos=cos(v1,v2);
    return (double) std::sqrt(1.0-_cos*_cos);
  }  
}

namespace aurostd {  // namespace aurostd
  template<class utype> double   // sin of angle between two vectors
  getsin(const xvector<utype>& v1,const xvector<utype>& v2) {      // for convasp
    return sin(v1,v2);
  }
}

namespace aurostd {
  template<class utype> double   // angle between two vectors in radiants !!!
  // namespace aurostd
  angle(const xvector<utype>& v1,const xvector<utype>& v2) {
    return (double) std::acos(cos(v1,v2));
  }  
}

namespace aurostd {
  template<class utype> double   // angle between two vectors in radiants !!!
  // namespace aurostd
  getangle(const xvector<utype>& v1,const xvector<utype>& v2) {     // for convasp
    return angle(v1,v2);
  }  
}

namespace aurostd {
  template<class utype> double   // angle between three vectors in radiants !!!
  // namespace aurostd
  angle(const xvector<utype>& v0,const xvector<utype>& v1,const xvector<utype>& v2) {
    return (double) std::acos(cos(v1-v0,v2-v0));
  }  
}

namespace aurostd {
  template<class utype> double   // angle between three vectors in radiants !!!
  // namespace aurostd
  getangle(const xvector<utype>& v0,const xvector<utype>& v1,const xvector<utype>& v2) {     // for convasp
    return angle(v1-v0,v2-v0);
  }  
}

namespace aurostd {
  template<class utype> bool  //determine if two vectors are collinear //CO 180409
    isCollinear(const xvector<utype>& v0,const xvector<utype>& v1,const utype& tol) {
      return abs(aurostd::angle(v0,v1))<tol;
    }
}

namespace aurostd {
  template<class utype> xvector<utype> //get centroid of data points //CO 180409
    getCentroid(const vector<xvector<utype> >& points) {
    xvector<utype> centroid;
    if(points.size()==0){return centroid;}
    centroid=points[0];
    for(uint i=1;i<points.size();i++){centroid+=points[i];}
    centroid/=points.size();
    return centroid;
    }
}

// ----------------------------------------------------------------------------
// ---------------------- trigonometric operations between two GENERAL xvectors<>
// TRIGONOMETRIC OPERATIONS BETWEEN TWO GENERAL VECTORS

//CO 180409
//given a vector, will give angles from x,y,z... axes (IN THAT ORDER)
//to get it in familar aximuthal and polar angles, need to swap indices
//I show an example below
namespace aurostd {
  template<class utype> xvector<double> 
    getGeneralAngles(const xvector<utype>& vec,const utype& tol) { //CO 180409
    //https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates
    xvector<double> angles(vec.urows-1,vec.lrows);
    for(int i=vec.lrows;i<vec.urows;i++){angles[i]=getGeneralAngle(vec,i,tol);}
    return angles;
  }

  template<class utype> double 
    getGeneralAngle(const xvector<utype>& _vec,int _i,const utype& tol) { //CO 180409
    //https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates

    //CO prefers lrows==0
    //force vec to have lrows==0, so code always works (robust)
    //don't worry, we shift the solution back later
    xvector<utype> vec=_vec;
    aurostd::shiftlrows(vec,0);
    int i=(_i-_vec.lrows);

    //this is a nice generalizable formulation
    //but it's a rotation from common azimuthal and polar angles we know in 3D and up
    //the angle i corresponds to angle off axis i, so x is really z for common 3D representation
    //so swap coordinates, easily seen in three coordinates:
    //x,y,z must become z,x,y
    //x,y,z (start)
    //z,y,x (swap 1,3)
    //z,x,y (swap 2,3)
    if(0){
      if(vec.rows>2){
        int j=0;
        for(int k=0;k<vec.rows-2;k++){aurostd::swap(vec,j++,vec.rows-1);}
        cerr << " " << vec << " " << endl; 
      }
    }
    
    //first check special case 1: that vec[i]!=0 and vec[k]==0 for all k>i
    bool special_case_1=abs(vec[i])>=tol;
    for(int k=i+1;k<vec.rows&&special_case_1;k++){special_case_1=(special_case_1 && abs(vec[k])<tol);}
    if(special_case_1){
      if(std::signbit(vec[i])){return pi;}
      else{return 0.0;}
    }

    double angle=0.0,denominator=0.0;
    for(int j=i;j<vec.rows;j++){denominator+=(vec[j]*vec[j]);}
    denominator=sqrt(denominator);
    if(abs(denominator)<tol){return angle;}  //special case 2: keep at 0
    angle=std::acos(vec[i]/denominator);
    if(i==vec.rows-2 && std::signbit(vec[i+1])){angle=2.0*pi-angle;}  //be wary, i<vec.rows
    return angle;
  }
}

namespace aurostd {
  template<class utype> xvector<utype>
    getGeneralNormal(const vector<xvector<utype> >& _directive_vectors){ //CO 180409
      bool LDEBUG=(FALSE || XHOST.DEBUG);
      string soliloquy="pflow::getGeneralNormal():";
      
      //tests of stupidity
      xvector<utype> dummy;
      if(!_directive_vectors.size()){return dummy;}
      int dim=_directive_vectors[0].rows;
      for(int i=1;i<(int)_directive_vectors.size();i++){if(_directive_vectors[i].rows!=dim){return dummy;}}
      int lrows=_directive_vectors[0].lrows;  //save for later
      for(int i=1;i<(int)_directive_vectors.size();i++){if(_directive_vectors[i].lrows!=lrows){return dummy;}}
      //in general, this is NOT needed, but it's safe to ensure exact normal (not ambiguous in direction)
      if(dim-1!=(int)_directive_vectors.size()){return dummy;}
      
      //CO prefers lrows==0
      //force directive vectors to have lrows==0, so code always works (robust)
      //don't worry, we shift the solution back later
      xvector<utype> directive_vector;
      vector<xvector<utype> > directive_vectors;
      for(int i=0;i<(int)_directive_vectors.size();i++){
        directive_vector=_directive_vectors[i];
        aurostd::shiftlrows(directive_vector,0);
        directive_vectors.push_back(directive_vector);
      }

      //normal is calculated by a generalized cross product, i.e., cofactor expansions
      //we simply generalize the method for cross-product in 3-space
      //take matrix of 2 vectors in last two rows of matrix
      // [ i j k ]
      // [ 1 2 3 ]
      // [ 4 5 6 ]
      //next we take determinant of submatrix (minor submatrix=minordet) for that dimension
      //submatrix of i:  det(2*6-3*5) forms coefficient for first dimension, etc.
      //don't forget to do alternating negative sign!
      //there's a few sources for this:
      //W. S. Massey, "Cross Products of Vectors in Higher Dimensional Euclidean Spaces" The American Mathematical Monthly, Vol. 90, No. 10 (Dec., 1983), pp. 697-701 (http://www.jstor.org/stable/2323537)
      //https://en.wikipedia.org/wiki/Cross_product#Multilinear_algebra
      //https://ef.gy/linear-algebra:normal-vectors-in-higher-dimensional-spaces
      //"it can be defined in a coordinate independent way as the Hodge dual of the wedge product of the arguments"
      xvector<utype> normal(dim-1,0);
      xmatrix<utype> mat(dim,dim-1,1,1); //must start at 1 to work with det(), minordet()
      for(int i=0;i<dim;i++){
        for(int j=0;j<dim-1;j++){
          mat(i+1,j+1)=directive_vectors[j][i];
        }
      }
      //mat is actually transpose of above
      // [ 1 4 i ]
      // [ 2 5 j ]
      // [ 3 6 k ]
      if(LDEBUG){cerr << soliloquy << " cross-product matrix:" << endl << mat << endl;}

      //since the i,j,k column is phony, we simply have to knock out the corresponding row
      //i.e., row-restricted cofactor
      //therefore, the note in wikipedia that says we need to take minordet(mat,i+1,dim) vs. minordet(mat,i+1,0)
      //is useless, they are the same for even/odd n, I checked.
      //so we keep as implemented
      //this is also verified in http://www.jstor.org/stable/2323537
      for(int i=0;i<dim;i++){normal[i]=std::pow(-1,i)*minordet(mat,i+1,0);}
      normal/=modulus(normal);    //normalize

      aurostd::shiftlrows(normal,lrows); //shift back to original!
      return normal;
    }
}

namespace aurostd {
  template<class utype> xvector<utype>
    pointLineIntersection(const xvector<utype>& a,const xvector<utype>& n,const xvector<utype>& p){ //CO 180520
      //https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line (vector formulation)
      //equation of line: x=a+t*n
      //a is point on line
      //n is direction of line
      //p is arbitrary point
      return p + ((a-p) - (aurostd::scalar_product((a-p),n)*n));
    }
  template<class utype> bool
    linePlaneIntersect(const xvector<utype>& p0,const xvector<utype>& n,const xvector<utype>& l0, const xvector<utype>& l,double& d,xvector<utype>& intersection){ //CO 180520
      bool LDEBUG=(FALSE || XHOST.DEBUG);
      string soliloquy="aurostd::linePlaneIntersect():";
      //https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection (algebraic form)
      //plane is defined as (p-p0)*n=0 for set of points p on plane
      //equation of line: p=d*l+l0
      double numerator=aurostd::scalar_product((p0-l0),n);
      double denominator=aurostd::scalar_product(l,n);
      d=0.0;
      if(LDEBUG){
        cerr << soliloquy << " p0=" << p0 << endl;
        cerr << soliloquy << " n=" << n << endl;
        cerr << soliloquy << " l0=" << l0 << endl;
        cerr << soliloquy << " l=" << l << endl;
        cerr << soliloquy << " numerator=" << numerator << endl;
        cerr << soliloquy << " denominator=" << denominator << endl;
      }
      if(aurostd::isequal(denominator,0.0,_ZERO_TOL_)){ //line and plane are parallel
        if(aurostd::isequal(numerator,0.0,_ZERO_TOL_)){intersection=l0;return true;} //line is contained in plane, so return back initial point
        else{return false;} //line and plane have no intersection
      }
      d=numerator/denominator;
      intersection=d*l+l0;
      if(LDEBUG){
        cerr << soliloquy << " d=" << d << endl;
        cerr << soliloquy << " intersection=" << intersection << endl;
      }
      return true;
    }
}

// ----------------------------------------------------------------------------
// ------------------------------------------------------- simple sort routines

namespace aurostd {
  template<class utype> xvector<utype>
  sort(const xvector<utype>& a) {  // function shellsort xvector<utype>
    return shellsort(a);  
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>
  shellsort(const xvector<utype>& a) {  // function shellsort xvector<utype>
    // returns a sorted vector !!
    // the input is not messed
    xvector<utype> c(a);
    // N^1.25 for random , N^1.5 for worst
    long i,j,inc=1,n=c.rows,as=c.lrows-1;
    utype v;
    do {
      inc*=3;
      inc++;
    } while (inc<=n);
    do {
      inc/=3;
      for(i=inc+1;i<=n;i++) {
	v=c[i+as];
	j=i;
	while (c[j-inc+as]>v) {
	  c[j+as]=c[j-inc+as];
	  j-=inc;
	  if(j<=inc) break;
	}
	c[j+as]=v;
      }
    } while (inc>1);
    return c;  
  }
}

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>
  heapsort(const xvector<utype>& a) {  // function shellsort xvector<utype>
    xvector<utype> c(a);
    // Nlog2N-ALWAYS
    long i,ir,j,l,n=c.rows,as=c.lrows-1;
    utype ra;
    if(n<2) return c;
    l=(n>>1)+1;
    ir=n;
    for(;;) {
      if(l>1) {
	ra=c[(--l)+as];
      } else {
	ra=c[ir+as];
	c[ir+as]=c[1+as];
	if(--ir==1) {
	  c[1+as]=ra;
	  break;
	}
      }
      i=l;
      j=l+l;
      while (j<=ir) {
	if(j<ir && c[j+as]<c[j+1+as]) j++;
	if(ra<c[j+as]) {
	  c[i+as]=c[j+as];
	  i=j;
	  j<<=1;
	} else j=ir+1;
      }
      c[i+as]=ra;
    }
    return c;  
  }
}

// ----------------------------------------------------------------------------
// ------------------------ SORT ARR following ARR order (from lower to higher)

#define _XQSORT_SWAP(a,b)        {temp=(a);(a)=(b);(b)=temp;}
#define _XQSORT_M 8              // below M qsort does shellsort
#define _XQSORT_NSTACK 100       // has to be 2log2(N) where N is the max

namespace aurostd {  // namespace aurostd
  template<class utype> xvector<utype>
  quicksort(const xvector<utype>& arr) {
    // Nlog2N for average, N^2 for worst
    long int i,j,k,l=1,ir=arr.rows,jstack=0,as=arr.lrows-1;
    xvector<int> istack(1,_XQSORT_NSTACK);
    utype a,temp;
    
    for(;;) {
      if(ir-l<_XQSORT_M) {
	for(j=l+1;j<=ir;j++) {
	  a=arr[j+as];
	  for(i=j-1;i>=1;i--) {
	    if(arr[i+as]<=a) break;
	    arr[i+1+as]=arr[i+as];
	  }
	  arr[i+1+as]=a;
	}
	if(jstack == 0) break;
	ir=istack[jstack--];
	l=istack[jstack--];
      } else {
	k=(l+ir)>>1;
	_XQSORT_SWAP(arr[k+as],arr[l+1+as]);
	if(arr[l+1+as] > arr[ir+as])
	  _XQSORT_SWAP(arr[l+1+as],arr[ir+as]);    
	if(arr[l+as] > arr[ir+as])
	  _XQSORT_SWAP(arr[l+as],arr[ir+as]);
	if(arr[l+1+as] > arr[l+as])
	  _XQSORT_SWAP(arr[l+1+as],arr[l+as]);
	i=l+1;
	j=ir;
	a=arr[l+as];
	for(;;) {
	  do i++; while (arr[i+as]<a);
	  do j--; while (arr[j+as]>a);
	  if(j<i) break;
	  _XQSORT_SWAP(arr[i+as],arr[j+as]);
	}
	arr[l+as]=arr[j+as];
	arr[j+as]=a;
	jstack += 2;
	if(jstack>_XQSORT_NSTACK) {
	  cerr << _AUROSTD_XLIBS_ERROR_ << " _XQSORT_NSTACK too small in sort." << endl;
	  exit(0);
	}
	if(ir-i+1 >= j-l) {
	  istack[jstack]=ir;
	  istack[jstack-1]=i;
	  ir=j-1;
	} else {
	  istack[jstack]=j-1;
	  istack[jstack-1]=l;
	  l=i;
	}
      }
    }
    return arr;
  }
}
#undef _XQSORT_M
#undef _XQSORT_NSTACK
#undef _XQSORT_SWAP

// ----------------------------------------------------------------------------
// ---------------- SORT ARR and BRR following ARR order (from lower to higher)

#define _XSORT_SWAPT(a,b,temp) {temp=(a);(a)=(b);(b)=temp;}
#define _XSORT_M 7
#define _XSORT_NSTACK 50

namespace aurostd {  // namespace aurostd
  template<class utype1,class utype2>                                 // function quicksort
  void quicksort2(unsigned long n, xvector<utype1>& arr, xvector<utype2>&  brr) {
    unsigned long i,ir=n,j,k,l=1;
    int jstack=0;
    utype1 a,atemp;
    utype2 b,btemp;
    
    xvector<int> istack(1,_XSORT_NSTACK);
    for(;;) {
      if(ir-l<_XSORT_M) {
	for(j=l+1;j<=ir;j++) {
	  a=arr[j];
	  b=brr[j];
	  for(i=j-1;i>=1;i--) {
	    if(arr[i] <=a) break;
	    arr[i+1]=arr[i];
	    brr[i+1]=brr[i];
	  }
	  arr[i+1]=a;
	  brr[i+1]=b;
	}
	if(!jstack) {
	  return;
	}
	ir=istack[jstack];
	l=istack[jstack-1];
	jstack -=2;
      } else {
	k=(l+ir)>> 1;
	_XSORT_SWAPT(arr[k],arr[l+1],atemp);
	_XSORT_SWAPT(brr[k],brr[l+1],btemp);
	if(arr[l+1]>arr[ir]) {
	  _XSORT_SWAPT(arr[l+1],arr[ir],atemp);
	  _XSORT_SWAPT(brr[l+1],brr[ir],btemp);
	}
	if(arr[l]>arr[ir]) {
	  _XSORT_SWAPT(arr[l],arr[ir],atemp);
	  _XSORT_SWAPT(brr[l],brr[ir],btemp);
	}
	if(arr[l+1]>arr[l]) {
	  _XSORT_SWAPT(arr[l+1],arr[l],atemp);
	  _XSORT_SWAPT(brr[l+1],brr[l],btemp);
	}
	i=l+1;
	j=ir;
	a=arr[l];
	b=brr[l];
	for(;;) {
	  do i++; while (arr[i]<a);
	  do j--; while (arr[j]>a);
	  if(j<i) break;
	  _XSORT_SWAPT(arr[i],arr[j],atemp);
	  _XSORT_SWAPT(brr[i],brr[j],btemp);
	}
	arr[l]=arr[j];arr[j]=a;
	brr[l]=brr[j];brr[j]=b;
	jstack +=2;
	if(jstack>_XSORT_NSTACK) {
	  cerr << _AUROSTD_XLIBS_ERROR_ << " _XSORT_NSTACK too small in sort2." << endl;
	  exit(0);
	}
	if(ir-i+1>=j-l) {
	  istack[jstack]=ir;
	  istack[jstack-1]=i;
	  ir=j-1;
	} else {
	  istack[jstack]=j-1;
	  istack[jstack-1]=l;
	  l=i;
	}
      }
    }
  }
}

// ----------------------------------------------------------------------------
// ---------------- SORT ARR,BRR,CRR following ARR order (from lower to higher)

namespace aurostd {  // namespace aurostd
  template<class utype1, class utype2, class utype3>                                 // function quicksort
  void quicksort3(unsigned long n, xvector<utype1>& arr, xvector<utype2>&  brr, xvector<utype3>&  crr) {
    unsigned long i,ir=n,j,k,l=1;
    int jstack=0;
    utype1 a,atemp;
    utype2 b,btemp;
    utype3 c,ctemp;
    
    xvector<int> istack(1,_XSORT_NSTACK);
    for(;;) {
      if(ir-l<_XSORT_M) {
	for(j=l+1;j<=ir;j++) {
	  a=arr[j];
	  b=brr[j];
	  c=crr[j];
	  for(i=j-1;i>=1;i--) {
	    if(arr[i] <=a) break;
	    arr[i+1]=arr[i];
	    brr[i+1]=brr[i];
	    crr[i+1]=crr[i];
	  }
	  arr[i+1]=a;
	  brr[i+1]=b;
	  crr[i+1]=c;
	}
	if(!jstack) {
	  return;
	}
	ir=istack[jstack];
	l=istack[jstack-1];
	jstack -=2;
      } else {
	k=(l+ir)>> 1;
	_XSORT_SWAPT(arr[k],arr[l+1],atemp);
	_XSORT_SWAPT(brr[k],brr[l+1],btemp);
	_XSORT_SWAPT(crr[k],crr[l+1],ctemp);
	if(arr[l+1]>arr[ir]) {
	  _XSORT_SWAPT(arr[l+1],arr[ir],atemp);
	  _XSORT_SWAPT(brr[l+1],brr[ir],btemp);
	  _XSORT_SWAPT(crr[l+1],crr[ir],ctemp);
	}
	if(arr[l]>arr[ir]) {
	  _XSORT_SWAPT(arr[l],arr[ir],atemp);
	  _XSORT_SWAPT(brr[l],brr[ir],btemp);
	  _XSORT_SWAPT(crr[l],crr[ir],ctemp);
	}
	if(arr[l+1]>arr[l]) {
	  _XSORT_SWAPT(arr[l+1],arr[l],atemp);
	  _XSORT_SWAPT(brr[l+1],brr[l],btemp);
	  _XSORT_SWAPT(crr[l+1],crr[l],ctemp);
	}
	i=l+1;
	j=ir;
	a=arr[l];
	b=brr[l];
	c=crr[l];
	for(;;) {
	  do i++; while (arr[i]<a);
	  do j--; while (arr[j]>a);
	  if(j<i) break;
	  _XSORT_SWAPT(arr[i],arr[j],atemp);
	  _XSORT_SWAPT(brr[i],brr[j],btemp);
	  _XSORT_SWAPT(crr[i],crr[j],ctemp);
	}
	arr[l]=arr[j];arr[j]=a;
	brr[l]=brr[j];brr[j]=b;
	crr[l]=crr[j];crr[j]=c;
	jstack +=2;
	if(jstack>_XSORT_NSTACK) {
	  cerr << _AUROSTD_XLIBS_ERROR_ << " _XSORT_NSTACK too small in sort2." << endl;
	  exit(0);
	}
	if(ir-i+1>=j-l) {
	  istack[jstack]=ir;
	  istack[jstack-1]=i;
	  ir=j-1;
	} else {
	  istack[jstack]=j-1;
	  istack[jstack-1]=l;
	  l=i;
	}
      }
    }
  }
}

// ----------------------------------------------------------------------------
// ------------ SORT ARR,BRR,CRR,DRR following ARR order (from lower to higher)

namespace aurostd {  // namespace aurostd
  template<class utype1, class utype2, class utype3, class utype4>                // function quicksort
  void quicksort4(unsigned long n, xvector<utype1>& arr, xvector<utype2>&  brr, xvector<utype3>&  crr, xvector<utype4>&  drr) {
    unsigned long i,ir=n,j,k,l=1;
    int jstack=0;
    utype1 a,atemp;
    utype2 b,btemp;
    utype3 c,ctemp;
    utype4 d,dtemp;
    
    xvector<int> istack(1,_XSORT_NSTACK);
    for(;;) {
      if(ir-l<_XSORT_M) {
	for(j=l+1;j<=ir;j++) {
	  a=arr[j];
	  b=brr[j];
	  c=crr[j];
	  d=drr[j];
	  for(i=j-1;i>=1;i--) {
	    if(arr[i] <=a) break;
	    arr[i+1]=arr[i];
	    brr[i+1]=brr[i];
	    crr[i+1]=crr[i];
	    drr[i+1]=drr[i];
	  }
	  arr[i+1]=a;
	  brr[i+1]=b;
	  crr[i+1]=c;
	  drr[i+1]=d;
	}
	if(!jstack) {
	  return;
	}
	ir=istack[jstack];
	l=istack[jstack-1];
	jstack -=2;
      } else {
	k=(l+ir)>> 1;
	_XSORT_SWAPT(arr[k],arr[l+1],atemp);
	_XSORT_SWAPT(brr[k],brr[l+1],btemp);
	_XSORT_SWAPT(crr[k],crr[l+1],ctemp);
	_XSORT_SWAPT(drr[k],drr[l+1],dtemp);
	if(arr[l+1]>arr[ir]) {
	  _XSORT_SWAPT(arr[l+1],arr[ir],atemp);
	  _XSORT_SWAPT(brr[l+1],brr[ir],btemp);
	  _XSORT_SWAPT(crr[l+1],crr[ir],ctemp);
	  _XSORT_SWAPT(drr[l+1],drr[ir],dtemp);
	}
	if(arr[l]>arr[ir]) {
	  _XSORT_SWAPT(arr[l],arr[ir],atemp);
	  _XSORT_SWAPT(brr[l],brr[ir],btemp);
	  _XSORT_SWAPT(crr[l],crr[ir],ctemp);
	  _XSORT_SWAPT(drr[l],drr[ir],dtemp);
	}
	if(arr[l+1]>arr[l]) {
	  _XSORT_SWAPT(arr[l+1],arr[l],atemp);
	  _XSORT_SWAPT(brr[l+1],brr[l],btemp);
	  _XSORT_SWAPT(crr[l+1],crr[l],ctemp);
	  _XSORT_SWAPT(drr[l+1],drr[l],dtemp);
	}
	i=l+1;
	j=ir;
	a=arr[l];
	b=brr[l];
	c=crr[l];
	d=drr[l];
	for(;;) {
	  do i++; while (arr[i]<a);
	  do j--; while (arr[j]>a);
	  if(j<i) break;
	  _XSORT_SWAPT(arr[i],arr[j],atemp);
	  _XSORT_SWAPT(brr[i],brr[j],btemp);
	  _XSORT_SWAPT(crr[i],crr[j],ctemp);
	  _XSORT_SWAPT(drr[i],drr[j],dtemp);
	}
	arr[l]=arr[j];arr[j]=a;
	brr[l]=brr[j];brr[j]=b;
	crr[l]=crr[j];crr[j]=c;
	drr[l]=drr[j];drr[j]=d;
	jstack +=2;
	if(jstack>_XSORT_NSTACK) {
	  cerr << _AUROSTD_XLIBS_ERROR_ << " _XSORT_NSTACK too small in sort2." << endl;
	  exit(0);
	}
	if(ir-i+1>=j-l) {
	  istack[jstack]=ir;
	  istack[jstack-1]=i;
	  ir=j-1;
	} else {
	  istack[jstack]=j-1;
	  istack[jstack-1]=l;
	  l=i;
	}
      }
    }
  }
}

#undef _XSORT_M
#undef _XSORT_NSTACK
#undef _XSORT_SWAPT

//CO190629 - START
namespace aurostd {  // namespace aurostd
  //sort by a particular index ONLY
  //using bool ascending MAY become slow for long containers, as it needs to check the bool with every comparison
  //split out into separate classes (compareVecElementAscending() vs. compareVecElementDescending()) for those cases
  //I don't believe this is a practical problem though
  template<class utype> 
    compareVecElement<utype>::compareVecElement(uint ind,bool ascending) : m_uindex_sort((uint)ind),m_iindex_sort((int)ind),m_ascending_sort(ascending) {} //CO190629
  template<class utype> 
    bool compareVecElement<utype>::operator() (const vector<utype>& a,const vector<utype>& b) { //CO190629
      if(a.size()!=b.size()){throw aurostd::xerror(_AFLOW_FILE_NAME_,"compareVecElement::operator()():","a.size()!=b.size()",_INDEX_MISMATCH_);}
      if(m_uindex_sort>=a.size()){throw aurostd::xerror(_AFLOW_FILE_NAME_,"compareVecElement::operator()():","index_sort>=a.size()",_INDEX_BOUNDS_);}
      if(m_ascending_sort){return a[m_uindex_sort]<b[m_uindex_sort];}
      return a[m_uindex_sort]>b[m_uindex_sort]; //descending sort
  }
  template<class utype> 
    bool compareVecElement<utype>::operator() (const xvector<utype>& a,const xvector<utype>& b) { //CO190629
      if(a.lrows!=b.lrows){throw aurostd::xerror(_AFLOW_FILE_NAME_,"compareVecElement::operator()():","a.lrows!=b.lrows",_INDEX_MISMATCH_);}
      if(a.rows!=b.rows){throw aurostd::xerror(_AFLOW_FILE_NAME_,"compareVecElement::operator()():","a.rows!=b.rows",_INDEX_MISMATCH_);}
      if(m_iindex_sort<a.lrows||m_iindex_sort>a.urows){throw aurostd::xerror(_AFLOW_FILE_NAME_,"compareVecElement::operator()():","index_sort<a.lrows||index_sort>a.urows",_INDEX_BOUNDS_);}
      if(m_ascending_sort){return a[m_iindex_sort]<b[m_iindex_sort];}
      return a[m_iindex_sort]>b[m_iindex_sort]; //descending sort
  }
  //sort by all indices in increasing order
  template<class utype>
    bool compareVecElements(const vector<utype>& a,const vector<utype>& b) { //CO190629
      if(a.size()!=b.size()){throw aurostd::xerror(_AFLOW_FILE_NAME_,"compareVecElements():","a.size()!=b.size()",_INDEX_MISMATCH_);}
      for(uint i=0;i<a.size();i++){if(a[i]!=b[i]){return a[i]<b[i];}}
      return false;
    }
  template<class utype>
    bool compareXVecElements(const aurostd::xvector<utype>& a,const aurostd::xvector<utype>& b) { //CO190629
      if(a.lrows!=b.lrows){throw aurostd::xerror(_AFLOW_FILE_NAME_,"compareXVecElements():","a.lrows!=b.lrows",_INDEX_MISMATCH_);}
      if(a.rows!=b.rows){throw aurostd::xerror(_AFLOW_FILE_NAME_,"compareXVecElements():","a.rows!=b.rows",_INDEX_MISMATCH_);}
      for(int i=a.lrows;i<=a.urows;i++){if(a[i]!=b[i]){return a[i]<b[i];}}
      return false;
    }
} // namespace aurostd
//CO190629 - STOP

// ----------------------------------------------------------------------------
// ----------------------------------------- STATS stuff

namespace aurostd {
  template<class utype> utype mean(const xvector<utype>& a){return sum(a)/a.rows;} //CO190520
  template<class utype> utype stddev(const xvector<utype>& a){ //CO190520
    utype avg=mean(a);
    utype sd=(utype)0,diff=(utype)0;
    for(int i=a.lrows;i<=a.urows;i++){diff=(a[i]-avg);sd+=diff*diff;}
    sd/=a.rows-1;
    sd=sqrt(sd);
    return sd;
  }
}

// ----------------------------------------------------------------------------
// quartiles - CO 171202
namespace aurostd {
  template<class utype> void 
    getQuartiles(const xvector<utype>& _a,utype& q1,utype& q2,utype& q3){ //CO 180409
    q1=q2=q3=(utype)AUROSTD_NAN;
    if(_a.rows<4){return;} //not enough points to do statistics (need at least 3 quartile)
    xvector<utype> a=sort(_a);  //unfortunate that we have to make a full copy here, but alas, we will
    shiftlrows(a,0);    //CO 180314 - even/odd specifications starting at 0
    //get first, second (median), and third quartiles
    int i1=a.rows/4+a.lrows;
    int i2=a.rows/2+a.lrows;
    int i3=i1+i2;
    if(a.rows%2==0){ //even, the harder of the two
      q1=(a[i1-1]+a[i1])/2.0;
      q2=(a[i2-1]+a[i2])/2.0; //not needed
      q3=(a[i3-1]+a[i3])/2.0;
    }else{  //odd, easy
      q1=a[i1];
      q2=a[i2];  //not needed
      q3=a[i3];
    }
  }
}

namespace aurostd {
  template<class utype> utype
    getMAD(const xvector<utype>& _a,utype median){  //absolute deviation around the median (MAD) //CO 180409
    //an EXCELLENT measure of spread of the data, robust to finding outliers
    //breakpoint = 50%!
    //see doi=10.1016/j.jesp.2013.03.013
    //ADDITIONAL NOTES:
    //"Absolute deviation from the median was (re-)discovered and 
    //popularized by Hampel (1974) who attributes the idea to Carl Friedrich Gauss (1777–1855). 
    //The median (M) is, like the mean, a measure of central tendency 
    //but offers the advantage of being very insensitive to the presence of outliers. 
    //One indicator of this insensitivity is the “breakdown point” (see, e.g., Donoho & Huber, 1983). 
    //The estimator's breakdown point is the maximum proportion of 
    //observations that can be contaminated (i.e., set to infinity) without forcing the estimator to 
    //result in a false value (infinite or null in the case of an estimator of scale). 
    //For example, when a single observation has an infinite value, the mean of all observations 
    //becomes infinite; hence the mean's breakdown point is 0. By contrast, the median value remains unchanged. 
    //The median becomes absurd only when more than 50% of the observations are infinite. 
    //With a breakdown point of 0.5, the median is the location estimator that has the highest 
    //breakdown point. Exactly the same can be said about the Median Absolute Deviation as an 
    //estimator of scale (see the formula below for a definition). 
    //Moreover, the MAD is totally immune to the sample size. 
    //These two properties have led Huber (1981) to describe the MAD as the 
    //“single most useful ancillary estimate of scale” (p. 107). 
    //It is for example more robust than the classical interquartile range (see Rousseeuw & Croux, 1993),
    //which has a breakdown point of 25% only.

    //if no median provided, find it
    if(median==(utype)AUROSTD_NAN){
      utype q1,q3;
      getQuartiles(_a,q1,median,q3);
    }
    xvector<utype> diff(_a); 
    for(int i=diff.lrows;i<=diff.urows;i++){diff[i]=abs(diff[i]-median);}
    utype q1,q2,q3;
    getQuartiles(diff,q1,q2,q3);
    return (utype)1.4826*q2;  //b = 1/Q(0.75), where Q(0.75) is the 0.75 quantile of that underlying distribution, this assumes NORMAL
  }
}
// ----------------------------------------------------------------------------
// ----------------------------------------- implementation for extra data type
#define DEBUG_CONVOLUTION 0
namespace aurostd {
  //CO190419 - convolution and moving average
  //see 'doc conv' in matlab
  //also see numerical recipes in C 2nd edition, page 538
  template<class utype> xvector<utype> convolution(const xvector<utype>& signal_input,const xvector<utype>& response_input,int SHAPE) {
    vector<uint> sum_counts;
    return convolution(signal_input,response_input,sum_counts,SHAPE);
  }
  template<class utype> xvector<utype> convolution(const xvector<utype>& signal_input,const xvector<utype>& response_input,vector<uint>& sum_counts,int SHAPE) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="aurostd::convolution():";
    if(LDEBUG){
      cerr << soliloquy << " signal_input=" << signal_input << endl;
      cerr << soliloquy << " response_input=" << response_input << endl;
    }
    if(signal_input.lrows!=response_input.lrows){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"signal_input.lrows!=response_input.lrows",_INDEX_MISMATCH_);}
    int lrows=signal_input.lrows; //fixed
    int size=signal_input.rows+response_input.rows-1;
    vector<uint> sum_counts_full((uint)size,0);
    xvector<utype> conv(size+(lrows-1),lrows);conv.reset();  //set it all to 0 //CO191110
    vector<int> ind_zero_padding; //keep indices that require zero-padding for 'valid'
    int ind=0;
    bool k_added=false;
    for(int k=conv.lrows;k<=conv.urows;k++){
      k_added=false;
      for(int j=signal_input.lrows;j<=signal_input.urows;j++){
        ind=k-j+1;
        if(j>=signal_input.lrows && j<=signal_input.urows &&
           ind>=response_input.lrows && ind<=response_input.urows){ //instead of zero-padding, we can check if the response_input index is valid
#if DEBUG_CONVOLUTION
          cerr << soliloquy << " k=" << k << " j=" << j << " ind=" << ind << endl;
#endif
          conv[k]+=signal_input[j]*response_input[ind];
          sum_counts_full[k]++;
        }else{ //keep k that require zero-padding (invalid indices), contains duplicates, but don't do work unless needed
          if(k_added==false){
            ind_zero_padding.push_back(k);
            k_added=true;
          }
        }
      }
    }
    if(LDEBUG){cerr << soliloquy << " full conv=" << conv << endl;}
    if(SHAPE==CONV_SHAPE_FULL){
      sum_counts.clear();for(uint i=0;i<sum_counts_full.size();i++){sum_counts.push_back(sum_counts_full[i]);}  //full copy
      return conv;
    }
    if(SHAPE==CONV_SHAPE_SAME){ //same - same size as signal_input, pick middle of full conv
      xvector<utype> conv_shape(signal_input.urows,signal_input.lrows);
      int ind1=((conv.rows-signal_input.rows)+1)/2+lrows;  //https://stackoverflow.com/questions/2745074/fast-ceiling-of-an-integer-division-in-c-c
      int ind2=lrows;
      if(LDEBUG){
        cerr << soliloquy << " ind1=" << ind1 << " conv[ind1]=" << conv[ind1] << endl;
        cerr << soliloquy << " ind2=" << ind2 << endl;
      }
      sum_counts.clear();
      for(int i=ind1;i<=ind1+signal_input.rows-1;i++){
        conv_shape[ind2++]=conv[i];
        sum_counts.push_back(sum_counts_full[i]);
      }
      if(LDEBUG){cerr << soliloquy << " same conv=" << conv_shape << endl;}
      return conv_shape;
    }
    else if(SHAPE==CONV_SHAPE_VALID){ //valid - only that section of conv that does NOT require zero-padding
      std::sort(ind_zero_padding.begin(),ind_zero_padding.end());ind_zero_padding.erase( std::unique( ind_zero_padding.begin(), ind_zero_padding.end() ), ind_zero_padding.end() );  //get unique only
      if(LDEBUG){cerr << soliloquy << " ind_zero_padding=" << aurostd::joinWDelimiter(ind_zero_padding,",") << endl;}
      size=conv.rows-ind_zero_padding.size();
      xvector<utype> conv_shape(size+(lrows-1),lrows);
      int ind2=lrows;
      sum_counts.clear();
      for(int i=conv.lrows;i<=conv.urows;i++){
        if(!aurostd::withinList(ind_zero_padding,i)){
          conv_shape[ind2++]=conv[i];
          sum_counts.push_back(sum_counts_full[i]);
        }
      }
      if(LDEBUG){cerr << soliloquy << " valid conv=" << conv_shape << endl;}
      return conv_shape;
    }
    else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"SHAPE specification unknown",_INPUT_UNKNOWN_);}
    return conv;
  }
  template<class utype> xvector<utype> moving_average(const xvector<utype>& signal_input,int window) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="aurostd::moving_average():";
    if(LDEBUG){
      cerr << soliloquy << " _signal_input=" << signal_input << endl;
      cerr << soliloquy << " window=" << window << endl;
    }
    //CO190622 - box_filter screws up edges
    //the averaging needs to be truncated for the edges
    //do NOT average beforehand
    //[CO190622 - box_filter screws up edges]xvector<utype> response_input=box_filter_xv<utype>(window,signal_input.lrows);  //note, padding response_input with 0s to make len same as signal_input will yield NO difference
    vector<uint> sum_counts;
    xvector<utype> response_input=ones_xv<utype>(window+(signal_input.lrows-1),signal_input.lrows);  //note, padding response_input with 0s to make len same as signal_input will yield NO difference
    xvector<utype> avg=convolution(signal_input,response_input,sum_counts,CONV_SHAPE_SAME);
    if(LDEBUG){
      cerr << soliloquy << " response_input=" << response_input << endl;
      cerr << soliloquy << " sum_counts=" << aurostd::joinWDelimiter(sum_counts,",") << endl;
    }
    for(int i=avg.lrows;i<=avg.urows;i++){avg[i]/=(utype)sum_counts[i-avg.lrows];}
    if(LDEBUG){cerr << soliloquy << " avg=" << avg << endl;}
    return avg;
  }
}

namespace aurostd {
  template<class utype> vector<int> getPeaks(const xvector<utype>& signal_input,uint smoothing_iterations,uint avg_window,int width_maximum,double significance_multiplier){  //CO190620
    xvector<utype> signal_smooth;
    return getPeaks(signal_input,signal_smooth,smoothing_iterations,avg_window,width_maximum,significance_multiplier);
  }
  template<class utype> vector<int> getPeaks(const xvector<utype>& signal_input,xvector<utype>& signal_smooth,uint smoothing_iterations,uint avg_window,int width_maximum,double significance_multiplier){  //CO190620
    //using method outlined here: https://dsp.stackexchange.com/questions/1302/peak-detection-approach
    //raw data is X
    //smooth data via moving average to get Y
    //stddev(X-Y) is sigma
    //detect peaks when (X-Y)>multiplier*sigma

    //smooth data
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="aurostd::getPeaks():";

    if(LDEBUG){
      cerr << soliloquy << " smoothing_iterations=" << smoothing_iterations << endl;
      cerr << soliloquy << " avg_window=" << avg_window << endl;
      cerr << soliloquy << " width_maximum=" << width_maximum << endl;
      cerr << soliloquy << " significance_multiplier=" << significance_multiplier << endl;
    }

    signal_smooth=signal_input;
    for(uint i=0;i<smoothing_iterations;i++){signal_smooth=aurostd::moving_average(signal_smooth,avg_window);}
   xvector<utype> diff=signal_input-signal_smooth;
    utype sigma=aurostd::stddev(diff);

    vector<int> peak_indices;
    bool local_maximum=false;
    bool significant=false;
    for(int i=signal_input.lrows;i<=signal_input.urows;i++){
      local_maximum=true;
      for(int j=1;j<=width_maximum&&local_maximum;j++){
        if(!((i-j)>=signal_input.lrows && (i+j)<=signal_input.urows && signal_input[i]>signal_input[i-j] && signal_input[i]>signal_input[i+j])){local_maximum=false;}
      }
      //[CO190620 - now width_maximum is a parameter]local_maximum=((i-1)>=signal_input.lrows && (i+1)<=signal_input.urows && signal_input[i]>signal_input[i-1] && signal_input[i]>signal_input[i+1]);
      significant=(diff[i]>significance_multiplier*sigma);
      if(local_maximum && significant){
        peak_indices.push_back(i);
        if(LDEBUG) {cerr << soliloquy << " PEAK[i=" << i << "]=" << signal_input[i] << endl;}
      }
    }
    return peak_indices;
  }
} // namespace aurostd

#endif  // _AUROSTD_XVECTOR_IMPLEMENTATIONS_


// ----------------------------------------------------------------------------
// --------------------------------------------------------------- constructors
// namespace aurostd {
//   // namespace aurostd
//   template<class utype>                                   // default constructor
//   xvector<utype>::xvector(void) {        
//     // allocate a xvector with subscript range [1..n]
//     int i;
//     lrows=1;
//     urows=_AUROSTD_XVECTOR_DEFAULT_SIZE_;
//     rows=urows-lrows+1;
// #ifdef _AUROSTD_XVECTOR_DEBUG_CONSTRUCTORS
//     printf("xxvector -> default constructor: lrows=%i, urows=%i,",lrows,urows);
// #endif
//     corpus=new utype[rows+XXEND];
//     if(!corpus) { 
//cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xvector<utype>: allocation failure in default constructor" << endl;
//exit(0);
//}
//     corpus+= -lrows+XXEND;
//     for(i=lrows;i<=urows;i++) corpus[i]=(utype) 0.0;  // clear
//     isfloat=aurostd::_isfloat((utype) 0);      
//     iscomplex=aurostd::_iscomplex((utype) 0);
//     size=(char) (sizeof(utype));
//     vsize=(long int) size*rows;
// #ifdef _AUROSTD_XVECTOR_DEBUG_CONSTRUCTORS
//     printf(" isfloat=%i, iscomplex=%i, sizeof=%i, vsize=%i\n",
// 	   isfloat,iscomplex,size,vsize);
// #endif
//   }
// }

// namespace aurostd {
//   // namespace aurostd
//   template<class utype>                                    // default constructor
//   xvector<utype>::xvector(int n) {        
//     // allocate a xvector with subscript range [1..n]
//     int i;
//     lrows=1;
//     urows=n;
//     rows=urows-lrows+1;
// #ifdef _AUROSTD_XVECTOR_DEBUG_CONSTRUCTORS
//     printf("xxvector -> default constructor: lrows=%i, urows=%i,",lrows,urows);
// #endif
//     corpus=new utype[rows+XXEND];
//     if(!corpus) {
//cerr << _AUROSTD_XLIBS_ERROR_ << "ERROR - aurostd::xvector<utype>: allocation failure in default constructor" << endl;
//exit(0);
//}
//     corpus+= -lrows+XXEND;
//     for(i=lrows;i<=urows;i++) corpus[i]=(utype) 0.0;  // clear
//     isfloat=aurostd::_isfloat((utype) 0);      
//     iscomplex=aurostd::_iscomplex((utype) 0);
//     size=(char) (sizeof(utype));
//     vsize=(long int) size*rows;
// #ifdef _AUROSTD_XVECTOR_DEBUG_CONSTRUCTORS
//     printf(" isfloat=%i, iscomplex=%i, sizeof=%i, vsize=%i\n",
// 	   isfloat,iscomplex,size,vsize);
// #endif
//   }
// }

// **************************************************************************
// *                                                                        *
// *             STEFANO CURTAROLO - Duke University 2003-2019              *
// *                                                                        *
// **************************************************************************
