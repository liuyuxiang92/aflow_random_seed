/*
    SymbolicC++ : An object oriented computer algebra system written in C++

    Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


// symmatrix.h

#ifndef SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX

#include "matrix.h"
//using namespace std; //DX20200625 - do not import entire namespace, now calling std functions when necessary (pair, bad_cast, list, ios, type_info, numeric_limits, and complex)

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX_FORWARD
#define SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX_FORWARD

namespace symbolic{ //DX20200625
class SymbolicMatrix;
} //namespace symbolic //DX20200625

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#ifndef SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX_DECLARE
#define SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX_DECLARE

namespace symbolic{ //DX20200625
class SymbolicMatrix
: public CloningSymbolicInterface, public Matrix<Symbolic>
{
 public: SymbolicMatrix(const SymbolicMatrix&);
         SymbolicMatrix(const Matrix<Symbolic>&);
         SymbolicMatrix(const std::list<std::list<Symbolic> >&);
         SymbolicMatrix(const string&,int,int);
         SymbolicMatrix(const char*,int,int);
         SymbolicMatrix(const Symbolic&,int,int);
         SymbolicMatrix(int,int);
         ~SymbolicMatrix();

         void print(ostream&) const;
         Symbolic subst(const Symbolic&,const Symbolic&,int &n) const;
         Simplified simplify() const;
         int compare(const Symbolic&) const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;
         Symbolic coeff(const Symbolic&) const;
         Expanded expand() const;
         int commute(const Symbolic&) const;
         PatternMatches match(const Symbolic&, const std::list<Symbolic>&) const;
         PatternMatches match_parts(const Symbolic&,
                                    const std::list<Symbolic>&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};
} //namespace symbolic //DX20200625

#endif
#endif

#define LIBSYMBOLICCPLUSPLUS

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX_DEFINE
#define SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX_DEFINE
#define SYMBOLIC_CPLUSPLUS_SYMBOLICMATRIX

namespace symbolic{ //DX20200625
SymbolicMatrix::SymbolicMatrix(const SymbolicMatrix &s)
: CloningSymbolicInterface(s), Matrix<Symbolic>(s) {}

SymbolicMatrix::SymbolicMatrix(const Matrix<Symbolic> &s)
: Matrix<Symbolic>(s)
{}

SymbolicMatrix::SymbolicMatrix(const std::list<std::list<Symbolic> > &sl)
{
 int cols = 0, k, l;
 std::list<Symbolic>::const_iterator j;
 std::list<std::list<Symbolic> >::const_iterator i;

 for(i=sl.begin();i!=sl.end();++i)
  if(int(i->size()) > cols) cols = i->size();

 Matrix<Symbolic>::resize(sl.size(),cols,Symbolic(0));
 for(k=0,i=sl.begin();i!=sl.end();++k,++i)
 {
  for(l=0,j=i->begin();j!=i->end();++l,++j)
   Matrix<Symbolic>::operator[](k)[l] = *j;
 }
}

SymbolicMatrix::SymbolicMatrix(const string &s,int n,int m)
 : Matrix<Symbolic>(n,m)
{
 for(int i=0;i<n;++i)
  for(int j=0;j<m;++j)
  {
   ostringstream os;
   if(n == 1 || m == 1)
    os << s << i + j;
   else
    os << s << "(" << i << "," << j << ")";
   Matrix<Symbolic>::operator[](i)[j] = Symbolic(os.str());
  }
}

SymbolicMatrix::SymbolicMatrix(const Symbolic &s,int n,int m)
 : Matrix<Symbolic>(n,m,s) {}

SymbolicMatrix::SymbolicMatrix(const char *s,int n,int m)
 : Matrix<Symbolic>(n,m)
{
 int i, j;
 for(i=0;i<n;++i)
  for(j=0;j<m;++j)
  {
   ostringstream os;
   if(n == 1 || m == 1)
    os << s << i + j;
   else
    os << s << "(" << i << "," << j << ")";
   Matrix<Symbolic>::operator[](i)[j] = Symbolic(os.str());
  }
}

SymbolicMatrix::SymbolicMatrix(int n,int m)
 : Matrix<Symbolic>(n,m,Symbolic(0)) {}

SymbolicMatrix::~SymbolicMatrix() {}

void SymbolicMatrix::print(ostream &o) const
{ o << endl << *this; }

Symbolic SymbolicMatrix::subst(const Symbolic &x,
                               const Symbolic &y,int &n) const
{
 if(*this == x) { ++n; return y; }

 SymbolicMatrix m(rows(),cols());
 int r, c;
 for(r = rows()-1;r>=0;r--)
  for(c = cols()-1;c>=0;c--)
   m[r][c] = Matrix<Symbolic>::operator[](r)[c].subst(x,y,n);

 return m;
}

Simplified SymbolicMatrix::simplify() const
{
 // single element matrix -> number
 if(rows() == 1 && cols() == 1)
  return Matrix<Symbolic>::operator[](0)[0].simplify();

 SymbolicMatrix m(rows(),cols());
 for(int r = rows()-1;r>=0;r--)
  for(int c = cols()-1;c>=0;c--)
   m[r][c] = Matrix<Symbolic>::operator[](r)[c].simplify();

 return m;
}

int SymbolicMatrix::compare(const Symbolic &s) const
{ 
 if(s.type() != type()) return 0;
 CastPtr<const SymbolicMatrix> m = s;
 return (*this) == (*m);
}

Symbolic SymbolicMatrix::df(const Symbolic &s) const
{
 SymbolicMatrix m(rows(),cols());
 for(int r = rows()-1;r>=0;r--)
  for(int c = cols()-1;c>=0;c--)
   m[r][c] = Matrix<Symbolic>::operator[](r)[c].df(s);

 return m;
}

Symbolic SymbolicMatrix::integrate(const Symbolic &s) const
{
 SymbolicMatrix m(rows(),cols());
 for(int r = rows()-1;r>=0;r--)
  for(int c = cols()-1;c>=0;c--)
   //DX20200625 - subst "::" with "symbolic::" - m[r][c] = ::integrate(Matrix<Symbolic>::operator[](r)[c],s);
   m[r][c] = symbolic::integrate(Matrix<Symbolic>::operator[](r)[c],s); //DX20200625 - subst "::" with "symbolic::"

 return m;
}

Symbolic SymbolicMatrix::coeff(const Symbolic &s) const
//DX 20190313 [OBSOLETE] { return 0; }
//DX 20190313 - need to use argument somehow, simple fix - START 
{ 
 if(s.type() == type()){ return 0; }
 return 0;
}
//DX 20190313 - need to use argument somehow, simple fix - END

Expanded SymbolicMatrix::expand() const
{
 SymbolicMatrix m(rows(),cols());
 for(int r = rows()-1;r>=0;r--)
  for(int c = cols()-1;c>=0;c--)
   m[r][c] = Matrix<Symbolic>::operator[](r)[c].expand();

 return m;
}

int SymbolicMatrix::commute(const Symbolic &s) const
{
 if(s.type() == typeid(SymbolicMatrix))
 {
  CastPtr<const SymbolicMatrix> m(s);
  return (*m) * (*this) - (*this) * (*m) == Matrix<Symbolic>(rows(),cols(),0);
 }
 return s.commute(*this);
}

PatternMatches
SymbolicMatrix::match(const Symbolic &s, const std::list<Symbolic> &p) const
{
 PatternMatches l;
 if(type() != s.type())
 {
  pattern_match_FALSE(l);
  return l;
 }

 CastPtr<SymbolicMatrix> sm(s);
 if(rows() != sm->rows() || cols() != sm->cols())
 {
  pattern_match_FALSE(l);
  return l;
 }

 for(int r = rows()-1;r>=0;r--)
  for(int c = cols()-1;c>=0;c--)
   if(r == rows()-1 && c == cols()-1)
    pattern_match_OR(l, Matrix<Symbolic>::operator[](r)[c].match(
                                                  sm->operator[](r)[c],p));
   else
    pattern_match_AND(l, Matrix<Symbolic>::operator[](r)[c].match(
                                                  sm->operator[](r)[c],p));
 return l;
}
 
PatternMatches
SymbolicMatrix::match_parts(const Symbolic &s, const std::list<Symbolic> &p) const
{
 PatternMatches l = s.match(*this, p);

 for(int r = rows()-1;r>=0;r--)
  for(int c = cols()-1;c>=0;c--)
   pattern_match_OR(l, Matrix<Symbolic>::operator[](r)[c].match_parts(s, p));
 return l;
}
} //namespace symbolic //DX20200625

#endif
#endif

#undef LIBSYMBOLICCPLUSPLUS

#endif
