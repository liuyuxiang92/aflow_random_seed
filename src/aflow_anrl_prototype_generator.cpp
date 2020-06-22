// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow DAVID HICKS - Duke University 2014-2020                 *
// *                                                                         *
// ***************************************************************************
// Written by David Hicks (DX) - 2020

#ifndef _AFLOW_ANRL_PROTOTYPE_GENERATOR_CPP
#define _AFLOW_ANRL_PROTOTYPE_GENERATOR_CPP
#include "aflow.h"
#include "aflow_anrl_prototype_generator.h"
#include "aflow_symmetry_spacegroup.h"
#include "aflow_compare_structure.h"
#include "aflow_pflow.h"

// *************************************************************************** 
// anrl::getLattice()
// *************************************************************************** 
namespace anrl {
  xmatrix<double> getLattice(const string& lattice_and_centering, 
      const char& space_group_letter, 
      const vector<double>& lattice_parameter_values, 
      uint mode){
    
    // Returns the lattice in the ANRL convention and populates with the 
    // relevant lattice parameters.
    // space_group_letter : needed to differentiate between the A, and C
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::getLattice():";
    stringstream message;

    if(LDEBUG){ cerr << function_name << " Lattice mode=" << mode; }

    // ---------------------------------------------------------------------------
    // triclinic crystal system
    if(lattice_and_centering == "aP"){
      return getTriclinicLattice(lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }
    
    // ---------------------------------------------------------------------------
    // monoclinic crystal system
    else if(lattice_and_centering == "mP" || lattice_and_centering == "mC"){
      return getMonoclinicLattice(lattice_and_centering, lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }
    
    // ---------------------------------------------------------------------------
    // orthorhombic crystal system
		else if(lattice_and_centering == "oP" || lattice_and_centering == "oC" || 
				lattice_and_centering == "oI" || lattice_and_centering == "oF"){
      return getOrthorhombicLattice(lattice_and_centering, space_group_letter, lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }
    
    // ---------------------------------------------------------------------------
    // tetragonal crystal system
    else if(lattice_and_centering == "tP" || lattice_and_centering == "tI"){
      return getTetragonalLattice(lattice_and_centering, lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }
    
    // ---------------------------------------------------------------------------
    // hexagonal crystal system (includes hex + rhl)
    else if(lattice_and_centering == "hP" || lattice_and_centering == "hR"){
      return getHexagonalLattice(lattice_and_centering, lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }
    
    // ---------------------------------------------------------------------------
    // cubic crystal system 
    else if(lattice_and_centering == "cP" || lattice_and_centering == "cF" || lattice_and_centering == "cI"){
      return getCubicLattice(lattice_and_centering, lattice_parameter_values, mode); // 0=primitive, 1=conventional
    }

    else{
      message << "Lattice type and centering are not possible: " << lattice_and_centering;
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_VALUE_ILLEGAL_);
    }
  }
}

// *************************************************************************** 
// anrl::getTriclinicLattice()
// *************************************************************************** 
namespace anrl {
  xmatrix<double> getTriclinicLattice(const vector<double>& lattice_parameter_values, 
      uint mode){
    
    // Returns the triclinic lattice in the ANRL convention and populates 
    // with the relevant lattice parameters.
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice
    // NOTE : the lattice_and_centering input is not required; triclinic 
    //        systems only have one centering option (P)

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::getTriclinicLattice():";
    stringstream message;

    // ---------------------------------------------------------------------------
    // check the number of inputs 
    if(lattice_parameter_values.size() != 6){
      message << "There needs to be 6 lattice parameters to build the triclinic lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
    }
    
    // ---------------------------------------------------------------------------
    // main variables 
    xmatrix<double> lattice;
    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);
    
    uint i=0;
    double a=lattice_parameter_values[i++];                  if(LDEBUG) { cerr << function_name << " a=" << a << endl;}
    double bovera=lattice_parameter_values[i++],b=bovera*a;  if(LDEBUG) { cerr << function_name << " b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=lattice_parameter_values[i++],c=covera*a;  if(LDEBUG) { cerr << function_name << " c=" << c << " (c/a=" << covera << ")" << endl;}
    double alpha=lattice_parameter_values[i++];              if(LDEBUG) { cerr << function_name << " alpha=" << alpha << endl;}
    double beta=lattice_parameter_values[i++];               if(LDEBUG) { cerr << function_name << " beta=" << beta << endl;}
    double gamma=lattice_parameter_values[i++];              if(LDEBUG) { cerr << function_name << " gamma=" << gamma << endl;}

    // ---------------------------------------------------------------------------
    // triclinic - tri (aP)
    if(mode ==0 || mode == 1){ // primitive and conventional cells are the same
    
      double cx=c*cos(deg2rad*beta);
      double cy=c*(cos(deg2rad*alpha)-cos(deg2rad*beta)*cos(deg2rad*gamma))/sin(deg2rad*gamma);
      double cz=sqrt(pow(c,2.0)-pow(cx,2.0)-pow(cy,2.0));

      a1=a*xn;
      a2=b*cos(deg2rad*gamma)*xn+b*sin(deg2rad*gamma)*yn;
      a3=cx*xn+cy*yn+cz*zn;

      // ---------------------------------------------------------------------------
      // build lattice 
      lattice(1,1)=a1(1);lattice(1,2)=a1(2);lattice(1,3)=a1(3);
      lattice(2,1)=a2(1);lattice(2,2)=a2(2);lattice(2,3)=a2(3);
      lattice(3,1)=a3(1);lattice(3,2)=a3(2);lattice(3,3)=a3(3);
	 
    }

    if(LDEBUG){ cerr << function_name << " lattice = " << lattice << endl; }

    return lattice;
  }
}

// *************************************************************************** 
// anrl::getMonoclinicLattice()
// *************************************************************************** 
namespace anrl {
  xmatrix<double> getMonoclinicLattice(const string& lattice_and_centering, 
      const vector<double>& lattice_parameter_values, 
      uint mode){
    
    // Returns the monoclinic lattice in the ANRL convention and populates 
    // with the relevant lattice parameters.
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice
    
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::getMonoclinicLattice():";
    stringstream message;
    
    // ---------------------------------------------------------------------------
    // check the number of inputs 
    if(lattice_parameter_values.size() != 4){
      message << "There needs to be 4 lattice parameters to build the monoclinic lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
    }
    
    // ---------------------------------------------------------------------------
    // main variables 
    xmatrix<double> lattice;
    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);
    
    uint i=0;
    double a=lattice_parameter_values[i++];                  if(LDEBUG) { cerr << function_name << " a=" << a << endl;}
    double bovera=lattice_parameter_values[i++],b=bovera*a;  if(LDEBUG) { cerr << function_name << " b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=lattice_parameter_values[i++],c=covera*a;  if(LDEBUG) { cerr << function_name << " c=" << c << " (c/a=" << covera << ")" << endl;}
    double beta=lattice_parameter_values[i++];               if(LDEBUG) { cerr << function_name << " beta=" << beta << endl;}

    // ---------------------------------------------------------------------------
    // primitive lattices
    if(mode == 0){
      // ---------------------------------------------------------------------------
      // simple monoclinic - mcl (mP)
      if(lattice_and_centering == "mP"){
        a1=a*xn;
        a2=b*yn;
        a3=c*cos(deg2rad*beta)*xn+c*sin(deg2rad*beta)*zn;
      } 
      
      // ---------------------------------------------------------------------------
      // base-centered monoclinic - mclc (mC)
      if(lattice_and_centering == "mC"){
        a1=(1.0/2.0)*a*xn-(1.0/2.0)*b*yn;
        a2=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
        a3=c*cos(deg2rad*beta)*xn+c*sin(deg2rad*beta)*zn;
      } 
    }

    // ---------------------------------------------------------------------------
    // conventional lattice
    else if(mode == 1){
      a1=a*xn;
      a2=b*yn;
      a3=c*cos(deg2rad*beta)*xn+c*sin(deg2rad*beta)*zn;
    }

    // ---------------------------------------------------------------------------
    // build lattice 
    lattice(1,1)=a1(1);lattice(1,2)=a1(2);lattice(1,3)=a1(3);
    lattice(2,1)=a2(1);lattice(2,2)=a2(2);lattice(2,3)=a2(3);
    lattice(3,1)=a3(1);lattice(3,2)=a3(2);lattice(3,3)=a3(3);
    
    if(LDEBUG){ cerr << function_name << " lattice = " << lattice << endl; }

    return lattice;
  }
}

// *************************************************************************** 
// anrl::getOrthorhombicLattice()
// *************************************************************************** 
namespace anrl {
  xmatrix<double> getOrthorhombicLattice(const string& lattice_and_centering,
      const char& space_group_letter,
      const vector<double>& lattice_parameter_values,
      uint mode){
    
    // Returns the orthorhombic lattice in the ANRL convention and populates 
    // with the relevant lattice parameters.
    // space_group_letter : needed to differentiate between the A, and C
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::getOrthorhombicLattice():";
    stringstream message;

    // ---------------------------------------------------------------------------
    // check the number of inputs 
    if(lattice_parameter_values.size() != 3){
      message << "There needs to be 3 lattice parameters to build the orthorhombic lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
    }
    
    // ---------------------------------------------------------------------------
    // main variables 
    xmatrix<double> lattice;
    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);
    
    uint i=0;
    double a=lattice_parameter_values[i++];                  if(LDEBUG) { cerr << function_name << " a=" << a << endl;}
    double bovera=lattice_parameter_values[i++],b=bovera*a;  if(LDEBUG) { cerr << function_name << " b=" << b << " (b/a=" << bovera << ")" << endl;}
    double covera=lattice_parameter_values[i++],c=covera*a;  if(LDEBUG) { cerr << function_name << " c=" << c << " (c/a=" << covera << ")" << endl;}

    // ---------------------------------------------------------------------------
    // primitive lattices
    if(mode == 0){
      // ---------------------------------------------------------------------------
      // simple orthorhombic - orc (oP)
      if(lattice_and_centering == "oP"){
        a1=a*xn;
        a2=b*yn;
        a3=c*zn;  
      } 
      
      // ---------------------------------------------------------------------------
      // base-centered orthorhombic - orcc (oC)
      if(lattice_and_centering == "oC"){
        // A-centered
        if(space_group_letter == 'A'){
          a1=a*xn;
          a2=(1.0/2.0)*b*yn-(1.0/2.0)*c*zn;
          a3=(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
        }
        // C-centered
        if(space_group_letter == 'C'){
          a1=(1.0/2.0)*a*xn-(1.0/2.0)*b*yn;
          a2=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
          a3=c*zn;
        }
      }
      
      // ---------------------------------------------------------------------------
      // body-centered orthorhombic - orci (oI)
      if(lattice_and_centering == "oI"){
        a1=-(1.0/2.0)*a*xn+(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
        a2=(1.0/2.0)*a*xn-(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
        a3=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn-(1.0/2.0)*c*zn;
      } 
      
      // ---------------------------------------------------------------------------
      // face-centered orthorhombic - orcf (oF)
      if(lattice_and_centering == "oF"){
        a1=(1.0/2.0)*b*yn+(1.0/2.0)*c*zn;
        a2=(1.0/2.0)*a*xn+(1.0/2.0)*c*zn;
        a3=(1.0/2.0)*a*xn+(1.0/2.0)*b*yn;
      } 
    }

    // ---------------------------------------------------------------------------
    // conventional lattice
    else if(mode == 1){
      a1=a*xn;
      a2=b*yn;
      a3=c*zn;  
    }

    // ---------------------------------------------------------------------------
    // build lattice 
    lattice(1,1)=a1(1);lattice(1,2)=a1(2);lattice(1,3)=a1(3);
    lattice(2,1)=a2(1);lattice(2,2)=a2(2);lattice(2,3)=a2(3);
    lattice(3,1)=a3(1);lattice(3,2)=a3(2);lattice(3,3)=a3(3);
    
    if(LDEBUG){ cerr << function_name << " lattice = " << lattice << endl; }

    return lattice;
  }
}

// *************************************************************************** 
// anrl::getTetragonaLattice()
// *************************************************************************** 
namespace anrl {
  xmatrix<double> getTetragonalLattice(const string& lattice_and_centering,
      const vector<double>& lattice_parameter_values,
      uint mode){
    
    // Returns the tetragonal lattice in the ANRL convention and populates 
    // with the relevant lattice parameters.
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::getTetragonalLattice():";
    stringstream message;

    // ---------------------------------------------------------------------------
    // check the number of inputs 
    if(lattice_parameter_values.size() != 2){
      message << "There needs to be 2 lattice parameters to build the tetragonal lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
    }
    
    // ---------------------------------------------------------------------------
    // main variables 
    xmatrix<double> lattice;
    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);
    
    uint i=0;
    double a=lattice_parameter_values[i++];                  if(LDEBUG) { cerr << function_name << " a=" << a << endl;}
    double covera=lattice_parameter_values[i++],c=covera*a;  if(LDEBUG) { cerr << function_name << " c=" << c << " (c/a=" << covera << ")" << endl;}

    // ---------------------------------------------------------------------------
    // primitive lattices
    if(mode == 0){
      // ---------------------------------------------------------------------------
      // simple tetragonal - tet (tP)
      if(lattice_and_centering == "tP"){
        a1=a*xn;
        a2=a*yn;
        a3=c*zn;
      } 
      
      // ---------------------------------------------------------------------------
      // body-centered tegtragonal - bct (tI)
      if(lattice_and_centering == "tI"){
        a1=-(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*c*zn;
        a2=(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*c*zn;
        a3=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn-(1.0/2.0)*c*zn;
      } 
    }
    // ---------------------------------------------------------------------------
    // conventional lattice
    else if(mode == 1){
      a1=a*xn;
      a2=a*yn;
      a3=c*zn;
    }

    // ---------------------------------------------------------------------------
    // build lattice 
    lattice(1,1)=a1(1);lattice(1,2)=a1(2);lattice(1,3)=a1(3);
    lattice(2,1)=a2(1);lattice(2,2)=a2(2);lattice(2,3)=a2(3);
    lattice(3,1)=a3(1);lattice(3,2)=a3(2);lattice(3,3)=a3(3);
    
    if(LDEBUG){ cerr << function_name << " lattice = " << lattice << endl; }

    return lattice;
  }
}

// *************************************************************************** 
// anrl::getHexagonalLattice()
// *************************************************************************** 
namespace anrl {
  xmatrix<double> getHexagonalLattice(const string& lattice_and_centering,
      const vector<double>& lattice_parameter_values,
      uint mode){
    
    // Returns the hexagonal/rhombohedral lattice in the ANRL convention and 
    // populates with the relevant lattice parameters.
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::getHexagonalLattice():";
    stringstream message;

    // ---------------------------------------------------------------------------
    // check number of inputs
    if(lattice_parameter_values.size() != 2){
      message << "There needs to be 2 lattice parameters to build the hexagonal/rhombohedral lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
    }
    
    // ---------------------------------------------------------------------------
    // main variables 
    xmatrix<double> lattice;
    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);
    
    uint i=0;
    double a=lattice_parameter_values[i++];                  if(LDEBUG) { cerr << function_name << " a=" << a << endl;}
    double covera=lattice_parameter_values[i++],c=covera*a;  if(LDEBUG) { cerr << function_name << " c=" << c << " (c/a=" << covera << ")" << endl;}

    // ---------------------------------------------------------------------------
    // primitive lattices
    if(mode == 0){
      // ---------------------------------------------------------------------------
      // hexagonal - hex (hP)
      if(lattice_and_centering == "hP"){
        a1=(1.0/2.0)*a*xn-(sqrt(3.0)/2.0)*a*yn;
        a2=(1.0/2.0)*a*xn+(sqrt(3.0)/2.0)*a*yn;
        a3=c*zn;
      } 
      
      // ---------------------------------------------------------------------------
      // rhombohedral - rhl (hR)
      if(lattice_and_centering == "hR"){
        a1=(1.0/2.0)*a*xn-(1.0/(2.0*sqrt(3.0)))*a*yn+(1.0/3.0)*c*zn;
        a2=(1.0/sqrt(3.0))*a*yn+(1.0/3.0)*c*zn;
        a3=-(1.0/2.0)*a*xn-(1.0/(2.0*sqrt(3.0)))*a*yn+(1.0/3.0)*c*zn;
      } 
    }
    // ---------------------------------------------------------------------------
    // conventional lattice
    else if(mode == 1){
      a1=(1.0/2.0)*a*xn-(sqrt(3.0)/2.0)*a*yn;
      a2=(1.0/2.0)*a*xn+(sqrt(3.0)/2.0)*a*yn;
      a3=c*zn;
    }

    // ---------------------------------------------------------------------------
    // build lattice
    lattice(1,1)=a1(1);lattice(1,2)=a1(2);lattice(1,3)=a1(3);
    lattice(2,1)=a2(1);lattice(2,2)=a2(2);lattice(2,3)=a2(3);
    lattice(3,1)=a3(1);lattice(3,2)=a3(2);lattice(3,3)=a3(3);
    
    if(LDEBUG){ cerr << function_name << " lattice = " << lattice << endl; }

    return lattice;
  }
}

// *************************************************************************** 
// anrl::getCubicLattice()
// *************************************************************************** 
namespace anrl {
  xmatrix<double> getCubicLattice(const string& lattice_and_centering,
      const vector<double>& lattice_parameter_values,
      uint mode){
    
    // Returns the cubic lattice in the ANRL convention and 
    // populates with the relevant lattice parameters.
    // lattice centering conventions
    // mode : specifies primitive(=0) or the conventional(=1) lattice

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::getCubicLattice():";
    stringstream message;

    // ---------------------------------------------------------------------------
    // check number of inputs
    if(lattice_parameter_values.size() != 1){
      message << "There needs to be 1 lattice parameters to build the cubic lattice (input size=" << lattice_parameter_values.size() << ")";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
    }
    
    // ---------------------------------------------------------------------------
    // main variables 
    xmatrix<double> lattice;
    xvector<double> xn(3);   xn(1)=1.0;xn(2)=0.0;xn(3)=0.0;
    xvector<double> yn(3);   yn(1)=0.0;yn(2)=1.0;yn(3)=0.0;
    xvector<double> zn(3);   zn(1)=0.0;zn(2)=0.0;zn(3)=1.0;
    xvector<double> a1(3),a2(3),a3(3);
    
    uint i=0;
    double a=lattice_parameter_values[i++];                  if(LDEBUG) { cerr << function_name << " a=" << a << endl;}

    // ---------------------------------------------------------------------------
    // primitive lattices
    if(mode == 0){
      // ---------------------------------------------------------------------------
      // simple cubic - cub (cP)
      if(lattice_and_centering == "cP"){
        a1=a*xn;
        a2=a*yn;
        a3=a*zn;
      } 
      
      // ---------------------------------------------------------------------------
      // face-centered cubic - fcc (cF)
      if(lattice_and_centering == "cF"){
        a1=(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
        a2=(1.0/2.0)*a*xn+(1.0/2.0)*a*zn;
        a3=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn;
      } 
      
      // ---------------------------------------------------------------------------
      // body-centered cubic - bcc (cI)
      if(lattice_and_centering == "cI"){
        a1=-(1.0/2.0)*a*xn+(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
        a2=(1.0/2.0)*a*xn-(1.0/2.0)*a*yn+(1.0/2.0)*a*zn;
        a3=(1.0/2.0)*a*xn+(1.0/2.0)*a*yn-(1.0/2.0)*a*zn;
      } 
    }
    
    // ---------------------------------------------------------------------------
    // conventional lattice
    else if(mode == 1){
      a1=a*xn;
      a2=a*yn;
      a3=a*zn;
    }

    // ---------------------------------------------------------------------------
    // build lattice 
    lattice(1,1)=a1(1);lattice(1,2)=a1(2);lattice(1,3)=a1(3);
    lattice(2,1)=a2(1);lattice(2,2)=a2(2);lattice(2,3)=a2(3);
    lattice(3,1)=a3(1);lattice(3,2)=a3(2);lattice(3,3)=a3(3);
    
    if(LDEBUG){ cerr << function_name << " lattice = " << lattice << endl; }

    return lattice;
  }
}

// *************************************************************************** 
// anrl::getAtomsFromWyckoff()
// *************************************************************************** 
namespace anrl {
  deque<_atom> getAtomsFromWyckoff(const vector<wyckoffsite_ITC>& Wyckoff_positions,
      const xmatrix<double>& lattice_conventional){
    
    // create atoms (deque<_atom>) from the Wyckoff positions by 
    // plugging in the parameter values into the Wyckoff equations 

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::getAtomsFromWyckoff():";
    stringstream message;

    // ---------------------------------------------------------------------------
    // variables 
    deque<_atom> atoms_conventional_cell;
    _atom atom_tmp;
    
    // ---------------------------------------------------------------------------
    // create an _atom for each Wyckoff position by plugging in the relevant parameters  
    for(uint i=0;i<Wyckoff_positions.size();i++){
      // get x, y, and z coordinates from the Wyckoff object
      string x_value_string = aurostd::utype2string<double>(Wyckoff_positions[i].coord(1));
      string y_value_string = aurostd::utype2string<double>(Wyckoff_positions[i].coord(2));
      string z_value_string = aurostd::utype2string<double>(Wyckoff_positions[i].coord(3));
      for(uint j=0;j<Wyckoff_positions[i].equations.size();j++){
        vector<string> coordinate_vstring = Wyckoff_positions[i].equations[j];
        xvector<double> coordinate(3);
        for(uint k=0;k<coordinate_vstring.size();k++){
          // substitute variable with value
          aurostd::StringSubst(coordinate_vstring[k],"x",x_value_string);
          aurostd::StringSubst(coordinate_vstring[k],"y",y_value_string);
          aurostd::StringSubst(coordinate_vstring[k],"z",z_value_string);
          // simplify string 
          vector<SYM::sdouble> component_tmp = SYM::simplify(coordinate_vstring[k]);
          string component_string = SYM::formatWyckoffPosition(component_tmp,false);
          // ensure the string is numeric
          if(!aurostd::isfloat(component_string)){
            message << "There are non-numeric characters in the string after variable-value substitution: component " << k << "=" << component_string << endl;
            throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_VALUE_ERROR_);
          }
          coordinate(k+1) = aurostd::string2utype<double>(component_string);
        }
        // store atom info
        atom_tmp.fpos = coordinate;
        atom_tmp.cpos = F2C(lattice_conventional,coordinate);
        atom_tmp.type = Wyckoff_positions[i].index; 
        atom_tmp.name = Wyckoff_positions[i].type;
        atoms_conventional_cell.push_back(atom_tmp);
      }
    }
    if(LDEBUG){
      cerr << function_name << " atoms in the conventional cell:" << endl;
      for(uint i=0;i<atoms_conventional_cell.size();i++){
        cerr << "atoms: " << atoms_conventional_cell[i] << " " << atoms_conventional_cell[i].name << endl;
      }
    }
    return atoms_conventional_cell;
  }
}

// *************************************************************************** 
// anrl::determineWyckoffVariables()
// *************************************************************************** 
namespace anrl {
	vector<string> determineWyckoffVariables(
			vector<wyckoffsite_ITC>& Wyckoff_positions){

    // Determines the variable coordinates in the Wyckoff positions and
    // returns a vector of Wyckoff variables (x, y, or z) that need to be 
    // specified. The Wyckoff positions should be ordered by Wyckoff letter
    // (alphabetic) in accordance with the ANRL convention.
    // The format is : x1, y1, z1, x2, y2, z2, x3, ...
    // Note: Wyckoff_positions are updated, i.e., the corresponding parameter
    // index is assigned so we can substitute values later

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::determineWyckoffVariables():";
    stringstream message;

    // ---------------------------------------------------------------------------
    // variables
    vector<string> Wyckoff_parameter_list;
    //uint w=0; // Wyckoff counter
    
    for(uint i=0;i<Wyckoff_positions.size();i++){
      bool contains_x=false; bool contains_y=false; bool contains_z=false;
      if(Wyckoff_positions[i].equations.size()>0){
        // ---------------------------------------------------------------------------
        // look at the representative Wyckoff position (index 0) and see if it 
        // x, y, and/or z
        for(uint j=0;j<Wyckoff_positions[i].equations[0].size();j++){
          if(aurostd::substring2bool(Wyckoff_positions[i].equations[0][j],"x")){ contains_x=true; }
          if(aurostd::substring2bool(Wyckoff_positions[i].equations[0][j],"y")){ contains_y=true; }
          if(aurostd::substring2bool(Wyckoff_positions[i].equations[0][j],"z")){ contains_z=true; }
        }
        
        // ---------------------------------------------------------------------------
        // check for x-coordinate 
				if(contains_x){
          string variable = "x"; 
					string variable_name = variable+aurostd::utype2string<uint>(i+1);
					Wyckoff_parameter_list.push_back(variable_name);
          Wyckoff_positions[i].parameter_index=i+1;
          //for(uint w=0;w<Wyckoff_positions[i].equations.size();w++){
          //  for(uint e=0;e<Wyckoff_positions[i].equations[w].size();e++){
          //    aurostd::StringSubst(Wyckoff_positions[i].equations[w][e],variable,variable_name);
          //  }
          //}
					//if(w<Wyckoff_parameter_values.size()){ Wyckoff_positions[i].coord(1) = Wyckoff_parameter_values[w++]; }
          // store parameter value
					//if(w<Wyckoff_parameter_values.size()){ Wyckoff_positions[i].coord(1) = Wyckoff_parameter_values[w++]; }
          //else{
          //  message << "There are too few input parameters; could not populate the x-coordinate for Wyckoff position " << i;
          //  throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
          //}
				}
        // ---------------------------------------------------------------------------
        // check for x-coordinate 
				if(contains_y){ 
          string variable = "y"; 
					string variable_name = variable+aurostd::utype2string<uint>(i+1);
					Wyckoff_parameter_list.push_back(variable_name);
          Wyckoff_positions[i].parameter_index=i+1;
          //for(uint w=0;w<Wyckoff_positions[i].equations.size();w++){
          //  for(uint e=0;e<Wyckoff_positions[i].equations[w].size();e++){
          //    aurostd::StringSubst(Wyckoff_positions[i].equations[w][e],variable,variable_name);
          //  }
          //}
					//if(w<Wyckoff_parameter_values.size()){ Wyckoff_positions[i].coord(2) = Wyckoff_parameter_values[w++]; }
          //else{
          //  message << "There are too few input parameters; could not populate the y-coordinate for Wyckoff position " << i;
          //  throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
          //}
				}
        // ---------------------------------------------------------------------------
        // check for x-coordinate 
				if(contains_z){ 
          string variable = "z"; 
					string variable_name = variable+aurostd::utype2string<uint>(i+1);
					Wyckoff_parameter_list.push_back(variable_name);
          Wyckoff_positions[i].parameter_index=i+1;
          //for(uint w=0;w<Wyckoff_positions[i].equations.size();w++){
          //  for(uint e=0;e<Wyckoff_positions[i].equations[w].size();e++){
          //    aurostd::StringSubst(Wyckoff_positions[i].equations[w][e],variable,variable_name);
          //  }
          //}
					//if(w<Wyckoff_parameter_values.size()){ Wyckoff_positions[i].coord(3) = Wyckoff_parameter_values[w++]; }
          //else{
          //  message << "There are too few input parameters; could not populate the z-coordinate for Wyckoff position " << i;
          //  throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
          //}
				}
      }
      else{
        message << "The equations for site " << i << "are not provided. Check symmetry.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_RUNTIME_ERROR_);
      }
    }

    if(LDEBUG){ cerr << function_name << " parameters=" << aurostd::joinWDelimiter(Wyckoff_parameter_list,",") << endl; }

    return Wyckoff_parameter_list;
  }
}

// *************************************************************************** 
// anrl::applyWyckoffValues()
// *************************************************************************** 
namespace anrl {
  void applyWyckoffValues(
      const vector<double>& Wyckoff_parameter_values,
      vector<wyckoffsite_ITC>& Wyckoff_positions){
    
    // Applies the Wyckoff parameter values to the Wyckoff object
    // i.e., the coord attribute, indicating the degree of freedom (x, y, or z) 
    // The Wyckoff positions should be ordered by Wyckoff letter
    // (alphabetic) in accordance with the ANRL convention.
    // The format is : x1, y1, z1, x2, y2, z2, x3, ...
    
    string function_name = XPID + "anrl::applyWyckoffValues():";
    stringstream message;
    
    // ---------------------------------------------------------------------------
    // variables
    vector<string> Wyckoff_parameter_list;
    uint w=0; // Wyckoff counter
    
    for(uint i=0;i<Wyckoff_positions.size();i++){
      bool contains_x=false; bool contains_y=false; bool contains_z=false;
      if(Wyckoff_positions[i].equations.size()>0){
        // ---------------------------------------------------------------------------
        // look at the representative Wyckoff position (index 0) and see if it 
        // x, y, and/or z
        for(uint j=0;j<Wyckoff_positions[i].equations[0].size();j++){
          if(aurostd::substring2bool(Wyckoff_positions[i].equations[0][j],"x")){ contains_x=true; }
          if(aurostd::substring2bool(Wyckoff_positions[i].equations[0][j],"y")){ contains_y=true; }
          if(aurostd::substring2bool(Wyckoff_positions[i].equations[0][j],"z")){ contains_z=true; }
        }
        
        // ---------------------------------------------------------------------------
        // check for x-coordinate 
				if(contains_x){ 
          // store parameter value
					if(w<Wyckoff_parameter_values.size()){ Wyckoff_positions[i].coord(1) = Wyckoff_parameter_values[w++]; }
          else{
            message << "There are too few input parameters; could not populate the x-coordinate for Wyckoff position " << i;
            throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
          }
				}
        // ---------------------------------------------------------------------------
        // check for x-coordinate 
				if(contains_y){ 
          // store parameter value
					if(w<Wyckoff_parameter_values.size()){ Wyckoff_positions[i].coord(2) = Wyckoff_parameter_values[w++]; }
          else{
            message << "There are too few input parameters; could not populate the y-coordinate for Wyckoff position " << i;
            throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
          }
				}
        // ---------------------------------------------------------------------------
        // check for x-coordinate 
				if(contains_z){ 
          // store parameter value
					if(w<Wyckoff_parameter_values.size()){ Wyckoff_positions[i].coord(3) = Wyckoff_parameter_values[w++]; }
          else{
            message << "There are too few input parameters; could not populate the z-coordinate for Wyckoff position " << i;
            throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
          }
				}
      }
      else{
        message << "The equations for site " << i << "are not provided. Check symmetry.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_RUNTIME_ERROR_);
      }
    }
  }
}

// *************************************************************************** 
// anrl::containsDuplicateWyckoffCoordinate()
// *************************************************************************** 
namespace anrl {
  bool containsDuplicateWyckoffCoordinate(const vector<wyckoffsite_ITC>& wyckoff_sites_ITC, bool already_ordered){
   
		// Checks if Wyckoff positions occur in the structure multiple times
    // Two cases:
    //   1) multiple instances of a fixed Wyckoff position (i.e., no variables) 
    //   2) multiple instances of a variable Wyckoff position with the same 
    //      parameters

    string function_name = XPID + "anrl::containsDuplicateWyckoffCoordinate():";
    stringstream message;

    // ---------------------------------------------------------------------------
    double _WYCKOFF_FRACTIONAL_TOL_ = 1e-6;
		
    vector<wyckoffsite_ITC> ordered_Wyckoff_sites = wyckoff_sites_ITC;
    if(!already_ordered){ std::sort(ordered_Wyckoff_sites.begin(), ordered_Wyckoff_sites.end(), sortWyckoffByLetter); } 

    for(uint i=0;i<ordered_Wyckoff_sites.size();i++){ 
      for(uint j=i+1;j<ordered_Wyckoff_sites.size();j++){ 
        if(ordered_Wyckoff_sites[i].letter == ordered_Wyckoff_sites[j].letter){
          // ---------------------------------------------------------------------------
          // case 1: no variables in representative Wyckoff positions (first equation)
          // means we should only have one instance of this Wyckoff position
          if(!aurostd::substring2bool(ordered_Wyckoff_sites[i].equations[0],"x") &&
              !aurostd::substring2bool(ordered_Wyckoff_sites[i].equations[0],"y") &&
              !aurostd::substring2bool(ordered_Wyckoff_sites[i].equations[0],"z")){
            message << "Contains multiple static (i.e., no variable) Wyckoff positions: " << ordered_Wyckoff_sites[i].letter;
            throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
            return true;
          }
          // ---------------------------------------------------------------------------
          // case 2: contains variables, but values are the same; this is a quick check 
          // (i.e., use fractional tol), more rigorous check done later
          else if(aurostd::isequal(ordered_Wyckoff_sites[i].coord(1),ordered_Wyckoff_sites[j].coord(1),_WYCKOFF_FRACTIONAL_TOL_) &&
              aurostd::isequal(ordered_Wyckoff_sites[i].coord(2),ordered_Wyckoff_sites[j].coord(2),_WYCKOFF_FRACTIONAL_TOL_) &&
              aurostd::isequal(ordered_Wyckoff_sites[i].coord(3),ordered_Wyckoff_sites[j].coord(3),_WYCKOFF_FRACTIONAL_TOL_)){
            message << "Contains duplicate Wyckoff letters with the same degrees of freedom: " << aurostd::joinWDelimiter(xvecDouble2vecString(ordered_Wyckoff_sites[i].coord),",");
            throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
            return true;
          }
        }
        // ---------------------------------------------------------------------------
        // since ordered by Wyckoff letter, we can skip the rest (break) and start
        // where we left off in the second loop (i=j)
        else{ i=j; break; }
      }
    }
    return false;
  }
}

// *************************************************************************** 
// anrl::getWyckoffSitesFromANRL()
// *************************************************************************** 
namespace anrl {
  vector<wyckoffsite_ITC> getWyckoffSitesFromANRL(
      const vector<string>& Wyckoff_tokens,
      const vector<string>& species,
      uint space_group_number,
      int setting){

    // Get Wyckoff positions/sites from the ANRL Wyckoff designation
    // i.e., given Wyckoff letter and number of times they are used 
    // (e.g., 4a, 2b, c) get the Wyckoff letter, multiplicity, site symmetry, 
    // and equations

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::getWyckoffSitesFromANRL():";
    stringstream message;

    vector<wyckoffsite_ITC> wyckoff_sites_ITC;
    for(uint i=0;i<Wyckoff_tokens.size();i++){

			if(LDEBUG){ cerr << function_name << " Wyckoff designation=" << Wyckoff_tokens[i] << endl; }

      uint Wyckoff_multiplication_factor = 1;
      stringstream ss_Wyckoff_letter, ss_factor;
        
      for(uint j=0;j<Wyckoff_tokens[i].size();j++){
        // ---------------------------------------------------------------------------
        // extract the prefactor (if it exists; prefactor=1 is not usually given)
        if(isdigit(Wyckoff_tokens[i][j])){
          ss_factor << Wyckoff_tokens[i][j];
          continue;
        }
        // ---------------------------------------------------------------------------
        // extract the Wyckoff letter 
        else{ 
          ss_Wyckoff_letter << Wyckoff_tokens[i][j]; 
          if(ss_factor.str().size()){ 
						Wyckoff_multiplication_factor=aurostd::string2utype<uint>(ss_factor.str()); 
          }
        }
				if(LDEBUG){ cerr << function_name << " " << Wyckoff_multiplication_factor << " Wyckoff position(s) with letter" << ss_Wyckoff_letter.str() << endl; }
        
				// ---------------------------------------------------------------------------
        // populates Wyckoff with corresponding letter, multiplicity, equations, etc. 
				wyckoffsite_ITC Wyckoff_tmp;
        Wyckoff_tmp.getWyckoffFromLetter(space_group_number, ss_Wyckoff_letter.str(), setting);
        Wyckoff_tmp.index = i; 
        Wyckoff_tmp.type = species[i]; 
				if(LDEBUG){ cerr << function_name << " extracted Wyckoff position:" << Wyckoff_tmp << endl; }
			
        // ---------------------------------------------------------------------------
        // store the Wyckoff position as indicated by the prefactor, 
        // e.g., 4b-> four Wyckoff positions with letter b
        for(uint m=0;m<Wyckoff_multiplication_factor;m++){
          wyckoff_sites_ITC.push_back(Wyckoff_tmp);
        }
        // reset
        Wyckoff_multiplication_factor=1;
        ss_Wyckoff_letter.str("");
        ss_factor.str("");
      }
    }

    return wyckoff_sites_ITC;
  }
}

// *************************************************************************** 
// anrl::extractPrototypeParameters()
// *************************************************************************** 
namespace anrl {
  string extractANRLPrototypeParameterValues(const string& label_anrl, 
      const string& number_id, 
      const string& variables, 
      bool& keep_anrl_lattice_parameter){ 

    // bool must be reference: toggles automatic volume scaling later

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::extractANRLPrototypeParameterValues():";
    stringstream message;

    // index 
    int choice = -1;
      
    // ---------------------------------------------------------------------------
    // only one degree of freedom
    if(variables == "a"){
      choice = 0; // only one possiblity
      if(number_id.size()!=0){ 
        message << label_anrl << " only has one degree of freedom (lattice parameter), i.e., no enumerated suffix necessary." << endl; 
        throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
      }
    }
    // ---------------------------------------------------------------------------
    // if number ID is given, i.e., 001, 002, etc.
    else if(number_id.size()!=0){
      choice = aurostd::string2utype<uint>(number_id) - 1; //number to index
			if(LDEBUG){ cerr << function_name << " id=" << number_id << " --> choice=" << choice << endl; }
    }

    // ---------------------------------------------------------------------------
    // extract parameters for a given label 
    vector<string> all_possible_vparameters = getANRLParameters(label_anrl, "", choice, keep_anrl_lattice_parameter);
    string parameters=all_possible_vparameters[0];
  
    // ---------------------------------------------------------------------------
    // extract parameters for a given label 
    vector<string> vparameters_library;
    aurostd::string2tokens(parameters,vparameters_library,",");
    
    if(vparameters_library.size()==0){
      message << "No parameters provided; add parameter values with --params=... or use tabulated enumeration suffix (see aflow --protos)" << endl;
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
    }

    return parameters;
  }
}

// *************************************************************************** 
// anrl::structureAndLabelConsistent()
// *************************************************************************** 
namespace anrl {
	bool structureAndLabelConsistent(const xstructure& _xstr, 
			const string& label_input, 
      string& label_and_params_calculated){

    // Checks if the created structure is consistent with the label;
    // it is possible that the provided parameters elevate the structure
    // to a higher symmetry 

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::structureAndLabelConsistent():";

    xstructure xstr = _xstr; //copy

    // ---------------------------------------------------------------------------
    // determine label from structure (reverse process) 
    label_and_params_calculated = structure2anrl(xstr, true); // true=calculate sym

    // cannot do a strict string comparison of labels, symmetry analysis may
    // change origin (i.e., Wyckoff letters); need to check if labels are  
    // isopointal (check SG and Wyckoff multiplicities and site symmetries)

    vector<string> label_fields;
    aurostd::string2tokens(label_input,label_fields,"_");

    // ---------------------------------------------------------------------------
    // check space groups
    uint space_group_in_label = aurostd::string2utype<uint>(label_fields[2]);
    if(!compare::matchableSpaceGroups(xstr.space_group_ITC, space_group_in_label)){
      if(LDEBUG){
        cerr << function_name << " the calculated and label-designated space groups are incommensurate: "
          << "calculated=" << xstr.space_group_ITC << " vs "
          << "label= " << space_group_in_label << endl;
      }
      return false;
    }

    // ---------------------------------------------------------------------------
    // check Wyckoff positions

    // get Wyckoff information from label and format to compare
    vector<vector<string> > Wyckoff_fields = compare::convertANRLWyckoffString2GroupedPositions(label_input);
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions_label;
    compare::groupWyckoffPositionsFromGroupedString(space_group_in_label,
        xstr.setting_ITC,
        Wyckoff_fields,
        grouped_Wyckoff_positions_label);

    // get Wyckoff information from xstructure and format to compare
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions_structure;
    compare::groupWyckoffPositions(xstr, grouped_Wyckoff_positions_structure);
    string Wyckoff_string_structure = anrl::groupedWyckoffPosition2ANRLString(grouped_Wyckoff_positions_structure, true);

    if(!compare::matchableWyckoffPositions(grouped_Wyckoff_positions_label,
          grouped_Wyckoff_positions_structure,
          true)){ // same_species=true
      if(LDEBUG){
        cerr << function_name << " the calculated and label-designated Wyckoff positions are incommensurate: "
          << "calculated=" << Wyckoff_string_structure << " vs "
          << "label= " << label_input << endl;
      }
      return false;
    }

    // ---------------------------------------------------------------------------
    // all tests passed; the structure and label are commensurate
    return true;
  }
}

// *************************************************************************** 
// anrl::PrototypeANRL_Generator()
// *************************************************************************** 
namespace anrl {
  xstructure PrototypeANRL_Generator(string& label,
      string& parameters,
      deque<string> &vatomX,
      deque<double> &vvolumeX,
      ostream& logstream,
      bool silence_logger){

    ofstream FileMESSAGE;
    
    xstructure prototype = PrototypeANRL_Generator(label, 
        parameters, 
        vatomX, 
        vvolumeX, 
        FileMESSAGE, 
        logstream,
        silence_logger);

    return prototype;
  }
}

namespace anrl {
  xstructure PrototypeANRL_Generator(string& label,
      string& parameters,
      deque<string> &vatomX,
      deque<double> &vvolumeX,
      ofstream& FileMESSAGE,
      ostream& logstream,
      bool silence_logger){

		// Returns a ANRL prototype structure based on the label and internal 
    // degrees of freedom.
    // The function is generic and will build ANY prototype as long as: 
    // 1) the label and parameters are valid (the function has many checks) AND
    // 2) the structure is a crystal (i.e., built from Wyckoff positions)
    // A symbolic representations of the crystal can be returned in terms of: 
    // lattice variables: a, b, c, alpha, beta, gamma AND
    // Wyckoff variables: x, y, and z 

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::PrototypeANRL_Generator()";
    stringstream message;

    xstructure str;

    // ---------------------------------------------------------------------------
    // determine print mode 
    uint print_mode = _PROTO_GENERATOR_GEOMETRY_FILE_; // no equations
    if(XHOST.vflag_pflow.flag("PROTO::EQUATIONS_ONLY")) {
      print_mode = _PROTO_GENERATOR_EQUATIONS_ONLY_; // equations only
      str.symbolic_math_representation_only=TRUE;  //DX20180618 print symbolic math representation only
      message << "Printing the symbolic equations only";
    }
    else if(XHOST.vflag_pflow.flag("PROTO::ADD_EQUATIONS")) {
      print_mode = _PROTO_GENERATOR_GEOMETRY_FILE_AND_EQUATIONS_; // equations + parameters
      str.constrained_symmetry_calculation=TRUE;  //DX20180618 appends information to geometry file for calculation
      message << "Printing the geometry file and the symbolic equations";
    }
    else{
      message << "Printing geometry file only";
    }
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, silence_logger, _LOGGER_MESSAGE_);
    
    // ---------------------------------------------------------------------------
    // declare variables
    vector<string> vproto,vlabel; 
    vector<uint>   vproto_nspecies,vproto_natoms,vproto_spacegroup,vproto_nunderscores,vproto_nparameters; 
    vector<string> vproto_Pearson_symbol,vproto_params,vproto_Strukturbericht,vproto_prototype,vproto_dialect;

    // ---------------------------------------------------------------------------
    // load existing ANRL labels
    anrl::PrototypeANRL_LoadList(vproto,vlabel,vproto_nspecies,vproto_natoms,
        vproto_spacegroup,vproto_nunderscores,vproto_nparameters,vproto_Pearson_symbol,vproto_params,
        vproto_Strukturbericht,vproto_prototype,vproto_dialect);

    vector<string> tokens;
    string label_anrl="";
    string number_id = ""; //for predefined anrls //DX20191207
    string label_permutations=""; deque<uint> vpermutation;
    
    // ---------------------------------------------------------------------------
    // search for label_permutations
    aurostd::string2tokens(label,tokens,".");
    if(LDEBUG) { cerr << function_name << ": tokens.size()=" << tokens.size() << endl;}

    if(tokens.size()==0) { label_anrl=label; }
    if(tokens.size()==1) { label_anrl=tokens.at(0); }
    if(tokens.size()==2) { label_anrl=tokens.at(0); label_permutations=tokens.at(1); }
    for(uint i=0;i<label_permutations.size();i++) vpermutation.push_back(aurostd::mod(label_permutations.at(i)-65,32));

    // ---------------------------------------------------------------------------
    // check if preset suffix is included with label, e.g., A_hR2_166_c-001 
    if(aurostd::substring2bool(label_anrl,"-")){
      tokens.clear();
      aurostd::string2tokens(label_anrl,tokens,"-");
      if(tokens.size()==2){
        label_anrl = tokens[0];
        number_id = tokens[1];
      }
    }
    message << "The input label is " << label_anrl;
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, silence_logger, _LOGGER_MESSAGE_);
    
    if(number_id.size()){
      message << "The preset parameters " << number_id << " will be extracted (if they exist)"; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, silence_logger, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // search if label already exists in the library
    bool found=FALSE;
    uint ifound=0;
    for(uint i=0;i<vlabel.size()&&!found;i++) {
      if(vlabel.at(i)==label_anrl) {  // FIX
        found=TRUE;
        ifound=i;
        message << "This prototype label exists in the AFLOW library (part 1 or 2) label=" << label_anrl << "; index=" << ifound;
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, silence_logger, _LOGGER_MESSAGE_);
      }
    }

    // ---------------------------------------------------------------------------
    // not found, new label
    if(!found) {
      message << "This label does not currently exist in the AFLOW library label=" << label_anrl;
      message << " Consider adding it to the library.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, silence_logger, _LOGGER_MESSAGE_);
    }
    
    // -------------------------------------------------------------------------
    // check if using original anrl lattice parameter value when using the 
    // preset parameter functionality
    bool keep_anrl_lattice_parameter = false;
    bool scale_volume_by_species = true;
    if(parameters=="use_anrl_lattice_param"){
      keep_anrl_lattice_parameter=true;
      scale_volume_by_species = false;
      parameters=""; // clear the hack
    }
  
    // -------------------------------------------------------------------------
    // if no parameters given 
    if(parameters.size()==0){
      // -------------------------------------------------------------------------
      // extract values from library 
      if(found){
        parameters = extractANRLPrototypeParameterValues(label_anrl,
            number_id,
            vproto_params[ifound],
            keep_anrl_lattice_parameter);
      }
      message << "Extracted the following parameters (internal degrees of freedom) from the AFLOW parameters= " << parameters;
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, silence_logger, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // split label into fields 
    aurostd::string2tokens(label_anrl,tokens,"_");
      
    message << "The AFLOW label has been partitioned into " << tokens.size() << " fields : " << aurostd::joinWDelimiter(tokens, " ");
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, silence_logger, _LOGGER_MESSAGE_);

    if(tokens.size()<4){ 
      message << "Number of fields in label is too small, should be 4 or more: label" << label_anrl;
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
    }
    
    // ---------------------------------------------------------------------------
    // check stoichometry
    string composition = tokens[0];
    vector<uint> stoichiometry = ::composition2stoichiometry(composition);
    vector<uint> reduced_stoich; aurostd::reduceByGCD(stoichiometry, reduced_stoich);
    if(!compare::sameStoichiometry(stoichiometry,reduced_stoich)){
      message << "The input stoichiometry (first field in label=" << composition << ") is not reduced, it should be: " << aurostd::joinWDelimiter(reduced_stoich,":");
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
    }
    
    // ---------------------------------------------------------------------------
    // iniialize species variables
    vector<string> species;
    XATOM_SplitAlloySpecies(composition, species);

    for(uint i=0;i<species.size();i++) { // number of species
      str.num_each_type.push_back(0);str.comp_each_type.push_back(0.0);
      str.species.push_back("");str.species_pp.push_back("");str.species_pp_type.push_back("");str.species_pp_version.push_back("");
      str.species_pp_ZVAL.push_back(0.0);
      str.species_pp_vLDAU.push_back(deque<double>());
      str.species_volume.push_back(0.0);
      str.species_mass.push_back(0.0);
    }

    // ---------------------------------------------------------------------------
    // get Pearson symbol
    string Pearson_symbol = tokens[1];
    char lattice_type = Pearson_symbol[0];
    char lattice_centering = Pearson_symbol[1];
    uint number_of_atoms_conventional = aurostd::string2utype<uint>(Pearson_symbol.substr(2,Pearson_symbol.size()));

    if(LDEBUG){ cerr << function_name << " # atoms in conventional cell (from Pearson): " << number_of_atoms_conventional << endl; } 

    // ---------------------------------------------------------------------------
    // get space group number
    uint space_group_number = aurostd::string2utype<uint>(tokens[2]);

    if(space_group_number < 1 || space_group_number > 230){
      message << "The space group number is invalid; it must be between 1-230: spacegroup=" << space_group_number;
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
    }
    
    // ---------------------------------------------------------------------------
    // check if Pearson symbol and space group match
    string lattice_and_centering_from_Pearson = Pearson_symbol.substr(0,2);
    string lattice_and_centering_from_sg = SYM::spacegroup2latticeAndCentering(space_group_number);
    if(lattice_and_centering_from_Pearson != lattice_and_centering_from_sg){
      message << "Pearson symbol and space group number are incommensurate; the lattice centerings do not match:"; 
      message << "Pearson=" << Pearson_symbol << " (centering=" << lattice_and_centering_from_Pearson << ") vs ";
      message << "SG=" << space_group_number << "(centering=" << lattice_and_centering_from_sg << ")";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
    }
    
    // ---------------------------------------------------------------------------
    // get space group information
    string space_group_symbol = GetSpaceGroupName(space_group_number);
    char space_group_letter = space_group_symbol[0];

    // ---------------------------------------------------------------------------
    // get Wyckoff positions (the remaining fields after the first three)
    vector<string> Wyckoff_tokens; Wyckoff_tokens.insert(Wyckoff_tokens.end(), tokens.begin()+3, tokens.end()); // get Wyckoff tokens
    
    // ---------------------------------------------------------------------------
    // check if number of Wyckoff positions match the number of species
    if(Wyckoff_tokens.size() != species.size()){
      message << "The number of species does not match the number of Wyckoff species: # species=" << species.size() << " vs ";
      message << "# Wyckoff species=" << Wyckoff_tokens.size() << " (input label=" << label << ")" << endl;
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
    }
    
    if(LDEBUG){
      cerr << function_name << " the Wyckoff sequences associated with each species are:" << endl;
      for(uint i=0;i<species.size();i++){ cerr << species[i] << ": " << Wyckoff_tokens[i] << endl; }
    }

    // ---------------------------------------------------------------------------
    // initialize ITC space group/Wyckoff position object
    uint setting=SG_SETTING_ANRL;
    
    vector<wyckoffsite_ITC> wyckoff_sites_ITC = anrl::getWyckoffSitesFromANRL(Wyckoff_tokens, 
        species, 
        space_group_number,
        setting);
    vector<uint> number_of_each_type = SYM::numberEachTypeFromWyckoff(wyckoff_sites_ITC);
    
    if(LDEBUG){
      cerr << function_name << " the Wyckoff positions are:" << endl;
      print(wyckoff_sites_ITC);
    }
    
    // ---------------------------------------------------------------------------
    // check Wyckoff positions and stoichiometry
    //vector<uint> Wyckoff_reduced_stoich = compare::gcdStoich(number_of_each_type);
    vector<uint> Wyckoff_reduced_stoich; aurostd::reduceByGCD(number_of_each_type, Wyckoff_reduced_stoich);
    if(!compare::sameStoichiometry(stoichiometry,Wyckoff_reduced_stoich)){
      message << "The input composition and Wyckoff positions yield different stoichiometries: composition=" << aurostd::joinWDelimiter(stoichiometry,":");  
      message << ", Wyckoff=" << aurostd::joinWDelimiter(Wyckoff_reduced_stoich,":") << endl;
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
    }
    
    // ---------------------------------------------------------------------------
    // check Wyckoff multiplicity and number of atoms in the conventional cell
    // (should make this into a function)
    uint Wyckoff_multiplicity_sum = 0;
    for(uint i=0;i<wyckoff_sites_ITC.size();i++){ Wyckoff_multiplicity_sum+=wyckoff_sites_ITC[i].multiplicity; }
    if(Wyckoff_multiplicity_sum!=number_of_atoms_conventional){
      message << "The sum of the Wyckoff multiplicity does not add up to the number of atoms in the conventional cell (from Pearson symbol); bad prototype label: ";
      message << " Wyckoff multiplicity sum = " << Wyckoff_multiplicity_sum;
      message << ", Pearson symbol = " << Pearson_symbol;
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // determine parameters
    vector<string> parameter_list, lattice_parameter_list, Wyckoff_parameter_list;
    vector<double> parameter_values, lattice_parameter_values, Wyckoff_parameter_values;
    
    // -------------------------------------------------------------------------
    // lattice parameters
    lattice_parameter_list = getANRLLatticeParameterString(lattice_type);
    
    // ---------------------------------------------------------------------------
    // reorder Wyckoff positions alphabetically by Wyckoff letter, then by species
		vector<wyckoffsite_ITC> ordered_Wyckoff_sites_ITC = wyckoff_sites_ITC;
    std::sort(ordered_Wyckoff_sites_ITC.begin(), ordered_Wyckoff_sites_ITC.end(), sortWyckoffByLetter); 
    
    if(LDEBUG){ 
      for(uint i=0;i<ordered_Wyckoff_sites_ITC.size();i++){
        cerr << function_name << "::Ordered Wyckoff site: " << ordered_Wyckoff_sites_ITC[i] << endl;
      }
    }
    
    // ---------------------------------------------------------------------------
    // determine degrees of freedom in Wyckoff positions 
    Wyckoff_parameter_list = determineWyckoffVariables(ordered_Wyckoff_sites_ITC); 
    
    
    // ---------------------------------------------------------------------------
    // combine parameter vectors
    parameter_list = lattice_parameter_list; 
    parameter_list.insert(parameter_list.end(), Wyckoff_parameter_list.begin(), Wyckoff_parameter_list.end());
    if(LDEBUG){ cerr << function_name << " parameters=" << aurostd::joinWDelimiter(parameter_list,",") << endl; }
    
    // -------------------------------------------------------------------------
    // if no parameters provided and more than one is needed
    if(parameters.size()==0 && parameter_list.size()!=1){
      message << "No parameters provided. Since this is a new prototype label with more than one degree of freedom,";
      message << "you must add parameter values with --params=..." << endl;
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
    }
    else if(parameters.size()==0 && parameter_list.size()==1){
      parameters = "1.0";
    }

    // ---------------------------------------------------------------------------
    // partition in parameter values
    vector<string> vparameters_temp;
    aurostd::string2tokens(parameters,vparameters_temp,",");
    vector<double> vparameters = aurostd::vectorstring2vectordouble(vparameters_temp);
    
    // ---------------------------------------------------------------------------
    // check for automatic volume scaling (i.e., first parameter is negative)
    if(vparameters[0]<=0.0){ //CO20181226 forget signbit, also include 0
      vparameters[0]=1.0; //fix
      vparameters_temp[0]="1.0"; //fix
      parameters=aurostd::joinWDelimiter(vparameters_temp,",");
    }

    if(vparameters.size() != parameter_list.size()){
      message << "The number of input parameters does not match the number required by the lattice and Wyckoff positions: ";
      message << "input parameters=" << parameters << " vs ";
      message << "parameter_list=" << aurostd::joinWDelimiter(parameter_list,",") << endl;
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_NUMBER_);
    }
    
    // ---------------------------------------------------------------------------
    // populate degree of freedom values 
    lattice_parameter_values.insert(lattice_parameter_values.end(), vparameters.begin(), vparameters.begin()+lattice_parameter_list.size());
    Wyckoff_parameter_values.insert(Wyckoff_parameter_values.end(), vparameters.begin()+lattice_parameter_list.size(), vparameters.end());
    
    anrl::applyWyckoffValues(Wyckoff_parameter_values, ordered_Wyckoff_sites_ITC); 
    
    // ---------------------------------------------------------------------------
    // check to ensure no duplicate 
    if(containsDuplicateWyckoffCoordinate(ordered_Wyckoff_sites_ITC,true)){
      message << "Contains duplicate Wyckoff letters with the same degrees of freedom.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_INPUT_ILLEGAL_);
    }

    // ---------------------------------------------------------------------------
    // generate lattice (in ANRL convention) based on symmetry and lattice
    // parameter values
    
    // primitive (0)
    xmatrix<double> lattice_primitive = getLattice(
        lattice_and_centering_from_Pearson, 
        space_group_letter, 
        lattice_parameter_values, 
        0);
    // conventional (1)
    xmatrix<double> lattice_conventional = getLattice(
        lattice_and_centering_from_Pearson, 
        space_group_letter, 
        lattice_parameter_values, 
        1);
   
    // ---------------------------------------------------------------------------
    // generate atoms based Wyckoff equations and Wyckoff parameter values
    deque<_atom> atoms_conventional_cell = getAtomsFromWyckoff(ordered_Wyckoff_sites_ITC,lattice_conventional);

    deque<_atom> atoms_primitive_cell;
    // special case: if using the rhombohedral setting, then the Wyckoff positions 
    // are already wrt to the primitive cell; no need to perform conversion 
    if(setting==SG_SETTING_ANRL && lattice_centering=='R'){
      atoms_primitive_cell = atoms_conventional_cell;
    }
    // generic case: convert conventional to primitive
    else{
      atoms_primitive_cell = foldAtomsInCell(atoms_conventional_cell,
          lattice_conventional,
          lattice_primitive,
          false,
          1e-6,
          false);
    }

    if(LDEBUG){
      cerr << function_name << " atoms in the primitive cell:" << endl;
      for(uint i=0;i<atoms_primitive_cell.size();i++){
        cerr << atoms_primitive_cell[i] << " " << atoms_primitive_cell[i].name << endl;
      }
    }
    
    // ---------------------------------------------------------------------------
    // check ratio between conventional and primitive atoms
    if(atoms_conventional_cell.size()%atoms_primitive_cell.size()!=0){
      message << "The ratio of atoms between the conventional cell and primitive cell is not an integer; check the tolerance: #conventional=" << atoms_conventional_cell.size() << " vs #primitive=" << atoms_primitive_cell.size();
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_RUNTIME_ERROR_);
    }
    uint ratio_calculated = atoms_conventional_cell.size()/atoms_primitive_cell.size();
    uint ratio_conventional2primitive = LATTICE::Conventional2PrimitiveRatio(lattice_centering);
    if(!aurostd::isequal(ratio_calculated, ratio_conventional2primitive)){
      if(!(ratio_calculated==1 && lattice_centering=='R' && setting==SG_SETTING_ANRL)){
        message << "The calculated ratio and the expected conventional2primtive ratio do not match: calculated=" << ratio_calculated << " vs expected=" << ratio_conventional2primitive;
        throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_RUNTIME_ERROR_);
      }
    }

    // ---------------------------------------------------------------------------
    // create xstructure
    str.iomode=IOVASP_AUTO;
    str.title=label+" params="+parameters+" SG="+aurostd::utype2string(space_group_number)+DOI_ANRL; //CO20190520
    str.scale=1.0;
    str.lattice = lattice_primitive;
    str.atoms = atoms_primitive_cell;

    // ---------------------------------------------------------------------------
    // add ANRL info to xstructure
    str.num_parameters = parameter_list.size(); 
    str.num_lattice_parameters = lattice_parameter_list.size(); 
    str.prototype_parameter_list = parameter_list; 
    str.prototype_parameter_values = parameter_values;
    str.setting_ITC=setting;
   
    // ---------------------------------------------------------------------------
    // convert RHL to HEX setting ([--hex] option)
    if(XHOST.vflag_pflow.flag("PROTO::HEX") && lattice_centering == 'R') {
      vector<double> vparameters; aurostd::string2tokens(parameters,vparameters,",");
      uint i=0;
      double a=vparameters.at(i++);
      double covera=vparameters.at(i++);
      double c=covera*a;
      str=rhl2hex(str,a,c);
    }

    for(uint iat=0;iat<str.atoms.size();iat++) {
      str.atoms.at(iat).name_is_given=TRUE;
      str.atoms.at(iat).number=iat;//iat;    // reference position for convasp
      str.atoms.at(iat).basis=iat;//iat;     // position in the basis
      if(print_mode!=_PROTO_GENERATOR_EQUATIONS_ONLY_){ //equations only //DX20180618
        str.atoms.at(iat).cpos=F2C(str.lattice,str.atoms.at(iat).fpos);
      }
      str.num_each_type.at(str.atoms.at(iat).type)++;
      //     str.comp_each_type.at(str.atoms.at(iat).type)+=1.0; inside code
      str.species.at(str.atoms.at(iat).type)=str.atoms.at(iat).name;	
    }
    
    // ---------------------------------------------------------------------------
    // symbolic representation of prototypes
    if(print_mode == _PROTO_GENERATOR_EQUATIONS_ONLY_ || print_mode == _PROTO_GENERATOR_GEOMETRY_FILE_AND_EQUATIONS_){ 
#ifdef COMPILE_SYMBOLIC
      // ---------------------------------------------------------------------------
      // get symbolic lattice
      Symbolic lattice_symbolic = SymbolicANRLPrimitiveLattices(lattice_and_centering_from_Pearson, space_group_letter);

      // ---------------------------------------------------------------------------
      // re-sort to alphabetic/type ordering 
      //vector<wyckoffsite_ITC> Wyckoff_sites_ordered_by_type = ordered_Wyckoff_sites_ITC;
      //std::sort(Wyckoff_sites_ordered_by_type.begin(), Wyckoff_sites_ordered_by_type.end(), sortWyckoffByType); 

      vector<SymbolicWyckoffSite> Wyckoff_sites_symbolic;
      //for(uint i=0;i<Wyckoff_sites_ordered_by_type.size();i++){
        //Wyckoff_sites_symbolic.push_back(initializeSymbolicWyckoffSite(Wyckoff_sites_ordered_by_type[i]));
      for(uint i=0;i<ordered_Wyckoff_sites_ITC.size();i++){
        Wyckoff_sites_symbolic.push_back(initializeSymbolicWyckoffSite(ordered_Wyckoff_sites_ITC[i]));
      }

      // ---------------------------------------------------------------------------
      // convert to symbolic equations 
      //vector<Symbolic> all_symbolic_equations;
      //for(uint i=0;i<Wyckoff_sites_ordered_by_type.size();i++){
      //  vector<Symbolic> symbolic_equation = equations2SymbolicEquations(Wyckoff_sites_ordered_by_type[i].equations);
      //  all_symbolic_equations.insert(all_symbolic_equations.end(), symbolic_equation.begin(), symbolic_equation.end());
      //}

      // ---------------------------------------------------------------------------
      // transform to ANRL primitive cell 
      //vector<Symbolic> reduced_symbolic_equation = convertEquations2FractionalEquations(lattice_and_centering_from_Pearson, lattice_symbolic, all_symbolic_equations);
      for(uint i=0;i<Wyckoff_sites_symbolic.size();i++){
        Wyckoff_sites_symbolic[i].equations = convertEquations2FractionalEquations(lattice_and_centering_from_Pearson, lattice_symbolic, Wyckoff_sites_symbolic[i].equations);
      }

      substituteVariableWithParameterDesignation(Wyckoff_sites_symbolic); 
      vector<Symbolic> symbolic_equations;
      for(uint i=0;i<Wyckoff_sites_symbolic.size();i++){
        symbolic_equations.insert(symbolic_equations.end(), Wyckoff_sites_symbolic[i].equations.begin(), Wyckoff_sites_symbolic[i].equations.end());
      }

      str.symbolic_math_lattice = symbolic::matrix2VectorVectorString(lattice_symbolic);
      addSymbolicEquation2Atoms(symbolic_equations, str.atoms);
      //addSymbolicEquation2Atoms(symbolic_equations, atoms_primitive_cell);
#else
      // ---------------------------------------------------------------------------
      // if the SymbolicC++ code is not compiled
      message << "The SymbolicC++ source code has not been compiled, symbolic equations cannot be printed";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_RUNTIME_ERROR_);
#endif
    }
    
    // ---------------------------------------------------------------------------
    // DONE
    if(print_mode!=_PROTO_GENERATOR_EQUATIONS_ONLY_){ //equations only //DX20180618
      xvector<double> data(6);
      data=Getabc_angles(str.lattice,DEGREES);
      str.a=data[1];str.b=data[2];str.c=data[3];str.alpha=data[4];str.beta=data[5];str.gamma=data[6];
      clear(str.origin);
    }
    //  if(vpflow.flag("STDPRIMCELL")) {cout << "EUREKA"<< endl;} //cout << GetStandardPrimitive(xstructure(cin,IOAFLOW_AUTO));

    // ---------------------------------------------------------------------------
    // NOW PLAY WITH PERMUTATIONS and ATOMX
    if(vpermutation.size()>0 || vatomX.size()>0) {
      if(LDEBUG) { cerr << function_name << " PERMUTATIONS" << endl;}
      if(LDEBUG) { cerr << function_name << " vpermutation.size()=" << vpermutation.size() << endl;}
      if(LDEBUG) { cerr << function_name << " vpermutation ="; for(uint i=0;i<vpermutation.size();i++) {cerr << " " << vpermutation.at(i);} cerr << endl;}     
      if(LDEBUG) { cerr << function_name << " ATOMX" << endl;}
      if(LDEBUG) { cerr << function_name << " vatomX.size()=" << vatomX.size() << endl;}
      if(LDEBUG) { cerr << function_name << " vatomX ="; for(uint i=0;i<vatomX.size();i++) {cerr << " " << vatomX.at(i);} cerr << endl;}
      if(print_mode!=_PROTO_GENERATOR_EQUATIONS_ONLY_){ //equations only //DX20180618
        std::deque<_atom> atoms;
        atoms=str.atoms;
        // STRIP ALL ATOMS
        while(str.atoms.size()>0) { str.RemoveAtom(0); }
        // ADD MODIFIED ATOMS
        for(uint i=0;i<atoms.size();i++) {
          uint type=atoms.at(i).type;
          if(vpermutation.size()>0)  { atoms.at(i).type=vpermutation.at(type); }  // PERMUTATIONS 
          if(vpermutation.size()>0 || vatomX.size()>0) { atoms.at(i).name=vatomX.at(atoms.at(i).type); }  // PERMUTATIONS AND ATOMX
          //	atoms.at(i).name=aurostd::mod(label_permutations.at(type)-65,32)+65;
          str.AddAtom(atoms.at(i));
          // DX20181205 - Volume scaling by atomic species - START
          // if a=1.0 for prototype (i.e., no scaling factor), use atomic species to get volume
          if(scale_volume_by_species==true){
            double volume=0.0;
            for(uint i=0;i<str.num_each_type.size();i++) {
              for(uint j=0;j<(uint)str.num_each_type[i];j++){
                volume+=vvolumeX[i];
                if(LDEBUG) { cerr << function_name << " volume=" << volume << "  (" << vvolumeX[i] << ")" << endl; }
              }
            }
            //[CO20190205 - OBSOLETE]str.scale=std::pow((double) (abs(volume)/det(str.lattice)),(double) 1.0/3.0);
            str.SetVolume(volume);  //CO20190205 - more robust
            str.neg_scale=TRUE;
          }
          //DX20181205 - Volume scaling by atomic species - END
        }
      }
      str.SpeciesPutAlphabetic();
    }

		// ---------------------------------------------------------------------------
		// fix title of geometry file 
    aurostd::StringSubst(str.title,label,aurostd::joinWDelimiter(str.species,"")+"/"+label);  //vlabel.at(ifound) //use label as we want permutations too //CO20181216
    if(scale_volume_by_species){aurostd::StringSubst(str.title," params=1.0"," params=-1");} //CO20181216
    
		// ---------------------------------------------------------------------------
    // if this is a new prototype (i.e., not in library), we should check the 
    // symmetry; it is possible that the provided parameters elevate the structure
    // to a higher symmetry 
    if(!found){
			string updated_label_and_params = "";
      if(!structureAndLabelConsistent(str, label_anrl, updated_label_and_params)){ 
        // if changes symmetry, give the appropriate label
        message << "The structure has a higher symmetry than indicated by the label. ";
        message << "The correct label and parameters for this structure are:" << endl;
        message << updated_label_and_params << endl;
        message << "Please feed this label and set of parameters into the prototype generator.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_RUNTIME_ERROR_);
      }
    } 

    return str;
  }
}


#ifdef COMPILE_SYMBOLIC
// *************************************************************************** 
// *************************************************************************** 
// Symbolic math functions leveraging the SymbolicC++ source code
// 
// SymbolicC++ is written by Y. Hardy, W.-H. Steeb, and T.K. Shi 
// SymbolicC++ is open-source via the GNU-GPL licensce (version 2)
// see http://symboliccpp.sourceforge.net/ for more information
//
// AFLOW implementation:
// The SymbolicC++ header files are included with AFLOW (compiled locally)
// Below are generic helper functions, e.g., convert string to Symbolic 
// representation and ANRL prototype specific functions
//
// Generic functions are in the "symbolic" namespace, and 
// ANRL specific functions are in the "anrl" namespace
//
// For more information on the AFLOW implementation, contact 
// David Hicks (DX) david.hicks@duke.edu
// *************************************************************************** 
// *************************************************************************** 

// ---------------------------------------------------------------------------
// GENERIC SYMBOLIC MATH FUNCTIONS
// (if more functions are added, perhaps make a new file: aflow_symbolic.cpp)

// *************************************************************************** 
// SymbolicWyckoffSite initializeSymbolicWyckoffSite() 
// *************************************************************************** 
SymbolicWyckoffSite initializeSymbolicWyckoffSite(const wyckoffsite_ITC& Wyckoff){

  // initialize symbolic Wyckoff site struct 

  SymbolicWyckoffSite Wyckoff_symbolic;
  Wyckoff_symbolic.coord = Wyckoff.coord;
  Wyckoff_symbolic.index = Wyckoff.index;
  Wyckoff_symbolic.type = Wyckoff.type;
  Wyckoff_symbolic.wyckoffSymbol = Wyckoff.wyckoffSymbol;
  Wyckoff_symbolic.letter = Wyckoff.letter;
  Wyckoff_symbolic.site_symmetry = Wyckoff.site_symmetry;
  Wyckoff_symbolic.multiplicity = Wyckoff.multiplicity;
  Wyckoff_symbolic.site_occupation = Wyckoff.site_occupation;
  Wyckoff_symbolic.equations = anrl::equations2SymbolicEquations(Wyckoff.equations);
  Wyckoff_symbolic.parameter_index = Wyckoff.parameter_index;
  return Wyckoff_symbolic;
}

// *************************************************************************** 
// substituteVariableWithParameterDesignation() 
// *************************************************************************** 
// vector<SymbolicWyckoffSite>
void substituteVariableWithParameterDesignation(vector<SymbolicWyckoffSite>& Wyckoff_sites_symbolic){
  for(uint i=0;i<Wyckoff_sites_symbolic.size();i++){
    substituteVariableWithParameterDesignation(Wyckoff_sites_symbolic[i]);
  }
}

// SymbolicWyckoffSite
void substituteVariableWithParameterDesignation(SymbolicWyckoffSite& Wyckoff_symbolic){
  for(uint i=0;i<Wyckoff_symbolic.equations.size();i++){
    Wyckoff_symbolic.equations[i]=Wyckoff_symbolic.equations[i].subst("x","x"+aurostd::utype2string<uint>(Wyckoff_symbolic.parameter_index));
    Wyckoff_symbolic.equations[i]=Wyckoff_symbolic.equations[i].subst("y","y"+aurostd::utype2string<uint>(Wyckoff_symbolic.parameter_index));
    Wyckoff_symbolic.equations[i]=Wyckoff_symbolic.equations[i].subst("z","z"+aurostd::utype2string<uint>(Wyckoff_symbolic.parameter_index));
  }
}

// *************************************************************************** 
// symbolic::string2symbolic
// *************************************************************************** 
namespace symbolic {
  Symbolic string2symbolic(const string& str){
    
    // Convert a string into symbolic math notation
    // The SYMBOLICC++ library cannot convert a string into a symbol so
    // we must do it ourselves
    // NOTE: this is preliminary, if this is to be used beyond Wyckoff 
    // positions, we need to add more operators/functions 
    // (e.g., sin, cos, exponentials, etc.)

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "symbolic::string2symbolic():";
    stringstream message;

    // ---------------------------------------------------------------------------
    // break up string into doubles and characters
    vector<SYM::sdouble> sdouble_temp = SYM::simplify(str);
    
    // ---------------------------------------------------------------------------
    // loop over characters and build symbolic expression 
    Symbolic out;
    for(uint i=0;i<sdouble_temp.size();i++){
      if(LDEBUG){
        cerr << function_name << " sdouble_temp: (dbl) " << sdouble_temp[i].dbl 
          << " (chr) " << sdouble_temp[i].chr << endl;
      }
      // ---------------------------------------------------------------------------
      // if contains char 
      if(sdouble_temp[i].chr != '\0'){
        stringstream ss_char; ss_char << sdouble_temp[i].chr;
        out += Symbolic(sdouble_temp[i].dbl)*Symbolic(ss_char.str()); // '*' is a big assumption, may need to generalize in the future
      }
      // ---------------------------------------------------------------------------
      // else contains only double
      else{
        out += Symbolic(sdouble_temp[i].dbl);
      }
    }

    if(LDEBUG){ cerr << function_name << " symbolic expression: out=" << out << endl; }
    return out;
  }
}

// *************************************************************************** 
// symbolic::isEqualSymbolic()
// *************************************************************************** 
namespace symbolic {
  bool isEqual(const Symbolic& a, const Symbolic& b){
    
    // need to be able to handle (1e-9)*x+(1e-6)*y+0
    // right now this relies on SymbolicC++'s implementation of ==
    // this may not be robust enough 

    bool VERBOSE=FALSE; // VERBOSE INSTEAD OF LDEBUG SINCE A NESTED FUNCTION

    Symbolic diff = a - b;

    if(VERBOSE){ 
      string function_name = XPID + "symbolic::isEqual():";
      cerr << function_name << " a-b=" << diff << endl;
    }

    return (diff==_SYMBOLIC_ZERO_);
  }
}

// *************************************************************************** 
// symbolic::isEqualSymbolicVector()
// *************************************************************************** 
namespace symbolic {
  bool isEqualVector(const Symbolic& a_vec, const Symbolic& b_vec){
    
    // perhaps check type (number, vector, matrix?) and have if statements?

    bool VERBOSE=FALSE; // VERBOSE INSTEAD OF LDEBUG SINCE A NESTED FUNCTION

    if(VERBOSE){ 
      string function_name = XPID + "symbolic::isEqualVector():";
      cerr << function_name << " a_vec-b_vec=" << (a_vec-b_vec) << endl;
    }

    for(uint i=0;i<3;i++){
      if(!isEqual(a_vec(i),b_vec(i))){ return false; }
    }

    return true;
  }
}

// *************************************************************************** 
// symbolic::matrix2VectorVectorString() 
// *************************************************************************** 
namespace symbolic {
  vector<vector<string> > matrix2VectorVectorString(const Symbolic& lattice){

    // convert symbolic matrix to vector<vector<string> >
    // perhaps check if type is matrix?

    vector<vector<string> > vvstring;

    for(uint i=0;i<3;i++){
      vector<string> row;
      for(uint j=0;j<3;j++){
        stringstream ss_component; ss_component << lattice.row(i)(j);
				row.push_back(ss_component.str());
      }
      vvstring.push_back(row);
    }

    return vvstring;
  }
}

// *************************************************************************** 
// symbolic::BringInCell() 
// *************************************************************************** 
namespace symbolic {
  Symbolic BringInCell(const Symbolic& vec_in, double tolerance, double upper_bound, double lower_bound){
    
    // bring symbolic math inside the cell
    // it is impossible to know what the variable (x, y, or z) will be 
    // to truly bring the coordinate in the cell, but this function brings 
    // the constant term in cell
    // uses the Symbolic.coeff(<variable>,<order>) function, where
    // <variable>: is the variable of interest
    // <order>: is the order/degree of the variable
    // TRICK: define a constant and provide order 1, i.e.,
    // Symbolic constant = Symbolic(1); equation.coeff(constant,1);
    // tune the upper and lower bound (e.g, 0 to 1 or -0.5 to 0.5)
    // NOTE: preference to lower bound (e.g., 0) vs upper bound (e.g., 1)

    Symbolic vec_out = vec_in;
    
    // ---------------------------------------------------------------------------
    // define a constant so we can determine the "coefficients of the constant"
    Symbolic constant = Symbolic(1); // 1 is a great constant

    for(uint i=0;i<3;i++){
      // ---------------------------------------------------------------------------
      // bring in cell based on constant term
      while(double(vec_out(i).coeff(constant,1))-upper_bound>=-tolerance){ vec_out(i) -= 1; }
      while(double(vec_out(i).coeff(constant,1)-lower_bound)<-tolerance){ vec_out(i) += 1; }
    }
    return vec_out;
  }
}


// ---------------------------------------------------------------------------
// ANRL SPECIFIC SYMBOLIC MATH FUNCTIONS

// *************************************************************************** 
// anrl::SymbolicANRLPrimitiveLattices()
// *************************************************************************** 
namespace anrl {
  Symbolic SymbolicANRLPrimitiveLattices(const string& lattice_and_centering, const char& space_group_letter){ 
    
    // Grab symbolic representation of primitive lattice in the ANRL convention.

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::SymbolicANRLPrimitiveLattices():";

    if(LDEBUG){
      cerr << function_name << " lattice and centering: " << lattice_and_centering << endl;
      cerr << function_name << " first character of space group symbol: " << space_group_letter << endl;
    }

    Symbolic lattice("L",3,3);

    // ---------------------------------------------------------------------------
    // create symbolic list of characters
    Symbolic a("a");
    Symbolic b("b");
    Symbolic c("c");
    Symbolic cx("cx");
    Symbolic cy("cy");
    Symbolic cz("cz");
    Symbolic beta("beta");
    Symbolic gamma("gamma");

    // ---------------------------------------------------------------------------
    // triclinic (aP)
    if(lattice_and_centering == "aP"){
	    lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
	              (b*cos(gamma), b*sin(gamma), _SYMBOLIC_ZERO_),
	              (cx, cy, cz));
    }

    // ---------------------------------------------------------------------------
    // simple monoclinic (mP)
    if(lattice_and_centering == "mP"){
      lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
                 (_SYMBOLIC_ZERO_, b, _SYMBOLIC_ZERO_),
                 (c*cos(beta), _SYMBOLIC_ZERO_, c*sin(beta)));
    }

    // ---------------------------------------------------------------------------
    // base-centered monoclinic (mC)
    if(lattice_and_centering == "mC"){
      lattice = ((1/2*a, -(1.0/2.0)*b, _SYMBOLIC_ZERO_),
	               (1.0/2.0*a, (1.0/2.0)*b, _SYMBOLIC_ZERO_),
	               (c*cos(beta), _SYMBOLIC_ZERO_ , c*sin(beta)));
    }

    // ---------------------------------------------------------------------------
    // simple orthorhombic (oP)
    if(lattice_and_centering == "oP"){
	    lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
	               (_SYMBOLIC_ZERO_, b, _SYMBOLIC_ZERO_),
                 (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }

    // ---------------------------------------------------------------------------
    // base-centered orthorhombic (oC)
    if(lattice_and_centering == "oC"){
	    if(space_group_letter == 'A'){
	      lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
	                 (_SYMBOLIC_ZERO_, (1.0/2.0)*b, -(1.0/2.0)*c),
	                 (_SYMBOLIC_ZERO_, (1.0/2.0)*b, (1.0/2.0)*c));
      }
	    if(space_group_letter == 'C'){
	      lattice = (((1.0/2.0)*a, -(1.0/2.0)*b, _SYMBOLIC_ZERO_),
	                 ((1.0/2.0)*a, (1.0/2.0)*b, _SYMBOLIC_ZERO_),
	                 (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
      }
    }

    // ---------------------------------------------------------------------------
    // body-centered orthorhombic (oI)
    if(lattice_and_centering == "oI"){
	    lattice = ((-(1.0/2.0)*a, (1.0/2.0)*b, (1.0/2.0)*c),
	               ((1.0/2.0)*a, -(1.0/2.0)*b, (1.0/2.0)*c),
	               ((1.0/2.0)*a, (1.0/2.0)*b, -(1.0/2.0)*c));
    }

    // ---------------------------------------------------------------------------
    // face-centered orthorhombic (oF)
    if(lattice_and_centering == "oF"){
	    lattice = ((_SYMBOLIC_ZERO_, (1.0/2.0)*b, (1.0/2.0)*c),
	               ((1.0/2.0)*a, _SYMBOLIC_ZERO_, (1.0/2.0)*c),
	               ((1.0/2.0)*a, (1.0/2.0)*b, _SYMBOLIC_ZERO_));
    }

    // ---------------------------------------------------------------------------
    // simple tetragonal (tP)
    if(lattice_and_centering == "tP"){
	    lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
	               (_SYMBOLIC_ZERO_, a, _SYMBOLIC_ZERO_),
	               (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }

    // ---------------------------------------------------------------------------
    // body-centered tetragonal (tI)
    if(lattice_and_centering == "tI"){
	    lattice = ((-(1.0/2.0)*a, (1.0/2.0)*a, (1.0/2.0)*c),
	               ((1.0/2.0)*a, -(1.0/2.0)*a, (1.0/2.0)*c),
	               ((1.0/2.0)*a, (1.0/2.0)*a, -(1.0/2.0)*c));
    }

    // ---------------------------------------------------------------------------
    // hexagonal/trigonal (hP)
    if(lattice_and_centering == "hP"){
	    lattice = (((1.0/2.0)*a, -(sqrt(3.0)/2.0)*a, _SYMBOLIC_ZERO_),
	               ((1.0/2.0)*a, (sqrt(3.0)/2.0)*a, _SYMBOLIC_ZERO_),
	               (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }

    // ---------------------------------------------------------------------------
    // rhombohedral (hR)
    if(lattice_and_centering == "hR"){
	    lattice = (((1.0/2.0)*a, -(1.0/(2.0*sqrt(3.0)))*a, (1.0/3.0)*c),
	               (_SYMBOLIC_ZERO_, (1.0/sqrt(3.0))*a, (1.0/3.0)*c),
	               (-(1.0/2.0)*a, -(1.0/(2.0*sqrt(3.0)))*a, (1.0/3.0)*c));
    }

    // ---------------------------------------------------------------------------
    // simple cubic (cP)
    if(lattice_and_centering == "cP"){
	    lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
	               (_SYMBOLIC_ZERO_, a, _SYMBOLIC_ZERO_),
	               (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, a));
    }

    // ---------------------------------------------------------------------------
    // body-centered cubic (cI)
    if(lattice_and_centering == "cI"){
	    lattice = ((-(1.0/2.0)*a, (1.0/2.0)*a, (1.0/2.0)*a),
	               ((1.0/2.0)*a, -(1.0/2.0)*a, (1.0/2.0)*a),
	               ((1.0/2.0)*a, (1.0/2.0)*a, -(1.0/2.0)*a));
    }
    
    // ---------------------------------------------------------------------------
    // face-centered cubic (cF)
    if(lattice_and_centering == "cF"){
	    lattice = ((_SYMBOLIC_ZERO_, (1.0/2.0)*a, (1.0/2.0)*a),
	               ((1.0/2.0)*a, _SYMBOLIC_ZERO_, (1.0/2.0)*a),
	               ((1.0/2.0)*a, (1.0/2.0)*a, _SYMBOLIC_ZERO_));
    }
		
		if(LDEBUG){
      cerr << function_name << " lattice: " << lattice << endl;
    }

    return lattice;
  }
}

// *************************************************************************** 
// anrl::equations2SymbolicEquations()
// *************************************************************************** 
namespace anrl {
  vector<Symbolic> equations2SymbolicEquations(const vector<vector<string> >& equations){
    
    // Convert equations (vector<vector<string> >) to symbolic notation.

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::equations2SymbolicEquations()";
    stringstream message;

    vector<Symbolic> symbolic_equations;

    for(uint i=0;i<equations.size();i++){
      Symbolic position("pos", 3);
      if(equations[i].size()!=3){
        message << "Equation " << i << " does not have 3 coordinates (problem with ITC library coordinates): " << aurostd::joinWDelimiter(equations[i],",");
        throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_GENERIC_ERROR_);
      }
      for(uint j=0;j<equations[i].size();j++){
        if(equations[i][j] == "0"){
          position(j) = _SYMBOLIC_ZERO_;
        }
        else{
          position(j) = symbolic::string2symbolic(equations[i][j]);
        }
      }
      symbolic_equations.push_back(position);
    }

    if(LDEBUG){
      for(uint i=0;i<symbolic_equations.size();i++){
        cerr << function_name << " equations (string): " << aurostd::joinWDelimiter(equations[i],",") << " --> equations (Symbolic): " << symbolic_equations[i] << endl;
			}
    }

    return symbolic_equations;
  }
}

// *************************************************************************** 
// anrl::cartesian2lattice()
// *************************************************************************** 
namespace anrl {
  Symbolic cartesian2lattice(const Symbolic& lattice, const Symbolic& cartesian_coordinate){ 
 
    // converts Cartesian coordinates to lattice coordinates
    // this is the symbolic math equivalent to C2F()

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::cartesian2lattice():";

    // ---------------------------------------------------------------------------
    // calculate volume (symbolic) 
    Symbolic volume = lattice.row(0)|((lattice.row(1)%(lattice.row(2))).transpose()); // | is dot product, % is cross product
    if(LDEBUG){
      cerr << function_name << " symbolic volume=" << volume << endl;
    }

    // ---------------------------------------------------------------------------
    // determine transformation matrix (vectors) 
    Symbolic b1 = (lattice.row(1)%(lattice.row(2)))/volume;
    Symbolic b2 = (lattice.row(2)%(lattice.row(0)))/volume;
    Symbolic b3 = (lattice.row(0)%(lattice.row(1)))/volume;
    
		if(LDEBUG){
      cerr << function_name << " transformation vectors b1=" << b1 << endl;
      cerr << function_name << " transformation vectors b2=" << b2 << endl;
      cerr << function_name << " transformation vectors b3=" << b3 << endl;
    }
    
    // ---------------------------------------------------------------------------
    // convert to lattice coordinate 
		Symbolic lattice_coordinate("latt_coord",3);
    lattice_coordinate(0)=(cartesian_coordinate|b1).simplify();
    lattice_coordinate(1)=(cartesian_coordinate|b2).simplify();
    lattice_coordinate(2)=(cartesian_coordinate|b3).simplify();
    
		if(LDEBUG){
      cerr << function_name << " cartesian coord=" << cartesian_coordinate << " --> " << lattice_coordinate << endl;
    }

    return lattice_coordinate;
  }
}

// *************************************************************************** 
// anrl::getXYZ2LatticeTransformation() 
// *************************************************************************** 
namespace anrl{
  Symbolic getXYZ2LatticeTransformation(const string& lattice_and_centering){
    
    // gets transformation matrix from XYZ coordinates (ITC equations) to 
    // the Cartesian coordinates with respect to the lattice vectors.
   
    Symbolic xyz2lattice("xyz2lattice",3,3);
    
    // ---------------------------------------------------------------------------
    // create symbolic list of characters
    Symbolic a("a");
    Symbolic b("b");
    Symbolic c("c");
    Symbolic cx("cx");
    Symbolic cy("cy");
    Symbolic cz("cz");
    Symbolic beta("beta");
    Symbolic gamma("gamma");

    // ---------------------------------------------------------------------------
    if(lattice_and_centering == "aP"){
	    xyz2lattice = ((a, b*cos(gamma), cx),
                     (_SYMBOLIC_ZERO_, b*sin(gamma), cy),
                     (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, cz));
	  } 
	                   
    else if(lattice_and_centering == "mP" || lattice_and_centering == "mC"){
	    xyz2lattice = ((a, _SYMBOLIC_ZERO_, c*cos(beta)),
                     (_SYMBOLIC_ZERO_, b, _SYMBOLIC_ZERO_),
                     (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c*sin(beta)));
    }
    else if(lattice_and_centering == "oP" || lattice_and_centering == "oC" || 
        lattice_and_centering == "oI" || lattice_and_centering == "oF"){
	    xyz2lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
                     (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
                     (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }
    else if(lattice_and_centering == "tP" || lattice_and_centering == "tI"){
	    xyz2lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
                     (_SYMBOLIC_ZERO_, a, _SYMBOLIC_ZERO_),
                     (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }
    else if(lattice_and_centering == "hP" || lattice_and_centering == "hR"){
	    xyz2lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
                     (_SYMBOLIC_ZERO_, a, _SYMBOLIC_ZERO_),
                     (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, c));
    }
    else if(lattice_and_centering == "cP" || lattice_and_centering == "cF" || lattice_and_centering == "cI"){
	    xyz2lattice = ((a, _SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_),
                     (_SYMBOLIC_ZERO_, a, _SYMBOLIC_ZERO_),
                     (_SYMBOLIC_ZERO_, _SYMBOLIC_ZERO_, a));
    }

    return xyz2lattice;
  }
}

// *************************************************************************** 
// anrl::getEquationsForCenteredLattices() 
// *************************************************************************** 
namespace anrl {
  vector<Symbolic> getEquationsForCenteredLattices(const string& lattice_and_centering, 
      const Symbolic& lattice, 
      const vector<Symbolic>& conventional_equations){ 
    
    // convert equations to lattice equations for centered lattice (C, I, F).

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl:getEquationsForCenteredLattices():";
    stringstream message;

    vector<Symbolic> lattice_equations;
    
		// ---------------------------------------------------------------------------
    // get symbolic transformation matrix for a particular lattice 
    Symbolic xyz2lattice = getXYZ2LatticeTransformation(lattice_and_centering);
    if(LDEBUG){ cerr << function_name << " symbolic transformation from xyz to lattice: " << xyz2lattice << endl; }

    // ---------------------------------------------------------------------------
    // transform symbolic coordinates 
    for(uint i=0;i<conventional_equations.size();i++){
      Symbolic vec = xyz2lattice*conventional_equations[i];
      lattice_equations.push_back(cartesian2lattice(lattice, vec).simplify());
    }

    // ---------------------------------------------------------------------------
    // [LDEBUG] print conventional vs centered lattice equations
		if(LDEBUG){
    	for(uint i=0;i<conventional_equations.size();i++){
        cerr << "conventional: " << conventional_equations[i] << " to centered lattice: " << lattice_equations[i] << endl;
      }
    }
   
    // ---------------------------------------------------------------------------
    // bring lattice equations in cell after transformation 
    for(uint i=0;i<lattice_equations.size();i++){
      lattice_equations[i] = symbolic::BringInCell(lattice_equations[i]);
    }
    
		// ---------------------------------------------------------------------------
    // [LDEBUG] print conventional vs centered lattice equations AFTER bring in cell
		if(LDEBUG){
    	for(uint i=0;i<conventional_equations.size();i++){
        cerr << "conventional: " << conventional_equations[i] << " to centered lattice: " << lattice_equations[i] << endl;
      }
    }

		// ---------------------------------------------------------------------------
    // remove duplicates after bring in cell
    vector<Symbolic> primitive_lattice_equations;
    for(uint i=0;i<lattice_equations.size();i++){
      bool is_unique_equation = true;
      for(uint j=0;j<primitive_lattice_equations.size();j++){
        if(symbolic::isEqualVector(primitive_lattice_equations[j],lattice_equations[i])){
          is_unique_equation = false;
          break;
        }
      }
      if(is_unique_equation){ primitive_lattice_equations.push_back(lattice_equations[i]); }
    }
		
		if(LDEBUG){
      cerr << "# equations before: " << lattice_equations.size() << " vs # equations after: " << primitive_lattice_equations.size() << endl;
    }

    return primitive_lattice_equations;
  }
}

// *************************************************************************** 
// anrl::getEquationsForCenteredLattices() 
// *************************************************************************** 
namespace anrl {
  vector<Symbolic> convertEquations2FractionalEquations(const string& lattice_and_centering, 
      const Symbolic& lattice, 
      const vector<Symbolic> conventional_equations){

    // convert equations to lattice equations (fractional).

    vector<Symbolic> transformed_equations;

    // ---------------------------------------------------------------------------
    // if C-, I-, F-centered, then the equations from ITC are xyz coordinates
    if(lattice_and_centering == "mC" || lattice_and_centering == "oC" || lattice_and_centering == "oI" ||
       lattice_and_centering == "oF" || lattice_and_centering == "tI" || lattice_and_centering == "cI" ||
       lattice_and_centering == "cF"){
      transformed_equations = getEquationsForCenteredLattices(lattice_and_centering, lattice, conventional_equations);  
    }

    // ---------------------------------------------------------------------------
    // if hexagonal/trigonal/rhombohedral, then 
    else if(lattice_and_centering == "hP" || lattice_and_centering == "hR"){
      transformed_equations = conventional_equations;
    }
    // ---------------------------------------------------------------------------
    // else 
    else{
      transformed_equations = conventional_equations;
    }

    return transformed_equations;
  }
}

// *************************************************************************** 
// symbolic::addSymbolicEquation2Atoms() 
// *************************************************************************** 
namespace anrl {
  void addSymbolicEquation2Atoms(const vector<Symbolic>& equations, deque<_atom>& atoms, bool isfpos){
   
    // add symbolic equations to atom.fpos_equation or atom.cpos_equation
    // converts from variables from Symbolic to string
    // DEFAULT: update atom.fpos_equation

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "anrl::addSymbolicEquation2Atoms():";
    stringstream message;

    // ---------------------------------------------------------------------------
    // ensure sizes of atoms and symbolic equations match 
		if(equations.size() != atoms.size()){
			message << "The number of equations and atoms do not match. Check tolerances. #equations=" << equations.size() << ", #atoms=" << atoms.size();
			throw aurostd::xerror(_AFLOW_FILE_NAME_,function_name,message,_RUNTIME_ERROR_);
		}

    // ---------------------------------------------------------------------------
    // ensure sizes of atoms and symbolic equations match 
		for(uint i=0;i<atoms.size();i++){
      vector<string> coordinate;
			for(uint j=0;j<3;j++){
				stringstream ss_pos; ss_pos << equations[i].row(j);
				coordinate.push_back(ss_pos.str());
      }
      if(isfpos){ atoms[i].fpos_equation = coordinate; }
      else{ atoms[i].cpos_equation = coordinate; }
		}

    // ---------------------------------------------------------------------------
    // [LDEBUG] print string (fpos/cpos) version of symbolic equation
		if(LDEBUG){
			if(isfpos){
			  for(uint i=0;i<atoms.size();i++){
          cerr << function_name << " fpos_equation=" << aurostd::joinWDelimiter(atoms[i].fpos_equation,",") << endl;
        } 
      }
			else{
			  for(uint i=0;i<atoms.size();i++){
          cerr << function_name << " cpos_equation=" << aurostd::joinWDelimiter(atoms[i].cpos_equation,",") << endl;
        } 
      }
    }
	}
}
#endif // COMPILE_SYMBOLIC

#endif // _AFLOW_ANRL_CPP

// Written by David Hicks (DX) - 2020
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow DAVID HICKS - Duke University 2014-2020                 *
// *                                                                         *
// ***************************************************************************
