// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow DAVID HICKS - Duke University 2014-2020                 *
// *                                                                         *
// ***************************************************************************
// Written by David Hicks (DX) - 2020

#ifndef _AFLOW_ANRL_H_
#define _AFLOW_ANRL_H_
#include "aflow.h"
#include "aflow_symmetry_spacegroup.h"
#include "aflow_compare_structure.h"
#include "aflow_symbolic.h"

// printing modes
#define _PROTO_GENERATOR_GEOMETRY_FILE_ 0
#define _PROTO_GENERATOR_EQUATIONS_ONLY_ 1
#define _PROTO_GENERATOR_GEOMETRY_FILE_AND_EQUATIONS_ 2

// toggle symbolic math
#define COMPILE_SYMBOLIC

// Symbolic math variables
#ifdef COMPILE_SYMBOLIC
#define _SYMBOLIC_ZERO_ symbolic::Symbolic(0)
#define _SYMBOLIC_TOL_ 1e-3 
#define _ANRL_LATTICE_VARIABLES_  "a,b,c,alpha,beta,gamma,cx,cy,cz"
#define _ANRL_TRIG_VARIABLES_  "sin,cos,tan,sec,csc,cot"
#define _ANRL_WYCKOFF_VARIABLES_  "x,y,z"
#endif
  
// ---------------------------------------------------------------------------
// below are functions limited to the aflow_anrl.cpp file, if you want to
// include functions into other parts of aflow, put the functions in aflow.h
// (since the SYMBOLICC++ is all header files, i.e., inline definitions, we
// cannot import the aflow_anrl.h file into other files)
namespace anrl{
  // ---------------------------------------------------------------------------
  // get lattice functions (primitive/conventional) - ITC/ANRL standard
  xmatrix<double> getLattice(const string& lattice_and_centering, const char& space_group_letter, const vector<double>& lattice_parameter_values, uint mode=0);
  xmatrix<double> getTriclinicLattice(const vector<double>& lattice_parameter_values, uint mode=0);
  xmatrix<double> getMonoclinicLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode=0);
  xmatrix<double> getOrthorhombicLattice(const string& lattice_and_centering, const char& space_group_letter, const vector<double>& lattice_parameter_values, uint mode=0);
  xmatrix<double> getTetragonalLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode=0);
  xmatrix<double> getHexagonalLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode=0);
  xmatrix<double> getCubicLattice(const string& lattice_and_centering, const vector<double>& lattice_parameter_values, uint mode=0);
  xstructure rhl2hex(const xstructure& str, double& a, double& c); 
  
  // ---------------------------------------------------------------------------
  // functions to determine atomic positions from Wyckoff and parameters
  deque<_atom> getAtomsFromWyckoff(const vector<wyckoffsite_ITC>& Wyckoff_positions, const xmatrix<double>& lattice_conventional);
  vector<string> determineWyckoffVariables(vector<wyckoffsite_ITC>& Wyckoff_positions);
  void applyWyckoffValues(const vector<double>& Wyckoff_parameter_values,vector<wyckoffsite_ITC>& Wyckoff_positions); //perhaps add this as a method to wyckoffsite_ITC?
  bool containsDuplicateWyckoffCoordinate(const vector<wyckoffsite_ITC>& wyckoff_sites_ITC, bool already_ordered=false);
  vector<wyckoffsite_ITC> getWyckoffSitesFromANRL(const vector<string>& Wyckoff_tokens, const vector<string>& species, uint space_group_number, int setting=SG_SETTING_ANRL); 
  
  // ---------------------------------------------------------------------------
  // checking functions
  bool structureAndLabelConsistent(const xstructure& _xstr, const string& label_input, string& label_and_params_calculated);
}

// Symbolic functions are defined in aflow_symbolic.cpp

#endif 
