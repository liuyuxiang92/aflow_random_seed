// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow DAVID HICKS - Duke University 2014-2020                 *
// *                                                                         *
// ***************************************************************************
// Written by David Hicks (DX) - 2020

#ifndef _AFLOW_AFLOW_PROTOTYPE_GENERATOR_H_
#define _AFLOW_AFLOW_PROTOTYPE_GENERATOR_H_
#include "aflow.h"
#include "aflow_symmetry_spacegroup.h"
#include "aflow_compare_structure.h"
#include "SYMBOLICC++/symbolicc++.h"

// printing modes
#define _PROTO_GENERATOR_GEOMETRY_FILE_ 0
#define _PROTO_GENERATOR_EQUATIONS_ONLY_ 1
#define _PROTO_GENERATOR_GEOMETRY_FILE_AND_EQUATIONS_ 2

// toggle symbolic math
#define COMPILE_SYMBOLIC

// Symbolic math variables
#ifdef COMPILE_SYMBOLIC
#define _SYMBOLIC_ZERO_ Symbolic(0)
#define _SYMBOLIC_TOL_ 1e-3 
#define _ANRL_LATTICE_VARIABLES_  "a,b,c,alpha,beta,gamma,cx,cy,cz"
#define _ANRL_TRIG_VARIABLES_  "sin,cos,tan,sec,csc,cot"
#define _ANRL_WYCKOFF_VARIABLES_  "x,y,z"
#endif

/*
namespace prototype {
  class Prototype{
    public:
      Prototype()                                                           //constructor operator
      ~Prototype()                                                          //destructor operator
      friend ostream& operator<<(ostream& oss, const Prototype& Prototype); // stringstream operator (printing)
      const Prototype& operator=(const Prototype& b);                       // assignment operator
      Prototype(const Prototype& b);                                        // copy constructor
    
      // label/params info
      string catalog;                                                       // prototype catalog 'anrl' or 'htqc'
      string label;                                                         // label (e.g., 201 or AB_cF8_225_a_b) 
      vector<string> parameter_list;                                        // list of degrees of freedom (a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,...)
      vector<double> parameter_values;                                      // values for degrees of freedom
      string parameter_set_id;                                              // parameter set enumertion (e.g., 001, 002, 003, etc.)
      string weblink;                                                       // link to the corresponding CrystalDatabase web page 
    
      // symmetry
      vector<uint> stoichiometry;                                           // reduced stoichiometry for prototype (e.g., equicompositional ternary=1:1:1)
      string Pearson_symbol;                                                // Pearson symbol
      uint space_group_number;                                              // space group number
      string space_group_symbol_H_M;                                        // space group symbol Hermann-Mauguin (optional or use AFLOW lookup table)
      string space_group_symbol_Hall;                                       // space group symbol Hall (optional or use AFLOW lookup table)
      string space_group_symbol_Schoenflies;                                // space group symbol Schoenflies (optional or use AFLOW lookup table)
      vector<vector<string> > Wyckoff_letters;                              // list of Wyckoff letters grouped by species ([[a,b],[c,d,e],[f,g,h,i],...])      
      vector<vector<string> > Wyckoff_site_symmetries;                      // list of Wyckoff site symmetries grouped by species ([mmm],[2mm,m2m],[mm2],...]) (optional, I can grab from look-up table)
      vector<vector<uint> > Wyckoff_multiplicities;                         // list of Wyckoff multiplicities grouped by species ([48],[24,24],[12,12,12][4,4,4,4],...]) (optional, I can grab from look-up table)

      // designations
      string prototype_material;                                            // common prototype material, e.g., NaCl
      string common_name;                                                   // common prototype name, e.g., half-Heusler
      string mineral_name;                                                  // mineral name, e.g., corundum
      string phase;                                                         // compound phase designation (alpha, beta, gamma, delta, etc.) (if applicable)
      string strukturbericht;                                               // Strukturbericht designation (if applicable)
      vector<string> similar_materials;                                     // list of similar compounds (if in same order as stoichometry we can easily decorate prototypes)
      vector<string> comments;                                              // noteworthy comments (included in ANRL document and webpage)
      string title;                                                         // title (for ANRL document/webpage)

      // citation info
      vector<string> references                                             // references (do we want DOIs, formatted refs?)
      vector<string> found_in_references                                    // references where the structure was first identified (e.g., Villars, Strukturbericht series)

      xstructure structure;
      vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;                               // Wyckoff positions grouped by site type
      
      vector<string> lattice_vectors
      vector<string> atom_positions;

      string latex_string;
      string html_string;
      
      //symbolic math

      //methods
      void generateXstructure();
      void generateSymbolicMath();
      void generateLatexPage();
      void generateHTMLPage();

      bool isNewPrototype();
      string findExistingPrototypeLabel();
      void loadExistingPrototypeData(const string& label);


  }
}
*/

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
  
  // ---------------------------------------------------------------------------
  // functions to determine atomic positions from Wyckoff and parameters
  string extractANRLPrototypeParameterValues(const string& label_anrl, const string& number_id, const string& variables, bool& keep_anrl_lattice_parameter); 
  deque<_atom> getAtomsFromWyckoff(const vector<wyckoffsite_ITC>& Wyckoff_positions, const xmatrix<double>& lattice_conventional);
  vector<string> determineWyckoffVariables(vector<wyckoffsite_ITC>& Wyckoff_positions);
  void applyWyckoffValues(const vector<double>& Wyckoff_parameter_values,vector<wyckoffsite_ITC>& Wyckoff_positions); //perhaps add this as a method to wyckoffsite_ITC?
  bool containsDuplicateWyckoffCoordinate(const vector<wyckoffsite_ITC>& wyckoff_sites_ITC, bool already_ordered=false);
  vector<wyckoffsite_ITC> getWyckoffSitesFromANRL(const vector<string>& Wyckoff_tokens, const vector<string>& species, uint space_group_number, int setting=SG_SETTING_ANRL); 
  
  // ---------------------------------------------------------------------------
  // checking functions
  bool structureAndLabelConsistent(const xstructure& _xstr, const string& label_input, string& label_and_params_calculated);
}
//SYMBOLIC

#ifdef COMPILE_SYMBOLIC
//SYMBOLIC

// ---------------------------------------------------------------------------
// mirrors wyckoffsite_ITC, requires SymbolicC++ 
// instead of adding Symbolic equations to wyckoffsite_ITC (don't want to put 
// large dependency in aflow.h)
struct SymbolicWyckoffSite {
  xvector<double> coord;
  uint index;
  string type;
  string wyckoffSymbol;
  string letter;
  string site_symmetry;
  uint multiplicity;
  double site_occupation;
  vector<Symbolic> equations;
  uint parameter_index;
};

SymbolicWyckoffSite initializeSymbolicWyckoffSite(const wyckoffsite_ITC& Wyckoff);
void substituteVariableWithParameterDesignation(vector<SymbolicWyckoffSite>& Wyckoff_sites_symbolic);
void substituteVariableWithParameterDesignation(SymbolicWyckoffSite& Wyckoff_symbolic);

// ---------------------------------------------------------------------------
// "pure" symbolic functions 
namespace symbolic {
  Symbolic string2symbolic(const string& str);
  bool isEqual(const Symbolic& a, const Symbolic& b);
  bool isEqualVector(const Symbolic& a_vec, const Symbolic& b_vec);
  vector<vector<string> > matrix2VectorVectorString(const Symbolic& lattice);
  Symbolic BringInCell(const Symbolic& vec_in, double tolerance=_ZERO_TOL_, double upper_bound=1.0, double lower_bound=0.0);
}

// ---------------------------------------------------------------------------
// ANRL symbolic functions 
namespace anrl {
  Symbolic SymbolicANRLPrimitiveLattices(const string& lattice_and_centering, const char& space_group_letter);
  vector<Symbolic> equations2SymbolicEquations(const vector<vector<string> >& equations);
  Symbolic cartesian2lattice(const Symbolic& lattice, const Symbolic& cartesian_coordinate); 
  Symbolic getXYZ2LatticeTransformation(const string& lattice_and_centering);
  vector<Symbolic> getEquationsForCenteredLattices(const string& lattice_and_centering, 
      const Symbolic& lattice, 
      const vector<Symbolic>& conventional_equations);  
  vector<Symbolic> convertEquations2FractionalEquations(const string& lattice_and_centering, 
      const Symbolic& lattice, 
      const vector<Symbolic> conventional_equations);
  void addSymbolicEquation2Atoms(const vector<Symbolic>& equations, deque<_atom>& atoms, bool isfpos=true);
}
#endif 

#endif 
