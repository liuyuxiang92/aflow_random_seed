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

   class xprototype {
         // label/params info
         string catalog;                                    // prototype catalog 'anrl' or 'htqc'
         string label;                                      // label (e.g., 201 or AB_cF8_225_a_b)
         vector<string> parameter_list;                     // list of degrees of freedom (a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,...)
         vector<double> parameter_values;                   // values for degrees of freedom
         string parameter_set_id;                           // parameter set enumertion (e.g., 001, 002, 003, etc.)
         string weblink;                                    // link to the corresponding CrystalDatabase web page

         // symmetry
         vector<uint> stoichiometry;                        // reduced stoichiometry for prototype (e.g., equicompositional ternary=1:1:1)
         string Pearson_symbol;                             // Pearson symbol
         uint space_group_number;                           // space group number
         string space_group_symbol_H_M;                     // space group symbol Hermann-Mauguin (optional or use AFLOW lookup table)
         string space_group_symbol_Hall;                    // space group symbol Hall (optional or use AFLOW lookup table)
         string space_group_symbol_Schoenflies;             // space group symbol Schoenflies (optional or use AFLOW lookup table)
         vector<vector<string> > Wyckoff_letters;           // list of Wyckoff letters grouped by species ([[a,b],[c,d,e],[f,g,h,i],...])
         vector<vector<string> > Wyckoff_site_symmetries;   // list of Wyckoff site symmetries grouped by species ([mmm],[2mm,m2m],[mm2],...]) (optional, I can grab from look-up table)
         vector<vector<uint> > Wyckoff_multiplicities;      // list of Wyckoff multiplicities grouped by species ([48],[24,24],[12,12,12][4,4,4,4],...]) (optional, I can grab from look-up table)

         // designations
         string prototype_material;                         // common prototype material, e.g., NaCl
         string common_name;                                // common prototype name, e.g., half-Heusler
         string mineral_name;                               // mineral name, e.g., corundum
         string phase;                                      // compound phase designation (alpha, beta, gamma, delta, etc.) (if applicable)
         string strukturbericht;                            // Strukturbericht designation (if applicable)
         vector<string> similar_materials;                  // list of similar compounds (if in same order as stoichiometry we can easily decorate prototypes)
         vector<string> comments;                           // noteworthy comments (included in ANRL document and webpage)
         string title;                                      // title (for ANRL document/webpage)

   };

// stannite example
_AFLOW_PROTOTYPE_ENCYCLOPEDIA_="http://aflow.org/CrystalDatabase/"

catalog="anrl";
label="A2BC4D_tI16_121_d_a_i_b";
parameter_list.push_back("a");parameter_list.push_back("c/a");parameter_list.push_back("x4");parameter_list.push_back("z4");
parameter_values.push_back(5.46);parameter_values.push_back(1.96428571429);parameter_values.push_back(0.245);parameter_values.push_back(0.132);
parameter_set_id="001";
weblink=_AFLOW_PROTOTYPE_ENCYCLOPEDIA_+"A2BC4D_tI16_121_d_a_i_b.html";

stoichiometry.push_back(2);stoichiometry.push_back(1);stoichiometry.push_back(4);stoichiometry.push_back(1);
Pearson_symbol="tI16";
space_group_number=121;
space_group_symbol_H_M="I-42m"; // or lookup table
space_group_symbol_Hall="I -4 2"; // or lookup table
space_group_symbol_Schoenflies="D_{2d}^{11}"; // or lookup table

vector<string> Wyckoff_tmp; 
// letters 
Wyckoff_tmp.push_back("d");Wyckoff_letters.push_back(Wyckoff_tmp);Wyckoff_tmp.clear();
Wyckoff_tmp.push_back("a");Wyckoff_letters.push_back(Wyckoff_tmp);Wyckoff_tmp.clear();
Wyckoff_tmp.push_back("i);Wyckoff_letters.push_back(Wyckoff_tmp);Wyckoff_tmp.clear();
Wyckoff_tmp.push_back("b");Wyckoff_letters.push_back(Wyckoff_tmp);Wyckoff_tmp.clear();
// site symmetries
Wyckoff_tmp.push_back("-4..");Wyckoff_site_symmetries.push_back(Wyckoff_tmp);Wyckoff_tmp.clear();
Wyckoff_tmp.push_back("-42m");Wyckoff_site_symmetries.push_back(Wyckoff_tmp);Wyckoff_tmp.clear();
Wyckoff_tmp.push_back("..m);Wyckoff_site_symmetries.push_back(Wyckoff_tmp);Wyckoff_tmp.clear();
Wyckoff_tmp.push_back("-42m");Wyckoff_site_symmetries.push_back(Wyckoff_tmp);Wyckoff_tmp.clear();
// site symmetries
vector<uint> Wyckoff_mult_tmp;
Wyckoff_mult_tmp.push_back(4);Wyckoff_multiplicities.push_back(Wyckoff_mult_tmp);Wyckoff_mult_tmp.clear();
Wyckoff_mult_tmp.push_back(2);Wyckoff_multiplicities.push_back(Wyckoff_mult_tmp);Wyckoff_mult_tmp.clear();
Wyckoff_mult_tmp.push_back(8);Wyckoff_multiplicities.push_back(Wyckoff_mult_tmp);Wyckoff_mult_tmp.clear();
Wyckoff_mult_tmp.push_back(2);Wyckoff_multiplicities.push_back(Wyckoff_mult_tmp);Wyckoff_mult_tmp.clear();

prototype_material="Cu2FeS4Sn";
common_name=NNN;
mineral_name="stannite";
phase=NNN;
strukturbericht="H2_{6}";
similar_materials.push_back("Cu2CdSe4Sn");similar_materials.push_back("CoCu2S4Sn");similar_materials.push_back("Cu2GeHgS4");similar_materials.push_back("Cu2HgS4Sn");similar_materials.push_back("Ag2FeS4Sn");
comments="If $c=2a$, $x=1/4$, and $z=3/8$, the atoms are on the sites of the diamond ($A4$) structure. If, in addition, the Cu, Fe, and Sn atoms are replaced by a single atom type, the crystal reduces to the zincblende ($B3$) structure."
title="Stannite (Cu2FeS4Sn, H2_{6}) Structure";

 
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
 /* 
  // ---------------------------------------------------------------------------
  // get existing prototype information
  uint PrototypeANRL_LoadList(vector<string>& vproto,
      vector<string>& vproto_label,
      vector<uint>& vproto_nspecies,
      vector<uint>& vproto_natoms,
      vector<uint>& vproto_spacegroup,
      vector<uint>& vproto_nunderscores,
      vector<uint>& vproto_nparameters,
      vector<string>& vproto_Pearson_symbol,
      vector<string>& vproto_params,
      vector<string>& vproto_Strukturbericht,
      vector<string>& vproto_prototype,
      vector<string>& vproto_dialect);
  vector<string> getANRLParameters(string anrl_label,
      string library="",
      int choice=-1,
      bool keep_original_lattice_parameter=false); //DX20181009 //DX20190227 - added keep_original_lattice_parameter
  bool vproto2tokens(string proto,
      string& label,
      uint& nspecies,
      uint& natoms,
      uint& spacegroup,
      uint& nunderscores,
      uint& nparameters,
      string& Pearson_symbol,
      string& params,
      string& Strukturbericht,
      string& prototype,
      string& dialect);
  
  // ---------------------------------------------------------------------------
  // map structure to label and internal degrees of freedom 
  string structure2anrl(istream& input, aurostd::xoption& vpflow);           // xoption
  string structure2anrl(xstructure& xstr, bool recalculate_symmetry=true);   // use default options //DX20191031 - added recalculate_symmetry
  string structure2anrl(xstructure& xstr, double tolerance);                 // specify symmetry tolerance //CO20190520 - removed pointers for bools and doubles, added const where possible
  string structure2anrl(xstructure& xstr, uint setting);                     // specify setting
  string structure2anrl(xstructure& xstr, double tolerance, uint setting, bool recalculate_symmetry=true);  // main function //CO20190520 - removed pointers for bools and doubles, added const where possible //DX20190829 - added recalculate_symmetry //DX20191031 - removed reference
  xstructure rhl2hex(xstructure& str, double& a, double& c); 
  
  // ---------------------------------------------------------------------------
  // helper functions to determine label and internal degrees of freedom 
  uint getANRLSettingChoice(int spacegroup); //DX20191031 - removed reference
  string groupedWyckoffPosition2ANRLString(const vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize);
  vector<string> getANRLLatticeParameterString(char& lattice_type);
  vector<double> getANRLLatticeParameterValuesFromWyccar(const vector<string>& wyccar_ITC, char lattice_type, char lattice_centering, uint setting); //DX20191031
  vector<double> getANRLLatticeParameterValuesFromABCAngles(const xstructure& xstr, char lattice_type, char lattice_centering, uint setting); //DX20191031
  vector<double> getANRLLatticeParameterValues(const vector<double>& all_lattice_parameters, char lattice_type, char lattice_centering, uint setting); //DX20191031
*/  
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
  //string extractANRLPrototypeParameterValues(const string& label_anrl, const string& number_id, const string& variables, bool& keep_anrl_lattice_parameter); 
  deque<_atom> getAtomsFromWyckoff(const vector<wyckoffsite_ITC>& Wyckoff_positions, const xmatrix<double>& lattice_conventional);
  vector<string> determineWyckoffVariables(vector<wyckoffsite_ITC>& Wyckoff_positions);
  void applyWyckoffValues(const vector<double>& Wyckoff_parameter_values,vector<wyckoffsite_ITC>& Wyckoff_positions); //perhaps add this as a method to wyckoffsite_ITC?
  bool containsDuplicateWyckoffCoordinate(const vector<wyckoffsite_ITC>& wyckoff_sites_ITC, bool already_ordered=false);
  vector<wyckoffsite_ITC> getWyckoffSitesFromANRL(const vector<string>& Wyckoff_tokens, const vector<string>& species, uint space_group_number, int setting=SG_SETTING_ANRL); 
  
  // ---------------------------------------------------------------------------
  // checking functions
  //vector<uint> extractStoichiometry(string& anrl_label);
  //bool PrototypeANRL_Consistency(uint vparameters_size,uint proto_nparameters,string proto_prototype,
  //    string proto_label,string proto_Strukturbericht,string proto_Pearson_symbol,
  //    uint proto_spacegroup, string proto_params, uint print_mode); //DX20180710 - added print_mode //DX20200207 - oss no longer needed
  bool structureAndLabelConsistent(const xstructure& _xstr, const string& label_input, string& label_and_params_calculated);
  /*
  // ---------------------------------------------------------------------------
  // generic prototype generator (main function)
  xstructure PrototypeANRL_Generator(string& label,
      string& parameters,
      deque<string> &vatomX,
      deque<double> &vvolumeX,
      ostream& logstream=cout,
      bool silence_logger=true); //DX20200528 - command line = no logger
  xstructure PrototypeANRL_Generator(string& label,
      string& parameters,
      deque<string> &vatomX,
      deque<double> &vvolumeX,
      ofstream& FileMESSAGE,
      ostream& logstream=cout,
      bool silence_logger=false); //DX20200528 - internal = logger
  
  // ---------------------------------------------------------------------------
  // [OLD] hard-coded generator (requires ANRL/ subdirectory)
  xstructure PrototypeANRL(ostream &oss,
      string label,
      string parameters,
      deque<string> &vatomX,
      deque<double> &vvolumeX,
      double volume_in,
      int mode,
      bool flip_option);*/
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
