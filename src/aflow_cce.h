// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2020              *
// *                                                                         *
// ***************************************************************************
// Written by Rico Friedrich, Corey Oses, and Marco Esters
// rico.friedrich@duke.edu

#ifndef _AFLOW_CCE_H_
#define _AFLOW_CCE_H_

static const string CCE_allowed_functionals = "PBE,LDA,SCAN,PBE+U_ICSD,exp"; // when adding a new functional also introduce new 'offset' in get_offset function needed for reading corrections from lookup table
static const string CCE_default_output_functionals= "PBE,LDA,SCAN,exp"; // corrections are given for these functionals if only a structure is given as input for the command line and web tools (i.e. --functionals= is not set)
static const string CCE_temperatures = "298.15,0"; // needs to be extended when adding new corrections for other temperatures
static const double _CCE_NN_DIST_TOL_ = 0.5; // 0.5 Ang tolerance between shortest and longest bonds for each cation-anion pair; works best up to now; in future maybe bonding could be explicitly determined via Bader analysis
static const double _CCE_NN_DIST_TOL_MULTI_ANION_ = 0.4; // 0.4 Ang tolerance between shortest and longest bonds for each bond when testing for multi-anion compound; it was found that the standard 0.5 Ang tol. is too large such that different anions appear to be bonded, which would prevent anions to be detected as such
static const double _CCE_OX_TOL_ = 0.001; // choose small finite value since sum of oxidation states might not be exactly zero due to numerics
static const double _CCE_SELF_DIST_TOL_ = 0.001; // distance tolerance in Ang for neighbor screening to savely exclude the cation itself having distance zero to itself
static const double _CCE_perox_cutoff_=1.6; // O-O bonds in peroxides for the studied examples are all shorter than 1.6 Ang
static const double _CCE_superox_cutoff_=1.4; // O-O bonds in superoxides for the studied examples are all shorter than 1.4 Ang
static const double _CCE_O2_molecule_cutoff_=1.2; // O-O bonds in the O2 molecule is about 1.21 Ang.
static const uint CCE_num_functionals_Bader=4; // Currently, Bader charges used to determine oxidation states are available for 4 functionals: PBE, LDA, SCAN, and PBE+U_ICSD and ONLY for oxides, see get_Bader_templates function

namespace cce {
  struct CCE_Variables {
    vector<double> dft_energies;
    vector<string> vfunctionals; // should be needed as long as output for corrected dft formation energies is based on vfunctionals
    vector<int> offset; // needed for reading corrections from lookup table for different functionals
    vector<string> vtemperatures;
    double standard_anion_charge;
    vector<double> electronegativities;
    vector<uint> multi_anion_atoms; // vector in which elements will be 1 for multi_anion atoms and 0 otherwise
    vector<double> oxidation_states;
    string anion_species;
    vector<double> cutoffs;
    vector<string> multi_anion_species; // vector storing all the multi anion species
    uint num_perox_bonds;
    uint num_superox_bonds;
    vector<uint> perox_indices; // vector in which elements will be 1 for peroxide O atoms and 0 otherwise; needed for correct setting of oxidation numbers below
    vector<uint> superox_indices; // vector in which elements will be 1 for superoxide O atoms and 0 otherwise; needed for correct setting of oxidation numbers below
    vector<uint> num_neighbors;
    vector<string> species_electronegativity_sorted;
    vector<int> num_pref_ox_states_electronegativity_sorted;
    vector< vector<double> > pref_ox_states_electronegativity_sorted;
    vector<int> num_all_ox_states_electronegativity_sorted;
    vector< vector<double> > all_ox_states_electronegativity_sorted;
    vector<vector<int> > cations_map;
    double oxidation_sum; // double because for superoxides O ox. number is -0.5
    vector<double> Bader_charges;
    vector<vector<double> > corrections_atom; // 1st dim. is number of functionals*2 (vfunctionals.size()*2) i.e. 298.15 & 0K corrections for each functional; 2nd dimension for corrections for each atom
    vector<vector<vector<double> > > multi_anion_corrections_atom; // 1st dim. for multi_anion_species, 2nd dim. for functionals and temperatures as above, 3rd dim. for corrections for each atom
    vector<double> perox_correction; // peroxide correction per cell for functionals and temperatures as above
    vector<double> superox_correction; // superoxide correction per cell for functionals and temperatures as above
    vector<double> cce_correction; // total correction per cell for functionals and temperatures as above
    vector<double> cce_form_energy_cell; // CCE formation enthalpy per cell for functionals and temperatures as above
  };

  // main CCE functions
  // for command line use, 
  // use inside AFLOW providing directory path or xstructure & functional string or flags and istream for web tool, 
  // and CCE core function called by all other main CCE functions
  void write_corrections(aurostd::xoption& flags, ostream& oss=std::cout);
  void write_corrections(xstructure& structure, aurostd::xoption& flags);
  void write_corrections(xstructure& structure, aurostd::xoption& flags, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void write_corrections(aurostd::xoption& flags, std::istream& ist); // ME20200213
  vector<double> calculate_corrections(const string& directory_path);
  vector<double> calculate_corrections(xstructure& structure, string& functional);
  void CCE_core(xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars);
  // read user input (from command line or directory path)
  xstructure read_structure(const string& structure_file, int=IOAFLOW_AUTO); // set xstructure mode argument only here and it is automoatically recognized in the main CCE cpp file
  xstructure read_structure(std::istream& ist);
  xstructure read_structure(xstructure& structure);
  void get_dft_form_energies_functionals(const string& dft_energies_input_str, const string& functionals_input_str, CCE_Variables& cce_vars);
  int get_offset(const string& functional);
  vector<double> get_oxidation_states(const string& oxidation_numbers_input_str, const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  string get_functional_from_aflow_in(const xstructure& structure, string& aflowin_file, string& outcar_file);
  // initialise flags and variables
  aurostd::xoption init_flags(); // ME20200213
  CCE_Variables init_variables(const xstructure&); // ME20200213
  // structural analysis
  string determine_anion_species(const xstructure& structure, CCE_Variables& cce_vars);
  vector<uint> check_for_multi_anion_system(xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, double tolerance=_CCE_NN_DIST_TOL_MULTI_ANION_, ostream& oss=std::cout);
  vector<uint> get_num_neighbors(xstructure& structure, double tolerance=_CCE_NN_DIST_TOL_);
  vector<uint> get_num_neighbors(xstructure& structure, const string& anion_species, double tolerance=_CCE_NN_DIST_TOL_);
  vector<uint> get_num_neighbors(xstructure& structure, const string& anion_species, xoption& cce_flags, CCE_Variables& cce_vars, double tolerance=_CCE_NN_DIST_TOL_, ostream& oss=std::cout);
  vector<double> get_dist_cutoffs(const xstructure& structure, double tolerance=_CCE_NN_DIST_TOL_);
  void check_per_super_oxides(xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  // determine oxidation numbers from electronegativities
  vector<double> get_oxidation_states_from_electronegativities(xstructure& structure);
  vector<double> get_oxidation_states_from_electronegativities(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void set_anion_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars);
  void sort_species_by_electronegativity(const xstructure& structure, CCE_Variables& cce_vars);
  void load_ox_states_templates_each_species(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars);
  void try_preferred_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars);
  void treat_SbO2_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void treat_Pb3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void treat_Ti_O_Magneli_phase_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void treat_Fe3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void treat_Mn3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void treat_Co3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void treat_alkali_sesquioxide_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  // following special cases only needed when determining oxidation states from Bader charges
  void treat_MnMoO4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void treat_Ca2Fe2O5_CaFe2O4_LDA_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void treat_FeTiO3_LDA_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void check_ox_nums_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss);
  void try_all_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars);
  void determine_cation_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars, const vector<vector<double> >& possible_ox_states); // ME Nov. 2019
  double print_oxidation_states_and_sum(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  double print_oxidation_states_and_get_sum(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  double get_oxidation_states_sum(CCE_Variables& cce_vars);
  double get_oxidation_states_sum(xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  // determine oxidation numbers from Bader charges (outdated)
  vector<double> get_oxidation_states_from_Bader(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  void get_system_name_functional_from_aflow_in(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, string& system_name, string& functional, ostream& oss=std::cout);
  vector<double> get_Bader_charges_from_Bader_file(const xstructure& structure, CCE_Variables& cce_vars, const string& Bader_file);
  vector<double> Bader_charges_to_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, string& functional);
  void general_attempt_fixing_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss=std::cout);
  // assign corrections
  void get_corrections(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, const string& considered_anion_species, const vector<uint>& num_neighbors, vector<vector<double> >& corrections_atom);
  void load_cation_corrections(const xstructure& structure, CCE_Variables& cce_vars, const string& corrections_line, vector<vector<double> >& corrections_atom, uint i);
  void set_anion_corrections(const xstructure& structure, CCE_Variables& cce_vars, vector<vector<double> >& corrections_atom, uint i);
  void check_get_per_super_ox_corrections(CCE_Variables& cce_vars);
  // apply corrections and get corrected formation enthalpies
  void check_apply_per_super_ox_corrections(CCE_Variables& cce_vars);
  //vector<double> get_formation_enthalpies(const vector<double>& cce_correction, CCE_Variables& cce_vars); // ME20200213
  // write output and citation
  string get_JSON(const xstructure& structure, const CCE_Variables& cce_vars); // ME20200213
  string write_output(const xstructure& structure, CCE_Variables& cce_vars, const vector<double>& cce_form_energy_cell);
  string write_test_output(CCE_Variables& cce_vars, const vector<double>& cce_form_energy_cell);
  string write_citation();
  // write user instructions
  string print_usage();
  // aflow_cce_data.cpp load corrections and other data
  string get_corrections_line(const string& cor_identifier);
  string get_Bader_templates(const string& element);
}

#endif  // _AFLOW_CCE_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2020              *
// *                                                                         *
// ***************************************************************************
