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

namespace cce {
  struct CCE_Variables {
    vector<double> dft_energies;
    vector<string> vfunctionals; // should be needed as long as output for corrected dft formation energies is based on vfunctionals
    vector<int> offset; // needed for reading corrections from lookup table for different functionals
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
  void CCE(aurostd::xoption& flags);
  vector<double> CCE_correct(const string& directory_path);
  vector<double> CCE_correct(xstructure& structure, string& functional);
  void CCE(aurostd::xoption& flags, std::istream& ist); // ME 200213
  void CCE_core(xstructure& structure, CCE_Variables& cce_vars, xoption& cce_flags);
  // read user input (from command line or directory path)
  xstructure CCE_read_structure(const string& structure_file, int=IOAFLOW_AUTO); // set xstructure mode argument only here and it is automoatically recognized in the main CCE cpp file
  void CCE_get_dft_form_energies_functionals(const string& dft_energies_input_str, const string& functionals_input_str, CCE_Variables& cce_vars);
  int CCE_get_offset(const string& functional);
  vector<double> CCE_get_oxidation_states(const string& oxidation_numbers_input_str, const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars);
  string CCE_get_functional_from_aflow_in(const xstructure& structure, string& aflowin_file);
  // initialise flags and variables
  aurostd::xoption CCE_init_flags(); // ME 200213
  CCE_Variables CCE_init_variables(const xstructure&); // ME 200213
  // structural analysis
  string CCE_determine_anion_species(const xstructure& structure, CCE_Variables& cce_vars);
  vector<uint> CCE_check_for_multi_anion_system(CCE_Variables& cce_vars, double tolerance, xstructure& structure, xoption& cce_flags);
  vector<double> CCE_get_dist_cutoffs(double tolerance, const xstructure& structure);
  vector<uint> CCE_get_num_neighbors(double tolerance, xstructure& structure);
  vector<uint> CCE_get_num_neighbors(const string& anion_species, double tolerance, xstructure& structure);
  vector<uint> CCE_get_num_neighbors(const string& anion_species, double tolerance, xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars);
  void CCE_check_per_super_oxides(xstructure& structure, CCE_Variables& cce_vars, xoption& cce_flags);
  // determine oxidation numbers from electronegativities
  vector<double> CCE_get_oxidation_states_from_electronegativities(xstructure& structure);
  vector<double> CCE_get_oxidation_states_from_electronegativities(const xstructure& structure, CCE_Variables& cce_vars, xoption& cce_flags);
  void CCE_set_anion_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars);
  void CCE_sort_species_by_electronegativity(const xstructure& structure, CCE_Variables& cce_vars);
  void CCE_load_ox_states_templates_each_species(const xstructure& structure, CCE_Variables& cce_vars, xoption& cce_flags);
  void CCE_try_preferred_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars);
  void CCE_treat_SbO2_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  void CCE_treat_Pb3O4_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  void CCE_treat_Ti_O_Magneli_phase_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  void CCE_treat_Fe3O4_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  void CCE_treat_Mn3O4_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  void CCE_treat_Co3O4_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  void CCE_treat_alkali_sesquioxide_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  // following special cases only needed when determining oxidation states from Bader charges
  void CCE_treat_MnMoO4_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  void CCE_treat_Ca2Fe2O5_CaFe2O4_LDA_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  void CCE_treat_FeTiO3_LDA_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  void CCE_try_all_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars);
  void CCE_determine_cation_oxidation_states(const vector<vector<double> >& possible_ox_states, const xstructure& structure, CCE_Variables& cce_vars); // ME Nov. 2019
  double CCE_print_oxidation_states_and_sum(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  double CCE_print_oxidation_states_and_get_sum(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  double CCE_get_oxidation_states_sum(CCE_Variables& cce_vars);
  double CCE_get_oxidation_states_sum(CCE_Variables& cce_vars, xoption& cce_flags);
  // determine oxidation numbers from Bader charges (outdated)
  vector<double> CCE_get_oxidation_states_from_Bader(xoption& cce_flags, const xstructure& structure, CCE_Variables& cce_vars);
  void CCE_get_system_name_functional_from_aflow_in(const xstructure& structure, string& system_name, xoption& cce_flags, string& functional, CCE_Variables& cce_vars);
  vector<double> CCE_get_Bader_charges_from_Bader_file(const string& Bader_file, const xstructure& structure, CCE_Variables& cce_vars);
  vector<double> CCE_Bader_charges_to_oxidation_states(xoption& cce_flags, const xstructure& structure, CCE_Variables& cce_vars, string& functional);
  void CCE_general_attempt_fixing_oxidation_states(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags);
  // assign corrections
  void CCE_get_corrections(xoption& cce_flags, const xstructure& structure, vector<vector<double> >& corrections_atom, CCE_Variables& cce_vars, const vector<uint>& num_neighbors, const string& considered_anion_species);
  void CCE_load_cation_corrections(CCE_Variables& cce_vars, const xstructure& structure, vector<vector<double> >& corrections_atom, const string& corrections_line, uint i);
  void CCE_set_anion_corrections(CCE_Variables& cce_vars, const xstructure& structure, vector<vector<double> >& corrections_atom, uint i);
  void CCE_check_get_per_super_ox_corrections(CCE_Variables& cce_vars);
  // apply corrections and get corrected formation enthalpies
  void CCE_check_apply_per_super_ox_corrections(CCE_Variables& cce_vars);
  vector<double> CCE_get_formation_enthalpies(const vector<double>& cce_correction, CCE_Variables& cce_vars); // ME 200213
  // write output and citation
  string CCE_get_JSON(const xstructure& structure, const CCE_Variables& cce_vars); // ME 200213
  string CCE_write_output(const xstructure& structure, CCE_Variables& cce_vars, const vector<double>& cce_form_energy_cell);
  string CCE_write_test_output(CCE_Variables& cce_vars, const vector<double>& cce_form_energy_cell);
  string CCE_write_citation();
  // write user instructions
  string CCE_print_usage();
  // aflow_cce_data.cpp load corrections and other data
  string CCE_get_corrections_line(const string& cor_identifier);
  string CCE_get_Bader_templates(const string& element);
}

#endif  // _AFLOW_CCE_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2020              *
// *                                                                         *
// ***************************************************************************
