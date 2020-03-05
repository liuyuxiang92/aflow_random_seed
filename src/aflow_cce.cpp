// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2020              *
// *                                                                         *
// ***************************************************************************
// Written by Rico Friedrich, Corey Oses, and Marco Esters
// rico.friedrich@duke.edu

#ifndef _AFLOW_CCE_CPP_
#define _AFLOW_CCE_CPP_

#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_cce.h"

using std::cout;
using std::cerr;
using std::endl;

#define CCE_DEBUG false
static const vector<string> CCE_vallowed_functionals={"PBE","LDA","SCAN","PBE+U_ICSD","exp"}; // when adding a new functional also introduce new 'offset' in CCE_get_offset function needed for reading corrections from lookup table
static const vector<string> CCE_vdefault_output_functionals={"PBE","LDA","SCAN","exp"}; // corrections are given for these functionals if only a structure is given as input for the command line and web tools (i.e. --functionals= is not set)
static const double _CCE_NN_DIST_TOL_ = 0.5; // 0.5 Ang tolerance between shortest and longest bonds for each cation-anion pair; works best up to now; in future maybe bonding could be explicitly determined via Bader analysis
static const double _CCE_NN_DIST_TOL_MULTI_ANION_ = 0.4; // 0.4 Ang tolerance between shortest and longest bonds for each bond when testing for multi-anion compound; it was found that the standard 0.5 Ang tol. is too large such that different anions appear to be bonded, which would prevent anions to be detected as such
static const double _CCE_OX_TOL_ = 0.001; // choose small finite value since sum of oxidation states might not be exactly zero due to numerics
static const double _CCE_SELF_DIST_TOL_ = 0.1; // distance tolerance in Ang for neighbor screening to savely exclude the cation itself having distance zero to itself

namespace cce {

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                           MAIN CCE FUNCTIONS                            //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // for command line use, 
  // use inside AFLOW providing directory path or xstructure & functional string or flags and istream for web tool, 
  // and CCE core function called by all other main CCE functions

  //CCE////////////////////////////////////////////////////////
  // main CCE function for command line use 
  // for reading input, analyzing structure, determining oxidation numbers, assigning corrections, 
  // calculating total corrections and corrected formation enthalpies, and writing output
  void CCE(aurostd::xoption& flags) {
    //bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);

    /************************************/
    // Print user instructions
    /************************************/
    // option to print user instructions and exit upon completion
    if(flags.flag("CCE_CORRECTION::USAGE")){
      cout << CCE_print_usage();
      return;
    }

    /************************************/
    // Read structure
    /************************************/
    //read structural data from structure file provided on command line
    xstructure structure=CCE_read_structure(flags.getattachedscheme("CCE_CORRECTION::POSCAR_PATH"));
    aurostd::xoption cce_flags = CCE_init_flags();
    if(!flags.getattachedscheme("CCE_CORRECTION::POSCAR_PATH").empty()){ // when POSCAR_PATH is provided the command line must have called the function instead of another part of AFLOW
      cce_flags.flag("COMMAND_LINE",TRUE);
    }
    CCE_Variables cce_vars = CCE_init_variables(structure);

    /********************************************************/
    // Read DFT formation energies and functionals if provided
    /********************************************************/
    CCE_get_dft_form_energies_functionals(flags.getattachedscheme("CCE_CORRECTION::DFT_FORMATION_ENERGIES"), flags.getattachedscheme("CCE_CORRECTION::FUNCTIONALS"), cce_vars); //provide precalc. DFT formation energies & corresponding functionals

    /********************************************************/
    // Read oxidation numbers if provided
    /********************************************************/
    cce_flags.flag("OX_NUMS_PROVIDED",FALSE);
    if(flags.flag("CCE_CORRECTION::OXIDATION_NUMBERS")){
      cce_vars.oxidation_states = CCE_get_oxidation_states(flags.getattachedscheme("CCE_CORRECTION::OXIDATION_NUMBERS"), structure, cce_flags, cce_vars);
      cce_flags.flag("OX_NUMS_PROVIDED",TRUE);
    }

    /********************************************************/
    // Get CCE corrections
    /********************************************************/
    // obtain total CCE corrections per cell from other CCE function only based on 
    // structure, oxidation numbers and functional information
    cce_flags.flag("CORRECTABLE",TRUE); // first assuming that formation energy of system IS correctable; will be set to not correctable if, for any atom, no correction can be identified
    CCE_core(structure, cce_vars, cce_flags);

    // CALCULATE CORRECTED FORMATION ENTHALPIES AT 298.15 AND 0K ###################################
    if (cce_flags.flag("CORRECTABLE")){
      //calculate CCE corrected DFT formation enthalpies if precalculated DFT formation energies are provided
      for(uint k=0,ksize=cce_vars.vfunctionals.size();k<ksize;k++){ // looping over and checking of vfunctionals is necessary to ensure correct correspondence between given formation energies [k] and corrections with respect to the functional
        if (cce_vars.vfunctionals[k] != "exp") {
          if(cce_vars.dft_energies.size()!=0){ // only when precalculated DFT formation energies are provided, calculate corrected values
            // for 298.15 K
            cce_vars.cce_form_energy_cell[2*k] = cce_vars.dft_energies[k] - cce_vars.cce_correction[2*k] ;
            // for 0 K
            cce_vars.cce_form_energy_cell[2*k+1] = cce_vars.dft_energies[k] - cce_vars.cce_correction[2*k+1] ;
          }
        } else { // for CCE@exp the formation enthalpy (only at 298.15 K) is directly calculated from the exp. formation enthalpies per bond times the number of cation-anion bonds and is not subtracted from a given value
          cce_vars.cce_form_energy_cell[2*k] = cce_vars.cce_correction[2*k] ;
        }
      }

      if (aurostd::toupper(flags.getattachedscheme("CCE_CORRECTION::PRINT")) == "JSON") {
        std::cout << CCE_get_JSON(structure, cce_vars) << std::endl;
      } else {
        // write CCE corrections & corrected formation enthalpies per cell and atom
        std::cout << CCE_write_output(structure, cce_vars, cce_vars.cce_form_energy_cell);
        // write CCE citation
        std::cout << CCE_write_citation();
      }
    }
  } // main CCE function for command line use



  //CCE_correct////////////////////////////////////////////////////////
  // main CCE function for calling inside AFLOW providing only directory path where data neded for correction (POSCAR.static & aflow.in) are located
  // for setting parameters, analyzing structure, determining oxidation numbers, assigning corrections,
  // calculating total corrections, converting correction vector, and returning corrections
  vector<double> CCE_correct(string directory_path) { 
    string soliloquy="cce::CCE_correct():";
    stringstream message;
    // get structure from POSCAR.static in directory
    string poscar_path=directory_path + "/CONTCAR.relax2";
    xstructure structure=CCE_read_structure(poscar_path); // AFLOW seems to automatically unzip and rezip zipped files so that only the file name without zipping extension needs to be given
    // get functional from aflow.in in directory
    string aflowin_file= directory_path + "/" + _AFLOWIN_;
    string functional=CCE_get_functional_from_aflow_in(structure, aflowin_file);
    if (functional.empty()) {
      message << "Functional cannot be determined from aflow.in. Corrections are available for PBE, LDA, SCAN, or PBE+U_ICSD.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    return CCE_correct(structure, functional);
  } // main CCE function for calling inside AFLOW with directory path



  //CCE_correct////////////////////////////////////////////////////////
  // main CCE function for calling inside AFLOW providing only structure (and functional) assuming that oxidation numbers will be obtained from Bader charges or electronegativities
  // for setting parameters, analyzing structure, determining oxidation numbers, assigning corrections,
  // calculating total corrections, converting correction vector, and returning corrections
  //vector<double> CCE(xstructure& structure) // OLD: functional will be automatically determined during Bader charge analysis for the current implementation, later when using e.g. electronegativities, it might be needed as input
  vector<double> CCE_correct(xstructure& structure, string functional) { // functional needed as input when determining oxidation numbers from electronegativities
    string soliloquy="cce::CCE_correct():";
    stringstream message;
    CCE_Variables cce_vars = CCE_init_variables(structure);
    aurostd::xoption cce_flags = CCE_init_flags();
    // determine functional
    if(functional!="exp"){
      functional=aurostd::toupper(functional);
    }
    if (!aurostd::withinList(CCE_vallowed_functionals, functional) || CCE_get_offset(functional) == -1) {
      message << "Unknown functional " << functional << ". Please choose PBE, LDA, or SCAN.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    cce_vars.vfunctionals.push_back(functional);
    cce_vars.offset.push_back(CCE_get_offset(functional));
    // run correction
    CCE_core(structure, cce_vars, cce_flags);
    // cce_vars.cce_corrections can be returned directly since there is always only one functional for this CCE function
    return cce_vars.cce_correction;
  } // main CCE function for calling inside AFLOW



  //CCE///////////////////////////////////////////////////////////////////////
  // ME 200213
  // For poscar2cce
  void CCE(aurostd::xoption& flags, std::istream& ist) {
    // initialize
    aurostd::xoption cce_flags = CCE_init_flags();
    if (!(aurostd::toupper(flags.getattachedscheme("CCE_CORRECTION::PRINT")) == "JSON")) {
      cce_flags.flag("COMMAND_LINE", true);
    }

    xstructure structure(ist);
    structure.ReScale(1.0);
    CCE_Variables cce_vars = CCE_init_variables(structure);

    // Get DFT energies 
    CCE_get_dft_form_energies_functionals(flags.getattachedscheme("CCE_CORRECTION::DFT_FORMATION_ENERGIES"), flags.getattachedscheme("CCE_CORRECTION::FUNCTIONALS"), cce_vars); //provide precalc. DFT formation energies & corresponding functionals

    // Get oxidation states
    if(flags.flag("CCE_CORRECTION::OXIDATION_NUMBERS")) {
      cce_vars.oxidation_states = CCE_get_oxidation_states(flags.getattachedscheme("CCE_CORRECTION::OXIDATION_NUMBERS"), structure, cce_flags, cce_vars);
      cce_flags.flag("OX_NUMS_PROVIDED",TRUE);
    }

    // Run correction
    CCE_core(structure, cce_vars, cce_flags);
    if (cce_flags.flag("CORRECTABLE")) {
      cce_vars.cce_form_energy_cell = CCE_get_formation_enthalpies(cce_vars.cce_correction, cce_vars);

      if (aurostd::toupper(flags.getattachedscheme("CCE_CORRECTION::PRINT")) == "JSON") {
        std::cout << CCE_get_JSON(structure, cce_vars) << std::endl;
      } else {
        // write CCE corrections & corrected formation enthalpies per cell and atom
        std::cout << CCE_write_output(structure, cce_vars, cce_vars.cce_form_energy_cell);
        // write CCE citation
        std::cout << CCE_write_citation();
      }
    }
  }



  //CCE_core////////////////////////////////////////////////////////
  // main CCE function core
  // analyzing structure, determining oxidation numbers, assigning corrections,
  // calculating total corrections, and returning full correction vector
  void CCE_core(xstructure& structure, CCE_Variables& cce_vars, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    // corrections_atom, (su-)perox_correction, cce_correction and cce_form_energy_cell can only be resized after vfunctionals.size() is known from previous main CCE functions calling CCE_core
    // vfunctionals.size()*2 for 298.15 and 0K corrections
    cce_vars.corrections_atom.clear();
    cce_vars.corrections_atom.resize(cce_vars.vfunctionals.size()*2, vector<double>(structure.atoms.size()));
    cce_vars.perox_correction.clear();
    cce_vars.perox_correction.resize(cce_vars.vfunctionals.size()*2);
    cce_vars.superox_correction.clear();
    cce_vars.superox_correction.resize(cce_vars.vfunctionals.size()*2);
    cce_vars.cce_correction.clear();
    cce_vars.cce_correction.resize(cce_vars.vfunctionals.size()*2);
    cce_vars.cce_form_energy_cell.clear();
    cce_vars.cce_form_energy_cell.resize(cce_vars.vfunctionals.size()*2);

    // DETERMINE NUMBER OF NEAREST ANION NEIGHBORS FOR EACH CATION: STRUCTURAL PART OF CORRECTION ##
    /********************************************************/
    // Determine anion species
    /********************************************************/
    // determine which species is the anion (for single anion systems)
    cce_vars.anion_species=CCE_determine_anion_species(structure, cce_vars);

    /********************************************************/
    // Check for multi anion system
    /********************************************************/
    // check whether it is a multi-anion system and set variables accordingly
    cce_flags.flag("MULTI_ANION_SYSTEM",FALSE);
    cce_flags.flag("O_MULTI_ANION_SPECIES",FALSE); // whether one of the multi anion species is O for which per/superoxide tests need to be made
    cce_vars.multi_anion_atoms=CCE_check_for_multi_anion_system(cce_vars, _CCE_NN_DIST_TOL_MULTI_ANION_, structure, cce_flags);
    // multi anion corrections can only be resized after number of multi anion species is known from check for multi anion system
    cce_vars.multi_anion_corrections_atom.clear();
    cce_vars.multi_anion_corrections_atom.resize(cce_vars.multi_anion_species.size(), vector<vector<double> >(cce_vars.vfunctionals.size()*2, vector<double>(structure.atoms.size())));

    /********************************************************/
    // Determine anion nearest neighbors for each cation
    /********************************************************/
    // apply species selective cutoffs to determine only nearest neighbors within respective cutoff
    cce_vars.num_neighbors=CCE_get_num_neighbors(cce_vars.anion_species, _CCE_NN_DIST_TOL_, structure, cce_flags, cce_vars);

    // determine anion nearest neighbors for cations bound to multi anion atoms if needed
    vector<vector<uint> > multi_anion_num_neighbors(cce_vars.multi_anion_species.size(), vector<uint>(structure.atoms.size()));
    if (cce_flags.flag("MULTI_ANION_SYSTEM")){
      for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
        if(LDEBUG){
          cerr  << "getting neighbors for multi anion species " << k << " (" << cce_vars.multi_anion_species[k] << ")" << endl;
        }
        multi_anion_num_neighbors[k]=CCE_get_num_neighbors(cce_vars.multi_anion_species[k], _CCE_NN_DIST_TOL_, structure, cce_flags, cce_vars);
      }
    }

    /********************************************************/
    // Check for per- and superoxide ions
    /********************************************************/
    // check for per- and superoxide ions based on O-O bond length and set number of (su-)peroxide bonds 
    // and indices accordingly if ions are detected
    if(cce_vars.anion_species == "O" || cce_flags.flag("O_MULTI_ANION_SPECIES")) {
      CCE_check_per_super_oxides(structure, cce_vars, cce_flags);
    }


    // DETERMINE CORRECTION FOR EACH CATION: OXIDATION STATE DEPENDENT PART OF COORECTION ##########
    /********************************************************/
    // Assign oxidation numbers
    /********************************************************/
    // determine oxidation numbers automatically from structure and Allen electronegativities if not provided on command line
    if(!cce_flags.flag("OX_NUMS_PROVIDED")) {
      cce_vars.oxidation_states=CCE_get_oxidation_states_from_electronegativities(structure, cce_vars, cce_flags);
      if(0) { // obtaining oxidation states from Bader charges is outdated but the functionality is kept mainly for test purposes
        cce_vars.oxidation_states=CCE_get_oxidation_states_from_Bader(cce_flags, structure, cce_vars);
      }
    }

    /********************************************************/
    // Assign corrections
    /********************************************************/
    // assigning corrections for each atom after oxidation numbers are determined if the system is correctable
    if (cce_flags.flag("CORRECTABLE")){
      CCE_get_corrections(cce_flags, structure, cce_vars.corrections_atom, cce_vars, cce_vars.num_neighbors, cce_vars.anion_species);
      // determine corrections for cations bound to multi anion atoms if needed
      if (cce_flags.flag("MULTI_ANION_SYSTEM")){
        for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
          if(LDEBUG){
            cerr  << endl;
            cerr  << "getting corrections for multi anion species " << k << " (" << cce_vars.multi_anion_species[k] << ")" << endl;
          }
          CCE_get_corrections(cce_flags, structure, cce_vars.multi_anion_corrections_atom[k], cce_vars, multi_anion_num_neighbors[k], cce_vars.multi_anion_species[k]);
        }
      }
      // load per- & superox. corrections if needed;
      if (cce_vars.num_perox_bonds > 0 || cce_vars.num_superox_bonds > 0){
        CCE_check_get_per_super_ox_corrections(cce_vars);
      }
    }


    // CALCULATE CCE CORRECTIONS AT 298.15 AND 0K ##################################################
    // calculate CCE correction if the system is correctable
    if (cce_flags.flag("CORRECTABLE")){
      // calculate total correction per cell only using cations although 
      // corrections for anions should be set to zero
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){
        if (structure.atoms[i].name != cce_vars.anion_species && cce_vars.multi_anion_atoms[i] != 1){ // exclude main anion species and multi anion atoms detected previously
          if (cce_vars.num_neighbors[i] > 0){ // are there actually bonds between the cation and the (main) anion species
            for (uint k = 0; k < cce_vars.vfunctionals.size(); k++) {
              // for exp this calculates already the CCE@exp formation enthalpy (only at 298.15 K) from the exp. formation enthalpies per bond
              // for 298.15 K
              cce_vars.cce_correction[2*k] += (cce_vars.num_neighbors[i]*cce_vars.corrections_atom[2*k][i]) ; // 2*k since to each functional belong 2 corrections for 298.15 and 0K
              if (cce_vars.vfunctionals[k] != "exp") {
                // for 0 K
                cce_vars.cce_correction[2*k+1] += (cce_vars.num_neighbors[i]*cce_vars.corrections_atom[2*k+1][i]) ;
              }
            }
          }
        }
      }
      // add correction for cations bound to multi anion atoms if needed
      if (cce_flags.flag("MULTI_ANION_SYSTEM")){
        for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
          if(LDEBUG){
            cerr  << endl;
            cerr  << "adding corrections for multi anion species " << k << " (" << cce_vars.multi_anion_species[k] << ")" << endl;
          }
          for(uint i=0,isize=structure.atoms.size();i<isize;i++){
            if (structure.atoms[i].name != cce_vars.anion_species && cce_vars.multi_anion_atoms[i] != 1){ // exclude main anion species and multi anion atoms detected previously
              if (multi_anion_num_neighbors[k][i] > 0){ // are there actually bonds between the cation and the anion species under consideration
                for (uint l = 0; l < cce_vars.vfunctionals.size(); l++) {
                  // for exp this calculates already the CCE@exp formation enthalpy from the exp. formation enthalpies per bond
                  // for 298.15 K
                  cce_vars.cce_correction[2*l] += (multi_anion_num_neighbors[k][i]*cce_vars.multi_anion_corrections_atom[k][2*l][i]) ;
                  if (cce_vars.vfunctionals[l] != "exp") {
                    // for 0 K
                    cce_vars.cce_correction[2*l+1] += (multi_anion_num_neighbors[k][i]*cce_vars.multi_anion_corrections_atom[k][2*l+1][i]) ;
                  }
                }
              }
            }
          }
        }
      }
      // add per- and superoxide corrections if needed
      if (cce_vars.num_perox_bonds > 0 || cce_vars.num_superox_bonds > 0){
        CCE_check_apply_per_super_ox_corrections(cce_vars);
      }
    }
  } // main CCE function core



  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                   READ USER INPUT (FROM COMMAND LINE)                   //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  //CCE_read_structure////////////////////////////////////////////////////////
  // read structural data from structure file provided on command line
  xstructure CCE_read_structure(const string& structure_file, int mode){ // first argument can be directly structure_file and not structure_file_path
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::CCE_read_structure():";
    stringstream message;
    //string structure_file = aurostd::file2string(structure_file_path); // first argument of CCE_read_structure_function does not need to be converted to string since it contains already the file content and not only the file name
    xstructure structure(structure_file, mode);
    structure.ReScale(1.0); // rescales scaling factor in second line of POSCAR to 1, needed for correct distances
    //let the program spit out what it thinks (input structure)
    if(LDEBUG){
      cerr << soliloquy << endl << "INPUT STRUCTURE:" << endl;
      cerr << structure << endl;
    }
    // check whether there are any atoms in the structure
    if (structure.atoms.size() == 0){
      message << "BAD NEWS: It seems there are no atoms in the structure file. Please adjust the structure file and rerun.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
    }
    // if species of atoms are not known like in VASP4 format, throw error
    for(uint k=0,ksize=structure.atoms.size();k<ksize;k++){
      if (structure.atoms[k].name == ""){
        message << "BAD NEWS: It seems you are providing a POSCAR without complete species information as input. This implementation requires a POSCAR in VASP5 format with the species information included. Please adjust the structure file and rerun.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
      }
    }
    // if there is only one species, it must be an elemental phase and is hence not correctable
    if (structure.species.size() == 1){
      message << "BAD NEWS: Only one species found. Enthalpies of elemental systems cannot be corrected with the CCE methodology.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
    }
    return structure;
  }

  //CCE_get_dft_form_energies_functionals////////////////////////////////////////////////////////
  // set the DFT formation energies and functionals according to the input and check consistency of functionals and formation energies
  void CCE_get_dft_form_energies_functionals(const string& dft_energies_input_str, const string& functionals_input_str, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::CCE_get_dft_form_energies_functionals():";
    stringstream message;
    // vectorize precalculated DFT formation energies if provided
    if(!dft_energies_input_str.empty()){ //if the DFT energies input string is not empty, convert it to vector of doubles
      aurostd::string2tokens<double>(dft_energies_input_str,cce_vars.dft_energies,",");
    } 
    // vectorize corresponding functionals
    if(!functionals_input_str.empty()){ // if functionals input string is not empty, convert it to vector of strings
      aurostd::string2tokens(functionals_input_str,cce_vars.vfunctionals,",");
    }
    //use PBE as default if only 1 DFT formation energy is provided
    ostream& oss = cout;
    ofstream FileMESSAGE;
    _aflags aflags;aflags.Directory=".";
    if(functionals_input_str.empty() && cce_vars.dft_energies.size() == 1){
      message << "Setting functionals=PBE since only 1 DFT formation energy is provided and PBE is the default functional!";
      pflow::logger(_AFLOW_FILE_NAME_,soliloquy, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
      cce_vars.vfunctionals.clear();
      cce_vars.vfunctionals.push_back("PBE");
    // otherwise if sizes of provided DFT formation energies and functionals do not match, throw error
    // if only functional argument is set corrections should only be returned for desired functional
    } else if(cce_vars.dft_energies.size()!=cce_vars.vfunctionals.size() && !dft_energies_input_str.empty() ){ 
      message << "BAD NEWS: The number of provided precalculated DFT formation energies and functionals must match.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
    }
    //let the program spit out what it thinks
    if(LDEBUG){
      bool roff = false;
      cerr << "INPUT DFT FORMATION ENERGIES & FUNCTIONALS:" << endl;
      cerr << soliloquy << " input dft formation energies=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(cce_vars.dft_energies,6,roff),",") << " (assumed to be in eV/cell)" << endl;
      cerr << soliloquy << " input functionals=" << functionals_input_str << endl;
    }
    // determine whether it is a functional for which corrections are available
    for(uint k=0,ksize=cce_vars.vfunctionals.size();k<ksize;k++){
      if(cce_vars.vfunctionals[k]!="exp"){
        cce_vars.vfunctionals[k]=aurostd::toupper(cce_vars.vfunctionals[k]);
      }
      if(LDEBUG){
        cerr << "cce_vars.vfunctionals[" << k << "]: " << cce_vars.vfunctionals[k] << endl;
      }
      if (!aurostd::withinList(CCE_vallowed_functionals, cce_vars.vfunctionals[k]) || CCE_get_offset(cce_vars.vfunctionals[k]) == -1) {
        message << "Unknown functional " << cce_vars.vfunctionals[k] << ". Please choose PBE, LDA, or SCAN.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
      }
      cce_vars.offset.push_back(CCE_get_offset(cce_vars.vfunctionals[k]));
    }
    // if only structure (and oxidation numbers) are provided, corrections should be given for PBE, LDA, and SCAN 
    // and the CCE@exp formation enthalpy from the exp. formation enthalpies per bond
    // ICSD correction should only be returned when explicitly asked for
    if(functionals_input_str.empty() && dft_energies_input_str.empty()){
      for(uint k=0,ksize=CCE_vdefault_output_functionals.size();k<ksize;k++){
        if (CCE_get_offset(CCE_vdefault_output_functionals[k]) == -1) {
          message << "Unknown functional " << cce_vars.vfunctionals[k] << ". Please choose PBE, LDA, or SCAN.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
        }
        cce_vars.vfunctionals.push_back(CCE_vdefault_output_functionals[k]); cce_vars.offset.push_back(CCE_get_offset(CCE_vdefault_output_functionals[k]));
      }
    }
    if(LDEBUG){
      cerr << "PBE: " << aurostd::withinList(cce_vars.vfunctionals, "PBE") << endl;
      cerr << "LDA: " << aurostd::withinList(cce_vars.vfunctionals, "LDA") << endl;
      cerr << "SCAN: " << aurostd::withinList(cce_vars.vfunctionals, "SCAN") << endl;
      cerr << "PBE+U_ICSD: " << aurostd::withinList(cce_vars.vfunctionals, "PBE+U_ICSD") << endl;
      cerr << "exp: " << aurostd::withinList(cce_vars.vfunctionals, "exp") << endl;
      cerr << endl;
    }
  }

  //CCE_get_offset////////////////////////////////////////////////////////
  // get offset needed for reading corrections from lookup table for different functionals
  int CCE_get_offset(const string& functional) {
  if (functional=="PBE")          {return 0;}
  if (functional=="LDA")          {return 2;}
  if (functional=="SCAN")         {return 4;}
  if (functional=="PBE+U_ICSD")   {return 6;}
  if (functional=="exp")          {return 8;}
  else {return -1;}
  }

  //CCE_get_oxidation_states////////////////////////////////////////////////////////
  // Retrieves the oxidation states of the material.
  vector<double> CCE_get_oxidation_states(const string& oxidation_numbers_input_str, const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars) {
    string soliloquy="cce::CCE_get_oxidation_states():";
    stringstream message;
    if(!oxidation_numbers_input_str.empty()){
      aurostd::string2tokens<double>(oxidation_numbers_input_str,cce_vars.oxidation_states,","); //if the oxidation numbers input string is not empty, convert it to vector of doubles
      //sizes of oxidation numbers and atoms must match
      if(cce_vars.oxidation_states.size()!=structure.atoms.size()){ 
        message << "BAD NEWS: The number of provided oxidation numbers does not match the number of atoms in the structure! Please correct and rerun.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
      }
      // print oxidation numbers
      if(cce_flags.flag("COMMAND_LINE")){
        cout << endl;
        cout << "INPUT OXIDATION NUMBERS:" << endl;
        for (uint k=0,ksize=cce_vars.oxidation_states.size();k<ksize;k++){
          cout << "Oxidation state of " << structure.atoms[k].name << " (atom[" << k << "]): " << cce_vars.oxidation_states[k] << endl;
        }
      }
      // get sum of oxidation numbers and validate (system should not be regarded correctable if sum over oxidation states is not zero)
      cce_vars.oxidation_sum = CCE_get_oxidation_states_sum(cce_vars, cce_flags); // double because for superoxides O ox. number is -0.5
      if(cce_flags.flag("COMMAND_LINE")){
        cout << endl;
      }
      if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
        string function = "cce::CCE_get_oxidation_states()";
        message << "BAD NEWS: The formation energy of this system is not correctable!"
        << " The oxidation numbers that you provided do not add up to zero!"
        << " Please correct and rerun.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _INPUT_ILLEGAL_);	
      }
    } else {
      message << "It seems you forgot to provide the oxidation numbers after \"--oxidation_numbers=\". Please add them or omit the option.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
    }
    return cce_vars.oxidation_states;
  }

  //CCE_get_functional_from_aflow_in////////////////////////////////////////////////////////
  // determine the functional from the aflow_in
  string CCE_get_functional_from_aflow_in(const xstructure& structure, string& aflowin_file) {
    string soliloquy="cce::CCE_get_functional_from_aflow_in():";
    stringstream message;
    string functional = "";
    string aflowIn = aurostd::RemoveComments(aurostd::file2string(aflowin_file));
    vector<string> vlines = aurostd::string2vectorstring(aflowIn);
    string line_a;
    bool ldau = false;
    bool ldau2 = false;
    // check whether it is a PBE calculation; PBE+U will be checked later
    for (uint i = 0; i < vlines.size(); i++) {
      line_a = aurostd::RemoveSpaces(vlines[i]);
      if (line_a.find("=potpaw_PBE") != string::npos || line_a.find("/potpaw_PBE") != string::npos){
        functional = "PBE";
      }
    }
    // check whether it is a DFT+U calculation with parameters as for the ICSD (PBE+U_ICSD calculation)
    // first check whether it is an LDAU2 calculation
    for (uint i = 0; i < vlines.size(); i++) { 
      line_a = aurostd::RemoveSpaces(vlines[i]);
      if (line_a.find("LDAU2=ON") != string::npos){ // string::npos is returned if string is not found
        ldau2 = true;
      }
    }
    // then read and check U parameters
    for (uint i = 0; i < vlines.size(); i++) { 
      line_a = aurostd::RemoveSpaces(vlines[i]);
      // new implementation checking Us explicitly
      if (line_a.find("LDAU_PARAMETERS=") != string::npos && ldau2){ // string::npos is returned if string is not found
        vector<string> ldau_line_vector;
        aurostd::string2tokens(line_a, ldau_line_vector, ";"); 
        // get species
        string species_part_1 = ldau_line_vector[0];
        vector<string> species_part_vector;
        aurostd::string2tokens(species_part_1, species_part_vector, "="); 
        string species_part = species_part_vector[1];
        vector<string> species_vector;
        aurostd::string2tokens(species_part, species_vector, ","); 
        if (species_vector.size() != structure.species.size()){
          message << "BAD NEWS: The number of species in the DFT+U settings differs from the total number of species for this structure. Please adapt and rerun.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
        }
        // get Us
        string Us_part = ldau_line_vector[2];
        vector<string> Us_string_vector;
        vector<double> Us_vector;
        aurostd::string2tokens(Us_part, Us_string_vector, ","); 
        for (uint k = 0; k < Us_string_vector.size(); k++) { 
          Us_vector.push_back(aurostd::string2utype<double>(Us_string_vector[k]));
        }
        if (species_vector.size() != Us_vector.size()){
          message << "BAD NEWS: The number of species in the DFT+U settings differs from the number of provided U values. Please adapt and rerun.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
        }
        // get standard Us for PBE+U_ICSD
        vector<double> standard_ICSD_Us_vector;
        bool LDAU=true;
        vector<string> vLDAUspecies;
        vector<uint> vLDAUtype;
        vector<int> vLDAUL;
        vector<double> vLDAUU;
        vector<double> vLDAUJ;
        for (uint k = 0; k < species_vector.size(); k++) { 
          AVASP_Get_LDAU_Parameters(species_vector[k],LDAU,vLDAUspecies,vLDAUtype,vLDAUL,vLDAUU,vLDAUJ);
          if (vLDAUtype[k] != 2 && vLDAUU[k] !=0){
            ostream& oss = cerr;
            ofstream FileMESSAGE;
            _aflags aflags;aflags.Directory=".";
            message << "BAD NEWS: It seems the standard DFT+U method for " << species_vector[k] << " was changed from Dudarev's approach (LDAU2=ON, vLDAUtype[" << k << "]=2) as used for the AFLOW ICSD database when obtaining the corrections to vLDAUtype[" << k << "]=" << vLDAUtype[k] << " now. If the standard U values have also been changed, then the corrections might have been constructed for other U values than used in this calculation and should not be applied. Please check this carefully!";
            pflow::logger(_AFLOW_FILE_NAME_,soliloquy, message, aflags, FileMESSAGE, oss, _LOGGER_WARNING_);
          }
          standard_ICSD_Us_vector.push_back(vLDAUU[k]);
        }
        // compare read and standard Us
        bool Us_disagree=false;
        for (uint k = 0; k < species_vector.size(); k++) { 
          if (Us_vector[k] != standard_ICSD_Us_vector[k]){
            Us_disagree=true;
            message << "BAD NEWS: For this DFT+U calculation with Dudarev's method the provided U value of " << Us_vector[k] << " eV for " << species_vector[k] << " does not match the standard value of " << standard_ICSD_Us_vector[k] << " eV. There are no corrections for this case.";
            throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
          }
        }
        if (!Us_disagree){
          functional = "PBE+U_ICSD";
        }
      }
    }
    // check whether it is a DFT+U calculation with different parameters than for PBE+U_ICSD
    for (uint i = 0; i < vlines.size(); i++) { 
      line_a = aurostd::RemoveSpaces(vlines[i]);
      if ((line_a.find("LDAU") != string::npos && line_a.find("=ON") != string::npos) || (line_a.find("LDAU") != string::npos && line_a.find("=ADIABATIC") != string::npos)){
        ldau = true;
        if (functional != "PBE+U_ICSD"){
          message << "BAD NEWS: It seems you are providing an aflow.in for a DFT+U calculation with different parameters than for the AFLOW ICSD database (Dudarev's approach, LDAU2=ON). There are no corrections for this case.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
        }
      }
    }
    // check whether it is an LDA or SCAN calculation
    if (functional != "PBE+U_ICSD"){
      if (!ldau){
        for (uint i = 0; i < vlines.size(); i++) {
          line_a = aurostd::RemoveSpaces(vlines[i]);
          if (line_a.find("=potpaw_LDA") != string::npos || line_a.find("/potpaw_LDA/") != string::npos){ // the first criterion seems to not find "=potpaw_LDA" when the line is commented by # due to RemoveComments above
            functional = "LDA";
            for (uint k = 0; k < vlines.size(); k++) { // it could still be a SCAN calculation with LDA PPs; allowing for that since to my experience right now the SCAN corrections from calculations with PBE PPs should be good for SCAN calculations with LDA PPs
              string line_b = aurostd::RemoveSpaces(vlines[k]);
              if (line_b.find("METAGGA=SCAN") != string::npos){
                functional = "SCAN";
              }
            }
          } else if (line_a.find("METAGGA=SCAN") != string::npos){
            functional = "SCAN";
          }
        }
      }
    }
    return functional;
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                      INITIALISE FLAGS AND VARIABLES                     //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  //CCE_init_flags////////////////////////////////////////////////////////////
  // ME 200213
  // Initializes the CCE flags to their default values.
  aurostd::xoption CCE_init_flags() {
    aurostd::xoption flags;
    flags.flag("COMMAND_LINE", false);
    flags.flag("CORRECTABLE", true);
    flags.flag("OX_NUMS_PROVIDED", false);
    flags.flag("MULTI_ANION_SYSTEM", false);
    flags.flag("O_MULTI_ANION_SPECIES", false); // whether one of the multi anion species is O for which per/superoxide tests need to be made
    return flags;
  }

  //CCE_init_variables////////////////////////////////////////////////////////
  // ME 200213
  // Initializes the CCE_variables struct.
  CCE_Variables CCE_init_variables(const xstructure& structure) {
    CCE_Variables cce_vars;
    cce_vars.dft_energies.clear();
    cce_vars.vfunctionals.clear();
    cce_vars.offset.clear();
    cce_vars.num_perox_bonds=0;
    cce_vars.num_superox_bonds=0;
    // clear and resize all vectors
    cce_vars.electronegativities.clear();
    cce_vars.electronegativities.resize(structure.species.size());
    cce_vars.multi_anion_atoms.clear();
    cce_vars.multi_anion_atoms.resize(structure.atoms.size());
    cce_vars.oxidation_states.clear();
    cce_vars.oxidation_states.resize(structure.atoms.size());
    cce_vars.multi_anion_species.clear();
    cce_vars.perox_indices.clear();
    cce_vars.perox_indices.resize(structure.atoms.size());
    cce_vars.superox_indices.clear();
    cce_vars.superox_indices.resize(structure.atoms.size());
    cce_vars.num_neighbors.clear();
    cce_vars.num_neighbors.resize(structure.atoms.size());
    cce_vars.num_pref_ox_states.clear();
    cce_vars.num_pref_ox_states.resize(structure.species.size());
    cce_vars.pref_ox_states_strings.clear();
    cce_vars.pref_ox_states_strings.resize(structure.species.size());
    cce_vars.num_all_ox_states.clear();
    cce_vars.num_all_ox_states.resize(structure.species.size());
    cce_vars.all_ox_states_strings.clear();
    cce_vars.all_ox_states_strings.resize(structure.species.size());
    cce_vars.species_electronegativity_sorted.clear();
    cce_vars.species_electronegativity_sorted.resize(structure.species.size());
    cce_vars.num_pref_ox_states_electronegativity_sorted.clear();
    cce_vars.num_pref_ox_states_electronegativity_sorted.resize(structure.species.size());
    cce_vars.pref_ox_states_strings_electronegativity_sorted.clear();
    cce_vars.pref_ox_states_strings_electronegativity_sorted.resize(structure.species.size());
    cce_vars.num_all_ox_states_electronegativity_sorted.clear();
    cce_vars.num_all_ox_states_electronegativity_sorted.resize(structure.species.size());
    cce_vars.all_ox_states_strings_electronegativity_sorted.clear();
    cce_vars.all_ox_states_strings_electronegativity_sorted.resize(structure.species.size());
    cce_vars.cations_map.clear();
    cce_vars.Bader_charges.clear();
    return cce_vars;
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                           STRUCTURAL ANALYSIS                           //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // counting the number of anion neighbours (multiple anion species possible) for each cation

  //CCE_determine_anion_species////////////////////////////////////////////////////////
  // determine anion species
  string CCE_determine_anion_species(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::CCE_determine_anion_species():";
    stringstream message;
    if(LDEBUG){
      cerr << "ANION SPECIES FROM ALLEN ELECTRONEGATIVITIES:" << endl;
    }
    double anion_electronegativity = 0;
    for(uint k=0,ksize=structure.species.size();k<ksize;k++){
      if (CCE_get_electronegativities_ox_nums(structure.species[k]) == "") {
        message << "VERY BAD NEWS: There is no known electronegativity value for " << structure.species[k] << ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
      } else{
        string electroneg_line = CCE_get_electronegativities_ox_nums(structure.species[k]);
        if(LDEBUG){
          cerr << "electronegativity line for species " << k << ": " << electroneg_line << endl;
        }
        vector<string> electroneg_tokens;
        aurostd::string2tokens(electroneg_line, electroneg_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one
        cce_vars.electronegativities[k] = aurostd::string2utype<double>(electroneg_tokens[0]);
        if(LDEBUG){
          cerr << "electronegativity of species " << k << " (" << structure.species[k] << "): " << cce_vars.electronegativities[k] << endl;
        }
        if (cce_vars.electronegativities[k] > anion_electronegativity) {
          anion_electronegativity = cce_vars.electronegativities[k];
          cce_vars.anion_species = structure.species[k];
        }
      }
    }
    // set anion charge and check whether it is negative
    string ox_nums_line = CCE_get_electronegativities_ox_nums(cce_vars.anion_species);
    vector<string> ox_nums_tokens_1;
    vector<string> ox_nums_tokens_2;
    aurostd::string2tokens(ox_nums_line, ox_nums_tokens_1, " "); // anion charge should be among all known oxidation states (last element of the line separated by spaces)
    aurostd::string2tokens(ox_nums_tokens_1.back(), ox_nums_tokens_2, ","); // anion charge should be last (most negative) oxidation number separated by commas
    cce_vars.standard_anion_charge = aurostd::string2utype<double>(ox_nums_tokens_2[ox_nums_tokens_2.size()-1]);
    if (cce_vars.standard_anion_charge > 0) {
      message << "VERY BAD NEWS: There is no known negative oxidation number for " << cce_vars.anion_species << " detected as anion species.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
    }
    if(LDEBUG){
      cerr << "anion electronegativity: " << anion_electronegativity << endl;
      cerr << "anion species: " << cce_vars.anion_species << endl;
      cerr << "anion charge: " << cce_vars.standard_anion_charge << endl;
      cerr << endl;
    }
    return cce_vars.anion_species;
  }

  //CCE_check_for_multi_anion_system////////////////////////////////////////////////////////
  // check whether it is a multi-anion system, i. e. whether atoms of another species than the the anion_species are only bound to atoms of lower electronegativity or of the same type
  vector<uint> CCE_check_for_multi_anion_system(CCE_Variables& cce_vars, double tolerance, xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::CCE_check_for_multi_anion_system():";
    stringstream message;
    if(LDEBUG){
      cerr << "CHECKING FOR MULTI-ANION SYSTEM:" << endl;
    }
    for ( uint i = 0; i < structure.atoms.size(); i++ ) { // initialize elements of vector to 0 
      cce_vars.multi_anion_atoms[i] = 0; 
    }
    cce_vars.cutoffs=CCE_get_dist_cutoffs(tolerance, structure);
    double cutoffs_max=aurostd::max(cce_vars.cutoffs);
    deque<deque<_atom> > neigh_mat;
    structure.GetStrNeighData(cutoffs_max,neigh_mat);
    for(uint i=0,isize=neigh_mat.size();i<isize;i++){ //same size as structure.atoms.size(); number of atoms in the structure (not determined by cutoff (or cutoffs_max))
      uint neighbors_count=0;
      uint multi_anion_count=0;
      for(uint j=0,jsize=neigh_mat[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by cutoffs_max
        const _atom& atom=neigh_mat[i][j]; // the atom object stands for the neighbors of each atom of the structure
        if (_CCE_SELF_DIST_TOL_ < AtomDist(structure.atoms[i],atom) && AtomDist(structure.atoms[i],atom) <= cce_vars.cutoffs[structure.atoms[i].type] ){ // distance must be larger than _CCE_SELF_DIST_TOL_ since GetStrNeighData includes also structure.atoms[i] itself as neighbor having distance zero to itself
          if (atom.name == cce_vars.anion_species){
            neighbors_count+=1;
          } else if (atom.name != cce_vars.anion_species && structure.atoms[i].name != cce_vars.anion_species){ // second condition set since the anion_species cannot be set as a multi-anion species again
            neighbors_count+=1;
            string electroneg_line_atom = CCE_get_electronegativities_ox_nums(structure.atoms[i].name);
            string electroneg_line_neighbor = CCE_get_electronegativities_ox_nums(atom.name);
            vector<string> electroneg_tokens;
            aurostd::string2tokens(electroneg_line_atom, electroneg_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one
            double electronegativity_atom = aurostd::string2utype<double>(electroneg_tokens[0]);
            aurostd::string2tokens(electroneg_line_neighbor, electroneg_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one
            double electronegativity_neighbor = aurostd::string2utype<double>(electroneg_tokens[0]);
            if(LDEBUG){
              cerr << "electronegativity of atom " << i << ": " << electronegativity_atom << endl;
              cerr << "electronegativity of neighbor " << j << ": " << electronegativity_neighbor << endl;
            }
            if(electronegativity_neighbor < electronegativity_atom || structure.atoms[i].name == atom.name){ // could be multi-anion atom if bound to only atoms of lower electroneg. or of same species
              multi_anion_count+=1;
            }
          }
        }
      }
      if(LDEBUG){
        cerr << "multi_anion_count for atom " << i << " (" << structure.atoms[i].name << "): " << multi_anion_count << endl;
        cerr << "neighbors_count for atom " << i << " (" << structure.atoms[i].name << "): " << neighbors_count << endl;
      }
      if (multi_anion_count == neighbors_count && structure.atoms[i].name != cce_vars.anion_species){ // anion_species should not be detected again as multi_anion_species
        if (!cce_flags.flag("MULTI_ANION_SYSTEM")){
          cce_flags.flag("MULTI_ANION_SYSTEM",TRUE);
          if(cce_flags.flag("COMMAND_LINE")){
            cout << endl;
            cout << "This system has been detected to be a multi-anion compound!" << endl;
          }
        }
        // set multi anion atoms
        cce_vars.multi_anion_atoms[i]=1;
        if(LDEBUG){
          cerr << "Atom " << i << " (" << structure.atoms[i].name << ") has been detected as a multi-anion atom." << endl;
        }
        // set multi anion species
        uint multi_anion_species_count=0;
        for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ // have to make sure that no anion species is added twice
          if(cce_vars.multi_anion_species[k] != structure.atoms[i].name){
            multi_anion_species_count+=1;
          }
        }
        if(multi_anion_species_count == cce_vars.multi_anion_species.size()){
          cce_vars.multi_anion_species.push_back(structure.atoms[i].name);
          if(LDEBUG){
            cerr << "New multi anion species: " << structure.atoms[i].name << endl;
          }
          if(structure.atoms[i].name == "O"){ // check whether one of the multi anion species is O for which per/superoxide test needs to be made
            cce_flags.flag("O_MULTI_ANION_SPECIES",TRUE);
            if(LDEBUG){
              cerr << "Oxygen was detected as multi anion species, i.e. system needs to be tested for (su-)peroxide ions." << endl;
            }
          }
        }
        // set multi anion oxidation numbers and check whether it is negative
        string ox_nums_line = CCE_get_electronegativities_ox_nums(structure.atoms[i].name);
        vector<string> ox_nums_tokens_1;
        vector<string> ox_nums_tokens_2;
        aurostd::string2tokens(ox_nums_line, ox_nums_tokens_1, " "); // anion charge should be among all known oxidation states (last element of the line separated by spaces)
        aurostd::string2tokens(ox_nums_tokens_1.back(), ox_nums_tokens_2, ","); // anion charge should be last (most negative) oxidation number separated by commas
        cce_vars.oxidation_states[i] = aurostd::string2utype<double>(ox_nums_tokens_2.back());
        if(LDEBUG){
          cerr << "Oxidation state for atom " << i << " (" << structure.atoms[i].name << ") has been set to: " << cce_vars.oxidation_states[i] << endl;
        }
        if (cce_vars.oxidation_states[i] > 0) {
          message << "VERY BAD NEWS: There is no known negative oxidation number for " << structure.atoms[i].name << " detected as multi anion species.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
        }
      }
    }
    if(LDEBUG){
      cerr << "number of multi anion species (1-total_number_anion_species since the main anion species with largest electronegativity is NOT counted as multi anion species): " << cce_vars.multi_anion_species.size() << endl;
      for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
        cout << "multi anion species " << k << ": " << cce_vars.multi_anion_species[k] << endl;
      }
      cout << endl;
    }
    return cce_vars.multi_anion_atoms;
  }

  //CCE_get_num_neighbors////////////////////////////////////////////////////////
  // trying to get the number of neighbors for each atom within respective species selective cutoff
  vector<uint> CCE_get_num_neighbors(double tolerance, xstructure& structure) {
    CCE_Variables cce_vars = CCE_init_variables(structure);
    string anion_species="";
    aurostd::xoption cce_flags = CCE_init_flags();
    return CCE_get_num_neighbors(anion_species, tolerance, structure, cce_flags, cce_vars);
  }

  //CCE_get_num_neighbors////////////////////////////////////////////////////////
  // trying to get the number of neighbors of a special anion_species for each atom within respective species selective cutoff
  vector<uint> CCE_get_num_neighbors(const string& anion_species, double tolerance, xstructure& structure) {
    CCE_Variables cce_vars = CCE_init_variables(structure);
    aurostd::xoption cce_flags = CCE_init_flags();
    return CCE_get_num_neighbors(anion_species, tolerance, structure, cce_flags, cce_vars);
  }

  //CCE_get_num_neighbors////////////////////////////////////////////////////////
  // trying to get the number of neighbors for each atom within respective species selective cutoff
  vector<uint> CCE_get_num_neighbors(const string& anion_species, double tolerance, xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars) { // anion_species here cannot be taken from cce_vars since function is also used to determine multi anion num_neighbors for which anion_species is the respective multi_anion_species
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << "STRUCTURAL ANALYSIS:" << endl;
    }
    cce_vars.cutoffs=CCE_get_dist_cutoffs(tolerance, structure);
    double cutoffs_max=aurostd::max(cce_vars.cutoffs);
    vector<uint> num_neighbors(structure.atoms.size());
    deque<deque<_atom> > neigh_mat;
    structure.GetStrNeighData(cutoffs_max,neigh_mat);
    for(uint i=0,isize=neigh_mat.size();i<isize;i++){ //same size as structure.atoms.size(); number of atoms in the structure (not determined by cutoff (or cutoffs_max))
      uint neighbors_count=0;
      bool warning = false;
      uint empty_line_count=0;
      for(uint j=0,jsize=neigh_mat[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by cutoffs_max
        const _atom& atom=neigh_mat[i][j];
        if (_CCE_SELF_DIST_TOL_ < AtomDist(structure.atoms[i],atom) && AtomDist(structure.atoms[i],atom) <= cce_vars.cutoffs[structure.atoms[i].type] ){ // distance must be larger than _CCE_SELF_DIST_TOL_ since GetStrNeighData includes also structure.atoms[i] itself as neighbor having distance zero to itself
          if (!anion_species.empty()){ // variable called anion type since function was developed for CCE for polar materials but it can be used to check for any atom type and only include those as neighbors
            // implement check whether each nearest neighbor is of the anion_species, otherwise throw warning; 
            if (atom.name == anion_species){
              neighbors_count+=1;
            } else if (atom.name != anion_species && structure.atoms[i].name != anion_species){ // second condition set since it is expected that the anion has predominantly other neighbors than its own type 
              if (!cce_flags.flag("MULTI_ANION_SYSTEM")){
                if (empty_line_count == 0){ // construction just to make sure that only one empty line is added at the beginning of the warning block
                  if(cce_flags.flag("COMMAND_LINE")){
                    cout << endl;
                  }
                  empty_line_count+=1;
                }
                warning = true;
                if(cce_flags.flag("COMMAND_LINE")){
                  cout << "WARNING: Not all nearest neighbors of " << structure.atoms[i].name << " (ATOM[" << i << "]) within the distance tolerance of " << tolerance << " Ang. are " << anion_species << ", there is also " << atom.name << endl;
                }
              }
            }
          } else {
            neighbors_count+=1;
          }
        }
      }
      if (!cce_flags.flag("MULTI_ANION_SYSTEM")){
        if (warning){
          if(cce_flags.flag("COMMAND_LINE")){
            cout << "WARNING: Not all nearest neighbors of " << structure.atoms[i].name << " (ATOM[" << i << "]) within the distance tolerance are " << anion_species << "!" << endl;
          }
        }
      }
      num_neighbors[i]=neighbors_count; // zero-based counting as for cutoffs array above
      if(LDEBUG){
        cerr << "number of " << anion_species << " nearest neighbors within " << tolerance << " Ang. tolerance of " << structure.atoms[i].name << " (ATOM[" << i << "]): " << num_neighbors[i] << endl;
        cerr << endl;
      }
    }
    return num_neighbors;
  }

  //CCE_get_dist_cutoffs////////////////////////////////////////////////////////
  // determine species selective nearest neighbor distances and then cutoffs accordingly
  vector<double> CCE_get_dist_cutoffs(double tolerance, const xstructure& structure) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    vector<double> cutoffs(structure.species.size());
    xmatrix<double> distsij=GetDistMatrix(structure); // gives matrix with nearest neighbor distances between all species pairs with species running over rows and columns
    vector<double> near_neigh_dist(structure.species.size());
    for(uint i=1,isize=distsij.rows+1;i<isize;i++){
      for(uint j=1,jsize=distsij.cols+1;j<jsize;j++){  
        if (j==1){
          near_neigh_dist[i-1]=distsij(i,j);
        }
        if (distsij(i,j) < near_neigh_dist[i-1]){
          near_neigh_dist[i-1]=distsij(i,j);
        }
      }
      if(LDEBUG){
        cerr << "nearest neighbor distance for species " << i << " is " << near_neigh_dist[i-1] << " Ang." << endl;
        cerr << "cutoff for the distance for species " << i << " is " << near_neigh_dist[i-1]+tolerance << " Ang." << endl; 
      }
      cutoffs[i-1]=near_neigh_dist[i-1]+tolerance; // -1 since counting of array elements starts from zero, NOT 1
    }
    return cutoffs;
  }

  //CCE_check_per_super_oxides////////////////////////////////////////////////////////
  // check whether the system contains per- or superoxide ions based on the O-O bond length and set variables accordingly
  void CCE_check_per_super_oxides(xstructure& structure, CCE_Variables& cce_vars, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::CCE_check_per_super_oxides():";
    stringstream message;
    if(LDEBUG){
      cerr << "CHECKING FOR (SU-)PEROXIDES:" << endl;
    }
    double perox_cutoff=1.6; // O-O bonds in peroxides for the studied examples are all shorter than 1.6 Ang
    double superox_cutoff=1.4; // O-O bonds in superoxides for the studied examples are all shorter than 1.4 Ang
    double O2_molecule_cutoff=1.2; // O-O bonds n the O2 molecule is about 1.21 Ang.
    uint perox_count=0;
    uint superox_count=0;
    for ( uint i = 0; i < structure.atoms.size(); i++ ) { // initialize elements of vectors to 0 
      cce_vars.perox_indices[i] = 0; 
      cce_vars.superox_indices[i] = 0; 
    }
    cce_vars.cutoffs=CCE_get_dist_cutoffs(_CCE_NN_DIST_TOL_, structure);
    double cutoffs_max=aurostd::max(cce_vars.cutoffs);
    deque<deque<_atom> > neigh_mat;
    structure.GetStrNeighData(cutoffs_max,neigh_mat);
    for(uint i=0,isize=neigh_mat.size();i<isize;i++){ //same size as structure.atoms.size(); number of atoms in the structure (not determined by cutoff (or cutoffs_max))
      if (structure.atoms[i].name == "O"){ // identify per- and superoxides by O-O bond length
        for(uint j=0,jsize=neigh_mat[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by the cutoffs_max
          const _atom& atom=neigh_mat[i][j];
          if (atom.name == "O"){
            if (_CCE_SELF_DIST_TOL_ < AtomDist(structure.atoms[i],atom) && AtomDist(structure.atoms[i],atom) <= O2_molecule_cutoff ){ // distance must be larger than _CCE_SELF_DIST_TOL_ to savely exclude the anion itself having distance zero to itself; if O-O bond is shorter than in O2 molecule (approx. 1.21 Ang) the result of the structural relaxation is most likely wrong 
              message << "THE DETERMINED OXYGEN-OXYGEN BOND LENGTH IS SHORTER THAN IN THE O2 MOLECULE; CHECK YOUR STRUCTURE! THE O-O BOND LENGTH IS: " << AtomDist(structure.atoms[i],atom) << " Ang.";
              throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
            } else if (O2_molecule_cutoff < AtomDist(structure.atoms[i],atom) && AtomDist(structure.atoms[i],atom) <= superox_cutoff ){
              if(LDEBUG){
                cerr << "WARNING: This should be a superoxide; the O-O bond length is: " << AtomDist(structure.atoms[i],atom) << " Ang." << endl;
              }
              superox_count+=1;
              cce_vars.superox_indices[i]=1;
            } else if (superox_cutoff < AtomDist(structure.atoms[i],atom) && AtomDist(structure.atoms[i],atom) <= perox_cutoff ){
              if(LDEBUG){
                cerr << "WARNING: This should be a peroxide; the O-O bond length is: " << AtomDist(structure.atoms[i],atom) << " Ang." << endl;
              }
              perox_count+=1;
              cce_vars.perox_indices[i]=1;
            }
          }
        }
      }
    }
    cce_vars.num_perox_bonds=perox_count/2; // needs to be divided by two due to double counting of O-O bonds when looping over all atoms of the structure
    cce_vars.num_superox_bonds=superox_count/2;
    // print out the number of per- & superoxide O-O bonds;
    if (cce_vars.num_perox_bonds > 0){
      if(cce_flags.flag("COMMAND_LINE")){
        cout << endl;
        cout << "WARNING: This should be a peroxide!" << endl;
        cout << "Number of peroxide O-O bonds in cell: " << cce_vars.num_perox_bonds << endl;
      }
    } else if (cce_vars.num_superox_bonds > 0) {
      if(cce_flags.flag("COMMAND_LINE")){
        cout << endl;
        cout << "WARNING: This should be a superoxide!" << endl;
        cout << "Number of superoxide O-O bonds in cell: " << cce_vars.num_superox_bonds << endl;
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //           DETERMINE OXIDATION STATES FROM ELECTRONEGATIVITIES           //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // determine oxidation numbers based on Allen electronegativities and treat mixed valence binary systems specially

  //CCE_get_oxidation_states_from_electronegativities////////////////////////////////////////////////////////
  // function overloading for below function to be able to use oxidation number determination independently of CCE
  vector<double> CCE_get_oxidation_states_from_electronegativities(xstructure& structure) {
    CCE_Variables cce_vars = CCE_init_variables(structure);
    cce_vars.anion_species=CCE_determine_anion_species(structure, cce_vars);
    aurostd::xoption cce_flags = CCE_init_flags();
    cce_vars.multi_anion_atoms=CCE_check_for_multi_anion_system(cce_vars, _CCE_NN_DIST_TOL_MULTI_ANION_, structure, cce_flags);
    // multi anion corrections can only be resized after number of multi anion species is known from check for multi anion system
    cce_vars.multi_anion_corrections_atom.clear();
    cce_vars.multi_anion_corrections_atom.resize(cce_vars.multi_anion_species.size(), vector<vector<double> >(cce_vars.vfunctionals.size()*2, vector<double>(structure.atoms.size())));
    cce_vars.num_neighbors=CCE_get_num_neighbors(cce_vars.anion_species, _CCE_NN_DIST_TOL_, structure, cce_flags, cce_vars);
    if(cce_vars.anion_species == "O" || cce_flags.flag("O_MULTI_ANION_SPECIES")) {
      CCE_check_per_super_oxides(structure, cce_vars, cce_flags);
    }
    cce_flags.flag("CORRECTABLE",TRUE);
    return CCE_get_oxidation_states_from_electronegativities(structure, cce_vars, cce_flags);
  }

  //CCE_get_oxidation_states_from_electronegativities////////////////////////////////////////////////////////
  // determine the oxidation numbers of the ions using preferred/all known oxidation numbers, electronegativities and structural information
  vector<double> CCE_get_oxidation_states_from_electronegativities(const xstructure& structure, CCE_Variables& cce_vars, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << "DETERMINATION OF OXIDATION NUMBERS FROM PREFERRED/ALL KNOWN OXIDATION STATES, ALLEN ELECTRONEGATIVITIES, AND STRUCTURE:" << endl;
    }
    // deal with anion charge in cell first
    CCE_set_anion_oxidation_states(structure, cce_vars);
    // treat cations
    if(LDEBUG){
      cerr << "CATION PART:" << endl;
    }
    // determine number of cation species (currently just all species minus one anion species, multi-anion atoms are dealt with separately)
    uint num_cation_species = structure.species.size()-1;
    // load needed info about preferred and all known oxidation states for all species
    cce_flags.flag("NO_PREF_OX_STATES",FALSE);
    cce_flags.flag("NO_OX_STATES",FALSE);
    cce_flags.flag("OX_STATES_DETERMINED",FALSE);
    CCE_load_ox_states_templates_each_species(structure, cce_vars, cce_flags);
    // sort species, electronegativities, number of preferred/all oxidation states, and preferred/all oxidation states strings ascending by electronegativity
    // the anion species should always be the last one with the highest electronegativity
    CCE_sort_species_pref_all_ox_states_by_electronegativity(structure, cce_vars);
    // ME Nov. 2019 for getting cations_map: a vector of vectors that lists for each cation species the atom numbers of the structure that are of this species (for Fe2ZnO4 there might be two Fe atoms at positions 0 and 1 in the structure)
    uint natoms = structure.atoms.size();
    cce_vars.cations_map.resize(num_cation_species);
    uint i = 0;
    for (uint at = 0; at < natoms; at++) {
      for (i = 0; i < num_cation_species; i++) {
        if (structure.atoms[at].name == cce_vars.species_electronegativity_sorted[i]) break; // if it finds that the atom belongs to the ith species sorted by electronegativities, break to insert it at the proper place for the ith species into cation_map
      }
      if (i < num_cation_species) cce_vars.cations_map[i].push_back(at);
    }
    // try to find proper oxidation states (making oxidation sum zero) by using preferred oxidation states
    if(!cce_flags.flag("NO_PREF_OX_STATES")){ // for He, Ne, Ar, and Xe there are no preferred ox. states
      CCE_try_preferred_oxidation_states(structure, cce_flags, cce_vars);
    }
    // if preferred oxidation states do not work, it could be a mixed valence (binary) system
    // treat them here as special cases
    if(!cce_flags.flag("OX_STATES_DETERMINED")){
      // for SbO2 the oxidation states are not identified properly, Bader analysis finds them automatically
      CCE_treat_SbO2_special_case(cce_vars, structure, cce_flags);
      // for Pb3O4 the oxidation states are not identified properly, Bader analysis finds them automatically
      CCE_treat_Pb3O4_special_case(cce_vars, structure, cce_flags);
      // Ti-O Magneli phases need to be treated specially since oxidation numbers are not recognized appropriately by above approach
      CCE_treat_Ti_O_Magneli_phase_special_case(cce_vars, structure, cce_flags);
      // for Fe3O4 in inverse spinel structure, the oxidation states are not identified properly
      CCE_treat_Fe3O4_special_case(cce_vars, structure, cce_flags);
      // for Mn3O4 in spinel structure, the oxidation states are not identified properly
      CCE_treat_Mn3O4_special_case(cce_vars, structure, cce_flags);
      // for Co3O4 in spinel structure, the oxidation states are not identified properly
      CCE_treat_Co3O4_special_case(cce_vars, structure, cce_flags);
      // alkali metal sesquioxides need to be treated specially since oxidation numbers and number of per- and superoxide bonds are not recognized appropriately
      CCE_treat_alkali_sesquioxide_special_case(cce_vars, structure, cce_flags);
      if(cce_flags.flag("OX_STATES_DETERMINED") && cce_flags.flag("COMMAND_LINE")){
        cout << endl;
      }
    }
    // if preferred oxidation states approach and mixed valence doesn't work, try to find proper oxidation states (making oxidation sum zero) by using all known oxidation states
    if(!cce_flags.flag("OX_STATES_DETERMINED")){
      if(!cce_flags.flag("NO_OX_STATES")){ // for He, Ne, and Ar there are no known oxidation states
        CCE_try_all_oxidation_states(structure, cce_vars);
        // print oxidation numbers and calculate sum
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
        }
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
        }
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << endl;
          cerr << "BAD NEWS: The determined oxidation numbers do not add up to zero!"  << endl;
          if(cce_flags.flag("COMMAND_LINE")){
            cerr << "The formation energy of this system is hence not correctable!"  << endl;
            cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
          }
          cerr << endl;
        }
      } else{ // error message that oxidation numbers cannot be determined since for at least one species there are no known oxidation numbers is already included in CCE_load_ox_states_templates_each_species function
        cce_flags.flag("CORRECTABLE",FALSE);
        if(cce_flags.flag("COMMAND_LINE")){
          cerr << "The formation energy of this system is hence not correctable!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        }
        cerr << endl;
      }
    }
    return cce_vars.oxidation_states;
  }

  //CCE_set_anion_oxidation_states////////////////////////////////////////////////////////
  // determine the oxidation numbers of the anions
  void CCE_set_anion_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << "ANION PART:" << endl;
    }
    double total_anion_charge=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].name == cce_vars.anion_species || cce_vars.multi_anion_atoms[i] == 1){
        if (structure.atoms[i].name == cce_vars.anion_species){ // for multi-anion atoms oxidation states have been assigned previously
          cce_vars.oxidation_states[i] = cce_vars.standard_anion_charge; // oxidation numbers for O are assumed to be -2 (other anions accordingly) and are corrected below if it is a per- or superoxide O atom as identified from the structural analysis in other functions
        }
        if (cce_vars.num_perox_bonds > 0){
          if (cce_vars.perox_indices[i]==1){
            cce_vars.oxidation_states[i]=-1;
          }
        }
        if (cce_vars.num_superox_bonds > 0){
          if (cce_vars.superox_indices[i]==1){
            cce_vars.oxidation_states[i]=-0.5;
          }
        }
        if(LDEBUG){
          if (cce_vars.multi_anion_atoms[i] == 1){ // for multi anion atoms oxidation states have been assigned previously
            cerr << "anion oxidation number for multi-anion ATOM[" << i << "] (" << structure.atoms[i].name << ") has been assigned previously to: " << cce_vars.oxidation_states[i] << endl;
          } else {
            cerr << "anion oxidation number for ATOM[" << i << "] (" << structure.atoms[i].name << "): " << cce_vars.oxidation_states[i] << endl;
          }
        }
        total_anion_charge += cce_vars.oxidation_states[i];
      }
    }
    if(LDEBUG){
      cerr << "Total anion charge in cell: " << total_anion_charge << endl;
    }
  }

  //CCE_load_ox_states_templates_each_species////////////////////////////////////////////////////////
  // load templates for preferred and other oxidation states
  void CCE_load_ox_states_templates_each_species(const xstructure& structure, CCE_Variables& cce_vars, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    for(uint i=0,isize=structure.species.size();i<isize;i++){ 
      string electroneg_line = CCE_get_electronegativities_ox_nums(structure.species[i]);
      if(LDEBUG){
        cerr << "electronegativity line for species " << i << "(" << structure.species[i] << "): " << electroneg_line << endl;
      }
      vector<string> electroneg_tokens;
      aurostd::string2tokens(electroneg_line, electroneg_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one
      // load preferred oxidation states for each species
      cce_vars.num_pref_ox_states[i] = aurostd::string2utype<uint>(electroneg_tokens[1]);
      if(cce_vars.num_pref_ox_states[i] > 0){
        cce_vars.pref_ox_states_strings[i] = electroneg_tokens[3];
        if(LDEBUG){
          cerr << "preferred oxidation states string of species " << i << " (" << structure.species[i] << "): " << cce_vars.pref_ox_states_strings[i] << endl;
        }
      } else{
        cce_flags.flag("NO_PREF_OX_STATES",TRUE);
        if(LDEBUG){
          cerr << endl;
          cerr << "BAD NEWS: There are no preferred oxidation states for species " << i << " (" << structure.species[i] <<")."  << endl;
          cerr << "Therefore the oxidation states cannot be determined on this basis." << endl;
        }
      }
      // load all oxidation states for each species
      cce_vars.num_all_ox_states[i] = aurostd::string2utype<uint>(electroneg_tokens[2]);
      if(cce_vars.num_all_ox_states[i] > 0){
        cce_vars.all_ox_states_strings[i] = electroneg_tokens[4];
        if(LDEBUG){
          cerr << "all oxidation states string of species " << i << " (" << structure.species[i] << "): " << cce_vars.all_ox_states_strings[i] << endl;
          cerr << endl;
        }
      } else{
        cce_flags.flag("NO_OX_STATES",TRUE);
        cerr << endl;
        cerr << "BAD NEWS: There are no known oxidation states for species " << i << " (" << structure.species[i] <<")."  << endl;
        cerr << "Therefore the oxidation states cannot be determined on this basis." << endl;
        cerr << endl;
      }
    }
  }

  //CCE_sort_species_pref_all_ox_states_by_electronegativity////////////////////////////////////////////////////////
  // sort species, electronegativities, number of preferred/all oxidation states, and preferred/all oxidation states strings ascending by electronegativity
  void CCE_sort_species_pref_all_ox_states_by_electronegativity(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    // using aurostd sort functions
    vector<double> electronegativities_sorted = cce_vars.electronegativities;
    for(uint j=0,jsize=structure.species.size();j<jsize;j++){ //loop over all species
     cce_vars.species_electronegativity_sorted[j] = structure.species[j];
    }
    cce_vars.num_pref_ox_states_electronegativity_sorted = cce_vars.num_pref_ox_states;
    cce_vars.pref_ox_states_strings_electronegativity_sorted = cce_vars.pref_ox_states_strings;
    cce_vars.num_all_ox_states_electronegativity_sorted = cce_vars.num_all_ox_states;
    cce_vars.all_ox_states_strings_electronegativity_sorted = cce_vars.all_ox_states_strings;
    // for the oxidation state algorithm the species, preferred and all known oxidation states of all species must be sorted by electronegativity
    aurostd::sort(electronegativities_sorted,cce_vars.species_electronegativity_sorted);
    electronegativities_sorted = cce_vars.electronegativities;
    aurostd::sort(electronegativities_sorted,cce_vars.num_pref_ox_states_electronegativity_sorted);
    electronegativities_sorted = cce_vars.electronegativities;
    aurostd::sort(electronegativities_sorted,cce_vars.pref_ox_states_strings_electronegativity_sorted);
    electronegativities_sorted = cce_vars.electronegativities;
    aurostd::sort(electronegativities_sorted,cce_vars.num_all_ox_states_electronegativity_sorted);
    electronegativities_sorted = cce_vars.electronegativities;
    aurostd::sort(electronegativities_sorted,cce_vars.all_ox_states_strings_electronegativity_sorted);
    for(uint j=0,jsize=structure.species.size();j<jsize;j++){ //loop over all species
      if(LDEBUG){
        cerr << "species_electronegativity_sorted[" << j << "]: " << cce_vars.species_electronegativity_sorted[j] << endl;
        cerr << "electronegativities_sorted[" << j << "]: " << electronegativities_sorted[j] << endl;
        cerr << "num_pref_ox_states_electronegativity_sorted[" << j << "]: " << cce_vars.num_pref_ox_states_electronegativity_sorted[j] << endl;
        cerr << "pref_ox_states_strings_electronegativity_sorted[" << j << "]: " << cce_vars.pref_ox_states_strings_electronegativity_sorted[j] << endl;
        cerr << "num_all_ox_states_electronegativity_sorted[" << j << "]: " << cce_vars.num_all_ox_states_electronegativity_sorted[j] << endl;
        cerr << "all_ox_states_strings_electronegativity_sorted[" << j << "]: " << cce_vars.all_ox_states_strings_electronegativity_sorted[j] << endl;
        cerr << endl;
      }
    }
  }

  //CCE_try_preferred_oxidation_states////////////////////////////////////////////////////////
  // try to determine the oxidation numbers using the preferred oxidation states for all cation species
  void CCE_try_preferred_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << "Trying preferred oxidation numbers:" << endl;
    }
    // determine possible_ox_states vector of vectors of doubles needed for Marco's function of my algorithm to determine oxidation numbers
    vector<vector<double> > pref_ox_states_electronegativity_sorted(structure.species.size());
    for(uint i=0,isize=structure.species.size();i<isize;i++){ 
      vector<string> pref_ox_states_tokens(cce_vars.num_pref_ox_states_electronegativity_sorted[i]);
      aurostd::string2tokens(cce_vars.pref_ox_states_strings_electronegativity_sorted[i], pref_ox_states_tokens, ",");
      for(uint j=0,jsize=cce_vars.num_pref_ox_states_electronegativity_sorted[i];j<jsize;j++){ 
        pref_ox_states_electronegativity_sorted[i].push_back(aurostd::string2utype<double>(pref_ox_states_tokens[j]));
      }
    }
    // use Marco's implementation of my algorithm to determine oxidation_numbers
    CCE_determine_cation_oxidation_states(pref_ox_states_electronegativity_sorted, structure, cce_vars);
    // print oxidation numbers and calculate sum
    cce_vars.oxidation_sum = CCE_get_oxidation_states_sum(cce_vars);
    // if sum of oxidation numbers is essentially zero, oxidation states should be regarded as determined correctly
    if (std::abs(cce_vars.oxidation_sum) <= _CCE_OX_TOL_) {
      cce_flags.flag("OX_STATES_DETERMINED",TRUE);
      if(cce_flags.flag("COMMAND_LINE")){
        cout << endl;
        CCE_print_oxidation_states_and_sum(cce_vars, structure);
        cout << endl;
      }
    }
  }

  //CCE_treat_SbO2_special_case////////////////////////////////////////////////////////
  // for SbO2 the oxidation states are not identified properly
  // with the actual formula Sb2O4 it is a mixed valence oxide with one Sb+3 with 4 Sb-O bonds and one Sb+5 with 6 Sb-O bonds
  void CCE_treat_SbO2_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::CCE_treat_SbO2_special_case():";
    stringstream message;
    if ( structure.species.size() == 2 && ((structure.species[0] == "O" && structure.species[1] == "Sb") || (structure.species[0] == "Sb" && structure.species[1] == "O")) ) {
      uint num_O_before_Sb;
      if ( structure.species[0] == "O" && structure.species[1] == "Sb" ) {
        num_O_before_Sb=4;
      } else {
        num_O_before_Sb=0;
      }
      double amount_Sb=0;
      double amount_O=0;
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].name == "Sb"){
          amount_Sb+=1;
        } else if (structure.atoms[i].name == "O"){
          amount_O+=1;
        }
      }
      if(LDEBUG){
        cerr << "number of Sb ions= " << amount_Sb << endl;
        cerr << "number of O ions= " << amount_O << endl;
      }
      double Sb_O_ratio;
      if ( amount_O != 0 ){
        Sb_O_ratio=amount_Sb/amount_O;
      } else {
        message << "Amount O determined to be ZERO. Please check your structure.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
      }
      if ( Sb_O_ratio == 0.5 ){
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
          cout << "WARNING: This system is identified as a mixed valence compound." << endl; 
          cout << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
          cout << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
          cout << "Sb2O4 with ratio of Sb/O= " << Sb_O_ratio << endl;
        }
        double num_formula_units_in_cell;
        num_formula_units_in_cell=structure.atoms.size()/6; // 6 for Sb2O4
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "Sb"){
            // taking the first Sb ions as +3 (without knowing whether they are actually +3); works only because 
            // the number of bonds for them will be adjusted to 4 as needed for Sb2O4 disregarding the actual 
            // number of Sb-O bonds (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (1+num_O_before_Sb)*num_formula_units_in_cell ){  // (1 Sb3+ ions per formula unit + 4*O listed before in alphabetic order) * number of formula units
              cce_vars.oxidation_states[i]=+3;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Sb+3 " << endl;
              }
              cce_vars.num_neighbors[i]=4; // for Sb2O4 the Sb3+ ions are 4-fold coordinated by oxygen
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "Modified number of neighbors of Sb (atom[" << i << "]) taken as Sb+3: " << cce_vars.num_neighbors[i] << endl;
              }
            } else {
              cce_vars.oxidation_states[i]=5;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Sb+5 " << endl;
              }
              cce_vars.num_neighbors[i]=6; // for Sb2O4 the Sb5+ ions are 6-fold coordinated by oxygen
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "Modified number of neighbors of Sb (atom[" << i << "]) taken as Sb+5: " << cce_vars.num_neighbors[i] << endl;
              }
            }
          }
        }
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) { // this case should never occur since for a Sb-O binary with Sb/O ratio 0.5 the above set oxidation states should always work
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        } else {
          cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
        }
      }
    }
  }

  //CCE_treat_Pb3O4_special_case////////////////////////////////////////////////////////
  // for Pb3O4 the oxidation states are not identified properly
  // https://en.wikipedia.org/wiki/Lead(II,IV)_oxide
  void CCE_treat_Pb3O4_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if ( structure.species.size() == 2 && ((structure.species[0] == "O" && structure.species[1] == "Pb") || (structure.species[0] == "Pb" && structure.species[1] == "O")) ) {
      uint num_O_before_Pb;
      if ( structure.species[0] == "O" && structure.species[1] == "Pb" ) {
        num_O_before_Pb=4;
      } else {
        num_O_before_Pb=0;
      }
      double amount_Pb=0;
      double amount_O=0;
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].name == "Pb"){
          amount_Pb+=1;
        } else if (structure.atoms[i].name == "O"){
          amount_O+=1;
        }
      }
      if(LDEBUG){
        cerr << "number of Pb ions= " << amount_Pb << endl;
        cerr << "number of O ions= " << amount_O << endl;
      }
      double Pb_O_ratio=amount_Pb/amount_O;
      if ( Pb_O_ratio == 0.75 ){
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
          cout << "WARNING: This system is identified as a mixed valence compound." << endl; 
          cout << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
          cout << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
          cout << "Pb3O4 with ratio of Pb/O= " << Pb_O_ratio << endl;
        }
        double num_formula_units_in_cell;
        num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Pb3O4
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "Pb"){
            // taking the first Pb ions as +4 (without knowing whether they are actually +4); works only because 
            // the number of bonds for them will be adjusted to 6 as needed for Pb3O4 disregarding the actual 
            // number of Pb-O bonds (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (1+num_O_before_Pb)*num_formula_units_in_cell ){  // (1 Pb4+ ions per formula unit + 4*O listed before in alphabetic order) * number of formula units
              cce_vars.oxidation_states[i]=+4;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Pb+4 " << endl;
              }
              cce_vars.num_neighbors[i]=6; // for Pb3O4 the Pb4+ ions are 6-fold coordinated by oxygen
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "Modified number of neighbors of Pb (atom[" << i << "]) taken as Pb+4: " << cce_vars.num_neighbors[i] << endl;
              }
            } else {
              cce_vars.oxidation_states[i]=+2;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Pb+2 " << endl;
              }
              cce_vars.num_neighbors[i]=3; // for Pb3O4 the Pb2+ ions are 3-fold coordinated by oxygen
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "Modified number of neighbors of Pb (atom[" << i << "]) taken as Pb+2: " << cce_vars.num_neighbors[i] << endl;
              }
            }
          }
        }
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) { // this case should never occur since for a Pb-O binary with Pb/O ratio 0.75 the above set oxidation states should always work
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        } else {
          cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
        }
      }
    }
  }

  //CCE_treat_Ti_O_Magneli_phase_special_case////////////////////////////////////////////////////////
  // treat Ti-O Magneli phases; there are always 2xTi+3 per formula unit and the rest is Ti+4; 
  // fortunately, both Ti+3 and Ti+4 have 6 Ti-O bonds, hence one only needs to know how many ions 
  // of the respective oxidation state there are, not which one is which
  void CCE_treat_Ti_O_Magneli_phase_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if ( structure.species.size() == 2 && ((structure.species[0] == "O" && structure.species[1] == "Ti") || (structure.species[0] == "Ti" && structure.species[1] == "O")) ) {
      uint num_O_before_Ti;
      if ( structure.species[0] == "O" && structure.species[1] == "Ti" ) {
        num_O_before_Ti=5;
      } else {
        num_O_before_Ti=0;
      }
      if(LDEBUG){
        cerr << "Ti-O system, Magneli for Ti_nO_(2n-1), i.e. Ti-O ratio= " << 3.0/5 << ", " << 4.0/7 << ", " << 5.0/9 << ", " << 6.0/11 << ", " << 7.0/13 << ", " << 8.0/15 << ", " << 9.0/17 << ", " << 10.0/19 << "..." << endl;
      }
      double amount_O=0;
      double amount_Ti=0;
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].name == "O"){
          amount_O+=1;
        } else if (structure.atoms[i].name == "Ti"){
          amount_Ti+=1;
        }
      }
      if(LDEBUG){
        cerr << "number of O ions= " << amount_O << endl;
        cerr << "number of Ti ions= " << amount_Ti << endl;
      }
      double Ti_O_ratio=amount_Ti/amount_O;
      if(LDEBUG){
        cerr << "ratio of Ti/O= " << Ti_O_ratio << endl;
      }
      double num_formula_units_in_cell;
      // check for Magneli composition Ti_(n)O_(2n-1)
      double n;
      bool magneli = false;
      for(n=3;n<101;n++){
        //cout << "n/(2*n-1)= " << n/(2*n-1) << endl;
        if ( Ti_O_ratio == n/(2*n-1) ){
          if(cce_flags.flag("COMMAND_LINE")){
            cout << endl;
            cout << "WARNING: This system is identified as a mixed valence compound." << endl; 
            cout << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
            cout << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
            cout << "n= " << n << " Magneli composition Ti_nO_(2n-1)" << endl;
          }
          magneli = true;
          num_formula_units_in_cell=amount_Ti/n;
          if(cce_flags.flag("COMMAND_LINE")){
            cout << "number of formula units in cell: " << num_formula_units_in_cell << endl;
          }
        }
      }
      if ( magneli == false){
        if(LDEBUG){
          cerr << "Not a Magneli composition." << endl;
        }
      }
      if (magneli){
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "Ti"){
            // taking the first Ti as +3 (without knowing whether they are actually +3) works only because 
            // for the Magneli phases both Ti+3 and Ti+4 have both always 6 Ti-O bonds
            if ( i < (2+num_O_before_Ti)*num_formula_units_in_cell ){
              cce_vars.oxidation_states[i]=+3;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Ti+3 " << endl;
              }
            } else {
              cce_vars.oxidation_states[i]=+4;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Ti+4 " << endl;
              }
            }
          }
        }
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) { // this case should never occur since for a Magneli phase the above set oxidation states should always work
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        } else {
          cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
        }
      }
    }
  }

  //CCE_treat_Fe3O4_special_case////////////////////////////////////////////////////////
  // for Fe3O4 in inverse spinel structure the oxidation states are not identified properly
  // According to Wikipedia the Fe2+ ions are octahedrally coordinated and the Fe3+ ions are evenly distributed 
  // between the octahedral and tetrahedral sites and one hence needs per formula unit 6x the Fe2+ correction 
  // and 4+6=10x the Fe3+ correction, see:
  // https://en.wikipedia.org/wiki/Iron(II,III)_oxide
  // for Co3O4 and Mn3O4 this is different; in the normal spinel structure Co2+/Mn2+ occupies only tetrahedral sites 
  // while Co3+/Mn3+ occupies octahedral sites; the correction hence needs to be different but at the present stage 
  // there are no corrections for Co3+ and Mn3+, see:
  // https://en.wikipedia.org/wiki/Cobalt(II,III)_oxide
  // https://en.wikipedia.org/wiki/Manganese(II,III)_oxide
  void CCE_treat_Fe3O4_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if ( structure.species.size() == 2 && ((structure.species[0] == "O" && structure.species[1] == "Fe") || (structure.species[0] == "Fe" && structure.species[1] == "O")) ) {
      uint num_O_before_Fe;
      if ( structure.species[0] == "O" && structure.species[1] == "Fe" ) {
        num_O_before_Fe=4;
      } else {
        num_O_before_Fe=0;
      }
      double amount_Fe=0;
      double amount_O=0;
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].name == "Fe"){
          amount_Fe+=1;
        } else if (structure.atoms[i].name == "O"){
          amount_O+=1;
        }
      }
      if(LDEBUG){
        cerr << "number of Fe ions= " << amount_Fe << endl;
        cerr << "number of O ions= " << amount_O << endl;
      }
      double Fe_O_ratio=amount_Fe/amount_O;
      if ( Fe_O_ratio == 0.75 ){
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
          cout << "WARNING: This system is identified as a mixed valence compound." << endl; 
          cout << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
          cout << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
          cout << "Fe3O4 with ratio of Fe/O= " << Fe_O_ratio << endl;
        }
        double num_formula_units_in_cell;
        num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Fe3O4
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "Fe"){
            // taking the first Fe ions as +2 (without knowing whether they are actually +2); works only because 
            // the number of bonds for them will be adjusted to 6 (octahedral) as needed for Fe3O4 disregarding the actual 
            // number of Fe-O bonds (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (1+num_O_before_Fe)*num_formula_units_in_cell ){  // 1 Fe2+ ions per formula unit * number of formula units
              cce_vars.oxidation_states[i]=+2;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Fe+2 " << endl;
              }
              cce_vars.num_neighbors[i]=6; // for Fe3O4 the Fe2+ ions are 6-fold coordinated by oxygen according to Wikipedia
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "Modified number of neighbors of Fe (atom[" << i << "]) taken as Fe+2: " << cce_vars.num_neighbors[i] << endl;
              }
            } else {
              cce_vars.oxidation_states[i]=+3;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Fe+3 " << endl;
              }
              cce_vars.num_neighbors[i]=5; // for Fe3O4 the Fe3+ ions are evenly 6- and 4-fold, so on average 5-fold (set here as a hack) coordinated by oxygen according to Wikipedia
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "Modified number of neighbors of Fe (atom[" << i << "]) taken as Fe+3 (average between even 6- and 4-fold coordination): " << cce_vars.num_neighbors[i] << endl;
              }
            }
          }
        }
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) { // this case should never occur since for a Fe-O binary with Fe/O ratio 0.75 the above set oxidation states should always work
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        } else {
          cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
        }
      }
    }
  }

  //CCE_treat_Mn3O4_special_case////////////////////////////////////////////////////////
  // for Co3O4 and Mn3O4 in the normal spinel structure Co2+/Mn2+ occupies only tetrahedral sites 
  // while Co3+/Mn3+ occupies octahedral sites; the correction hence needs to be different than for Fe3O4
  // but at the present stage there are no corrections for Co3+ and Mn3+, see:
  // https://en.wikipedia.org/wiki/Cobalt(II,III)_oxide
  // https://en.wikipedia.org/wiki/Manganese(II,III)_oxide
  void CCE_treat_Mn3O4_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if ( structure.species.size() == 2 && ((structure.species[0] == "O" && structure.species[1] == "Mn") || (structure.species[0] == "Mn" && structure.species[1] == "O")) ) {
      uint num_O_before_Mn;
      if ( structure.species[0] == "O" && structure.species[1] == "Mn" ) {
        num_O_before_Mn=4;
      } else {
        num_O_before_Mn=0;
      }
      double amount_Mn=0;
      double amount_O=0;
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].name == "Mn"){
          amount_Mn+=1;
        } else if (structure.atoms[i].name == "O"){
          amount_O+=1;
        }
      }
      if(LDEBUG){
        cerr << "number of Mn ions= " << amount_Mn << endl;
        cerr << "number of O ions= " << amount_O << endl;
      }
      double Mn_O_ratio=amount_Mn/amount_O;
      if ( Mn_O_ratio == 0.75 ){
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
          cout << "WARNING: This system is identified as a mixed valence compound." << endl; 
          cout << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
          cout << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
          cout << "Mn3O4 with ratio of Mn/O= " << Mn_O_ratio << endl;
        }
        double num_formula_units_in_cell;
        num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Mn3O4
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "Mn"){
            // taking the first Mn ions as +2 (without knowing whether they are actually +2); works only because 
            // the number of bonds for them will be adjusted to 4 (tetrahedral) as needed for Mn3O4 disregarding the actual 
            // number of Mn-O bonds (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (1+num_O_before_Mn)*num_formula_units_in_cell ){  // 1 Mn2+ ions per formula unit * number of formula units
              cce_vars.oxidation_states[i]=+2;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Mn+2 " << endl;
              }
              cce_vars.num_neighbors[i]=4; // for Mn3O4 the Mn2+ ions are 4-fold coordinated by oxygen according to Wikipedia
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "Modified number of neighbors of Mn (atom[" << i << "]) taken as Mn+2: " << cce_vars.num_neighbors[i] << endl;
              }
            } else {
              cce_vars.oxidation_states[i]=+3;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Mn+3 " << endl;
              }
              cce_vars.num_neighbors[i]=6; // for Mn3O4 the Mn3+ ions are 6-fold coordinated by oxygen according to Wikipedia
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "Modified number of neighbors of Mn (atom[" << i << "]) taken as Mn+3: " << cce_vars.num_neighbors[i] << endl;
              }
            }
          }
        }
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) { // this case should never occur since for a Mn-O binary with Mn/O ratio 0.75 the above set oxidation states should always work
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        } else {
          cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
        }
      }
    }
  }

  //CCE_treat_Co3O4_special_case////////////////////////////////////////////////////////
  // for Co3O4 and Mn3O4 in the normal spinel structure Co2+/Mn2+ occupies only tetrahedral sites 
  // while Co3+/Mn3+ occupies octahedral sites; the correction hence needs to be different than for Fe3O4
  // but at the present stage there are no corrections for Co3+ and Mn3+, see:
  // https://en.wikipedia.org/wiki/Cobalt(II,III)_oxide
  // https://en.wikipedia.org/wiki/Manganese(II,III)_oxide
  void CCE_treat_Co3O4_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if ( structure.species.size() == 2 && ((structure.species[0] == "O" && structure.species[1] == "Co") || (structure.species[0] == "Co" && structure.species[1] == "O")) ) {
      uint num_O_before_Co;
      if ( structure.species[0] == "O" && structure.species[1] == "Co" ) {
        num_O_before_Co=4;
      } else {
        num_O_before_Co=0;
      }
      double amount_Co=0;
      double amount_O=0;
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].name == "Co"){
          amount_Co+=1;
        } else if (structure.atoms[i].name == "O"){
          amount_O+=1;
        }
      }
      if(LDEBUG){
        cerr << "number of Co ions= " << amount_Co << endl;
        cerr << "number of O ions= " << amount_O << endl;
      }
      double Co_O_ratio=amount_Co/amount_O;
      if ( Co_O_ratio == 0.75 ){
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
          cout << "WARNING: This system is identified as a mixed valence compound." << endl; 
          cout << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
          cout << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
          cout << "Co3O4 with ratio of Co/O= " << Co_O_ratio << endl;
        }
        double num_formula_units_in_cell;
        num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Co3O4
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "Co"){
            // taking the first Co ions as +2 (without knowing whether they are actually +2); works only because 
            // the number of bonds for them will be adjusted to 4 (tetrahedral) as needed for Co3O4 disregarding the actual 
            // number of Co-O bonds (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (1+num_O_before_Co)*num_formula_units_in_cell ){  // 1 Co2+ ions per formula unit * number of formula units
              cce_vars.oxidation_states[i]=+2;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Co+2 " << endl;
              }
              cce_vars.num_neighbors[i]=4; // for Co3O4 the Co2+ ions are 4-fold coordinated by oxygen according to Wikipedia
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "Modified number of neighbors of Co (atom[" << i << "]) taken as Co+2: " << cce_vars.num_neighbors[i] << endl;
              }
            } else {
              cce_vars.oxidation_states[i]=+3;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Co+3 " << endl;
              }
              cce_vars.num_neighbors[i]=6; // for Co3O4 the Co3+ ions are 6-fold coordinated by oxygen according to Wikipedia
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "Modified number of neighbors of Co (atom[" << i << "]) taken as Co+3: " << cce_vars.num_neighbors[i] << endl;
              }
            }
          }
        }
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) { // this case should never occur since for a Co-O binary with Co/O ratio 0.75 the above set oxidation states should always work
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        } else {
          cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
        }
      }
    }
  }

  //CCE_treat_alkali_sesquioxide_special_case////////////////////////////////////////////////////////
  // for alkali metal sesquioxides (X2O3) the oxidation states and corrections might not be found correctly
  // since these systems are formulated to contain both per- and superoxide ions
  // see: https://en.wikipedia.org/wiki/Sesquioxide
  // the only known example to date for which also a structure was found on Springer Materials and the AFLOW ICSD
  // is Rb2O3 formulated as [(Rb+)4(O2−2)(O2-1)2] (for Cs4O6 there is only a single entry in a scanned pdf on Springer Materials)
  // this function will treat it as an exceptional case for all possible alkali metal sesquioxides
  // the oxidation numbers and per- as well as superoxide corrections will be adjusted accordingly
  void CCE_treat_alkali_sesquioxide_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string alkali_metals = "Li,Na,K,Rb,Cs,Fr";
    vector<string> valkali_metals;
    aurostd::string2tokens(alkali_metals, valkali_metals, ",");
    if ( structure.species.size() == 2 && ((structure.species[0] == "O" && aurostd::withinList(valkali_metals, structure.species[1])) || (aurostd::withinList(valkali_metals, structure.species[0]) && structure.species[1] == "O")) ) { // check whether it is a binary alkali metal oxide
      if(LDEBUG){
        cerr << "This is a binary alkali metal oxide,checking whether it is an alkali metal sesquioxide..." << endl;
      }
      uint num_alkali_before_O; // num cations before O not O before cations since setting oxidation states of anions below, not for cations as in other cases
      string alkali_metal;
      if ( structure.species[0] == "O" && aurostd::withinList(valkali_metals, structure.species[1]) ) {
        num_alkali_before_O=0;
        alkali_metal=structure.species[1];
      } else {
        num_alkali_before_O=2;
        alkali_metal=structure.species[0];
      }
      double amount_O=0;
      double amount_alkali=0;
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].name == "O"){
          amount_O+=1;
        } else {
          amount_alkali+=1;
        }
      }
      if(LDEBUG){
        cerr << "number of O ions= " << amount_O << endl;
        cerr << "number of alkali (" << alkali_metal << ") ions= " << amount_alkali << endl;
      }
      double O_alkali_ratio=amount_O/amount_alkali;
      if(LDEBUG){
        cerr << "ratio of O/" << alkali_metal << "= " << O_alkali_ratio << endl;
      }
      // check for sesqui-composition alkali_metal2O3
      if ( O_alkali_ratio == 1.5 ){
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
          cout << "WARNING: This system is identified as an alkali metal sesquioxide (formally) containing both per- and superoxide ions." << endl; 
          cout << "The oxidation numbers and the number of per- and superoxide bonds will be set as known for this special case." << endl;
          cout << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each oxidation state occurs should be correct." << endl;
          cout << alkali_metal << "2O3 with ratio of O/" << alkali_metal << "= " << O_alkali_ratio << endl;
        }
        double num_formula_units_in_cell;
        num_formula_units_in_cell=structure.atoms.size()/5; // 5 for alkali_metal2O3
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "O"){
            // taking the first O ions as -0.5 (superox. Os, without knowing whether they are actually -0.5); works only because 
            // the number of superoxide bonds will be adjusted to 1 per formula unit as needed for alkali metal sesquioxides
            // taking last O ions as -1 (perox. Os, without knowing whether they are actually -1); works only because 
            // the number of peroxide bonds will be adjusted to 0.5 per formula unit as needed for alkali metal sesquioxides
            // (this is a hack since I know how many bonds there should be for each ion type)
            if ( i < (2+num_alkali_before_O)*num_formula_units_in_cell ){  // 2 superoxide O-atoms per formula unit * number of formula units
              cce_vars.oxidation_states[i]=-0.5;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to O-0.5 " << endl;
              }
              cce_vars.superox_indices[i]=1;
            } else {
              cce_vars.oxidation_states[i]=-1;
              if(cce_flags.flag("COMMAND_LINE")){
                cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to O-1 " << endl;
              }
              cce_vars.perox_indices[i]=1;
            }
          }
        }
        cce_vars.num_superox_bonds=1*num_formula_units_in_cell; // 1 superox. bonds per formula unit * num. formula units
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "Number of superoxide bonds set to: " << cce_vars.num_superox_bonds << endl;
        }
        cce_vars.num_perox_bonds=0.5*num_formula_units_in_cell; // 0.5 perox. bonds per formula unit * num. formula units; this should neve yield a non-integer since there should be at least one complete peroxide ion per cell
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "Number of peroxide bonds set to: " << cce_vars.num_perox_bonds << endl;
        }
        if(cce_flags.flag("COMMAND_LINE")){
          cout << endl;
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) { // this case should never occur since for an alkali metal sesquioxide set oxidation states should always work
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        } else {
          cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
        }
      }
    }
  }

  //CCE_try_all_oxidation_states////////////////////////////////////////////////////////
  // try to determine the oxidation numbers using all known oxidation states for all cation species
  void CCE_try_all_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << "Trying all known oxidation numbers:" << endl;
    }
    // determine possible_ox_states vector of vectors of doubles needed for Marco's function of my algorithm to determine oxidation numbers
    vector<vector<double> > all_ox_states_electronegativity_sorted(structure.species.size());
    for(uint i=0,isize=structure.species.size();i<isize;i++){ 
      vector<string> all_ox_states_tokens(cce_vars.num_all_ox_states_electronegativity_sorted[i]);
      aurostd::string2tokens(cce_vars.all_ox_states_strings_electronegativity_sorted[i], all_ox_states_tokens, ",");
      for(uint j=0,jsize=cce_vars.num_all_ox_states_electronegativity_sorted[i];j<jsize;j++){ 
        all_ox_states_electronegativity_sorted[i].push_back(aurostd::string2utype<double>(all_ox_states_tokens[j]));
      }
    }
    // use Marco's implementation of my algorithm to determine oxidation_numbers
    CCE_determine_cation_oxidation_states(all_ox_states_electronegativity_sorted, structure, cce_vars);
  }

  //CCE_determine_cation_oxidation_states////////////////////////////////////////////////////////
  // ME Nov. 2019
  // for avoiding recursion algorithm to determine cation oxidation_numbers
  // determine the cation oxidation numbers by using possible oxidation states for each species 
  // which can be either the preferred or all known oxidation numbers
  void CCE_determine_cation_oxidation_states(const vector<vector<double> >& possible_ox_states, const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    uint num_cation_species = cce_vars.cations_map.size(); // the number of cation species
    uint natoms = cce_vars.oxidation_states.size();

    // Initialize
    for (uint i = 0; i < num_cation_species; i++) { // loop over all cations
      for (uint j = 0; j < cce_vars.cations_map[i].size(); j++) { // loop over atoms of the ith cation type that is given by the second index of cation_map
        if (cce_vars.multi_anion_atoms[cce_vars.cations_map[i][j]] != 1){ // exclude atoms that have been identified as multi anion atoms previously
          if(LDEBUG){
            cerr << "Updating oxidation number for atom " << cce_vars.cations_map[i][j] << " (" << structure.atoms[cce_vars.cations_map[i][j]].name << ") to " << possible_ox_states[i][0] << endl;
          }
          cce_vars.oxidation_states[cce_vars.cations_map[i][j]] = possible_ox_states[i][0]; // possible_ox_states (either preferred or all) for all species should be electronegativity sorted
        }
      }
    }
    if(LDEBUG){
      for(uint n=0,nsize=structure.atoms.size();n<nsize;n++){ 
        cerr << "chosen oxdiation state for atom " << n << " (" << structure.atoms[n].name << "): " << cce_vars.oxidation_states[n] << endl;
      }
    }
    // check
    double total_ox = 0.0;
    for (uint at = 0; at < natoms; at++) total_ox += cce_vars.oxidation_states[at];
    bool iszero = aurostd::isequal(total_ox, 0.0, _CCE_OX_TOL_);
    // indicate whether these oxidation states were successful or not
    if (iszero) {
      if(LDEBUG){
        cerr << "Oxidation numbers successfully determined since oxidation sum is " << total_ox << endl;
      }
    } else{
      if(LDEBUG){
        cerr << "No successful determination of oxidation numbers since oxidation sum is " << total_ox << endl;
      }
    }

    // Cycle through oxidation states
    if (!iszero) { // maybe the oxidation numbers from the initialization work already, then nothing else needs to be tried
      uint maxsize = 0; // maximum number of preferred or all oxidation states for a cation species in the system
      for (uint j = 0; j < num_cation_species; j++) {
        if (possible_ox_states[j].size() > maxsize) maxsize = possible_ox_states[j].size();
      }

      uint n = 0;
      uint i = 1; // can start at 1 since 0th pref./all ox states has been used during initialization
      while (!iszero && (i < maxsize)) {
        for (uint j = num_cation_species; j > 0; j--) { // have to loop downward since species are sorted ascending by electronegativity and oxidation number should be changed for less electronegative species first
          n = possible_ox_states[j-1].size(); // j-1 since species loop goes downward and j_min should be comapred to be > 0 for uint, how many pref./all ox. states there are for the jth cation species
          if (i < n) { // the index i looping over all possible ox states for species j must be smaller than the maximum number of pref./all ox states for this species 
            for (uint k = 0; k < cce_vars.cations_map[j-1].size(); k++) {
              if (cce_vars.multi_anion_atoms[cce_vars.cations_map[j-1][k]] != 1){ // exclude atoms that have been identified as multi anion atoms previously
                if(LDEBUG){
                  cerr << "Updating oxidation number for atom " << cce_vars.cations_map[j-1][k] << "(" << structure.atoms[cce_vars.cations_map[j-1][k]].name << ") to " << possible_ox_states[j-1][i] << endl;
                }
                cce_vars.oxidation_states[cce_vars.cations_map[j-1][k]] = possible_ox_states[j-1][i]; // i goes over all preferred/all oxidation states
              }
            }
            if(LDEBUG){
              for(uint n=0,nsize=structure.atoms.size();n<nsize;n++){ 
                cerr << "chosen oxidation state for atom " << n << " (" << structure.atoms[n].name << "): " << cce_vars.oxidation_states[n] << endl;
              }
            }
            // check
            total_ox = 0.0;
            for (uint at = 0; at < natoms; at++) total_ox += cce_vars.oxidation_states[at];
            iszero = aurostd::isequal(total_ox, 0.0, _CCE_OX_TOL_);
            // indicate whether these oxidation states were successful or not
            if (iszero) {
              if(LDEBUG){
                cerr << "Oxidation numbers successfully determined since oxidation sum is " << total_ox << endl;
              }
              break; // BREAK for loop over possible oxidation states k for species j-1 upon success (sum over oxidation numbers is zero)
            } else{
              if(LDEBUG){
                cerr << "No successful determination of oxidation numbers since oxidation sum is " << total_ox << endl;
              }
            }
          }
        }
        i++;
      }
    }
  }

  //CCE_print_oxidation_states_and_sum////////////////////////////////////////////////////////
  // print out previously determined oxidation numbers and sum
  void CCE_print_oxidation_states_and_sum(CCE_Variables& cce_vars, const xstructure& structure) {
    aurostd::xoption cce_flags = CCE_init_flags();
    cce_flags.flag("COMMAND_LINE",TRUE);
    // double returned by following function call is not needed
    CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
  }

  //CCE_print_oxidation_states_and_get_sum////////////////////////////////////////////////////////
  // print out previously determined oxidation numbers and get sum
  double CCE_print_oxidation_states_and_get_sum(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    // print oxidation numbers
    if(cce_flags.flag("COMMAND_LINE")){
      cout << "OXIDATION NUMBERS:" << endl;
      for (uint k=0,ksize=structure.atoms.size();k<ksize;k++){
        cout << "Oxidation state of " << structure.atoms[k].name << " (atom[" << k << "]): " << cce_vars.oxidation_states[k] << endl;
      }
      cout << "CHECK whether this is what you are expecting!" << endl;
    }
    // calculate and print sum of oxidation numbers
    cce_vars.oxidation_sum = CCE_get_oxidation_states_sum(cce_vars, cce_flags);
    return cce_vars.oxidation_sum;
  }

  //CCE_get_oxidation_states_sum////////////////////////////////////////////////////////
  // Calculate and return the sum of all oxidation states without writing output.
  double CCE_get_oxidation_states_sum(CCE_Variables& cce_vars) {
    aurostd::xoption cce_flags = CCE_init_flags();
    cce_flags.flag("COMMAND_LINE",FALSE);
    return CCE_get_oxidation_states_sum(cce_vars, cce_flags);
  }

  //CCE_get_oxidation_states_sum////////////////////////////////////////////////////////
  // Calculate and return the sum of all oxidation states.
  double CCE_get_oxidation_states_sum(CCE_Variables& cce_vars, xoption& cce_flags) {
    cce_vars.oxidation_sum=0;
    for (uint k=0,ksize=cce_vars.oxidation_states.size();k<ksize;k++){
      cce_vars.oxidation_sum += cce_vars.oxidation_states[k];
    }
    if(cce_flags.flag("COMMAND_LINE")){
      cout << "Sum over all oxidation numbers is (should be ZERO): " << cce_vars.oxidation_sum << endl;
    }
    return cce_vars.oxidation_sum;
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //              DETERMINE OXIDATION STATES FROM BADER CHARGES              //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // determine oxidation numbers based on Bader charges (outdated)
  // functions for treating known special cases not included for the electronegativity based implementation are included here

  //CCE_get_oxidation_states_from_Bader////////////////////////////////////////////////////////
  // determine the oxidation numbers of the ions by an analysis of the Bader charges of the system and handle known exceptional cases explicitly
  vector<double> CCE_get_oxidation_states_from_Bader(xoption& cce_flags, const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << "DETERMINATION OF OXIDATION NUMBERS FROM BADER CHARGES:" << endl;
    }
    // in this part functional dependent settings will be primarily evaluated with "if, else if" conditions, 
    // since the calculation (for getting the Bader charges) could have only been done with one functional and not more than one
    // oxidation numbers should not be functional dependent since results are for a specific calculation, i.e. PBE or LDA or SCAN or PBE+U_ICSD and not more than one
    // for a given structure there should be only one correct assignment of oxidation numbers
    string system_name;
    string functional;
    // check whether aflow.in exists
    if (aurostd::FileExist(_AFLOWIN_)) {
      // get system name to identify Bader charges file and functional to distinguish corrections to be loaded from aflow.in
      CCE_get_system_name_functional_from_aflow_in(system_name, cce_flags, functional, cce_vars);
      if (cce_flags.flag("CORRECTABLE")){
        // check whether Bader file exists
        string Bader_file = system_name + "_abader.out";
        if (aurostd::FileExist(Bader_file) || aurostd::EFileExist(Bader_file)) {
          // determine Bader charges from Bader file and store in array/vector Bader_charges; O Bader charges will be included although they might not be used
          cce_vars.Bader_charges = CCE_get_Bader_charges_from_Bader_file(Bader_file, structure, cce_vars);
          // determine oxidation numbers from Bader charges
          cce_vars.oxidation_states = CCE_Bader_charges_to_oxidation_states(cce_flags, structure, cce_vars, functional);
          // check whether sum of oxidation numbers is equal zero and try to correct for known special cases (for some cases even if it is zero)
          if (cce_flags.flag("CORRECTABLE")){
            // print oxidation numbers and calculate sum
            cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
            // SECOND point to deal with special cases known from (binary+ternary) oxides
            // Ti-O Magneli phases need to be treated specially since oxidation numbers are not recognized appropriately for all functionals
            cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
            CCE_treat_Ti_O_Magneli_phase_special_case(cce_vars, structure, cce_flags);
            // for Fe3O4 in inverse spinel structure, the oxidation states are not identified properly via the Bader charges.
            CCE_treat_Fe3O4_special_case(cce_vars, structure, cce_flags);
            // for Mn3O4 in spinel structure, the oxidation states are not identified properly
            CCE_treat_Mn3O4_special_case(cce_vars, structure, cce_flags);
            // for Co3O4 in spinel structure, the oxidation states are not identified properly
            CCE_treat_Co3O4_special_case(cce_vars, structure, cce_flags);
            // MnMoO4 needs to be treated specially since oxidation numbers are not recognized appropriately
            CCE_treat_MnMoO4_special_case(cce_vars, structure, cce_flags);
            if (functional == "LDA") {
              // Ca2Fe2O5 & CaFe2O4 need to be treated specially for LDA since oxidation numbers are not recognized appropriately
              CCE_treat_Ca2Fe2O5_CaFe2O4_LDA_special_case(cce_vars, structure, cce_flags);
              // FeTiO3 needs to be treated specially for LDA since oxidation numbers are not recognized appropriately
              CCE_treat_FeTiO3_LDA_special_case(cce_vars, structure, cce_flags);
            }
            if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
              // general scheme to repair wrong oxidation numbers based on changing oxidation state of several ions known to be problematic
              CCE_general_attempt_fixing_oxidation_states(cce_vars, structure, cce_flags);
              // print oxidation numbers and calculate sum
              cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
              // system should not be regarded correctable if sum over oxidation states is not zero
              if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
                cce_flags.flag("CORRECTABLE",FALSE);
                cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
                cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
              }
            }
            if(cce_flags.flag("COMMAND_LINE")){
              cout << endl; // empty line to separate formation enthalpy values from the rest
            }
          }
        } else {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "Bader file " << Bader_file << " (or xz, bz2, gz version) not found. A Bader file is required to determine the oxidation numbers." << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
          cerr << endl;
        }
      }
    } else {
      cce_flags.flag("CORRECTABLE",FALSE);
      cerr << "aflow.in file not found. An aflow.in file is required to identify the functional and to find the Bader charges file to determine the oxidation numbers." << endl;
      cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
      cerr << endl;
    }
    return cce_vars.oxidation_states;
  }

  //CCE_get_system_name_functional_from_aflow_in////////////////////////////////////////////////////////
  // determine the system name and functional from the aflow_in
  void CCE_get_system_name_functional_from_aflow_in(string& system_name, xoption& cce_flags, string& functional, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string aflowIn = aurostd::RemoveComments(aurostd::file2string(_AFLOWIN_));
    vector<string> vlines = aurostd::string2vectorstring(aflowIn);
    string line_a;
    // define local functional variable to check whether Bader data is for different functional than asked for in vfunctionals and throw warning accordingly
    functional="PBE";
    bool ldau = false;
    for (uint i = 0; i < vlines.size(); i++) { 
      line_a = aurostd::RemoveSpaces(vlines[i]);
      if (line_a.find("[AFLOW]SYSTEM=") != string::npos){ // string::npos is returned if string is not found
        system_name = aurostd::RemoveSubStringFirst(line_a, "[AFLOW]SYSTEM=");
        system_name=aurostd::RemoveSubString(system_name,"\r"); // remove carriage return characters at end of string in case they exist (maybe due to sending aflow.in via email); leads to problems with string addition
        if (system_name.find("_ICSD_") != string::npos){ // needs maybe more reliable check based on DFT+U parameters in the future
          functional="PBE+U_ICSD";
        }
      }
      if ((line_a.find("LDAU") != string::npos && line_a.find("=ON") != string::npos) || (line_a.find("LDAU") != string::npos && line_a.find("=ADIABATIC") != string::npos)){
        ldau = true;
        if (functional != "PBE+U_ICSD"){
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << endl;
          cerr << "BAD NEWS: It seems you are providing an aflow.in for a DFT+U calculation with different paramters than for the AFLOW ICSD database. There are no corrections for this case." << endl;
          cerr << endl;
        }
      }
    }
    if (functional != "PBE+U_ICSD"){
      if (!ldau){
        for (uint i = 0; i < vlines.size(); i++) {
          line_a = aurostd::RemoveSpaces(vlines[i]);
          if (line_a.find("=potpaw_LDA") != string::npos || line_a.find("/potpaw_LDA/") != string::npos){ // the first criterion seems to not find "=potpaw_LDA" when the line is commented by # due to RemoveComments above
            functional="LDA";
          } else if (line_a.find("METAGGA=SCAN") != string::npos){
            functional="SCAN";
          }
        }
      }
    }
    if(LDEBUG){
      cerr << "functional determined from aflow.in: " << functional << endl;
      cerr << "PBE: " << aurostd::withinList(cce_vars.vfunctionals, "PBE") << endl;
      cerr << "LDA: " << aurostd::withinList(cce_vars.vfunctionals, "LDA") << endl;
      cerr << "SCAN: " << aurostd::withinList(cce_vars.vfunctionals, "SCAN") << endl;
      cerr << "PBE+U_ICSD: " << aurostd::withinList(cce_vars.vfunctionals, "PBE+U_ICSD") << endl;
      cerr << "exp: " << aurostd::withinList(cce_vars.vfunctionals, "exp") << endl;
    }
    // if functional determined from aflow.in is different from the ones given by the input options, 
    // throw warning that oxidation numbers are only determined on the basis of a specific functional
    // also when only the formation enthalpy from the exprimental values per bond shall be obtained,
    // a warning should be given from which functional the data are used to determine the ox nums.
    if(cce_flags.flag("COMMAND_LINE")){
      cout << endl;
      if(!(cce_vars.vfunctionals.size() == 1 && cce_vars.vfunctionals[0] == functional)){
        if (functional == "LDA") {
          cout << "WARNING: The oxidation numbers are only determined on the basis of an " << functional << " calculation." << endl;
        } else {
          cout << "WARNING: The oxidation numbers are only determined on the basis of a " << functional << " calculation." << endl;
        }
      }
    }
  }

  //CCE_get_Bader_charges_from_Bader_file////////////////////////////////////////////////////////
  // determine the Bader charges from the Bader file
  vector<double> CCE_get_Bader_charges_from_Bader_file(const string& Bader_file, const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    vector<string> vlines; 
    aurostd::efile2vectorstring(Bader_file,vlines);
    vector<string> tokens;
    for (uint i = 1; i < vlines.size(); i++) {
      aurostd::string2tokens(vlines[i], tokens, " ");
      cce_vars.Bader_charges.push_back(aurostd::string2utype<double>(tokens[2]));
    }
    for (uint k=0,ksize=structure.atoms.size();k<ksize;k++){ // there should always be as many Bader charges as atoms in the structure
      if(LDEBUG){
        cerr << "Bader_charges[" << k << "]: " << cce_vars.Bader_charges[k] << endl;
      }
    }
    return cce_vars.Bader_charges;
  }

  //CCE_Bader_charges_to_oxidation_states////////////////////////////////////////////////////////
  // compare Bader charges to template values from fitting set and set oxidation numbers accordingly
  vector<double> CCE_Bader_charges_to_oxidation_states(xoption& cce_flags, const xstructure& structure, CCE_Variables& cce_vars, string& functional) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::CCE_Bader_charges_to_oxidation_states():";
    stringstream message;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].name != cce_vars.anion_species){
        string Bader_templ_line;
        if (CCE_get_Bader_templates(structure.atoms[i].name) == "") {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "VERY BAD NEWS: There is no correction for " << structure.atoms[i].name << " (ATOM[" << i << "])" << " since this species was not included in the set for deducing corrections!"  << endl;
          cerr << endl;
        } else {
          Bader_templ_line=CCE_get_Bader_templates(structure.atoms[i].name);
          if(LDEBUG){
            cerr << "Bader templates: " << Bader_templ_line << endl;
          }
          vector<string> Bader_tokens;
          aurostd::string2tokens(Bader_templ_line, Bader_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. Bader_tokens[1] is not a space but the oxidation number
          uint num_ox_states = aurostd::string2utype<uint>(Bader_tokens[0]);
          if(LDEBUG){
            cerr << "number of oxidation states for which Bader charges are available: " << num_ox_states << endl;
          }
          double Bader_template=0; // Bader cation charge from the binary oxide from which the correction was deduced
          double Bader_tolerance=0.26; // tolerance with which the Bader charge of each atom is allowed to deviate from any of the template values to still assign the correction (and oxidation number)
          double Bader_deviation=0.5; // deviation of the Bader charge of the cation from the template value; initial value safely larger than Bader tolerance for identifying oxidation number so that if value is not changed, no correction is made
          double Bader_deviation_0=Bader_deviation; // initial Bader deviation for later check whether any corrections were found, i. e. Bader deviation was changed
          double Bader_deviation_min=17; // initialize with high value; used later to output smallest Bader deviation if no correction according to the Bader tolerance was identified
          for(uint n=0;n<num_ox_states;n++){ //loop over all oxidation number for which Bader charges are available and read the Bader charges from the respective positions
            // here the functional determined from the aflow.in is decisive; exp should not occur
            vector<string> vfunctional_aflow_in;
            vfunctional_aflow_in.push_back(functional);
            if (CCE_get_offset(functional) == -1) {
              message << "Unknown functional " << functional << ". Please choose PBE, LDA, or SCAN.";
              throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
            }
            int offset=CCE_get_offset(functional); // define also local offset variable since offset must correspond to functional detected from aflow.in
            for (uint i = 0; i < vfunctional_aflow_in.size(); i++) {
              Bader_template = aurostd::string2utype<double>(Bader_tokens[offset/2+2+5*n]);
              if(LDEBUG){
                cerr << "Bader_template " << vfunctional_aflow_in[i] << ": " << Bader_template << endl;
              }
            }
            if ( std::abs(Bader_template-cce_vars.Bader_charges[i]) < Bader_tolerance && std::abs(Bader_template-cce_vars.Bader_charges[i]) < Bader_deviation ){ // must be compared to Bader_deviation to load correction for Bader charge closest to template and not to be overloaded by other corrections for later tested Bader charges (oxidation states)
              Bader_deviation= std::abs(Bader_template-cce_vars.Bader_charges[i]);
              if(LDEBUG){
                cerr << "Bader_deviation: " << Bader_deviation << endl;
              }
              cce_vars.oxidation_states[i]= aurostd::string2utype<double>(Bader_tokens[1+5*n]);
            } else { // only if element is found but no corrections can be assigned because above conditions are not met, update Bader_deviation_min
              if ( std::abs(Bader_template-cce_vars.Bader_charges[i]) < Bader_deviation_min ) {
                Bader_deviation_min= std::abs(Bader_template-cce_vars.Bader_charges[i]);
              }
            }
          }
          if ( Bader_deviation == Bader_deviation_0 ){ // only if no oxidation state was found, i. e. Bader deviation was not changed but there are corrections (Bader charges) for this species (CCE_get_Bader_templates returns non-empty string), display the following error message after handling the special cases 
            // FIRST point to deal with special cases known from (ternary) oxides; other cases will be dealt with when looking at the oxidation numbers and checking the sum
            // since it can happen that the Bader charge for W in the compound (especially for W+6) is too far away 
	    // from the Bader charge template, here W+6 is set and later the sum over the oxidation states will still be checked
            if (structure.atoms[i].name == "W") {
              cce_vars.oxidation_states[i]=+6;
            // since it can happen that the Bader charge for Pb in the compound (especially for Pb+2) is too far away 
	    // from the Bader charge template, here Pb+2 is set and later the sum over the oxidation states will still be checked
            } else if (structure.atoms[i].name == "Pb") {
              cce_vars.oxidation_states[i]=+2;
            // since it can happen that the Bader charge for Ag in the compound (especially for Ag+1) is too far away 
	    // from the Bader charge template, here Ag+1 is set and later the sum over the oxidation states will still be checked
            } else if (structure.atoms[i].name == "Ag") {
              cce_vars.oxidation_states[i]=+1;
            } else {
              cce_flags.flag("CORRECTABLE",FALSE);
              cerr << "BAD NEWS: The oxidation number (and hence the correction) for " << structure.atoms[i].name << " (ATOM[" << i << "])" << " cannot be identified from the Bader charges!"  << endl;
              cerr << "The deviation of the Bader charge from the closest tested template value is: " << Bader_deviation_min << " electrons. This is larger than the tolerance: " << Bader_tolerance  << " electrons." << endl;
              // list all oxidation states of the element for which corrections are available
              string ox_nums_avail="";
              string separator=", ";
              for(uint n=0;n<num_ox_states;n++){ 
                if (n<num_ox_states-1){
                  ox_nums_avail+= Bader_tokens[1+5*n] + separator ; 
                } else if (n==num_ox_states-1){
                  ox_nums_avail+= Bader_tokens[1+5*n]; 
                }
              }
              cerr << "Corrections for " << structure.atoms[i].name << " coordinated by " << cce_vars.anion_species << " are available for oxidation states: " << ox_nums_avail << endl;
              cerr << "If the desired oxidation state is listed but it is just not correctly determined from the Bader charges," << endl;
              cerr << "you might want to consider supplying the oxidation numbers manually by using the option --oxidation_numbers=." << endl;
              cerr << endl;
            }
          }
        }
      } else if (structure.atoms[i].name == cce_vars.anion_species) {
        cce_vars.oxidation_states[i] = cce_vars.standard_anion_charge; // oxidation numbers for O are assumed to be -2 (other anions accordingly) and are corrected below if it is a per- or superoxide O atom as identified from the structural analysis in other functions
        if (cce_vars.num_perox_bonds > 0){
          for (uint j=0,jsize=structure.atoms.size();j<jsize;j++){
            if (cce_vars.perox_indices[j]==1 && j == i){
              cce_vars.oxidation_states[i]=-1;
            }
          }
        }
        if (cce_vars.num_superox_bonds > 0){
          for (uint j=0,jsize=structure.atoms.size();j<jsize;j++){
            if (cce_vars.superox_indices[j]==1 && j == i){
              cce_vars.oxidation_states[i]=-0.5;
            }
          }
        }
      }
    }
    return cce_vars.oxidation_states;
  }

  //CCE_treat_MnMoO4_special_case////////////////////////////////////////////////////////
  // for MnMoO4 the oxidation numbers are determined to be +4 for both Mn and Mo from the Bader charges; 
  // it should be Mn+2 & Mo+6; however since the sum of the falsely determined oxidation numbers 
  // is accidentally 0 it needs to be corrected individually
  void CCE_treat_MnMoO4_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if ( structure.species.size() == 3 && structure.species[0] == "Mn" && structure.species[1] == "Mo" && structure.species[2] == "O" ) {
      double amount_Mn=0;
      double amount_Mo=0;
      double amount_O=0;
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].name == "Mn"){
          amount_Mn+=1;
        } else if (structure.atoms[i].name == "Mo"){
          amount_Mo+=1;
        } else if (structure.atoms[i].name == "O"){
          amount_O+=1;
        }
      }
      if(LDEBUG){
        cerr << "number of Mn ions= " << amount_Mn << endl;
        cerr << "number of Mo ions= " << amount_Mo << endl;
        cerr << "number of O ions= " << amount_O << endl;
      }
      if (amount_Mn/amount_Mo == 1 && amount_Mn/amount_O == 0.25 && amount_Mo/amount_O == 0.25) {
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "MnMoO4 special treatment since sum over oxdiation states is zero but individual oxidation numbers are wrong!!!" << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "Mn"){
            cce_vars.oxidation_states[i]=+2;
            if(cce_flags.flag("COMMAND_LINE")){
             cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Mn+2 " << endl;
            }
          }
          if (structure.atoms[i].name == "Mo"){
            cce_vars.oxidation_states[i]=+6;
            if(cce_flags.flag("COMMAND_LINE")){
              cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Mo+6 " << endl;
            }
          }
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        }
      }
    }
  }

  //CCE_treat_Ca2Fe2O5_CaFe2O4_LDA_special_case////////////////////////////////////////////////////////
  // for Ca2Fe2O5 and CaFe2O4 for LDA the oxidation numbers of Fe are not correctly determined 
  // to be Fe+3 but are partly Fe+2 which will be corrected here
  void CCE_treat_Ca2Fe2O5_CaFe2O4_LDA_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if ( structure.species.size() == 3 && structure.species[0] == "Ca" && structure.species[1] == "Fe" && structure.species[2] == "O" ) {
      double amount_Ca=0;
      double amount_Fe=0;
      double amount_O=0;
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].name == "Ca"){
          amount_Ca+=1;
        } else if (structure.atoms[i].name == "Fe"){
          amount_Fe+=1;
        } else if (structure.atoms[i].name == "O"){
          amount_O+=1;
        }
      }
      if(LDEBUG){
        cerr << "number of Ca ions= " << amount_Ca << endl;
        cerr << "number of Fe ions= " << amount_Fe << endl;
        cerr << "number of O ions= " << amount_O << endl;
      }
      //making sure it is Ca2Fe2O5
      if (amount_Ca/amount_Fe == 1 && amount_Ca/amount_O == 0.4 && amount_Fe/amount_O == 0.4) {
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "Ca2Fe2O5 special treatment for LDA since oxidation numbers for Fe, which should be Fe+3, are not correctly determined from Bader charges for all Fe!!!" << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "Fe"){
            cce_vars.oxidation_states[i]=+3;
            if(cce_flags.flag("COMMAND_LINE")){
              cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Fe+3 " << endl;
            }
          }
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        }
      //making sure it is CaFe2O4
      } else if (amount_Ca/amount_Fe == 0.5 && amount_Ca/amount_O == 0.25 && amount_Fe/amount_O == 0.5) {
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "CaFe2O4 special treatment for LDA since oxidation numbers for Fe, which should be Fe+3, are not correctly determined from Bader charges for all Fe!!!" << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "Fe"){
            cce_vars.oxidation_states[i]=+3;
            if(cce_flags.flag("COMMAND_LINE")){
              cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Fe+3 " << endl;
            }
          }
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        }
      }
    }
  }

  //CCE_treat_FeTiO3_LDA_special_case////////////////////////////////////////////////////////
  // for FeTiO3 for LDA the oxidation numbers of Ti are not correctly determined to be Ti+4 and using 
  // the general fixes would modify both the Ti AND the Fe oxidation numbers resulting again 
  // in non-zero oxidation number sum, which is fixed here
  void CCE_treat_FeTiO3_LDA_special_case(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if ( structure.species.size() == 3 && structure.species[0] == "Fe" && structure.species[1] == "O" && structure.species[2] == "Ti" ) {
      double amount_Fe=0;
      double amount_O=0;
      double amount_Ti=0;
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].name == "Fe"){
          amount_Fe+=1;
        } else if (structure.atoms[i].name == "O"){
          amount_O+=1;
        } else if (structure.atoms[i].name == "Ti"){
          amount_Ti+=1;
        }
      }
      if(LDEBUG){
        cerr << "number of Fe ions= " << amount_Fe << endl;
        cerr << "number of O ions= " << amount_O << endl;
        cerr << "number of Ti ions= " << amount_Ti << endl;
      }
      //making sure it is FeTiO3
      if (amount_Fe/amount_Ti == 1 && amount_Fe/amount_O == 1.0/3 && amount_Ti/amount_O == 1.0/3) {
        if(cce_flags.flag("COMMAND_LINE")){
          cout << "FeTiO3 special treatment for LDA since oxidation numbers for Ti, which should be Ti+4, are not correctly determined from Bader charges and using other fixing would also change Fe oxidation numbers!!!" << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].name == "Ti"){
            cce_vars.oxidation_states[i]=+4;
            if(cce_flags.flag("COMMAND_LINE")){
              cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Ti+4 " << endl;
            }
          }
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = CCE_print_oxidation_states_and_get_sum(cce_vars, structure, cce_flags);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        }
      }
    }
  }

  //CCE_general_attempt_fixing_oxidation_states////////////////////////////////////////////////////////
  // for some ions known to be problematic for the Bader analysis (e.g. V+5, Fe+2, Fe+3, Ti+4) 
  // a brute force ansatz can be implemented to just change their oxidation state and later check 
  // whether it fixes the oxidation number sum rule
  void CCE_general_attempt_fixing_oxidation_states(CCE_Variables& cce_vars, const xstructure& structure, xoption& cce_flags) {
    if(cce_flags.flag("COMMAND_LINE")){
      cout << "The sum over all oxidation numbers for all atoms of the system is NOT zero, trying to repair that based on known problematic cases (Ti, V, Fe). This might or might not work." << endl;
    }
    // repairing only by considering non mixed valence oxides
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ 
      if (structure.atoms[i].name == "Ti"){
        if ( cce_vars.oxidation_states[i] == 3){
          cce_vars.oxidation_states[i]=+4;
          if(cce_flags.flag("COMMAND_LINE")){
	    cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Ti+4 " << endl;
          }
        }
      }
      if (structure.atoms[i].name == "Fe"){
        if ( cce_vars.oxidation_states[i] == 2){
          cce_vars.oxidation_states[i]=+3;
          if(cce_flags.flag("COMMAND_LINE")){
	    cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Fe+3 " << endl;
          }
        } else if ( cce_vars.oxidation_states[i] == 3){
          cce_vars.oxidation_states[i]=+2;
          if(cce_flags.flag("COMMAND_LINE")){
	    cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to Fe+2 " << endl;
          }
        }
      }
      if (structure.atoms[i].name == "V"){
        if ( cce_vars.oxidation_states[i] == 4){
          cce_vars.oxidation_states[i]=+5;
          if(cce_flags.flag("COMMAND_LINE")){
	    cout << "setting oxidation state of " << structure.atoms[i].name << " (atom[" << i << "]) to V+5 " << endl;
          }
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                            ASSIGN CORRECTIONS                           //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
 
  //CCE_get_corrections////////////////////////////////////////////////////////
  // determine the corrections after oxidation numbers are known (from input or Bader)
  void CCE_get_corrections(xoption& cce_flags, const xstructure& structure, vector<vector<double> >& corrections_atom, CCE_Variables& cce_vars, const vector<uint>& num_neighbors, const string& considered_anion_species) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if(LDEBUG){
      cerr << "ASSIGNMENT OF CORRECTIONS:" << endl;
    }
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].name != cce_vars.anion_species && cce_vars.multi_anion_atoms[i] != 1){ // exclude main anion species and multi anion atoms detected previously
        string corrections_line;
        if (num_neighbors[i] > 0){ // are there actually bonds between the cation and the anion species under consideration (to load only needed corrections in multi-anion systems)
          if ( CCE_get_corrections_line(structure.atoms[i].name + "_+" + aurostd::utype2string<double>(cce_vars.oxidation_states[i]) + "_" + considered_anion_species) == "") { // the considered anion species can be the main anion species or a multi anion species
            cce_flags.flag("CORRECTABLE",FALSE);
            cerr << "BAD NEWS: No correction available for " << structure.atoms[i].name << " (ATOM[" << i << "])" << " in oxidation state " << cce_vars.oxidation_states[i] << " when coordinated by " << considered_anion_species << "." << endl;
            //checking for which oxidation states corrections are available and throw out errors accordingly
            uint ox_nums_count=0;
            vector<uint> ox_nums_avail_vec;
            for(uint k=0,ksize=12;k<ksize;k++){ // larger than ox. num. +12 should not occur
              if ( CCE_get_corrections_line(structure.atoms[i].name + "_+" + aurostd::utype2string<uint>(k) + "_" + considered_anion_species) != "") {
                ox_nums_count+=1;
                ox_nums_avail_vec.push_back(k);
              }
            }
            if ( ox_nums_count == 0) {
              cerr << "Currently no corrections available for " << structure.atoms[i].name << " when coordinated by "<< considered_anion_species << "." << endl;
            } else {
              // list all oxidation states of the element for which corrections are available
              string ox_nums_avail="";
              string separator=", ";
              for(uint n=0;n<ox_nums_count;n++){ 
                if (n<ox_nums_count-1){
                  ox_nums_avail+= "+" + aurostd::utype2string<uint>(ox_nums_avail_vec[n]) + separator ; 
                } else if (n==ox_nums_count-1){
                  ox_nums_avail+= "+" + aurostd::utype2string<uint>(ox_nums_avail_vec[n]); 
                }
              }
              cerr << "Corrections for " << structure.atoms[i].name << " coordinated by " << considered_anion_species << " are available for oxidation states: " << ox_nums_avail << endl;
              cerr << "If the desired oxidation state is listed but it is just not correctly determined," << endl;
              cerr << "you might want to consider supplying the oxidation numbers manually by using the option --oxidation_numbers=." << endl;
            }
            cerr << endl;
          } else {
            corrections_line=CCE_get_corrections_line(structure.atoms[i].name + "_+" + aurostd::utype2string<double>(cce_vars.oxidation_states[i]) + "_" + considered_anion_species);
            if(LDEBUG){
              cerr << "corrections line: " << corrections_line << endl;
            }
            // load cation corrections
            CCE_load_cation_corrections(cce_vars, structure, corrections_atom, corrections_line, i);
          }
        }
      } else if (structure.atoms[i].name == cce_vars.anion_species || cce_vars.multi_anion_atoms[i] == 1) {
        // set anion corrections (to zero)
        CCE_set_anion_corrections(cce_vars, structure, corrections_atom, i);
      }
    }
  }

  //CCE_load_cation_corrections////////////////////////////////////////////////////////
  // load corrections for the cations from the corrections line according to the functional
  // here for every functional there will be a separate "if" evaluation, since when only a structure (+ oxidation numbers) 
  // are provided as input, one might want to get the corrections for more than one functional
  void CCE_load_cation_corrections(CCE_Variables& cce_vars, const xstructure& structure, vector<vector<double> >& corrections_atom, const string& corrections_line, uint i) { // also provide increment for loop of CCE_get_corrections function
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    vector<string> corrections_tokens;
    aurostd::string2tokens(corrections_line, corrections_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. corrections_tokens[1] is not a space but the oxidation number
    for (uint k = 0; k < cce_vars.vfunctionals.size(); k++) {
      corrections_atom[2*k][i]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+2]); // 2*k since to each functional belong 2 corrections for 298.15 and 0K
      if(LDEBUG){
        cerr << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].name << " (atom[" << i << "]) for 298.15K: " << corrections_atom[2*k][i] << " eV/bond" << endl;
      }
      if (cce_vars.vfunctionals[k] != "exp") {
        corrections_atom[2*k+1][i]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+3]);
        if(LDEBUG){
          cerr << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].name << " (atom[" << i << "]) for 0K: " << corrections_atom[2*k+1][i] << " eV/bond" << endl;
        }
      }
    }
  }

  //CCE_set_anion_corrections////////////////////////////////////////////////////////
  // set corrections for anion atoms to zero since only for cations number of bonds with anions are needed; 
  // per- & superoxides will be dealt with in other function
  void CCE_set_anion_corrections(CCE_Variables& cce_vars, const xstructure& structure, vector<vector<double> >& corrections_atom, uint i) { // also provide increment for loop of CCE_get_corrections function
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    for (uint k = 0; k < cce_vars.vfunctionals.size(); k++) {
      corrections_atom[2*k][i]=0; // 2*k since to each functional belong 2 corrections for 298.15 and 0K
      if(LDEBUG){
        cerr << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].name << " (atom[" << i << "]) for 298.15K: " << corrections_atom[2*k][i] << " eV/bond" << endl;
      }
      if (cce_vars.vfunctionals[k] != "exp") {
        corrections_atom[2*k+1][i]=0;
        if(LDEBUG){
          cerr << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].name << " (atom[" << i << "]) for 0K: " << corrections_atom[2*k+1][i] << " eV/bond" << endl;
        }
      }
    }
  }

  //CCE_check_get_per_super_ox_corrections////////////////////////////////////////////////////////
  // check whether corrections for per- or superoxide ions are needed and assign accordingly
  void CCE_check_get_per_super_ox_corrections(CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    // check for per- and superoxides should always be kept separate since a compound could have both per- and superoxide ions (alkali-metal sesquioxides)
    if (cce_vars.num_perox_bonds > 0){
      string corrections_line;
      corrections_line=CCE_get_corrections_line("O2_-2_O");
      if(LDEBUG){
        cerr << "corrections line peroxides: " << corrections_line << endl;
      }
      vector<string> corrections_tokens;
      aurostd::string2tokens(corrections_line, corrections_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. corrections_tokens[1] is not a space but the oxidation number
      for (uint k = 0; k < cce_vars.vfunctionals.size(); k++) {
        cce_vars.perox_correction[2*k]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+2]);
        if(LDEBUG){
          cerr << cce_vars.vfunctionals[k] << " peroxide correction for 298.15K: " << cce_vars.perox_correction[2*k] << " eV/bond" << endl;
        }
        if (cce_vars.vfunctionals[k] != "exp") {
          cce_vars.perox_correction[2*k+1]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+3]);
          if(LDEBUG){
            cerr << cce_vars.vfunctionals[k] << " peroxide correction for 0K: " << cce_vars.perox_correction[2*k+1] << " eV/bond" << endl;
          }
        }
      }
    } 
    if (cce_vars.num_superox_bonds > 0) {
      string corrections_line;
      corrections_line=CCE_get_corrections_line("O2_-1_O");
      if(LDEBUG){
        cerr << "corrections line superoxides: " << corrections_line << endl;
      }
      vector<string> corrections_tokens;
      aurostd::string2tokens(corrections_line, corrections_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. corrections_tokens[1] is not a space but the oxidation number
      for (uint k = 0; k < cce_vars.vfunctionals.size(); k++) {
        cce_vars.superox_correction[2*k]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+2]);
        if(LDEBUG){
          cerr << cce_vars.vfunctionals[k] << " superoxide correction for 298.15K: " << cce_vars.superox_correction[2*k] << " eV/bond" << endl;
        }
        if (cce_vars.vfunctionals[k] != "exp") {
          cce_vars.superox_correction[2*k+1]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+3]);
          if(LDEBUG){
            cerr << cce_vars.vfunctionals[k] << " superoxide correction for 0K: " << cce_vars.superox_correction[2*k+1] << " eV/bond" << endl;
          }
        }
      }
    }
    if(LDEBUG){
      cerr << endl;
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //         APPLY CORRECTIONS AND GET CORRECTED FORMATION ENTHALPIES        //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
 
  //CCE_check_apply_per_super_ox_corrections////////////////////////////////////////////////////////
  // check whether corrections for per- or superoxide ions are needed and apply accordingly
  void CCE_check_apply_per_super_ox_corrections(CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if (cce_vars.num_perox_bonds > 0){
      for (uint k = 0; k < cce_vars.vfunctionals.size(); k++) {
        cce_vars.cce_correction[2*k] += (cce_vars.num_perox_bonds * cce_vars.perox_correction[2*k]) ;
        if(LDEBUG){
          cerr << cce_vars.vfunctionals[k] << " peroxide correction for 298.15K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_perox_bonds * cce_vars.perox_correction[2*k]) << " eV" << endl;
        }
        if (cce_vars.vfunctionals[k] != "exp") {
          cce_vars.cce_correction[2*k+1] += (cce_vars.num_perox_bonds * cce_vars.perox_correction[2*k+1]) ;
          if(LDEBUG){
            cerr << cce_vars.vfunctionals[k] << " peroxide correction for 0K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_perox_bonds * cce_vars.perox_correction[2*k+1]) << " eV" << endl;
            cerr << endl;
          }
        }
      }
    }
    if (cce_vars.num_superox_bonds > 0){
      for (uint k = 0; k < cce_vars.vfunctionals.size(); k++) {
        cce_vars.cce_correction[2*k] += (cce_vars.num_superox_bonds * cce_vars.superox_correction[2*k]) ;
        if(LDEBUG){
          cerr << cce_vars.vfunctionals[k] << " superoxide correction for 298.15K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_superox_bonds * cce_vars.superox_correction[2*k]) << " eV" << endl;
        }
        if (cce_vars.vfunctionals[k] != "exp") {
          cce_vars.cce_correction[2*k+1] += (cce_vars.num_superox_bonds * cce_vars.superox_correction[2*k+1]) ;
          if(LDEBUG){
            cerr << cce_vars.vfunctionals[k] << " superoxide correction for 0K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_superox_bonds * cce_vars.superox_correction[2*k+1]) << " eV" << endl ;
            cerr << endl;
          }
        }
      }
    }
  }

  //CCE_get_formation_enthalpies//////////////////////////////////////////////
  // ME 200213
  // Returns the formation enthalpy per cell for each functional
  vector<double> CCE_get_formation_enthalpies(const vector<double>& cce_correction, CCE_Variables& cce_vars) {
    vector<double> formation_enthalpies(cce_correction.size(), 0.0);
    if (cce_vars.dft_energies.size() > 0) {
      for (uint i = 0; i < cce_vars.vfunctionals.size(); i++) {
        if (cce_vars.vfunctionals[i] == "exp" ) continue;
        // for 298.15 K
        formation_enthalpies[2*i] = cce_vars.dft_energies[i] - cce_correction[2*i];
        // for 0 K
        formation_enthalpies[2*i+1] = cce_vars.dft_energies[i] - cce_correction[2*i+1];
      }
    }
    // treat exp. separately since formation enthalpy is equal to sum of corrections in this case
    for (uint i = 0; i < cce_vars.vfunctionals.size(); i++) {
      if (cce_vars.vfunctionals[i] == "exp") {
        formation_enthalpies[2*i] = cce_correction[2*i];
      }
    }
    return formation_enthalpies;
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                        WRITE OUTPUT AND CITATION                        //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
 
  //CCE_get_JSON/////////////////////////////////////////////////////////////
  // ME 200213
  // Returns CCE results in JSON format
  string CCE_get_JSON(const xstructure& structure, const CCE_Variables& cce_vars) {
    stringstream json;
    bool print_Hf = (cce_vars.dft_energies.size() > 0);
    uint nfuncs = cce_vars.vfunctionals.size();
    uint natoms = structure.atoms.size();

    json << "{";
    json << "\"oxidation_states\":";
    json << "[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(cce_vars.oxidation_states),",") << "],";
    json << "\"CCE\":{";
    for (uint i = 0; i < nfuncs; i++) {
      json << "\"" << cce_vars.vfunctionals[i] << "\":{";
      json << "\"298.15K\":{";
      json << "\"cce_correction_cell\":" << cce_vars.cce_correction[2*i] << ",";
      json << "\"cce_correction_atom\":" << (cce_vars.cce_correction[2*i]/natoms);
      if (print_Hf) {
        json << ",";
        json << "\"formation_enthalpy_cell\":" << cce_vars.cce_form_energy_cell[2*i] << ",";
        json << "\"formation_enthalpy_atom\":" << (cce_vars.cce_form_energy_cell[2*i]/natoms);
      }
      json << "}";
      if (cce_vars.vfunctionals[i] != "exp") {
        json << ",\"0K\":{";
        json << "\"cce_correction_cell\":" << cce_vars.cce_correction[2*i+1] << ",";
        json << "\"cce_correction_atom\":" << (cce_vars.cce_correction[2*i+1]/natoms);
        if (print_Hf) {
          json << ",";
          json << "\"formation_enthalpy_cell\":" << cce_vars.cce_form_energy_cell[2*i+1] << ",";
          json << "\"formation_enthalpy_atom\":" << (cce_vars.cce_form_energy_cell[2*i+1]/natoms);
        }
        json << "}";
      }
      json << "}";
      if (i < nfuncs - 1) json << ",";
    }
    json << "}";
    json << "}";
    return json.str();
  }

  //CCE_write_output////////////////////////////////////////////////////////
  // write CCE corrections and corrected formation enthalpies if precalculated DFT values are provided
  string CCE_write_output(const xstructure& structure, CCE_Variables& cce_vars, const vector<double>& cce_form_energy_cell) {
    stringstream output;
    // print out CCE corrections per cell and atom for functionals selected
    if (!(cce_vars.vfunctionals.size() == 1 && cce_vars.vfunctionals[0] == "exp")){ // if only exp is set as functional CCE CORRECTIONS: should not be written
      output << "CCE CORRECTIONS (to be subtracted from precalculated DFT formation energies):" << endl;
    }
    for (uint k = 0; k < cce_vars.vfunctionals.size(); k++) {
      if (cce_vars.vfunctionals[k] != "exp") {
        output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k] << " eV/cell //" << "CCE@" << cce_vars.vfunctionals[k] << " correction for 298.15K." << endl;
        output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k+1] << " eV/cell //" << "CCE@" << cce_vars.vfunctionals[k] << " correction for 0K." << endl;
        output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k]/structure.atoms.size() << " eV/atom //" << "CCE@" << cce_vars.vfunctionals[k] << " correction for 298.15K." << endl;
        output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k+1]/structure.atoms.size() << " eV/atom //" << "CCE@" << cce_vars.vfunctionals[k] << " correction for 0K." << endl;
        output << endl;
      }
    }
    // exp result should always be written at the end, hence write only after writing output for other functionals
    for (uint k = 0; k < cce_vars.vfunctionals.size(); k++) {
      if (cce_vars.vfunctionals[k] == "exp" && cce_vars.dft_energies.size()==0) { // second condition for that if precalc. form. energies are given and asking for exp., exp. result is not written twice
        output << "CCE FORMATION ENTHALPIES:" << endl;
        output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k] << " eV/cell //" << "CCE@exp formation enthalpy at 298.15K from exp. formation enthalpies per bond." << endl;
        output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k]/structure.atoms.size() << " eV/atom //" << "CCE@exp formation enthalpy at 298.15K from exp. formation enthalpies per bond." << endl;
        output << "Note that this provides A ROUGH GUESS with an estimated average accuracy of only about 250 meV/atom (from test for ternary oxides)!" << endl;
        output << endl;
      }
    }
    // print CCE formation enthalpies per cell and atom for functionals selected 
    // if precalculated DFT values are provided
    if(cce_vars.dft_energies.size()!=0){ 
      output << "CCE FORMATION ENTHALPIES:" << endl;
      for (uint k = 0; k < cce_vars.vfunctionals.size(); k++) {
        if (cce_vars.vfunctionals[k] != "exp") {
          output << cce_form_energy_cell[2*k] << " eV/cell //" << "CCE@" << cce_vars.vfunctionals[k] << " formation enthalpy at 298.15K." << endl;
          output << cce_form_energy_cell[2*k+1] << " eV/cell //" << "CCE@" << cce_vars.vfunctionals[k] << " formation enthalpy at 0K." << endl;
          output << cce_form_energy_cell[2*k]/structure.atoms.size() << " eV/atom //" << "CCE@" << cce_vars.vfunctionals[k] << " formation enthalpy at 298.15K." << endl;
          output << cce_form_energy_cell[2*k+1]/structure.atoms.size() << " eV/atom //" << "CCE@" << cce_vars.vfunctionals[k] << " formation enthalpy at 0K." << endl;
          output << endl;
        } else if (cce_vars.vfunctionals[k] == "exp") {
          output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k] << " eV/cell //" << "CCE@exp formation enthalpy at 298.15K from exp. formation enthalpies per bond." << endl;
          output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k]/structure.atoms.size() << " eV/atom //" << "CCE@exp formation enthalpy at 298.15K from exp. formation enthalpies per bond." << endl;
          output << "Note that this provides A ROUGH GUESS with an estimated average accuracy of only about 250 meV/atom (from test for ternary oxides)!" << endl;
          output << endl;
        } 
      }
    }
    return output.str();
  }

  //CCE_write_citation////////////////////////////////////////////////////////
  // write citation information at end of output
  string CCE_write_citation() {
    stringstream oss;
    oss << "#########################################################################################################" << endl;
    oss << "When you obtain results using the CCE methodology and/or this implementation," << endl; 
    oss << "please cite the following article:" << endl;
    oss << "Friedrich et al., Coordination corrected ab initio formation enthalpies, npj Comput. Mater. 5, 59 (2019);" << endl; 
    oss << "https://doi.org/10.1038/s41524-019-0192-1" << endl;
    oss << "#########################################################################################################" << endl;
    oss << endl;
    return oss.str();
  }

  //CCE_print_usage////////////////////////////////////////////////////////
  // For printing user instructions if no additional input is provided.
  string CCE_print_usage() {
    stringstream oss;
    oss << endl;
    oss << "Written by Rico Friedrich, Corey Oses, and Marco Esters, 2019-2020" << endl;
    oss << endl;
    oss << "USER INSTRUCTIONS:" << endl;
    oss << endl;
    oss << "(i) GENERAL INFORMATION:" << endl;
    oss << "Implementation to obtain corrected DFT formation enthalpies based on the coordination corrected" << endl; 
    oss << "enthalpies (CCE) methodology described in:" << endl; 
    oss << "Friedrich et al., Coordination corrected ab initio formation enthalpies, npj Comput. Mater. 5, 59 (2019);" << endl;
    oss << "https://doi.org/10.1038/s41524-019-0192-1" << endl;
    oss << "Please cite this article when using this method and/or this implementation." << endl;
    oss << "The corrections depend on the number of cation-anion bonds and on the cation oxidation state." << endl;
    oss << "More general information and a list of elements and oxidation states for which corrections are available" << endl; 
    oss << "can be found in README_AFLOW_CCE.TXT." << endl;
    oss << endl;
    oss << "(ii) AVAILABLE OPTIONS:" << endl;
    oss << "--cce                       Prints these user instructions." << endl;
    oss << "--cce=POSCAR_FILE_PATH      Provide the path to the structure file in VASP5 POSCAR format." << endl; 
    oss << "--dft_formation_energies=   Provide a comma separated list of precalculated DFT formation energies," << endl; 
    oss << "                            they are assumed to be: (i) negative for compounds lower in energy" << endl; 
    oss << "                            than the elements, (ii) in eV/cell. Currently, corrections are available" << endl; 
    oss << "                            for PBE, LDA and SCAN." << endl;
    oss << "--functionals=              Provide a comma separated list of functionals for which corrections" << endl;
    oss << "                            should be returned. If used together with --dft_formation energies," << endl;
    oss << "                            the functionals must be in the same sequence as the DFT formation" << endl;
    oss << "                            energies they correspond to. Available functionals are:" << endl;
    oss << "                            (i) PBE, (ii) LDA or (iii) SCAN. Default: PBE (if only one DFT formation" << endl; 
    oss << "                            energy is provided)." << endl;
    oss << "--oxidation_numbers=        Provide as a comma separated list the oxidation numbers. It is" << endl;
    oss << "                            assumed that: (i) one is provided for each atom of the structure and" << endl; 
    oss << "                            (ii) they are in the same sequence as the corresponding atoms in the" << endl;
    oss << "                            provided POSCAR file." << endl;
    oss << "--poscar2cce < POSCAR       Determines the CCE corrections for the structure in file POSCAR (must" << endl;
    oss << "                            be in VASP5 POSCAR format)." << endl;
    oss << endl;
    oss << "(iii) EXAMPLE INPUT STRUCTURE FOR ROCKSALT MgO:" << endl;
    oss << "Mg1O1   [FCC,FCC,cF8] (STD_PRIM doi:10.1  [FCC,FCC,cF8] (STD_PRIM doi:10.1016/j.commatsci.2010.05.010)" << endl;
    oss << "1.224745" << endl;
    oss << "   0.00000000000000   1.73568248770103   1.73568248770103" << endl;
    oss << "   1.73568248770103   0.00000000000000   1.73568248770103" << endl;
    oss << "   1.73568248770103   1.73568248770103   0.00000000000000" << endl;
    oss << "Mg O" << endl;
    oss << "1 1" << endl;
    oss << "Direct(2) [A1B1]" << endl;
    oss << "   0.00000000000000   0.00000000000000   0.00000000000000  Mg" << endl;
    oss << "   0.50000000000000   0.50000000000000   0.50000000000000  O" << endl;
    oss << endl;
    oss << "(iv) EXAMPLE COMMANDS:" << endl;
    oss << "Assuming that AFLOW is in your PATH and you saved the above example structure file for MgO" << endl; 
    oss << "in the current directory as POSCAR, the following commands can be executed:" << endl;
    oss << endl;
    oss << "aflow --cce=POSCAR --dft_formation_energies=-5.434,-6.220,-6.249 --functionals=PBE,LDA,SCAN" << endl;
    oss << "This will give you the CCE corrections and CCE formation enthalpies for PBE, LDA, and SCAN for MgO." << endl;
    oss << endl;
    oss << "aflow --cce=POSCAR --dft_formation_energies=-6.220 --functionals=LDA" << endl;
    oss << "This gives you only the CCE corrections and CCE formation enthalpies for LDA." << endl;
    oss << endl;
    oss << "aflow --cce=POSCAR --dft_formation_energies=-5.434" << endl;
    oss << "This gives you the CCE corrections and CCE formation enthalpies for PBE with a warning that" << endl; 
    oss << "PBE is assumed as functional." << endl;
    oss << endl;
    oss << "aflow --cce=POSCAR" << endl;
    oss << "This gives you the CCE corrections for PBE, LDA, and SCAN and a rough guess of the formation" << endl;
    oss << "enthalpy based on experimental formation enthalpies per bond." << endl;
    oss << endl;
    oss << "aflow --cce=POSCAR --oxidation_numbers=2,-2" << endl;
    oss << "Oxidation numbers for each atom can also be provided as input." << endl;
    oss << endl;
    return oss.str();
  }
} // namespace cce

#endif // _AFLOW_CCE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow RICO FRIEDRICH - Duke University 2018-2020              *
// *                                                                         *
// ***************************************************************************
