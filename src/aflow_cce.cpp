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
#include "aflowlib.h"
#include "aflow_pflow.h"
#include "aflow_cce.h"

using std::cout;
using std::cerr;
using std::endl;

#define CCE_DEBUG false

namespace cce {

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                           MAIN CCE FUNCTIONS                            //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // for command line use, 
  // use inside AFLOW providing directory path or xstructure & functional string or flags and istream for web tool, 
  // and CCE core function called by all other main CCE functions

  //write_corrections////////////////////////////////////////////////////////
  // main CCE function for command line use 
  // for reading input, analyzing structure, determining oxidation numbers, assigning corrections, 
  // calculating total corrections and corrected formation enthalpies, and writing output
  void write_corrections(aurostd::xoption& flags, ostream& oss) {

    /************************************/
    // Print user instructions
    /************************************/
    // option to print user instructions and exit upon completion
    if(flags.flag("CCE_CORRECTION::USAGE")){
      oss << print_usage();
      return;
    }

    /************************************/
    // Read structure
    /************************************/
    //read structural data from structure file provided on command line
    xstructure structure=read_structure(flags.getattachedscheme("CCE_CORRECTION::POSCAR_PATH"));

    return write_corrections(structure, flags);
  }

  void write_corrections(xstructure& structure, aurostd::xoption& flags) {
    aurostd::xoption cce_flags = init_flags();
    if(!(aurostd::toupper(flags.getattachedscheme("CCE_CORRECTION::PRINT")) == "JSON") && !aurostd::toupper(flags.flag("CCE_CORRECTION::TEST"))){ // additional output like oxidation numbers should not be provided for tests and json
      cce_flags.flag("COMMAND_LINE",TRUE);
    }
    CCE_Variables cce_vars = init_variables(structure);

    /********************************************************/
    // Read DFT formation energies and functionals if provided
    /********************************************************/
    get_dft_form_energies_functionals(flags.getattachedscheme("CCE_CORRECTION::DFT_FORMATION_ENERGIES"), flags.getattachedscheme("CCE_CORRECTION::FUNCTIONALS"), cce_vars); //provide precalc. DFT formation energies & corresponding functionals

    /********************************************************/
    // Read oxidation numbers if provided
    /********************************************************/
    if(flags.flag("CCE_CORRECTION::OXIDATION_NUMBERS")){
      cce_vars.oxidation_states = get_oxidation_states(flags.getattachedscheme("CCE_CORRECTION::OXIDATION_NUMBERS"), structure, cce_flags, cce_vars);
      cce_flags.flag("OX_NUMS_PROVIDED",TRUE);
    }

    return write_corrections(structure, flags, cce_flags, cce_vars);
  }

  void write_corrections(xstructure& structure, aurostd::xoption& flags, aurostd::xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    /********************************************************/
    // Get CCE corrections
    /********************************************************/
    // obtain total CCE corrections per cell from CCE core function only based on
    // structure, oxidation numbers and functional information
    CCE_core(structure, cce_flags, cce_vars);

    // CALCULATE CORRECTED FORMATION ENTHALPIES AT 298.15 AND 0K ###################################
    if (cce_flags.flag("CORRECTABLE")){
      //calculate CCE corrected DFT formation enthalpies if precalculated DFT formation energies are provided
      uint num_funcs=cce_vars.vfunctionals.size();
      for(uint k=0,ksize=num_funcs;k<ksize;k++){ // looping over and checking of vfunctionals is necessary to ensure correct correspondence between given formation energies [k] and corrections with respect to the functional
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
        oss << get_JSON(structure, cce_vars) << std::endl;
      } else if (aurostd::toupper(flags.flag("CCE_CORRECTION::TEST"))) {
        oss << write_test_output(cce_vars, cce_vars.cce_form_energy_cell) << std::endl;
      } else {
        // write CCE corrections & corrected formation enthalpies per cell and atom
        oss << write_output(structure, cce_vars, cce_vars.cce_form_energy_cell);
        // write CCE citation
        oss << write_citation();
      }
    }
  } // main CCE function for command line use



  //write_corrections///////////////////////////////////////////////////////////////////////
  // ME20200213
  // For poscar2cce
  void write_corrections(aurostd::xoption& flags, std::istream& ist) {
    // read structure
    xstructure structure=read_structure(ist);

    return write_corrections(structure, flags);
  }



  //calculate_corrections////////////////////////////////////////////////////////
  // main CCE function for calling inside AFLOW providing only directory path where data neded for correction (POSCAR.static & aflow.in) are located
  // for setting parameters, analyzing structure, determining oxidation numbers, assigning corrections,
  // calculating total corrections, converting correction vector, and returning corrections
  vector<double> calculate_corrections(const string& directory_path) {
    string soliloquy="cce::calculate_corrections():";
    stringstream message;
    // get structure
    string poscar;
    try { // try most relaxed CONTCAR first
      poscar=aflowlib::vaspfile2stringstream(directory_path,"CONTCAR");
    } catch (aurostd::xerror e) { // if that doesn't work try POSCAR
      poscar=aflowlib::vaspfile2stringstream(directory_path,"POSCAR");
    }
    xstructure structure=read_structure(poscar); // AFLOW seems to automatically unzip and rezip zipped files so that only the file name without zipping extension needs to be given
    // get functional from aflow.in in directory
    string aflowin_file= directory_path + "/" + _AFLOWIN_;
    string outcar_file= directory_path + "/" + "OUTCAR.relax1";
    string functional=get_functional_from_aflow_in(structure, aflowin_file, outcar_file);
    if (functional.empty()) {
      message << "Functional cannot be determined from aflow.in. Corrections are available for PBE, LDA, SCAN, or PBE+U_ICSD.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    return calculate_corrections(structure, functional);
  } // main CCE function for calling inside AFLOW with directory path



  //calculate_corrections////////////////////////////////////////////////////////
  // main CCE function for calling inside AFLOW providing only structure (and functional) assuming that oxidation numbers will be obtained from Bader charges or electronegativities
  // for setting parameters, analyzing structure, determining oxidation numbers, assigning corrections,
  // calculating total corrections, converting correction vector, and returning corrections
  //vector<double> CCE(xstructure& structure) // OLD: functional will be automatically determined during Bader charge analysis for the current implementation, later when using e.g. electronegativities, it might be needed as input
  vector<double> calculate_corrections(xstructure& structure, string& functional) { // functional needed as input when determining oxidation numbers from electronegativities
    string soliloquy="cce::calculate_corrections():";
    stringstream message;
    CCE_Variables cce_vars = init_variables(structure);
    aurostd::xoption cce_flags = init_flags();
    // determine functional
    vector<string> CCE_vallowed_functionals;
    aurostd::string2tokens(CCE_allowed_functionals, CCE_vallowed_functionals, ",");
    if(functional!="exp"){
      functional=aurostd::toupper(functional);
    }
    if (!aurostd::withinList(CCE_vallowed_functionals, functional) || get_offset(functional) == -1) {
      message << "Unknown functional " << functional << ". Please choose PBE, LDA, or SCAN.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    cce_vars.vfunctionals.push_back(functional);
    cce_vars.offset.push_back(get_offset(functional));
    // run main CCE function to determine correction
    CCE_core(structure, cce_flags, cce_vars);
    // cce_vars.cce_corrections can be returned directly since there is always only one functional for this CCE function
    return cce_vars.cce_correction;
  } // main CCE function for calling inside AFLOW



  //CCE_core////////////////////////////////////////////////////////
  // main CCE function core called by all other main CCE functions
  // analyzing structure, determining oxidation numbers, assigning corrections,
  // and calculating total corrections
  void CCE_core(xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::CCE_core():";
    // corrections_atom, (su-)perox_correction, cce_correction and cce_form_energy_cell can only be resized after vfunctionals.size() is known from previous main CCE functions calling CCE_core
    // vfunctionals.size()*2 for 298.15 and 0K corrections
    cce_vars.corrections_atom.resize(cce_vars.vfunctionals.size()*2, vector<double>(structure.atoms.size()));
    cce_vars.perox_correction.resize(cce_vars.vfunctionals.size()*2);
    cce_vars.superox_correction.resize(cce_vars.vfunctionals.size()*2);
    cce_vars.cce_correction.resize(cce_vars.vfunctionals.size()*2);
    cce_vars.cce_form_energy_cell.resize(cce_vars.vfunctionals.size()*2);

    // DETERMINE NUMBER OF NEAREST ANION NEIGHBORS FOR EACH CATION: STRUCTURAL PART OF CORRECTION ##
    /********************************************************/
    // Determine anion species
    /********************************************************/
    // determine which species is the anion (for single anion systems)
    cce_vars.anion_species=determine_anion_species(structure, cce_vars);

    /********************************************************/
    // Check for multi anion system
    /********************************************************/
    // check whether it is a multi-anion system and set variables accordingly
    cce_flags.flag("MULTI_ANION_SYSTEM",FALSE);
    cce_flags.flag("O_MULTI_ANION_SPECIES",FALSE); // whether one of the multi anion species is O for which per/superoxide tests need to be made
    cce_vars.multi_anion_atoms=check_for_multi_anion_system(structure, cce_flags, cce_vars);
    // multi anion corrections can only be resized after number of multi anion species is known from check for multi anion system
    cce_vars.multi_anion_corrections_atom.resize(cce_vars.multi_anion_species.size(), vector<vector<double> >(cce_vars.vfunctionals.size()*2, vector<double>(structure.atoms.size())));

    /********************************************************/
    // Determine anion nearest neighbors for each cation
    /********************************************************/
    // apply species selective cutoffs to determine only nearest neighbors within respective cutoff
    cce_vars.num_neighbors=get_num_neighbors(structure, cce_vars.anion_species, cce_flags, cce_vars);

    // determine anion nearest neighbors for cations bound to multi anion atoms if needed
    vector<vector<uint> > multi_anion_num_neighbors(cce_vars.multi_anion_species.size(), vector<uint>(structure.atoms.size()));
    if (cce_flags.flag("MULTI_ANION_SYSTEM")){
      for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
        if(LDEBUG){
          cerr << soliloquy << "getting neighbors for multi anion species " << k << " (" << cce_vars.multi_anion_species[k] << ")" << endl;
        }
        multi_anion_num_neighbors[k]=get_num_neighbors(structure, cce_vars.multi_anion_species[k], cce_flags, cce_vars);
      }
    }

    /********************************************************/
    // Check for per- and superoxide ions
    /********************************************************/
    // check for per- and superoxide ions based on O-O bond length and set number of (su-)peroxide bonds 
    // and indices accordingly if ions are detected
    if(cce_vars.anion_species == "O" || cce_flags.flag("O_MULTI_ANION_SPECIES")) {
      check_per_super_oxides(structure, cce_flags, cce_vars);
    }


    // DETERMINE CORRECTION FOR EACH CATION: OXIDATION STATE DEPENDENT PART OF COORECTION ##########
    /********************************************************/
    // Assign oxidation numbers
    /********************************************************/
    // determine oxidation numbers automatically from structure and Allen electronegativities if not provided on command line
    if(!cce_flags.flag("OX_NUMS_PROVIDED")) {
      if(DEFAULT_CCE_OX_METHOD == "ELECTRONEGATIVITY_ALLEN") {
        cce_vars.oxidation_states=get_oxidation_states_from_electronegativities(structure, cce_flags, cce_vars);
      } else if(DEFAULT_CCE_OX_METHOD == "BADER") { // obtaining oxidation states from Bader charges is outdated but the functionality is kept mainly for test purposes
        cce_vars.oxidation_states=get_oxidation_states_from_Bader(structure, cce_flags, cce_vars);
      }
    }

    /********************************************************/
    // Assign corrections
    /********************************************************/
    // assigning corrections for each atom after oxidation numbers are determined if the system is correctable
    if (cce_flags.flag("CORRECTABLE")){
      get_corrections(structure, cce_flags, cce_vars, cce_vars.anion_species, cce_vars.num_neighbors, cce_vars.corrections_atom);
      // determine corrections for cations bound to multi anion atoms if needed
      if (cce_flags.flag("MULTI_ANION_SYSTEM")){
        for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
          if(LDEBUG){
            cerr << endl;
            cerr << soliloquy << "getting corrections for multi anion species " << k << " (" << cce_vars.multi_anion_species[k] << ")" << endl;
          }
          get_corrections(structure, cce_flags, cce_vars, cce_vars.multi_anion_species[k], multi_anion_num_neighbors[k], cce_vars.multi_anion_corrections_atom[k]);
        }
      }
      // load per- & superox. corrections if needed;
      if (cce_vars.num_perox_bonds > 0 || cce_vars.num_superox_bonds > 0){
        check_get_per_super_ox_corrections(cce_vars);
      }
    }


    // CALCULATE CCE CORRECTIONS AT 298.15 AND 0K ##################################################
    // calculate CCE correction if the system is correctable
    if (cce_flags.flag("CORRECTABLE")){
      // calculate total correction per cell only using cations although 
      // corrections for anions should be set to zero
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){
        if (structure.atoms[i].cleanname != cce_vars.anion_species && cce_vars.multi_anion_atoms[i] != 1){ // exclude main anion species and multi anion atoms detected previously
          if (cce_vars.num_neighbors[i] > 0){ // are there actually bonds between the cation and the (main) anion species
            uint num_funcs=cce_vars.vfunctionals.size();
            for (uint k = 0; k < num_funcs; k++) {
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
            cerr << endl;
            cerr << soliloquy << "adding corrections for multi anion species " << k << " (" << cce_vars.multi_anion_species[k] << ")" << endl;
          }
          for(uint i=0,isize=structure.atoms.size();i<isize;i++){
            if (structure.atoms[i].cleanname != cce_vars.anion_species && cce_vars.multi_anion_atoms[i] != 1){ // exclude main anion species and multi anion atoms detected previously
              if (multi_anion_num_neighbors[k][i] > 0){ // are there actually bonds between the cation and the anion species under consideration
                uint num_funcs=cce_vars.vfunctionals.size();
                for (uint l = 0; l < num_funcs; l++) {
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
        check_apply_per_super_ox_corrections(cce_vars);
      }
    }
  } // main CCE function core



  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                   READ USER INPUT (FROM COMMAND LINE)                   //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  //read_structure////////////////////////////////////////////////////////
  // read structural data from structure file provided on command line
  xstructure read_structure(const string& structure_file, int mode){ // first argument can be directly structure_file and not structure_file_path
    string soliloquy="cce::read_structure():";
    stringstream message;
    //string structure_file = aurostd::file2string(structure_file_path); // first argument of read_structure_function does not need to be converted to string since it contains already the file content and not only the file name
    xstructure structure(structure_file, mode);
    return read_structure(structure);
  }

  //read_structure////////////////////////////////////////////////////////
  // read structural data from istream
  xstructure read_structure(std::istream& ist){
    string soliloquy="cce::read_structure():";
    xstructure structure(ist);
    return read_structure(structure);
  }

  //read_structure////////////////////////////////////////////////////////
  xstructure read_structure(xstructure& structure){
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::read_structure():";
    stringstream message;
    structure.ReScale(1.0); // rescales scaling factor in second line of POSCAR to 1, needed for correct distances
    //let the program spit out what it thinks (input structure)
    if(LDEBUG){
      cerr << soliloquy << endl << "INPUT STRUCTURE:" << endl;
      cerr << soliloquy << structure << endl;
    }
    // check whether there are any atoms in the structure
    if (structure.atoms.size() == 0){
      message << "BAD NEWS: It seems there are no atoms in the structure file. Please adjust the structure file and rerun.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
    }
    // if species of atoms are not known like in VASP4 format, throw error
    for(uint k=0,ksize=structure.atoms.size();k<ksize;k++){
      if (structure.atoms[k].cleanname == ""){
        message << "BAD NEWS: It seems you are providing a POSCAR without complete species information as input. This implementation requires a structure in VASP POSCAR format with the species information included. Please adjust the structure file and rerun.";
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

  //get_dft_form_energies_functionals////////////////////////////////////////////////////////
  // set the DFT formation energies and functionals according to the input and check consistency of functionals and formation energies
  void get_dft_form_energies_functionals(const string& dft_energies_input_str, const string& functionals_input_str, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::get_dft_form_energies_functionals():";
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
      cerr << soliloquy << "INPUT DFT FORMATION ENERGIES & FUNCTIONALS:" << endl;
      cerr << soliloquy << " input dft formation energies=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(cce_vars.dft_energies,6,roff),",") << " (assumed to be in eV/cell)" << endl;
      cerr << soliloquy << " input functionals=" << functionals_input_str << endl;
    }
    // determine whether it is a functional for which corrections are available
    vector<string> CCE_vallowed_functionals;
    aurostd::string2tokens(CCE_allowed_functionals, CCE_vallowed_functionals, ",");
    uint num_funcs=cce_vars.vfunctionals.size();
    for(uint k=0,ksize=num_funcs;k<ksize;k++){
      if(cce_vars.vfunctionals[k]!="exp"){
        cce_vars.vfunctionals[k]=aurostd::toupper(cce_vars.vfunctionals[k]);
      }
      if(LDEBUG){
        cerr << soliloquy << "cce_vars.vfunctionals[" << k << "]: " << cce_vars.vfunctionals[k] << endl;
      }
      if (!aurostd::withinList(CCE_vallowed_functionals, cce_vars.vfunctionals[k]) || get_offset(cce_vars.vfunctionals[k]) == -1) {
        message << "Unknown functional " << cce_vars.vfunctionals[k] << ". Please choose PBE, LDA, or SCAN.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
      }
      cce_vars.offset.push_back(get_offset(cce_vars.vfunctionals[k]));
    }
    // if only structure (and oxidation numbers) are provided, corrections should be given for PBE, LDA, and SCAN 
    // and the CCE@exp formation enthalpy from the exp. formation enthalpies per bond
    // ICSD correction should only be returned when explicitly asked for
    if(functionals_input_str.empty() && dft_energies_input_str.empty()){
      vector<string> CCE_vdefault_output_functionals;
      aurostd::string2tokens(CCE_default_output_functionals, CCE_vdefault_output_functionals, ",");
      for(uint k=0,ksize=CCE_vdefault_output_functionals.size();k<ksize;k++){
        if (get_offset(CCE_vdefault_output_functionals[k]) == -1) {
          message << "Unknown functional " << cce_vars.vfunctionals[k] << ". Please choose PBE, LDA, or SCAN.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
        }
        cce_vars.vfunctionals.push_back(CCE_vdefault_output_functionals[k]); cce_vars.offset.push_back(get_offset(CCE_vdefault_output_functionals[k]));
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "PBE: " << aurostd::withinList(cce_vars.vfunctionals, "PBE") << endl;
      cerr << soliloquy << "LDA: " << aurostd::withinList(cce_vars.vfunctionals, "LDA") << endl;
      cerr << soliloquy << "SCAN: " << aurostd::withinList(cce_vars.vfunctionals, "SCAN") << endl;
      cerr << soliloquy << "PBE+U_ICSD: " << aurostd::withinList(cce_vars.vfunctionals, "PBE+U_ICSD") << endl;
      cerr << soliloquy << "exp: " << aurostd::withinList(cce_vars.vfunctionals, "exp") << endl;
      cerr << endl;
    }
  }

  //get_offset////////////////////////////////////////////////////////
  // get offset needed for reading corrections from lookup table for different functionals
  int get_offset(const string& functional) {
  if (functional=="PBE")          {return 0;}
  if (functional=="LDA")          {return 2;}
  if (functional=="SCAN")         {return 4;}
  if (functional=="PBE+U_ICSD")   {return 6;}
  if (functional=="exp")          {return 8;}
  else {return -1;}
  }

  //get_oxidation_states////////////////////////////////////////////////////////
  // Retrieves the oxidation states of the material.
  vector<double> get_oxidation_states(const string& oxidation_numbers_input_str, const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    string soliloquy="cce::get_oxidation_states():";
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
        oss << endl;
        oss << "INPUT OXIDATION NUMBERS:" << endl;
        for (uint k=0,ksize=cce_vars.oxidation_states.size();k<ksize;k++){
          oss << "Oxidation state of " << structure.atoms[k].cleanname << " (atom[" << k << "]): " << cce_vars.oxidation_states[k] << endl;
        }
      }
      // get sum of oxidation numbers and validate (system should not be regarded correctable if sum over oxidation states is not zero)
      cce_vars.oxidation_sum = get_oxidation_states_sum(cce_flags, cce_vars); // double because for superoxides O ox. number is -0.5
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
      }
      if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
        string function = "cce::get_oxidation_states()";
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

  //get_functional_from_aflow_in////////////////////////////////////////////////////////
  // determine the functional from the aflow_in
  string get_functional_from_aflow_in(const xstructure& structure, string& aflowin_file, string& outcar_file) {
    string soliloquy="cce::get_functional_from_aflow_in():";
    stringstream message;
    string functional = "";
    string aflowIn = aurostd::RemoveComments(aurostd::file2string(aflowin_file));
    vector<string> vlines = aurostd::string2vectorstring(aflowIn);
    string line_a = "";
    bool pbe = false;
    bool ldau = false;
    bool pbe_u_icsd = false;
    bool ldau2 = false;
    bool lda = false;
    bool scan = false;
    _kflags kflags;
    _aflags aflags;
    _vflags vflags=KBIN::VASP_Get_Vflags_from_AflowIN(aflowIn,aflags,kflags);
    //_xvasp xvasp;xvasp.clear();
    //KBIN::readModulesFromAflowIn(aflowIn, kflags, xvasp);  // ME20181027
    //cout << "xvasp.str.species_pp: " << xvasp.str.species_pp << endl;
    //cout << "xvasp.str: " << xvasp.str << endl;
    if (aurostd::FileExist(outcar_file) || aurostd::EFileExist(outcar_file)) {
      xOUTCAR outcar;
      outcar.GetPropertiesFile(outcar_file);
      uint pbe_count=0;
      uint lda_count=0;
      for(uint i=0;i<outcar.species_pp_type.size();i++) {
        if (outcar.species_pp_type.at(i).find("PAW_PBE") != string::npos){
          pbe_count+=1;
        } else if (outcar.species_pp_type.at(i).find("PAW_LDA") != string::npos){
          lda_count+=1;
        }
      }
      if (outcar.species_pp_type.size() == pbe_count){
        pbe = true;
      } else if (outcar.species_pp_type.size() == lda_count){
        lda = true;
      }
    } else {
      for (uint i = 0; i < vlines.size(); i++) {
        line_a = aurostd::RemoveSpaces(vlines[i]);
        // check whether it is a PBE calculation; PBE+U will be checked later
        if (line_a.find("=potpaw_PBE") != string::npos || line_a.find("potpaw_PBE/") != string::npos){
          pbe = true;
        }
        // check whether it is an LDA
        if (line_a.find("=potpaw_LDA") != string::npos || line_a.find("potpaw_LDA/") != string::npos){
          lda = true;
        }
      }
    }
    //check for SCAN
    if (vflags.KBIN_VASP_FORCE_OPTION_METAGGA.xscheme == "SCAN"){
      scan = true;
    }
    // check whether it is a DFT+U calculation with parameters as for the ICSD (PBE+U_ICSD calculation)
    // new implementation checking Us explicitly
    if (vflags.KBIN_VASP_LDAU_PARAMETERS != "" && vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry && !vflags.KBIN_VASP_FORCE_OPTION_LDAU_ADIABATIC.isentry && !vflags.KBIN_VASP_FORCE_OPTION_LDAU_CUTOFF.isentry){
      ldau2=true;
      vector<string> ldau_params_vector;
      aurostd::string2tokens(vflags.KBIN_VASP_LDAU_PARAMETERS, ldau_params_vector, ";");
      // get species
      string species_part = ldau_params_vector[0];
      vector<string> species_vector;
      aurostd::string2tokens(species_part, species_vector, ",");
      if (species_vector.size() != structure.species.size()){
        message << "BAD NEWS: The number of species in the DFT+U settings differs from the total number of species for this structure. Please adapt and rerun.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
      }
      // get Us
      string Us_part = ldau_params_vector[2];
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
        pbe_u_icsd = true;
      }
    }
    // check whether it is a DFT+U calculation with different parameters than for PBE+U_ICSD
    if ((vflags.KBIN_VASP_FORCE_OPTION_LDAU1.isentry || vflags.KBIN_VASP_FORCE_OPTION_LDAU2.isentry) && !pbe_u_icsd){
      ldau = true;
      message << "BAD NEWS: It seems you are providing an aflow.in for a DFT+U calculation with different parameters than for the AFLOW ICSD database (Dudarev's approach, LDAU2=ON). There are no corrections for this case.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
    }

    if (pbe){
      if (pbe_u_icsd){
        functional = "PBE+U_ICSD";
      } else if (scan){ //SCAN calculations are usually done with PBE PPs
        functional = "SCAN";
      } else if (!ldau && !ldau2){
        functional = "PBE";
      }
    } else if (lda){
      if (scan){
        functional = "SCAN"; //SCAN calcs. with LDA PPs should be okay
      } else if (!ldau && !ldau2){
        functional = "LDA";
      }
    } else if (scan){
      if (!ldau && !ldau2){
        functional = "SCAN";
      }
    }
    //cout << "functional: " << functional << endl;
    return functional;
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                      INITIALISE FLAGS AND VARIABLES                     //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////

  //init_flags////////////////////////////////////////////////////////////
  // ME20200213
  // Initializes the CCE flags to their default values.
  aurostd::xoption init_flags() {
    aurostd::xoption flags;
    flags.flag("COMMAND_LINE", false);
    flags.flag("CORRECTABLE", true); // first assuming that formation energy of system IS correctable; will be set to not correctable if, for any atom, no correction can be identified
    flags.flag("OX_NUMS_PROVIDED", false);
    flags.flag("MULTI_ANION_SYSTEM", false);
    flags.flag("O_MULTI_ANION_SPECIES", false); // whether one of the multi anion species is O for which per/superoxide tests need to be made
    return flags;
  }

  //init_variables////////////////////////////////////////////////////////
  // ME20200213
  // Initializes the variables struct.
  CCE_Variables init_variables(const xstructure& structure) {
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
    cce_vars.species_electronegativity_sorted.clear();
    cce_vars.species_electronegativity_sorted.resize(structure.species.size());
    cce_vars.num_pref_ox_states_electronegativity_sorted.clear();
    cce_vars.num_pref_ox_states_electronegativity_sorted.resize(structure.species.size());
    cce_vars.num_all_ox_states_electronegativity_sorted.clear();
    cce_vars.num_all_ox_states_electronegativity_sorted.resize(structure.species.size());
    cce_vars.pref_ox_states_electronegativity_sorted.clear();
    cce_vars.pref_ox_states_electronegativity_sorted.resize(structure.species.size());
    cce_vars.all_ox_states_electronegativity_sorted.clear();
    cce_vars.all_ox_states_electronegativity_sorted.resize(structure.species.size());
    cce_vars.cations_map.clear();
    cce_vars.Bader_charges.clear();
    cce_vars.corrections_atom.clear();
    cce_vars.multi_anion_corrections_atom.clear();
    cce_vars.perox_correction.clear();
    cce_vars.superox_correction.clear();
    cce_vars.cce_correction.clear();
    cce_vars.cce_form_energy_cell.clear();
    return cce_vars;
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                           STRUCTURAL ANALYSIS                           //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // counting the number of anion neighbours (multiple anion species possible) for each cation

  //determine_anion_species////////////////////////////////////////////////////////
  // determine anion species
  string determine_anion_species(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::determine_anion_species():";
    stringstream message;
    if(LDEBUG){
      cerr << soliloquy << "ANION SPECIES FROM ALLEN ELECTRONEGATIVITIES:" << endl;
    }
    _atom atom;
    uint z = 0;

    double anion_electronegativity = 0;
    for(uint k=0,ksize=structure.species.size();k<ksize;k++){
      z = GetAtomNumber(KBIN::VASP_PseudoPotential_CleanName(structure.species[k]));
      xelement element(z);
      if (element.electronegativity_Allen == NNN) {
        message << "VERY BAD NEWS: There is no known electronegativity value for " << KBIN::VASP_PseudoPotential_CleanName(structure.species[k]) << ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
      } else{
        cce_vars.electronegativities[k] = element.electronegativity_Allen;
        if(LDEBUG){
          cerr << soliloquy << "electronegativity of species " << k << " (" << KBIN::VASP_PseudoPotential_CleanName(structure.species[k]) << "): " << cce_vars.electronegativities[k] << endl;
        }
        if (cce_vars.electronegativities[k] > anion_electronegativity) {
          anion_electronegativity = cce_vars.electronegativities[k];
          cce_vars.anion_species = KBIN::VASP_PseudoPotential_CleanName(structure.species[k]);
        }
      }
    }
    // set anion charge and check whether it is negative
    z = GetAtomNumber(cce_vars.anion_species);
    xelement element(z);
    cce_vars.standard_anion_charge = element.oxidation_states[element.oxidation_states.size()-1];
    if (cce_vars.standard_anion_charge > 0) {
      message << "VERY BAD NEWS: There is no known negative oxidation number for " << cce_vars.anion_species << " detected as anion species.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
    }
    if(LDEBUG){
      cerr << soliloquy << "anion electronegativity: " << anion_electronegativity << endl;
      cerr << soliloquy << "anion species: " << cce_vars.anion_species << endl;
      cerr << soliloquy << "anion charge: " << cce_vars.standard_anion_charge << endl;
      cerr << endl;
    }
    return cce_vars.anion_species;
  }

  //check_for_multi_anion_system////////////////////////////////////////////////////////
  // check whether it is a multi-anion system, i. e. whether atoms of another species than the the anion_species are only bound to atoms of lower electronegativity or of the same type
  vector<uint> check_for_multi_anion_system(xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, double tolerance, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::check_for_multi_anion_system():";
    stringstream message;
    if(LDEBUG){
      cerr << soliloquy << "CHECKING FOR MULTI-ANION SYSTEM:" << endl;
    }
    for ( uint i = 0; i < structure.atoms.size(); i++ ) { // initialize elements of vector to 0 
      cce_vars.multi_anion_atoms[i] = 0; 
    }
    cce_vars.cutoffs=get_dist_cutoffs(structure, tolerance);
    double cutoffs_max=aurostd::max(cce_vars.cutoffs);
    deque<deque<_atom> > neigh_mat;
    structure.GetStrNeighData(cutoffs_max,neigh_mat);
    uint z = 0;
    for(uint i=0,isize=neigh_mat.size();i<isize;i++){ //same size as structure.atoms.size(); number of atoms in the structure (not determined by cutoff (or cutoffs_max))
      uint neighbors_count=0;
      uint multi_anion_count=0;
      for(uint j=0,jsize=neigh_mat[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by cutoffs_max
        const _atom& atom=neigh_mat[i][j]; // the atom object stands for the neighbors of each atom of the structure
        if (_CCE_SELF_DIST_TOL_ < AtomDist(structure.atoms[i],atom) && AtomDist(structure.atoms[i],atom) <= cce_vars.cutoffs[structure.atoms[i].type] ){ // distance must be larger than _CCE_SELF_DIST_TOL_ since GetStrNeighData includes also structure.atoms[i] itself as neighbor having distance zero to itself
          if (atom.cleanname == cce_vars.anion_species){
            neighbors_count+=1;
          } else if (atom.cleanname != cce_vars.anion_species && structure.atoms[i].cleanname != cce_vars.anion_species){ // second condition set since the anion_species cannot be set as a multi-anion species again
            neighbors_count+=1;
            z = GetAtomNumber(structure.atoms[i].cleanname);
            xelement atom_element(z);
            z = GetAtomNumber(atom.cleanname);
            xelement neigh_element(z);
            double electronegativity_atom = atom_element.electronegativity_Allen;
            double electronegativity_neighbor = neigh_element.electronegativity_Allen;
            if(LDEBUG){
              cerr << soliloquy << "electronegativity of atom " << i << ": " << electronegativity_atom << endl;
              cerr << soliloquy << "electronegativity of neighbor " << j << ": " << electronegativity_neighbor << endl;
            }
            if(electronegativity_neighbor < electronegativity_atom || structure.atoms[i].cleanname == atom.cleanname){ // could be multi-anion atom if bound to only atoms of lower electroneg. or of same species
              multi_anion_count+=1;
            }
          }
        }
      }
      if(LDEBUG){
        cerr << soliloquy << "multi_anion_count for atom " << i << " (" << structure.atoms[i].cleanname << "): " << multi_anion_count << endl;
        cerr << soliloquy << "neighbors_count for atom " << i << " (" << structure.atoms[i].cleanname << "): " << neighbors_count << endl;
      }
      if (multi_anion_count == neighbors_count && structure.atoms[i].cleanname != cce_vars.anion_species){ // anion_species should not be detected again as multi_anion_species
        if (!cce_flags.flag("MULTI_ANION_SYSTEM")){
          cce_flags.flag("MULTI_ANION_SYSTEM",TRUE);
          if(cce_flags.flag("COMMAND_LINE")){
            oss << endl;
            oss << "This system has been detected to be a multi-anion compound!" << endl;
          }
        }
        // set multi anion atoms
        cce_vars.multi_anion_atoms[i]=1;
        if(LDEBUG){
          cerr << soliloquy << "Atom " << i << " (" << structure.atoms[i].cleanname << ") has been detected as a multi-anion atom." << endl;
        }
        // set multi anion species
        uint multi_anion_species_count=0;
        for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ // have to make sure that no anion species is added twice
          if(cce_vars.multi_anion_species[k] != structure.atoms[i].cleanname){
            multi_anion_species_count+=1;
          }
        }
        if(multi_anion_species_count == cce_vars.multi_anion_species.size()){
          cce_vars.multi_anion_species.push_back(structure.atoms[i].cleanname);
          if(LDEBUG){
            cerr << soliloquy << "New multi anion species: " << structure.atoms[i].cleanname << endl;
          }
          if(structure.atoms[i].cleanname == "O"){ // check whether one of the multi anion species is O for which per/superoxide test needs to be made
            cce_flags.flag("O_MULTI_ANION_SPECIES",TRUE);
            if(LDEBUG){
              cerr << soliloquy << "Oxygen was detected as multi anion species, i.e. system needs to be tested for (su-)peroxide ions." << endl;
            }
          }
        }
        // set multi anion oxidation numbers and check whether it is negative
        _atom atom;
        z = GetAtomNumber(structure.atoms[i].cleanname);
        xelement atom_element(z);
        cce_vars.oxidation_states[i] = atom_element.oxidation_states[atom_element.oxidation_states.size()-1];
        if(LDEBUG){
          cerr << soliloquy << "Oxidation state for atom " << i << " (" << structure.atoms[i].cleanname << ") has been set to: " << cce_vars.oxidation_states[i] << endl;
        }
        if (cce_vars.oxidation_states[i] > 0) {
          message << "VERY BAD NEWS: There is no known negative oxidation number for " << structure.atoms[i].cleanname << " detected as multi anion species.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
        }
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of multi anion species (1-total_number_anion_species since the main anion species with largest electronegativity is NOT counted as multi anion species): " << cce_vars.multi_anion_species.size() << endl;
      for(uint k=0,ksize=cce_vars.multi_anion_species.size();k<ksize;k++){ 
        cerr << soliloquy << "multi anion species " << k << ": " << cce_vars.multi_anion_species[k] << endl;
      }
      cerr << endl;
    }
    return cce_vars.multi_anion_atoms;
  }

  //get_num_neighbors////////////////////////////////////////////////////////
  // trying to get the number of neighbors for each atom within respective species selective cutoff
  vector<uint> get_num_neighbors(xstructure& structure, double tolerance) {
    CCE_Variables cce_vars = init_variables(structure);
    string anion_species="";
    aurostd::xoption cce_flags = init_flags();
    return get_num_neighbors(structure, anion_species, cce_flags, cce_vars, tolerance);
  }

  //get_num_neighbors////////////////////////////////////////////////////////
  // trying to get the number of neighbors of a special anion_species for each atom within respective species selective cutoff
  vector<uint> get_num_neighbors(xstructure& structure, const string& anion_species, double tolerance) {
    CCE_Variables cce_vars = init_variables(structure);
    aurostd::xoption cce_flags = init_flags();
    return get_num_neighbors(structure, anion_species, cce_flags, cce_vars, tolerance);
  }

  //get_num_neighbors////////////////////////////////////////////////////////
  // trying to get the number of neighbors for each atom within respective species selective cutoff
  vector<uint> get_num_neighbors(xstructure& structure, const string& anion_species, xoption& cce_flags, CCE_Variables& cce_vars, double tolerance, ostream& oss) { // anion_species here cannot be taken from cce_vars since function is also used to determine multi anion num_neighbors for which anion_species is the respective multi_anion_species
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::get_num_neighbors():";
    if(LDEBUG){
      cerr << soliloquy << "STRUCTURAL ANALYSIS:" << endl;
    }
    cce_vars.cutoffs=get_dist_cutoffs(structure, tolerance);
    double cutoffs_max=aurostd::max(cce_vars.cutoffs);
    vector<uint> num_neighbors(structure.atoms.size());
    deque<deque<_atom> > neigh_mat;
    structure.GetStrNeighData(cutoffs_max,neigh_mat);
    for(uint i=0,isize=neigh_mat.size();i<isize;i++){ //same size as structure.atoms.size(); number of atoms in the structure (not determined by cutoff (or cutoffs_max))
      uint neighbors_count=0;
      bool warning = false;
      bool insert_empty_line=true;
      for(uint j=0,jsize=neigh_mat[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by cutoffs_max
        const _atom& atom=neigh_mat[i][j];
        if (_CCE_SELF_DIST_TOL_ < AtomDist(structure.atoms[i],atom) && AtomDist(structure.atoms[i],atom) <= cce_vars.cutoffs[structure.atoms[i].type] ){ // distance must be larger than _CCE_SELF_DIST_TOL_ since GetStrNeighData includes also structure.atoms[i] itself as neighbor having distance zero to itself
          if (!anion_species.empty()){ // variable called anion type since function was developed for CCE for polar materials but it can be used to check for any atom type and only include those as neighbors
            // implement check whether each nearest neighbor is of the anion_species, otherwise throw warning; 
            if (atom.cleanname == anion_species){
              neighbors_count+=1;
            } else if (atom.cleanname != anion_species && structure.atoms[i].cleanname != anion_species){ // second condition set since it is expected that the anion has predominantly other neighbors than its own type 
              if (!cce_flags.flag("MULTI_ANION_SYSTEM")){
                if (insert_empty_line && cce_flags.flag("COMMAND_LINE")) { oss << endl; insert_empty_line = false; } // construction just to make sure that only one empty line is added at the beginning of the warning block
                warning = true;
                if(cce_flags.flag("COMMAND_LINE")){
                  oss << "WARNING: Not all nearest neighbors of " << structure.atoms[i].cleanname << " (ATOM[" << i << "]) within the distance tolerance of " << tolerance << " Ang. are " << anion_species << ", there is also " << atom.cleanname << endl;
                }
              }
            }
          } else {
            neighbors_count+=1;
          }
        }
      }
      if (!cce_flags.flag("MULTI_ANION_SYSTEM") && warning && cce_flags.flag("COMMAND_LINE")){
        oss << "WARNING: Not all nearest neighbors of " << structure.atoms[i].cleanname << " (ATOM[" << i << "]) within the distance tolerance are " << anion_species << "!" << endl;
      }
      num_neighbors[i]=neighbors_count; // zero-based counting as for cutoffs array above
      if(LDEBUG){
        cerr << soliloquy << "number of " << anion_species << " nearest neighbors within " << tolerance << " Ang. tolerance of " << structure.atoms[i].cleanname << " (ATOM[" << i << "]): " << num_neighbors[i] << endl;
        cerr << endl;
      }
    }
    return num_neighbors;
  }

  //get_dist_cutoffs////////////////////////////////////////////////////////
  // determine species selective nearest neighbor distances and then cutoffs accordingly
  vector<double> get_dist_cutoffs(const xstructure& structure, double tolerance) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::get_dist_cutoffs():";
    vector<double> cutoffs(structure.species.size());
    xmatrix<double> distsij=GetDistMatrix(structure); // gives matrix with nearest neighbor distances between all species pairs with species running over rows and columns
    vector<double> near_neigh_dist(structure.species.size());
    for(int i=1,isize=distsij.rows+1;i<isize;i++){
      for(int j=1,jsize=distsij.cols+1;j<jsize;j++){
        if (j==1){
          near_neigh_dist[i-1]=distsij(i,j);
        }
        if (distsij(i,j) < near_neigh_dist[i-1]){
          near_neigh_dist[i-1]=distsij(i,j);
        }
      }
      if(LDEBUG){
        cerr << soliloquy << "nearest neighbor distance for species " << i << " is " << near_neigh_dist[i-1] << " Ang." << endl;
        cerr << soliloquy << "cutoff for the distance for species " << i << " is " << near_neigh_dist[i-1]+tolerance << " Ang." << endl;
      }
      cutoffs[i-1]=near_neigh_dist[i-1]+tolerance; // -1 since counting of array elements starts from zero, NOT 1
    }
    return cutoffs;
  }

  //check_per_super_oxides////////////////////////////////////////////////////////
  // check whether the system contains per- or superoxide ions based on the O-O bond length and set variables accordingly
  void check_per_super_oxides(xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::check_per_super_oxides():";
    stringstream message;
    if(LDEBUG){
      cerr << soliloquy << "CHECKING FOR (SU-)PEROXIDES:" << endl;
    }
    uint perox_count=0;
    uint superox_count=0;
    for ( uint i = 0; i < structure.atoms.size(); i++ ) { // initialize elements of vectors to 0 
      cce_vars.perox_indices[i] = 0; 
      cce_vars.superox_indices[i] = 0; 
    }
    cce_vars.cutoffs=get_dist_cutoffs(structure);
    double cutoffs_max=aurostd::max(cce_vars.cutoffs);
    deque<deque<_atom> > neigh_mat;
    structure.GetStrNeighData(cutoffs_max,neigh_mat);
    for(uint i=0,isize=neigh_mat.size();i<isize;i++){ //same size as structure.atoms.size(); number of atoms in the structure (not determined by cutoff (or cutoffs_max))
      if (structure.atoms[i].cleanname == "O"){ // identify per- and superoxides by O-O bond length
        for(uint j=0,jsize=neigh_mat[i].size();j<jsize;j++){  //number of nearest neighbors within cutoff of atom i; number of neighbors of each atom i determined by the cutoffs_max
          const _atom& atom=neigh_mat[i][j];
          if (atom.cleanname == "O"){
            if (_CCE_SELF_DIST_TOL_ < AtomDist(structure.atoms[i],atom) && AtomDist(structure.atoms[i],atom) <= _CCE_O2_molecule_cutoff_ ){ // distance must be larger than _CCE_SELF_DIST_TOL_ to savely exclude the anion itself having distance zero to itself; if O-O bond is shorter than in O2 molecule (approx. 1.21 Ang) the result of the structural relaxation is most likely wrong
              message << "THE DETERMINED OXYGEN-OXYGEN BOND LENGTH IS SHORTER THAN IN THE O2 MOLECULE; CHECK YOUR STRUCTURE! THE O-O BOND LENGTH IS: " << AtomDist(structure.atoms[i],atom) << " Ang.";
              throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
            } else if (_CCE_O2_molecule_cutoff_ < AtomDist(structure.atoms[i],atom) && AtomDist(structure.atoms[i],atom) <= _CCE_superox_cutoff_ ){
              if(LDEBUG){
                cerr << soliloquy << "WARNING: This should be a superoxide; the O-O bond length is: " << AtomDist(structure.atoms[i],atom) << " Ang." << endl;
              }
              superox_count+=1;
              cce_vars.superox_indices[i]=1;
            } else if (_CCE_superox_cutoff_ < AtomDist(structure.atoms[i],atom) && AtomDist(structure.atoms[i],atom) <= _CCE_perox_cutoff_ ){
              if(LDEBUG){
                cerr << soliloquy << "WARNING: This should be a peroxide; the O-O bond length is: " << AtomDist(structure.atoms[i],atom) << " Ang." << endl;
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
        oss << endl;
        oss << "WARNING: This should be a peroxide!" << endl;
        oss << "Number of peroxide O-O bonds in cell: " << cce_vars.num_perox_bonds << endl;
      }
    } else if (cce_vars.num_superox_bonds > 0) {
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
        oss << "WARNING: This should be a superoxide!" << endl;
        oss << "Number of superoxide O-O bonds in cell: " << cce_vars.num_superox_bonds << endl;
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //           DETERMINE OXIDATION STATES FROM ELECTRONEGATIVITIES           //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  // determine oxidation numbers based on Allen electronegativities and treat mixed valence binary systems specially

  //get_oxidation_states_from_electronegativities////////////////////////////////////////////////////////
  // function overloading for below function to be able to use oxidation number determination independently of CCE
  vector<double> get_oxidation_states_from_electronegativities(xstructure& structure) {
    CCE_Variables cce_vars = init_variables(structure);
    cce_vars.anion_species=determine_anion_species(structure, cce_vars);
    aurostd::xoption cce_flags = init_flags();
    cce_vars.multi_anion_atoms=check_for_multi_anion_system(structure, cce_flags, cce_vars);
    // multi anion corrections can only be resized after number of multi anion species is known from check for multi anion system
    cce_vars.multi_anion_corrections_atom.resize(cce_vars.multi_anion_species.size(), vector<vector<double> >(cce_vars.vfunctionals.size()*2, vector<double>(structure.atoms.size())));
    cce_vars.num_neighbors=get_num_neighbors(structure, cce_vars.anion_species, cce_flags, cce_vars);
    if(cce_vars.anion_species == "O" || cce_flags.flag("O_MULTI_ANION_SPECIES")) {
      check_per_super_oxides(structure, cce_flags, cce_vars);
    }
    cce_flags.flag("CORRECTABLE",TRUE);
    return get_oxidation_states_from_electronegativities(structure, cce_flags, cce_vars);
  }

  //get_oxidation_states_from_electronegativities////////////////////////////////////////////////////////
  // determine the oxidation numbers of the ions using preferred/all known oxidation numbers, electronegativities and structural information
  vector<double> get_oxidation_states_from_electronegativities(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::get_oxidation_states_from_electronegativities():";
    if(LDEBUG){
      cerr << soliloquy << "DETERMINATION OF OXIDATION NUMBERS FROM PREFERRED/ALL KNOWN OXIDATION STATES, ALLEN ELECTRONEGATIVITIES, AND STRUCTURE:" << endl;
    }
    // deal with anion charge in cell first
    set_anion_oxidation_states(structure, cce_vars);
    // treat cations
    if(LDEBUG){
      cerr << soliloquy << "CATION PART:" << endl;
    }
    // determine number of cation species (currently just all species minus one anion species, multi-anion atoms are dealt with separately)
    uint num_cation_species = structure.species.size()-1;
    // sort species ascending by electronegativity
    // preferred and all known oxidation states will be automatically sorted by electronegativity during subsequent loading
    // the anion species should always be the last one with the highest electronegativity
    sort_species_by_electronegativity(structure, cce_vars);
    // load needed info about preferred and all known oxidation states for all species (sorted by electronegativity)
    cce_flags.flag("NO_PREF_OX_STATES",FALSE);
    cce_flags.flag("NO_OX_STATES",FALSE);
    cce_flags.flag("OX_STATES_DETERMINED",FALSE);
    load_ox_states_templates_each_species(structure, cce_flags, cce_vars);
    // ME Nov. 2019 for getting cations_map: a vector of vectors that lists for each cation species the atom numbers of the structure that are of this species (for Fe2ZnO4 there might be two Fe atoms at positions 0 and 1 in the structure)
    uint natoms = structure.atoms.size();
    cce_vars.cations_map.resize(num_cation_species);
    uint i = 0;
    for (uint at = 0; at < natoms; at++) {
      for (i = 0; i < num_cation_species; i++) {
        if (structure.atoms[at].cleanname == cce_vars.species_electronegativity_sorted[i]) break; // if it finds that the atom belongs to the ith species sorted by electronegativities, break to insert it at the proper place for the ith species into cation_map
      }
      if (i < num_cation_species) cce_vars.cations_map[i].push_back(at);
    }
    // try to find proper oxidation states (making oxidation sum zero) by using preferred oxidation states
    if(!cce_flags.flag("NO_PREF_OX_STATES")){ // for He, Ne, Ar, and Xe there are no preferred ox. states
      try_preferred_oxidation_states(structure, cce_flags, cce_vars);
    }
    // if preferred oxidation states do not work, it could be a mixed valence (binary) system
    // treat them here as special cases
    if(!cce_flags.flag("OX_STATES_DETERMINED")){
      // for SbO2 the oxidation states are not identified properly, Bader analysis finds them automatically
      treat_SbO2_special_case(structure, cce_flags, cce_vars);
      // for Pb3O4 the oxidation states are not identified properly, Bader analysis finds them automatically
      treat_Pb3O4_special_case(structure, cce_flags, cce_vars);
      // Ti-O Magneli phases need to be treated specially since oxidation numbers are not recognized appropriately by above approach
      treat_Ti_O_Magneli_phase_special_case(structure, cce_flags, cce_vars);
      // for Fe3O4 in inverse spinel structure, the oxidation states are not identified properly
      treat_Fe3O4_special_case(structure, cce_flags, cce_vars);
      // for Mn3O4 in spinel structure, the oxidation states are not identified properly
      treat_Mn3O4_special_case(structure, cce_flags, cce_vars);
      // for Co3O4 in spinel structure, the oxidation states are not identified properly
      treat_Co3O4_special_case(structure, cce_flags, cce_vars);
      // alkali metal sesquioxides need to be treated specially since oxidation numbers and number of per- and superoxide bonds are not recognized appropriately
      treat_alkali_sesquioxide_special_case(structure, cce_flags, cce_vars);
      if(cce_flags.flag("OX_STATES_DETERMINED") && cce_flags.flag("COMMAND_LINE")){
        oss << endl;
      }
    }
    // if preferred oxidation states approach and mixed valence doesn't work, try to find proper oxidation states (making oxidation sum zero) by using all known oxidation states
    if(!cce_flags.flag("OX_STATES_DETERMINED")){
      if(!cce_flags.flag("NO_OX_STATES")){ // for He, Ne, and Ar there are no known oxidation states
        try_all_oxidation_states(structure, cce_vars);
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = print_oxidation_states_and_sum(structure, cce_flags, cce_vars);
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
      } else{ // error message that oxidation numbers cannot be determined since for at least one species there are no known oxidation numbers is already included in load_ox_states_templates_each_species function
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

  //set_anion_oxidation_states////////////////////////////////////////////////////////
  // determine the oxidation numbers of the anions
  void set_anion_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::set_anion_oxidation_states():";
    if(LDEBUG){
      cerr << soliloquy << soliloquy << "ANION PART:" << endl;
    }
    double total_anion_charge=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == cce_vars.anion_species || cce_vars.multi_anion_atoms[i] == 1){
        if (structure.atoms[i].cleanname == cce_vars.anion_species){ // for multi-anion atoms oxidation states have been assigned previously
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
            cerr << soliloquy << "anion oxidation number for multi-anion ATOM[" << i << "] (" << structure.atoms[i].cleanname << ") has been assigned previously to: " << cce_vars.oxidation_states[i] << endl;
          } else {
            cerr << soliloquy << "anion oxidation number for ATOM[" << i << "] (" << structure.atoms[i].cleanname << "): " << cce_vars.oxidation_states[i] << endl;
          }
        }
        total_anion_charge += cce_vars.oxidation_states[i];
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "Total anion charge in cell: " << total_anion_charge << endl;
    }
  }

  //sort_species_by_electronegativity////////////////////////////////////////////////////////
  // sort species ascending by electronegativity
  void sort_species_by_electronegativity(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::sort_species_by_electronegativity():";
    // using aurostd sort functions
    vector<double> electronegativities_sorted = cce_vars.electronegativities;
    for(uint j=0,jsize=structure.species.size();j<jsize;j++){ //loop over all species
     cce_vars.species_electronegativity_sorted[j] = KBIN::VASP_PseudoPotential_CleanName(structure.species[j]);
    }
    // for the oxidation state algorithm the species must be sorted by electronegativity and the preferred 
    // and all known oxidation states will be automatically sorted by electronegativity in the subsequent loading
    // sort species by electronegativity
    aurostd::sort(electronegativities_sorted,cce_vars.species_electronegativity_sorted);
    for(uint j=0,jsize=structure.species.size();j<jsize;j++){ //loop over all species
      if(LDEBUG){
        cerr << soliloquy << "species_electronegativity_sorted[" << j << "]: " << cce_vars.species_electronegativity_sorted[j] << endl;
        cerr << soliloquy << "electronegativities_sorted[" << j << "]: " << electronegativities_sorted[j] << endl;
      }
    }
  }

  //load_ox_states_templates_each_species////////////////////////////////////////////////////////
  // load templates for preferred and other oxidation states
  void load_ox_states_templates_each_species(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::load_ox_states_templates_each_species():";
    _atom atom;
    uint z = 0;
    for(uint i=0,isize=structure.species.size();i<isize;i++){ 
      z = GetAtomNumber(KBIN::VASP_PseudoPotential_CleanName(cce_vars.species_electronegativity_sorted[i]));
      xelement element(z);
      // load preferred oxidation states for each species
      cce_vars.num_pref_ox_states_electronegativity_sorted[i] = element.oxidation_states_preferred.size();
      if(element.oxidation_states_preferred[0] != NNN){
        for (uint k=0,ksize=cce_vars.num_pref_ox_states_electronegativity_sorted[i];k<ksize;k++) {
          cce_vars.pref_ox_states_electronegativity_sorted[i].push_back(element.oxidation_states_preferred[k]);
        }
        if(LDEBUG){
          cerr << soliloquy << "num_pref_ox_states_electronegativity_sorted[" << i << "]: " << cce_vars.num_pref_ox_states_electronegativity_sorted[i] << endl;
          for (uint k=0,ksize=cce_vars.num_pref_ox_states_electronegativity_sorted[i];k<ksize;k++) {
            cerr << soliloquy << "preferred oxidation state " << k << " of species " << i << " (" << KBIN::VASP_PseudoPotential_CleanName(cce_vars.species_electronegativity_sorted[i]) << "): " <<  cce_vars.pref_ox_states_electronegativity_sorted[i][k] << endl;
          }
        }
      } else{
        cce_vars.num_pref_ox_states_electronegativity_sorted[i] = 0;
        cce_flags.flag("NO_PREF_OX_STATES",TRUE);
        if(LDEBUG){
          cerr << endl;
          cerr << soliloquy << "BAD NEWS: There are no preferred oxidation states for species " << i << " (" << KBIN::VASP_PseudoPotential_CleanName(cce_vars.species_electronegativity_sorted[i]) << ")."  << endl;
          cerr << soliloquy << "Therefore the oxidation states cannot be determined on this basis." << endl;
        }
      }
      // load all oxidation states for each species
      cce_vars.num_all_ox_states_electronegativity_sorted[i] = element.oxidation_states.size();
      if(element.oxidation_states[0] != NNN){
        for (uint k=0,ksize=cce_vars.num_all_ox_states_electronegativity_sorted[i];k<ksize;k++) {
          cce_vars.all_ox_states_electronegativity_sorted[i].push_back(element.oxidation_states[k]);
        }
        if(LDEBUG){
          cerr << soliloquy << "num_all_ox_states_electronegativity_sorted[" << i << "]: " << cce_vars.num_all_ox_states_electronegativity_sorted[i] << endl;
          for (uint k=0,ksize=cce_vars.num_all_ox_states_electronegativity_sorted[i];k<ksize;k++) {
            cerr << soliloquy << "all oxidation state " << k << " of species " << i << " (" << KBIN::VASP_PseudoPotential_CleanName(cce_vars.species_electronegativity_sorted[i]) << "): " << cce_vars.all_ox_states_electronegativity_sorted[i][k] << endl;
          }
          cerr << endl;
        }
      } else{
        cce_vars.num_all_ox_states_electronegativity_sorted[i] = 0;
        cce_flags.flag("NO_OX_STATES",TRUE);
        cerr << endl;
        cerr << "BAD NEWS: There are no known oxidation states for species " << i << " (" << KBIN::VASP_PseudoPotential_CleanName(cce_vars.species_electronegativity_sorted[i]) << ")."  << endl;
        cerr << "Therefore the oxidation states cannot be determined on this basis." << endl;
        cerr << endl;
      }
    }
  }

  //try_preferred_oxidation_states////////////////////////////////////////////////////////
  // try to determine the oxidation numbers using the preferred oxidation states for all cation species
  void try_preferred_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::try_preferred_oxidation_states():";
    if(LDEBUG){
      cerr << soliloquy << "Trying preferred oxidation numbers:" << endl;
    }
    // use Marco's implementation of my algorithm to determine oxidation_numbers
    determine_cation_oxidation_states(structure, cce_vars, cce_vars.pref_ox_states_electronegativity_sorted);
    // print oxidation numbers and calculate sum
    cce_vars.oxidation_sum = get_oxidation_states_sum(cce_vars);
    // if sum of oxidation numbers is essentially zero, oxidation states should be regarded as determined correctly
    if (std::abs(cce_vars.oxidation_sum) <= _CCE_OX_TOL_) {
      cce_flags.flag("OX_STATES_DETERMINED",TRUE);
      if(cce_flags.flag("COMMAND_LINE")){
        print_oxidation_states_and_sum(structure, cce_flags, cce_vars);
      }
    }
  }

  //treat_SbO2_special_case////////////////////////////////////////////////////////
  // for SbO2 the oxidation states are not identified properly
  // with the actual formula Sb2O4 it is a mixed valence oxide with one Sb+3 with 4 Sb-O bonds and one Sb+5 with 6 Sb-O bonds
  void treat_SbO2_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::treat_SbO2_special_case():";
    stringstream message;
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Sb") || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Sb" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_Sb = 0;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Sb" ) {
      num_O_before_Sb=4;
    } else {
      num_O_before_Sb=0;
    }
    double amount_Sb=0;
    double amount_O=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == "Sb"){
        amount_Sb+=1;
      } else if (structure.atoms[i].cleanname == "O"){
        amount_O+=1;
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of Sb ions= " << amount_Sb << endl;
      cerr << soliloquy << "number of O ions= " << amount_O << endl;
    }
    double Sb_O_ratio;
    if ( amount_O != 0 ){
      Sb_O_ratio=amount_Sb/amount_O;
    } else {
      message << "Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    if ( aurostd::isequal(Sb_O_ratio,0.5) ){
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
        oss << "WARNING: This system is identified as a mixed valence compound." << endl;
        oss << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
        oss << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
        oss << "Sb2O4 with ratio of Sb/O= " << Sb_O_ratio << endl;
      }
      uint num_formula_units_in_cell = 0;
      num_formula_units_in_cell=structure.atoms.size()/6; // 6 for Sb2O4
      if(cce_flags.flag("COMMAND_LINE")){
        oss << "number of formula units in cell: " << num_formula_units_in_cell << endl;
      }
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].cleanname == "Sb"){
          // taking the first Sb ions as +3 (without knowing whether they are actually +3); works only because 
          // the number of bonds for them will be adjusted to 4 as needed for Sb2O4 disregarding the actual 
          // number of Sb-O bonds (this is a hack since I know how many bonds there should be for each ion type)
          if ( i < (1+num_O_before_Sb)*num_formula_units_in_cell ){  // (1 Sb3+ ions per formula unit + 4*O listed before in alphabetic order) * number of formula units
            cce_vars.oxidation_states[i]=+3;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Sb+3 " << endl;
            }
            cce_vars.num_neighbors[i]=4; // for Sb2O4 the Sb3+ ions are 4-fold coordinated by oxygen
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "Modified number of neighbors of Sb (atom[" << i << "]) taken as Sb+3: " << cce_vars.num_neighbors[i] << endl;
            }
          } else {
            cce_vars.oxidation_states[i]=5;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Sb+5 " << endl;
            }
            cce_vars.num_neighbors[i]=6; // for Sb2O4 the Sb5+ ions are 6-fold coordinated by oxygen
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "Modified number of neighbors of Sb (atom[" << i << "]) taken as Sb+5: " << cce_vars.num_neighbors[i] << endl;
            }
          }
        }
      }
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
      }
      // print oxidation numbers and calculate sum
      cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
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

  //treat_Pb3O4_special_case////////////////////////////////////////////////////////
  // for Pb3O4 the oxidation states are not identified properly
  // https://en.wikipedia.org/wiki/Lead(II,IV)_oxide
  void treat_Pb3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    string soliloquy="cce::treat_Pb3O4_special_case():";
    stringstream message;
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Pb") || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Pb" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_Pb = 0;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Pb" ) {
      num_O_before_Pb=4;
    } else {
      num_O_before_Pb=0;
    }
    double amount_Pb=0;
    double amount_O=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == "Pb"){
        amount_Pb+=1;
      } else if (structure.atoms[i].cleanname == "O"){
        amount_O+=1;
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of Pb ions= " << amount_Pb << endl;
      cerr << soliloquy << "number of O ions= " << amount_O << endl;
    }
    double Pb_O_ratio = 0.0;
    if ( amount_O != 0 ){
      Pb_O_ratio=amount_Pb/amount_O;
    } else {
      message << "Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    if ( aurostd::isequal(Pb_O_ratio,0.75) ){
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
        oss << "WARNING: This system is identified as a mixed valence compound." << endl;
        oss << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
        oss << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
        oss << "Pb3O4 with ratio of Pb/O= " << Pb_O_ratio << endl;
      }
      uint num_formula_units_in_cell;
      num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Pb3O4
      if(cce_flags.flag("COMMAND_LINE")){
        oss << "number of formula units in cell: " << num_formula_units_in_cell << endl;
      }
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].cleanname == "Pb"){
          // taking the first Pb ions as +4 (without knowing whether they are actually +4); works only because 
          // the number of bonds for them will be adjusted to 6 as needed for Pb3O4 disregarding the actual 
          // number of Pb-O bonds (this is a hack since I know how many bonds there should be for each ion type)
          if ( i < (1+num_O_before_Pb)*num_formula_units_in_cell ){  // (1 Pb4+ ions per formula unit + 4*O listed before in alphabetic order) * number of formula units
            cce_vars.oxidation_states[i]=+4;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Pb+4 " << endl;
            }
            cce_vars.num_neighbors[i]=6; // for Pb3O4 the Pb4+ ions are 6-fold coordinated by oxygen
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "Modified number of neighbors of Pb (atom[" << i << "]) taken as Pb+4: " << cce_vars.num_neighbors[i] << endl;
            }
          } else {
            cce_vars.oxidation_states[i]=+2;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Pb+2 " << endl;
            }
            cce_vars.num_neighbors[i]=3; // for Pb3O4 the Pb2+ ions are 3-fold coordinated by oxygen
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "Modified number of neighbors of Pb (atom[" << i << "]) taken as Pb+2: " << cce_vars.num_neighbors[i] << endl;
            }
          }
        }
      }
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
      }
      // print oxidation numbers and calculate sum
      cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
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

  //treat_Ti_O_Magneli_phase_special_case////////////////////////////////////////////////////////
  // treat Ti-O Magneli phases; there are always 2xTi+3 per formula unit and the rest is Ti+4; 
  // fortunately, both Ti+3 and Ti+4 have 6 Ti-O bonds, hence one only needs to know how many ions 
  // of the respective oxidation state there are, not which one is which
  void treat_Ti_O_Magneli_phase_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::treat_Ti_O_Magneli_phase_special_case():";
    stringstream message;
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Ti") || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Ti" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_Ti = 0;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Ti" ) {
      num_O_before_Ti=5;
    } else {
      num_O_before_Ti=0;
    }
    if(LDEBUG){
      cerr << soliloquy << "Ti-O system, Magneli for Ti_nO_(2n-1), i.e. Ti-O ratio= " << 3.0/5 << ", " << 4.0/7 << ", " << 5.0/9 << ", " << 6.0/11 << ", " << 7.0/13 << ", " << 8.0/15 << ", " << 9.0/17 << ", " << 10.0/19 << "..." << endl;
    }
    double amount_O=0;
    double amount_Ti=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == "O"){
        amount_O+=1;
      } else if (structure.atoms[i].cleanname == "Ti"){
        amount_Ti+=1;
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of O ions= " << amount_O << endl;
      cerr << soliloquy << "number of Ti ions= " << amount_Ti << endl;
    }
    double Ti_O_ratio;
    if ( amount_O != 0 ){
      Ti_O_ratio=amount_Ti/amount_O;
    } else {
      message << "Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    if(LDEBUG){
      cerr << soliloquy << "ratio of Ti/O= " << Ti_O_ratio << endl;
    }
    uint num_formula_units_in_cell;
    // check for Magneli composition Ti_(n)O_(2n-1)
    double n;
    bool magneli = false;
    for(n=3;n<101;n++){
      //oss << "n/(2*n-1)= " << n/(2*n-1) << endl;
      if ( aurostd::isequal(Ti_O_ratio,n/(2*n-1)) ){
        if(cce_flags.flag("COMMAND_LINE")){
          oss << endl;
          oss << "WARNING: This system is identified as a mixed valence compound." << endl;
          oss << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
          oss << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
          oss << "n= " << n << " Magneli composition Ti_nO_(2n-1)" << endl;
        }
        magneli = true;
        num_formula_units_in_cell=amount_Ti/n;
        if(cce_flags.flag("COMMAND_LINE")){
          oss << "number of formula units in cell: " << num_formula_units_in_cell << endl;
        }
      }
    }
    if ( magneli == false){
      if(LDEBUG){
        cerr << soliloquy << "Not a Magneli composition." << endl;
      }
    }
    if (magneli){
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].cleanname == "Ti"){
          // taking the first Ti as +3 (without knowing whether they are actually +3) works only because 
          // for the Magneli phases both Ti+3 and Ti+4 have both always 6 Ti-O bonds
          if ( i < (2+num_O_before_Ti)*num_formula_units_in_cell ){
            cce_vars.oxidation_states[i]=+3;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Ti+3 " << endl;
            }
          } else {
            cce_vars.oxidation_states[i]=+4;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Ti+4 " << endl;
            }
          }
        }
      }
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
      }
      // print oxidation numbers and calculate sum
      cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
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

  //treat_Fe3O4_special_case////////////////////////////////////////////////////////
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
  void treat_Fe3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::treat_Fe3O4_special_case():";
    stringstream message;
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Fe") || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Fe" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_Fe;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Fe" ) {
      num_O_before_Fe=4;
    } else {
      num_O_before_Fe=0;
    }
    double amount_Fe=0;
    double amount_O=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == "Fe"){
        amount_Fe+=1;
      } else if (structure.atoms[i].cleanname == "O"){
        amount_O+=1;
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of Fe ions= " << amount_Fe << endl;
      cerr << soliloquy << "number of O ions= " << amount_O << endl;
    }
    double Fe_O_ratio;
    if ( amount_O != 0 ){
      Fe_O_ratio=amount_Fe/amount_O;
    } else {
      message << "Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    if ( aurostd::isequal(Fe_O_ratio,0.75) ){
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
        oss << "WARNING: This system is identified as a mixed valence compound." << endl;
        oss << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
        oss << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
        oss << "Fe3O4 with ratio of Fe/O= " << Fe_O_ratio << endl;
      }
      uint num_formula_units_in_cell;
      num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Fe3O4
      if(cce_flags.flag("COMMAND_LINE")){
        oss << "number of formula units in cell: " << num_formula_units_in_cell << endl;
      }
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].cleanname == "Fe"){
          // taking the first Fe ions as +2 (without knowing whether they are actually +2); works only because 
          // the number of bonds for them will be adjusted to 6 (octahedral) as needed for Fe3O4 disregarding the actual 
          // number of Fe-O bonds (this is a hack since I know how many bonds there should be for each ion type)
          if ( i < (1+num_O_before_Fe)*num_formula_units_in_cell ){  // 1 Fe2+ ions per formula unit * number of formula units
            cce_vars.oxidation_states[i]=+2;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+2 " << endl;
            }
            cce_vars.num_neighbors[i]=6; // for Fe3O4 the Fe2+ ions are 6-fold coordinated by oxygen according to Wikipedia
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "Modified number of neighbors of Fe (atom[" << i << "]) taken as Fe+2: " << cce_vars.num_neighbors[i] << endl;
            }
          } else {
            cce_vars.oxidation_states[i]=+3;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+3 " << endl;
            }
            cce_vars.num_neighbors[i]=5; // for Fe3O4 the Fe3+ ions are evenly 6- and 4-fold, so on average 5-fold (set here as a hack) coordinated by oxygen according to Wikipedia
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "Modified number of neighbors of Fe (atom[" << i << "]) taken as Fe+3 (average between even 6- and 4-fold coordination): " << cce_vars.num_neighbors[i] << endl;
            }
          }
        }
      }
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
      }
      // print oxidation numbers and calculate sum
      cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
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

  //treat_Mn3O4_special_case////////////////////////////////////////////////////////
  // for Co3O4 and Mn3O4 in the normal spinel structure Co2+/Mn2+ occupies only tetrahedral sites 
  // while Co3+/Mn3+ occupies octahedral sites; the correction hence needs to be different than for Fe3O4
  // but at the present stage there are no corrections for Co3+ and Mn3+, see:
  // https://en.wikipedia.org/wiki/Cobalt(II,III)_oxide
  // https://en.wikipedia.org/wiki/Manganese(II,III)_oxide
  void treat_Mn3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::treat_Mn3O4_special_case():";
    stringstream message;
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Mn") || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Mn" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_Mn;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Mn" ) {
      num_O_before_Mn=4;
    } else {
      num_O_before_Mn=0;
    }
    double amount_Mn=0;
    double amount_O=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == "Mn"){
        amount_Mn+=1;
      } else if (structure.atoms[i].cleanname == "O"){
        amount_O+=1;
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of Mn ions= " << amount_Mn << endl;
      cerr << soliloquy << "number of O ions= " << amount_O << endl;
    }
    double Mn_O_ratio;
    if ( amount_O != 0 ){
      Mn_O_ratio=amount_Mn/amount_O;
    } else {
      message << "Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    if ( aurostd::isequal(Mn_O_ratio,0.75) ){
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
        oss << "WARNING: This system is identified as a mixed valence compound." << endl;
        oss << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
        oss << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
        oss << "Mn3O4 with ratio of Mn/O= " << Mn_O_ratio << endl;
      }
      uint num_formula_units_in_cell;
      num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Mn3O4
      if(cce_flags.flag("COMMAND_LINE")){
        oss << "number of formula units in cell: " << num_formula_units_in_cell << endl;
      }
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].cleanname == "Mn"){
          // taking the first Mn ions as +2 (without knowing whether they are actually +2); works only because 
          // the number of bonds for them will be adjusted to 4 (tetrahedral) as needed for Mn3O4 disregarding the actual 
          // number of Mn-O bonds (this is a hack since I know how many bonds there should be for each ion type)
          if ( i < (1+num_O_before_Mn)*num_formula_units_in_cell ){  // 1 Mn2+ ions per formula unit * number of formula units
            cce_vars.oxidation_states[i]=+2;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Mn+2 " << endl;
            }
            cce_vars.num_neighbors[i]=4; // for Mn3O4 the Mn2+ ions are 4-fold coordinated by oxygen according to Wikipedia
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "Modified number of neighbors of Mn (atom[" << i << "]) taken as Mn+2: " << cce_vars.num_neighbors[i] << endl;
            }
          } else {
            cce_vars.oxidation_states[i]=+3;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Mn+3 " << endl;
            }
            cce_vars.num_neighbors[i]=6; // for Mn3O4 the Mn3+ ions are 6-fold coordinated by oxygen according to Wikipedia
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "Modified number of neighbors of Mn (atom[" << i << "]) taken as Mn+3: " << cce_vars.num_neighbors[i] << endl;
            }
          }
        }
      }
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
      }
      // print oxidation numbers and calculate sum
      cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
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

  //treat_Co3O4_special_case////////////////////////////////////////////////////////
  // for Co3O4 and Mn3O4 in the normal spinel structure Co2+/Mn2+ occupies only tetrahedral sites 
  // while Co3+/Mn3+ occupies octahedral sites; the correction hence needs to be different than for Fe3O4
  // but at the present stage there are no corrections for Co3+ and Mn3+, see:
  // https://en.wikipedia.org/wiki/Cobalt(II,III)_oxide
  // https://en.wikipedia.org/wiki/Manganese(II,III)_oxide
  void treat_Co3O4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::treat_Co3O4_special_case():";
    stringstream message;
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Co") || (KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Co" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;}
    uint num_O_before_Co;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Co" ) {
      num_O_before_Co=4;
    } else {
      num_O_before_Co=0;
    }
    double amount_Co=0;
    double amount_O=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == "Co"){
        amount_Co+=1;
      } else if (structure.atoms[i].cleanname == "O"){
        amount_O+=1;
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of Co ions= " << amount_Co << endl;
      cerr << soliloquy << "number of O ions= " << amount_O << endl;
    }
    double Co_O_ratio;
    if ( amount_O != 0 ){
      Co_O_ratio=amount_Co/amount_O;
    } else {
      message << "Amount of O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    if ( aurostd::isequal(Co_O_ratio,0.75) ){
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
        oss << "WARNING: This system is identified as a mixed valence compound." << endl;
        oss << "The oxidation numbers and the number of cation-anion bonds will be set as known for this special case." << endl;
        oss << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each cation oxidation state occurs should be correct." << endl;
        oss << "Co3O4 with ratio of Co/O= " << Co_O_ratio << endl;
      }
      uint num_formula_units_in_cell;
      num_formula_units_in_cell=structure.atoms.size()/7; // 7 for Co3O4
      if(cce_flags.flag("COMMAND_LINE")){
        oss << "number of formula units in cell: " << num_formula_units_in_cell << endl;
      }
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].cleanname == "Co"){
          // taking the first Co ions as +2 (without knowing whether they are actually +2); works only because 
          // the number of bonds for them will be adjusted to 4 (tetrahedral) as needed for Co3O4 disregarding the actual 
          // number of Co-O bonds (this is a hack since I know how many bonds there should be for each ion type)
          if ( i < (1+num_O_before_Co)*num_formula_units_in_cell ){  // 1 Co2+ ions per formula unit * number of formula units
            cce_vars.oxidation_states[i]=+2;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Co+2 " << endl;
            }
            cce_vars.num_neighbors[i]=4; // for Co3O4 the Co2+ ions are 4-fold coordinated by oxygen according to Wikipedia
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "Modified number of neighbors of Co (atom[" << i << "]) taken as Co+2: " << cce_vars.num_neighbors[i] << endl;
            }
          } else {
            cce_vars.oxidation_states[i]=+3;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Co+3 " << endl;
            }
            cce_vars.num_neighbors[i]=6; // for Co3O4 the Co3+ ions are 6-fold coordinated by oxygen according to Wikipedia
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "Modified number of neighbors of Co (atom[" << i << "]) taken as Co+3: " << cce_vars.num_neighbors[i] << endl;
            }
          }
        }
      }
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
      }
      // print oxidation numbers and calculate sum
      cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
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

  //treat_alkali_sesquioxide_special_case////////////////////////////////////////////////////////
  // for alkali metal sesquioxides (X2O3) the oxidation states and corrections might not be found correctly
  // since these systems are formulated to contain both per- and superoxide ions
  // see: https://en.wikipedia.org/wiki/Sesquioxide
  // the only known example to date for which also a structure was found on Springer Materials and the AFLOW ICSD
  // is Rb2O3 formulated as [(Rb+)4(O2−2)(O2-1)2] (for Cs4O6 there is only a single entry in a scanned pdf on Springer Materials)
  // this function will treat it as an exceptional case for all possible alkali metal sesquioxides
  // the oxidation numbers and per- as well as superoxide corrections will be adjusted accordingly
  void treat_alkali_sesquioxide_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::treat_alkali_sesquioxide_special_case():";
    stringstream message;
    string alkali_metals = "Li,Na,K,Rb,Cs,Fr";
    vector<string> valkali_metals;
    aurostd::string2tokens(alkali_metals, valkali_metals, ",");
    if (! ( structure.species.size() == 2 && ((KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && aurostd::withinList(valkali_metals, KBIN::VASP_PseudoPotential_CleanName(structure.species[1]))) || (aurostd::withinList(valkali_metals, KBIN::VASP_PseudoPotential_CleanName(structure.species[0])) && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O")) )) {return;} // check whether it is a binary alkali metal oxide
    if(LDEBUG){
      cerr << soliloquy << "This is a binary alkali metal oxide, checking whether it is an alkali metal sesquioxide..." << endl;
    }
    uint num_alkali_before_O; // num cations before O not O before cations since setting oxidation states of anions below, not for cations as in other cases
    string alkali_metal;
    if ( KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "O" && aurostd::withinList(valkali_metals, KBIN::VASP_PseudoPotential_CleanName(structure.species[1])) ) {
      num_alkali_before_O=0;
      alkali_metal=KBIN::VASP_PseudoPotential_CleanName(structure.species[1]);
    } else {
      num_alkali_before_O=2;
      alkali_metal=KBIN::VASP_PseudoPotential_CleanName(structure.species[0]);
    }
    double amount_O=0;
    double amount_alkali=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == "O"){
        amount_O+=1;
      } else {
        amount_alkali+=1;
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of O ions= " << amount_O << endl;
      cerr << soliloquy << "number of alkali (" << alkali_metal << ") ions= " << amount_alkali << endl;
    }
    double O_alkali_ratio;
    if ( amount_alkali != 0 ){
      O_alkali_ratio=amount_O/amount_alkali;
    } else {
      message << "Amount of alkali atoms determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
    if(LDEBUG){
      cerr << soliloquy << "ratio of O/" << alkali_metal << "= " << O_alkali_ratio << endl;
    }
    // check for sesqui-composition alkali_metal2O3
    if ( aurostd::isequal(O_alkali_ratio,1.5) ){
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
        oss << "WARNING: This system is identified as an alkali metal sesquioxide (formally) containing both per- and superoxide ions." << endl;
        oss << "The oxidation numbers and the number of per- and superoxide bonds will be set as known for this special case." << endl;
        oss << "The individual oxidation numbers might therefore not be assigned to the correct atoms." << endl; //, but at least how often each oxidation state occurs should be correct." << endl;
        oss << alkali_metal << "2O3 with ratio of O/" << alkali_metal << "= " << O_alkali_ratio << endl;
      }
      uint num_formula_units_in_cell;
      num_formula_units_in_cell=structure.atoms.size()/5; // 5 for alkali_metal2O3
      if(cce_flags.flag("COMMAND_LINE")){
        oss << "number of formula units in cell: " << num_formula_units_in_cell << endl;
      }
      for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
        if (structure.atoms[i].cleanname == "O"){
          // taking the first O ions as -0.5 (superox. Os, without knowing whether they are actually -0.5); works only because 
          // the number of superoxide bonds will be adjusted to 1 per formula unit as needed for alkali metal sesquioxides
          // taking last O ions as -1 (perox. Os, without knowing whether they are actually -1); works only because 
          // the number of peroxide bonds will be adjusted to 0.5 per formula unit as needed for alkali metal sesquioxides
          // (this is a hack since I know how many bonds there should be for each ion type)
          if ( i < (2+num_alkali_before_O)*num_formula_units_in_cell ){  // 2 superoxide O-atoms per formula unit * number of formula units
            cce_vars.oxidation_states[i]=-0.5;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to O-0.5 " << endl;
            }
            cce_vars.superox_indices[i]=1;
          } else {
            cce_vars.oxidation_states[i]=-1;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to O-1 " << endl;
            }
            cce_vars.perox_indices[i]=1;
          }
        }
      }
      cce_vars.num_superox_bonds=1*num_formula_units_in_cell; // 1 superox. bonds per formula unit * num. formula units
      if(cce_flags.flag("COMMAND_LINE")){
        oss << "Number of superoxide bonds set to: " << cce_vars.num_superox_bonds << endl;
      }
      cce_vars.num_perox_bonds=0.5*num_formula_units_in_cell; // 0.5 perox. bonds per formula unit * num. formula units; this should neve yield a non-integer since there should be at least one complete peroxide ion per cell
      if(cce_flags.flag("COMMAND_LINE")){
        oss << "Number of peroxide bonds set to: " << cce_vars.num_perox_bonds << endl;
      }
      if(cce_flags.flag("COMMAND_LINE")){
        oss << endl;
      }
      // print oxidation numbers and calculate sum
      cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
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

  // following special cases only needed when determining oxidation states from Bader charges

  //treat_MnMoO4_special_case////////////////////////////////////////////////////////
  // for MnMoO4 the oxidation numbers are determined to be +4 for both Mn and Mo from the Bader charges; 
  // it should be Mn+2 & Mo+6; however since the sum of the falsely determined oxidation numbers 
  // is accidentally 0 it needs to be corrected individually
  void treat_MnMoO4_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::treat_MnMoO4_special_case():";
    stringstream message;
    if (! ( structure.species.size() == 3 && KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Mn" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Mo" && KBIN::VASP_PseudoPotential_CleanName(structure.species[2]) == "O" )) {return;}
    double amount_Mn=0;
    double amount_Mo=0;
    double amount_O=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == "Mn"){
        amount_Mn+=1;
      } else if (structure.atoms[i].cleanname == "Mo"){
        amount_Mo+=1;
      } else if (structure.atoms[i].cleanname == "O"){
        amount_O+=1;
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of Mn ions= " << amount_Mn << endl;
      cerr << soliloquy << "number of Mo ions= " << amount_Mo << endl;
      cerr << soliloquy << "number of O ions= " << amount_O << endl;
    }
    if (amount_Mo != 0 && amount_O != 0) {
      if (aurostd::isequal(amount_Mn/amount_Mo,1.0) && aurostd::isequal(amount_Mn/amount_O,0.25) && aurostd::isequal(amount_Mo/amount_O,0.25)) {
        if(cce_flags.flag("COMMAND_LINE")){
          oss << "MnMoO4 special treatment since sum over oxdiation states is zero but individual oxidation numbers are wrong!!!" << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Mn"){
            cce_vars.oxidation_states[i]=+2;
            if(cce_flags.flag("COMMAND_LINE")){
             oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Mn+2 " << endl;
            }
          }
          if (structure.atoms[i].cleanname == "Mo"){
            cce_vars.oxidation_states[i]=+6;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Mo+6 " << endl;
            }
          }
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        }
      }
    } else {
      message << "Amount of Mo or O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
  }

  //treat_Ca2Fe2O5_CaFe2O4_LDA_special_case////////////////////////////////////////////////////////
  // for Ca2Fe2O5 and CaFe2O4 for LDA the oxidation numbers of Fe are not correctly determined 
  // to be Fe+3 but are partly Fe+2 which will be corrected here
  void treat_Ca2Fe2O5_CaFe2O4_LDA_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::treat_Ca2Fe2O5_CaFe2O4_LDA_special_case():";
    stringstream message;
    if (! ( structure.species.size() == 3 && KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Ca" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "Fe" && KBIN::VASP_PseudoPotential_CleanName(structure.species[2]) == "O" )) {return;}
    double amount_Ca=0;
    double amount_Fe=0;
    double amount_O=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == "Ca"){
        amount_Ca+=1;
      } else if (structure.atoms[i].cleanname == "Fe"){
        amount_Fe+=1;
      } else if (structure.atoms[i].cleanname == "O"){
        amount_O+=1;
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of Ca ions= " << amount_Ca << endl;
      cerr << soliloquy << "number of Fe ions= " << amount_Fe << endl;
      cerr << soliloquy << "number of O ions= " << amount_O << endl;
    }
    //making sure it is Ca2Fe2O5
    if (amount_Fe != 0 && amount_O != 0) {
      if (aurostd::isequal(amount_Ca/amount_Fe,1.0) && aurostd::isequal(amount_Ca/amount_O,0.4) && aurostd::isequal(amount_Fe/amount_O,0.4)) {
        if(cce_flags.flag("COMMAND_LINE")){
          oss << "Ca2Fe2O5 special treatment for LDA since oxidation numbers for Fe, which should be Fe+3, are not correctly determined from Bader charges for all Fe!!!" << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Fe"){
            cce_vars.oxidation_states[i]=+3;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+3 " << endl;
            }
          }
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        }
      //making sure it is CaFe2O4
      } else if (amount_Ca/amount_Fe == 0.5 && amount_Ca/amount_O == 0.25 && amount_Fe/amount_O == 0.5) {
        if(cce_flags.flag("COMMAND_LINE")){
          oss << "CaFe2O4 special treatment for LDA since oxidation numbers for Fe, which should be Fe+3, are not correctly determined from Bader charges for all Fe!!!" << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Fe"){
            cce_vars.oxidation_states[i]=+3;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+3 " << endl;
            }
          }
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        }
      }
    } else {
      message << "Amount of Fe or O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
  }

  //treat_FeTiO3_LDA_special_case////////////////////////////////////////////////////////
  // for FeTiO3 for LDA the oxidation numbers of Ti are not correctly determined to be Ti+4 and using 
  // the general fixes would modify both the Ti AND the Fe oxidation numbers resulting again 
  // in non-zero oxidation number sum, which is fixed here
  void treat_FeTiO3_LDA_special_case(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::treat_FeTiO3_LDA_special_case():";
    stringstream message;
    if (! ( structure.species.size() == 3 && KBIN::VASP_PseudoPotential_CleanName(structure.species[0]) == "Fe" && KBIN::VASP_PseudoPotential_CleanName(structure.species[1]) == "O" && KBIN::VASP_PseudoPotential_CleanName(structure.species[2]) == "Ti" )) {return;}
    double amount_Fe=0;
    double amount_O=0;
    double amount_Ti=0;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname == "Fe"){
        amount_Fe+=1;
      } else if (structure.atoms[i].cleanname == "O"){
        amount_O+=1;
      } else if (structure.atoms[i].cleanname == "Ti"){
        amount_Ti+=1;
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "number of Fe ions= " << amount_Fe << endl;
      cerr << soliloquy << "number of O ions= " << amount_O << endl;
      cerr << soliloquy << "number of Ti ions= " << amount_Ti << endl;
    }
    //making sure it is FeTiO3
    if (amount_Ti != 0 && amount_O != 0) {
      if (aurostd::isequal(amount_Fe/amount_Ti,1.0) && aurostd::isequal(amount_Fe/amount_O,1.0/3) && aurostd::isequal(amount_Ti/amount_O,1.0/3)) {
        if(cce_flags.flag("COMMAND_LINE")){
          oss << "FeTiO3 special treatment for LDA since oxidation numbers for Ti, which should be Ti+4, are not correctly determined from Bader charges and using other fixing would also change Fe oxidation numbers!!!" << endl;
        }
        for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
          if (structure.atoms[i].cleanname == "Ti"){
            cce_vars.oxidation_states[i]=+4;
            if(cce_flags.flag("COMMAND_LINE")){
              oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Ti+4 " << endl;
            }
          }
        }
        // print oxidation numbers and calculate sum
        cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
        // system should not be regarded correctable if sum over oxidation states is not zero
        if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
          cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
        }
      }
    } else {
      message << "Amount of Ti or O determined to be ZERO. Please check your structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_VALUE_ILLEGAL_);
    }
  }

  //try_all_oxidation_states////////////////////////////////////////////////////////
  // try to determine the oxidation numbers using all known oxidation states for all cation species
  void try_all_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::try_all_oxidation_states():";
    if(LDEBUG){
      cerr << soliloquy << "Trying all known oxidation numbers:" << endl;
    }
    // use Marco's implementation of my algorithm to determine oxidation_numbers
    determine_cation_oxidation_states(structure, cce_vars, cce_vars.all_ox_states_electronegativity_sorted);
  }

  //determine_cation_oxidation_states////////////////////////////////////////////////////////
  // ME Nov. 2019
  // for avoiding recursion algorithm to determine cation oxidation_numbers
  // determine the cation oxidation numbers by using possible oxidation states for each species 
  // which can be either the preferred or all known oxidation numbers
  void determine_cation_oxidation_states(const xstructure& structure, CCE_Variables& cce_vars, const vector<vector<double> >& possible_ox_states) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::determine_cation_oxidation_states():";
    uint num_cation_species = cce_vars.cations_map.size(); // the number of cation species
    uint natoms = cce_vars.oxidation_states.size();

    // Initialize
    for (uint i = 0; i < num_cation_species; i++) { // loop over all cations
      for (uint j = 0; j < cce_vars.cations_map[i].size(); j++) { // loop over atoms of the ith cation type that is given by the second index of cation_map
        if (cce_vars.multi_anion_atoms[cce_vars.cations_map[i][j]] != 1){ // exclude atoms that have been identified as multi anion atoms previously
          if(LDEBUG){
            cerr << soliloquy << "Updating oxidation number for atom " << cce_vars.cations_map[i][j] << " (" << structure.atoms[cce_vars.cations_map[i][j]].cleanname << ") to " << possible_ox_states[i][0] << endl;
          }
          cce_vars.oxidation_states[cce_vars.cations_map[i][j]] = possible_ox_states[i][0]; // possible_ox_states (either preferred or all) for all species should be electronegativity sorted
        }
      }
    }
    if(LDEBUG){
      for(uint n=0,nsize=structure.atoms.size();n<nsize;n++){ 
        cerr << soliloquy << "chosen oxidation state for atom " << n << " (" << structure.atoms[n].cleanname << "): " << cce_vars.oxidation_states[n] << endl;
      }
    }
    // check
    double total_ox = 0.0;
    for (uint at = 0; at < natoms; at++) total_ox += cce_vars.oxidation_states[at];
    bool iszero = aurostd::isequal(total_ox, 0.0, _CCE_OX_TOL_);
    // indicate whether these oxidation states were successful or not
    if (iszero) {
      if(LDEBUG){
        cerr << soliloquy << "Oxidation numbers successfully determined since oxidation sum is " << total_ox << endl;
      }
    } else{
      if(LDEBUG){
        cerr << soliloquy << "No successful determination of oxidation numbers since oxidation sum is " << total_ox << endl;
      }
    }

    // Cycle through oxidation states
    if (!iszero) { // maybe the oxidation numbers from the initialization work already, then nothing else needs to be done
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
                  cerr << soliloquy << "Updating oxidation number for atom " << cce_vars.cations_map[j-1][k] << "(" << structure.atoms[cce_vars.cations_map[j-1][k]].cleanname << ") to " << possible_ox_states[j-1][i] << endl;
                }
                cce_vars.oxidation_states[cce_vars.cations_map[j-1][k]] = possible_ox_states[j-1][i]; // i goes over all preferred/all oxidation states
              }
            }
            if(LDEBUG){
              for(uint n=0,nsize=structure.atoms.size();n<nsize;n++){ 
                cerr << soliloquy << "chosen oxidation state for atom " << n << " (" << structure.atoms[n].cleanname << "): " << cce_vars.oxidation_states[n] << endl;
              }
            }
            // check
            total_ox = 0.0;
            for (uint at = 0; at < natoms; at++) total_ox += cce_vars.oxidation_states[at];
            iszero = aurostd::isequal(total_ox, 0.0, _CCE_OX_TOL_);
            // indicate whether these oxidation states were successful or not
            if (iszero) {
              if(LDEBUG){
                cerr << soliloquy << "Oxidation numbers successfully determined since oxidation sum is " << total_ox << endl;
              }
              break; // BREAK for loop over possible oxidation states k for species j-1 upon success (sum over oxidation numbers is zero)
            } else{
              if(LDEBUG){
                cerr << soliloquy << "No successful determination of oxidation numbers since oxidation sum is " << total_ox << endl;
              }
            }
          }
        }
        i++;
      }
    }
  }

  //print_oxidation_states_and_sum////////////////////////////////////////////////////////
  // print out previously determined oxidation numbers and sum
  double print_oxidation_states_and_sum(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    if(cce_flags.flag("COMMAND_LINE")){
      oss << endl;
    }
    cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
    if(cce_flags.flag("COMMAND_LINE")){
      oss << endl;
    }
    return cce_vars.oxidation_sum;
  }

  //print_oxidation_states_and_get_sum////////////////////////////////////////////////////////
  // print out previously determined oxidation numbers and get sum
  double print_oxidation_states_and_get_sum(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    // print oxidation numbers
    if(cce_flags.flag("COMMAND_LINE")){
      oss << "OXIDATION NUMBERS:" << endl;
      for (uint k=0,ksize=structure.atoms.size();k<ksize;k++){
        oss << "Oxidation state of " << structure.atoms[k].cleanname << " (atom[" << k << "]): " << cce_vars.oxidation_states[k] << endl;
      }
      oss << "CHECK whether this is what you are expecting!" << endl;
    }
    // calculate and print sum of oxidation numbers
    cce_vars.oxidation_sum = get_oxidation_states_sum(cce_flags, cce_vars);
    return cce_vars.oxidation_sum;
  }

  //get_oxidation_states_sum////////////////////////////////////////////////////////
  // Calculate and return the sum of all oxidation states without writing output.
  double get_oxidation_states_sum(CCE_Variables& cce_vars) {
    aurostd::xoption cce_flags = init_flags();
    cce_flags.flag("COMMAND_LINE",FALSE);
    return get_oxidation_states_sum(cce_flags, cce_vars);
  }

  //get_oxidation_states_sum////////////////////////////////////////////////////////
  // Calculate and return the sum of all oxidation states.
  double get_oxidation_states_sum(xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    cce_vars.oxidation_sum=0;
    for (uint k=0,ksize=cce_vars.oxidation_states.size();k<ksize;k++){
      cce_vars.oxidation_sum += cce_vars.oxidation_states[k];
    }
    if(cce_flags.flag("COMMAND_LINE")){
      oss << "Sum over all oxidation numbers is (should be ZERO): " << cce_vars.oxidation_sum << endl;
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

  //get_oxidation_states_from_Bader////////////////////////////////////////////////////////
  // determine the oxidation numbers of the ions by an analysis of the Bader charges of the system and handle known exceptional cases explicitly
  vector<double> get_oxidation_states_from_Bader(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::get_oxidation_states_from_Bader():";
    if(LDEBUG){
      cerr << soliloquy << "DETERMINATION OF OXIDATION NUMBERS FROM BADER CHARGES:" << endl;
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
      get_system_name_functional_from_aflow_in(structure, cce_flags, cce_vars, system_name, functional);
      if (cce_flags.flag("CORRECTABLE")){
        // check whether Bader file exists
        string Bader_file = system_name + "_abader.out";
        if (aurostd::FileExist(Bader_file) || aurostd::EFileExist(Bader_file)) {
          // determine Bader charges from Bader file and store in array/vector Bader_charges; O Bader charges will be included although they might not be used
          cce_vars.Bader_charges = get_Bader_charges_from_Bader_file(structure, cce_vars, Bader_file);
          // determine oxidation numbers from Bader charges
          cce_vars.oxidation_states = Bader_charges_to_oxidation_states(structure, cce_flags, cce_vars, functional);
          // check whether sum of oxidation numbers is equal zero and try to correct for known special cases (for some cases even if it is zero)
          if (cce_flags.flag("CORRECTABLE")){
            // print oxidation numbers and calculate sum
            cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
            // SECOND point to deal with special cases known from (binary+ternary) oxides
            // Ti-O Magneli phases need to be treated specially since oxidation numbers are not recognized appropriately for all functionals
            cce_flags.flag("OX_STATES_DETERMINED",TRUE); // needed for algorithm determining oxidation numbers from electronegativities
            treat_Ti_O_Magneli_phase_special_case(structure, cce_flags, cce_vars);
            // for Fe3O4 in inverse spinel structure, the oxidation states are not identified properly via the Bader charges.
            treat_Fe3O4_special_case(structure, cce_flags, cce_vars);
            // for Mn3O4 in spinel structure, the oxidation states are not identified properly
            treat_Mn3O4_special_case(structure, cce_flags, cce_vars);
            // for Co3O4 in spinel structure, the oxidation states are not identified properly
            treat_Co3O4_special_case(structure, cce_flags, cce_vars);
            // MnMoO4 needs to be treated specially since oxidation numbers are not recognized appropriately
            treat_MnMoO4_special_case(structure, cce_flags, cce_vars);
            if (functional == "LDA") {
              // Ca2Fe2O5 & CaFe2O4 need to be treated specially for LDA since oxidation numbers are not recognized appropriately
              treat_Ca2Fe2O5_CaFe2O4_LDA_special_case(structure, cce_flags, cce_vars);
              // FeTiO3 needs to be treated specially for LDA since oxidation numbers are not recognized appropriately
              treat_FeTiO3_LDA_special_case(structure, cce_flags, cce_vars);
            }
            if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
              // general scheme to repair wrong oxidation numbers based on changing oxidation state of several ions known to be problematic
              general_attempt_fixing_oxidation_states(structure, cce_flags, cce_vars);
              // print oxidation numbers and calculate sum
              cce_vars.oxidation_sum = print_oxidation_states_and_get_sum(structure, cce_flags, cce_vars);
              // system should not be regarded correctable if sum over oxidation states is not zero
              if (std::abs(cce_vars.oxidation_sum) > _CCE_OX_TOL_) {
                cce_flags.flag("CORRECTABLE",FALSE);
                cerr << "BAD NEWS: The formation energy of this system is not correctable! The determined and fixed oxidation numbers do not add up to zero!"  << endl;
                cerr << "You can also provide oxidation numbers as a comma separated list as input via the option --oxidation_numbers=." << endl;
              }
            }
            if(cce_flags.flag("COMMAND_LINE")){
              oss << endl; // empty line to separate formation enthalpy values from the rest
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

  //get_system_name_functional_from_aflow_in////////////////////////////////////////////////////////
  // determine the system name and functional from the aflow_in
  void get_system_name_functional_from_aflow_in(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, string& system_name, string& functional, ostream& oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::get_system_name_functional_from_aflow_in():";
    string outcar="OUTCAR.relax1";
    functional=get_functional_from_aflow_in(structure, _AFLOWIN_,outcar);
    string aflowIn = aurostd::RemoveComments(aurostd::file2string(_AFLOWIN_));
    vector<string> vlines = aurostd::string2vectorstring(aflowIn);
    string line_a;
    for (uint i = 0; i < vlines.size(); i++) { 
      line_a = aurostd::RemoveSpaces(vlines[i]);
      if (line_a.find("[AFLOW]SYSTEM=") != string::npos){ // string::npos is returned if string is not found
        system_name = aurostd::RemoveSubStringFirst(line_a, "[AFLOW]SYSTEM=");
        system_name=aurostd::RemoveSubString(system_name,"\r"); // remove carriage return characters at end of string in case they exist (maybe due to sending aflow.in via email); leads to problems with string addition
      }
    }
    if(LDEBUG){
      cerr << soliloquy << "functional determined from aflow.in: " << functional << endl;
      cerr << soliloquy << "PBE: " << aurostd::withinList(cce_vars.vfunctionals, "PBE") << endl;
      cerr << soliloquy << "LDA: " << aurostd::withinList(cce_vars.vfunctionals, "LDA") << endl;
      cerr << soliloquy << "SCAN: " << aurostd::withinList(cce_vars.vfunctionals, "SCAN") << endl;
      cerr << soliloquy << "PBE+U_ICSD: " << aurostd::withinList(cce_vars.vfunctionals, "PBE+U_ICSD") << endl;
      cerr << soliloquy << "exp: " << aurostd::withinList(cce_vars.vfunctionals, "exp") << endl;
    }
    // if functional determined from aflow.in is different from the ones given by the input options, 
    // throw warning that oxidation numbers are only determined on the basis of a specific functional
    // also when only the formation enthalpy from the exprimental values per bond shall be obtained,
    // a warning should be given from which functional the data are used to determine the ox nums.
    if(cce_flags.flag("COMMAND_LINE")){
      oss << endl;
      if(!(cce_vars.vfunctionals.size() == 1 && cce_vars.vfunctionals[0] == functional)){
        oss << "WARNING: The oxidation numbers are only determined on the basis of a" << (functional == "LDA"?"n ":" ") << functional << " calculation." << endl;
      }
    }
  }

  //get_Bader_charges_from_Bader_file////////////////////////////////////////////////////////
  // determine the Bader charges from the Bader file
  vector<double> get_Bader_charges_from_Bader_file(const xstructure& structure, CCE_Variables& cce_vars, const string& Bader_file) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::get_Bader_charges_from_Bader_file():";
    vector<string> vlines; 
    aurostd::efile2vectorstring(Bader_file,vlines);
    vector<string> tokens;
    for (uint i = 1; i < vlines.size(); i++) {
      aurostd::string2tokens(vlines[i], tokens, " ");
      cce_vars.Bader_charges.push_back(aurostd::string2utype<double>(tokens[2]));
    }
    for (uint k=0,ksize=structure.atoms.size();k<ksize;k++){ // there should always be as many Bader charges as atoms in the structure
      if(LDEBUG){
        cerr << soliloquy << "Bader_charges[" << k << "]: " << cce_vars.Bader_charges[k] << endl;
      }
    }
    return cce_vars.Bader_charges;
  }

  //Bader_charges_to_oxidation_states////////////////////////////////////////////////////////
  // compare Bader charges to template values from fitting set and set oxidation numbers accordingly
  vector<double> Bader_charges_to_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, string& functional) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::Bader_charges_to_oxidation_states():";
    stringstream message;
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname != cce_vars.anion_species){
        string Bader_templ_line;
        if (get_Bader_templates(structure.atoms[i].cleanname) == "") {
          cce_flags.flag("CORRECTABLE",FALSE);
          cerr << "VERY BAD NEWS: There is no correction for " << structure.atoms[i].cleanname << " (ATOM[" << i << "])" << " since this species was not included in the set for deducing corrections!"  << endl;
          cerr << endl;
        } else {
          Bader_templ_line=get_Bader_templates(structure.atoms[i].cleanname);
          if(LDEBUG){
            cerr << soliloquy << "Bader templates: " << Bader_templ_line << endl;
          }
          vector<string> Bader_tokens;
          aurostd::string2tokens(Bader_templ_line, Bader_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. Bader_tokens[1] is not a space but the oxidation number
          uint num_ox_states = aurostd::string2utype<uint>(Bader_tokens[0]);
          if(LDEBUG){
            cerr << soliloquy << "number of oxidation states for which Bader charges are available: " << num_ox_states << endl;
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
            if (get_offset(functional) == -1) {
              message << "Unknown functional " << functional << ". Please choose PBE, LDA, or SCAN.";
              throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,message,_INPUT_ILLEGAL_);
            }
            int offset=get_offset(functional); // define also local offset variable since offset must correspond to functional detected from aflow.in
            for (uint i = 0; i < vfunctional_aflow_in.size(); i++) {
              Bader_template = aurostd::string2utype<double>(Bader_tokens[offset/2+2+5*n]);
              if(LDEBUG){
                cerr << soliloquy << "Bader_template " << vfunctional_aflow_in[i] << ": " << Bader_template << endl;
              }
            }
            if ( std::abs(Bader_template-cce_vars.Bader_charges[i]) < Bader_tolerance && std::abs(Bader_template-cce_vars.Bader_charges[i]) < Bader_deviation ){ // must be compared to Bader_deviation to load correction for Bader charge closest to template and not to be overloaded by other corrections for later tested Bader charges (oxidation states)
              Bader_deviation= std::abs(Bader_template-cce_vars.Bader_charges[i]);
              if(LDEBUG){
                cerr << soliloquy << "Bader_deviation: " << Bader_deviation << endl;
              }
              cce_vars.oxidation_states[i]= aurostd::string2utype<double>(Bader_tokens[1+5*n]);
            } else { // only if element is found but no corrections can be assigned because above conditions are not met, update Bader_deviation_min
              if ( std::abs(Bader_template-cce_vars.Bader_charges[i]) < Bader_deviation_min ) {
                Bader_deviation_min= std::abs(Bader_template-cce_vars.Bader_charges[i]);
              }
            }
          }
          if ( Bader_deviation == Bader_deviation_0 ){ // only if no oxidation state was found, i. e. Bader deviation was not changed but there are corrections (Bader charges) for this species (get_Bader_templates returns non-empty string), display the following error message after handling the special cases
            // FIRST point to deal with special cases known from (ternary) oxides; other cases will be dealt with when looking at the oxidation numbers and checking the sum
            // since it can happen that the Bader charge for W in the compound (especially for W+6) is too far away 
	    // from the Bader charge template, here W+6 is set and later the sum over the oxidation states will still be checked
            if (structure.atoms[i].cleanname == "W") {
              cce_vars.oxidation_states[i]=+6;
            // since it can happen that the Bader charge for Pb in the compound (especially for Pb+2) is too far away 
	    // from the Bader charge template, here Pb+2 is set and later the sum over the oxidation states will still be checked
            } else if (structure.atoms[i].cleanname == "Pb") {
              cce_vars.oxidation_states[i]=+2;
            // since it can happen that the Bader charge for Ag in the compound (especially for Ag+1) is too far away 
            // from the Bader charge template, here Ag+1 is set and later the sum over the oxidation states will still be checked
            } else if (structure.atoms[i].cleanname == "Ag") {
              cce_vars.oxidation_states[i]=+1;
            } else {
              cce_flags.flag("CORRECTABLE",FALSE);
              cerr << "BAD NEWS: The oxidation number (and hence the correction) for " << structure.atoms[i].cleanname << " (ATOM[" << i << "])" << " cannot be identified from the Bader charges!"  << endl;
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
              cerr << "Corrections for " << structure.atoms[i].cleanname << " coordinated by " << cce_vars.anion_species << " are available for oxidation states: " << ox_nums_avail << endl;
              cerr << "If the desired oxidation state is listed but it is just not correctly determined from the Bader charges," << endl;
              cerr << "you might want to consider supplying the oxidation numbers manually by using the option --oxidation_numbers=." << endl;
              cerr << endl;
            }
          }
        }
      } else if (structure.atoms[i].cleanname == cce_vars.anion_species) {
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

  //general_attempt_fixing_oxidation_states////////////////////////////////////////////////////////
  // for some ions known to be problematic for the Bader analysis (e.g. V+5, Fe+2, Fe+3, Ti+4) 
  // a brute force ansatz can be implemented to just change their oxidation state and later check 
  // whether it fixes the oxidation number sum rule
  void general_attempt_fixing_oxidation_states(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, ostream& oss) {
    if(cce_flags.flag("COMMAND_LINE")){
      oss << "The sum over all oxidation numbers for all atoms of the system is NOT zero, trying to repair that based on known problematic cases (Ti, V, Fe). This may or may not work." << endl;
    }
    // repairing only by considering non mixed valence oxides
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ 
      if (structure.atoms[i].cleanname == "Ti"){
        if ( cce_vars.oxidation_states[i] == 3){
          cce_vars.oxidation_states[i]=+4;
          if(cce_flags.flag("COMMAND_LINE")){
            oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Ti+4 " << endl;
          }
        }
      }
      if (structure.atoms[i].cleanname == "Fe"){
        if ( cce_vars.oxidation_states[i] == 2){
          cce_vars.oxidation_states[i]=+3;
          if(cce_flags.flag("COMMAND_LINE")){
            oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+3 " << endl;
          }
        } else if ( cce_vars.oxidation_states[i] == 3){
          cce_vars.oxidation_states[i]=+2;
          if(cce_flags.flag("COMMAND_LINE")){
            oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to Fe+2 " << endl;
          }
        }
      }
      if (structure.atoms[i].cleanname == "V"){
        if ( cce_vars.oxidation_states[i] == 4){
          cce_vars.oxidation_states[i]=+5;
          if(cce_flags.flag("COMMAND_LINE")){
            oss << "setting oxidation state of " << structure.atoms[i].cleanname << " (atom[" << i << "]) to V+5 " << endl;
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
 
  //get_corrections////////////////////////////////////////////////////////
  // determine the corrections after oxidation numbers are known (from input or Bader)
  void get_corrections(const xstructure& structure, xoption& cce_flags, CCE_Variables& cce_vars, const string& considered_anion_species, const vector<uint>& num_neighbors, vector<vector<double> >& corrections_atom) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::get_corrections():";
    if(LDEBUG){
      cerr << soliloquy << "ASSIGNMENT OF CORRECTIONS:" << endl;
    }
    for(uint i=0,isize=structure.atoms.size();i<isize;i++){ //loop over all atoms in structure
      if (structure.atoms[i].cleanname != cce_vars.anion_species && cce_vars.multi_anion_atoms[i] != 1){ // exclude main anion species and multi anion atoms detected previously
        string corrections_line;
        if (num_neighbors[i] > 0){ // are there actually bonds between the cation and the anion species under consideration (to load only needed corrections in multi-anion systems)
          if ( get_corrections_line(structure.atoms[i].cleanname + "_+" + aurostd::utype2string<double>(cce_vars.oxidation_states[i]) + "_" + considered_anion_species) == "") { // the considered anion species can be the main anion species or a multi anion species
            cce_flags.flag("CORRECTABLE",FALSE);
            cerr << "BAD NEWS: No correction available for " << structure.atoms[i].cleanname << " (ATOM[" << i << "])" << " in oxidation state " << cce_vars.oxidation_states[i] << " when coordinated by " << considered_anion_species << "." << endl;
            //checking for which oxidation states corrections are available and throw out errors accordingly
            uint ox_nums_count=0;
            vector<uint> ox_nums_avail_vec;
            for(uint k=0,ksize=12;k<ksize;k++){ // larger than ox. num. +12 should not occur
              if ( get_corrections_line(structure.atoms[i].cleanname + "_+" + aurostd::utype2string<uint>(k) + "_" + considered_anion_species) != "") {
                ox_nums_count+=1;
                ox_nums_avail_vec.push_back(k);
              }
            }
            if ( ox_nums_count == 0) {
              cerr << "Currently no corrections available for " << structure.atoms[i].cleanname << " when coordinated by "<< considered_anion_species << "." << endl;
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
              cerr << "Corrections for " << structure.atoms[i].cleanname << " coordinated by " << considered_anion_species << " are available for oxidation states: " << ox_nums_avail << endl;
              cerr << "If the desired oxidation state is listed but it is just not correctly determined," << endl;
              cerr << "you might want to consider supplying the oxidation numbers manually by using the option --oxidation_numbers=." << endl;
            }
            cerr << endl;
          } else {
            corrections_line=get_corrections_line(structure.atoms[i].cleanname + "_+" + aurostd::utype2string<double>(cce_vars.oxidation_states[i]) + "_" + considered_anion_species);
            if(LDEBUG){
              cerr << soliloquy << "corrections line: " << corrections_line << endl;
            }
            // load cation corrections
            load_cation_corrections(structure, cce_vars, corrections_line, corrections_atom, i);
          }
        }
      } else if (structure.atoms[i].cleanname == cce_vars.anion_species || cce_vars.multi_anion_atoms[i] == 1) {
        // set anion corrections (to zero)
        set_anion_corrections(structure, cce_vars, corrections_atom, i);
      }
    }
  }

  //load_cation_corrections////////////////////////////////////////////////////////
  // load corrections for the cations from the corrections line according to the functional
  // here for every functional there will be a separate "if" evaluation, since when only a structure (+ oxidation numbers) 
  // are provided as input, one might want to get the corrections for more than one functional
  void load_cation_corrections(const xstructure& structure, CCE_Variables& cce_vars, const string& corrections_line, vector<vector<double> >& corrections_atom, uint i) { // also provide increment for loop of get_corrections function
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::load_cation_corrections():";
    vector<string> corrections_tokens;
    aurostd::string2tokens(corrections_line, corrections_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. corrections_tokens[1] is not a space but the oxidation number
    uint num_funcs=cce_vars.vfunctionals.size();
    for (uint k = 0; k < num_funcs; k++) {
      corrections_atom[2*k][i]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+2]); // 2*k since to each functional belong 2 corrections for 298.15 and 0K
      if(LDEBUG){
        cerr << soliloquy << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].cleanname << " (atom[" << i << "]) for 298.15K: " << corrections_atom[2*k][i] << " eV/bond" << endl;
      }
      if (cce_vars.vfunctionals[k] != "exp") {
        corrections_atom[2*k+1][i]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+3]);
        if(LDEBUG){
          cerr << soliloquy << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].cleanname << " (atom[" << i << "]) for 0K: " << corrections_atom[2*k+1][i] << " eV/bond" << endl;
        }
      }
    }
  }

  //set_anion_corrections////////////////////////////////////////////////////////
  // set corrections for anion atoms to zero since only for cations number of bonds with anions are needed; 
  // per- & superoxides will be dealt with in other function
  void set_anion_corrections(const xstructure& structure, CCE_Variables& cce_vars, vector<vector<double> >& corrections_atom, uint i) { // also provide increment for loop of get_corrections function
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::set_anion_corrections():";
    uint num_funcs=cce_vars.vfunctionals.size();
    for (uint k = 0; k < num_funcs; k++) {
      corrections_atom[2*k][i]=0; // 2*k since to each functional belong 2 corrections for 298.15 and 0K
      if(LDEBUG){
        cerr << soliloquy << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].cleanname << " (atom[" << i << "]) for 298.15K: " << corrections_atom[2*k][i] << " eV/bond" << endl;
      }
      if (cce_vars.vfunctionals[k] != "exp") {
        corrections_atom[2*k+1][i]=0;
        if(LDEBUG){
          cerr << soliloquy << cce_vars.vfunctionals[k] << " correction for " << structure.atoms[i].cleanname << " (atom[" << i << "]) for 0K: " << corrections_atom[2*k+1][i] << " eV/bond" << endl;
        }
      }
    }
  }

  //check_get_per_super_ox_corrections////////////////////////////////////////////////////////
  // check whether corrections for per- or superoxide ions are needed and assign accordingly
  void check_get_per_super_ox_corrections(CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::check_get_per_super_ox_corrections():";
    // check for per- and superoxides should always be kept separate since a compound could have both per- and superoxide ions (alkali-metal sesquioxides)
    if (cce_vars.num_perox_bonds > 0){
      string corrections_line;
      corrections_line=get_corrections_line("O2_-2_O");
      if(LDEBUG){
        cerr << soliloquy << "corrections line peroxides: " << corrections_line << endl;
      }
      vector<string> corrections_tokens;
      aurostd::string2tokens(corrections_line, corrections_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. corrections_tokens[1] is not a space but the oxidation number
      uint num_funcs=cce_vars.vfunctionals.size();
      for (uint k = 0; k < num_funcs; k++) {
        cce_vars.perox_correction[2*k]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+2]);
        if(LDEBUG){
          cerr << soliloquy << cce_vars.vfunctionals[k] << " peroxide correction for 298.15K: " << cce_vars.perox_correction[2*k] << " eV/bond" << endl;
        }
        if (cce_vars.vfunctionals[k] != "exp") {
          cce_vars.perox_correction[2*k+1]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+3]);
          if(LDEBUG){
            cerr << soliloquy << cce_vars.vfunctionals[k] << " peroxide correction for 0K: " << cce_vars.perox_correction[2*k+1] << " eV/bond" << endl;
          }
        }
      }
    } 
    if (cce_vars.num_superox_bonds > 0) {
      string corrections_line;
      corrections_line=get_corrections_line("O2_-1_O");
      if(LDEBUG){
        cerr << soliloquy << "corrections line superoxides: " << corrections_line << endl;
      }
      vector<string> corrections_tokens;
      aurostd::string2tokens(corrections_line, corrections_tokens, " "); // seems to automatically reduce the number of multiple spaces in a row to one so that e.g. corrections_tokens[1] is not a space but the oxidation number
      uint num_funcs=cce_vars.vfunctionals.size();
      for (uint k = 0; k < num_funcs; k++) {
        cce_vars.superox_correction[2*k]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+2]);
        if(LDEBUG){
          cerr << soliloquy << cce_vars.vfunctionals[k] << " superoxide correction for 298.15K: " << cce_vars.superox_correction[2*k] << " eV/bond" << endl;
        }
        if (cce_vars.vfunctionals[k] != "exp") {
          cce_vars.superox_correction[2*k+1]= aurostd::string2utype<double>(corrections_tokens[cce_vars.offset[k]+3]);
          if(LDEBUG){
            cerr << soliloquy << cce_vars.vfunctionals[k] << " superoxide correction for 0K: " << cce_vars.superox_correction[2*k+1] << " eV/bond" << endl;
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
 
  //check_apply_per_super_ox_corrections////////////////////////////////////////////////////////
  // check whether corrections for per- or superoxide ions are needed and apply accordingly
  void check_apply_per_super_ox_corrections(CCE_Variables& cce_vars) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || CCE_DEBUG);
    string soliloquy="cce::check_apply_per_super_ox_corrections():";
    if (cce_vars.num_perox_bonds > 0){
      uint num_funcs=cce_vars.vfunctionals.size();
      for (uint k = 0; k < num_funcs; k++) {
        cce_vars.cce_correction[2*k] += (cce_vars.num_perox_bonds * cce_vars.perox_correction[2*k]) ;
        if(LDEBUG){
          cerr << soliloquy << cce_vars.vfunctionals[k] << " peroxide correction for 298.15K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_perox_bonds * cce_vars.perox_correction[2*k]) << " eV" << endl;
        }
        if (cce_vars.vfunctionals[k] != "exp") {
          cce_vars.cce_correction[2*k+1] += (cce_vars.num_perox_bonds * cce_vars.perox_correction[2*k+1]) ;
          if(LDEBUG){
            cerr << soliloquy << cce_vars.vfunctionals[k] << " peroxide correction for 0K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_perox_bonds * cce_vars.perox_correction[2*k+1]) << " eV" << endl;
            cerr << endl;
          }
        }
      }
    }
    if (cce_vars.num_superox_bonds > 0){
      uint num_funcs=cce_vars.vfunctionals.size();
      for (uint k = 0; k < num_funcs; k++) {
        cce_vars.cce_correction[2*k] += (cce_vars.num_superox_bonds * cce_vars.superox_correction[2*k]) ;
        if(LDEBUG){
          cerr << soliloquy << cce_vars.vfunctionals[k] << " superoxide correction for 298.15K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_superox_bonds * cce_vars.superox_correction[2*k]) << " eV" << endl;
        }
        if (cce_vars.vfunctionals[k] != "exp") {
          cce_vars.cce_correction[2*k+1] += (cce_vars.num_superox_bonds * cce_vars.superox_correction[2*k+1]) ;
          if(LDEBUG){
            cerr << soliloquy << cce_vars.vfunctionals[k] << " superoxide correction for 0K per cell: " << std::setprecision(3) << std::fixed << (cce_vars.num_superox_bonds * cce_vars.superox_correction[2*k+1]) << " eV" << endl;
            cerr << endl;
          }
        }
      }
    }
  }

  ////get_formation_enthalpies//////////////////////////////////////////////
  //// ME20200213
  //// Returns the formation enthalpy per cell for each functional
  //vector<double> get_formation_enthalpies(const vector<double>& cce_correction, CCE_Variables& cce_vars) {
  //  vector<double> formation_enthalpies(cce_correction.size(), 0.0);
  //  uint num_funcs=cce_vars.vfunctionals.size();
  //  for (uint i = 0; i < num_funcs; i++) {
  //    if (cce_vars.vfunctionals[i] != "exp" ) {
  //      if (cce_vars.dft_energies.size() > 0) {
  //        // for 298.15 K
  //        formation_enthalpies[2*i] = cce_vars.dft_energies[i] - cce_correction[2*i];
  //        // for 0 K
  //        formation_enthalpies[2*i+1] = cce_vars.dft_energies[i] - cce_correction[2*i+1];
  //      }
  //    } else { // treat exp. separately since formation enthalpy is equal to sum of corrections in this case
  //      formation_enthalpies[2*i] = cce_correction[2*i];
  //    }
  //  }
  //  return formation_enthalpies;
  //}

  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                        WRITE OUTPUT AND CITATION                        //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
 
  //get_JSON/////////////////////////////////////////////////////////////
  // ME20200213
  // Returns CCE results in JSON format
  string get_JSON(const xstructure& structure, const CCE_Variables& cce_vars) {
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
    json << "},";
    json << "\"publication\":";
    json << "\"Friedrich et al., Coordination corrected ab initio formation enthalpies, npj Comput. Mater. 5, 59 (2019). ";
    json << "https://doi.org/10.1038/s41524-019-0192-1\"";
    json << "}";
    return json.str();
  }

  //write_output////////////////////////////////////////////////////////
  // write CCE corrections and corrected formation enthalpies if precalculated DFT values are provided
  string write_output(const xstructure& structure, CCE_Variables& cce_vars, const vector<double>& cce_form_energy_cell) {
    stringstream output;
    // print out CCE corrections per cell and atom for functionals selected
    if (!(cce_vars.vfunctionals.size() == 1 && cce_vars.vfunctionals[0] == "exp")){ // if only exp is set as functional CCE CORRECTIONS: should not be written
      output << "CCE CORRECTIONS (to be subtracted from precalculated DFT formation energies):" << endl;
    }
    uint num_funcs=cce_vars.vfunctionals.size();
    for (uint k = 0; k < num_funcs; k++) {
      if (cce_vars.vfunctionals[k] != "exp") {
        output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k] << " eV/cell //" << "CCE@" << cce_vars.vfunctionals[k] << " correction for 298.15K." << endl;
        output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k+1] << " eV/cell //" << "CCE@" << cce_vars.vfunctionals[k] << " correction for 0K." << endl;
        output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k]/structure.atoms.size() << " eV/atom //" << "CCE@" << cce_vars.vfunctionals[k] << " correction for 298.15K." << endl;
        output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k+1]/structure.atoms.size() << " eV/atom //" << "CCE@" << cce_vars.vfunctionals[k] << " correction for 0K." << endl;
        output << endl;
      }
    }
    // exp result should always be written at the end, hence write only after writing output for other functionals
    for (uint k = 0; k < num_funcs; k++) {
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
      uint num_funcs=cce_vars.vfunctionals.size();
      for (uint k = 0; k < num_funcs; k++) {
        if (cce_vars.vfunctionals[k] != "exp") {
          output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k] << " eV/cell //" << "CCE@" << cce_vars.vfunctionals[k] << " formation enthalpy at 298.15K." << endl;
          output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k+1] << " eV/cell //" << "CCE@" << cce_vars.vfunctionals[k] << " formation enthalpy at 0K." << endl;
          output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k]/structure.atoms.size() << " eV/atom //" << "CCE@" << cce_vars.vfunctionals[k] << " formation enthalpy at 298.15K." << endl;
          output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k+1]/structure.atoms.size() << " eV/atom //" << "CCE@" << cce_vars.vfunctionals[k] << " formation enthalpy at 0K." << endl;
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

  //write_test_output////////////////////////////////////////////////////////
  // write CCE corrections and corrected formation enthalpies for testing
  string write_test_output(CCE_Variables& cce_vars, const vector<double>& cce_form_energy_cell) {
    stringstream output;
    uint num_funcs=cce_vars.vfunctionals.size();
    for (uint k = 0; k < num_funcs; k++) {
      if (cce_vars.vfunctionals[k] != "exp") {
        if (k<num_funcs-1 || cce_vars.dft_energies.size()!=0) {
          output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k] << ",";
          output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k+1] << ",";
        } else {
          output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k] << ",";
          output << std::setprecision(3) << std::fixed << cce_vars.cce_correction[2*k+1];
        }
      }
    }
    // exp part
    for (uint k = 0; k < num_funcs; k++) {
      if (cce_vars.vfunctionals[k] == "exp" && cce_vars.dft_energies.size()==0) { // second condition for that if precalc. form. energies are given and asking for exp., exp. result is not written twice
        output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k];
      }
    }
    // corrected DFT formation energies
    if(cce_vars.dft_energies.size()!=0){
      uint num_funcs=cce_vars.vfunctionals.size();
      for (uint k = 0; k < num_funcs; k++) {
        if (cce_vars.vfunctionals[k] != "exp") {
          if (k<num_funcs-1) {
            output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k] << ",";
            output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k+1] << ",";
          } else {
            output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k] << ",";
            output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k+1];
          }
        } else if (cce_vars.vfunctionals[k] == "exp") {
          if (k<num_funcs-1) {
            output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k] << ",";
          } else {
            output << std::setprecision(3) << std::fixed << cce_form_energy_cell[2*k];
          }
        }
      }
    }
    return output.str();
  }

  //write_citation////////////////////////////////////////////////////////
  // write citation information at end of output
  string write_citation() {
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

  //print_usage////////////////////////////////////////////////////////
  // For printing user instructions if no additional input is provided.
  string print_usage() {
    stringstream oss;
    oss << endl;
    oss << "Written by Rico Friedrich, Corey Oses, and Marco Esters, 2018-2020" << endl;
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
    oss << "--cce/--cce=POSCAR --usage  Prints these user instructions." << endl;
    oss << "--cce=POSCAR_FILE_PATH      Provide the path to the structure file in VASP POSCAR format." << endl;
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
    oss << "                            be in VASP POSCAR format)." << endl;
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
