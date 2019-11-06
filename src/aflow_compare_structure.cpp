// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                                                                         *
// ***************************************************************************
// AFLOW-XTAL-MATCH (compare crystal structures)
// Written by David Hicks (david.hicks@duke.edu) 
// Contributors: Carlo De Santo

#include<fstream>
#include<iostream>
#include<vector>
#include<string>
#include<exception>
#include<algorithm>
#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_compare_structure.h"
#include "aflow_symmetry_spacegroup.h"

#undef AFLOW_COMPARE_MULTITHREADS_ENABLE

#if GCC_VERSION >= 40400   // added two zeros
#define AFLOW_COMPARE_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of AFLOW-COMPARE will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

// ***************************************************************************
// AFLOW COMMANDS 
// ***************************************************************************

// ======== Single Comparison ======== //

  // aflow --compare_material=POSCAR_1,POSCAR_2
  //       Description: Compares POSCAR_1 to POSCAR_2 if the ratios are 
  //                    commensurate and types of atoms are the same 
  //                    (i.e. same material). 

  // aflow --compare_structure=POSCAR_1,POSCAR_2
  //       Description: Compares POSCAR_1 to POSCAR_2 (no requirement on type 
  //                    of atoms, just stoichiometry ratios). 

// ======== Directory Comparison ======== //

  // aflow --compare_material_directory|--compare_material_dir
  //       Description: Determines the unique materials (same atomic species) 
  //                     within a given directory. Returns a JSON and TXT file 
  //                     with the results. Default Directory: "."

  // aflow --compare_structure_directory|--compare_structure_dir
  //       Description: Determines the unique structure prototypes within a 
  //                    given directory. Returns a JSON and TXT file with the 
  //                    results. Default Directory: "."

  // OPTIONS:
  //       -D "PATH":   User can specify a specific directory to compare. 
  //                    Output will be placed there also.

// ***************************************************************************
// OPTIONS FOR ALL AFLOW COMMANDS
// ***************************************************************************
// --np=xx          : Number of processors. Algorithm is thread-friendly. 
//                    (Default: 8 processors)

// ***************************************************************************
// OVERVIEW OF ALGORITHM
// ***************************************************************************
// This algorithm takes two crystal structures and compares them to one another
// and determines their level of similarity as described by H. Burzlaff 
// (see his paper for more information: Acta Cryst., A53, 217-224 (1997)).

// Steps:
//   1) Scale volumes of two structures to be commensurate
//   2) Determine the least frequently occuring atom (LFA)
//   3) Shift the least frequently occuring atom (LFA) for 
//      both structures; these will be used to indicate our lattice
//   4) Create a 3x3x3 supercell of structure2
//   5) Search for all possible quadruplets of LFA atoms in the 
//      supercell and see if we can match it with structure1 
//      using Burzlaff's criteria
//   6) Once possible quadruplet/lattice is found, check contents
//      of lattice (i.e. atoms). Check we can have a one-to-one 
//      mapping
//   7) Of the best match (smallest misfit (mis)):
//      If mis <= 0.1:
//        Structures similar. Print out new representation of 
//        structure2 and the figure of misfit
//      else if 0.1 < mis <=0.2:  
//        Structures in the same family (possible symmetric
//        group-subgroup relation). Print out new representation 
//        of structure2 and the figure of misfit
//      else mis>0.2:
//        Structures are not the same. Print "No match"

// ***************************************************************************

// ***************************************************************************
// pflow::CompareStructures - Prepares comparison from command line input 
// ***************************************************************************
namespace pflow {
  string compareStructures(aurostd::xoption& vpflow){ 
    bool LDEBUG=(false || XHOST.DEBUG);
    ostringstream oss;
    ofstream FileMESSAGE; //DX 20190319 - added FileMESSAGE
    ostream& logstream = cout; //DX 20190424
    stringstream message; //DX 20190424
    string directory = "";

    string function_name = "pflow::compareStuctures()";
   
    string usage_material="aflow --compare_material=POSCAR1,POSCAR2";
    string usage_structure="aflow --compare_structure=POSCAR1,POSCAR2";
    //[THREADS]string options_single="--np=xx (default 8),--print";
    string options_single="--print"; //fast
 
    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //Defalut=1
    if(vpflow.flag("COMPARE_STRUCTURE::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_STRUCTURE::NP"));
    }

    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool same_species=false;
    // material comparisons (find duplicates)
    if(vpflow.flag("COMPARE_MATERIAL")){
      same_species=true;
    }
    // structure comparisons (find structure prototypes)
    else if (vpflow.flag("COMPARE_STRUCTURE")){ 
      same_species=false;
    }

    // ---------------------------------------------------------------------------
    // FLAG: grab input structures to compare
    vector<string> vinput;
    //DX 20190425 [OBSOLETE]    if(vpflow.flag("COMPARE_STRUCTURE::STRUCTURES_1")) {
    //DX 20190425 [OBSOLETE]      vinput.push_back(vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURES_1"));
    //DX 20190425 [OBSOLETE]    }
    //DX 20190425 [OBSOLETE]    if(vpflow.flag("COMPARE_STRUCTURE::STRUCTURES_2")) {
    //DX 20190425 [OBSOLETE]      vinput.push_back(vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURES_2"));
    //DX 20190425 [OBSOLETE]    }
    if(vpflow.flag("COMPARE_STRUCTURE::STRUCTURE_LIST")){ //DX 20190425
      aurostd::string2tokens(vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURE_LIST"),vinput,","); //DX 20190425
    } //DX 20190425
 
    // ---------------------------------------------------------------------------
    // FLAG: consider magnetic structure in comparison 
    bool magnetic_comparison=false;
    vector<string> magmoms_for_systems;
    if(vpflow.flag("COMPARE_STRUCTURE::MAGNETIC")){
      magnetic_comparison=true;
      string magnetic_info=vpflow.getattachedscheme("COMPARE_STRUCTURE::MAGNETIC");
      aurostd::string2tokens(magnetic_info,magmoms_for_systems,":");
      message << "OPTIONS: Including magnetic moment information in comparisons. Magnetic input detected for " << magmoms_for_systems.size() << " systems.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
 
    // ---------------------------------------------------------------------------
    // FLAG: print mapping details
    bool print=false;
    if(vpflow.flag("COMPARE_STRUCTURE::PRINT")) {
      print=true;
    }
    
    // ---------------------------------------------------------------------------
    // optimize match (default: false)
    bool optimize_match=false;
    if(vpflow.flag("COMPARE_STRUCTURE::OPTIMIZE_MATCH")) {
      optimize_match=true;
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: scale volume
    bool scale_volume=true;
    if(vpflow.flag("COMPARE_STRUCTURE::NO_SCALE_VOLUME")) {
      scale_volume=false;
    }

    //DX 20190424 - START 
    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    bool ignore_Wyckoff=false;
    if(vpflow.flag("COMPARE_STRUCTURE::IGNORE_WYCKOFF")) {
      ignore_Wyckoff=true;
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    bool ignore_symmetry=false;
    if(vpflow.flag("COMPARE_STRUCTURE::IGNORE_SYMMETRY")) {
      ignore_symmetry=true;
      ignore_Wyckoff=true;
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      }
    
    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    bool ignore_environment=false;
    if(vpflow.flag("COMPARE_STRUCTURE::IGNORE_ENVIRONMENT_ANALYSIS")) {
      ignore_environment=true;
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: remove duplicate compounds (useful for non-biased statistics)
    bool remove_duplicate_compounds=false;
    if(vpflow.flag("COMPARE_STRUCTURE::REMOVE_DUPLICATE_COMPOUNDS")) {
      remove_duplicate_compounds=true;
      message << "OPTIONS: Remove duplicate compounds first, useful for non-biased prototype statistics."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      }
    
    // ---------------------------------------------------------------------------
    // FLAG: order structures by smallest ICSD number if designation is in title 
    bool ICSD_comparison=false;
    if(vpflow.flag("COMPARE_STRUCTURE::ICSD_COMPARISON")) {
      ICSD_comparison=true;
      message << "OPTIONS: Running on ICSD structures; use oldest ICSD number as representative prototype.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    //DX 20190424 - END 

    // ---------------------------------------------------------------------------
    // FLAG: match unique structures to the AFLOW prototypes 
    bool add_matching_aflow_protos=false;
    if(vpflow.flag("COMPARE_STRUCTURE::MATCH_TO_AFLOW_PROTOS")) {
      add_matching_aflow_protos=true;
      message << "OPTIONS: Compare unique structures to the AFLOW prototypes."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: get AFLOW ANRL designation for unique structures
    bool get_aflow_prototype_designation=false;
    if(vpflow.flag("COMPARE_STRUCTURE::ADD_AFLOW_PROTOTYPE_DESIGNATION")) {
      get_aflow_prototype_designation=true;
      message << "OPTIONS: Cast unique structures into AFLOW standard designation."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: do not calculate unique permutations 
    bool calculate_unique_permutations=true;
    if(vpflow.flag("COMPARE_STRUCTURE::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS")) {
      calculate_unique_permutations=false;
      message << "OPTIONS: Do not calculate unique permutations."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // DX TODO: allow for comparing multiple inputs
    bool multiple_comparisons=false;
    bool single_comparison_round=false;
    bool clean_unmatched=true; //DX 20190504
    bool store_comparison_logs = false; //DX 20190822 - add log bool

    // ---------------------------------------------------------------------------
    // check input structures
    //DX 20190424 - START
    if(LDEBUG) cerr << "pflow::compareStructures: Number of structures: " << vinput.size() << endl;
    vector<StructurePrototype> all_structures = compare::loadStructuresFromStructureList(vinput, magmoms_for_systems, same_species,  FileMESSAGE);
    if(LDEBUG) cerr << "pflow::compareStructures: Number of successfully loaded structures: " << all_structures.size() << endl;

    // check if more than two structures, i.e., perform multiple comparisons
    if(all_structures.size()>2){ multiple_comparisons=true;}

    // check if at least two structures
    if(all_structures.size()<2){ 
      cerr << "pflow::compareStructures: ERROR: There must be at least two structures to compare." << endl;
      return oss.str();
    }

    //DX 20190424 - END

    //DX 20190425 [OBSOLETE] double final_misfit=-1.0;
    //DX 20190425 [OBSOLETE] if(LDEBUG) cerr << "pflow::compareStructures: begin" << endl;
    //DX 20190425 [OBSOLETE] if(LDEBUG) cerr << "pflow::compareStructures: vinput.size()=" << vinput.size() << endl;

    //DX 20190425 [OBSOLETE] // ensure there are two structures 
    //DX 20190425 [OBSOLETE] if(vinput.size()!=2) {      cerr << "pflow::compareStructures: ERROR vinput.size()!=2" << endl;
    //DX 20190425 [OBSOLETE]   if(vpflow.flag("COMPARE_MATERIAL")) {
    //DX 20190425 [OBSOLETE]     init::ErrorOption(cout,options_single,"pflow::compareStructures()",usage_material);
    //DX 20190425 [OBSOLETE]   }
    //DX 20190425 [OBSOLETE]   else if(vpflow.flag("COMPARE_STRUCTURE")) {
    //DX 20190425 [OBSOLETE]     init::ErrorOption(cout,options_single,"pflow::compareStructures()",usage_structure);
    //DX 20190425 [OBSOLETE]   }
    //DX 20190425 [OBSOLETE]   return oss.str();
    //DX 20190425 [OBSOLETE] }

    //DX 20190425 [OBSOLETE] if(LDEBUG) cerr << "pflow::compareStructures: loading vinput.at(0)=" << vinput.at(0) << endl;

    //DX 20190425 [OBSOLETE] // check and load structure 1 
    //DX 20190425 [OBSOLETE] if(!aurostd::FileExist(vinput.at(0))){
    //DX 20190425 [OBSOLETE]   oss << "pflow::compareStructures: ERROR file vinput.at(0)=" << vinput.at(0) << "  not found" << endl;
    //DX 20190425 [OBSOLETE]   return oss.str();
    //DX 20190425 [OBSOLETE] }
    //DX 20190425 [OBSOLETE] stringstream sss1;
    //DX 20190425 [OBSOLETE] aurostd::efile2stringstream(vinput.at(0),sss1);
    //DX 20190425 [OBSOLETE] xstructure xstr1(sss1);  

    //DX 20190425 [OBSOLETE] // check and load structure 2 
    //DX 20190425 [OBSOLETE] if(LDEBUG) cerr << "pflow::compareStructures: loading vinput.at(1)=" << vinput.at(1) << endl;
    //DX 20190425 [OBSOLETE] if(!aurostd::FileExist(vinput.at(1))){
    //DX 20190425 [OBSOLETE]   cerr << "pflow::compareStructures: ERROR file vinput.at(1)=" << vinput.at(1) << "  not found" << endl;
    //DX 20190425 [OBSOLETE]   return oss.str();
    //DX 20190425 [OBSOLETE] }
    //DX 20190425 [OBSOLETE] stringstream sss2;
    //DX 20190425 [OBSOLETE] aurostd::efile2stringstream(vinput.at(1),sss2);
    //DX 20190425 [OBSOLETE] xstructure xstr2(sss2);  

    // ---------------------------------------------------------------------------
    // compare structures
    if(!multiple_comparisons){
      double final_misfit=-1.0; //DX 20190424
      store_comparison_logs = true; //DX 20190822 - add log bool
      // call main comparison function
      // DX 20190424 [OBSOLETE] compare::aflowCompareStructure(num_proc,xstr1,xstr2,same_species, scale_volume, optimize_match, oss,final_misfit);
      compare::aflowCompareStructure(num_proc,all_structures[0].representative_structure,all_structures[1].representative_structure,same_species, scale_volume, optimize_match, oss,final_misfit); //DX 2010424
      if(print==true){
        // return mapping details
        return oss.str();
      }
      else {
        // return abbreviated results (i.e., misfit value along with match, same family, or no match text
        oss.str("");
        oss.clear();
        oss << final_misfit << " : ";
        if(final_misfit <=0.1 && (final_misfit+1.0)> 1e-3){
          oss << "MATCH" << endl;
        }
        else if(final_misfit > 0.1 && final_misfit <= 0.2){
          oss << "SAME FAMILY" << endl;
        }
        else if(final_misfit > 0.2 || (final_misfit+1.0) < 1e-3){ 
          oss << "NOT A MATCH" << endl;
        }
      }
    }
    else {
      compare::compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, remove_duplicate_compounds, calculate_unique_permutations, add_matching_aflow_protos, get_aflow_prototype_designation, ICSD_comparison, store_comparison_logs); //DX 20190504 -added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
    }
    return oss.str();
  }
}

// ***************************************************************************
// pflow::comparePermutations()
// ***************************************************************************
namespace pflow {
  string comparePermutations(istream& input, aurostd::xoption& vpflow){
    ostringstream oss;
    ofstream FileMESSAGE; //DX 20190319 - added FileMESSAGE
  
    string usage="aflow --compare_permutation<POSCAR";
    string options="[--usage] [--np=<number>] [--print_misfit]";

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE_PERMUTATION::USAGE")) {
      stringstream ss_usage;
      init::ErrorOption(ss_usage,vpflow.getattachedscheme("COMPARE_PERMUTATION"),"pflow::comparePermutations()",aurostd::liststring2string(usage,options));
      return ss_usage.str();
    }

    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE_PERMUTATION::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_PERMUTATION::NP"));
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: print misfit values for duplicate permutations 
    bool print_misfit=false; //defalut=false
    if(vpflow.flag("COMPARE_PERMUTATION::PRINT")) { print_misfit=true; }
    
    //DX 20190504 - START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    string format = "both";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      format = "txt";
    }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      format = "json";
    }
    //DX 20190504 - END
    
    // ---------------------------------------------------------------------------
    // FLAG: optimize match (default: false)
    // bool optimize_match=false; // permutation comparisons do not need to have a best match; let's save time
    // ---------------------------------------------------------------------------
    // optimize match (default: false)
    bool optimize_match=false;
    if(vpflow.flag("COMPARE_PERMUTATION::OPTIMIZE_MATCH")) {
      optimize_match=true;
    }

    // ---------------------------------------------------------------------------
    // load structure
    xstructure xstr(input,IOAFLOW_AUTO);
    
    // ---------------------------------------------------------------------------
    // calculate unique/duplicate permutations 
    vector<string> unique_permutations = compare::getUniquePermutations(xstr, num_proc, optimize_match, print_misfit, oss, FileMESSAGE); //DX 20190319 - added FileMESSAGE

    return oss.str();
  }
}

// ***************************************************************************
// compare::getUniquePermutations()
// ***************************************************************************
namespace compare{
  vector<string> getUniquePermutations(xstructure& xstr){
    uint num_proc=1;
    return compare::getUniquePermutations(xstr, num_proc);
  }
}

namespace compare{
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc){
    bool optimize_match=false;
    return compare::getUniquePermutations(xstr, num_proc, optimize_match);
  }
}

namespace compare{
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc, bool& optimize_match){
    bool print_misfit=false;
    ostringstream oss;
    ofstream FileMESSAGE; //DX 20190319 - added FileMESSAGE
    return compare::getUniquePermutations(xstr, num_proc, optimize_match, print_misfit, oss, FileMESSAGE); //DX 20190319 - added FileMESSAGE
  }
}

namespace compare{
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc, bool& optimize_match, bool& print_misfit, ostream& oss, ofstream& FileMESSAGE){ //DX 20190319 - added FileMESSAGE
    
    vector<string> unique_permutations;
    stringstream ss_output; //DX 20190506

    //DX 20190506 - START
    // ---------------------------------------------------------------------------
    // print format 
    string format = "text";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      format = "text";
    }
    if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      format = "json";
    }
    //DX 20190506 - END

    // ---------------------------------------------------------------------------
    // permutation comparisons must compare the same species 
    bool same_species=true; 

    // ---------------------------------------------------------------------------
    // quick check: check if any sites have the same number of atoms; if not, then no need to try comparing
    if(!print_misfit){
      if(!arePermutationsComparableViaStoichiometry(xstr)){ //DX 20190624 - put into function
        vector<uint> reduced_stoichiometry = gcdStoich(xstr.num_each_type); //DX 20190508
        unique_permutations = generatePermutationString(reduced_stoichiometry); //DX 20190508
        if(format=="text"){ //DX 20190506
          ss_output << "Unique permutations (" << unique_permutations.size() << "): " << endl; 
          ss_output << " " << aurostd::joinWDelimiter(unique_permutations,"\n ") << endl;
        }
        if(format=="json"){ //DX 20190506
          ss_output << "{\"unique_permutations\":["; 
          ss_output << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(unique_permutations,"\""),",") << "]";
          ss_output << "]}" << endl;
        }
        oss << ss_output.str();
        return unique_permutations;
      }
    }
    
    // ---------------------------------------------------------------------------
    // load input structure
    StructurePrototype structure;
    structure.representative_structure = xstr;
    structure.representative_structure_name = "input geometry";
    structure.stoichiometry = compare::getStoichiometry(xstr,true);
    structure.elements = compare::getElements(xstr);
    // update xstructure species
    if(structure.representative_structure.species.size()==0){
      deque<string> deque_species; for(uint j=0;j<structure.elements.size();j++){deque_species.push_back(structure.elements[j]);}
      structure.representative_structure.SetSpecies(deque_species);
      structure.representative_structure.SpeciesPutAlphabetic();
    }
    structure.representative_structure_generated = true; 
    structure.representative_structure_from = "input"; 

    // ---------------------------------------------------------------------------
    // get the unique permutations for the structure
    vector<StructurePrototype> final_permutations = compare::comparePermutations(structure,num_proc,optimize_match,oss,FileMESSAGE); //DX 20190319 - added FileMESSAGE

    // ---------------------------------------------------------------------------
    // print results
    if(format=="text"){ //DX 20190506
      ss_output << "Unique permutations (" << final_permutations.size() << "): " << endl; 

      for(uint j=0;j<final_permutations.size();j++){
        ss_output << " " << final_permutations[j].representative_structure_name;
        for (uint k=0;k<final_permutations[j].duplicate_structures_names.size();k++){
          ss_output << " = " << final_permutations[j].duplicate_structures_names[k];
        }
        ss_output << endl;
      }
    }
    //DX 20190506 - START
    else if(format=="json"){
      stringstream sscontent_json;
      vector<string> vcontent_json;
      sscontent_json << "\"unique_permutations\":[";
      for(uint j=0;j<final_permutations.size();j++){
        stringstream sstmp;
        vector<string> equivalent_permutations;
        equivalent_permutations.push_back(final_permutations[j].representative_structure_name);
        equivalent_permutations.insert(equivalent_permutations.end(),final_permutations[j].duplicate_structures_names.begin(),final_permutations[j].duplicate_structures_names.end());
        sstmp << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(equivalent_permutations,"\""),",") << "]";
        vcontent_json.push_back(sstmp.str()); sstmp.str("");
      }
      sscontent_json << aurostd::joinWDelimiter(vcontent_json,",");
      sscontent_json << "]";
      ss_output << "{" << sscontent_json.str() << "}" << endl;
    }
    //DX 20190506 - END

    // ---------------------------------------------------------------------------
    // print misfit results
    if(print_misfit){
      if(format=="text"){ //DX 20190506
        ss_output << "Misfit values: " << endl; 
        stringstream ss_text;
        compare::printResults(ss_text, same_species, final_permutations, "text");
        ss_output << ss_text.str();
      }
      else if(format=="json"){ //DX 20190506
        ss_output.str(""); // need to clear content abbreviated content from above
        stringstream ss_json;
        compare::printResults(ss_json, same_species, final_permutations, "json");
        ss_output << ss_json.str() << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // store unique permutations in vector
    for(uint j=0;j<final_permutations.size();j++){
      unique_permutations.push_back(final_permutations[j].representative_structure_name);
    }

    // update oss
    oss << ss_output.str();
    return unique_permutations;
  }
}

// ***************************************************************************
// pflow::compareMultipleStructures()
// ***************************************************************************
namespace pflow {
  string compareMultipleStructures(aurostd::xoption& vpflow){ //DX 20190425 - changed name, more general
    
    // This function compares multiple structures (i.e., more than two).

    string function_name = "pflow::compareMultipleStructures()";
    bool LDEBUG=(false || XHOST.DEBUG);
    ostringstream oss;
    ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE_STRUCTURE::USAGE")) {
      stringstream ss_usage;
      // material-type comparisons
      if(vpflow.flag("COMPARE_MATERIAL_DIRECTORY")){
        string usage_material_comparison="aflow --compare_materials -D <dir_path>";
        string options_material_comparison="[--usage] [--np=|--num_proc=<number>] [--optimize_match] [--no_scale_volume] [--ignore_symmetry] [--ignore_Wyckoff]";
        init::ErrorOption(ss_usage,vpflow.getattachedscheme("COMPARE_STRUCTURE"),"pflow::compareMultipleStructures()",aurostd::liststring2string(usage_material_comparison,options_material_comparison));
      }
      // structure-type comparisons
      else if(vpflow.flag("COMPARE_STRUCTURE_DIRECTORY")){
        string usage_structure_comparison="aflow --compare_structures -D <dir_path>";
        string options_structure_comparison="[--usage] [--np=|--num_proc=<number>] [--optimize_match] [--no_scale_volume] [--ignore_symmetry] [--ignore_Wyckoff] [--remove_duplicates|--remove_duplicate_compounds]";
        init::ErrorOption(ss_usage,vpflow.getattachedscheme("COMPARE_STRUCTURE"),"pflow::compareMultipleStructures()",aurostd::liststring2string(usage_structure_comparison,options_structure_comparison));
      }
      return ss_usage.str();
    }
    
    // ---------------------------------------------------------------------------
    // distinguish structures coming from directory or file
    string structures_from = ""; // "structure_list", "directory" or "file"

    // from list appended to command
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURE_LIST").empty()){
      structures_from = "structure_list";
    }
    // from directory
    //DX 20190424 [OBSOLETE] if(vpflow.flag("COMPARE_MATERIAL_DIRECTORY") || vpflow.flag("COMPARE_STRUCTURE_DIRECTORY")){
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::DIRECTORY").empty()){
      structures_from = "directory";
    }
    // from file
    //DX 20190424 [OBSOLETE] if(vpflow.flag("COMPARE_MATERIAL_FILE") || vpflow.flag("COMPARE_STRUCTURE_FILE")){
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::FILE").empty()){
      structures_from = "file";
    }

    // ---------------------------------------------------------------------------
    // FLAG: directory of structures to compare
    vector<string> file_list; //DX 20190424
    string directory=".";
    string filename="";
    //DX 20190424 - START
    if(structures_from=="structure_list") {
      aurostd::string2tokens(vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURE_LIST"),file_list,",");
      message << "List of files to compare: " << aurostd::joinWDelimiter(file_list,",");
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    //DX 20190424 - END
    else if(structures_from=="directory") {
      directory=vpflow.getattachedscheme("COMPARE_STRUCTURE::DIRECTORY");
      if(!aurostd::FileExist(directory)) {
        message << "Unable to locate directory: " << directory << ". Exiting." << endl;
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
        return oss.str();
      }
      message << "Comparison directory: " << directory;
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    // ---------------------------------------------------------------------------
    // FLAG: file of structures to compare
    else if(structures_from=="file") {
      filename=vpflow.getattachedscheme("COMPARE_STRUCTURE::FILE");
      if(!aurostd::FileExist(filename)) {
        message << "Unable to locate file: " << filename << ". Exiting." << endl;
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
        return oss.str();
      }
      message << "Comparison file: " << filename;
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    else {
      message << "Need to specify location of structures to compare: -D <directory> or -F=<filename>. Exiting." << endl;
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE_STRUCTURE::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_STRUCTURE::NP"));
      message << "OPTIONS: Using multiple threads; np = " << num_proc << ".";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
  
    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool same_species=false;
    if(vpflow.flag("COMPARE_MATERIAL")){  // material comparisons (find duplicates) //DX 20190429 - generalized
      same_species=true;
      message << "OPTIONS: Performing material type comparisons (comparing alike atomic species).";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    else if (vpflow.flag("COMPARE_STRUCTURE")){ // structure comparisons (find structure prototypes) //DX 20190425 - generalized
      same_species=false;
      message << "OPTIONS: Performing structure type comparisons (any atomic species).";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
   
    // ---------------------------------------------------------------------------
    // FLAG: consider magnetic structure in comparison 
    bool magnetic_comparison=false;
    vector<string> magmoms_for_systems;
    if(vpflow.flag("COMPARE_STRUCTURE::MAGNETIC")){
      magnetic_comparison=true;
      string magnetic_info=vpflow.getattachedscheme("COMPARE_STRUCTURE::MAGNETIC");
      aurostd::string2tokens(magnetic_info,magmoms_for_systems,":");
      message << "OPTIONS: Including magnetic moment information in comparisons. Magnetic input detected for " << magmoms_for_systems.size() << " systems.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
 
    // ---------------------------------------------------------------------------
    // FLAG: ICSD comparison - structure with minimum ICSD number as representative prototype
    // in general: smaller ICSD number = older = more reliable
    bool ICSD_comparison=false;
    if(vpflow.flag("COMPARE_STRUCTURE::ICSD_COMPARISON")) {
      ICSD_comparison=true;
      message << "OPTIONS: Running on ICSD structures; use oldest ICSD number as representative prototype.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // optimize match (default: false)
    bool optimize_match=false;
    if(vpflow.flag("COMPARE_STRUCTURE::OPTIMIZE_MATCH")) {
      optimize_match=true;
      message << "OPTIONS: Finding optimal match; exploring all possible lattices and origins to find the best match (note: this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    bool scale_volume=true;
    if(vpflow.flag("COMPARE_STRUCTURE::NO_SCALE_VOLUME")) {
      scale_volume=false;
      message << "OPTIONS: Suppressing volume scaling; useful for distinguishing structures at different pressures."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    bool ignore_Wyckoff=false;
    if(vpflow.flag("COMPARE_STRUCTURE::IGNORE_WYCKOFF")) {
      ignore_Wyckoff=true;
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    bool ignore_symmetry=false;
    if(vpflow.flag("COMPARE_STRUCTURE::IGNORE_SYMMETRY")) {
      ignore_symmetry=true;
      ignore_Wyckoff=true;
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    bool ignore_environment=false;
    if(vpflow.flag("COMPARE_STRUCTURE::IGNORE_ENVIRONMENT_ANALYSIS")) {
      ignore_environment=true;
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: remove duplicate compounds (useful for non-biased statistics)
    bool remove_duplicate_compounds=false;
    if(vpflow.flag("COMPARE_STRUCTURE::REMOVE_DUPLICATE_COMPOUNDS")) {
      remove_duplicate_compounds=true;
      message << "OPTIONS: Remove duplicate compounds first, useful for non-biased prototype statistics."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
   
    // ---------------------------------------------------------------------------
    // FLAG: match unique structures to the AFLOW prototypes 
    bool add_matching_aflow_protos=false;
    if(vpflow.flag("COMPARE_STRUCTURE::MATCH_TO_AFLOW_PROTOS")) {
      add_matching_aflow_protos=true;
      message << "OPTIONS: Compare unique structures to the AFLOW prototypes."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: get AFLOW ANRL designation for unique structures
    bool get_aflow_prototype_designation=false;
    if(vpflow.flag("COMPARE_STRUCTURE::ADD_AFLOW_PROTOTYPE_DESIGNATION")) {
      get_aflow_prototype_designation=true;
      message << "OPTIONS: Cast unique structures into AFLOW standard designation."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: do not calculate unique permutations 
    bool calculate_unique_permutations=true;
    if(vpflow.flag("COMPARE_STRUCTURE::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS")) {
      calculate_unique_permutations=false;
      message << "OPTIONS: Do not calculate unique permutations."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
   
    //DX 20190425 - added print and screen only flag - START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    string format = "both";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      format = "txt";
    }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      format = "json";
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: print format
    bool screen_only = false;
    if(vpflow.flag("COMPARE_STRUCTURE::SCREEN_ONLY")) {
      screen_only=true;
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: print mapping details (for two-structure comparison only)
    bool print=false;
    if(vpflow.flag("COMPARE_STRUCTURE::PRINT")) {
      print=true;
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: do not remove unmatched structures from the StructurePrototype Object
    // keeps results of each comparison
    bool clean_unmatched=true;
    if(vpflow.flag("COMPARE_STRUCTURE::KEEP_UNMATCHED")) {
      clean_unmatched=false;
    }
    
    // ---------------------------------------------------------------------------
    // single round of comparisons
    bool single_comparison_round = false; //compare all structures until matched or exhausted all comparisons

    bool store_comparison_logs = false; //DX 20190822 - add log bool

    // ---------------------------------------------------------------------------
    // check if two-structure comparison
    if(file_list.size()==2){
      single_comparison_round = true;
      clean_unmatched = false;
      store_comparison_logs = true; //DX 20190822 - add log bool
    }
    //DX 20190425 - added print and screen only flag - END

    //BETA // ===== FLAG: REMOVE DUPLICATE COMPOUNDS LAST ===== //
    //BETA bool remove_duplicate_compounds_last=false;
    //BETA if(vpflow.flag("COMPARE_STRUCTURE::REMOVE_DUPLICATE_COMPOUNDS_LAST")) {
    //BETA   remove_duplicate_compounds_last=true;
    //BETA   message << "OPTIONS: Remove duplicate compounds after performing structure-type comparisons, useful for non-biased prototype statistics."; 
    //BETA   pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    //BETA }
    
    // ---------------------------------------------------------------------------
    // load structures 
    vector<StructurePrototype> final_prototypes;
    if(structures_from=="structure_list") {
      final_prototypes = compare::compareStructuresFromStructureList(file_list, magmoms_for_systems, oss, FileMESSAGE, num_proc, same_species, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, remove_duplicate_compounds, calculate_unique_permutations, add_matching_aflow_protos, get_aflow_prototype_designation, ICSD_comparison, store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
    }
    else if(structures_from=="directory") {
      final_prototypes = compare::compareStructuresFromDirectory(directory, magmoms_for_systems, oss, FileMESSAGE, num_proc, same_species, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, remove_duplicate_compounds, calculate_unique_permutations, add_matching_aflow_protos, get_aflow_prototype_designation, ICSD_comparison, store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
      //all_structures = compare::loadStructuresFromDirectory(directory, same_species);
    }
    if(structures_from=="file") {
      final_prototypes = compare::compareStructuresFromFile(filename, magmoms_for_systems, oss, FileMESSAGE, num_proc, same_species, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, remove_duplicate_compounds, calculate_unique_permutations, add_matching_aflow_protos, get_aflow_prototype_designation, ICSD_comparison, store_comparison_logs);  //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
      //all_structures = compare::loadStructuresFromDirectory(directory, same_species);
    }

    // ---------------------------------------------------------------------------
    // prepare JSON output
    stringstream ss_json;
    compare::printResults(ss_json, same_species, final_prototypes, "json");
  
    // ---------------------------------------------------------------------------
    // prepare TEXT (.out) output
    stringstream ss_out;
    compare::printResults(ss_out, same_species, final_prototypes, "txt");
  
    //DX 20190429 - added screen only option - START
    // ---------------------------------------------------------------------------
    // write results to screen and return immediately (do not write file)
    if(screen_only){
      if(format=="json"){ return ss_json.str(); }
      // default is txt
      else { return ss_out.str(); }
    }
    //DX 20190429 - added screen only option - END

    //DX 20190429 - added format options - START
    // ---------------------------------------------------------------------------
    // if only two comparisons and text only, print mismatch information 
    if(file_list.size()==2){
      // return abbreviated results (i.e., misfit value along with match, same family, or no match text
      double final_misfit = -1.0;
      if(final_prototypes[0].misfits.size()==1){
        final_misfit =  final_prototypes[0].misfits[0];
      }
      message << final_misfit << " : ";
      if(final_misfit <=0.1 && (final_misfit+1.0)> 1e-3){
        message << "MATCH" << endl;
      }
      else if(final_misfit > 0.1 && final_misfit <= 0.2){
        message << "SAME FAMILY" << endl;
      }
      else if(final_misfit > 0.2 || (final_misfit+1.0) < 1e-3){ 
        message << "NOT A MATCH" << endl;
      }
      if(XHOST.QUIET){
        oss << message.str();
      }
      else {
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
      }
      if(print){
        if(final_prototypes[0].misfits.size()==1){
          oss << final_prototypes[0].duplicate_comparison_logs[0];
        }
      }
      return oss.str();
    }
  
    // ---------------------------------------------------------------------------
    // write results to files
    if(same_species==true){
      if(format=="json"){
        aurostd::stringstream2file(ss_json,directory+"/material_comparison_output.json");
        message << "RESULTS: See " << directory << "/material_comparison_output.json" << " for list of unique/duplicate materials.";
      }
      else if(format=="txt"){
        aurostd::stringstream2file(ss_out,directory+"/material_comparison_output.out");
        message << "RESULTS: See " << directory << "/material_comparison_output.out" << " for list of unique/duplicate materials.";
      }
      else if(format=="both"){
        aurostd::stringstream2file(ss_json,directory+"/material_comparison_output.json");
        aurostd::stringstream2file(ss_out,directory+"/material_comparison_output.out");
        message << "RESULTS: See " << directory << "/material_comparison_output.out" << " or " << directory << "/material_comparison_output.json" << " for list of unique/duplicate materials.";
      }
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }
    else if(same_species==false){
      if(format=="json"){
        aurostd::stringstream2file(ss_json,directory+"/structure_comparison_output.json");
        message << "RESULTS: See " << directory << "/structure_comparison_output.json" << " for list of unique/duplicate structures.";
      }
      else if(format=="txt"){
        aurostd::stringstream2file(ss_out,directory+"/structure_comparison_output.out");
        message << "RESULTS: See " << directory << "/structure_comparison_output.out" << " for list of unique/duplicate structures.";
      }
      else if(format=="both"){
        aurostd::stringstream2file(ss_json,directory+"/structure_comparison_output.json");
        aurostd::stringstream2file(ss_out,directory+"/structure_comparison_output.out");
        message << "RESULTS: See " << directory << "/structure_comparison_output.out" << " or " << directory << "/structure_comparison_output.json" << " for list of unique/duplicate structures.";
      }
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }
    //DX 20190429 - added format options - END

    return oss.str();
  }
}

//DX 20190314 - added new function - START
// ***************************************************************************
// pflow::getMatchingPrototype - returns corresponding prototype label
// ***************************************************************************
namespace pflow {
  vector<string> getMatchingPrototypes(xstructure& xstr, string& catalog){ 
    
    // Returns the matching prototype label, if any exists

    aurostd::xoption vpflow;
    vpflow.flag("COMPARE2PROTOTYPES",TRUE);
    
    // ---------------------------------------------------------------------------
    // specify catalog
    vpflow.flag("COMPARE2PROTOTYPES::CATALOG",TRUE); //DX 20190329 - need to make scheme before attaching, otherwise it doesn't work
    vpflow.push_attached("COMPARE2PROTOTYPES::CATALOG",catalog); 
    
    // ---------------------------------------------------------------------------
    // do not calculate unique permutations
    vpflow.flag("COMPARE2PROTOTYPES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS",TRUE);
   
    // ---------------------------------------------------------------------------
    // quiet output
    bool original_quiet = XHOST.QUIET;
    XHOST.QUIET=TRUE;

    // ---------------------------------------------------------------------------
    // compare structure to AFLOW prototypes 
    vector<StructurePrototype> prototypes = compare2prototypes(xstr,vpflow);

    // ---------------------------------------------------------------------------
    // global quiet back to default
    XHOST.QUIET=original_quiet;

    return prototypes[0].duplicate_structures_names; // duplicates names are prototype labels 
  }
}
//DX 20190314 - added new function - START

//DX 20190314 - added overloads for compare2prototypes - START
// ***************************************************************************
// pflow::compare2prototypes - identifies corresponding protos 
// ***************************************************************************
namespace pflow {
  vector<StructurePrototype> compare2prototypes(istream& input, aurostd::xoption& vpflow){ 
    
    // ---------------------------------------------------------------------------
    // load input structure
    xstructure xstr(input,IOAFLOW_AUTO);
  
    return compare2prototypes(xstr, vpflow);
  }
}

// ***************************************************************************
// pflow::printMatchingPrototypes - returns list of matching structures 
// ***************************************************************************
namespace pflow {
  string printMatchingPrototypes(istream& input, aurostd::xoption& vpflow){ 
    
    // ---------------------------------------------------------------------------
    // load input structure
    xstructure xstr(input,IOAFLOW_AUTO);

    return printMatchingPrototypes(xstr,vpflow);
  }
}

// ***************************************************************************
// pflow::printMatchingPrototypes - returns list of matching structures 
// ***************************************************************************
namespace pflow {
  string printMatchingPrototypes(xstructure& xstr, aurostd::xoption& vpflow){ 
    
    //DX 20190425 - added print flag - START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    string format = "both";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      format = "txt";
    }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      format = "json";
    }
    //DX 20190425 - added print flag - END
    
    //DX 20190425 - added screen only flag - START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    bool screen_only = false;
    if(vpflow.flag("COMPARE2PROTOTYPES::SCREEN_ONLY")) {
      screen_only=true;
    }
    //DX 20190425 - added screen only flag - END
    
    vector<StructurePrototype> prototypes = compare2prototypes(xstr,vpflow);
    
    // ---------------------------------------------------------------------------
    // print results 
    //DX 20190509 [OBSOLETE-moved down] stringstream ss_out;
    bool same_species = false; //default for prototypes
    stringstream ss_json;
    compare::printResults(ss_json, same_species, prototypes, "json");
    
    stringstream ss_out;
    compare::printResults(ss_out, same_species, prototypes, "txt");
    
    //DX 20190429 - added screen only option - START
    // ---------------------------------------------------------------------------
    // write results to screen and return immediately (do not write file)
    if(screen_only){
      if(format=="json"){ return ss_json.str(); }
      // default is txt
      else { return ss_out.str(); }
    }
    //DX 20190429 - added screen only option - END
  
    if(format=="json"){ return ss_json.str(); }
    return ss_out.str();
  }
}

//DX 20190314 - added overloads for compare2prototypes - END
// ***************************************************************************
// pflow::compare2prototypes - identifies corresponding protos 
// ***************************************************************************
namespace pflow {
  vector<StructurePrototype> compare2prototypes(xstructure& xstr, aurostd::xoption& vpflow){ 
    bool LDEBUG=(false || XHOST.DEBUG);
    
    string function_name = "pflow::compare2prototypes()";
    ostringstream oss;
    ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;

    string directory="";

    string usage="aflow --compare2protos|--compare2prototypes < POSCAR";
    string options="";

    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE2PROTOTYPES::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE2PROTOTYPES::NP"));
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: catalog (htqc, anrl, or all)
    string catalog="all";
    if(vpflow.flag("COMPARE2PROTOTYPES::CATALOG")) {
      catalog=aurostd::tolower(vpflow.getattachedscheme("COMPARE2PROTOTYPES::CATALOG"));
      if(catalog!="htqc" && catalog!="anrl" && catalog!="all"){
        message << "Catalog/library can only be htqc, anrl, or all.";     
        throw aurostd::xerror(function_name,message,_INPUT_ILLEGAL_); //DX 20191031 - exit to xerror
      }
      message << "OPTIONS: Catalog/library (htqc, anrl, or all): " << catalog << endl; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // same spaces 
    bool same_species=false; //compare2prototype: by definition, want to compare backbone structure, i.e., ignore species
    
    // ---------------------------------------------------------------------------
    // optimize match (default: false)
    bool optimize_match=false; //compare2prototype: list similar structures; no need for best match
    
    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    bool scale_volume=true;
    if(vpflow.flag("COMPARE2PROTOTYPES::NO_SCALE_VOLUME")) {
      scale_volume=false;
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    bool ignore_Wyckoff=false;
    if(vpflow.flag("COMPARE2PROTOTYPES::IGNORE_WYCKOFF")) {
      ignore_Wyckoff=true;
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    bool ignore_symmetry=false;
    if(vpflow.flag("COMPARE2PROTOTYPES::IGNORE_SYMMETRY")) {
      ignore_symmetry=true;
      ignore_Wyckoff=true;
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    bool ignore_environment=false;
    if(vpflow.flag("COMPARE2PROTOTYPES::IGNORE_ENVIRONMENT_ANALYSIS")) {
      ignore_environment=true;
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: do not calculate unique permutations 
    bool calculate_unique_permutations=true;
    if(vpflow.flag("COMPARE2PROTOTYPES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS")) {
      calculate_unique_permutations=false;
      message << "OPTIONS: Do not calculate unique permutations."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ===== Single round of comparisons ===== //
    bool single_comparison_round = true; //perform only first comparison
    bool clean_unmatched = true; // clean unmatched structures from object //DX 20190506
    bool store_comparison_logs = false; //DX 20190822 - add log bool

    vector<StructurePrototype> all_structures;
   
    //DX 20190314 [OBSOLETE] // ---------------------------------------------------------------------------
    //DX 20190314 [OBSOLETE] // load input structure
    //DX 20190314 [OBSOLETE] xstructure xstr(input,IOAFLOW_AUTO);

    StructurePrototype input_structure;
    input_structure.representative_structure = xstr;
    input_structure.representative_structure_name = "input geometry";
    input_structure.stoichiometry = compare::getStoichiometry(xstr,true); //true preserves the stoich order for the structure
    input_structure.elements = compare::getElements(xstr);
    input_structure.representative_structure_compound = compare::getCompoundName(xstr);
    input_structure.representative_structure_generated = true;
    stringstream ss_input; ss_input << xstr;
    input_structure.representative_structure_from = ss_input.str(); 
    all_structures.push_back(input_structure);

    // ---------------------------------------------------------------------------
    // symmetry
    if(!ignore_symmetry && xstr.space_group_ITC==0){ //DX 20190829 - don't recalculate symmetry if already calcualted
      message << "Calculating the symmetry of the input structure.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      uint one_proc=1;
      compare::calculateSymmetries(all_structures,one_proc); 
      message << "Symmetry calculated.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }
    else if(!ignore_symmetry && xstr.space_group_ITC!=0){
      for(uint i=0;i<all_structures.size();i++){
        all_structures[i].space_group = all_structures[i].representative_structure.space_group_ITC;
        vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
        compare::groupWyckoffPositions(all_structures[i].representative_structure, grouped_Wyckoff_positions);
        all_structures[i].grouped_Wyckoff_positions=grouped_Wyckoff_positions;
      }
    }
    
    if(LDEBUG) {
      cerr << function_name << ": Wyckoff positions of input structure:" << endl;
      for(uint i=0;i<all_structures[0].grouped_Wyckoff_positions.size();i++){
        cerr << all_structures[0].grouped_Wyckoff_positions[i] << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // variables to pass to GetPrototypes functions
    vector<uint> stoichiometry = all_structures[0].stoichiometry;
	  std::sort(stoichiometry.begin(),stoichiometry.end()); // order stoichiometry, so we can match to AFLOW prototypes more easily
    uint space_group_num = all_structures[0].space_group;
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions = all_structures[0].grouped_Wyckoff_positions;
    
    if(LDEBUG) {
      cerr << function_name << ": Wyckoff positions of input structure:" << endl;
      for(uint i=0;i<all_structures[0].grouped_Wyckoff_positions.size();i++){
        cerr << all_structures[0].grouped_Wyckoff_positions[i] << endl;
      }
    }

    vector<string> vlabel;
    vector<uint> prototype_space_groups;

    // ---------------------------------------------------------------------------
    // find prototypes based on stoichiometry only
    if(ignore_symmetry){
      message << "Load prototypes with the same stoichiometry as the input.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      vlabel = aflowlib::GetPrototypesByStoichiometry(stoichiometry, catalog);
    }
    // find prototypes based on stoichiometry and space group only
    else if(!ignore_symmetry && ignore_Wyckoff){
      message << "Load prototypes with the same stoichiometry and space group as the input.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      vector<GroupedWyckoffPosition> empty_Wyckoff_positions;
      vlabel = aflowlib::GetPrototypesBySymmetry(stoichiometry, space_group_num, empty_Wyckoff_positions, prototype_space_groups, 2, catalog);
    }
    // find prototypes based on stoichiometry, space group, and Wyckoff positions only (recommended/default)
    else {
      message << "Load prototypes with the same stoichiometry, space group, and Wyckoff positions as the input.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      vlabel = aflowlib::GetPrototypesBySymmetry(stoichiometry, space_group_num, grouped_Wyckoff_positions, prototype_space_groups, 2, catalog);
    }
    message << "Potential compatible prototypes: " << vlabel.size() << " (" << aurostd::joinWDelimiter(vlabel,",") << ").";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    if(LDEBUG) {
      cerr << function_name << ": Wyckoff positions of input structure:" << endl;
      for(uint i=0;i<all_structures[0].grouped_Wyckoff_positions.size();i++){
        cerr << all_structures[0].grouped_Wyckoff_positions[i] << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // load compatible aflow prototypes
      
    // DX 20190830 - to avoid multiple threads being spun-up (here and in aflow_xproto.cpp), turn of aflow_pthreads
    uint uint_backup=AFLOW_PTHREADS::MAX_PTHREADS;
    AFLOW_PTHREADS::MAX_PTHREADS=1;

    compare::addAFLOWPrototypes2StructurePrototypeVector(all_structures, vlabel); 
    if(LDEBUG) {
      cerr << function_name << ": Wyckoff positions of input structure:" << endl;
      for(uint i=0;i<all_structures[0].grouped_Wyckoff_positions.size();i++){
        cerr << all_structures[0].grouped_Wyckoff_positions[i] << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // group into objects based on stoichiometry and symmetry (Pearson and space group)
    message << "Grouping sets of comparisons.";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> comparison_schemes = compare::groupStructurePrototypes(all_structures, same_species, ignore_symmetry, ignore_Wyckoff, ignore_environment, false); //DX 20190731 - add ignore_environment //DX 20190829 - false for duplicates_removed 
    message << "Number of comparison groups: " << comparison_schemes.size() << ".";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
   
    // ---------------------------------------------------------------------------
    // compare structures 
    message << "Running comparisons ...";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes = compare::runComparisonScheme(num_proc, comparison_schemes, same_species, false, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, false, store_comparison_logs, oss, FileMESSAGE);  //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190731 - added ignore_symmetry/Wyckoff/environment //DX 20190822 - add log bool //DX 20190930 - false=duplicates_removed

    message << "Comparisons complete ...";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    
    AFLOW_PTHREADS::MAX_PTHREADS = uint_backup; // DX 20190830 - set back to original setting
    
    // ---------------------------------------------------------------------------
    // return if there are no similar structures
    if(final_prototypes.size()==0){ return final_prototypes; } // DX 20190314 originally : return oss.str()

    comparison_schemes.clear();

    // ---------------------------------------------------------------------------
    // get unique permutations 
    if(calculate_unique_permutations){
      message << "Identifying unique permutations for representative structures ...";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

      for(uint i=0;i<final_prototypes.size();i++){
        // check if xstructure is generated; if not, make it
        if(!final_prototypes[i].representative_structure_generated){
          if(!compare::generateStructure(final_prototypes[i].representative_structure_name,final_prototypes[i].representative_structure_from,final_prototypes[i].representative_structure,oss)){
            message << "Could not generate structure (" << final_prototypes[i].representative_structure_name << ").";
            throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_); //DX 20191031 - exit to xerror
          }
        }        
        if(LDEBUG){ //DX 20190601 - added LDEBUG
          cerr << "Finding unique permutations for " << final_prototypes[i].representative_structure_name << ".";
        }        
        vector<StructurePrototype> final_permutations = compare::comparePermutations(final_prototypes[i],num_proc,optimize_match,oss,FileMESSAGE); //DX 20190319 - added FileMESSAGE
        for(uint j=0;j<final_permutations.size();j++){
          final_prototypes[i].unique_permutations.push_back(final_permutations[j].representative_structure_name);
        }
      }
      message << "Unique permutations found.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }

    // ---------------------------------------------------------------------------
    // print results 
    //DX 20190314 -moved to overload [OBSOLETE] stringstream ss_out;
    //DX 20190314 -moved to overload [OBSOLETE] compare::printResults(ss_out, same_species, final_prototypes);
    //DX 20190314 -moved to overload [OBSOLETE]
    //DX 20190314 -moved to overload [OBSOLETE] oss << ss_out.str();
    //DX 20190314 -moved to overload [OBSOLETE] return oss.str();
    return final_prototypes; //DX 20190314 - new return type
  }
}

// ***************************************************************************
// pflow::compare2database - compares database 
// ***************************************************************************
namespace pflow {
  string compare2database(istream& input, aurostd::xoption& vpflow){
    bool LDEBUG=(false || XHOST.DEBUG);
    
    string function_name = "pflow::compare2database()";
    string directory = "";
    ostringstream oss;
    ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;

    string usage="aflow --compare2database < POSCAR";
    string options="";

    vector<string> tokens,sub_tokens;
    vector<string> matchbook; //aflux - filter/get properties
    vector<string> schema; //aflux - get metadata (e.g., units)
    vector<string> property_units;

    bool same_species = true;
    bool single_comparison_round = true;
    bool store_comparison_logs = false; //DX 20190822 - add log bool

    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE2DATABASE::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE2DATABASE::NP"));
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: optimize match
    bool optimize_match=false; //false
    if(vpflow.flag("COMPARE2DATABASE::OPTIMIZE_MATCH")) {
      optimize_match=true;
      message << "OPTIONS: Finding optimal match; exploring all possible lattices and origins to find the best match (note: this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    bool scale_volume=true;
    if(vpflow.flag("COMPARE2DATABASE::NO_SCALE_VOLUME")) {
      scale_volume=false;
      message << "OPTIONS: Suppressing volume scaling; useful for distinguishing structures at different pressures."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    bool ignore_Wyckoff=false;
    if(vpflow.flag("COMPARE2DATABASE::IGNORE_WYCKOFF")) {
      ignore_Wyckoff=true;
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    bool ignore_symmetry=false;
    if(vpflow.flag("COMPARE2DATABASE::IGNORE_SYMMETRY")) {
      ignore_symmetry=true;
      ignore_Wyckoff=true;
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    bool ignore_environment=false;
    if(vpflow.flag("COMPARE2DATABASE::IGNORE_ENVIRONMENT_ANALYSIS")) {
      ignore_environment=true;
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool structure_comparison=false;
    if(vpflow.flag("COMPARE2DATABASE::STRUCTURE")) {
      structure_comparison=true;
      same_species = false;
      message << "OPTIONS: Structure-type comparison, i.e., ignore atomic species."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: property list to extract from database (using AFLUX)
    vector<string> property_list;
    if(vpflow.flag("COMPARE2DATABASE::PROPERTY_LIST")) {
      aurostd::string2tokens(vpflow.getattachedscheme("COMPARE2DATABASE::PROPERTY_LIST"),property_list,",");

      // put properties in schema and matchbook for AFLUX call
      schema.push_back("schema("+vpflow.getattachedscheme("COMPARE2DATABASE::PROPERTY_LIST")+")"); //to get units
      matchbook.insert(matchbook.end(), property_list.begin(), property_list.end());

      message << "OPTIONS: Extracting the following properties: " << aurostd::joinWDelimiter(property_list,", "); 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: specify the geometry file to grab (orig, relax1, relax2, static, bands, POSCAR, CONTCAR)
    // DX TODO 
    string geometry_file = "";
    if(vpflow.flag("COMPARE2DATABASE::GEOMETRY_FILE")) {
      geometry_file = vpflow.getattachedscheme("COMPARE2DATABASE::GEOMETRY_FILE");
      message << "OPTIONS: Structure type (POSCAR.orig, POSCAR.relax1, POSCAR.relax2, CONTCAR.relax1, ...): " << geometry_file << endl; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: catalog (icsd, lib1, lib2, lib3, ...)
    string catalog = "";
    string catalog_summons = "";
    if(vpflow.flag("COMPARE2DATABASE::CATALOG")) {
      catalog = aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::CATALOG")); //DX 20190329 -- added tolower
      if(catalog != "all"){ //DX 20190329 - added if-statement since AFLUX doesn't use "all"
        catalog_summons = "catalog(\'" + catalog + "\')";
        matchbook.push_back(catalog_summons);
      } //DX 20190329 - added if-statement since AFLUX doesn't use "all"
      message << "OPTIONS: Catalog/library (icsd, lib1, lib2, lib3, ...): " << catalog << endl; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    //DX 20190425 - added print flag - START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    string format = "both";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      format = "txt";
    }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      format = "json";
    }
    //DX 20190425 - added print flag - END
    
    //DX 20190508 - added keep unmatched - START
    // ---------------------------------------------------------------------------
    // FLAG: do not remove unmatched structures from the StructurePrototype Object
    // keeps results of each comparison
    bool clean_unmatched=true;
    if(vpflow.flag("COMPARE2DATABASE::KEEP_UNMATCHED")) {
      clean_unmatched=false;
    }
    //DX 20190508 - added keep unmatched - END
    
    //DX 20190425 - added screen only flag - START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    bool screen_only = false;
    if(vpflow.flag("COMPARE2DATABASE::SCREEN_ONLY")) {
      screen_only=true;
    }
    //DX 20190425 - added screen only flag - END
    
    // ---------------------------------------------------------------------------
    // load input structure
    xstructure xstr(input,IOAFLOW_AUTO);

    //DX 20190329 - added species check - START
    // check if fake names for same species comparison
    if(xstr.species[0]=="A" && !structure_comparison){
      message << "Atomic species are missing for the input structure. Cannot compare to database materials without species.";     
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
      return "";
    }
    //DX 20190329 - added species check - END

    // ---------------------------------------------------------------------------
    // calculate symmetry
    int space_group_number = xstr.SpaceGroup_ITC();
    
    // ---------------------------------------------------------------------------
    // get stoichiometries
    vector<uint> stoichiometry = compare::getStoichiometry(xstr,same_species);
    uint stoichiometry_sum = aurostd::sum(stoichiometry);
    vector<double> normalized_stoichiometry;
    for(uint i=0;i<stoichiometry.size();i++){normalized_stoichiometry.push_back((double)stoichiometry[i]/(double)stoichiometry_sum);}

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get entries with compatible space groups, i.e., same or enantiomorph
    if(!ignore_symmetry){
      string space_group_summons = "";
      // check if enantiomorphic space group
      int enantiomorph_space_group_number = SYM::getEnantiomorphSpaceGroupNumber(space_group_number);
      if(space_group_number == enantiomorph_space_group_number){
        // relaxed: need to match last in string, i.e., "*,<sg_symbol> <sg_number>" (comma necessary or we may grab the orig symmetry)
        space_group_summons = "sg2(*%27," + GetSpaceGroupName(space_group_number) + "%20%23" + aurostd::utype2string<int>(space_group_number) + "%27)";
      }
      else { // need to get enantiomorph too
        space_group_summons = "sg2(*%27," + GetSpaceGroupName(space_group_number) + "%20%23" + aurostd::utype2string<int>(space_group_number) + "%27";
        space_group_summons += ":*%27," + GetSpaceGroupName(enantiomorph_space_group_number) + "%20%23" + aurostd::utype2string<int>(enantiomorph_space_group_number) + "%27)";
      }
      matchbook.push_back(space_group_summons);
    }
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
    compare::groupWyckoffPositions(xstr, grouped_Wyckoff_positions);

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get aurl for entry
    string aurl = "aurl";
    matchbook.push_back(aurl);
    
    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get species and number of species
    string species_summons = "";
    if(!structure_comparison){
      species_summons = "species(" + aurostd::joinWDelimiter(xstr.species,",") + ")";
    }
    string nspecies_summons = "nspecies(" + aurostd::utype2string<uint>(xstr.num_each_type.size()) + ")";
    matchbook.push_back(species_summons);
    matchbook.push_back(nspecies_summons);

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: format AFLUX output 
    string aflux_format = "format(aflow)";
    string paging = "paging(0)";
    matchbook.push_back(aflux_format);
    matchbook.push_back(paging);

    // ---------------------------------------------------------------------------
    // construct aflux summons, i.e., combine matchbook
    string Summons = aurostd::joinWDelimiter(matchbook,",");
    message << "AFLUX matchbook request: " << Summons;
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // call AFLUX 
    string response = aflowlib::AFLUXCall(Summons);

    message << "Number of entries returned: " << aurostd::string2tokens(response,tokens,"\n");
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    if(LDEBUG) {cerr << function_name << "::AFLUX response:" << endl << response << endl;}
   
    // ---------------------------------------------------------------------------
    // extract properties from AFLUX response
    vector<vector<std::pair<string,string> > > properties_response = aflowlib::getPropertiesFromAFLUXResponse(response);
    if(LDEBUG) {
      for(uint i=0;i<properties_response.size();i++){
        for(uint j=0;j<properties_response[i].size();j++){
          cerr << properties_response[i][j].first << " = " << properties_response[i][j].second << ", ";
        }
        cerr << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // extract aurl, auid, and compound type from properties variable
    vector<string> auids, aurls, compounds;
		for(uint i=0;i<properties_response.size();i++){
			for(uint j=0;j<properties_response[i].size();j++){
				if(properties_response[i][j].first=="aurl"){
					aurls.push_back(properties_response[i][j].second);
				}
				if(properties_response[i][j].first=="auid"){
					auids.push_back(properties_response[i][j].second);
				}
				if(properties_response[i][j].first=="compound"){
					compounds.push_back(properties_response[i][j].second);
				}
			}
		}
    //cerr << "==============================" << endl;
    //::print(auids);
    //::print(aurls);
    //::print(compounds);
    
    // ---------------------------------------------------------------------------
    // get AFLUX schema, i.e., metadata (for the units)
    if(schema.size()>0){
      schema.push_back(aflux_format);
      schema.push_back(paging);

      // call AFLUX to get schema
      response = aflowlib::AFLUXCall(schema);
      vector<vector<std::pair<string,string> > > schema_response = aflowlib::getPropertiesFromAFLUXResponse(response);

			// extract units
			for(uint i=0;i<schema_response.size();i++){
				bool units_found = false;
				for(uint j=0;j<schema_response[i].size();j++){
					if(schema_response[i][j].first=="units"){
						property_units.push_back(schema_response[i][j].second);
						units_found=true;
						break;
					}
				}
				if(!units_found){
					property_units.push_back("");
				}
			}
      if(LDEBUG) {
        for(uint i=0;i<property_units.size();i++){ cerr << function_name << ": units for " << property_list[i] << ": " << property_units[i] << endl; }
      }
    }

    message << "Total number of candidate structures from database: " << auids.size();
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Loading structures ... ";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);


    vector<StructurePrototype> all_structures;
   
    // ---------------------------------------------------------------------------
    // store input structure 
    StructurePrototype input_structure;
    input_structure.representative_structure = xstr;
    input_structure.representative_structure.ReScale(1.0); //DX 20191105
    input_structure.representative_structure_name = "input geometry";
    input_structure.stoichiometry = compare::getStoichiometry(xstr,same_species);
    input_structure.elements = compare::getElements(xstr);
    input_structure.representative_structure_compound = compare::getCompoundName(xstr);
    //DX 20191105 [MOVED LATER - SAME AS SYMMETRY] input_structure.LFA_environments= compare::computeLFAEnvironment(input_structure.representative_structure); //DX 20190711
    input_structure.representative_structure_generated = true; 
    stringstream ss_input; ss_input << xstr;
    input_structure.representative_structure_from = ss_input.str(); 
    input_structure.property_names = property_list;
    input_structure.property_units = property_units;
    all_structures.push_back(input_structure);

    // ---------------------------------------------------------------------------
    // load and store entries from the database 
    for(uint i=0; i<auids.size(); i++){
      // first, get stoichiometry from entry
      vector<string> species; vector<double> natoms;
      XATOM_SplitAlloySpecies(compounds[i], species, natoms);
      vector<uint> tmp_stoich;
      for(uint j=0;j<natoms.size();j++){
        if(aurostd::isinteger(natoms[j])){
          tmp_stoich.push_back((uint)aurostd::nint(natoms[j]));
        }
        else {
          message << "Expected natoms in " << auids[i] << " to be an integer.";     
          throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_); //DX 20191031 - exit to xerror
        }
      }

      vector<uint> tmp_reduced_stoich = compare::gcdStoich(tmp_stoich);
      //DX 20190402 - need to sort if ignoring species - START
      if(!same_species){
        for(uint i=0; i<tmp_reduced_stoich.size(); i++){
	        std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
        }
      }
      //DX 20190402 - need to sort if ignoring species - END
      // second, check if stoichiometries are compatible
      // note: do not include in AFLUX matchbook, we would need to specify a range of compatible stoichs (could be expensive)
      // instead: filter on stoichiometry after recieving AFLUX response
      if(compare::sameStoichiometry(stoichiometry,tmp_reduced_stoich)){
        aflowlib::_aflowlib_entry entry; entry.auid=auids[i]; entry.aurl=aurls[i]; 
        if(!loadXstructures(entry,FileMESSAGE,oss)){ cerr << "WARNING::Could not load structure (auid=" << entry.auid << ") ... skipping..." << endl; continue;}
        if(entry.vstr.size()==1){
          // store entry from database
          StructurePrototype tmp;
          deque<string> deque_species; for(uint j=0;j<species.size();j++){deque_species.push_back(species[j]);}
          entry.vstr[0].SetSpecies(deque_species);
          tmp.representative_structure = entry.vstr[0];
          tmp.representative_structure.ReScale(1.0); //DX 20191105
          tmp.representative_structure_name=entry.getPathAURL(FileMESSAGE,oss,false); //DX 20190321 - changed to false, i.e., do not load from common
          tmp.representative_structure.directory=tmp.representative_structure_name; //DX 20190718 - update xstructure.directoryr
          tmp.representative_structure_generated=true;
          tmp.representative_structure_from="aurl";
          tmp.stoichiometry=tmp_reduced_stoich;
          tmp.representative_structure_compound = compare::getCompoundName(entry.vstr[0]); //DX 20190430 - added
          //DX 20191105 [MOVED LATER - SAME AS SYMMETRY] tmp.LFA_environments= compare::computeLFAEnvironment(tmp.representative_structure); //DX 20190711
          tmp.elements=species;
          tmp.number_of_atoms = entry.vstr[0].atoms.size(); //DX 20191031
          tmp.number_types = entry.vstr[0].num_each_type.size(); //DX 20191031
          // store any properties 
          for(uint l=0;l<properties_response[i].size();l++){
            bool property_requested = false;
            for(uint m=0;m<property_list.size();m++){
              if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
            }
            if(property_requested){
              tmp.representative_structure_properties.push_back(properties_response[i][l].second);
            }
          }
          if(LDEBUG) {
            cerr << "pflow::compareStructureDirectory() Found structure: " << tmp.representative_structure_name << endl;
          }
          all_structures.push_back(tmp);
        }
        else {
          message << "More structures loaded than anticipated for auid=" << auids[i] << ".";     
          throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_); //DX 20191031 - exit to xerror
        }
      }
    }
    message << "Total number of candidate structures loaded: " << all_structures.size(); //DX 20190403
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_); //DX 20190403

    // ---------------------------------------------------------------------------
    // calculate symmetry of database entries (need Wyckoff positions, but database is not sufficiently populated)
    // in the meantime, we calculate
    if(!ignore_symmetry){
      message << "Calculating the symmetry of the structures.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

      compare::calculateSymmetries(all_structures,num_proc); 

      message << "Symmetries calculated.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }  
    else {
      for(uint i=0; i<all_structures.size(); i++){
        all_structures[i].Pearson = "xX";
        all_structures[i].space_group = 0;
        vector<GroupedWyckoffPosition> tmp;
        all_structures[i].grouped_Wyckoff_positions = tmp;
      }
    }
    
    // ---------------------------------------------------------------------------
    // calculate LFA environments of  database entries 
    if(!ignore_environment){
      message << "Calculating the environments of the structures.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    
      compare::calculateLFAEnvironments(all_structures,num_proc); 
    
      message << "Environments calculated.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }  
    
    // ---------------------------------------------------------------------------
    // group into objects based on stoichiometry and symmetry (Pearson and space group)
    message << "Grouping sets of comparisons.";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> comparison_schemes = compare::groupStructurePrototypes(all_structures, same_species, ignore_symmetry, ignore_Wyckoff, ignore_environment, false); //DX 20190731 - add ignore_environment //DX 20190829 - false for duplicates_removed
    //cerr << "number of schemes: " << comparison_schemes.size() << endl;
    //cerr << "comparison_schemes: " << comparison_schemes.size() << endl;
    //cerr << "property names size: " << comparison_schemes[0].property_names.size() << endl;

    // ---------------------------------------------------------------------------
    // only compare entries to the input representation, the rest are extraneous comparisons
    vector<StructurePrototype> input_structure_comparison_scheme_only; input_structure_comparison_scheme_only.push_back(comparison_schemes[0]);
    message << "Number of structures to compare to input structure: " << input_structure_comparison_scheme_only[0].duplicate_structures_names.size();
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // compare structures
    message << "Running comparisons...";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes = compare::runComparisonScheme(num_proc, input_structure_comparison_scheme_only, same_species, false, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, false, store_comparison_logs, oss, FileMESSAGE); //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched //DX 20190731 - added ignore_symmetry/Wyckoff/environment //DX 20190822 - add log bool //DX 20190830 - false=duplicates_removed
    
    // ---------------------------------------------------------------------------
    // return if there are no similar structures
    if(final_prototypes.size()==0){
      return oss.str();
    }

    comparison_schemes.clear();

    // ---------------------------------------------------------------------------
    // print results 
    stringstream ss_out;
    compare::printResults(ss_out, same_species, final_prototypes);
    stringstream ss_json;
    compare::printResults(ss_json, same_species, final_prototypes, "json");
     
    // DEBUG oss << ss_out.str();
    message << "Number of structures in database matching with the input structure: " << final_prototypes[0].duplicate_structures.size() << "." << endl;
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    
    //DX 20190429 - added screen only option - START
    // ---------------------------------------------------------------------------
    // write results to screen and return immediately (do not write file)
    if(screen_only){
      if(format=="json"){ return ss_json.str(); }
      // default is txt
      else { return ss_out.str(); }
    }
    //DX 20190429 - added screen only option - END
    
    // ---------------------------------------------------------------------------
    // write results to files
    if(!structure_comparison){  
      if(format=="json"){
        aurostd::stringstream2file(ss_json,"database_material_comparison_output.json");
        message << "RESULTS: See database_material_comparison_output.json for list of unique/duplicate materials in database.";
      }
      else if(format=="txt"){
        aurostd::stringstream2file(ss_out,"database_material_comparison_output.out");
        message << "RESULTS: See database_material_comparison_output.out for list of unique/duplicate materials in database.";
      }
      else if(format=="both"){
        aurostd::stringstream2file(ss_json,"database_material_comparison_output.json");
        aurostd::stringstream2file(ss_out,"database_material_comparison_output.out");
        message << "RESULTS: See database_material_comparison_output.out or database_material_comparison_output.json for list of unique/duplicate materials in database.";
      }
    }
    else {
      if(format=="json"){
        aurostd::stringstream2file(ss_json,"database_structure_comparison_output.json");
        message << "RESULTS: See database_structure_comparison_output.json for list of unique/duplicate structures in database.";
      }
      else if(format=="txt"){
        aurostd::stringstream2file(ss_out,"database_structure_comparison_output.out");
        message << "RESULTS: See database_structure_comparison_output.out for list of unique/duplicate structures in database.";
      }
      else if(format=="both"){
        aurostd::stringstream2file(ss_json,"database_structure_comparison_output.json");
        aurostd::stringstream2file(ss_out,"database_structure_comparison_output.out");
        message << "RESULTS: See database_structure_comparison_output.out or database_structure_comparison_output.json for list of unique/duplicate structures in database.";
      }
    }
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);

    return oss.str();

  }
}

//DX - COMPARE DATABASE ENTRIES - START
// ***************************************************************************
// pflow::compareDatabaseEntries - compares database entries
// ***************************************************************************
namespace pflow {
  string compareDatabaseEntries(aurostd::xoption& vpflow){
    bool LDEBUG=(false || XHOST.DEBUG);
    
    string function_name = "pflow::compareDatabaseEntries()";
    string directory = ".";
    ostringstream oss;
    ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;

    string usage="aflow --compare_database_entries < POSCAR";
    string options="";

    vector<string> tokens,sub_tokens;
    vector<string> matchbook; //aflux - filter/get properties
    vector<string> schema; //aflux - get metadata (e.g., units)
    vector<string> property_units;

    bool same_species = true;
    bool single_comparison_round = false;
    bool store_comparison_logs = false; //DX 20190822 - add log bool

    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::NP"));
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: optimize match
    bool optimize_match=false; //false
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::OPTIMIZE_MATCH")) {
      optimize_match=true;
      message << "OPTIONS: Finding optimal match; exploring all possible lattices and origins to find the best match (note: this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    bool scale_volume=true;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::NO_SCALE_VOLUME")) {
      scale_volume=false;
      message << "OPTIONS: Suppressing volume scaling; useful for distinguishing structures at different pressures."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    bool ignore_Wyckoff=false;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::IGNORE_WYCKOFF")) {
      ignore_Wyckoff=true;
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    bool ignore_symmetry=false;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::IGNORE_SYMMETRY")) {
      ignore_symmetry=true;
      ignore_Wyckoff=true;
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    bool ignore_environment=false;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::IGNORE_ENVIRONMENT_ANALYSIS")) {
      ignore_environment=true;
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: specify number of species 
    uint arity=0; //Defalut=0 : all 
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::ARITY")) {
      arity=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::ARITY"));
    }
   
    // ---------------------------------------------------------------------------
    // FLAG: specify alloy systems
    vector<string> species;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::ALLOY")){
      string alloy_string = vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::ALLOY");
      // split by comma
      if(aurostd::substring2bool(alloy_string,",")){
        aurostd::string2tokens(alloy_string,species,",");
      }
      // split by colon
      else if(aurostd::substring2bool(alloy_string,":")){
        aurostd::string2tokens(alloy_string,species,":");
      }
      // split by alloy species (no delimiter)
      else{
        XATOM_SplitAlloySpecies(alloy_string, species);
      }
    }
    
    //TODO
    // ===== FLAG: STOICHIOMETRY ===== //
    //string alloy_string = vpflow.getattachedscheme("COMPARE_ALLOY::ALLOY");
    //TODO
    
    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool structure_comparison=false;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::STRUCTURE")) {
      structure_comparison=true;
      same_species = false;
      message << "OPTIONS: Structure-type comparison, i.e., ignore atomic species."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: do not remove unmatched structures from the StructurePrototype Object
    // keeps results of each comparison
    bool clean_unmatched=true;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::KEEP_UNMATCHED")) {
      clean_unmatched=false;
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: remove duplicate compounds (useful for non-biased statistics)
    bool remove_duplicate_compounds=false;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::REMOVE_DUPLICATE_COMPOUNDS")) {
      remove_duplicate_compounds=true;
      message << "OPTIONS: Remove duplicate compounds first, useful for non-biased prototype statistics."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: match unique structures to the AFLOW prototypes 
    bool add_matching_aflow_protos=false;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::MATCH_TO_AFLOW_PROTOS")) {
      add_matching_aflow_protos=true;
      message << "OPTIONS: Compare unique structures to the AFLOW prototypes."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: get AFLOW ANRL designation for unique structures
    bool get_aflow_prototype_designation=false;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::ADD_AFLOW_PROTOTYPE_DESIGNATION")) {
      get_aflow_prototype_designation=true;
      message << "OPTIONS: Cast unique structures into AFLOW standard designation."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: do not calculate unique permutations 
    bool calculate_unique_permutations=true;
    if(vpflow.flag("COMPARE_DATABSE_ENTRIES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS")) {
      calculate_unique_permutations=false;
      message << "OPTIONS: Do not calculate unique permutations."; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: property list to extract from database (using AFLUX)
    vector<string> property_list;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::PROPERTY_LIST")) {
      aurostd::string2tokens(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::PROPERTY_LIST"),property_list,",");

      // put properties in schema and matchbook for AFLUX call
      schema.push_back("schema("+vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::PROPERTY_LIST")+")"); //to get units
      matchbook.insert(matchbook.end(), property_list.begin(), property_list.end());

      message << "OPTIONS: Extracting the following properties: " << aurostd::joinWDelimiter(property_list,", "); 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: specify the geometry file to grab (orig, relax1, relax2, static, bands, POSCAR, CONTCAR)
    // DX TODO 
    string geometry_file = "";
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::GEOMETRY_FILE")) {
      geometry_file = vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::GEOMETRY_FILE");
      message << "OPTIONS: Structure type (POSCAR.orig, POSCAR.relax1, POSCAR.relax2, CONTCAR.relax1, ...): " << geometry_file << endl; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    
    // ---------------------------------------------------------------------------
    // FLAG: catalog (icsd, lib1, lib2, lib3, ...)
    string catalog = "";
    string catalog_summons = "";
    bool ICSD_comparison = false;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::CATALOG")) {
      catalog = aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::CATALOG")); //DX 20190718
      catalog_summons = "catalog(\'" + catalog + "\')";
      matchbook.push_back(catalog_summons);
      message << "OPTIONS: Catalog/library (icsd, lib1, lib2, lib3, ...): " << catalog << endl; 
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

      if(catalog=="" || catalog=="icsd" || catalog=="all"){
        ICSD_comparison=true;
      }
    }
    
    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get aurl for entry
    string aurl = "aurl";
    matchbook.push_back(aurl);
    
    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get species and number of species
    string species_summons = "";
    if(species.size()!=0){
      species_summons = "species(" + aurostd::joinWDelimiter(species,",") + ")";
    }

    string nspecies_summons = "";
    if(arity!=0){
      nspecies_summons = "nspecies(" + aurostd::utype2string<uint>(arity) + ")";
    }

    matchbook.push_back(species_summons);
    matchbook.push_back(nspecies_summons);

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: format AFLUX output 
    string format = "format(aflow)";
    string paging = "paging(0)";
    matchbook.push_back(format);
    matchbook.push_back(paging);

    // ---------------------------------------------------------------------------
    // construct aflux summons, i.e., combine matchbook
    string Summons = aurostd::joinWDelimiter(matchbook,",");
    message << "AFLUX matchbook request: " << Summons;
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // call AFLUX 
    string response = aflowlib::AFLUXCall(Summons);

    message << "Number of entries returned: " << aurostd::string2tokens(response,tokens,"\n");
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    if(LDEBUG){cerr << function_name << "::AFLUX response:" << endl << response << endl;}
   
    // ---------------------------------------------------------------------------
    // extract properties from AFLUX response
    vector<vector<std::pair<string,string> > > properties_response = aflowlib::getPropertiesFromAFLUXResponse(response);
    if(LDEBUG){
      for(uint i=0;i<properties_response.size();i++){
        for(uint j=0;j<properties_response[i].size();j++){
          cerr << properties_response[i][j].first << " = " << properties_response[i][j].second << ", ";
        }
        cerr << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // extract aurl, auid, and compound type from properties variable
    vector<string> auids, aurls, compounds;
    for(uint i=0;i<properties_response.size();i++){
      for(uint j=0;j<properties_response[i].size();j++){
        if(properties_response[i][j].first=="aurl"){
          aurls.push_back(properties_response[i][j].second);
        }
        if(properties_response[i][j].first=="auid"){
          auids.push_back(properties_response[i][j].second);
        }
        if(properties_response[i][j].first=="compound"){
          compounds.push_back(properties_response[i][j].second);
        }
      }
    }
    //cerr << "==============================" << endl;
    //::print(auids);
    //::print(aurls);
    //::print(compounds);
    
    // ---------------------------------------------------------------------------
    // get AFLUX schema, i.e., metadata (for the units)
    if(schema.size()>0){
      schema.push_back(format);
      schema.push_back(paging);

      // call AFLUX to get schema
      response = aflowlib::AFLUXCall(schema);
      vector<vector<std::pair<string,string> > > schema_response = aflowlib::getPropertiesFromAFLUXResponse(response);

      // extract units
      for(uint i=0;i<schema_response.size();i++){
        bool units_found = false;
        for(uint j=0;j<schema_response[i].size();j++){
          if(schema_response[i][j].first=="units"){
            property_units.push_back(schema_response[i][j].second);
            units_found=true;
            break;
          }
        }
        if(!units_found){
          property_units.push_back("");
        }
      }
      if(LDEBUG){
        for(uint i=0;i<property_units.size();i++){ cerr << function_name << ": units for " << property_list[i] << ": " << property_units[i] << endl; }
      }
    }

    message << "Total number of candidate structures from database: " << auids.size();
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Loading structures ..." << auids.size();
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
    message << "Splitting into threads...";     
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    // ---------------------------------------------------------------------------
    // distribute threads via indices
    uint number_of_structures = auids.size();
    uint num_threads = aurostd::min(num_proc,number_of_structures); // cannot have more threads than structures
    vector<uint> start_indices, end_indices;
    compare::splitTaskIntoThreads(number_of_structures,num_threads,start_indices,end_indices);
    message << "Done. Split into threads.";     
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // initialize vector of objects 
    vector<StructurePrototype> all_structures;
    for(uint i=0; i<auids.size(); i++){
      StructurePrototype tmp;
      // first, get stoichiometry from entry
      vector<string> species; vector<double> natoms;
      XATOM_SplitAlloySpecies(compounds[i], species, natoms);
      vector<uint> tmp_stoich;
      for(uint j=0;j<natoms.size();j++){
        if(aurostd::isinteger(natoms[j])){
          tmp_stoich.push_back((uint)aurostd::nint(natoms[j]));
        }
        else{
          message << "Expected natoms in " << auids[i] << " to be an integer.";     
          throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_); //DX 20191031 - exit to xerror
          //DX 20191031 pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
          //DX 20191031 exit(1);
        }
      }

      vector<uint> tmp_reduced_stoich = compare::gcdStoich(tmp_stoich);
      //DX 20190402 - need to sort if ignoring species - START
      if(!same_species){
        for(uint i=0; i<tmp_reduced_stoich.size(); i++){
	        std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
        }
      }
      tmp.stoichiometry=tmp_reduced_stoich;
      tmp.elements=species;
      tmp.representative_structure_name=aurls[i];
      tmp.representative_structure_from="aurl";
      all_structures.push_back(tmp);
    }
    message << "Finished initializing StructurePrototype object, now spawn threads.";     
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    
    // ---------------------------------------------------------------------------
    // Run threads 
    vector<std::thread> threads;
    for(uint n=0; n<num_threads; n++){
	    threads.push_back(std::thread(compare::generateStructuresInRange,std::ref(all_structures),std::ref(oss),start_indices[n],end_indices[n]));
    }
    // ---------------------------------------------------------------------------
    // Join threads
	  for(uint t=0;t<num_threads;t++){
      threads[t].join();
	  }
    message << "Threads complete. " << all_structures.size() << " structures. Adding properties.";     
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
   
    // ---------------------------------------------------------------------------
    // add properties information 
    for(uint i=0;i<all_structures.size();i++){
      all_structures[i].property_names = property_list; //DX 20190326
      all_structures[i].property_units = property_units; //DX 20190326
      // store any properties 
      for(uint l=0;l<properties_response[i].size();l++){
        bool property_requested = false;
        for(uint m=0;m<property_list.size();m++){
          if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
        }
        if(property_requested){
          all_structures[i].representative_structure_properties.push_back(properties_response[i][l].second);
        }
      }
    }

    for(uint i=0;i<all_structures.size();i++){
      if(all_structures[i].representative_structure_generated){
        deque<string> deque_species; for(uint j=0;j<all_structures[i].elements.size();j++){deque_species.push_back(all_structures[i].elements[j]);}
        all_structures[i].representative_structure.SetSpecies(deque_species);
        all_structures[i].representative_structure_compound = compare::getCompoundName(all_structures[i].representative_structure);
      }
    }

    message << "Properties added, now removing non-generated structures" << endl;     
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    // ---------------------------------------------------------------------------
    // remove non-generated structures
    compare::removeNonGeneratedStructures(all_structures);
          
    for(uint i=0;i<all_structures.size();i++){
      all_structures[i].representative_structure.ReScale(1.0); //DX 20191105
      //DX 20191105 [MOVED LATER - SAME AS SYMMETRY] all_structures[i].LFA_environments= compare::computeLFAEnvironment(all_structures[i].representative_structure); //DX 20190711
    }
         
#else
    // ---------------------------------------------------------------------------
    // load and store entries from the database 
    vector<StructurePrototype> all_structures;
    for(uint i=0; i<auids.size(); i++){
      // first, get stoichiometry from entry
      vector<string> species; vector<double> natoms;
      XATOM_SplitAlloySpecies(compounds[i], species, natoms);
      vector<uint> tmp_stoich;
      for(uint j=0;j<natoms.size();j++){
        if(aurostd::isinteger(natoms[j])){
          tmp_stoich.push_back((uint)aurostd::nint(natoms[j]));
        }
        else{
          message << "Expected natoms in " << auids[i] << " to be an integer.";     
          throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_); //DX 20191031 - exit to xerror
        }
      }

      vector<uint> tmp_reduced_stoich = compare::gcdStoich(tmp_stoich);
      //DX 20190402 - need to sort if ignoring species - START
      if(!same_species){
        for(uint i=0; i<tmp_reduced_stoich.size(); i++){
          std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
        }
      }
      //DX 20190402 - need to sort if ignoring species - END
      // second, check if stoichiometries are compatible
      // note: do not include in AFLUX matchbook, we would need to specify a range of compatible stoichs (could be expensive)
      // instead: filter on stoichiometry after recieving AFLUX response
      //if(compare::sameStoichiometry(stoichiometry,tmp_reduced_stoich)){
      aflowlib::_aflowlib_entry entry; entry.auid=auids[i]; entry.aurl=aurls[i]; 
      if(!loadXstructures(entry,FileMESSAGE,oss)){ cerr << "WARNING::Could not load structure (auid=" << entry.auid << ") ... skipping..." << endl; continue;}
      if(entry.vstr.size()==1){
        // store entry from database
        StructurePrototype tmp;
        deque<string> deque_species; for(uint j=0;j<species.size();j++){deque_species.push_back(species[j]);}
        entry.vstr[0].SetSpecies(deque_species);
        tmp.representative_structure = entry.vstr[0];
        tmp.representative_structure.ReScale(1.0); //DX 20191105
        tmp.representative_structure_name=entry.getPathAURL(FileMESSAGE,oss,false); //DX 20190321 - changed to false, i.e., do not load from common
        tmp.representative_structure.directory=tmp.representative_structure_name; //DX 20190718 - update xstructure.directoryr
        tmp.representative_structure_generated=true;
        tmp.representative_structure_from="aurl";
        tmp.stoichiometry=tmp_reduced_stoich;
        tmp.elements=species;
        tmp.number_of_atoms = entry.vstr[0].atoms.size(); //DX 20191031
        tmp.number_types = entry.vstr[0].num_each_type.size(); //DX 20191031
        tmp.representative_structure_compound = compare::getCompoundName(entry.vstr[0]);
        //DX 20191105 [MOVED LATER - SAME AS SYMMETRY] tmp.LFA_environments= compare::computeLFAEnvironment(tmp.representative_structure); //DX 20190711
        tmp.property_names = property_list; //DX 20190326
        tmp.property_units = property_units; //DX 20190326
        // store any properties 
        for(uint l=0;l<properties_response[i].size();l++){
          bool property_requested = false;
          for(uint m=0;m<property_list.size();m++){
            if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
          }
          if(property_requested){
            tmp.representative_structure_properties.push_back(properties_response[i][l].second);
          }
        }
        if(LDEBUG){
          cerr << "pflow::compareStructureDirectory() Found structure: " << tmp.representative_structure_name << endl;
        }
        all_structures.push_back(tmp);
      }
      else{
        message << "More structures loaded than anticipated for auid=" << auids[i] << ".";     
        throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_); //DX 20191031 - exit to xerror
      }
      //}
    }
#endif

    message << "Total number of candidate structures loaded: " << all_structures.size(); //DX 20190403
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_); //DX 20190403
   
    vector<StructurePrototype> final_prototypes = compare::compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, remove_duplicate_compounds, calculate_unique_permutations, add_matching_aflow_protos, get_aflow_prototype_designation, ICSD_comparison, store_comparison_logs); //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool

    // ---------------------------------------------------------------------------
    // print results 
    stringstream ss_out;
    compare::printResults(ss_out, same_species, final_prototypes);
    stringstream ss_json;
    compare::printResults(ss_json, same_species, final_prototypes, "json");

    // ---------------------------------------------------------------------------
    // write results to files
    if(!structure_comparison){  
      aurostd::stringstream2file(ss_json,directory+"/database_entries_material_comparison_output.json");
      aurostd::stringstream2file(ss_out,directory+"/database_entries_material_comparison_output.out");
      message << "RESULTS: See database_entries_material_comparison_output.out or database_entries_material_comparison_output.json for list of unique/duplicate materials in database.";
    }
    else{
      aurostd::stringstream2file(ss_json,directory+"/database_entries_structure_comparison_output.json");
      aurostd::stringstream2file(ss_out,directory+"/database_entries_structure_comparison_output.out");
      message << "RESULTS: See database_entries_structure_comparison_output.out or database_entries_structure_comparison_output.json for list of unique/duplicate structures in database.";
    }
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);

    return oss.str();
    
  }
}

//DX - COMPARE DATABASE ENTRIES - END

//DX 20190424 - START
// ***************************************************************************
// compare::compareStructuresFromStructureList()
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compareStructuresFromStructureList(vector<string>& filenames, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species,
    bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds, 
    bool calculate_unique_permutations, bool add_matching_aflow_protos, bool get_aflow_prototype_designation, bool ICSD_comparison, bool store_comparison_logs){ //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
    
    // ---------------------------------------------------------------------------
    // directory to write results
    string directory = "."; // for now this is fixed

    // ---------------------------------------------------------------------------
    // load structures appended to command
    vector<StructurePrototype> all_structures = compare::loadStructuresFromStructureList(filenames, magmoms_for_systems, same_species, FileMESSAGE); //DX 20190319 - added FileMESSAGE 
    
    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    return compare::compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, remove_duplicate_compounds, calculate_unique_permutations, add_matching_aflow_protos, get_aflow_prototype_designation, ICSD_comparison, store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool

  }
}
//DX 20190424 - END

// ***************************************************************************
// compare::compareStructuresFromDirectory()
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compareStructuresFromDirectory(string& directory, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species,
    bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds, 
    bool calculate_unique_permutations, bool add_matching_aflow_protos, bool get_aflow_prototype_designation, bool ICSD_comparison, bool store_comparison_logs){ //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool

    // ---------------------------------------------------------------------------
    // load structures in directory 
    vector<StructurePrototype> all_structures = compare::loadStructuresFromDirectory(directory, magmoms_for_systems, same_species, FileMESSAGE); //DX 20190319 - added FileMESSAGE
    
    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    return compare::compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, remove_duplicate_compounds, calculate_unique_permutations, add_matching_aflow_protos, get_aflow_prototype_designation, ICSD_comparison, store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool

  }
}

// ***************************************************************************
// compare::compareStructuresFromFile()
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compareStructuresFromFile(string& filename, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species,
    bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds, 
    bool calculate_unique_permutations, bool add_matching_aflow_protos, bool get_aflow_prototype_designation, bool ICSD_comparison, bool store_comparison_logs){ //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool

    // ---------------------------------------------------------------------------
    // directory to write results
    string directory = "."; // for now this is fixed

    // ---------------------------------------------------------------------------
    // load structures in file
    vector<StructurePrototype> all_structures = compare::loadStructuresFromFile(filename, magmoms_for_systems, same_species, FileMESSAGE); //DX 20190319 - added FileMESSAGE
    
    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    return compare::compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, remove_duplicate_compounds, calculate_unique_permutations, add_matching_aflow_protos, get_aflow_prototype_designation, ICSD_comparison, store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool

  }
}

// ***************************************************************************
// compare::compareMultipleStructures()
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& all_structures, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, string& directory){ //DX 20190319 - added FileMESSAGE
    bool scale_volume = true;
    bool optimize_match = false;
    bool ignore_symmetry = false;
    bool ignore_Wyckoff = false;
    bool ignore_environment = false; //DX 20190731
    bool single_comparison_round = false;
    bool store_comparison_logs = false; //DX 20190822 - add log bool
    bool clean_unmatched = true; //DX 20190506
    bool remove_duplicate_compounds = false;
    bool add_matching_aflow_protos = false; //DX 20190724
    bool get_aflow_prototype_designation = false; //DX 20190724
    bool calculate_unique_permutations = true; //DX 20190724
    bool ICSD_comparison = false;
    if(!same_species){
      remove_duplicate_compounds=true;
    }
    
    return compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, remove_duplicate_compounds, calculate_unique_permutations, add_matching_aflow_protos, get_aflow_prototype_designation, ICSD_comparison, store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
  }
}

namespace compare {
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& all_structures, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, string& directory,
    bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds, 
    bool calculate_unique_permutations, bool add_matching_aflow_protos, bool get_aflow_prototype_designation, bool ICSD_comparison, bool store_comparison_logs){ //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
    
    string function_name = "compare::compareMultipleStructures()";
    bool LDEBUG=(false || XHOST.DEBUG);
    ostream& logstream = cout;
    stringstream message;
    //DX 20190319 [OBSOLETE] ofstream FileMESSAGE;

    message << "Total number of structures to compare: " << all_structures.size();
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // calculate symmetries of structures
    // if already calculated, do not recalculate
    bool all_symmetries_calculated = true;
    for(uint i=0;i<all_structures.size();i++){ all_symmetries_calculated*=all_structures[i].isSymmetryCalculated(); }

    if(!ignore_symmetry && !all_symmetries_calculated){
      message << "Calculating the symmetry of the structures.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      compare::calculateSymmetries(all_structures,num_proc); 
      message << "Symmetries calculated.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }  
    else if(!all_symmetries_calculated){
      for(uint i=0; i<all_structures.size(); i++){
        all_structures[i].Pearson = "xX";
        all_structures[i].space_group = 0;
        vector<GroupedWyckoffPosition> tmp;
        all_structures[i].grouped_Wyckoff_positions = tmp;
      }
    }
    
    // ---------------------------------------------------------------------------
    // calculate LFA environments of  database entries 
    // if already calculated, do not recalculate
    bool all_environments_calculated = true;
    for(uint i=0;i<all_structures.size();i++){ all_environments_calculated*=all_structures[i].isLFAEnvironmentCalculated(); }
    if(!ignore_environment && !all_environments_calculated){
      message << "Calculating the environments of the structures.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    
      compare::calculateLFAEnvironments(all_structures,num_proc); 
    
      message << "Environments calculated.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }
   
    // ---------------------------------------------------------------------------
    // remove duplicate compounds first; uses recursion of this function
    if(remove_duplicate_compounds){
      bool tmp_same_species = true;
      bool tmp_scale_volume = false;
      bool tmp_remove_duplicates = false; // to avoid an infinite recursive loop
      bool tmp_calculate_unique_permutations = false; //DX 20190724
      bool tmp_add_matching_aflow_protos = false; //DX 20190724
      bool tmp_get_aflow_prototype_designation = false; //DX 20190724

      message << "Comparing to remove duplicate materials.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      vector<StructurePrototype> unique_compounds = compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, tmp_same_species, directory, 
          tmp_scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, 
          clean_unmatched, tmp_remove_duplicates, tmp_calculate_unique_permutations, tmp_add_matching_aflow_protos, tmp_get_aflow_prototype_designation, ICSD_comparison, store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool

      message << "Duplicate materials removed.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

      // include duplicate compounds count in object
      for(uint i=0;i<unique_compounds.size();i++){
        unique_compounds[i].number_compounds_matching_representative = unique_compounds[i].numberOfComparisons();
      }
    
      // prepare JSON output
      stringstream ss_json_remove_duplicates;
      compare::printResults(ss_json_remove_duplicates, tmp_same_species, unique_compounds, "json");
      
      // prepare TEXT (.out) output
      stringstream ss_out_remove_duplicates;
      compare::printResults(ss_out_remove_duplicates, tmp_same_species, unique_compounds, "txt");
      
     
      // write results to files
      if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
        aurostd::stringstream2file(ss_json_remove_duplicates,directory+"/duplicate_compounds_output.json");
        message << "RESULTS: See " << directory << "/duplicate_compounds_output.json for list of unique/duplicate structures.";
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
      }
      else if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
        aurostd::stringstream2file(ss_out_remove_duplicates,directory+"/duplicate_compounds_output.out");
        message << "RESULTS: See " << directory << "/duplicate_compounds_output.out for list of unique/duplicate structures.";
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
      }
      else {
        aurostd::stringstream2file(ss_json_remove_duplicates,directory+"/duplicate_compounds_output.json");
        aurostd::stringstream2file(ss_out_remove_duplicates,directory+"/duplicate_compounds_output.out");
        message << "RESULTS: See " << directory << "/duplicate_compounds_output.out" << " or " << directory << "/duplicate_compounds_output.json" << " for list of unique/duplicate structures.";
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
      }

      // overwrite original all_structures variable
      all_structures = unique_compounds;

      // remove all duplicate information from all_structures
      bool remove_duplicate_count = false; // but, keep the duplicate count information
      for(uint i=0;i<all_structures.size();i++){
        all_structures[i].removeDuplicates(remove_duplicate_count);
      }

    }
    

    // ---------------------------------------------------------------------------
    // group structures based on stoichiometry and symmetry (unless ignoring symmetry/Wyckoff)
    message << "Grouping sets of comparisons.";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    vector<StructurePrototype> comparison_schemes = compare::groupStructurePrototypes(all_structures, same_species, ignore_symmetry, ignore_Wyckoff, ignore_environment, false); //DX 20190731 - add ignore_environment //DX 20190829 - added remove_duplicate_compounds

    message << "Number of comparison groups: " << comparison_schemes.size() << ".";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
   
    // ---------------------------------------------------------------------------
    // if ICSD comparison, make structure with minimum ICSD number the representative structure
    if(ICSD_comparison){
      compare::representativePrototypeForICSDRuns(comparison_schemes);
    }

    // ---------------------------------------------------------------------------
    // compare structures 
    message << "Running comparisons ...";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes = compare::runComparisonScheme(num_proc, comparison_schemes, same_species, false, scale_volume, optimize_match, 
        ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, ICSD_comparison, store_comparison_logs, oss, FileMESSAGE);  //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched //DX 20190731 - added ignore_symmetry/Wyckoff/environment //DX 20190822 - add log bool
    
    // ---------------------------------------------------------------------------
    // return if there are no similar structures
    if(final_prototypes.size()==0){
      //return oss.str();
      return final_prototypes;
    }
    comparison_schemes.clear();

    // ---------------------------------------------------------------------------
    // combine prototypes regardless of having different space groups
    // BETA TESTING (perhaps we shouldn't do this) - compare::checkPrototypes(num_proc,same_species,final_prototypes);
 
    message << "Number of unique prototypes: " << final_prototypes.size() << " (out of " << all_structures.size() << " structures).";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);

    // ---------------------------------------------------------------------------
    // get unique permutations of prototype (representative) structures
    if(!same_species && calculate_unique_permutations){ 
      message << "Determining the unique permuations for each prototype.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      // find unique permutations of prototype
      for(uint i=0;i<final_prototypes.size();i++){
        if(arePermutationsComparableViaStoichiometry(final_prototypes[i].stoichiometry) && 
           arePermutationsComparableViaSymmetry(final_prototypes[i].grouped_Wyckoff_positions)){
          // check if xstructure is generated; if not, make it
          if(!final_prototypes[i].representative_structure_generated){
            if(!compare::generateStructure(final_prototypes[i].representative_structure_name,final_prototypes[i].representative_structure_from,final_prototypes[i].representative_structure,oss)){
              message << "Could not generate structure (" << final_prototypes[i].representative_structure_name << ").";
              throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_); //DX 20191031 - exit to xerror
            }
          }
          vector<StructurePrototype> final_permutations = compare::comparePermutations(final_prototypes[i],num_proc,optimize_match,oss,FileMESSAGE);
         
          // store permutation results in main StructurePrototype object
          for(uint j=0;j<final_permutations.size();j++){
            final_prototypes[i].unique_permutations.push_back(final_permutations[j].representative_structure_name);
          }
          final_permutations.clear(); //DX 20190624
        }
        else{
          vector<string> unique_permutations = generatePermutationString(final_prototypes[i].stoichiometry); //DX 20190508
          // store permutation results in main StructurePrototype object
          for(uint j=0;j<unique_permutations.size();j++){
            final_prototypes[i].unique_permutations.push_back(unique_permutations[j]);
          }
        }
      }
    }
   
    // ---------------------------------------------------------------------------
    // for structure-type comparisons, remove duplicate compounds (avoid biased duplicate statistics)
    if(same_species==false && remove_duplicate_compounds){
    //DX 20190220 [BETA]  message << "Performing comparisons to removing duplicate compounds";
    //DX 20190220 [BETA]  pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    //DX 20190220 [BETA]
    //DX 20190220 [BETA]  vector<StructurePrototype> duplicate_compound_comparisons = compare::compareDuplicateCompounds(final_prototypes, num_proc, ICSD_comparison, oss);
    //DX 20190220 [BETA]  
    //DX 20190220 [BETA]  message << "Removing duplicates from final comparison list";
    //DX 20190220 [BETA]  pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    //DX 20190220 [BETA]  // remove duplicates compounds from final comparison results
    //DX 20190220 [BETA]  compare::removeDuplicateCompounds(final_prototypes, duplicate_compound_comparisons);
    //DX 20190220 [BETA]
    //DX 20190220 [BETA]  // prepare JSON output
    //DX 20190220 [BETA]  stringstream ss_json_remove_duplicates;
    //DX 20190220 [BETA]  compare::printResults(ss_json_remove_duplicates, same_species, final_prototypes, "json");
    //DX 20190220 [BETA]  
    //DX 20190220 [BETA]  // prepare TEXT (.out) output
    //DX 20190220 [BETA]  stringstream ss_out_remove_duplicates;
    //DX 20190220 [BETA]  compare::printResults(ss_out_remove_duplicates, same_species, final_prototypes, "txt");
    //DX 20190220 [BETA]  
    //DX 20190220 [BETA]  // write results to files
    //DX 20190220 [BETA]  aurostd::stringstream2file(ss_json_remove_duplicates,directory+"/structure_comparison_no_duplicate_compounds_output.json");
    //DX 20190220 [BETA]  aurostd::stringstream2file(ss_out_remove_duplicates,directory+"/structure_comparison_no_duplicate_compounds_output.out");
    //DX 20190220 [BETA]  message << "RESULTS: See " << directory << "/structure_comparison_no_duplicate_compounds_output.out" << " or " << directory << "/structure_comparison_no_duplicate_compounds_output.json" << " for list of unique/duplicate structures.";
    //DX 20190220 [BETA]  pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }

    if(get_aflow_prototype_designation){
      message << "Determining the AFLOW standard designation.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
      // ---------------------------------------------------------------------------
      // split task into threads 
      uint number_of_structures = final_prototypes.size();
      uint number_of_threads = aurostd::min(num_proc,number_of_structures); // cannot have more threads than structures
      vector<uint> start_indices, end_indices;
      splitTaskIntoThreads(number_of_structures, number_of_threads, start_indices, end_indices);
      
			// ---------------------------------------------------------------------------
      // [THREADED] determine AFLOW standard designation 
      vector<std::thread> threads;
      for(uint n=0; n<number_of_threads; n++){
        threads.push_back(std::thread(compare::getPrototypeDesignationsInRange,std::ref(final_prototypes),start_indices[n], end_indices[n]));
			}
      for(uint t=0;t<threads.size();t++){
        threads[t].join();
	    }
      
      // ---------------------------------------------------------------------------
      // update once all are collected (safer) 
      for(uint i=0;i<final_prototypes.size();i++){
        final_prototypes[i].aflow_label = final_prototypes[i].representative_structure.prototype;
        final_prototypes[i].aflow_parameter_list = final_prototypes[i].representative_structure.prototype_parameter_list;
        final_prototypes[i].aflow_parameter_values = final_prototypes[i].representative_structure.prototype_parameter_values;
      }
#else
      // ---------------------------------------------------------------------------
      // [NON-THREADED] determine AFLOW standard designation 
      for(uint i=0;i<final_prototypes.size();i++){
        anrl::structure2anrl(final_prototypes[i].representative_structure,false); //DX 20190829 - false for recalculate_symmetry
        final_prototypes[i].aflow_label = final_prototypes[i].representative_structure.prototype;
        final_prototypes[i].aflow_parameter_list = final_prototypes[i].representative_structure.prototype_parameter_list;
        final_prototypes[i].aflow_parameter_values = final_prototypes[i].representative_structure.prototype_parameter_values;
      }
#endif
    }
      

    // ---------------------------------------------------------------------------
    // for testing/development; in case the subsequent analyses fails, checkpoint file 
    bool store_checkpoint=false;
    if(store_checkpoint){
      stringstream ss_json;
      compare::printResults(ss_json, same_species, final_prototypes, "json");
      stringstream ss_out;
      compare::printResults(ss_out, same_species, final_prototypes, "text");
      aurostd::stringstream2file(ss_json,directory+"/structure_comparison_output.json");
      aurostd::stringstream2file(ss_out,directory+"/structure_comparison_output.out");
      message << "RESULTS: See [tmp]" << directory << "/structure_comparison_output.out" << " or " << directory << "/structure_comparison_output.json" << " for list of unique/duplicate structures.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }

    if(add_matching_aflow_protos){
      message << "Determining if representative structures map to any of the AFLOW prototypes.";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

      aurostd::xoption vpflow_protos;
      vpflow_protos.flag("COMPARE2PROTOTYPES",TRUE);

      // ---------------------------------------------------------------------------
      // specify catalog
      vpflow_protos.flag("COMPARE2PROTOTYPES::CATALOG",TRUE);
      vpflow_protos.push_attached("COMPARE2PROTOTYPES::CATALOG","all"); 
      
      // ---------------------------------------------------------------------------
      // specify number of processors
      vpflow_protos.flag("COMPARE2PROTOTYPES::NP",TRUE);
      vpflow_protos.push_attached("COMPARE2PROTOTYPES::NP",aurostd::utype2string<uint>(num_proc)); 
    
			// ---------------------------------------------------------------------------
    	// do not calculate unique permutations since this was already done
    	vpflow_protos.flag("COMPARE2PROTOTYPES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS",TRUE);

      // ---------------------------------------------------------------------------
      // match to AFLOW prototypes 
      for(uint i=0;i<final_prototypes.size();i++){
        vector<StructurePrototype> matching_protos = pflow::compare2prototypes(final_prototypes[i].representative_structure, vpflow_protos);
        final_prototypes[i].matching_aflow_prototypes = matching_protos[0].duplicate_structures_names;
      }
    }
   
    return final_prototypes;
  }
}

// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, const bool &same_species) {
    ostringstream oss;
    uint num_proc=1;
    double final_misfit=-1;
    bool scale_volume=true; //default is true
    bool optimize_match=false; //default is false
    return aflowCompareStructure(num_proc, xstr1, xstr2, same_species, scale_volume, optimize_match, oss, final_misfit);
  }
}

namespace compare {
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, const bool& same_species, const bool& scale_volume, const bool& optimize_match) {
    ostringstream oss;
    uint num_proc = 1;
    double final_misfit = -1;
    return aflowCompareStructure(num_proc, xstr1, xstr2, same_species, scale_volume, optimize_match, oss, final_misfit);
  }
}


// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  double aflowCompareStructureMisfit(const xstructure& xstr1, const xstructure& xstr2, const bool &same_species) {
    ostringstream oss;
    uint num_proc=1;
    double final_misfit=-1;
    bool scale_volume=true; //default is true
    bool optimize_match=false; //default is false
    aflowCompareStructure(num_proc, xstr1, xstr2, same_species, scale_volume, optimize_match, oss, final_misfit);
    return final_misfit;
  }
}

// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const uint& num_proc, const xstructure& xstr1, const xstructure& xstr2, 
      const bool &same_species, const bool& scale_volume, const bool& optimize_match, 
      ostream& oss, double& final_misfit) {

    // This is the main comparison function, which  compares two crystal structures
    // and determines their level of similarity based on the idea discussed 
    // in H. Burzlaff's paper (Acta Cryst., A53, 217-224 (1997)).

    bool LDEBUG=(false || XHOST.DEBUG);

    oss << "==================================================================================" << endl;

    // ---------------------------------------------------------------------------
    // prepare structures (swap structure order if necessary, fix lattices, rescale, etc.)
    xstructure xstr_base, xstr_test;
    int mode = 1; // 0: expand structure with more atoms (robust using 3x3x3 method); 1: make xstr1 structure with smaller lattice point radius (faster using dims)
    // expand larger system to ensure we find a commensurate unit cell between each structure
    if(mode==0){
      if(xstr1.atoms.size()>xstr2.atoms.size()){
        xstr_base = xstr2;
        xstr_test = xstr1;
        if(LDEBUG) {
          cerr << "compare::aflowCompareStructure: WARNING: Swapping order of xstructure 1 and 2 since 1 is larger than the other." << endl;
        }
      }
      else { xstr_base = xstr1; xstr_test = xstr2; }
    }
    else if(mode==1){
      xvector<double> abc_angles_q1=Getabc_angles(xstr1.lattice,DEGREES);
      double radius_xstr1 = aurostd::max(abc_angles_q1(1),abc_angles_q1(2),abc_angles_q1(3)); 
      xvector<double> abc_angles_q2=Getabc_angles(xstr2.lattice,DEGREES);
      double radius_xstr2 = aurostd::max(abc_angles_q2(1),abc_angles_q2(2),abc_angles_q2(3)); 
      if(radius_xstr2<radius_xstr1){ 
        xstr_base = xstr2; 
        xstr_test = xstr1;
        if(LDEBUG) {
          cerr << "compare::aflowCompareStructure: WARNING: Swapping order of xstructure 1 and 2 since 2 has smaller lattice search radius than the other." << endl;
        }
      }
      else { xstr_base = xstr1; xstr_test = xstr2; }
    }
    else { xstr_base = xstr1; xstr_test = xstr2; }

    xstr_base.FixLattices();
    xstr_test.FixLattices();
    xstr_base.ReScale(1.0);
    xstr_test.ReScale(1.0);
    xstr_base.BringInCell(); //DX 20190329 - need to ensure incell; otherwise supercell expansion breaks
    xstr_test.BringInCell(); //DX 20190329 - need to ensure incell; otherwise supercell expansion breaks

    // ---------------------------------------------------------------------------
    // clean atom names (remove pseudopotential information)
    for(uint i=0;i<xstr_base.species.size();i++){ xstr_base.species[i]=KBIN::VASP_PseudoPotential_CleanName(xstr_base.species[i]); }
    for(uint i=0;i<xstr_test.species.size();i++){ xstr_test.species[i]=KBIN::VASP_PseudoPotential_CleanName(xstr_test.species[i]); }

    for(uint i=0;i<xstr_base.atoms.size();i++){ xstr_base.atoms[i].name=KBIN::VASP_PseudoPotential_CleanName(xstr_base.atoms[i].name); }
    for(uint i=0;i<xstr_test.atoms.size();i++){ xstr_test.atoms[i].name=KBIN::VASP_PseudoPotential_CleanName(xstr_test.atoms[i].name); }

    // ---------------------------------------------------------------------------
    // standardize structure (not default) 
    // below is no longer necessary, algorithm handles supercells/conventional/prim
    bool primitivize=false;
    bool niggli=false;
    if(primitivize){
      xstr_base=GetStandardPrimitive(xstr_base);
      xstr_test=GetStandardPrimitive(xstr_test);
    }
    if(niggli){
      xstr_base.NiggliUnitCellForm();
      xstr_test.NiggliUnitCellForm();
    }

    // ---------------------------------------------------------------------------
    // determine if structures are matchable (same species and/or same stoichiometry)
    int type_match=0;
    bool criteria_met = false;
    if(same_species == true){
      type_match=2;
      // if atoms are not labeled in either structure; assign fake names
      if(xstr_base.atoms.at(0).name == "" || xstr_test.atoms.at(0).name == ""){ 
        if(LDEBUG) {cerr << "compare:: " << "Atoms not labeled ... Assigning fake names." << endl;}
        fakeAtomsName(xstr_base);
        fakeAtomsName(xstr_test);
      }
    }
    else if(same_species == false){
      type_match=1;
    }
    if(matchableSpecies(xstr_base,xstr_test,same_species)==true){
      criteria_met = true;
    }
    //cerr << "type_match: " << type_match << endl; 

    // ---------------------------------------------------------------------------
    // standardize structure (not used) 
    if(criteria_met == true){
      oss << "=========================================================" << endl; 

      oss << "STRUCTURE 1: " << endl;  
      oss << xstr_base << endl;
      //cerr << xstr_base << endl;

      oss << "=========================================================" << endl;

      oss << "STRUCTURE 2: " << endl;
      oss << xstr_test << endl;	
      //cerr << xstr_test << endl;

      oss << "=========================================================" << endl;

      // ---------------------------------------------------------------------------
      // comparison types
      // type_match determines how the atoms should be matched
      //  1: assigns fake names to atoms (allows for structural comparison regardless of type of atom)
      //  2: uses the names given in POSCAR (allows for structural comparison of material; type of atom necessary)

      if(same_species==true){
        type_match=2;
      }
      if(same_species==false){
        type_match=1;
      }	

      // ---------------------------------------------------------------------------
      // variables
      uint i=0;
      xvector<double> origin;
      xstructure proto;           
      vector<xstructure> vprotos,vprotos_tmp;
      vector<vector<uint> > IM1, IM2;
      vector<vector<double> > vmin_dists;
      vector<uint> im1, im2;
      vector<string> PAIR1, PAIR2;
      double minMis=1;

      oss<<"-------------------------------------------------------"<<endl;

      // ---------------------------------------------------------------------------
      // normalize scaling factors 
      if(LDEBUG) {cerr << "compare:: " << "Scale structures."<<endl;} 
      // structures should already be scaled to the same scaling factor, below may be redundant
      rescaleStructure(xstr_base,xstr_test);

      // ---------------------------------------------------------------------------
      // scale volumes of structures
      if(scale_volume==true){atomicNumberDensity(xstr_base, xstr_test);}

      // ---------------------------------------------------------------------------
      // assign fake atom names 
      if(type_match==1){
        fakeAtomsName(xstr_base);
        fakeAtomsName(xstr_test);
      }

      // OBSOLETE THIS PRINTS OUT XSTRUCTURES WITH ATOM ZERO SHIFTED TO ORIGIN...
      // OBSOLETE oss<<"========================================================="<<endl;
      // OBSOLETE oss << xstr_base << endl;
      // OBSOLETE oss<<"========================================================="<<endl;
      // OBSOLETE oss << xstr_test << endl;		

      // ---------------------------------------------------------------------------
      // assign fake atom names 
      printParameters(xstr_base,oss);
      printParameters(xstr_test,oss);

      oss << "========================================================="<<endl;    
      oss << "QUADRUPLETS METHOD" << endl;

      xmatrix<double> q_base=xstr_base.lattice;; 
      
      // ---------------------------------------------------------------------------
      // compare structures
      if(LDEBUG) {cerr << "compare:: " << "WAIT... Computing quadruplets..."<<endl;} 
      // creates the threads for checking quadruplets (lattices)
      //DX 20190530 - OLD threadGeneration(num_proc,q_base,xstr_test,vprotos,xstr_base,type_match,optimize_match,minMis,oss);
      latticeAndOriginSearch(xstr_base,xstr_test,num_proc,q_base,vprotos,minMis,type_match,optimize_match,oss); //DX 20190530

      if(LDEBUG) {cerr << "compare:: " << "Total # of possible matching representations: " << vprotos.size() << endl;}	
      final_misfit=minMis;

      // ---------------------------------------------------------------------------
      // find matches
      // note, this is done in threadGeneration(), so it is redundant, we can save time  

      // This first match finder is based on the best fitting between each atoms. This means that 
      // the routine looks for the closest atoms in order to find the match.
      // Can happen that one atom is the best matching for more than one atom in the second structure;
      // in this case the match is cancelled (cleanMatch)

      //DX [OBSOLETE - best-match shift is expensive]      for(i=0; i<vprotos.size(); i++){
      //DX [OBSOLETE - best-match shift is expensive]        //cerr << "xstr_base " << xstr_base << endl;
      //DX [OBSOLETE - best-match shift is expensive]        //cerr << "vprotos[i] " << vprotos[i] << endl;
      //DX [OBSOLETE - best-match shift is expensive]        //cerr << "orig: " << endl;
      //DX [OBSOLETE - best-match shift is expensive]        //for(uint j=0;j<im1.size();j++){
      //DX [OBSOLETE - best-match shift is expensive]        //  cerr << im1[j] << " == " << im2[j] << endl; 
      //DX [OBSOLETE - best-match shift is expensive]        //}
      //DX [OBSOLETE - best-match shift is expensive]        im1.clear(); im2.clear();
      //DX [OBSOLETE - best-match shift is expensive]        //cerr << "after: " << endl;
      //DX [OBSOLETE - best-match shift is expensive]        vector<double> min_dists;
      //DX [OBSOLETE - best-match shift is expensive]        //findMatch(xstr_base,vprotos.at(i),im1,im2);
      //DX [OBSOLETE - best-match shift is expensive]        //cerr << "find new match" << endl;
      //DX [OBSOLETE - best-match shift is expensive]        findMatch(xstr_base,vprotos.at(i),im1,im2,min_dists,type_match);
      //DX [OBSOLETE - best-match shift is expensive]        //cerr << "im1.size(): " << im1.size() << endl;
      //DX [OBSOLETE - best-match shift is expensive]        //for(uint j=0;j<im1.size();j++){
      //DX [OBSOLETE - best-match shift is expensive]        //  cerr << im1[j] << " == " << im2[j] << " (" << min_dists[j] << ")" << endl; 
      //DX [OBSOLETE - best-match shift is expensive]        //}
      //DX [OBSOLETE - best-match shift is expensive]        if(cleanMatch(im2)==false && cleanMatch(im1)==false){
      //DX [OBSOLETE - best-match shift is expensive]          //cerr << "cleanMatch" << endl;
      //DX [OBSOLETE - best-match shift is expensive]          vprotos_tmp.push_back(vprotos.at(i));
      //DX [OBSOLETE - best-match shift is expensive]          IM1.push_back(im1);
      //DX [OBSOLETE - best-match shift is expensive]          IM2.push_back(im2);
      //DX [OBSOLETE - best-match shift is expensive]          vmin_dists.push_back(min_dists);
      //DX [OBSOLETE - best-match shift is expensive]        }
      //DX [OBSOLETE - best-match shift is expensive]      }
      //DX [OBSOLETE - best-match shift is expensive]      if(LDEBUG) {cerr << "compare:: " << "Number of matching representations: "<< vprotos_tmp.size() << endl;}
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]      vprotos.clear();
      //DX [OBSOLETE - best-match shift is expensive]      vprotos=vprotos_tmp;
      //DX [OBSOLETE - best-match shift is expensive]      vprotos_tmp.clear();
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]      if(vprotos.size()!=0){
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]        vector<vector<uint> > auxstr_base,auxstr_test;
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]        auxstr_base=IM1;
      //DX [OBSOLETE - best-match shift is expensive]        auxstr_test=IM2;	
      //DX [OBSOLETE - best-match shift is expensive]        IM1.clear();
      //DX [OBSOLETE - best-match shift is expensive]        IM2.clear();
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]        // sameAtomType allow to check that the atoms matched are of the same type.
      //DX [OBSOLETE - best-match shift is expensive]        // (can happen that certain transformation find matches between atoms of different type)
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]        for(i=0; i<vprotos.size(); i++){
      //DX [OBSOLETE - best-match shift is expensive]          //cerr << "Same atom" << endl;
      //DX [OBSOLETE - best-match shift is expensive]          //if(sameAtomType(xstr_base,vprotos.at(i),auxstr_base.at(i),auxstr_test.at(i),type_match)==true){
      //DX [OBSOLETE - best-match shift is expensive]          //  cerr << "IN Same atom" << endl;
      //DX [OBSOLETE - best-match shift is expensive]            vprotos_tmp.push_back(vprotos.at(i));
      //DX [OBSOLETE - best-match shift is expensive]            IM1.push_back(auxstr_base.at(i));
      //DX [OBSOLETE - best-match shift is expensive]            IM2.push_back(auxstr_test.at(i));
      //DX [OBSOLETE - best-match shift is expensive]          //}
      //DX [OBSOLETE - best-match shift is expensive]        }
      //DX [OBSOLETE - best-match shift is expensive]        if(LDEBUG) {cerr << "compare:: " << "Number of valid matches with the same type: " << vprotos_tmp.size() << endl;}
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]        vprotos.clear();
      //DX [OBSOLETE - best-match shift is expensive]        vprotos=vprotos_tmp;
      //DX [OBSOLETE - best-match shift is expensive]        vprotos_tmp.clear();
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]      }
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]      // ---------------------------------------------------------------------------
      //DX [OBSOLETE - best-match shift is expensive]      // calculate misfit values for all matching structures 
      //DX [OBSOLETE - best-match shift is expensive]      if(vprotos.size()!=0){
      //DX [OBSOLETE - best-match shift is expensive]        //  follows the computation of the figure of misfit described by Burzlaff: 
      //DX [OBSOLETE - best-match shift is expensive]        //  Burzlaff H., Malinovsky Y. (1996), "A Procedure for the Classification of Non-Organic Crystal structures."
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]        xstructure xstr_base_tmp = xstr_base;
      //DX [OBSOLETE - best-match shift is expensive]        vector<double> diag_sum1,diag_sum2,diag_diff1,diag_diff2;
      //DX [OBSOLETE - best-match shift is expensive]        double scale, lattdev;
      //DX [OBSOLETE - best-match shift is expensive]        vector<double> vLattDevs;
      //DX [OBSOLETE - best-match shift is expensive]        vector<double> vCoordDevs, vfails;
      //DX [OBSOLETE - best-match shift is expensive]        double coorddev=1e9, fail_figure=1e9;
      //DX [OBSOLETE - best-match shift is expensive]        double mis=1e9;
      //DX [OBSOLETE - best-match shift is expensive]        vector<double> misfits;
      //DX [OBSOLETE - best-match shift is expensive]        double min=1e9;
      //DX [OBSOLETE - best-match shift is expensive]        xstructure better;
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]        for(i=0; i<vprotos.size(); i++){
      //DX [OBSOLETE - best-match shift is expensive]          xstructure proto = vprotos[i];
      //DX [OBSOLETE - best-match shift is expensive]          diag_sum1.clear();	diag_sum2.clear();
      //DX [OBSOLETE - best-match shift is expensive]          diag_diff1.clear();	diag_diff2.clear();
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          scale=xstr_base.Volume()/vprotos.at(i).Volume();
      //DX [OBSOLETE - best-match shift is expensive]          scale=pow(scale,0.3333);
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          cellDiagonal(xstr_base,diag_sum1,diag_diff1,1);
      //DX [OBSOLETE - best-match shift is expensive]          cellDiagonal(vprotos.at(i),diag_sum2,diag_diff2,scale);
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          lattdev=latticeDeviation(diag_sum1,diag_sum2,diag_diff1,diag_diff2);
      //DX [OBSOLETE - best-match shift is expensive]          //cerr << "lattdev: " << lattdev << endl;
      //DX [OBSOLETE - best-match shift is expensive]          vLattDevs.push_back(lattdev);
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          vector<double> all_nn1 = computeNearestNeighbors(xstr_base); 
      //DX [OBSOLETE - best-match shift is expensive]          vector<double> all_nn_proto = computeNearestNeighbors(proto);
      //DX [OBSOLETE - best-match shift is expensive]          coordinateDeviation(xstr_base_tmp,proto,all_nn1,all_nn_proto,IM1.at(i),IM2.at(i),vmin_dists[i],coorddev,fail_figure);
      //DX [OBSOLETE - best-match shift is expensive]          //cerr << "coorddev: " << coorddev << endl;
      //DX [OBSOLETE - best-match shift is expensive]          //cerr << "fail_figure: " << fail_figure << endl;
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          vCoordDevs.push_back(coorddev);
      //DX [OBSOLETE - best-match shift is expensive]          vfails.push_back(fail_figure);
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          mis=computeMisfit(lattdev,coorddev,fail_figure);
      //DX [OBSOLETE - best-match shift is expensive]          min=mis;
      //DX [OBSOLETE - best-match shift is expensive]	        if(LDEBUG) {cerr << "compare:: " << "misfit: " << mis << "   (lattice deviation: " << lattdev << "  coordinate displacement: " 
      //DX [OBSOLETE - best-match shift is expensive]                          << coorddev << "  figure of fail: " << fail_figure << ")" << endl;}
      //DX [OBSOLETE - best-match shift is expensive]          misfits.push_back(mis);
      //DX [OBSOLETE - best-match shift is expensive]          vprotos_tmp.push_back(vprotos.at(i));
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          /*
      //DX [OBSOLETE - best-match shift is expensive]          //If between 0.1 and 0.2 will try shifting method to see if we can obtain a figure of misfit under 0.1
      //DX [OBSOLETE - best-match shift is expensive]          if(minMis <= 0.2 && minMis > 0.1){	
      //DX [OBSOLETE - best-match shift is expensive]            if(fail_figure!=0 && std::isnan(misfits.at(i))==false){
      //DX [OBSOLETE - best-match shift is expensive]              // This part is fundamental because it allows us to correct 
      //DX [OBSOLETE - best-match shift is expensive]              // a trivial simplification done during the transformation: 
      //DX [OBSOLETE - best-match shift is expensive]              // When the structure has been rotated, then it is shifted 
      //DX [OBSOLETE - best-match shift is expensive]              // to each atom and brought in the cell to look for matchings. 
      //DX [OBSOLETE - best-match shift is expensive]              // The atom shifted to the origin will coincide perfectly with 
      //DX [OBSOLETE - best-match shift is expensive]              // the atom in the origin for the reference structure. 
      //DX [OBSOLETE - best-match shift is expensive]              // However, this can lead to a matching failure between other 
      //DX [OBSOLETE - best-match shift is expensive]              // pairs of atoms in a position different from the origin 
      //DX [OBSOLETE - best-match shift is expensive]              // (see definition of failure on the paper). This routine 
      //DX [OBSOLETE - best-match shift is expensive]              // aims to take the structure from the atom in the origin and 
      //DX [OBSOLETE - best-match shift is expensive]              // move the entire structures little by little around this 
      //DX [OBSOLETE - best-match shift is expensive]              // position. The rigid translation of all the atoms allow us 
      //DX [OBSOLETE - best-match shift is expensive]              // to compute many different figures of misfit with the 
      //DX [OBSOLETE - best-match shift is expensive]              // possibility that some failures disapper returning in the 
      //DX [OBSOLETE - best-match shift is expensive]              // allowed tolerances.	
      //DX [OBSOLETE - best-match shift is expensive]              if(LDEBUG) {cerr << "compare:: " << "Attempting shift of structure since the minimum misfit is just above the similarity threshold..." << endl;}
      //DX [OBSOLETE - best-match shift is expensive]              proto=vprotos.at(i);	
      //DX [OBSOLETE - best-match shift is expensive]              for(j=0; j<proto.atoms.size(); j++){
      //DX [OBSOLETE - best-match shift is expensive]                if(proto.atoms.at(j).cpos==origin){
      //DX [OBSOLETE - best-match shift is expensive]                  delta=0.01*shortestDistance(proto,j);		
      //DX [OBSOLETE - best-match shift is expensive]                  inc=0.2*delta;
      //DX [OBSOLETE - best-match shift is expensive]                }
      //DX [OBSOLETE - best-match shift is expensive]              }
      //DX [OBSOLETE - best-match shift is expensive]              //cerr << "delta: " << delta << endl;
      //DX [OBSOLETE - best-match shift is expensive]              for(double j=-delta; j<=delta; j=j+inc){
      //DX [OBSOLETE - best-match shift is expensive]                //cerr << "j: " << j << endl;
      //DX [OBSOLETE - best-match shift is expensive]                for(double k=-delta; k<=delta; k=k+inc){
      //DX [OBSOLETE - best-match shift is expensive]                  for(double w=-delta; w<=delta; w=w+inc){
      //DX [OBSOLETE - best-match shift is expensive]                    proto=vprotos.at(i);
      //DX [OBSOLETE - best-match shift is expensive]                    for(uint iat=0; iat<proto.atoms.size(); iat++){
      //DX [OBSOLETE - best-match shift is expensive]                      proto.atoms.at(iat).cpos(1)+=j;
      //DX [OBSOLETE - best-match shift is expensive]                      proto.atoms.at(iat).cpos(2)+=k;
      //DX [OBSOLETE - best-match shift is expensive]                      proto.atoms.at(iat).cpos(3)+=w;
      //DX [OBSOLETE - best-match shift is expensive]                    }	
      //DX [OBSOLETE - best-match shift is expensive]                    	
      //DX [OBSOLETE - best-match shift is expensive]                    vector<double> min_dists; 
      //DX [OBSOLETE - best-match shift is expensive]                    findMatch(xstr_base,proto,im1,im2,min_dists,type_match);
      //DX [OBSOLETE - best-match shift is expensive]                    //coordinateDeviation(xstr_base,proto,IM1.at(i),IM2.at(i),coorddev,fail_figure);
      //DX [OBSOLETE - best-match shift is expensive]                    coordinateDeviation(xstr_base,proto,im1,im2,coorddev,fail_figure);
      //DX [OBSOLETE - best-match shift is expensive]                    //cerr << "lattdev: " << lattdev << endl;
      //DX [OBSOLETE - best-match shift is expensive]                    //cerr << "coorddev: " << coorddev << endl;
      //DX [OBSOLETE - best-match shift is expensive]                    mis_tmp=computeMisfit(lattdev,coorddev,fail_figure);
      //DX [OBSOLETE - best-match shift is expensive]                    //cerr << "mis_tmp: " << mis_tmp << endl;
      //DX [OBSOLETE - best-match shift is expensive]		                if(mis_tmp<min){
      //DX [OBSOLETE - best-match shift is expensive]                      min=mis_tmp;
      //DX [OBSOLETE - best-match shift is expensive]                      coorddev_tmp=coorddev;
      //DX [OBSOLETE - best-match shift is expensive]                      fail_figure_tmp=fail_figure;
      //DX [OBSOLETE - best-match shift is expensive]                      better=proto;
      //DX [OBSOLETE - best-match shift is expensive]                    }
      //DX [OBSOLETE - best-match shift is expensive]                  }
      //DX [OBSOLETE - best-match shift is expensive]                }
      //DX [OBSOLETE - best-match shift is expensive]              }		
      //DX [OBSOLETE - best-match shift is expensive]              if(min<mis){
      //DX [OBSOLETE - best-match shift is expensive]                vprotos_tmp.at(i)=better;
      //DX [OBSOLETE - best-match shift is expensive]                vCoordDevs.at(i)=coorddev_tmp;
      //DX [OBSOLETE - best-match shift is expensive]                vfails.at(i)=fail_figure_tmp;
      //DX [OBSOLETE - best-match shift is expensive]                misfits.at(i)=min;
      //DX [OBSOLETE - best-match shift is expensive]              }
      //DX [OBSOLETE - best-match shift is expensive]            }
      //DX [OBSOLETE - best-match shift is expensive]          }
      //DX [OBSOLETE - best-match shift is expensive]          */
      //DX [OBSOLETE - best-match shift is expensive]        }
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]        vprotos.clear();
      //DX [OBSOLETE - best-match shift is expensive]        vprotos=vprotos_tmp;
      //DX [OBSOLETE - best-match shift is expensive]        vprotos_tmp.clear();
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]        if(LDEBUG) {cerr << "compare:: " << "Number of Misfits Computed:	"<<vprotos.size()<<endl;}
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]        // ---------------------------------------------------------------------------
      //DX [OBSOLETE - best-match shift is expensive]        // print results 
      //DX [OBSOLETE - best-match shift is expensive]        if(vprotos.size()!=0){
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          uint min_index = 0; //DX 5/14/18 - added initialization
      //DX [OBSOLETE - best-match shift is expensive]          int flag=0;
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          for(i=0; i<misfits.size(); i++){
      //DX [OBSOLETE - best-match shift is expensive]            if(flag==0){
      //DX [OBSOLETE - best-match shift is expensive]              min=misfits.at(i);
      //DX [OBSOLETE - best-match shift is expensive]              min_index=i;
      //DX [OBSOLETE - best-match shift is expensive]              flag=1;
      //DX [OBSOLETE - best-match shift is expensive]            }
      //DX [OBSOLETE - best-match shift is expensive]            else {
      //DX [OBSOLETE - best-match shift is expensive]              if(misfits.at(i)<min){
      //DX [OBSOLETE - best-match shift is expensive]                min=misfits.at(i);
      //DX [OBSOLETE - best-match shift is expensive]                min_index=i;
      //DX [OBSOLETE - best-match shift is expensive]              }   
      //DX [OBSOLETE - best-match shift is expensive]            }   
      //DX [OBSOLETE - best-match shift is expensive]          }	   
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          oss << endl <<"**************************** RESULTS ****************************"<<endl;
      //DX [OBSOLETE - best-match shift is expensive]          final_misfit=misfits.at(min_index);
      //DX [OBSOLETE - best-match shift is expensive]          if(misfits.at(min_index)<0.1){
      //DX [OBSOLETE - best-match shift is expensive]            oss << endl <<"MISFIT" <<":			" << misfits.at(min_index)<<"  STRUCTURES ARE COMPATIBLE" << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss <<"----------------------------------------------------"<<endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << "Figure of Deviation:	"<< vLattDevs.at(min_index) << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << "Figure of Displacement:	"<<vCoordDevs.at(min_index) << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << "Figure of Failure:	"<<vfails.at(min_index) << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss <<"----------------------------------------------------"<<endl;
      //DX [OBSOLETE - best-match shift is expensive]            printMatch(IM1.at(min_index),IM2.at(min_index),vprotos.at(min_index),xstr_base,oss);
      //DX [OBSOLETE - best-match shift is expensive]            oss <<"----------------------------------------------------"<<endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << "FINAL - REFERENCE STRUCTURE: " << endl;	
      //DX [OBSOLETE - best-match shift is expensive]            oss << xstr_base << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss <<"----------------------------------------------------"<<endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << "FINAL - MAPPED STRUCTURE: " << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << vprotos.at(min_index);
      //DX [OBSOLETE - best-match shift is expensive]          }
      //DX [OBSOLETE - best-match shift is expensive]          else {
      //DX [OBSOLETE - best-match shift is expensive]            if(misfits.at(min_index)<0.2){
      //DX [OBSOLETE - best-match shift is expensive]              oss << endl <<"MISFIT" <<":                    " << misfits.at(min_index)<<"  STRUCTURES ARE IN THE SAME FAMILY" << endl;
      //DX [OBSOLETE - best-match shift is expensive]            }
      //DX [OBSOLETE - best-match shift is expensive]            else {
      //DX [OBSOLETE - best-match shift is expensive]              oss << endl <<"MISFIT" <<":			" << misfits.at(min_index)<<"  STRUCTURES ARE INCOMPATIBLE (No match found)" << endl;
      //DX [OBSOLETE - best-match shift is expensive]            }
      //DX [OBSOLETE - best-match shift is expensive]            oss <<"----------------------------------------------------"<<endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << "Figure of Deviation:	"<< vLattDevs.at(min_index) << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << "Figure of Displacement:	"<<vCoordDevs.at(min_index) << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << "Figure of Failure:	"<<vfails.at(min_index) << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss <<"----------------------------------------------------"<<endl;
      //DX [OBSOLETE - best-match shift is expensive]            printMatch(IM1.at(min_index),IM2.at(min_index),vprotos.at(min_index),xstr_base,oss);
      //DX [OBSOLETE - best-match shift is expensive]            oss <<"----------------------------------------------------"<<endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << "FINAL - REFERENCE STRUCTURE: " << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << xstr_base << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss <<"----------------------------------------------------"<<endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << "FINAL - MAPPED STRUCTURE: " << endl;
      //DX [OBSOLETE - best-match shift is expensive]            oss << vprotos.at(min_index) << endl;
      //DX [OBSOLETE - best-match shift is expensive]          } 
      //DX [OBSOLETE - best-match shift is expensive]
      //DX [OBSOLETE - best-match shift is expensive]          //---------------------------------------------------------------------------//
      //DX [OBSOLETE - best-match shift is expensive]        }
      //DX [OBSOLETE - best-match shift is expensive]        else {
      //DX [OBSOLETE - best-match shift is expensive]          oss << "[ERROR]: No match found!" << endl;
      //DX [OBSOLETE - best-match shift is expensive]        }
      //DX [OBSOLETE - best-match shift is expensive]      }
      //DX [OBSOLETE - best-match shift is expensive]      else { 
      //DX [OBSOLETE - best-match shift is expensive]        oss << "[ERROR]: No match found!" << endl;
      //DX [OBSOLETE - best-match shift is expensive]      }
      //DX [OBSOLETE - best-match shift is expensive]      oss << endl << "*********************  THE END - FINE  **********************" << endl << endl;
      if(final_misfit<0.1 && !((final_misfit+1.0)<1e-3)){
        return true;
      }
      else {
        return false;
      }
    } //end of criteria_met
    return false;
  }
} //end of compare namespace

// AFLOW-XTAL-MATCH (compare crystal structures)
// Written by David Hicks (david.hicks@duke.edu) 
// Contributors: Carlo De Santo
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                                                                         *
// ***************************************************************************
