// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow DAVID HICKS - Duke University 2014-2020                 *
// *                                                                         *
// ***************************************************************************
// AFLOW-XtalMatch (compare crystal structures)
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
// compare::comparePermutations()
// ***************************************************************************
namespace compare {
  string comparePermutations(istream& input, const aurostd::xoption& vpflow){
    ostringstream results_ss;
    ostringstream oss;
    ofstream FileMESSAGE; //DX20190319 - added FileMESSAGE

    string usage="aflow --compare_permutation<POSCAR";
    string options="[--usage] [--np=<number>] [--print_misfit]";

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE_PERMUTATION::USAGE")) {
      //[CO20200624 - OBSOLETE]stringstream ss_usage;
      //[CO20200624 - OBSOLETE]init::ErrorOption(ss_usage,vpflow.getattachedscheme("COMPARE_PERMUTATION"),"compare::comparePermutations()",aurostd::liststring2string(usage,options));
      //[CO20200624 - OBSOLETE]return ss_usage.str();
      init::ErrorOption(vpflow.getattachedscheme("COMPARE_PERMUTATION"),"compare::comparePermutations()",aurostd::liststring2string(usage,options));
    }
    
    // ---------------------------------------------------------------------------
    // instantiate XtalFinder calculator 
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        1,
        oss);
    
    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions("permutation"); //DX20200103
    xtal_finder.getOptions(vpflow, comparison_options);

    /*
    // ---------------------------------------------------------------------------
    // FLAG: misfit threshold //DX20201119
    double misfit_match_threshold = DEFAULT_XTALFINDER_MISFIT_MATCH;
    double misfit_family_threshold = DEFAULT_XTALFINDER_MISFIT_FAMILY;
    if(vpflow.flag("COMPARE_PERMUTATION::MISFIT_MATCH")) {
      misfit_match_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE_PERMUTATION::MISFIT_MATCH"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_MATCH",vpflow.getattachedscheme("COMPARE_PERMUTATION::MISFIT_MATCH"));
    }
    if(vpflow.flag("COMPARE_PERMUTATION::MISFIT_FAMILY")) {
      misfit_family_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE_PERMUTATION::MISFIT_FAMILY"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_FAMILY",vpflow.getattachedscheme("COMPARE_PERMUTATION::MISFIT_FAMILY"));
    }
    // match threshold must be less than family threshold
    if(misfit_match_threshold>misfit_family_threshold){
      message << "Matching misfit threshold must be less than the same family threshold:"
        << " misfit_match_threshold: " << misfit_match_threshold
        << " misfit_family_threshold: " << misfit_family_threshold;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ILLEGAL_);
    }
    message << "Misfit theshold for matched structures: " << misfit_match_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_MATCH << ")"; 
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Misfit theshold for structures in the same family: " << misfit_family_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_FAMILY << ")";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    */

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

    //DX20190504 START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    string format = "text";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      format = "text";
    }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      format = "json";
    }
    //DX20190504 END

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
    //vector<string> unique_permutations = compare::getUniquePermutations(xstr, num_proc, optimize_match, print_misfit, oss, FileMESSAGE); //DX20190319 - added FileMESSAGE
    vector<string> unique_permutations = xtal_finder.getUniquePermutations(xstr, num_proc, optimize_match, print_misfit, results_ss, comparison_options); //DX20190319 - added FileMESSAGE

    return results_ss.str();
  }
}

/*
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
    ofstream FileMESSAGE; //DX20190319 - added FileMESSAGE
    return compare::getUniquePermutations(xstr, num_proc, optimize_match, print_misfit, oss, FileMESSAGE); //DX20190319 - added FileMESSAGE
  }
}

namespace compare{
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc, bool& optimize_match, bool& print_misfit, ostream& oss, ofstream& FileMESSAGE){ //DX20190319 - added FileMESSAGE

    vector<string> unique_permutations;
    stringstream ss_output; //DX20190506

    //DX20190506 START
    // ---------------------------------------------------------------------------
    // print format 
    string format = "text";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      format = "text";
    }
    if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      format = "json";
    }
    //DX20190506 END

    // ---------------------------------------------------------------------------
    // permutation comparisons must compare the same species 
    bool same_species=true; 

    // ---------------------------------------------------------------------------
    // quick check: check if any sites have the same number of atoms; if not, then no need to try comparing
    if(!print_misfit){
      if(!arePermutationsComparableViaComposition(xstr)){ //DX20190624 - put into function
        //DX20191125 [OBSOLETE] vector<uint> reduced_stoichiometry = gcdStoich(xstr.num_each_type); //DX20190508
        deque<int> reduced_stoichiometry; aurostd::reduceByGCD(xstr.num_each_type, reduced_stoichiometry); //DX20191125
        deque<uint> reduced_stoichiometry_uint; for(uint i=0;i<reduced_stoichiometry.size(); i++){ reduced_stoichiometry_uint.push_back((uint)reduced_stoichiometry[i]); } //DX20191125
        generatePermutationString(reduced_stoichiometry_uint, unique_permutations); //DX20190508
        if(format=="text"){ //DX20190506
          ss_output << "Unique atom decorations (" << unique_permutations.size() << "): " << endl; 
          ss_output << " " << aurostd::joinWDelimiter(unique_permutations,"\n ") << endl;
        }
        if(format=="json"){ //DX20190506
          ss_output << "{\"atom_decorations_equivalent\":["; 
          ss_output << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(unique_permutations,"\""),",") << "]"; //DX20191125 - Vec to Dec
          ss_output << "]}" << endl;
        }
        oss << ss_output.str();
        return unique_permutations;
      }
    }

    // ---------------------------------------------------------------------------
    // load input structure
    StructurePrototype structure;
    structure.structure_representative = xstr;
    structure.structure_representative_name = "input geometry";
    structure.stoichiometry = xstr.GetReducedComposition(false);
    structure.elements = xstr.GetElements(true,true); // true: clean names
    // update xstructure species
    if(structure.structure_representative.species.size()==0){
      deque<string> deque_species; for(uint j=0;j<structure.elements.size();j++){deque_species.push_back(structure.elements[j]);}
      structure.structure_representative.SetSpecies(deque_species);
      structure.structure_representative.SpeciesPutAlphabetic();
    }
    structure.structure_representative_generated = true; 
    structure.structure_representative_source = "input";
    structure.structure_representative_relaxation_step = 0; //DX20200429 input is assumed to be unrelaxed

    // ---------------------------------------------------------------------------
    // get the unique atom decorations for the structure
    vector<StructurePrototype> final_permutations = compare::comparePermutations(structure,num_proc,optimize_match,oss,FileMESSAGE); //DX20190319 - added FileMESSAGE

    // ---------------------------------------------------------------------------
    // print results
    if(format=="text"){ //DX20190506
      ss_output << "Unique atom decorations (" << final_permutations.size() << "): " << endl; 

      for(uint j=0;j<final_permutations.size();j++){
        ss_output << " " << final_permutations[j].structure_representative_name;
        for (uint k=0;k<final_permutations[j].structures_duplicate_names.size();k++){
          ss_output << " = " << final_permutations[j].structures_duplicate_names[k];
        }
        ss_output << endl;
      }
    }
    //DX20190506 START
    else if(format=="json"){
      stringstream sscontent_json;
      vector<string> vcontent_json;
      sscontent_json << "\"atom_decorations_equivalent\":[";
      for(uint j=0;j<final_permutations.size();j++){
        stringstream sstmp;
        vector<string> equivalent_permutations;
        equivalent_permutations.push_back(final_permutations[j].structure_representative_name);
        equivalent_permutations.insert(equivalent_permutations.end(),final_permutations[j].structures_duplicate_names.begin(),final_permutations[j].structures_duplicate_names.end());
        sstmp << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(equivalent_permutations,"\""),",") << "]";
        vcontent_json.push_back(sstmp.str()); sstmp.str("");
      }
      sscontent_json << aurostd::joinWDelimiter(vcontent_json,",");
      sscontent_json << "]";
      ss_output << "{" << sscontent_json.str() << "}" << endl;
    }
    //DX20190506 END

    // ---------------------------------------------------------------------------
    // print misfit results
    if(print_misfit){
      if(format=="text"){ //DX20190506
        ss_output << "Misfit values: " << endl; 
        stringstream ss_text;
        compare::printResults(ss_text, same_species, final_permutations, "text");
        ss_output << ss_text.str();
      }
      else if(format=="json"){ //DX20190506
        ss_output.str(""); // need to clear content abbreviated content from above
        stringstream ss_json;
        compare::printResults(ss_json, same_species, final_permutations, "json");
        ss_output << ss_json.str() << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // store unique atom decorations in vector
    for(uint j=0;j<final_permutations.size();j++){
      unique_permutations.push_back(final_permutations[j].structure_representative_name);
    }

    // update oss
    oss << ss_output.str();
    return unique_permutations;
  }
}
*/

// ***************************************************************************
// XtalFinderCalculator::getUniquePermutations()
// ***************************************************************************
  vector<string> XtalFinderCalculator::getUniquePermutations(xstructure& xstr){
    uint num_proc=1;
    return getUniquePermutations(xstr, num_proc);
  }

  vector<string> XtalFinderCalculator::getUniquePermutations(xstructure& xstr, uint num_proc){
    bool optimize_match=false;
    return getUniquePermutations(xstr, num_proc, optimize_match);
  }

  vector<string> XtalFinderCalculator::getUniquePermutations(xstructure& xstr, uint num_proc, bool optimize_match){
    bool print_misfit=false;
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions("permutation"); //DX20200103
    stringstream oss;
    return getUniquePermutations(xstr, num_proc, optimize_match, print_misfit, oss, comparison_options); //DX20190319 - added FileMESSAGE
  }

  vector<string> XtalFinderCalculator::getUniquePermutations(xstructure& xstr, uint num_proc, bool optimize_match, bool print_misfit, ostream& oss, aurostd::xoption& comparison_options){ //DX20190319 - added FileMESSAGE

    vector<string> unique_permutations;
    stringstream ss_output; //DX20190506

    //DX20190506 START
    // ---------------------------------------------------------------------------
    // print format 
    string format = "text";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      format = "text";
    }
    if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      format = "json";
    }
    //DX20190506 END

    // ---------------------------------------------------------------------------
    // permutation comparisons must compare the same species 
    bool same_species=true; 

    // ---------------------------------------------------------------------------
    // quick check: check if any sites have the same number of atoms; if not, then no need to try comparing
    if(!print_misfit){
      if(!compare::arePermutationsComparableViaComposition(xstr)){ //DX20190624 - put into function
        //DX20191125 [OBSOLETE] vector<uint> reduced_stoichiometry = gcdStoich(xstr.num_each_type); //DX20190508
        deque<int> reduced_stoichiometry; aurostd::reduceByGCD(xstr.num_each_type, reduced_stoichiometry); //DX20191125
        deque<uint> reduced_stoichiometry_uint; for(uint i=0;i<reduced_stoichiometry.size(); i++){ reduced_stoichiometry_uint.push_back((uint)reduced_stoichiometry[i]); } //DX20191125
        compare::generatePermutationString(reduced_stoichiometry_uint, unique_permutations); //DX20190508
        //if(format=="text"){ //DX20190506
        //  ss_output << "Unique atom decorations (" << unique_permutations.size() << "): " << endl; 
        //  ss_output << " " << aurostd::joinWDelimiter(unique_permutations,"\n ") << endl;
        //}
        //if(format=="json"){ //DX20190506
        //  ss_output << "{\"atom_decorations_equivalent\":["; 
        //  ss_output << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(unique_permutations,"\""),",") << "]"; //DX20191125 - Vec to Dec
        //  ss_output << "]}" << endl;
        //}
        //oss << ss_output.str();
        return unique_permutations;
      }
    }

    // ---------------------------------------------------------------------------
    // load input structure
    stringstream xstr_ss; xstr_ss << xstr;
    //addStructureToContainer(xstr, "input geometry", xstr_ss.str(), 0, same_species);
    addStructureToContainer(xstr, "input geometry", xstr_ss.str(), 0, false); // false: so we can find decorations on systems without atoms
    StructurePrototype structure;
    uint container_index = 0;
    addToStructureRepresentative(structure,container_index);
   
    /*
    StructurePrototype structure;
    structure.structure_representative = xstr;
    structure.structure_representative_name = "input geometry";
    structure.stoichiometry = xstr.GetReducedComposition(false);
    structure.elements = xstr.GetElements(true,true); // true: clean names
    // update xstructure species
    if(structure.structure_representative.species.size()==0){
      deque<string> deque_species; for(uint j=0;j<structure.elements.size();j++){deque_species.push_back(structure.elements[j]);}
      structure.structure_representative.SetSpecies(deque_species);
      structure.structure_representative.SpeciesPutAlphabetic();
    }
    structure.structure_representative_generated = true; 
    structure.structure_representative_source = "input";
    structure.structure_representative_relaxation_step = 0; //DX20200429 input is assumed to be unrelaxed
    */

    // ---------------------------------------------------------------------------
    // get the unique atom decorations for the structure
    XtalFinderCalculator xtal_finder_permutations;
    vector<StructurePrototype> final_permutations = xtal_finder_permutations.comparePermutations(structure,num_proc,optimize_match); //DX20190319 - added FileMESSAGE

    // ---------------------------------------------------------------------------
    // print results
    if(format=="text"){ //DX20190506
      ss_output << "Unique atom decorations (" << final_permutations.size() << "): " << endl; 
      for(uint j=0;j<final_permutations.size();j++){
        ss_output << " " << final_permutations[j].structure_representative_struct->name;
        for (uint k=0;k<final_permutations[j].structures_duplicate_struct.size();k++){
          ss_output << " = " << final_permutations[j].structures_duplicate_struct[k]->name;
        }
        ss_output << endl;
      }
    }
    //DX20190506 START
    else if(format=="json"){
      stringstream sscontent_json;
      vector<string> vcontent_json;
      sscontent_json << "\"atom_decorations_equivalent\":[";
      for(uint j=0;j<final_permutations.size();j++){
        stringstream sstmp;
        vector<string> equivalent_permutations;
        equivalent_permutations.push_back(final_permutations[j].structure_representative_struct->name);
        for(uint k=0;k<final_permutations[j].structures_duplicate_struct.size();k++){
          equivalent_permutations.push_back(final_permutations[j].structures_duplicate_struct[k]->name);
        }
        sstmp << "[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(equivalent_permutations,"\""),",") << "]";
        vcontent_json.push_back(sstmp.str()); sstmp.str("");
      }
      sscontent_json << aurostd::joinWDelimiter(vcontent_json,",");
      sscontent_json << "]";
      ss_output << "{" << sscontent_json.str() << "}" << endl;
    }
    //DX20190506 END

    // ---------------------------------------------------------------------------
    // print misfit results
    if(print_misfit){
      if(format=="text"){ //DX20190506
        ss_output << "Misfit values: " << endl; 
        stringstream ss_text;
        printResults(ss_text, same_species, final_permutations, "text");
        ss_output << ss_text.str();
      }
      else if(format=="json"){ //DX20190506
        ss_output.str(""); // need to clear content abbreviated content from above
        stringstream ss_json;
        printResults(ss_json, same_species, final_permutations, "json");
        ss_output << ss_json.str() << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // store unique atom decorations in vector
    for(uint j=0;j<final_permutations.size();j++){
      unique_permutations.push_back(final_permutations[j].structure_representative_struct->name);
    }

    // update oss
    oss << ss_output.str();
    return unique_permutations;
  }

// ***************************************************************************
// compare::compareMultipleStructures()
// ***************************************************************************
namespace compare {
  string compareMultipleStructures(const aurostd::xoption& vpflow, ostream& logstream){ //DX20190425 - changed name, more general

    // This function compares multiple structures (i.e., more than two).

    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = XPID + "compare::compareMultipleStructures():";
    ostringstream oss;
    //DX20200103 ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;

    if(LDEBUG){cerr << function_name << " BEGIN" << endl;}  //CO20200508
    
    // ---------------------------------------------------------------------------
    // instantiate XtalFinder calculator 
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        1,
        logstream);

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103
    xtal_finder.getOptions(vpflow, comparison_options);

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE_STRUCTURE::USAGE")) {
      //[CO20200624 - OBSOLETE]stringstream ss_usage;
      // material-type comparisons
      if(vpflow.flag("COMPARE_MATERIAL_DIRECTORY")){
        string usage_material_comparison="aflow --compare_materials -D <dir_path>";
        string options_material_comparison="[--usage] [--np=|--num_proc=<number>] [--optimize_match] [--no_scale_volume] [--ignore_symmetry] [--ignore_Wyckoff]";
        //[CO20200624 - OBSOLETE]init::ErrorOption(ss_usage,vpflow.getattachedscheme("COMPARE_STRUCTURE"),"compare::compareMultipleStructures()",aurostd::liststring2string(usage_material_comparison,options_material_comparison));
        init::ErrorOption(vpflow.getattachedscheme("COMPARE_STRUCTURE"),"compare::compareMultipleStructures()",aurostd::liststring2string(usage_material_comparison,options_material_comparison));
      }
      // structure-type comparisons
      else if(vpflow.flag("COMPARE_STRUCTURE_DIRECTORY")){
        string usage_structure_comparison="aflow --compare_structures -D <dir_path>";
        string options_structure_comparison="[--usage] [--np=|--num_proc=<number>] [--optimize_match] [--no_scale_volume] [--ignore_symmetry] [--ignore_Wyckoff] [--remove_duplicates|--remove_duplicate_compounds]";
        //[CO20200624 - OBSOLETE]init::ErrorOption(ss_usage,vpflow.getattachedscheme("COMPARE_STRUCTURE"),"compare::compareMultipleStructures()",aurostd::liststring2string(usage_structure_comparison,options_structure_comparison));
        init::ErrorOption(vpflow.getattachedscheme("COMPARE_STRUCTURE"),"compare::compareMultipleStructures()",aurostd::liststring2string(usage_structure_comparison,options_structure_comparison));
      }
      //[CO20200624 - OBSOLETE]return ss_usage.str();
      return "";
    }

    // ---------------------------------------------------------------------------
    // distinguish structures coming from directory or file
    string structures_source = ""; // "structure_list", "directory" or "file"

    // from list appended to command
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURE_LIST").empty()){
      structures_source = "structure_list";
    }
    // from directory
    //DX20190424 [OBSOLETE] if(vpflow.flag("COMPARE_MATERIAL_DIRECTORY") || vpflow.flag("COMPARE_STRUCTURE_DIRECTORY"))
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::DIRECTORY").empty())
    { //CO20200106 - patching for auto-indenting
      structures_source = "directory";
    }
    // from file
    //DX20190424 [OBSOLETE] if(vpflow.flag("COMPARE_MATERIAL_FILE") || vpflow.flag("COMPARE_STRUCTURE_FILE"))
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::FILE").empty())
    { //CO20200106 - patching for auto-indenting
      structures_source = "file";
    }

    // ---------------------------------------------------------------------------
    // FLAG: directory of structures to compare
    vector<string> file_list; //DX20190424
    string directory=".";
    string filename="";
    //DX20190424 START
    if(structures_source=="structure_list") {
      aurostd::string2tokens(vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURE_LIST"),file_list,",");
      message << "List of files to compare: " << aurostd::joinWDelimiter(file_list,",");
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    //DX20190424 END
    else if(structures_source=="directory") {
      directory=vpflow.getattachedscheme("COMPARE_STRUCTURE::DIRECTORY");
      if(!aurostd::FileExist(directory)) {
        message << "Unable to locate directory: " << directory << "." << endl;
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
        return oss.str();
      }
      message << "Comparison directory: " << directory;
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    // ---------------------------------------------------------------------------
    // FLAG: file of structures to compare
    else if(structures_source=="file") {
      filename=vpflow.getattachedscheme("COMPARE_STRUCTURE::FILE");
      if(!aurostd::FileExist(filename)) {
        message << "Unable to locate file: " << filename << "." << endl;
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
        return oss.str();
      }
      message << "Comparison file: " << filename;
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    else {
      message << "Need to specify location of structures to compare: -D <directory> or -F=<filename>." << endl;
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
    }

    /*
    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE_STRUCTURE::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_STRUCTURE::NP"));
      message << "OPTIONS: Using multiple threads; np = " << num_proc << ".";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
*/
    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool same_species=false;
    if(vpflow.flag("COMPARE_MATERIAL")){  // material comparisons (find duplicates) //DX20190429 - generalized
      same_species=true;
      message << "OPTIONS: Performing material type comparisons (comparing alike atomic species).";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    else if (vpflow.flag("COMPARE_STRUCTURE")){ // structure comparisons (find structure prototypes) //DX20190425 - generalized
      same_species=false;
      message << "OPTIONS: Performing structure type comparisons (any atomic species).";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
/*
    // ---------------------------------------------------------------------------
    // FLAG: misfit threshold //DX20201119
    double misfit_match_threshold = DEFAULT_XTALFINDER_MISFIT_MATCH;
    double misfit_family_threshold = DEFAULT_XTALFINDER_MISFIT_FAMILY;
    if(vpflow.flag("COMPARE_STRUCTURE::MISFIT_MATCH")) {
      misfit_match_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE_STRUCTURE::MISFIT_MATCH"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_MATCH",vpflow.getattachedscheme("COMPARE_STRUCTURE::MISFIT_MATCH"));
    }
    if(vpflow.flag("COMPARE_STRUCTURE::MISFIT_FAMILY")) {
      misfit_family_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE_STRUCTURE::MISFIT_FAMILY"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_FAMILY",vpflow.getattachedscheme("COMPARE_STRUCTURE::MISFIT_FAMILY"));
    }
    // match threshold must be less than family threshold
    if(misfit_match_threshold>misfit_family_threshold){
      message << "Matching misfit threshold must be less than the same family threshold:"
        << " misfit_match_threshold: " << misfit_match_threshold
        << " misfit_family_threshold: " << misfit_family_threshold;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ILLEGAL_);
    }
    message << "Misfit theshold for matched structures: " << misfit_match_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_MATCH << ")"; 
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Misfit theshold for structures in the same family: " << misfit_family_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_FAMILY << ")";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
*/
    // ---------------------------------------------------------------------------
    // FLAG: consider magnetic structure in comparison 
    bool magnetic_comparison=false;
    vector<string> magmoms_for_systems;
    if(vpflow.flag("COMPARE::MAGNETIC")){
      magnetic_comparison=true;
      string magnetic_info=vpflow.getattachedscheme("COMPARE::MAGNETIC");
      aurostd::string2tokens(magnetic_info,magmoms_for_systems,":");
      message << "OPTIONS: Including magnetic moment information in comparisons. Magnetic input detected for " << magmoms_for_systems.size() << " systems.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    if(magnetic_comparison){} //CO20200508 - keep it busy
/*
    // ---------------------------------------------------------------------------
    // FLAG: ICSD comparison - structure with minimum ICSD number as representative prototype
    // in general: smaller ICSD number = older = more reliable
    if(vpflow.flag("COMPARE_STRUCTURE::ICSD_COMPARISON")) {
      comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON",TRUE);
      message << "OPTIONS: Running on ICSD structures; use oldest ICSD number as representative prototype.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // optimize match (default: false)
    if(vpflow.flag("COMPARE_STRUCTURE::OPTIMIZE_MATCH")) {
      comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",TRUE);
      message << "OPTIONS: Finding optimal match; exploring all possible lattices and origins to find the best match (note: this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    if(vpflow.flag("COMPARE_STRUCTURE::NO_SCALE_VOLUME")) {
      comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
      message << "OPTIONS: Suppressing volume scaling; useful for distinguishing structures at different pressures."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    if(vpflow.flag("COMPARE_STRUCTURE::IGNORE_WYCKOFF")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    if(vpflow.flag("COMPARE_STRUCTURE::IGNORE_SYMMETRY")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    if(vpflow.flag("COMPARE_STRUCTURE::IGNORE_ENVIRONMENT_ANALYSIS")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",TRUE);
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: remove duplicate compounds (useful for non-biased statistics)
    if(vpflow.flag("COMPARE_STRUCTURE::REMOVE_DUPLICATE_COMPOUNDS")) {
      comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",TRUE);
      message << "OPTIONS: Remove duplicate compounds first, useful for non-biased prototype statistics."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: match unique structures to the AFLOW prototypes 
    if(vpflow.flag("COMPARE_STRUCTURE::MATCH_TO_AFLOW_PROTOS")) {
      comparison_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS",TRUE);
      message << "OPTIONS: Compare unique structures to the AFLOW prototypes."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: get AFLOW ANRL designation for unique structures
    if(vpflow.flag("COMPARE_STRUCTURE::ADD_AFLOW_PROTOTYPE_DESIGNATION")) {
      comparison_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION",TRUE);
      message << "OPTIONS: Cast unique structures into AFLOW standard designation."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: do not calculate unique atom decorations 
    if(vpflow.flag("COMPARE_STRUCTURE::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS")) {
      comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",FALSE);
      message << "OPTIONS: Do not calculate unique atom decorations."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: primitivize structures //DX20201005
    if(vpflow.flag("COMPARE_STRUCTURE::PRIMITIVIZE")) {
      comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE",TRUE);
      message << "OPTIONS: Converting all structures to a primitive representation.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: Minkowski reduction //DX20201005
    if(vpflow.flag("COMPARE_STRUCTURE::MINKOWSKI")) {
      comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI",TRUE);
      message << "OPTIONS: Performing Minkowski lattice reduction on all structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: Niggli reduction //DX20201005
    if(vpflow.flag("COMPARE_STRUCTURE::NIGGLI")) {
      comparison_options.flag("COMPARISON_OPTIONS::NIGGLI",TRUE);
      message << "OPTIONS: Performing Niggli lattice reduction on all structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
*/
    //DX20190425 - added print and screen only flag - START
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
/*
    // ---------------------------------------------------------------------------
    // FLAG: do not remove unmatched structures from the StructurePrototype Object
    // keeps results of each comparison
    if(vpflow.flag("COMPARE_STRUCTURE::KEEP_UNMATCHED")) {
      comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE);
    }
*/
    // ---------------------------------------------------------------------------
    // check if two-structure comparison
    if(file_list.size()==2){
      comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::STORE_COMPARISON_LOGS",TRUE);  //DX20200113 - fixed typo
    }
    //DX20190425 - added print and screen only flag - END

    //BETA // ===== FLAG: REMOVE DUPLICATE COMPOUNDS LAST ===== //
    //BETA bool remove_duplicate_compounds_last=false;
    //BETA if(vpflow.flag("COMPARE_STRUCTURE::REMOVE_DUPLICATE_COMPOUNDS_LAST")) {
    //BETA   remove_duplicate_compounds_last=true;
    //BETA   message << "OPTIONS: Remove duplicate compounds after performing structure-type comparisons, useful for non-biased prototype statistics."; 
    //BETA   pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    //BETA }


    // ---------------------------------------------------------------------------
    // load structures 
    vector<StructurePrototype> final_prototypes;
    if(structures_source=="structure_list") {
      //final_prototypes = compare::compareStructuresFromStructureList(file_list, magmoms_for_systems, oss, FileMESSAGE, num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
      final_prototypes = xtal_finder.compareStructuresFromStructureList(file_list, magmoms_for_systems, xtal_finder.num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
    }
    else if(structures_source=="directory") {
      //final_prototypes = compare::compareStructuresFromDirectory(directory, magmoms_for_systems, oss, FileMESSAGE, num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
      final_prototypes = xtal_finder.compareStructuresFromDirectory(directory, magmoms_for_systems, xtal_finder.num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
    }
    if(structures_source=="file") {
      //final_prototypes = compare::compareStructuresFromFile(filename, magmoms_for_systems, oss, FileMESSAGE, num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
      final_prototypes = xtal_finder.compareStructuresFromFile(filename, magmoms_for_systems, xtal_finder.num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
    }

    // ---------------------------------------------------------------------------
    // prepare JSON output
    stringstream ss_json;
    xtal_finder.printResults(ss_json, same_species, final_prototypes, "json");

    // ---------------------------------------------------------------------------
    // prepare TEXT (.out) output
    stringstream ss_out;
    xtal_finder.printResults(ss_out, same_species, final_prototypes, "txt");

    //DX20190429 - added screen only option - START
    // ---------------------------------------------------------------------------
    // write results to screen and return immediately (do not write file)
    if(screen_only){
      if(format=="json"){ return ss_json.str(); }
      // default is txt
      else { return ss_out.str(); }
    }
    //DX20190429 - added screen only option - END
    //DX20190429 - added format options - START
    // ---------------------------------------------------------------------------
    // if only two comparisons and text only, print mismatch information 
    if(file_list.size()==2){
      // return abbreviated results (i.e., misfit value along with match, same family, or no match text
      double final_misfit = AUROSTD_MAX_DOUBLE;
      if(final_prototypes[0].structure_misfits_duplicate.size()==1){
        final_misfit =  final_prototypes[0].structure_misfits_duplicate[0].misfit;
      }
      if(final_misfit <= xtal_finder.misfit_match && (final_misfit+1.0)> 1e-3){
        message << final_misfit << " : " << "MATCH" << endl;
      }
      else if(final_misfit > xtal_finder.misfit_match && final_misfit <= xtal_finder.misfit_family){
        message << final_misfit << " : " << "SAME FAMILY" << endl;
      }
      else if(final_misfit > xtal_finder.misfit_family && final_misfit <= 1.0){ 
        message << final_misfit << " : " << "NOT A MATCH" << endl;
      }
      else if(aurostd::isequal(final_misfit,AUROSTD_MAX_DOUBLE) || (final_misfit+1.0) < 1e-3){ 
        message << "UNMATCHABLE" << endl;
      }
      if(XHOST.QUIET){
        oss << message.str();
      }
      else {
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
      }
      if(print){
        if(final_prototypes[0].structure_misfits_duplicate.size()==1){
          oss << xtal_finder.printStructureMappingResults(final_prototypes[0].structure_misfits_duplicate[0],
              final_prototypes[0].structure_representative_struct->structure,
              final_prototypes[0].structures_duplicate_struct[0]->structure);
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
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
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
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }
    //DX20190429 - added format options - END
    return oss.str();
  }
}

// ***************************************************************************
// compare::getIsopointalPrototypes - returns corresponding prototype label
// ***************************************************************************
namespace compare {
  string isopointalPrototypes(istream& input, const aurostd::xoption& vpflow){ 

    string function_name = "compare::IsopointalPrototypes():";
    string usage="aflow --isopointal_prototypes|--get_isopointal_prototypes < POSCAR";
    string options="";

    // ---------------------------------------------------------------------------
    // load input structure
    xstructure xstr(input,IOAFLOW_AUTO);

    // ---------------------------------------------------------------------------
    // FLAG: catalog (htqc, anrl, or all)
    string catalog="all";
    if(vpflow.flag("ISOPOINTAL_PROTOTYPES::CATALOG")) {
      catalog=aurostd::tolower(vpflow.getattachedscheme("ISOPOINTAL_PROTOTYPES::CATALOG"));
      if(catalog!="htqc" && catalog!="anrl" && catalog!="all"){
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name, "Catalog/library can only be htqc, anrl, or all.",_INPUT_ILLEGAL_);
      }
    }

    // ---------------------------------------------------------------------------
    // get isopointal structures 
    // (calculates symmetry of input structure and grabs symmetrically similar prototypes)
    vector<string> isopointal_prototypes = compare::getIsopointalPrototypes(xstr, catalog);

    if(isopointal_prototypes.size()==0){
      return "no isopointal prototypes in AFLOW";
    }

    return aurostd::joinWDelimiter(isopointal_prototypes,",");
  }
}

// ***************************************************************************
// compare::getIsopointalPrototypes - returns corresponding prototype label
// ***************************************************************************
namespace compare {
  vector<string> getIsopointalPrototypes(xstructure& xstr, string& catalog){ 

    string function_name = "compare::getIsopointalPrototypes():";

    // ---------------------------------------------------------------------------
    // stoichiometry
    vector<uint> stoichiometry = xstr.GetReducedComposition(true);

    // ---------------------------------------------------------------------------
    // symmetry
    if(xstr.space_group_ITC<1 || xstr.space_group_ITC>230){ // don't recalculate symmetry if already calculated 
      double use_tol = SYM::defaultTolerance(xstr); //DX20200821
      xstr.SpaceGroup_ITC(use_tol, SG_SETTING_ANRL); //DX20200821 - added ANRL setting
    }
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
    compare::groupWyckoffPositions(xstr, grouped_Wyckoff_positions);

    // ---------------------------------------------------------------------------
    // extract isopointal prototypes from AFLOW 
    vector<string> vlabel;
    vector<uint> prototype_space_groups;
    vlabel = aflowlib::GetPrototypesBySymmetry(stoichiometry, xstr.space_group_ITC, grouped_Wyckoff_positions, prototype_space_groups, SG_SETTING_ANRL, catalog);

    return vlabel;

  }
}

/*
//DX20190314 - added new function - START
// ***************************************************************************
// compare::getMatchingPrototype - returns corresponding prototype label
// ***************************************************************************
namespace compare {
  vector<string> getMatchingPrototypes(xstructure& xstr, string& catalog){ 

    // Returns the matching prototype label, if any exists

    aurostd::xoption vpflow;
    vpflow.flag("COMPARE2PROTOTYPES",TRUE);

    // ---------------------------------------------------------------------------
    // specify catalog
    vpflow.flag("COMPARE2PROTOTYPES::CATALOG",TRUE); //DX20190329 - need to make scheme before attaching, otherwise it doesn't work
    vpflow.push_attached("COMPARE2PROTOTYPES::CATALOG",catalog); 

    // ---------------------------------------------------------------------------
    // do not calculate unique atom decorations
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

    return prototypes[0].structures_duplicate_names; // duplicates names are prototype labels 
  }
}
//DX20190314 - added new function - START
*/

//DX20190314 - added new function - START
// ***************************************************************************
// compare::getMatchingPrototype - returns corresponding prototype label
// ***************************************************************************
namespace compare {
  vector<string> getMatchingPrototypes(xstructure& xstr, string& catalog){ 

    // Returns the matching prototype label, if any exists

    aurostd::xoption vpflow;
    vpflow.flag("COMPARE2PROTOTYPES",TRUE);

    // ---------------------------------------------------------------------------
    // specify catalog
    vpflow.flag("COMPARE2PROTOTYPES::CATALOG",TRUE); //DX20190329 - need to make scheme before attaching, otherwise it doesn't work
    vpflow.push_attached("COMPARE2PROTOTYPES::CATALOG",catalog); 

    // ---------------------------------------------------------------------------
    // do not calculate unique atom decorations
    vpflow.flag("COMPARE2PROTOTYPES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS",TRUE);

    // ---------------------------------------------------------------------------
    // quiet output
    bool original_quiet = XHOST.QUIET;
    XHOST.QUIET=TRUE;

    // ---------------------------------------------------------------------------
    // compare structure to AFLOW prototypes
    XtalFinderCalculator xtal_finder;
    vector<StructurePrototype> prototypes = xtal_finder.compare2prototypes(xstr,vpflow);

    // ---------------------------------------------------------------------------
    // global quiet back to default
    XHOST.QUIET=original_quiet;

    //return prototypes[0].structures_duplicate_names; // duplicates names are prototype labels 
    vector<string> matching_prototypes;
    for(uint i=0;i<prototypes[0].structures_duplicate_struct.size();i++){
      matching_prototypes.push_back(prototypes[0].structures_duplicate_struct[i]->name);
    }
    return matching_prototypes; // duplicates names are prototype labels 
  }
}
//DX20190314 - added new function - START


//DX20190314 - added overloads for compare2prototypes - START
// ***************************************************************************
// compare::compare2prototypes - identifies corresponding protos 
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compare2prototypes(istream& input, const aurostd::xoption& vpflow, ostream& logstream){ 
    ofstream FileMESSAGE;
    return compare2prototypes(input, vpflow, FileMESSAGE, logstream);
  }
  vector<StructurePrototype> compare2prototypes(istream& input, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){ 

    // ---------------------------------------------------------------------------
    // load input structure
    xstructure xstr(input,IOAFLOW_AUTO);

    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        1,
        logstream);
    return xtal_finder.compare2prototypes(xstr,vpflow);
    //return compare2prototypes(xstr, vpflow, FileMESSAGE, logstream);
  }
}

// ***************************************************************************
// compare::printMatchingPrototypes - returns list of matching structures 
// ***************************************************************************
namespace compare {
  string printMatchingPrototypes(istream& input, const aurostd::xoption& vpflow){ 

    // ---------------------------------------------------------------------------
    // load input structure
    xstructure xstr(input,IOAFLOW_AUTO);

    XtalFinderCalculator xtal_finder;
    return xtal_finder.printMatchingPrototypes(xstr,vpflow);

  }
}

/*
// ***************************************************************************
// compare::printMatchingPrototypes - returns list of matching structures 
// ***************************************************************************
namespace compare {
  string printMatchingPrototypes(xstructure& xstr, const aurostd::xoption& vpflow){ 

    //DX20190425 - added print flag - START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    string format = "both";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      format = "txt";
    }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      format = "json";
    }
    //DX20190425 - added print flag - END

    //DX20190425 - added screen only flag - START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    bool screen_only = false;
    if(vpflow.flag("COMPARE2PROTOTYPES::SCREEN_ONLY")) {
      screen_only=true;
    }
    //DX20190425 - added screen only flag - END

    vector<StructurePrototype> prototypes = compare2prototypes(xstr,vpflow);

    // ---------------------------------------------------------------------------
    // print results 
    //DX20190509 [OBSOLETE-moved down] stringstream ss_out;
    bool same_species = false; //default for prototypes
    stringstream ss_json;
    compare::printResults(ss_json, same_species, prototypes, "json");

    stringstream ss_out;
    compare::printResults(ss_out, same_species, prototypes, "txt");

    //DX20190429 - added screen only option - START
    // ---------------------------------------------------------------------------
    // write results to screen and return immediately (do not write file)
    if(screen_only){
      if(format=="json"){ return ss_json.str(); }
      // default is txt
      else { return ss_out.str(); }
    }
    //DX20190429 - added screen only option - END

    if(format=="json"){ return ss_json.str(); }
    return ss_out.str();
  }
}
*/

// ***************************************************************************
// XtalFinderCalculator::printMatchingPrototypes - returns list of matching structures 
// ***************************************************************************
  string XtalFinderCalculator::printMatchingPrototypes(xstructure& xstr, const aurostd::xoption& vpflow){ 

    //DX20190425 - added print flag - START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    string format = "both";
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      format = "txt";
    }
    else if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      format = "json";
    }
    //DX20190425 - added print flag - END

    //DX20190425 - added screen only flag - START
    // ---------------------------------------------------------------------------
    // FLAG: print format
    bool screen_only = false;
    if(vpflow.flag("COMPARE2PROTOTYPES::SCREEN_ONLY")) {
      screen_only=true;
    }
    //DX20190425 - added screen only flag - END

    vector<StructurePrototype> prototypes = compare2prototypes(xstr,vpflow);

    // ---------------------------------------------------------------------------
    // print results 
    //DX20190509 [OBSOLETE-moved down] stringstream ss_out;
    bool same_species = false; //default for prototypes
    stringstream ss_json;
    printResults(ss_json, same_species, prototypes, "json");

    stringstream ss_out;
    printResults(ss_out, same_species, prototypes, "txt");

    //DX20190429 - added screen only option - START
    // ---------------------------------------------------------------------------
    // write results to screen and return immediately (do not write file)
    if(screen_only){
      if(format=="json"){ return ss_json.str(); }
      // default is txt
      else { return ss_out.str(); }
    }
    //DX20190429 - added screen only option - END

    if(format=="json"){ return ss_json.str(); }
    return ss_out.str();
  }

/*
//DX20190314 - added overloads for compare2prototypes - END
// ***************************************************************************
// compare::compare2prototypes - identifies corresponding protos 
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compare2prototypes(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream){ 
    ofstream FileMESSAGE;
    return compare2prototypes(xstrIN, vpflow, FileMESSAGE, logstream);
  }
  vector<StructurePrototype> compare2prototypes(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){ 
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    string function_name = XPID + "compare::compare2prototypes():";
    ostringstream oss;
    //DX20200103 ostream& logstream = cout;
    bool quiet = false;
    stringstream message;
    //DX20200226 [OBSOLETE] ofstream FileMESSAGE;

    string directory="";

    string usage="aflow --compare2protos|--compare2prototypes < POSCAR";
    string options="";

    xstructure xstr = xstrIN; //DX20200226 - copy

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103

    // ---------------------------------------------------------------------------
    // FLAG: misfit threshold //DX20201119
    double misfit_match_threshold = DEFAULT_XTALFINDER_MISFIT_MATCH;
    double misfit_family_threshold = DEFAULT_XTALFINDER_MISFIT_FAMILY;
    if(vpflow.flag("COMPARE2PROTOTYPES::MISFIT_MATCH")) {
      misfit_match_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE2PROTOTYPES::MISFIT_MATCH"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_MATCH",vpflow.getattachedscheme("COMPARE2PROTOTYPES::MISFIT_MATCH"));
    }
    if(vpflow.flag("COMPARE2PROTOTYPES::MISFIT_FAMILY")) {
      misfit_family_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE2PROTOTYPES::MISFIT_FAMILY"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_FAMILY",vpflow.getattachedscheme("COMPARE2PROTOTYPES::MISFIT_FAMILY"));
    }
    // match threshold must be less than family threshold
    if(misfit_match_threshold>misfit_family_threshold){
      message << "Matching misfit threshold must be less than the same family threshold:"
        << " misfit_match_threshold: " << misfit_match_threshold
        << " misfit_family_threshold: " << misfit_family_threshold;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ILLEGAL_);
    }
    message << "Misfit theshold for matched structures: " << misfit_match_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_MATCH << ")"; 
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Misfit theshold for structures in the same family: " << misfit_family_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_FAMILY << ")";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

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
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ILLEGAL_);
      }
      message << "OPTIONS: Catalog/library (htqc, anrl, or all): " << catalog << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // same spaces 
    bool same_species=false; //compare2prototype: by definition, want to compare backbone structure, i.e., ignore species

    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    if(vpflow.flag("COMPARE2PROTOTYPES::NO_SCALE_VOLUME")) {
      comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    if(vpflow.flag("COMPARE2PROTOTYPES::IGNORE_WYCKOFF")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    if(vpflow.flag("COMPARE2PROTOTYPES::IGNORE_SYMMETRY")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    if(vpflow.flag("COMPARE2PROTOTYPES::IGNORE_ENVIRONMENT_ANALYSIS")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",TRUE);
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: do not calculate unique atom decorations
    if(vpflow.flag("COMPARE2PROTOTYPES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS")) {
      comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",FALSE);
      message << "OPTIONS: Do not calculate unique atom decorations."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // single round of comparisons 
    comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);

    vector<StructurePrototype> all_structures;

    //DX20190314 [OBSOLETE] // ---------------------------------------------------------------------------
    //DX20190314 [OBSOLETE] // load input structure
    //DX20190314 [OBSOLETE] xstructure xstr(input,IOAFLOW_AUTO);

    StructurePrototype input_structure;
    input_structure.structure_representative = xstr;
    input_structure.structure_representative_name = "input geometry";
    input_structure.stoichiometry = xstr.GetReducedComposition(); //preserves the stoich order for the structure
    input_structure.elements = xstr.GetElements(true,true); // true: clean names
    input_structure.structure_representative_compound = pflow::prettyPrintCompound(input_structure.elements,input_structure.stoichiometry,no_vrt,false,txt_ft);
    input_structure.structure_representative_generated = true;
    stringstream ss_input; ss_input << xstr;
    input_structure.structure_representative_source = ss_input.str();
    input_structure.structure_representative_relaxation_step = 0; //DX20200429 input is assumed to be unrelaxed
    all_structures.push_back(input_structure);

    // ---------------------------------------------------------------------------
    // symmetry
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && (xstr.space_group_ITC<1 || xstr.space_group_ITC>230)){ //DX20190829 - don't recalculate symmetry if already calculated //DX20191220 - put range instead of ==0
      message << "Calculating the symmetry of the input structure.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      uint one_proc=1;
      compare::calculateSymmetries(all_structures,one_proc); 
      message << "Symmetry calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }
    else if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && xstr.space_group_ITC>=1 && xstr.space_group_ITC<=230){ //DX20191220 - put range instead of !=0
      for(uint i=0;i<all_structures.size();i++){
        all_structures[i].space_group = all_structures[i].structure_representative.space_group_ITC;
        vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
        compare::groupWyckoffPositions(all_structures[i].structure_representative, grouped_Wyckoff_positions);
        all_structures[i].grouped_Wyckoff_positions=grouped_Wyckoff_positions;
      }
    }

    if(LDEBUG) {
      cerr << function_name << " Wyckoff positions of input structure:" << endl;
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
      cerr << function_name << " Wyckoff positions of input structure:" << endl;
      for(uint i=0;i<all_structures[0].grouped_Wyckoff_positions.size();i++){
        cerr << all_structures[0].grouped_Wyckoff_positions[i] << endl;
      }
    }

    vector<string> vlabel;
    vector<uint> prototype_space_groups;

    // ---------------------------------------------------------------------------
    // find prototypes based on stoichiometry only
    if(comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY")){
      message << "Load prototypes with the same stoichiometry as the input.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      vlabel = aflowlib::GetPrototypesByStoichiometry(stoichiometry, catalog);
    }
    // find prototypes based on stoichiometry and space group only
    else if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF")){
      message << "Load prototypes with the same stoichiometry and space group as the input.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      vector<GroupedWyckoffPosition> empty_Wyckoff_positions;
      vlabel = aflowlib::GetPrototypesBySymmetry(stoichiometry, space_group_num, empty_Wyckoff_positions, prototype_space_groups, SG_SETTING_ANRL, catalog);
    }
    // find prototypes based on stoichiometry, space group, and Wyckoff positions only (recommended/default)
    else {
      message << "Load prototypes with the same stoichiometry, space group, and Wyckoff positions as the input.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      vlabel = aflowlib::GetPrototypesBySymmetry(stoichiometry, space_group_num, grouped_Wyckoff_positions, prototype_space_groups, SG_SETTING_ANRL, catalog);
    }
    message << "Potential compatible prototypes: " << vlabel.size() << " (" << aurostd::joinWDelimiter(vlabel,",") << ").";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    if(LDEBUG) {
      cerr << function_name << " Wyckoff positions of input structure:" << endl;
      for(uint i=0;i<all_structures[0].grouped_Wyckoff_positions.size();i++){
        cerr << all_structures[0].grouped_Wyckoff_positions[i] << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // load compatible aflow prototypes

    //DX20190830 - to avoid multiple threads being spun-up (here and in aflow_xproto.cpp), turn of aflow_pthreads
    uint uint_backup=AFLOW_PTHREADS::MAX_PTHREADS;
    AFLOW_PTHREADS::MAX_PTHREADS=1;

    compare::addAFLOWPrototypes2StructurePrototypeVector(all_structures, vlabel); 
    if(LDEBUG) {
      cerr << function_name << " Wyckoff positions of input structure:" << endl;
      for(uint i=0;i<all_structures[0].grouped_Wyckoff_positions.size();i++){
        cerr << all_structures[0].grouped_Wyckoff_positions[i] << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // group into objects based on stoichiometry and symmetry (Pearson and space group)
    message << "Grouping sets of comparisons.";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> comparison_schemes = compare::groupStructurePrototypes(all_structures, 
        same_species, 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320 - added environment angles
        false); //DX20200103 - condensed booleans to xoptions 
    message << "Number of comparison groups: " << comparison_schemes.size() << ".";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // compare structures 
    message << "Running comparisons ...";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes = compare::runComparisonScheme(comparison_schemes, same_species, num_proc, comparison_options, oss, FileMESSAGE, quiet, logstream); //DX20200103 - condensed booleans to xoptions 

    message << "Comparisons complete ...";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);

    AFLOW_PTHREADS::MAX_PTHREADS = uint_backup; //DX20190830 - set back to original setting

    // ---------------------------------------------------------------------------
    // return if there are no similar structures
    if(final_prototypes.size()==0){ return final_prototypes; } //DX20190314 originally : return oss.str()

    comparison_schemes.clear();

    // ---------------------------------------------------------------------------
    // get unique atom decorations
    if(comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS")){
      message << "Identifying unique atom decorations for representative structures ...";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

      for(uint i=0;i<final_prototypes.size();i++){
        // check if xstructure is generated; if not, make it
        if(!final_prototypes[i].structure_representative_generated){
          if(!compare::generateStructure(final_prototypes[i].structure_representative_name,final_prototypes[i].structure_representative_source,final_prototypes[i].structure_representative_relaxation_step,final_prototypes[i].structure_representative,oss)){ //DX20200429
            message << "Could not generate structure (" << final_prototypes[i].structure_representative_name << ").";
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
          }
        }        
        if(LDEBUG){ //DX20190601 - added LDEBUG
          cerr << "Finding unique atom decorations for " << final_prototypes[i].structure_representative_name << ".";
        }        
        vector<StructurePrototype> final_permutations = compare::comparePermutations(final_prototypes[i],num_proc,comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH"),oss,FileMESSAGE); //DX20200103 - condensed booleans to xoptions
        for(uint j=0;j<final_permutations.size();j++){
          vector<string> tmp_permutations; 
          tmp_permutations.push_back(final_permutations[j].structure_representative_name); //push back representative permutation
          for(uint d=0;d<final_permutations[j].structures_duplicate_names.size();d++){ tmp_permutations.push_back(final_permutations[j].structures_duplicate_names[d]); } //push back equivalent permutations
          final_prototypes[i].atom_decorations_equivalent.push_back(tmp_permutations);
        }
      }
      message << "Unique atom decorations found.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }

    // ---------------------------------------------------------------------------
    // print results 
    //DX20190314 -moved to overload [OBSOLETE] stringstream ss_out;
    //DX20190314 -moved to overload [OBSOLETE] compare::printResults(ss_out, same_species, final_prototypes);
    //DX20190314 -moved to overload [OBSOLETE]
    //DX20190314 -moved to overload [OBSOLETE] oss << ss_out.str();
    //DX20190314 -moved to overload [OBSOLETE] return oss.str();
    return final_prototypes; //DX20190314 - new return type
  }
}
*/

// ***************************************************************************
// XtalFinderCalculator::compare2prototypes - identifies corresponding protos 
// ***************************************************************************
  vector<StructurePrototype> XtalFinderCalculator::compare2prototypes(
      const xstructure& xstrIN,
      const aurostd::xoption& vpflow){ 
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    string function_name = XPID + "XtalFinderCalculator::compare2prototypes():";
    stringstream message;
    bool quiet = false;

    string directory="";

    string usage="aflow --compare2protos|--compare2prototypes < POSCAR";
    string options="";

    xstructure xstr = xstrIN; //DX20200226 - copy

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103
    getOptions(vpflow, comparison_options);
    /*
    // ---------------------------------------------------------------------------
    // FLAG: misfit threshold //DX20201119
    //double misfit_match_threshold = DEFAULT_XTALFINDER_MISFIT_MATCH;
    //double misfit_family_threshold = DEFAULT_XTALFINDER_MISFIT_FAMILY;
    if(vpflow.flag("COMPARE2PROTOTYPES::MISFIT_MATCH")) {
      misfit_match = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE2PROTOTYPES::MISFIT_MATCH"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_MATCH",vpflow.getattachedscheme("COMPARE2PROTOTYPES::MISFIT_MATCH"));
    }
    if(vpflow.flag("COMPARE2PROTOTYPES::MISFIT_FAMILY")) {
      misfit_family = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE2PROTOTYPES::MISFIT_FAMILY"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_FAMILY",vpflow.getattachedscheme("COMPARE2PROTOTYPES::MISFIT_FAMILY"));
    }
    // match threshold must be less than family threshold
    if(misfit_match>misfit_family){
      message << "Matching misfit threshold must be less than the same family threshold:"
        << " misfit_match_threshold: " << misfit_match
        << " misfit_family_threshold: " << misfit_family;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ILLEGAL_);
    }
    message << "Misfit theshold for matched structures: " << misfit_match << " (default: " << DEFAULT_XTALFINDER_MISFIT_MATCH << ")"; 
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Misfit theshold for structures in the same family: " << misfit_family << " (default: " << DEFAULT_XTALFINDER_MISFIT_FAMILY << ")";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE2PROTOTYPES::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE2PROTOTYPES::NP"));
    }
*/
    // ---------------------------------------------------------------------------
    // FLAG: catalog (htqc, anrl, or all)
    string catalog="all";
    if(vpflow.flag("COMPARE2PROTOTYPES::CATALOG")) {
      catalog=aurostd::tolower(vpflow.getattachedscheme("COMPARE2PROTOTYPES::CATALOG"));
      if(catalog!="htqc" && catalog!="anrl" && catalog!="all"){
        message << "Catalog/library can only be htqc, anrl, or all.";     
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ILLEGAL_);
      }
      message << "OPTIONS: Catalog/library (htqc, anrl, or all): " << catalog << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // same spaces 
    bool same_species=false; //compare2prototype: by definition, want to compare backbone structure, i.e., ignore species

    /*
    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    if(vpflow.flag("COMPARE2PROTOTYPES::NO_SCALE_VOLUME")) {
      comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    if(vpflow.flag("COMPARE2PROTOTYPES::IGNORE_WYCKOFF")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    if(vpflow.flag("COMPARE2PROTOTYPES::IGNORE_SYMMETRY")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    if(vpflow.flag("COMPARE2PROTOTYPES::IGNORE_ENVIRONMENT_ANALYSIS")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",TRUE);
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: do not calculate unique atom decorations
    if(vpflow.flag("COMPARE2PROTOTYPES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS")) {
      comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",FALSE);
      message << "OPTIONS: Do not calculate unique atom decorations."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }
*/
    // ---------------------------------------------------------------------------
    // single round of comparisons 
    comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);

    //
    // add structure to container
    stringstream ss_input; ss_input << xstr;
    addStructureToContainer(xstr, "input geometry", ss_input.str(), 0, false);

    vector<StructurePrototype> all_structures;

    //DX20190314 [OBSOLETE] // ---------------------------------------------------------------------------
    //DX20190314 [OBSOLETE] // load input structure
    //DX20190314 [OBSOLETE] xstructure xstr(input,IOAFLOW_AUTO);

    // ---------------------------------------------------------------------------
    // symmetry
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && (xstr.space_group_ITC<1 || xstr.space_group_ITC>230)){ //DX20190829 - don't recalculate symmetry if already calculated //DX20191220 - put range instead of ==0
      message << "Calculating the symmetry of the input structure.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      uint one_proc=1;
      calculateSymmetries(one_proc); 
      message << "Symmetry calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    }
    else if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && xstr.space_group_ITC>=1 && xstr.space_group_ITC<=230){ //DX20191220 - put range instead of !=0
      for(uint i=0;i<structure_containers.size();i++){
        structure_containers[i].Pearson = "xX";
        structure_containers[i].space_group = 0;
        vector<GroupedWyckoffPosition> vGWyckoffPos_tmp;
        structure_containers[i].grouped_Wyckoff_positions = vGWyckoffPos_tmp;
      }
    }

    if(LDEBUG) {
      cerr << function_name << " Wyckoff positions of input structure:" << endl;
      for(uint i=0;i<structure_containers[0].grouped_Wyckoff_positions.size();i++){
        cerr << structure_containers[0].grouped_Wyckoff_positions[i] << endl;
      }
    }

    // ---------------------------------------------------------------------------
    // variables to pass to GetPrototypes functions
    vector<uint> stoichiometry = structure_containers[0].stoichiometry;
    std::sort(stoichiometry.begin(),stoichiometry.end()); // order stoichiometry, so we can match to AFLOW prototypes more easily
    uint space_group_num = structure_containers[0].space_group;
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions = structure_containers[0].grouped_Wyckoff_positions;

    vector<string> vlabel;
    vector<uint> prototype_space_groups;

    // ---------------------------------------------------------------------------
    // find prototypes based on stoichiometry only
    if(comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY")){
      message << "Load prototypes with the same stoichiometry as the input.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      vlabel = aflowlib::GetPrototypesByStoichiometry(stoichiometry, catalog);
    }
    // find prototypes based on stoichiometry and space group only
    else if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF")){
      message << "Load prototypes with the same stoichiometry and space group as the input.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      vector<GroupedWyckoffPosition> empty_Wyckoff_positions;
      vlabel = aflowlib::GetPrototypesBySymmetry(stoichiometry, space_group_num, empty_Wyckoff_positions, prototype_space_groups, SG_SETTING_ANRL, catalog);
    }
    // find prototypes based on stoichiometry, space group, and Wyckoff positions only (recommended/default)
    else {
      message << "Load prototypes with the same stoichiometry, space group, and Wyckoff positions as the input.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      vlabel = aflowlib::GetPrototypesBySymmetry(stoichiometry, space_group_num, grouped_Wyckoff_positions, prototype_space_groups, SG_SETTING_ANRL, catalog);
    }
    message << "Potential compatible prototypes: " << vlabel.size() << " (" << aurostd::joinWDelimiter(vlabel,",") << ").";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // load compatible aflow prototypes

    //DX20190830 - to avoid multiple threads being spun-up (here and in aflow_xproto.cpp), turn of aflow_pthreads
    uint uint_backup=AFLOW_PTHREADS::MAX_PTHREADS;
    AFLOW_PTHREADS::MAX_PTHREADS=1;

    //compare::addAFLOWPrototypes2StructurePrototypeVector(all_structures, vlabel); 
    addAFLOWPrototypes2container(vlabel); 

    // ---------------------------------------------------------------------------
    // group into objects based on stoichiometry and symmetry (Pearson and space group)
    message << "Grouping sets of comparisons.";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    vector<StructurePrototype> comparison_schemes = groupStructurePrototypes( 
        same_species, 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320 - added environment angles
        false); //DX20200103 - condensed booleans to xoptions
    message << "Number of comparison groups: " << comparison_schemes.size() << ".";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // compare structures 
    message << "Running comparisons ...";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes = runComparisonScheme(comparison_schemes, same_species, num_proc, comparison_options, quiet); //DX20200103 - condensed booleans to xoptions 

    message << "Comparisons complete ...";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);

    AFLOW_PTHREADS::MAX_PTHREADS = uint_backup; //DX20190830 - set back to original setting

    // ---------------------------------------------------------------------------
    // return if there are no similar structures
    if(final_prototypes.size()==0){ return final_prototypes; } //DX20190314 originally : return oss.str()

    comparison_schemes.clear();

    // ---------------------------------------------------------------------------
    // get unique atom decorations
    if(comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS")){
      message << "Identifying unique atom decorations for representative structures ...";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

      for(uint i=0;i<final_prototypes.size();i++){
        // check if xstructure is generated; if not, make it
        if(!final_prototypes[i].structure_representative_struct->is_structure_generated){
          if(!compare::generateStructure(final_prototypes[i].structure_representative_struct->name,final_prototypes[i].structure_representative_struct->source,final_prototypes[i].structure_representative_struct->relaxation_step,final_prototypes[i].structure_representative_struct->structure,*p_oss)){ //DX20200429
            message << "Could not generate structure (" << final_prototypes[i].structure_representative_struct->name << ").";
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
          }
        }        
        if(LDEBUG){ //DX20190601 - added LDEBUG
          cerr << "Finding unique atom decorations for " << final_prototypes[i].structure_representative_struct->name << ".";
        }        
        XtalFinderCalculator xtal_finder_permutations;
        vector<StructurePrototype> final_permutations = xtal_finder_permutations.comparePermutations(final_prototypes[i],num_proc,comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH"));
        for(uint j=0;j<final_permutations.size();j++){
          vector<string> tmp_permutations; 
          tmp_permutations.push_back(final_permutations[j].structure_representative_struct->name); //push back representative permutation
          for(uint d=0;d<final_permutations[j].structures_duplicate_struct.size();d++){ tmp_permutations.push_back(final_permutations[j].structures_duplicate_struct[d]->name); } //push back equivalent permutations
          final_prototypes[i].atom_decorations_equivalent.push_back(tmp_permutations);
        }
      }
      message << "Unique atom decorations found.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    }

    // ---------------------------------------------------------------------------
    // print results 
    //DX20190314 -moved to overload [OBSOLETE] stringstream ss_out;
    //DX20190314 -moved to overload [OBSOLETE] compare::printResults(ss_out, same_species, final_prototypes);
    //DX20190314 -moved to overload [OBSOLETE]
    //DX20190314 -moved to overload [OBSOLETE] oss << ss_out.str();
    //DX20190314 -moved to overload [OBSOLETE] return oss.str();
    return final_prototypes; //DX20190314 - new return type
  }

// ***************************************************************************
// compare::isMatchingStructureInDatabase - boolean if match found in database 
// ***************************************************************************
namespace compare {
  // load input structure
  bool isMatchingStructureInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream){ ofstream FileMESSAGE; return isMatchingStructureInDatabase(xstrIN, vpflow, FileMESSAGE, logstream); }
  bool isMatchingStructureInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){

    // ---------------------------------------------------------------------------
    // main compare2database() function
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        1,
        logstream);
    vector<StructurePrototype> final_prototypes = xtal_finder.compare2database(xstrIN, vpflow);

    // ---------------------------------------------------------------------------
    // safety against bad input geometry files
    if(final_prototypes.empty()){
      string function_name = "compare::isMatchingStructureInDatabase():";
      stringstream message;
      message << "The input geometry file is invalid (could not be read, corrupt, etc.); it could not be compared to the database.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_FILE_ERROR_);
    }

    // ---------------------------------------------------------------------------
    // database contains input structure
    //if(final_prototypes[0].structures_duplicate_names.size()){
    if(final_prototypes[0].structures_duplicate_struct.size()){
      return true;
    }

    // ---------------------------------------------------------------------------
    // database DOESN'T contain equivalent structure in input
    return false;

  }
}

// ***************************************************************************
// compare::matchingStructuresInDatabase - returns matching structures in database 
// ***************************************************************************
namespace compare {
  // load input structure
  vector<matching_structure> matchingStructuresInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream){ ofstream FileMESSAGE; return matchingStructuresInDatabase(xstrIN, vpflow, FileMESSAGE, logstream); }
  vector<matching_structure> matchingStructuresInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){

    // ---------------------------------------------------------------------------
    // main compare2database() function
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        1,
        logstream);
    vector<StructurePrototype> final_prototypes = xtal_finder.compare2database(xstrIN, vpflow);


    vector<matching_structure> matched_database_structures;

    // ---------------------------------------------------------------------------
    // database DOESN'T contain equivalent structure to input
    //if(!final_prototypes[0].structures_duplicate_names.size()){
    if(!final_prototypes[0].structures_duplicate_struct.size()){
      return matched_database_structures;
    }

    // ---------------------------------------------------------------------------
    // return equivalent structures to input
    //for(uint i=0;i<final_prototypes[0].structures_duplicate_names.size();i++){
    for(uint i=0;i<final_prototypes[0].structures_duplicate_struct.size();i++){
      matching_structure database_entry;
      database_entry.name = final_prototypes[0].structures_duplicate_struct[i]->name;
      database_entry.misfit = final_prototypes[0].structure_misfits_duplicate[i].misfit;
      matched_database_structures.push_back(database_entry);
    }

    return matched_database_structures;

  }
}

/*
// ***************************************************************************
// compare::compare2database - compares database 
// ***************************************************************************
namespace compare {
  // load input structure
  vector<StructurePrototype> compare2database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream){ ofstream FileMESSAGE; return compare2database(xstrIN, vpflow, FileMESSAGE, logstream); }
  vector<StructurePrototype> compare2database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    vector<StructurePrototype> final_prototypes; //DX20200225
    xstructure xstr = xstrIN; //copy //DX20200225

    string function_name = XPID + "compare::compare2database():";
    string directory = "";
    ostringstream oss;
    stringstream message;
    bool quiet = false;
    //DX20200103 ostream& logstream = cout;
    //DX20200225 ofstream FileMESSAGE;

    string usage="aflow --compare2database < POSCAR";
    string options="";

    vector<string> tokens,sub_tokens;
    vector<string> matchbook; //aflux - filter/get properties
    vector<string> schema; //get metadata of properties (e.g., units)
    vector<string> property_units;

    bool same_species = true;

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103

    // ---------------------------------------------------------------------------
    // FLAG: misfit threshold //DX20201119
    double misfit_match_threshold = DEFAULT_XTALFINDER_MISFIT_MATCH;
    double misfit_family_threshold = DEFAULT_XTALFINDER_MISFIT_FAMILY;
    if(vpflow.flag("COMPARE2DATABASE::MISFIT_MATCH")) {
      misfit_match_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE2DATABASE::MISFIT_MATCH"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_MATCH",vpflow.getattachedscheme("COMPARE2DATABASE::MISFIT_MATCH"));
    }
    if(vpflow.flag("COMPARE2DATABASE::MISFIT_FAMILY")) {
      misfit_family_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE2DATABASE::MISFIT_FAMILY"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_FAMILY",vpflow.getattachedscheme("COMPARE2DATABASE::MISFIT_FAMILY"));
    }
    // match threshold must be less than family threshold
    if(misfit_match_threshold>misfit_family_threshold){
      message << "Matching misfit threshold must be less than the same family threshold:"
        << " misfit_match_threshold: " << misfit_match_threshold
        << " misfit_family_threshold: " << misfit_family_threshold;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ILLEGAL_);
    }
    message << "Misfit theshold for matched structures: " << misfit_match_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_MATCH << ")"; 
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Misfit theshold for structures in the same family: " << misfit_family_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_FAMILY << ")";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // single round of comparisons 
    comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);

    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE2DATABASE::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE2DATABASE::NP"));
    }

    // ---------------------------------------------------------------------------
    // FLAG: optimize match
    if(vpflow.flag("COMPARE2DATABASE::OPTIMIZE_MATCH")) {
      comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",TRUE); //DX20200225 - was missing
      message << "OPTIONS: Finding optimal match; exploring all possible lattices and origins to find the best match (note: this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    if(vpflow.flag("COMPARE2DATABASE::NO_SCALE_VOLUME")) {
      comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
      message << "OPTIONS: Suppressing volume scaling; useful for distinguishing structures at different pressures."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    if(vpflow.flag("COMPARE2DATABASE::IGNORE_WYCKOFF")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    if(vpflow.flag("COMPARE2DATABASE::IGNORE_SYMMETRY")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    if(vpflow.flag("COMPARE2DATABASE::IGNORE_ENVIRONMENT_ANALYSIS")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",TRUE);
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: primitivize structures //DX20201005
    if(vpflow.flag("COMPARE2DATABASE::PRIMITIVIZE")) {
      comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE",TRUE);
      message << "OPTIONS: Converting all structures to a primitive representation.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: Minkowski reduction //DX20201005
    if(vpflow.flag("COMPARE2DATABASE::MINKOWSKI")) {
      comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI",TRUE);
      message << "OPTIONS: Performing Minkowski lattice reduction on all structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: Niggli reduction //DX20201005
    if(vpflow.flag("COMPARE2DATABASE::NIGGLI")) {
      comparison_options.flag("COMPARISON_OPTIONS::NIGGLI",TRUE);
      message << "OPTIONS: Performing Niggli lattice reduction on all structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool structure_comparison=false;
    if(vpflow.flag("COMPARE2DATABASE::STRUCTURE")) {
      structure_comparison=true;
      same_species = false;
      message << "OPTIONS: Structure-type comparison, i.e., ignore atomic species."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
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
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify the geometry file to grab (orig, relax1, relax2, static, bands, POSCAR, CONTCAR)
    //DX TODO 
    string geometry_file = "";
    if(vpflow.flag("COMPARE2DATABASE::GEOMETRY_FILE")) {
      geometry_file = vpflow.getattachedscheme("COMPARE2DATABASE::GEOMETRY_FILE");
      message << "OPTIONS: Structure type (POSCAR.orig, POSCAR.relax1, POSCAR.relax2, CONTCAR.relax1, ...): " << geometry_file << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify the relaxation step to grab (orig, relax1, relax2, static, bands, POSCAR, CONTCAR)
    //DX TODO 
    uint relaxation_step = _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_;
    bool load_most_relaxed_structure_only = true; 
    if(vpflow.flag("COMPARE2DATABASE::RELAXATION_STEP")) {
      if(aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP")), "orig") || 
          vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP") == "0"){
        relaxation_step = _COMPARE_DATABASE_GEOMETRY_ORIGINAL_;
        load_most_relaxed_structure_only = false;
      }
      if(aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP")), "relax1") || 
          aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP")), "middle_relax") || 
          vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP") == "1"){
        relaxation_step = _COMPARE_DATABASE_GEOMETRY_RELAX1_;
        load_most_relaxed_structure_only = false;
      }
      message << "OPTIONS: Relaxation step (0=original, 1=relax1, 2=most_relaxed): " << relaxation_step << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // should not grab properties and compare structures other than the most relaxed structure 
    if(property_list.size() && relaxation_step != _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_){
      string relaxation_name = "";
      if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_ORIGINAL_){ relaxation_name = "original"; }
      else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_RELAX1_){ relaxation_name = "relax1"; }
      message << "The " << relaxation_name << " structures will be extracted; the properties will not correspond to these structures. Proceed with caution."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: catalog (icsd, lib1, lib2, lib3, ...)
    string catalog = "";
    string catalog_summons = "";
    if(vpflow.flag("COMPARE2DATABASE::CATALOG")) {
      catalog = aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::CATALOG")); //DX20190329 -- added tolower
      if(catalog != "all"){ //DX20190329 - added if-statement since AFLUX doesn't use "all"
        catalog_summons = "catalog(\'" + catalog + "\')";
        matchbook.push_back(catalog_summons);
      } //DX20190329 - added if-statement since AFLUX doesn't use "all"
      message << "OPTIONS: Catalog/library (icsd, lib1, lib2, lib3, ...): " << catalog << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    //DX20190508 - added keep unmatched - START
    // ---------------------------------------------------------------------------
    // FLAG: do not remove unmatched structures from the StructurePrototype Object
    // keeps results of each comparison
    if(vpflow.flag("COMPARE2DATABASE::KEEP_UNMATCHED")) {
      comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE);
    }
    //DX20190508 - added keep unmatched - END

    // ---------------------------------------------------------------------------
    // fix species (remove pseudopotentials, etc.) 
    string species_str = aurostd::joinWDelimiter(xstr.species, ""); //DX20200212 
    vector<string> vspecies = aurostd::getElements(species_str); //DX20200212
    xstr.species = aurostd::vector2deque(vspecies); //DX20200212 - needed to perform material comparisons with database entries
    xstr.SetSpecies(xstr.species);

    //DX20190329 - added species check - START
    // check if fake names for same species comparison
    if(LDEBUG){cerr << function_name << " input structure species=" << aurostd::joinWDelimiter(vspecies,",") << endl;}
    if(vspecies[0]=="A" && !structure_comparison){
      message << "Atomic species are missing for the input structure. Cannot compare to database materials without species.";     
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_ERROR_);
      return final_prototypes; //empty //DX20200225
    }
    //DX20190329 - added species check - END

    // ---------------------------------------------------------------------------
    // calculate symmetry
    int space_group_number = xstr.SpaceGroup_ITC();

    // ---------------------------------------------------------------------------
    // get stoichiometries
    vector<uint> stoichiometry = xstr.GetReducedComposition(!same_species);
    uint stoichiometry_sum = aurostd::sum(stoichiometry);
    vector<double> normalized_stoichiometry;
    for(uint i=0;i<stoichiometry.size();i++){normalized_stoichiometry.push_back((double)stoichiometry[i]/(double)stoichiometry_sum);}

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get entries with compatible space groups, i.e., same or enantiomorph
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY")){
      string space_group_summons = aflowlib::getSpaceGroupAFLUXSummons(space_group_number,relaxation_step); //DX20200929 - consolidated formatting to a function
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
      species_summons = "species(" + aurostd::joinWDelimiter(vspecies,",") + ")"; //DX20200212
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
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // call AFLUX 
    string response = aflowlib::AFLUXCall(Summons);

    message << "Number of entries returned: " << aurostd::string2tokens(response,tokens,"\n");
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    if(LDEBUG) {cerr << function_name << " AFLUX response:" << endl << response << endl;}

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
    // get schema from xoptions, i.e., metadata (for the units)
    string schema_unit = "";
    for(uint p=0;p<property_list.size();p++){
      schema_unit = "SCHEMA::UNIT:" + aurostd::toupper(property_list[p]);
      property_units.push_back(XHOST.vschema.getattachedscheme(schema_unit));
    }

    //DX20191230 [OBSOLETE - FROM AFLUX]     // ---------------------------------------------------------------------------
    //DX20191230 [OBSOLETE - FROM AFLUX]     // get AFLUX schema, i.e., metadata (for the units)
    //DX20191230 [OBSOLETE - FROM AFLUX]     if(schema.size()>0){
    //DX20191230 [OBSOLETE - FROM AFLUX]       schema.push_back(aflux_format);
    //DX20191230 [OBSOLETE - FROM AFLUX]       schema.push_back(paging);
    //DX20191230 [OBSOLETE - FROM AFLUX] 
    //DX20191230 [OBSOLETE - FROM AFLUX]       // call AFLUX to get schema
    //DX20191230 [OBSOLETE - FROM AFLUX]       response = aflowlib::AFLUXCall(schema);
    //DX20191230 [OBSOLETE - FROM AFLUX]       vector<vector<std::pair<string,string> > > schema_response = aflowlib::getPropertiesFromAFLUXResponse(response);
    //DX20191230 [OBSOLETE - FROM AFLUX] 
    //DX20191230 [OBSOLETE - FROM AFLUX] 			// extract units
    //DX20191230 [OBSOLETE - FROM AFLUX] 			for(uint i=0;i<schema_response.size();i++){
    //DX20191230 [OBSOLETE - FROM AFLUX] 				bool units_found = false;
    //DX20191230 [OBSOLETE - FROM AFLUX] 				for(uint j=0;j<schema_response[i].size();j++){
    //DX20191230 [OBSOLETE - FROM AFLUX] 					if(schema_response[i][j].first=="units"){
    //DX20191230 [OBSOLETE - FROM AFLUX] 						property_units.push_back(schema_response[i][j].second);
    //DX20191230 [OBSOLETE - FROM AFLUX] 						units_found=true;
    //DX20191230 [OBSOLETE - FROM AFLUX] 						break;
    //DX20191230 [OBSOLETE - FROM AFLUX] 					}
    //DX20191230 [OBSOLETE - FROM AFLUX] 				}
    //DX20191230 [OBSOLETE - FROM AFLUX] 				if(!units_found){
    //DX20191230 [OBSOLETE - FROM AFLUX] 					property_units.push_back("");
    //DX20191230 [OBSOLETE - FROM AFLUX] 				}
    //DX20191230 [OBSOLETE - FROM AFLUX] 			}
    //DX20191230 [OBSOLETE - FROM AFLUX]       if(LDEBUG) {
    //DX20191230 [OBSOLETE - FROM AFLUX]         for(uint i=0;i<property_units.size();i++){ cerr << function_name << ": units for " << property_list[i] << ": " << property_units[i] << endl; }
    //DX20191230 [OBSOLETE - FROM AFLUX]       }
    //DX20191230 [OBSOLETE - FROM AFLUX]     }

    message << "Total number of candidate structures from database: " << auids.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Loading structures ... ";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);


    vector<StructurePrototype> all_structures;

    // ---------------------------------------------------------------------------
    // store input structure 
    StructurePrototype input_structure;
    input_structure.structure_representative = xstr;
    input_structure.natoms = xstr.atoms.size(); //DX20200421
    input_structure.ntypes = xstr.num_each_type.size(); //DX20200421
    input_structure.structure_representative.ReScale(1.0); //DX20191105
    input_structure.structure_representative.BringInCell(); //DX20200707
    input_structure.structure_representative_name = "input geometry";
    input_structure.stoichiometry = xstr.GetReducedComposition(!same_species);
    input_structure.elements = xstr.GetElements(true,true); // true: clean names
    input_structure.structure_representative_compound = pflow::prettyPrintCompound(input_structure.elements,input_structure.stoichiometry,no_vrt,false,txt_ft);
    //DX20191105 [MOVED LATER - SAME AS SYMMETRY] input_structure.LFA_environments= compare::computeLFAEnvironment(input_structure.structure_representative); //DX20190711
    input_structure.structure_representative_generated = true; 
    stringstream ss_input; ss_input << xstr;
    input_structure.structure_representative_source = ss_input.str();
    input_structure.structure_representative_relaxation_step = 0; //DX20200429 - input assumed to be unrelaxed
    input_structure.property_names = property_list;
    input_structure.property_units = property_units;
    all_structures.push_back(input_structure);

    // ---------------------------------------------------------------------------
    // load and store entries from the database 
    for(uint i=0; i<auids.size(); i++){
      // first, get stoichiometry from entry
      //DX20191106 [OBSOLETE - switch to getElements] vector<string> species; vector<double> natoms;
      //DX20191106 [OBSOLETE - switch to getElements] XATOM_SplitAlloySpecies(compounds[i], species, natoms);
      vector<double> vcomposition;
      vector<string> species = aurostd::getElements(compounds[i], vcomposition);
      if(LDEBUG){cerr << function_name << " species=" << aurostd::joinWDelimiter(species,",") << endl;}
      vector<uint> tmp_stoich;
      //DX20191106 [OBSOLETE - switch to getElements] for(uint j=0;j<natoms.size();j++)
      for(uint j=0;j<vcomposition.size();j++) //DX20191106
      { //CO20200106 - patching for auto-indenting
        if(aurostd::isinteger(vcomposition[j])){
          tmp_stoich.push_back((uint)aurostd::nint(vcomposition[j]));
        }
        else {
          message << "Expected natoms in " << auids[i] << " to be an integer.";     
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
      }

      //DX20191125 [OBSOLETE] vector<uint> tmp_reduced_stoich = compare::gcdStoich(tmp_stoich);
      vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
      //DX20190402 - need to sort if ignoring species - START
      if(!same_species){
        //DX20191125 [OBSOLETE - REDUNDANT] for(uint i=0; i<tmp_reduced_stoich.size(); i++){
        std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
        //DX20191125 [OBSOLETE - REDUNDANT] }
      }
      //DX20190402 - need to sort if ignoring species - END
      // second, check if stoichiometries are compatible
      // note: do not include in AFLUX matchbook, we would need to specify a range of compatible stoichs (could be expensive)
      // instead: filter on stoichiometry after recieving AFLUX response
      if(compare::sameStoichiometry(stoichiometry,tmp_reduced_stoich)){
        aflowlib::_aflowlib_entry entry; entry.auid=auids[i]; entry.aurl=aurls[i]; 
        vector<string> structure_files;
        if(!pflow::loadXstructures(entry,structure_files,FileMESSAGE,oss,load_most_relaxed_structure_only)){
          pflow::logger(_AFLOW_FILE_NAME_, function_name, "Could not load structure (auid="+entry.auid+") ... skipping...", FileMESSAGE, logstream, _LOGGER_WARNING_);
          continue;
        }
        bool found_structure = false;
        uint structure_index = 0;
        if(load_most_relaxed_structure_only && entry.vstr.size()==1){
          found_structure = true;
          structure_index = 0;
        }
        else if(!load_most_relaxed_structure_only){
          if(entry.vstr.size()==3 && structure_files.size()==3){
            if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_ORIGINAL_ && 
                (structure_files[0] == "POSCAR.orig" || 
                 structure_files[0] == "POSCAR.relax1")){ 
              structure_index = 0; 
              found_structure = true; 
              if(LDEBUG){cerr << function_name << " loaded original structure: " << structure_files[0] << endl;}
            }
            else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_RELAX1_ && 
                (structure_files[1] == "POSCAR.relax2" || 
                 structure_files[1] == "CONTCAR.relax1")){ 
              structure_index = 1; 
              found_structure = true; 
              if(LDEBUG){cerr << function_name << " loaded relax1 structure: " << structure_files[1] << endl;}
            }
          }
        }
        if(found_structure){
          // store entry from database
          StructurePrototype str_proto_tmp;
          deque<string> deque_species; for(uint j=0;j<species.size();j++){deque_species.push_back(species[j]);}
          entry.vstr[structure_index].SetSpecies(deque_species);
          str_proto_tmp.structure_representative = entry.vstr[structure_index];
          str_proto_tmp.natoms = entry.vstr[structure_index].atoms.size(); //DX20200421
          str_proto_tmp.ntypes = entry.vstr[structure_index].num_each_type.size(); //DX20200421
          str_proto_tmp.structure_representative.ReScale(1.0); //DX20191105
          str_proto_tmp.structure_representative.BringInCell(); //DX20200707
          str_proto_tmp.structure_representative_name=entry.getPathAURL(FileMESSAGE,oss,false); //DX20190321 - changed to false, i.e., do not load from common
          str_proto_tmp.structure_representative.directory=str_proto_tmp.structure_representative_name; //DX20190718 - update xstructure.directoryr
          str_proto_tmp.structure_representative_generated=true;
          str_proto_tmp.structure_representative_source="aurl";
          str_proto_tmp.structure_representative_relaxation_step=relaxation_step; //DX20200429
          str_proto_tmp.stoichiometry=tmp_reduced_stoich;
          str_proto_tmp.elements=species; //DX20200903 - needs to be before prettyPrintCompound()
          str_proto_tmp.structure_representative_compound = pflow::prettyPrintCompound(str_proto_tmp.elements,str_proto_tmp.stoichiometry,no_vrt,false,txt_ft);
          //DX20191105 [MOVED LATER - SAME AS SYMMETRY] str_proto_tmp.LFA_environments= compare::computeLFAEnvironment(str_proto_tmp.structure_representative); //DX20190711
          str_proto_tmp.natoms = entry.vstr[structure_index].atoms.size(); //DX20191031
          str_proto_tmp.ntypes = entry.vstr[structure_index].num_each_type.size(); //DX20191031
          // store any properties 
          for(uint l=0;l<properties_response[i].size();l++){
            bool property_requested = false;
            for(uint m=0;m<property_list.size();m++){
              if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
            }
            if(property_requested){
              str_proto_tmp.properties_structure_representative.push_back(properties_response[i][l].second);
            }
          }
          if(LDEBUG) {
            cerr << XPID << "compare::compareStructureDirectory() Found structure: " << str_proto_tmp.structure_representative_name << endl;
          }
          all_structures.push_back(str_proto_tmp);
        }
        else {
          message << "More structures loaded than anticipated for auid=" << auids[i] << " (# structures=" << entry.vstr.size() << ").";     
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
      }
    }
    message << "Total number of candidate structures loaded: " << all_structures.size(); //DX20190403
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_); //DX20190403

    // ---------------------------------------------------------------------------
    // calculate symmetry of database entries (need Wyckoff positions, but database is not sufficiently populated)
    // in the meantime, we calculate
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY")){
      message << "Calculating the symmetry of the structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

      compare::calculateSymmetries(all_structures,num_proc); 

      message << "Symmetries calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
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
    // calculate LFA environments of database entries 
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS")){
      message << "Calculating the environments of the structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

      compare::calculateLFAEnvironments(all_structures,num_proc); 

      message << "Environments calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }  

    // ---------------------------------------------------------------------------
    // group into objects based on stoichiometry and symmetry (Pearson and space group)
    message << "Grouping sets of comparisons.";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> comparison_schemes = compare::groupStructurePrototypes(all_structures, 
        same_species, 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320 - added environment angles
        false); //DX20200103 - condensed booleans to xoptions 
    //cerr << "number of schemes: " << comparison_schemes.size() << endl;
    //cerr << "comparison_schemes: " << comparison_schemes.size() << endl;
    //cerr << "property names size: " << comparison_schemes[0].property_names.size() << endl;

    // ---------------------------------------------------------------------------
    // only compare entries to the input representation, the rest are extraneous comparisons
    vector<StructurePrototype> input_structure_comparison_scheme_only; input_structure_comparison_scheme_only.push_back(comparison_schemes[0]);
    message << "Number of structures to compare to input structure: " << input_structure_comparison_scheme_only[0].structures_duplicate_names.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // compare structures
    message << "Running comparisons...";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    final_prototypes = compare::runComparisonScheme(input_structure_comparison_scheme_only, same_species, num_proc, comparison_options, oss, FileMESSAGE, quiet, logstream); //DX20200103 - condensed booleans to xoptions //DX20200225 - declare variable at top of function

    return final_prototypes;
  }
}
*/

// ***************************************************************************
// XtalFinderCalculator::compare2database - compares database 
// ***************************************************************************
  // load input structure
  vector<StructurePrototype> XtalFinderCalculator::compare2database(
      const xstructure& xstrIN, const aurostd::xoption& vpflow){
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    vector<StructurePrototype> final_prototypes; //DX20200225
    xstructure xstr = xstrIN; //copy //DX20200225

    string function_name = XPID + "XtalFinderCalculator::compare2database():";
    string directory = "";
    stringstream message;
    bool quiet = false;

    string usage="aflow --compare2database < POSCAR";
    string options="";

    vector<string> tokens,sub_tokens;
    vector<string> matchbook; //aflux - filter/get properties
    vector<string> schema; //get metadata of properties (e.g., units)
    vector<string> property_units;

    bool same_species = true;

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103
    
    // ---------------------------------------------------------------------------
    // get options from vpflow/command line
    getOptions(vpflow,comparison_options); //DX20200103

    // ---------------------------------------------------------------------------
    // single round of comparisons 
    comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);

    /*
    // ---------------------------------------------------------------------------
    // FLAG: misfit threshold //DX20201119
    //double misfit_match_threshold = DEFAULT_XTALFINDER_MISFIT_MATCH;
    //double misfit_family_threshold = DEFAULT_XTALFINDER_MISFIT_FAMILY;
    if(vpflow.flag("COMPARE2DATABASE::MISFIT_MATCH")) {
      misfit_match = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE2DATABASE::MISFIT_MATCH"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_MATCH",vpflow.getattachedscheme("COMPARE2DATABASE::MISFIT_MATCH"));
    }
    if(vpflow.flag("COMPARE2DATABASE::MISFIT_FAMILY")) {
      misfit_family = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE2DATABASE::MISFIT_FAMILY"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_FAMILY",vpflow.getattachedscheme("COMPARE2DATABASE::MISFIT_FAMILY"));
    }
    // match threshold must be less than family threshold
    if(misfit_match>misfit_family){
      message << "Matching misfit threshold must be less than the same family threshold:"
        << " misfit_match_threshold: " << misfit_match
        << " misfit_family_threshold: " << misfit_family;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ILLEGAL_);
    }
    message << "Misfit theshold for matched structures: " << misfit_match << " (default: " << DEFAULT_XTALFINDER_MISFIT_MATCH << ")"; 
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Misfit theshold for structures in the same family: " << misfit_family << " (default: " << DEFAULT_XTALFINDER_MISFIT_FAMILY << ")";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // single round of comparisons 
    comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);

    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE2DATABASE::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE2DATABASE::NP"));
    }

    // ---------------------------------------------------------------------------
    // FLAG: optimize match
    if(vpflow.flag("COMPARE2DATABASE::OPTIMIZE_MATCH")) {
      comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",TRUE); //DX20200225 - was missing
      message << "OPTIONS: Finding optimal match; exploring all possible lattices and origins to find the best match (note: this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    if(vpflow.flag("COMPARE2DATABASE::NO_SCALE_VOLUME")) {
      comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
      message << "OPTIONS: Suppressing volume scaling; useful for distinguishing structures at different pressures."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    if(vpflow.flag("COMPARE2DATABASE::IGNORE_WYCKOFF")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    if(vpflow.flag("COMPARE2DATABASE::IGNORE_SYMMETRY")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    if(vpflow.flag("COMPARE2DATABASE::IGNORE_ENVIRONMENT_ANALYSIS")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",TRUE);
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: primitivize structures //DX20201005
    if(vpflow.flag("COMPARE2DATABASE::PRIMITIVIZE")) {
      comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE",TRUE);
      message << "OPTIONS: Converting all structures to a primitive representation.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: Minkowski reduction //DX20201005
    if(vpflow.flag("COMPARE2DATABASE::MINKOWSKI")) {
      comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI",TRUE);
      message << "OPTIONS: Performing Minkowski lattice reduction on all structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: Niggli reduction //DX20201005
    if(vpflow.flag("COMPARE2DATABASE::NIGGLI")) {
      comparison_options.flag("COMPARISON_OPTIONS::NIGGLI",TRUE);
      message << "OPTIONS: Performing Niggli lattice reduction on all structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }
    */

    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool structure_comparison=false;
    if(vpflow.flag("COMPARE2DATABASE::STRUCTURE")) {
      structure_comparison=true;
      same_species = false;
      message << "OPTIONS: Structure-type comparison, i.e., ignore atomic species."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
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
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify the geometry file to grab (orig, relax1, relax2, static, bands, POSCAR, CONTCAR)
    string geometry_file = "";
    if(vpflow.flag("COMPARE2DATABASE::GEOMETRY_FILE")) {
      geometry_file = vpflow.getattachedscheme("COMPARE2DATABASE::GEOMETRY_FILE");
      message << "OPTIONS: Structure type (POSCAR.orig, POSCAR.relax1, POSCAR.relax2, CONTCAR.relax1, ...): " << geometry_file << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }
    
    uint relaxation_step = _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_;
    bool load_most_relaxed_structure_only = true; 
    if(vpflow.flag("COMPARE2DATABASE::RELAXATION_STEP")) {
      if(aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP")), "orig") || 
          vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP") == "0"){
        relaxation_step = _COMPARE_DATABASE_GEOMETRY_ORIGINAL_;
        load_most_relaxed_structure_only = false;
      }
      if(aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP")), "relax1") || 
          aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP")), "middle_relax") || 
          vpflow.getattachedscheme("COMPARE2DATABASE::RELAXATION_STEP") == "1"){
        relaxation_step = _COMPARE_DATABASE_GEOMETRY_RELAX1_;
        load_most_relaxed_structure_only = false;
      }
      message << "OPTIONS: Relaxation step (0=original, 1=relax1, 2=most_relaxed): " << relaxation_step << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // should not grab properties and compare structures other than the most relaxed structure 
    if(property_list.size() && relaxation_step != _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_){
      string relaxation_name = "";
      if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_ORIGINAL_){ relaxation_name = "original"; }
      else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_RELAX1_){ relaxation_name = "relax1"; }
      message << "The " << relaxation_name << " structures will be extracted; the properties will not correspond to these structures. Proceed with caution."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: catalog (icsd, lib1, lib2, lib3, ...)
    string catalog = "";
    string catalog_summons = "";
    if(vpflow.flag("COMPARE2DATABASE::CATALOG")) {
      catalog = aurostd::tolower(vpflow.getattachedscheme("COMPARE2DATABASE::CATALOG")); //DX20190329 -- added tolower
      if(catalog != "all"){ //DX20190329 - added if-statement since AFLUX doesn't use "all"
        catalog_summons = "catalog(\'" + catalog + "\')";
        matchbook.push_back(catalog_summons);
      } //DX20190329 - added if-statement since AFLUX doesn't use "all"
      message << "OPTIONS: Catalog/library (icsd, lib1, lib2, lib3, ...): " << catalog << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }

    //DX20190508 - added keep unmatched - START
    // ---------------------------------------------------------------------------
    // FLAG: do not remove unmatched structures from the StructurePrototype Object
    // keeps results of each comparison
    if(vpflow.flag("COMPARE2DATABASE::KEEP_UNMATCHED")) {
      comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE);
    }
    //DX20190508 - added keep unmatched - END

    // ---------------------------------------------------------------------------
    // fix species (remove pseudopotentials, etc.) 
    string species_str = aurostd::joinWDelimiter(xstr.species, ""); //DX20200212 
    vector<string> vspecies = aurostd::getElements(species_str); //DX20200212
    xstr.species = aurostd::vector2deque(vspecies); //DX20200212 - needed to perform material comparisons with database entries
    xstr.SetSpecies(xstr.species);

    //DX20190329 - added species check - START
    // check if fake names for same species comparison
    if(LDEBUG){cerr << function_name << " input structure species=" << aurostd::joinWDelimiter(vspecies,",") << endl;}
    if(vspecies[0]=="A" && !structure_comparison){
      message << "Atomic species are missing for the input structure. Cannot compare to database materials without species.";     
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
      return final_prototypes; //empty //DX20200225
    }
    //DX20190329 - added species check - END
    
    // ---------------------------------------------------------------------------
    // add structure to container
    stringstream ss_input; ss_input << xstr;
    addStructureToContainer(xstr, "input geometry", ss_input.str(), 0, false);
    
    // ---------------------------------------------------------------------------
    // symmetry
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && (xstr.space_group_ITC<1 || xstr.space_group_ITC>230)){ //DX20190829 - don't recalculate symmetry if already calculated //DX20191220 - put range instead of ==0
      message << "Calculating the symmetry of the input structure.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      uint one_proc=1;
      calculateSymmetries(one_proc); 
      message << "Symmetry calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    }
    else if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && xstr.space_group_ITC>=1 && xstr.space_group_ITC<=230){ //DX20191220 - put range instead of !=0
      for(uint i=0;i<structure_containers.size();i++){
        structure_containers[i].Pearson = "xX";
        structure_containers[i].space_group = 0;
        vector<GroupedWyckoffPosition> vGWyckoffPos_tmp;
        structure_containers[i].grouped_Wyckoff_positions = vGWyckoffPos_tmp;
      }
    }

    if(LDEBUG) {
      cerr << function_name << " Wyckoff positions of input structure:" << endl;
      for(uint i=0;i<structure_containers[0].grouped_Wyckoff_positions.size();i++){
        cerr << structure_containers[0].grouped_Wyckoff_positions[i] << endl;
      }
    }

    //// ---------------------------------------------------------------------------
    //// calculate symmetry
    //int space_group_number = xstr.SpaceGroup_ITC();

    // ---------------------------------------------------------------------------
    // get stoichiometries
    vector<uint> stoichiometry = structure_containers[0].stoichiometry; // xstr.GetReducedComposition(!same_species);
    uint stoichiometry_sum = aurostd::sum(stoichiometry);
    vector<double> normalized_stoichiometry;
    for(uint i=0;i<stoichiometry.size();i++){normalized_stoichiometry.push_back((double)stoichiometry[i]/(double)stoichiometry_sum);}

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get entries with compatible space groups, i.e., same or enantiomorph
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY")){
      //string space_group_summons = aflowlib::getSpaceGroupAFLUXSummons(space_group_number,relaxation_step); //DX20200929 - consolidated formatting to a function
      string space_group_summons = aflowlib::getSpaceGroupAFLUXSummons(structure_containers[0].space_group,relaxation_step); //DX20200929 - consolidated formatting to a function
      matchbook.push_back(space_group_summons);
    }
    //vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
    //compare::groupWyckoffPositions(xstr, grouped_Wyckoff_positions);

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get aurl for entry
    string aurl = "aurl";
    matchbook.push_back(aurl);

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get species and number of species
    string species_summons = "";
    if(!structure_comparison){
      species_summons = "species(" + aurostd::joinWDelimiter(vspecies,",") + ")"; //DX20200212
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
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // call AFLUX 
    string response = aflowlib::AFLUXCall(Summons);

    message << "Number of entries returned: " << aurostd::string2tokens(response,tokens,"\n");
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    if(LDEBUG) {cerr << function_name << " AFLUX response:" << endl << response << endl;}

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
    // get schema from xoptions, i.e., metadata (for the units)
    string schema_unit = "";
    for(uint p=0;p<property_list.size();p++){
      schema_unit = "SCHEMA::UNIT:" + aurostd::toupper(property_list[p]);
      property_units.push_back(XHOST.vschema.getattachedscheme(schema_unit));
    }

    //DX20191230 [OBSOLETE - FROM AFLUX]     // ---------------------------------------------------------------------------
    //DX20191230 [OBSOLETE - FROM AFLUX]     // get AFLUX schema, i.e., metadata (for the units)
    //DX20191230 [OBSOLETE - FROM AFLUX]     if(schema.size()>0){
    //DX20191230 [OBSOLETE - FROM AFLUX]       schema.push_back(aflux_format);
    //DX20191230 [OBSOLETE - FROM AFLUX]       schema.push_back(paging);
    //DX20191230 [OBSOLETE - FROM AFLUX] 
    //DX20191230 [OBSOLETE - FROM AFLUX]       // call AFLUX to get schema
    //DX20191230 [OBSOLETE - FROM AFLUX]       response = aflowlib::AFLUXCall(schema);
    //DX20191230 [OBSOLETE - FROM AFLUX]       vector<vector<std::pair<string,string> > > schema_response = aflowlib::getPropertiesFromAFLUXResponse(response);
    //DX20191230 [OBSOLETE - FROM AFLUX] 
    //DX20191230 [OBSOLETE - FROM AFLUX] 			// extract units
    //DX20191230 [OBSOLETE - FROM AFLUX] 			for(uint i=0;i<schema_response.size();i++){
    //DX20191230 [OBSOLETE - FROM AFLUX] 				bool units_found = false;
    //DX20191230 [OBSOLETE - FROM AFLUX] 				for(uint j=0;j<schema_response[i].size();j++){
    //DX20191230 [OBSOLETE - FROM AFLUX] 					if(schema_response[i][j].first=="units"){
    //DX20191230 [OBSOLETE - FROM AFLUX] 						property_units.push_back(schema_response[i][j].second);
    //DX20191230 [OBSOLETE - FROM AFLUX] 						units_found=true;
    //DX20191230 [OBSOLETE - FROM AFLUX] 						break;
    //DX20191230 [OBSOLETE - FROM AFLUX] 					}
    //DX20191230 [OBSOLETE - FROM AFLUX] 				}
    //DX20191230 [OBSOLETE - FROM AFLUX] 				if(!units_found){
    //DX20191230 [OBSOLETE - FROM AFLUX] 					property_units.push_back("");
    //DX20191230 [OBSOLETE - FROM AFLUX] 				}
    //DX20191230 [OBSOLETE - FROM AFLUX] 			}
    //DX20191230 [OBSOLETE - FROM AFLUX]       if(LDEBUG) {
    //DX20191230 [OBSOLETE - FROM AFLUX]         for(uint i=0;i<property_units.size();i++){ cerr << function_name << ": units for " << property_list[i] << ": " << property_units[i] << endl; }
    //DX20191230 [OBSOLETE - FROM AFLUX]       }
    //DX20191230 [OBSOLETE - FROM AFLUX]     }

    message << "Total number of candidate structures from database: " << auids.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    message << "Loading structures ... ";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);


    vector<StructurePrototype> all_structures;

    /*
    // ---------------------------------------------------------------------------
    // store input structure 
    StructurePrototype input_structure;
    input_structure.structure_representative = xstr;
    input_structure.natoms = xstr.atoms.size(); //DX20200421
    input_structure.ntypes = xstr.num_each_type.size(); //DX20200421
    input_structure.structure_representative.ReScale(1.0); //DX20191105
    input_structure.structure_representative.BringInCell(); //DX20200707
    input_structure.structure_representative_name = "input geometry";
    input_structure.stoichiometry = xstr.GetReducedComposition(!same_species);
    input_structure.elements = xstr.GetElements(true,true); // true: clean names
    input_structure.structure_representative_compound = pflow::prettyPrintCompound(input_structure.elements,input_structure.stoichiometry,no_vrt,false,txt_ft);
    //DX20191105 [MOVED LATER - SAME AS SYMMETRY] input_structure.LFA_environments= compare::computeLFAEnvironment(input_structure.structure_representative); //DX20190711
    input_structure.structure_representative_generated = true; 
    stringstream ss_input; ss_input << xstr;
    input_structure.structure_representative_source = ss_input.str();
    input_structure.structure_representative_relaxation_step = 0; //DX20200429 - input assumed to be unrelaxed
    input_structure.property_names = property_list;
    input_structure.property_units = property_units;
    all_structures.push_back(input_structure);
*/


    // ---------------------------------------------------------------------------
    // load and store entries from the database 
    for(uint i=0; i<auids.size(); i++){
      // first, get stoichiometry from entry
      //DX20191106 [OBSOLETE - switch to getElements] vector<string> species; vector<double> natoms;
      //DX20191106 [OBSOLETE - switch to getElements] XATOM_SplitAlloySpecies(compounds[i], species, natoms);
      vector<double> vcomposition;
      vector<string> species = aurostd::getElements(compounds[i], vcomposition);
      if(LDEBUG){cerr << function_name << " species=" << aurostd::joinWDelimiter(species,",") << endl;}
      vector<uint> tmp_stoich;
      //DX20191106 [OBSOLETE - switch to getElements] for(uint j=0;j<natoms.size();j++)
      for(uint j=0;j<vcomposition.size();j++) //DX20191106
      { //CO20200106 - patching for auto-indenting
        if(aurostd::isinteger(vcomposition[j])){
          tmp_stoich.push_back((uint)aurostd::nint(vcomposition[j]));
        }
        else {
          message << "Expected natoms in " << auids[i] << " to be an integer.";     
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
      }

      //DX20191125 [OBSOLETE] vector<uint> tmp_reduced_stoich = compare::gcdStoich(tmp_stoich);
      vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
      //DX20190402 - need to sort if ignoring species - START
      if(!same_species){
        //DX20191125 [OBSOLETE - REDUNDANT] for(uint i=0; i<tmp_reduced_stoich.size(); i++){
        std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
        //DX20191125 [OBSOLETE - REDUNDANT] }
      }
      //DX20190402 - need to sort if ignoring species - END
      // second, check if stoichiometries are compatible
      // note: do not include in AFLUX matchbook, we would need to specify a range of compatible stoichs (could be expensive)
      // instead: filter on stoichiometry after recieving AFLUX response
      if(compare::sameStoichiometry(stoichiometry,tmp_reduced_stoich)){
        aflowlib::_aflowlib_entry entry; entry.auid=auids[i]; entry.aurl=aurls[i]; 
        vector<string> structure_files;
        if(!pflow::loadXstructures(entry,structure_files,*p_FileMESSAGE,*p_oss,load_most_relaxed_structure_only)){
          pflow::logger(_AFLOW_FILE_NAME_, function_name, "Could not load structure (auid="+entry.auid+") ... skipping...", *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
          continue;
        }
        bool found_structure = false;
        uint structure_index = 0;
        if(load_most_relaxed_structure_only && entry.vstr.size()==1){
          found_structure = true;
          structure_index = 0;
        }
        else if(!load_most_relaxed_structure_only){
          if(entry.vstr.size()==3 && structure_files.size()==3){
            if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_ORIGINAL_ && 
                (structure_files[0] == "POSCAR.orig" || 
                 structure_files[0] == "POSCAR.relax1")){ 
              structure_index = 0; 
              found_structure = true; 
              if(LDEBUG){cerr << function_name << " loaded original structure: " << structure_files[0] << endl;}
            }
            else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_RELAX1_ && 
                (structure_files[1] == "POSCAR.relax2" || 
                 structure_files[1] == "CONTCAR.relax1")){ 
              structure_index = 1; 
              found_structure = true; 
              if(LDEBUG){cerr << function_name << " loaded relax1 structure: " << structure_files[1] << endl;}
            }
          }
        }
        if(found_structure){
          // store entry from database
          deque<string> deque_species; for(uint j=0;j<species.size();j++){deque_species.push_back(species[j]);}
          entry.vstr[structure_index].SetSpecies(deque_species);
          string str_path = entry.getPathAURL(*p_FileMESSAGE,*p_oss,false);
          addStructureToContainer(entry.vstr[structure_index], str_path, "aurl", relaxation_step, same_species);
          /*
          StructurePrototype str_proto_tmp;
          deque<string> deque_species; for(uint j=0;j<species.size();j++){deque_species.push_back(species[j]);}
          entry.vstr[structure_index].SetSpecies(deque_species);
          str_proto_tmp.structure_representative = entry.vstr[structure_index];
          str_proto_tmp.natoms = entry.vstr[structure_index].atoms.size(); //DX20200421
          str_proto_tmp.ntypes = entry.vstr[structure_index].num_each_type.size(); //DX20200421
          str_proto_tmp.structure_representative.ReScale(1.0); //DX20191105
          str_proto_tmp.structure_representative.BringInCell(); //DX20200707
          str_proto_tmp.structure_representative_name=entry.getPathAURL(*p_FileMESSAGE,oss,false); //DX20190321 - changed to false, i.e., do not load from common
          str_proto_tmp.structure_representative.directory=str_proto_tmp.structure_representative_name; //DX20190718 - update xstructure.directoryr
          str_proto_tmp.structure_representative_generated=true;
          str_proto_tmp.structure_representative_source="aurl";
          str_proto_tmp.structure_representative_relaxation_step=relaxation_step; //DX20200429
          str_proto_tmp.stoichiometry=tmp_reduced_stoich;
          str_proto_tmp.elements=species; //DX20200903 - needs to be before prettyPrintCompound()
          str_proto_tmp.structure_representative_compound = pflow::prettyPrintCompound(str_proto_tmp.elements,str_proto_tmp.stoichiometry,no_vrt,false,txt_ft);
          //DX20191105 [MOVED LATER - SAME AS SYMMETRY] str_proto_tmp.LFA_environments= compare::computeLFAEnvironment(str_proto_tmp.structure_representative); //DX20190711
          str_proto_tmp.natoms = entry.vstr[structure_index].atoms.size(); //DX20191031
          str_proto_tmp.ntypes = entry.vstr[structure_index].num_each_type.size(); //DX20191031
          */
          // store any properties 
          for(uint l=0;l<properties_response[i].size();l++){
            bool property_requested = false;
            for(uint m=0;m<property_list.size();m++){
              if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
            }
            if(property_requested){
              //str_proto_tmp.properties_structure_representative.push_back(properties_response[i][l].second);
              structure_containers.back().properties.push_back(properties_response[i][l].second);
              structure_containers.back().properties_names = property_list;
              structure_containers.back().properties_units = property_units;
              //structure_containers.back().properties_types = property_types;
            }
          }
          //if(LDEBUG) {
          //  cerr << XPID << "compare::compareStructureDirectory() Found structure: " << str_proto_tmp.structure_representative_name << endl;
          //}
          //all_structures.push_back(str_proto_tmp);
        }
        else {
          message << "More structures loaded than anticipated for auid=" << auids[i] << " (# structures=" << entry.vstr.size() << ").";     
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
      }
    }
    message << "Total number of candidate structures loaded: " << all_structures.size(); //DX20190403
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_); //DX20190403

    /*
    // ---------------------------------------------------------------------------
    // calculate symmetry of database entries (need Wyckoff positions, but database is not sufficiently populated)
    // in the meantime, we calculate
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY")){
      message << "Calculating the symmetry of the structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      calculateSymmetries(one_proc); 
      message << "Symmetries calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    }  
    else {
      for(uint i=0;i<structure_containers.size();i++){
        structure_containers[i].Pearson = "xX";
        structure_containers[i].space_group = 0;
        vector<GroupedWyckoffPosition> vGWyckoffPos_tmp;
        structure_containers[i].grouped_Wyckoff_positions = vGWyckoffPos_tmp;
      }
    }

    // ---------------------------------------------------------------------------
    // calculate LFA environments of database entries 
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS")){
      message << "Calculating the environments of the structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      calculateLFAEnvironments(num_proc); 
      message << "Environments calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    } 
    

    // ---------------------------------------------------------------------------
    // group into objects based on stoichiometry and symmetry (Pearson and space group)
    message << "Grouping sets of comparisons.";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    vector<StructurePrototype> comparison_schemes = groupStructurePrototypes(all_structures, 
        same_species, 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320 - added environment angles
        false); //DX20200103 - condensed booleans to xoptions 
    //cerr << "number of schemes: " << comparison_schemes.size() << endl;
    //cerr << "comparison_schemes: " << comparison_schemes.size() << endl;
    //cerr << "property names size: " << comparison_schemes[0].property_names.size() << endl;

    // ---------------------------------------------------------------------------
    // only compare entries to the input representation, the rest are extraneous comparisons
    vector<StructurePrototype> input_structure_comparison_scheme_only; input_structure_comparison_scheme_only.push_back(comparison_schemes[0]);
    message << "Number of structures to compare to input structure: " << input_structure_comparison_scheme_only[0].structures_duplicate_names.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // compare structures
    message << "Running comparisons...";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    final_prototypes = compare::runComparisonScheme(input_structure_comparison_scheme_only, same_species, num_proc, comparison_options, oss, *p_FileMESSAGE, quiet, *p_oss); //DX20200103 - condensed booleans to xoptions //DX20200225 - declare variable at top of function
    return final_prototypes;
*/
    return compareMultipleStructures(num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions 
  }

// ***************************************************************************
// compare::compare2database - compares database 
// ***************************************************************************
namespace compare {
  // load input structure
  string printCompare2Database(istream& input, const aurostd::xoption& vpflow, ostream& logstream){xstructure xstr(input,IOAFLOW_AUTO); ofstream FileMESSAGE; return printCompare2Database(xstr,vpflow,FileMESSAGE,logstream);}  //DX20200225
  string printCompare2Database(istream& input, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){xstructure xstr(input,IOAFLOW_AUTO);return printCompare2Database(xstr,vpflow,FileMESSAGE,logstream);}  //CO20200225
  string printCompare2Database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = "compare::printCompare2Database():";
    stringstream message;
    ostringstream oss;

    if(LDEBUG){cerr << function_name << " BEGIN" << endl;}  //CO20200508

    // ---------------------------------------------------------------------------
    // main compare2database() function
    //vector<StructurePrototype> final_prototypes = compare::compare2database(xstrIN, vpflow, FileMESSAGE, logstream);
    XtalFinderCalculator xtal_finder_database(DEFAULT_XTALFINDER_MISFIT_MATCH,DEFAULT_XTALFINDER_MISFIT_FAMILY,FileMESSAGE,1,logstream);
    vector<StructurePrototype> final_prototypes = xtal_finder_database.compare2database(xstrIN, vpflow);

    // ---------------------------------------------------------------------------
    // return if there are no similar structures
    if(final_prototypes.size()==0){
      return oss.str();
    }

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
    if(vpflow.flag("COMPARE2DATABASE::SCREEN_ONLY")) {
      screen_only=true;
    }

    // ---------------------------------------------------------------------------
    // FLAG: type of comparison (material-type or structure-type)
    bool same_species = true;
    if(vpflow.flag("COMPARE2DATABASE::STRUCTURE")) {
      same_species = false;
    }

    // ---------------------------------------------------------------------------
    // print results 
    stringstream ss_out;
    xtal_finder_database.printResults(ss_out, same_species, final_prototypes, "text");
    stringstream ss_json;
    xtal_finder_database.printResults(ss_json, same_species, final_prototypes, "json");

    // DEBUG oss << ss_out.str();
    message << "Number of structures in database matching with the input structure: " << final_prototypes[0].structures_duplicate_struct.size() << "." << endl;
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    //DX20190429 - added screen only option - START
    // ---------------------------------------------------------------------------
    // write results to screen and return immediately (do not write file)
    if(screen_only){
      if(format=="json"){ return ss_json.str(); }
      // default is txt
      else { return ss_out.str(); }
    }
    //DX20190429 - added screen only option - END

    // ---------------------------------------------------------------------------
    // write results to files
    if(same_species){  
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
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);

    return oss.str();

  }
}

/*
//DX - COMPARE DATABASE ENTRIES - START
// ***************************************************************************
// compare::compareDatabaseEntries - compares database entries
// ***************************************************************************
namespace compare {
  string compareDatabaseEntries(const aurostd::xoption& vpflow, ostream& logstream){ //DX20191125 - added ofstream overload and added ostream as input
    ofstream FileMESSAGE;
    return compareDatabaseEntries(vpflow, FileMESSAGE, logstream);
  }

  string compareDatabaseEntries(const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){ //DX20191125 - added ofstream and ostream
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    string function_name = XPID + "compare::compareDatabaseEntries():";
    string directory = ".";
    stringstream message;
    stringstream oss;

    string usage="aflow --compare_database_entries < POSCAR";
    string options="";

    vector<string> tokens,sub_tokens;
    vector<string> matchbook; //aflux - filter/get properties
    vector<string> schema; //get metadata of properties (e.g., units)
    vector<string> property_units;

    bool same_species = true;

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103

    // ---------------------------------------------------------------------------
    // FLAG: misfit threshold //DX20201119
    double misfit_match_threshold = DEFAULT_XTALFINDER_MISFIT_MATCH;
    double misfit_family_threshold = DEFAULT_XTALFINDER_MISFIT_FAMILY;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::MISFIT_MATCH")) {
      misfit_match_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::MISFIT_MATCH"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_MATCH",vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::MISFIT_MATCH"));
    }
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::MISFIT_FAMILY")) {
      misfit_family_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::MISFIT_FAMILY"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_FAMILY",vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::MISFIT_FAMILY"));
    }
    // match threshold must be less than family threshold
    if(misfit_match_threshold>misfit_family_threshold){
      message << "Matching misfit threshold must be less than the same family threshold:"
        << " misfit_match_threshold: " << misfit_match_threshold
        << " misfit_family_threshold: " << misfit_family_threshold;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ILLEGAL_);
    }
    message << "Misfit theshold for matched structures: " << misfit_match_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_MATCH << ")"; 
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Misfit theshold for structures in the same family: " << misfit_family_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_FAMILY << ")";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::NP"));
    }

    // ---------------------------------------------------------------------------
    // FLAG: optimize match
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::OPTIMIZE_MATCH")) {
      comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",TRUE);
      message << "OPTIONS: Finding optimal match; exploring all possible lattices and origins to find the best match (note: this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::NO_SCALE_VOLUME")) {
      comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
      message << "OPTIONS: Suppressing volume scaling; useful for distinguishing structures at different pressures."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::IGNORE_WYCKOFF")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::IGNORE_SYMMETRY")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::IGNORE_ENVIRONMENT_ANALYSIS")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",TRUE);
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: primitivize structures //DX20201005
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::PRIMITIVIZE")) {
      comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE",TRUE);
      message << "OPTIONS: Converting all structures to a primitive representation.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: Minkowski reduction //DX20201005
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::MINKOWSKI")) {
      comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI",TRUE);
      message << "OPTIONS: Performing Minkowski lattice reduction on all structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: Niggli reduction //DX20201005
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::NIGGLI")) {
      comparison_options.flag("COMPARISON_OPTIONS::NIGGLI",TRUE);
      message << "OPTIONS: Performing Niggli lattice reduction on all structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
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
        //DX20191106 [OBSOLETE - switch to getElements] XATOM_SplitAlloySpecies(alloy_string, species);
        species = aurostd::getElements(alloy_string); //DX20191106
      }
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify space group
    vector<uint> space_groups;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::SPACE_GROUP")){
      string space_group_input = vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::SPACE_GROUP");
      vector<string> space_group_strings;
      aurostd::string2tokens(space_group_input,space_group_strings,",");

      uint sg_tmp = 0;
      for(uint i=0;i<space_group_strings.size();i++){
        sg_tmp = aurostd::string2utype<uint>(space_group_strings[i]);
        if(sg_tmp<1 || sg_tmp>230){
          message << "Invalid space group requested: " << space_group_strings[i] << ". Please check input.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ERROR_);
        }
        space_groups.push_back(sg_tmp);
      }
      message << "OPTIONS: Requesting the following space groups: " << space_group_input;
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

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
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: do not remove unmatched structures from the StructurePrototype Object
    // keeps results of each comparison
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::KEEP_UNMATCHED")) {
      comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE);
    }

    // ---------------------------------------------------------------------------
    // FLAG: remove duplicate compounds (useful for non-biased statistics)
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::REMOVE_DUPLICATE_COMPOUNDS")) {
      comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",TRUE);
      message << "OPTIONS: Remove duplicate compounds first, useful for non-biased prototype statistics."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: match unique structures to the AFLOW prototypes 
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::MATCH_TO_AFLOW_PROTOS")) {
      comparison_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS",TRUE);
      message << "OPTIONS: Compare unique structures to the AFLOW prototypes."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: get AFLOW ANRL designation for unique structures
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::ADD_AFLOW_PROTOTYPE_DESIGNATION")) {
      comparison_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION",TRUE);
      message << "OPTIONS: Cast unique structures into AFLOW standard designation."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: consider magnetic structure in comparison //DX20191212
    bool magnetic_comparison=false;
    vector<string> magmoms_for_systems;
    if(vpflow.flag("COMPARE::MAGNETIC")){
      magnetic_comparison=true;
      message << "OPTIONS: Including magnetic moment information in comparisons.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    if(magnetic_comparison){} //CO20200508 - keep it busy

    // ---------------------------------------------------------------------------
    // FLAG: do not calculate unique atom decorations 
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS")) {
      comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",FALSE);
      message << "OPTIONS: Do not calculate unique atom decorations."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
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
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify the geometry file to grab (orig, relax1, relax2, static, bands, POSCAR, CONTCAR)
    //DX TODO 
    string geometry_file = "";
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::GEOMETRY_FILE")) {
      geometry_file = vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::GEOMETRY_FILE");
      message << "OPTIONS: Structure type (POSCAR.orig, POSCAR.relax1, POSCAR.relax2, CONTCAR.relax1, ...): " << geometry_file << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify the relaxation step to grab (orig, relax1, relax2, static, bands, POSCAR, CONTCAR)
    uint relaxation_step = _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_;
    bool load_most_relaxed_structure_only = true; 
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")) {
      if(aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")), "orig") || 
          vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP") == "0" ){
        relaxation_step = _COMPARE_DATABASE_GEOMETRY_ORIGINAL_;
        load_most_relaxed_structure_only = false;
      }
      if(aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")), "relax1") || 
          aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")), "middle_relax") || 
          vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP") == "1"){
        relaxation_step = _COMPARE_DATABASE_GEOMETRY_RELAX1_;
        load_most_relaxed_structure_only = false;
      }
      message << "OPTIONS: Relaxation step (0=original, 1=relax1, 2=most_relaxed): " << relaxation_step << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    if(load_most_relaxed_structure_only){} //CO20200508 - keep it busy

    // ---------------------------------------------------------------------------
    // FLAG: catalog (icsd, lib1, lib2, lib3, ...)
    string catalog = "";
    string catalog_summons = "";
    //DX20200331 [OBSOLETE] bool ICSD_comparison = false;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::CATALOG")) {
      catalog = aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::CATALOG")); //DX20190718
      catalog_summons = "catalog(\'" + catalog + "\')";
      matchbook.push_back(catalog_summons);
      message << "OPTIONS: Catalog/library (icsd, lib1, lib2, lib3, ...): " << catalog << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    if(catalog=="" || catalog=="icsd" || catalog=="all"){ //DX20191108 - needs to be outside of loop
      comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON",TRUE);
      //DX20200331 [OBSOLETE] ICSD_comparison=true;
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

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get space group(s) //DX20200929
    string sg_summons = "";
    if(space_groups.size()!=0){
      sg_summons = aflowlib::getSpaceGroupAFLUXSummons(space_groups, relaxation_step);
      if(LDEBUG){ cerr << "Space group summons: " << sg_summons << endl; }
    }

    matchbook.push_back(species_summons);
    matchbook.push_back(nspecies_summons);
    matchbook.push_back(sg_summons);

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
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // call AFLUX 
    string response = aflowlib::AFLUXCall(Summons);

    message << "Number of entries returned: " << aurostd::string2tokens(response,tokens,"\n");
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

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
    // get schema from xoptions, i.e., metadata (for the units) (DX20191230)
    string schema_unit = "";
    for(uint p=0;p<property_list.size();p++){
      schema_unit = "SCHEMA::UNIT:" + aurostd::toupper(property_list[p]);
      property_units.push_back(XHOST.vschema.getattachedscheme(schema_unit));
    }

    //DX20191230 [OBSOLETE - FROM AFLUX]    // ---------------------------------------------------------------------------
    //DX20191230 [OBSOLETE - FROM AFLUX]    // get AFLUX schema, i.e., metadata (for the units)
    //DX20191230 [OBSOLETE - FROM AFLUX]    if(schema.size()>0){
    //DX20191230 [OBSOLETE - FROM AFLUX]      schema.push_back(format);
    //DX20191230 [OBSOLETE - FROM AFLUX]      schema.push_back(paging);
    //DX20191230 [OBSOLETE - FROM AFLUX]
    //DX20191230 [OBSOLETE - FROM AFLUX]      // call AFLUX to get schema
    //DX20191230 [OBSOLETE - FROM AFLUX]      response = aflowlib::AFLUXCall(schema);
    //DX20191230 [OBSOLETE - FROM AFLUX]      vector<vector<std::pair<string,string> > > schema_response = aflowlib::getPropertiesFromAFLUXResponse(response);
    //DX20191230 [OBSOLETE - FROM AFLUX]
    //DX20191230 [OBSOLETE - FROM AFLUX]      // extract units
    //DX20191230 [OBSOLETE - FROM AFLUX]      for(uint i=0;i<schema_response.size();i++){
    //DX20191230 [OBSOLETE - FROM AFLUX]        bool units_found = false;
    //DX20191230 [OBSOLETE - FROM AFLUX]        for(uint j=0;j<schema_response[i].size();j++){
    //DX20191230 [OBSOLETE - FROM AFLUX]          if(schema_response[i][j].first=="units"){
    //DX20191230 [OBSOLETE - FROM AFLUX]            property_units.push_back(schema_response[i][j].second);
    //DX20191230 [OBSOLETE - FROM AFLUX]            units_found=true;
    //DX20191230 [OBSOLETE - FROM AFLUX]            break;
    //DX20191230 [OBSOLETE - FROM AFLUX]          }
    //DX20191230 [OBSOLETE - FROM AFLUX]        }
    //DX20191230 [OBSOLETE - FROM AFLUX]        if(!units_found){
    //DX20191230 [OBSOLETE - FROM AFLUX]          property_units.push_back("");
    //DX20191230 [OBSOLETE - FROM AFLUX]        }
    //DX20191230 [OBSOLETE - FROM AFLUX]      }
    //DX20191230 [OBSOLETE - FROM AFLUX]      if(LDEBUG){
    //DX20191230 [OBSOLETE - FROM AFLUX]        for(uint i=0;i<property_units.size();i++){ cerr << function_name << ": units for " << property_list[i] << ": " << property_units[i] << endl; }
    //DX20191230 [OBSOLETE - FROM AFLUX]      }
    //DX20191230 [OBSOLETE - FROM AFLUX]    }

    message << "Total number of candidate structures from database: " << auids.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Loading structures ..." << auids.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
    message << "Splitting into threads...";     
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    // ---------------------------------------------------------------------------
    // distribute threads via indices
    uint number_of_structures = auids.size();
    //DX20191210 [URL TIMEOUT] uint num_threads = aurostd::min(num_proc,number_of_structures); // cannot have more threads than structures
    uint num_threads = 1; //DX20191210 [URL TIMEOUT]
    //DX20191107 [switching to getThreadDistribution] - vector<uint> start_indices, end_indices;
    //DX20191107 [switching to getThreadDistribution] - compare::splitTaskIntoThreads(number_of_structures,num_threads,start_indices,end_indices);
    vector<vector<int> > thread_distribution = getThreadDistribution(number_of_structures, num_threads); //DX20191107 
    message << "Done. Split into threads.";     
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // initialize vector of objects 
    vector<StructurePrototype> all_structures;
    for(uint i=0; i<auids.size(); i++){
      StructurePrototype str_proto_tmp;
      // first, get stoichiometry from entry
      //DX20191106 [OBSOLETE - switch to getElements] vector<string> species; vector<double> natoms;
      //DX20191106 [OBSOLETE - switch to getElements] XATOM_SplitAlloySpecies(compounds[i], species, natoms);
      vector<double> vcomposition;
      vector<string> species = aurostd::getElements(compounds[i], vcomposition);
      vector<uint> tmp_stoich;
      //DX20191106 [OBSOLETE - switch to getElements] for(uint j=0;j<natoms.size();j++)
      for(uint j=0;j<vcomposition.size();j++) //DX20191106
      { //CO20200106 - patching for auto-indenting
        if(aurostd::isinteger(vcomposition[j])){
          tmp_stoich.push_back((uint)aurostd::nint(vcomposition[j]));
        }
        else{
          message << "Expected natoms in " << auids[i] << " to be an integer.";     
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
      }

      //DX20191125 [OBSOLETE] vector<uint> tmp_reduced_stoich = compare::gcdStoich(tmp_stoich);
      vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
      //DX20190402 - need to sort if ignoring species - START
      if(!same_species){
        //DX20191125 [OBSOLETE - REDUNDANT]for(uint i=0; i<tmp_reduced_stoich.size(); i++){
        std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
        //DX20191125 [OBSOLETE - REDUNDANT] }
      }
      str_proto_tmp.stoichiometry=tmp_reduced_stoich;
      str_proto_tmp.elements=species;
      str_proto_tmp.structure_representative_name=aurls[i];
      str_proto_tmp.structure_representative_source="aurl";
      str_proto_tmp.structure_representative_relaxation_step=relaxation_step; //DX20200429
      all_structures.push_back(str_proto_tmp);
    }
    message << "Finished initializing StructurePrototype object, now spawn threads.";     
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // Run threads (now using pointer to thread)
    vector<std::thread*> threads;
    for(uint n=0; n<num_threads; n++){
      //DX20191107 [OBOSLETE - switch to getThreadDistribution convention] threads.push_back(std::thread(compare::generateStructures,std::ref(all_structures),std::ref(oss),start_indices[n],end_indices[n]));
      threads.push_back(new std::thread(&compare::generateStructures, std::ref(all_structures),std::ref(logstream),thread_distribution[n][0],thread_distribution[n][1])); //DX20191107 //DX20200225 - oss to logstream
    }
    // ---------------------------------------------------------------------------
    // Join threads
    for(uint t=0;t<num_threads;t++){
      threads[t]->join();
      delete threads[t];
    }
    // ---------------------------------------------------------------------------
    // try to regenerate via one thread; sometimes multithreading does not load 
    // all structures (timeout issue)
    if(num_threads!=1){
      compare::generateStructures(all_structures);
    }
    message << "Threads complete. " << all_structures.size() << " structures. Adding properties.";     
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // add properties information 
    for(uint i=0;i<all_structures.size();i++){
      all_structures[i].property_names = property_list; //DX20190326
      all_structures[i].property_units = property_units; //DX20190326
      // store any properties 
      for(uint l=0;l<properties_response[i].size();l++){
        bool property_requested = false;
        for(uint m=0;m<property_list.size();m++){
          if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
        }
        if(property_requested){
          all_structures[i].properties_structure_representative.push_back(properties_response[i][l].second);
        }
      }
    }

    for(uint i=0;i<all_structures.size();i++){
      if(all_structures[i].structure_representative_generated){
        deque<string> deque_species; for(uint j=0;j<all_structures[i].elements.size();j++){deque_species.push_back(all_structures[i].elements[j]);}
        all_structures[i].structure_representative.SetSpecies(deque_species);
        all_structures[i].structure_representative_compound = pflow::prettyPrintCompound(all_structures[i].elements,all_structures[i].stoichiometry,no_vrt,false,txt_ft);
      }
    }

    message << "Properties added, now removing non-generated structures" << endl;     
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    // ---------------------------------------------------------------------------
    // remove non-generated structures
    compare::removeNonGeneratedStructures(all_structures);

    for(uint i=0;i<all_structures.size();i++){
      all_structures[i].structure_representative.ReScale(1.0); //DX20191105
      all_structures[i].structure_representative.BringInCell(); //DX20200707
      //DX20191105 [MOVED LATER - SAME AS SYMMETRY] all_structures[i].LFA_environments= compare::computeLFAEnvironment(all_structures[i].structure_representative); //DX20190711
    }

#else
    // ---------------------------------------------------------------------------
    // load and store entries from the database 
    vector<StructurePrototype> all_structures;
    for(uint i=0; i<auids.size(); i++){
      // first, get stoichiometry from entry
      //DX20191106 [OBSOLETE - switch to getElements] vector<string> species; vector<double> natoms;
      //DX20191106 [OBSOLETE - switch to getElements] XATOM_SplitAlloySpecies(compounds[i], species, natoms);
      vector<double> vcomposition;
      vector<string> species = aurostd::getElements(compounds[i], vcomposition);
      vector<uint> tmp_stoich;
      //DX20191106 [OBSOLETE - switch to getElements] for(uint j=0;j<natoms.size();j++){
      for(uint j=0;j<vcomposition.size();j++){ //DX20191106
        if(aurostd::isinteger(vcomposition[j])){
          tmp_stoich.push_back((uint)aurostd::nint(vcomposition[j]));
        }
        else{
          message << "Expected natoms in " << auids[i] << " to be an integer.";     
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
      }

      //DX20191125 [OBSOLETE} vector<uint> tmp_reduced_stoich = compare::gcdStoich(tmp_stoich);
      vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
      //DX20190402 - need to sort if ignoring species - START
      if(!same_species){
        //DX20191125 [OBSOLETE - REDUNDANT] for(uint i=0; i<tmp_reduced_stoich.size(); i++){
        std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
        //DX20191125 [OBSOLETE - REDUNDANT] }
      }
      //DX20190402 - need to sort if ignoring species - END
      // second, check if stoichiometries are compatible
      // note: do not include in AFLUX matchbook, we would need to specify a range of compatible stoichs (could be expensive)
      // instead: filter on stoichiometry after recieving AFLUX response
      //if(compare::sameStoichiometry(stoichiometry,tmp_reduced_stoich)){
      aflowlib::_aflowlib_entry entry; entry.auid=auids[i]; entry.aurl=aurls[i]; 
      vector<string> structure_files;
      if(!pflow::loadXstructures(entry,structure_files,FileMESSAGE,oss,load_most_relaxed_structure_only)){
        pflow::logger(_AFLOW_FILE_NAME_, function_name, "Could not load structure (auid="+entry.auid+") ... skipping...", FileMESSAGE, logstream, _LOGGER_WARNING_);
        continue;
      }
      //DX20200225 - added compare to particular geometry files - START
      bool found_structure = false;
      uint structure_index = 0;
      if(load_most_relaxed_structure_only && entry.vstr.size()==1){
        found_structure = true;
        structure_index = 0;
      }
      else if(!load_most_relaxed_structure_only){
        if(entry.vstr.size()==3 && structure_files.size()==3){
          if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_ORIGINAL_ && 
              (structure_files[0] == "POSCAR.orig" || 
               structure_files[0] == "POSCAR.relax1")){ 
            structure_index = 0; 
            found_structure = true; 
            if(LDEBUG){cerr << function_name << " loaded original structure: " << structure_files[0] << endl;}
          }
          else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_RELAX1_ && 
              (structure_files[1] == "POSCAR.relax2" || 
               structure_files[1] == "CONTCAR.relax1")){ 
            structure_index = 1; 
            found_structure = true; 
            if(LDEBUG){cerr << function_name << " loaded relax1 structure: " << structure_files[1] << endl;}
          }
        }
      }
      //DX20200225 - added compare to particular geometry files - END
      if(found_structure){
        // store entry from database
        StructurePrototype str_proto_tmp;
        deque<string> deque_species; for(uint j=0;j<species.size();j++){deque_species.push_back(species[j]);}
        entry.vstr[structure_index].SetSpecies(deque_species);
        str_proto_tmp.structure_representative = entry.vstr[structure_index];
        str_proto_tmp.structure_representative.ReScale(1.0); //DX20191105
        str_proto_tmp.structure_representative.BringInCell(); //DX20200707
        str_proto_tmp.structure_representative_name=entry.getPathAURL(FileMESSAGE,oss,false); //DX20190321 - changed to false, i.e., do not load from common
        str_proto_tmp.structure_representative.directory=str_proto_tmp.structure_representative_name; //DX20190718 - update xstructure.directory
        str_proto_tmp.structure_representative_generated=true;
        str_proto_tmp.structure_representative_source="aurl";
        str_proto_tmp.structure_representative_relaxation_step=relaxation_step; //DX20200429
        str_proto_tmp.stoichiometry=tmp_reduced_stoich;
        str_proto_tmp.elements=species;
        str_proto_tmp.natoms = entry.vstr[structure_index].atoms.size(); //DX20191031
        str_proto_tmp.ntypes = entry.vstr[structure_index].num_each_type.size(); //DX20191031
        str_proto_tmp.structure_representative_compound = pflow::prettyPrintCompound(str_proto_tmp.elements,str_proto_tmp.stoichiometry,no_vrt,false,txt_ft);
        //DX20191105 [MOVED LATER - SAME AS SYMMETRY] str_proto_tmp.LFA_environments= compare::computeLFAEnvironment(tmp.structure_representative); //DX20190711
        str_proto_tmp.property_names = property_list; //DX20190326
        str_proto_tmp.property_units = property_units; //DX20190326
        // store any properties 
        for(uint l=0;l<properties_response[i].size();l++){
          bool property_requested = false;
          for(uint m=0;m<property_list.size();m++){
            if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
          }
          if(property_requested){
            str_proto_tmp.properties_structure_representative.push_back(properties_response[i][l].second);
          }
        }
        if(LDEBUG){
          cerr << XPID << "compare::compareStructureDirectory() Found structure: " << str_proto_tmp.structure_representative_name << endl;
        }
        all_structures.push_back(str_proto_tmp);
      }
      else{
        message << "More structures loaded than anticipated for auid=" << auids[i] << ".";     
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
      }
      //}
    }
#endif

    message << "Total number of candidate structures loaded: " << all_structures.size(); //DX20190403
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_); //DX20190403

    vector<StructurePrototype> final_prototypes = compare::compareMultipleStructures(all_structures, 
        oss, 
        FileMESSAGE, 
        num_proc, 
        same_species, 
        directory, 
        comparison_options); //DX20200103 - condensed booleans to xoptions 

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
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);

    return oss.str();

  }
}

//DX - COMPARE DATABASE ENTRIES - END
*/

//DX - COMPARE DATABASE ENTRIES - START
// ***************************************************************************
// compare::compareDatabaseEntries - compares database entries
// ***************************************************************************
namespace compare {
  string compareDatabaseEntries(const aurostd::xoption& vpflow, ostream& logstream){ //DX20191125 - added ofstream overload and added ostream as input
    ofstream FileMESSAGE;
    return compareDatabaseEntries(vpflow, FileMESSAGE, logstream);
  }

  string compareDatabaseEntries(const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){ //DX20191125 - added ofstream and ostream
    bool LDEBUG=(FALSE || XHOST.DEBUG);

    string function_name = XPID + "compare::compareDatabaseEntries():";
    string directory = ".";
    stringstream message;
    stringstream oss;

    string usage="aflow --compare_database_entries < POSCAR";
    string options="";

    vector<string> tokens,sub_tokens;
    vector<string> matchbook; //aflux - filter/get properties
    vector<string> schema; //get metadata of properties (e.g., units)
    vector<string> property_units;

    bool same_species = true;
    
    // ---------------------------------------------------------------------------
    // instantiate XtalFinder calculator 
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        1,
        logstream);
    
    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103
    xtal_finder.getOptions(vpflow, comparison_options);

    /*
    // ---------------------------------------------------------------------------
    // FLAG: misfit threshold //DX20201119
    double misfit_match_threshold = DEFAULT_XTALFINDER_MISFIT_MATCH;
    double misfit_family_threshold = DEFAULT_XTALFINDER_MISFIT_FAMILY;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::MISFIT_MATCH")) {
      misfit_match_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::MISFIT_MATCH"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_MATCH",vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::MISFIT_MATCH"));
    }
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::MISFIT_FAMILY")) {
      misfit_family_threshold = aurostd::string2utype<double>(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::MISFIT_FAMILY"));
      comparison_options.push_attached("COMPARISON_OPTIONS::MISFIT_FAMILY",vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::MISFIT_FAMILY"));
    }
    // match threshold must be less than family threshold
    if(misfit_match_threshold>misfit_family_threshold){
      message << "Matching misfit threshold must be less than the same family threshold:"
        << " misfit_match_threshold: " << misfit_match_threshold
        << " misfit_family_threshold: " << misfit_family_threshold;
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ILLEGAL_);
    }
    message << "Misfit theshold for matched structures: " << misfit_match_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_MATCH << ")"; 
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Misfit theshold for structures in the same family: " << misfit_family_threshold << " (default: " << DEFAULT_XTALFINDER_MISFIT_FAMILY << ")";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // FLAG: number of processors (multithreading) 
    uint num_proc=1; //defalut=1
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::NP")) {
      num_proc=aurostd::string2utype<uint>(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::NP"));
    }

    // ---------------------------------------------------------------------------
    // FLAG: optimize match
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::OPTIMIZE_MATCH")) {
      comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",TRUE);
      message << "OPTIONS: Finding optimal match; exploring all possible lattices and origins to find the best match (note: this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: no volume scaling
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::NO_SCALE_VOLUME")) {
      comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
      message << "OPTIONS: Suppressing volume scaling; useful for distinguishing structures at different pressures."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore Wyckoff positions
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::IGNORE_WYCKOFF")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring Wyckoff positions when grouping comparisons, but will group by space group (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore symmetry
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::IGNORE_SYMMETRY")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF",TRUE);
      message << "OPTIONS: Ignoring symmetry when grouping comparisons, i.e., do not group by space group and Wyckoff positions (note: do not use for making prototypes; this will slow down the comparisons)."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: ignore LFA environment analysis
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::IGNORE_ENVIRONMENT_ANALYSIS")) {
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS",TRUE);
      message << "OPTIONS: Ignoring LFA environment analysis when grouping comparisons."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: primitivize structures //DX20201005
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::PRIMITIVIZE")) {
      comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE",TRUE);
      message << "OPTIONS: Converting all structures to a primitive representation.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: Minkowski reduction //DX20201005
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::MINKOWSKI")) {
      comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI",TRUE);
      message << "OPTIONS: Performing Minkowski lattice reduction on all structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: Niggli reduction //DX20201005
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::NIGGLI")) {
      comparison_options.flag("COMPARISON_OPTIONS::NIGGLI",TRUE);
      message << "OPTIONS: Performing Niggli lattice reduction on all structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
*/
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
        //DX20191106 [OBSOLETE - switch to getElements] XATOM_SplitAlloySpecies(alloy_string, species);
        species = aurostd::getElements(alloy_string); //DX20191106
      }
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify space group
    vector<uint> space_groups;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::SPACE_GROUP")){
      string space_group_input = vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::SPACE_GROUP");
      vector<string> space_group_strings;
      aurostd::string2tokens(space_group_input,space_group_strings,",");

      uint sg_tmp = 0;
      for(uint i=0;i<space_group_strings.size();i++){
        sg_tmp = aurostd::string2utype<uint>(space_group_strings[i]);
        if(sg_tmp<1 || sg_tmp>230){
          message << "Invalid space group requested: " << space_group_strings[i] << ". Please check input.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_INPUT_ERROR_);
        }
        space_groups.push_back(sg_tmp);
      }
      message << "OPTIONS: Requesting the following space groups: " << space_group_input;
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

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
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
/*
    // ---------------------------------------------------------------------------
    // FLAG: do not remove unmatched structures from the StructurePrototype Object
    // keeps results of each comparison
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::KEEP_UNMATCHED")) {
      comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE);
    }

    // ---------------------------------------------------------------------------
    // FLAG: remove duplicate compounds (useful for non-biased statistics)
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::REMOVE_DUPLICATE_COMPOUNDS")) {
      comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",TRUE);
      message << "OPTIONS: Remove duplicate compounds first, useful for non-biased prototype statistics."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: match unique structures to the AFLOW prototypes 
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::MATCH_TO_AFLOW_PROTOS")) {
      comparison_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS",TRUE);
      message << "OPTIONS: Compare unique structures to the AFLOW prototypes."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: get AFLOW ANRL designation for unique structures
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::ADD_AFLOW_PROTOTYPE_DESIGNATION")) {
      comparison_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION",TRUE);
      message << "OPTIONS: Cast unique structures into AFLOW standard designation."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
*/
    // ---------------------------------------------------------------------------
    // FLAG: consider magnetic structure in comparison //DX20191212
    bool magnetic_comparison=false;
    vector<string> magmoms_for_systems;
    if(vpflow.flag("COMPARE::MAGNETIC")){
      magnetic_comparison=true;
      message << "OPTIONS: Including magnetic moment information in comparisons.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    if(magnetic_comparison){} //CO20200508 - keep it busy
/*
    // ---------------------------------------------------------------------------
    // FLAG: do not calculate unique atom decorations 
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS")) {
      comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",FALSE);
      message << "OPTIONS: Do not calculate unique atom decorations."; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
*/
    // ---------------------------------------------------------------------------
    // FLAG: property list to extract from database (using AFLUX)
    vector<string> property_list;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::PROPERTY_LIST")) {
      aurostd::string2tokens(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::PROPERTY_LIST"),property_list,",");

      // put properties in schema and matchbook for AFLUX call
      schema.push_back("schema("+vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::PROPERTY_LIST")+")"); //to get units
      matchbook.insert(matchbook.end(), property_list.begin(), property_list.end());

      message << "OPTIONS: Extracting the following properties: " << aurostd::joinWDelimiter(property_list,", "); 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify the geometry file to grab (orig, relax1, relax2, static, bands, POSCAR, CONTCAR)
    //DX TODO 
    string geometry_file = "";
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::GEOMETRY_FILE")) {
      geometry_file = vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::GEOMETRY_FILE");
      message << "OPTIONS: Structure type (POSCAR.orig, POSCAR.relax1, POSCAR.relax2, CONTCAR.relax1, ...): " << geometry_file << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // ---------------------------------------------------------------------------
    // FLAG: specify the relaxation step to grab (orig, relax1, relax2, static, bands, POSCAR, CONTCAR)
    uint relaxation_step = _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_;
    bool load_most_relaxed_structure_only = true; 
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")) {
      if(aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")), "orig") || 
          vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP") == "0" ){
        relaxation_step = _COMPARE_DATABASE_GEOMETRY_ORIGINAL_;
        load_most_relaxed_structure_only = false;
      }
      if(aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")), "relax1") || 
          aurostd::substring2bool(aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP")), "middle_relax") || 
          vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::RELAXATION_STEP") == "1"){
        relaxation_step = _COMPARE_DATABASE_GEOMETRY_RELAX1_;
        load_most_relaxed_structure_only = false;
      }
      message << "OPTIONS: Relaxation step (0=original, 1=relax1, 2=most_relaxed): " << relaxation_step << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    if(load_most_relaxed_structure_only){} //CO20200508 - keep it busy

    // ---------------------------------------------------------------------------
    // FLAG: catalog (icsd, lib1, lib2, lib3, ...)
    string catalog = "";
    string catalog_summons = "";
    //DX20200331 [OBSOLETE] bool ICSD_comparison = false;
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::CATALOG")) {
      catalog = aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::CATALOG")); //DX20190718
      catalog_summons = "catalog(\'" + catalog + "\')";
      matchbook.push_back(catalog_summons);
      message << "OPTIONS: Catalog/library (icsd, lib1, lib2, lib3, ...): " << catalog << endl; 
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    if(catalog=="" || catalog=="icsd" || catalog=="all"){ //DX20191108 - needs to be outside of loop
      comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON",TRUE);
      //DX20200331 [OBSOLETE] ICSD_comparison=true;
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

    // ---------------------------------------------------------------------------
    // AFLUX matchbook preparations: get space group(s) //DX20200929
    string sg_summons = "";
    if(space_groups.size()!=0){
      sg_summons = aflowlib::getSpaceGroupAFLUXSummons(space_groups, relaxation_step);
      if(LDEBUG){ cerr << "Space group summons: " << sg_summons << endl; }
    }

    matchbook.push_back(species_summons);
    matchbook.push_back(nspecies_summons);
    matchbook.push_back(sg_summons);

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
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // call AFLUX 
    string response = aflowlib::AFLUXCall(Summons);

    message << "Number of entries returned: " << aurostd::string2tokens(response,tokens,"\n");
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

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
    // get schema from xoptions, i.e., metadata (for the units) (DX20191230)
    string schema_unit = "";
    for(uint p=0;p<property_list.size();p++){
      schema_unit = "SCHEMA::UNIT:" + aurostd::toupper(property_list[p]);
      property_units.push_back(XHOST.vschema.getattachedscheme(schema_unit));
    }

    //DX20191230 [OBSOLETE - FROM AFLUX]    // ---------------------------------------------------------------------------
    //DX20191230 [OBSOLETE - FROM AFLUX]    // get AFLUX schema, i.e., metadata (for the units)
    //DX20191230 [OBSOLETE - FROM AFLUX]    if(schema.size()>0){
    //DX20191230 [OBSOLETE - FROM AFLUX]      schema.push_back(format);
    //DX20191230 [OBSOLETE - FROM AFLUX]      schema.push_back(paging);
    //DX20191230 [OBSOLETE - FROM AFLUX]
    //DX20191230 [OBSOLETE - FROM AFLUX]      // call AFLUX to get schema
    //DX20191230 [OBSOLETE - FROM AFLUX]      response = aflowlib::AFLUXCall(schema);
    //DX20191230 [OBSOLETE - FROM AFLUX]      vector<vector<std::pair<string,string> > > schema_response = aflowlib::getPropertiesFromAFLUXResponse(response);
    //DX20191230 [OBSOLETE - FROM AFLUX]
    //DX20191230 [OBSOLETE - FROM AFLUX]      // extract units
    //DX20191230 [OBSOLETE - FROM AFLUX]      for(uint i=0;i<schema_response.size();i++){
    //DX20191230 [OBSOLETE - FROM AFLUX]        bool units_found = false;
    //DX20191230 [OBSOLETE - FROM AFLUX]        for(uint j=0;j<schema_response[i].size();j++){
    //DX20191230 [OBSOLETE - FROM AFLUX]          if(schema_response[i][j].first=="units"){
    //DX20191230 [OBSOLETE - FROM AFLUX]            property_units.push_back(schema_response[i][j].second);
    //DX20191230 [OBSOLETE - FROM AFLUX]            units_found=true;
    //DX20191230 [OBSOLETE - FROM AFLUX]            break;
    //DX20191230 [OBSOLETE - FROM AFLUX]          }
    //DX20191230 [OBSOLETE - FROM AFLUX]        }
    //DX20191230 [OBSOLETE - FROM AFLUX]        if(!units_found){
    //DX20191230 [OBSOLETE - FROM AFLUX]          property_units.push_back("");
    //DX20191230 [OBSOLETE - FROM AFLUX]        }
    //DX20191230 [OBSOLETE - FROM AFLUX]      }
    //DX20191230 [OBSOLETE - FROM AFLUX]      if(LDEBUG){
    //DX20191230 [OBSOLETE - FROM AFLUX]        for(uint i=0;i<property_units.size();i++){ cerr << function_name << ": units for " << property_list[i] << ": " << property_units[i] << endl; }
    //DX20191230 [OBSOLETE - FROM AFLUX]      }
    //DX20191230 [OBSOLETE - FROM AFLUX]    }

    message << "Total number of candidate structures from database: " << auids.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Loading structures ..." << auids.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    
    // ---------------------------------------------------------------------------
    // load and store entries from the database 
    for(uint i=0; i<auids.size(); i++){
      // first, get stoichiometry from entry
      //DX20191106 [OBSOLETE - switch to getElements] vector<string> species; vector<double> natoms;
      //DX20191106 [OBSOLETE - switch to getElements] XATOM_SplitAlloySpecies(compounds[i], species, natoms);
      vector<double> vcomposition;
      vector<string> species = aurostd::getElements(compounds[i], vcomposition);
      if(LDEBUG){cerr << function_name << " species=" << aurostd::joinWDelimiter(species,",") << endl;}
      vector<uint> tmp_stoich;
      //DX20191106 [OBSOLETE - switch to getElements] for(uint j=0;j<natoms.size();j++)
      for(uint j=0;j<vcomposition.size();j++) //DX20191106
      { //CO20200106 - patching for auto-indenting
        if(aurostd::isinteger(vcomposition[j])){
          tmp_stoich.push_back((uint)aurostd::nint(vcomposition[j]));
        }
        else {
          message << "Expected natoms in " << auids[i] << " to be an integer.";     
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
      }

      //DX20191125 [OBSOLETE] vector<uint> tmp_reduced_stoich = compare::gcdStoich(tmp_stoich);
      vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
      //DX20190402 - need to sort if ignoring species - START
      if(!same_species){
        //DX20191125 [OBSOLETE - REDUNDANT] for(uint i=0; i<tmp_reduced_stoich.size(); i++){
        std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
        //DX20191125 [OBSOLETE - REDUNDANT] }
      }
      //DX20190402 - need to sort if ignoring species - END
      // second, check if stoichiometries are compatible
      // note: do not include in AFLUX matchbook, we would need to specify a range of compatible stoichs (could be expensive)
      // instead: filter on stoichiometry after recieving AFLUX response
      //if(compare::sameStoichiometry(stoichiometry,tmp_reduced_stoich)){
        aflowlib::_aflowlib_entry entry; entry.auid=auids[i]; entry.aurl=aurls[i]; 
        vector<string> structure_files;
        if(!pflow::loadXstructures(entry,structure_files,FileMESSAGE,logstream,load_most_relaxed_structure_only)){
          pflow::logger(_AFLOW_FILE_NAME_, function_name, "Could not load structure (auid="+entry.auid+") ... skipping...", FileMESSAGE, logstream, _LOGGER_WARNING_);
          continue;
        }
        bool found_structure = false;
        uint structure_index = 0;
        if(load_most_relaxed_structure_only && entry.vstr.size()==1){
          found_structure = true;
          structure_index = 0;
        }
        else if(!load_most_relaxed_structure_only){
          if(entry.vstr.size()==3 && structure_files.size()==3){
            if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_ORIGINAL_ && 
                (structure_files[0] == "POSCAR.orig" || 
                 structure_files[0] == "POSCAR.relax1")){ 
              structure_index = 0; 
              found_structure = true; 
              if(LDEBUG){cerr << function_name << " loaded original structure: " << structure_files[0] << endl;}
            }
            else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_RELAX1_ && 
                (structure_files[1] == "POSCAR.relax2" || 
                 structure_files[1] == "CONTCAR.relax1")){ 
              structure_index = 1; 
              found_structure = true; 
              if(LDEBUG){cerr << function_name << " loaded relax1 structure: " << structure_files[1] << endl;}
            }
          }
        }
        if(found_structure){
          // store entry from database
          deque<string> deque_species; for(uint j=0;j<species.size();j++){deque_species.push_back(species[j]);}
          entry.vstr[structure_index].SetSpecies(deque_species);
          string str_path = entry.getPathAURL(FileMESSAGE,logstream,false);
          xtal_finder.addStructureToContainer(entry.vstr[structure_index], str_path, "aurl", relaxation_step, same_species);
          /*
          StructurePrototype str_proto_tmp;
          deque<string> deque_species; for(uint j=0;j<species.size();j++){deque_species.push_back(species[j]);}
          entry.vstr[structure_index].SetSpecies(deque_species);
          str_proto_tmp.structure_representative = entry.vstr[structure_index];
          str_proto_tmp.natoms = entry.vstr[structure_index].atoms.size(); //DX20200421
          str_proto_tmp.ntypes = entry.vstr[structure_index].num_each_type.size(); //DX20200421
          str_proto_tmp.structure_representative.ReScale(1.0); //DX20191105
          str_proto_tmp.structure_representative.BringInCell(); //DX20200707
          str_proto_tmp.structure_representative_name=entry.getPathAURL(FileMESSAGE,oss,false); //DX20190321 - changed to false, i.e., do not load from common
          str_proto_tmp.structure_representative.directory=str_proto_tmp.structure_representative_name; //DX20190718 - update xstructure.directoryr
          str_proto_tmp.structure_representative_generated=true;
          str_proto_tmp.structure_representative_source="aurl";
          str_proto_tmp.structure_representative_relaxation_step=relaxation_step; //DX20200429
          str_proto_tmp.stoichiometry=tmp_reduced_stoich;
          str_proto_tmp.elements=species; //DX20200903 - needs to be before prettyPrintCompound()
          str_proto_tmp.structure_representative_compound = pflow::prettyPrintCompound(str_proto_tmp.elements,str_proto_tmp.stoichiometry,no_vrt,false,txt_ft);
          //DX20191105 [MOVED LATER - SAME AS SYMMETRY] str_proto_tmp.LFA_environments= compare::computeLFAEnvironment(str_proto_tmp.structure_representative); //DX20190711
          str_proto_tmp.natoms = entry.vstr[structure_index].atoms.size(); //DX20191031
          str_proto_tmp.ntypes = entry.vstr[structure_index].num_each_type.size(); //DX20191031
          */
          // store any properties 
          for(uint l=0;l<properties_response[i].size();l++){
            bool property_requested = false;
            for(uint m=0;m<property_list.size();m++){
              if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
            }
            if(property_requested){
              //str_proto_tmp.properties_structure_representative.push_back(properties_response[i][l].second);
              xtal_finder.structure_containers.back().properties.push_back(properties_response[i][l].second);
              xtal_finder.structure_containers.back().properties_names = property_list;
              xtal_finder.structure_containers.back().properties_units = property_units;
              //structure_containers.back().properties_types = property_types;
            }
          }
          //if(LDEBUG) {
          //  cerr << XPID << "compare::compareStructureDirectory() Found structure: " << str_proto_tmp.structure_representative_name << endl;
          //}
          //all_structures.push_back(str_proto_tmp);
        }
        else {
          message << "More structures loaded than anticipated for auid=" << auids[i] << " (# structures=" << entry.vstr.size() << ").";     
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
      //}
    }

    /*
#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
    message << "Splitting into threads...";     
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    // ---------------------------------------------------------------------------
    // distribute threads via indices
    uint number_of_structures = auids.size();
    //DX20191210 [URL TIMEOUT] uint num_threads = aurostd::min(num_proc,number_of_structures); // cannot have more threads than structures
    uint num_threads = 1; //DX20191210 [URL TIMEOUT]
    //DX20191107 [switching to getThreadDistribution] - vector<uint> start_indices, end_indices;
    //DX20191107 [switching to getThreadDistribution] - compare::splitTaskIntoThreads(number_of_structures,num_threads,start_indices,end_indices);
    vector<vector<int> > thread_distribution = getThreadDistribution(number_of_structures, num_threads); //DX20191107 
    message << "Done. Split into threads.";     
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // initialize vector of objects 
    vector<StructurePrototype> all_structures;
    for(uint i=0; i<auids.size(); i++){
      StructurePrototype str_proto_tmp;
      // first, get stoichiometry from entry
      //DX20191106 [OBSOLETE - switch to getElements] vector<string> species; vector<double> natoms;
      //DX20191106 [OBSOLETE - switch to getElements] XATOM_SplitAlloySpecies(compounds[i], species, natoms);
      vector<double> vcomposition;
      vector<string> species = aurostd::getElements(compounds[i], vcomposition);
      vector<uint> tmp_stoich;
      //DX20191106 [OBSOLETE - switch to getElements] for(uint j=0;j<natoms.size();j++)
      for(uint j=0;j<vcomposition.size();j++) //DX20191106
      { //CO20200106 - patching for auto-indenting
        if(aurostd::isinteger(vcomposition[j])){
          tmp_stoich.push_back((uint)aurostd::nint(vcomposition[j]));
        }
        else{
          message << "Expected natoms in " << auids[i] << " to be an integer.";     
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
      }

      //DX20191125 [OBSOLETE] vector<uint> tmp_reduced_stoich = compare::gcdStoich(tmp_stoich);
      vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
      //DX20190402 - need to sort if ignoring species - START
      if(!same_species){
        //DX20191125 [OBSOLETE - REDUNDANT]for(uint i=0; i<tmp_reduced_stoich.size(); i++){
        std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
        //DX20191125 [OBSOLETE - REDUNDANT] }
      }
      str_proto_tmp.stoichiometry=tmp_reduced_stoich;
      str_proto_tmp.elements=species;
      str_proto_tmp.structure_representative_name=aurls[i];
      str_proto_tmp.structure_representative_source="aurl";
      str_proto_tmp.structure_representative_relaxation_step=relaxation_step; //DX20200429
      all_structures.push_back(str_proto_tmp);
    }
    message << "Finished initializing StructurePrototype object, now spawn threads.";     
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // Run threads (now using pointer to thread)
    vector<std::thread*> threads;
    for(uint n=0; n<num_threads; n++){
      //DX20191107 [OBOSLETE - switch to getThreadDistribution convention] threads.push_back(std::thread(compare::generateStructures,std::ref(all_structures),std::ref(oss),start_indices[n],end_indices[n]));
      threads.push_back(new std::thread(&compare::generateStructures, std::ref(all_structures),std::ref(logstream),thread_distribution[n][0],thread_distribution[n][1])); //DX20191107 //DX20200225 - oss to logstream
    }
    // ---------------------------------------------------------------------------
    // Join threads
    for(uint t=0;t<num_threads;t++){
      threads[t]->join();
      delete threads[t];
    }
    // ---------------------------------------------------------------------------
    // try to regenerate via one thread; sometimes multithreading does not load 
    // all structures (timeout issue)
    if(num_threads!=1){
      compare::generateStructures(all_structures);
    }
    message << "Threads complete. " << all_structures.size() << " structures. Adding properties.";     
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // add properties information 
    for(uint i=0;i<all_structures.size();i++){
      all_structures[i].property_names = property_list; //DX20190326
      all_structures[i].property_units = property_units; //DX20190326
      // store any properties 
      for(uint l=0;l<properties_response[i].size();l++){
        bool property_requested = false;
        for(uint m=0;m<property_list.size();m++){
          if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
        }
        if(property_requested){
          all_structures[i].properties_structure_representative.push_back(properties_response[i][l].second);
        }
      }
    }

    for(uint i=0;i<all_structures.size();i++){
      if(all_structures[i].structure_representative_generated){
        deque<string> deque_species; for(uint j=0;j<all_structures[i].elements.size();j++){deque_species.push_back(all_structures[i].elements[j]);}
        all_structures[i].structure_representative.SetSpecies(deque_species);
        all_structures[i].structure_representative_compound = pflow::prettyPrintCompound(all_structures[i].elements,all_structures[i].stoichiometry,no_vrt,false,txt_ft);
      }
    }

    message << "Properties added, now removing non-generated structures" << endl;     
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    // ---------------------------------------------------------------------------
    // remove non-generated structures
    compare::removeNonGeneratedStructures(all_structures);

    for(uint i=0;i<all_structures.size();i++){
      all_structures[i].structure_representative.ReScale(1.0); //DX20191105
      all_structures[i].structure_representative.BringInCell(); //DX20200707
      //DX20191105 [MOVED LATER - SAME AS SYMMETRY] all_structures[i].LFA_environments= compare::computeLFAEnvironment(all_structures[i].structure_representative); //DX20190711
    }

#else
    // ---------------------------------------------------------------------------
    // load and store entries from the database 
    vector<StructurePrototype> all_structures;
    for(uint i=0; i<auids.size(); i++){
      // first, get stoichiometry from entry
      //DX20191106 [OBSOLETE - switch to getElements] vector<string> species; vector<double> natoms;
      //DX20191106 [OBSOLETE - switch to getElements] XATOM_SplitAlloySpecies(compounds[i], species, natoms);
      vector<double> vcomposition;
      vector<string> species = aurostd::getElements(compounds[i], vcomposition);
      vector<uint> tmp_stoich;
      //DX20191106 [OBSOLETE - switch to getElements] for(uint j=0;j<natoms.size();j++){
      for(uint j=0;j<vcomposition.size();j++){ //DX20191106
        if(aurostd::isinteger(vcomposition[j])){
          tmp_stoich.push_back((uint)aurostd::nint(vcomposition[j]));
        }
        else{
          message << "Expected natoms in " << auids[i] << " to be an integer.";     
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
      }

      //DX20191125 [OBSOLETE} vector<uint> tmp_reduced_stoich = compare::gcdStoich(tmp_stoich);
      vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
      //DX20190402 - need to sort if ignoring species - START
      if(!same_species){
        //DX20191125 [OBSOLETE - REDUNDANT] for(uint i=0; i<tmp_reduced_stoich.size(); i++){
        std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
        //DX20191125 [OBSOLETE - REDUNDANT] }
      }
      //DX20190402 - need to sort if ignoring species - END
      // second, check if stoichiometries are compatible
      // note: do not include in AFLUX matchbook, we would need to specify a range of compatible stoichs (could be expensive)
      // instead: filter on stoichiometry after recieving AFLUX response
      //if(compare::sameStoichiometry(stoichiometry,tmp_reduced_stoich)){
      aflowlib::_aflowlib_entry entry; entry.auid=auids[i]; entry.aurl=aurls[i]; 
      vector<string> structure_files;
      if(!pflow::loadXstructures(entry,structure_files,FileMESSAGE,oss,load_most_relaxed_structure_only)){
        pflow::logger(_AFLOW_FILE_NAME_, function_name, "Could not load structure (auid="+entry.auid+") ... skipping...", FileMESSAGE, logstream, _LOGGER_WARNING_);
        continue;
      }
      //DX20200225 - added compare to particular geometry files - START
      bool found_structure = false;
      uint structure_index = 0;
      if(load_most_relaxed_structure_only && entry.vstr.size()==1){
        found_structure = true;
        structure_index = 0;
      }
      else if(!load_most_relaxed_structure_only){
        if(entry.vstr.size()==3 && structure_files.size()==3){
          if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_ORIGINAL_ && 
              (structure_files[0] == "POSCAR.orig" || 
               structure_files[0] == "POSCAR.relax1")){ 
            structure_index = 0; 
            found_structure = true; 
            if(LDEBUG){cerr << function_name << " loaded original structure: " << structure_files[0] << endl;}
          }
          else if(relaxation_step == _COMPARE_DATABASE_GEOMETRY_RELAX1_ && 
              (structure_files[1] == "POSCAR.relax2" || 
               structure_files[1] == "CONTCAR.relax1")){ 
            structure_index = 1; 
            found_structure = true; 
            if(LDEBUG){cerr << function_name << " loaded relax1 structure: " << structure_files[1] << endl;}
          }
        }
      }
      //DX20200225 - added compare to particular geometry files - END
      if(found_structure){
        // store entry from database
        StructurePrototype str_proto_tmp;
        deque<string> deque_species; for(uint j=0;j<species.size();j++){deque_species.push_back(species[j]);}
        entry.vstr[structure_index].SetSpecies(deque_species);
        str_proto_tmp.structure_representative = entry.vstr[structure_index];
        str_proto_tmp.structure_representative.ReScale(1.0); //DX20191105
        str_proto_tmp.structure_representative.BringInCell(); //DX20200707
        str_proto_tmp.structure_representative_name=entry.getPathAURL(FileMESSAGE,oss,false); //DX20190321 - changed to false, i.e., do not load from common
        str_proto_tmp.structure_representative.directory=str_proto_tmp.structure_representative_name; //DX20190718 - update xstructure.directory
        str_proto_tmp.structure_representative_generated=true;
        str_proto_tmp.structure_representative_source="aurl";
        str_proto_tmp.structure_representative_relaxation_step=relaxation_step; //DX20200429
        str_proto_tmp.stoichiometry=tmp_reduced_stoich;
        str_proto_tmp.elements=species;
        str_proto_tmp.natoms = entry.vstr[structure_index].atoms.size(); //DX20191031
        str_proto_tmp.ntypes = entry.vstr[structure_index].num_each_type.size(); //DX20191031
        str_proto_tmp.structure_representative_compound = pflow::prettyPrintCompound(str_proto_tmp.elements,str_proto_tmp.stoichiometry,no_vrt,false,txt_ft);
        //DX20191105 [MOVED LATER - SAME AS SYMMETRY] str_proto_tmp.LFA_environments= compare::computeLFAEnvironment(tmp.structure_representative); //DX20190711
        str_proto_tmp.property_names = property_list; //DX20190326
        str_proto_tmp.property_units = property_units; //DX20190326
        // store any properties 
        for(uint l=0;l<properties_response[i].size();l++){
          bool property_requested = false;
          for(uint m=0;m<property_list.size();m++){
            if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
          }
          if(property_requested){
            str_proto_tmp.properties_structure_representative.push_back(properties_response[i][l].second);
          }
        }
        if(LDEBUG){
          cerr << XPID << "compare::compareStructureDirectory() Found structure: " << str_proto_tmp.structure_representative_name << endl;
        }
        all_structures.push_back(str_proto_tmp);
      }
      else{
        message << "More structures loaded than anticipated for auid=" << auids[i] << ".";     
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
      }
      //}
    }
#endif
*/
    message << "Total number of candidate structures loaded: " << xtal_finder.structure_containers.size(); //DX20190403
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_); //DX20190403

    vector<StructurePrototype> final_prototypes = xtal_finder.compareMultipleStructures(xtal_finder.num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions 

    // ---------------------------------------------------------------------------
    // print results 
    stringstream ss_out;
    xtal_finder.printResults(ss_out, same_species, final_prototypes, "text");
    stringstream ss_json;
    xtal_finder.printResults(ss_json, same_species, final_prototypes, "json");

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
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);

    return oss.str();

  }
}
//DX - COMPARE DATABASE ENTRIES - END
/*
//DX20190424 START
// ***************************************************************************
// compare::compareStructuresFromStructureList()
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compareStructuresFromStructureList(vector<string>& filenames, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, const aurostd::xoption& comparison_options){ //DX20200103 - condensed booleans to xoptions

    // ---------------------------------------------------------------------------
    // directory to write results
    string directory = "."; // for now this is fixed

    // ---------------------------------------------------------------------------
    // load structures appended to command
    vector<StructurePrototype> all_structures = compare::loadStructuresFromStructureList(filenames, magmoms_for_systems, same_species, FileMESSAGE); //DX20190319 - added FileMESSAGE 

    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    return compare::compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions 

  }
}
//DX20190424 END
*/

//DX20190424 START
// ***************************************************************************
// XtalFinderCalculator::compareStructuresFromStructureList()
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareStructuresFromStructureList(const vector<string>& filenames, vector<string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options){ //DX20200103 - condensed booleans to xoptions

    // ---------------------------------------------------------------------------
    // directory to write results
    string directory = "."; // for now this is fixed

    // ---------------------------------------------------------------------------
    // load structures appended to command
    loadStructuresFromStructureList(
        filenames,
        magmoms_for_systems,
        same_species);

    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    return compareMultipleStructures(num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions 

}
//DX20190424 END

/*
// ***************************************************************************
// compare::compareStructuresFromDirectory()
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compareStructuresFromDirectory(string& directory, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, const aurostd::xoption& comparison_options){ //DX20200103 - condensed booleans to xoptions

    // ---------------------------------------------------------------------------
    // load structures in directory 
    vector<StructurePrototype> all_structures = compare::loadStructuresFromDirectory(directory, magmoms_for_systems, same_species, FileMESSAGE); //DX20190319 - added FileMESSAGE

    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    return compare::compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions 

  }
}
*/

// ***************************************************************************
// XtalFinderCalculator::compareStructuresFromDirectory()
// ***************************************************************************
  vector<StructurePrototype> XtalFinderCalculator::compareStructuresFromDirectory(const string& directory, vector<string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options){

    // ---------------------------------------------------------------------------
    // load structures in directory 
    loadStructuresFromDirectory(
        directory,
        magmoms_for_systems,
        same_species);

    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    return compareMultipleStructures(num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions 

  }

/*
// ***************************************************************************
// compare::compareStructuresFromFile()
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compareStructuresFromFile(string& filename, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, const aurostd::xoption& comparison_options){ //DX20200103 - condensed booleans to xoptions

    // ---------------------------------------------------------------------------
    // directory to write results
    string directory = "."; // for now this is fixed

    // ---------------------------------------------------------------------------
    // load structures in file
    vector<StructurePrototype> all_structures = compare::loadStructuresFromFile(filename, magmoms_for_systems, same_species, FileMESSAGE); //DX20190319 - added FileMESSAGE

    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    return compare::compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions

  }
}
*/

// ***************************************************************************
// XtalFinderCalculator::compareStructuresFromFile()
// ***************************************************************************
  vector<StructurePrototype> XtalFinderCalculator::compareStructuresFromFile(const string& filename, vector<string>& magmoms_for_systems, uint num_proc, bool same_species, const aurostd::xoption& comparison_options){ //DX20200103 - condensed booleans to xoptions

    // ---------------------------------------------------------------------------
    // directory to write results
    string directory = "."; // for now this is fixed

    // ---------------------------------------------------------------------------
    // load structures in file
    loadStructuresFromFile(filename,
        magmoms_for_systems,
        same_species);

    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    return compareMultipleStructures(num_proc, same_species, directory, comparison_options);

  }

//TO DO
/*
// ***************************************************************************
// compare::getUniqueEntries() //DX20201111
// ***************************************************************************
namespace compare {
  vector<aflowlib::_aflowlib_entry> getUniqueEntries(vector<aflowlib::_aflowlib_entry>& entries, uint num_proc, bool same_species, bool scale_volume, bool optimize_match){
    ostringstream oss;
    ofstream FileMESSAGE;
    return getUniqueEntries(entries, oss, FileMESSAGE, num_proc, same_species, scale_volume, optimize_match);
  }
}

namespace compare {
  vector<aflowlib::_aflowlib_entry> getUniqueEntries(vector<aflowlib::_aflowlib_entry>& entries, ostream& oss, ofstream& FileMESSAGE, uint num_proc, bool same_species, bool scale_volume, bool optimize_match){

    // ---------------------------------------------------------------------------
    // load structures from aflowlib entries
    vector<string> magmoms_for_systems; //DX20201111 - not included for now
    vector<StructurePrototype> all_structures = compare::loadStructuresFromAflowlibEntries(entries, magmoms_for_systems, same_species, FileMESSAGE); //DX20190319 - added FileMESSAGE

    // ---------------------------------------------------------------------------
    // directory to write results
    string directory = "."; // for now this is fixed

    // ---------------------------------------------------------------------------
    // default comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions();
    comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",scale_volume);
    comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",optimize_match);

    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    vector<StructurePrototype> grouped_structures = compare::compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions

    // ---------------------------------------------------------------------------
    // filter out duplicate aflowlib entries
    vector<aflowlib::_aflowlib_entry> entries_unique;
    for(uint i=0;i<grouped_structures.size();i++){
      for(uint j=0;j<entries.size();j++){
        if(grouped_structures[i].structure_representative_name==entries[j].auid){
          entries_unique.push_back(entries[j]);
          break;
        }
      }
    }

    return entries_unique;

  }
}
*/

/*
// ***************************************************************************
// compare::compareMultipleStructures()
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& all_structures, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, string& directory, ostream& logstream){ //DX20190319 - added FileMESSAGE //DX20200226 - added logstream

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); 
    if(!same_species){
      comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",TRUE);
    }

    return compareMultipleStructures(all_structures, oss, FileMESSAGE, num_proc, same_species, directory, comparison_options, logstream);   //CO20200508

  }
}

namespace compare {
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& all_structures, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, string& directory, const aurostd::xoption& comparison_options, ostream& logstream){ //DX20200103 - condensed booleans to xoptions //DX20200226 - added logstream

    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = XPID + "compare::compareMultipleStructures():";
    //DX20200226 ostream& logstream = cout;
    stringstream message;
    bool quiet = false;
    //DX20190319 [OBSOLETE] ofstream FileMESSAGE;

    if(LDEBUG){cerr << function_name << " BEGIN" << endl;}  //CO20200508

    message << "Total number of structures to compare: " << all_structures.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // convert to structures to certain representations //DX20201006
    // conversion type(s) is indicated in the comparison_options flag
    if(comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE") ||
        comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI") ||
        comparison_options.flag("COMPARISON_OPTIONS::NIGGLI")){
      message << "Converting structures standard representation (primitive, Minkowski, and/or Niggli).";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      compare::convertStructures(all_structures,comparison_options,num_proc);
      message << "All structures converted.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }

    // ---------------------------------------------------------------------------
    // calculate symmetries of structures
    // if already calculated, do not recalculate
    bool all_symmetries_calculated = true;
    //for(uint i=0;i<all_structures.size();i++){ all_symmetries_calculated*=all_structures[i].isSymmetryCalculated(); } //DX20200810 - gcc-10 warnings
    for(uint i=0;i<all_structures.size();i++){ all_symmetries_calculated = (all_symmetries_calculated&&all_structures[i].isSymmetryCalculated()); } //DX20200810

    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && !all_symmetries_calculated){
      message << "Calculating the symmetry of the structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      compare::calculateSymmetries(all_structures,num_proc); 
      message << "Symmetries calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }  
    else if(!all_symmetries_calculated){
      for(uint i=0; i<all_structures.size(); i++){
        all_structures[i].Pearson = "xX";
        all_structures[i].space_group = 0;
        vector<GroupedWyckoffPosition> vGWyckoffPos_tmp;
        all_structures[i].grouped_Wyckoff_positions = vGWyckoffPos_tmp;
      }
    }

    // ---------------------------------------------------------------------------
    // calculate LFA environments of  database entries 
    // if already calculated, do not recalculate
    bool all_environments_calculated = true;
    //for(uint i=0;i<all_structures.size();i++){ all_environments_calculated*=all_structures[i].isLFAEnvironmentCalculated(); } //DX20200925 - gcc-10 warnings
    for(uint i=0;i<all_structures.size();i++){ all_environments_calculated = (all_environments_calculated&&all_structures[i].isLFAEnvironmentCalculated()); } //DX20200925 - gcc-10 warnings
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS") && !all_environments_calculated){
      message << "Calculating the environments of the structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

      compare::calculateLFAEnvironments(all_structures,num_proc); 

      message << "Environments calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }

    // ---------------------------------------------------------------------------
    // remove duplicate compounds first; uses recursion of this function
    if(comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS")){
      bool tmp_same_species = true;
      aurostd::xoption remove_duplicates_options = comparison_options;
      remove_duplicates_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
      remove_duplicates_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",FALSE);
      remove_duplicates_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",FALSE);
      remove_duplicates_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS",FALSE);
      remove_duplicates_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION",FALSE);

      message << "Comparing to remove duplicate materials.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      vector<StructurePrototype> unique_compounds = compareMultipleStructures(all_structures, 
          oss, 
          FileMESSAGE, 
          num_proc, 
          tmp_same_species, 
          directory,
          remove_duplicates_options); //DX20200103 - condensed booleans to xoptions

      message << "Duplicate materials removed.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

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
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
      }
      else if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
        aurostd::stringstream2file(ss_out_remove_duplicates,directory+"/duplicate_compounds_output.out");
        message << "RESULTS: See " << directory << "/duplicate_compounds_output.out for list of unique/duplicate structures.";
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
      }
      else {
        aurostd::stringstream2file(ss_json_remove_duplicates,directory+"/duplicate_compounds_output.json");
        aurostd::stringstream2file(ss_out_remove_duplicates,directory+"/duplicate_compounds_output.out");
        message << "RESULTS: See " << directory << "/duplicate_compounds_output.out" << " or " << directory << "/duplicate_compounds_output.json" << " for list of unique/duplicate structures.";
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
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
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    vector<StructurePrototype> comparison_schemes = compare::groupStructurePrototypes(all_structures, 
        same_species, 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320 - added environment angles
        false); //DX20200103 - condensed booleans to xoptions

    message << "Number of comparison groups: " << comparison_schemes.size() << ".";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // if ICSD comparison, make structure with minimum ICSD number the representative structure
    if(comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON")){
      compare::representativePrototypeForICSDRuns(comparison_schemes);
    }

    // ---------------------------------------------------------------------------
    // compare structures 
    message << "Running comparisons ...";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes = compare::runComparisonScheme(comparison_schemes, same_species, num_proc, comparison_options, oss, FileMESSAGE, quiet, logstream); //DX20200103 - condensed booleans to xoptions 

    // ---------------------------------------------------------------------------
    // return if there are no similar structures
    if(final_prototypes.size()==0){
      return final_prototypes;
    }
    comparison_schemes.clear();

    // ---------------------------------------------------------------------------
    // combine prototypes regardless of having different space groups
    // BETA TESTING (perhaps we shouldn't do this) - compare::checkPrototypes(num_proc,same_species,final_prototypes);

    message << "Number of unique prototypes: " << final_prototypes.size() << " (out of " << all_structures.size() << " structures).";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);

    // ---------------------------------------------------------------------------
    // get unique atom decorations prototype (representative) structures
    if(!same_species && comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS")){ 
      message << "Determining the unique atom decorations for each prototype.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      // find unique atom decorations of prototype
      for(uint i=0;i<final_prototypes.size();i++){
        if(arePermutationsComparableViaComposition(final_prototypes[i].stoichiometry) && 
            arePermutationsComparableViaSymmetry(final_prototypes[i].grouped_Wyckoff_positions)){
          // check if xstructure is generated; if not, make it
          if(!final_prototypes[i].structure_representative_generated){
            if(!compare::generateStructure(final_prototypes[i].structure_representative_name,final_prototypes[i].structure_representative_source,final_prototypes[i].structure_representative_relaxation_step,final_prototypes[i].structure_representative,oss)){ //DX20200429
              message << "Could not generate structure (" << final_prototypes[i].structure_representative_name << ").";
              throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
            }
          }
          vector<StructurePrototype> final_permutations = compare::comparePermutations(final_prototypes[i],num_proc,comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH"),oss,FileMESSAGE); //DX20200103 - condensed booleans to xoptions

          // store permutation results in main StructurePrototype object
          for(uint j=0;j<final_permutations.size();j++){
            vector<string> permutations_tmp; 
            permutations_tmp.push_back(final_permutations[j].structure_representative_name); //push back representative permutation
            for(uint d=0;d<final_permutations[j].structures_duplicate_names.size();d++){ permutations_tmp.push_back(final_permutations[j].structures_duplicate_names[d]); } //push back equivalent permutations
            final_prototypes[i].atom_decorations_equivalent.push_back(permutations_tmp);
          }
          final_permutations.clear(); //DX20190624
        }
        else{
          vector<string> unique_permutations; generatePermutationString(final_prototypes[i].stoichiometry, unique_permutations); //DX20191125
          // store permutation results in main StructurePrototype object
          for(uint j=0;j<unique_permutations.size();j++){
            vector<string> permutations_tmp; permutations_tmp.push_back(unique_permutations[j]);
            final_prototypes[i].atom_decorations_equivalent.push_back(permutations_tmp);
          }
        }
      }
    }

    // ---------------------------------------------------------------------------
    // for structure-type comparisons, remove duplicate compounds (avoid biased duplicate statistics)
    if(same_species==false && comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS")){
      //DX20190220 [BETA]  message << "Performing comparisons to removing duplicate compounds";
      //DX20190220 [BETA]  pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      //DX20190220 [BETA]
      //DX20190220 [BETA]  vector<StructurePrototype> duplicate_compound_comparisons = compare::compareDuplicateCompounds(final_prototypes, num_proc, comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISONS"), oss);
      //DX20190220 [BETA]  
      //DX20190220 [BETA]  message << "Removing duplicates from final comparison list";
      //DX20190220 [BETA]  pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      //DX20190220 [BETA]  // remove duplicates compounds from final comparison results
      //DX20190220 [BETA]  compare::removeDuplicateCompounds(final_prototypes, duplicate_compound_comparisons);
      //DX20190220 [BETA]
      //DX20190220 [BETA]  // prepare JSON output
      //DX20190220 [BETA]  stringstream ss_json_remove_duplicates;
      //DX20190220 [BETA]  compare::printResults(ss_json_remove_duplicates, same_species, final_prototypes, "json");
      //DX20190220 [BETA]  
      //DX20190220 [BETA]  // prepare TEXT (.out) output
      //DX20190220 [BETA]  stringstream ss_out_remove_duplicates;
      //DX20190220 [BETA]  compare::printResults(ss_out_remove_duplicates, same_species, final_prototypes, "txt");
      //DX20190220 [BETA]  
      //DX20190220 [BETA]  // write results to files
      //DX20190220 [BETA]  aurostd::stringstream2file(ss_json_remove_duplicates,directory+"/structure_comparison_no_duplicate_compounds_output.json");
      //DX20190220 [BETA]  aurostd::stringstream2file(ss_out_remove_duplicates,directory+"/structure_comparison_no_duplicate_compounds_output.out");
      //DX20190220 [BETA]  message << "RESULTS: See " << directory << "/structure_comparison_no_duplicate_compounds_output.out" << " or " << directory << "/structure_comparison_no_duplicate_compounds_output.json" << " for list of unique/duplicate structures.";
      //DX20190220 [BETA]  pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }

    if(comparison_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION")){
      message << "Determining the AFLOW standard designation.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
      // ---------------------------------------------------------------------------
      // split task into threads 
      uint number_of_structures = final_prototypes.size();
      uint number_of_threads = aurostd::min(num_proc,number_of_structures); // cannot have more threads than structures
      //DX20191107 [switching to getThreadDistribution] vector<uint> start_indices, end_indices;
      //DX20191107 [switching to getThreadDistribution] splitTaskIntoThreads(number_of_structures, number_of_threads, start_indices, end_indices);
      vector<vector<int> > thread_distribution = getThreadDistribution(number_of_structures, number_of_threads); //DX20191107 

      // ---------------------------------------------------------------------------
      // [THREADED] determine AFLOW standard designation 
      vector<std::thread*> threads;
      for(uint n=0; n<number_of_threads; n++){
        //DX20191107 [switching to getThreadDistribution convention] threads.push_back(std::thread(compare::getPrototypeDesignations,std::ref(final_prototypes),start_indices[n], end_indices[n]));
        threads.push_back(new std::thread(&compare::getPrototypeDesignations,std::ref(final_prototypes),thread_distribution[n][0], thread_distribution[n][1])); //DX20191107
      }
      for(uint t=0;t<threads.size();t++){
        threads[t]->join();
        delete threads[t];
      }

      // ---------------------------------------------------------------------------
      // update once all are collected (safer) 
      for(uint i=0;i<final_prototypes.size();i++){
        final_prototypes[i].aflow_label = final_prototypes[i].structure_representative.prototype;
        final_prototypes[i].aflow_parameter_list = final_prototypes[i].structure_representative.prototype_parameter_list;
        final_prototypes[i].aflow_parameter_values = final_prototypes[i].structure_representative.prototype_parameter_values;
      }
#else
      // ---------------------------------------------------------------------------
      // [NON-THREADED] determine AFLOW standard designation 
      for(uint i=0;i<final_prototypes.size();i++){
        anrl::structure2anrl(final_prototypes[i].structure_representative,false); //DX20190829 - false for recalculate_symmetry
        final_prototypes[i].aflow_label = final_prototypes[i].structure_representative.prototype;
        final_prototypes[i].aflow_parameter_list = final_prototypes[i].structure_representative.prototype_parameter_list;
        final_prototypes[i].aflow_parameter_values = final_prototypes[i].structure_representative.prototype_parameter_values;
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
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
    }

    if(comparison_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS")){
      message << "Determining if representative structures map to any of the AFLOW prototypes.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

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
      // do not calculate unique atom decorations since this was already done
      vpflow_protos.flag("COMPARE2PROTOTYPES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS",TRUE);

      // ---------------------------------------------------------------------------
      // match to AFLOW prototypes 
      for(uint i=0;i<final_prototypes.size();i++){
        vector<StructurePrototype> matching_protos = compare::compare2prototypes(final_prototypes[i].structure_representative, vpflow_protos);
        final_prototypes[i].matching_aflow_prototypes = matching_protos[0].structures_duplicate_names;
      }
    }

    return final_prototypes;
  }
}
*/

// ***************************************************************************
// compare::compareMultipleStructures()
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareMultipleStructures(
    uint num_proc, bool same_species,
    const string& directory){ //DX20190319 - added FileMESSAGE //DX20200226 - added logstream

    // ---------------------------------------------------------------------------
    // create xoptions to contain all comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); 
    if(!same_species){
      comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",TRUE);
    }

    return compareMultipleStructures(num_proc, same_species, directory, comparison_options);   //CO20200508

  }

vector<StructurePrototype> XtalFinderCalculator::compareMultipleStructures(
    uint num_proc,
    bool same_species,
    const string& directory,
    const aurostd::xoption& comparison_options){ //DX20200103 - condensed booleans to xoptions //DX20200226 - added logstream

    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = XPID + "XtalFinderCalculator::compareMultipleStructures():";
    //DX20200226 ostream& logstream = cout;
    stringstream message;
    bool quiet = false;
    //DX20190319 [OBSOLETE] ofstream FileMESSAGE;

    if(LDEBUG){cerr << function_name << " BEGIN" << endl;}  //CO20200508

    message << "Total number of structures to compare: " << structure_containers.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // convert to structures to certain representations //DX20201006
    // conversion type(s) is indicated in the comparison_options flag
    if(comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE") ||
        comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI") ||
        comparison_options.flag("COMPARISON_OPTIONS::NIGGLI")){
      message << "Converting structures standard representation (primitive, Minkowski, and/or Niggli).";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      //compare::convertStructures(all_structures,comparison_options,num_proc);
      convertStructures(comparison_options,num_proc);
      message << "All structures converted.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    }

    // ---------------------------------------------------------------------------
    // calculate symmetries of structures
    // if already calculated, do not recalculate
    bool all_symmetries_calculated = true;
    //for(uint i=0;i<all_structures.size();i++){ all_symmetries_calculated*=all_structures[i].isSymmetryCalculated(); } //DX20200810 - gcc-10 warnings
    //for(uint i=0;i<all_structures.size();i++){ all_symmetries_calculated = (all_symmetries_calculated&&all_structures[i].isSymmetryCalculated()); } //DX20200810
    for(uint i=0;i<structure_containers.size();i++){ all_symmetries_calculated = (all_symmetries_calculated&&isSymmetryCalculated(structure_containers[i])); } //DX20200810

    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && !all_symmetries_calculated){
      message << "Calculating the symmetry of the structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      calculateSymmetries(num_proc); 
      message << "Symmetries calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    }  
    else if(!all_symmetries_calculated){
      for(uint i=0; i<structure_containers.size(); i++){
        structure_containers[i].Pearson = "xX";
        structure_containers[i].space_group = 0;
        vector<GroupedWyckoffPosition> vGWyckoffPos_tmp;
        structure_containers[i].grouped_Wyckoff_positions = vGWyckoffPos_tmp;
      }
    }

    // ---------------------------------------------------------------------------
    // calculate LFA environments of  database entries 
    // if already calculated, do not recalculate
    bool all_environments_calculated = true;
    //for(uint i=0;i<all_structures.size();i++){ all_environments_calculated*=all_structures[i].isLFAEnvironmentCalculated(); } //DX20200925 - gcc-10 warnings
    //for(uint i=0;i<all_structures.size();i++){ all_environments_calculated = (all_environments_calculated&&all_structures[i].isLFAEnvironmentCalculated()); } //DX20200925 - gcc-10 warnings
    for(uint i=0;i<structure_containers.size();i++){ all_environments_calculated = (all_environments_calculated&&isLFAEnvironmentCalculated(structure_containers[i])); } //DX20200810
    if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS") && !all_environments_calculated){
      message << "Calculating the environments of the structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

      calculateLFAEnvironments(num_proc); 

      message << "Environments calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    }
    
    // ---------------------------------------------------------------------------
    // get nearest neighbor info, can perhaps only calculate this if necessary
    // (i.e., if we know we will perform a comparison)
    // if already calculated, do not recalculate
    bool all_neighbors_calculated = true;
    //for(uint i=0;i<all_structures.size();i++){ all_neighbors_calculated = (all_neighbors_calculated&&(all_structures[i].structure_representative_struct->nearest_neighbor_distances.size()!=0)); } //DX20200925 - gcc-10 warnings
    for(uint i=0;i<structure_containers.size();i++){ all_neighbors_calculated = (all_neighbors_calculated&&(structure_containers[i].nearest_neighbor_distances.size()!=0)); } //DX20200925 - gcc-10 warnings
    if(!all_neighbors_calculated){
      message << "Calculating the nearest neighbors of all the structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

      getNearestNeighbors(num_proc); 

      message << "Nearest neighbors information calculated.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    }

    // ---------------------------------------------------------------------------
    // remove duplicate compounds first; uses recursion of this function
    if(comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS")){
      bool tmp_same_species = true;
      aurostd::xoption remove_duplicates_options = comparison_options;
      remove_duplicates_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",FALSE);
      remove_duplicates_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",FALSE);
      remove_duplicates_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS",FALSE);
      remove_duplicates_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS",FALSE);
      remove_duplicates_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION",FALSE);

      message << "Comparing to remove duplicate materials.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      vector<StructurePrototype> unique_compounds = compareMultipleStructures(
          num_proc, 
          tmp_same_species, 
          directory,
          remove_duplicates_options); //DX20200103 - condensed booleans to xoptions

      message << "Duplicate materials removed.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

      // include duplicate compounds count in object
      for(uint i=0;i<unique_compounds.size();i++){
        unique_compounds[i].structure_representative_struct->number_compounds_matching_structure = unique_compounds[i].numberOfComparisons();
      }

      // prepare JSON output
      stringstream ss_json_remove_duplicates;
      printResults(ss_json_remove_duplicates, tmp_same_species, unique_compounds, "json");

      // prepare TEXT (.out) output
      stringstream ss_out_remove_duplicates;
      printResults(ss_out_remove_duplicates, tmp_same_species, unique_compounds, "txt");

      // write results to files
      if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
        aurostd::stringstream2file(ss_json_remove_duplicates,directory+"/duplicate_compounds_output.json");
        message << "RESULTS: See " << directory << "/duplicate_compounds_output.json for list of unique/duplicate structures.";
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
      }
      else if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
        aurostd::stringstream2file(ss_out_remove_duplicates,directory+"/duplicate_compounds_output.out");
        message << "RESULTS: See " << directory << "/duplicate_compounds_output.out for list of unique/duplicate structures.";
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
      }
      else {
        aurostd::stringstream2file(ss_json_remove_duplicates,directory+"/duplicate_compounds_output.json");
        aurostd::stringstream2file(ss_out_remove_duplicates,directory+"/duplicate_compounds_output.out");
        message << "RESULTS: See " << directory << "/duplicate_compounds_output.out" << " or " << directory << "/duplicate_compounds_output.json" << " for list of unique/duplicate structures.";
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
      }

      // ---------------------------------------------------------------------------
      // remove duplicates in structure container
      vector<string> duplicates_name;
      for(uint i=0;i<unique_compounds.size();i++){
        for(uint j=0;j<unique_compounds[i].structures_duplicate_struct.size();j++){
          duplicates_name.push_back(unique_compounds[i].structures_duplicate_struct[j]->name);
        }
      }
      for(uint i=0;i<duplicates_name.size();i++){
        removeStructureFromContainerByName(duplicates_name[i]);
      }
    
    }

    // ---------------------------------------------------------------------------
    // group structures based on stoichiometry and symmetry (unless ignoring symmetry/Wyckoff)
    message << "Grouping sets of comparisons.";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    vector<StructurePrototype> comparison_schemes = groupStructurePrototypes( 
        same_species, 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"), 
        comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320 - added environment angles
        false); //DX20200103 - condensed booleans to xoptions

    message << "Number of comparison groups: " << comparison_schemes.size() << ".";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    // ---------------------------------------------------------------------------
    // if ICSD comparison, make structure with minimum ICSD number the representative structure
    if(comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON")){
      representativePrototypeForICSDRunsNEW(comparison_schemes);
    }

    // ---------------------------------------------------------------------------
    // compare structures 
    message << "Running comparisons ...";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes = runComparisonScheme(comparison_schemes, same_species, num_proc, comparison_options, quiet); //DX20200103 - condensed booleans to xoptions 

    //cerr << "DONE:" << endl;
    // ---------------------------------------------------------------------------
    // return if there are no similar structures
    if(final_prototypes.size()==0){
      return final_prototypes;
    }
    comparison_schemes.clear();

    // ---------------------------------------------------------------------------
    // combine prototypes regardless of having different space groups
    // BETA TESTING (perhaps we shouldn't do this) - compare::checkPrototypes(num_proc,same_species,final_prototypes);

    //message << "Number of unique prototypes: " << final_prototypes.size() << " (out of " << all_structures.size() << " structures).";
    message << "Number of unique prototypes: " << final_prototypes.size() << " (out of " << structure_containers.size() << " structures).";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);

    
    // ---------------------------------------------------------------------------
    // get unique atom decorations prototype (representative) structures
    if(!same_species && comparison_options.flag("COMPARISON_OPTIONS::CALCULATE_UNIQUE_PERMUTATIONS")){ 
      message << "Determining the unique atom decorations for each prototype.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      // find unique atom decorations of prototype
      for(uint i=0;i<final_prototypes.size();i++){
        if(compare::arePermutationsComparableViaComposition(final_prototypes[i].stoichiometry) && 
           (comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") || 
            compare::arePermutationsComparableViaSymmetry(final_prototypes[i].grouped_Wyckoff_positions))){
          // check if xstructure is generated; if not, make it
          if(!final_prototypes[i].structure_representative_struct->is_structure_generated){
            if(!compare::generateStructure(final_prototypes[i].structure_representative_struct->name,final_prototypes[i].structure_representative_struct->source,final_prototypes[i].structure_representative_struct->relaxation_step,final_prototypes[i].structure_representative_struct->structure,*p_oss)){ //DX20200429
              message << "Could not generate structure (" << final_prototypes[i].structure_representative_struct->name << ").";
              throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
            }
          }
          XtalFinderCalculator xtal_finder_permutations;
          vector<StructurePrototype> final_permutations = xtal_finder_permutations.comparePermutations(final_prototypes[i],num_proc,comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH"));
          // store permutation results in main StructurePrototype object
          for(uint j=0;j<final_permutations.size();j++){
            vector<string> permutations_tmp; 
            permutations_tmp.push_back(final_permutations[j].structure_representative_struct->name); //push back representative permutation
            for(uint d=0;d<final_permutations[j].structures_duplicate_struct.size();d++){ permutations_tmp.push_back(final_permutations[j].structures_duplicate_struct[d]->name); } //push back equivalent permutations
            final_prototypes[i].atom_decorations_equivalent.push_back(permutations_tmp);
          }
          final_permutations.clear(); //DX20190624
        }
        else{
          vector<string> unique_permutations; compare::generatePermutationString(final_prototypes[i].stoichiometry, unique_permutations); //DX20191125
          // store permutation results in main StructurePrototype object
          for(uint j=0;j<unique_permutations.size();j++){
            vector<string> permutations_tmp; permutations_tmp.push_back(unique_permutations[j]);
            final_prototypes[i].atom_decorations_equivalent.push_back(permutations_tmp);
          }
        }
      }
    }

    // ---------------------------------------------------------------------------
    // for structure-type comparisons, remove duplicate compounds (avoid biased duplicate statistics)
    if(same_species==false && comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS")){
      //DX20190220 [BETA]  message << "Performing comparisons to removing duplicate compounds";
      //DX20190220 [BETA]  pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      //DX20190220 [BETA]
      //DX20190220 [BETA]  vector<StructurePrototype> duplicate_compound_comparisons = compare::compareDuplicateCompounds(final_prototypes, num_proc, comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISONS"), oss);
      //DX20190220 [BETA]  
      //DX20190220 [BETA]  message << "Removing duplicates from final comparison list";
      //DX20190220 [BETA]  pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
      //DX20190220 [BETA]  // remove duplicates compounds from final comparison results
      //DX20190220 [BETA]  compare::removeDuplicateCompounds(final_prototypes, duplicate_compound_comparisons);
      //DX20190220 [BETA]
      //DX20190220 [BETA]  // prepare JSON output
      //DX20190220 [BETA]  stringstream ss_json_remove_duplicates;
      //DX20190220 [BETA]  compare::printResults(ss_json_remove_duplicates, same_species, final_prototypes, "json");
      //DX20190220 [BETA]  
      //DX20190220 [BETA]  // prepare TEXT (.out) output
      //DX20190220 [BETA]  stringstream ss_out_remove_duplicates;
      //DX20190220 [BETA]  compare::printResults(ss_out_remove_duplicates, same_species, final_prototypes, "txt");
      //DX20190220 [BETA]  
      //DX20190220 [BETA]  // write results to files
      //DX20190220 [BETA]  aurostd::stringstream2file(ss_json_remove_duplicates,directory+"/structure_comparison_no_duplicate_compounds_output.json");
      //DX20190220 [BETA]  aurostd::stringstream2file(ss_out_remove_duplicates,directory+"/structure_comparison_no_duplicate_compounds_output.out");
      //DX20190220 [BETA]  message << "RESULTS: See " << directory << "/structure_comparison_no_duplicate_compounds_output.out" << " or " << directory << "/structure_comparison_no_duplicate_compounds_output.json" << " for list of unique/duplicate structures.";
      //DX20190220 [BETA]  pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    }

    if(comparison_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION")){
      message << "Determining the AFLOW standard designation.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
      // ---------------------------------------------------------------------------
      // split task into threads 
      uint number_of_structures = final_prototypes.size();
      uint number_of_threads = aurostd::min(num_proc,number_of_structures); // cannot have more threads than structures
      //DX20191107 [switching to getThreadDistribution] vector<uint> start_indices, end_indices;
      //DX20191107 [switching to getThreadDistribution] splitTaskIntoThreads(number_of_structures, number_of_threads, start_indices, end_indices);
      vector<vector<int> > thread_distribution = getThreadDistribution(number_of_structures, number_of_threads); //DX20191107 

      // ---------------------------------------------------------------------------
      // [THREADED] determine AFLOW standard designation 
      vector<std::thread*> threads;
      for(uint n=0; n<number_of_threads; n++){
        //DX20191107 [switching to getThreadDistribution convention] threads.push_back(std::thread(compare::getPrototypeDesignations,std::ref(final_prototypes),start_indices[n], end_indices[n]));
        threads.push_back(new std::thread(&XtalFinderCalculator::getPrototypeDesignations,this,std::ref(final_prototypes),thread_distribution[n][0], thread_distribution[n][1])); //DX20191107
      }
      for(uint t=0;t<threads.size();t++){
        threads[t]->join();
        delete threads[t];
      }

      // ---------------------------------------------------------------------------
      // update once all are collected (safer) 
      for(uint i=0;i<final_prototypes.size();i++){
        final_prototypes[i].aflow_label = final_prototypes[i].structure_representative_struct->structure.prototype;
        final_prototypes[i].aflow_parameter_list = final_prototypes[i].structure_representative_struct->structure.prototype_parameter_list;
        final_prototypes[i].aflow_parameter_values = final_prototypes[i].structure_representative_struct->structure.prototype_parameter_values;
      }
#else
      // ---------------------------------------------------------------------------
      // [NON-THREADED] determine AFLOW standard designation 
      for(uint i=0;i<final_prototypes.size();i++){
        anrl::structure2anrl(final_prototypes[i].structure_representative_struct->structure,false); //DX20190829 - false for recalculate_symmetry
        final_prototypes[i].aflow_label = final_prototypes[i].structure_representative_struct->structure.prototype;
        final_prototypes[i].aflow_parameter_list = final_prototypes[i].structure_representative_struct->structure.prototype_parameter_list;
        final_prototypes[i].aflow_parameter_values = final_prototypes[i].structure_representative_struct->structure.prototype_parameter_values;
      }
#endif
    }


    // ---------------------------------------------------------------------------
    // for testing/development; in case the subsequent analyses fails, checkpoint file 
    bool store_checkpoint=false;
    if(store_checkpoint){
      stringstream ss_json;
      printResults(ss_json, same_species, final_prototypes, "json");
      stringstream ss_out;
      printResults(ss_out, same_species, final_prototypes, "text");
      aurostd::stringstream2file(ss_json,directory+"/structure_comparison_output.json");
      aurostd::stringstream2file(ss_out,directory+"/structure_comparison_output.out");
      message << "RESULTS: See [tmp]" << directory << "/structure_comparison_output.out" << " or " << directory << "/structure_comparison_output.json" << " for list of unique/duplicate structures.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    }

    if(comparison_options.flag("COMPARISON_OPTIONS::MATCH_TO_AFLOW_PROTOS")){
      message << "Determining if representative structures map to any of the AFLOW prototypes.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

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
      // do not calculate unique atom decorations since this was already done
      vpflow_protos.flag("COMPARE2PROTOTYPES::DO_NOT_CALCULATE_UNIQUE_PERMUTATIONS",TRUE);

      // ---------------------------------------------------------------------------
      // match to AFLOW prototypes 
      for(uint i=0;i<final_prototypes.size();i++){
        XtalFinderCalculator xtal_finder_protos(misfit_match,misfit_family,*p_FileMESSAGE,num_proc,*p_oss);
        //vector<StructurePrototype> matching_protos = compare::compare2prototypes(final_prototypes[i].structure_representative, vpflow_protos);
        vector<StructurePrototype> matching_protos = xtal_finder_protos.compare2prototypes(final_prototypes[i].structure_representative_struct->structure, vpflow_protos);
        for(uint j=0;j<matching_protos[0].structures_duplicate_struct.size();j++){
          final_prototypes[i].matching_aflow_prototypes.push_back(matching_protos[0].structures_duplicate_struct[j]->name);
        }
      }
    }

    //vector<StructurePrototype> final_prototypes;
    return final_prototypes;
  }


// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species) { //DX20191108 - remove const & from bools //DX20191122 - move ostream to end
    uint num_proc=1;
    double final_misfit=AUROSTD_MAX_DOUBLE;
    bool scale_volume=true; //default is true
    bool optimize_match=false; //default is false
    ostringstream comparison_log; //DX20191202 
    return aflowCompareStructure(num_proc, xstr1, xstr2, same_species, scale_volume, optimize_match, final_misfit, comparison_log); //DX20191122 - move ostream to end
  }
}

namespace compare {
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool optimize_match) { //DX20191108 - remove const & from bools //DX20191122 - move ostream to end
    uint num_proc = 1;
    double final_misfit = AUROSTD_MAX_DOUBLE;
    ostringstream comparison_log; //DX20191202 
    return aflowCompareStructure(num_proc, xstr1, xstr2, same_species, scale_volume, optimize_match, final_misfit, comparison_log); //DX20191122 - move ostream to end and add default
  }
}


// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  double aflowCompareStructureMisfit(const xstructure& xstr1, const xstructure& xstr2, bool same_species) { //DX20191108 - remove const & from bools
    uint num_proc=1;
    double final_misfit=AUROSTD_MAX_DOUBLE;
    bool scale_volume=true; //default is true
    bool optimize_match=false; //default is false
    ostringstream comparison_log; //DX20191202 
    aflowCompareStructure(num_proc, xstr1, xstr2, same_species, scale_volume, optimize_match, final_misfit, comparison_log); //DX20191122 - move ostream to end
    return final_misfit;
  }
}

// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const uint& num_proc, const xstructure& xstr1, const xstructure& xstr2, 
      bool same_species, bool scale_volume, bool optimize_match, 
      double& final_misfit, ostream& comparison_log) {

    structure_misfit final_misfit_info = compare::initialize_misfit_struct(); //DX20191218
    

    return aflowCompareStructure(num_proc, xstr1, xstr2, same_species, scale_volume, optimize_match, final_misfit, final_misfit_info, comparison_log); //DX20191122 - move ostream to end
  }
}

// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const uint& num_proc, const xstructure& xstr1, const xstructure& xstr2, 
      bool same_species, bool scale_volume, bool optimize_match, 
      double& final_misfit, structure_misfit& final_misfit_info, ostream& comparison_log) { //DX20191108 - remove const & from bools //DX20191122 - move ostream to end and add default

    _structure_representative str_rep = compare::initializeStructureRepresentativeStruct(xstr1); 
    _structure_representative str_matched = compare::initializeStructureRepresentativeStruct(xstr2); 
    
    // clean structures
    str_rep.structure.FixLattices();
    str_matched.structure.FixLattices();
    str_rep.structure.ReScale(1.0);
    str_matched.structure.ReScale(1.0);
    str_rep.structure.BringInCell(); //DX20190329 - need to ensure incell; otherwise supercell expansion breaks
    str_matched.structure.BringInCell(); //DX20190329 - need to ensure incell; otherwise supercell expansion breaks

    // ---------------------------------------------------------------------------
    // clean atom names (remove pseudopotential information)
    for(uint i=0;i<str_rep.structure.species.size();i++){ str_rep.structure.species[i]=KBIN::VASP_PseudoPotential_CleanName(str_rep.structure.species[i]); }
    for(uint i=0;i<str_matched.structure.species.size();i++){ str_matched.structure.species[i]=KBIN::VASP_PseudoPotential_CleanName(str_matched.structure.species[i]); }

    for(uint i=0;i<str_rep.structure.atoms.size();i++){ str_rep.structure.atoms[i].name=KBIN::VASP_PseudoPotential_CleanName(str_rep.structure.atoms[i].name); }
    for(uint i=0;i<str_matched.structure.atoms.size();i++){ str_matched.structure.atoms[i].name=KBIN::VASP_PseudoPotential_CleanName(str_matched.structure.atoms[i].name); }

    XtalFinderCalculator xtal_finder;
    xtal_finder.compareStructures(str_rep,str_matched,final_misfit_info,same_species,scale_volume,optimize_match);

    final_misfit = final_misfit_info.misfit;
   /* 
    cerr << "str_matched.name: " << str_matched.structure << endl;;
    cerr << "basis transformation: " << str_matched.basis_transformation << endl;;
    cerr << "rotation: " << str_matched.rotation << endl;;

    print(str_matched.atom_map);
    print(str_matched.basis_map);
    print(str_matched.distances_mapped);
    */
    return(final_misfit<=xtal_finder.misfit_match);

  }
}

// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
void XtalFinderCalculator::compareStructures(
    _structure_representative& str_rep,
    _structure_representative& str_matched,
    structure_misfit& match_info,
    bool same_species,
    bool scale_volume,
    bool optimize_match) {

    // This is the main comparison function, which  compares two crystal structures
    // and determines their level of similarity based on the idea discussed 
    // in H. Burzlaff's paper (Acta Cryst., A53, 217-224 (1997)).

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "XtalFinderCalculator::compareStructures():";

    //DX20201130 [COMPARISON LOG] comparison_log << "==================================================================================" << endl;

    // ---------------------------------------------------------------------------
    // prepare structures (swap structure order if necessary, fix lattices, rescale, etc.)
    //xstructure xstr_base, xstr_test;
    //xstr_base = xstr1; xstr_test = xstr2; //DX20201020

    // clean structures
    //DX20201210 [OBSOLETE] - accounted for earlier - str_rep.structure.FixLattices();
    //DX20201210 [OBSOLETE] - accounted for earlier - str_matched.structure.FixLattices();
    //DX20201210 [OBSOLETE] - accounted for earlier - str_rep.structure.ReScale(1.0);
    //DX20201210 [OBSOLETE] - accounted for earlier - str_matched.structure.ReScale(1.0);
    //DX20201210 [OBSOLETE] - accounted for earlier - str_rep.structure.BringInCell(); //DX20190329 - need to ensure incell; otherwise supercell expansion breaks
    //DX20201210 [OBSOLETE] - accounted for earlier - str_matched.structure.BringInCell(); //DX20190329 - need to ensure incell; otherwise supercell expansion breaks

    // ---------------------------------------------------------------------------
    // clean atom names (remove pseudopotential information)
    //DX20201210 [OBSOLETE] - accounted for earlier - for(uint i=0;i<str_rep.structure.species.size();i++){ str_rep.structure.species[i]=KBIN::VASP_PseudoPotential_CleanName(str_rep.structure.species[i]); }
    //DX20201210 [OBSOLETE] - accounted for earlier - for(uint i=0;i<str_matched.structure.species.size();i++){ str_matched.structure.species[i]=KBIN::VASP_PseudoPotential_CleanName(str_matched.structure.species[i]); }

    //DX20201210 [OBSOLETE] - accounted for earlier - for(uint i=0;i<str_rep.structure.atoms.size();i++){ str_rep.structure.atoms[i].name=KBIN::VASP_PseudoPotential_CleanName(str_rep.structure.atoms[i].name); }
    //DX20201210 [OBSOLETE] - accounted for earlier - for(uint i=0;i<str_matched.structure.atoms.size();i++){ str_matched.structure.atoms[i].name=KBIN::VASP_PseudoPotential_CleanName(str_matched.structure.atoms[i].name); }

    // ---------------------------------------------------------------------------
    // NOW DONE BEFOREHAND standardize structure (not default) 
    // below is no longer necessary, algorithm handles supercells/conventional/prim
    //bool primitivize=false;
    //bool niggli=false;
    //if(primitivize){
    //  str_rep.structure=GetStandardPrimitive(str_rep.structure);
    //  str_matched.structure=GetStandardPrimitive(str_matched.structure);
    //}
    //if(niggli){
    //  str_rep.structure.NiggliUnitCellForm();
    //  str_matched.structure.NiggliUnitCellForm();
    //}

    // ---------------------------------------------------------------------------
    // determine minimum interatomic distances of structures (resolution of atoms) //DX20200623
    // //DX20200715 - may need to rescale this if the structures are being scaled later....
    //if(str_rep.structure.dist_nn_min==AUROSTD_NAN){ cerr << "yo: " << endl; str_rep.structure.dist_nn_min=SYM::minimumDistance(str_rep.structure); }
    //if(str_matched.structure.dist_nn_min==AUROSTD_NAN){ cerr << "yo: " <<endl; str_matched.structure.dist_nn_min=SYM::minimumDistance(str_matched.structure); }
    if(str_rep.structure.dist_nn_min==AUROSTD_NAN){ 
      if(str_rep.nearest_neighbor_distances.size()){ str_rep.structure.dist_nn_min=aurostd::min(str_rep.nearest_neighbor_distances); }
      else{ str_rep.structure.dist_nn_min=SYM::minimumDistance(str_rep.structure); }
    }
    if(str_matched.structure.dist_nn_min==AUROSTD_NAN){ 
      if(str_matched.nearest_neighbor_distances.size()){ str_matched.structure.dist_nn_min=aurostd::min(str_matched.nearest_neighbor_distances); }
      else{ str_matched.structure.dist_nn_min=SYM::minimumDistance(str_matched.structure); }
    }
    
    // ---------------------------------------------------------------------------
    // determine if structures are matchable (same species and/or same stoichiometry)
    if(same_species == true){
      // if atoms are not labeled in either structure; assign fake names
      if(str_rep.structure.atoms.at(0).name == "" || str_matched.structure.atoms.at(0).name == ""){ 
        if(LDEBUG) {cerr << function_name << " Atoms are not labeled. Assigning fake element names." << endl;}
        str_rep.structure.DecorateWithFakeElements();
        str_matched.structure.DecorateWithFakeElements();
      }
    }
    if(compare::matchableSpecies(str_rep.structure,str_matched.structure,same_species)==true){
      //DX20201130 [COMPARISON LOG] comparison_log << "=========================================================" << endl; 

      //DX20201130 [COMPARISON LOG] comparison_log << "STRUCTURE 1: " << endl;  
      //DX20201130 [COMPARISON LOG] comparison_log << str_rep.structure << endl;
      //cerr << str_rep.structure << endl;

      //DX20201130 [COMPARISON LOG] comparison_log << "=========================================================" << endl;

      //DX20201130 [COMPARISON LOG] comparison_log << "STRUCTURE 2: " << endl;
      //DX20201130 [COMPARISON LOG] comparison_log << str_matched.structure << endl;	
      //cerr << str_matched.structure << endl;

      //DX20201130 [COMPARISON LOG] comparison_log << "=========================================================" << endl;

      // ---------------------------------------------------------------------------
      // comparison types
      // type_match determines how the atoms should be matched
      //  1: assigns fake names to atoms (allows for structural comparison regardless of type of atom)
      //  2: uses the names given in POSCAR (allows for structural comparison of material; type of atom necessary)

      //if(same_species==true){ type_match=2; }
      //if(same_species==false){ type_match=1; }	

      // ---------------------------------------------------------------------------
      // variables
      //[CO20200508 - OBSOLETE]uint i=0;
      //xvector<double> origin;
      //xstructure proto;           
      //vector<xstructure> vprotos,vprotos_tmp;
      //vector<vector<uint> > IM1, IM2;
      //vector<vector<double> > vmin_dists;
      //vector<uint> im1, im2;
      //vector<string> PAIR1, PAIR2;
      //structure_misfit min_misfit_info = compare::initialize_misfit_struct();
      //bool magnetic_analysis = (str_rep.structure.atoms[0].spin_is_given || str_rep.structure.atoms[0].noncoll_spin_is_given);
      //min_misfit_info.is_magnetic_misfit=(magnetic_analysis && _CALCULATE_MAGNETIC_MISFIT_); //DX20191218

      //DX20201130 [COMPARISON LOG] comparison_log<<"-------------------------------------------------------"<<endl;

      // ---------------------------------------------------------------------------
      // normalize scaling factors 
      //if(LDEBUG) {cerr << "compare:: " << "Scale structures."<<endl;} 
      // structures should already be scaled to the same scaling factor, below may be redundant
      //DX20201210 [OBSOLETE?] compare::rescaleStructure(str_rep.structure,str_matched.structure);

      // ---------------------------------------------------------------------------
      // scale volumes of structures
      // MOVED INTO LATTICE SEARCH if(scale_volume==true){compare::atomicNumberDensity(str_rep.structure, str_matched.structure, match_info.rescale_factor);}

      // ---------------------------------------------------------------------------
      // assign fake atom names 
      //if(type_match==1){
      //  str_rep.structure.DecorateWithFakeElements();
      //  str_matched.structure.DecorateWithFakeElements();
      //}

      // OBSOLETE THIS PRINTS OUT XSTRUCTURES WITH ATOM ZERO SHIFTED TO ORIGIN...
      // OBSOLETE comparison_log<<"========================================================="<<endl;
      // OBSOLETE comparison_log << str_rep.structure << endl;
      // OBSOLETE comparison_log<<"========================================================="<<endl;
      // OBSOLETE comparison_log << str_matched.structure << endl;		

      // ---------------------------------------------------------------------------
      // print lattice parameters
      //DX20201130 [COMPARISON LOG] printParameters(str_rep.structure,comparison_log);
      //DX20201130 [COMPARISON LOG] printParameters(str_matched.structure,comparison_log);

      //DX20201130 [COMPARISON LOG] comparison_log << "========================================================="<<endl;    
      //DX20201130 [COMPARISON LOG] comparison_log << "QUADRUPLETS METHOD" << endl;

      // ---------------------------------------------------------------------------
      // compare structures
      if(LDEBUG) {cerr << function_name << " Searching for new representation of test structure ..."<<endl;} 
      // creates the threads for checking quadruplets (lattices)
      //DX20190530 - OLD threadGeneration(num_proc,q_base,str_matched.structure,vprotos,str_rep.structure,type_match,optimize_match,minMis,comparison_log);
      //latticeAndOriginSearch(str_rep.structure,str_matched.structure,num_proc,q_base,vprotos,min_misfit_info,type_match,optimize_match,scale_volume,comparison_log); //DX20190530 //DX20200422 - scale_volume added
      //latticeSearch(str_rep.structure,str_matched.structure,num_proc,q_base,vprotos,min_misfit_info,type_match,optimize_match,scale_volume,comparison_log); //DX20190530 //DX20200422 - scale_volume added

      //_structure_representative str_representative = compare::initializeStructureRepresentativeStruct(str_rep.structure); 
      //structure_matched str_matched = compare::initializeStructureMatched(str_matched.structure); 

      //stringstream comparison_log; ////DX20201130 [COMPARISON LOG]
      //compare::latticeSearch(str_rep,str_matched,match_info,type_match,optimize_match,scale_volume,num_proc,comparison_log); //DX20190530 //DX20200422 - scale_volume added
      latticeSearch(str_rep,str_matched,match_info,same_species,optimize_match,scale_volume,num_proc); //DX20190530 //DX20200422 - scale_volume added


      //if(LDEBUG) {cerr << "compare:: " << "Total # of possible matching representations: " << vprotos.size() << endl;}	
      //ORIG final_misfit=min_misfit_info.misfit;
      //ORIG final_misfit_info=min_misfit_info; //DX20191218
      //DXOBS final_misfit=str_matched.misfit_info.misfit;
      //DXOBSfinal_misfit_info=str_matched.misfit_info; //DX20191218

    }
  }
//} //end of compare namespace

// AFLOW-XtalMatch (compare crystal structures)
// Written by David Hicks (david.hicks@duke.edu) 
// Contributors: Carlo De Santo
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow DAVID HICKS - Duke University 2014-2020                 *
// *                                                                         *
// ***************************************************************************
