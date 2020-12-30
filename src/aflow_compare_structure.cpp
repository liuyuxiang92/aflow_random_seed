// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow DAVID HICKS - Duke University 2014-2020                 *
// *                                                                         *
// ***************************************************************************
// AFLOW-XtalFinder (compare crystal structures)
// Written by David Hicks (david.hicks@duke.edu)
// Contributors: Carlo De Santo

#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_compare_structure.h"
#include "aflow_symmetry_spacegroup.h"

#undef AFLOW_COMPARE_MULTITHREADS_ENABLE

#if GCC_VERSION >= 40400   // added two zeros
#define AFLOW_COMPARE_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of AFLOW-XtalFinder will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

// ***************************************************************************
// SEE README_AFLOW_COMPARE.TXT FOR THE FULL LIST OF AFLOW COMMANDS
// ***************************************************************************

// ***************************************************************************
// compare::compareAtomDecorations()
// ***************************************************************************
namespace compare {
  string compareAtomDecorations(istream& input, const aurostd::xoption& vpflow){
    ostringstream results_ss;
    ostringstream oss;
    ofstream FileMESSAGE; //DX20190319 - added FileMESSAGE

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE::USAGE")) {
      string usage="aflow --compare_atom_decorations|--unique_atom_decorations [GENERAL_COMPARISON_OPTIONS] < file";
      string options_function_string = "unique_atom_decorations_options: [--print_misfit]";

      vector<string> options, options_general, options_function;
      aurostd::string2tokens(GENERAL_OPTIONS_LIST,options_general," ");
      aurostd::string2tokens(options_function_string,options_function," ");
      options.push_back(usage);
      options.insert(options.end(), options_general.begin(), options_general.end());
      options.insert(options.end(), options_function.begin(), options_function.end());

      init::ErrorOption("--usage","compare::compareAtomDecorations()",options);
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
    // load structure
    xstructure xstr(input,IOAFLOW_AUTO);

    // ---------------------------------------------------------------------------
    // calculate unique/duplicate permutations
    vector<string> unique_permutations = xtal_finder.getUniquePermutations(xstr, xtal_finder.num_proc, print_misfit, results_ss, comparison_options);

    return results_ss.str();
  }
}


// ***************************************************************************
// XtalFinderCalculator::getUniquePermutations() //DX20201201
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
  aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions("permutation");
  comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",optimize_match);
  stringstream oss;
  return getUniquePermutations(xstr, num_proc, print_misfit, oss, comparison_options);
}

vector<string> XtalFinderCalculator::getUniquePermutations(
    xstructure& xstr,
    uint num_proc,
    bool print_misfit,
    ostream& oss,
    aurostd::xoption& comparison_options){

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
      deque<int> reduced_stoichiometry; aurostd::reduceByGCD(xstr.num_each_type, reduced_stoichiometry); //DX20191125
      deque<uint> reduced_stoichiometry_uint; for(uint i=0;i<reduced_stoichiometry.size(); i++){ reduced_stoichiometry_uint.push_back((uint)reduced_stoichiometry[i]); } //DX20191125
      unique_permutations = getSpeciesPermutedStrings(reduced_stoichiometry_uint); //DX20191125
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
  stringstream xstr_ss; xstr_ss << xstr;
  addStructure2container(xstr, "input geometry", xstr_ss.str(), 0, false); // false: so we can find decorations on systems without atoms
  StructurePrototype structure;
  uint container_index = 0;
  setStructureAsRepresentative(structure,container_index);

  // ---------------------------------------------------------------------------
  // get the unique atom decorations for the structure
  XtalFinderCalculator xtal_finder_permutations;
  xtal_finder_permutations.misfit_match = misfit_match; //copy misfit_match
  xtal_finder_permutations.misfit_family = misfit_family; //copy misfit_family
  xtal_finder_permutations.num_proc = num_proc; //copy num_proc
  vector<StructurePrototype> final_permutations = xtal_finder_permutations.compareAtomDecorations(structure,num_proc,comparison_options);

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

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
    string function_name = XPID + "compare::compareMultipleStructures():";
    ostringstream oss;
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
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions();
    // get options from command line and pass to comparison options //DX20201218
    xtal_finder.getOptions(vpflow, comparison_options);

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE::USAGE")) {
      string usage = "";
      // material-type comparisons
      if(vpflow.flag("COMPARE_MATERIAL")){
        usage="aflow --compare_materials=str1,str2,str3,... | aflow --compare_materials -D <dir_path> | aflow --compare_materials -F=<filename>";
      }
      // structure-type comparisons
      else if(vpflow.flag("COMPARE_STRUCTURE")){
        usage="aflow --compare_structures=str1,str2,str3,... | aflow --compare_structures -D <dir_path> | aflow --compare_structures -F=<filename>";
      }
      vector<string> options, options_general;
      aurostd::string2tokens(GENERAL_OPTIONS_LIST,options_general," ");
      options.push_back(usage);
      options.insert(options.end(), options_general.begin(), options_general.end());
      init::ErrorOption("--usage","compare::compareMultipleStructures()",options);
    }

    // ---------------------------------------------------------------------------
    // distinguish structures coming from directory or file
    string structures_source = ""; // "structure_list", "directory" or "file"

    // from list appended to command
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::STRUCTURE_LIST").empty()){
      structures_source = "structure_list";
    }
    // from directory
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::DIRECTORY").empty()){
      structures_source = "directory";
    }
    // from file
    if(!vpflow.getattachedscheme("COMPARE_STRUCTURE::FILE").empty()){
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

    // ---------------------------------------------------------------------------
    // check if two-structure comparison
    if(file_list.size()==2){
      comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);
      comparison_options.flag("COMPARISON_OPTIONS::CLEAN_UNMATCHED",FALSE);
      comparison_options.flag("COMPARISON_OPTIONS::STORE_COMPARISON_LOGS",TRUE);  //DX20200113 - fixed typo
    }

    // ---------------------------------------------------------------------------
    // load structures
    vector<StructurePrototype> final_prototypes;
    if(structures_source=="structure_list") {
      final_prototypes = xtal_finder.compareStructuresFromStructureList(file_list, magmoms_for_systems, xtal_finder.num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
    }
    else if(structures_source=="directory") {
      final_prototypes = xtal_finder.compareStructuresFromDirectory(directory, magmoms_for_systems, xtal_finder.num_proc, same_species, comparison_options); //DX20200103 - condensed booleans to xoptions
    }
    if(structures_source=="file") {
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
    // write results to files //DX20201229 - consolidated into functions
    if(format=="json"){
      xtal_finder.writeComparisonOutputFile(ss_json, directory, "JSON", "compare_input", same_species);
    }
    else if(format=="text"){
      xtal_finder.writeComparisonOutputFile(ss_out, directory, "TEXT", "compare_input", same_species);
    }
    else if(format=="both"){
      xtal_finder.writeComparisonOutputFile(ss_json, directory, "JSON", "compare_input", same_species);
      xtal_finder.writeComparisonOutputFile(ss_out, directory, "TEXT", "compare_input", same_species);
    }

    return oss.str();
  }
}

// ***************************************************************************
// compare::getIsopointalPrototypes - returns corresponding prototype label
// ***************************************************************************
namespace compare {
  string isopointalPrototypes(istream& input, const aurostd::xoption& vpflow){

    string function_name = "compare::isopointalPrototypes():";

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("ISOPOINTAL_PROTOTYPES::USAGE")) {
      string usage="aflow --get_isopointal_prototypes|--isopointal_prototypes|--get_same_symmetry_prototypes < file";
      string options_function_string = "options: [--catalog=aflow|htqc|all]";
      vector<string> options, options_function;
      aurostd::string2tokens(options_function_string,options_function," ");
      options.push_back(usage);
      options.insert(options.end(), options_function.begin(), options_function.end());

      init::ErrorOption("--usage","compare::isopointalPrototypes()",options);
    }

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
    XtalFinderCalculator xtal_finder; //DX20201218
    vector<StructurePrototype> prototypes = xtal_finder.compare2prototypes(xstr,vpflow);

    // ---------------------------------------------------------------------------
    // global quiet back to default
    XHOST.QUIET=original_quiet;

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
  }
}

// ***************************************************************************
// compare::printMatchingPrototypes - returns list of matching structures
// ***************************************************************************
namespace compare {
  string printMatchingPrototypes(istream& input, const aurostd::xoption& vpflow){

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE::USAGE")) {
      string usage="aflow --compare2prototypes|--compare2protos [GENERAL_COMPARISON_OPTIONS] [COMPARE2PROTOTYPES_OPTIONS] < file";
      string options_function_string = "compare2protos_options: [--catalog=aflow|htqc|all]";

      vector<string> options, options_general, options_function;
      aurostd::string2tokens(GENERAL_OPTIONS_LIST,options_general," ");
      aurostd::string2tokens(options_function_string,options_function," ");
      options.push_back(usage);
      options.insert(options.end(), options_general.begin(), options_general.end());
      options.insert(options.end(), options_function.begin(), options_function.end());

      init::ErrorOption("--usage","compare::printMatchingPrototypes()",options);
    }

    // ---------------------------------------------------------------------------
    // load input structure
    xstructure xstr(input,IOAFLOW_AUTO);

    XtalFinderCalculator xtal_finder;
    return xtal_finder.printMatchingPrototypes(xstr,vpflow);

  }
}

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

// ***************************************************************************
// XtalFinderCalculator::compare2prototypes - identifies corresponding protos
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compare2prototypes(
    const xstructure& xstrIN,
    const aurostd::xoption& vpflow){

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  string function_name = XPID + "XtalFinderCalculator::compare2prototypes():";
  stringstream message;
  bool quiet = false;

  xstructure xstr = xstrIN; //DX20200226 - copy

  // ---------------------------------------------------------------------------
  // create xoptions to contain all comparison options
  aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103
  getOptions(vpflow, comparison_options);

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

  // ---------------------------------------------------------------------------
  // single round of comparisons
  comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);

  // ---------------------------------------------------------------------------
  // add structure to container
  stringstream ss_input; ss_input << xstr;
  addStructure2container(xstr, "input geometry", ss_input.str(), 0, false);

  vector<StructurePrototype> all_structures;

  // ---------------------------------------------------------------------------
  // symmetry
  if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && (xstr.space_group_ITC<1 || xstr.space_group_ITC>230)){ //DX20190829 - don't recalculate symmetry if already calculated //DX20191220 - put range instead of ==0
    calculateSymmetries(1);  //1: one structure -> one processor
  }
  else if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && xstr.space_group_ITC>=1 && xstr.space_group_ITC<=230){ //DX20191220 - put range instead of !=0
    setSymmetryPlaceholders();
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
  vector<StructurePrototype> comparison_schemes = groupStructurePrototypes(
      same_species,
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320 - added environment angles
      false); //DX20200103 - condensed booleans to xoptions

  // ---------------------------------------------------------------------------
  // compare structures
  vector<StructurePrototype> final_prototypes = runComparisonScheme(comparison_schemes, same_species, num_proc, comparison_options, quiet); //DX20200103 - condensed booleans to xoptions

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
      vector<StructurePrototype> final_permutations = xtal_finder_permutations.compareAtomDecorations(final_prototypes[i],num_proc,comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH"));
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
    // instantiate class
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        1,
        logstream);

    // ---------------------------------------------------------------------------
    // main compare2database() function
    vector<StructurePrototype> final_prototypes = xtal_finder.compare2database(xstrIN, vpflow);

    vector<matching_structure> matched_database_structures;

    // ---------------------------------------------------------------------------
    // database DOESN'T contain equivalent structure to input
    if(!final_prototypes[0].structures_duplicate_struct.size()){
      return matched_database_structures;
    }

    // ---------------------------------------------------------------------------
    // return equivalent structures to input
    for(uint i=0;i<final_prototypes[0].structures_duplicate_struct.size();i++){
      matching_structure database_entry;
      database_entry.name = final_prototypes[0].structures_duplicate_struct[i]->name;
      database_entry.misfit = final_prototypes[0].structure_misfits_duplicate[i].misfit;
      matched_database_structures.push_back(database_entry);
    }

    return matched_database_structures;
  }
}

// ***************************************************************************
// XtalFinderCalculator::compare2database - compares database
// ***************************************************************************
// load input structure
vector<StructurePrototype> XtalFinderCalculator::compare2database(
    const xstructure& xstrIN, const aurostd::xoption& vpflow){
  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

  string function_name = XPID + "XtalFinderCalculator::compare2database():";
  string directory = "";
  stringstream message;

  vector<StructurePrototype> final_prototypes; //DX20200225
  xstructure xstr = xstrIN; //copy //DX20200225

  vector<string> tokens,sub_tokens;
  vector<string> matchbook; //aflux - filter/get properties
  vector<string> schema; //get metadata of properties (e.g., units)
  vector<string> property_units;

  bool same_species = true;

  // ---------------------------------------------------------------------------
  // create xoptions to contain all comparison options
  aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions(); //DX20200103
  // get options from vpflow/command line
  getOptions(vpflow,comparison_options); //DX20200103

  // ---------------------------------------------------------------------------
  // single round of comparisons
  comparison_options.flag("COMPARISON_OPTIONS::SINGLE_COMPARISON_ROUND",TRUE);

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
  addStructure2container(xstr, "input geometry", ss_input.str(), 0, false);

  // ---------------------------------------------------------------------------
  // symmetry
  if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && (xstr.space_group_ITC<1 || xstr.space_group_ITC>230)){ //DX20190829 - don't recalculate symmetry if already calculated //DX20191220 - put range instead of ==0
    calculateSymmetries(1);  //1: one structure -> one processor
  }
  else if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && xstr.space_group_ITC>=1 && xstr.space_group_ITC<=230){ //DX20191220 - put range instead of !=0
    setSymmetryPlaceholders();
  }

  if(LDEBUG) {
    cerr << function_name << " Wyckoff positions of input structure:" << endl;
    for(uint i=0;i<structure_containers[0].grouped_Wyckoff_positions.size();i++){
      cerr << structure_containers[0].grouped_Wyckoff_positions[i] << endl;
    }
  }

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
  // (will switch over to CO's AFLUX+AFLOW code when available)
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

  message << "Loading structures ... ";
  pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  vector<StructurePrototype> all_structures;

  // ---------------------------------------------------------------------------
  // load and store entries from the database
  for(uint i=0; i<auids.size(); i++){
    // first, get stoichiometry from entry
    vector<double> vcomposition;
    vector<string> species = aurostd::getElements(compounds[i], vcomposition);
    if(LDEBUG){cerr << function_name << " species=" << aurostd::joinWDelimiter(species,",") << endl;}
    vector<uint> tmp_stoich;
    for(uint j=0;j<vcomposition.size();j++){
      if(aurostd::isinteger(vcomposition[j])){
        tmp_stoich.push_back((uint)aurostd::nint(vcomposition[j]));
      }
      else {
        message << "Expected natoms in " << auids[i] << " to be an integer.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
      }
    }

    vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
    //DX20190402 - need to sort if ignoring species - START
    if(!same_species){
      std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
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
        addStructure2container(entry.vstr[structure_index], str_path, "aurl", relaxation_step, same_species);
        // store any properties
        for(uint l=0;l<properties_response[i].size();l++){
          bool property_requested = false;
          for(uint m=0;m<property_list.size();m++){
            if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
          }
          if(property_requested){
            structure_containers.back().properties.push_back(properties_response[i][l].second);
            structure_containers.back().properties_names = property_list;
            structure_containers.back().properties_units = property_units;
          }
        }
      }
      else {
        message << "More structures loaded than anticipated for auid=" << auids[i] << " (# structures=" << entry.vstr.size() << ").";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
      }
    }
  }
  message << "Total number of candidate structures loaded: " << structure_containers.size(); //DX20190403
  pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_); //DX20190403

  /*
    // ---------------------------------------------------------------------------
    // only compare entries to the input representation, the rest are extraneous comparisons
    vector<StructurePrototype> input_structure_comparison_scheme_only; input_structure_comparison_scheme_only.push_back(comparison_schemes[0]);
    message << "Number of structures to compare to input structure: " << input_structure_comparison_scheme_only[0].structures_duplicate_names.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
*/
  return compareMultipleStructures(num_proc, same_species, directory, comparison_options); //DX20200103 - condensed booleans to xoptions
}

// ***************************************************************************
// compare::compare2database - compares database
// ***************************************************************************
namespace compare {
  // load input structure
  string printCompare2Database(istream& input, const aurostd::xoption& vpflow, ostream& logstream){

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE::USAGE")) {
      string usage="aflow --compare2database [GENERAL_COMPARISON_OPTIONS] [COMPARE2DATABASE_OPTIONS] < file";
      string options_function_string = "compare2database_options: [--catalog=lib1|lib2|lib3|lib4|lib6|lib7|icsd] [--properties=enthalpy_atom,natoms,...] [--relaxation_step=original|relax1|most_relaxed]";

      vector<string> options, options_general, options_function;
      aurostd::string2tokens(GENERAL_OPTIONS_LIST,options_general," ");
      aurostd::string2tokens(options_function_string,options_function," ");
      options.push_back(usage);
      options.insert(options.end(), options_general.begin(), options_general.end());
      options.insert(options.end(), options_function.begin(), options_function.end());

      init::ErrorOption("--usage","compare::printCompare2Database()",options);
    }

    xstructure xstr(input,IOAFLOW_AUTO);
    ofstream FileMESSAGE;
    return printCompare2Database(xstr,vpflow,FileMESSAGE,logstream);
  }  //DX20200225
}

namespace compare {
  string printCompare2Database(istream& input, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){
    xstructure xstr(input,IOAFLOW_AUTO);
    return printCompare2Database(xstr,vpflow,FileMESSAGE,logstream);
  }  //CO20200225
}

namespace compare {
  string printCompare2Database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream){

    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
    string function_name = "compare::printCompare2Database():";
    stringstream message;
    ostringstream oss;

    if(LDEBUG){cerr << function_name << " BEGIN" << endl;}  //CO20200508

    // ---------------------------------------------------------------------------
    // main compare2database() function
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
    // write results to files //DX20201229 - consolidated into functions
    string directory = aurostd::getPWD();
    if(format=="json"){
      xtal_finder_database.writeComparisonOutputFile(ss_json, directory, "JSON", "compare2database", same_species);
    }
    else if(format=="text"){
      xtal_finder_database.writeComparisonOutputFile(ss_out, directory, "TEXT", "compare2database", same_species);
    }
    else if(format=="both"){
      xtal_finder_database.writeComparisonOutputFile(ss_json, directory, "JSON", "compare2database", same_species);
      xtal_finder_database.writeComparisonOutputFile(ss_out, directory, "TEXT", "compare2database", same_species);
    }

    return oss.str();

  }
}

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
    bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);

    string function_name = XPID + "compare::compareDatabaseEntries():";
    string directory = aurostd::getPWD(); //DX20201229
    stringstream message;
    stringstream oss;

    // ---------------------------------------------------------------------------
    // FLAG: usage
    if(vpflow.flag("COMPARE::USAGE")) {
      string usage="aflow --compare_database_entries [GENERAL_COMPARISON_OPTIONS] [COMPARE_DATABASE_ENTRIES_OPTIONS] < file";
      string options_function_string = "compare_database_entries_options: [--alloy=AgAlMn...] [--nspecies=3] [--catalog=lib1|lib2|lib3|lib4|lib6|lib7|icsd] [--properties=enthalpy_atom,natoms,...] [--relaxation_step=original|relax1|most_relaxed] [--space_group=225,186,227,...] [--stoichiometry=1:2:3:...]";

      vector<string> options, options_general, options_function;
      aurostd::string2tokens(GENERAL_OPTIONS_LIST,options_general," ");
      aurostd::string2tokens(options_function_string,options_function," ");
      options.push_back(usage);
      options.insert(options.end(), options_general.begin(), options_general.end());
      options.insert(options.end(), options_function.begin(), options_function.end());

      init::ErrorOption("--usage","compare::compareDatabaseEntries()",options);
    }

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
    string catalog = "", catalog_summons = "";
    if(vpflow.flag("COMPARE_DATABASE_ENTRIES::CATALOG")) {
      catalog = aurostd::tolower(vpflow.getattachedscheme("COMPARE_DATABASE_ENTRIES::CATALOG")); //DX20190718
      catalog_summons = "catalog(\'" + catalog + "\')";
      matchbook.push_back(catalog_summons);
      message << "OPTIONS: Catalog/library (icsd, lib1, lib2, lib3, ...): " << catalog << endl;
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }
    if(catalog=="" || catalog=="icsd" || catalog=="all"){ //DX20191108 - needs to be outside of loop
      comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON",TRUE);
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

    message << "Total number of candidate structures from database: " << auids.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    message << "Loading structures ..." << auids.size();
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    // ---------------------------------------------------------------------------
    // load and store entries from the database
    for(uint i=0; i<auids.size(); i++){
      // first, get stoichiometry from entry
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

      vector<uint> tmp_reduced_stoich; aurostd::reduceByGCD(tmp_stoich, tmp_reduced_stoich); //DX20191125
      //DX20190402 - need to sort if ignoring species - START
      if(!same_species){
        std::sort(tmp_reduced_stoich.begin(),tmp_reduced_stoich.end());
      }
      //DX20190402 - need to sort if ignoring species - END
      // second, check if stoichiometries are compatible
      // note: do not include in AFLUX matchbook, we would need to specify a range of compatible stoichs (could be expensive)
      // instead: filter on stoichiometry after recieving AFLUX response
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
          xtal_finder.addStructure2container(entry.vstr[structure_index], str_path, "aurl", relaxation_step, same_species);
          // store any properties
          for(uint l=0;l<properties_response[i].size();l++){
            bool property_requested = false;
            for(uint m=0;m<property_list.size();m++){
              if(properties_response[i][l].first == property_list[m]){ property_requested=true; break;}
            }
            if(property_requested){
              xtal_finder.structure_containers.back().properties.push_back(properties_response[i][l].second);
              xtal_finder.structure_containers.back().properties_names = property_list;
              xtal_finder.structure_containers.back().properties_units = property_units;
            }
          }
        }
        else {
          message << "More structures loaded than anticipated for auid=" << auids[i] << " (# structures=" << entry.vstr.size() << ").";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
        }
    }

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
    // write results to files //DX20201229 - consolidated into functions
    xtal_finder.writeComparisonOutputFile(ss_out, directory, "TEXT", "compare_database_entries", !structure_comparison);
    xtal_finder.writeComparisonOutputFile(ss_json, directory, "JSON", "compare_database_entries", !structure_comparison);

    return oss.str();

  }
}
//DX - COMPARE DATABASE ENTRIES - END

// ***************************************************************************
// XtalFinderCalculator::compareStructuresFromStructureList() //DX20201201
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareStructuresFromStructureList(
    const vector<string>& filenames,
    vector<string>& magmoms_for_systems,
    uint num_proc,
    bool same_species,
    const aurostd::xoption& comparison_options){

  // ---------------------------------------------------------------------------
  // directory to write results
  string directory = aurostd::getPWD();

  // ---------------------------------------------------------------------------
  // load structures appended to command
  loadStructuresFromStructureList(
      filenames,
      magmoms_for_systems,
      same_species);

  // ---------------------------------------------------------------------------
  // compare structures returns vector<StructureProtoype> of unique/duplicate info
  return compareMultipleStructures(num_proc, same_species, directory, comparison_options);

}

// ***************************************************************************
// XtalFinderCalculator::compareStructuresFromDirectory() //DX20201201
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareStructuresFromDirectory(
    const string& directory,
    vector<string>& magmoms_for_systems,
    uint num_proc,
    bool same_species,
    const aurostd::xoption& comparison_options){

  // ---------------------------------------------------------------------------
  // load structures in directory
  loadStructuresFromDirectory(
      directory,
      magmoms_for_systems,
      same_species);

  // ---------------------------------------------------------------------------
  // compare structures returns vector<StructureProtoype> of unique/duplicate info
  return compareMultipleStructures(num_proc, same_species, directory, comparison_options);

}

// ***************************************************************************
// XtalFinderCalculator::compareStructuresFromFile() //DX20201201
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareStructuresFromFile(
    const string& filename,
    vector<string>& magmoms_for_systems,
    uint num_proc,
    bool same_species,
    const aurostd::xoption& comparison_options){

  // ---------------------------------------------------------------------------
  // directory to write results
  string directory = aurostd::getPWD();

  // ---------------------------------------------------------------------------
  // load structures in file
  loadStructuresFromFile(filename,
      magmoms_for_systems,
      same_species);

  // ---------------------------------------------------------------------------
  // compare structures returns vector<StructureProtoype> of unique/duplicate info
  return compareMultipleStructures(num_proc, same_species, directory, comparison_options);

}

// ***************************************************************************
// compare::getUniqueEntries() //DX20201111
// ***************************************************************************
// Get structurally unique aflowlib entries. Helper function to GFA-code.
// If needed, this function can be extended to read in some of the
// structural properties from the database (i.e., symmetry) and suppress the
// on-the-fly analysis (potential speed-up).
namespace compare {
  vector<aflowlib::_aflowlib_entry> getUniqueEntries(vector<aflowlib::_aflowlib_entry>& entries, uint num_proc, bool same_species, bool scale_volume, bool optimize_match){
    ostringstream oss;
    ofstream FileMESSAGE;
    return getUniqueEntries(entries, oss, FileMESSAGE, num_proc, same_species, scale_volume, optimize_match);
  }
}

namespace compare {
  vector<aflowlib::_aflowlib_entry> getUniqueEntries(vector<aflowlib::_aflowlib_entry>& entries,
      ostream& oss,
      ofstream& FileMESSAGE,
      uint num_proc,
      bool same_species,
      bool scale_volume,
      bool optimize_match){

    // ---------------------------------------------------------------------------
    // instantiate XtalFinder calculator
    XtalFinderCalculator xtal_finder(
        DEFAULT_XTALFINDER_MISFIT_MATCH,
        DEFAULT_XTALFINDER_MISFIT_FAMILY,
        FileMESSAGE,
        num_proc,
        oss);

    // ---------------------------------------------------------------------------
    // load structures from aflowlib entries
    vector<string> magmoms_for_systems; //DX20201111 - not included for now
    xtal_finder.loadStructuresFromAflowlibEntries(
        entries,
        magmoms_for_systems,
        same_species);

    // ---------------------------------------------------------------------------
    // directory to write results
    string directory = aurostd::getPWD();

    // ---------------------------------------------------------------------------
    // default comparison options
    aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions();
    comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",scale_volume);
    comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH",optimize_match);

    // ---------------------------------------------------------------------------
    // compare structures returns vector<StructureProtoype> of unique/duplicate info
    vector<StructurePrototype> grouped_structures = xtal_finder.compareMultipleStructures(
        xtal_finder.num_proc,
        same_species,
        directory,
        comparison_options);

    // ---------------------------------------------------------------------------
    // filter out duplicate aflowlib entries
    vector<aflowlib::_aflowlib_entry> entries_unique;
    for(uint i=0;i<grouped_structures.size();i++){
      for(uint j=0;j<entries.size();j++){
        if(grouped_structures[i].structure_representative_struct->name==entries[j].auid){
          entries_unique.push_back(entries[j]);
          break;
        }
      }
    }

    return entries_unique;
  }
}

// ***************************************************************************
// compare::compareMultipleStructures() //DX20201201
// ***************************************************************************
vector<StructurePrototype> XtalFinderCalculator::compareMultipleStructures(
    uint num_proc, bool same_species,
    const string& directory){

  // ---------------------------------------------------------------------------
  // create xoptions to contain all comparison options
  aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions();
  if(!same_species){
    comparison_options.flag("COMPARISON_OPTIONS::REMOVE_DUPLICATE_COMPOUNDS",TRUE);
  }

  return compareMultipleStructures(num_proc, same_species, directory, comparison_options);

}

vector<StructurePrototype> XtalFinderCalculator::compareMultipleStructures(
    uint num_proc,
    bool same_species,
    const string& directory,
    const aurostd::xoption& comparison_options){

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  string function_name = XPID + "XtalFinderCalculator::compareMultipleStructures():";
  stringstream message;
  bool quiet = false;

  if(LDEBUG){cerr << function_name << " BEGIN" << endl;}  //CO20200508

  message << "Total number of structures to compare: " << structure_containers.size();
  pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

  // ---------------------------------------------------------------------------
  // convert to structures to certain representations //DX20201006
  // conversion type(s) is indicated in the comparison_options flag
  if(comparison_options.flag("COMPARISON_OPTIONS::PRIMITIVIZE") ||
      comparison_options.flag("COMPARISON_OPTIONS::MINKOWSKI") ||
      comparison_options.flag("COMPARISON_OPTIONS::NIGGLI")){
    convertStructures(comparison_options,num_proc);
  }

  // ---------------------------------------------------------------------------
  // calculate symmetries of structures
  // if already calculated, do not recalculate
  bool all_symmetries_calculated = true;
  for(uint i=0;i<structure_containers.size();i++){ all_symmetries_calculated = (all_symmetries_calculated&&isSymmetryCalculated(structure_containers[i])); } //DX20200810

  if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY") && !all_symmetries_calculated){
    calculateSymmetries(num_proc);
  }
  else if(!all_symmetries_calculated){
    setSymmetryPlaceholders();
  }

  // ---------------------------------------------------------------------------
  // calculate LFA environments of database entries
  // if already calculated, do not recalculate
  bool all_environments_calculated = true;
  for(uint i=0;i<structure_containers.size();i++){ all_environments_calculated = (all_environments_calculated&&isLFAEnvironmentCalculated(structure_containers[i])); } //DX20200810
  if(!comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS") && !all_environments_calculated){
    calculateLFAEnvironments(num_proc);
  }

  // ---------------------------------------------------------------------------
  // get nearest neighbor info, can perhaps only calculate this if necessary
  // (i.e., if we know we will perform a comparison)
  // if already calculated, do not recalculate //DX20201201
  bool all_neighbors_calculated = true;
  for(uint i=0;i<structure_containers.size();i++){ all_neighbors_calculated = (all_neighbors_calculated&&(structure_containers[i].nearest_neighbor_distances.size()!=0)); } //DX20200925 - gcc-10 warnings
  if(!all_neighbors_calculated){
    getNearestNeighbors(num_proc);
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
        remove_duplicates_options);

    message << "Duplicate materials removed.";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    // include duplicate compounds count in object
    for(uint i=0;i<unique_compounds.size();i++){
      unique_compounds[i].structure_representative_struct->number_compounds_matching_structure = unique_compounds[i].numberOfComparisons();
    }

    // ---------------------------------------------------------------------------
    // write results to files //DX20201229 - consolidated into functions
    stringstream ss_json_remove_duplicates, ss_out_remove_duplicates;

    if(XHOST.vflag_control.flag("PRINT_MODE::JSON")){
      printResults(ss_json_remove_duplicates, tmp_same_species, unique_compounds, "json");
      writeComparisonOutputFile(ss_json_remove_duplicates, directory, "JSON", "duplicate_compounds", true);
    }
    else if(XHOST.vflag_control.flag("PRINT_MODE::TXT")){
      printResults(ss_out_remove_duplicates, tmp_same_species, unique_compounds, "txt");
      writeComparisonOutputFile(ss_out_remove_duplicates, directory, "TEXT", "duplicate_compounds", true);
    }
    else{
      printResults(ss_json_remove_duplicates, tmp_same_species, unique_compounds, "json");
      writeComparisonOutputFile(ss_json_remove_duplicates, directory, "JSON", "duplicate_compounds", true);
      printResults(ss_out_remove_duplicates, tmp_same_species, unique_compounds, "txt");
      writeComparisonOutputFile(ss_out_remove_duplicates, directory, "TEXT", "duplicate_compounds", true);
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
  vector<StructurePrototype> comparison_schemes = groupStructurePrototypes(
      same_species,
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_SYMMETRY"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_WYCKOFF"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANALYSIS"),
      comparison_options.flag("COMPARISON_OPTIONS::IGNORE_ENVIRONMENT_ANGLES"), //DX20200320 - added environment angles
      false); //DX20200103 - condensed booleans to xoptions

  // ---------------------------------------------------------------------------
  // if ICSD comparison, make structure with minimum ICSD number the representative structure
  if(comparison_options.flag("COMPARISON_OPTIONS::ICSD_COMPARISON")){
    representativePrototypeForICSDRunsNEW(comparison_schemes);
  }

  // ---------------------------------------------------------------------------
  // compare structures
  vector<StructurePrototype> final_prototypes = runComparisonScheme(comparison_schemes, same_species, num_proc, comparison_options, quiet); //DX20200103 - condensed booleans to xoptions

  // ---------------------------------------------------------------------------
  // return if there are no similar structures
  if(final_prototypes.size()==0){ return final_prototypes; }
  comparison_schemes.clear();

  // ---------------------------------------------------------------------------
  // combine prototypes regardless of having different space groups
  // BETA TESTING (perhaps we shouldn't do this) - compare::checkPrototypes(num_proc,same_species,final_prototypes);

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
          if(!compare::generateStructure(
                final_prototypes[i].structure_representative_struct->name,
                final_prototypes[i].structure_representative_struct->source,
                final_prototypes[i].structure_representative_struct->relaxation_step,
                final_prototypes[i].structure_representative_struct->structure,*p_oss)){ //DX20200429
            message << "Could not generate structure (" << final_prototypes[i].structure_representative_struct->name << ").";
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function_name,message,_RUNTIME_ERROR_);
          }
        }
        XtalFinderCalculator xtal_finder_permutations;
        vector<StructurePrototype> final_permutations = xtal_finder_permutations.compareAtomDecorations(final_prototypes[i],num_proc,comparison_options.flag("COMPARISON_OPTIONS::OPTIMIZE_MATCH"));
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
        XtalFinderCalculator xtal_finder_permutations;
        vector<string> unique_permutations = xtal_finder_permutations.getSpeciesPermutedStrings(final_prototypes[i].stoichiometry); //DX20191125
        // store permutation results in main StructurePrototype object
        for(uint j=0;j<unique_permutations.size();j++){
          vector<string> permutations_tmp; permutations_tmp.push_back(unique_permutations[j]);
          final_prototypes[i].atom_decorations_equivalent.push_back(permutations_tmp);
        }
      }
    }
  }

  if(comparison_options.flag("COMPARISON_OPTIONS::ADD_AFLOW_PROTOTYPE_DESIGNATION")){
    // SEPARATE FUNCTION????
    message << "Determining the AFLOW standard designation.";
    pflow::logger(_AFLOW_FILE_NAME_, function_name, message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
    // ---------------------------------------------------------------------------
    // split task into threads
    uint number_of_structures = final_prototypes.size();
    uint number_of_threads = aurostd::min(num_proc,number_of_structures); // cannot have more threads than structures
    vector<vector<int> > thread_distribution = getThreadDistribution(number_of_structures, number_of_threads); //DX20191107

    // ---------------------------------------------------------------------------
    // [THREADED] determine AFLOW standard designation
    vector<std::thread*> threads;
    for(uint n=0; n<number_of_threads; n++){
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
    stringstream ss_json, ss_out;
    printResults(ss_json, same_species, final_prototypes, "json");
    writeComparisonOutputFile(ss_json, directory, "JSON", "compare_input", same_species);
    printResults(ss_out, same_species, final_prototypes, "text");
    writeComparisonOutputFile(ss_out, directory, "TEXT", "compare_input", same_species);
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
      vector<StructurePrototype> matching_protos = xtal_finder_protos.compare2prototypes(final_prototypes[i].structure_representative_struct->structure, vpflow_protos);
      for(uint j=0;j<matching_protos[0].structures_duplicate_struct.size();j++){
        final_prototypes[i].matching_aflow_prototypes.push_back(matching_protos[0].structures_duplicate_struct[j]->name);
      }
    }
  }

  return final_prototypes;
}

// ***************************************************************************
// XtalFinderCalculator::groupSimilarXstructures() //DX20201229
// ***************************************************************************
vector<vector<uint> > XtalFinderCalculator::groupSimilarXstructures(
    const vector<xstructure>& vxstrs,
    bool same_species,
    bool scale_volume) {

  // Compares a vector of xstructures and groups the indices of the
  // xstructures by structural similarity

  uint nstrs= vxstrs.size();

  stringstream structure_name;
  string source = "input";
  uint relaxation_step = 0; // unknowable from input
  string directory = aurostd::getPWD();

  // add structures to containers
  for(uint i=0;i<nstrs;i++){
    structure_name.str("");
    structure_name << i; // name by index
    addStructure2container(vxstrs[i], structure_name.str(), source, relaxation_step, same_species);
  }

  // ---------------------------------------------------------------------------
  // set comparison options
  aurostd::xoption comparison_options = compare::loadDefaultComparisonOptions();
  comparison_options.flag("COMPARISON_OPTIONS::SCALE_VOLUME",scale_volume);

  // ---------------------------------------------------------------------------
  // compare structures returns vector<StructureProtoype> of unique/duplicate info
  vector<StructurePrototype> grouped_structures = compareMultipleStructures(
      num_proc,
      same_species,
      directory,
      comparison_options); //DX20200103 - condensed booleans to xoptions

  // ---------------------------------------------------------------------------
  // group indices based on structural similarity
  vector<vector<uint> > grouped_indices;
  vector<uint> structure_indices_equivalent;
  for(uint i=0;i<grouped_structures.size();i++){
    structure_indices_equivalent.clear();
    structure_indices_equivalent.push_back(aurostd::string2utype<uint>(grouped_structures[i].structure_representative_struct->name));
    for(uint j=0;j<grouped_structures[i].structures_duplicate_struct.size();j++){
      structure_indices_equivalent.push_back(aurostd::string2utype<uint>(grouped_structures[i].structures_duplicate_struct[j]->name));
    }
    grouped_indices.push_back(structure_indices_equivalent);
  }

  return grouped_indices;
}

// ***************************************************************************
// compare::structuresMatch() //DX20201228 - renamed
// ***************************************************************************
namespace compare {
  bool structuresMatch(const xstructure& xstr1,
      const xstructure& xstr2,
      bool same_species,
      uint num_proc) {

    double final_misfit = AUROSTD_MAX_DOUBLE;
    bool scale_volume=true; //default is true
    bool optimize_match=false; //default is false
    return aflowCompareStructure(xstr1, xstr2, same_species, scale_volume, optimize_match, final_misfit, num_proc); //DX20191122 - move ostream to end
  }
}

namespace compare {
  bool structuresMatch(const xstructure& xstr1,
      const xstructure& xstr2,
      bool same_species,
      bool scale_volume,
      bool optimize_match,
      uint num_proc) {

    double final_misfit = AUROSTD_MAX_DOUBLE;
    return aflowCompareStructure(xstr1, xstr2, same_species, scale_volume, optimize_match, final_misfit, num_proc); //DX20191122 - move ostream to end and add default
  }
}


// ***************************************************************************
// compare::getMisfitBetweenStructures() //DX20201228 - renamed
// ***************************************************************************
namespace compare {
  double getMisfitBetweenStructures(const xstructure& xstr1,
      const xstructure& xstr2,
      bool same_species,
      uint num_proc) {

    double final_misfit=AUROSTD_MAX_DOUBLE;
    bool scale_volume=true; //default is true
    bool optimize_match=false; //default is false
    aflowCompareStructure(xstr1, xstr2, same_species, scale_volume, optimize_match, final_misfit, num_proc); //DX20191122 - move ostream to end
    return final_misfit;
  }
}

// ***************************************************************************
// compare::getTransformationBetweenStructures() //DX20201229
// ***************************************************************************
namespace compare {
  structure_misfit getTransformationBetweenStructures(const xstructure& xstr1, const xstructure& xstr2, bool same_species, uint num_proc) { //DX20191108 - remove const & from bools
    double final_misfit=AUROSTD_MAX_DOUBLE;
    bool scale_volume=true; //default is true
    bool optimize_match=false; //default is false
    structure_misfit final_misfit_info = compare::initialize_misfit_struct();
    aflowCompareStructure(xstr1, xstr2, same_species, scale_volume, optimize_match, final_misfit, final_misfit_info, num_proc); //DX20191122 - move ostream to end
    return final_misfit_info;
  }
}

// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const xstructure& xstr1,
      const xstructure& xstr2,
      bool same_species,
      bool scale_volume,
      bool optimize_match,
      double& final_misfit,
      uint num_proc) {

    structure_misfit final_misfit_info = compare::initialize_misfit_struct(); //DX20191218

    return aflowCompareStructure(xstr1, xstr2, same_species, scale_volume, optimize_match, final_misfit, final_misfit_info, num_proc); //DX20191122 - move ostream to end
  }
}

// ***************************************************************************
// compare::aflowCompareStructure - MAIN FUNCTION
// ***************************************************************************
namespace compare {
  bool aflowCompareStructure(const xstructure& xstr1,
      const xstructure& xstr2,
      bool same_species,
      bool scale_volume,
      bool optimize_match,
      double& final_misfit,
      structure_misfit& final_misfit_info,
      uint num_proc) { //DX20191108 - remove const & from bools //DX20191122 - move ostream to end and add default

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

    XtalFinderCalculator xtal_finder(num_proc);
    xtal_finder.compareStructures(str_rep,str_matched,final_misfit_info,same_species,scale_volume,optimize_match);

    final_misfit = final_misfit_info.misfit;
    return(final_misfit<=xtal_finder.misfit_match);

  }
}

// ***************************************************************************
// XtalFinderCalculator::compareStructures() - MAIN COMPARISON FUNCTION
// ***************************************************************************
void XtalFinderCalculator::compareStructures(
    _structure_representative& str_rep,
    _structure_representative& str_matched,
    structure_misfit& match_info,
    bool same_species,
    bool scale_volume,
    bool optimize_match) {

  // This is the main comparison function that compares two crystal structures
  // and determines their similarity level based on the idea discussed
  // in H. Burzlaff's paper (Acta Cryst., A53, 217-224 (1997)).

  bool LDEBUG=(FALSE || XHOST.DEBUG || _DEBUG_COMPARE_);
  string function_name = XPID + "XtalFinderCalculator::compareStructures():";

  // ---------------------------------------------------------------------------
  // determine minimum interatomic distances of structures (resolution of atoms) //DX20200623
  // //DX20200715 - may need to rescale this if the structures are being scaled later....
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
    if(LDEBUG) {cerr << function_name << " Searching for new representation of test structure ..."<<endl;}
    latticeSearch(str_rep,str_matched,match_info,same_species,optimize_match,scale_volume,num_proc); //DX20190530 //DX20200422 - scale_volume added
  }
}

// AFLOW-XtalFinder (compare crystal structures)
// Written by David Hicks (david.hicks@duke.edu)
// Contributors: Carlo De Santo
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow DAVID HICKS - Duke University 2014-2020                 *
// *                                                                         *
// ***************************************************************************
