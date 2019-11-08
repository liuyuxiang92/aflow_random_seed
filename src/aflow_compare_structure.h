// ***************************************************************************
// AFLOW_COMPARE_STRUCTURE
// ***************************************************************************

#include<iostream>
#include "aflow_pflow.h"
#include "aflow.h"
#include "math.h"
#include<vector>
#include<string>
//#include <stdatomic.h> // Needed to communicate between threads
//#include<atomic>

//DX 20190314 [OBSOLETE] using aurostd::isequal;
//DX 20190314 [OBSOLETE] using aurostd::string2utype;
//DX 20190314 [OBSOLETE] using aurostd::FileExist;

#ifndef __AFLOW_COMPARE_STRUCTURE_H_
#define __AFLOW_COMPARE_STRUCTURE_H_

// Create the version of GCC, we will uset it for multithread parts of code,
// to check if the current compiling version of gcc is able to compile the
// std::thead features
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

#define JSON_MODE 0

//// ===== GroupedWyckoffPosition Class ===== //
//class GroupedWyckoffPosition{
//  public:
//    GroupedWyckoffPosition();
//    ~GroupedWyckoffPosition();
//    friend ostream& operator<<(ostream& oss, const GroupedWyckoffPosition& GroupedWyckoffPosition);
//    const GroupedWyckoffPosition& operator=(const GroupedWyckoffPosition& b);
//    GroupedWyckoffPosition(const GroupedWyckoffPosition& b);
//    uint type;
//    string element;
//    vector<string> site_symmetries;
//    vector<uint> multiplicities;
//  private:
//    void free();
//    void copy(const GroupedWyckoffPosition& b);
//};

// ===== AtomEnvironment Class ===== //
class AtomEnvironment{
  public:
    AtomEnvironment();                                                                      // constructor operator
    ~AtomEnvironment();                                                                     // destructor operator
    friend ostream& operator<<(ostream& oss, const AtomEnvironment& AtomEnvironment);       // stringstream operator (printing)
    const AtomEnvironment& operator=(const AtomEnvironment& b);                             // assignment operator
    AtomEnvironment(const AtomEnvironment& b);                                              // copy constructor
    string center_element;                                                                  // species/element at center of environment                                                                   
    uint center_type;                                                                       // type (uint) at center of environment
    vector<string> neighbor_elements;                                                       // species/element of atoms neighboring center atom
    vector<uint> neighbor_types;                                                            // types (uint) of atoms neighboring center atom
    vector<double> neighbor_distances;                                                      // distances to atoms neighboring atoms (typically put in a bin with small tolerance threshold)                                             
    vector<uint> neighbor_frequencies;                                                      // frequency of neighboring distance                                              
    vector<vector<xvector<double> > > neighbor_coordinates;                                 // coordinates of atoms neighboring atoms (center is assumed to be zero,i.e. coord=neighbor-origin)
    //functions
  private:
    void free();                                                                            // free operator
    void copy(const AtomEnvironment& b);                                                    // copy constructor
};


// ===== StructurePrototype Class ===== //
class StructurePrototype{
  public:
    StructurePrototype();                                                                   // constructor operator
    ~StructurePrototype();                                                                  // destructor operator
    friend ostream& operator<<(ostream& oss, const StructurePrototype& StructurePrototype); // stringstream operator (printing)
    const StructurePrototype& operator=(const StructurePrototype& b);                       // assignment operator
    StructurePrototype(const StructurePrototype& b);                                        // copy constructor
    int iomode;                                                                             // mode for printing
    string structure_representative_name;                                                   // name of representative structure
    string structure_representative_compound;                                               // compound name of representative structure (w/reduced stoichometry), e.g., Ag1Br2
    xstructure structure_representative;                                                    // xstructure of representative structure
    bool structure_representative_generated;                                                // boolean indicating if xstructure is generated
    string structure_representative_from;                                                   // string indicating where structure came from, i.e., input, auid, aflow protos, etc.
    uint number_compounds_matching_representative;                                          // number of compounds that match with the representative structure
    int number_of_types;                                                                    // number of types in prototype
    vector<string> elements;                                                                // list of elements exhibiting in this protoype (from representative and duplicate structures)
    vector<uint> stoichiometry;                                                             // reduced stoichiometry of prototype
    uint number_of_atoms;                                                                   // number of atoms in the prototype (from the representative structure; not necessarily reduced)
    vector<string> unique_permutations;                                                     // unique permutations of prototype
    string Pearson;                                                                         // Pearson symbol of prototype
    uint space_group;                                                                       // space group number of prototype
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;                               // Wyckoff positions grouped by site type
    vector<string> wyckoff_site_symmetry;                                                   // vector of Wyckoff site symmetries of prototype
    vector<int> wyckoff_multiplicity;                                                       // vector of Wyckoff multiplicities of prototype
    vector<int> wyckoff_letter;                                                             // vector of Wyckoff letters of prototype
    vector<AtomEnvironment> environments_LFA;                                               // vector of LFA atom environments
    string aflow_label;                                                                     // AFLOW label designation  
    vector<string> aflow_parameter_list;                                                    // vector of strings corresponding to AFLOW parameter variables 
    vector<double> aflow_parameter_values;                                                  // vector of doubles corresponding to AFLOW parameter values
    vector<string> matching_aflow_prototypes;                                               // vector of strings indicating matching AFLOW prototype labels
    vector<string> structures_duplicate_names;                                              // vector of names for duplicate structures
    vector<string> structures_duplicate_compounds;                                          // vector of compounds names for duplicate structures (w/reduced stoichometry), e.g., Ag1Br2
    vector<xstructure> structures_duplicate;                                                // vector of xstructures for duplicate structures
    vector<bool> structures_duplicate_generated;                                            // vector of booleans indicating if xstructure is generated for duplicate structures
    vector<string> structures_duplicate_from;                                               // vector of strings indicating where structure came from, i.e., input, auid, aflow protos, etc. 
    vector<vector<GroupedWyckoffPosition> > structures_duplicate_grouped_Wyckoff_positions; // Wyckoff positions grouped by site type
    vector<uint> number_compounds_matching_duplicate;                                       // vector of number of compounds that match with the duplicate structures
    vector<string> duplicate_comparison_logs; //DX 20190506                                 // vector of comparison logs for duplicate structures, CAREFUL: very long string //DX 20190506
    vector<string> structures_family_names;                                                 // vector of names for structures within the same family
    vector<xstructure> structures_family;                                                   // vector of xstructures for structures within the same family
    vector<bool> structures_family_generated;                                               // vector of booleans indicating if xstructure is generated for structures within the same family
    vector<string> structures_family_from;                                                  // vector of strings indicating where structure came from, i.e., input, auid, aflow protos, etc.
    vector<vector<GroupedWyckoffPosition> > structures_family_grouped_Wyckoff_positions;    // Wyckoff positions grouped by site type
    vector<uint> number_compounds_matching_family;                                          // vector of number of compounds that match with the same family structures
    vector<string> family_comparison_logs; //DX 20190506                                    // vector of comparison logs for same family structures, CAREFUL: very long string //DX 20190506
    vector<double> misfits_duplicate;                                                       // vector of misfit values when comparing the duplicate to the representative structure
    vector<double> misfits_family;                                                          // vector of misfit values when comparing the same family structure to the representative structure  
    vector<string> property_names;                                                          // vector of property names (if using AFLUX)
    vector<string> property_units;                                                          // vector of property units (if using AFLUX)
    vector<string> properties_structure_representative;                                     // vector of property values for the representative structure (if using AFLUX)
    vector<vector<string> > properties_structures_duplicate;                                // vector of property values for the duplicate structures (if using AFLUX) 
    vector<vector<string> > properties_structures_family; //DX 20190506                     // vector of property values for the family structures (if using AFLUX) //DX 20190425 
    // functions
    uint numberOfDuplicates() const; //DX 20190506                                          // return the number of duplicate structures for this prototype (i.e., checks misfit value)
    uint numberOfComparisons(); //DX 20181221                                               // return the number of comparisons for this prototype 
    bool isSymmetryCalculated(); //DX 20190228
    bool isLFAEnvironmentCalculated(); //DX 20191105
    bool calculateSymmetry(); //DX 20190118                                                 // calculate space group symmetry and populate symmetry attributes for StructurePrototype
    bool addStructurePrototypeAsDuplicate(StructurePrototype& b);                           // combine structure with another StructurePrototype object as a potential duplicate
    void putDuplicateAsFamily(uint index, bool keep_generated=false);                       // make duplicate structure a same family structure in the same StructurePrototypeObject //DX 20190814 
    bool copyPrototypeInformation(StructurePrototype& b);                                   // copy prototype information from one StructurePrototype object to another
    bool putDuplicateAsRepresentative(StructurePrototype& b, uint& index);                  // make duplicate structure a representative structure in a new StructurePrototype object
    //bool putRepresentativeAsDuplicate(StructurePrototype& b);                               // make representative structure a duplicate structure in a new StructurePrototype object //DX 20190730
    bool copyDuplicate(StructurePrototype& b, uint& index, bool copy_misfit=false);         // copy duplicate info to another StructurePrototype object
    bool removeNonDuplicate(uint& index);                                                   // remove non-duplicate structures 
    bool removeDuplicates(bool remove_duplicate_count);                                     // remove duplicate structure information
  private:
    void free();                                                                            // free operator
    void copy(const StructurePrototype& b);                                                 // copy constructor
};


namespace compare{

  // ===== Main functions ===== //
  vector<StructurePrototype> compareStructuresFromStructureList(vector<string>& filenames, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species,
      bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds,
      bool calculate_unique_permutations, bool add_matching_aflow_protos, bool get_aflow_prototype_designation, bool ICSD_comparison, bool store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
  vector<StructurePrototype> compareStructuresFromDirectory(string& directory, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species,
      bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds, 
      bool calculate_unique_permutations, bool add_matching_aflow_protos, bool get_aflow_prototype_designation, bool ICSD_comparison, bool store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
  vector<StructurePrototype> compareStructuresFromFile(string& filename, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species,
      bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds, 
      bool calculate_unique_permutations, bool add_matching_aflow_protos, bool get_aflow_prototype_designation, bool ICSD_comparison, bool store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190506 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& all_structures, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, string& directory); //DX 20190319 - added FileMESSAGE
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& all_structures, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, string& directory,
      bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds, 
      bool calculate_unique_permutations, bool add_matching_aflow_protos, bool get_aflow_prototype_designation, bool ICSD_comparison, bool store_comparison_logs); //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched //DX 20190724 - added add_matching_aflow_protos, get_aflow_prototype_designation, calculate_unique_permutations, ignore_environment //DX 20190822 - add log bool
  bool aflowCompareStructure(const uint& num_proc, const xstructure& xstr1, const xstructure& xstr2, 
      bool same_species, bool scale_volume, bool fast_match, ostream& oss, double& final_misfit); //Main function //DX 20191108 - remove const & from bools
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species); //Overloaded, returns true (match), false (no match) //DX 20191108 - remove const & from bools
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool fast_match);  //DX 20191108 - remove const & from bools
  double aflowCompareStructureMisfit(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool fast_match); //Overloaded, returns misfit value //DX 20191108 - remove const & from bools

  // permutaion comparisons
  vector<string> getUniquePermutations(xstructure& xstr);
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc); // add number of threads
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc, bool& optimize_match); // add number of threads and optimize matches
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc, bool& optimize_match, bool& print_misfit, ostream& oss, ofstream& FileMESSAGE); // full function //DX 20190319 - added FileMESSAGE


  //  string CompareStructures(aurostd::xoption& vpflow);
  //  string CompareStructureDirectory(aurostd::xoption& vpflow);

  // ===== Compare Directory Functions ===== //
  vector<StructurePrototype> loadStructuresFromStructureList(vector<string>& filenames, vector<string>& magmoms_for_systems, bool same_species, ofstream& FileMESSAGE); //DX 20190424 //DX 20190801 - added vector<string>& magmoms_for_systems
  vector<StructurePrototype> loadStructuresFromDirectory(string& directory, vector<string>& magmoms_for_systems, bool same_species, ofstream& FileMESSAGE); //DX 20190319 - added FileMESSAGE //DX 20190801 - added vector<string>& magmoms_for_systems,
  vector<StructurePrototype> loadStructuresFromFile(string& directory, vector<string>& magmoms_for_systems, bool same_species, ofstream& FileMESSAGE); //DX 20190319 - added FileMESSAGE //DX 20190801 - added vector<string>& magmoms_for_systems,
  void generateStructures(vector<StructurePrototype>& structures, ostream& oss); //DX 20191105
  void generateStructuresInRange(vector<StructurePrototype>& structures, ostream& oss, uint start_index, uint end_index); //DX 20191105
  bool generateStructure(string& structure_name, string& structure_from, xstructure& structure, ostream& oss);
  void removeNonGeneratedStructures(vector<StructurePrototype>& structures); //DX 20191105
  vector<uint> getStoichiometry(string& compositions, const bool& same_species);
  bool addAFLOWPrototypes2StructurePrototypeVector(vector<StructurePrototype>& all_structures, vector<string>& vlabel);
  string getCompoundName(xstructure& xstr, bool remove_ones=false);
  string getCompoundName(vector<string>& elements, vector<uint>& stoichiometry, bool remove_ones=false);
  vector<uint> getStoichiometry(const xstructure& xstr, const bool& same_species);
  vector<string> getElements(xstructure& xstr);
  vector<uint> gcdStoich(const vector<uint>& numbers); //DX 20181009
  vector<uint> gcdStoich(const deque<int>& numbers);
  //DX 20191108 [OBSOLETE - switching to getThreadDistribution] bool prepareSymmetryThreads(vector<xstructure>& vxstrs, uint& num_proc,
  //DX 20191108 [OBSOLETE - switching to getThreadDistribution]     vector<uint>& start_indices, vector<uint>& end_indices);
  //DX 20191108 [OBSOLETE - switching to getThreadDistribution] bool prepareSymmetryThreads(uint& number_of_structures, uint& num_proc,
  //DX 20191108 [OBSOLETE - switching to getThreadDistribution]     vector<uint>& start_indices, vector<uint>& end_indices);
  //DX 20191108 [OBSOLETE - switching to getThreadDistribution] bool splitTaskIntoThreads(uint& number_of_tasks, uint& num_proc,              
  //DX 20191108 [OBSOLETE - switching to getThreadDistribution]     vector<uint>& start_indices, vector<uint>& end_indices); //DX 20190530 - renamed
  void calculateSymmetries(vector<xstructure>& vxstrs, vector<string>& vpearsons, vector<uint>& vsgroups, 
      vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions, uint& num_proc);
  void calculateSpaceGroups(vector<StructurePrototype>& structures, uint start_index, uint end_index); //DX 20191107 - remove & from uint
  void calculateLFAEnvironmentsInSetRange(vector<StructurePrototype>& structures, uint start_index, uint end_index);  //DX 20191105
  void calculateLFAEnvironments(vector<StructurePrototype>& structures, uint num_proc);  //DX 20191105
  void calculateSymmetries(vector<StructurePrototype>& structures, uint& num_proc);  //DX 20190118
  void calculateSymmetry(xstructure& xstr, vector<string>& vpearsons, vector<uint>& vsgroups,
      vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions);
  bool groupWyckoffPositions(xstructure& xstr, vector<GroupedWyckoffPosition>& grouped_positions);
  bool groupWyckoffPositionsFromGroupedString(uint& space_group_number, uint& setting, vector<vector<string> >& grouped_Wyckoff_string, vector<GroupedWyckoffPosition>& grouped_positions);
  string printWyckoffString(const vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize=false);
  vector<GroupedWyckoffPosition> sortSiteSymmetryOfGroupedWyckoffPositions(const vector<GroupedWyckoffPosition>& grouped_Wyckoffs); //DX 20190219
  bool matchableWyckoffPositions(const vector<GroupedWyckoffPosition>& temp_grouped_Wyckoffs,
      const vector<GroupedWyckoffPosition>& representative_grouped_Wyckoffs,
      bool same_species);
  bool matchableWyckoffPositionSet(vector<vector<vector<string> > > grouped_possible_Wyckoff_letters,
      vector<vector<string> > grouped_Wyckoff_letters);

  vector<vector<string> > convertANRLWyckoffString2GroupedPositions(string label);
  vector<vector<string> > convertWyckoffString2GroupedPositions(string Wyckoff_letter_string);
  bool sameStoichiometry(const vector<uint>& stoich1, const vector<uint>& stoich2);
  bool matchableSpaceGroups(uint space_group_1, uint space_group_2);
  bool matchableEnantiomorphicSpaceGroups(uint space_group_1, uint space_group_2);
  bool filterPrototypes(uint& species_count, string& reduced_stoichiometry, uint& space_group,
      vector<vector<vector<string> > >& grouped_possible_Wyckoff_letters,
      vector<string>& prototype_labels, vector<uint>& species_counts,
      vector<uint>& space_groups);
  /*vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& comparison_schemes, 
    vector<xstructure>& vxstrs, const bool& same_species, 
    const vector< vector<string> >& vvelements,
    vector< vector<uint> >& vstoichs, vector<string>& vpearsons, 
    vector<uint>& vsgroups, 
    vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
    const string& directory, const vector<string>& vfiles, 
    vector<bool>& vstructures_generated,
    vector<string>& vstructures_from,
    const bool& ignore_symmetry, const bool& ignore_Wyckoff,
    const bool& structures_generated);
    void createStructurePrototypes(vector<StructurePrototype>& comparison_schemes, 
    vector<xstructure>& vxstrs, const bool& same_species, 
    const vector< vector<string> >& vvelements,
    vector< vector<uint> >& vstoichs, vector<string>& vpearsons, 
    vector<uint>& vsgroups, 
    vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
    vector<string>& property_names,
    vector<string>& property_units,
    vector<vector<string> >& property_values,
    const string& directory, const vector<string>& vfiles,
    vector<bool>& vstructures_generated,
    vector<string>& vstructures_from,
    const bool& ignore_symmetry, const bool& ignore_Wyckoff,
    const bool& structures_generated);
    */
  void createStructurePrototypes(vector<StructurePrototype>& comparison_schemes,
      vector<xstructure>& vxstrs, const bool& same_species,
      const vector< vector<string> >& vvelements,
      vector< vector<uint> >& vstoichs, vector<string>& vpearsons,
      vector<uint>& vsgroups,
      vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
      const string& directory, const vector<string>& vfiles,
      vector<bool>& vstructures_generated,
      vector<string>& vstructures_from,
      const bool& ignore_symmetry, const bool& ignore_Wyckoff);
  void createStructurePrototypes(vector<StructurePrototype>& comparison_schemes,
      vector<xstructure>& vxstrs, const bool& same_species,
      const vector< vector<string> >& vvelements,
      vector< vector<uint> >& vstoichs, vector<string>& vpearsons,
      vector<uint>& vsgroups,
      vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
      vector<string>& property_names,
      vector<string>& property_units,
      vector<vector<string> >& property_values,
      const string& directory, const vector<string>& vfiles,
      vector<bool>& vstructures_generated,
      vector<string>& vstructures_from,
      const bool& ignore_symmetry, const bool& ignore_Wyckoff);
  bool structuresCompatible(const StructurePrototype& structure1,
      const StructurePrototype& structure2, bool same_species, 
      bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment); //DX 20190730
  vector<StructurePrototype> groupStructurePrototypes(vector<StructurePrototype>& structures, 
      bool same_species, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool duplicates_removed); //DX 20190731 - remove const and & //DX 20190830 - added duplicates_removed

  vector<StructurePrototype> checkForBetterMatches(vector<StructurePrototype>& prototype_schemes, 
      ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool check_for_better_matches, bool same_species,
      bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, 
      bool ignore_environment, bool clean_unmatched, bool ICSD_comparison, bool quiet=false); //DX 20190730 
  vector<StructurePrototype> compareDuplicateCompounds(vector<StructurePrototype>& prototype_schemes, uint& num_proc, 
      bool& ICSD_comparison, ostringstream& oss);
  vector<StructurePrototype> createComparisonSchemeForDuplicateCompounds(StructurePrototype& prototype_scheme); 
  void removeDuplicateCompounds(vector<StructurePrototype>& final_prototypes, vector<StructurePrototype>& duplicate_compound_comparisons);

  bool representativePrototypeForICSDRuns(vector<StructurePrototype>& comparison_schemes);
  bool splitComparisonIntoThreads(vector<StructurePrototype>& comparison_schemes, uint& num_proc,
      vector<std::pair<uint,uint> >& start_indices,
      vector<std::pair<uint,uint> >& end_indices);
  void runComparisonThreads(vector<StructurePrototype>& comparison_schemes, 
      std::pair<uint,uint>& start_indices,
      std::pair<uint,uint>& end_indices,
      ostream& oss,
      bool same_species, 
      bool scale_volume, bool optimize_match,
      bool store_comparison_logs); //DX 20190822 - added comparison log bool
  //OLD vector<StructurePrototype> runComparisonScheme(uint& num_proc, vector<StructurePrototype>& comparison_schemes, const bool& same_species, 
  //OLD                                                const bool& scale_volume, const bool& optimize_match, const bool& single_comparison_round, const bool& structures_generated, const bool& ICSD_comparison, ostringstream& oss);
  vector<StructurePrototype> runComparisonScheme(uint num_proc, vector<StructurePrototype>& comparison_schemes, 
      bool same_species, bool check_other_grouping, bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, 
      bool ignore_environment, bool single_comparison_round, bool clean_unmatched, bool ICSD_comparison, bool store_comparison_logs, 
      ostream& oss, ofstream& FileMESSAGE, bool quiet=false); //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched //DX 20190731 - removed const and &, added ignore_symmetry/Wyckoff/environment //DX 20190822 - add log bool
  vector<std::pair<uint,uint> > calculateDivisors(const int& number);
  bool checkNumberOfGroupings(vector<StructurePrototype>& comparison_schemes, uint number);
  void createStructurePermutations(vector<StructurePrototype>& comparison_schemes, const vector<vector<string> >& name_order,
      vector<vector<GroupedWyckoffPosition> >& permutation_grouped_Wyckoff_positions,
      const vector<xstructure>& vxstrs, const bool& same_species);
  bool makeRepresentativeEvenPermutation(vector<StructurePrototype>& comparison_schemes, const vector<vector<string> >& name_order);  
  int numberMismatches(const vector<StructurePrototype> comparison_schemes);
  void appendStructurePrototypes(vector<StructurePrototype>& comparison_schemes, 
      vector<StructurePrototype>& final_prototypes,
      bool clean_unmatched, //DX 20190506
      bool quiet=false);
  void checkPrototypes(const uint& num_proc, const bool& same_species, vector<StructurePrototype>& final_prototypes);
  void printResults(ostream& ss_out, const bool& same_species, const vector<StructurePrototype>& final_prototypes, string mode="txt");
  void printStructureMappingResults(ostream& oss, 
      const xstructure& xstr_reference,
      const xstructure& xstr_transformed,
      const double misfit,
      const double lattice_deviation,
      const double coordinate_displacement,
      const double failure,
      const double magnetic_displacement,
      const double magnetic_failure,
      const vector<uint>& matching_indices_1, 
      const vector<uint>& matching_indices_2,
      const vector<double>& minimum_distances,
      bool magnetic_analysis,
      string mode="TEXT"); //DX 20190802

  bool SVD(xmatrix<double>& A);
  string findICSDName(string& name);
  string findMinimumICSDEntry(vector<string>& ICSD_entries);
  bool groupSameRatios(vector<int>& stoich, vector<int>& unique_stoich, vector<vector<int> >& type_index);
  //vector<vector<int> > generatePermutations(uint& num_elements, vector<int>& indices);
  vector<StructurePrototype> comparePermutations(StructurePrototype& structure, uint& num_proc, bool& optmize_match, ostream& oss, ofstream& FileMESSAGE); //DX 20190319 - added FileMESSAGE
  vector<StructurePrototype> generatePermutationStructures(StructurePrototype& structure);
  vector<string> generatePermutationString(vector<uint>& stoichiometry); //DX 20190508
  bool generatePermutations(uint& num_elements, vector<uint>& indices, vector<string>& names, vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, vector<vector<uint> >& permutations, vector<vector<string> >&name_order, vector<vector<GroupedWyckoffPosition> >& permutation_grouped_Wyckoff_positions);
  bool arePermutationsComparableViaStoichiometry(const xstructure& xstr); //DX 20190624 
  bool arePermutationsComparableViaStoichiometry(vector<uint>& stoichiometry, bool reduce_stoichiometry=false); //DX 20190624
  bool arePermutationsComparableViaSymmetry(vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions); //DX 20190624
  bool makePermutations(StructurePrototype& structure, vector<vector<string> >& name_order, vector<StructurePrototype>& permutation_structures);
  bool makePermutations(xstructure& xstr, vector<vector<string> >& name_order, vector<xstructure>& xstr_permutations);


  // ===== All other functions ===== //
  bool matchableSpecies(const xstructure& xstr1, const xstructure& xstr2, const bool& same_species);
  xstructure GetSuperCell3x3x3(const xstructure& aa, const xmatrix<double> &supercell);
  bool sameSpecies(const xstructure& x1, const xstructure& x2, const bool& display); 
  void rescaleStructure(xstructure& x1, xstructure& x2);
  void atomicNumberDensity(xstructure& xstr1, xstructure& xstr2);
  vector<string> fakeElements(const uint& number_of_species);
  void fakeAtomsName(xstructure& x1);
  void printParameters(xstructure& xstr, ostream& oss);
  string leastFrequentAtom(const xstructure& xstr);
  vector<string> getLeastFrequentAtomSpecies(const xstructure& xstr);
  bool sortBySecondPair(const std::pair<string,uint>& a, const std::pair<string,uint>& b);
  vector<string> sortSpeciesByFrequency(const xstructure& xstr);
  bool checkTolerance(xvector<double> d1, xmatrix<double> Q2);
  bool checkABCTolerance(xvector<double> d1, xvector<double> d2);
  bool checkAngleTolerance(xvector<double> d1, xvector<double> d2);
  void resetLatticeDimensions(const xmatrix<double>& lattice, double radius, xvector<int>& dims,
      vector<xvector<double> >& l1, vector<xvector<double> >& l2, 
      vector<xvector<double> >& l3, vector<int>& a_index, 
      vector<int>& b_index, vector<int>& c_index); //DX 20190705
  void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
      double& min_dist, uint& frequency, vector<xvector<double> >& coordinates); //DX 20191105
  void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
      double& min_dist, uint& frequency, vector<xvector<double> >& coordinates, double radius); //DX 20191105
  void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice, xvector<int>& dims,
      vector<xvector<double> >& l1, vector<xvector<double> >& l2, vector<xvector<double> >& l3, 
      vector<int>& a_index, vector<int>& b_index, vector<int>& c_index, 
      double& min_dist, uint& frequency, vector<xvector<double> >& coordinates,
      double radius); //DX 20191105
  void minimumCoordinationShell(const xstructure& xstr, uint center_index, 
      double& min_dist, uint& frequency, vector<xvector<double> >& coordinates); //DX 20191105
  void minimumCoordinationShell(const xstructure& xstr, uint center_index, 
      double& min_dist, uint& frequency, vector<xvector<double> >& coordinates, const string& type); //DX 20191105
  xvector<double> centroid_with_PBC(const xstructure& xstr);
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, const xmatrix<double>& lattice);
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, vector<double>& weights,
      const xmatrix<double>& lattice);
  //bool findMatch(const xstructure& xstr1, const xstructure& PROTO,vector<uint>& im1, vector<uint>& im2, vector<double>& min_dists, const int& type_match);
  bool findMatch(const deque<_atom>& xstr1_atoms, const deque<_atom>& PROTO_atoms,
      const xmatrix<double>& PROTO_lattice,
      vector<uint>& im1, vector<uint>& im2, vector<double>& min_dists, 
      const int& type_match); //DX 20190716
  void clusterize(const xstructure& xstr1, const vector<uint>& im1, vector<string>& TYPE1,
      vector<uint>& QTA1, vector<vector<uint> >& I1);
  bool sameAtomType(const xstructure& xstr1, const xstructure& xstr2, const vector<uint>& im1,
      const vector<uint>& im2, const int& type_match);
  bool cleanMatch(const vector<uint>& im1);
  void cellDiagonal(xstructure& xstr, vector<double>& diag_sum, vector<double>& diag_diff, const double& scale);
  void cellDiagonal(xmatrix<double>& lattice, vector<double>& diag_sum, vector<double>& diag_diff, const double& scale);
  double latticeDeviation(const vector<double>& diag_sum1,const vector<double>& diag_sum2, 
      const vector<double>& diag_diff1,const vector<double>& diag_diff2); 
  vector<AtomEnvironment> computeLFAEnvironment(const xstructure& xstr, bool unique_only=true); //DX 20190711
  bool compatibleEnvironmentSets(const vector<AtomEnvironment>& env_set1, 
      const vector<AtomEnvironment>& env_set2, bool same_species, bool exact_match); //DX 20190711
  bool compatibleEnvironments(const AtomEnvironment& env_1, 
      const AtomEnvironment& env_2, bool same_species, bool compare_frequency, bool exact_match); //DX 20190711
  bool compatibleEnvironments(const AtomEnvironment& env_1, 
      const AtomEnvironment& env_2, vector<vector<string> > & matched_species, 
      bool same_species, bool compare_frequency, bool exact_match); //DX 20190711
  vector<vector<double> > getAnglesBetweenMixedSpeciesEnvironments(const vector<vector<xvector<double> > >& neighbor_coordinates); //DX 20190715
  bool compatibleNearestNeighborTypesEnvironments(const vector<vector<double> >& nn_lfa_with_types_1,
      const vector<vector<double> >& nn_lfa_with_types_2,
      int type_match);
  vector<AtomEnvironment> getUniqueTypesAtomEnvironmentForLFA(const xstructure& xstr, const string lfa, const vector<string>& LFA); //DX 20190711
  double shortestDistanceRestrictType(const xstructure& xstr, const uint& k, string lfa); //DX 20190710
  double shortestDistanceRestrictType(const xstructure& xstr, const uint& k, uint& frequency, vector<xvector<double> >& coordinates, string lfa); //DX 20190710
  vector<double> computeNearestNeighbors(xstructure& xstr);
  double shortestDistance(const xstructure& xstr, const uint& k);
  void coordinateDeviation(const xstructure& xstr1, const xstructure& xstr2,
      const vector<double>& all_nn1, const vector<double>& all_nn_proto,
      const vector<uint>& indexMatch1, const vector<uint>& indexMatch2, vector<double>& min_dists,
      double& cd, double& fail_figure);
  void magneticDeviation(const xstructure& xstr1, const xstructure& xstr2, 
      const vector<uint>& indexMatch1, const vector<uint>& indexMatch2,
      double& magnetic_deviation, double& magnetic_fail); //DX 20190801
  double computeMisfit(const double& dev, const double& dis, const double& fail);
  double computeMagneticMisfit(const double dev, const double dis, const double fail, const double mag_dis, const double mag_fail); //DX 20190801
  void printMatch(const vector<uint>& indexMatch1, const vector<uint>& indexMatch2,
      const vector<double>& distances,
      const xstructure& PROTO, const xstructure& xstr1, ostream& oss); //DX 20190802 - added distances
  xvector<double> bringCoordinateInCell(xvector<double>& coord);
  double checkLatticeDeviation(double& xstr1_vol, xmatrix<double>& q2,vector<double>& D1,vector<double>& F1);
  bool quadrupletPeriodic(const xmatrix<double>& quad, const xstructure& lfa_supercell, const int& i, 
      const int& j, const int& k, const int& w);
  xstructure GetLFASupercell(const xstructure& xstr, const xvector<int>& dims, const string& lfa_name); //DX 20190530
  bool atomInCell(const _atom& atom); //DX 20190717 
  bool vectorPeriodic(const xvector<double>& vec, const xstructure& lfa_supercell, const int& i, 
      const int& j);
  // [OBSOLETE - DX 20190717] void threadGeneration(const uint& num_proc,xmatrix<double>& q1, xstructure& xstr2, 
  // [OBSOLETE - DX 20190717]                       vector<xstructure> &vprotos, xstructure &xstr1, const int& type_match, 
  // [OBSOLETE - DX 20190717]                       const bool& optimize_match, double& minMis, ostream& oss);
  void latticeAndOriginSearch(xstructure& xstr1, xstructure& xstr2, const uint& num_proc,xmatrix<double>& q1, vector<xstructure> &vprotos, 
      double& minMis, int type_match, bool optimize_match, ostream& oss); //DX 20190530
  void quadrupletSearch(const xmatrix<double>& q1, const xstructure& xstr_LFA_supercell,
      const xstructure& xstr2,
      vector<xvector<double> >& lattice_vecs, vector<vector<uint> >& ij_index);
  //  void quadrupletSearch(std::atomic_bool &misfit_in_threshold_found, const string& lfa,
  //                        const xmatrix<double>& q1,
  //                        const xstructure& xstr_LFA_supercell,
  //                        vector<xstructure>& vprotos, xstructure& xstr1, const xstructure& xstr2,
  //                        const int& type_match, double& possible_minMis, ostream& oss,
  //                        vector<xvector<double> >& lattice_vecs,
  //                        vector<vector<uint> >& ij_index); // 1 bracket
  bool buildSimilarLattices(vector<xvector<double> >& translation_vectors, xmatrix<double>& q1, double& xstr1_vol, double& abs_det_q1, xvector<double>& abc_angles_q1, vector<xmatrix<double> >& lattices, vector<xmatrix<double> >& clattices, vector<double>& latt_devs, const bool& optimize_match);
  //  bool structureSearch(const string& lfa,
  //                        const vector<double>& all_nn1,
  //			                  const xstructure& xstr_supercell, //DX 20190530 - added "_supercell"; more descriptive 
  //			                  vector<xstructure>& vprotos, xstructure& xstr1, 
  //                        const int& type_match, double& possible_minMis,
  //                        vector<xmatrix<double> >& lattices,
  //                        vector<xmatrix<double> >& clattices,
  //                        vector<double>& latt_devs,
  //			                  const bool& optimize_match,
  //                        const uint& start_index, const uint& end_index); //DX 20190530 - new 
  bool structureSearch(
      const xstructure& xstr1, 
      const xstructure& xstr_supercell, //DX 20190530 - added "_supercell"; more descriptive 
      const vector<double>& all_nn1,
      const string& lfa, 
      const int type_match, 
      const vector<xmatrix<double> >& lattices,
      const vector<xmatrix<double> >& clattices, 
      const vector<double>& latt_devs, 
      const uint start_index, const uint end_index,
      double& min_misfit, double& min_latt_dev, 
      double& min_coordinate_dis, double& min_failure,
      double& min_magnetic_dis, double& min_magnetic_failure,
      vector<uint>& index_match_1, vector<uint>& index_match_2,
      vector<double>& min_distances,
      vector<xstructure>& vprotos,
      bool optimize_match);  //DX 20190802 - new input format
  // [OBSOLETE - DX 20190717]bool structureSearch(const string& lfa,
  // [OBSOLETE - DX 20190717]                      const vector<double>& all_nn1,
  // [OBSOLETE - DX 20190717]                      const xstructure& xstr,
  // [OBSOLETE - DX 20190717]                      vector<xstructure>& vprotos, xstructure& xstr1, const xstructure& xstr2,
  // [OBSOLETE - DX 20190717]                      const int& type_match, double& possible_minMis,
  // [OBSOLETE - DX 20190717]                      vector<xmatrix<double> >& lattices,
  // [OBSOLETE - DX 20190717]                      vector<xmatrix<double> >& clattices,
  // [OBSOLETE - DX 20190717]                      vector<double>& latt_devs,
  // [OBSOLETE - DX 20190717]                      const bool& optimize_match);

  void getPrototypeDesignations(vector<StructurePrototype>& prototypes); //DX 20190725
  void getPrototypeDesignationsInRange(vector<StructurePrototype>& prototypes, uint start_index, uint end_index); //DX 20190725

} // end namespace

#endif 

// ***************************************************************************
// AFLOW_COMPARE_STRUCTURE
// ***************************************************************************

