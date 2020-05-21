// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow DAVID HICKS - Duke University 2014-2020                 *
// *                                                                         *
// ***************************************************************************
// AFLOW-XtalMatch (compare crystal structures) - Functions
// Written by David Hicks (david.hicks@duke.edu) 
// Contributors: Carlo De Santo

#include<iostream>
#include "aflow_pflow.h"
#include "aflow.h"
//#include <math.h> //CO20191217 - this is a STANDARD library  //ME20200514 - not needed
#include<vector>
#include<string>
//#include <stdatomic.h> // Needed to communicate between threads
//#include<atomic>

//DX20190314 [OBSOLETE] using aurostd::isequal;
//DX20190314 [OBSOLETE] using aurostd::string2utype;
//DX20190314 [OBSOLETE] using aurostd::FileExist;

#ifndef __AFLOW_COMPARE_STRUCTURE_H_
#define __AFLOW_COMPARE_STRUCTURE_H_

// Create the version of GCC, we will uset it for multithread parts of code,
// to check if the current compiling version of gcc is able to compile the
// std::thead features
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

#define JSON_MODE 0
#define _CALCULATE_MAGNETIC_MISFIT_ false
#define _SPIN_TOL_ 0.1

#define _COMPARE_DATABASE_GEOMETRY_ORIGINAL_ 0
#define _COMPARE_DATABASE_GEOMETRY_RELAX1_ 1
#define _COMPARE_DATABASE_GEOMETRY_MOST_RELAXED_ 2

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
//DX20191120 [MOVED TO aflow.h]

//DX20191212 - added 
struct structure_misfit {
  bool is_magnetic_misfit;              // boolean indicating if a magnetic system and using magnetic as misfit
  double misfit;                        // Burzlaff's structural misfit (=1-(1-lattice_deviation)(1-coordinate_displacement)(1-failure))
  double lattice_deviation;             // Burzlaff's lattice deviation, captures differences between lattices
  double coordinate_displacement;       // Burzlaff's coordinate displacement; captures differences between atom positions (relatively close together)
  double failure;                       // Burzlaff's figure of failure; captures differences between atom positions (significantly far apart)
  double magnetic_misfit;               //DX's magnetic misfit (inspired by Burzlaff's misfit; =1-(1-magnetic_displacement)(1-magnetic_failure))
  double magnetic_displacement;         //DX's magnetic displacement; captures differences between magnetic moment magnitude (and angle for non-collinear)
  double magnetic_failure;              //DX's magnetic failure; captures spin flip differences 
};
namespace compare{
  structure_misfit initialize_misfit_struct(bool magnetic=false);
}

//DX20200225 - temp struct; working on more robust scheme
struct matching_structure {
  string name;
  double misfit;
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
    string structure_representative_source;                                                 // string indicating where structure came from, i.e., input, auid, aflow protos, etc.
    uint structure_representative_relaxation_step;                                          // number specifying relaxation step of the representative structure (0=original, 1=one relaxation, 2=most relaxed)
    uint number_compounds_matching_representative;                                          // number of compounds that match with the representative structure
    int ntypes;                                                                             // number of types in prototype
    vector<string> elements;                                                                // list of elements exhibiting in this protoype (from representative and duplicate structures)
    vector<uint> stoichiometry;                                                             // reduced stoichiometry of prototype
    uint natoms;                                                                            // number of atoms in the prototype (from the representative structure; not necessarily reduced)
    //DX20191111 [OBSOLETE] vector<string> atom_decorations_unique;                         // unique atom decorations (permutations) of prototype
    vector<vector<string> > atom_decorations_equivalent;                                    // equivalent atom decorations (permutations) of prototype
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
    vector<string> structures_duplicate_source;                                             // vector of strings indicating where structure came from, i.e., input, auid, aflow protos, etc.
    vector<uint> structures_duplicate_relaxation_step;                                      // vector of numbers specifying relaxation steps of the duplicate structures (0=original, 1=one relaxation, 2=most relaxed)
    vector<vector<GroupedWyckoffPosition> > structures_duplicate_grouped_Wyckoff_positions; // Wyckoff positions grouped by site type
    vector<uint> number_compounds_matching_duplicate;                                       // vector of number of compounds that match with the duplicate structures
    vector<string> duplicate_comparison_logs; //DX20190506                                  // vector of comparison logs for duplicate structures, CAREFUL: very long string //DX20190506
    vector<string> structures_family_names;                                                 // vector of names for structures within the same family
    vector<xstructure> structures_family;                                                   // vector of xstructures for structures within the same family
    vector<bool> structures_family_generated;                                               // vector of booleans indicating if xstructure is generated for structures within the same family
    vector<string> structures_family_source;                                                // vector of strings indicating where structure came from, i.e., input, auid, aflow protos, etc.
    vector<uint> structures_family_relaxation_step;                                         // vector of numbers specifying relaxation steps of the same family structures (0=original, 1=one relaxation, 2=most relaxed)
    vector<vector<GroupedWyckoffPosition> > structures_family_grouped_Wyckoff_positions;    // Wyckoff positions grouped by site type
    vector<uint> number_compounds_matching_family;                                          // vector of number of compounds that match with the same family structures
    vector<string> family_comparison_logs; //DX20190506                                     // vector of comparison logs for same family structures, CAREFUL: very long string //DX20190506
    vector<structure_misfit> structure_misfits_duplicate;                                   // vector of C++ structs containing broken-down structural misfit information between the duplicate and representative structures //DX20191217
    vector<structure_misfit> structure_misfits_family;                                      // vector of C++ structs containing broken-down structural misfit information between the family and representative structures //DX20191217
    vector<string> property_names;                                                          // vector of property names (if using AFLUX)
    vector<string> property_units;                                                          // vector of property units (if using AFLUX)
    vector<string> properties_structure_representative;                                     // vector of property values for the representative structure (if using AFLUX)
    vector<vector<string> > properties_structures_duplicate;                                // vector of property values for the duplicate structures (if using AFLUX)
    vector<vector<string> > properties_structures_family; //DX20190506                      // vector of property values for the family structures (if using AFLUX) //DX20190425
    // functions
    uint numberOfDuplicates() const; //DX20190506                                           // return the number of duplicate structures for this prototype (i.e., checks misfit value)
    uint numberOfComparisons(); //DX20181221                                                // return the number of comparisons for this prototype
    bool isSymmetryCalculated(); //DX20190228
    bool isLFAEnvironmentCalculated(); //DX20191105
    bool calculateSymmetry(); //DX20190118                                                  // calculate space group symmetry and populate symmetry attributes for StructurePrototype
    bool addStructurePrototypeAsDuplicate(StructurePrototype& b);                           // combine structure with another StructurePrototype object as a potential duplicate
    void putDuplicateAsFamily(uint index, bool keep_generated=false);                       // make duplicate structure a same family structure in the same StructurePrototypeObject //DX20190814 
    bool copyPrototypeInformation(StructurePrototype& b);                                   // copy prototype information from one StructurePrototype object to another
    bool putDuplicateAsRepresentative(StructurePrototype& b, uint& index);                  // make duplicate structure a representative structure in a new StructurePrototype object
    //bool putRepresentativeAsDuplicate(StructurePrototype& b);                             // make representative structure a duplicate structure in a new StructurePrototype object //DX20190730
    bool copyDuplicate(StructurePrototype& b, uint& index, bool copy_misfit=false);         // copy duplicate info to another StructurePrototype object
    bool removeNonDuplicate(uint& index);                                                   // remove non-duplicate structures 
    bool removeDuplicates(bool remove_duplicate_count);                                     // remove duplicate structure information
  private:
    void free();                                                                            // free operator
    void copy(const StructurePrototype& b);                                                 // copy constructor
};


namespace compare{

  // ===== Main functions ===== //
  string compareMultipleStructures(const aurostd::xoption& vpflow, ostream& logstream=cout); //DX //DX20190425
  vector<StructurePrototype> compareStructuresFromStructureList(vector<string>& filenames, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, const aurostd::xoption& comparison_options); //DX20200103 - condensed bools to xoptions
  vector<StructurePrototype> compareStructuresFromDirectory(string& directory, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, const aurostd::xoption& comparison_options); //DX20200103 - condensed bools to xoptions
  vector<StructurePrototype> compareStructuresFromFile(string& filename, vector<string>& magmoms_for_systems, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, const aurostd::xoption& comparison_options); //DX20200103 - condensed bools to xoptions
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& all_structures, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, string& directory, ostream& logstream=cout); //DX20190319 - added FileMESSAGE
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& all_structures, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, string& directory, const aurostd::xoption& comparison_options, ostream& logstream=cout); //DX20200103 - condensed bools to xoptions 
  bool aflowCompareStructure(const uint& num_proc, const xstructure& xstr1, const xstructure& xstr2, 
      bool same_species, bool scale_volume, bool optimize_match, double& final_misfit, ostream& oss); //Main function //DX20191108 - remove const & from bools //DX20191122 - move ostream to end and add default
  bool aflowCompareStructure(const uint& num_proc, const xstructure& xstr1, const xstructure& xstr2, 
      bool same_species, bool scale_volume, bool optimize_match, double& final_misfit, structure_misfit& final_misfit_info, ostream& oss); //Main function //DX20191108 - remove const & from bools //DX20191122 - move ostream to end and add default //DX20191210 - added lattice_dev, coordinate_dis, and failure
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species); //Overco, returns true (match), false (no match) //DX20191108 - remove const & from bools
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool scale_volume, bool optmize_match);  //DX20191108 - remove const & from bools
  double aflowCompareStructureMisfit(const xstructure& xstr1, const xstructure& xstr2, bool same_species, bool optimize_match); //Overloaded, returns misfit value //DX20191108 - remove const & from bools

  // ---------------------------------------------------------------------------
  // comparisons to AFLOW database 
  vector<StructurePrototype> compare2database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20200225
  vector<StructurePrototype> compare2database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20200225
  bool isMatchingStructureInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20200225
  bool isMatchingStructureInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20200225
  vector<matching_structure> matchingStructuresInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20200225
  vector<matching_structure> matchingStructuresInDatabase(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20200225
  string printCompare2Database(istream& input, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20200225
  string printCompare2Database(istream& input, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout);  //CO20200225
  string printCompare2Database(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20200225

  // ---------------------------------------------------------------------------
  // comparisons between entries in AFLOW database 
  string compareDatabaseEntries(const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20191125
  string compareDatabaseEntries(const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20191125

  // ---------------------------------------------------------------------------
  // comparisons to AFLOW prototype library 
  vector<StructurePrototype> compare2prototypes(istream& input, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20181004 //DX20190314 - changed return value
  vector<StructurePrototype> compare2prototypes(istream& input, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20181004 //DX20190314 - changed return value
  vector<StructurePrototype> compare2prototypes(const xstructure& xstrIN, const aurostd::xoption& vpflow, ostream& logstream=cout); //DX20190314 - overloaded 
  vector<StructurePrototype> compare2prototypes(const xstructure& xstrIN, const aurostd::xoption& vpflow, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20190314 - overloaded 
  string printMatchingPrototypes(xstructure& xstr, const aurostd::xoption& vpflow); //DX20190314 
  string printMatchingPrototypes(istream& cin, const aurostd::xoption& vpflow); //DX20190314 
  vector<string> getMatchingPrototypes(xstructure& xstr, string& catalog); //DX20190314 

  // ---------------------------------------------------------------------------
  // permutaion comparisons
  string comparePermutations(istream& input, const aurostd::xoption& vpflow); //DX //DX20181004
  vector<string> getUniquePermutations(xstructure& xstr);
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc); // add number of threads
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc, bool& optimize_match); // add number of threads and optimize matches
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc, bool& optimize_match, bool& print_misfit, ostream& oss, ofstream& FileMESSAGE); // full function //DX20190319 - added FileMESSAGE

  // ---------------------------------------------------------------------------
  // isopointal AFLOW prototype functions 
  string isopointalPrototypes(istream& input, const aurostd::xoption& vpflow); //DX20200131 
  vector<string> getIsopointalPrototypes(xstructure& xstr, string& catalog); //DX20200131 

  //  string CompareStructures(aurostd::xoption& vpflow);
  //  string CompareStructureDirectory(aurostd::xoption& vpflow);

  aurostd::xoption loadDefaultComparisonOptions(string mode=""); //DX20200103
  // ===== Compare Directory Functions ===== //
  vector<StructurePrototype> loadStructuresFromStructureList(const vector<string>& filenames, const vector<string>& magmoms_for_systems, bool same_species, ostream& logstream=cout); //DX20191122
  vector<StructurePrototype> loadStructuresFromStructureList(const vector<string>& filenames, const vector<string>& magmoms_for_systems, bool same_species, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20190424 //DX20190801 - added vector<string>& magmoms_for_systems //DX20191122 - added ostream and consts
  vector<StructurePrototype> loadStructuresFromDirectory(const string& directory, const vector<string>& magmoms_for_systems, bool same_species, ostream& logstream=cout); //DX20191122 
  vector<StructurePrototype> loadStructuresFromDirectory(const string& directory, const vector<string>& magmoms_for_systems, bool same_species, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20190319 - added FileMESSAGE //DX20190801 - added vector<string>& magmoms_for_systems //DX20191122 - added ostream and consts
  vector<StructurePrototype> loadStructuresFromFile(const string& directory, const vector<string>& magmoms_for_systems, bool same_species, ostream& logstream=cout); //DX20191122
  vector<StructurePrototype> loadStructuresFromFile(const string& directory, const vector<string>& magmoms_for_systems, bool same_species, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20190319 - added FileMESSAGE //DX20190801 - added vector<string>& magmoms_for_systems, //DX20191122 - added ostream and consts
  void generateStructures(vector<StructurePrototype>& structures, ostream& oss=cout, uint start_index=0, uint end_index=AUROSTD_MAX_UINT); //DX20191122
  bool generateStructure(string& structure_name, string& structure_source, uint relaxation_step, xstructure& structure, ostream& oss); //DX20200429 - added relaxation_step
  void removeNonGeneratedStructures(vector<StructurePrototype>& structures); //DX20191105
  vector<uint> getStoichiometry(string& compositions, const bool& same_species);
  bool addAFLOWPrototypes2StructurePrototypeVector(vector<StructurePrototype>& all_structures, vector<string>& vlabel);
  string getCompoundName(xstructure& xstr, bool remove_ones=false);
  string getCompoundName(vector<string>& elements, vector<uint>& stoichiometry, bool remove_ones=false);
  vector<uint> getStoichiometry(const xstructure& xstr, const bool& same_species);
  vector<string> getElements(xstructure& xstr);
  //DX20191125 [OBSOLETE - USING AUROSTD VERSION] vector<uint> gcdStoich(const vector<uint>& numbers); //DX20181009
  //DX20191125 [OBSOLETE - USING AUROSTD VERSION] vector<uint> gcdStoich(const deque<int>& numbers);
  //DX20191108 [OBSOLETE - switching to getThreadDistribution] bool prepareSymmetryThreads(vector<xstructure>& vxstrs, uint& num_proc,
  //DX20191108 [OBSOLETE - switching to getThreadDistribution]     vector<uint>& start_indices, vector<uint>& end_indices);
  //DX20191108 [OBSOLETE - switching to getThreadDistribution] bool prepareSymmetryThreads(uint& number_of_structures, uint& num_proc,
  //DX20191108 [OBSOLETE - switching to getThreadDistribution]     vector<uint>& start_indices, vector<uint>& end_indices);
  //DX20191108 [OBSOLETE - switching to getThreadDistribution] bool splitTaskIntoThreads(uint& number_of_tasks, uint& num_proc,              
  //DX20191108 [OBSOLETE - switching to getThreadDistribution]     vector<uint>& start_indices, vector<uint>& end_indices); //DX20190530 - renamed
  void calculateSymmetries(vector<xstructure>& vxstrs, vector<string>& vpearsons, vector<uint>& vsgroups, 
      vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions, uint& num_proc);
  void calculateSpaceGroups(vector<StructurePrototype>& structures, uint start_index=0, uint end_index=AUROSTD_MAX_UINT, uint setting=0); //DX20191230 added setting
  void calculateLFAEnvironments(vector<StructurePrototype>& structures, uint num_proc);  //DX20191105
  void calculateSymmetries(vector<StructurePrototype>& structures, uint& num_proc);  //DX20190118
  void calculateSymmetry(xstructure& xstr, vector<string>& vpearsons, vector<uint>& vsgroups,
      vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions);
  bool groupWyckoffPositions(xstructure& xstr, vector<GroupedWyckoffPosition>& grouped_positions);
  bool groupWyckoffPositionsFromGroupedString(uint& space_group_number, uint& setting, vector<vector<string> >& grouped_Wyckoff_string, vector<GroupedWyckoffPosition>& grouped_positions);
  string printWyckoffString(const vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize=false);
  vector<GroupedWyckoffPosition> sortSiteSymmetryOfGroupedWyckoffPositions(const vector<GroupedWyckoffPosition>& grouped_Wyckoffs); //DX20190219
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
  //vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& comparison_schemes, 
  //vector<xstructure>& vxstrs, const bool& same_species, 
  //const vector< vector<string> >& vvelements,
  //vector< vector<uint> >& vstoichs, vector<string>& vpearsons, 
  //vector<uint>& vsgroups, 
  //vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
  //const string& directory, const vector<string>& vfiles, 
  //vector<bool>& vstructures_generated,
  //vector<string>& vstructures_source,
  //const bool& ignore_symmetry, const bool& ignore_Wyckoff,
  //const bool& structures_generated);
  //void createStructurePrototypes(vector<StructurePrototype>& comparison_schemes, 
  //vector<xstructure>& vxstrs, const bool& same_species, 
  //const vector< vector<string> >& vvelements,
  //vector< vector<uint> >& vstoichs, vector<string>& vpearsons, 
  //vector<uint>& vsgroups, 
  //vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
  //vector<string>& property_names,
  //vector<string>& property_units,
  //vector<vector<string> >& property_values,
  //const string& directory, const vector<string>& vfiles,
  //vector<bool>& vstructures_generated,
  //vector<string>& vstructures_source,
  //const bool& ignore_symmetry, const bool& ignore_Wyckoff,
  //const bool& structures_generated);
  void createStructurePrototypes(vector<StructurePrototype>& comparison_schemes,
      vector<xstructure>& vxstrs, const bool& same_species,
      const vector< vector<string> >& vvelements,
      vector< vector<uint> >& vstoichs, vector<string>& vpearsons,
      vector<uint>& vsgroups,
      vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
      const string& directory, const vector<string>& vfiles,
      vector<bool>& vstructures_generated,
      vector<string>& vstructures_source,
      vector<uint>& vstructures_relaxation_step,
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
      vector<string>& vstructures_source,
      vector<uint>& vstructures_relaxation_step,
      const bool& ignore_symmetry, const bool& ignore_Wyckoff);
  bool structuresCompatible(const StructurePrototype& structure1,
      const StructurePrototype& structure2, bool same_species, 
      bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool ignore_environment_angles); //DX20190730 //DX20200320 - added environment angles
  vector<StructurePrototype> groupStructurePrototypes(vector<StructurePrototype>& structures, 
      bool same_species, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment, bool ignore_environment_angles, bool duplicates_removed); //DX20190731 - remove const and & //DX20190830 - added duplicates_removed //DX20200320 - added environment angles

  vector<StructurePrototype> checkForBetterMatches(vector<StructurePrototype>& prototype_schemes, 
      ostream& oss, uint& num_proc, bool check_for_better_matches, bool same_species,
      const aurostd::xoption& comparison_options, bool quiet=false, ostream& logstream=cout); //DX20200103 - condensed bools to xoptions
  vector<StructurePrototype> checkForBetterMatches(vector<StructurePrototype>& prototype_schemes, 
      ostream& oss, uint& num_proc, bool check_for_better_matches, bool same_species,
      const aurostd::xoption& comparison_options, ofstream& FileMESSAGE, bool quiet=false, ostream& logstream=cout); //DX20200103 - condensed bools to xoptions
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
      bool store_comparison_logs); //DX20190822 - added comparison log bool
  //OLD vector<StructurePrototype> runComparisonScheme(uint& num_proc, vector<StructurePrototype>& comparison_schemes, const bool& same_species, 
  //OLD                                                const bool& scale_volume, const bool& optimize_match, const bool& single_comparison_round, const bool& structures_generated, const bool& ICSD_comparison, ostringstream& oss);
  vector<StructurePrototype> runComparisonScheme(vector<StructurePrototype>& comparison_schemes, 
      bool same_species, uint num_proc, const aurostd::xoption& comparison_options,
      ostream& oss, bool quiet=false, ostream& logstream=cout); //DX20200103 - condensed bools to xoptions
  vector<StructurePrototype> runComparisonScheme(vector<StructurePrototype>& comparison_schemes, 
      bool same_species, uint num_proc, const aurostd::xoption& comparison_options, 
      ostream& oss, ofstream& FileMESSAGE, bool quiet=false, ostream& logstream=cout); //DX20200103 - condensed bools to xoptions
  vector<std::pair<uint,uint> > calculateDivisors(const int& number);
  bool checkNumberOfGroupings(vector<StructurePrototype>& comparison_schemes, uint number);
  void createStructurePermutations(vector<StructurePrototype>& comparison_schemes, const vector<vector<string> >& name_order,
      vector<vector<GroupedWyckoffPosition> >& permutation_grouped_Wyckoff_positions,
      const vector<xstructure>& vxstrs, const bool& same_species);
  bool makeRepresentativeEvenPermutation(vector<StructurePrototype>& comparison_schemes, const vector<vector<string> >& name_order);  
  int numberMismatches(const vector<StructurePrototype> comparison_schemes);
  void appendStructurePrototypes(vector<StructurePrototype>& comparison_schemes, 
      vector<StructurePrototype>& final_prototypes,
      bool clean_unmatched, //DX20190506
      bool quiet=false,
      ostream& logstream=cout);
  void appendStructurePrototypes(vector<StructurePrototype>& comparison_schemes, 
      vector<StructurePrototype>& final_prototypes,
      bool clean_unmatched, 
      ofstream& FileMESSAGE,
      bool quiet=false,
      ostream& logstream=cout); //DX20191125
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
      string mode="TEXT"); //DX20190802

  bool SVD(xmatrix<double>& A);
  string findICSDName(string& name);
  string findMinimumICSDEntry(vector<string>& ICSD_entries);
  bool groupSameRatios(vector<int>& stoich, vector<int>& unique_stoich, vector<vector<int> >& type_index);
  //vector<vector<int> > generatePermutations(uint& num_elements, vector<int>& indices);
  vector<StructurePrototype> comparePermutations(StructurePrototype& structure, uint& num_proc, bool optmize_match, ostream& oss, ostream& logstream=cout); //DX20191125 
  vector<StructurePrototype> comparePermutations(StructurePrototype& structure, uint& num_proc, bool optmize_match, ostream& oss, ofstream& FileMESSAGE, ostream& logstream=cout); //DX20190319 - added FileMESSAGE //DX20191125 - added ostream
  vector<StructurePrototype> generatePermutationStructures(StructurePrototype& structure);
  void generatePermutationString(const deque<uint>& stoichiometry, vector<string>& permutation); //DX20190508 //DX20191125 - changed from vector to deque
  void generatePermutationString(const vector<uint>& stoichiometry, vector<string>& permutation); //DX20190508
  bool generatePermutations(uint& num_elements, vector<uint>& indices, vector<string>& names, vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, vector<vector<uint> >& permutations, vector<vector<string> >&name_order, vector<vector<GroupedWyckoffPosition> >& permutation_grouped_Wyckoff_positions);
  bool arePermutationsComparableViaComposition(const xstructure& xstr); //DX20190624 
  bool arePermutationsComparableViaComposition(vector<uint>& composition, bool reduce_composition=false); //DX20190624
  bool arePermutationsComparableViaSymmetry(vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions); //DX20190624
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
  //DX20191122 [MOVED TO XATOM] void resetLatticeDimensions(const xmatrix<double>& lattice, double radius, xvector<int>& dims,
  //DX20191122 [MOVED TO XATOM]     vector<xvector<double> >& l1, vector<xvector<double> >& l2, 
  //DX20191122 [MOVED TO XATOM]     vector<xvector<double> >& l3, vector<int>& a_index, 
  //DX20191122 [MOVED TO XATOM]     vector<int>& b_index, vector<int>& c_index); //DX20190705
  //DX20191122 [MOVED TO XATOM] void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
  //DX20191122 [MOVED TO XATOM]     double& min_dist, uint& frequency, vector<xvector<double> >& coordinates); //DX20191105
  //DX20191122 [MOVED TO XATOM] void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
  //DX20191122 [MOVED TO XATOM]     double& min_dist, uint& frequency, vector<xvector<double> >& coordinates, double radius); //DX20191105
  //DX20191122 [MOVED TO XATOM] void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice, xvector<int>& dims,
  //DX20191122 [MOVED TO XATOM]     vector<xvector<double> >& l1, vector<xvector<double> >& l2, vector<xvector<double> >& l3, 
  //DX20191122 [MOVED TO XATOM]     vector<int>& a_index, vector<int>& b_index, vector<int>& c_index, 
  //DX20191122 [MOVED TO XATOM]     double& min_dist, uint& frequency, vector<xvector<double> >& coordinates,
  //DX20191122 [MOVED TO XATOM]     double radius); //DX20191105
  //DX20191122 [MOVED TO XATOM] void minimumCoordinationShell(const xstructure& xstr, uint center_index, 
  //DX20191122 [MOVED TO XATOM]     double& min_dist, uint& frequency, vector<xvector<double> >& coordinates); //DX20191105
  //DX20191122 [MOVED TO XATOM] void minimumCoordinationShell(const xstructure& xstr, uint center_index, 
  //DX20191122 [MOVED TO XATOM]     double& min_dist, uint& frequency, vector<xvector<double> >& coordinates, const string& type); //DX20191105
  xvector<double> centroid_with_PBC(const xstructure& xstr);
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, const xmatrix<double>& lattice);
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, vector<double>& weights,
      const xmatrix<double>& lattice);
  //bool findMatch(const xstructure& xstr1, const xstructure& PROTO,vector<uint>& im1, vector<uint>& im2, vector<double>& min_dists, const int& type_match);
  bool findMatch(const deque<_atom>& xstr1_atoms, const deque<_atom>& PROTO_atoms,
      const xmatrix<double>& PROTO_lattice,
      vector<uint>& mapping_index_str1, vector<uint>& mapping_index_str2, vector<double>& min_dists,
      const int& type_match); //DX20190716
  void clusterize(const xstructure& xstr1, const vector<uint>& im1, vector<string>& TYPE1,
      vector<uint>& QTA1, vector<vector<uint> >& I1);
  bool sameAtomType(const xstructure& xstr1, const xstructure& xstr2, const vector<uint>& im1,
      const vector<uint>& im2, const int& type_match);
  bool cleanMatch(const vector<uint>& im1);
  void cellDiagonal(xstructure& xstr, vector<double>& diag_sum, vector<double>& diag_diff, const double& scale);
  void cellDiagonal(xmatrix<double>& lattice, vector<double>& diag_sum, vector<double>& diag_diff, const double& scale);
  double latticeDeviation(const vector<double>& diag_sum1,const vector<double>& diag_sum2, 
      const vector<double>& diag_diff1,const vector<double>& diag_diff2); 
  void computeLFAEnvironments(vector<StructurePrototype>& structures, uint start_index=0, uint end_index=AUROSTD_MAX_UINT);  //DX20191122
  vector<AtomEnvironment> computeLFAEnvironment(const xstructure& xstr, bool unique_only=true); //DX20190711
  bool compatibleEnvironmentSets(const vector<AtomEnvironment>& env_set1, 
      const vector<AtomEnvironment>& env_set2, bool same_species, bool ignore_environment_angles, bool exact_match); //DX20190711 //DX20200320 - added environment angles
  bool compatibleEnvironments(const AtomEnvironment& env_1, 
      const AtomEnvironment& env_2, bool same_species, bool ignore_environment_angles, bool exact_match); //DX20190711 //DX20200320 - added environment angles
  bool compatibleEnvironments(const AtomEnvironment& env_1, 
      const AtomEnvironment& env_2, vector<vector<string> > & matched_species, 
      bool same_species, bool ignore_environment_angles, bool exact_match); //DX20190711 //DX20200320 - added environment angles
  vector<vector<double> > getAnglesBetweenMixedSpeciesEnvironments(const vector<vector<xvector<double> > >& neighbor_coordinates); //DX20190715
  bool compatibleNearestNeighborTypesEnvironments(const vector<vector<double> >& nn_lfa_with_types_1,
      const vector<vector<double> >& nn_lfa_with_types_2,
      int type_match);
  vector<AtomEnvironment> getUniqueTypesAtomEnvironmentForLFA(const xstructure& xstr, const string lfa, const vector<string>& LFA); //DX20190711
  double shortestDistanceRestrictType(const xstructure& xstr, const uint& k, string lfa); //DX20190710
  double shortestDistanceRestrictType(const xstructure& xstr, const uint& k, uint& frequency, vector<xvector<double> >& coordinates, string lfa); //DX20190710
  vector<double> computeNearestNeighbors(xstructure& xstr);
  double shortestDistance(const xstructure& xstr, const uint& k);
  void coordinateDeviation(const xstructure& xstr1, const xstructure& xstr2,
      const vector<double>& all_nn1, const vector<double>& all_nn_proto,
      const vector<uint>& indexMatch1, const vector<uint>& indexMatch2, vector<double>& min_dists,
      double& cd, double& fail_figure);
  void magneticDeviation(const xstructure& xstr1, const xstructure& xstr2, 
      const vector<uint>& indexMatch1, const vector<uint>& indexMatch2,
      double& magnetic_deviation, double& magnetic_fail); //DX20190801
  double computeMisfit(const double& dev, const double& dis, const double& fail);
  double computeMagneticMisfit(const double dev, const double dis, const double fail, const double mag_dis, const double mag_fail); //DX20190801
  void printMatch(const vector<uint>& indexMatch1, const vector<uint>& indexMatch2,
      const vector<double>& distances,
      const xstructure& PROTO, const xstructure& xstr1, ostream& oss); //DX20190802 - added distances
  xvector<double> bringCoordinateInCell(xvector<double>& coord);
  double checkLatticeDeviation(double& xstr1_vol, xmatrix<double>& q2,vector<double>& D1,vector<double>& F1);
  bool quadrupletPeriodic(const xmatrix<double>& quad, const xstructure& lfa_supercell, const int& i, 
      const int& j, const int& k, const int& w);
  xstructure GetLFASupercell(const xstructure& xstr, const xvector<int>& dims, const string& lfa_name); //DX20190530
  //DX20191125 [OBSOLETE - MOVED TO XATOM] bool atomInCell(const _atom& atom); //DX20190717 
  bool vectorPeriodic(const xvector<double>& vec, const xstructure& lfa_supercell, const int& i, 
      const int& j);
  // [OBSOLETE - DX20190717] void threadGeneration(const uint& num_proc,xmatrix<double>& q1, xstructure& xstr2, 
  // [OBSOLETE - DX20190717]                       vector<xstructure> &vprotos, xstructure &xstr1, const int& type_match, 
  // [OBSOLETE - DX20190717]                       const bool& optimize_match, double& minMis, ostream& oss);
  //void latticeAndOriginSearch(xstructure& xstr1, xstructure& xstr2, const uint& num_proc,xmatrix<double>& q1, vector<xstructure> &vprotos, 
  //    double& minMis, int type_match, bool optimize_match, ostream& oss); //DX20190530
  void latticeAndOriginSearch(xstructure& xstr1, xstructure& xstr2, const uint& num_proc,xmatrix<double>& q1, vector<xstructure> &vprotos, 
      structure_misfit& min_misfit_info, int type_match, bool optimize_match, bool scale_volume, ostream& oss); //DX20190530 //DX20191210 - added lattice_dev, coordinate_dis, and failure //DX20200422 - added scale volume
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
  bool buildSimilarLattices(vector<xvector<double> >& translation_vectors, xmatrix<double>& q1, double& xstr1_vol, double& abs_det_q1, xvector<double>& abc_angles_q1, vector<xmatrix<double> >& lattices, vector<xmatrix<double> >& clattices, vector<double>& latt_devs, bool optimize_match, bool scale_volume); //DX20200422 - scale volume
  //  bool structureSearch(const string& lfa,
  //                        const vector<double>& all_nn1,
  //			                  const xstructure& xstr_supercell, //DX20190530 - added "_supercell"; more descriptive 
  //			                  vector<xstructure>& vprotos, xstructure& xstr1, 
  //                        const int& type_match, double& possible_minMis,
  //                        vector<xmatrix<double> >& lattices,
  //                        vector<xmatrix<double> >& clattices,
  //                        vector<double>& latt_devs,
  //			                  const bool& optimize_match,
  //                        const uint& start_index, const uint& end_index); //DX20190530 - new 
  bool structureSearch(
      const xstructure& xstr1, 
      const xstructure& xstr_supercell, //DX20190530 - added "_supercell"; more descriptive 
      const vector<double>& all_nn1,
      const string& lfa, 
      const int type_match, 
      const vector<xmatrix<double> >& lattices,
      const vector<xmatrix<double> >& clattices, 
      const vector<double>& latt_devs, 
      const uint start_index, const uint end_index,
      structure_misfit& min_misfit_info,
      //double& min_misfit, double& min_latt_dev, 
      //double& min_coordinate_dis, double& min_failure,
      //double& min_magnetic_dis, double& min_magnetic_failure,
      vector<uint>& index_match_1, vector<uint>& index_match_2,
      vector<double>& min_distances,
      vector<xstructure>& vprotos,
      bool optimize_match);  //DX20190802 - new input format
  // [OBSOLETE - DX20190717]bool structureSearch(const string& lfa,
  // [OBSOLETE - DX20190717]                      const vector<double>& all_nn1,
  // [OBSOLETE - DX20190717]                      const xstructure& xstr,
  // [OBSOLETE - DX20190717]                      vector<xstructure>& vprotos, xstructure& xstr1, const xstructure& xstr2,
  // [OBSOLETE - DX20190717]                      const int& type_match, double& possible_minMis,
  // [OBSOLETE - DX20190717]                      vector<xmatrix<double> >& lattices,
  // [OBSOLETE - DX20190717]                      vector<xmatrix<double> >& clattices,
  // [OBSOLETE - DX20190717]                      vector<double>& latt_devs,
  // [OBSOLETE - DX20190717]                      const bool& optimize_match);

  void getPrototypeDesignations(vector<StructurePrototype>& prototypes, uint start_index=0, uint end_index=AUROSTD_MAX_UINT); //DX20191122

} // end namespace

#endif 

// ***************************************************************************
// AFLOW_COMPARE_STRUCTURE
// ***************************************************************************

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow DAVID HICKS - Duke University 2014-2020                 *
// *                                                                         *
// ***************************************************************************
