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
//    void Free();
//    void Copy(const GroupedWyckoffPosition& b);
//};

// ===== StructurePrototype Class ===== //
class StructurePrototype{
  public:
    StructurePrototype();                                                                   // constructor operator
    ~StructurePrototype();                                                                  // destructor operator
    friend ostream& operator<<(ostream& oss, const StructurePrototype& StructurePrototype); // stringstream operator (printing)
    const StructurePrototype& operator=(const StructurePrototype& b);                       // assignment operator
    StructurePrototype(const StructurePrototype& b);                                        // copy constructor
    int iomode;                                                                             // mode for printing
    string representative_structure_name;                                                   // name of representative structure
    string representative_structure_compound;                                               // compound name of representative structure (w/reduced stoichometry), e.g., Ag1Br2
    xstructure representative_structure;                                                    // xstructure of representative structure
    bool representative_structure_generated;                                                // boolean indicating if xstructure is generated
    string representative_structure_from;                                                   // string indicating where structure came from, i.e., input, auid, aflow protos, etc.
    uint number_compounds_matching_representative;                                          // number of compounds that match with the representative structure
    int number_types;                                                                       // number of types in prototype
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
    vector<string> duplicate_structures_names;                                              // vector of names for duplicate structures
    vector<string> duplicate_structures_compounds;                                          // vector of compounds names for duplicate structures (w/reduced stoichometry), e.g., Ag1Br2
    vector<xstructure> duplicate_structures;                                                // vector of xstructures for duplicate structures
    vector<bool> duplicate_structures_generated;                                            // vector of booleans indicating if xstructure is generated for duplicate structures
    vector<string> duplicate_structures_from;                                               // vector of strings indicating where structure came from, i.e., input, auid, aflow protos, etc. 
    vector<uint> number_compounds_matching_duplicate;                                       // vector of number of compounds that match with the duplicate structures
    vector<string> duplicate_comparison_logs; //DX 20190506                                              // vector of comparison logs for duplicate structures, CAREFUL: very long string //DX 20190506
    vector<string> family_structures_names;                                                 // vector of names for structures within the same family
    vector<xstructure> family_structures;                                                   // vector of xstructures for structures within the same family
    vector<bool> family_structures_generated;                                               // vector of booleans indicating if xstructure is generated for structures within the same family
    vector<string> family_structures_from;                                                  // vector of strings indicating where structure came from, i.e., input, auid, aflow protos, etc.
    vector<uint> number_compounds_matching_family;                                          // vector of number of compounds that match with the same family structures
    vector<string> family_comparison_logs; //DX 20190506                                                 // vector of comparison logs for same family structures, CAREFUL: very long string //DX 20190506
    vector<double> misfits;                                                                 // vector of misfit values when comparing the duplicate to the representative structure
    vector<double> family_misfits;                                                          // vector of misfit values when comparing the same family structure to the representative structure  
    vector<string> property_names;                                                          // vector of property names (if using AFLUX)
    vector<string> property_units;                                                          // vector of property units (if using AFLUX)
    vector<string> representative_structure_properties;                                     // vector of property values for the representative structure (if using AFLUX)
    vector<vector<string> > duplicate_structures_properties;                                // vector of property values for the duplicate structures (if using AFLUX) 
    vector<vector<string> > family_structures_properties; //DX 20190506                                   // vector of property values for the family structures (if using AFLUX) //DX 20190425 
    // functions
    uint numberOfDuplicates() const; //DX 20190506                                                // return the number of duplicate structures for this prototype (i.e., checks misfit value)
    uint numberOfComparisons(); //DX 20181221                                               // return the number of comparisons for this prototype 
    bool isSymmetryCalculated(); //DX 20190228
    bool calculateSymmetry(); //DX 20190118                                                 // calculate space group symmetry and populate symmetry attributes for StructurePrototype
    bool addStructurePrototypeAsDuplicate(StructurePrototype& b);                           // combine structure with another StructurePrototype object as a potential duplicate
    bool copyPrototypeInformation(StructurePrototype& b);                                   // copy prototype information from one StructurePrototype object to another
    bool putDuplicateAsRepresentative(StructurePrototype& b, uint& index);                  // make duplicate structure a representative structure in a new StructurePrototype object
    bool copyDuplicate(StructurePrototype& b, uint& index);                                 // copy duplicate info to another StructurePrototype object
    bool removeNonDuplicate(uint& index);                                                   // remove non-duplicate structures 
    bool removeDuplicates(bool remove_duplicate_count);                                     // remove duplicate structure information
  private:
    void Free();                                                                            // free operator
    void Copy(const StructurePrototype& b);                                                 // copy constructor
};


namespace compare{

  // ===== Main functions ===== //
  vector<StructurePrototype> compareStructuresFromStructureList(vector<string>& filenames, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species,
    bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds, bool ICSD_comparison); //DX 20190424 //DX 20190504 - added clean unmatched
  vector<StructurePrototype> compareStructuresFromDirectory(string& directory, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species,
    bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds, bool ICSD_comparison); //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched
  vector<StructurePrototype> compareStructuresFromFile(string& filename, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species,
    bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool single_comparison_round, bool clean_unmatched, bool remove_duplicate_compounds, bool ICSD_comparison); //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& all_structures, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, string& directory); //DX 20190319 - added FileMESSAGE
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& all_structures, ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool same_species, string& directory,
    bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, bool single_comparison_round, bool clean_ummatched, bool remove_duplicate_compounds, bool ICSD_comparison); //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched
  bool aflowCompareStructure(const uint& num_proc, const xstructure& xstr1, const xstructure& xstr2, 
                               const bool& same_species, const bool& scale_volume, const bool& fast_match, ostream& oss, double& final_misfit); //Main function
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, const bool &same_species); //Overloaded, returns true (match), false (no match)
  bool aflowCompareStructure(const xstructure& xstr1, const xstructure& xstr2, const bool& same_species, const bool& scale_volume, const bool& fast_match); 
  double aflowCompareStructureMisfit(const xstructure& xstr1, const xstructure& xstr2, const bool &same_species, const bool& fast_match); //Overloaded, returns misfit value

  // permutaion comparisons
  vector<string> getUniquePermutations(xstructure& xstr);
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc); // add number of threads
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc, bool& optimize_match); // add number of threads and optimize matches
  vector<string> getUniquePermutations(xstructure& xstr, uint& num_proc, bool& optimize_match, bool& print_misfit, ostream& oss, ofstream& FileMESSAGE); // full function //DX 20190319 - added FileMESSAGE


//  string CompareStructures(aurostd::xoption& vpflow);
//  string CompareStructureDirectory(aurostd::xoption& vpflow);

  // ===== Compare Directory Functions ===== //
  vector<StructurePrototype> loadStructuresFromStructureList(vector<string>& filenames, bool& same_species, ofstream& FileMESSAGE); //DX 20190424 
  vector<StructurePrototype> loadStructuresFromDirectory(string& directory, bool& same_species, ofstream& FileMESSAGE); //DX 20190319 - added FileMESSAGE
  vector<StructurePrototype> loadStructuresFromFile(string& directory, bool& same_species, ofstream& FileMESSAGE); //DX 20190319 - added FileMESSAGE
  bool generateStructure(string& structure_name, string& structure_from, xstructure& structure, ostream& oss);
  vector<uint> getStoichiometry(string& compositions, const bool& same_species);
  bool addAFLOWPrototypes2StructurePrototypeVector(vector<StructurePrototype>& all_structures, vector<string>& vlabel);
  string getCompoundName(xstructure& xstr, bool remove_ones=false);
  string getCompoundName(vector<string>& elements, vector<uint>& stoichiometry, bool remove_ones=false);
  vector<uint> getStoichiometry(const xstructure& xstr, const bool& same_species);
  vector<string> getElements(xstructure& xstr);
  vector<uint> gcdStoich(const vector<uint>& numbers); //DX 20181009
  vector<uint> gcdStoich(const deque<int>& numbers);
  vector<vector<xstructure> > prepareSymmetryThreads(vector<xstructure>& vxstrs, uint& num_proc);
  bool prepareSymmetryThreads(vector<xstructure>& vxstrs, uint& num_proc,
                              vector<uint>& start_indices, vector<uint>& end_indices);
  bool prepareSymmetryThreads(uint& number_of_structures, uint& num_proc,
                              vector<uint>& start_indices, vector<uint>& end_indices);
  void calculateSymmetries(vector<xstructure>& vxstrs, vector<string>& vpearsons, vector<uint>& vsgroups, 
                         vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions, uint& num_proc);
  void calculateSpaceGroups(vector<StructurePrototype>& structures, uint& start_index, uint& end_index);
  void calculateSymmetries(vector<StructurePrototype>& structures, uint& num_proc);  //DX 20190118
  void calculateSymmetry(xstructure& xstr, vector<string>& vpearsons, vector<uint>& vsgroups,
                         vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions);
  bool groupWyckoffPositions(xstructure& xstr, vector<GroupedWyckoffPosition>& grouped_positions);
  bool groupWyckoffPositionsFromGroupedString(uint& space_group_number, uint& setting, vector<vector<string> >& grouped_Wyckoff_string, vector<GroupedWyckoffPosition>& grouped_positions);
  string printWyckoffString(const vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize=false);
  vector<GroupedWyckoffPosition> sortSiteSymmetryOfGroupedWyckoffPositions(vector<GroupedWyckoffPosition>& grouped_Wyckoffs); //DX 20190219
  bool matchableWyckoffPositions(vector<GroupedWyckoffPosition>& temp_grouped_Wyckoffs,
                                 vector<GroupedWyckoffPosition>& representative_grouped_Wyckoffs,
                                 const bool& same_species);
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
  vector<StructurePrototype> groupStructurePrototypes(vector<StructurePrototype>& structures, 
				 const bool& same_species, 
         const bool& ignore_symmetry, const bool& ignore_Wyckoff);
  
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
      const bool& same_species, 
      const bool& scale_volume, const bool& optimize_match, 
      ostream& oss);
  //OLD vector<StructurePrototype> runComparisonScheme(uint& num_proc, vector<StructurePrototype>& comparison_schemes, const bool& same_species, 
  //OLD                                                const bool& scale_volume, const bool& optimize_match, const bool& single_comparison_round, const bool& structures_generated, const bool& ICSD_comparison, ostringstream& oss);
  vector<StructurePrototype> runComparisonScheme(uint& num_proc, vector<StructurePrototype>& comparison_schemes, const bool& same_species, 
                                                 const bool& scale_volume, const bool& optimize_match, const bool& single_comparison_round, const bool& clean_unmatched, const bool& ICSD_comparison, ostream& oss, ofstream& FileMESSAGE, bool quiet=false); //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched
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

  bool SVD(xmatrix<double>& A);
  string findICSDName(string& name);
  string findMinimumICSDEntry(vector<string>& ICSD_entries);
  bool groupSameRatios(vector<int>& stoich, vector<int>& unique_stoich, vector<vector<int> >& type_index);
  //vector<vector<int> > generatePermutations(uint& num_elements, vector<int>& indices);
  vector<StructurePrototype> comparePermutations(StructurePrototype& structure, uint& num_proc, bool& optmize_match, ostream& oss, ofstream& FileMESSAGE); //DX 20190319 - added FileMESSAGE
  vector<StructurePrototype> generatePermutationStructures(StructurePrototype& structure);
  vector<string> generatePermutationString(vector<uint>& stoichiometry); //DX 20190508
  bool generatePermutations(uint& num_elements, vector<uint>& indices, vector<string>& names, vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, vector<vector<uint> >& permutations, vector<vector<string> >&name_order, vector<vector<GroupedWyckoffPosition> >& permutation_grouped_Wyckoff_positions);
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
  vector<string> leastFrequentAtom2(const xstructure& xstr);
  bool checkTolerance(xvector<double> d1, xmatrix<double> Q2);
  bool checkABCTolerance(xvector<double> d1, xvector<double> d2);
  bool checkAngleTolerance(xvector<double> d1, xvector<double> d2);
  xvector<double> centroid_with_PBC(const xstructure& xstr);
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, const xmatrix<double>& lattice);
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, vector<double>& weights,
                                    const xmatrix<double>& lattice);
  bool findMatch(const xstructure& xstr1, const xstructure& PROTO,vector<uint>& im1, vector<uint>& im2, vector<double>& min_dists, const int& type_match);
  void clusterize(const xstructure& xstr1, const vector<uint>& im1, vector<string>& TYPE1,
                  vector<uint>& QTA1, vector<vector<uint> >& I1);
  bool sameAtomType(const xstructure& xstr1, const xstructure& xstr2, const vector<uint>& im1,
                    const vector<uint>& im2, const int& type_match);
  bool cleanMatch(const vector<uint>& im1);
  void cellDiagonal(xstructure& xstr, vector<double>& diag_sum, vector<double>& diag_diff, const double& scale);
  void cellDiagonal(xmatrix<double>& lattice, vector<double>& diag_sum, vector<double>& diag_diff, const double& scale);
  double latticeDeviation(const vector<double>& diag_sum1,const vector<double>& diag_sum2, 
                          const vector<double>& diag_diff1,const vector<double>& diag_diff2); 
  vector<double> computeNearestNeighbors(xstructure& xstr);
  double shortestDistance(const xstructure& xstr, const uint& k);
  void coordinateDeviation(const xstructure& xstr1, const xstructure& xstr2,
                      const vector<double>& all_nn1, const vector<double>& all_nn_proto,
                      const vector<uint>& indexMatch1, const vector<uint>& indexMatch2, vector<double>& min_dists,
                      double& cd, double& fail_figure);
  double computeMisfit(const double& dev, const double& dis, const double& fail);
  void printMatch(const vector<uint>& indexMatch1, const vector<uint>& indexMatch2,
                  const xstructure& PROTO, const xstructure& xstr1, ostream& oss);
  xvector<double> bringCoordinateInCell(xvector<double>& coord);
  double checkLatticeDeviation(double& xstr1_vol, xmatrix<double>& q2,vector<double>& D1,vector<double>& F1);
  bool quadrupletPeriodic(const xmatrix<double>& quad, const xstructure& lfa_supercell, const int& i, 
                       const int& j, const int& k, const int& w);
  bool vectorPeriodic(const xvector<double>& vec, const xstructure& lfa_supercell, const int& i, 
                       const int& j);
  void threadGeneration(const uint& num_proc,xmatrix<double>& q1, xstructure& xstr2, 
                        vector<xstructure> &vprotos, xstructure &xstr1, const int& type_match, 
                        const bool& optimize_match, double& minMis, ostream& oss);
  void quadrupletSearch(const xmatrix<double>& q1, const xstructure& xstr_LFA_supercell,
                        vector<xvector<double> >& lattice_vecs, vector<vector<uint> >& ij_index);
//  void quadrupletSearch(std::atomic_bool &misfit_in_threshold_found, const string& lfa,
//                        const xmatrix<double>& q1,
//                        const xstructure& xstr_LFA_supercell,
//                        vector<xstructure>& vprotos, xstructure& xstr1, const xstructure& xstr2,
//                        const int& type_match, double& possible_minMis, ostream& oss,
//                        vector<xvector<double> >& lattice_vecs,
//                        vector<vector<uint> >& ij_index); // 1 bracket
  bool buildSimilarLattices(vector<xvector<double> >& translation_vectors, xmatrix<double>& q1, double& xstr1_vol, double& abs_det_q1, xvector<double>& abc_angles_q1, vector<xmatrix<double> >& lattices, vector<xmatrix<double> >& clattices, vector<double>& latt_devs, const bool& optimize_match);
  bool structureSearch(const string& lfa,
                        const vector<double>& all_nn1,
                        const xstructure& xstr,
                        vector<xstructure>& vprotos, xstructure& xstr1, const xstructure& xstr2,
                        const int& type_match, double& possible_minMis,
                        vector<xmatrix<double> >& lattices,
                        vector<xmatrix<double> >& clattices,
                        vector<double>& latt_devs,
                        const bool& optimize_match);

} // end namespace

#endif 

// ***************************************************************************
// AFLOW_COMPARE_STRUCTURE
// ***************************************************************************

