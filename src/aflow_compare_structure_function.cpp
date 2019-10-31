// ***************************************************************************
// 			AFLOW Compare Structure - Functions
//		David Hicks (d.hicks@duke.edu) and Carlo de Santo
// ***************************************************************************
#include<fstream>
#include<iostream>
#include<vector>
#include<string>
#include "aflow.h"
#include "aflow_pflow.h"
#include "aflow_compare_structure.h"
#include "aflow_symmetry_spacegroup.h"

#undef AFLOW_COMPARE_MULTITHREADS_ENABLE

#if GCC_VERSION >= 40400   // added two zeros
#define AFLOW_COMPARE_MULTITHREADS_ENABLE 1
#include <thread>
//DX 20190510 [OBSOLETE] #include <atomic>
#else
#warning "The multithread parts of AFLOW-COMPARE will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

// ***************************************************************************
// Prototype Class 
// ***************************************************************************
// ===== Constructor ===== //
StructurePrototype::StructurePrototype(){ 
  iomode=JSON_MODE;
  representative_structure_name="";
  representative_structure_compound="";
  representative_structure.Clear();
  representative_structure_generated=false;
  representative_structure_from="";
  number_compounds_matching_representative=0;
  number_types=0;
  elements.clear();
  stoichiometry.clear();
  number_of_atoms=0;
  unique_permutations.clear();
  Pearson="";
  space_group=0;
  grouped_Wyckoff_positions.clear();
  wyckoff_site_symmetry.clear();
  wyckoff_multiplicity.clear();
  wyckoff_letter.clear();
  duplicate_structures_names.clear();
  duplicate_structures_compounds.clear();
  duplicate_structures.clear();
  duplicate_structures_generated.clear();
  duplicate_structures_from.clear();
  number_compounds_matching_duplicate.clear();
  duplicate_comparison_logs.clear(); //DX 20190506
  family_structures_names.clear();
  family_structures.clear();
  family_structures_generated.clear();
  family_structures_from.clear();
  number_compounds_matching_family.clear();
  family_comparison_logs.clear(); //DX 20190506
  misfits.clear();
  family_misfits.clear();
  property_names.clear();
  property_units.clear();
  representative_structure_properties.clear();
  duplicate_structures_properties.clear();
  family_structures_properties.clear(); //DX 20190425
}

// ===== Free  ===== //
void StructurePrototype::free(){
}

// ===== Destructor ===== //
StructurePrototype::~StructurePrototype(){ 
  representative_structure.Clear();
  elements.clear();
  stoichiometry.clear();
  unique_permutations.clear();
  grouped_Wyckoff_positions.clear();
  wyckoff_site_symmetry.clear();
  wyckoff_multiplicity.clear();
  wyckoff_letter.clear();
  duplicate_structures_names.clear();
  duplicate_structures_compounds.clear();
  duplicate_structures.clear();
  duplicate_structures_generated.clear();
  duplicate_structures_from.clear();
  number_compounds_matching_duplicate.clear();
  duplicate_comparison_logs.clear(); //DX 20190506
  family_structures_names.clear();
  family_structures.clear();
  family_structures_generated.clear();
  family_structures_from.clear();
  number_compounds_matching_family.clear();
  family_comparison_logs.clear(); //DX 20190506
  misfits.clear();
  family_misfits.clear();
  property_names.clear();
  property_units.clear();
  representative_structure_properties.clear();
  duplicate_structures_properties.clear();
  family_structures_properties.clear(); //DX 20190425
  free();
}

// ===== Copy Constructor ===== //
StructurePrototype::StructurePrototype(const StructurePrototype& b){
  copy(b);
}

// ===== Copy Constructor Function ===== //
void StructurePrototype::copy(const StructurePrototype& b) {
  if(this != &b){
    iomode=b.iomode;
    representative_structure_name=b.representative_structure_name; 
    representative_structure_compound=b.representative_structure_compound; 
    representative_structure=b.representative_structure; 
    representative_structure_generated=b.representative_structure_generated; 
    representative_structure_from=b.representative_structure_from; 
    number_compounds_matching_representative=b.number_compounds_matching_representative;
    number_types=b.number_types;
    elements=b.elements;
    stoichiometry=b.stoichiometry;
    number_of_atoms=b.number_of_atoms;
    unique_permutations=b.unique_permutations;
    Pearson=b.Pearson;
    space_group=b.space_group;
    grouped_Wyckoff_positions=b.grouped_Wyckoff_positions;
    wyckoff_site_symmetry=b.wyckoff_site_symmetry;
    wyckoff_multiplicity=b.wyckoff_multiplicity;
    wyckoff_letter=b.wyckoff_letter;
    duplicate_structures_names=b.duplicate_structures_names;
    duplicate_structures_compounds=b.duplicate_structures_compounds;
    duplicate_structures=b.duplicate_structures;
    duplicate_structures_generated=b.duplicate_structures_generated;
    duplicate_structures_from=b.duplicate_structures_from;
    number_compounds_matching_duplicate=b.number_compounds_matching_duplicate;
    duplicate_comparison_logs=b.duplicate_comparison_logs; //DX 20190506
    family_structures_names=b.family_structures_names;
    family_structures=b.family_structures;
    family_structures_generated=b.family_structures_generated;
    family_structures_from=b.family_structures_from;
    number_compounds_matching_family=b.number_compounds_matching_family;
    family_comparison_logs=b.family_comparison_logs; //DX 20190506
    misfits=b.misfits;
    family_misfits=b.family_misfits;
    property_names=b.property_names;
    property_units=b.property_units;
    representative_structure_properties=b.representative_structure_properties;
    duplicate_structures_properties=b.duplicate_structures_properties;
    family_structures_properties=b.family_structures_properties; //DX 20190425
  }
}

// ===== Assignment Operator (operator=) ===== //
// ***************************************************************************
// StructurePrototype::operator=
// ***************************************************************************
const StructurePrototype& StructurePrototype::operator=(const StructurePrototype& b){
  if(this!=&b){
    free();
    copy(b);
  }
  return *this;
}

// ===== Output Handler ===== //
// ***************************************************************************
// StructurePrototype::operator<< 
// ***************************************************************************
ostream& operator<<(ostream& oss, const StructurePrototype& StructurePrototype){
  
  // operator<< for the StruturePrototype object (looks like a JSON)

  //if(StructurePrototype.iomode!=JSON_MODE){ //A safeguard until we construct more output schemes.
  //  StructurePrototype.iomode=JSON_MODE;
  //}

  if(StructurePrototype.iomode==JSON_MODE){
    string eendl="";
    bool roff=true; //round off
    stringstream sscontent_json;
    vector<string> vcontent_json, tmp;
    
    // representative_structure 
    sscontent_json << "\"representative_structure\":\"" << StructurePrototype.representative_structure_name << "\"" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // number_compounds_matching_representative 
    sscontent_json << "\"number_compounds_matching_representative\":" << StructurePrototype.number_compounds_matching_representative << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // number_types 
    sscontent_json << "\"number_types\":" << StructurePrototype.number_types << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // elements 
    sscontent_json << "\"elements\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.elements,"\""),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    // stoichiometry 
    sscontent_json << "\"stoichiometry\":[" << aurostd::joinWDelimiter(StructurePrototype.stoichiometry,",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // number_of_atoms 
    sscontent_json << "\"number_of_atoms\":" << StructurePrototype.number_of_atoms << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // unique_permutations
    if(StructurePrototype.unique_permutations.size()!=0){ //DX 20190425 - only print if calculated
    sscontent_json << "\"unique_permutations\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.unique_permutations,"\""),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    } //DX 20190425
    
    //DX 20190425 [OBSOLETE] // Pearson
    //DX 20190425 [OBSOLETE] sscontent_json << "\"Pearson\":\"" << StructurePrototype.Pearson << "\"" << eendl;
    //DX 20190425 [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // space_group
    sscontent_json << "\"space_group\":" << StructurePrototype.space_group << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // grouped_Wyckoff_positions
    sscontent_json << "\"grouped_Wyckoff_positions\":[";
    tmp.clear();
    for(uint i=0;i<StructurePrototype.grouped_Wyckoff_positions.size();i++){
      stringstream ss_tmp;
      ss_tmp << StructurePrototype.grouped_Wyckoff_positions[i];
      tmp.push_back(ss_tmp.str());
    }
    sscontent_json << aurostd::joinWDelimiter(tmp,",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

    //DX 20190425 [OBSOLETE] // Wyckoff_multiplicities
    //DX 20190425 [OBSOLETE] sscontent_json << "\"Wyckoff_multiplicities\":[" << aurostd::joinWDelimiter(StructurePrototype.wyckoff_multiplicity,",") << "]" << eendl;
    //DX 20190425 [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    //DX 20190425 [OBSOLETE] // Wyckoff_site_symmetries
    //DX 20190425 [OBSOLETE] sscontent_json << "\"Wyckoff_site_symmetries\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.wyckoff_site_symmetry,"\""),",") << "]" << eendl;
    //DX 20190425 [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    //DX 20190425 [OBSOLETE] // Wyckoff_letters
    //DX 20190425 [OBSOLETE] sscontent_json << "\"Wyckoff_letters\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.wyckoff_site_symmetry,"\""),",") << "]" << eendl;
    //DX 20190425 [OBSOLETE] vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // duplicate_structures
    sscontent_json << "\"duplicate_structures\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.duplicate_structures_names,"\""),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // misfits
    sscontent_json << "\"misfits\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(StructurePrototype.misfits,8,roff),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // number_compounds_matching_duplicate
    sscontent_json << "\"number_compounds_matching_duplicate\":[" << aurostd::joinWDelimiter(StructurePrototype.number_compounds_matching_duplicate,",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // family_structures
    sscontent_json << "\"family_structures\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.family_structures_names,"\""),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // family_misfits
    sscontent_json << "\"family_misfits\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(StructurePrototype.family_misfits,8,roff),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // number_compounds_matching_family
    sscontent_json << "\"number_compounds_matching_family\":[" << aurostd::joinWDelimiter(StructurePrototype.number_compounds_matching_family,",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    if(StructurePrototype.property_names.size()!=0){ //DX 20190425 - only print if requested
    // property_names
    sscontent_json << "\"property_names\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.property_names,"\""),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // property_units
    sscontent_json << "\"property_units\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.property_units,"\""),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // representative_structure_properties
    sscontent_json << "\"representative_structure_properties\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.representative_structure_properties,"\""),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // duplicate_structure_properties
      sscontent_json << "\"duplicate_structures_properties\":[";
    tmp.clear();
    //DX 20190326 - should be duplicate for(uint i=0;i<StructurePrototype.representative_structure_properties.size();i++){
    for(uint i=0;i<StructurePrototype.duplicate_structures_properties.size();i++){ // DX 20190326 - correctly changed to duplicate_sturctures_properties
      tmp.push_back(aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.duplicate_structures_properties[i],"\""),","));
    }
    sscontent_json << aurostd::joinWDelimiter(aurostd::wrapVecEntries(tmp,"[","]"),",");
    sscontent_json << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      
      //DX 20190425 - START
      // family_structure_properties
      sscontent_json << "\"family_structures_properties\":[";
      tmp.clear();
      for(uint i=0;i<StructurePrototype.family_structures_properties.size();i++){
        tmp.push_back(aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.family_structures_properties[i],"\""),","));
      }
      sscontent_json << aurostd::joinWDelimiter(aurostd::wrapVecEntries(tmp,"[","]"),",");
      sscontent_json << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
      //DX 20190425 - END
    } //DX 20190425
    
    // Put into json StructurePrototype object
    oss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
    vcontent_json.clear();
  }
  return oss;
}

// ***************************************************************************
// StructurePrototype::numberOfDuplicatess() 
// ***************************************************************************
uint StructurePrototype::numberOfDuplicates() const {

  // Count the number of duplicate structures
  
  uint number_of_duplicates = 0;
  for(uint i=0;i<misfits.size();i++){
    if(misfits[i]<=0.1 && (misfits[i]+1.0)>1e-3 ){
      number_of_duplicates+=1;
    }
  }

  return number_of_duplicates;
}

// ***************************************************************************
// StructurePrototype::numberOfComparisons() 
// ***************************************************************************
uint StructurePrototype::numberOfComparisons(){

  // Count the number of comparisons

  return duplicate_structures_names.size(); 
}

// ***************************************************************************
// StructurePrototype::isSymmetryCalculated() 
// ***************************************************************************
bool StructurePrototype::isSymmetryCalculated(){

  // Check if the space group symmetry is calculated

  if(space_group==0){return false;}

  return true; 
}

// ***************************************************************************
// StructurePrototype::calculateSymmetry() 
// ***************************************************************************
bool StructurePrototype::calculateSymmetry(){

  // Calculate the symmetry of the representative structure
  // (i.e., space group and Wyckoff positions)
  // The Pearson symbol calculation takes more time, so I do not calculate it

  Pearson = "";
  space_group = representative_structure.SpaceGroup_ITC();
  vector<GroupedWyckoffPosition> tmp_grouped_Wyckoff_positions; 
  compare::groupWyckoffPositions(representative_structure, tmp_grouped_Wyckoff_positions);
  grouped_Wyckoff_positions = tmp_grouped_Wyckoff_positions;
  return true;
}

// ***************************************************************************
// StructurePrototype::addStructurePrototypeAsDuplicate() 
// ***************************************************************************
bool StructurePrototype::addStructurePrototypeAsDuplicate(StructurePrototype& b){

  // Add representative structure as a duplicate structure in this object
  
  // add structure info
  duplicate_structures_names.push_back(b.representative_structure_name);
  duplicate_structures_generated.push_back(b.representative_structure_generated);
 
  // only add xstructure if it has been generated
  if(b.representative_structure_generated){
    duplicate_structures.push_back(b.representative_structure);
    duplicate_structures_compounds.push_back(compare::getCompoundName(b.representative_structure)); //DX 20190111 - added compound, e.g., Ag1Br2
  }
  else if(!b.representative_structure_compound.empty()){
    duplicate_structures_compounds.push_back(b.representative_structure_compound); //DX 20190111 - added compound, e.g., Ag1Br2
  }
  else {
    duplicate_structures_compounds.push_back(""); //DX 20190111 - added compound, e.g., Ag1Br2
  }
  duplicate_structures_from.push_back(b.representative_structure_from);
  number_compounds_matching_duplicate.push_back(b.number_compounds_matching_representative); //DX 20190228 - added duplicate compound count

  // add structure properties
  if(b.representative_structure_properties.size()!=0){
    duplicate_structures_properties.push_back(b.representative_structure_properties); //DX 20181218 - added property_values
  }

  // signifies comparison has not been performed yet
  misfits.push_back(-1.0);

  return true;
}

// ***************************************************************************
// StructurePrototype::copyPrototypeInformation() 
// ***************************************************************************
bool StructurePrototype::copyPrototypeInformation(StructurePrototype& b){

  // Copy prototype information to new StructurePrototype object
  // (i.e., stoichiometry, number of types, Pearson symbol, space group,
  // and Wyckoff positions)

  //elements=comparison_schemes[i].elements;
  stoichiometry=b.stoichiometry;
  number_types=b.stoichiometry.size();
  Pearson=b.Pearson;
  space_group=b.space_group;
  grouped_Wyckoff_positions=b.grouped_Wyckoff_positions;
  return true;
}

// ***************************************************************************
// StructurePrototype::putDuplicateAsRepresentative() 
// ***************************************************************************
bool StructurePrototype::putDuplicateAsRepresentative(StructurePrototype& b, uint& index){

  // Make structure in duplicate the new representative structure
  // if structure (b.duplicate_structure[index]) does not match with the 
  // original representative structure (i.e., misfit>1.0),
  // make a new object with this structure as the representative 
  
  // load structure info
  representative_structure_name=b.duplicate_structures_names[index];
  // number_types=b.duplicate_structures[index].num_each_type.size();
  // elements??? don't want to propogate elements which do not apply
  representative_structure_compound=b.duplicate_structures_compounds[index]; //DX 20190111 - added compound, e.g., Ag1Br2
  representative_structure_generated=b.duplicate_structures_generated[index];
  representative_structure_from=b.duplicate_structures_from[index];
  number_compounds_matching_representative=b.number_compounds_matching_duplicate[index]; //DX 20190228 - added duplicate compound count
  if(b.duplicate_structures_generated[index]){
    number_of_atoms=b.duplicate_structures[index].atoms.size();
    representative_structure=b.duplicate_structures[index];
  }
  
  // load property info
  if(b.property_names.size()!=0){
    property_names=b.property_names;
    property_units=b.property_units;
    representative_structure_properties=b.duplicate_structures_properties[index];
  }
  return true;
}

// ***************************************************************************
// StructurePrototype::copyDuplicate() 
// ***************************************************************************
bool StructurePrototype::copyDuplicate(StructurePrototype& b, uint& index){

  // Copy structure at index as a potential duplicate in this object 

  duplicate_structures_names.push_back(b.duplicate_structures_names[index]);
  duplicate_structures_compounds.push_back(b.duplicate_structures_compounds[index]);
  duplicate_structures_generated.push_back(b.duplicate_structures_generated[index]);
  duplicate_structures_from.push_back(b.duplicate_structures_from[index]);
  number_compounds_matching_duplicate.push_back(b.number_compounds_matching_duplicate[index]); //DX 20190228 - added duplicate compound count
  if(b.duplicate_structures_generated[index]){
    duplicate_structures.push_back(b.duplicate_structures[index]);
  }
  // signifies comparison has not been performed yet
  misfits.push_back(-1.0);
  
  // load property info
  if(property_names.size()!=0){
    duplicate_structures_properties.push_back(b.duplicate_structures_properties[index]);
  }

  return true;
}

// ***************************************************************************
// StructurePrototype::removeNonDuplicate() 
// ***************************************************************************
bool StructurePrototype::removeNonDuplicate(uint& index){

  // If structure (b.duplicate_structure[index]) does not match with the 
  // original representative structure (i.e., misfit>1.0),
  // remove from scheme

  // remove structure information
  duplicate_structures_names.erase(duplicate_structures_names.begin()+index);
  duplicate_structures_compounds.erase(duplicate_structures_compounds.begin()+index); //DX 20190111 - added compound, e.g., Ag1Br2
  // DX - may need to be careful here.  If we have a mix of generated and non-generated, we may have difficulties
  if(duplicate_structures_generated[index]){
    duplicate_structures.erase(duplicate_structures.begin()+index);
  }
  duplicate_structures_generated.erase(duplicate_structures_generated.begin()+index);
  duplicate_structures_from.erase(duplicate_structures_from.begin()+index);
  number_compounds_matching_duplicate.erase(number_compounds_matching_duplicate.begin()+index);

  // remove misfit
  misfits.erase(misfits.begin()+index);

  //DX 20190504 - START
  // remove comparison log
  if(duplicate_comparison_logs.size()!=0){
    duplicate_comparison_logs.erase(duplicate_comparison_logs.begin()+index);
  }
  //DX 20190504 - END

  // remove properties
  if(property_names.size()!=0){
    duplicate_structures_properties.erase(duplicate_structures_properties.begin()+index);
  }
  return true;
} 

// ***************************************************************************
// StructurePrototype::clearDuplicateInformation() 
// ***************************************************************************
bool StructurePrototype::removeDuplicates(bool remove_duplicate_count){

  // Remove duplicate structure information from object
  
  // remove structure information
  duplicate_structures_names.clear();
  duplicate_structures_compounds.clear(); //DX 20190111 - added compound, e.g., Ag1Br2
  duplicate_structures.clear();
  duplicate_structures_generated.clear();
  duplicate_structures_from.clear();

  // may want to keep this information
  if(remove_duplicate_count){number_compounds_matching_duplicate.clear();}

  // remove comparison logs 
  duplicate_comparison_logs.clear(); //DX 20190506

  // remove misfit
  misfits.clear();

  // remove properties
  duplicate_structures_properties.clear();

  return true;
}

// ***************************************************************************
// END:: Prototype Class
// ***************************************************************************

// ***************************************************************************
// START:: GroupedWyckoffPosition Class
// ***************************************************************************
// ===== Constructor ===== //
GroupedWyckoffPosition::GroupedWyckoffPosition(){
  type = 0;
  element = "";
  site_symmetries.clear();
  multiplicities.clear();
  letters.clear();
}

// ===== Free  ===== //
void GroupedWyckoffPosition::free(){
}

// ===== Destructor  ===== //
GroupedWyckoffPosition::~GroupedWyckoffPosition(){ 
  site_symmetries.clear();
  multiplicities.clear();
  letters.clear();
}

// ===== Copy Constructor ===== //
GroupedWyckoffPosition::GroupedWyckoffPosition(const GroupedWyckoffPosition& b){
  copy(b);
}

// ===== Copy Constructor Function ===== //
void GroupedWyckoffPosition::copy(const GroupedWyckoffPosition& b) {
  if(this != &b){
    type=b.type;
    element=b.element;
    site_symmetries=b.site_symmetries;
    multiplicities=b.multiplicities;
    letters=b.letters;
  }
}

// ===== Assignment Operator (operator=) ===== //
const GroupedWyckoffPosition& GroupedWyckoffPosition::operator=(const GroupedWyckoffPosition& b){
  if(this!=&b){
    free();
    copy(b);
  }
  return *this;
}

// ===== Output Handler ===== //
ostream& operator<<(ostream& oss, const GroupedWyckoffPosition& GroupedWyckoffPosition){

  // operator<< for the GroupedWyckoffPosition object (looks like a JSON)

  //if(StructurePrototype.iomode!=JSON_MODE){ //A safeguard until we construct more output schemes.
  //  StructurePrototype.iomode=JSON_MODE;
  //}

  string eendl="";
  stringstream sscontent_json;
  vector<string> vcontent_json;
  
  // type 
  sscontent_json << "\"type\":" << GroupedWyckoffPosition.type << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // element
  sscontent_json << "\"element\":\"" << GroupedWyckoffPosition.element << "\"" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // site_symmetries
  sscontent_json << "\"site_symmetries\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(GroupedWyckoffPosition.site_symmetries,"\""),",") << "]" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // multiplicities
  sscontent_json << "\"multiplicities\":[" << aurostd::joinWDelimiter(GroupedWyckoffPosition.multiplicities,",") << "]" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // letters
  sscontent_json << "\"letters\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(GroupedWyckoffPosition.letters,"\""),",") << "]" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // Put into json GroupedWyckoffPosition object
  oss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
  vcontent_json.clear();
  
  return oss;
}

// ===== Less than operator<, for sortin ===== //
bool GroupedWyckoffPosition::operator<(const GroupedWyckoffPosition& b) const{       // operator<

  // allows sorting vector by element name

  if(element<b.element){
    return true;
  }
  return false;
}

// ***************************************************************************
// END:: GroupedWyckoffPosition Class
// ***************************************************************************

// ***************************************************************************


// ***************************************************************************
// *                                                                         *
// *                             FUNCTIONS                                   *
// *                                                                         *
// ***************************************************************************

// ***************************************************************************
// loadStructuresFromDirectory() 
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> loadStructuresFromDirectory(string& directory, bool& same_species, ofstream& FileMESSAGE){ //DX 20190319 - added FileMESSAGE

    // load all structures from a directory into a vector of StructurePrototype 
    // objects
    string function_name = "compare:loadStructuresFromDirectory()";
    
    bool LDEBUG=(false || XHOST.DEBUG);
    ostream& logstream = cout;
    stringstream message;
    //DX 20190319 [OBSOLETE] ofstream FileMESSAGE;
    
    vector<StructurePrototype> all_structures;
    vector<string> vfiles;
    aurostd::DirectoryLS(directory, vfiles);
    std::sort(vfiles.begin(),vfiles.end()); //CO 180830

    message << "Loading " << vfiles.size() << " files in directory ... ";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    for(uint i=0; i<vfiles.size(); i++){
      if(LDEBUG) {cerr << "compare:: " << i << "/" << vfiles.size() << " " << vfiles[i] << endl;}
      if(vfiles[i].find("comparison_output.json") != std::string::npos || 
         vfiles[i].find("comparison_output.out") != std::string::npos || 
         vfiles[i].find("comparison_no_duplicate_compounds_output.json") != std::string::npos || 
         vfiles[i].find("comparison_no_duplicate_compounds_output.out") != std::string::npos || 
         vfiles[i].find("duplicate_compounds_output.json") != std::string::npos || 
         vfiles[i].find("duplicate_compounds_output.out") != std::string::npos || 
         vfiles[i].find("nohup.out") != std::string::npos){
        vfiles.erase(vfiles.begin()+i);
        i--;
      }
      else {
        StructurePrototype structure_tmp;
        stringstream sss1;
        aurostd::efile2stringstream(directory+"/"+vfiles[i],sss1);
        xstructure xstr1(sss1);
        structure_tmp.representative_structure = xstr1;
        structure_tmp.representative_structure_name = directory+"/"+vfiles[i];
        structure_tmp.stoichiometry = compare::getStoichiometry(xstr1,same_species);
        structure_tmp.elements = compare::getElements(xstr1);
        structure_tmp.number_of_atoms = xstr1.atoms.size(); //DX 20190425
        structure_tmp.number_types = xstr1.num_each_type.size(); //DX 20190425
	structure_tmp.representative_structure_compound = compare::getCompoundName(xstr1); //remove ones is true  //DX 20190311 //DX 20190313 - use xstr1
        // update xstructure species
        if(structure_tmp.representative_structure.species.size()==0){
          deque<string> deque_species; for(uint j=0;j<structure_tmp.elements.size();j++){deque_species.push_back(structure_tmp.elements[j]);}
          structure_tmp.representative_structure.SetSpecies(deque_species);
          structure_tmp.representative_structure.SpeciesPutAlphabetic();
        }
        // check if fake names for same species comparison
        if(structure_tmp.representative_structure.species[0]=="A" && same_species){
          message << "Atomic species are missing for " << structure_tmp.representative_structure_name << " cannot perform material comparison; skipping structure.";     
          pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);
          continue;
        }
        structure_tmp.representative_structure_generated = true; 
        structure_tmp.representative_structure_from = "file"; 
        if(LDEBUG) {
          cerr << "compare::loadStructureFromDirectory() Found structure: " << directory+"/"+vfiles[i] << endl;
        }
        all_structures.push_back(structure_tmp);
      }
    }
    return all_structures; 
  }
}

// ***************************************************************************
// loadStructuresFromFile() 
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> loadStructuresFromFile(string& filename, bool& same_species, ofstream& FileMESSAGE){ //DX 20190319 - added FileMESSAGE

    // load all structures from a file into a vector of StructurePrototype object
    // useful for reading in aflow.in relaxation steps or pocc structures

    string function_name = "compare:loadStructuresFromFile()";
    
    bool LDEBUG=(false || XHOST.DEBUG);
    ostream& logstream = cout;
    stringstream message;
    //DX [OBSOLETE] ofstream FileMESSAGE;
    
    vector<StructurePrototype> all_structures;
   
    // ---------------------------------------------------------------------------
    // file to stringstream
    stringstream input_file;
    aurostd::efile2stringstream(filename, input_file);

    // ---------------------------------------------------------------------------
    // tokenize stringstream by newline
    vector<string> lines;
    aurostd::string2tokens(input_file.str(),lines,"\n");
    
    // ---------------------------------------------------------------------------
    // structure delimiters 
    string START="[VASP_POSCAR_MODE_EXPLICIT]START";
    string STOP="[VASP_POSCAR_MODE_EXPLICIT]STOP";

    // ---------------------------------------------------------------------------
    // used to find the total number of structures
    vector<string> start_string; 
    aurostd::substring2strings(input_file.str(),start_string,START);
    
    message << "Loading " << start_string.size() << " structures in file ... ";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);

    bool structure_lines = false;
    uint structure_count = 0;
    stringstream geometry;geometry.clear();geometry.str(std::string());
    for(uint i=0;i<lines.size();i++){
      if(aurostd::substring2bool(lines[i],START)){ 
        geometry.clear();geometry.str("");
        structure_lines = true; 
        structure_count+=1;
      }
      else if(structure_lines && !aurostd::substring2bool(lines[i],STOP)){
        geometry << lines[i] << endl;
      }
      else if(aurostd::substring2bool(lines[i],STOP)){ 
        structure_lines = false; 
        if(LDEBUG) {cerr << "compare:: loading " << structure_count << "/" << start_string.size() << endl;}
        StructurePrototype structure_tmp;
        xstructure xstr1(geometry);
        structure_tmp.representative_structure = xstr1;
        stringstream designation; designation << "file structure # " << structure_count << "/" << start_string.size();
        structure_tmp.representative_structure_name = designation.str();
        structure_tmp.stoichiometry = compare::getStoichiometry(xstr1,same_species);
        structure_tmp.elements = compare::getElements(xstr1);
        structure_tmp.number_of_atoms = xstr1.atoms.size(); //DX 20190425
        structure_tmp.number_types = xstr1.num_each_type.size(); //DX 20190425
	structure_tmp.representative_structure_compound = compare::getCompoundName(xstr1); //remove ones is true  //DX 20190311 //DX 20190313 - use xstr
        // update xstructure species
        if(structure_tmp.representative_structure.species.size()==0){
          deque<string> deque_species; for(uint j=0;j<structure_tmp.elements.size();j++){deque_species.push_back(structure_tmp.elements[j]);}
          structure_tmp.representative_structure.SetSpecies(deque_species);
          structure_tmp.representative_structure.SpeciesPutAlphabetic();
        }
        // check if fake names for same species comparison
        if(structure_tmp.representative_structure.species[0]=="A" && same_species){
          message << "Atomic species are missing for " << structure_tmp.representative_structure_name << " cannot perform material comparison; skipping strucutre.";     
          pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);
          continue;
        }
        structure_tmp.representative_structure_generated = true; 
        structure_tmp.representative_structure_from = input_file.str(); 
        if(LDEBUG) {
          cerr << "compare::loadStructureFromFile(): loaded structure " << i << endl;
        }
        all_structures.push_back(structure_tmp);
      }
    }
    return all_structures; 
  }
}

//DX 20190424 - START
// ***************************************************************************
// loadStructuresFromStructureList() 
// ***************************************************************************
namespace compare {
  vector<StructurePrototype> loadStructuresFromStructureList(vector<string>& filenames, bool& same_species, ofstream& FileMESSAGE){

    // load all structures from a vector of filenames into a vector of StructurePrototype object

    string function_name = "compare:loadStructuresFromStructureList()";
    
    bool LDEBUG=(false || XHOST.DEBUG);
    ostream& logstream = cout;
    stringstream message;
    //DX [OBSOLETE] ofstream FileMESSAGE;
    
    vector<StructurePrototype> all_structures;
  
    // ---------------------------------------------------------------------------
    // read in files
    for(uint i=0;i<filenames.size();i++){
      StructurePrototype structure_tmp;
      if(!aurostd::FileExist(filenames[i])){
        cerr << function_name << ": ERROR " << filenames[i] << " file not found" << endl;
        exit(1);
      }
      stringstream sss;
      aurostd::efile2stringstream(filenames[i],sss);
      xstructure xstr(sss);  
      structure_tmp.representative_structure = xstr;
      structure_tmp.representative_structure_name = filenames[i];
      structure_tmp.stoichiometry = compare::getStoichiometry(xstr,same_species);
      structure_tmp.elements = compare::getElements(xstr);
      structure_tmp.number_of_atoms = xstr.atoms.size(); //DX 20190425
      structure_tmp.number_types = xstr.num_each_type.size(); //DX 20190425
	    structure_tmp.representative_structure_compound = compare::getCompoundName(xstr); //remove ones is true  //DX 20190311 //DX 20190313 - use xstr
      // update xstructure species
      if(structure_tmp.representative_structure.species.size()==0){
        deque<string> deque_species; for(uint j=0;j<structure_tmp.elements.size();j++){deque_species.push_back(structure_tmp.elements[j]);}
        structure_tmp.representative_structure.SetSpecies(deque_species);
        structure_tmp.representative_structure.SpeciesPutAlphabetic();
      }
      // check if fake names for same species comparison
      if(structure_tmp.representative_structure.species[0]=="A" && same_species){
        message << "Atomic species are missing for " << structure_tmp.representative_structure_name << " cannot perform material comparison; skipping strucutre.";     
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);
        continue;
      }
      structure_tmp.representative_structure_generated = true; 
      structure_tmp.representative_structure_from = "file"; 
      if(LDEBUG) {
        cerr << function_name << ": loaded structure " << i << endl;
      }
      all_structures.push_back(structure_tmp);
    }
    return all_structures;
  }
}
//DX 20190424 - END

// ***************************************************************************
// generateStructure 
// ***************************************************************************
namespace compare {
  bool generateStructure(string& structure_name, string& structure_from, xstructure& structure, ostream& oss){

    // generate the xstructure object, having this separate function allows us to load structures in 
    // a threaded environment

    // there are three modes of generating xstructures
    // 1) AFLOW prototypes
    // 2) AURL
    // 3) input (cin)

    string function_name = "compare::generateStructure()";
    ofstream FileMESSAGE;
    vector<string> tokens;

    if(structure_from=="aflow_prototypes"){
      // htqc or anrl
      structure = aflowlib::PrototypeLibraries(cout,structure_name,"",2); 
    }
    else if(structure_from=="aurl"){
      aflowlib::_aflowlib_entry entry; entry.aurl = structure_name; 
      //DX 20190326 - need to put url path, i.e., structure name, [OBSOLETE] if(!pflow::loadXstructures(entry,FileMESSAGE,oss,true,structure_name,true)){ cerr << function_name << "WARNING::Could not load structure via aurl..." << endl; return false;}
      if(!pflow::loadXstructures(entry,FileMESSAGE,oss,true,structure_name,true)){ cerr << function_name << "WARNING::Could not load structure via aurl..." << endl; return false;} //DX 20190326
      if(entry.vstr.size()==1){
        structure = entry.vstr[0];
      }
      else {
        cerr << function_name << "::WARNING: More structures loaded than anticipated." << endl;
        return false;
      }
    }
    else if(structure_from=="file"){
      stringstream sss;
      aurostd::efile2stringstream(structure_name,sss);
      xstructure xstr(sss);
      structure = xstr;
    }
    else if(aurostd::substring2bool(structure_name,"file structure ")){
      aurostd::string2tokens(structure_name,tokens,"#");
      string structure_designation = aurostd::RemoveWhiteSpaces(tokens[1]);
      uint structure_number=0; uint number_of_structures=0;
      tokens.clear(); aurostd::string2tokens(structure_designation,tokens,"/");
      structure_number=aurostd::string2utype<uint>(tokens[0]);
      number_of_structures=aurostd::string2utype<uint>(tokens[1]);
    
      // ---------------------------------------------------------------------------
      // tokenize stringstream by newline
      vector<string> lines;
      aurostd::string2tokens(structure_from,lines,"\n");
    
      // ---------------------------------------------------------------------------
      // structure delimiters 
      string START="[VASP_POSCAR_MODE_EXPLICIT]START";
      string STOP="[VASP_POSCAR_MODE_EXPLICIT]STOP";

      bool structure_lines = false;
      uint structure_count = 0;
      stringstream geometry;geometry.clear();geometry.str(std::string());
      for(uint i=0;i<lines.size();i++){
        if(aurostd::substring2bool(lines[i],START)){ 
          stringstream geometry;geometry.clear();geometry.str(std::string());
          structure_lines = true; 
          structure_count+=1;
        }
        else if(structure_lines && structure_count==structure_number && !aurostd::substring2bool(lines[i],STOP)){
          geometry << lines[i] << endl;
        }
        else if(aurostd::substring2bool(lines[i],STOP) && structure_lines){ 
          structure_lines = false; 
          xstructure xstr(geometry);
          structure = xstr;
          break;
        }
      }
    }
    else if(structure_name=="input geometry"){
      stringstream sss; sss << structure_from;
      xstructure xstr(sss);
      structure = xstr;
    }
    else {
      cerr << function_name << "::WARNING: Structure location (from=" << structure_from << ") is not specified correctly for " << structure_name << " (i.e., input, aflow_prototype, aurl, etc.)." << endl;
      return false;
    }
    return true;
  }
}
  
/*
// ***************************************************************************
// SVD Decomposition 
// ***************************************************************************
namespace compare{
  bool SVD(xmatrix<double> A){
    xmatrix<double> ATA = aurostd::trasp(A)*A;
    //to bidiagonal
    xmatrix<double> test = ATA;
    xmatrix<double> Q = pflow::generalHouseHolderQRDecomposition(test);
    //cerr << "A: " << A << endl;
    //cerr << "ATA: " << ATA << endl;
    //cerr << "Q: " << Q << endl;
    //cerr << "R: " << test << endl;
    return true;
  }
}
*/

// ***************************************************************************
// Find ICSD name - Find ICSD name 
// ***************************************************************************
namespace compare{
  string findICSDName(string& name){

    // Find ICSD substring within path name
    // In order for this to work, the following ICSD name format must be 
    // present in the string (e.g., Ag1_ICSD_#####)

    string ICSD_substring = "";
    bool ICSD_substring_found = false;
    if(aurostd::substring2bool(name,"/")){
      vector<string> tokens;
      aurostd::string2tokens(name,tokens,"/");
      for(uint i=0;i<tokens.size();i++){
        if(aurostd::substring2bool(tokens[i],"_ICSD_")){
          ICSD_substring = tokens[i];
          ICSD_substring_found = true;
        }
      }
    }
    else {
      cerr << aurostd::substring2bool(name,"_ICSD_") << endl;
      if(aurostd::substring2bool(name,"_ICSD_")){ 
        ICSD_substring = name;
        ICSD_substring_found = true;
      }
    }
    if(!ICSD_substring_found){
      cerr << "compare::findICSDName: WARNING: Could not find ICSD substring in name.  representative prototype will not necessarily be the minimum ICSD number." << endl; 
      cerr << "string: " << name << endl;
    }
    return ICSD_substring;
  }
}

// ***************************************************************************
// Find Minimum ICSD Entry - Find the minimum ICSD number in a set of names 
// ***************************************************************************
namespace compare{
  string findMinimumICSDEntry(vector<string>& ICSD_entries){

    // Identify the structure with the minimum ICSD number (to use 
    // as the representative structure in the comparison)

    string min_ICSD = "";
    int min_num = 1e9;
    for(uint i=0;i<ICSD_entries.size();i++){
      vector<string> tokens;
      aurostd::string2tokens(ICSD_entries[i],tokens,"_"); 
      string num_string = tokens[tokens.size()-1];
      for(uint j=0;j<num_string.size();j++){
        if(isalpha(num_string[j])){
          num_string.erase(num_string.begin()+j);
        }
      }
      int num = aurostd::string2utype<int>(tokens[tokens.size()-1]); // ICSD number is aways of the form (Ag1_ICSD_######)
      if(num < min_num){
        min_ICSD = ICSD_entries[i];
        min_num = num;
      }
    }
    return min_ICSD;
  }
}

// ***************************************************************************
// groupSameRatios
// ***************************************************************************
namespace compare{
  bool groupSameRatios(vector<int>& stoich, vector<int>& unique_stoich, vector<vector<int> >& type_index){

    // map two sets of stoichiometries with one another 

    for(uint i=0;i<stoich.size();i++){
      bool stoich_stored = false;      
      for(uint j=0;j<unique_stoich.size();j++){     
        if(stoich[i] == unique_stoich[j]){
          stoich_stored = true;
          type_index[j].push_back(i);
          break;
        }
      }
      if(!stoich_stored){
        unique_stoich.push_back(stoich[i]);
        vector<int> tmp; tmp.push_back(i);
        type_index.push_back(tmp);
      }
    }
    return true;
  }
}

// ***************************************************************************
// generatePermutations 
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> comparePermutations(StructurePrototype& structure, uint& num_proc, bool& optimize_match, ostream& oss, ofstream& FileMESSAGE){ //DX 20190319 - added FileMESSAGE

    bool same_species = true; // permutation comparisons must compare the same species 
    bool scale_volume = true; // permutations are generated from the same structure, so they will have the same volume anyway
    //bool optimize_match = false;
    bool ignore_symmetry = false; // duplicate permutations should have the same symmetry
    bool ignore_Wyckoff = false; // duplicate permutations should have the same Wyckoff positions
    bool single_comparison_round = false; // compare all permutations until matched or exhausted all comparisons
    bool clean_unmatched = true; // remove unmatched structures from object //DX 20190504 
    bool quiet = true; //true 

    // ===== Call main permutation comparison function ===== //
    vector<StructurePrototype> permutation_comparisons;

    vector<uint> stoichiometry = compare::getStoichiometry(structure.representative_structure,true);

    // Calculate symmetry
    if(structure.space_group ==0){
      structure.calculateSymmetry();
    }

    //cerr << "generating all permutations" << endl;
    // generate all permuations structures
    vector<StructurePrototype> permutation_structures = compare::generatePermutationStructures(structure);

    //cerr << "store naming" << endl;
    vector<vector<string> > name_order;
    for(uint i=0;i<permutation_structures.size();i++){
      vector<string> vtmp; 
      for(uint j=0;j<permutation_structures[i].representative_structure_name.size();j++){
        stringstream ss_tmp; ss_tmp << permutation_structures[i].representative_structure_name[j];
        vtmp.push_back(ss_tmp.str());
      }
      name_order.push_back(vtmp);
    }

    // group comparable permutations
    permutation_comparisons = compare::groupStructurePrototypes(permutation_structures, same_species, ignore_symmetry, ignore_Wyckoff);
    
    // ensure the representative stucture is an even permutation
    compare::makeRepresentativeEvenPermutation(permutation_comparisons, name_order);

    // compare permutations
    vector<StructurePrototype> final_permutations = compare::runComparisonScheme(num_proc, permutation_comparisons, same_species, scale_volume, optimize_match, single_comparison_round, clean_unmatched, false,oss,FileMESSAGE,quiet); //DX 20190319 - added FileMESSAGE 
    if(!compare::checkNumberOfGroupings(final_permutations, name_order.size())){
      cerr << "compare::comparePermutations():: ERROR: Compared groupings of permutations do not follow number theory. Contact David Hicks (d.hicks@duke.edu). Exiting." << endl;
      exit(1);
    }

    //DEBUG cerr << "final_permuations:" << final_permutations.size() << endl;
    //DEBUG for(uint j=0;j<final_permutations.size();j++){
    //DEBUG   cerr << "j: " << j << endl;
    //DEBUG   cerr << final_permutations[j] << endl;
    //DEBUG }
    return final_permutations;
  }
}

// ***************************************************************************
// generatePermutationStructures
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> generatePermutationStructures(StructurePrototype& structure){
   
    vector<StructurePrototype> permutation_structures;
 
    vector<string> names = fakeElements(structure.stoichiometry.size());
    vector<uint> indices = structure.stoichiometry;
    vector<vector<string> > name_order;  
    uint num_elements = structure.stoichiometry.size();
    
    bool is_symmetry_calculated = structure.isSymmetryCalculated(); //DX 20190508 - check if Wyckoff positions determined
 
    // Permutation algorithm based on Heap's algorithm (https://en.wikipedia.org/wiki/Heap%27s_algorithm)
    vector<uint> new_indices;
    for(uint i=0;i<num_elements;i++){new_indices.push_back(0);}
 
    name_order.push_back(names); 
  
    //DX 20190508 [OBSOLETE] vector<GroupedWyckoffPosition> grouped_Wyckoff_positions = structure.grouped_Wyckoff_positions;
    //DX 20190508 - add option to generate permutation and ignore Wyckoff - START
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions; 
    vector<vector<GroupedWyckoffPosition> > permutation_grouped_Wyckoff_positions;

    if(is_symmetry_calculated){
      grouped_Wyckoff_positions = structure.grouped_Wyckoff_positions;
    for(uint i=0;i<grouped_Wyckoff_positions.size();i++){grouped_Wyckoff_positions[i].element=names[i];}
    //DX 20190508 [OBSOLETE] vector<vector<GroupedWyckoffPosition> > permutation_grouped_Wyckoff_positions;
    permutation_grouped_Wyckoff_positions.push_back(grouped_Wyckoff_positions);
    }
    //DX 20190508 - add option to generate permutation and ignore Wyckoff - END

    uint i=0;
    while(i<num_elements){
      if(new_indices[i] < i){
        //LDEBUBcerr << "orig: ";
        //for(uint n=0;n<indices.size();n++){cerr << indices[n] << " ";}
        if(i%2==0){
          //LDEBUBcerr << "even: swapping " << indices[0] << " and " << indices[i] << endl;
          int swap1 = indices[0];
          int swap2 = indices[i];
          indices[0]=swap2; indices[i]=swap1;
          string swap_name1 = names[0];
          string swap_name2 = names[i];
          //GroupedWyckoffPosition swap_position1 = grouped_Wyckoff_positions[0]; swap_position1.element = names[i];
          //GroupedWyckoffPosition swap_position2 = grouped_Wyckoff_positions[i]; swap_position2.element = names[0];
          names[0]=swap_name2; names[i]=swap_name1;
          if(is_symmetry_calculated){ //DX 20190508
          grouped_Wyckoff_positions[0].element=swap_name2; grouped_Wyckoff_positions[i].element=swap_name1;
        }
        }
        else {
          //LDEBUBcerr << "odd: swapping " << indices[new_indices[i]] << " and " << indices[i] << endl;
          int swap1 = indices[new_indices[i]];
          int swap2 = indices[i];
          indices[new_indices[i]]=swap2; indices[i]=swap1;
          string swap_name1 = names[new_indices[i]];
          string swap_name2 = names[i];
          //GroupedWyckoffPosition swap_position1 = grouped_Wyckoff_positions[new_indices[i]]; swap_position1.element = names[i];
          //GroupedWyckoffPosition swap_position2 = grouped_Wyckoff_positions[i]; swap_position2.element = names[new_indices[i]];
          names[new_indices[i]]=swap_name2; names[i]=swap_name1;
          if(is_symmetry_calculated){ //DX 20190508
          grouped_Wyckoff_positions[new_indices[i]].element=swap_name2; grouped_Wyckoff_positions[i].element=swap_name1;
        }
        }
        //LDEBUBcerr << "storing: ";
        //for(uint n=0;n<indices.size();n++){cerr << indices[n] << " ";}
        //cerr << endl;
        //permutations.push_back(indices);
        name_order.push_back(names);
        //DX 20190508 - check for Wyckoff - START
        if(is_symmetry_calculated){
        vector<GroupedWyckoffPosition> permuted_grouped_Wyckoff_positions = grouped_Wyckoff_positions;
        std::sort(permuted_grouped_Wyckoff_positions.begin(), permuted_grouped_Wyckoff_positions.end()); //sort Wyckoff positions alphabetically by element
        permutation_grouped_Wyckoff_positions.push_back(permuted_grouped_Wyckoff_positions);
        }
        else {
          permutation_grouped_Wyckoff_positions.push_back(grouped_Wyckoff_positions); // push back empty
        }
        //DX 20190508 - check for Wyckoff - END
        new_indices[i]++;
        i=0;
      } 
      else {
        //LDEBUBcerr << "moving on" << endl;
        new_indices[i]=0;
        i++;
      }   
    }

    // create permuted structure    
    for(uint i=0;i<name_order.size();i++){
      xstructure xstr_tmp = structure.representative_structure;
      deque<string> species; 
      for(uint j=0;j<name_order[i].size();j++){species.push_back(name_order[i][j]);}
      xstr_tmp.SetSpecies(species);
      xstr_tmp.SpeciesPutAlphabetic();
      deque<int> sizes = SYM::arrange_atoms(xstr_tmp.atoms);
      //LDEBUGfor(uint j=0;j<sizes.size();j++){cerr << "sizes[j]: " << sizes[j] << endl;}
      xstr_tmp = pflow::SetNumEachType(xstr_tmp, sizes);
      //if (xstr_out.num_each_type.size() != names.size()){
      //  xstr_out = pflow::SetAllAtomNames(xstr_out, in_names);
      //}

      StructurePrototype tmp;
      tmp.representative_structure = xstr_tmp;
      tmp.representative_structure_name = aurostd::joinWDelimiter(species,"");
      tmp.representative_structure_generated = true;
      tmp.representative_structure_from = "permutation";
      tmp.copyPrototypeInformation(structure);
      tmp.grouped_Wyckoff_positions = permutation_grouped_Wyckoff_positions[i];
      permutation_structures.push_back(tmp);
    }
   
    return permutation_structures;
  }
}

//DX 20190508 - added permutation string function - START
// ***************************************************************************
// generatePermutationString
// ***************************************************************************
namespace compare{
  vector<string> generatePermutationString(vector<uint>& stoichiometry){
   
    vector<StructurePrototype> permutation_structures;
 
    vector<string> names = fakeElements(stoichiometry.size());
    vector<uint> indices = stoichiometry;
    vector<vector<string> > name_order;  
    uint num_elements = stoichiometry.size();
    
    // Permutation algorithm based on Heap's algorithm (https://en.wikipedia.org/wiki/Heap%27s_algorithm)
    vector<uint> new_indices;
    for(uint i=0;i<num_elements;i++){new_indices.push_back(0);}

    name_order.push_back(names); 
 
    uint i=0;
    while(i<num_elements){
      if(new_indices[i] < i){
        if(i%2==0){
          int swap1 = indices[0];
          int swap2 = indices[i];
          indices[0]=swap2; indices[i]=swap1;
          string swap_name1 = names[0];
          string swap_name2 = names[i];
          names[0]=swap_name2; names[i]=swap_name1;
        }
        else {
          int swap1 = indices[new_indices[i]];
          int swap2 = indices[i];
          indices[new_indices[i]]=swap2; indices[i]=swap1;
          string swap_name1 = names[new_indices[i]];
          string swap_name2 = names[i];
          names[new_indices[i]]=swap_name2; names[i]=swap_name1;
        }
        name_order.push_back(names);
        new_indices[i]++;
        i=0;
      } 
      else {
        new_indices[i]=0;
        i++;
      }   
    }

    vector<string> permuted_string;
    for(uint i=0;i<name_order.size();i++){ permuted_string.push_back(aurostd::joinWDelimiter(name_order[i],"")); }

    return permuted_string; 
  }
}
//DX 20190508 - added permutation string function - END

// ***************************************************************************
// generatePermutations 
// ***************************************************************************
namespace compare{
  bool generatePermutations(uint& num_elements, vector<uint>& indices, vector<string>& names, vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions, vector<vector<uint> >& permutations, vector<vector<string> >&name_order, vector<vector<GroupedWyckoffPosition> >& permutation_grouped_Wyckoff_positions){

    // Permutation algorithm based on Heap's algorithm (https://en.wikipedia.org/wiki/Heap%27s_algorithm)

    //vector<vector<int> > permutations;
    vector<uint> new_indices;
    for(uint i=0;i<indices.size();i++){new_indices.push_back(0);}
   
    permutations.push_back(indices);
    name_order.push_back(names);    
    for(uint i=0;i<grouped_Wyckoff_positions.size();i++){grouped_Wyckoff_positions[i].element=names[i];}
    permutation_grouped_Wyckoff_positions.push_back(grouped_Wyckoff_positions);

    uint i=0;
    while(i<num_elements){
      if(new_indices[i] < i){
        //LDEBUBcerr << "orig: ";
        //for(uint n=0;n<indices.size();n++){cerr << indices[n] << " ";}
        if(i%2==0){
          //LDEBUBcerr << "even: swapping " << indices[0] << " and " << indices[i] << endl;
          int swap1 = indices[0];
          int swap2 = indices[i];
          indices[0]=swap2; indices[i]=swap1;
          string swap_name1 = names[0];
          string swap_name2 = names[i];
          //GroupedWyckoffPosition swap_position1 = grouped_Wyckoff_positions[0]; swap_position1.element = names[i];
          //GroupedWyckoffPosition swap_position2 = grouped_Wyckoff_positions[i]; swap_position2.element = names[0];
          names[0]=swap_name2; names[i]=swap_name1;
          grouped_Wyckoff_positions[0].element=swap_name2; grouped_Wyckoff_positions[i].element=swap_name1;
        }
        else {
          //LDEBUBcerr << "odd: swapping " << indices[new_indices[i]] << " and " << indices[i] << endl;
          int swap1 = indices[new_indices[i]];
          int swap2 = indices[i];
          indices[new_indices[i]]=swap2; indices[i]=swap1;
          string swap_name1 = names[new_indices[i]];
          string swap_name2 = names[i];
          //GroupedWyckoffPosition swap_position1 = grouped_Wyckoff_positions[new_indices[i]]; swap_position1.element = names[i];
          //GroupedWyckoffPosition swap_position2 = grouped_Wyckoff_positions[i]; swap_position2.element = names[new_indices[i]];
          names[new_indices[i]]=swap_name2; names[i]=swap_name1;
          grouped_Wyckoff_positions[new_indices[i]].element=swap_name2; grouped_Wyckoff_positions[i].element=swap_name1;
        }
        //LDEBUBcerr << "storing: ";
        //for(uint n=0;n<indices.size();n++){cerr << indices[n] << " ";}
        //cerr << endl;
        permutations.push_back(indices);
        name_order.push_back(names);
        permutation_grouped_Wyckoff_positions.push_back(grouped_Wyckoff_positions);
        new_indices[i]++;
        i=0;
      } 
      else {
        //LDEBUBcerr << "moving on" << endl;
        new_indices[i]=0;
        i++;
      }   
    }
    return true;
  }
}
//ABOVE FUNCTION MAY NEED A DOUBLE CHECK

// ***************************************************************************
// makePermutations 
// ***************************************************************************
namespace compare{
  bool makePermutations(StructurePrototype& structure, vector<vector<string> >& name_order, vector<StructurePrototype>& permutation_structures){

    // make vector<StructurePrototype> of permutations
    
    //vector<int> unique_stoich;
    //vector<vector<int> > type_index;
    //groupSameRatios(stoich,unique_stoich,type_index);

    for(uint i=0;i<name_order.size();i++){
      xstructure xstr_tmp = structure.representative_structure;
      deque<string> species; 
      for(uint j=0;j<name_order[i].size();j++){cerr << name_order[i][j] << endl; species.push_back(name_order[i][j]);}
      xstr_tmp.SetSpecies(species);
      xstr_tmp.SpeciesPutAlphabetic();
      deque<int> sizes = SYM::arrange_atoms(xstr_tmp.atoms);
      //LDEBUGfor(uint j=0;j<sizes.size();j++){cerr << "sizes[j]: " << sizes[j] << endl;}
      xstr_tmp = pflow::SetNumEachType(xstr_tmp, sizes);
      //if (xstr_out.num_each_type.size() != names.size()){
      //  xstr_out = pflow::SetAllAtomNames(xstr_out, in_names);
      //}

      StructurePrototype tmp;
      tmp.representative_structure = xstr_tmp;
      tmp.representative_structure_name = aurostd::joinWDelimiter(species,"");
      tmp.representative_structure_generated = true;
      tmp.representative_structure_from = "permutation";
      tmp.copyPrototypeInformation(structure);
      permutation_structures.push_back(tmp);
    }
   
    return true;
  }
}



// ***************************************************************************
// makePermutations 
// ***************************************************************************
namespace compare{
  bool makePermutations(xstructure& xstr, vector<vector<string> >& name_order, vector<xstructure>& xstr_permutations){

    // make xstructures of permutations
    
    //vector<int> unique_stoich;
    //vector<vector<int> > type_index;
    //groupSameRatios(stoich,unique_stoich,type_index);

    for(uint i=0;i<name_order.size();i++){
      xstructure xstr_tmp = xstr;
      deque<string> species; 
      for(uint j=0;j<name_order[i].size();j++){cerr << name_order[i][j] << endl; species.push_back(name_order[i][j]);}
      xstr_tmp.SetSpecies(species);
      xstr_tmp.SpeciesPutAlphabetic();
      deque<int> sizes = SYM::arrange_atoms(xstr_tmp.atoms);
      //LDEBUGfor(uint j=0;j<sizes.size();j++){cerr << "sizes[j]: " << sizes[j] << endl;}
      xstr_tmp = pflow::SetNumEachType(xstr_tmp, sizes);
      //if (xstr_out.num_each_type.size() != names.size()){
      //  xstr_out = pflow::SetAllAtomNames(xstr_out, in_names);
      //}
      xstr_permutations.push_back(xstr_tmp);
    }   
    return true;
  }
}

// ***************************************************************************
//  Get Stoichiometry - Obtain stoichiometries from composition string
// ***************************************************************************
namespace compare{
  vector<uint> getStoichiometry(string& composition, const bool& same_species){

    // Obtains the least common multiple representation of the stoichiometry.

    //cerr << "composition string: " << composition << endl;

    vector<uint> stoich;
    if(composition.size()==1){
       stoich.push_back(1);
    }
    else {
      bool is_previous_alpha = false;
      vector<uint> stoichiometry;
      stringstream tmp;
      for(uint i=0;i<composition.size();i++){
        if(isalpha(composition[i])){
          if(is_previous_alpha){
            stoichiometry.push_back(1);
          }
          else if(i!=0 && tmp.str().size() == 0){
            stoichiometry.push_back(aurostd::string2utype<uint>(tmp.str()));
            tmp.str("");
          }
        }
        else if(isdigit(composition[i])){
          tmp << composition[i];
        }    
        is_previous_alpha = isalpha(composition[i]);
      }
      if(is_previous_alpha){
        stoichiometry.push_back(1);
      }
      else if(tmp.str().size() == 0){
        stoichiometry.push_back(aurostd::string2utype<uint>(tmp.str()));
        tmp.str("");
      }
      //::print(stoichiometry);
      stoich=gcdStoich(stoichiometry);
    }
    // If a structure prototype comparison (not material type), ensure 
    // stoichiometries are in numerical order for comparison.  Else, 
    // leave in position indicating atomic species count.
    if(same_species==false){
      for(uint i=0; i<stoich.size(); i++){
	      std::sort(stoich.begin(),stoich.end());
      }
    }
    return stoich;
  }
}

// ***************************************************************************
//  addAFLOWPrototypes2StructurePrototypeVector() 
// ***************************************************************************
namespace compare{
  bool addAFLOWPrototypes2StructurePrototypeVector(vector<StructurePrototype>& all_structures, vector<string>& vlabel){ 

    // Add the AFLOW prototypes to the vector of StructurePrototype objects
    // Note: The AFLOW labels should already be filtered to the relevant
    // strutures for comparison
    
    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "addAFLOWPrototypes2StructurePrototypeVector()";

    for(uint i=0;i<vlabel.size();i++){
      if(LDEBUG) { cerr << function_name << ": Storing AFLOW prototype information for " << vlabel[i] << endl; }

      // anrl prototypes
      vector<string> tokens;
      if(aurostd::string2tokens(vlabel[i],tokens,"_")>=4) {
        vector<string> vparameter_values = anrl::getANRLParameters(vlabel[i],"all");
        for(uint j=0;j<vparameter_values.size();j++){
					StructurePrototype structure_tmp;
					//if one degree of freedom, then no number scheme required
					if(aurostd::string2tokens(vparameter_values[j],tokens,",")==1){
						structure_tmp.representative_structure_name = vlabel[i];
					}
					//if multiple degrees of freedom, then number scheme is required
					else {
						stringstream tmp; tmp << vlabel[i] << "-" << std::setw(3) << std::setfill('0') << j+1;
						structure_tmp.representative_structure_name = tmp.str();
    			}
					structure_tmp.stoichiometry = all_structures[0].stoichiometry;
					//structure_tmp.space_group = prototype_space_groups[i]; 
					structure_tmp.space_group = all_structures[0].space_group; // same as representative structure (either will be the same, or we are forcing it to be for the ignore_symmetry/ignore_Wyckoff run)
					structure_tmp.grouped_Wyckoff_positions = all_structures[0].grouped_Wyckoff_positions;
					structure_tmp.elements = compare::fakeElements(all_structures[0].stoichiometry.size());
				  structure_tmp.representative_structure_compound = compare::getCompoundName(structure_tmp.elements,structure_tmp.stoichiometry,true); //remove ones is true 
					structure_tmp.representative_structure_generated = false; 
					structure_tmp.representative_structure_from = "aflow_prototypes"; 
					all_structures.push_back(structure_tmp);
			  }
			}

			// htqc prototypes
			else {
				StructurePrototype structure_tmp;
				structure_tmp.representative_structure_name = vlabel[i];
				structure_tmp.stoichiometry = all_structures[0].stoichiometry;
				//structure_tmp.space_group = prototype_space_groups[i];
				structure_tmp.space_group = all_structures[0].space_group; // same as representative structure (either will be the same, or we are forcing it to be for the ignore_symmetry/ignore_Wyckoff run)
				structure_tmp.grouped_Wyckoff_positions = all_structures[0].grouped_Wyckoff_positions;
				vector<string> elements;
			  structure_tmp.elements = compare::fakeElements(all_structures[0].stoichiometry.size());
				structure_tmp.representative_structure_compound = compare::getCompoundName(structure_tmp.elements,structure_tmp.stoichiometry,true); //remove ones is true 
				structure_tmp.representative_structure_generated = false; 
				structure_tmp.representative_structure_from = "aflow_prototypes"; 
				all_structures.push_back(structure_tmp);
			}
    }
    return true;
  }
} 

// ***************************************************************************
//  Get Stoichiometry - Obtain stoichiometries from xstructure
// ***************************************************************************
namespace compare{
  string getCompoundName(xstructure& xstr, bool remove_ones){
    
    // Obtains compound name from xstructure in the reduced stoichiometric form.

    vector<uint> stoichiometry = compare::getStoichiometry(xstr,true);
    vector<string> elements = compare::getElements(xstr);

    return getCompoundName(elements,stoichiometry,remove_ones);
  }
}

// ***************************************************************************
//  Get Stoichiometry - Obtain stoichiometries from xstructure
// ***************************************************************************
namespace compare{
  string getCompoundName(vector<string>& elements, vector<uint>& stoichiometry, bool remove_ones){
    
    // Obtains compound name from the elements and stoichiometry vectors.
    // Note: Does not guarantee reduced stoichiometric form.

    stringstream compound_name;
    if(stoichiometry.size() == elements.size()){
      for(uint i=0;i<stoichiometry.size();i++){
        if(remove_ones && stoichiometry[i]==1){
          compound_name << elements[i];
        }
        else {
          compound_name << elements[i] << stoichiometry[i];
        }
      }
    }
    else {
      cerr << "compare::getCompoundName():WARNING: Size of elements != stoichiometry. Not returning reduced compound name" << endl;
      return "";
    }
    return compound_name.str();
  }
}

// ***************************************************************************
//  Get Stoichiometry - Obtain stoichiometries from xstructure
// ***************************************************************************
namespace compare{
  vector<uint> getStoichiometry(const xstructure& xstr, const bool& same_species){
    
    // Obtains the least common multiple representation of the stoichiometry.

    vector<uint> stoich;
    if(xstr.species.size()==1){
       stoich.push_back(1);
    }
    else {
      stoich=gcdStoich(xstr.num_each_type);
    }
    // If a structure prototype comparison (not material type), ensure 
    // stoichiometries are in numerical order for comparison.  Else, 
    // leave in position indicating atomic species count.
    if(same_species==false){
      for(uint i=0; i<stoich.size(); i++){
	      std::sort(stoich.begin(),stoich.end());
      }
    }
    return stoich;
  }
}

// ***************************************************************************
//  Get Elements 
// ***************************************************************************
namespace compare{
  vector<string> getElements(xstructure& xstr){
    
    // Obtains the elements in the xstructure.

    bool LDEBUG=(false || XHOST.DEBUG);
    vector<string> velements;
    // If atoms in poscar not labeled in either POSCAR; assign fake names
    if (xstr.atoms[0].name == ""){
      if(LDEBUG) {cerr << "compare::getElements():" << "WARNING!!!!! Atoms not labeled ... Assigning Fake names" << endl;}
      fakeAtomsName(xstr);
    }

    string prev_element="";
    for(uint i=0; i<xstr.atoms.size(); i++){
      if(KBIN::VASP_PseudoPotential_CleanName(xstr.atoms[i].name) != prev_element){ //DX 20190329 - remove pseudopotential info
	      velements.push_back(KBIN::VASP_PseudoPotential_CleanName(xstr.atoms[i].name)); //DX 20190329 - remove pseudopotential info
	      prev_element=KBIN::VASP_PseudoPotential_CleanName(xstr.atoms[i].name); //DX 20190329 - remove pseudopotential info
      }
    }
    return velements;
  }
}

// ***************************************************************************
// gcdStoich - Euler's Greatest Common Divisor Algorithm
// ***************************************************************************
namespace compare{
  vector<uint> gcdStoich(const vector<uint>& numbers){

    // This is Euler's Greates Common Divisor Algorithm.  It is used to determine 
    // the least common multiple representation for the stoichiometry.
    // vector version

    deque<int> number_deque; 
    for(uint i=0;i<numbers.size();i++){number_deque.push_back((int)numbers[i]);}
    return gcdStoich(number_deque);
  }
}

namespace compare{
  vector<uint> gcdStoich(const deque<int>& numbers){

    // This is Euler's Greates Common Divisor Algorithm.  It is used to determine 
    // the least common multiple representation for the stoichiometry.
    // deque version

    int global_GCD = 0; //DX 5/14/18 - added initialization
    int GCD = 0; //DX 5/14/18 - added initialization
    vector<uint> reduced_numbers;
    // Find min number first
    int min=0;
    for(uint i=0; i<numbers.size(); i++){
      if(i==0){
	min=numbers[i];
      }
      else {
	if(numbers[i]<min){
	   min=numbers[i];
	}
      }
    }
    bool found_GCD=true;
    for(uint i=0; i<numbers.size(); i++){
      if(numbers[i]%min != 0){
	found_GCD=false;
	break;
      }
    }
    if(found_GCD==true){
      global_GCD=min;
    }
    else if(found_GCD==false){
      int remainder=1000;
      int divisor=min;
      for(uint i=0; i<numbers.size(); i++){
	int num=numbers[i];
	while(remainder != 0){
	  remainder=(num%divisor);
	  if(remainder==0){
	    GCD=divisor;
	  }
	  else {
	    num=divisor;
	    divisor=remainder;
	  }
	}
	divisor=GCD;
	if(i==0){
	  global_GCD=GCD;
	}
	else if(GCD < global_GCD){
	   global_GCD=GCD;
	}
	remainder=100;
      }
    }
    for(uint i=0; i<numbers.size(); i++){
       reduced_numbers.push_back((uint)(numbers[i]/global_GCD));
       if(numbers[i]%global_GCD){
	  cerr << "compare::ERROR - Error in GCD procedure. Exiting. (David Hicks: d.hicks@duke.edu)" << endl;
	  exit(1);
       }
    }
    return reduced_numbers;
  }
}

// ***************************************************************************
// prepareSymmetryThreads - 
// ***************************************************************************
namespace compare{
  vector<vector<xstructure> > prepareSymmetryThreads(vector<xstructure>& vxstrs, uint& num_proc){

    // Split xstructures into a vector, i.e., to be used in different threads for 
    // calculating the symmetry (space group and Wyckoff positions)
    // More efficient method is to split into indices; see other funtion with the same name (below)

	  vector<vector<xstructure> > vxstrs_split;
	  uint num_per_thread = vxstrs.size()/num_proc;
	  uint residual = vxstrs.size()%num_proc;
	  bool accounted_for_residual=false;
	  if(residual!=0){num_per_thread+=1;}
    
    uint count = 0;
    uint thread_count = 0;
    vector<xstructure> tmp;
    for(uint l=0; l<vxstrs.size(); l++){
      tmp.push_back(vxstrs[l]);
      count+=1;
      if(count == num_per_thread && thread_count<num_proc-1){
        thread_count+=1;
        vxstrs_split.push_back(tmp);
        tmp.clear();
        count = 0;
      }
      else if(thread_count==num_proc-1 && l==vxstrs.size()-1){
        thread_count+=1;
        vxstrs_split.push_back(tmp);
        tmp.clear();
        count = 0;
      }
      if(!accounted_for_residual && residual!=0 && thread_count==residual){
        accounted_for_residual=true;
        num_per_thread=num_per_thread-1;
      }
    }

    uint recovered=0;
    //Need the following safety in case the number of threads is greater than the number of xstructures to test
    uint num_of_threads=0;
    if(vxstrs_split.size()>=num_proc){
      num_of_threads=num_proc;
    }
    else if(vxstrs_split.size()<num_proc){
      num_of_threads=vxstrs_split.size();
    }
    for(uint n=0; n<num_of_threads; n++){
      for(uint h=0;h<vxstrs_split[n].size();h++){
        recovered+=1;
      }
    } 
    if(recovered != vxstrs.size()){
      cerr << "compare::prepareSymmetryThreads(): The splitting of jobs failed...not all were accounted for: " << recovered << " != " << vxstrs.size() << endl;
      exit(1);
    }
    //DEBUG for(uint i=0;i<vxstrs_split.size();i++){
    //DEBUG   cerr << "num of xstrs for thread: " << i << " = " << vxstrs_split[i].size() << endl;
    //DEBUG }
    return vxstrs_split;
  }
}

// ***************************************************************************
// prepareSymmetryThreads - 
// ***************************************************************************
namespace compare{
  bool prepareSymmetryThreads(vector<xstructure>& vxstrs, uint& num_proc,
                              vector<uint>& start_indices, vector<uint>& end_indices){

    // Split xstructures via indices, i.e., to be used in different threads for 
    // calculating the symmetry (space group and Wyckoff positions)

	  vector<vector<xstructure> > vxstrs_split;
	  uint num_per_thread = vxstrs.size()/num_proc;
	  uint residual = vxstrs.size()%num_proc;
	  bool accounted_for_residual=false;
	  if(residual!=0){num_per_thread+=1;}
    
    uint count = 0;
    uint thread_count = 0;
    vector<xstructure> tmp;
    uint tmp_start_index=0;
    for(uint l=0; l<vxstrs.size(); l++){
      count+=1;
      if(count == num_per_thread && thread_count<num_proc-1){
        thread_count+=1;
        start_indices.push_back(tmp_start_index);
        end_indices.push_back(l);
        tmp_start_index=l+1;
        count = 0;
      }
      else if(thread_count==num_proc-1 && l==vxstrs.size()-1){
        thread_count+=1;
        start_indices.push_back(tmp_start_index);
        end_indices.push_back(l);
        tmp_start_index=l+1;
        count = 0;
      }
      if(!accounted_for_residual && residual!=0 && thread_count==residual){
        accounted_for_residual=true;
        num_per_thread=num_per_thread-1;
      }
    }

    //Need the following safety in case the number of threads is greater than the number of structures to test
    uint recovered=0;
    uint num_of_threads=0;
    if(start_indices.size()>=num_proc){
      num_of_threads=num_proc;
    }
    else if(start_indices.size()<num_proc){
      num_of_threads=start_indices.size();
    }
    for(uint n=0; n<num_of_threads; n++){
      for(uint i=start_indices[n];i<=end_indices[n];i++){
        recovered+=1;
      }
    } 
    if(recovered != vxstrs.size()){
      cerr << "compare::prepareSymmetryThreads(): The splitting of jobs failed...not all were accounted for: " << recovered << " != " << vxstrs.size() << endl;
      exit(1);
    }
    //DEBUG for(uint i=0;i<vxstrs_split.size();i++){
    //DEBUG   cerr << "num of xstrs for thread: " << i << " = " << vxstrs_split[i].size() << endl;
    //DEBUG }
    return true;
  }
}


// ***************************************************************************
// prepareComparisonThreads - 
// ***************************************************************************
namespace compare{
  bool splitComparisonIntoThreads(vector<StructurePrototype>& comparison_schemes, uint& num_proc,
                              vector<std::pair<uint,uint> >& start_indices,
                              vector<std::pair<uint,uint> >& end_indices){

    // split comparisons into threads by indicating comparison indices
    
    bool LDEBUG=(false || XHOST.DEBUG);
    uint number_of_comparisons = 0;
    for(uint i=0;i<comparison_schemes.size();i++){ number_of_comparisons += comparison_schemes[i].numberOfComparisons(); }
    //DEBUG cerr << "# of comparisons: " << number_of_comparisons << endl;

    if(number_of_comparisons==0){
      if(LDEBUG) {
        cerr << "Number of comparisons is zero, no need to split into threads." << endl;
      } 
      return true;
    }   
 
	  uint num_per_thread = number_of_comparisons/num_proc;
	  uint residual = number_of_comparisons%num_proc;
	  bool accounted_for_residual=false;
	  if(residual!=0){num_per_thread+=1;}
   
    if(LDEBUG) {
      cerr << "Number of comparisons per thread: " << num_per_thread << endl;
    }

    uint tmp =0;

    uint count = 0;
    uint thread_count = 0;
    std::pair<uint,uint> tmp_start, tmp_end;
    std::pair<uint,uint> indices;
    for(uint i=0;i<comparison_schemes.size();i++){
      //DEBUG cerr << "splitting comparison indices: i: " << i << "/" << comparison_schemes.size() << endl;
      for(uint j=0;j<comparison_schemes[i].duplicate_structures_names.size();j++){
        indices.first=i, indices.second=j;
        count+=1;
        tmp+=1;
        if(count == num_per_thread && thread_count<num_proc-1){
          thread_count+=1;
          start_indices.push_back(tmp_start);
          //update tmp_start
          if(j+1>=comparison_schemes[i].duplicate_structures_names.size()-1 && i+1<comparison_schemes.size()-1){
            tmp_start.first=i+1; tmp_start.second=0; 
            tmp_end.first=i+1; tmp_end.second=0;
          }
          else {
            tmp_start.first=i; tmp_start.second=j+1;
            tmp_end.first=i; tmp_end.second=j+1;
          }
          end_indices.push_back(tmp_end);
          count = 0;
        }
        else if(thread_count==num_proc-1 && i==comparison_schemes.size()-1 && j==comparison_schemes[i].duplicate_structures_names.size()-1){
          thread_count+=1;
          start_indices.push_back(tmp_start);
          //update tmp_start
          if(j+1>=comparison_schemes[i].duplicate_structures_names.size()-1 && i+1<comparison_schemes.size()-1){
            tmp_start.first=i+1; tmp_start.second=0; 
            tmp_end.first=i+1; tmp_end.second=0;
          }
          else {
            tmp_start.first=i; tmp_start.second=j+1;
            tmp_end.first=i; tmp_end.second=j+1;
          }
          end_indices.push_back(tmp_end);
          count = 0;
        }
        if(!accounted_for_residual && residual!=0 && thread_count==residual){
          accounted_for_residual=true;
          num_per_thread=num_per_thread-1;
        }
      }
    }
    // need to add last if it was not already added, using "indices" variable to populate last indices
    if(count){
      start_indices.push_back(tmp_start);
      //tmp_end.first=indices.first; tmp_end.second=indices.second; //+1 
      tmp_end.first=comparison_schemes.size()-1; tmp_end.second=comparison_schemes[tmp_end.first].duplicate_structures_names.size();
      end_indices.push_back(tmp_end);
    }

    //cerr << "number of indices indicators: " << start_indices.size() << endl;
    //cerr << "number gone through: " << tmp << endl;
    uint recovered=0;
    //Need the following safety in case the number of threads is greater than the number of xstructures to test
    uint num_of_threads=0;
    if(start_indices.size()>=num_proc){
      num_of_threads=num_proc;
    }
    else if(start_indices.size()<num_proc){
      num_of_threads=start_indices.size();
    }

    for(uint n=0; n<num_of_threads; n++){
      //cerr << start_indices[n].first << "," << start_indices[n].second << " - " << end_indices[n].first << "," << end_indices[n].second << endl;
      uint i_min=start_indices[n].first; uint i_max=end_indices[n].first;
      uint j_min=0; uint j_max=0;
      for(uint i=0;i<comparison_schemes.size();i++){
        // to loop properly
        if(i==i_min){
          j_min=start_indices[n].second;
          if(i==i_max){j_max=end_indices[n].second;}
          else {j_max=comparison_schemes[i].duplicate_structures_names.size();} //-1 since in loop: j<=j_max
        }
        else if(i==i_max){j_min=0; j_max=end_indices[n].second;}
        else {j_min=0; j_max=comparison_schemes[i].duplicate_structures_names.size();} //-1 since in loop: j<=j_max
        for(uint j=0;j<comparison_schemes[i].duplicate_structures_names.size();j++){
          if(i>=i_min && j>=j_min &&
             i<=i_max && j<j_max){
            recovered+=1;
          }
        }
      }
    } 
    if(recovered != number_of_comparisons){
      cerr << "compare::splitComparisonsIntoThreads(): The splitting of jobs failed...not all were accounted for: " << recovered << " != " << number_of_comparisons << endl;
      exit(1);
    }
    return true;
  }
}

// ***************************************************************************
// calculateSymmetries - Calculate Symmetries
// ***************************************************************************
namespace compare{
  void calculateSymmetries(vector<xstructure>& vxstrs, vector<string>& vpearsons, vector<uint>& vsgroups, 
                         vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions, uint& num_proc){

    // Calculates the symmetry (space group and Wyckoff positions) of each structure 
    // and stores it in the relevant vector
    // Same as the StructurePrototype version, just not as concise

    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name="compare::calculateSymmetries()";
    if(LDEBUG) {cerr << function_name << ": Number of threads=" << num_proc << endl;}

#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
    // THREADED VERSION - START
    
    // Distribute threads via indices
    vector<uint> start_indices, end_indices;
    prepareSymmetryThreads(vxstrs,num_proc,start_indices,end_indices);

    // Run threads 
    vector<std::thread> threads;
    for(uint n=0; n<num_proc; n++){
	    threads.push_back(std::thread(SYM::calculateSpaceGroupsInSetRange,std::ref(vxstrs),std::ref(start_indices[n]),std::ref(end_indices[n])));
    }
    // Join threads 
	  for(uint t=0;t<num_proc;t++){
      threads[t].join();
	  }
    // THREADED VERSION - END

#else
    // NONTHREADS - START
    for(uint i=0; i<vxstrs.size(); i++){
      compare::calculateSymmetry(vxstrs[i],vpearsons,vsgroups,vgrouped_Wyckoff_positions);
    }   
    // NONTHREADS - END

#endif

    // Populate symmetry vectors
    for(uint i=0;i<vxstrs.size();i++){
      vpearsons.push_back("xX");
      vsgroups.push_back(vxstrs[i].space_group_ITC);
      vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
      groupWyckoffPositions(vxstrs[i], grouped_Wyckoff_positions);
      vgrouped_Wyckoff_positions.push_back(grouped_Wyckoff_positions);
    }
  }
}

// ***************************************************************************
// prepareSymmetryThreads - 
// ***************************************************************************
namespace compare{
  bool prepareSymmetryThreads(uint& number_of_structures, uint& num_proc,
                              vector<uint>& start_indices, vector<uint>& end_indices){

    // Split number of structures via indices, i.e., to be used in different threads for 
    // calculating the symmetry (space group and Wyckoff positions)
    // Could be generalized for any type of job splitting

	  uint num_per_thread = number_of_structures/num_proc;
	  uint residual = number_of_structures%num_proc;
	  bool accounted_for_residual=false;
	  if(residual!=0){num_per_thread+=1;}
    
    uint count = 0;
    uint thread_count = 0;
    uint tmp_start_index=0;
    for(uint l=0; l<number_of_structures; l++){
      count+=1;
      if(count == num_per_thread && thread_count<num_proc-1){
        thread_count+=1;
        start_indices.push_back(tmp_start_index);
        end_indices.push_back(l);
        tmp_start_index=l+1;
        count = 0;
      }
      else if(thread_count==num_proc-1 && l==number_of_structures-1){
        thread_count+=1;
        start_indices.push_back(tmp_start_index);
        end_indices.push_back(l);
        tmp_start_index=l+1;
        count = 0;
      }
      if(!accounted_for_residual && residual!=0 && thread_count==residual){
        accounted_for_residual=true;
        num_per_thread=num_per_thread-1;
      }
    }

    //Need the following safety in case the number of threads is greater than the number of structures to test
    uint recovered=0;
    uint num_of_threads=0;
    if(start_indices.size()>=num_proc){
      num_of_threads=num_proc;
    }
    else if(start_indices.size()<num_proc){
      num_of_threads=start_indices.size();
    }
    for(uint n=0; n<num_of_threads; n++){
      for(uint i=start_indices[n];i<=end_indices[n];i++){
        recovered+=1;
      }
    } 
    if(recovered != number_of_structures){
      cerr << "compare::prepareSymmetryThreads(): The splitting of jobs failed...not all were accounted for: " << recovered << " != " << number_of_structures << endl;
      exit(1);
    }
    return true;
  }
}

// ***************************************************************************
// calculateSpaceGroupsInSetRange
// ***************************************************************************
namespace compare {
  void calculateSpaceGroupsInSetRange(vector<StructurePrototype>& structures, uint& start_index, uint& end_index){
   
    // Calculates the space group and Wyckoff positions for the representative 
    // structure in the StructurePrototype object
    // Mirrors SYM::calculateSpaceGroupsInSetRange(), but is specific for 
    // StructurePrototype objects, as opposed to xstructures

    for(uint i=start_index;i<=end_index;i++){
      structures[i].space_group = structures[i].representative_structure.SpaceGroup_ITC();
      vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
      groupWyckoffPositions(structures[i].representative_structure, grouped_Wyckoff_positions);
      structures[i].grouped_Wyckoff_positions=grouped_Wyckoff_positions;
    }
  }
}

// ***************************************************************************
// calculateSymmetries - Calculate Symmetries
// ***************************************************************************
namespace compare{
  void calculateSymmetries(vector<StructurePrototype>& structures, uint& num_proc){

    // Calculates the symmetry (space group and Wyckoff positions) of each structure 
    // and stores it in the StructurePrototype object
    // Same as the vector<xstructure> version, just more concise
    
    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name="compare::calculateSymmetries()";
    if(LDEBUG) {cerr << function_name << ": Number of threads=" << num_proc << endl;}

#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
    // THREADED VERISON - START
    
    // Distribute threads via indices
    uint number_of_structures = structures.size();
    uint num_threads = aurostd::min(num_proc,number_of_structures); // cannot have more threads than structures
    vector<uint> start_indices, end_indices;
    prepareSymmetryThreads(number_of_structures,num_threads,start_indices,end_indices);

    // Run threads 
    vector<std::thread> threads;
    for(uint n=0; n<num_threads; n++){
	    threads.push_back(std::thread(compare::calculateSpaceGroupsInSetRange,std::ref(structures),std::ref(start_indices[n]),std::ref(end_indices[n])));
    }
    // Join threads
	  for(uint t=0;t<num_threads;t++){
      threads[t].join();
	  }
    // THREADED VERISON - END

#else
    // NON-THREADED VERSION - START
    for(uint i=0; i<structures.size(); i++){
      structures[i].calculateSymmetry();
    }   
    // NON-THREADED VERSION - END

#endif

  }
}

// ***************************************************************************
// calculateSymmetry - Calculate Symmetry
// ***************************************************************************
namespace compare{
  void calculateSymmetry(xstructure& xstr, vector<string>& vpearsons, vector<uint>& vsgroups, 
                         vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions){

    // Calculates the Pearson symbol and space group of the xstructure 
    // and stores into the relevant vector
    // There is an equivalent method for the StructurePrototype object

    //xstr.GetLatticeType(); //slow; consider a different method -> (SpaceGroup_ITC) finds this
    //vpearsons.push_back(xstr.pearson_symbol);
    vpearsons.push_back("");
    vsgroups.push_back(xstr.SpaceGroup_ITC());
    vector<GroupedWyckoffPosition> grouped_Wyckoff_positions; 
    groupWyckoffPositions(xstr, grouped_Wyckoff_positions);
    vgrouped_Wyckoff_positions.push_back(grouped_Wyckoff_positions);

    //for(uint i=0;i<xstr.wyckoff_sites_ITC.size();i++){
    // cerr << xstr.wyckoff_sites_ITC[i].wyckoffSymbol << endl;
    //}
  }
}

// ***************************************************************************
// groupWyckoffPositions
// ***************************************************************************
namespace compare{
  bool groupWyckoffPositions(xstructure& xstr, vector<GroupedWyckoffPosition>& grouped_positions){

    // Groups the Wyckoff positions via species 
    // Obtains information from xstructure
    // Assumes xstr.SpaceGroup_ITC() has been called, otherwise, this will fail

    uint type_count = 0; //DX 20190425 - add type indicator

    for(uint i=0;i<xstr.wyckoff_sites_ITC.size();i++){
      //DX 20191030 [OBOSLETE] vector<string> tokens;
      //DX 20191030 [OBOSLETE] aurostd::string2tokens(xstr.wyckoff_sites_ITC[i].wyckoffSymbol,tokens," ");  
      //DX 20191030 [OBOSLETE] uint multiplicity = aurostd::string2utype<uint>(tokens[0]);     
      //DX 20191030 [OBOSLETE] string letter = aurostd::string2utype<string>(tokens[1]);     
      //DX 20191030 [OBOSLETE] string site_symmetry = aurostd::string2utype<string>(tokens[2]);     
      uint multiplicity = xstr.wyckoff_sites_ITC[i].multiplicity; //DX 20191031 
      string letter = xstr.wyckoff_sites_ITC[i].letter; //DX 20191031
      string site_symmetry = xstr.wyckoff_sites_ITC[i].site_symmetry; //DX 20191031

      bool element_found = false;
      uint element_index = 0;
      for(uint j=0;j<grouped_positions.size();j++){
        if(KBIN::VASP_PseudoPotential_CleanName(xstr.wyckoff_sites_ITC[i].type) == KBIN::VASP_PseudoPotential_CleanName(grouped_positions[j].element)){ //DX 20190329 - remove pseudopotential info   
          element_found = true;
          element_index = j;
          break;
        }
      }
      if(element_found == false){
        GroupedWyckoffPosition tmp;
        tmp.type = type_count; //DX 20190425 - added type   
        tmp.element = KBIN::VASP_PseudoPotential_CleanName(xstr.wyckoff_sites_ITC[i].type); //DX 20190329 - remove pseudopotential info   
        tmp.site_symmetries.push_back(site_symmetry);          
        tmp.multiplicities.push_back(multiplicity);          
        tmp.letters.push_back(letter); //DX 20190208 - add Wyckoff letters      
        grouped_positions.push_back(tmp);
        type_count++; //DX 20190425
      }
      else {
        grouped_positions[element_index].site_symmetries.push_back(site_symmetry);
        grouped_positions[element_index].multiplicities.push_back(multiplicity);
        grouped_positions[element_index].letters.push_back(letter); //DX 20190208 - add Wyckoff letters
      }
    }

    //cerr << "xstr: " << xstr << endl;
    //for(uint j=0;j<grouped_positions.size();j++){
    //  cerr << "grouped wyckoffs: " << grouped_positions[j] << endl;
    //}
    return true;
  }
}

// ***************************************************************************
// groupWyckoffPositions
// ***************************************************************************
namespace compare{
  bool groupWyckoffPositionsFromGroupedString(uint& space_group_number, uint& setting, vector<vector<string> >& grouped_Wyckoff_string, vector<GroupedWyckoffPosition>& grouped_positions){
    
    // Groups the Wyckoff positions via species 
    // Obtains information from the string of the following form: a,f,g;b,c;a,a
    // (i.e., ";" separates by species, and "," splits the Wyckoff letters for that species)

    stringstream axis_cell;
    axis_cell.str(std::string());
    axis_cell << setting; //DX 20180806 - use setting
    SymmetryInformationITC ITC_sym_info; //DX 20190215
    ITC_sym_info.initsgs(axis_cell.str()); //DX 20190215
    //DX 20190215 [OBSOLETE] SYM::initsgs(axis_cell.str());
    //DX 20190215 [OBSOLETE] using SYM::gl_sgs;
    string spacegroupstring = ITC_sym_info.gl_sgs[space_group_number - 1]; //DX 20190215
    for(uint i=0;i<grouped_Wyckoff_string.size();i++){
      GroupedWyckoffPosition tmp;
      tmp.element = aurostd::utype2string<uint>(i);   
      for(uint j=0;j<grouped_Wyckoff_string[i].size();j++){
        string Wyckoff_letter = grouped_Wyckoff_string[i][j];
        uint Wyckoff_multiplicity = 0;
        string Wyckoff_site_symmetry = "";
        vector<string> positions;
        SYM::get_Wyckoff_from_letter(spacegroupstring, Wyckoff_letter, Wyckoff_multiplicity, Wyckoff_site_symmetry, positions);
        tmp.site_symmetries.push_back(Wyckoff_site_symmetry);          
        tmp.multiplicities.push_back(Wyckoff_multiplicity);          
        tmp.letters.push_back(Wyckoff_letter);          
      }
      grouped_positions.push_back(tmp);
    }
    return true;
  }
}

// ***************************************************************************
// groupWyckoffPositions
// ***************************************************************************
namespace compare{
  string printWyckoffString(const vector<GroupedWyckoffPosition>& grouped_positions, bool alphabetize){
    
    // Prints the grouped Wyckoff positions as a string of the following form: a,f,g;b,c;a,a
    // (i.e., ";" separates by species, and "," splits the Wyckoff letters for that species)

    vector<string> all_Wyckoff_sets;
    for(uint i=0;i<grouped_positions.size();i++){
      vector<string> Wyckoff_set;
      for(uint j=0;j<grouped_positions[i].letters.size();j++){
        Wyckoff_set.push_back(grouped_positions[i].letters[j]);
      }
      if(alphabetize){
        std::sort(Wyckoff_set.begin(),Wyckoff_set.end());
      }
      all_Wyckoff_sets.push_back(aurostd::joinWDelimiter(Wyckoff_set,","));
    }
    return aurostd::joinWDelimiter(all_Wyckoff_sets,";");
  }
}

// ***************************************************************************
// sortSiteSymmetryOfGroupedWyckoffPositions
// ***************************************************************************
namespace compare{
  vector<GroupedWyckoffPosition> sortSiteSymmetryOfGroupedWyckoffPositions(vector<GroupedWyckoffPosition>& grouped_Wyckoffs){

    // Sort they Wyckoff site symmetry for each Wyckoff position.
    // It is possible to match Wyckoff positions with different ordering of the 
    // site symmetry depending on the cell choice.
    // To account for this, we need to split the site symmetry into a vector 
    // containing the symmetries along the primary, secondary, and tertiary 
    // directions.  Then, we sort the vector (alphabetic will suffice, 
    // just need a standard comparison method).
    // NOTE: Sorting is only for comparing site symmetry strings to negate
    // choice of directions. Sortng may yield "site symmetries" that are not physical.

    vector<GroupedWyckoffPosition> sorted_site_symmetry_Wyckoff_positions;

    for(uint i=0;i<grouped_Wyckoffs.size();i++){
      sorted_site_symmetry_Wyckoff_positions.push_back(grouped_Wyckoffs[i]);
      for(uint j=0;j<grouped_Wyckoffs[i].site_symmetries.size();j++){
        vector<string> split_site_symmetry = SYM::splitSiteSymmetry(grouped_Wyckoffs[i].site_symmetries[j]);
        std::sort(split_site_symmetry.begin(),split_site_symmetry.end());
        sorted_site_symmetry_Wyckoff_positions[i].site_symmetries[j] = aurostd::joinWDelimiter(split_site_symmetry,"");
      }
    }
  
    return sorted_site_symmetry_Wyckoff_positions;
  }
}

// ***************************************************************************
// matchableWyckoffPositions
// ***************************************************************************
namespace compare{
  bool matchableWyckoffPositions(vector<GroupedWyckoffPosition>& temp_grouped_Wyckoffs,
                                 vector<GroupedWyckoffPosition>& representative_grouped_Wyckoffs, 
                                 const bool& same_species){

    // Determines if two sets of grouped Wyckoff positions are commensurate 
    // with one another, i.e., checks if the Wyckoff multiplicities and site 
    // symmetries are the same.  Comparing the same species is optional.

    bool LDEBUG=(false || XHOST.DEBUG);
    // quick check: are the number of Wyckoff positions the same; cannot match otherwise
    if(temp_grouped_Wyckoffs.size() != representative_grouped_Wyckoffs.size()){
      if(LDEBUG) {
        cerr << "compare::matchableWyckoffPositions(): # of Wyckoff positions does not match (" 
        << temp_grouped_Wyckoffs.size() << " vs " << representative_grouped_Wyckoffs.size() << endl;
      }
      return false;
    }

    //sort site symmetries to account for different cell choices
    vector<GroupedWyckoffPosition> sorted_temp_grouped_Wyckoffs = compare::sortSiteSymmetryOfGroupedWyckoffPositions(temp_grouped_Wyckoffs);
    vector<GroupedWyckoffPosition> sorted_representative_grouped_Wyckoffs = compare::sortSiteSymmetryOfGroupedWyckoffPositions(representative_grouped_Wyckoffs);

    vector<vector<bool> > found_matches;
    for(uint i=0;i<sorted_temp_grouped_Wyckoffs.size();i++){
      vector<bool> tmp;
      for(uint m=0;m<sorted_temp_grouped_Wyckoffs[i].multiplicities.size();m++){
        tmp.push_back(false);
      }
      found_matches.push_back(tmp);
    }

    for(uint i=0;i<sorted_temp_grouped_Wyckoffs.size();i++){
      for(uint j=0;j<sorted_representative_grouped_Wyckoffs.size();j++){
        if(same_species && sorted_temp_grouped_Wyckoffs[i].element == sorted_representative_grouped_Wyckoffs[j].element &&
           sorted_temp_grouped_Wyckoffs[i].multiplicities.size() == sorted_representative_grouped_Wyckoffs[j].multiplicities.size()){
          uint match_counts = 0;
          for(uint m=0;m<sorted_temp_grouped_Wyckoffs[i].multiplicities.size();m++){
            for(uint n=0;n<sorted_representative_grouped_Wyckoffs[j].multiplicities.size();n++){
              if(sorted_temp_grouped_Wyckoffs[i].multiplicities[m] == sorted_representative_grouped_Wyckoffs[j].multiplicities[n] &&
                 sorted_temp_grouped_Wyckoffs[i].site_symmetries[m] == sorted_representative_grouped_Wyckoffs[j].site_symmetries[n]){
                found_matches[i][m] = true;
                match_counts++;
              }
            }
          }
        }
        else if(!same_species && sorted_temp_grouped_Wyckoffs[i].multiplicities.size() == sorted_representative_grouped_Wyckoffs[j].multiplicities.size()){
          uint match_counts = 0;
          for(uint m=0;m<sorted_temp_grouped_Wyckoffs[i].multiplicities.size();m++){
            for(uint n=0;n<sorted_representative_grouped_Wyckoffs[j].multiplicities.size();n++){
              if(sorted_temp_grouped_Wyckoffs[i].multiplicities[m] == sorted_representative_grouped_Wyckoffs[j].multiplicities[n] &&
                 sorted_temp_grouped_Wyckoffs[i].site_symmetries[m] == sorted_representative_grouped_Wyckoffs[j].site_symmetries[n] &&
                 !found_matches[i][m]){ // and not already matched
                found_matches[i][m] = true;
                match_counts++;
              }
            }
          } 
          // if any match, all need to match; otherwise the Wyckoff positions are not matchable
          if(match_counts>0 && match_counts != sorted_temp_grouped_Wyckoffs[i].multiplicities.size()){ 
            vector<bool> tmp; for(uint m=0;m<sorted_temp_grouped_Wyckoffs[i].multiplicities.size();m++){tmp.push_back(false);}
            //cerr << "match_counts do not match the multiplicity: " << match_counts << " vs " << sorted_temp_grouped_Wyckoffs[i].multiplicities.size() << endl;
            found_matches[i] = tmp;
          }
        }
      }
    }

    for(uint i=0;i<found_matches.size();i++){
      for(uint j=0;j<found_matches[i].size();j++){
        if(found_matches[i][j] == false){
          //cerr << "could not match!!!: " << i << " " << j << endl;
          return false;
        }
      }
    }
    return true;
  }
}

// ***************************************************************************
// convertWyckoffString2GroupedPositions - 
// ***************************************************************************
namespace compare{
  vector<vector<string> > convertANRLWyckoffString2GroupedPositions(string label){
    
    // Converts the ANRL Wyckoff string to a grouped Wyckoff position object
    // ANRL Wyckoff string example: a2b_4bcd_e2f3g

    vector<vector<string> > anrl_grouped_Wyckoff_letters;
    vector<string> anrl_Wyckoff_set, tmp; 
    aurostd::string2tokens(label, tmp, "_");
    for(uint i=0;i<tmp.size();i++){ if(i>2){ anrl_Wyckoff_set.push_back(tmp[i]); }}
    //::print(anrl_Wyckoff_set);
    for(uint i=0;i<anrl_Wyckoff_set.size();i++){
      vector<string> Wyckoff_letters_for_species;
      uint Wyckoff_count = 1;
      bool is_previous_digit = false;
      for(uint j=0;j<anrl_Wyckoff_set[i].size();j++){
        if(isdigit(anrl_Wyckoff_set[i][j])){
          stringstream ss_tmp; 
          if(is_previous_digit){
            ss_tmp << Wyckoff_count << anrl_Wyckoff_set[i][j]; 
          }
          else {
            ss_tmp << anrl_Wyckoff_set[i][j];
          }
          Wyckoff_count = aurostd::string2utype<uint>(ss_tmp.str());
          is_previous_digit=true;
          continue;
        }
        else {
          for(uint c=0;c<Wyckoff_count;c++){
            stringstream ss_tmp; ss_tmp << anrl_Wyckoff_set[i][j];
            Wyckoff_letters_for_species.push_back(ss_tmp.str());
          }
          Wyckoff_count = 1; //reset
          is_previous_digit=false;
        }
      }
      anrl_grouped_Wyckoff_letters.push_back(Wyckoff_letters_for_species);
    }
    //for(uint i=0;i<anrl_grouped_Wyckoff_letters.size();i++){
    //  ::print(anrl_grouped_Wyckoff_letters[i]);
    //}
    return anrl_grouped_Wyckoff_letters;
  }
}

// ***************************************************************************
// convertWyckoffString2GroupedPositions - 
// ***************************************************************************
namespace compare{
  vector<vector<string> > convertWyckoffString2GroupedPositions(string Wyckoff_letter_string){
    
    // Groups the Wyckoff positions via species 
    // Obtains information from the string of the following form: a,f,g;b,c;a,a
    // (i.e., ";" separates by species, and "," splits the Wyckoff letters for that species)

    vector<vector<string> > grouped_Wyckoff_letters;
    vector<string> Wyckoff_set, token; 
    aurostd::string2tokens(Wyckoff_letter_string, token, ";");
    for(uint i=0;i<token.size();i++){ Wyckoff_set.push_back(token[i]); }
    for(uint i=0;i<Wyckoff_set.size();i++){
      aurostd::string2tokens(Wyckoff_set[i], token, ",");
      vector<string> Wyckoff_letters_for_species;
      for(uint j=0;j<token.size();j++){ Wyckoff_letters_for_species.push_back(token[j]); }
      grouped_Wyckoff_letters.push_back(Wyckoff_letters_for_species);
    }
    return grouped_Wyckoff_letters;
  }
}

// ***************************************************************************
// WyckoffsMatchable 
// ***************************************************************************
namespace compare{
  bool matchableWyckoffPositionSet(vector<vector<vector<string> > > grouped_possible_Wyckoff_letters,
                                   vector<vector<string> > grouped_Wyckoff_letters){
    bool LDEBUG=(false || XHOST.DEBUG);
    bool all_Wyckoffs_matched = false;
    //quick check
    if(grouped_possible_Wyckoff_letters.size() != grouped_Wyckoff_letters.size()){
      //cerr << "sizes do not match" << endl;
      return false;
    }

    //check if each set has the same number of Wyckoff letters for a given species
    vector<int> number_of_letters_1, number_of_letters_2;
    for(uint i=0;i<grouped_possible_Wyckoff_letters.size();i++){
      number_of_letters_1.push_back((int)grouped_possible_Wyckoff_letters[i].size());
      number_of_letters_2.push_back((int)grouped_Wyckoff_letters[i].size());
    }
    aurostd::sort(number_of_letters_1);
    aurostd::sort(number_of_letters_2);
    
    for(uint i=0;i<number_of_letters_1.size();i++){
      if(number_of_letters_1[i]!=number_of_letters_2[i]){
        if(LDEBUG) {cerr << "matchableWyckoffPositionSet(): Number of Wyckoff letters does not match between structures." << endl;} 
        return false;
      }
    }

    //identify matchable sets (i.e., same number of Wyckoff letters)
    //vector<bool> group_matched_once; for(uint i=0;i<grouped_possible_Wyckoff_letters.size();i++){group_matched_once.push_back(false);}
    //uint matches = 0;

    vector<std::pair<uint,uint> > matchable_indices;
    uint number_of_matched_sets = 0;
    for(uint m=0;m<grouped_Wyckoff_letters.size();m++){
    //for(uint i=0;i<grouped_possible_Wyckoff_letters.size();i++){
      bool matched_set = false;
      for(uint i=0;i<grouped_possible_Wyckoff_letters.size();i++){
      //for(uint m=0;m<grouped_Wyckoff_letters.size();m++){
        if(grouped_possible_Wyckoff_letters[i].size()==grouped_Wyckoff_letters[m].size()){
          uint matched_Wyckoffs = 0;
          vector<uint> index_matched;
          for(uint j=0;j<grouped_possible_Wyckoff_letters[i].size();j++){
            bool Wyckoff_matched = false;
            for(uint k=0;k<grouped_possible_Wyckoff_letters[i][j].size();k++){
              for(uint n=0;n<grouped_Wyckoff_letters[m].size();n++){
                bool index_used = false;
                for(uint t=0;t<index_matched.size();t++){if(n==index_matched[t]){index_used=true;}}
                if(index_used){break;}
                if(grouped_possible_Wyckoff_letters[i][j][k] == grouped_Wyckoff_letters[m][n]){
                  index_matched.push_back(n);
                  Wyckoff_matched = true;
                  matched_Wyckoffs++;
                  break;
                }
              }
              if(Wyckoff_matched==true){break;}
            }
          }
          if(matched_Wyckoffs==grouped_possible_Wyckoff_letters[i].size() && matched_Wyckoffs==grouped_Wyckoff_letters[m].size()){
            matched_set=true;
          }
        }
        if(matched_set==true){
          number_of_matched_sets++;
          matched_set=false;
        }
      }
    }
    if(number_of_matched_sets==grouped_possible_Wyckoff_letters.size()){
      //cerr << "possible match!" << endl;
      return true;
    }
    cerr << number_of_matched_sets << " == " << grouped_possible_Wyckoff_letters.size() << endl;
                
    //cerr << "all_Wyckoffs_matched: " << all_Wyckoffs_matched << endl;
    return all_Wyckoffs_matched;
  }
}

// ***************************************************************************
// matchable space groups 
// ***************************************************************************
namespace compare{
  bool matchableSpaceGroups(uint space_group_1, uint space_group_2){

    // Determines if two space groups are commensurate
    // The space groups must either be the same or enantiomorphic pairs

    if(space_group_1 == space_group_2){return true;}
    else {
      return matchableEnantiomorphicSpaceGroups(space_group_1, space_group_2);
    }
    return false;
  }
}

// ***************************************************************************
// matchable enantiomorphs 
// ***************************************************************************
namespace compare{
  bool matchableEnantiomorphicSpaceGroups(uint space_group_1, uint space_group_2){
    
    // Check if the space group has an enantimorphic pair
    // If it does have an enantiomorphic pair, and it is equal to the 
    // second space group, then return true
    // If it does not, the function below returns the input space group and 
    // compares it to the second space group 

    return (SYM::getEnantiomorphSpaceGroupNumber(space_group_1)==space_group_2);
  }
}

// ***************************************************************************
// filterPrototypes - Filter posssible prototypes by stoichiometry and symmetry
// ***************************************************************************
namespace compare{
  bool filterPrototypes(uint& species_count, string& reduced_stoichiometry, uint& space_group,
                        vector<vector<vector<string> > >& grouped_possible_Wyckoff_letters,
                        vector<string>& prototype_labels, vector<uint>& species_counts, 
                        vector<uint>& space_groups){
    vector<string> candidate_prototype_labels;
    vector<uint> candidate_prototype_species_counts;
    vector<string> candidate_prototype_stoichiometries;
    vector<uint> candidate_prototype_space_groups;

    // Filters AFLOW prototypes based on species count/stoich/spacegroup/Wyckoff positions
    // OBSOLETE: This functionality has been moved to aflow_xproto.cpp and has been improved


    //cerr << "structure: " << endl;
    //cerr << "species count: " << species_count << endl;
    //cerr << "reduced stoichiometry: " << reduced_stoichiometry << endl;
    //cerr << "space group: " << space_group << endl;

    for(uint i=0;i<prototype_labels.size();i++){ 
      if(species_count == species_counts[i]){
        //if(space_group == space_groups[i]){ // only 230, filter before stoich which can have virtually infinite
        if(matchableSpaceGroups(space_group,space_groups[i])){ // only 230, filter before stoich which can have virtually infinite
          vector<uint> anrl_stoichiometry = anrl::extractStoichiometry(prototype_labels[i]);
          string anrl_stoich_string = aurostd::joinWDelimiter(anrl_stoichiometry,":"); 
          if(reduced_stoichiometry == anrl_stoich_string){
            vector<vector<string> > anrl_grouped_Wyckoff_letters = convertANRLWyckoffString2GroupedPositions(prototype_labels[i]);
            bool all_Wyckoffs_matched = true;
            for(uint j=0;j<grouped_possible_Wyckoff_letters.size();j++){
              bool Wyckoff_matched = false;
              for(uint k=0;k<grouped_possible_Wyckoff_letters[j].size();k++){
                for(uint l=0;l<grouped_possible_Wyckoff_letters[j][k].size();l++){
                  if(grouped_possible_Wyckoff_letters[j][k][l] == anrl_grouped_Wyckoff_letters[j][k]){
                    Wyckoff_matched = true;
                    break;
                  }
                }
                if(Wyckoff_matched==true){ break; }
              }
              if(Wyckoff_matched == false){
                //could not match Wyckoff letters
                all_Wyckoffs_matched = false;
                break;
              }
            }
            if(all_Wyckoffs_matched){
              candidate_prototype_labels.push_back(prototype_labels[i]);
              candidate_prototype_space_groups.push_back(space_groups[i]);
            }
          }
        }
      }
    }
    prototype_labels = candidate_prototype_labels;
    space_groups = candidate_prototype_space_groups;
    return true;
  }
} 


// ***************************************************************************
// createStructurePrototypes - Group structures by Pearson symbol, then space group
// ***************************************************************************
namespace compare{
  void createStructurePrototypes(vector<StructurePrototype>& comparison_schemes, 
				 vector<xstructure>& vxstrs, const bool& same_species, 
				 const vector< vector<string> >& vvelements,
				 vector< vector<uint> >& vstoichs, vector<string>& vpearsons, 
				 vector<uint>& vsgroups, 
         vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
         const string& directory, const vector<string>& vfiles, 
         vector<bool>& vstructures_generated,
         vector<string>& vstructures_from,
         const bool& ignore_symmetry, const bool& ignore_Wyckoff){
    vector<string> property_names, property_units;
    vector<vector<string> > property_values;
    
    // Populates the structure information into the StructurePrototype object.
    // It groups structure based on their stoichiometry and pearson symbol and 
    // space group. A "representative" structure is chosen and will be compared to the 
    // possible "duplicates". The misfit values are set to -1.0 until compared.
    // Overloaded: In case the properties of the material are not given.
    // OBSOLETE: Created a cleaner function which requires less input arguments,
    // i.e., groupStructurePrototypes()   

    createStructurePrototypes(comparison_schemes, vxstrs, same_species, 
                                     vvelements, vstoichs, vpearsons, vsgroups, 
                                     vgrouped_Wyckoff_positions, property_names, property_units, property_values, 
                                     directory, vfiles,
                                     vstructures_generated, vstructures_from,
                                     ignore_symmetry, ignore_Wyckoff);
    
  }
}

// ***************************************************************************
// createStructurePrototypes - Group structures by Pearson symbol, then space group
// ***************************************************************************
namespace compare{
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
         const bool& ignore_symmetry, const bool& ignore_Wyckoff){

    // Populates the structure information into the StructurePrototype object.
    // It groups structure based on their stoichiometry and pearson symbol and 
    // space group. A "representative" structure is chosen and will be compared to the 
    // possible "duplicates". The misfit values are set to -1.0 until compared.
    // OBSOLETE: Created a cleaner function which requires less input arguments,
    // i.e., groupStructurePrototypes()   
 
    //cerr << "vstoichs: " << vstoichs.size() << endl;
    //cerr << "vstructures_generated: " << vstructures_generated.size() << endl;
    //cerr << "vstructures_from: " << vstructures_from.size() << endl;
    bool LDEBUG=(false || XHOST.DEBUG);

    // First, separate by stoichiometry
    for(uint i=0;i<vstoichs.size(); i++){
      bool scheme_created=false;
      if(i==0){
        StructurePrototype tmp;
        if(directory==""){
          tmp.representative_structure_name = vfiles[i];
        }
        else {
          tmp.representative_structure_name = directory+"/"+vfiles[i];  aurostd::StringSubst(tmp.representative_structure_name,"//","/"); //DX 20181003
        }
        tmp.representative_structure_generated=vstructures_generated[i];
        tmp.representative_structure_from=vstructures_from[i];
        //tmp.number_types=vxstrs[i].num_each_type.size();
        tmp.elements=vvelements[i];
        tmp.stoichiometry=vstoichs[i];
        //tmp.number_of_atoms=vxstrs[i].atoms.size();
        tmp.Pearson=vpearsons[i];
        tmp.space_group=vsgroups[i];
        tmp.grouped_Wyckoff_positions=vgrouped_Wyckoff_positions[i];
        if(property_names.size()!=0){
          tmp.property_names=property_names; //DX 20181218 - added property_names
          tmp.property_units=property_units; //DX 20181218 - added property_units
          tmp.representative_structure_properties=property_values[i]; //DX 20181218 - added property_values
        }
        if(vstructures_generated[i]){
          tmp.representative_structure=vxstrs[i];
          tmp.number_types=vxstrs[i].num_each_type.size();
          tmp.number_of_atoms=vxstrs[i].atoms.size();
          tmp.representative_structure_compound=getCompoundName(vxstrs[i]); //DX 20190111 - added compound, e.g., Ag1Br2
        }
        comparison_schemes.push_back(tmp);
      }
      else {
        for(uint j=0; j<comparison_schemes.size(); j++){
          bool same_material_stoich=false;
          ostringstream tmp;
          tmp.clear();
          if(same_species==true && 
             matchableSpecies(vxstrs[i],comparison_schemes[j].representative_structure,same_species)==true){
            same_material_stoich=true;
          }
          else if(same_species==false){
            same_material_stoich=true;
          }
          if(same_material_stoich==true && vstoichs[i] == comparison_schemes[j].stoichiometry && 
             ((ignore_symmetry && ignore_Wyckoff) ||    
              (!ignore_symmetry && ignore_Wyckoff &&
               vpearsons[i] == comparison_schemes[j].Pearson && 
               //vsgroups[i] == comparison_schemes[j].space_group) || 
               matchableSpaceGroups(vsgroups[i],comparison_schemes[j].space_group)) || 
              (!ignore_symmetry && !ignore_Wyckoff &&
             vpearsons[i] == comparison_schemes[j].Pearson && 
             //vsgroups[i] == comparison_schemes[j].space_group &&
             matchableSpaceGroups(vsgroups[i],comparison_schemes[j].space_group) &&
               matchableWyckoffPositions(vgrouped_Wyckoff_positions[i], comparison_schemes[j].grouped_Wyckoff_positions,same_species)))){
          //DX ORIG if(same_material_stoich==true && 
          //DX ORIG   vstoichs[i] == comparison_schemes[j].stoichiometry && 
          //DX ORIG   vpearsons[i] == comparison_schemes[j].Pearson && 
          //DX ORIG   vsgroups[i] == comparison_schemes[j].space_group &&
          //DX ORIG   matchableWyckoffPositions(vgrouped_Wyckoff_positions[i], comparison_schemes[j].grouped_Wyckoff_positions,same_species)){
            if(same_species==false){
              for(uint e=0;e<vvelements[i].size();e++){
                bool already_in=false;
                for(uint f=0;f<comparison_schemes[j].elements.size();f++){
                  if(vvelements[i][e]==comparison_schemes[j].elements[f]){
                    already_in=true;
                    break;
                  }
                }
                if(already_in==false){
                  comparison_schemes[j].elements.push_back(vvelements[i][e]);
                }
              }
            }
            string duplicate_name = "";
            if(directory==""){
              duplicate_name = vfiles[i];
            }
            else {
              duplicate_name = directory+"/"+vfiles[i];  aurostd::StringSubst(duplicate_name,"//","/"); //DX 20181003
            }
            comparison_schemes[j].duplicate_structures_names.push_back(duplicate_name);
            comparison_schemes[j].duplicate_structures_generated.push_back(vstructures_generated[i]);
            comparison_schemes[j].duplicate_structures_from.push_back(vstructures_from[i]);
            //cerr << "adding to " << j << " (name): " << comparison_schemes[j].duplicate_structures_names.size() << endl;
            //cerr << "adding to " << j << " (gen): " << comparison_schemes[j].duplicate_structures_generated.size() << endl;
            //cerr << "adding to " << j << " (from): " << comparison_schemes[j].duplicate_structures_from.size() << endl;
            if(vstructures_generated[i]){
              comparison_schemes[j].duplicate_structures.push_back(vxstrs[i]);
              comparison_schemes[j].duplicate_structures_compounds.push_back(getCompoundName(vxstrs[i])); //DX 20190111 - added compound, e.g., Ag1Br2
            }
            comparison_schemes[j].misfits.push_back(-1.0);
            if(property_names.size()!=0){
              comparison_schemes[j].duplicate_structures_properties.push_back(property_values[i]); //DX 20181218 - added property_values
            }
            scheme_created=true;
            break;
          }
          //cerr << "!!!!!!!!!!!!!!!!!!!!UNMATCHABLE: " << comparison_schemes[j].representative_structure_name << " and " << directory+"/"+vfiles[i] <<  endl;
        }
        if(scheme_created==false){
          StructurePrototype tmp;
          if(directory==""){
            tmp.representative_structure_name = vfiles[i];
          }
          else {
            tmp.representative_structure_name = directory+"/"+vfiles[i];  aurostd::StringSubst(tmp.representative_structure_name,"//","/"); //DX 20181003
          }
          tmp.representative_structure_generated=vstructures_generated[i];
          tmp.representative_structure_from=vstructures_from[i];
          //tmp.number_types=vxstrs[i].num_each_type.size();
          tmp.elements=vvelements[i];
          tmp.stoichiometry=vstoichs[i];
          //tmp.number_of_atoms=vxstrs[i].atoms.size();
          tmp.Pearson=vpearsons[i];
          tmp.space_group=vsgroups[i];
          tmp.grouped_Wyckoff_positions=vgrouped_Wyckoff_positions[i];
          if(property_names.size()!=0){
            tmp.property_names=property_names; //DX 20181218 - added property_names
            tmp.property_units=property_units; //DX 20181218 - added property_units
            tmp.representative_structure_properties=property_values[i]; //DX 20181218 - added property_values
          }
          if(vstructures_generated[i]){
            tmp.representative_structure=vxstrs[i];
            tmp.number_types=vxstrs[i].num_each_type.size();
            tmp.number_of_atoms=vxstrs[i].atoms.size();
            tmp.representative_structure_compound=getCompoundName(vxstrs[i]); //DX 20190111 - added compound, e.g., Ag1Br2
          }
          comparison_schemes.push_back(tmp);
        }
      }
    }
    if(LDEBUG) {
      cerr << "Prepared comparison sets: " << endl;
      stringstream ss_test;
      compare::printResults(ss_test, same_species, comparison_schemes);
      cerr << ss_test.str() << endl;
    }
    //for(uint i=0;i<comparison_schemes.size();i++){
    //      cerr << i << "duplicate_structures.size(): " << comparison_schemes[i].duplicate_structures.size() << endl;
    //      cerr << i << "duplicate_structures_generated.size(): " << comparison_schemes[i].duplicate_structures_generated.size() << endl;
    //      cerr << i << "duplicate_structures_from.size(): " << comparison_schemes[i].duplicate_structures_from.size() << endl;
    //}
  }
}

// ***************************************************************************
// createStructurePrototypes - Group structures by Pearson symbol, then space group
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> groupStructurePrototypes(vector<StructurePrototype>& structures, 
				 const bool& same_species, 
         const bool& ignore_symmetry, const bool& ignore_Wyckoff){

    // Populates the structure information into the StructurePrototype object.
    // It groups structure based on their stoichiometry, space group, and Wyckoff positions. 
    // A "representative" structure is chosen and will be compared to the 
    // possible "duplicates". The misfit values are set to -1.0 until compared.
    
    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "groupStructurePrototypes()";
   
    // variable to store structure sets to compare 
    vector<StructurePrototype> comparison_schemes;

    if(LDEBUG) {cerr << function_name << ": Number of structures to group: " << structures.size() << endl;}

    // Loop over structures. 
    // Group structures that have comparable by stoichiometry and symmetry 
    // Optional booleans control certain grouping requirements:
    //   same_species    : groups structures comprised of the same species
    //   ignore_symmetry : ignores space group when grouping (i.e., can group 
    //                     structures with different space groups)
    //   ignore_Wyckoff : ignores Wyckoff positions when grouping (i.e., can group 
    //                     structures with the same space group number, but different Wyckoff positions)

    for(uint i=0;i<structures.size(); i++){
      bool scheme_created=false;
      for(uint j=0; j<comparison_schemes.size(); j++){
        bool same_material_stoich=false;
        ostringstream tmp;
        tmp.clear();
        //DX 20190430 - this may take longer, use compound if(same_species==true && matchableSpecies(structures[i].representative_structure,comparison_schemes[j].representative_structure,same_species)==true){
        if(same_species==true && structures[i].representative_structure_compound==comparison_schemes[j].representative_structure_compound){ //DX 20190430 - quicker
          same_material_stoich=true;
        }
        else if(same_species==false){
          same_material_stoich=true;
        }
        if(same_material_stoich==true && structures[i].stoichiometry==comparison_schemes[j].stoichiometry && 
           ((ignore_symmetry && ignore_Wyckoff) ||    
            (!ignore_symmetry && ignore_Wyckoff &&
             structures[i].Pearson == comparison_schemes[j].Pearson && 
             matchableSpaceGroups(structures[i].space_group,comparison_schemes[j].space_group)) || 
            (!ignore_symmetry && !ignore_Wyckoff &&
           structures[i].Pearson == comparison_schemes[j].Pearson && 
           matchableSpaceGroups(structures[i].space_group,comparison_schemes[j].space_group) &&
             matchableWyckoffPositions(structures[i].grouped_Wyckoff_positions, comparison_schemes[j].grouped_Wyckoff_positions,same_species)))){
          if(same_species==false){
            for(uint e=0;e<structures[i].elements.size();e++){
              bool already_in=false;
              for(uint f=0;f<comparison_schemes[j].elements.size();f++){
                if(structures[i].elements[e]==comparison_schemes[j].elements[f]){
                  already_in=true;
                  break;
                }
              }
              if(already_in==false){
                comparison_schemes[j].elements.push_back(structures[i].elements[e]);
              }
            }
          }
          comparison_schemes[j].addStructurePrototypeAsDuplicate(structures[i]);
          scheme_created=true;
          break;
        }
      }
      if(scheme_created==false){
        StructurePrototype tmp = structures[i];
        comparison_schemes.push_back(tmp);
      }
    }
    if(LDEBUG) {
      cerr << "Prepared comparison sets: " << endl;
      stringstream ss_test;
      compare::printResults(ss_test, same_species, comparison_schemes);
      cerr << ss_test.str() << endl;
    }
    // DEBUG for(uint i=0;i<comparison_schemes.size();i++){
    // DEBUG  cerr << i << "duplicate_structures.size(): " << comparison_schemes[i].duplicate_structures.size() << endl;
    // DEBUG  cerr << i << "duplicate_structures_names.size(): " << comparison_schemes[i].duplicate_structures_names.size() << endl;
    // DEBUG  cerr << i << "duplicate_structures_generated.size(): " << comparison_schemes[i].duplicate_structures_generated.size() << endl;
    // DEBUG  cerr << i << "duplicate_structures_from.size(): " << comparison_schemes[i].duplicate_structures_from.size() << endl;
    // DEBUG }
    return comparison_schemes;
  }
}

/*
// ***************************************************************************
// removeDuplicateCompounds 
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& comparison_schemes, 
				 vector<xstructure>& vxstrs, const bool& same_species, 
				 const vector< vector<string> >& vvelements,
				 vector< vector<uint> >& vstoichs, vector<string>& vpearsons, 
				 vector<uint>& vsgroups, 
         vector<vector<GroupedWyckoffPosition> >& vgrouped_Wyckoff_positions,
         const string& directory, const vector<string>& vfiles, 
         vector<bool>& vstructures_generated,
         vector<string>& vstructures_from,
         const bool& ignore_symmetry, const bool& ignore_Wyckoff,
         const bool& structures_generated){
    vector<string> property_names, property_units;
    vector<vector<string> > property_values;
    return compareMultipleStructures(comparison_schemes, vxstrs, same_species, 
                                     vvelements, vstoichs, vpearsons, vsgroups, 
                                     vgrouped_Wyckoff_positions, property_names, property_units, property_values, 
                                     directory, vfiles,
                                     vstructures_generated, vstructures_from,
                                     ignore_symmetry, ignore_Wyckoff, structures_generated);
  }
}

// ***************************************************************************
// 
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> compareMultipleStructures(vector<StructurePrototype>& comparison_schemes, 
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
         const bool& structures_generated){

    string function_name = "compare::compareMultipleStructures()";
    ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;

    message << "Grouping sets of comparisons.";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    // === Organize into objects based on stoichiometry and symmetry (Pearson and space group)

    cerr << "vstructures_generated: " << vstructures_generated.size() << endl;
    cerr << "vstructures_from: " << vstructures_from.size() << endl;
    compare::createStructurePrototypes(comparison_schemes, vxstrs, same_species, 
                                       vvelements, vstoichs, vpearsons, vsgroups, 
                                       vgrouped_Wyckoff_positions, directory, vfiles,
                                       vstructures_generated, vstructures_from,
                                       ignore_symmetry, ignore_Wyckoff, structures_generated);
   
    // === If an ICSD comparison, make minimum ICSD number as the representative prototype === // 
    if(ICSD_comparison){
      compare::representativePrototypeForICSDRuns(comparison_schemes);
    }

    message << "Running comparisons ...";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes = compare::runComparisonScheme(num_proc, comparison_schemes, same_species, scale_volume, optimize_match, single_comparison_round, structures_generated, ICSD_comparison, oss); 
    
    if(final_prototypes.size()==0){
      return oss.str();
    }
    comparison_schemes.clear();
 
    // ========== Check final_prototypes ========== //
    // It is possible that two prototypes are the same regardless 
    // of having different space groups.
    //DX - BETA TESTING - compare::checkPrototypes(num_proc,same_species,final_prototypes);
 
    message << "Number of unique prototypes: " << final_prototypes.size() << " (out of " << vxstrs.size() << " structures).";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_COMPLETE_);
   
    return final_prototypes; 
  }
}
*/

// ***************************************************************************
// removeDuplicateCompounds 
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> compareDuplicateCompounds(vector<StructurePrototype>& prototype_schemes, uint& num_proc, 
                                                      bool& ICSD_comparison, ostringstream& oss){
    string function_name = "compare::compareDuplicateCompounds()";
    ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;

    vector<StructurePrototype> duplicate_check_schemes;
    for(uint i=0;i<prototype_schemes.size();i++){
      if(prototype_schemes[i].duplicate_structures_names.size()>0){ //if 0, none to check; 
        vector<StructurePrototype> tmp = compare::createComparisonSchemeForDuplicateCompounds(prototype_schemes[i]);
        duplicate_check_schemes.insert(duplicate_check_schemes.end(), tmp.begin(), tmp.end());
      }
    }
    
    if(ICSD_comparison){
      compare::representativePrototypeForICSDRuns(duplicate_check_schemes);
    }

    // set options
    bool same_species=true;
    bool scale_volume=false;
    bool optimize_match=false;
    bool single_comparison_round=false;
    bool clean_unmatched=true; //DX 20190504
    //bool structures_generated=false;

    message << "Running comparisons to remove duplicate compounds ...";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes_reduced = compare::runComparisonScheme(num_proc, duplicate_check_schemes, same_species, scale_volume, optimize_match, single_comparison_round, clean_unmatched, ICSD_comparison, oss, FileMESSAGE); //DX 20190319 - added FileMESSAGE 

    return final_prototypes_reduced;

  }
}

// ***************************************************************************
// createComparisonSchemeForDuplicateCompounds 
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> createComparisonSchemeForDuplicateCompounds(StructurePrototype& prototype_scheme){ 
                                                      
    // Group structures by similar compounds.
    // The incoming prototype is assumed to be grouped by stoichiometry, space group, and Wyckoff positions,  
    // so there is not need to check.

    vector<StructurePrototype> duplicate_check_schemes;
    for(uint i=0;i<prototype_scheme.duplicate_structures_names.size();i++){
      bool scheme_created=false;
      //cerr << i << endl;
      if(i==0){
        StructurePrototype tmp = prototype_scheme;
        tmp.representative_structure_name = prototype_scheme.duplicate_structures_names[i];
        tmp.representative_structure_compound = prototype_scheme.duplicate_structures_compounds[i];
        if(prototype_scheme.duplicate_structures_properties.size()>0){
          tmp.representative_structure_properties = prototype_scheme.duplicate_structures_properties[i];
        }
        else {
          tmp.representative_structure_properties.clear();
        }
        tmp.representative_structure.Clear(); //clear xstructure;
        tmp.representative_structure_generated=false;
        tmp.representative_structure_from=prototype_scheme.duplicate_structures_from[i];
        tmp.duplicate_structures_names.clear(); tmp.duplicate_structures_compounds.clear(); tmp.duplicate_structures.clear(); tmp.misfits.clear(); tmp.duplicate_structures_properties.clear();
        tmp.duplicate_structures_generated.clear(); tmp.duplicate_structures_from.clear();
        tmp.family_structures_names.clear(); tmp.family_structures.clear(); tmp.family_misfits.clear(); tmp.family_structures_properties.clear(); //DX 20190425 - added properties
        tmp.family_structures_generated.clear(); tmp.family_structures_from.clear();
        tmp.elements.clear();
        duplicate_check_schemes.push_back(tmp);
      }
      else {
        for(uint j=0; j<duplicate_check_schemes.size(); j++){
          if(prototype_scheme.duplicate_structures_compounds[i] == duplicate_check_schemes[j].representative_structure_compound){
            scheme_created=true;
            duplicate_check_schemes[j].duplicate_structures_names.push_back(prototype_scheme.duplicate_structures_names[i]);
            duplicate_check_schemes[j].duplicate_structures_compounds.push_back(prototype_scheme.duplicate_structures_compounds[i]);
            xstructure tmp_xstr;
            duplicate_check_schemes[j].duplicate_structures.push_back(tmp_xstr);
            duplicate_check_schemes[j].duplicate_structures_generated.push_back(false);
            duplicate_check_schemes[j].duplicate_structures_from.push_back(prototype_scheme.duplicate_structures_from[i]);
            if(prototype_scheme.duplicate_structures_properties.size()>0){
              duplicate_check_schemes[j].duplicate_structures_properties.push_back(prototype_scheme.duplicate_structures_properties[i]);
            }
            else {
              duplicate_check_schemes[j].duplicate_structures_properties.clear();
            }
            duplicate_check_schemes[j].misfits.push_back(-1.0);
            break;
          }
        }
        if(!scheme_created){
          StructurePrototype tmp = prototype_scheme;
          tmp.representative_structure_name = prototype_scheme.duplicate_structures_names[i];
          tmp.representative_structure_compound = prototype_scheme.duplicate_structures_compounds[i];
          if(prototype_scheme.duplicate_structures_properties.size()>0){
            tmp.representative_structure_properties = prototype_scheme.duplicate_structures_properties[i];
          }
          else {
            tmp.representative_structure_properties.clear();
          }
          tmp.representative_structure.Clear(); //clear xstructure;
          tmp.representative_structure_generated=false;
          tmp.representative_structure_from=prototype_scheme.duplicate_structures_from[i];
          tmp.duplicate_structures_names.clear(); tmp.duplicate_structures_compounds.clear(); tmp.duplicate_structures.clear(); tmp.misfits.clear(); tmp.duplicate_structures_properties.clear();
          tmp.duplicate_structures_generated.clear(); tmp.duplicate_structures_from.clear();
          tmp.family_structures_names.clear(); tmp.family_structures.clear(); tmp.family_misfits.clear(); tmp.family_structures_properties.clear(); //DX 20190425 - added properties
          tmp.family_structures_generated.clear(); tmp.family_structures_from.clear();
          duplicate_check_schemes.push_back(tmp);
        }
      }
    }
    return duplicate_check_schemes;
  }
}

// ***************************************************************************
// removeDuplicateCompounds 
// ***************************************************************************
namespace compare{
  void removeDuplicateCompounds(vector<StructurePrototype>& final_prototypes, vector<StructurePrototype>& duplicate_compound_comparisons){

    // Remove the duplicate compounds from the final prototype list

    for(uint i=0;i<final_prototypes.size();i++){
      for(uint j=0;j<final_prototypes[i].duplicate_structures_names.size();j++){
        bool remove_structure = true;
        for(uint k=0;k<duplicate_compound_comparisons.size();k++){
          if(final_prototypes[i].duplicate_structures_names[j] == duplicate_compound_comparisons[k].representative_structure_name){
            remove_structure = false;
            break;
          }
        }
        if(remove_structure){
          final_prototypes[i].duplicate_structures_names.erase(final_prototypes[i].duplicate_structures_names.begin()+j);
          final_prototypes[i].duplicate_structures_compounds.erase(final_prototypes[i].duplicate_structures_compounds.begin()+j);
          final_prototypes[i].duplicate_structures_generated.erase(final_prototypes[i].duplicate_structures_generated.begin()+j);
          final_prototypes[i].duplicate_structures_from.erase(final_prototypes[i].duplicate_structures_from.begin()+j);
          if(final_prototypes[i].property_names.size()!=0){
            final_prototypes[i].duplicate_structures_properties.erase(final_prototypes[i].duplicate_structures_properties.begin()+j);
          }
          final_prototypes[i].misfits.erase(final_prototypes[i].misfits.begin()+j);
        }
      }
    }
  }
}

// ***************************************************************************
// representativePrototypeForICSDRuns: Determine representative prototype for ICSD runs
// ***************************************************************************
namespace compare{
  bool representativePrototypeForICSDRuns(vector<StructurePrototype>& comparison_schemes){

    // For ICSD comparisons; make the structure with the smallest ICSD number the 
    // representative structure (normally the oldest).

    for(uint i=0;i<comparison_schemes.size();i++){
      if(comparison_schemes[i].duplicate_structures_names.size()){
        vector<string> ICSD_entries;
        ICSD_entries.push_back(findICSDName(comparison_schemes[i].representative_structure_name));
        for(uint j=0;j<comparison_schemes[i].duplicate_structures_names.size();j++){
          ICSD_entries.push_back(findICSDName(comparison_schemes[i].duplicate_structures_names[j]));
        }
        string min_ICSD_entry = findMinimumICSDEntry(ICSD_entries);
        if(!aurostd::substring2bool(comparison_schemes[i].representative_structure_name,min_ICSD_entry)){
          for(uint j=0;j<comparison_schemes[i].duplicate_structures_names.size();j++){
            if(aurostd::substring2bool(comparison_schemes[i].duplicate_structures_names[j],min_ICSD_entry)){
              string old_representative_ID = comparison_schemes[i].representative_structure_name;
              string old_representative_compound = comparison_schemes[i].representative_structure_compound;
              bool old_representative_generated = comparison_schemes[i].representative_structure_generated;
              string old_representative_from = comparison_schemes[i].representative_structure_from;
              uint old_representative_duplicate_count = comparison_schemes[i].number_compounds_matching_representative;
              vector<string> old_representative_properties = comparison_schemes[i].representative_structure_properties;
              comparison_schemes[i].representative_structure_name = comparison_schemes[i].duplicate_structures_names[j];  
              comparison_schemes[i].representative_structure_compound = comparison_schemes[i].duplicate_structures_compounds[j];  
              comparison_schemes[i].representative_structure_generated = comparison_schemes[i].duplicate_structures_generated[j];  
              comparison_schemes[i].representative_structure_from = comparison_schemes[i].duplicate_structures_from[j]; 
              comparison_schemes[i].number_compounds_matching_representative = comparison_schemes[i].number_compounds_matching_duplicate[j]; 
              if(old_representative_generated){
                xstructure old_representative_xstr = comparison_schemes[i].representative_structure;
                if(comparison_schemes[i].duplicate_structures_generated[j]){
                  comparison_schemes[i].representative_structure = comparison_schemes[i].duplicate_structures[j]; 
                  comparison_schemes[i].duplicate_structures[j] = old_representative_xstr;
                  comparison_schemes[i].duplicate_structures_generated[j] = old_representative_generated;
                }
                else {
                  comparison_schemes[i].representative_structure.Clear();
                  comparison_schemes[i].duplicate_structures_generated[j]=false; //cannot guarantee the rest of the vector is generated; may populate wrong index
                 //DX 20190304 - should I also clear out vector<xvstructure>?
                }
              }
              else {
                comparison_schemes[i].duplicate_structures_generated[j] = old_representative_generated;
              }
              if(comparison_schemes[i].duplicate_structures_properties.size()>0){
                comparison_schemes[i].representative_structure_properties = comparison_schemes[i].duplicate_structures_properties[j];
              }
              else {
                comparison_schemes[i].representative_structure_properties.clear();
              }
              comparison_schemes[i].duplicate_structures_names[j] = old_representative_ID;
              //DX 20190304 - moved into if statment up above - comparison_schemes[i].duplicate_structures[j] = old_representative_xstr;
              comparison_schemes[i].duplicate_structures_compounds[j] = old_representative_compound;
              //DX 20190304 - moved into if statment up above - comparison_schemes[i].duplicate_structures_generated[j] = old_representative_generated;
              comparison_schemes[i].duplicate_structures_from[j] = old_representative_from;
              comparison_schemes[i].number_compounds_matching_duplicate[j] = old_representative_duplicate_count;
              if(comparison_schemes[i].duplicate_structures_properties.size()>0){
                comparison_schemes[i].duplicate_structures_properties[j] = old_representative_properties;
              }
              else {
                comparison_schemes[i].duplicate_structures_properties.clear();
              }
              break;
            }  
          }
        }
      }
    }
    return true;
  }
}

// ***************************************************************************
// runComparisonThreads: Runs comparison threads
// ***************************************************************************
namespace compare{
  void runComparisonThreads(vector<StructurePrototype>& comparison_schemes, 
      std::pair<uint,uint>& start_indices,
      std::pair<uint,uint>& end_indices,
      const bool& same_species, 
      const bool& scale_volume, const bool& optimize_match, 
      ostream& oss){
 
    // Run comparison thread
    // If the xstructure is not generated, it will generate a local copy for 
    // the single comparison only (prevents overwriting in the comparisons) 
 
    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "compare::runComparisonThreads()";

    uint i_min=start_indices.first; uint i_max=end_indices.first;
    uint j_min=0; uint j_max=0;

    for(uint i=i_min; i<=i_max; i++){
      xstructure representative_structure;

      // to loop properly
      if(i==i_min){
        j_min=start_indices.second;
        if(i==i_max){j_max=end_indices.second;}
        else {j_max=comparison_schemes[i].duplicate_structures_names.size();} //-1 since in loop: j<=j_max
      }
      else if(i==i_max){j_min=0; j_max=end_indices.second;}
      else {j_min=0; j_max=comparison_schemes[i].duplicate_structures_names.size();} //-1 since in loop: j<=j_max
     
      for(uint j=j_min; j<j_max; j++){

        // get representative structure 
        if(!comparison_schemes[i].representative_structure_generated){
          if(!generateStructure(comparison_schemes[i].representative_structure_name,comparison_schemes[i].representative_structure_from,representative_structure,oss)){
            cerr << function_name << "ERROR: Could not generate structure (" << comparison_schemes[i].representative_structure_name << ")." << endl;
            exit(1);
          }
        }
        else {
          representative_structure = comparison_schemes[i].representative_structure;
        }
        //if(LDEBUG) { cerr << function_name << ": Loaded representative structure = " << comparison_schemes[i].representative_structure_name << endl; }
      
        // get prototype structure 
        xstructure duplicate_structure;
        if(!comparison_schemes[i].duplicate_structures_generated[j]){
          if(!generateStructure(comparison_schemes[i].duplicate_structures_names[j],comparison_schemes[i].duplicate_structures_from[j],duplicate_structure,oss)){
            cerr << function_name << "ERROR: Could not generate structure (" << comparison_schemes[i].duplicate_structures_names[j] << ")." << endl;
            exit(1);
          }
        }
        else {
          duplicate_structure = comparison_schemes[i].duplicate_structures[j];
        }
        if(LDEBUG) { cerr << function_name << ": Loaded duplicate structure = " << comparison_schemes[i].duplicate_structures_names[j] << endl; }
        
        // call the main comparison function
        ostringstream tmp_oss; tmp_oss.clear();
        double final_misfit=-1.0;
        if(LDEBUG) { cerr << function_name << ": Comparing " << comparison_schemes[i].representative_structure_name << " and " << comparison_schemes[i].duplicate_structures_names[j] <<  endl; }
        compare::aflowCompareStructure(1, representative_structure, //num_proc -> 1 for threads (not sure how it behaves otherwise)
                                       duplicate_structure,
                                       same_species, scale_volume, optimize_match, tmp_oss, final_misfit);
        
        // store the figure of misfit
        if(LDEBUG) { cerr << function_name << ": Comparison complete, misfit = " << final_misfit << "." << endl; }
        comparison_schemes[i].misfits[j]=final_misfit;
        comparison_schemes[i].duplicate_comparison_logs.push_back(tmp_oss.str()); //DX 20190506
      }
    }
  }
}

// ***************************************************************************
// runComparisonScheme: Runs comparisons automatically 
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> runComparisonScheme(uint& num_proc, vector<StructurePrototype>& comparison_schemes, const bool& same_species, const bool& scale_volume, const bool& optimize_match, const bool& single_comparison_round, const bool& clean_unmatched, const bool& ICSD_comparison, ostream& oss, ofstream& FileMESSAGE, bool quiet){ //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched
    
    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "compare::runComparisonScheme()";

    ostream& logstream = cout;
    stringstream message;
    //DX 20190319 [OBSOLETE] ofstream FileMESSAGE;

    // print initial grouped sets of comparisons
    if(LDEBUG) {
      cerr << function_name << ": Number of comparion sets: " << comparison_schemes.size() << endl;
      stringstream ss_test;
      compare::printResults(ss_test, same_species, comparison_schemes);
      cerr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      cerr << ss_test.str() << endl;
    }

    uint number_of_comparisons = 0;
    vector<std::pair<uint,uint> > start_indices, end_indices;
    
    uint num_comparison_threads = 1;
#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
    // THREADED VERSION - START
    if(LDEBUG) { cerr << function_name << ": Threaded version." << endl; } 

    // split into threads
    number_of_comparisons = 0;
    for(uint i=0;i<comparison_schemes.size();i++){ number_of_comparisons += comparison_schemes[i].numberOfComparisons(); }
  
    num_comparison_threads = aurostd::min(num_proc,number_of_comparisons);
    splitComparisonIntoThreads(comparison_schemes, num_comparison_threads, start_indices, end_indices);

    // run threads
    vector<std::thread> threads;
    for(uint n=0; n<num_comparison_threads; n++){
	    threads.push_back(std::thread(compare::runComparisonThreads,std::ref(comparison_schemes),std::ref(start_indices[n]),std::ref(end_indices[n]),
                                    std::ref(same_species),std::ref(scale_volume),std::ref(optimize_match),
                                    std::ref(oss)));
    }        

    // join threads
    for(uint t=0;t<threads.size();t++){
	    threads[t].join();
    }

    // THREADED VERSION - END

#else
   // NON-THREADED VERISON - START
    if(LDEBUG) { cerr << function_name << ": Non-threaded version." << endl; } 

    for(uint i=0; i<comparison_schemes.size(); i++){
      xstructure representative_structure;
      for(uint j=0; j<comparison_schemes[i].duplicate_structures_names.size(); j++){
        if(j==0){
          if(!comparison_schemes[i].representative_structure_generated){
            if(!generateStructure(comparison_schemes[i].representative_structure_name,comparison_schemes[i].representative_structure_from,representative_structure,oss)){
              cerr << function_name << "ERROR: Could not generate structure (" << comparison_schemes[i].representative_structure_name << ")." << endl;
              exit(1);
            }
          }
          else {
            representative_structure = comparison_schemes[i].representative_structure;
          }
        }
        ostringstream tmp_oss;
        tmp_oss.clear();
        double final_misfit=-1.0;
        xstructure duplicate_structure;
        if(!comparison_schemes[i].duplicate_structures_generated[j]){
          if(!generateStructure(comparison_schemes[i].duplicate_structures_names[j],comparison_schemes[i].duplicate_structures_from[j],duplicate_structure,oss)){
            cerr << function_name << "ERROR: Could not generate structure (" << comparison_schemes[i].duplicate_structures_names[j] << ")." << endl;
            exit(1);
          }
        }
        else {
          duplicate_structure = comparison_schemes[i].duplicate_structures[j];
        }
        // Call the main comparison function
        compare::aflowCompareStructure(num_proc, representative_structure,
                                       duplicate_structure,
                                       same_species, scale_volume, optimize_match, tmp_oss, final_misfit);
        // Store the figure of misfit
        comparison_schemes[i].misfits[j]=final_misfit;
      }
    }
   //SINGLE THREAD - END
   
#endif

    // count the number of mismatches (i.e. mis > 0.1)
    int num_mismatches_orig=compare::numberMismatches(comparison_schemes);
    int num_mismatches=num_mismatches_orig;

    //DX 20190504 - added clean unmatched option - START
    if(!clean_unmatched && single_comparison_round){
      return comparison_schemes;
    }
    //DX 20190504 - added clean unmatched option - END

    if(num_mismatches > 0 && !single_comparison_round && !quiet){
      message << "Number of unmatched structures: " << num_mismatches << ". Continuing comparisons ...";
      pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    }

    // create new object for comparisons
    vector<StructurePrototype> final_prototypes;

    // regroup comparisons based on misfit value
    if(num_mismatches==0 && !single_comparison_round){
      compare::appendStructurePrototypes(comparison_schemes, final_prototypes, clean_unmatched, quiet); //DX 20190506 -- added clean_unmatched
    }
    
    // Loop: continue comparison until all strucutures are matched or all comparisons schemes exhaused
    while(num_mismatches!=0){
      // regroup comparisons based on misfit value
      compare::appendStructurePrototypes(comparison_schemes, final_prototypes, clean_unmatched, quiet); //DX 20190506 -- added clean unmatched

      // exit if only one round of comparison is requested
      if(single_comparison_round==true){return final_prototypes;}

      // reorder structures so minimum ICSD is the representative structure
      if(ICSD_comparison){ compare::representativePrototypeForICSDRuns(comparison_schemes); }

      // split into threads
      number_of_comparisons=0;
      for(uint i=0;i<comparison_schemes.size();i++){ number_of_comparisons += comparison_schemes[i].numberOfComparisons(); }
      num_comparison_threads = aurostd::min(num_proc,number_of_comparisons);

      if(number_of_comparisons>0){
        if(LDEBUG) { cerr << function_name << ": Number of comparisons is not zero... " << number_of_comparisons << endl; }
#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
				// THREADED VERISON - START

        // split into threads
				start_indices.clear(); end_indices.clear();
				splitComparisonIntoThreads(comparison_schemes, num_comparison_threads, start_indices, end_indices);
				threads.clear();
        // run threads
				for(uint n=0; n<num_comparison_threads; n++){
					threads.push_back(std::thread(compare::runComparisonThreads,std::ref(comparison_schemes),std::ref(start_indices[n]),std::ref(end_indices[n]),
																			std::ref(same_species),std::ref(scale_volume),std::ref(optimize_match),
																			std::ref(oss)));
				}        
        // join threads
				for(uint t=0;t<threads.size();t++){
					threads[t].join();
				}
				// THREADED VERISON - END
#else
      //SINGLE THREAD - START

			start_indices.clear(); end_indices.clear();
      uint single_thread=1;
      splitComparisonIntoThreads(comparison_schemes, single_thread, start_indices, end_indices);
			compare::runComparisonThreads(comparison_schemes,start_indices[0],end_indices[0],
																		same_species,scale_volume,optimize_match,oss);
      //SINGLE THREAD - END
#endif
      }

      // update number of mismatches
      num_mismatches_orig=num_mismatches;
      num_mismatches=compare::numberMismatches(comparison_schemes);

      if(num_mismatches > 0 && !quiet){
        message << "Number of unmatched structures: " << num_mismatches << ". Continuing comparisons ...";
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      }

      // ensure while loop is not uncontrolled
      if(num_mismatches>=num_mismatches_orig){
         oss << "compare::ERROR - Number of mismatches is increasing ... impossible (bug in comparision framework). "
              << "Contact David Hicks (d.hicks@duke.edu)." << endl;
         final_prototypes.clear();
         return final_prototypes;
      }
    }
    // end of while loop

    // append new prototype groupings
    compare::appendStructurePrototypes(comparison_schemes, final_prototypes, clean_unmatched, quiet); //DX 20190303 - added to properly clear xstructure and generated structure variable //DX 20190504 - added clean unmatched
    //DX ORIG 20190303 - final_prototypes.insert(final_prototypes.end(),comparison_schemes.begin(),comparison_schemes.end());
    return final_prototypes;
  }
}

// ***************************************************************************
// calculateDivisors 
// ***************************************************************************
namespace compare{
  vector<std::pair<uint,uint> > calculateDivisors(const uint& number){

    // Determine the divisors for a number 
    // Used to check if the number of permutations follows number theory

    vector<std::pair<uint,uint> > divisor_pairs;
    std::pair<uint,uint> tmp_pair; 
    for(uint i=1; i<=number; i++){
      if(number%i==0){
        tmp_pair.first = i; tmp_pair.second = number/i;
        divisor_pairs.push_back(tmp_pair);
      }
    }
    return divisor_pairs;   
  }
}

// ***************************************************************************
// calculateDivisors 
// ***************************************************************************
namespace compare{
  bool checkNumberOfGroupings(vector<StructurePrototype>& comparison_schemes, uint number){

    // Check if the number of comparisons for the permutation comparisons are commensurate.  
    // The number of unqiue permutations must be divisors of the total number of permutations,
    // e.g., 6 permutations -> 1 unique in a group of 6; 2 unique in groups of 3;
    // 3 unique in groups of 2; or 6 unique in groups of 1.
    // If this is not the case, then there is something wrong with the comparison result.

    vector<std::pair<uint,uint> > divisor_pairs = calculateDivisors(number);
    uint divisor_index = 0;
    for(uint i=0;i<divisor_pairs.size();i++){
      if(comparison_schemes.size() == divisor_pairs[i].first){
        divisor_index = i;
        break;
      }
    }
    uint num_consistent = 0;
    for(uint j=0; j<comparison_schemes.size(); j++){
      if(comparison_schemes[j].duplicate_structures_names.size()+1 == divisor_pairs[divisor_index].second){ //+1 to include representative
        num_consistent++;
      }
    }
    if(num_consistent != comparison_schemes.size()){
      return false;
    }
    return true;
  }
}

// ***************************************************************************
// createStructurePermutations - Group structures by Pearson symbol, then space group
// ***************************************************************************
namespace compare{
  void createStructurePermutations(vector<StructurePrototype>& comparison_schemes, const vector<vector<string> >& name_order,
                                   vector<vector<GroupedWyckoffPosition> >& permutation_grouped_Wyckoff_positions,
				   const vector<xstructure>& vxstrs, const bool& same_species){

    // Populates the structure information into the StructurePrototype object.
    // It groups structure based on their stoichiometry and pearson symbol and 
    // space group. A "representative" structure is chosen and will be compared to the 
    // possible "duplicates". The misfit values are set to -1.0 until compared.

    // First, separate by stoichiometry
    //LDEBUG cerr << "name_order.size(): " << name_order.size() << endl;
    for(uint i=0;i<name_order.size(); i++){
      bool scheme_created=false;
      string name = "";
      for(uint j=0;j<name_order[i].size();j++){name+=name_order[i][j];}
      if(i==0){
        StructurePrototype tmp;
        tmp.representative_structure_name=name;
        tmp.representative_structure=vxstrs[i];
        tmp.representative_structure_generated=true;
        tmp.representative_structure_from="input";
        tmp.number_types=vxstrs[i].num_each_type.size();
        tmp.number_of_atoms=vxstrs[i].atoms.size();
        tmp.grouped_Wyckoff_positions=permutation_grouped_Wyckoff_positions[i];
        //cerr << "tmp: " << tmp << endl;
        comparison_schemes.push_back(tmp);
      }
      else {
        //cerr << comparison_schemes.size() << endl;
        for(uint j=0; j<comparison_schemes.size(); j++){
          bool same_material_stoich=false;
          ostringstream tmp;
          tmp.clear();
          if(same_species==true && matchableSpecies(vxstrs[i],comparison_schemes[j].representative_structure,same_species)==true){
            same_material_stoich=true;
          }
          else if(same_species==false){
            same_material_stoich=true;
          }
          if(same_material_stoich==true &&
             matchableWyckoffPositions(permutation_grouped_Wyckoff_positions[i], comparison_schemes[j].grouped_Wyckoff_positions,same_species)){
            comparison_schemes[j].duplicate_structures_names.push_back(name);
            comparison_schemes[j].duplicate_structures.push_back(vxstrs[i]);
            comparison_schemes[j].duplicate_structures_generated.push_back(true);
            comparison_schemes[j].duplicate_structures_from.push_back("input");
            comparison_schemes[j].misfits.push_back(-1.0);
            scheme_created=true;
            break;
          }
        }
        if(scheme_created==false){
          StructurePrototype tmp;
          tmp.representative_structure_name=name;
          tmp.representative_structure_generated=true;
          tmp.representative_structure_from="input";
          tmp.representative_structure=vxstrs[i];
          tmp.number_types=vxstrs[i].num_each_type.size();
          tmp.number_of_atoms=vxstrs[i].atoms.size();
          tmp.grouped_Wyckoff_positions=permutation_grouped_Wyckoff_positions[i];
          comparison_schemes.push_back(tmp);
        }
      }
    }
    // Check number of sets and elements in set are consistent with number theory
    // (Groupings/sets must be a divisor of the total number of permutations)
    if(!checkNumberOfGroupings(comparison_schemes,name_order.size())){
      cerr << "ERROR: Initial groupings of permutations do not follow number theory. Contact David Hicks (d.hicks@duke.edu). Exiting." << endl;
      exit(1);
    }
  }
}

// ***************************************************************************
// makeRepresentativeEvenPermutation
// ***************************************************************************
namespace compare{
  bool makeRepresentativeEvenPermutation(vector<StructurePrototype>& comparison_schemes, const vector<vector<string> >& name_order){ 
    // Make sure the even permutation is the representative.  
    // If there are multiple even permutations in a given set of comparisons,
    // default to the mininum even permutation. (DX, may want to change default)

    for(uint i=0; i<comparison_schemes.size(); i++){	
      //Find representative permutation number
      uint representative_permutation_num = 0;
      for(uint j=0; j<name_order.size(); j++){			
        string name="";
        for(uint k=0;k<name_order[j].size();k++){name+=name_order[j][k];}
        if(comparison_schemes[i].representative_structure_name == name){
          representative_permutation_num = j;
          break;
        }
      }
      //Find proto (potential duplicate) permutation number
      uint min_even_duplicate_permutation_num = 1e9;
      uint duplicate_permutation_num = 0;
      uint min_duplicate_index = 0;
      for(uint p=0; p<comparison_schemes[i].duplicate_structures_names.size(); p++){
        for(uint j=0; j<name_order.size(); j++){			
          string name="";
          for(uint k=0;k<name_order[j].size();k++){name+=name_order[j][k];}
          if(comparison_schemes[i].duplicate_structures_names[p] == name){
            duplicate_permutation_num = j;
            // Check if proto is minimum/even permutation
            if(duplicate_permutation_num < min_even_duplicate_permutation_num && duplicate_permutation_num%2 == 0){ 
              min_duplicate_index = p;     
              min_even_duplicate_permutation_num = duplicate_permutation_num;
            }
          }
        }
      }
      // If representative permutation is already even, only swap representative and proto if proto permutation is less and even
      // Else replace automatically if proto permutation is not 1e9
      if((representative_permutation_num%2 == 0 && min_even_duplicate_permutation_num < representative_permutation_num) ||
         (representative_permutation_num%2 != 0 && min_even_duplicate_permutation_num != 1e9)){
        comparison_schemes[i].duplicate_structures_names.push_back(comparison_schemes[i].representative_structure_name);  
        comparison_schemes[i].duplicate_structures.push_back(comparison_schemes[i].representative_structure);  
        comparison_schemes[i].duplicate_structures_generated.push_back(comparison_schemes[i].representative_structure_generated);  
        comparison_schemes[i].duplicate_structures_from.push_back(comparison_schemes[i].representative_structure_from);  
        comparison_schemes[i].representative_structure_name = comparison_schemes[i].duplicate_structures_names[min_duplicate_index];  
        comparison_schemes[i].representative_structure = comparison_schemes[i].duplicate_structures[min_duplicate_index];  
        comparison_schemes[i].representative_structure_generated = comparison_schemes[i].duplicate_structures_generated[min_duplicate_index];  
        comparison_schemes[i].representative_structure_from = comparison_schemes[i].duplicate_structures_from[min_duplicate_index];  
        comparison_schemes[i].duplicate_structures_names.erase(comparison_schemes[i].duplicate_structures_names.begin()+min_duplicate_index);
        comparison_schemes[i].duplicate_structures.erase(comparison_schemes[i].duplicate_structures.begin()+min_duplicate_index);
        comparison_schemes[i].duplicate_structures_generated.erase(comparison_schemes[i].duplicate_structures_generated.begin()+min_duplicate_index);
        comparison_schemes[i].duplicate_structures_from.erase(comparison_schemes[i].duplicate_structures_from.begin()+min_duplicate_index);
      }
    }
    return true;
  } 
}

// ***************************************************************************
// NumberMismaches - Count the number of non-matches
// ***************************************************************************
namespace compare{
  int numberMismatches(const vector<StructurePrototype> comparison_schemes){

    // Count the number of comparisons that have a misfit greater than 0.1
    // i.e., not matching

    int num_mismatches=0;
    for(uint i=0; i<comparison_schemes.size(); i++){
      for(uint j=0; j<comparison_schemes[i].misfits.size(); j++){
        if(comparison_schemes[i].misfits[j] > 0.1 || comparison_schemes[i].misfits[j] == -1.0 ){
          num_mismatches+=1;
        }
      }
    }
    return num_mismatches;
  }
}

// ***************************************************************************
// appendStructurePrototypes - Create new structure prototypes after comparisons
// ***************************************************************************
namespace compare{
  void appendStructurePrototypes(vector<StructurePrototype>& comparison_schemes, 
				 vector<StructurePrototype>& final_prototypes,
         bool clean_unmatched, //DX 20190506
         bool quiet){

    // This cleans the StrucuturePrototype objects by removing all the mismatches.
    // Then, it takes the mismatches and makes them into new StructurePrototype objects
    // to be compared.
      
    //LDEBUG stringstream ss_test;
    //LDEBUG compare::printResults(ss_test, true, comparison_schemes);
    //LDEBUG cerr << ss_test.str() << endl;

    ostringstream oss;
    ostream& logstream = cout;
    stringstream message;
    ofstream FileMESSAGE;
    string function_name = "compare::appendStructurePrototypes()";

    vector<StructurePrototype> tmp_list;
    for(uint i=0; i<comparison_schemes.size(); i++){
      bool first_mismatch=true;
      for(uint j=0; j<comparison_schemes[i].misfits.size(); j++){
        //DX TEST if(comparison_schemes[i].misfits[j] > 0.1 || comparison_schemes[i].misfits[j] == -1.0 ){
        if(comparison_schemes[i].misfits[j] > 0.1 || comparison_schemes[i].misfits[j]+1.0 <_ZERO_TOL_){ //DX 20190301 - more robust
          // First, store any family prototype information
          if(comparison_schemes[i].misfits[j] > 0.1 && comparison_schemes[i].misfits[j] <= 0.2){
            comparison_schemes[i].family_structures_names.push_back(comparison_schemes[i].duplicate_structures_names[j]);
            comparison_schemes[i].family_misfits.push_back(comparison_schemes[i].misfits[j]);
            //DX 20190424 - store properties of family structures - START
            if(comparison_schemes[i].property_names.size()!=0){
              comparison_schemes[i].family_structures_properties.push_back(comparison_schemes[i].duplicate_structures_properties[j]);
            }
            //DX 20190424 - store properties of family structures - END
          }
          // Take first mismatch and make as the representative structure in the new object
          if(first_mismatch==true){
            StructurePrototype tmp;
            tmp.copyPrototypeInformation(comparison_schemes[i]);
            tmp.putDuplicateAsRepresentative(comparison_schemes[i],j);
            tmp_list.push_back(tmp);
            if(clean_unmatched){ comparison_schemes[i].removeNonDuplicate(j); j--; } //DX 20190504 - put in if-statement
            first_mismatch=false;
          }
          // If not the first mismatch, add as a proto structure in the new object
          else if(first_mismatch==false){
            tmp_list.back().copyDuplicate(comparison_schemes[i],j);
            if(clean_unmatched){ comparison_schemes[i].removeNonDuplicate(j); j--; } //DX 20190504 - put in if-statement
          }
        }
        //DX 20181220 - if they are matched, then we should delete the xstructure, since we no longer need the structure (save memory)
        else if(comparison_schemes[i].misfits[j] <= 0.1 && !std::signbit(comparison_schemes[i].misfits[j])){
          //comparison_schemes[i].duplicate_structures.erase(comparison_schemes[i].duplicate_structures.begin()+j);
          // check if the structure is generated first
          if(comparison_schemes[i].duplicate_structures_generated[j]){comparison_schemes[i].duplicate_structures[j].Clear(); comparison_schemes[i].duplicate_structures_generated[j]=false; } //DX 20190303 - update generated flag
        }
      }
      //DX 20181220 - can clear representative structure since we will no longer compare it (save memory)
      comparison_schemes[i].representative_structure.Clear();
      comparison_schemes[i].representative_structure_generated=false;

      // if not quiet, print the comparison results to the screen 
      // (useful for long comparison times or if the program terminates early)
      if(!quiet){
        message << "Identified unique prototype: " << endl;
        message << "   prototype=" << comparison_schemes[i].representative_structure_name << endl;
        if(comparison_schemes[i].duplicate_structures_names.size()==0){
          message << "   No duplicates. " << endl;
        }
        else {
          message << "   " << setw(80) << std::left << "List of duplicates"
                  << setw(15) << std::left << "misfit value" << endl;
          message << "   " << setw(80) << std::left 
                  << "-----------------------------------------------------------------------------------------------" << endl;
          for(uint d=0;d<comparison_schemes[i].duplicate_structures_names.size();d++){
            message << "   " << setw(80) << std::left << comparison_schemes[i].duplicate_structures_names[d]
                    << setw(15) << std::left << comparison_schemes[i].misfits[d] << endl;
          }
        }
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
      }
      
      // Store finished (already compared) schemes in final_prototypes
      final_prototypes.push_back(comparison_schemes[i]);
      comparison_schemes.erase(comparison_schemes.begin()+i);
      i--;
    }
    // Store newly generated schemes (not compared yet) into comparison_schemes
    comparison_schemes=tmp_list;
  }
}

// ***************************************************************************
// checkPrototypes - Ensure prototypes of different SG are compared
// ***************************************************************************
namespace compare{
  void checkPrototypes(const uint& num_proc, const bool& same_species, 
		       vector<StructurePrototype>& final_prototypes){

    // Checks to see if prototypes of different space groups are similar. 
    // If they are, then we combine the StructurePrototype objects into one.
    // When combining, we keep the "representative" prototype as the one with a higher 
    // symmetry (i.e. higher space group).
    // This is an optional function; we may not want to do this when 
    // comparing prototypes or comparing material properties 

    for(uint i=0;i<final_prototypes.size();i++){
      vector<int> store_indices;
      vector<double> store_misfits;
      int min_index=-1;
      double min_misfit=-1.0;
      for(uint j=i;j<final_prototypes.size();j++){
        ostringstream tmp_oss;
        tmp_oss.clear();
        if(
           // If same_species==true
           (same_species==true && 
            matchableSpecies(final_prototypes[i].representative_structure,final_prototypes[j].representative_structure,same_species)==true &&
            final_prototypes[i].stoichiometry==final_prototypes[j].stoichiometry && 
            final_prototypes[i].Pearson==final_prototypes[j].Pearson &&
            //final_prototypes[i].space_group!=final_prototypes[j].space_group) ||
            !matchableSpaceGroups(final_prototypes[i].space_group,final_prototypes[j].space_group)) ||
           // If same_species==false
           (same_species==false && 
            final_prototypes[i].stoichiometry==final_prototypes[j].stoichiometry &&
            final_prototypes[i].Pearson==final_prototypes[j].Pearson && 
            //final_prototypes[i].space_group!=final_prototypes[j].space_group)
            !matchableSpaceGroups(final_prototypes[i].space_group,final_prototypes[j].space_group))
          ){
          double final_misfit=-1.0;
          bool scale_volume=true; //default is true
          bool optimize_match=false; //default is false
          aflowCompareStructure(num_proc, final_prototypes[i].representative_structure, 
					 final_prototypes[j].representative_structure, same_species, 
					 scale_volume, optimize_match, tmp_oss, final_misfit);
          if(final_misfit < min_misfit){
             min_misfit=final_misfit;
             min_index=j;
          }
        }
      }
      // If one prototype is similar to another, add to one with higher space group
      if(min_misfit!=-1.0){
        int sg_ind=-1;
        uint other_ind=-1;
        if(final_prototypes[i].space_group > final_prototypes[min_index].space_group){
          sg_ind=i;
          other_ind=min_index;
        }
        else {
          sg_ind=min_index;
          other_ind=i;
        }
        // Transfer info to prototype with higher space group
        final_prototypes[sg_ind].duplicate_structures_names.push_back(final_prototypes[other_ind].representative_structure_name);
        final_prototypes[sg_ind].duplicate_structures_names.insert(final_prototypes[sg_ind].duplicate_structures_names.end(),
				final_prototypes[other_ind].duplicate_structures_names.begin(),
				final_prototypes[other_ind].duplicate_structures_names.end());
        final_prototypes[sg_ind].duplicate_structures.push_back(final_prototypes[other_ind].representative_structure);
        final_prototypes[sg_ind].duplicate_structures.insert(final_prototypes[sg_ind].duplicate_structures.end(),
				final_prototypes[other_ind].duplicate_structures.begin(),
				final_prototypes[other_ind].duplicate_structures.end());
        final_prototypes[sg_ind].duplicate_structures_compounds.push_back(final_prototypes[other_ind].representative_structure_compound);
        final_prototypes[sg_ind].duplicate_structures_compounds.insert(final_prototypes[sg_ind].duplicate_structures_compounds.end(),
				final_prototypes[other_ind].duplicate_structures_compounds.begin(),
				final_prototypes[other_ind].duplicate_structures_compounds.end());
        final_prototypes[sg_ind].duplicate_structures_generated.push_back(final_prototypes[other_ind].representative_structure_generated);
        final_prototypes[sg_ind].duplicate_structures_generated.insert(final_prototypes[sg_ind].duplicate_structures_generated.end(),
				final_prototypes[other_ind].duplicate_structures_generated.begin(),
				final_prototypes[other_ind].duplicate_structures_generated.end());
        final_prototypes[sg_ind].duplicate_structures_from.push_back(final_prototypes[other_ind].representative_structure_from);
        final_prototypes[sg_ind].duplicate_structures_from.insert(final_prototypes[sg_ind].duplicate_structures_from.end(),
				final_prototypes[other_ind].duplicate_structures_from.begin(),
				final_prototypes[other_ind].duplicate_structures_from.end());
        // Delete the prototype with the lower space group
        final_prototypes.erase(final_prototypes.begin()+other_ind);
        // If the index deleted was less than the initial loop (i), then need to reduce iterator
        if(other_ind<=i){
          i--;
        }
      }
    }
  }
}

// ***************************************************************************
// printResults - Displays results for .txt file
// ***************************************************************************
namespace compare{
  void printResults(ostream& ss_out, const bool& same_species, 
			 const vector<StructurePrototype>& final_prototypes,
       string mode){

    // Print the comparison results in either a JSON or TXT format
    // In general, the structure along with the misfit value is printed
    // If material properties are provided, then the properties will be displayed 
    // next to the misfit value

    // JSON MODE
    if(aurostd::tolower(mode)=="json"){ //case insensitive
      ss_out << "[" << endl;
      for(uint j=0; j<final_prototypes.size(); j++){
         ss_out << final_prototypes[j];
         if(j!=final_prototypes.size()-1){
           ss_out << "," << endl;
         }
      }
      ss_out << endl << "]" << endl;
    }

    if(aurostd::tolower(mode)=="txt" || aurostd::tolower(mode)=="text"){ //case insensitive
      // TXT MODE
      int indent_spacing = 2;
      int structure_spacing = 80; // structure name spacing
      int misfit_spacing = 15;
      int property_spacing = 35;

      // Displays comparison information in a TXT file.

      for(uint j=0; j<final_prototypes.size(); j++){
        //ss_out << std::string(indent_spacing+structure_spacing+misfit_spacing, '=') << endl;
        ss_out << std::string(indent_spacing+structure_spacing+misfit_spacing, '=');
        if(final_prototypes[j].property_names.size()!=0){
          ss_out << std::string(final_prototypes[j].property_names.size()*property_spacing, '=');
        }
        ss_out << endl;
        ss_out << "# ";
        if(same_species==true){
          //DX 2019090311 [OBSOLETE] for(uint k=0;k<final_prototypes[j].elements.size();k++){
          //DX 2019090311 [OBSOLETE]   ss_out << final_prototypes[j].elements[k] << final_prototypes[j].stoichiometry[k];
          //DX 2019090311 [OBSOLETE] }
          //LDEBUG cerr << j << " compound: " << final_prototypes[j].representative_structure_compound; //DX 20190311
          ss_out << final_prototypes[j].representative_structure_compound; //DX 20190311
          ss_out << "  SG=#" << final_prototypes[j].space_group;
          // ORIG ss_out << "  Wyckoffs=" << compare::printWyckoffString(final_prototypes[j].grouped_Wyckoff_positions,true) << endl;
          ss_out << "  Wyckoffs=" << compare::printWyckoffString(final_prototypes[j].grouped_Wyckoff_positions,true); //DX 20190228 - remove count
          uint number_of_duplicates = final_prototypes[j].numberOfDuplicates(); //DX 20190506 - made function
          ss_out << "  duplicate_compounds=" << number_of_duplicates << endl; //DX 20190228 - add count
        if(final_prototypes[j].representative_structure_properties.size()!=0){
          ss_out << "  " << setw(structure_spacing) << std::left << "structure";
          ss_out << setw(misfit_spacing) << std::right << "misfit";
          for(uint l=0;l<final_prototypes[j].property_names.size();l++){
            if(final_prototypes[j].property_units[l].size()!=0){
              ss_out << setw(property_spacing) << std::right 
                << final_prototypes[j].property_names[l]+"("+final_prototypes[j].property_units[l]+")";
            }
            else {ss_out << setw(property_spacing) << std::right << final_prototypes[j].property_names[l];}
          }
          ss_out << endl;
        }
        if(final_prototypes[j].representative_structure_properties.size()!=0){
        ss_out << std::string(indent_spacing+structure_spacing+misfit_spacing, '-');
        if(final_prototypes[j].property_names.size()!=0){
          ss_out << std::string(final_prototypes[j].representative_structure_properties.size()*property_spacing, '-');
          ss_out << std::string(final_prototypes[j].property_names.size()*property_spacing, '-');
        }
        ss_out << endl;
        }
          ss_out << "  " << setw(structure_spacing) << std::left << "prototype="+final_prototypes[j].representative_structure_name;
          if(final_prototypes[j].representative_structure_properties.size()!=0){
            ss_out << setw(misfit_spacing) << std::right << "-";
            for(uint l=0;l<final_prototypes[j].representative_structure_properties.size();l++){
              ss_out << setw(property_spacing) << std::right << final_prototypes[j].representative_structure_properties[l];
            }
          }
          //ss_out << endl;
        }
        else if(same_species==false){
          for(uint k=0;k<final_prototypes[j].stoichiometry.size();k++){
            if(k==0){
              ss_out << final_prototypes[j].stoichiometry[k];
            }
            else {
              ss_out << ":" << final_prototypes[j].stoichiometry[k];
            }
          }
          ss_out << "  SG=#" << final_prototypes[j].space_group;
          ss_out << "  Wyckoffs=" << compare::printWyckoffString(final_prototypes[j].grouped_Wyckoff_positions,true);;
          uint number_of_duplicates = final_prototypes[j].numberOfDuplicates(); //DX 20190506 - made function
          ss_out << "  duplicate_structures=" << number_of_duplicates; //DX 20190228 - add count
          uint number_duplicate_compounds = 0;
          for(uint k=0;k<final_prototypes[j].number_compounds_matching_duplicate.size();k++){
            number_duplicate_compounds+=final_prototypes[j].number_compounds_matching_duplicate[k];
          }
          number_duplicate_compounds+= number_of_duplicates+final_prototypes[j].number_compounds_matching_representative; //DX 20190321 - need to update variable, otherwise may not enter if statement
          if(number_duplicate_compounds!=0){
            ss_out << "  duplicate_compounds=" << number_duplicate_compounds; //DX 20190228 - add count
          }
          ss_out << endl;
          ss_out << "  " << setw(structure_spacing) << std::left << "prototype="+final_prototypes[j].representative_structure_name;
          // perhaps add which permutations are duplicates
          if(final_prototypes[j].unique_permutations.size()!=0){
            ss_out << endl << "  " << setw(structure_spacing) << std::left << "unique permutations="+aurostd::joinWDelimiter(final_prototypes[j].unique_permutations,",");
          }
        }
        ss_out << endl;
        ss_out << std::string(indent_spacing+structure_spacing+misfit_spacing, '-');
        if(final_prototypes[j].property_names.size()!=0){
          ss_out << std::string(final_prototypes[j].property_names.size()*property_spacing, '-');
        }
        ss_out << endl;
        if(final_prototypes[j].duplicate_structures_names.size()!=0 && final_prototypes[j].representative_structure_properties.size()==0){
          ss_out << "  " << setw(structure_spacing) << std::left << "list of duplicates";
          ss_out << setw(misfit_spacing) << std::right << "misfit";
        }
        else if(final_prototypes[j].duplicate_structures_names.size()==0){
          ss_out << "  " << setw(structure_spacing) << std::left << "no duplicates";
        }
        if(final_prototypes[j].representative_structure_properties.size()==0){
          for(uint l=0;l<final_prototypes[j].property_names.size();l++){
            if(final_prototypes[j].property_units[l].size()!=0){
              ss_out << setw(property_spacing) << std::right 
                << final_prototypes[j].property_names[l]+"("+final_prototypes[j].property_units[l]+")";
            }
            else {ss_out << setw(property_spacing) << std::right << final_prototypes[j].property_names[l];}
          }
          ss_out << endl;
          ss_out << std::string(indent_spacing+structure_spacing+misfit_spacing, '-');
          ss_out << std::string(final_prototypes[j].property_names.size()*property_spacing, '-');
          ss_out << endl;
        }
        if(final_prototypes[j].duplicate_structures_names.size()>0){
        for(uint k=0;k<final_prototypes[j].duplicate_structures_names.size();k++){
           ss_out << "  " << setw(structure_spacing) << std::left << final_prototypes[j].duplicate_structures_names[k] 
                  << setw(misfit_spacing) << std::right << final_prototypes[j].misfits[k];
           if(final_prototypes[j].property_names.size()!=0){
             for(uint l=0;l<final_prototypes[j].duplicate_structures_properties[k].size();l++){
               ss_out << setw(property_spacing) << std::right << final_prototypes[j].duplicate_structures_properties[k][l];
             }   
           }     
           ss_out << endl;
        }
        }
        else {
          ss_out << endl;
        }
      }
    }
  }
}

// ***************************************************************************
// sameStoichiometry
// ***************************************************************************
namespace compare{
  bool sameStoichiometry(const vector<uint>& stoich1, const vector<uint>& stoich2){

    // Determine if two stoichiometries are equivalent
    // Stoichiometries must be in the same order to match

    // quick check
    if(stoich1.size()!=stoich2.size()){return false;}
    
    for(uint i=0;i<stoich1.size();i++){
      if(stoich1[i]!=stoich2[i]){
        return false;
      }
    }
    return true;
  }
}

// ***************************************************************************
// Matchable species
// ***************************************************************************
namespace compare{
  bool matchableSpecies(const xstructure& xstr1, const xstructure& xstr2, 
			const bool& same_species){

    // Determine if it is possible to match up species based on the number of atoms
    // (i.e, reduced stoichiometries are equal)

    bool LDEBUG=(false || XHOST.DEBUG);
    vector<uint> stoich1;
    vector<uint> stoich2;
    bool matchable=true;
    if(xstr1.species.size()==xstr2.species.size()){
      if(xstr1.species.size()==1){
     	  stoich1.push_back(1); stoich2.push_back(1);
      }
      else {
        stoich1=gcdStoich(xstr1.num_each_type);
        stoich2=gcdStoich(xstr2.num_each_type);
      }
      uint matches=0;
      // Check if we can match to same species (atoms and stoichs)
      if(same_species==true){
        bool commensurate=false;
        for(uint i=0; i<stoich1.size(); i++){
          for(uint j=0; j<stoich2.size(); j++){
            if(stoich1[i]==stoich2[j] && KBIN::VASP_PseudoPotential_CleanName(xstr1.species[i])==KBIN::VASP_PseudoPotential_CleanName(xstr2.species[j])){ //DX 20190329 - remove pseudopotential information
              //cerr << "matching: " << stoich1[i] << "==" << stoich2[j] << " && " << xstr1.species[i] << "==" << xstr2.species[j] << endl;
              matches++;
              commensurate=true;
              break;
            }
	  }
          if(commensurate==false){
            matchable=false;
            break;
          }
	      }
      }
      // Check if we can match stoichs only
      else {
        for(uint i=0; i<stoich1.size(); i++){
          std::sort(stoich1.begin(),stoich1.end());
          std::sort(stoich2.begin(),stoich2.end());
        }
        for(uint i=0; i<stoich1.size(); i++){
          if(stoich1[i]!=stoich2[i]){
            matchable=false;
            break;
          }
          else {
            matches++;
          }
        }
      }
      if(matchable==true && matches==stoich1.size()){
        //cerr << "match found" << endl;
	      return true;
      }
      else {
	      return false;
      }
    }
    else {
      if(LDEBUG) {cerr << "compare:: " << "NUMBER OF TYPES OF ATOMIC SPECIES:   xstr1:  " << xstr1.num_each_type.size() << " " << xstr1.title << endl << xstr1
	                    << "           xstr2:  " << xstr2.num_each_type.size() << " " << xstr2.title << endl << xstr2 << endl;}
      if(LDEBUG) {cerr << "compare:: " << "NUMBER OF TYPES OF ATOMIC SPECIES IS NOT THE SAME...QUITTING..." << endl;}
      return false;
    }
  }
}

// ***************************************************************************
// Same Species 	
// ***************************************************************************
namespace compare{
  bool sameSpecies(const xstructure& xstr1, const xstructure& xstr2, const bool& display){

    // Determine if the structures have the same types and counts of species

    bool LDEBUG=(false || XHOST.DEBUG);

    // Check number of types
    if(xstr1.num_each_type.size() != xstr2.num_each_type.size()){
      // Display counts
      if(display==true){
        if(LDEBUG) {
          cerr << "compare::NUMBER OF TYPES OF ATOMIC SPECIES:   xstr1:  " << xstr1.num_each_type.size() 
               << " 		xstr2:  " << xstr2.num_each_type.size() << endl;
          cerr << "compare::NUMBER OF TYPES OF ATOMIC SPECIES IS NOT THE SAME..." << endl;
        }
        return false;
      }
    }

    // Check counts and species
    for(uint i=0;i<xstr1.num_each_type.size();i++){
      bool matched = false;
      for(uint j=0;j<xstr2.num_each_type.size();j++){
        // DX IS SPECIES CHECK TOO STRICT? if(xstr1.num_each_type[i] == xstr2.num_each_type[j] &&
        // DX IS SPECIES CHECK TOO STRICT?   xstr1.species[i] == xstr2.species[j]){
        if(xstr1.num_each_type[i] == xstr2.num_each_type[j]){
          matched = true;
          break;
        }
      }
      if(matched == false){
        if(display==true){ 
          if(LDEBUG) {
            cerr << "compare::WARNING:: TYPE OF ATOMIC SPECIES OR NUMBER PER TYPE ARE NOT THE SAME..." << endl;  
          }
        }
        return false;
      }
    }
    if(display==true){
      if(LDEBUG) {
	      cerr << "compare::NUMBER AND TYPE OF ATOMIC SPECIES ARE THE SAME...PROCEEDING..." << endl;
      }
    }
    return true;
  }
}

// ***************************************************************************
// Rescale Structure
// ***************************************************************************
namespace compare{
  void rescaleStructure(xstructure& xstr1, xstructure& xstr2){

    // If the scale factor is different, the two structures are rescaled to 1.00 

    if(abs(xstr1.scale-xstr2.scale)>0.001){
       xstr1.ReScale(1.0);
       xstr2.ReScale(1.0);
       xstr1.FixLattices();
       xstr2.FixLattices();
    }
  }
}

// ***************************************************************************
// Atomic Number Desnity
// ***************************************************************************
namespace compare{
  void atomicNumberDensity(xstructure& xstr1, xstructure& xstr2) { 

    // To compare structure with different volumes
    // we rescale the cell so that the volume divided by
    // the number of atoms is the same.
    // Update the Cartesian coordinates after scaling (DX 20181003)

    double scale;
    //cerr << xstr1.Volume()/xstr1.atoms.size() << " vs " << xstr2.Volume()/xstr2.atoms.size() << endl;
    scale=(xstr1.Volume()/xstr1.atoms.size())/(xstr2.Volume()/xstr2.atoms.size());
    xstr2.InflateVolume(scale);
    // update Cartesian coordinates
    for(uint i=0; i<xstr2.atoms.size(); i++){
      xstr2.atoms[i].cpos=F2C(xstr2.lattice,xstr2.atoms[i].fpos);
    }
  }
}

// ***************************************************************************
// Fake Atoms Name
// ***************************************************************************
namespace compare{
  vector<string> fakeElements(const uint& number_of_species){

    // Return vector of fake letters

    vector<string> elements;
   
    string letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    for(uint i=0;i<number_of_species;i++){
      stringstream ss_letter; ss_letter << letters[i]; // cannot type cast char to string directly
      elements.push_back(ss_letter.str());
    }
    
    return elements;
  }
}

// ***************************************************************************
// Fake Atoms Name
// ***************************************************************************
namespace compare{
  void fakeAtomsName(xstructure& xstr){

    // Assign a fake letter to each atom type. In case of materials with more 
    // than 26 species it is necessary to add more characters to this string

    string letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    int iat=0;

    for(uint i=0; i<xstr.num_each_type.size(); i++){
      xstr.species[i]=letters[i];
      for(int j=0; j<xstr.num_each_type[i]; j++){
        xstr.atoms[iat].name=letters[i];
	iat++;
      }
    }
  }
}

// ***************************************************************************
// Print Structure Parameters
// ***************************************************************************
namespace compare{
  void printParameters(xstructure& xstr, ostream& oss) { 

    // Print lattice parameters and volume for the xstructure

    oss << "========================================================" << endl;
    oss << xstr.title << endl;
    oss << "Parameters:" << endl;
    oss << "a:     " << xstr.a << endl;
    oss << "b:     " << xstr.b << endl;
    oss << "c:     " << xstr.c << endl;
    oss << "alpha: " << xstr.alpha << endl;
    oss << "beta:  " << xstr.beta << endl;
    oss << "gamma: " << xstr.gamma << endl;
    oss << "volume:" << xstr.GetVolume() << endl;
  }
}

// ***************************************************************************
// Least Frequent Atom
// ***************************************************************************
namespace compare{
  string leastFrequentAtom(const xstructure& xstr) {

    // Least frequent atom can be exploited in case we want to reduce the number 
    // of quadruplets to compute the rotational matrix by considering only the 
    // LFA species. It differs from other leastfrequentatom2 function because it 
    // finds the first LFA atom in xstructure 1. We do not need to find more the 
    // one LFA species in xstructure 1 because we will find all the LFAs for 
    // xstructure 2 and if the structures are a match we are guaranteed to find 
    // a one-to-one correspondence from one of the LFAs in xstructure 2.

    int flag=0, leastFrequentAtomCount=0;
    string leastFrequentAtomType;


    for(uint i=0; i<xstr.num_each_type.size(); i++){
      if(flag==0){
        // Number of occurrences of the LFA
        leastFrequentAtomCount = xstr.num_each_type[i];   
        // LFA Species 
        leastFrequentAtomType = KBIN::VASP_PseudoPotential_CleanName(xstr.species[i]); //DX 20190329 - remove pseudopotential info        
        flag=1;
      }
      else {
        if(leastFrequentAtomCount>xstr.num_each_type[i]){
          leastFrequentAtomCount = xstr.num_each_type[i];
          leastFrequentAtomType = KBIN::VASP_PseudoPotential_CleanName(xstr.species[i]); //DX 20190329 - remove pseudopotential info
        }   
      }   
    }   
    return leastFrequentAtomType;
  }
}

// ***************************************************************************
// Least Frequent Atom 2
// ***************************************************************************
namespace compare{
  vector<string> leastFrequentAtom2(const xstructure& xstr) {

    // This least frequent atom function finds all possible least frequent atoms 
    // for an xstructure and stores them in a vector. All of these LFAs are used 
    // in the quadruplet search.
    // We may not need to search over multiple LFAs during the quadruplet search. 
    // If a match is not found for one LFA, it won't be found for another since 
    // we need to map all atoms in one structure to the other structure. We will 
    // leave this implementation in for now, but may speed up the quadruplet 
    // search if we consider only one LFA.

    int flag=0, leastFrequentAtomCount=0;
    vector<string> leastFrequentAtomType;

    for(uint i=0; i<xstr.num_each_type.size(); i++){
      if(flag==0){
        // Number of occurrences of the LFA
        leastFrequentAtomCount = xstr.num_each_type[i];   
        // LFA Species 
        leastFrequentAtomType.push_back(KBIN::VASP_PseudoPotential_CleanName(xstr.species[i])); //DX 20190329 - remove pseudopotential info  
        flag=1;
      }
      else {
        if(leastFrequentAtomCount>xstr.num_each_type[i]){
          leastFrequentAtomCount = xstr.num_each_type[i];
          leastFrequentAtomType.clear();
          leastFrequentAtomType.push_back(KBIN::VASP_PseudoPotential_CleanName(xstr.species[i])); //DX 20190329 - remove pseudopotential info
        }
        if(leastFrequentAtomCount==xstr.num_each_type[i] && leastFrequentAtomType[0]!=KBIN::VASP_PseudoPotential_CleanName(xstr.species[i])){ //DX 20190329 - remove pseudopotential info 
          // Added the statement after '&&' (above); ensures no double counting from previous if statement
          leastFrequentAtomType.push_back(KBIN::VASP_PseudoPotential_CleanName(xstr.species[i])); //DX 20190329 - remove pseudopotential info
        }
      }
    }
    return leastFrequentAtomType;
  }
}

// ***************************************************************************
// Check Tolerances
// ***************************************************************************
namespace compare{
  bool checkTolerance(xvector<double> d1, xvector<double> d2){

    // Look for 2 corresponding reference frames, check that
    // the length of the 3 vectors and the angles between them are within
    // a given tolerance (in this case 30% of the reference value has been set).

    double tol_length=0.3, tol_angle=0.3;

    if( abs(d1(1)-d2(1)) < tol_length*abs(d1(1)) &&
        abs(d1(2)-d2(2)) < tol_length*abs(d1(2)) &&
        abs(d1(3)-d2(3)) < tol_length*abs(d1(3)) &&
        abs(d1(4)-d2(4)) < tol_angle*abs(d1(4)) &&
        abs(d1(5)-d2(5)) < tol_angle*abs(d1(5)) &&
        abs(d1(6)-d2(6)) < tol_angle*abs(d1(6))
      ){
      return false;
    }
    else {	
      return true;
    }
  }
}

// ***************************************************************************
// Check ABC Tolerances
// ***************************************************************************
namespace compare{
  bool checkABCTolerance(xvector<double> d1, xvector<double> d2){

    // Similar to checkTolerance, but it only looks at the length of the lattice 
    // vectors (for screening structures; makes faster)

    double tol_length=0.3;

    for(uint i=1;i<4;i++){
      for(uint j=1;j<4;j++){
	if(j!=i){
	  for(uint k=1;k<4;k++){
	    if(k!=i && k!=j){
	      if(abs(d1(i)-d2(1)) < tol_length*abs(d1(i)) &&
		 abs(d1(j)-d2(2)) < tol_length*abs(d1(j)) &&
		 abs(d1(k)-d2(3)) < tol_length*abs(d1(k))
		){
		return false;
	      }
	    }
	  }
	}
      }
    }
    return true;
  }
}

// ***************************************************************************
// Check Angle Tolerances
// ***************************************************************************
namespace compare{
  bool checkAngleTolerance(xvector<double> d1, xvector<double> d2){

    // Similar to checkTolerance, but it only looks at the angles betweeen the 
    // lattice vectors (for screening structures; makes faster)

    double tol_angle=0.3;

    for(uint i=4;i<7;i++){
      for(uint j=4;j<7;j++){
	if(j!=i){
	  for(uint k=4;k<7;k++){
	    if(k!=i && k!=j){	
	      if(abs(d1(i)-d2(4)) < tol_angle*abs(d1(4)) &&
		 abs(d1(j)-d2(5)) < tol_angle*abs(d1(5)) &&
		 abs(d1(k)-d2(6)) < tol_angle*abs(d1(6))
		){
		return false;
	      }
	    }
	  }
	}
      }
    }
    return true;
  }
}

// ***************************************************************************
// Find centroid for system with periodic boundary conditions
// ***************************************************************************
namespace compare{
  xvector<double> centroid_with_PBC(const xstructure& xstr){

    // Calculate the "best" centroid in a system with periodic boundary conditions.
    // This is based on the algorithm proposed in: 
    // https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
    // Used to find the best origin/centroid for a crystal structure.
    // Use of this for finding the best origin is still in beta testing; 
    // the method has some issues

    vector<xvector<double> > coordinates;
    for(uint i=0;i<xstr.atoms.size();i++){
      coordinates.push_back(xstr.atoms[i].cpos); //or cpos
    }
    return centroid_with_PBC(coordinates,xstr.lattice);
  }
}

// ***************************************************************************
// Find centroid for system with periodic boundary conditions
// ***************************************************************************
namespace compare{
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, const xmatrix<double>& lattice){

    // Calculate the "best" centroid in a system with periodic boundary conditions.
    // This is based on the algorithm proposed in: 
    // https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
    // If there are no weights (geometric center), then the weights are set to 1

    vector<double> weights;
    for(uint i=0;i<coordinates.size();i++){
      weights.push_back(1.0);
    }
    return centroid_with_PBC(coordinates,weights,lattice);
  }
}

// ***************************************************************************
// Find centroid for system with periodic boundary conditions
// ***************************************************************************
namespace compare{
  xvector<double> centroid_with_PBC(vector<xvector<double> >& coordinates, vector<double>& weights, 
				    const xmatrix<double>& lattice){

    // Calculate the "best" centroid in a system with periodic boundary conditions.
    // This is based on the algorithm proposed in: 
    // https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions

    xvector<double> centroid;
    for(uint i=1;i<4;i++){
      double zi_avg = 0.0;
      double zeta_avg = 0.0;
      double theta_avg =0.0;
      for(uint j=0;j<coordinates.size();j++){
        double theta = coordinates[j][i]*(2.0*Pi_r)/(aurostd::modulus(lattice(i)));
	double zi = std::cos(theta)*weights[j];
	double zeta = std::sin(theta)*weights[j];
	zi_avg += zi/coordinates.size();
	zeta_avg += zeta/coordinates.size();
      }
      theta_avg = std::atan2(-zeta_avg,-zi_avg);
      centroid(i) = theta_avg*(aurostd::modulus(lattice(i))/(2.0*Pi_r));
    }
    return centroid;
  } 
}


// ***************************************************************************
// Find Matches
// ***************************************************************************
namespace compare{
  bool findMatch(const xstructure& xstr1, const xstructure& PROTO, 
		    vector<uint>& im1, vector<uint>& im2, vector<double>& min_dists, 
                    const int& type_match) {

    // In order to find the best matchings the routine computes 
    // the difference between one atom and all the others, 
    // building a matrix of differences of coordinates.
    // Then, it checks which atoms have the best matching
    // with another atom in the second structure.

    // A | 1    A1, A2     I can check which is the best matching
    // B | 2 -> B1, B2  -> for the atom 1 and 2 in the structure
    // C |      C1, C2     with A,B,C,D 
    // D |      D1, D2    

    uint j=0,k=0;
    int i1=0,i2=0;                                  //Indices corresponding atoms

    vector<double> vdiffs;                      //Difference btwn atoms coords
    vector<std::pair<xvector<double>,xvector<double> > > min_positions;                      //Store sets of Cartesian coords which minimize distance
    vector<vector<double> > all_vdiffs;         //For all the atoms

    vector<uint> im1_tmp;
    vector<uint> im2_tmp;
    vector<string> im1_name;
    vector<string> im2_name;
    im1.clear();
    im2.clear();
    vdiffs.clear();
    all_vdiffs.clear();

    xmatrix<double> lattice=PROTO.lattice;

    double tmp = 1e9;
    uint i1_real=0;
    uint i2_real=0;
    string i1_name = "";
    string i2_name = "";

    //DX 20190226 [BETA] xvector<double> best_centroid1 = centroid_with_PBC(xstr1); 
    //DX 20190226 [BETA] xvector<double> best_centroid2 = centroid_with_PBC(PROTO); 
    
    vector<xvector<double> > l1, l2, l3;
    vector<int> a_index, b_index, c_index;
    for(int a=-2;a<=2;a++){l1.push_back(a*lattice(1));a_index.push_back(a);} //DX calc once and store
    for(int b=-2;b<=2;b++){l2.push_back(b*lattice(2));b_index.push_back(b);} //DX calc once and store
    for(int c=-2;c<=2;c++){l3.push_back(c*lattice(3));c_index.push_back(c);} //DX calc once and store

    for(j=0;j<xstr1.atoms.size();j++){
      std::pair<xvector<double>,xvector<double> > tmp_pair;
      xvector<double> tmp_xvec = xstr1.atoms[j].cpos;
      tmp_pair.first = tmp_xvec;
      vdiffs.clear();
      double match_dist=1e9;
      for(k=0;k<PROTO.atoms.size();k++){
        // speed increase, if same species map, don't try to map atoms of different types
        if(type_match==2 && xstr1.atoms[j].name!=PROTO.atoms[k].name){
          continue;
        }
        double dist=1e9;
        xvector<double> min_xvec;
        xvector<int> abc;
        xvector<double> incell_dist = xstr1.atoms[j].cpos-PROTO.atoms[k].cpos;
        // Need to find the min distance; thus check distance between neighboring cells to find true minimum.
        //DX - running vector in each loop saves computations; fewer duplicate operations
        for(uint m=0;m<l1.size();m++){
          xvector<double> a_component = incell_dist + l1[m];    // DX : coord1-coord2+a*lattice(1)
          for(uint n=0;n<l2.size();n++){
            xvector<double> ab_component = a_component + l2[n]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
            for(uint p=0;p<l3.size();p++){
              tmp_xvec = ab_component + l3[p];                       // DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
              tmp=aurostd::modulus(tmp_xvec);
        //DX 20190226 [OBSOLETE -SLOWER] for(int a=-2;a<=2;a++){
        //DX 20190226 [OBSOLETE -SLOWER]  for(int b=-2;b<=2;b++){
        //DX 20190226 [OBSOLETE -SLOWER]    for(int c=-2;c<=2;c++){
        //DX 20190226 [OBSOLETE -SLOWER]      tmp_xvec = PROTO.atoms[k].cpos+a*lattice(1)+b*lattice(2)+c*lattice(3);
        //DX 20190226 [OBSOLETE -SLOWER]      tmp = aurostd::modulus(xstr1.atoms[j].cpos-tmp_xvec); 
              if(tmp < dist){
                //DX 20190226 [OBSOLETE -SLOWER]abc(1)=a;
                //DX 20190226 [OBSOLETE -SLOWER]abc(2)=b;
                //DX 20190226 [OBSOLETE -SLOWER]abc(3)=c;
                abc(1)=a_index[m];
                abc(2)=b_index[n];
                abc(3)=c_index[p];
                i1 = j;
                i2 = k;
                //cerr << tmp << " vs " << dist << " - " << j << ", " << k << endl;
                dist = tmp;
                min_xvec = tmp_xvec;
              }
              // DXdist=aurostd::min(dist,aurostd::modulus(PROTO.atoms[j].cpos-xstr1.atoms[k].cpos+a*klattice(1)+b*klattice(2)+c*klattice(3)));
            }
          }
        }
        //cerr << "match_dist: " << match_dist << endl;
        if(dist<match_dist){
          i1_real=i1;
          i2_real=i2;
          match_dist = dist;
          i1_name = xstr1.atoms[i1_real].name;
          i2_name = PROTO.atoms[i2_real].name;
          tmp_pair.second = min_xvec;
        }
        vdiffs.push_back(dist);
      }
      // Check if same species match
      if(type_match == 2){ // same species
        if(i1_name != i2_name){
          //cerr << "WARNING: Matching species are not the same type, throwing out match (same species comparison)" << endl;
          return false;
        }
      }
      // Check for one-to-one mappings 
      for(uint i=0;i<im1_name.size();i++){
        // Check if i1 index has multiple mappings
        if(i1_real == im1[i]){
          //cerr << "WARNING: Used the same index for matching in i1! (" << i1_real << " == " << im1[i] << ")"<< endl;
          //cerr << match_dist << " vs " << min_dists[i] << endl;
          return false;
        }
        // Check if i2 index has multiple mappings
        if(i2_real == im2[i]){
          //cerr << "WARNING: Used the same index for matching in i2! (" << i2_real << " == " << im2[i] << ")"<< endl;
          //cerr << match_dist << " vs " << min_dists[i] << endl;
          return false;
        }
        // Check if types are not consistently mapped to a single type
        if(i1_name == im1_name[i]){
          if(i2_name != im2_name[i]){
            //cerr << "WARNING: Matching one type of atom to more than one type! (" << i1_name << " == " << i2_name << " | " << im1_name[i] << " == " << im2_name[i] << ")" <<  endl;
            //for(uint j=0;j<im1_name.size();j++){
            //  cerr << im1[j] << " == " << im2[j] << " | " << xstr1.atoms[im1[j]].cpos << " == " << PROTO.atoms[im2[j]].cpos << " (" << min_dists[i] << ") | " << im1_name[j] << " == " << im2_name[j] << endl;
            //}
            //cerr << i1_real << " == " << i2_real << " | " << xstr1.atoms[i1_real].cpos << " == " << PROTO.atoms[i2_real].cpos << " (" << match_dist << ") | " << i1_name << " == " << i2_name << endl;            
            return false;
          }
        }
      }
      im1_tmp.push_back(i1_real);
      im2_tmp.push_back(i2_real);
      im1.push_back(i1_real);
      im2.push_back(i2_real);
      im1_name.push_back(i1_name);
      im2_name.push_back(i2_name);
      min_dists.push_back(match_dist);
      all_vdiffs.push_back(vdiffs);
      min_positions.push_back(tmp_pair);
    }
    return true;
  }
}

// ***************************************************************************
// clusterize
// ***************************************************************************
namespace compare{
  void clusterize(const xstructure& xstr1, const vector<uint>& im1, vector<string>& TYPE1, 
		  vector<uint>& QTA1, vector<vector<uint> >& I1){

    // This function builds clusters/vectors of atoms of the same type. It is necessary when
    // we want to check correspondences between atoms of the same type or between atoms
    //of a specific species.

    uint i,j;
    int flag;
    vector<uint> i1;

    i1.clear();

    for(i=0; i<im1.size(); i++){
      i1.clear();
      if(i==0){
	TYPE1.push_back(xstr1.atoms[im1[i]].name);
	QTA1.push_back(1);
	i1.push_back(im1[i]);
	I1.push_back(i1);
      }
      else {
	flag=0;
	for(j=0; j<TYPE1.size(); j++){
	  if(xstr1.atoms[im1[i]].name==TYPE1[j]){
	    QTA1[j]++;
	    i1=I1[j];
	    i1.push_back(im1[i]);
	    I1[j]=i1;
	    flag=1;
	    break;
	  }
	}
	if(flag==0){
	  TYPE1.push_back(xstr1.atoms[im1[i]].name);
	  QTA1.push_back(1);
	  i1.push_back(im1[i]);
	  I1.push_back(i1);
	}
      }
    }
  }
}

// ***************************************************************************
// Same Atom Type
// ***************************************************************************
namespace compare{
  bool sameAtomType(const xstructure& xstr1, const xstructure& xstr2, const vector<uint>& im1, 
		    const vector<uint>& im2, const int& type_match){

    // Checks if mapping indices vectors are possible

    uint i,j,k,w,z;
    vector<string> TYPE1, TYPE2;
    vector<uint> QTA1,QTA2;
    vector<vector<uint> > I1,I2;

    vector<int> flag,checkType;

    clusterize(xstr1,im1,TYPE1,QTA1,I1);
    clusterize(xstr2,im2,TYPE2,QTA2,I2);

    checkType.clear();

    for(i=0; i<TYPE1.size(); i++){
      for(j=0; j<TYPE2.size(); j++){
	if(QTA1[i]==QTA2[j]){	
	  flag.clear();
	  for(w=0; w<I1[i].size(); w++){
	    for(z=0; z<I2[j].size(); z++){
	      for(k=0; k<im1.size(); k++){
		//Means we want the same atomic species to be matched up.
		if(type_match==2){           
		  if(I1[i][w]==im1[k] && I2[j][z]==im2[k] && 
		     xstr1.atoms[im1[k]].name == xstr2.atoms[im2[k]].name
		    ) 
		    flag.push_back(1);
		}
		else {
		  if(I1[i][w]==im1[k] && I2[j][z]==im2[k]) flag.push_back(1);
		}
	      }
	    }
	  }
	  if(flag.size()==QTA1[i])
	    checkType.push_back(1);	
	}
      }
    }
    if(checkType.size()==TYPE1.size()) return true;
    else return false;
  }
}

// ***************************************************************************
// Clean Match
// ***************************************************************************
namespace compare{
  bool cleanMatch(const vector<uint>& im1) {
    
    // The result of the findMatch function is a pair of vectors containing the indices of the matched
    // atoms of the structures. This function allows to delate the matchings where the same
    // index appears more the one time.

    uint i,j;

    for(i=0; i<im1.size(); i++){
      for(j=0; j<im1.size(); j++){
	if(i!=j){
	  if(im1[i]==im1[j]){
	    return true;
          }
        }
      }	
    }
    return false;
  }
}

// ***************************************************************************
// Cell Diagonal
// ***************************************************************************
namespace compare{
  void cellDiagonal(xstructure& xstr, vector<double>& diag_sum, 
		    vector<double>& diag_diff, const double& scale) {
    xmatrix<double> lattice = xstr.lattice;
    cellDiagonal(lattice, diag_sum, diag_diff, scale);
  }
}
namespace compare{
  void cellDiagonal(xmatrix<double>& lattice, vector<double>& diag_sum, 
		    vector<double>& diag_diff, const double& scale) {

    // The cell diagonals are represented by the shortest and longest diagonals
    // (in the case of cubic lattice the 2 diagonals are equal)
    // They are obtained as the vectorial sum and difference of the lattice
    // basis vectors.

    xvector<double> origin;

    diag_sum.push_back(abs(distance(origin,(lattice(1)+lattice(2))*scale)));
    diag_sum.push_back(abs(distance(origin,(lattice(1)+lattice(3))*scale)));
    diag_sum.push_back(abs(distance(origin,(lattice(2)+lattice(3))*scale)));
    diag_diff.push_back(abs(distance(lattice(1),lattice(2))*scale));
    diag_diff.push_back(abs(distance(lattice(1),lattice(3))*scale));
    diag_diff.push_back(abs(distance(lattice(2),lattice(3))*scale));
  }
}

// ***************************************************************************
// Lattice Deviation
// ***************************************************************************
namespace compare{
  double latticeDeviation(const vector<double>& diag_sum1,const vector<double>& diag_sum2, 
			  const vector<double>& diag_diff1,const vector<double>& diag_diff2) {

    // Lattice deviation is computed as the deviation of the 2 diagonals of each face
    // normalized on the sum of the diagonals of the faces of the reference structure.
    // The images of the diagonals of the mapped structure must be rescaled such that
    // the volume of the image of its unit cell of equals the volume of the unit
    // cell of the reference one.

    uint i;
    double d;
    vector<double> dev;

    for(i=0; i<3; i++) 
      dev.push_back((abs(diag_sum1[i]-diag_sum2[i])+abs(diag_diff1[i]-diag_diff2[i]))/(diag_sum2[i]+diag_diff2[i]));

    d=1;
    for(i=0;i<dev.size();i++) 
      d=d*(1-dev[i]);
    d=1-d;

    return d;
  }
}

// ***************************************************************************
// Compute nearest neighbors 
// ***************************************************************************
namespace compare{
  vector<double> computeNearestNeighbors(xstructure& xstr){

    // Determine the nearest neighbor distances centered on each atom
    // of the structure (needed for misfit calculation)

    vector<double> all_nn_distances;
    double nn = 1e9;
    for(uint i=0;i<xstr.atoms.size();i++){
      nn = shortestDistance(xstr,i);
      all_nn_distances.push_back(nn);
    }
    return all_nn_distances;
  }
}

// ***************************************************************************
// Shortest Distance from one atom
// ***************************************************************************
namespace compare{
  double shortestDistance(const xstructure& xstr, const uint& k) {

    // Find the minimum interatomic distance in the structure 

    //double min;
    double min_dist=1e9;
    xmatrix<double> lattice = xstr.lattice; //NEW

    //DX speed increase
    //perhaps can speed up even more, since the lattice doesn't change for the xstr...
    vector<xvector<double> > l1, l2, l3;
    vector<int> a_index, b_index, c_index;
    l1.clear(); l2.clear(); l3.clear();
    a_index.clear(); b_index.clear(); c_index.clear();
    for(int a=-2;a<=2;a++){l1.push_back(a*lattice(1));a_index.push_back(a);} //DX calc once and store, faster
    for(int b=-2;b<=2;b++){l2.push_back(b*lattice(2));b_index.push_back(b);} //DX calc once and store, faster
    for(int c=-2;c<=2;c++){l3.push_back(c*lattice(3));c_index.push_back(c);} //DX calc once and store, faster
    
    xvector<double> tmp;

    for(uint ii=0; ii<xstr.atoms.size(); ii++){
      if(ii!=k){
        xvector<double> incell_dist = xstr.atoms[k].cpos-xstr.atoms[ii].cpos;
        //DX 4/23/18 - running vector in each loop saves computations; fewer duplicate operations
        for(uint m=0;m<l1.size();m++){
          xvector<double> a_component = incell_dist + l1[m];    // DX : coord1-coord2+a*lattice(1)
          for(uint n=0;n<l2.size();n++){
            xvector<double> ab_component = a_component + l2[n]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
            for(uint p=0;p<l3.size();p++){
              tmp = ab_component + l3[p];                       // DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
              min_dist=aurostd::min(min_dist,aurostd::modulus(tmp));
            }
          }
        }
      }
    }
    //DX 20190226 [OBSOLETE - SLOWER] for(uint ii=0; ii<xstr.atoms.size(); ii++){
    //DX 20190226 [OBSOLETE - SLOWER]   if(ii!=k){
	  //DX 20190226 [OBSOLETE - SLOWER]     for(int a=-2;a<=2;a++){
	  //DX 20190226 [OBSOLETE - SLOWER]       for(int b=-2;b<=2;b++){
	  //DX 20190226 [OBSOLETE - SLOWER]         for(int c=-2;c<=2;c++){
    //DX 20190226 [OBSOLETE - SLOWER]           double tmp = aurostd::modulus(xstr.atoms[k].cpos-xstr.atoms[ii].cpos+a*lattice(1)+b*lattice(2)+c*lattice(3));
	  //DX 20190226 [OBSOLETE - SLOWER]           min_dist=aurostd::min(min_dist,tmp);
	  //DX 20190226 [OBSOLETE - SLOWER]         }
	  //DX 20190226 [OBSOLETE - SLOWER]       }
	  //DX 20190226 [OBSOLETE - SLOWER]     }      
    //DX 20190226 [OBSOLETE - SLOWER]   }
    //DX 20190226 [OBSOLETE - SLOWER] }

    return min_dist;
  }
}

// ***************************************************************************
// Coordinates Deviation
// ***************************************************************************
namespace compare{
  void coordinateDeviation(const xstructure& xstr1, const xstructure& xstr2, 
                      const vector<double>& all_nn1, const vector<double>& all_nn_proto,
		      const vector<uint>& indexMatch1, const vector<uint>& indexMatch2, vector<double>& min_dists,
		      double& cd, double& fail_figure) {

    // Compute the coordinates deviation by looking at each pair of atoms from
    // the reference and mapped structure

    uint j;
    double num=0, den=0, nfail=0;
    double dd, nn1, nn2; //dd=delta distance, nn=nearest neighbour 
    int  fail1, fail2;
    xmatrix<double> klattice = xstr1.lattice;
    for(j=0; j<indexMatch1.size(); j++){
      //nn1=shortestDistance(xstr1,indexMatch1[j]);
      //nn2=shortestDistance(xstr2, indexMatch2[j]);
      nn1 = all_nn1[indexMatch1[j]];
      nn2 = all_nn_proto[indexMatch2[j]];
      dd = min_dists[j];
//      dd=1e9;
//      // Need to find the min distance; thus check distance between neighboring cells to find true minimum.
//      for(int a=-2;a<=2;a++){
//	for(int b=-2;b<=2;b++){
//	  for(int c=-2;c<=2;c++){
//	    dd=aurostd::min(dd,modulus(xstr1.atoms.at(indexMatch1[j]).cpos-xstr2.atoms.at(indexMatch2[j]).cpos+a*klattice(1)+b*klattice(2)+c*klattice(3)));
//	  }
//	}
//      }

      //DX [OBSOLETE] if(dd<=0.5*nn1) fail1=0;
      //DX [OBSOLETE] if(dd>0.5*nn1) fail1=1;
      //DX [OBSOLETE] if(dd<=0.5*nn2) fail2=0;
      //DX [OBSOLETE] if(dd>0.5*nn2) fail2=1;
      //DX 20190226 - below is a bit faster than above
      if(dd<=0.5*nn1){ fail1=0; }
      else { fail1=1; }
      if(dd<=0.5*nn2) fail2=0;
      else { fail2=1; }

      if(fail1==0){
	num=num+dd;
	den=den+nn1;
      }   
      if(fail2==0){
	num=num+dd;
	den=den+nn2;
      }   
      if(fail1==1) nfail++;
      if(fail2==1) nfail++;
    }   

    if(den==0) cd=1;
    else cd=num/den;
    
    //Consider unmatched atoms
    int flag=0;
    for(uint i=0; i<xstr1.atoms.size();i++){
      flag=0;
      for(uint k=0; k<indexMatch1.size();k++){
	if(i==indexMatch1[k]){
	  flag=1;
	}
      }
      if(flag==0){	//Meaning this atom does not have a match; increase the figure of failure
	nfail++;
      }
    }
    for(uint i=0; i<xstr2.atoms.size();i++){
      flag=0;
      for(uint k=0; k<indexMatch2.size();k++){
	if(i==indexMatch2[k]){
	  flag=1;
	}
      }
      if(flag==0){        //Meaning this atom does not have a match; increase the figure of failure
	nfail++;
      }
    }
    
    fail_figure=(nfail/(xstr1.atoms.size()+xstr2.atoms.size()));
  }   
}

// ***************************************************************************
// Compute Misfit
// ***************************************************************************
namespace compare{
  double computeMisfit(const double& dev, const double& dis, const double& fail) {

    // Combines differences between all aspects of crystal structure into a 
    // figure of misfit. (See Burzlaff)

    double mis;

    mis=1-((1-dev)*(1-dis)*(1-fail));
    return mis;
  }   
}

// ***************************************************************************
// Print Matching Between Atoms
// ***************************************************************************
namespace compare{
  void printMatch(const vector<uint>& indexMatch1, const vector<uint>& indexMatch2, 
		  const xstructure& PROTO, const xstructure& xstr1, ostream& oss) {

    // With this function we print the atoms matched in the previous function
    // whose indices are stored in the indexMatch vector.
    // This allows us to call them directly.
    // Lastly, the translational term, specific of each mapping, is printed.

    uint i,j;

    oss << "              Reference                               Mapped"<<endl;
    for(j=0; j< indexMatch1.size(); j++){
      oss << indexMatch1[j]<<"-"<<indexMatch2[j]<<"    ";
      oss << xstr1.atoms[indexMatch1[j]].cpos << " " << xstr1.atoms[indexMatch1[j]].name;
      oss << "       ";
      oss << PROTO.atoms[indexMatch2[j]].cpos << " " << PROTO.atoms[indexMatch2[j]].name << endl;
    }

    int flag=0;

    oss <<"----------------------------------------------------"<<endl;
    oss << "Missing Atoms in Reference structure:"<< endl;
    for(j=0; j<xstr1.atoms.size(); j++){
      flag=0;
      for(i=0; i<indexMatch1.size(); i++){
	if(j==indexMatch1[i])
	  flag=1;
      }
      if(flag==0){
	oss << "# "<< j << "   " << xstr1.atoms[j].cpos << "   " << xstr1.atoms[j].name << endl;
      }
    }

    oss << "Missing Atoms in Mapped structure:" << endl;
    for(j=0; j<PROTO.atoms.size(); j++){
      flag=0;
      for(i=0; i<indexMatch2.size(); i++){
	if(indexMatch2[i]==j) flag=1;
      }
      if(flag==0) oss << "# "<< j << "   " << PROTO.atoms[j].cpos << "   " << PROTO.atoms[j].name << endl;
    }
  }
}

// ***************************************************************************
// Bring Coordinate in the cell (Similar to xstructure.BringInCell())
// ***************************************************************************
namespace compare{
  xvector<double> bringCoordinateInCell(xvector<double>& coord){

    // This function brings a coordinate back in the unit cell.  
    // It is not an xstructure attribute, so we cannot use xstructure.BringInCell().

    double tol=1e-6;
    for(uint i=1;i<4;i++){
      if(coord(i)<-tol){
	coord(i)+=1.0;
      }
      else if(coord(i)>=1.0-tol){
	coord(i)-=1.0;
      }
    }
    return coord;
  }
}

// ***************************************************************************
// Determine if vector is Periodic
// ***************************************************************************
namespace compare{
  bool vectorPeriodic(const xvector<double>& vec, const xstructure& lfa_supercell, 
		       const int& i, const int& j){

    // Once we have a possible quadruplet (lattice), we need to make sure that this 
    // choice of the primitive cell preserves the periodicity o the lattice. 
    // Therefore, we check that each of the quadruplet atoms maps onto another atom 
    // in the supercell. Helpful analogy: Lattice periodicty vs crystal periodicity. 
    // The quadruplets form the lattice and in this function we check for lattice
    // periodicity. The misfit criteria checks the crystal periodicity.

    double tolerance = 0.01; // Hundredth of an angstrom
    deque<_atom> atoms = lfa_supercell.atoms;
    xmatrix<double> lattice = lfa_supercell.lattice;
    xmatrix<double> f2c = trasp(lattice);
    //DX 20190619 [OBSOLETE] xmatrix<double> c2f = inverse(trasp(lattice));
    bool skew = false;

    vector<int> ind(2); ind[0]=i, ind[1]=j;
    xvector<double> tmp;

    uint count=0;
    deque<_atom> transformed;
    deque<uint> index_to_check;

    // ===== Check if applying the symmetry element along with internal translation maps to another atom ===== //
    for(uint d=0;d<atoms.size();d++){
      _atom tmp;
      tmp.type = atoms[d].type;
      tmp.cpos = atoms[d].cpos+vec;
      tmp.fpos = C2F(lattice,tmp.cpos);
      if(SYM::MapAtom(atoms,tmp,true,lattice,f2c,skew,tolerance)){ //DX 20190619 - lattice and f2c as input
        transformed.push_back(tmp);
        index_to_check.push_back(d);
        count++;
      }
      else {
        return false;
      }
    }
    if(count == atoms.size()){
      return true;
    }
    return false;
  }
}

// ***************************************************************************
// Determine if LFA Quadruplet is Periodic
// ***************************************************************************
namespace compare{
  bool quadrupletPeriodic(const xmatrix<double>& quad, const xstructure& lfa_supercell, 
		       const int& i, const int& j, const int& k, const int& w){

    // Once we have a possible quadruplet (lattice), we need to make sure that this 
    // choice of the primitive cell preserves the periodicity o the lattice. 
    // Therefore, we check that each of the quadruplet atoms maps onto another atom 
    // in the supercell. Helpful analogy: Lattice periodicty vs crystal periodicity. 
    // The quadruplets form the lattice and in this function we check for lattice
    // periodicity. The misfit criteria checks the crystal periodicity.

    vector<int> ind(4); ind[0]=i, ind[1]=j; ind[2]=k; ind[3]=w;
    vector<xvector<double> > latt_vecs(3); 
    latt_vecs[0]=quad(1); latt_vecs[1]=quad(2); latt_vecs[2]=quad(3);
    xvector<double> tmp;

    for(uint b=0; b<ind.size();b++){
      for(uint c=0;c<latt_vecs.size();c++){
	tmp=lfa_supercell.atoms.at(ind[b]).cpos+latt_vecs[c];
	bool match_found=false;
	for(uint d=0;d<lfa_supercell.atoms.size();d++){
	  if(abs(lfa_supercell.atoms.at(d).cpos(1)-tmp(1))<0.01 && 
	     abs(lfa_supercell.atoms.at(d).cpos(2)-tmp(2))<0.01 && 
	     abs(lfa_supercell.atoms.at(d).cpos(3)-tmp(3))<0.01){ // Less than hundredth of Angstrom?
	    match_found=true;
	    break;
	  }
	}
	if(match_found==false){
	  xvector<double> tmp_frac=C2F(lfa_supercell.lattice,tmp);
	  tmp_frac=bringCoordinateInCell(tmp_frac);
	  xstructure lfa_supercell_tmp = lfa_supercell;
	  for(uint f=0;f<lfa_supercell_tmp.atoms.size();f++){
	    lfa_supercell_tmp.atoms.at(f).fpos=C2F(lfa_supercell_tmp.lattice,lfa_supercell_tmp.atoms.at(f).cpos);
	    if(abs(lfa_supercell_tmp.atoms.at(f).fpos(1)-tmp_frac(1))<0.01 && 
	       abs(lfa_supercell_tmp.atoms.at(f).fpos(2)-tmp_frac(2))<0.01 && 
	       abs(lfa_supercell_tmp.atoms.at(f).fpos(3)-tmp_frac(3))<0.01){
	      match_found=true;
	      break;
	    }
	  }
	  if(match_found==false){
	    return false;
	  }
	}
      }
    }
    return true;
  }
}

// ***************************************************************************
// Thread Generation (For parallel processing of quadruplets)
// ***************************************************************************
namespace compare{
  void threadGeneration(const uint& num_proc,xmatrix<double>& q1, xstructure& xstr2, 
			vector<xstructure> &vprotos, xstructure &xstr1, const int& type_match, 
			const bool& optimize_match, double& minMis, ostream& oss){ 

    // This function creates the supercell of the second structure and begins 
    // the quadruplets search within a supercell. Due to the costly nature of 
    // this algorithm, the quadruplet search is parallelized. The splitting 
    // of computation of done here

    bool LDEBUG=(false || XHOST.DEBUG);
    bool test_one_lfa_only = false; //DX 20190318
    if(type_match==2){ test_one_lfa_only=true;} //DX 20190318


    xstructure xstr1_tmp = xstr1;

    //cerr << "LFA" << endl;
    vector<string> LFA_str1=leastFrequentAtom2(xstr1);
    vector<string> LFA_str2=leastFrequentAtom2(xstr2);
    string lfa, lfa_str1;

    //cerr << "SUPERCELL" << endl;
    //DX 20190319 - START
    vector<int> sc2pcMap, pc2scMap;
    bool get_symmetry=false;
    bool get_full_basis=false;
    bool force_supercell_matrix=true;
    xstructure xstr=GetSuperCell(xstr2,3,0,0,0,3,0,0,0,3,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix); //DX 20190319 - use supercell matrix in expansion

    uint y=0;
    uint x=0;

    // DX TEST
    // Consider all LFAs
    for(y=0;y<LFA_str2.size();y++){
    for(x=0;x<LFA_str1.size();x++){
    // DX TEST
      lfa_str1=LFA_str1[x];
      lfa=LFA_str2[y];
      if(type_match == 2 && lfa_str1 != lfa){
        continue;
      }
      oss << "===> LFA: "<<lfa<<endl;

      if(LDEBUG) {
        cerr << "===> LFA_1: " << lfa_str1 <<endl;
        cerr << "===> LFA: "<<lfa<<endl;
      }
      
      // DX NEW so we don't need to do this in an inner loop
      for(uint a=0;a<xstr.atoms.size();a++){
        if(xstr.atoms[a].name == lfa){
          xstr.ShifOriginToAtom(a);
          break;
        }
      }
      // DX NEW

      //cerr << "xstr1 centroid: " << endl;

      //cerr << "SHIFT" << endl;
      // NEED TO SHIFT origin of xstr1_tmp to one of the LFA (this was missing before and caused ICSD_102428.BCA, and CBA to not match, but they should
      for(uint i=0;i<xstr1_tmp.atoms.size();i++){
        if(xstr1_tmp.atoms[i].name==lfa_str1){
          xstr1_tmp.ShifOriginToAtom(i);
          xstr1_tmp.BringInCell(1e-10);
          break;
        }
      }

      // NEED TO SHIFT origin of xstr2 to one of the LFA
      for(uint i=0;i<xstr2.atoms.size();i++){
        if(xstr2.atoms[i].name==lfa){
          xstr2.ShifOriginToAtom(i);
          xstr2.BringInCell(1e-10);
          break;
        }
      }

      // When checking the quadruplets/lattice, we only need to generate a 
      // LFA supercell (i.e. take out the other atoms).  This greatly reduces 
      // the time of computation (don't need to scan through unnecessary atoms) 
      xstructure xstr_LFA_only=xstr2;
      xstr_LFA_only.ClearSymmetry(); //DX 20181022
      uint num_atoms=xstr_LFA_only.atoms.size();
      for(uint i=0;i<num_atoms;i++){
	if(xstr_LFA_only.atoms[i].name!=lfa){
	  xstr_LFA_only.RemoveAtom(i);
	  num_atoms--;
	  i--;
	}
      }
      //cerr << "LFA SUPERELL" << endl;
      xstructure xstr_LFA_supercell=GetSuperCell(xstr_LFA_only,3,0,0,0,3,0,0,0,3,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix); //DX 20190319 - use supercell matrix in expansion
      //cerr << "created LFA supercell" << endl;
      // Determines the number of LFAs in the supercell.
      int num_LFAs=-1; //-1 as default value 
      for(uint q=0; q<xstr.num_each_type.size();q++){
	if(xstr.species[q] == lfa){ 
	  num_LFAs= xstr.num_each_type[q];
	  if(LDEBUG) {cerr << "compare:: " << "Number of LFAs in supercell: " << xstr.species[q] << "= " << num_LFAs << endl;}
	}
      }

      // === THREAD PREPARATION FOR PARALLEL PROCESSING OF QUADRUPLET SEARCH === //

      // DECLARATION OF ATOMIC BOOL SECTION: Allows the threads to communicate. 
      // This is useful for stopping the threads if the misfit falls below
      // the compatible misfit criterion (mis<0.1) in any of the threads. 
      // [OBSOLETE] std::atomic_bool misfit_in_threshold_found (false);

      //[NONTHREADS]bool misfit_in_threshold_found=false;
   
      vector<xstructure> xstr1_for_thread;
      vector<double> possible_minMis;
      vector<vector<xstructure> > vvprotos;
      //vector<std::thread> threads;
      //DX NEW - used to be done in one of the inner loops in structureSearch 

      // compute xstr1 information once only (perhaps we can use this in the directory scheme!!!!!!!! 
      // and really only calculate once; but not that expensive, may be more expensive to store --memory!)
      vector<double> D1,F1;
      cellDiagonal(xstr1_tmp,D1,F1,1);
      xstr1_tmp.lattice=GetClat(xstr1_tmp.a,xstr1_tmp.b,xstr1_tmp.c,xstr1_tmp.alpha,xstr1_tmp.beta,xstr1_tmp.gamma);
      for(uint iat=0; iat<xstr1_tmp.atoms.size(); iat++){
        xstr1_tmp.atoms[iat].cpos=F2C(xstr1_tmp.lattice,xstr1_tmp.atoms[iat].fpos);
      }
      vector<double> all_nn1 = computeNearestNeighbors(xstr1_tmp);
      // DX NEW
      //for(uint n=0; n<num_proc; n++){
	    //  vector<xstructure> vprotos_tmp;
	    //  vvprotos.push_back(vprotos_tmp);
	    //  xstr1_for_thread.push_back(xstr1_tmp);
	    //  possible_minMis.push_back(1.0);
      //}
      
      vector<xmatrix<double> > lattices;
      vector<xmatrix<double> > clattices;
      vector<vector<uint> > ij_index;
      vector<double> latt_devs;

      //vector<std::thread> threads0;
      vector<xvector<double> > translation_vectors;
	    quadrupletSearch(q1,xstr_LFA_supercell,translation_vectors,ij_index);
      //cerr << "FINDING TRANSLATION VECTORS: " << endl;
      double abs_det_q1=abs(det(q1));
      xvector<double> abc_angles_q1=Getabc_angles(q1,DEGREES);

      buildSimilarLattices(translation_vectors, q1, abs_det_q1, abs_det_q1, abc_angles_q1, lattices, clattices, latt_devs, optimize_match);
      if(LDEBUG) {cerr << "pflow::threadGeneration: Number of lattices to compare: " << lattices.size() << endl;}

      if(lattices.size()>0){
  for(uint n=0; n<num_proc; n++){
	  vector<xstructure> vprotos_tmp;
	  vvprotos.push_back(vprotos_tmp);
	  xstr1_for_thread.push_back(xstr1_tmp);
	  possible_minMis.push_back(1.0);
  }
#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
	vector<std::thread> threads1;
	vector<vector<xmatrix<double> > > lattices_split;
	vector<vector<xmatrix<double> > > clattices_split;
	vector<vector<double> > latt_devs_split;
	uint num_per_thread = lattices.size()/num_proc;
	uint residual = lattices.size()%num_proc;
	bool accounted_for_residual=false;
	if(residual!=0){
	  num_per_thread+=1;
	}
	uint count = 0;
	uint thread_count = 0;
	vector<xmatrix<double> > latt_tmp, clatt_tmp;
	vector<double> tmp_dev;
	for(uint l=0; l<lattices.size(); l++){
	  latt_tmp.push_back(lattices[l]); clatt_tmp.push_back(clattices[l]); tmp_dev.push_back(latt_devs[l]);
	  count+=1;
	  if(count == num_per_thread && thread_count<num_proc-1){
	    thread_count+=1;
	    lattices_split.push_back(latt_tmp);
	    clattices_split.push_back(clatt_tmp);
	    latt_devs_split.push_back(tmp_dev);
	    latt_tmp.clear(); clatt_tmp.clear(); tmp_dev.clear();
	    count = 0;
	  }
	  else if(thread_count==num_proc-1 && l==lattices.size()-1){
	    thread_count+=1;
	    lattices_split.push_back(latt_tmp);
	    clattices_split.push_back(clatt_tmp);
	    latt_devs_split.push_back(tmp_dev);
	    latt_tmp.clear(); clatt_tmp.clear(); tmp_dev.clear();
	    count = 0;
	  }
	  if(!accounted_for_residual && residual!=0 && thread_count==residual){
	    accounted_for_residual=true;
	    num_per_thread=num_per_thread-1;
	  }
	}
	uint recovered=0;
	//Need the following safety in case the number of threads is greater than the number of lattices to test
	uint num_of_threads=0;
	if(lattices_split.size()>=num_proc){
	  num_of_threads=num_proc;
	}
	else if(lattices_split.size()<num_proc){
	  num_of_threads=lattices_split.size();
	}
	for(uint n=0; n<num_of_threads; n++){
	  for(uint h=0;h<lattices_split[n].size();h++){
	    recovered+=1;
	    //cerr << "recovered: " << recovered << " - " << lattices_split[n][h] << endl;
	  }
	} 
	if(recovered != lattices.size()){
	  cerr << "The splitting of jobs failed...not all were accounted for: " << recovered << " != " << lattices.size() << endl;
	  exit(1);
	}
  if(LDEBUG) {cerr << "pflow::threadGeneration: Performing structure search on " << lattices.size() << " lattices ..." << endl;}
	for(uint n=0; n<num_of_threads; n++){
	  threads1.push_back(std::thread(structureSearch,lfa,all_nn1,xstr,
			    std::ref(vvprotos[n]),std::ref(xstr1_for_thread[n]),xstr2,type_match,std::ref(possible_minMis[n]),
			    std::ref(lattices_split[n]),std::ref(clattices_split[n]),std::ref(latt_devs_split[n]),
                            optimize_match));
	}         
	for(uint t=0;t<threads1.size();t++){
	  threads1[t].join();
	}
#else
  uint n=0;
	structureSearch(lfa,all_nn1,xstr,vvprotos[n],xstr1_for_thread[n],xstr2,type_match,possible_minMis[n],
	   	        lattices,clattices,latt_devs,optimize_match);
#endif
      }
      //cerr << "========== possible_minMis.size(): " << possible_minMis.size() << endl;
      for(uint p=0;p<possible_minMis.size();p++){
	if(p==0 && y==0 && x==0){ // DX 2/8/17 - need to add x==0 ortherwise matches can be overwritten
	  minMis=possible_minMis[p];
	  xstr1=xstr1_for_thread[p];
	  vprotos=vvprotos[p];
	}
	else {
	  if(possible_minMis[p]<=minMis){
	    minMis=possible_minMis[p];
	    xstr1=xstr1_for_thread[p];
	    vprotos=vvprotos[p];
	  }
	}
        //cerr << "minMis: " << minMis << endl;
      }
      //break if(minMis<=0.1) break;

    // //DX 20190226 - fast return, no need to check other LFAs if a match is found - START
    if(minMis<0.1 && !optimize_match){
      return;
    }
    // //DX 20190226 - fast return, no need to check other LFAs if a match is found - END
    // DX 20190318 - START
    if(test_one_lfa_only){
        return;
    }
    // DX 20190318 - START
    } 
    } 
  }  
}

// ***************************************************************************
// Quadruplet Search
// ***************************************************************************
namespace compare{
  void quadrupletSearch(const xmatrix<double>& q1, const xstructure& xstr_LFA_supercell, 
			vector<xvector<double> >& lattice_vecs, vector<vector<uint> >& ij_index){

    // This function scans through the possible quadruplets (sets of 4 LFA atoms) i
    // to find a lattice which is commensurate with the reference structure (xstr1). 
    // This function is parallelized since it is the time-limiting function.

    bool LDEBUG=(false || XHOST.DEBUG);

    double min_q1_a = aurostd::modulus(q1(1))-aurostd::modulus(q1(1))*0.1;
    double max_q1_a = aurostd::modulus(q1(1))+aurostd::modulus(q1(1))*0.1;
    double min_q1_b = aurostd::modulus(q1(2))-aurostd::modulus(q1(2))*0.1;
    double max_q1_b = aurostd::modulus(q1(2))+aurostd::modulus(q1(2))*0.1;
    double min_q1_c = aurostd::modulus(q1(3))-aurostd::modulus(q1(3))*0.1;
    double max_q1_c = aurostd::modulus(q1(3))+aurostd::modulus(q1(3))*0.1;

    if(LDEBUG) {
      cerr << "compare::quadrupletSearch: Lattice parameters: " << aurostd::modulus(q1(1)) << ", " << aurostd::modulus(q1(2)) << ", " << aurostd::modulus(q1(3)) << endl;
      cerr << "compare::quadrupletSearch: Modulus search range for lattice vector a: " << min_q1_a << " - " << max_q1_a << endl;
      cerr << "compare::quadrupletSearch: Modulus search range for lattice vector b: " << min_q1_b << " - " << max_q1_b << endl;
      cerr << "compare::quadrupletSearch: Modulus search range for lattice vector c: " << min_q1_c << " - " << max_q1_c << endl;
    }

    xvector<double> tmp_vec;
    double tmp_mod = 0.0;
    
    // Search all possible vectors with modulus comparable to one of the lattice vectors 
    for(uint i=0; i<xstr_LFA_supercell.atoms.size(); i++){
      for(uint j=i+1; j<xstr_LFA_supercell.atoms.size(); j++){ //upper triangular
	tmp_vec = xstr_LFA_supercell.atoms[j].cpos-xstr_LFA_supercell.atoms[i].cpos;
        tmp_mod = aurostd::modulus(tmp_vec);
        if((tmp_mod <= max_q1_a && tmp_mod >= min_q1_a) || 
           (tmp_mod <= max_q1_b && tmp_mod >= min_q1_b) || 
           (tmp_mod <= max_q1_c && tmp_mod >= min_q1_c)){ 
          bool vec_stored = false;
          for(uint p=0;p<lattice_vecs.size();p++){
            if(identical(lattice_vecs[p],tmp_vec,1e-3)){ //DX 20190318 - changed from -10 to -3
              vec_stored = true;
              break;
            }
          }       
          if(vec_stored == false){
            lattice_vecs.push_back(tmp_vec);
            // Store indices of atoms comprising the vector
            vector<uint> ij;
            ij.push_back(i); ij.push_back(j);
            ij_index.push_back(ij); 
            // Store negative (may not be needed)
            //lattice_vecs.push_back(-tmp_vec);
            //vector<uint> ji;
            //ji.push_back(j); ji.push_back(i);
            //ij_index.push_back(ji); 
          }
        }
      }
    }

    if(LDEBUG) {
      cerr << "compare::quadrupletSearch: Number of potential lattice vectors: " << lattice_vecs.size() << endl;
    } 
    // Removing non-periodic lattice vectors
    vector<xvector<double> > lattice_vecs_periodic;
    for(uint i=0;i<lattice_vecs.size();i++){
      if(vectorPeriodic(lattice_vecs[i],xstr_LFA_supercell,ij_index[i][0],ij_index[i][1])){
        lattice_vecs_periodic.push_back(lattice_vecs[i]);
        //DX 20190318 [OBSOLETE] lattice_vecs_periodic.push_back(-lattice_vecs[i]);
      }
    }
    vector<xvector<double> > final_lattice_vecs = lattice_vecs_periodic; //DX 20190320
    //DX 20190318 - only store negative if not a duplicate - START
    for(uint i=0;i<lattice_vecs_periodic.size();i++){ //DX 20190320
      xvector<double> tmp = -lattice_vecs_periodic[i];
      bool vec_stored = false;
      for(uint p=0;p<final_lattice_vecs.size();p++){
        if(identical(final_lattice_vecs[p],tmp,1e-3)){
          vec_stored = true;
          break;
        }
      }       
      if(!vec_stored){
        final_lattice_vecs.push_back(tmp);
      }
    }
    //DX 20190318 - only store negative if not a duplicate - END
    lattice_vecs = final_lattice_vecs; //DX 20190320
    if(LDEBUG) {
      cerr << "compare::quadrupletSearch: Number of lattice vectors (preserves periodicity): " << lattice_vecs.size() << endl;
      for(uint i=0;i<lattice_vecs.size();i++){
        cerr << "compare::quadrupletSearch: lattice vector " << i << ": " << lattice_vecs[i] << " (" << aurostd::modulus(lattice_vecs[i]) << ")" << endl; 
      }
    } 
  }
}

// ***************************************************************************
// Build All Lattices
// ***************************************************************************
namespace compare{
  bool buildSimilarLattices(vector<xvector<double> >& translation_vectors, xmatrix<double>& q1, double& xstr1_vol, double& abs_det_q1, 
                            xvector<double>& abc_angles_q1, vector<xmatrix<double> >& lattices, vector<xmatrix<double> >& clattices, 
                            vector<double>& latt_devs, const bool& optimize_match){

    bool LDEBUG=(false || XHOST.DEBUG);

    vector<double> D1,F1;
    cellDiagonal(q1,D1,F1,1);

    double tol_vol=0.1;
    double det_tol=tol_vol*abs_det_q1;

    double tol_a=abc_angles_q1[1]*0.3;
    double tol_b=abc_angles_q1[2]*0.3;
    double tol_c=abc_angles_q1[3]*0.3;

    if(LDEBUG) {
      cerr << "compare::buildSimilarLattices: Tolerance for a (Angstroms): " << tol_a << endl;
      cerr << "compare::buildSimilarLattices: Tolerance for b (Angstroms): " << tol_b << endl;
      cerr << "compare::buildSimilarLattices: Tolerance for c (Angstroms): " << tol_c << endl;
      cerr << "compare::buildSimilarLattices: Tolerance for alpha (degrees): " << abc_angles_q1[4]*0.3 << endl;
      cerr << "compare::buildSimilarLattices: Tolerance for beta (degrees): " << abc_angles_q1[5]*0.3 << endl;
      cerr << "compare::buildSimilarLattices: Tolerance for gamma (degrees): " << abc_angles_q1[6]*0.3 << endl;
    }

    xmatrix<double> tmp_lattice(3,3);
    xmatrix<double> tmp_clatt(3,3);

    int store=0;
    vector<double> translations_mod;
    for(uint i=0;i<translation_vectors.size();i++){
      translations_mod.push_back(aurostd::modulus(translation_vectors[i]));
    }
    if(LDEBUG) {
      cerr << "buildSimilarLattices:: Number of lattice vectors: " << translation_vectors.size() << endl;
    }

    // Build all possible unit cells with combinations of lattice vectors (order matters, hence not upper triangular)
    for(uint i=0;i<translation_vectors.size();i++){
      if(abs(translations_mod[i]-abc_angles_q1[1])<tol_a){ //check a
	for(uint j=0;j<translation_vectors.size();j++){
	  if(j!=i){
	    if(abs(translations_mod[j]-abc_angles_q1[2])<tol_b){ // check b
	      for(uint k=0;k<translation_vectors.size();k++){
		if(k!=i && k!=j){
		  if(abs(translations_mod[k]-abc_angles_q1[3])<tol_c){ //check c
		    tmp_lattice = SYM::xvec2xmat(translation_vectors[i],translation_vectors[j],translation_vectors[k]);         
		    if(abs(abs_det_q1-abs(det(tmp_lattice))) < det_tol){ //check determinant/volume
		      xvector<double> abc_angles_q2=Getabc_angles(tmp_lattice,DEGREES);
		      if(checkTolerance(abc_angles_q1,abc_angles_q2)==false){
// DX 20190318 only calc if lattice dev is small - START
//			    tmp_clatt=GetClat(abc_angles_q2[1],abc_angles_q2[2],abc_angles_q2[3],abc_angles_q2[4],abc_angles_q2[5],abc_angles_q2[6]);
//			    bool unique = true;
//			    for(uint t=0;t<lattices.size();t++){
//			      if(identical(tmp_lattice,lattices[t],1e-10)){
//			    // DX TEST - CANNOT DO: Eliminates potential matches for(uint t=0;t<clattices.size();t++){
//			    // DX TEST - CANNOT DO: Eliminates potential matches if(identical(tmp_clatt,clattices[t],1e-10)){
//			    unique=false;
//			    break;
//			  }
//			}
// DX 20190318 only calc if lattice dev is small - END
			double tmp_latt_dev = checkLatticeDeviation(xstr1_vol,tmp_lattice,D1,F1);
                        // Time-saver, keep lattices that have a deviation smaller than Burzlaff's same-family requirement
                        // otherwise, there is no possible way that it could match with anything or be in the same-family
			// DX 20190318 [OBSOLETE] if(unique && !optimize_match && tmp_latt_dev <= 0.1){ //fast match doesn't care about finding same family information
			if(!optimize_match && tmp_latt_dev <= 0.1){ //fast match doesn't care about finding same family information //DX 20190318 - removed unique since it doesn't exist yet
			bool unique = true;
			for(uint t=0;t<lattices.size();t++){
			  if(identical(tmp_lattice,lattices[t],1e-10)){
			    // DX TEST - CANNOT DO: Eliminates potential matches for(uint t=0;t<clattices.size();t++){
			    // DX TEST - CANNOT DO: Eliminates potential matches if(identical(tmp_clatt,clattices[t],1e-10)){
			    unique=false;
			    break;
			  }
			}
                        //DX 20190318 - now store calc and store clattice - START
                          if(unique){
			  lattices.push_back(tmp_lattice); // stores original original orientation
			  tmp_clatt=GetClat(abc_angles_q2[1],abc_angles_q2[2],abc_angles_q2[3],abc_angles_q2[4],abc_angles_q2[5],abc_angles_q2[6]);
			  clattices.push_back(tmp_clatt); // store Cartesian lattice, alignes with XYZ coordinates
			  latt_devs.push_back(tmp_latt_dev);
			  store++;
                        }
                        }
                        //DX 20190318 - now store calc and store clattice - END
			// DX 20190318 [OBSOLETE] else if(unique && optimize_match && tmp_latt_dev <= 0.2){
			else if(optimize_match && tmp_latt_dev <= 0.2){ //DX 20190318 - removed unique since it doesn't exist yet
			            bool unique = true;
			            for(uint t=0;t<lattices.size();t++){
			              if(identical(tmp_lattice,lattices[t],1e-10)){
			    // DX TEST - CANNOT DO: Eliminates potential matches for(uint t=0;t<clattices.size();t++){
			    // DX TEST - CANNOT DO: Eliminates potential matches if(identical(tmp_clatt,clattices[t],1e-10)){
			    unique=false;
			    break;
			  }
			}
                        //DX 20190318 - now store calc and store clattice - START
                          if(unique){
			  lattices.push_back(tmp_lattice); // stores original original orientation
			  tmp_clatt=GetClat(abc_angles_q2[1],abc_angles_q2[2],abc_angles_q2[3],abc_angles_q2[4],abc_angles_q2[5],abc_angles_q2[6]);
			  clattices.push_back(tmp_clatt); // store Cartesian lattice, alignes with XYZ coordinates
			  latt_devs.push_back(tmp_latt_dev);
			  store++;
                          }
                        //DX 20190318 - now store calc and store clattice - END
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    return true;
  }
}

// ***************************************************************************
// checkLatticeDeviation
// ***************************************************************************
namespace compare{
  double checkLatticeDeviation(double& xstr1_vol, xmatrix<double>& q2,vector<double>& D1,vector<double>& F1){
    double scale=xstr1_vol/(aurostd::abs(aurostd::det(q2)));
    scale=pow(scale,0.3333);
    vector<double> D2,F2;
    cellDiagonal(q2,D2,F2,scale);
    double latt_dev=latticeDeviation(D1,D2,F1,F2);
    return latt_dev;
  }
}

// ***************************************************************************
// Internal structure
// ***************************************************************************
namespace compare{
  bool structureSearch(const string& lfa, 
                        const vector<double>& all_nn1, 
			const xstructure& xstr, 
			vector<xstructure>& vprotos, xstructure& xstr1, const xstructure& xstr2, 
			const int& type_match, double& possible_minMis,
			vector<xmatrix<double> >& lattices,
			vector<xmatrix<double> >& clattices, 
			vector<double>& latt_devs, 
			const bool& optimize_match){ 


    bool LDEBUG=(false || XHOST.DEBUG);

    double mis=1;  
    xstructure proto;
    int flag=0;
    xstructure xstr2_tmp = xstr2;

    for(uint p=0;p<lattices.size();p++){
      if(LDEBUG) {
        cerr << "compare::structureSearch: Trying lattice " << p << endl;
      }
      proto=xstr;
      proto.lattice=lattices[p];

      // Transform
      for(uint iat=0;iat<proto.atoms.size();iat++){
	proto.atoms[iat].fpos=C2F(proto.lattice,proto.atoms[iat].cpos);
      }
      proto.lattice=clattices[p];
	for(uint iat=0;iat<proto.atoms.size();iat++){
	proto.atoms[iat].cpos=F2C(proto.lattice,proto.atoms[iat].fpos);
      }
      xstructure proto_new;
      proto_new.title=proto.title;
      proto_new.lattice=clattices[p];

      // DX NEW - START =======================
      xmatrix<double> f2c = trasp(proto.lattice);
      xmatrix<double> c2f = aurostd::inverse(trasp(proto.lattice));
      bool skew = false;
      double tol=0.01;
      deque<_atom> new_basis;
      for(uint j=0;j<proto.atoms.size();j++){
	if(new_basis.size()==0){
	  proto.atoms[j].fpos = BringInCell(proto.atoms[j].fpos,1e-10);
	  proto.atoms[j].cpos = f2c*proto.atoms[j].fpos;
	  new_basis.push_back(proto.atoms[j]);
	  //proto_new.AddAtom(proto.atoms[j]);
	}
	else {
	  bool duplicate_lattice_point=false;
	  for(uint a=0; a<new_basis.size(); a++){
	    xvector<double> tmp = BringInCell(proto.atoms[j].fpos,1e-10);
	    if(SYM::MapAtom(new_basis[a].fpos,tmp,proto.lattice,proto.f2c,skew,tol)){ //DX 20190619 - lattice and f2c as input
	      duplicate_lattice_point=true;
	      break;
	    }
	  }
	  if(duplicate_lattice_point==false){
	    proto.atoms[j].fpos = BringInCell(proto.atoms[j].fpos,1e-10);
	    proto.atoms[j].cpos = f2c*proto.atoms[j].fpos;
	    new_basis.push_back(proto.atoms[j]);
	    //proto_new.AddAtom(proto.atoms[j]);
	  }
	}
      }
      proto_new.atoms = new_basis;
      proto_new.BringInCell(1e-10); 
      proto_new.FixLattices();
      proto_new.SpeciesPutAlphabetic();
      deque<int> sizes = SYM::arrange_atoms(new_basis);
      proto_new = pflow::SetNumEachType(proto_new, sizes);
      proto = proto_new;
      if(sameSpecies(proto,xstr1,false)){
	vector<double> all_nn_proto;
	bool all_nn_calculated = false;
	for(uint iat=0; iat<proto.atoms.size();iat++){
	  if(proto.atoms[iat].name==lfa){
	    proto.ShifOriginToAtom(iat);
	    proto.BringInCell(1e-10);
	    vector<uint> im1, im2;
	    vector<double> min_dists;
	    if(findMatch(xstr1,proto,im1,im2,min_dists,type_match)){;
	      double cd, f;
	      // Only calculate the NN for the proto if we found suitable matches.  
	      // Only calculate once, nothing changes between shifts to origin (affine)
	      if(!all_nn_calculated){
                all_nn_proto = computeNearestNeighbors(proto);
		if(LDEBUG) {
		  cerr << "compare::structureSearch: Nearest neighbors:" << endl;
		  for(uint a=0;a<all_nn_proto.size();a++){
		    cerr << "compare::structureSearch: Nearest neighbor distance from " << a << " atom: " << all_nn_proto[a] << endl;
		  }
		}
		all_nn_calculated = true;
              }
	      coordinateDeviation(xstr1,proto,all_nn1,all_nn_proto,im1,im2,min_dists,cd,f);
	      mis=computeMisfit(latt_devs[p],cd,f);
	      //if(LDEBUG) {
	      //  cerr << "mis,latt_dev,cd,f: " << mis << ", " <<latt_devs[p] << ", " << cd << ", " << f <<  endl;
	      //}
	      if(flag==0){
	        flag=1;
		//cerr << "storing: " << proto << endl;
		vprotos.push_back(proto);
		possible_minMis=mis;
	      }
	      else {
	        if(mis<possible_minMis){
		  vprotos.clear();
		  possible_minMis=mis;
		  //cerr << "storing: " << proto << endl;
		  vprotos.push_back(proto); //to here
		}
	      }
              // If we want to simply find a match and not find the best match, we can exit early
	      if(mis<0.1 && !optimize_match) {
                return true;
	        //DEBUGGING
		//cerr <<"Winning combo: "<<i<<","<<j<<","<<k<<","<<w<<endl;
		//cerr << "proto.lattice: " << proto.lattice << endl;
		//cerr << "lattice(1): " << modulus(proto.lattice(1)) << endl;
		//cerr << "lattice(2): " << modulus(proto.lattice(2)) << endl;
		//cerr << "lattice(3): " << modulus(proto.lattice(3)) << endl;
	      }
	    }
	  }
	}
      }// end of if protos.size()...
    }
    return true;
  }
}

//---------------------------------------------------------------

// ***************************************************************************
//                                 END - FINE
// 			AFLOW Compare Structure - Functions
//		David Hicks (d.hicks@duke.edu) and Carlo de Santo
// ***************************************************************************
