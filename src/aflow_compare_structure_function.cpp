// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                                                                         *
// ***************************************************************************
// AFLOW-XTAL-MATCH (compare crystal structures) - Functions
// Written by David Hicks (david.hicks@duke.edu) 
// Contributors: Carlo De Santo

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
  structure_representative_name="";
  structure_representative_compound="";
  structure_representative.Clear();
  structure_representative_generated=false;
  structure_representative_from="";
  number_compounds_matching_representative=0;
  number_of_types=0;
  elements.clear();
  stoichiometry.clear();
  number_of_atoms=0;
  atom_decorations_equivalent.clear();
  Pearson="";
  space_group=0;
  grouped_Wyckoff_positions.clear();
  wyckoff_site_symmetry.clear();
  wyckoff_multiplicity.clear();
  wyckoff_letter.clear();
  aflow_label=""; //DX 20190724
  aflow_parameter_list.clear(); //DX 20190724
  aflow_parameter_values.clear(); //DX 20190724
  matching_aflow_prototypes.clear(); //DX 20190724
  environments_LFA.clear(); //DX 20190711
  structures_duplicate_names.clear();
  structures_duplicate_compounds.clear();
  structures_duplicate.clear();
  structures_duplicate_generated.clear();
  structures_duplicate_from.clear();
  structures_duplicate_grouped_Wyckoff_positions.clear(); //DX 20190813
  number_compounds_matching_duplicate.clear();
  duplicate_comparison_logs.clear(); //DX 20190506
  structures_family_names.clear();
  structures_family.clear();
  structures_family_generated.clear();
  structures_family_from.clear();
  structures_family_grouped_Wyckoff_positions.clear(); //DX 20190813
  number_compounds_matching_family.clear();
  family_comparison_logs.clear(); //DX 20190506
  misfits_duplicate.clear();
  misfits_family.clear();
  property_names.clear();
  property_units.clear();
  properties_structure_representative.clear();
  properties_structures_duplicate.clear();
  properties_structures_family.clear(); //DX 20190425
}

// ===== Free  ===== //
void StructurePrototype::free(){
}

// ===== Destructor ===== //
StructurePrototype::~StructurePrototype(){ 
  structure_representative.Clear();
  elements.clear();
  stoichiometry.clear();
  atom_decorations_equivalent.clear();
  grouped_Wyckoff_positions.clear();
  wyckoff_site_symmetry.clear();
  wyckoff_multiplicity.clear();
  wyckoff_letter.clear();
  aflow_parameter_list.clear(); //DX 20190724
  aflow_parameter_values.clear(); //DX 20190724
  matching_aflow_prototypes.clear(); //DX 20190724
  environments_LFA.clear(); //DX 20190711
  structures_duplicate_names.clear();
  structures_duplicate_compounds.clear();
  structures_duplicate.clear();
  structures_duplicate_generated.clear();
  structures_duplicate_from.clear();
  structures_duplicate_grouped_Wyckoff_positions.clear(); //DX 20190813
  number_compounds_matching_duplicate.clear();
  duplicate_comparison_logs.clear(); //DX 20190506
  structures_family_names.clear();
  structures_family.clear();
  structures_family_generated.clear();
  structures_family_from.clear();
  structures_family_grouped_Wyckoff_positions.clear(); //DX 20190813
  number_compounds_matching_family.clear();
  family_comparison_logs.clear(); //DX 20190506
  misfits_duplicate.clear();
  misfits_family.clear();
  property_names.clear();
  property_units.clear();
  properties_structure_representative.clear();
  properties_structures_duplicate.clear();
  properties_structures_family.clear(); //DX 20190425
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
    structure_representative_name=b.structure_representative_name; 
    structure_representative_compound=b.structure_representative_compound; 
    structure_representative=b.structure_representative; 
    structure_representative_generated=b.structure_representative_generated; 
    structure_representative_from=b.structure_representative_from; 
    number_compounds_matching_representative=b.number_compounds_matching_representative;
    number_of_types=b.number_of_types;
    elements=b.elements;
    stoichiometry=b.stoichiometry;
    number_of_atoms=b.number_of_atoms;
    atom_decorations_equivalent=b.atom_decorations_equivalent;
    Pearson=b.Pearson;
    space_group=b.space_group;
    grouped_Wyckoff_positions=b.grouped_Wyckoff_positions;
    wyckoff_site_symmetry=b.wyckoff_site_symmetry;
    wyckoff_multiplicity=b.wyckoff_multiplicity;
    wyckoff_letter=b.wyckoff_letter;
    aflow_label=b.aflow_label; //DX 20190724
    aflow_parameter_list=b.aflow_parameter_list; //DX 20190724
    aflow_parameter_values=b.aflow_parameter_values; //DX 20190724
    matching_aflow_prototypes=b.matching_aflow_prototypes; //DX 20190724
    environments_LFA=b.environments_LFA; //DX 20190711
    structures_duplicate_names=b.structures_duplicate_names;
    structures_duplicate_compounds=b.structures_duplicate_compounds;
    structures_duplicate=b.structures_duplicate;
    structures_duplicate_generated=b.structures_duplicate_generated;
    structures_duplicate_from=b.structures_duplicate_from;
    structures_duplicate_grouped_Wyckoff_positions=b.structures_duplicate_grouped_Wyckoff_positions; //DX 20190813
    number_compounds_matching_duplicate=b.number_compounds_matching_duplicate;
    duplicate_comparison_logs=b.duplicate_comparison_logs; //DX 20190506
    structures_family_names=b.structures_family_names;
    structures_family=b.structures_family;
    structures_family_generated=b.structures_family_generated;
    structures_family_from=b.structures_family_from;
    structures_family_grouped_Wyckoff_positions=b.structures_family_grouped_Wyckoff_positions; //DX 20190813
    number_compounds_matching_family=b.number_compounds_matching_family;
    family_comparison_logs=b.family_comparison_logs; //DX 20190506
    misfits_duplicate=b.misfits_duplicate;
    misfits_family=b.misfits_family;
    property_names=b.property_names;
    property_units=b.property_units;
    properties_structure_representative=b.properties_structure_representative;
    properties_structures_duplicate=b.properties_structures_duplicate;
    properties_structures_family=b.properties_structures_family; //DX 20190425
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
    
    // structure_representative 
    sscontent_json << "\"structure_representative\":\"" << StructurePrototype.structure_representative_name << "\"" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // number_compounds_matching_representative 
    sscontent_json << "\"number_compounds_matching_representative\":" << StructurePrototype.number_compounds_matching_representative << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // number_of_types 
    sscontent_json << "\"number_of_types\":" << StructurePrototype.number_of_types << eendl;
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
    
    // atom_decorations_equivalent
    if(StructurePrototype.atom_decorations_equivalent.size()!=0){ //DX 20190425 - only print if calculated
      //DX 20191111 [OBSOLETE] sscontent_json << "\"atom_decorations_equivalent\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.atom_decorations_equivalent,"\""),",") << "]" << eendl;
      sscontent_json << "\"atom_decorations_equivalent\":[";
      tmp.clear();
      for(uint i=0;i<StructurePrototype.atom_decorations_equivalent.size();i++){
        tmp.push_back("["+aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.atom_decorations_equivalent[i],"\""),",")+"]");
      }
      sscontent_json << aurostd::joinWDelimiter(tmp,",") << "]" << eendl;
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
    
    if(StructurePrototype.aflow_label.size()!=0){
      // aflow_label 
      sscontent_json << "\"aflow_label\":\"" << StructurePrototype.aflow_label << "\"" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // aflow_parameter_list
      sscontent_json << "\"aflow_parameter_list\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.aflow_parameter_list,"\""),",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // aflow_parameter_values
      sscontent_json << "\"aflow_parameter_values\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(StructurePrototype.aflow_parameter_values,8,roff),",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }

    // matching_aflow_prototypes
    if(StructurePrototype.matching_aflow_prototypes.size()!=0){
      sscontent_json << "\"matching_aflow_prototypes\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.matching_aflow_prototypes,"\""),",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    }

    // environments_LFA
    sscontent_json << "\"environments_LFA\":[";
    tmp.clear();
    for(uint i=0;i<StructurePrototype.environments_LFA.size();i++){
      stringstream ss_tmp;
      ss_tmp << StructurePrototype.environments_LFA[i];
      tmp.push_back(ss_tmp.str());
    }
    sscontent_json << aurostd::joinWDelimiter(tmp,",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // structures_duplicate
    sscontent_json << "\"structures_duplicate\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.structures_duplicate_names,"\""),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // misfits_duplicate
    sscontent_json << "\"misfits_duplicate\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(StructurePrototype.misfits_duplicate,8,roff),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // number_compounds_matching_duplicate
    sscontent_json << "\"number_compounds_matching_duplicate\":[" << aurostd::joinWDelimiter(StructurePrototype.number_compounds_matching_duplicate,",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // structures_family
    sscontent_json << "\"structures_family\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.structures_family_names,"\""),",") << "]" << eendl;
    vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
    
    // misfits_family
    sscontent_json << "\"misfits_family\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(StructurePrototype.misfits_family,8,roff),",") << "]" << eendl;
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

      // properties_structure_representative
      sscontent_json << "\"properties_structure_representative\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.properties_structure_representative,"\""),",") << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      // duplicate_structure_properties
      sscontent_json << "\"properties_structures_duplicate\":[";
      tmp.clear();
      //DX 20190326 - should be duplicate for(uint i=0;i<StructurePrototype.properties_structure_representative.size();i++){
      for(uint i=0;i<StructurePrototype.properties_structures_duplicate.size();i++){ // DX 20190326 - correctly changed to duplicate_sturctures_properties
        tmp.push_back(aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.properties_structures_duplicate[i],"\""),","));
      }
      sscontent_json << aurostd::joinWDelimiter(aurostd::wrapVecEntries(tmp,"[","]"),",");
      sscontent_json << "]" << eendl;
      vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

      //DX 20190425 - START
      // family_structure_properties
      sscontent_json << "\"properties_structures_family\":[";
      tmp.clear();
      for(uint i=0;i<StructurePrototype.properties_structures_family.size();i++){
        tmp.push_back(aurostd::joinWDelimiter(aurostd::wrapVecEntries(StructurePrototype.properties_structures_family[i],"\""),","));
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
  for(uint i=0;i<misfits_duplicate.size();i++){
    if(misfits_duplicate[i]<=0.1 && (misfits_duplicate[i]+1.0)>1e-3 ){
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

  return structures_duplicate_names.size(); 
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
// StructurePrototype::isSymmetryCalculated() 
// ***************************************************************************
bool StructurePrototype::isLFAEnvironmentCalculated(){

  // Check if the LFA environment is calculated

  if(environments_LFA.size()==0){return false;}

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
  space_group = structure_representative.SpaceGroup_ITC();
  vector<GroupedWyckoffPosition> tmp_grouped_Wyckoff_positions; 
  compare::groupWyckoffPositions(structure_representative, tmp_grouped_Wyckoff_positions);
  grouped_Wyckoff_positions = tmp_grouped_Wyckoff_positions;
  return true;
}

// ***************************************************************************
// StructurePrototype::addStructurePrototypeAsDuplicate() 
// ***************************************************************************
bool StructurePrototype::addStructurePrototypeAsDuplicate(StructurePrototype& b){

  // Add representative structure as a duplicate structure in this object
  
  // add structure info
  structures_duplicate_names.push_back(b.structure_representative_name);
  structures_duplicate_generated.push_back(b.structure_representative_generated);
 
  // only add xstructure if it has been generated
  if(b.structure_representative_generated){
    structures_duplicate.push_back(b.structure_representative);
    structures_duplicate_compounds.push_back(compare::getCompoundName(b.structure_representative)); //DX 20190111 - added compound, e.g., Ag1Br2
  }
  else if(!b.structure_representative_compound.empty()){
    structures_duplicate_compounds.push_back(b.structure_representative_compound); //DX 20190111 - added compound, e.g., Ag1Br2
  }
  else {
    structures_duplicate_compounds.push_back(""); //DX 20190111 - added compound, e.g., Ag1Br2
  }
  structures_duplicate_from.push_back(b.structure_representative_from);
  // add Wyckoff positions
  structures_duplicate_grouped_Wyckoff_positions.push_back(b.grouped_Wyckoff_positions); //DX 20190814 - added Wyckoff information for duplicates
  number_compounds_matching_duplicate.push_back(b.number_compounds_matching_representative); //DX 20190228 - added duplicate compound count

  // add structure properties
  if(b.properties_structure_representative.size()!=0){
    properties_structures_duplicate.push_back(b.properties_structure_representative); //DX 20181218 - added property_values
  }

  // signifies comparison has not been performed yet
  misfits_duplicate.push_back(-1.0);

  return true;
}

// ***************************************************************************
// StructurePrototype::addStructurePrototypeAsFamilyStructure() 
// ***************************************************************************
void StructurePrototype::putDuplicateAsFamily(uint index, bool keep_generated){

  // Move duplicate structure information to the same family attributes in this 
  // StructurePrototype object
            
  // add structure info
  structures_family_names.push_back(structures_duplicate_names[index]);
  structures_family_grouped_Wyckoff_positions.push_back(structures_duplicate_grouped_Wyckoff_positions[index]);
  
  // add misfit
  misfits_family.push_back(misfits_duplicate[index]);
  if(duplicate_comparison_logs.size()!=0){
    family_comparison_logs.push_back(duplicate_comparison_logs[index]);
  }

  // generated structure info
  structures_family_from.push_back(structures_duplicate_from[index]);
  if(keep_generated && structures_duplicate_generated[index]){
    structures_family.push_back(structures_duplicate[index]);
    structures_family_generated.push_back(structures_duplicate_generated[index]);
  }
  else{ //default
    structures_family_generated.push_back(false);
  }

  // store properties of family structures
  if(property_names.size()!=0){
    properties_structures_family.push_back(properties_structures_duplicate[index]);
  }
  
  return;
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
  number_of_types=b.stoichiometry.size();
  Pearson=b.Pearson;
  space_group=b.space_group;
  grouped_Wyckoff_positions=b.grouped_Wyckoff_positions;
  environments_LFA=b.environments_LFA; //DX 20190711
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
  // NOTE: what about elements???: No, don't want to propogate elements which do not apply to 
  // new prototype system
  
  // load structure info
  structure_representative_name=b.structures_duplicate_names[index];
  // number_of_types=b.structures_duplicate[index].num_each_type.size();
  structure_representative_compound=b.structures_duplicate_compounds[index]; //DX 20190111 - added compound, e.g., Ag1Br2
  elements = pflow::getElements(structure_representative_compound); //DX 20191003 - add elements of this compound only!
  structure_representative_generated=b.structures_duplicate_generated[index];
  structure_representative_from=b.structures_duplicate_from[index];
  grouped_Wyckoff_positions=b.structures_duplicate_grouped_Wyckoff_positions[index]; //DX 20190814
  number_compounds_matching_representative=b.number_compounds_matching_duplicate[index]; //DX 20190228 - added duplicate compound count
  if(b.structures_duplicate_generated[index]){
    number_of_atoms=b.structures_duplicate[index].atoms.size();
    structure_representative=b.structures_duplicate[index];
  }
  
  // load property info
  if(b.property_names.size()!=0){
    property_names=b.property_names;
    property_units=b.property_units;
    properties_structure_representative=b.properties_structures_duplicate[index];
  }
  return true;
}

/*
// ***************************************************************************
// StructurePrototype::putRepresentativeAsDuplicate() 
// ***************************************************************************
bool StructurePrototype::putRepresentativeAsDuplicate(StructurePrototype& b){

  // Make structure in representative the new duplicate structure
  
  // load structure info
  structures_duplicate_names.push_back(b.structure_representative_name);
  structures_duplicate_compounds.push_back(b.structure_representative_compound); //DX 20190111 - added compound, e.g., Ag1Br2
  structures_duplicate_generated.push_back(b.structure_representative_generated);
  structures_duplicate_from.push_back(b.structure_representative_from);
  number_compounds_matching_duplicate.push_back(b.number_compounds_matching_representative); //DX 20190228 - added duplicate compound count
  if(b.structure_representative_generated){
    structures_duplicate.push_back(b.structure_representative);
  }
  
  // load property info
  if(b.property_names.size()!=0){
    properties_structures_duplicate.push_back(b.properties_structure_representative);
  }

  return true;
}
*/

// ***************************************************************************
// StructurePrototype::copyDuplicate() 
// ***************************************************************************
bool StructurePrototype::copyDuplicate(StructurePrototype& b, uint& index, bool copy_misfit){ //DX 20190730 - added copy_misfit option

  // Copy structure at index as a potential duplicate in this object 

  structures_duplicate_names.push_back(b.structures_duplicate_names[index]);
  structures_duplicate_compounds.push_back(b.structures_duplicate_compounds[index]);
  structures_duplicate_generated.push_back(b.structures_duplicate_generated[index]);
  structures_duplicate_from.push_back(b.structures_duplicate_from[index]);
  structures_duplicate_grouped_Wyckoff_positions.push_back(b.structures_duplicate_grouped_Wyckoff_positions[index]);
  number_compounds_matching_duplicate.push_back(b.number_compounds_matching_duplicate[index]); //DX 20190228 - added duplicate compound count
  if(b.structures_duplicate_generated[index]){
    structures_duplicate.push_back(b.structures_duplicate[index]);
  }

  if(!copy_misfit){
  // signifies comparison has not been performed yet
  misfits_duplicate.push_back(-1.0);
  }
  else{ //DX 20190730
    misfits_duplicate.push_back(b.misfits_duplicate[index]);
  }
  
  // load property info
  if(property_names.size()!=0){
    properties_structures_duplicate.push_back(b.properties_structures_duplicate[index]);
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
  structures_duplicate_names.erase(structures_duplicate_names.begin()+index);
  structures_duplicate_compounds.erase(structures_duplicate_compounds.begin()+index); //DX 20190111 - added compound, e.g., Ag1Br2
  // DX - may need to be careful here.  If we have a mix of generated and non-generated, we may have difficulties
  if(structures_duplicate_generated[index]){
    structures_duplicate.erase(structures_duplicate.begin()+index);
  }
  structures_duplicate_generated.erase(structures_duplicate_generated.begin()+index);
  structures_duplicate_from.erase(structures_duplicate_from.begin()+index);
  number_compounds_matching_duplicate.erase(number_compounds_matching_duplicate.begin()+index);
  structures_duplicate_grouped_Wyckoff_positions.erase(structures_duplicate_grouped_Wyckoff_positions.begin()+index); //DX 20190814

  // remove misfit
  misfits_duplicate.erase(misfits_duplicate.begin()+index);

  //DX 20190504 - START
  // remove comparison log
  if(duplicate_comparison_logs.size()!=0){
    duplicate_comparison_logs.erase(duplicate_comparison_logs.begin()+index);
  }
  //DX 20190504 - END

  // remove properties
  if(property_names.size()!=0){
    properties_structures_duplicate.erase(properties_structures_duplicate.begin()+index);
  }
  return true;
} 

// ***************************************************************************
// StructurePrototype::clearDuplicateInformation() 
// ***************************************************************************
bool StructurePrototype::removeDuplicates(bool remove_duplicate_count){

  // Remove duplicate structure information from object
  
  // remove structure information
  structures_duplicate_names.clear();
  structures_duplicate_compounds.clear(); //DX 20190111 - added compound, e.g., Ag1Br2
  structures_duplicate.clear();
  structures_duplicate_generated.clear();
  structures_duplicate_from.clear();

  // may want to keep this information
  if(remove_duplicate_count){number_compounds_matching_duplicate.clear();}

  // remove comparison logs 
  duplicate_comparison_logs.clear(); //DX 20190506

  // remove misfit
  misfits_duplicate.clear();

  // remove properties
  properties_structures_duplicate.clear();

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
// AtomEnvironment Class 
// ***************************************************************************
// ===== Constructor ===== //
AtomEnvironment::AtomEnvironment(){
  center_element="";
  center_type=0;
  neighbor_elements.clear();
  neighbor_types.clear();
  neighbor_distances.clear();
  neighbor_frequencies.clear();
  neighbor_coordinates.clear();
}

// ===== Free  ===== //
void AtomEnvironment::free(){
}

// ===== Destructor ===== //
AtomEnvironment::~AtomEnvironment(){
  neighbor_elements.clear();
  neighbor_types.clear();
  neighbor_distances.clear();
  neighbor_frequencies.clear();
  neighbor_coordinates.clear();
  free();
}

// ===== Copy Constructor ===== //
AtomEnvironment::AtomEnvironment(const AtomEnvironment& b){
  copy(b);
}

// ===== Copy Constructor Function ===== //
void AtomEnvironment::copy(const AtomEnvironment& b) {
  if(this != &b){
    center_element=b.center_element;
    center_type=b.center_type;
    neighbor_elements=b.neighbor_elements;
    neighbor_types=b.neighbor_types;
    neighbor_distances=b.neighbor_distances;
    neighbor_frequencies=b.neighbor_frequencies;
    neighbor_coordinates=b.neighbor_coordinates;
  }
}

// ===== Assignment Operator (operator=) ===== //
// ***************************************************************************
// AtomEnvironment::operator=
// ***************************************************************************
const AtomEnvironment& AtomEnvironment::operator=(const AtomEnvironment& b){
  if(this!=&b){
    free();
    copy(b);
  }
  return *this;
}

// ===== Output Handler ===== //
// ***************************************************************************
// AtomEnvironment::operator<< 
// ***************************************************************************
ostream& operator<<(ostream& oss, const AtomEnvironment& AtomEnvironment){
  
  // operator<< for the StruturePrototype object (looks like a JSON)

  //if(AtomEnvironment.iomode!=JSON_MODE){ //A safeguard until we construct more output schemes.
  //  AtomEnvironment.iomode=JSON_MODE;
  //}

  string eendl="";
  bool roff=true; //round off
  stringstream sscontent_json;
  vector<string> vcontent_json, tmp;
  
  // center_element 
  sscontent_json << "\"center_element\":\"" << AtomEnvironment.center_element << "\"" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // center_type
  sscontent_json << "\"center_type\":" << AtomEnvironment.center_type << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // neighbor_elements
  sscontent_json << "\"neighbor_elements\":[" << aurostd::joinWDelimiter(aurostd::wrapVecEntries(AtomEnvironment.neighbor_elements,"\""),",") << "]" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // neighbor_distances
  sscontent_json << "\"neighbor_distances\":[" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(AtomEnvironment.neighbor_distances,8,true,1e-6),",") << "]" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // neighbor_frequencies
  sscontent_json << "\"neighbor_frequencies\":[" << aurostd::joinWDelimiter(AtomEnvironment.neighbor_frequencies,",") << "]" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");
  
  // neighbor_coordinates
  sscontent_json << "\"neighbor_coordinates\":[";
  vector<string> coordinate_sets;
  for(uint i=0;i<AtomEnvironment.neighbor_coordinates.size();i++){
    vector<string> coordinates;
    for(uint j=0;j<AtomEnvironment.neighbor_coordinates[i].size();j++){
      stringstream ss_tmp; ss_tmp << "[" << aurostd::joinWDelimiter(aurostd::xvecDouble2vecString(AtomEnvironment.neighbor_coordinates[i][j],8,true,1e-6),",") << "]";
      coordinates.push_back(ss_tmp.str());
      ss_tmp.clear();
    }
    coordinate_sets.push_back("["+aurostd::joinWDelimiter(coordinates,",")+"]");
  }
  sscontent_json << aurostd::joinWDelimiter(coordinate_sets,",") << "]" << eendl;
  vcontent_json.push_back(sscontent_json.str()); sscontent_json.str("");

  // Put into json AtomEnvironment object
  oss << "{" << aurostd::joinWDelimiter(vcontent_json,",")  << "}";
  vcontent_json.clear();

  return oss;
}


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
  vector<StructurePrototype> loadStructuresFromDirectory(string& directory, vector<string>& magmoms_for_systems, 
  bool same_species, ofstream& FileMESSAGE){ //DX 20190319 - added FileMESSAGE

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
        message << "Ignoring file=" << vfiles[i] << endl;
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);
        vfiles.erase(vfiles.begin()+i);
        i--;
      }
      else {
        StructurePrototype structure_tmp;
        stringstream sss1;
        aurostd::efile2stringstream(directory+"/"+vfiles[i],sss1);
        xstructure xstr1; //DX 20190718
        try { xstr1 = sss1; } //DX 20190718
        catch(aurostd::xerror& excpt) { message << "Could not load structure " << vfiles[i] << "...skipping structure"; pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);  continue; } //DX 20190718
        xstr1.directory = directory+"/"+vfiles[i]; //DX 20190718 - need to pass in since passing structure as a string
        if(magmoms_for_systems.size()==vfiles.size()){
          try { pflow::ProcessAndAddSpinToXstructure(xstr1, magmoms_for_systems[i]); } //DX 20190801
          catch(aurostd::xerror& excpt) { message << "Magnetic information could not be loaded (" << magmoms_for_systems[i] << "...skipping structure"; pflow::logger(function_name, message, FileMESSAGE, _LOGGER_WARNING_); continue; } //DX 20190801
        }
        structure_tmp.structure_representative = xstr1;
        structure_tmp.structure_representative.ReScale(1.0); //DX 20190715
        structure_tmp.structure_representative_name = directory+"/"+vfiles[i];
        structure_tmp.stoichiometry = compare::getStoichiometry(xstr1,same_species);
        structure_tmp.elements = compare::getElements(xstr1);
        structure_tmp.number_of_atoms = xstr1.atoms.size(); //DX 20190425
        structure_tmp.number_of_types = xstr1.num_each_type.size(); //DX 20190425
        structure_tmp.structure_representative_compound = compare::getCompoundName(xstr1); //remove ones is true  //DX 20190311 //DX 20190313 - use xstr1
        // update xstructure species
        if(structure_tmp.structure_representative.species.size()==0){
          deque<string> deque_species; for(uint j=0;j<structure_tmp.elements.size();j++){deque_species.push_back(structure_tmp.elements[j]);}
          structure_tmp.structure_representative.SetSpecies(deque_species);
          structure_tmp.structure_representative.SpeciesPutAlphabetic();
        }
        // clean species
        else{
          for(uint s=0;s<structure_tmp.structure_representative.species.size();s++){structure_tmp.structure_representative.species[s]=KBIN::VASP_PseudoPotential_CleanName(structure_tmp.structure_representative.species[s]); } //DX 20190711
          structure_tmp.structure_representative.SetSpecies(structure_tmp.structure_representative.species);
        }
        // check if fake names for same species comparison
        if(structure_tmp.structure_representative.species[0]=="A" && same_species){
          message << "Atomic species are missing for " << structure_tmp.structure_representative_name << " cannot perform material comparison; skipping structure.";     
          pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);
          continue;
        }
        //DX 20191105 [MOVED LATER - SAME AS SYMMETRY] structure_tmp.environments_LFA= compare::computeLFAEnvironment(xstr1); //DX 20190711
        structure_tmp.structure_representative_generated = true; 
        structure_tmp.structure_representative_from = "file"; 
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
  vector<StructurePrototype> loadStructuresFromFile(string& filename, vector<string>& magmoms_for_systems, 
  bool same_species, ofstream& FileMESSAGE){ //DX 20190319 - added FileMESSAGE

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
        xstructure xstr1; //DX 20190718
        try { xstr1 = geometry; } //DX 20190718
        catch(aurostd::xerror& excpt) { message << "Could not load structure " << structure_count << "/" << start_string.size() << "...skipping structure"; pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);  continue; } //DX 20190718
        stringstream designation; designation << "file structure # " << structure_count << "/" << start_string.size();
        xstr1.directory = designation.str(); //DX 20190718 - need to pass in since passing structure as a string
        if(magmoms_for_systems.size()==lines.size()){
          try { pflow::ProcessAndAddSpinToXstructure(xstr1, magmoms_for_systems[i]); } //DX 20190801
          catch(aurostd::xerror& excpt) { message << "Magnetic information could not be loaded (" << magmoms_for_systems[i] << "...skipping structure"; pflow::logger(function_name, message, FileMESSAGE, _LOGGER_WARNING_); continue; } //DX 20190801
        }
        structure_tmp.structure_representative = xstr1;
        structure_tmp.structure_representative.ReScale(1.0); //DX 20190715
        structure_tmp.structure_representative_name = designation.str();
        structure_tmp.stoichiometry = compare::getStoichiometry(xstr1,same_species);
        structure_tmp.elements = compare::getElements(xstr1);
        structure_tmp.number_of_atoms = xstr1.atoms.size(); //DX 20190425
        structure_tmp.number_of_types = xstr1.num_each_type.size(); //DX 20190425
        structure_tmp.structure_representative_compound = compare::getCompoundName(xstr1); //remove ones is true  //DX 20190311 //DX 20190313 - use xstr
        // update xstructure species
        if(structure_tmp.structure_representative.species.size()==0){
          deque<string> deque_species; for(uint j=0;j<structure_tmp.elements.size();j++){deque_species.push_back(structure_tmp.elements[j]);}
          structure_tmp.structure_representative.SetSpecies(deque_species);
          structure_tmp.structure_representative.SpeciesPutAlphabetic();
        }
        // clean species
        else{
          for(uint s=0;s<structure_tmp.structure_representative.species.size();s++){structure_tmp.structure_representative.species[s]=KBIN::VASP_PseudoPotential_CleanName(structure_tmp.structure_representative.species[s]); } //DX 20190711
          structure_tmp.structure_representative.SetSpecies(structure_tmp.structure_representative.species);
        }
        // check if fake names for same species comparison
        if(structure_tmp.structure_representative.species[0]=="A" && same_species){
          message << "Atomic species are missing for " << structure_tmp.structure_representative_name << " cannot perform material comparison; skipping strucutre.";     
          pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);
          continue;
        }
        //DX 20191105 [MOVED LATER - SAME AS SYMMETRY] structure_tmp.environments_LFA=compare::computeLFAEnvironment(structure_tmp.structure_representative); //DX 20190711
        structure_tmp.structure_representative_generated = true; 
        structure_tmp.structure_representative_from = input_file.str(); 
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
  vector<StructurePrototype> loadStructuresFromStructureList(vector<string>& filenames, vector<string>& magmoms_for_systems, 
  bool same_species, ofstream& FileMESSAGE){

    // load all structures from a vector of filenames into a vector of StructurePrototype object

    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "compare:loadStructuresFromStructureList()";
    stringstream message;
    ostream& logstream = cout;
    //DX [OBSOLETE] ofstream FileMESSAGE;
    
    vector<StructurePrototype> all_structures;
  
    // ---------------------------------------------------------------------------
    // read in files
    for(uint i=0;i<filenames.size();i++){
      StructurePrototype structure_tmp;
      if(!aurostd::FileExist(filenames[i])){
        message << filenames[i] << " file not found.";
        throw aurostd::xerror(function_name,message,_FILE_NOT_FOUND_);
      }
      stringstream sss;
      aurostd::efile2stringstream(filenames[i],sss);
      xstructure xstr(sss);  
      xstr.directory = filenames[i]; //DX 20190718 - need to pass in since passing structure as a string
      if(magmoms_for_systems.size()==filenames.size()){
        try { pflow::ProcessAndAddSpinToXstructure(xstr, magmoms_for_systems[i]); } //DX 20190801
        catch(aurostd::xerror& excpt) { message << "Magnetic information could not be loaded (" << magmoms_for_systems[i] << "."; throw aurostd::xerror(function_name, message, _INPUT_ERROR_); } //DX 20190801
      }
      structure_tmp.structure_representative = xstr;
      structure_tmp.structure_representative.ReScale(1.0); //DX 20190715
      structure_tmp.structure_representative_name = filenames[i];
      structure_tmp.stoichiometry = compare::getStoichiometry(xstr,same_species);
      structure_tmp.elements = compare::getElements(xstr);
      structure_tmp.number_of_atoms = xstr.atoms.size(); //DX 20190425
      structure_tmp.number_of_types = xstr.num_each_type.size(); //DX 20190425
	    structure_tmp.structure_representative_compound = compare::getCompoundName(xstr); //remove ones is true  //DX 20190311 //DX 20190313 - use xstr
      // update xstructure species
      if(structure_tmp.structure_representative.species.size()==0){
        deque<string> deque_species; for(uint j=0;j<structure_tmp.elements.size();j++){deque_species.push_back(structure_tmp.elements[j]);}
        structure_tmp.structure_representative.SetSpecies(deque_species);
        structure_tmp.structure_representative.SpeciesPutAlphabetic();
      }
      // clean species
      else{
        for(uint s=0;s<structure_tmp.structure_representative.species.size();s++){structure_tmp.structure_representative.species[s]=KBIN::VASP_PseudoPotential_CleanName(structure_tmp.structure_representative.species[s]); } //DX 20190711
        structure_tmp.structure_representative.SetSpecies(structure_tmp.structure_representative.species);
      }
      // check if fake names for same species comparison
      if(structure_tmp.structure_representative.species[0]=="A" && same_species){
        message << "Atomic species are missing for " << structure_tmp.structure_representative_name << " cannot perform material comparison; skipping strucutre.";     
        pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);
        continue;
      }
      //DX 20191105 [MOVED LATER - SAME AS SYMMETRY] structure_tmp.environments_LFA=compare::computeLFAEnvironment(structure_tmp.structure_representative); //DX 20190711
      structure_tmp.structure_representative_generated = true; 
      structure_tmp.structure_representative_from = "file"; 
      if(LDEBUG) {
        cerr << function_name << ": loaded structure " << i << endl;
      }
      all_structures.push_back(structure_tmp);
    }
    return all_structures;
  }
}
//DX 20190424 - END

//DX 20191105 - load multiple structures (useful for multithreaded structure loading) - START
// ***************************************************************************
// generateStructures()
// ***************************************************************************
namespace compare {
  void generateStructures(vector<StructurePrototype>& structures, ostream& oss){

    uint start_index = 0;
    uint end_index = structures.size(); // since we use i<=end_index

    generateStructuresInRange(structures, oss, start_index, end_index);
  }

  void generateStructuresInRange(vector<StructurePrototype>& structures, ostream& oss, uint start_index, uint end_index){
    
    for(uint i=start_index;i<end_index;i++){
      structures[i].structure_representative_generated = generateStructure(structures[i].structure_representative_name, structures[i].structure_representative_from,
          structures[i].structure_representative, oss);
    }

  }
}
//DX 20191105 - load multiple structures (useful for multithreaded structure loading) - END 

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

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = "compare::generateStructure()";
    ofstream FileMESSAGE;
    vector<string> tokens;

    if(LDEBUG){
      cerr << function_name << ": generating structure: " << structure_name << " from " << structure_from << endl;
    }

    // ---------------------------------------------------------------------------
    // load from AFLOW prototypes
    if(structure_from=="aflow_prototypes"){
      // htqc or anrl
      structure = aflowlib::PrototypeLibraries(cout,structure_name,"",2); 
    }
    // ---------------------------------------------------------------------------
    // load from AURL
    else if(structure_from=="aurl"){
      aflowlib::_aflowlib_entry entry; entry.aurl = structure_name; 
      //DX 20190326 - need to put url path, i.e., structure name, [OBSOLETE] if(!pflow::loadXstructures(entry,FileMESSAGE,oss,true,structure_name,true)){ cerr << function_name << "WARNING::Could not load structure via aurl..." << endl; return false;}
      //DX ORIG B4 20191105 - if(!pflow::loadXstructures(entry,FileMESSAGE,oss,true,structure_name,true)){ cerr << function_name << "WARNING::Could not load structure via aurl..." << endl; return false;} //DX 20190326
      if(!pflow::loadXstructures(entry,FileMESSAGE,oss)){ cerr << "WARNING::Could not load structure (aurl=" << entry.aurl << ") ... skipping..." << endl; return false;} //DX 20191105
      if(entry.vstr.size()==1){
        structure = entry.vstr[0];
      }
      else {
        cerr << function_name << "::WARNING: More structures loaded than anticipated." << endl;
        return false;
      }
    }
    // ---------------------------------------------------------------------------
    // load from file
    else if(structure_from=="file"){
      stringstream sss;
      aurostd::efile2stringstream(structure_name,sss);
      xstructure xstr(sss);
      structure = xstr;
    }
    // ---------------------------------------------------------------------------
    // load from file containing list of structures (like aflow.in or pocc)
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
    // ---------------------------------------------------------------------------
    // load from input (istream, e.g., from 'cat' or redirect '<') 
    else if(structure_name=="input geometry"){
      stringstream sss; sss << structure_from;
      xstructure xstr(sss);
      structure = xstr;
    }
    // ---------------------------------------------------------------------------
    // load permutation 
    else if(aurostd::substring2bool(structure_from, "permutation of: ")){
      //cerr << "permutation generator: " << structure_from << endl;
      //cerr << "permutation of: " << structure_name << endl;
      string tmp_from = structure_from;
      stringstream sss; sss << aurostd::StringSubst(tmp_from, "permutation of: ", "");
      xstructure xstr(sss);
      deque<string> species; 
      for(uint j=0;j<structure_name.size();j++){stringstream ss_site; ss_site << structure_name[j]; species.push_back(ss_site.str());}
      xstr.SetSpecies(species);
      xstr.SpeciesPutAlphabetic();
      deque<int> sizes = SYM::arrange_atoms(xstr.atoms);
      xstr = pflow::SetNumEachType(xstr, sizes);
      structure = xstr;
    }
    // ---------------------------------------------------------------------------
    // load-type not accounted for 
    else {
      cerr << function_name << "::WARNING: Structure location (from=" << structure_from << ") is not specified correctly for " << structure_name << " (i.e., input, aflow_prototype, aurl, etc.)." << endl;
      return false;
    }

    return true;
  }
}
  
//DX 20191105 - remove non-generated structures - START
// ***************************************************************************
// removeNonGeneratedStructures 
// ***************************************************************************
namespace compare {
	void removeNonGeneratedStructures(vector<StructurePrototype>& structures){

		// remove structures that have not been generated
		// typically used if there is an error in loading a structure

		for(uint i=0;i<structures.size();i++){
			if(!structures[i].structure_representative_generated){
				structures.erase(structures.begin()+i);
				i--;
			}
		}
	}
}
//DX 20191105 - remove non-generated structures - END
  
//DX 20191108 [OBOSLETE] // ***************************************************************************
//DX 20191108 [OBOSLETE] // SVD Decomposition 
//DX 20191108 [OBOSLETE] // ***************************************************************************
//DX 20191108 [OBOSLETE] namespace compare{
//DX 20191108 [OBOSLETE]   bool SVD(xmatrix<double> A){
//DX 20191108 [OBOSLETE]     xmatrix<double> ATA = aurostd::trasp(A)*A;
//DX 20191108 [OBOSLETE]     //to bidiagonal
//DX 20191108 [OBOSLETE]     xmatrix<double> test = ATA;
//DX 20191108 [OBOSLETE]     xmatrix<double> Q = pflow::generalHouseHolderQRDecomposition(test);
//DX 20191108 [OBOSLETE]     //cerr << "A: " << A << endl;
//DX 20191108 [OBOSLETE]     //cerr << "ATA: " << ATA << endl;
//DX 20191108 [OBOSLETE]     //cerr << "Q: " << Q << endl;
//DX 20191108 [OBOSLETE]     //cerr << "R: " << test << endl;
//DX 20191108 [OBOSLETE]     return true;
//DX 20191108 [OBOSLETE]   }
//DX 20191108 [OBOSLETE] }

// ***************************************************************************
// Find ICSD name - Find ICSD name 
// ***************************************************************************
namespace compare{
  string findICSDName(string& name){

    // Find ICSD substring within path name
    // In order for this to work, the following ICSD name format must be 
    // present in the string (e.g., Ag1_ICSD_#####)

    bool LDEBUG=(false || XHOST.DEBUG);

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
      if(aurostd::substring2bool(name,"_ICSD_")){ 
        ICSD_substring = name;
        ICSD_substring_found = true;
      }
    }
    if(!ICSD_substring_found){
      if(LDEBUG){
        cerr << "compare::findICSDName: WARNING: Could not find ICSD substring in name.  representative prototype will not necessarily be the minimum ICSD number." << endl; 
        cerr << "string: " << name << endl;
      }
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
      if(ICSD_entries[i].empty()){ continue; } //DX 20191108 - if not an ICSD, skip
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

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    bool VERBOSE=false;
    string function_name = "compare::comparePermutations()";
    stringstream message;
    ostream& logstream = cout;

    // ---------------------------------------------------------------------------
    // fixed options for permutation comparisons 
    bool same_species = true; // permutation comparisons must compare the same species 
    bool scale_volume = true; // permutations are generated from the same structure, so they will have the same volume anyway
    bool ignore_symmetry = false; // duplicate permutations should have the same symmetry
    bool ignore_Wyckoff = false; // duplicate permutations should have the same Wyckoff positions
    bool ignore_environment = false; // duplicate permutations should have the same environment
    bool single_comparison_round = false; // compare all permutations until matched or exhausted all comparisons
    bool clean_unmatched = true; // remove unmatched structures from object //DX 20190504 
    bool store_comparison_logs = false; // do not store comparison logs //DX 20190822
    bool check_other_grouping = false; // DX 20190830
    bool quiet = true; //true 

    vector<StructurePrototype> final_permutations;

    // ---------------------------------------------------------------------------
    // get stoichiometry
    vector<uint> stoichiometry = compare::getStoichiometry(structure.structure_representative,true);

    // ---------------------------------------------------------------------------
    // calculate symmetry (if not already calculated)
    if(structure.space_group ==0){
      structure.calculateSymmetry();
    }

    // ---------------------------------------------------------------------------
    // generate all permuations structures
    vector<StructurePrototype> permutation_structures = compare::generatePermutationStructures(structure);

    //cerr << "store naming" << endl;
    vector<vector<string> > name_order;
    for(uint i=0;i<permutation_structures.size();i++){
      vector<string> vtmp; 
      for(uint j=0;j<permutation_structures[i].structure_representative_name.size();j++){
        stringstream ss_tmp; ss_tmp << permutation_structures[i].structure_representative_name[j];
        vtmp.push_back(ss_tmp.str());
      }
      name_order.push_back(vtmp);
    }

    // ---------------------------------------------------------------------------
    // loop over grouping modes
    // mode=0: use LFA environment to filter
    // mode=1: if incommensurate groupings, then ignore LFA environment in grouping
    for(uint mode=0;mode<2;mode++){
      if(mode==0){
        if(!quiet || LDEBUG){
          message << "Considering environment analysis in grouping permutations (mode=0)." << endl;
          pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
        }
      }
      if(mode==1){ 
        ignore_environment=true; 
        if(!quiet || LDEBUG){
          message << "Could not find commensurate pemutations when grouping via environment. Ignoring environment analysis in grouping permutations (mode=1)." << endl;
          pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
        }
      }
      final_permutations.clear();

      // ---------------------------------------------------------------------------
      // group comparable permutations
      vector<StructurePrototype> permutation_comparisons;
      permutation_comparisons = compare::groupStructurePrototypes(permutation_structures, same_species, ignore_symmetry, ignore_Wyckoff, ignore_environment, false); //DX 20190731 - added ignore_environment //DX 20190829 - false for duplicates_removed

      // ---------------------------------------------------------------------------
      // ensure the representative stucture is an even permutation
      compare::makeRepresentativeEvenPermutation(permutation_comparisons, name_order);

      if(VERBOSE){ for(uint i=0;i<permutation_comparisons.size();i++){ cerr << "Initial permutation groupings: " << permutation_comparisons[i] << endl; } }

      // ---------------------------------------------------------------------------
      // compare permutations
      final_permutations = compare::runComparisonScheme(num_proc, permutation_comparisons, same_species, check_other_grouping, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, false,store_comparison_logs,oss,FileMESSAGE,quiet); //DX 20190319 - added FileMESSAGE  //DX 20190731 - added ignore_symmetry/Wyckoff/environment //DX 20190822 - add log bool 
      
      // ---------------------------------------------------------------------------
      // check if matched permutations are physically possible
      if(!compare::checkNumberOfGroupings(final_permutations, name_order.size())){
        if(!quiet || LDEBUG){
          message << "Compared groupings of permutations do not follow number theory (# unique=" << final_permutations.size() << " vs # total=" << name_order.size() << ")" << endl;
          // comprehensive output
          if(LDEBUG){ 
            for(uint i=0;i<final_permutations.size();i++){ message << final_permutations[i] << endl; } 
          }
          // minimal output
          else{ 
            for(uint i=0;i<final_permutations.size();i++){ message << final_permutations[i].structure_representative_name << " = " << aurostd::joinWDelimiter(final_permutations[i].structures_duplicate_names,",") << " (misfits_duplicate=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(final_permutations[i].misfits_duplicate,8,true),",")  << ")" << endl; } 
          }
          message << "Trying to check if duplicates match better with other representative structures ... " << endl;
          pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_WARNING_);
        }

        // ---------------------------------------------------------------------------
        // check if better matchings; perhaps matched structures would have smaller misfits if matched to different representatives
        final_permutations = compare::checkForBetterMatches(final_permutations, oss, FileMESSAGE, num_proc, true, same_species, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, clean_unmatched, false, quiet); 

        // ---------------------------------------------------------------------------
        // check if NEW matched permutations are physically possible
        if(!compare::checkNumberOfGroupings(final_permutations, name_order.size())){
          message << "Compared groupings of permutations do not follow number theory (# unique=" << final_permutations.size() << " vs # total=" << name_order.size() << ")" << endl;
          // comprehensive output
          if(LDEBUG){ 
            for(uint i=0;i<final_permutations.size();i++){ message << final_permutations[i] << endl; } 
          }
          // minimal output
          else{ 
            for(uint i=0;i<final_permutations.size();i++){ message << final_permutations[i].structure_representative_name << " = " << aurostd::joinWDelimiter(final_permutations[i].structures_duplicate_names,",") << " (misfits_duplicate=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(final_permutations[i].misfits_duplicate,8,true),",") << ")" << endl; } 
          }
          if(mode==1){  // exhausted checks
            message << "Please contact David Hicks (david.hicks@duke.edu) and provide the corresponding example." << endl;
            throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_);
          }
        }
        else{ break; }
      }
      else{ break; }
    }

    if(VERBOSE){ for(uint i=0;i<final_permutations.size();i++){ cerr << "Final permutation groupings: " << final_permutations[i] << endl; } }

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
      xstructure xstr_tmp = structure.structure_representative;
      deque<string> species; 
      for(uint j=0;j<name_order[i].size();j++){species.push_back(name_order[i][j]);}
      xstr_tmp.SetSpecies(species);
      //DX TEST xstr_tmp.species_pp = species; //for vasp5 20190731
      xstr_tmp.species = species; //DX 20190813
      xstr_tmp.SpeciesPutAlphabetic();
      xstr_tmp.species_pp = xstr_tmp.species; //for vasp5 20190731, after ordered
      deque<int> sizes = SYM::arrange_atoms(xstr_tmp.atoms);
      //LDEBUGfor(uint j=0;j<sizes.size();j++){cerr << "sizes[j]: " << sizes[j] << endl;}
      xstr_tmp = pflow::SetNumEachType(xstr_tmp, sizes);
      //if (xstr_out.num_each_type.size() != names.size()){
      //  xstr_out = pflow::SetAllAtomNames(xstr_out, in_names);
      //}
      StructurePrototype tmp;
      tmp.structure_representative = xstr_tmp;
      tmp.structure_representative.ReScale(1.0); //DX 20190715
      tmp.structure_representative_name = aurostd::joinWDelimiter(species,"");
      tmp.structure_representative_generated = true;
      //DX 20190730 - ORIG - tmp.structure_representative_from = "permutation";
      stringstream ss_str; ss_str << "permutation of: " << structure.structure_representative; //DX 20190730
      tmp.structure_representative_from = ss_str.str(); //DX 20190730
      tmp.copyPrototypeInformation(structure);
      tmp.stoichiometry = compare::getStoichiometry(tmp.structure_representative, true); //DX 20190529 - need updated stoich
      // compound name is always alphabetic; order swapping is dictated by stoichiometry
      //tmp.structure_representative_compound = compare::getCompoundName(name_order[0],tmp.stoichiometry,false); //remove ones is true //DX 20190529 - need compound name to group later;
      tmp.structure_representative_compound = compare::getCompoundName(tmp.structure_representative); //remove ones is true  //DX 20190311 //DX 20190313 - use xstr1
      tmp.environments_LFA=compare::computeLFAEnvironment(tmp.structure_representative); //DX 20190711
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
// arePermutationsComparableViaStoichiometry() 
// ***************************************************************************
namespace compare{
  bool arePermutationsComparableViaStoichiometry(const xstructure& xstr){ 
    
    // check if permutations of a structure are possible via stoichiometry
    // i.e., check if stoichiometric ratio value occurs more than once
    // xstructure version

    bool matchable_sites=false;
    for(uint i=0;i<xstr.num_each_type.size();i++){
      for(uint j=i+1;j<xstr.num_each_type.size();j++){
        if(xstr.num_each_type[i]==xstr.num_each_type[j]){ matchable_sites=true; break; }
      }
      if(matchable_sites){ break; }
    }
    return matchable_sites;
  }
}

// ***************************************************************************
// arePermutationsComparableViaStoichiometry() 
// ***************************************************************************
namespace compare{
  bool arePermutationsComparableViaStoichiometry(vector<uint>& stoichiometry, bool reduce_stoichiometry){ 
    
    // check if permutations of a structure are possible via stoichiometry
    // i.e., check if stoichiometric ratio value occurs more than once
    // vector<uint> stoichiometry version

    // ---------------------------------------------------------------------------
    // reduce stoichiometry first if necessary 
    vector<uint> tmp_stoich;
    if(reduce_stoichiometry){ tmp_stoich = compare::gcdStoich(stoichiometry); }
    else{ tmp_stoich = stoichiometry; } // assuming it is already reduced

    // ---------------------------------------------------------------------------
    // check if stoichiometries occur more than once
    bool matchable_sites=false;
    for(uint i=0;i<tmp_stoich.size();i++){
      for(uint j=i+1;j<tmp_stoich.size();j++){
        if(tmp_stoich[i]==tmp_stoich[j]){ matchable_sites=true; break; }
      }
      if(matchable_sites){ break; }
    }
    return matchable_sites;
  }
}

// ***************************************************************************
// arePermutationsComparableViaSymmetry() 
// ***************************************************************************
namespace compare{
  bool arePermutationsComparableViaSymmetry(vector<GroupedWyckoffPosition>& grouped_Wyckoff_positions){
    
    // check if permutations of a structure are possible via symmetry (Wyckoff positions)
    // i.e., check if there are matchable Wyckoff positions in a given structure
    
    // ---------------------------------------------------------------------------
    // re-write site symmetries to account for cell choice differences 
    // note: the resulting site symmetries may NOT be physical (this is just a 
    // speed up; fast way to compare)
    vector<GroupedWyckoffPosition> sorted_temp_grouped_Wyckoffs = compare::sortSiteSymmetryOfGroupedWyckoffPositions(grouped_Wyckoff_positions);

    // ---------------------------------------------------------------------------
    // check if Wyckoff site symmetries and multiplicites occur more than once
    for(uint i=0;i<sorted_temp_grouped_Wyckoffs.size();i++){
      string set_1_site_symmetry = aurostd::joinWDelimiter(sorted_temp_grouped_Wyckoffs[i].site_symmetries,",");
      string set_1_multiplicity = aurostd::joinWDelimiter(sorted_temp_grouped_Wyckoffs[i].multiplicities,",");
      for(uint j=i+1;j<sorted_temp_grouped_Wyckoffs.size();j++){
        string set_2_site_symmetry = aurostd::joinWDelimiter(sorted_temp_grouped_Wyckoffs[j].site_symmetries,",");
        string set_2_multiplicity = aurostd::joinWDelimiter(sorted_temp_grouped_Wyckoffs[j].multiplicities,",");
        if(set_1_site_symmetry == set_2_site_symmetry && set_1_multiplicity == set_2_multiplicity){
          return true;
        }
      }
    }
    return false;
  }
}
  
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
      xstructure xstr_tmp = structure.structure_representative;
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
      tmp.structure_representative = xstr_tmp;
      tmp.structure_representative.ReScale(1.0); //DX 20190715
      tmp.structure_representative_name = aurostd::joinWDelimiter(species,"");
      tmp.environments_LFA=compare::computeLFAEnvironment(tmp.structure_representative); //DX 20190711
      tmp.structure_representative_generated = true;
      tmp.structure_representative_from = "permutation";
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
            structure_tmp.structure_representative_name = vlabel[i];
          }
          //if multiple degrees of freedom, then number scheme is required
          else {
            stringstream tmp; tmp << vlabel[i] << "-" << std::setw(3) << std::setfill('0') << j+1;
            structure_tmp.structure_representative_name = tmp.str();
          }
          structure_tmp.stoichiometry = all_structures[0].stoichiometry;
          //structure_tmp.space_group = prototype_space_groups[i]; 
          structure_tmp.space_group = all_structures[0].space_group; // same as representative structure (either will be the same, or we are forcing it to be for the ignore_symmetry/ignore_Wyckoff run)
          structure_tmp.grouped_Wyckoff_positions = all_structures[0].grouped_Wyckoff_positions;
          structure_tmp.elements = compare::fakeElements(all_structures[0].stoichiometry.size());
          structure_tmp.structure_representative_compound = compare::getCompoundName(structure_tmp.elements,structure_tmp.stoichiometry,true); //remove ones is true 
          structure_tmp.structure_representative_generated = false; 
          structure_tmp.structure_representative_from = "aflow_prototypes"; 
          all_structures.push_back(structure_tmp);
        }
      }

      // htqc prototypes
      else {
        StructurePrototype structure_tmp;
        structure_tmp.structure_representative_name = vlabel[i];
        structure_tmp.stoichiometry = all_structures[0].stoichiometry;
        //structure_tmp.space_group = prototype_space_groups[i];
        structure_tmp.space_group = all_structures[0].space_group; // same as representative structure (either will be the same, or we are forcing it to be for the ignore_symmetry/ignore_Wyckoff run)
        structure_tmp.grouped_Wyckoff_positions = all_structures[0].grouped_Wyckoff_positions;
        vector<string> elements;
        structure_tmp.elements = compare::fakeElements(all_structures[0].stoichiometry.size());
        structure_tmp.structure_representative_compound = compare::getCompoundName(structure_tmp.elements,structure_tmp.stoichiometry,true); //remove ones is true 
        structure_tmp.structure_representative_generated = false; 
        structure_tmp.structure_representative_from = "aflow_prototypes"; 
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

    string function_name = "compare::gcdStoich():";
    stringstream message;

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
        message << "Error in GCD procedure. Contact David Hicks (david.hicks@duke.edu)";
        throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_);
      }
    }
    return reduced_numbers;
  }
}

//DX 20191108 [OBSOLETE - switching to getThreadDistribution] // ***************************************************************************
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] // prepareSymmetryThreads - 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] // ***************************************************************************
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] namespace compare{
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]   bool prepareSymmetryThreads(vector<xstructure>& vxstrs, uint& num_proc,
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       vector<uint>& start_indices, vector<uint>& end_indices){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     // Split xstructures via indices, i.e., to be used in different threads for 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     // calculating the symmetry (space group and Wyckoff positions)
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     string function_name = "compare::prepareSymmetryThreads()";
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     stringstream message;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     vector<vector<xstructure> > vxstrs_split;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint num_per_thread = vxstrs.size()/num_proc;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint residual = vxstrs.size()%num_proc;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     bool accounted_for_residual=false;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     if(residual!=0){num_per_thread+=1;}
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint count = 0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint thread_count = 0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     vector<xstructure> tmp;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint tmp_start_index=0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     for(uint l=0; l<vxstrs.size(); l++){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       count+=1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       if(count == num_per_thread && thread_count<num_proc-1){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         thread_count+=1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         start_indices.push_back(tmp_start_index);
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         end_indices.push_back(l);
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         tmp_start_index=l+1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         count = 0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       else if(thread_count==num_proc-1 && l==vxstrs.size()-1){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         thread_count+=1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         start_indices.push_back(tmp_start_index);
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         end_indices.push_back(l);
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         tmp_start_index=l+1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         count = 0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       if(!accounted_for_residual && residual!=0 && thread_count==residual){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         accounted_for_residual=true;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         num_per_thread=num_per_thread-1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     //Need the following safety in case the number of threads is greater than the number of structures to test
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint recovered=0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint num_of_threads=0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     if(start_indices.size()>=num_proc){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       num_of_threads=num_proc;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     else if(start_indices.size()<num_proc){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       num_of_threads=start_indices.size();
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     for(uint n=0; n<num_of_threads; n++){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       for(uint i=start_indices[n];i<=end_indices[n];i++){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         recovered+=1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     } 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     if(recovered != vxstrs.size()){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       message << "The splitting of jobs failed...not all were accounted for: " << recovered << " != " << vxstrs.size();
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_); //DX 20190717 - exit to xerror
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     //DEBUG for(uint i=0;i<vxstrs_split.size();i++){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     //DEBUG   cerr << "num of xstrs for thread: " << i << " = " << vxstrs_split[i].size() << endl;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     //DEBUG }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     return true;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]   }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] }


// ***************************************************************************
// splitComparisonIntoThreads()
// ***************************************************************************
namespace compare{
  bool splitComparisonIntoThreads(vector<StructurePrototype>& comparison_schemes, uint& num_proc,
                              vector<std::pair<uint,uint> >& start_indices,
                              vector<std::pair<uint,uint> >& end_indices){

    // ---------------------------------------------------------------------------
    // split comparisons into threads via indices
    string function_name = "compare::splitComparisonIntoThreads()";
    stringstream message;
    bool LDEBUG=(false || XHOST.DEBUG);
    
    bool safety_check=false; // safety check if split incorrectly

    uint number_of_comparisons = 0;
    for(uint i=0;i<comparison_schemes.size();i++){ number_of_comparisons += comparison_schemes[i].numberOfComparisons(); }
    //DEBUG cerr << "# of comparisons: " << number_of_comparisons << endl;

    if(number_of_comparisons==0){
      if(LDEBUG) {
        cerr << function_name << ": Number of comparisons is zero, no need to split into threads." << endl;
      } 
      return true;
    }   
 
	  uint num_per_thread = number_of_comparisons/num_proc;
	  uint residual = number_of_comparisons%num_proc;
	  bool accounted_for_residual=false;
	  if(residual!=0){num_per_thread+=1;}
   
    if(LDEBUG) {
      cerr << function_name << ": Number of comparisons per thread: " << num_per_thread << endl;
    }

    uint tmp =0;

    uint count = 0;
    uint thread_count = 0;
    std::pair<uint,uint> tmp_start, tmp_end;
    std::pair<uint,uint> indices;
    for(uint i=0;i<comparison_schemes.size();i++){
      //DEBUG cerr << "splitting comparison indices: i: " << i << "/" << comparison_schemes.size() << endl;
      for(uint j=0;j<comparison_schemes[i].structures_duplicate_names.size();j++){
        indices.first=i, indices.second=j;
        count+=1;
        tmp+=1;
        if(count == num_per_thread && thread_count<num_proc-1){
          thread_count+=1;
          start_indices.push_back(tmp_start);
          //update tmp_start
          if(j+1>=comparison_schemes[i].structures_duplicate_names.size()-1 && i+1<comparison_schemes.size()-1){
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
        else if(thread_count==num_proc-1 && i==comparison_schemes.size()-1 && j==comparison_schemes[i].structures_duplicate_names.size()-1){
          thread_count+=1;
          start_indices.push_back(tmp_start);
          //update tmp_start
          if(j+1>=comparison_schemes[i].structures_duplicate_names.size()-1 && i+1<comparison_schemes.size()-1){
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
      tmp_end.first=comparison_schemes.size()-1; tmp_end.second=comparison_schemes[tmp_end.first].structures_duplicate_names.size();
      end_indices.push_back(tmp_end);
    }

    // ---------------------------------------------------------------------------
    // check if split correctly
    // in case the number of threads is greater than the number of xstructures to test
    // put inside an if-statement (on 20190715) to save time; the function has been well-tested
    if(safety_check){
      uint recovered=0;
      uint num_of_threads=0;
      if(start_indices.size()>=num_proc){
        num_of_threads=num_proc;
      }
      else if(start_indices.size()<num_proc){
        num_of_threads=start_indices.size();
      }

      for(uint n=0; n<num_of_threads; n++){
        uint i_min=start_indices[n].first; uint i_max=end_indices[n].first;
        uint j_min=0; uint j_max=0;
        for(uint i=0;i<comparison_schemes.size();i++){
          // to loop properly
          if(i==i_min){
            j_min=start_indices[n].second;
            if(i==i_max){j_max=end_indices[n].second;}
            else {j_max=comparison_schemes[i].structures_duplicate_names.size();} //-1 since in loop: j<=j_max
          }
          else if(i==i_max){j_min=0; j_max=end_indices[n].second;}
          else {j_min=0; j_max=comparison_schemes[i].structures_duplicate_names.size();} //-1 since in loop: j<=j_max
          for(uint j=0;j<comparison_schemes[i].structures_duplicate_names.size();j++){
            if(i>=i_min && j>=j_min &&
                i<=i_max && j<j_max){
              recovered+=1;
            }
          }
        }
      } 
      if(recovered != number_of_comparisons){
        message << "The splitting of jobs failed...not all were accounted for: " << recovered << " != " << number_of_comparisons;
        throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_); //DX 20190717 - exit to xerror
      }
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
    //DX 20191107 [switching to getThreadDistribution] - vector<uint> start_indices, end_indices;
    //DX 20191107 [switching to getThreadDistribution] - prepareSymmetryThreads(vxstrs,num_proc,start_indices,end_indices);
    uint number_of_structures = vxstrs.size(); //DX 20191107
    vector<vector<int> > thread_distribution = getThreadDistribution(number_of_structures, num_proc); //DX 20191107 

    // Run threads 
    vector<std::thread*> threads;
    for(uint n=0; n<num_proc; n++){
      //DX 20191107 [switching to getThreadDistribution] - threads.push_back(std::thread(SYM::calculateSpaceGroupsInSetRange,std::ref(vxstrs),std::ref(start_indices[n]),std::ref(end_indices[n])));
      threads.push_back(new std::thread(&SYM::calculateSpaceGroupsInSetRange,std::ref(vxstrs),thread_distribution[n][0],thread_distribution[n][1]));
    }
    // Join threads 
    for(uint t=0;t<num_proc;t++){
      threads[t]->join();
      delete threads[t];
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

//DX 20191108 [OBSOLETE - switching to getThreadDistribution] // ***************************************************************************
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] // splitTaskIntoThreads - 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] // ***************************************************************************
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] namespace compare{
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]   bool splitTaskIntoThreads(uint& number_of_tasks, uint& num_proc,
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       vector<uint>& start_indices, vector<uint>& end_indices){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     // Split number of tasks via indices, e.g., to be used in different threads for 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     // calculating the symmetry (space group and Wyckoff positions)
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     // It is generalized for any type of job splitting
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     // ---------------------------------------------------------------------------
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     // split comparisons into threads via indices
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     string function_name = "compare::splitTaskIntoThread()";
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     stringstream message;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     bool safety_check=false; // safety check if split incorrectly
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint num_per_thread = number_of_tasks/num_proc;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint residual = number_of_tasks%num_proc;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     bool accounted_for_residual=false;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     if(residual!=0){num_per_thread+=1;}
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint count = 0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint thread_count = 0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     uint tmp_start_index=0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     for(uint l=0; l<number_of_tasks; l++){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       count+=1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       if(count == num_per_thread && thread_count<num_proc-1){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         thread_count+=1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         start_indices.push_back(tmp_start_index);
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         end_indices.push_back(l);
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         tmp_start_index=l+1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         count = 0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       else if(thread_count==num_proc-1 && l==number_of_tasks-1){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         thread_count+=1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         start_indices.push_back(tmp_start_index);
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         end_indices.push_back(l);
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         tmp_start_index=l+1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         count = 0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       if(!accounted_for_residual && residual!=0 && thread_count==residual){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         accounted_for_residual=true;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         num_per_thread=num_per_thread-1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     // ---------------------------------------------------------------------------
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     // check if split correctly
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     // put inside an if-statement (on 20190715) to save time; the function has been well-tested
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     if(safety_check){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       uint recovered=0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       uint num_of_threads=0;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       if(start_indices.size()>=num_proc){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         num_of_threads=num_proc;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       else if(start_indices.size()<num_proc){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         num_of_threads=start_indices.size();
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       for(uint n=0; n<num_of_threads; n++){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         for(uint i=start_indices[n];i<=end_indices[n];i++){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]           recovered+=1;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       } 
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       if(recovered != number_of_tasks){
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         message << "The splitting of jobs failed...not all were accounted for: " << recovered << " != " << number_of_tasks;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]         throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_);
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]       }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]     return true;
//DX 20191108 [OBSOLETE - switching to getThreadDistribution]   }
//DX 20191108 [OBSOLETE - switching to getThreadDistribution] }

// ***************************************************************************
// calculateSpaceGroupsInSetRange
// ***************************************************************************
namespace compare {
  void calculateSpaceGroups(vector<StructurePrototype>& structures, uint start_index, uint end_index){ //DX 20191108 - removed & from uint
   
    // Calculates the space group and Wyckoff positions for the representative 
    // structure in the StructurePrototype object
    // Mirrors SYM::calculateSpaceGroupsInSetRange(), but is specific for 
    // StructurePrototype objects, as opposed to xstructures

    for(uint i=start_index;i<end_index;i++){ //DX 20191107 - switching convention <= vs <
      structures[i].space_group = structures[i].structure_representative.SpaceGroup_ITC();
      vector<GroupedWyckoffPosition> grouped_Wyckoff_positions;
      groupWyckoffPositions(structures[i].structure_representative, grouped_Wyckoff_positions);
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
    //DX 20191107 [switching to getThreadDistribution] - vector<uint> start_indices, end_indices;
    //DX 20191107 [switching to getThreadDistribution] - splitTaskIntoThreads(number_of_structures,num_threads,start_indices,end_indices); //DX 20190530 - renamed
    vector<vector<int> > thread_distribution = getThreadDistribution(number_of_structures, num_threads); //DX 20191107 

    // Run threads (DX 20191108 thread pointer)
    vector<std::thread*> threads;
    for(uint n=0; n<num_threads; n++){
      //DX 20191107 [switching to getThreadDistribution] - threads.push_back(std::thread(compare::calculateSpaceGroupsInSetRange,std::ref(structures),std::ref(start_indices[n]),std::ref(end_indices[n])));
      threads.push_back(new std::thread(&compare::calculateSpaceGroups,std::ref(structures),thread_distribution[n][0],thread_distribution[n][1])); //DX 20191107 [switching to getThreadDistribution] 
    }
    // Join threads
    for(uint t=0;t<num_threads;t++){
      threads[t]->join();
      delete threads[t];
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
// calculateLFAEnvironmentsInSetRange
// ***************************************************************************
namespace compare {
  void calculateLFAEnvironmentsInSetRange(vector<StructurePrototype>& structures, uint start_index, uint end_index){
   
    // Calculates the LFA environments for a structure and
    // stores it in the StructurePrototype object

    for(uint i=start_index;i<end_index;i++){ //DX 20191107 switching end index convention <= vs <
      structures[i].environments_LFA = compare::computeLFAEnvironment(structures[i].structure_representative);
    }
  }
}

// ***************************************************************************
// calculateLFAEnvironments - Calculate LFA environments
// ***************************************************************************
namespace compare{
  void calculateLFAEnvironments(vector<StructurePrototype>& structures, uint num_proc){

    // Calculates the LFA environments for a structure and
    // stores it in the StructurePrototype object
    
    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name="compare::calculateLFAEnvironments()";
    if(LDEBUG) {cerr << function_name << ": Number of threads=" << num_proc << endl;}

#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
    // THREADED VERISON - START
    
    // Distribute threads via indices
    uint number_of_structures = structures.size();
    uint num_threads = aurostd::min(num_proc,number_of_structures); // cannot have more threads than structures
    //DX 20191107 [switching to getThreadDistribution] - vector<uint> start_indices, end_indices;
    //DX 20191107 [switching to getThreadDistribution] - splitTaskIntoThreads(number_of_structures,num_threads,start_indices,end_indices); //DX 20190530 - renamed
    vector<vector<int> > thread_distribution = getThreadDistribution(number_of_structures, num_threads); //DX 20191107 

    // Run threads 
    vector<std::thread*> threads;
    for(uint n=0; n<num_threads; n++){
      //DX 20191107 [switching to getThreadDistribution] - threads.push_back(std::thread(compare::calculateLFAEnvironmentsInSetRange,std::ref(structures),start_indices[n],end_indices[n]));
      threads.push_back(new std::thread(&compare::calculateLFAEnvironmentsInSetRange,std::ref(structures),thread_distribution[n][0],thread_distribution[n][1])); //DX 20191107 [switching to getThreadDistribution] -
    }
    // Join threads
    for(uint t=0;t<num_threads;t++){
      threads[t]->join();
      delete threads[t];
    }
    // THREADED VERISON - END

#else
    // NON-THREADED VERSION - START
    for(uint i=0; i<structures.size(); i++){
      structures[i].environments_LFA = compare::computeLFAEnvironment(structures[i].structure_representative);
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
  vector<GroupedWyckoffPosition> sortSiteSymmetryOfGroupedWyckoffPositions(const vector<GroupedWyckoffPosition>& grouped_Wyckoffs){

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
  bool matchableWyckoffPositions(const vector<GroupedWyckoffPosition>& temp_grouped_Wyckoffs,
                                 const vector<GroupedWyckoffPosition>& representative_grouped_Wyckoffs, 
                                 bool same_species){

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
          tmp.structure_representative_name = vfiles[i];
        }
        else {
          tmp.structure_representative_name = directory+"/"+vfiles[i];  aurostd::StringSubst(tmp.structure_representative_name,"//","/"); //DX 20181003
        }
        tmp.structure_representative_generated=vstructures_generated[i];
        tmp.structure_representative_from=vstructures_from[i];
        //tmp.number_of_types=vxstrs[i].num_each_type.size();
        tmp.elements=vvelements[i];
        tmp.stoichiometry=vstoichs[i];
        //tmp.number_of_atoms=vxstrs[i].atoms.size();
        tmp.Pearson=vpearsons[i];
        tmp.space_group=vsgroups[i];
        tmp.grouped_Wyckoff_positions=vgrouped_Wyckoff_positions[i];
        if(property_names.size()!=0){
          tmp.property_names=property_names; //DX 20181218 - added property_names
          tmp.property_units=property_units; //DX 20181218 - added property_units
          tmp.properties_structure_representative=property_values[i]; //DX 20181218 - added property_values
        }
        if(vstructures_generated[i]){
          tmp.structure_representative=vxstrs[i];
          tmp.number_of_types=vxstrs[i].num_each_type.size();
          tmp.number_of_atoms=vxstrs[i].atoms.size();
          tmp.structure_representative_compound=getCompoundName(vxstrs[i]); //DX 20190111 - added compound, e.g., Ag1Br2
        }
        comparison_schemes.push_back(tmp);
      }
      else {
        for(uint j=0; j<comparison_schemes.size(); j++){
          bool same_material_stoich=false;
          ostringstream tmp;
          tmp.clear();
          if(same_species==true && 
             matchableSpecies(vxstrs[i],comparison_schemes[j].structure_representative,same_species)==true){
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
            comparison_schemes[j].structures_duplicate_names.push_back(duplicate_name);
            comparison_schemes[j].structures_duplicate_generated.push_back(vstructures_generated[i]);
            comparison_schemes[j].structures_duplicate_from.push_back(vstructures_from[i]);
            //cerr << "adding to " << j << " (name): " << comparison_schemes[j].structures_duplicate_names.size() << endl;
            //cerr << "adding to " << j << " (gen): " << comparison_schemes[j].structures_duplicate_generated.size() << endl;
            //cerr << "adding to " << j << " (from): " << comparison_schemes[j].structures_duplicate_from.size() << endl;
            if(vstructures_generated[i]){
              comparison_schemes[j].structures_duplicate.push_back(vxstrs[i]);
              comparison_schemes[j].structures_duplicate_compounds.push_back(getCompoundName(vxstrs[i])); //DX 20190111 - added compound, e.g., Ag1Br2
            }
            comparison_schemes[j].misfits_duplicate.push_back(-1.0);
            if(property_names.size()!=0){
              comparison_schemes[j].properties_structures_duplicate.push_back(property_values[i]); //DX 20181218 - added property_values
            }
            scheme_created=true;
            break;
          }
          //cerr << "!!!!!!!!!!!!!!!!!!!!UNMATCHABLE: " << comparison_schemes[j].structure_representative_name << " and " << directory+"/"+vfiles[i] <<  endl;
        }
        if(scheme_created==false){
          StructurePrototype tmp;
          if(directory==""){
            tmp.structure_representative_name = vfiles[i];
          }
          else {
            tmp.structure_representative_name = directory+"/"+vfiles[i];  aurostd::StringSubst(tmp.structure_representative_name,"//","/"); //DX 20181003
          }
          tmp.structure_representative_generated=vstructures_generated[i];
          tmp.structure_representative_from=vstructures_from[i];
          //tmp.number_of_types=vxstrs[i].num_each_type.size();
          tmp.elements=vvelements[i];
          tmp.stoichiometry=vstoichs[i];
          //tmp.number_of_atoms=vxstrs[i].atoms.size();
          tmp.Pearson=vpearsons[i];
          tmp.space_group=vsgroups[i];
          tmp.grouped_Wyckoff_positions=vgrouped_Wyckoff_positions[i];
          if(property_names.size()!=0){
            tmp.property_names=property_names; //DX 20181218 - added property_names
            tmp.property_units=property_units; //DX 20181218 - added property_units
            tmp.properties_structure_representative=property_values[i]; //DX 20181218 - added property_values
          }
          if(vstructures_generated[i]){
            tmp.structure_representative=vxstrs[i];
            tmp.number_of_types=vxstrs[i].num_each_type.size();
            tmp.number_of_atoms=vxstrs[i].atoms.size();
            tmp.structure_representative_compound=getCompoundName(vxstrs[i]); //DX 20190111 - added compound, e.g., Ag1Br2
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
    //      cerr << i << "structures_duplicate.size(): " << comparison_schemes[i].structures_duplicate.size() << endl;
    //      cerr << i << "structures_duplicate_generated.size(): " << comparison_schemes[i].structures_duplicate_generated.size() << endl;
    //      cerr << i << "structures_duplicate_from.size(): " << comparison_schemes[i].structures_duplicate_from.size() << endl;
    //}
  }
}
        
// ***************************************************************************
// structuresCompatible - check compatiblity of structures 
// ***************************************************************************
namespace compare{
  bool structuresCompatible(const StructurePrototype& structure1,
    const StructurePrototype& structure2, bool same_species,  
    bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment,
    bool duplicates_removed){ //DX 20190829 - added duplicates_removed

    // ---------------------------------------------------------------------------
    // check if species/stoichiometries are compatible
    //DX 20190430 - this may take longer, use compound if(same_species==true && matchableSpecies(structures[i].structure_representative,comparison_schemes[j].structure_representative,same_species)==true){
    if(same_species==true && structure1.structure_representative_compound!=structure2.structure_representative_compound){ //DX 20190430 - quicker //DX 20190702 - changed to "!=" and "false" for speed increase
      return false;
    }
    else if(same_species==false && structure1.stoichiometry!=structure2.stoichiometry){ //DX 20190702 - changed to "!=" and "false" for speed increase
      return false;
    }
    // if already removed duplicate compounds, then structures were already compared, so don't compare again
    else if(same_species==false && duplicates_removed && structure1.structure_representative_compound==structure2.structure_representative_compound){
      return false;
    }
    // ---------------------------------------------------------------------------
    // check if LFA environments are compatible - DX 20190711
    if(!ignore_environment && !compatibleEnvironmentSets(structure1.environments_LFA,structure2.environments_LFA,same_species,false)){
      return false;
    }
    // ---------------------------------------------------------------------------
    // check symmetry (if applicable) 
    //DX 20190702 - checking stoich is redundant for compound checking - if(same_material_stoich==true && structures[i].stoichiometry==comparison_schemes[j].stoichiometry && 
    //DX 20190702 [OBSOLETE] if(same_material_stoich==true &&  //DX 20190702 - moved stoichiometry up
    if(((ignore_symmetry && ignore_Wyckoff) ||    
        (!ignore_symmetry && ignore_Wyckoff &&
         structure1.Pearson == structure2.Pearson && 
         matchableSpaceGroups(structure1.space_group,structure2.space_group)) || 
        (!ignore_symmetry && !ignore_Wyckoff &&
         structure1.Pearson == structure2.Pearson && 
         matchableSpaceGroups(structure1.space_group,structure2.space_group) &&
         matchableWyckoffPositions(structure1.grouped_Wyckoff_positions,structure2.grouped_Wyckoff_positions,same_species)))){
       return true;
    }
    return false;
  }
}

// ***************************************************************************
// createStructurePrototypes - Group structures by Pearson symbol, then space group
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> groupStructurePrototypes(vector<StructurePrototype>& structures, 
				 bool same_species, bool ignore_symmetry, bool ignore_Wyckoff, bool ignore_environment,
         bool duplicates_removed){ //DX 20190829 - added duplicates_removed

    // Populates the structure information into the StructurePrototype object.
    // It groups structure based on their stoichiometry, space group, and Wyckoff positions. 
    // A "representative" structure is chosen and will be compared to the 
    // possible "duplicates". The misfit values are set to -1.0 until compared.
    
    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "compare::groupStructurePrototypes()";
   
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
    //   ignore_environment : ignores LFA environment analysis 

    for(uint i=0;i<structures.size(); i++){
      bool scheme_created=false;
      for(uint j=0; j<comparison_schemes.size(); j++){
        //bool same_material_stoich=false;
        ostringstream tmp;
        tmp.clear();

        if(structuresCompatible(structures[i], comparison_schemes[j], same_species, ignore_symmetry, ignore_Wyckoff, ignore_environment, duplicates_removed)){ //DX 20190829 - added duplicates_removed
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
      cerr << function_name << ": Prepared comparison sets: " << endl;
      stringstream ss_test;
      compare::printResults(ss_test, same_species, comparison_schemes);
      cerr << ss_test.str() << endl;
    }
    // DEBUG for(uint i=0;i<comparison_schemes.size();i++){
    // DEBUG  cerr << i << "structures_duplicate.size(): " << comparison_schemes[i].structures_duplicate.size() << endl;
    // DEBUG  cerr << i << "structures_duplicate_names.size(): " << comparison_schemes[i].structures_duplicate_names.size() << endl;
    // DEBUG  cerr << i << "structures_duplicate_generated.size(): " << comparison_schemes[i].structures_duplicate_generated.size() << endl;
    // DEBUG  cerr << i << "structures_duplicate_from.size(): " << comparison_schemes[i].structures_duplicate_from.size() << endl;
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
// checkForBetterMatches 
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> checkForBetterMatches(vector<StructurePrototype>& prototype_schemes, 
    ostream& oss, ofstream& FileMESSAGE, uint& num_proc, bool check_for_better_matches, bool same_species,
    bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, 
    bool ignore_environment, bool clean_unmatched, bool ICSD_comparison, bool quiet){ 

    // this function checks if other groups based on 


    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "compare::checkForBetterMatches()";
    stringstream message;
    ostream& logstream = cout;

    double misfit_min = 0.01; // used for check_for_better_matches : this is quite strict; if too expensive, make more loose
    double misfit_max = 0.1; // used for check_for_better_matches : otherwise we will compare same family structures which have already been moved (if using !clean_unmatched)
    bool single_comparison_round = true; // always true for this function
    bool store_comparison_logs = false; //DX 20190822 - add log bool
    bool duplicates_removed = false; //DX 20190829 - should have no affect in this function

    // ---------------------------------------------------------------------------
    // check which structures could potentially match to others based on misfit
    vector<StructurePrototype> comparison_groups;
  
    for(uint i=0;i<prototype_schemes.size();i++){
      for(uint j=0;j<prototype_schemes[i].structures_duplicate_names.size();j++){
        if((check_for_better_matches && prototype_schemes[i].misfits_duplicate[j] > misfit_min && prototype_schemes[i].misfits_duplicate[j] < misfit_max) || // find better match
           (!check_for_better_matches && (prototype_schemes[i].misfits_duplicate[j] > 0.1 || aurostd::isequal(prototype_schemes[i].misfits_duplicate[j],1.0,1e-6) || aurostd::isequal(prototype_schemes[i].misfits_duplicate[j],-1.0,1e-6)))){   // find a match
          StructurePrototype tmp;
          bool found_new_match=false;
          // ---------------------------------------------------------------------------
          // check for other compatible representative structures 
          // start_index=i+1 : (only need to search forward for better matches, due to appendStructurePrototypes() scheme)
          uint start_index = 0;
          if(check_for_better_matches){ start_index = i+1; }
          for(uint k=start_index;k<prototype_schemes.size();k++){
            if(structuresCompatible(prototype_schemes[i], prototype_schemes[k], same_species, ignore_symmetry, ignore_Wyckoff, ignore_environment, duplicates_removed)){ // can check based on representatives; duplicate info matches its representative info //DX 20190829 - added duplicates_removed
              if(!quiet || LDEBUG){
                message << "Found potential match for " << prototype_schemes[i].structures_duplicate_names[j] << ": " << prototype_schemes[k].structure_representative_name; 
                pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
              }

              // ---------------------------------------------------------------------------
              // reverse the convention, single structurePrototype object to find best match
              // i.e., duplicate -> representative and representatives -> duplicates 
              if(!found_new_match){
                tmp.copyPrototypeInformation(prototype_schemes[i]);
                tmp.putDuplicateAsRepresentative(prototype_schemes[i],j);
                // store the current match so we can check fast if it matches to any other
                tmp.addStructurePrototypeAsDuplicate(prototype_schemes[i]); // store the current match structure
                //tmp.misfits_duplicate.push_back(prototype_schemes[i].misfits_duplicate[j]); // store the current match misfit 
                found_new_match=true;
              }
              if(k!=j){
              tmp.addStructurePrototypeAsDuplicate(prototype_schemes[k]);
              }
            }
          }
          if(found_new_match){
            comparison_groups.push_back(tmp);
          }
        }
      }
    }

    // ---------------------------------------------------------------------------
    // compare structures 
    vector<StructurePrototype> other_matches_schemes = compare::runComparisonScheme(num_proc, comparison_groups, same_species, duplicates_removed, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, false, ICSD_comparison, store_comparison_logs, oss, FileMESSAGE, quiet);  //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched //DX 20190731 - added ignore_symmetry/Wyckoff/environment //DX 20190822 - add log bool

    // ---------------------------------------------------------------------------
    // check if there are any better matches and reorganize if necessary
    // the original match is stored in the first position
    for(uint i=0;i<other_matches_schemes.size();i++){
      double min_misfit = aurostd::abs(other_matches_schemes[i].misfits_duplicate[0]); // put first as min //abs to turn -1 into 1 for comparison
      uint min_index = 0;
      for(uint j=1;j<other_matches_schemes[i].misfits_duplicate.size();j++){
        if(other_matches_schemes[i].misfits_duplicate[j]<min_misfit && aurostd::isdifferent(other_matches_schemes[i].misfits_duplicate[j],-1.0,1e-6)){
          min_misfit=other_matches_schemes[i].misfits_duplicate[j];
          min_index=j;
        }
      }
      // ---------------------------------------------------------------------------
      // move if found better match
      if(min_index!=0){
        for(uint j=0;j<prototype_schemes.size();j++){
          // add structure to its better matching representative
          if(prototype_schemes[j].structure_representative_name == other_matches_schemes[i].structures_duplicate_names[min_index]){
            if(!quiet || LDEBUG){
              message << other_matches_schemes[i].structure_representative_name << " matches better with " << prototype_schemes[j].structure_representative_name; 
              pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
            }
            prototype_schemes[j].addStructurePrototypeAsDuplicate(other_matches_schemes[i]);
            prototype_schemes[j].misfits_duplicate.back()=other_matches_schemes[i].misfits_duplicate[min_index];
          }
          // remove from old representative
          if(clean_unmatched && prototype_schemes[j].structure_representative_name == other_matches_schemes[i].structures_duplicate_names[0]){
            for(uint k=0;k<prototype_schemes[j].structures_duplicate_names.size();k++){
              if(prototype_schemes[j].structures_duplicate_names[k] == other_matches_schemes[i].structure_representative_name){
                if(!quiet || LDEBUG){
                  message << "removing " << other_matches_schemes[i].structure_representative_name << " from " << prototype_schemes[j].structure_representative_name << " set"; 
                  pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
                }
                prototype_schemes[j].removeNonDuplicate(k);
                break;
              }
            }
          }
        }
      }
      else{
        if(!quiet){
          message << other_matches_schemes[i].structure_representative_name << " matches better with original set " << other_matches_schemes[i].structures_duplicate_names[0];
          pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
        }
      }
    }
    
    if(LDEBUG){
      for(uint i=0;i<prototype_schemes.size();i++){ cerr << function_name << " prototype_schemes[i]: " << prototype_schemes[i] << endl; }
    }
    return prototype_schemes;
  }
}

// ***************************************************************************
// compareDuplicateCompounds()
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
      if(prototype_schemes[i].structures_duplicate_names.size()>0){ //if 0, none to check; 
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
    bool ignore_symmetry=false;
    bool ignore_Wyckoff=false;
    bool ignore_environment=false;
    bool optimize_match=false;
    bool single_comparison_round=false;
    bool clean_unmatched=true; //DX 20190504
    bool store_comparison_logs=false; //DX 20190822
    bool check_other_groupings=false; //DX 20190830
    //bool structures_generated=false;

    message << "Running comparisons to remove duplicate compounds ...";
    pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
    vector<StructurePrototype> final_prototypes_reduced = compare::runComparisonScheme(num_proc, duplicate_check_schemes, same_species, check_other_groupings, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, ignore_environment, single_comparison_round, clean_unmatched, ICSD_comparison, store_comparison_logs, oss, FileMESSAGE); //DX 20190319 - added FileMESSAGE //DX 20190731 - //DX 20190731 - added ignore_symmetry/Wyckoff/environment //DX 20190822 - add log bool

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
    for(uint i=0;i<prototype_scheme.structures_duplicate_names.size();i++){
      bool scheme_created=false;
      //cerr << i << endl;
      if(i==0){
        StructurePrototype tmp = prototype_scheme;
        tmp.structure_representative_name = prototype_scheme.structures_duplicate_names[i];
        tmp.structure_representative_compound = prototype_scheme.structures_duplicate_compounds[i];
        if(prototype_scheme.properties_structures_duplicate.size()>0){
          tmp.properties_structure_representative = prototype_scheme.properties_structures_duplicate[i];
        }
        else {
          tmp.properties_structure_representative.clear();
        }
        tmp.structure_representative.Clear(); //clear xstructure;
        tmp.structure_representative_generated=false;
        tmp.structure_representative_from=prototype_scheme.structures_duplicate_from[i];
        tmp.structures_duplicate_names.clear(); tmp.structures_duplicate_compounds.clear(); tmp.structures_duplicate.clear(); tmp.misfits_duplicate.clear(); tmp.properties_structures_duplicate.clear();
        tmp.structures_duplicate_generated.clear(); tmp.structures_duplicate_from.clear(); tmp.structures_duplicate_grouped_Wyckoff_positions.clear(); //DX 20190814 - added Wyckoff positions
        tmp.structures_family_names.clear(); tmp.structures_family.clear(); tmp.misfits_family.clear(); tmp.properties_structures_family.clear(); //DX 20190425 - added properties
        tmp.structures_family_generated.clear(); tmp.structures_family_from.clear();
        tmp.elements.clear();
        duplicate_check_schemes.push_back(tmp);
      }
      else {
        for(uint j=0; j<duplicate_check_schemes.size(); j++){
          if(prototype_scheme.structures_duplicate_compounds[i] == duplicate_check_schemes[j].structure_representative_compound){
            scheme_created=true;
            duplicate_check_schemes[j].structures_duplicate_names.push_back(prototype_scheme.structures_duplicate_names[i]);
            duplicate_check_schemes[j].structures_duplicate_compounds.push_back(prototype_scheme.structures_duplicate_compounds[i]);
            xstructure tmp_xstr;
            duplicate_check_schemes[j].structures_duplicate.push_back(tmp_xstr);
            duplicate_check_schemes[j].structures_duplicate_generated.push_back(false);
            duplicate_check_schemes[j].structures_duplicate_from.push_back(prototype_scheme.structures_duplicate_from[i]);
            duplicate_check_schemes[j].structures_duplicate_grouped_Wyckoff_positions.push_back(prototype_scheme.structures_duplicate_grouped_Wyckoff_positions[i]);
            if(prototype_scheme.properties_structures_duplicate.size()>0){
              duplicate_check_schemes[j].properties_structures_duplicate.push_back(prototype_scheme.properties_structures_duplicate[i]);
            }
            else {
              duplicate_check_schemes[j].properties_structures_duplicate.clear();
            }
            duplicate_check_schemes[j].misfits_duplicate.push_back(-1.0);
            break;
          }
        }
        if(!scheme_created){
          StructurePrototype tmp = prototype_scheme;
          tmp.structure_representative_name = prototype_scheme.structures_duplicate_names[i];
          tmp.structure_representative_compound = prototype_scheme.structures_duplicate_compounds[i];
          if(prototype_scheme.properties_structures_duplicate.size()>0){
            tmp.properties_structure_representative = prototype_scheme.properties_structures_duplicate[i];
          }
          else {
            tmp.properties_structure_representative.clear();
          }
          tmp.structure_representative.Clear(); //clear xstructure;
          tmp.structure_representative_generated=false;
          tmp.structure_representative_from=prototype_scheme.structures_duplicate_from[i];
          tmp.structures_duplicate_names.clear(); tmp.structures_duplicate_compounds.clear(); tmp.structures_duplicate.clear(); tmp.misfits_duplicate.clear(); tmp.properties_structures_duplicate.clear();
          tmp.structures_duplicate_generated.clear(); tmp.structures_duplicate_from.clear(); tmp.structures_duplicate_grouped_Wyckoff_positions.clear(); //DX 20190814 - added Wyckoff positions
          tmp.structures_family_names.clear(); tmp.structures_family.clear(); tmp.misfits_family.clear(); tmp.properties_structures_family.clear(); //DX 20190425 - added properties
          tmp.structures_family_generated.clear(); tmp.structures_family_from.clear();
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
      for(uint j=0;j<final_prototypes[i].structures_duplicate_names.size();j++){
        bool remove_structure = true;
        for(uint k=0;k<duplicate_compound_comparisons.size();k++){
          if(final_prototypes[i].structures_duplicate_names[j] == duplicate_compound_comparisons[k].structure_representative_name){
            remove_structure = false;
            break;
          }
        }
        if(remove_structure){
          final_prototypes[i].structures_duplicate_names.erase(final_prototypes[i].structures_duplicate_names.begin()+j);
          final_prototypes[i].structures_duplicate_compounds.erase(final_prototypes[i].structures_duplicate_compounds.begin()+j);
          final_prototypes[i].structures_duplicate_generated.erase(final_prototypes[i].structures_duplicate_generated.begin()+j);
          final_prototypes[i].structures_duplicate_from.erase(final_prototypes[i].structures_duplicate_from.begin()+j);
          if(final_prototypes[i].property_names.size()!=0){
            final_prototypes[i].properties_structures_duplicate.erase(final_prototypes[i].properties_structures_duplicate.begin()+j);
          }
          final_prototypes[i].misfits_duplicate.erase(final_prototypes[i].misfits_duplicate.begin()+j);
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
      if(comparison_schemes[i].structures_duplicate_names.size()){
        vector<string> ICSD_entries;
        ICSD_entries.push_back(findICSDName(comparison_schemes[i].structure_representative_name));
        for(uint j=0;j<comparison_schemes[i].structures_duplicate_names.size();j++){
          ICSD_entries.push_back(findICSDName(comparison_schemes[i].structures_duplicate_names[j]));
        }
        string min_ICSD_entry = findMinimumICSDEntry(ICSD_entries);
        if(!aurostd::substring2bool(comparison_schemes[i].structure_representative_name,min_ICSD_entry) && !min_ICSD_entry.empty()){ //DX 20191108 - add not empty case
          for(uint j=0;j<comparison_schemes[i].structures_duplicate_names.size();j++){
            if(aurostd::substring2bool(comparison_schemes[i].structures_duplicate_names[j],min_ICSD_entry)){
              string old_representative_ID = comparison_schemes[i].structure_representative_name;
              string old_representative_compound = comparison_schemes[i].structure_representative_compound;
              bool old_representative_generated = comparison_schemes[i].structure_representative_generated;
              string old_representative_from = comparison_schemes[i].structure_representative_from;
              vector<GroupedWyckoffPosition> old_Wyckoff_positions = comparison_schemes[i].grouped_Wyckoff_positions; //DX 20190813 - need to update; otherwise, the compound name and Wyckoff positions may not match; especially for structure-type comparisons
              uint old_representative_duplicate_count = comparison_schemes[i].number_compounds_matching_representative;
              vector<string> old_representative_properties = comparison_schemes[i].properties_structure_representative;
              comparison_schemes[i].structure_representative_name = comparison_schemes[i].structures_duplicate_names[j];  
              comparison_schemes[i].structure_representative_compound = comparison_schemes[i].structures_duplicate_compounds[j];  
              comparison_schemes[i].structure_representative_generated = comparison_schemes[i].structures_duplicate_generated[j];  
              comparison_schemes[i].structure_representative_from = comparison_schemes[i].structures_duplicate_from[j]; 
              comparison_schemes[i].grouped_Wyckoff_positions = comparison_schemes[i].structures_duplicate_grouped_Wyckoff_positions[j]; //DX 20190813 - need to update; otherwise, the compound name and Wyckoff positions may not match; especially for structure-type comparisons
              comparison_schemes[i].number_compounds_matching_representative = comparison_schemes[i].number_compounds_matching_duplicate[j]; 
              if(old_representative_generated){
                xstructure old_representative_xstr = comparison_schemes[i].structure_representative;
                if(comparison_schemes[i].structures_duplicate_generated[j]){
                  comparison_schemes[i].structure_representative = comparison_schemes[i].structures_duplicate[j]; 
                  comparison_schemes[i].structures_duplicate[j] = old_representative_xstr;
                  comparison_schemes[i].structures_duplicate_generated[j] = old_representative_generated;
                }
                else {
                  comparison_schemes[i].structure_representative.Clear();
                  comparison_schemes[i].structures_duplicate_generated[j]=false; //cannot guarantee the rest of the vector is generated; may populate wrong index
                 //DX 20190304 - should I also clear out vector<xvstructure>?
                }
              }
              else {
                comparison_schemes[i].structures_duplicate_generated[j] = old_representative_generated;
              }
              if(comparison_schemes[i].properties_structures_duplicate.size()>0){
                comparison_schemes[i].properties_structure_representative = comparison_schemes[i].properties_structures_duplicate[j];
              }
              else {
                comparison_schemes[i].properties_structure_representative.clear();
              }
              comparison_schemes[i].structures_duplicate_names[j] = old_representative_ID;
              //DX 20190304 - moved into if statment up above - comparison_schemes[i].structures_duplicate[j] = old_representative_xstr;
              comparison_schemes[i].structures_duplicate_compounds[j] = old_representative_compound;
              //DX 20190304 - moved into if statment up above - comparison_schemes[i].structures_duplicate_generated[j] = old_representative_generated;
              comparison_schemes[i].structures_duplicate_from[j] = old_representative_from;
              comparison_schemes[i].structures_duplicate_grouped_Wyckoff_positions[j] = old_Wyckoff_positions; //DX 20190813
              comparison_schemes[i].number_compounds_matching_duplicate[j] = old_representative_duplicate_count;
              if(comparison_schemes[i].properties_structures_duplicate.size()>0){
                comparison_schemes[i].properties_structures_duplicate[j] = old_representative_properties;
              }
              else {
                comparison_schemes[i].properties_structures_duplicate.clear();
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
      ostream& oss,
      bool same_species, 
      bool scale_volume, bool optimize_match, 
      bool store_comparison_logs){ 
 
    // Run comparison thread
    // If the xstructure is not generated, it will generate a local copy for 
    // the single comparison only (prevents overwriting in the comparisons) 
 
    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "compare::runComparisonThreads()";
    stringstream message;
    //bool store_comparison_logs = false; //DX 20190624
   
    //// check if only one comparison, then we store the comparison logs //DX 20190802
    //uint number_of_comparisons = 0;
    //for(uint i=0;i<comparison_schemes.size();i++){ number_of_comparisons += comparison_schemes[i].numberOfComparisons(); }
    //if(number_of_comparisons==1){store_comparison_logs=true;}

    uint i_min=start_indices.first; uint i_max=end_indices.first;
    uint j_min=0; uint j_max=0;

    for(uint i=i_min; i<=i_max; i++){
      xstructure structure_representative;

      // to loop properly
      if(i==i_min){
        j_min=start_indices.second;
        if(i==i_max){j_max=end_indices.second;}
        else {j_max=comparison_schemes[i].structures_duplicate_names.size();} //-1 since in loop: j<=j_max
      }
      else if(i==i_max){j_min=0; j_max=end_indices.second;}
      else {j_min=0; j_max=comparison_schemes[i].structures_duplicate_names.size();} //-1 since in loop: j<=j_max
     
      for(uint j=j_min; j<j_max; j++){

        // get representative structure 
        if(!comparison_schemes[i].structure_representative_generated){
          if(!generateStructure(comparison_schemes[i].structure_representative_name,comparison_schemes[i].structure_representative_from,structure_representative,oss)){
            message << "Could not generate representative structure (" << comparison_schemes[i].structure_representative_name << ").";
            throw aurostd::xerror(function_name,message,_INPUT_ERROR_); //DX 20190717 - exit to xerror
          }
        }
        else {
          structure_representative = comparison_schemes[i].structure_representative;
        }
        //if(LDEBUG) { cerr << function_name << ": Loaded representative structure = " << comparison_schemes[i].structure_representative_name << endl; }
      
        // get prototype structure 
        xstructure duplicate_structure;
        if(!comparison_schemes[i].structures_duplicate_generated[j]){
          if(!generateStructure(comparison_schemes[i].structures_duplicate_names[j],comparison_schemes[i].structures_duplicate_from[j],duplicate_structure,oss)){
            message << "Could not generate duplicate structure (" << comparison_schemes[i].structures_duplicate_names[j] << ").";
            throw aurostd::xerror(function_name,message,_INPUT_ERROR_); //DX 20190717 - exit to xerror
          }
        }
        else {
          duplicate_structure = comparison_schemes[i].structures_duplicate[j];
        }
        if(LDEBUG) { cerr << function_name << ": Loaded duplicate structure = " << comparison_schemes[i].structures_duplicate_names[j] << endl; }
        
        // call the main comparison function
        ostringstream tmp_oss; tmp_oss.clear();
        double final_misfit=-1.0;
        if(LDEBUG) { cerr << function_name << ": Comparing " << comparison_schemes[i].structure_representative_name << " and " << comparison_schemes[i].structures_duplicate_names[j] <<  endl; }
        compare::aflowCompareStructure(1, structure_representative, //num_proc -> 1 for threads (not sure how it behaves otherwise)
                                       duplicate_structure,
                                       same_species, scale_volume, optimize_match, tmp_oss, final_misfit);
        
        // store the figure of misfit
        if(LDEBUG) { cerr << function_name << ": Comparison complete, misfit = " << final_misfit << "." << endl; }
        comparison_schemes[i].misfits_duplicate[j]=final_misfit;
        if(store_comparison_logs){comparison_schemes[i].duplicate_comparison_logs.push_back(tmp_oss.str());} //DX 20190506
      }
    }
  }
}

// ***************************************************************************
// runComparisonScheme: Runs comparisons automatically 
// ***************************************************************************
namespace compare{
  vector<StructurePrototype> runComparisonScheme(uint num_proc, vector<StructurePrototype>& comparison_schemes, 
      bool same_species, bool check_other_grouping, bool scale_volume, bool optimize_match, bool ignore_symmetry, bool ignore_Wyckoff, 
      bool ignore_environment, bool single_comparison_round, bool clean_unmatched, bool ICSD_comparison, bool store_comparison_logs, 
      ostream& oss, ofstream& FileMESSAGE, bool quiet){ //DX 20190319 - added FileMESSAGE //DX 20190504 - added clean unmatched //DX 20190731 - removed const and &, added ignore_symmetry/Wyckoff/environment //DX 20190822 - added logs bool //DX 20190829 - added check_other_grouping

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

    // run threads (DX 20191108 - thread pointer)
    vector<std::thread*> threads;
    for(uint n=0; n<num_comparison_threads; n++){
      threads.push_back(new std::thread(&compare::runComparisonThreads,std::ref(comparison_schemes),std::ref(start_indices[n]),std::ref(end_indices[n]),
            std::ref(oss),same_species,scale_volume,optimize_match,store_comparison_logs));
    }        

    // join threads
    for(uint t=0;t<threads.size();t++){
      threads[t]->join();
      delete threads[t];
    }

    // THREADED VERSION - END

#else
    // NON-THREADED VERISON - START
    if(LDEBUG) { cerr << function_name << ": Non-threaded version." << endl; } 

    for(uint i=0; i<comparison_schemes.size(); i++){
      xstructure structure_representative;
      for(uint j=0; j<comparison_schemes[i].structures_duplicate_names.size(); j++){
        if(j==0){
          if(!comparison_schemes[i].structure_representative_generated){
            if(!generateStructure(comparison_schemes[i].structure_representative_name,comparison_schemes[i].structure_representative_from,structure_representative,oss)){
              message << "Could not generate representative structure (" << comparison_schemes[i].structure_representative_name << ").";
              throw aurostd::xerror(function_name,message,_INPUT_ERROR_); //DX 20190717 - exit to xerror
            }
          }
          else {
            structure_representative = comparison_schemes[i].structure_representative;
          }
        }
        ostringstream tmp_oss;
        tmp_oss.clear();
        double final_misfit=-1.0;
        xstructure duplicate_structure;
        if(!comparison_schemes[i].structures_duplicate_generated[j]){
          if(!generateStructure(comparison_schemes[i].structures_duplicate_names[j],comparison_schemes[i].structures_duplicate_from[j],duplicate_structure,oss)){
            message << "Could not generate duplicate structure (" << comparison_schemes[i].structures_duplicate_names[j] << ").";
            throw aurostd::xerror(function_name,message,_INPUT_ERROR_); //DX 20190717 - exit to xerror
          }
        }
        else {
          duplicate_structure = comparison_schemes[i].structures_duplicate[j];
        }
        // Call the main comparison function
        compare::aflowCompareStructure(num_proc, structure_representative,
            duplicate_structure,
            same_species, scale_volume, optimize_match, tmp_oss, final_misfit);
        // Store the figure of misfit
        comparison_schemes[i].misfits_duplicate[j]=final_misfit;
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

    // for this first iteration: LFA may have incorrectly grouped, so we need to check if the structures belong to other groups
    // OR after removing duplicate compounds these compounds remain separate, so for !same_species comparisons we need to check 
    // if they should match with other groups
    if(!ignore_environment || check_other_grouping){
      comparison_schemes = compare::checkForBetterMatches(comparison_schemes, oss, FileMESSAGE, num_proc, false, same_species, scale_volume, optimize_match, ignore_symmetry, ignore_Wyckoff, true, clean_unmatched, false, quiet); 
    }

    // regroup comparisons based on misfit value
    if(num_mismatches==0 && !single_comparison_round){
      compare::appendStructurePrototypes(comparison_schemes, final_prototypes, clean_unmatched, quiet); //DX 20190506 -- added clean_unmatched
    }

    // Loop: continue comparison until all strucutures are matched or all comparisons schemes exhaused
    while(num_mismatches!=0){
      // regroup comparisons based on misfit value
      compare::appendStructurePrototypes(comparison_schemes, final_prototypes, clean_unmatched, quiet); //DX 20190506 -- added clean unmatched

      // return if only one round of comparison is requested
      if(single_comparison_round==true){return final_prototypes;}

      // reorder structures so minimum ICSD is the representative structure
      if(ICSD_comparison){ compare::representativePrototypeForICSDRuns(comparison_schemes); }

      // split into threads
      number_of_comparisons=0;
      for(uint i=0;i<comparison_schemes.size();i++){ number_of_comparisons += comparison_schemes[i].numberOfComparisons(); }
      num_comparison_threads = aurostd::min(num_proc,number_of_comparisons);

      if(number_of_comparisons>0){
        if(!quiet){
          message << "Continuing comparisons to match " << num_mismatches << " structures ...";
          pflow::logger(function_name, message, FileMESSAGE, logstream, _LOGGER_MESSAGE_);
        }
        if(LDEBUG) { cerr << function_name << ": Number of comparisons is not zero... " << number_of_comparisons << endl; }
#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
        // THREADED VERISON - START

        // split into threads
        start_indices.clear(); end_indices.clear();
        splitComparisonIntoThreads(comparison_schemes, num_comparison_threads, start_indices, end_indices);
        threads.clear();
        // run threads
        for(uint n=0; n<num_comparison_threads; n++){
          threads.push_back(new std::thread(&compare::runComparisonThreads,std::ref(comparison_schemes),std::ref(start_indices[n]),std::ref(end_indices[n]),
                std::ref(oss),same_species,scale_volume,optimize_match,store_comparison_logs));
        }        
        // join threads
        for(uint t=0;t<threads.size();t++){
          threads[t]->join();
          delete threads[t];
        }
        // THREADED VERISON - END
#else
        //SINGLE THREAD - START

        start_indices.clear(); end_indices.clear();
        uint single_thread=1;
        splitComparisonIntoThreads(comparison_schemes, single_thread, start_indices, end_indices);
        compare::runComparisonThreads(comparison_schemes,start_indices[0],end_indices[0],
            oss,same_species,scale_volume,optimize_match,store_comparison_logs);
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
          << "Contact David Hicks (david.hicks@duke.edu)." << endl;
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
      if(comparison_schemes[j].structures_duplicate_names.size()+1 == divisor_pairs[divisor_index].second){ //+1 to include representative
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

    string function_name = "compare::createStructurePermutations()";
    stringstream message;

    // ---------------------------------------------------------------------------
    // First, separate by stoichiometry
    for(uint i=0;i<name_order.size(); i++){
      bool scheme_created=false;
      string name = "";
      for(uint j=0;j<name_order[i].size();j++){name+=name_order[i][j];}
      if(i==0){
        StructurePrototype tmp;
        tmp.structure_representative_name=name;
        tmp.structure_representative=vxstrs[i];
        tmp.structure_representative_generated=true;
        tmp.structure_representative_from="input";
        tmp.number_of_types=vxstrs[i].num_each_type.size();
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
          if(same_species==true && matchableSpecies(vxstrs[i],comparison_schemes[j].structure_representative,same_species)==true){
            same_material_stoich=true;
          }
          else if(same_species==false){
            same_material_stoich=true;
          }
          if(same_material_stoich==true &&
             matchableWyckoffPositions(permutation_grouped_Wyckoff_positions[i], comparison_schemes[j].grouped_Wyckoff_positions,same_species)){
            comparison_schemes[j].structures_duplicate_names.push_back(name);
            comparison_schemes[j].structures_duplicate.push_back(vxstrs[i]);
            comparison_schemes[j].structures_duplicate_generated.push_back(true);
            comparison_schemes[j].structures_duplicate_from.push_back("input");
            comparison_schemes[j].misfits_duplicate.push_back(-1.0);
            scheme_created=true;
            break;
          }
        }
        if(scheme_created==false){
          StructurePrototype tmp;
          tmp.structure_representative_name=name;
          tmp.structure_representative_generated=true;
          tmp.structure_representative_from="input";
          tmp.structure_representative=vxstrs[i];
          tmp.number_of_types=vxstrs[i].num_each_type.size();
          tmp.number_of_atoms=vxstrs[i].atoms.size();
          tmp.grouped_Wyckoff_positions=permutation_grouped_Wyckoff_positions[i];
          comparison_schemes.push_back(tmp);
        }
      }
    }
    // ---------------------------------------------------------------------------
    // check number of sets and elements in set are consistent with number theory
    // (Groupings/sets must be a divisor of the total number of permutations)
    if(!checkNumberOfGroupings(comparison_schemes,name_order.size())){
      message << "Initial groupings of permutations do not follow number theory." << endl; 
      for(uint i=0;i<comparison_schemes.size();i++){ //DX 20190601 - added more info
        message << comparison_schemes[i] << endl;
      }
      message << "Please contact David Hicks (david.hicks@duke.edu) and provide the corresponding example." << endl;
      throw aurostd::xerror(function_name,message,_RUNTIME_ERROR_);
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
        if(comparison_schemes[i].structure_representative_name == name){
          representative_permutation_num = j;
          break;
        }
      }
      //Find proto (potential duplicate) permutation number
      uint min_even_duplicate_permutation_num = 1e9;
      uint duplicate_permutation_num = 0;
      uint min_duplicate_index = 0;
      for(uint p=0; p<comparison_schemes[i].structures_duplicate_names.size(); p++){
        for(uint j=0; j<name_order.size(); j++){			
          string name="";
          for(uint k=0;k<name_order[j].size();k++){name+=name_order[j][k];}
          if(comparison_schemes[i].structures_duplicate_names[p] == name){
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
        comparison_schemes[i].structures_duplicate_names.push_back(comparison_schemes[i].structure_representative_name);  
        comparison_schemes[i].structures_duplicate.push_back(comparison_schemes[i].structure_representative);  
        comparison_schemes[i].structures_duplicate_generated.push_back(comparison_schemes[i].structure_representative_generated);  
        comparison_schemes[i].structures_duplicate_from.push_back(comparison_schemes[i].structure_representative_from);  
        comparison_schemes[i].structure_representative_name = comparison_schemes[i].structures_duplicate_names[min_duplicate_index];  
        comparison_schemes[i].structure_representative = comparison_schemes[i].structures_duplicate[min_duplicate_index];  
        comparison_schemes[i].structure_representative_generated = comparison_schemes[i].structures_duplicate_generated[min_duplicate_index];  
        comparison_schemes[i].structure_representative_from = comparison_schemes[i].structures_duplicate_from[min_duplicate_index];  
        comparison_schemes[i].structures_duplicate_names.erase(comparison_schemes[i].structures_duplicate_names.begin()+min_duplicate_index);
        comparison_schemes[i].structures_duplicate.erase(comparison_schemes[i].structures_duplicate.begin()+min_duplicate_index);
        comparison_schemes[i].structures_duplicate_generated.erase(comparison_schemes[i].structures_duplicate_generated.begin()+min_duplicate_index);
        comparison_schemes[i].structures_duplicate_from.erase(comparison_schemes[i].structures_duplicate_from.begin()+min_duplicate_index);
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
      for(uint j=0; j<comparison_schemes[i].misfits_duplicate.size(); j++){
        if(comparison_schemes[i].misfits_duplicate[j] > 0.1 || aurostd::isequal(comparison_schemes[i].misfits_duplicate[j],-1.0,1e-6)){
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
      for(uint j=0; j<comparison_schemes[i].misfits_duplicate.size(); j++){
        //DX TEST if(comparison_schemes[i].misfits_duplicate[j] > 0.1 || comparison_schemes[i].misfits_duplicate[j] == -1.0 ){
        if(comparison_schemes[i].misfits_duplicate[j] > 0.1 || comparison_schemes[i].misfits_duplicate[j]+1.0 <_ZERO_TOL_){ //DX 20190301 - more robust
          // First, store any family prototype information
          if(comparison_schemes[i].misfits_duplicate[j] > 0.1 && comparison_schemes[i].misfits_duplicate[j] <= 0.2){
            comparison_schemes[i].putDuplicateAsFamily(j); //DX 20190814 - consolidated below into single function
            //DX 20190814 [OBSOLETE] - comparison_schemes[i].structures_family_names.push_back(comparison_schemes[i].structures_duplicate_names[j]);
            //DX 20190814 [OBSOLETE] - comparison_schemes[i].misfits_family.push_back(comparison_schemes[i].misfits_duplicate[j]);
            //DX 20190814 [OBSOLETE] - comparison_schemes[i].structures_family_grouped_Wyckoff_positions.push_back(comparison_schemes[i].structures_duplicate_grouped_Wyckoff_positions[j]); //DX 20190814 
            //DX 20190814 [OBSOLETE] - //DX 20190424 - store properties of family structures - START
            //DX 20190814 [OBSOLETE] - if(comparison_schemes[i].property_names.size()!=0){
            //DX 20190814 [OBSOLETE] -   comparison_schemes[i].properties_structures_family.push_back(comparison_schemes[i].properties_structures_duplicate[j]);
            //DX 20190814 [OBSOLETE] - }
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
        else if(comparison_schemes[i].misfits_duplicate[j] <= 0.1 && !std::signbit(comparison_schemes[i].misfits_duplicate[j])){
          //comparison_schemes[i].structures_duplicate.erase(comparison_schemes[i].structures_duplicate.begin()+j);
          // check if the structure is generated first
          if(comparison_schemes[i].structures_duplicate_generated[j]){comparison_schemes[i].structures_duplicate[j].Clear(); comparison_schemes[i].structures_duplicate_generated[j]=false; } //DX 20190303 - update generated flag
        }
      }
      //DX 20181220 - can clear representative structure since we will no longer compare it (save memory)
      //DX 20190521 [BREAKS WITH AURL MODE IN PARALLEL] comparison_schemes[i].structure_representative.Clear();
      //DX 20190521 [BREAKS WITH AURL MODE IN PARALLEL] comparison_schemes[i].structure_representative_generated=false;

      // if not quiet, print the comparison results to the screen 
      // (useful for long comparison times or if the program terminates early)
      if(!quiet){
        message << "Identified unique prototype: " << endl;
        message << "   prototype=" << comparison_schemes[i].structure_representative_name << endl;
        if(comparison_schemes[i].structures_duplicate_names.size()==0){
          message << "   No duplicates. " << endl;
        }
        else {
          message << "   " << setw(80) << std::left << "List of duplicates"
                  << setw(15) << std::left << "misfit value" << endl;
          message << "   " << setw(80) << std::left 
                  << "-----------------------------------------------------------------------------------------------" << endl;
          for(uint d=0;d<comparison_schemes[i].structures_duplicate_names.size();d++){
            message << "   " << setw(80) << std::left << comparison_schemes[i].structures_duplicate_names[d]
                    << setw(15) << std::left << comparison_schemes[i].misfits_duplicate[d] << endl;
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
            matchableSpecies(final_prototypes[i].structure_representative,final_prototypes[j].structure_representative,same_species)==true &&
            final_prototypes[i].stoichiometry==final_prototypes[j].stoichiometry && 
            final_prototypes[i].Pearson==final_prototypes[j].Pearson &&
            !matchableSpaceGroups(final_prototypes[i].space_group,final_prototypes[j].space_group)) ||
           // If same_species==false
           (same_species==false && 
            final_prototypes[i].stoichiometry==final_prototypes[j].stoichiometry &&
            final_prototypes[i].Pearson==final_prototypes[j].Pearson && 
            !matchableSpaceGroups(final_prototypes[i].space_group,final_prototypes[j].space_group))
          ){
          double final_misfit=-1.0;
          bool scale_volume=true; //default is true
          bool optimize_match=false; //default is false
          aflowCompareStructure(num_proc, final_prototypes[i].structure_representative, 
					 final_prototypes[j].structure_representative, same_species, 
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
        final_prototypes[sg_ind].structures_duplicate_names.push_back(final_prototypes[other_ind].structure_representative_name);
        final_prototypes[sg_ind].structures_duplicate_names.insert(final_prototypes[sg_ind].structures_duplicate_names.end(),
            final_prototypes[other_ind].structures_duplicate_names.begin(),
            final_prototypes[other_ind].structures_duplicate_names.end());
        final_prototypes[sg_ind].structures_duplicate.push_back(final_prototypes[other_ind].structure_representative);
        final_prototypes[sg_ind].structures_duplicate.insert(final_prototypes[sg_ind].structures_duplicate.end(),
            final_prototypes[other_ind].structures_duplicate.begin(),
            final_prototypes[other_ind].structures_duplicate.end());
        final_prototypes[sg_ind].structures_duplicate_compounds.push_back(final_prototypes[other_ind].structure_representative_compound);
        final_prototypes[sg_ind].structures_duplicate_compounds.insert(final_prototypes[sg_ind].structures_duplicate_compounds.end(),
            final_prototypes[other_ind].structures_duplicate_compounds.begin(),
            final_prototypes[other_ind].structures_duplicate_compounds.end());
        final_prototypes[sg_ind].structures_duplicate_generated.push_back(final_prototypes[other_ind].structure_representative_generated);
        final_prototypes[sg_ind].structures_duplicate_generated.insert(final_prototypes[sg_ind].structures_duplicate_generated.end(),
            final_prototypes[other_ind].structures_duplicate_generated.begin(),
            final_prototypes[other_ind].structures_duplicate_generated.end());
        final_prototypes[sg_ind].structures_duplicate_from.push_back(final_prototypes[other_ind].structure_representative_from);
        final_prototypes[sg_ind].structures_duplicate_from.insert(final_prototypes[sg_ind].structures_duplicate_from.end(),
            final_prototypes[other_ind].structures_duplicate_from.begin(),
            final_prototypes[other_ind].structures_duplicate_from.end());
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

    bool roff=true; //round off

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
          //LDEBUG cerr << j << " compound: " << final_prototypes[j].structure_representative_compound; //DX 20190311
          ss_out << final_prototypes[j].structure_representative_compound; //DX 20190311
          ss_out << "  SG=#" << final_prototypes[j].space_group;
          // ORIG ss_out << "  Wyckoffs=" << compare::printWyckoffString(final_prototypes[j].grouped_Wyckoff_positions,true) << endl;
          ss_out << "  Wyckoffs=" << compare::printWyckoffString(final_prototypes[j].grouped_Wyckoff_positions,true); //DX 20190228 - remove count
          uint number_of_duplicates = final_prototypes[j].numberOfDuplicates(); //DX 20190506 - made function
          ss_out << "  duplicate_compounds=" << number_of_duplicates << endl; //DX 20190228 - add count
          if(final_prototypes[j].aflow_label.size()!=0){
            ss_out << "  aflow_label=" << final_prototypes[j].aflow_label << endl; 
            ss_out << "  aflow_parameter_list=" << aurostd::joinWDelimiter(final_prototypes[j].aflow_parameter_list,",") << endl; 
            ss_out << "  aflow_parameter_values=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(final_prototypes[j].aflow_parameter_values,8,roff),",") << endl;
          } 
          if(final_prototypes[j].matching_aflow_prototypes.size()!=0){
            ss_out << "  matching_aflow_prototypes=" << aurostd::joinWDelimiter(final_prototypes[j].matching_aflow_prototypes,",") << endl; 
          } 
          if(final_prototypes[j].properties_structure_representative.size()!=0){
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
          if(final_prototypes[j].properties_structure_representative.size()!=0){
            ss_out << std::string(indent_spacing+structure_spacing+misfit_spacing, '-');
            if(final_prototypes[j].property_names.size()!=0){
              ss_out << std::string(final_prototypes[j].properties_structure_representative.size()*property_spacing, '-');
              ss_out << std::string(final_prototypes[j].property_names.size()*property_spacing, '-');
            }
            ss_out << endl;
          }
          ss_out << "  " << setw(structure_spacing) << std::left << "prototype="+final_prototypes[j].structure_representative_name;
          if(final_prototypes[j].properties_structure_representative.size()!=0){
            ss_out << setw(misfit_spacing) << std::right << "-";
            for(uint l=0;l<final_prototypes[j].properties_structure_representative.size();l++){
              ss_out << setw(property_spacing) << std::right << final_prototypes[j].properties_structure_representative[l];
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
          ss_out << "  Wyckoffs=" << compare::printWyckoffString(final_prototypes[j].grouped_Wyckoff_positions,true);
          uint number_of_duplicates = final_prototypes[j].numberOfDuplicates(); //DX 20190506 - made function
          ss_out << "  structures_duplicate=" << number_of_duplicates; //DX 20190228 - add count
          uint number_duplicate_compounds = 0;
          for(uint k=0;k<final_prototypes[j].number_compounds_matching_duplicate.size();k++){
            number_duplicate_compounds+=final_prototypes[j].number_compounds_matching_duplicate[k];
          }
          number_duplicate_compounds+= number_of_duplicates+final_prototypes[j].number_compounds_matching_representative; //DX 20190321 - need to update variable, otherwise may not enter if statement
          if(number_duplicate_compounds!=0){
            ss_out << "  duplicate_compounds=" << number_duplicate_compounds; //DX 20190228 - add count
          }
          ss_out << endl;
          if(final_prototypes[j].aflow_label.size()!=0){
            ss_out << "  aflow_label=" << final_prototypes[j].aflow_label << endl; 
            ss_out << "  aflow_parameter_list=" << aurostd::joinWDelimiter(final_prototypes[j].aflow_parameter_list,",") << endl; 
            ss_out << "  aflow_parameter_values=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(final_prototypes[j].aflow_parameter_values,8,roff),",") << endl;
          }
          if(final_prototypes[j].matching_aflow_prototypes.size()!=0){
            ss_out << "  matching_aflow_prototypes=" << aurostd::joinWDelimiter(final_prototypes[j].matching_aflow_prototypes,",") << endl; 
          } 
          ss_out << "  " << setw(structure_spacing) << std::left << "prototype="+final_prototypes[j].structure_representative_name;
          // perhaps add which permutations are duplicates
          if(final_prototypes[j].atom_decorations_equivalent.size()!=0){
            //ss_out << endl << "  " << setw(structure_spacing) << std::left << "unique atom decorations="+aurostd::joinWDelimiter(final_prototypes[j].atom_decorations_equivalent,",");
            vector<string> unique_decorations;
            for(uint d=0;d<final_prototypes[j].atom_decorations_equivalent.size();d++){ unique_decorations.push_back(final_prototypes[j].atom_decorations_equivalent[d][0]); }
            ss_out << endl << "  " << setw(structure_spacing) << std::left << "unique atom decorations="+aurostd::joinWDelimiter(unique_decorations,",");
          }
        }
        ss_out << endl;
        ss_out << std::string(indent_spacing+structure_spacing+misfit_spacing, '-');
        if(final_prototypes[j].property_names.size()!=0){
          ss_out << std::string(final_prototypes[j].property_names.size()*property_spacing, '-');
        }
        ss_out << endl;
        if(final_prototypes[j].structures_duplicate_names.size()!=0 && final_prototypes[j].properties_structure_representative.size()==0){
          ss_out << "  " << setw(structure_spacing) << std::left << "list of duplicates";
          ss_out << setw(misfit_spacing) << std::right << "misfit";
        }
        else if(final_prototypes[j].structures_duplicate_names.size()==0){
          ss_out << "  " << setw(structure_spacing) << std::left << "no duplicates";
        }
        if(final_prototypes[j].properties_structure_representative.size()==0){
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
        if(final_prototypes[j].structures_duplicate_names.size()>0){
          for(uint k=0;k<final_prototypes[j].structures_duplicate_names.size();k++){
            ss_out << "  " << setw(structure_spacing) << std::left << final_prototypes[j].structures_duplicate_names[k] 
              << setw(misfit_spacing) << std::right << final_prototypes[j].misfits_duplicate[k];
            if(final_prototypes[j].property_names.size()!=0){
              for(uint l=0;l<final_prototypes[j].properties_structures_duplicate[k].size();l++){
                ss_out << setw(property_spacing) << std::right << final_prototypes[j].properties_structures_duplicate[k][l];
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
// printStructureMappingResults - Displays comprehensive mapping information
// ***************************************************************************
namespace compare{
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
    string mode){

    if(aurostd::toupper(mode) == "TEXT" || aurostd::toupper(mode) == "TXT"){
      oss << endl <<"**************************** RESULTS ****************************"<<endl;
      if(misfit<=0.1){
        oss << endl <<"MISFIT:			" << misfit << "  STRUCTURES ARE COMPATIBLE" << endl;
      }
      else if(misfit<=0.2){
        oss << endl <<"MISFIT:      " << misfit <<"  STRUCTURES ARE IN THE SAME FAMILY" << endl;
      }
      else {
        oss << endl <<"MISFIT:			" << misfit <<"  STRUCTURES ARE INCOMPATIBLE (No match found)" << endl;
      }
      oss <<"----------------------------------------------------"<<endl;
      oss << "Figure of Deviation:	" << lattice_deviation << endl;
      oss << "Figure of Displacement:	" << coordinate_displacement << endl;
      oss << "Figure of Failure:	" << failure << endl;
      if(magnetic_analysis){
        oss << "Figure of Magnetic Displacement:	" << magnetic_displacement << endl;
        oss << "Figure of Magnetic Failure:	" << magnetic_failure << endl;
      }
      oss <<"----------------------------------------------------"<<endl;
      printMatch(matching_indices_1,matching_indices_2,minimum_distances,xstr_transformed,xstr_reference,oss);
      oss <<"----------------------------------------------------"<<endl;
      oss << "FINAL - REFERENCE STRUCTURE: " << endl;	
      oss << xstr_reference << endl;
      oss <<"----------------------------------------------------"<<endl;
      oss << "FINAL - MAPPED STRUCTURE: " << endl;
      oss << xstr_transformed;
      oss << endl << "*********************  THE END - FINE  **********************" << endl << endl;
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
    string function_name = "compare::sameSpecies()";

    // ---------------------------------------------------------------------------
    // Check number of types
    if(xstr1.num_each_type.size() != xstr2.num_each_type.size()){
      // Display counts
      if(display && LDEBUG) { //DX 20190702 - condense if-statements
        cerr << function_name << ": Number of element types are not the same." 
          << " xstr 1: " << xstr1.num_each_type.size() 
          << " and xstr2: " << xstr2.num_each_type.size() << endl;
      }
      return false; //DX 20190702 - bug fix, should not be in if-statement
    }

    // ---------------------------------------------------------------------------
    // Check counts and species
    // sort, then check (to check if matchable) = fast
    deque<int> xstr1_num_each_type = xstr1.num_each_type;
    deque<int> xstr2_num_each_type = xstr2.num_each_type;

    std::sort(xstr1_num_each_type.begin(), xstr1_num_each_type.end());
    std::sort(xstr2_num_each_type.begin(), xstr2_num_each_type.end());

    if(xstr1_num_each_type!=xstr2_num_each_type){
      if(display && LDEBUG) { cerr << function_name << ": Number of each type of element are incompatible." << endl; }
      return false;
    }
    if(display && LDEBUG) { cerr << function_name << ": Number of each type of element are compatible; proceeding." << endl; }

    //DX 20190702 [OBSOLETE - not robust and slow] for(uint i=0;i<xstr1.num_each_type.size();i++){
    //DX 20190702 [OBSOLETE - not robust and slow]   bool matched = false;
    //DX 20190702 [OBSOLETE - not robust and slow]   for(uint j=0;j<xstr2.num_each_type.size();j++){
    //DX 20190702 [OBSOLETE - not robust and slow]     // DX IS SPECIES CHECK TOO STRICT? if(xstr1.num_each_type[i] == xstr2.num_each_type[j] &&
    //DX 20190702 [OBSOLETE - not robust and slow]     // DX IS SPECIES CHECK TOO STRICT?   xstr1.species[i] == xstr2.species[j]){
    //DX 20190702 [OBSOLETE - not robust and slow]     if(xstr1.num_each_type[i] == xstr2.num_each_type[j]){
    //DX 20190702 [OBSOLETE - not robust and slow]       matched = true;
    //DX 20190702 [OBSOLETE - not robust and slow]       break;
    //DX 20190702 [OBSOLETE - not robust and slow]     }
    //DX 20190702 [OBSOLETE - not robust and slow]   }
    //DX 20190702 [OBSOLETE - not robust and slow]   if(matched == false){
    //DX 20190702 [OBSOLETE - not robust and slow]     if(display==true){ 
    //DX 20190702 [OBSOLETE - not robust and slow]       if(LDEBUG) {
    //DX 20190702 [OBSOLETE - not robust and slow]         cerr << "compare::WARNING:: TYPE OF ATOMIC SPECIES OR NUMBER PER TYPE ARE NOT THE SAME..." << endl;  
    //DX 20190702 [OBSOLETE - not robust and slow]       }
    //DX 20190702 [OBSOLETE - not robust and slow]     }
    //DX 20190702 [OBSOLETE - not robust and slow]     return false;
    //DX 20190702 [OBSOLETE - not robust and slow]   }
    //DX 20190702 [OBSOLETE - not robust and slow] }
    //DX 20190702 [OBSOLETE - not robust and slow] if(display==true){
    //DX 20190702 [OBSOLETE - not robust and slow]   if(LDEBUG) {
	  //DX 20190702 [OBSOLETE - not robust and slow]     cerr << "compare::NUMBER AND TYPE OF ATOMIC SPECIES ARE THE SAME...PROCEEDING..." << endl;
    //DX 20190702 [OBSOLETE - not robust and slow]   }
    //DX 20190702 [OBSOLETE - not robust and slow] }
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
// getLeastFreqentAtomSpecies
// ***************************************************************************
namespace compare{
  vector<string> getLeastFrequentAtomSpecies(const xstructure& xstr) {

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
// sortBySecondPair 
// ***************************************************************************
namespace compare{
  bool sortBySecondPair(const std::pair<string,uint>& a, const std::pair<string,uint>& b) {
    return (a.second<b.second);
  }
}

// ***************************************************************************
// sortSpeciesByFrequency
// ***************************************************************************
namespace compare{
  vector<string> sortSpeciesByFrequency(const xstructure& xstr) {

    vector<std::pair<string,uint> > species_and_counts;

    for(uint i=0;i<xstr.num_each_type.size();i++){
      std::pair<string,uint> tmp_pair;
      tmp_pair.first=xstr.species[i]; tmp_pair.second=xstr.num_each_type[i];
      species_and_counts.push_back(tmp_pair);
    }

    std::stable_sort(species_and_counts.begin(),species_and_counts.end(),compare::sortBySecondPair);

    vector<string> reordered_species;
    for(uint i=0;i<species_and_counts.size();i++){ reordered_species.push_back(species_and_counts[i].first); }
    
    return reordered_species;
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

    // ---------------------------------------------------------------------------
    // switch between relative and absolute tolerance
    bool relative = false;

    // ---------------------------------------------------------------------------
    // relative
    if(relative){
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
    // ---------------------------------------------------------------------------
    // absolute
    else{
      double tol_length=1.0, tol_angle=5; // 1 Angstrom; 5 degrees

      if( abs(d1(1)-d2(1)) < tol_length &&
          abs(d1(2)-d2(2)) < tol_length &&
          abs(d1(3)-d2(3)) < tol_length &&
          abs(d1(4)-d2(4)) < tol_angle &&
          abs(d1(5)-d2(5)) < tol_angle &&
          abs(d1(6)-d2(6)) < tol_angle
        ){
        return false;
      }
      else {	
        return true;
      }
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
// Reset dims for RadiusSphereLattice() 
// ***************************************************************************
namespace compare{
  void resetLatticeDimensions(const xmatrix<double>& lattice, double radius, xvector<int>& dims,
      vector<xvector<double> >& l1, vector<xvector<double> >& l2, 
      vector<xvector<double> >& l3, vector<int>& a_index, 
      vector<int>& b_index, vector<int>& c_index){

    // resets the lattice dimensions (dims) based on radius
    // generates lattice vectors (l1,l2,l3) right away = speed increase
    // stores dimension indices (a_index,b_index,c_index)
    // new dims explore order : zeroth cell to max dims = speed increase 
    // (can break early if match is found)
    // DX create function date: 20190705

    // ---------------------------------------------------------------------------
    // get new dimensions based on radius
    if(radius<=_ZERO_TOL_){ dims[1]=1; dims[2]=1; dims[3]=1; }
    else{ dims=LatticeDimensionSphere(lattice,radius); }
    //cerr << "using dims: " << dims << endl; 
    // ---------------------------------------------------------------------------
    // clear old 
    l1.clear(); l2.clear(); l3.clear();
    a_index.clear(); b_index.clear(); c_index.clear();

    // ---------------------------------------------------------------------------
    // [NEW] - go from zeroth cell out
    // more likely to find match close to origin, why start so far away
    
    // push back zeroth cell : dims[1]=dims[2]=dims[3]=0
    l1.push_back(0*lattice(1));a_index.push_back(0);
    l2.push_back(0*lattice(2));b_index.push_back(0);
    l3.push_back(0*lattice(3));c_index.push_back(0);

    // push back 1,-1,2,-2,...dims,-dims
    for(int a=1;a<=dims[1];a++){l1.push_back(a*lattice(1));a_index.push_back(a); l1.push_back(-a*lattice(1));a_index.push_back(-a);}
    for(int b=1;b<=dims[2];b++){l2.push_back(b*lattice(2));b_index.push_back(b); l2.push_back(-b*lattice(2));b_index.push_back(-b);}
    for(int c=1;c<=dims[3];c++){l3.push_back(c*lattice(3));c_index.push_back(c); l3.push_back(-c*lattice(3));c_index.push_back(-c);}

  }
}

// ***************************************************************************
// minimumCoordinationShellLatticeOnly() 
// ***************************************************************************
namespace compare{
  void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
      double& min_dist, uint& frequency, vector<xvector<double> >& coordinates){
    
    // determine the minimum coordination shell of the lattice
    // i.e., find the set of closest neighbors to the origin
    // (overload: uses lattice radius) 

    // ---------------------------------------------------------------------------
    // determine necessary search radius
    double radius=RadiusSphereLattice(lattice);
  
    minimumCoordinationShellLatticeOnly(lattice, min_dist, frequency, coordinates, radius);
  }
}

// ***************************************************************************
// minimumCoordinationShellLatticeOnly() 
// ***************************************************************************
namespace compare{
  void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice,
      double& min_dist, uint& frequency, vector<xvector<double> >& coordinates, double radius){
    
    // determine the minimum coordination shell of the lattice
    // i.e., find the set of closest neighbors to the origin
    // (overload: instantiates lattice dimension information) 

    // ---------------------------------------------------------------------------
    // instantiate lattice vectors 
    vector<xvector<double> > l1, l2, l3; 
    vector<int> a_index, b_index, c_index;
    xvector<int> dims(3); dims[1]=dims[2]=dims[3]=0; // declare/reset
    resetLatticeDimensions(lattice,radius,dims,l1,l2,l3,a_index,b_index,c_index);
  
    minimumCoordinationShellLatticeOnly(lattice, dims, l1, l2, l3, 
        a_index, b_index, c_index, 
        min_dist, frequency, coordinates, radius);
  }
}

// ***************************************************************************
// minimumCoordinationShellLatticeOnly() 
// ***************************************************************************
namespace compare{
  void minimumCoordinationShellLatticeOnly(const xmatrix<double>& lattice, xvector<int>& dims,
      vector<xvector<double> >& l1, vector<xvector<double> >& l2, vector<xvector<double> >& l3, 
      vector<int>& a_index, vector<int>& b_index, vector<int>& c_index, 
      double& min_dist, uint& frequency, vector<xvector<double> >& coordinates,
      double radius){

    // determine the minimum coordination shell environment of the lattice
    // i.e., find the set of closest neighbors to the origin
    // stores l1, l2, l3, a_index, b_index, and c_index for external use
    // optional "radius" as enables more control over search space 
    // (and potential speed up, may not need to search as far as the lattice radius)
    
    xvector<double> tmp;

    // ---------------------------------------------------------------------------
    // reset lattice dimensions 
    resetLatticeDimensions(lattice,radius,dims,l1,l2,l3,a_index,b_index,c_index);

    double relative_tolerance = 10.0; // coordination shell thickness is ten percent of minimum distance

    // ---------------------------------------------------------------------------
    // loop through lattice vectors (stored before-hand in l1,l2,l3)
    for(uint m=0;m<l1.size();m++){
      xvector<double> a_component = l1[m];                  // DX : coord1-coord2+a*lattice(1)
      for(uint n=0;n<l2.size();n++){
        xvector<double> ab_component = a_component + l2[n]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
        for(uint p=0;p<l3.size();p++){
          if(!(m==0 && n==0 && p==0)){
            tmp = ab_component + l3[p];                     // DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
            double tmp_mod = aurostd::modulus(tmp);
            // ---------------------------------------------------------------------------
            // if found a new minimum distance and update coordination/frequency and coordinate 
            if(tmp_mod<min_dist){
              // ---------------------------------------------------------------------------
              // if new distance is close to the original it is the same coordination shell (add to coordination)
              // otherwise, reset coordination shell
              // DX - FIXED TOL (bad for undecorated prototypes) - if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
              if(aurostd::isequal(tmp_mod,min_dist,(min_dist/relative_tolerance))){ frequency+=1; coordinates.push_back(tmp); } // tenth of min_dist
              else{ frequency=1; coordinates.clear(); coordinates.push_back(tmp); } //initialize
              min_dist=tmp_mod;
              // ---------------------------------------------------------------------------
              // diminishing dims: if minimum distance changed, then we may not need to search as far
              // reset loop and search again based on new minimum distance
              if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
                resetLatticeDimensions(lattice,min_dist,dims,l1,l2,l3,a_index,b_index,c_index);
                m=n=p=0;
                frequency=0; //reset
              }
            }
            // DX - FIXED TOL (bad for undecorated prototypes) - else if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
            else if(aurostd::isequal(tmp_mod,min_dist,(min_dist/relative_tolerance))){ frequency+=1; coordinates.push_back(tmp); } // tenth of min dist
          }
        }
      }
    }
  }
}

// ***************************************************************************
// minimumCoordinationShell() 
// ***************************************************************************
namespace compare{
  void minimumCoordinationShell(const xstructure& xstr, uint center_index, 
      double& min_dist, uint& frequency, vector<xvector<double> >& coordinates){

    string type = "";
  
    minimumCoordinationShell(xstr, center_index, min_dist, frequency, coordinates, type);
  }
}

// ***************************************************************************
// minimumCoordinationShell() 
// ***************************************************************************
namespace compare{
  void minimumCoordinationShell(const xstructure& xstr, uint center_index, 
      double& min_dist, uint& frequency, vector<xvector<double> >& coordinates, const string& type){

    // determine the minimum coordination shell environment
    // "type" enables the search of environments by certain elements/types only
    // (e.g., find the neighborhood of oxygen atoms surrounding a magnesium center)
    
    xvector<double> tmp;

    // ---------------------------------------------------------------------------
    // instantiate lattice vectors 
    vector<xvector<double> > l1, l2, l3; 
    vector<int> a_index, b_index, c_index;
    xvector<int> dims(3); dims[1]=dims[2]=dims[3]=0; // declare/reset
    //resetLatticeDimensions(lattice,radius,dims,l1,l2,l3,a_index,b_index,c_index);

    double relative_tolerance = 10.0; // coordination shell thickness is ten percent of minimum distance

    for(uint ii=0; ii<xstr.atoms.size(); ii++){
      // ---------------------------------------------------------------------------
      // if atom ii is not environment center, find minimum distance between center atom ii's images 
      if(ii!=center_index && (xstr.atoms[ii].name == type || type == "")){ //DX 20191105 - added type=="" 
        xvector<double> incell_dist = xstr.atoms[center_index].cpos-xstr.atoms[ii].cpos;
        double incell_mod = aurostd::modulus(incell_dist);
        if(!(dims[1]==1 && dims[2]==1 && dims[3]==1) && incell_mod!=1e9){
          resetLatticeDimensions(xstr.lattice,incell_mod,dims,l1,l2,l3,a_index,b_index,c_index);
        }
        //DX 4/23/18 - running vector in each loop saves computations; fewer duplicate operations
        for(uint m=0;m<l1.size();m++){
          xvector<double> a_component = incell_dist + l1[m];    // DX : coord1-coord2+a*lattice(1)
          for(uint n=0;n<l2.size();n++){
            xvector<double> ab_component = a_component + l2[n]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
            for(uint p=0;p<l3.size();p++){
              tmp = ab_component + l3[p];                       // DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
              double tmp_mod = aurostd::modulus(tmp);
              if(tmp_mod<min_dist){
                //DX - FIXED TOL (bad for undecorated prototypes) - if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
                if(aurostd::isequal(tmp_mod,min_dist,(min_dist/relative_tolerance))){ frequency+=1; coordinates.push_back(tmp); } // tenth of min_dist
                else{ frequency=1; coordinates.clear(); coordinates.push_back(tmp); } //initialize
                min_dist=tmp_mod;
              }
              //DX - FIXED TOL (bad for undecorated prototypes) else if(aurostd::isequal(tmp_mod,min_dist,0.5)){ frequency+=1; } // within half an Angstrom
              else if(aurostd::isequal(tmp_mod,min_dist,(min_dist/relative_tolerance))){ frequency+=1; coordinates.push_back(tmp); } // tenth of min_dist
            }
          }
        }
      }
      // ---------------------------------------------------------------------------
      // if atom is environment center check its images, but only need to search as 
      // far as min_dist or lattice_radius (whichever is smaller)
      else if(ii==center_index && (xstr.atoms[ii].name == type || type == "")){ //DX 20191105 - added type==""
        double lattice_radius=RadiusSphereLattice(xstr.lattice);
        double search_radius=min(lattice_radius,min_dist);

        // ---------------------------------------------------------------------------
        // use variant that stores the lattice dimension information so it can be 
        // updated for the "minimumCoordinationShell" function
        minimumCoordinationShellLatticeOnly(xstr.lattice, dims, l1, l2, l3, 
            a_index, b_index, c_index, 
            min_dist, frequency, coordinates, search_radius);
      }
    }
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
  bool findMatch(const deque<_atom>& xstr1_atoms, const deque<_atom>& PROTO_atoms,
        const xmatrix<double>& PROTO_lattice,
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

    bool LDEBUG=(false || XHOST.DEBUG);
    bool VERBOSE=false;

    string function_name = "compare::findMatch()";

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

    //xmatrix<double> lattice=PROTO.lattice;
    xmatrix<double> lattice=PROTO_lattice;

    double tmp = 1e9;
    uint i1_real=0;
    uint i2_real=0;
    string i1_name = "";
    string i2_name = "";

    //DX 20190226 [BETA] xvector<double> best_centroid1 = centroid_with_PBC(xstr1); 
    //DX 20190226 [BETA] xvector<double> best_centroid2 = centroid_with_PBC(PROTO); 
    
    vector<xvector<double> > l1, l2, l3;
    vector<int> a_index, b_index, c_index;
    //double radius=RadiusSphereLattice(lattice); //DX 20190701 - use robust method
    xvector<int> dims(3); //DX 20190701 - use robust method
    //resetLatticeDimensions(lattice,radius,dims,l1,l2,l3,a_index,b_index,c_index);

    for(j=0;j<xstr1_atoms.size();j++){
      //cerr << "xstr1.atoms[j]: " << xstr1.atoms[j] << endl;
      std::pair<xvector<double>,xvector<double> > tmp_pair;
      xvector<double> tmp_xvec = xstr1_atoms[j].cpos;
      tmp_pair.first = tmp_xvec;
      vdiffs.clear();
      double match_dist=1e9;
      double prev_match_dist=0; //to avoid recalculating dims if nothing has changed
      dims[1]=dims[2]=dims[3]=0; //reset
      for(k=0;k<PROTO_atoms.size();k++){
        //cerr << "PROTO.atoms[k]: " << PROTO.atoms[k] << endl;
        //cerr << "[orig] dims: " << dims << endl;
        //DX 20190701 - use diminishing dims - START
        if(match_dist<prev_match_dist){
          if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
            resetLatticeDimensions(lattice,match_dist,dims,l1,l2,l3,a_index,b_index,c_index);
            prev_match_dist=match_dist;
          }
        }
        double dist=1e9;
        xvector<double> min_xvec;
        xvector<double> incell_dist = xstr1_atoms[j].cpos-PROTO_atoms[k].cpos;
        double incell_mod = aurostd::modulus(incell_dist);
        //DX 20190701 - use diminishing dims - START
        if(incell_mod < match_dist){
          if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
            resetLatticeDimensions(lattice,incell_mod,dims,l1,l2,l3,a_index,b_index,c_index);
            prev_match_dist=incell_mod;
          }
        }
        // Need to find the min distance; thus check distance between neighboring cells to find true minimum.
        //DX - running vector in each loop saves computations; fewer duplicate operations
        if(incell_mod>0.25){
          for(uint m=0;m<l1.size();m++){
            xvector<double> a_component = incell_dist + l1[m];    // DX : coord1-coord2+a*lattice(1)
            for(uint n=0;n<l2.size();n++){
              xvector<double> ab_component = a_component + l2[n]; // DX : coord1-coord2+a*lattice(1) + (b*lattice(2))
              for(uint p=0;p<l3.size();p++){
                tmp_xvec = ab_component + l3[p];                       // DX : coord1-coord2+a*lattice(1) + (b*lattice(2)) + (c*lattice(3))
                tmp=aurostd::modulus(tmp_xvec);
                if(tmp < dist){
                  i1 = j;
                  i2 = k;
                  dist = tmp;
                  min_xvec = tmp_xvec;
                }
                if(dist<0.25){ break; }
              }
              if(dist<0.25){ break; }
            }
            if(dist<0.25){ break; }
          }
        }
        else{
          if(incell_mod < dist){
            i1 = j;
            i2 = k;
            dist = incell_mod;
            min_xvec = incell_dist;
          }
        }
        //cerr << "match_dist: " << match_dist << endl;
        if(dist<match_dist){
          i1_real=i1;
          i2_real=i2;
          match_dist = dist;
          i1_name = xstr1_atoms[i1_real].name;
          i2_name = PROTO_atoms[i2_real].name;
          tmp_pair.second = min_xvec;
        }
        vdiffs.push_back(dist);
        //DX 20190701 - speed increase, not possible to match to anything else if less than quarter of an Angstrom
        // note this will truncate vdiffs, so if we need it, then do not use the break below
        if(dist<0.25){
          break;
        }
      }
      if(VERBOSE){
        cerr << function_name << ": set of distances for original structure's " << j << "th atom: " << endl; 
        for(uint d=0;d<vdiffs.size();d++){
          cerr << d << ": " << vdiffs[d] << endl;
        }
      }
      // Check if same species match
      if(type_match == 2){ // same species
        if(i1_name != i2_name){
          if(VERBOSE){
            cerr << function_name << ":WARNING: Matching species are not the same type, throwing out match (same species comparison)" << endl;
          }
          return false;
        }
      }
      // Check for one-to-one mappings 
      for(uint i=0;i<im1_name.size();i++){
        // Check if i1 index has multiple mappings
        if(i1_real == im1[i]){
          if(LDEBUG){
            cerr << "WARNING: Used the same index for matching in i1! (" << i1_real << " == " << im1[i] << ")"<< endl;
            cerr << "                                             i2! (" << i2_real << " == " << im2[i] << ")"<< endl;
            cerr << match_dist << " vs " << min_dists[i] << endl;
          }
          return false;
        }
        // Check if i2 index has multiple mappings
        if(i2_real == im2[i]){
          if(LDEBUG){
            cerr << "WARNING: Used the same index for matching in i2! (" << i2_real << " == " << im2[i] << ")"<< endl;
            cerr << "                                             i1! (" << i1_real << " == " << im1[i] << ")"<< endl;
            cerr << match_dist << " vs " << min_dists[i] << endl;
          }
          return false;
        }
        // Check if types are not consistently mapped to a single type
        if(i1_name == im1_name[i]){
          if(i2_name != im2_name[i]){
            if(LDEBUG){
              cerr << "WARNING: Matching one type of atom to more than one type! (" << i1_name << " == " << i2_name << " | " << im1_name[i] << " == " << im2_name[i] << ")" <<  endl;
              for(uint j=0;j<im1_name.size();j++){
            //  cerr << im1[j] << " == " << im2[j] << " | " << xstr1.atoms[im1[j]].cpos << " == " << PROTO.atoms[im2[j]].cpos << " (" << min_dists[i] << ") | " << im1_name[j] << " == " << im2_name[j] << endl;
              }
            //cerr << i1_real << " == " << i2_real << " | " << xstr1.atoms[i1_real].cpos << " == " << PROTO.atoms[i2_real].cpos << " (" << match_dist << ") | " << i1_name << " == " << i2_name << endl;            
            }
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
      //DX 20191112 [ORIG] dev.push_back((abs(diag_sum1[i]-diag_sum2[i])+abs(diag_diff1[i]-diag_diff2[i]))/(diag_sum2[i]+diag_diff2[i]));
      dev.push_back((abs(diag_sum2[i]-diag_sum1[i])+abs(diag_diff2[i]-diag_diff1[i]))/(diag_sum1[i]+diag_diff1[i])); //DX 20191112 - should be compared to the reference, (1) not (2)?

    d=1;
    for(i=0;i<dev.size();i++) 
      d=d*(1-dev[i]);
    d=1-d;

    return d;
  }
}

// ***************************************************************************
// compute LFA environment 
// ***************************************************************************
namespace compare{
  vector<AtomEnvironment> computeLFAEnvironment(const xstructure& xstr, bool unique_only){
    
    // computes an abridged atomic environment around all LFA atom centers
    // only determines the closest distance fore each atom type (i.e., one for each species)
    // this is a quick way to look at the atom environment before trying to compare
    // atom environment is invariant of unit cell representation/origin choice

    // ---------------------------------------------------------------------------
    // determine all LFA atoms in the structure (could be more than one) 
    vector<string> LFAs=getLeastFrequentAtomSpecies(xstr);

    // ---------------------------------------------------------------------------
    // compute all LFA environments, looping through each LFA type 
    vector<AtomEnvironment> all_environments_LFA;

    for(uint i=0;i<LFAs.size();i++){
      vector<AtomEnvironment> environments_LFA = getUniqueTypesAtomEnvironmentForLFA(xstr, LFAs[i], LFAs);
      // ---------------------------------------------------------------------------
      // may have non-primitive cell, but we only want unique/smallest set of information (fast)
      if(unique_only){
        for(uint j=0;j<environments_LFA.size();j++){
          bool duplicate = false;
          for(uint k=0;k<all_environments_LFA.size();k++){
            if(compatibleEnvironments(environments_LFA[j],all_environments_LFA[k],true,true,true)){
              duplicate = true;
              break;
            }
          }
          if(!duplicate){
            all_environments_LFA.push_back(environments_LFA[j]);
          }
        }
      }
      // ---------------------------------------------------------------------------
      // primitivized cell, push back all 
      else{
        all_environments_LFA.insert(all_environments_LFA.end(),environments_LFA.begin(),environments_LFA.end());
      }
    }

    return all_environments_LFA;
  }
}

// ***************************************************************************
// compute LFA environment 
// ***************************************************************************
namespace compare{
  bool compatibleEnvironmentSets(const vector<AtomEnvironment>& env_set1, 
      const vector<AtomEnvironment>& env_set2, bool same_species, bool exact_match){
 
    // determines if set of LFA environments are similar
    // same_species : requires the same atom decorations/types
    // exact_match  : signals if an exact match (remove duplicates) or conversely 
    //                if it possible to match later via structure comparison
    // if this used elsewhere, it may need to be modified since it only considers
    // the special case of comparing LFA environments

    // ---------------------------------------------------------------------------
    // if one is not calculated then we need to assume they may be compatible 
    if(env_set1.size()==0 || env_set2.size()==0){
      return true;
    }

    // ---------------------------------------------------------------------------
    // check if atom environment set has only one atom 
    // (signal more comprehensive environment comparison)
    bool only_single_LFA_atoms = true;
    string first_lfa_element = env_set1[0].center_element; //safe to acess element since I checked earlier
    for(uint i=1;i<env_set1.size();i++){ // start after first 1
      if(first_lfa_element ==env_set1[i].center_element){
        only_single_LFA_atoms = false;
        break;
      }
    }
    bool compare_frequency = only_single_LFA_atoms;
    // ---------------------------------------------------------------------------
    // if same species
    if(same_species){
      for(uint i=0;i<env_set1.size();i++){
        bool matched_set = false;
        vector<vector<string> > matched_species;
        for(uint j=0;j<env_set2.size();j++){
          matched_set = compatibleEnvironments(env_set1[i],env_set2[j],matched_species,same_species,compare_frequency,exact_match);
          if(matched_set) { break; }
        }
        if(!matched_set){ return false; }
      }
    }
    // ---------------------------------------------------------------------------
    // any species
    else{
      // this is a bit more complicated
      // need to think about this; come back later
      // for now return true so we do a robust match
      return true; 
    }
    return true;
  }
}

// ***************************************************************************
// compute LFA environment 
// ***************************************************************************
namespace compare{
  bool compatibleEnvironments(const AtomEnvironment& env_1, 
      const AtomEnvironment& env_2, bool same_species, bool compare_frequency, 
      bool exact_match){

    vector<vector<string> > matched_species;

    return compatibleEnvironments(env_1,env_2,matched_species,same_species,compare_frequency,exact_match);
  }
}

// ***************************************************************************
// compute LFA environment 
// ***************************************************************************
namespace compare{
  bool compatibleEnvironments(const AtomEnvironment& env_1, 
      const AtomEnvironment& env_2, vector<vector<string> > & matched_species, 
      bool same_species, bool compare_frequency, bool exact_match){

    // determines if set of LFA environments are similar
    // same_species     : requires the same atom decorations/types
    // exact_match      : signals if an exact match (remove duplicates) or conversely 
    //                    if it possible to match later via structure comparison
   
    bool LDEBUG=(false || XHOST.DEBUG);
    bool VERBOSE=false;
    string function_name = "compare::compatibleEnvironments()";

    double _TOL_EXACT_MATCH_ = 0.01; // hundredth of an Angstrom, perhaps put in header?
    double _TOL_RELATIVE_MATCH_ = 0.10; // ten percent, perhaps put in header? //DX 20190724 - changed from 0.25 to 0.1
    double _TOL_LOOSE_MATCH_ = aurostd::min(env_1.neighbor_distances)/2.0; // ten percent, perhaps put in header? //DX 20190724 - changed from 0.25 to 0.1

    // ---------------------------------------------------------------------------
    // check for element for center first (fast)
    if(same_species){
      if(env_1.center_element!=env_2.center_element){ return false; }
    }

    // ---------------------------------------------------------------------------
    // check neighboring sites
    for(uint i=0;i<env_1.neighbor_elements.size();i++){
      bool match_found = false;
      vector<string> species;
      for(uint j=0;j<env_2.neighbor_elements.size();j++){
        // ---------------------------------------------------------------------------
        // check same species, if applicable
        if(same_species && env_1.neighbor_elements[i]!=env_2.neighbor_elements[j]){ continue; }
        
        // ---------------------------------------------------------------------------
        // check frequency of distance
        // this is sensitive to tolerance of cutoff; use with caution
        // DX - THIS IS TOO SENSITIVE - if(compare_frequency && env_1.neighbor_frequencies[i]!=env_2.neighbor_frequencies[j]){ continue; }
        
        // ---------------------------------------------------------------------------
        // exact match 
        if(exact_match && aurostd::abs(env_1.neighbor_distances[i]-env_2.neighbor_distances[j])<_TOL_EXACT_MATCH_){
           match_found = true; species.push_back(env_2.neighbor_elements[j]);
        }

        // ---------------------------------------------------------------------------
        // relative match 
        else if(!exact_match && 
                //aurostd::abs(env_1.neighbor_distances[i]-env_2.neighbor_distances[j])/(env_1.neighbor_distances[i]+env_2.neighbor_distances[j])<_TOL_RELATIVE_MATCH_){ //DX 20190730 - too strict
                aurostd::abs(env_1.neighbor_distances[i]-env_2.neighbor_distances[j])<_TOL_LOOSE_MATCH_){ //DX 20190730
            match_found = true; species.push_back(env_2.neighbor_elements[j]);
        }
        
      }
      if(!match_found){ return false; }
      matched_species.push_back(species);
    }

    // ---------------------------------------------------------------------------
    // check relationship between LFA environments
    // only check for on one atom center (cheaper)
    // check this if all checks have passed thus far
    // this is sensitive to tolerance of cutoff; use with caution
    if(compare_frequency && same_species){ //same species for now
      vector<vector<double> > angles_sets_1 = getAnglesBetweenMixedSpeciesEnvironments(env_1.neighbor_coordinates);
      vector<vector<double> > angles_sets_2 = getAnglesBetweenMixedSpeciesEnvironments(env_2.neighbor_coordinates);
  
      for(uint i=0;i<angles_sets_1.size();i++){
        // ---------------------------------------------------------------------------
        // use soft-cutoff; if number of coordinates is not equal, then check that 
        // the smallest set matches; a better alternative then the hard-cutoff freqency match
        // case 1) if set 1 < set 2 
        if(angles_sets_1[i].size()<=angles_sets_2[i].size()){
          for(uint j=0;j<angles_sets_1[i].size();j++){
            bool matched=false;
            for(uint k=0;k<angles_sets_2[i].size();k++){
              //if(aurostd::isequal(angles_sets_1[i][j],angles_sets_2[i][k],10.0)){ //equal within 10 degrees
              if(aurostd::isequal(angles_sets_1[i][j],angles_sets_2[i][k],20.0)){ //equal within 10 degrees
                matched=true;
                break;
              }
              if(VERBOSE){cerr << function_name << " angles dont match: " << angles_sets_1[i][j] << " vs " << angles_sets_2[i][k] << " | diff=" << angles_sets_1[i][j]-angles_sets_2[i][k] << endl;}
            }
            if(!matched){ matched_species.clear(); return false; }
          }
        }
        // ---------------------------------------------------------------------------
        // case 2) if set 1 > set 2 
        else if(angles_sets_1[i].size()>angles_sets_2[i].size()){
          for(uint j=0;j<angles_sets_2[i].size();j++){
            bool matched=false;
            for(uint k=0;k<angles_sets_1[i].size();k++){
              if(aurostd::isequal(angles_sets_2[i][j],angles_sets_1[i][k],10.0)){ //equal within 10 degrees
                matched=true;
                break;
              }
            }
            if(!matched){ matched_species.clear(); return false; }
          }
        }
      }
    }
  
    if(LDEBUG){ cerr << function_name << " environments are compatible." << endl; }

    return true;
  }
}

// ***************************************************************************
//  
// ***************************************************************************
namespace compare{
  vector<vector<double> > getAnglesBetweenMixedSpeciesEnvironments(const vector<vector<xvector<double> > >& neighbor_coordinates){
      
    vector<vector<double> > angles_sets;
    for(uint j=1;j<neighbor_coordinates.size();j++){
      vector<double> angles;
      for(uint i=0;i<neighbor_coordinates[0].size();i++){
        for(uint k=0;k<neighbor_coordinates[j].size();k++){
          angles.push_back(aurostd::angle(neighbor_coordinates[0][i],neighbor_coordinates[j][k])*rad2deg);
        }
      }
      angles_sets.push_back(angles);
    }
    return angles_sets;
  }
}

// ***************************************************************************
// compatible nearest neighbor environments
// ***************************************************************************
namespace compare{
  bool compatibleNearestNeighborTypesEnvironments(const vector<vector<double> >& nn_lfa_with_types_1,
      const vector<vector<double> >& nn_lfa_with_types_2,
      int type_match){

    // Determine if LFA nearest neighbor distances are compatible
    // i.e., quick filter for determining if we can match atoms
    // hinges on alphabetic, perhaps make more robust

    bool VERBOSE=false;

    if(VERBOSE){
      for(uint i=0;i<nn_lfa_with_types_1.size();i++){
        cerr <<"1 " << i << ": " << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(nn_lfa_with_types_1[i]),",") << endl;
      } 
      for(uint i=0;i<nn_lfa_with_types_2.size();i++){
        cerr <<"2 " << i << ": " << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(nn_lfa_with_types_2[i]),",") << endl;
      } 
    }

    // ---------------------------------------------------------------------------
    // if same species comparison
    if(type_match==2){
      if(nn_lfa_with_types_1.size()==nn_lfa_with_types_2.size()){
        for(uint i=0;i<nn_lfa_with_types_1.size();i++){
          bool matched_set=false;
          for(uint j=0;j<nn_lfa_with_types_2.size();j++){
            bool matched = true;
            for(uint ii=0;ii<nn_lfa_with_types_1[i].size();ii++){
              if(aurostd::abs(nn_lfa_with_types_1[i][ii]-nn_lfa_with_types_2[j][ii])/(nn_lfa_with_types_1[i][ii]+nn_lfa_with_types_2[j][ii])>0.1){ // less than 10% relative error
                matched=false;
                break;
              }
            }
            if(matched){
              matched_set=true;
            }
          }
          if(!matched_set){
            return false;
          }
        }
      }
      else{
        // come back here to deal with comparing supercells
        // for now, putting true and forcing comparison
        return true;
      }
    }
    // ---------------------------------------------------------------------------
    // if not same species
    else{
      // come back here later (this is a bit more complicated)
      return true;
    }

    return true;
  }
}

// ***************************************************************************
// get nearest neighbors 
// ***************************************************************************
namespace compare{
  vector<AtomEnvironment> getUniqueTypesAtomEnvironmentForLFA(const xstructure& xstr, const string lfa, 
      const vector<string>& LFAs){

    // Calculates the nearest neighbor distance to a particular atom
    // hinges on alphabetic, perhaps make more robust

    vector<AtomEnvironment> environments_LFA;

    for(uint i=0;i<xstr.atoms.size();i++){
      if(xstr.atoms[i].name == lfa){
        AtomEnvironment LFA_env; 
        LFA_env.center_element = xstr.atoms[i].name;
        LFA_env.center_type = xstr.atoms[i].type;

        for(uint j=0;j<LFAs.size();j++){
        //DX ORIG for(uint j=0;j<xstr.species.size();j++){
          //if(xstr.species[j]!=lfa){
            uint frequency = 0;
            vector<xvector<double> > coordinates;
            //DX ORIG LFA_env.neighbor_elements.push_back(xstr.species[j]);
            LFA_env.neighbor_elements.push_back(LFAs[j]); //TEST 
            LFA_env.neighbor_types.push_back(j);
            //DX ORIG LFA_env.neighbor_distances.push_back(shortestDistanceRestrictType(xstr,i,frequency,coordinates,xstr.species[j]));
            LFA_env.neighbor_distances.push_back(shortestDistanceRestrictType(xstr,i,frequency,coordinates,LFAs[j])); //TEST
            LFA_env.neighbor_frequencies.push_back(frequency);
            LFA_env.neighbor_coordinates.push_back(coordinates);
          //}
        }
        environments_LFA.push_back(LFA_env);
      }
    }
    return environments_LFA;
  }
}

// ***************************************************************************
// Shortest Distance from one atom
// ***************************************************************************
namespace compare{
  double shortestDistanceRestrictType(const xstructure& xstr, const uint& k, string type) {
    // Find the minimum interatomic distance from a central atom to 
    // a particular element/type 
    // (overload)
    
    uint frequency = 0;
    vector<xvector<double> > coordinates;
    return shortestDistanceRestrictType(xstr,k,frequency,coordinates,type);
  }
}

namespace compare{
  double shortestDistanceRestrictType(const xstructure& xstr, const uint& k, 
      uint& frequency, vector<xvector<double> >& coordinates, string type) {

    // Find the minimum interatomic distance from a central atom to 
    // a particular element/type and store frequency/coordination and coordinates

    // ---------------------------------------------------------------------------
    // instantiate variables
    double min_dist=1e9;
    minimumCoordinationShell(xstr, k, min_dist, frequency, coordinates, type);

    return min_dist;
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

    // Find the minimum interatomic distance in the structure to atom k
    // (perhaps integrate with SYM::minimumDistance())

    double min_dist=1e9;
    double prev_min_dist=0; //DX 20190716
    xmatrix<double> lattice = xstr.lattice; //NEW

    //DX speed increase
    //perhaps can speed up even more, since the lattice doesn't change for the xstr...
    vector<xvector<double> > l1, l2, l3;
    vector<int> a_index, b_index, c_index;
    xvector<int> dims(3); //DX 20190710 - use robust method
    dims[1]=dims[2]=dims[3]=0; //reset
    
    xvector<double> tmp;

    for(uint ii=0; ii<xstr.atoms.size(); ii++){
      if(ii!=k){
        if(min_dist<prev_min_dist){
          if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
            resetLatticeDimensions(lattice,min_dist,dims,l1,l2,l3,a_index,b_index,c_index);
            prev_min_dist=min_dist;
          }
        }
        xvector<double> incell_dist = xstr.atoms[k].cpos-xstr.atoms[ii].cpos;
        double incell_mod = aurostd::modulus(incell_dist);
        if(incell_mod<min_dist){
          if(!(dims[1]==1 && dims[2]==1 && dims[3]==1)){
            resetLatticeDimensions(lattice,incell_mod,dims,l1,l2,l3,a_index,b_index,c_index);
          }
          prev_min_dist=incell_mod;
        }
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
// Magnetic Deviation (beta)
// ***************************************************************************
namespace compare{
  void magneticDeviation(const xstructure& xstr1, const xstructure& xstr2, 
		const vector<uint>& indexMatch1, const vector<uint>& indexMatch2,
    double& magnetic_deviation, double& magnetic_fail) {

    // BETA functionality: Determine magnetic deviation between magnetic structures

    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "compare::magneticDeviation()";
    double projection_sum = 0.0;
    double mag_fail = 0.0;
    bool is_non_collinear = xstr1.atoms[0].noncoll_spin_is_given;

    for(uint j=0; j<indexMatch1.size(); j++){
      
      double magnetic_projection_str1 = 0.0;
      double magnetic_projection_str2 = 0.0;

      // ---------------------------------------------------------------------------
      // collinear 
      if(!is_non_collinear){
        //if(aurostd::isequal(xstr1.atoms[indexMatch1[j]].spin,0.0,1e-6)){
        //  projection_sum = projection_sum + (
        //}
        double dot_product = xstr1.atoms[indexMatch1[j]].spin*xstr2.atoms[indexMatch2[j]].spin;
        double cos_theta = dot_product/(aurostd::abs(xstr1.atoms[indexMatch1[j]].spin)*aurostd::abs(xstr2.atoms[indexMatch2[j]].spin));
        magnetic_projection_str1 = 1.0-(aurostd::abs(xstr2.atoms[indexMatch2[j]].spin)*cos_theta/aurostd::abs(xstr1.atoms[indexMatch1[j]].spin));        
        magnetic_projection_str2 = 1.0-(aurostd::abs(xstr1.atoms[indexMatch1[j]].spin)*cos_theta/aurostd::abs(xstr2.atoms[indexMatch2[j]].spin));        
        if(LDEBUG){
          cerr << function_name << " matching xstr1 (" << indexMatch1[j] << ") mag=" << xstr1.atoms[indexMatch1[j]].spin << " to xstr2 (" << indexMatch2[j] << ") mag=" << xstr2.atoms[indexMatch2[j]].spin << endl;
          cerr << function_name << " dot_product: " << xstr1.atoms[indexMatch1[j]].spin*xstr2.atoms[indexMatch2[j]].spin << endl; 
          cerr << function_name << " cos_theta: " << cos_theta << endl; 
          cerr << function_name << " magnetic_projection_str1: " << magnetic_projection_str1 << endl;
          cerr << function_name << " magnetic_projection_str2: " << magnetic_projection_str2 << endl;
        }
        if(aurostd::isequal(cos_theta,-1.0,1e-6)){
          mag_fail += 1.0; //str1
          mag_fail += 1.0; //str2
        }
        else if(aurostd::isequal(dot_product,0.0,1e-6)){
          projection_sum = projection_sum + 0.0;
        }
        else{
          projection_sum = projection_sum + (1.0-magnetic_projection_str1) + (1.0-magnetic_projection_str2);
        }

      }
      // ---------------------------------------------------------------------------
      // non-collinear 
      else if(is_non_collinear){
        double projection = aurostd::abs(aurostd::scalar_product(xstr1.atoms[indexMatch1[j]].noncoll_spin,xstr2.atoms[indexMatch2[j]].noncoll_spin));
        double mod_spin_1 = aurostd::modulus(xstr1.atoms[indexMatch1[j]].noncoll_spin);
        double mod_spin_2 = aurostd::modulus(xstr2.atoms[indexMatch2[j]].noncoll_spin);
        magnetic_projection_str1 = projection/(mod_spin_1*mod_spin_2);
        magnetic_projection_str2 = projection/(mod_spin_1*mod_spin_2);
        projection_sum = projection_sum + (1.0-magnetic_projection_str1) + (1.0-magnetic_projection_str2);
      
      }
      if(LDEBUG){
        cerr << "projection_sum: " << projection_sum << endl;
      }
    }

    magnetic_deviation = (projection_sum/(xstr1.atoms.size()+xstr2.atoms.size()-mag_fail));
    magnetic_fail = (mag_fail/(xstr1.atoms.size()+xstr2.atoms.size()));
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
// Compute Magnetic Misfit (beta)
// ***************************************************************************
namespace compare{
  double computeMagneticMisfit(const double dev, const double dis, const double fail, const double mag_dis, const double mag_fail) {

    // Combines differences between all aspects of crystal structure into a 
    // figure of misfit. (See Burzlaff)
    // It also adds in a new magnetic contribution (created by us), in the 
    // same spirit as the Burzlaff criteria
    // BETA FUNCTIONALITY

    double mis;

    mis=1-((1-dev)*(1-dis)*(1-fail)*(1-mag_dis)*(1-mag_fail));
    return mis;
  }   
}

// ***************************************************************************
// Print Matching Between Atoms
// ***************************************************************************
namespace compare{
  void printMatch(const vector<uint>& indexMatch1, const vector<uint>& indexMatch2, 
      const vector<double>& distances,
      const xstructure& PROTO, const xstructure& xstr1, ostream& oss) { //DX 20190802 - added distances

    // With this function we print the atoms matched in the previous function
    // whose indices are stored in the indexMatch vector.
    // This allows us to call them directly.
    // Lastly, the translational term, specific of each mapping, is printed.

    uint i,j;

    oss << "              Reference                               Mapped                               Distances"<<endl;
    for(j=0; j< indexMatch1.size(); j++){
      oss << indexMatch1[j]<<"-"<<indexMatch2[j]<<"    ";
      oss << xstr1.atoms[indexMatch1[j]].cpos << " " << xstr1.atoms[indexMatch1[j]].name;
      oss << "       ";
      oss << PROTO.atoms[indexMatch2[j]].cpos << " " << PROTO.atoms[indexMatch2[j]].name;
      oss << "       ";
      oss << distances[j] << endl; //DX 20190802
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
// atomInCell() 
// ***************************************************************************
namespace compare{
  bool atomInCell(const _atom& atom){ 

    // check if the atom is in the unit cell based on fractional coordinates
    // this alone is not robust; this should be used in tandem with 
    // SYM::MapAtom() to account for periodic boundary conditions
    // filtering with this first with soft cutoffs is much faster, 
    // especially if there are many atoms to check (e.g., 20,000)
    
    for(uint f=1;f<4;f++){
      if(atom.fpos[f] > 1.05 || atom.fpos[f] < -0.05){ //soft cutoff, use hard cutoff later
        return false;
      }
    }
    return true;
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
    //DX 20190619 xmatrix<double> c2f = inverse(trasp(lattice));
    bool skew = false;

    //vector<int> ind(2); ind[0]=i, ind[1]=j;
    //xvector<double> tmp;

    uint count=0;
    //deque<_atom> transformed;
    //deque<uint> index_to_check;

    // ===== Check if applying the symmetry element along with internal translation maps to another atom ===== //
    for(uint d=0;d<atoms.size();d++){
      //DX 20190702 - no name information: _atom tmp;
      //DX 20190702 - use assignment op to get type - tmp.type = atoms[d].type;
      _atom tmp = atoms[d]; //copy names, types, etc. //DX 20190702
      tmp.cpos = atoms[d].cpos+vec;
      tmp.fpos = C2F(lattice,tmp.cpos);
      if(SYM::MapAtom(atoms,tmp,true,lattice,f2c,skew,tolerance)){ //DX 20190619 - removed c2f
        //transformed.push_back(tmp);
        //index_to_check.push_back(d);
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
// GetLFASupercell
// ***************************************************************************
namespace compare{
  xstructure GetLFASupercell(const xstructure& xstr, const xvector<int>& dims, const string& lfa_name){

    // build a supercell comprised only of LFA atoms
    // to speed up translation vector search

    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "compare::GetLFASupercell()";

    // ---------------------------------------------------------------------------
    // remove all atoms that are not of the LFA type 
    xstructure xstr_LFA_only=xstr;
    xstr_LFA_only.ClearSymmetry(); //DX 20181022
    uint num_atoms=xstr_LFA_only.atoms.size();
    for(uint i=0;i<num_atoms;i++){
	    if(xstr_LFA_only.atoms[i].name!=lfa_name){
	      xstr_LFA_only.RemoveAtom(i);
	      num_atoms--;
	      i--;
	    }
    }

    // ---------------------------------------------------------------------------
    // shift an LFA atom to origin
    for(uint a=0;a<xstr_LFA_only.atoms.size();a++){
      if(xstr_LFA_only.atoms[a].name == lfa_name){
        xstr_LFA_only.ShifOriginToAtom(a);
        break;
      }
    }

    /*
    // ---------------------------------------------------------------------------
    // create supercell (fast) 
    vector<int> sc2pcMap, pc2scMap;
    bool get_symmetry=false;
    bool get_full_basis=false;
    bool force_supercell_matrix=true;
    xstructure xstr_LFA_supercell=GetSuperCell(xstr_LFA_only,3,0,0,0,3,0,0,0,3,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix); //DX 20190319 - use supercell matrix in expansion
    */

    // ---------------------------------------------------------------------------
    // create supercell (fast/robust)
    xstructure xstr_LFA_supercell=xstr_LFA_only;
    GenerateGridAtoms(xstr_LFA_supercell,dims(1),dims(2),dims(3));

    // ---------------------------------------------------------------------------
    // update atoms
    xstr_LFA_supercell.atoms = xstr_LFA_supercell.grid_atoms;
    xstr_LFA_supercell.grid_atoms.clear();
	  //xstr_LFA_supercell = pflow::SetNumEachType(xstr_LFA_supercell, sizes);

	  if(LDEBUG){cerr << function_name << ": Number of LFAs in supercell: " << xstr_LFA_supercell.atoms.size() << endl;}
      
    return xstr_LFA_supercell;
  }
}

// ***************************************************************************
// latticeAndOriginSearch
// ***************************************************************************
namespace compare{
  void latticeAndOriginSearch(xstructure& xstr1, xstructure& xstr2, 
    const uint& num_proc,xmatrix<double>& q1, vector<xstructure> &vprotos, 
    double& minMis, int type_match, bool optimize_match, ostream& oss){

    // Performs lattice and origin search

    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "compare::latticeAndOriginSearch()";

    bool test_one_lfa_only = false; //DX 20190318
    //DX - SPEED UP BUT NOT ROBUST - if(type_match==2){ test_one_lfa_only=true;} //DX 20190318

    bool magnetic_analysis = (xstr1.atoms[0].spin_is_given || xstr1.atoms[0].noncoll_spin_is_given);
    double min_latt_dev=0.0; double min_coordinate_dis=0.0; double min_failure=0.0;
    double min_magnetic_dis=0.0; double min_magnetic_failure=0.0;
    vector<uint> matching_indices_1, matching_indices_2;
    vector<double> minimum_distances;

    // ---------------------------------------------------------------------------
    // determine least-frequently occuring atom type (LFA) for each structure
    // (there may be more than one)
    vector<string> LFA_str1=getLeastFrequentAtomSpecies(xstr1);
    vector<string> LFA_str2=getLeastFrequentAtomSpecies(xstr2);
    string lfa_str1=LFA_str1[0]; //initialize
    string lfa_str2=LFA_str2[0]; //initialize

    double abs_det_q1=abs(det(q1)); // volume
    xvector<double> abc_angles_q1=Getabc_angles(q1,DEGREES); // lattice parameters
   
    // ---------------------------------------------------------------------------
    // determine supercell size via search radius/dims
    double search_radius = aurostd::max(abc_angles_q1(1),abc_angles_q1(2),abc_angles_q1(3));
    xvector<int> dims = LatticeDimensionSphere(xstr2.lattice,search_radius);

    if(LDEBUG){cerr << function_name << ": lattice search radius: " << search_radius << endl;}
    if(LDEBUG){cerr << function_name << ": lattice dims : " << dims << endl;}

    // ---------------------------------------------------------------------------
    // peform supercell expansion on LFA atoms in structure2 
    xstructure xstr_LFA_supercell = compare::GetLFASupercell(xstr2, dims, lfa_str2);
   
    // ---------------------------------------------------------------------------
    // find possible translation vectors 
    vector<xvector<double> > translation_vectors;
    vector<vector<uint> > ij_index;
	  quadrupletSearch(q1,xstr_LFA_supercell,xstr2,translation_vectors,ij_index); //DX 20190701 - added xstr2

    // ---------------------------------------------------------------------------
    // build possible lattices
    vector<xmatrix<double> > lattices;
    vector<xmatrix<double> > clattices;
    vector<double> latt_devs;
    buildSimilarLattices(translation_vectors, q1, abs_det_q1, abs_det_q1, abc_angles_q1, lattices, clattices, latt_devs, optimize_match);
    if(LDEBUG){cerr << function_name << ": Number of lattices to compare: " << lattices.size() << endl;}
      
    if(lattices.size()>0){
    
      xstructure xstr1_tmp = xstr1;
      // ---------------------------------------------------------------------------
      // calculate attributes of structure 1 (volume, lattice parameters, nearest neighbor distances, etc.)
      vector<double> D1,F1;
      cellDiagonal(xstr1_tmp,D1,F1,1); // cell diagonals
      // convert to clattice representation
      xstr1_tmp.lattice=GetClat(xstr1_tmp.a,xstr1_tmp.b,xstr1_tmp.c,xstr1_tmp.alpha,xstr1_tmp.beta,xstr1_tmp.gamma);
      for(uint iat=0; iat<xstr1_tmp.atoms.size(); iat++){
        xstr1_tmp.atoms[iat].cpos=F2C(xstr1_tmp.lattice,xstr1_tmp.atoms[iat].fpos);
      }
      vector<double> all_nn1 = computeNearestNeighbors(xstr1_tmp); // nearest neighbor distances (invariant of origin shifts) 
      
      // ---------------------------------------------------------------------------
      // peform expansion on structure2
      // wait until we confirm there are similar lattices, otherwise we build it for nothing (i.e. unnecessary cost)
      xstructure xstr_supercell=xstr2;
      GenerateGridAtoms(xstr_supercell,dims(1),dims(2),dims(3));

      // ---------------------------------------------------------------------------
      // update atoms
      xstr_supercell.atoms = xstr_supercell.grid_atoms;
      xstr_supercell.grid_atoms.clear();

#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
      // ---------------------------------------------------------------------------
      // split task into threads 
      uint number_of_lattices = lattices.size();
      uint number_of_threads = aurostd::min(num_proc,number_of_lattices); // cannot have more threads than lattices
      //DX 20191107 [switching to getThreadDistribution] - vector<uint> start_indices, end_indices;
      //DX 20191107 [switching to getThreadDistribution] - splitTaskIntoThreads(number_of_lattices, number_of_threads, start_indices, end_indices);
      vector<vector<int> > thread_distribution = getThreadDistribution(number_of_lattices, number_of_threads); //DX 20191107 
#endif

      /*
      // ---------------------------------------------------------------------------
      // identify rotation between original and new structure
      // [ODD BEHAVIOR - REVISIT]
      for(uint p=0;p<lattices.size();p++){
        cerr << "lattice:" << lattices[p] << endl;
        cerr << "clattice:" << clattices[p] << endl;
        xmatrix<double> lattice_metric_tensor = MetricTensor(lattices[p]);
        xmatrix<double> clattice_metric_tensor = MetricTensor(clattices[p]);
        cerr << "lattice metric tensor: " << lattice_metric_tensor << endl;
        cerr << "clattice metric tensor: " << clattice_metric_tensor << endl;
        //xmatrix<double> rotation = lattices[p]*aurostd::inverse(clattices[p]);
        xmatrix<double> rotation1 = aurostd::inverse(lattices[p])*clattices[p];
        cerr << "rotation1: " << rotation1 << endl; 
        xstructure xstr_rot1 = xstr2;
        xstr_rot1 = Rotate(xstr2,trasp(rotation1));
        xmatrix<double> rotation2 = trasp(clattices[p])*trasp(aurostd::inverse(lattices[p]));
        cerr << "rotation2: " << rotation2 << endl; 
        xstructure xstr_rot2 = xstr2;
        xstr_rot2 = Rotate(xstr2,rotation2);
        xmatrix<double> rotation3 = trasp(clattices[p])*trasp(aurostd::inverse(lattices[p]));
        cerr << "rotation3: " << rotation3 << endl; 
        xstructure xstr_rot3 = xstr2;
        xstr_rot3 = GetLTCell(rotation3,xstr2);
        xmatrix<double> rotation4 = trasp(rotation1); 
        cerr << "rotation4: " << rotation4 << endl; 
        xstructure xstr_rot4 = xstr2;
        xstr_rot4 = Rotate(xstr2,rotation4);

        cerr << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
        cerr << "orig xstr2: " << xstr2 << endl;
        cerr << "-----------------------------------------------------------------" << endl;
        cerr << "rotated xstr2 - 1: " << xstr_rot1 << endl;
        cerr << "-----------------------------------------------------------------" << endl;
        cerr << "rotated xstr2 - 2: " << xstr_rot2 << endl;
        cerr << "-----------------------------------------------------------------" << endl;
        cerr << "rotated xstr2 - 3: " << xstr_rot3 << endl;
        cerr << "-----------------------------------------------------------------" << endl;
        cerr << "rotated xstr2 - 4: " << xstr_rot4 << endl;
      }
      */

      // ---------------------------------------------------------------------------
      // test origin shifts 
      for(uint y=0;y<LFA_str2.size();y++){
        for(uint x=0;x<LFA_str1.size();x++){
          lfa_str1=LFA_str1[x];
          lfa_str2=LFA_str2[y];
          if(type_match == 2 && lfa_str1 != lfa_str2){ continue;}

          oss << "===> LFA (structure 1): " << lfa_str1 << endl;
          oss << "===> LFA (structure 2): " << lfa_str2 << endl;

          if(LDEBUG){
            cerr << function_name << ": LFA (structure 1): " << lfa_str1 << endl;
            cerr << function_name << ": LFA (structure 2): " << lfa_str2 << endl;
          }
        
          // ---------------------------------------------------------------------------
          // shift supercell to LFA
          // cerr << "xstr_supercell.atoms.size(): " << xstr_supercell.atoms.size() << endl;
          // for(uint a=0;a<xstr_supercell.atoms.size();a++){
          //   if(xstr_supercell.atoms[a].name == lfa_str2){
          //     xstr_supercell.ShifOriginToAtom(a);
          //     break;
          //   }
          // }

          // ---------------------------------------------------------------------------
          // shift representative structure to LFA
          // NEED TO SHIFT origin of xstr1_tmp to one of the LFA (this was missing before and caused ICSD_102428.BCA, and CBA to not match, but they should
          for(uint i=0;i<xstr1_tmp.atoms.size();i++){
            if(xstr1_tmp.atoms[i].name==lfa_str1){
              xstr1_tmp.ShifOriginToAtom(i);
              xstr1_tmp.BringInCell(1e-10);
              break;
            }
          }

          // ---------------------------------------------------------------------------
          // create vector of variables for each thread 
          vector<xstructure> xstr1_for_thread;
          vector<double> possible_minMis, possible_minLattDev, possible_minCoordDis, possible_minFail, possible_minMagneticDis, possible_minMagneticFail;
          vector<vector<uint> > possible_matching_indices_1, possible_matching_indices_2;
          vector<vector<double> > possible_minimum_distances;
          vector<vector<xstructure> > vvprotos;
          for(uint n=0; n<num_proc; n++){
            vector<xstructure> vprotos_tmp;
            vvprotos.push_back(vprotos_tmp);
            xstr1_for_thread.push_back(xstr1_tmp);
            possible_minMis.push_back(1.0);
            possible_minLattDev.push_back(0.0);
            possible_minCoordDis.push_back(0.0);
            possible_minFail.push_back(0.0);
            possible_minMagneticDis.push_back(0.0);
            possible_minMagneticFail.push_back(0.0);
            vector<uint> tmp_indices;
            possible_matching_indices_1.push_back(tmp_indices);
            possible_matching_indices_2.push_back(tmp_indices);
            vector<double> tmp_distances;
            possible_minimum_distances.push_back(tmp_distances);
          }
 
#ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
          // ---------------------------------------------------------------------------
          // threaded (DX 20191107 thread pointer) 
          vector<std::thread*> threads;
          if(LDEBUG){cerr << function_name << ": Searching for possible matching structures [THREADED VERSION]" << endl;}
          for(uint n=0; n<number_of_threads; n++){
            threads.push_back(new std::thread(structureSearch,
                  std::ref(xstr1_for_thread[n]),
                  std::ref(xstr_supercell),
                  std::ref(all_nn1),
                  std::ref(lfa_str2),
                  type_match,
                  std::ref(lattices),std::ref(clattices),std::ref(latt_devs),
                  //DX 20191107 [switching to getThreadDistribution] - start_indices[n], end_indices[n],
                  thread_distribution[n][0], thread_distribution[n][1],
                  std::ref(possible_minMis[n]),std::ref(possible_minLattDev[n]),
                  std::ref(possible_minCoordDis[n]),std::ref(possible_minFail[n]),
                  std::ref(possible_minMagneticDis[n]),std::ref(possible_minMagneticFail[n]),
                  std::ref(possible_matching_indices_1[n]),std::ref(possible_matching_indices_2[n]),
                  std::ref(possible_minimum_distances[n]),std::ref(vvprotos[n]),
                  optimize_match));
          }         
          for(uint t=0;t<threads.size();t++){
            threads[t]->join();
            delete threads[t];
          }
#else
          // ---------------------------------------------------------------------------
          // non-threaded 
          uint n=0;
          uint start_index=0;
          uint end_index=lattices.size();  //DX 20191107 switching end point convention
          if(LDEBUG){cerr << function_name << ": Searching for possible matching structures [NON-THREADED VERSION]" << endl;}
          //structureSearch(lfa_str2,all_nn1,xstr_supercell,vvprotos[n],xstr1_for_thread[n],type_match,possible_minMis[n],
          //                lattices,clattices,latt_devs,optimize_match,start_index,end_index);
          structureSearch(
              xstr1_for_thread[n],
              xstr_supercell,
              all_nn1,
              lfa_str2,
              type_match,
              lattices,clattices,latt_devs,
              start_index, end_index,
              possible_minMis[n],possible_minLattDev[n],
              possible_minCoordDis[n],possible_minFail[n],
              possible_minMagneticDis[n],possible_minMagneticFail[n],
              possible_matching_indices_1[n],possible_matching_indices_2[n],
              possible_minimum_distances[n],vvprotos[n],
              optimize_match);
#endif
      
          // ---------------------------------------------------------------------------
          // collect misfits and matching structure representations
          for(uint p=0;p<possible_minMis.size();p++){
            if(p==0 && y==0 && x==0){ // DX 2/8/17 - need to add x==0 ortherwise matches can be overwritten
              minMis=possible_minMis[p];
              min_latt_dev=possible_minLattDev[p];
              min_coordinate_dis=possible_minCoordDis[p];
              min_failure=possible_minFail[p];
              min_magnetic_dis=possible_minMagneticDis[p];
              min_magnetic_failure=possible_minMagneticFail[p];
              matching_indices_1=possible_matching_indices_1[p];
              matching_indices_2=possible_matching_indices_2[p];
              minimum_distances=possible_minimum_distances[p];
              xstr1=xstr1_for_thread[p];
              vprotos=vvprotos[p];
            }
            else {
              if(possible_minMis[p]<=minMis){
                minMis=possible_minMis[p];
                min_latt_dev=possible_minLattDev[p];
                min_coordinate_dis=possible_minCoordDis[p];
                min_failure=possible_minFail[p];
                min_magnetic_dis=possible_minMagneticDis[p];
                min_magnetic_failure=possible_minMagneticFail[p];
                matching_indices_1=possible_matching_indices_1[p];
                matching_indices_2=possible_matching_indices_2[p];
                minimum_distances=possible_minimum_distances[p];
                xstr1=xstr1_for_thread[p];
                vprotos=vvprotos[p];
              }
            }
          }

          // ---------------------------------------------------------------------------
          // quick return if found a match
          if(minMis<0.1 && !optimize_match){
            if(LDEBUG){cerr << function_name << ": Found match (misfit = " << minMis << ")! Terminating search early." << endl;}
            printStructureMappingResults(oss,xstr1,vprotos[0],minMis,min_latt_dev,min_coordinate_dis,min_failure,min_magnetic_dis,min_magnetic_failure,
                matching_indices_1,matching_indices_2,minimum_distances,magnetic_analysis);
            return;
          }

          // ---------------------------------------------------------------------------
          // quick return if testing only one LFA set
          //DX 20190702 - can i do this: if(!optimize_match && minMis==1){ test_one_lfa_only=true;}
          if(!optimize_match && minMis==1 && type_match==2){ test_one_lfa_only=true;} //DX 20190809 - need type match here; otherwise we may miss structure-type matches
          if(test_one_lfa_only){
            if(LDEBUG){cerr << function_name << ": No match found. Searched only one LFA set. Terminating search early." << endl;}
            return;
          }
        } 
      } 
    }  
    if(minMis!=1.0 && vprotos.size()>0){
      printStructureMappingResults(oss,xstr1,vprotos[0],minMis,min_latt_dev,min_coordinate_dis,min_failure,min_magnetic_dis,min_magnetic_failure,
          matching_indices_1,matching_indices_2,minimum_distances,magnetic_analysis);
    }
  }
}

// [OBSOLETE - DX 20190717] // ***************************************************************************
// [OBSOLETE - DX 20190717] // Thread Generation (For parallel processing of quadruplets)
// [OBSOLETE - DX 20190717] // ***************************************************************************
// [OBSOLETE - DX 20190717] namespace compare{
// [OBSOLETE - DX 20190717]   void threadGeneration(const uint& num_proc,xmatrix<double>& q1, xstructure& xstr2, 
// [OBSOLETE - DX 20190717] 			vector<xstructure> &vprotos, xstructure &xstr1, const int& type_match, 
// [OBSOLETE - DX 20190717] 			const bool& optimize_match, double& minMis, ostream& oss){ 
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]     // This function creates the supercell of the second structure and begins 
// [OBSOLETE - DX 20190717]     // the quadruplets search within a supercell. Due to the costly nature of 
// [OBSOLETE - DX 20190717]     // this algorithm, the quadruplet search is parallelized. The splitting 
// [OBSOLETE - DX 20190717]     // of computation of done here
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]     bool LDEBUG=(false || XHOST.DEBUG);
// [OBSOLETE - DX 20190717]     bool test_one_lfa_only = false; //DX 20190318
// [OBSOLETE - DX 20190717]     //DX if(type_match==2){ test_one_lfa_only=true;} //DX 20190318 - need to comment out for permutation matching
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]     xstructure xstr1_tmp = xstr1;
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]     //cerr << "LFA" << endl;
// [OBSOLETE - DX 20190717]     vector<string> LFA_str1=getLeastFrequentAtomSpecies(xstr1);
// [OBSOLETE - DX 20190717]     vector<string> LFA_str2=getLeastFrequentAtomSpecies(xstr2);
// [OBSOLETE - DX 20190717]     string lfa, lfa_str1;
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]     //cerr << "SUPERCELL" << endl;
// [OBSOLETE - DX 20190717]     //DX 20190319 - START
// [OBSOLETE - DX 20190717]     vector<int> sc2pcMap, pc2scMap;
// [OBSOLETE - DX 20190717]     bool get_symmetry=false;
// [OBSOLETE - DX 20190717]     bool get_full_basis=false;
// [OBSOLETE - DX 20190717]     bool force_supercell_matrix=true;
// [OBSOLETE - DX 20190717]     xstructure xstr=GetSuperCell(xstr2,3,0,0,0,3,0,0,0,3,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix); //DX 20190319 - use supercell matrix in expansion
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]     uint y=0;
// [OBSOLETE - DX 20190717]     uint x=0;
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]     // DX TEST
// [OBSOLETE - DX 20190717]     // Consider all LFAs
// [OBSOLETE - DX 20190717]     for(y=0;y<LFA_str2.size();y++){
// [OBSOLETE - DX 20190717]       for(x=0;x<LFA_str1.size();x++){
// [OBSOLETE - DX 20190717]         // DX TEST
// [OBSOLETE - DX 20190717]         lfa_str1=LFA_str1[x];
// [OBSOLETE - DX 20190717]         lfa=LFA_str2[y];
// [OBSOLETE - DX 20190717]         if(type_match == 2 && lfa_str1 != lfa){
// [OBSOLETE - DX 20190717]           continue;
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717]         oss << "===> LFA: "<<lfa<<endl;
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         if(LDEBUG) {
// [OBSOLETE - DX 20190717]           cerr << "===> LFA_1: " << lfa_str1 <<endl;
// [OBSOLETE - DX 20190717]           cerr << "===> LFA: "<<lfa<<endl;
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         // DX NEW so we don't need to do this in an inner loop
// [OBSOLETE - DX 20190717]         for(uint a=0;a<xstr.atoms.size();a++){
// [OBSOLETE - DX 20190717]           if(xstr.atoms[a].name == lfa){
// [OBSOLETE - DX 20190717]             xstr.ShifOriginToAtom(a);
// [OBSOLETE - DX 20190717]             break;
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717]         // DX NEW
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         //cerr << "xstr1 centroid: " << endl;
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         //cerr << "SHIFT" << endl;
// [OBSOLETE - DX 20190717]         // NEED TO SHIFT origin of xstr1_tmp to one of the LFA (this was missing before and caused ICSD_102428.BCA, and CBA to not match, but they should
// [OBSOLETE - DX 20190717]         for(uint i=0;i<xstr1_tmp.atoms.size();i++){
// [OBSOLETE - DX 20190717]           if(xstr1_tmp.atoms[i].name==lfa_str1){
// [OBSOLETE - DX 20190717]             xstr1_tmp.ShifOriginToAtom(i);
// [OBSOLETE - DX 20190717]             xstr1_tmp.BringInCell(1e-10);
// [OBSOLETE - DX 20190717]             break;
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         // NEED TO SHIFT origin of xstr2 to one of the LFA
// [OBSOLETE - DX 20190717]         for(uint i=0;i<xstr2.atoms.size();i++){
// [OBSOLETE - DX 20190717]           if(xstr2.atoms[i].name==lfa){
// [OBSOLETE - DX 20190717]             xstr2.ShifOriginToAtom(i);
// [OBSOLETE - DX 20190717]             xstr2.BringInCell(1e-10);
// [OBSOLETE - DX 20190717]             break;
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         // When checking the quadruplets/lattice, we only need to generate a 
// [OBSOLETE - DX 20190717]         // LFA supercell (i.e. take out the other atoms).  This greatly reduces 
// [OBSOLETE - DX 20190717]         // the time of computation (don't need to scan through unnecessary atoms) 
// [OBSOLETE - DX 20190717]         xstructure xstr_LFA_only=xstr2;
// [OBSOLETE - DX 20190717]         xstr_LFA_only.ClearSymmetry(); //DX 20181022
// [OBSOLETE - DX 20190717]         uint num_atoms=xstr_LFA_only.atoms.size();
// [OBSOLETE - DX 20190717]         for(uint i=0;i<num_atoms;i++){
// [OBSOLETE - DX 20190717]           if(xstr_LFA_only.atoms[i].name!=lfa){
// [OBSOLETE - DX 20190717]             xstr_LFA_only.RemoveAtom(i);
// [OBSOLETE - DX 20190717]             num_atoms--;
// [OBSOLETE - DX 20190717]             i--;
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717]         //cerr << "LFA SUPERELL" << endl;
// [OBSOLETE - DX 20190717]         xstructure xstr_LFA_supercell=GetSuperCell(xstr_LFA_only,3,0,0,0,3,0,0,0,3,sc2pcMap,pc2scMap,get_symmetry,get_full_basis,force_supercell_matrix); //DX 20190319 - use supercell matrix in expansion
// [OBSOLETE - DX 20190717]         //cerr << "created LFA supercell" << endl;
// [OBSOLETE - DX 20190717]         // Determines the number of LFAs in the supercell.
// [OBSOLETE - DX 20190717]         int num_LFAs=-1; //-1 as default value 
// [OBSOLETE - DX 20190717]         for(uint q=0; q<xstr.num_each_type.size();q++){
// [OBSOLETE - DX 20190717]           if(xstr.species[q] == lfa){ 
// [OBSOLETE - DX 20190717]             num_LFAs= xstr.num_each_type[q];
// [OBSOLETE - DX 20190717]             if(LDEBUG) {cerr << "compare:: " << "Number of LFAs in supercell: " << xstr.species[q] << "= " << num_LFAs << endl;}
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         // === THREAD PREPARATION FOR PARALLEL PROCESSING OF QUADRUPLET SEARCH === //
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         // DECLARATION OF ATOMIC BOOL SECTION: Allows the threads to communicate. 
// [OBSOLETE - DX 20190717]         // This is useful for stopping the threads if the misfit falls below
// [OBSOLETE - DX 20190717]         // the compatible misfit criterion (mis<0.1) in any of the threads. 
// [OBSOLETE - DX 20190717]         // [OBSOLETE] std::atomic_bool misfit_in_threshold_found (false);
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         //[NONTHREADS]bool misfit_in_threshold_found=false;
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         vector<xstructure> xstr1_for_thread;
// [OBSOLETE - DX 20190717]         vector<double> possible_minMis;
// [OBSOLETE - DX 20190717]         vector<vector<xstructure> > vvprotos;
// [OBSOLETE - DX 20190717]         //vector<std::thread> threads;
// [OBSOLETE - DX 20190717]         //DX NEW - used to be done in one of the inner loops in structureSearch 
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         // compute xstr1 information once only (perhaps we can use this in the directory scheme!!!!!!!! 
// [OBSOLETE - DX 20190717]         // and really only calculate once; but not that expensive, may be more expensive to store --memory!)
// [OBSOLETE - DX 20190717]         vector<double> D1,F1;
// [OBSOLETE - DX 20190717]         cellDiagonal(xstr1_tmp,D1,F1,1);
// [OBSOLETE - DX 20190717]         xstr1_tmp.lattice=GetClat(xstr1_tmp.a,xstr1_tmp.b,xstr1_tmp.c,xstr1_tmp.alpha,xstr1_tmp.beta,xstr1_tmp.gamma);
// [OBSOLETE - DX 20190717]         for(uint iat=0; iat<xstr1_tmp.atoms.size(); iat++){
// [OBSOLETE - DX 20190717]           xstr1_tmp.atoms[iat].cpos=F2C(xstr1_tmp.lattice,xstr1_tmp.atoms[iat].fpos);
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717]         vector<double> all_nn1 = computeNearestNeighbors(xstr1_tmp);
// [OBSOLETE - DX 20190717]         // DX NEW
// [OBSOLETE - DX 20190717]         //for(uint n=0; n<num_proc; n++){
// [OBSOLETE - DX 20190717]         //  vector<xstructure> vprotos_tmp;
// [OBSOLETE - DX 20190717]         //  vvprotos.push_back(vprotos_tmp);
// [OBSOLETE - DX 20190717]         //  xstr1_for_thread.push_back(xstr1_tmp);
// [OBSOLETE - DX 20190717]         //  possible_minMis.push_back(1.0);
// [OBSOLETE - DX 20190717]         //}
// [OBSOLETE - DX 20190717]       
// [OBSOLETE - DX 20190717]         vector<xmatrix<double> > lattices;
// [OBSOLETE - DX 20190717]         vector<xmatrix<double> > clattices;
// [OBSOLETE - DX 20190717]         vector<vector<uint> > ij_index;
// [OBSOLETE - DX 20190717]         vector<double> latt_devs;
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         //vector<std::thread> threads0;
// [OBSOLETE - DX 20190717]         vector<xvector<double> > translation_vectors;
// [OBSOLETE - DX 20190717]         quadrupletSearch(q1,xstr_LFA_supercell,xstr2,translation_vectors,ij_index); //DX 20190701 - added xstr2
// [OBSOLETE - DX 20190717]         //cerr << "FINDING TRANSLATION VECTORS: " << endl;
// [OBSOLETE - DX 20190717]         double abs_det_q1=abs(det(q1));
// [OBSOLETE - DX 20190717]         xvector<double> abc_angles_q1=Getabc_angles(q1,DEGREES);
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         buildSimilarLattices(translation_vectors, q1, abs_det_q1, abs_det_q1, abc_angles_q1, lattices, clattices, latt_devs, optimize_match);
// [OBSOLETE - DX 20190717]         if(LDEBUG) {cerr << "pflow::threadGeneration: Number of lattices to compare: " << lattices.size() << endl;}
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         if(lattices.size()>0){
// [OBSOLETE - DX 20190717]           for(uint n=0; n<num_proc; n++){
// [OBSOLETE - DX 20190717]             vector<xstructure> vprotos_tmp;
// [OBSOLETE - DX 20190717]             vvprotos.push_back(vprotos_tmp);
// [OBSOLETE - DX 20190717]             xstr1_for_thread.push_back(xstr1_tmp);
// [OBSOLETE - DX 20190717]             possible_minMis.push_back(1.0);
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717] #ifdef AFLOW_COMPARE_MULTITHREADS_ENABLE
// [OBSOLETE - DX 20190717]           vector<std::thread> threads1;
// [OBSOLETE - DX 20190717]           vector<vector<xmatrix<double> > > lattices_split;
// [OBSOLETE - DX 20190717]           vector<vector<xmatrix<double> > > clattices_split;
// [OBSOLETE - DX 20190717]           vector<vector<double> > latt_devs_split;
// [OBSOLETE - DX 20190717]           uint num_per_thread = lattices.size()/num_proc;
// [OBSOLETE - DX 20190717]           uint residual = lattices.size()%num_proc;
// [OBSOLETE - DX 20190717]           bool accounted_for_residual=false;
// [OBSOLETE - DX 20190717]           if(residual!=0){
// [OBSOLETE - DX 20190717]             num_per_thread+=1;
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]           uint count = 0;
// [OBSOLETE - DX 20190717]           uint thread_count = 0;
// [OBSOLETE - DX 20190717]           vector<xmatrix<double> > latt_tmp, clatt_tmp;
// [OBSOLETE - DX 20190717]           vector<double> tmp_dev;
// [OBSOLETE - DX 20190717]           for(uint l=0; l<lattices.size(); l++){
// [OBSOLETE - DX 20190717]             latt_tmp.push_back(lattices[l]); clatt_tmp.push_back(clattices[l]); tmp_dev.push_back(latt_devs[l]);
// [OBSOLETE - DX 20190717]             count+=1;
// [OBSOLETE - DX 20190717]             if(count == num_per_thread && thread_count<num_proc-1){
// [OBSOLETE - DX 20190717]               thread_count+=1;
// [OBSOLETE - DX 20190717]               lattices_split.push_back(latt_tmp);
// [OBSOLETE - DX 20190717]               clattices_split.push_back(clatt_tmp);
// [OBSOLETE - DX 20190717]               latt_devs_split.push_back(tmp_dev);
// [OBSOLETE - DX 20190717]               latt_tmp.clear(); clatt_tmp.clear(); tmp_dev.clear();
// [OBSOLETE - DX 20190717]               count = 0;
// [OBSOLETE - DX 20190717]             }
// [OBSOLETE - DX 20190717]             else if(thread_count==num_proc-1 && l==lattices.size()-1){
// [OBSOLETE - DX 20190717]               thread_count+=1;
// [OBSOLETE - DX 20190717]               lattices_split.push_back(latt_tmp);
// [OBSOLETE - DX 20190717]               clattices_split.push_back(clatt_tmp);
// [OBSOLETE - DX 20190717]               latt_devs_split.push_back(tmp_dev);
// [OBSOLETE - DX 20190717]               latt_tmp.clear(); clatt_tmp.clear(); tmp_dev.clear();
// [OBSOLETE - DX 20190717]               count = 0;
// [OBSOLETE - DX 20190717]             }
// [OBSOLETE - DX 20190717]             if(!accounted_for_residual && residual!=0 && thread_count==residual){
// [OBSOLETE - DX 20190717]               accounted_for_residual=true;
// [OBSOLETE - DX 20190717]               num_per_thread=num_per_thread-1;
// [OBSOLETE - DX 20190717]             }
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]           uint recovered=0;
// [OBSOLETE - DX 20190717]           //Need the following safety in case the number of threads is greater than the number of lattices to test
// [OBSOLETE - DX 20190717]           uint num_of_threads=0;
// [OBSOLETE - DX 20190717]           if(lattices_split.size()>=num_proc){
// [OBSOLETE - DX 20190717]             num_of_threads=num_proc;
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]           else if(lattices_split.size()<num_proc){
// [OBSOLETE - DX 20190717]             num_of_threads=lattices_split.size();
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]           for(uint n=0; n<num_of_threads; n++){
// [OBSOLETE - DX 20190717]             for(uint h=0;h<lattices_split[n].size();h++){
// [OBSOLETE - DX 20190717]               recovered+=1;
// [OBSOLETE - DX 20190717]               //cerr << "recovered: " << recovered << " - " << lattices_split[n][h] << endl;
// [OBSOLETE - DX 20190717]             }
// [OBSOLETE - DX 20190717]           } 
// [OBSOLETE - DX 20190717]           if(recovered != lattices.size()){
// [OBSOLETE - DX 20190717]             cerr << "The splitting of jobs failed...not all were accounted for: " << recovered << " != " << lattices.size() << endl;
// [OBSOLETE - DX 20190717]             exit(1);
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]           if(LDEBUG) {cerr << "pflow::threadGeneration: Performing structure search on " << lattices.size() << " lattices ..." << endl;}
// [OBSOLETE - DX 20190717]           /*
// [OBSOLETE - DX 20190717]              for(uint n=0; n<num_of_threads; n++){
// [OBSOLETE - DX 20190717]              threads1.push_back(std::thread(structureSearch,lfa,all_nn1,xstr,
// [OBSOLETE - DX 20190717]              std::ref(vvprotos[n]),std::ref(xstr1_for_thread[n]),xstr2,type_match,std::ref(possible_minMis[n]),
// [OBSOLETE - DX 20190717]              std::ref(lattices_split[n]),std::ref(clattices_split[n]),std::ref(latt_devs_split[n]),
// [OBSOLETE - DX 20190717]              optimize_match));
// [OBSOLETE - DX 20190717]              }         
// [OBSOLETE - DX 20190717]              for(uint t=0;t<threads1.size();t++){
// [OBSOLETE - DX 20190717]              threads1[t].join();
// [OBSOLETE - DX 20190717]              }*/
// [OBSOLETE - DX 20190717] #else
// [OBSOLETE - DX 20190717]           uint n=0;
// [OBSOLETE - DX 20190717]           structureSearch(lfa,all_nn1,xstr,vvprotos[n],xstr1_for_thread[n],xstr2,type_match,possible_minMis[n],
// [OBSOLETE - DX 20190717]               lattices,clattices,latt_devs,optimize_match);
// [OBSOLETE - DX 20190717] #endif
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717]         //cerr << "========== possible_minMis.size(): " << possible_minMis.size() << endl;
// [OBSOLETE - DX 20190717]         for(uint p=0;p<possible_minMis.size();p++){
// [OBSOLETE - DX 20190717]           if(p==0 && y==0 && x==0){ // DX 2/8/17 - need to add x==0 ortherwise matches can be overwritten
// [OBSOLETE - DX 20190717]             minMis=possible_minMis[p];
// [OBSOLETE - DX 20190717]             xstr1=xstr1_for_thread[p];
// [OBSOLETE - DX 20190717]             vprotos=vvprotos[p];
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]           else {
// [OBSOLETE - DX 20190717]             if(possible_minMis[p]<=minMis){
// [OBSOLETE - DX 20190717]               minMis=possible_minMis[p];
// [OBSOLETE - DX 20190717]               xstr1=xstr1_for_thread[p];
// [OBSOLETE - DX 20190717]               vprotos=vvprotos[p];
// [OBSOLETE - DX 20190717]             }
// [OBSOLETE - DX 20190717]           }
// [OBSOLETE - DX 20190717]           //cerr << "minMis: " << minMis << endl;
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717]         //break if(minMis<=0.1) break;
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]         // //DX 20190226 - fast return, no need to check other LFAs if a match is found - START
// [OBSOLETE - DX 20190717]         if(minMis<0.1 && !optimize_match){
// [OBSOLETE - DX 20190717]           return;
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717]         // //DX 20190226 - fast return, no need to check other LFAs if a match is found - END
// [OBSOLETE - DX 20190717]         // DX 20190318 - START
// [OBSOLETE - DX 20190717]         if(test_one_lfa_only){
// [OBSOLETE - DX 20190717]           return;
// [OBSOLETE - DX 20190717]         }
// [OBSOLETE - DX 20190717]         // DX 20190318 - START
// [OBSOLETE - DX 20190717]       } 
// [OBSOLETE - DX 20190717]     } 
// [OBSOLETE - DX 20190717]   }  
// [OBSOLETE - DX 20190717] }

// ***************************************************************************
// structureSearch
// ***************************************************************************
namespace compare{
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
      bool optimize_match){ 

    bool LDEBUG=(false || XHOST.DEBUG);
    bool VERBOSE=false;

    double mis=1;  
    double mag_dis=0.0; double mag_fail=0.0;
    xstructure proto;
    //xstructure xstr2_tmp = xstr2;
   

    vector<string> species_str1=sortSpeciesByFrequency(xstr1);
    deque<_atom> xstr1_atoms;
    for(uint i=0;i<species_str1.size();i++){
      for(uint j=0;j<xstr1.atoms.size();j++){
        if(species_str1[i]==xstr1.atoms[j].name){
          xstr1_atoms.push_back(xstr1.atoms[j]);
        }
      }
    }

    for(uint p=start_index;p<end_index;p++){ //DX 20191107 switching end index convention <= vs <
      if(LDEBUG){
        cerr << "compare::structureSearch: Trying lattice " << p << endl;
        cerr << "lattice=" << lattices[p] << endl;
      }

      // ---------------------------------------------------------------------------
      // lattice rotation (beta) 
      xmatrix<double> lattice_rotation = aurostd::inverse(lattices[p])*clattices[p];
      if(LDEBUG){
        cerr << "rotation between lattice " << lattices[p] << endl;
        cerr << "and clattice: " << clattices[p] << endl;
        cerr << "is " << lattice_rotation << endl;
        cerr << "validate (lattice * rotation =? clattice) : " << lattices[p]*lattice_rotation << endl;
      }

      // ---------------------------------------------------------------------------
      // make smaller lattice the new lattice in the supercell structure 
      // note: lattices[p] are oriented wrt to supercell (it has to be), otherwise could break periodicity
      proto=xstr_supercell; //DX 20190530 - added "_supercell"; more descriptive
      proto.lattice=lattices[p];

      // ---------------------------------------------------------------------------
      // C2F - (i.e., will provide fractional coordinates wrt to new lattice) 
      // AND remove all atoms outside unit cell based on fractional coordinates
      // speed increase: ensure this is in cell before computing F2C 
      // (don't calculate unnecessary matrix-vector multiplication)
      // Note: C2F (done later) changes lattice to one that is aligned with Cartesian directions (a along +X, etc.) 
      //       this is like rotating the global coordinates, therefore, fpos does not change
      deque<_atom> new_basis_2;
      for(uint iat=0;iat<proto.atoms.size();iat++){
	      proto.atoms[iat].fpos=C2F(proto.lattice,proto.atoms[iat].cpos);
        if(atomInCell(proto.atoms[iat])){
          new_basis_2.push_back(proto.atoms[iat]);
        }
      }
      
      xstructure proto_new;
      proto_new.title=proto.title;
      proto_new.lattice=clattices[p];
      
      // DX NEW - START =======================
      xmatrix<double> f2c = trasp(proto_new.lattice); //DX 20190717
      xmatrix<double> c2f = aurostd::inverse(trasp(proto_new.lattice)); //DX 20190717
      //DX 20190717 [OBSOLETE] xmatrix<double> f2c = trasp(proto.lattice);
      //DX 20190717 [OBSOLETE] xmatrix<double> c2f = aurostd::inverse(trasp(proto.lattice));
      bool skew = false;
      double tol=0.01;
      
      deque<_atom> new_basis;
      for(uint j=0;j<new_basis_2.size();j++){
        bool duplicate_lattice_point=false;
        for(uint a=0; a<new_basis.size(); a++){
          xvector<double> tmp = BringInCell(new_basis_2[j].fpos,1e-10);
          if(SYM::MapAtom(new_basis[a].fpos,tmp,proto_new.lattice,f2c,skew,tol)){
            duplicate_lattice_point=true;
            break;
          }
        }
        if(duplicate_lattice_point==false){
          new_basis_2[j].fpos = BringInCell(new_basis_2[j].fpos,1e-10);
          new_basis_2[j].cpos = f2c*new_basis_2[j].fpos;
          new_basis.push_back(new_basis_2[j]);
        }
      }
      std::stable_sort(new_basis.begin(),new_basis.end(),sortAtomsNames); //DX 20190709 - need to sort now
      proto_new.atoms = new_basis;
      proto_new.BringInCell(1e-10); 
      proto_new.FixLattices();
      proto_new.SpeciesPutAlphabetic();
      deque<int> sizes = SYM::arrange_atoms(new_basis);
      proto_new = pflow::SetNumEachType(proto_new, sizes);
      proto_new.species = proto.species; //DX 20190718
      proto = proto_new;
      
      if(sameSpecies(proto,xstr1,false)){
        vector<string> species_str2=sortSpeciesByFrequency(proto);
        vector<double> all_nn_proto;
        bool all_nn_calculated = false;
        for(uint iat=0; iat<proto.atoms.size();iat++){
          if(proto.atoms[iat].name==lfa){
            proto.ShifOriginToAtom(iat);
            proto.BringInCell(1e-10);
            if(VERBOSE){
              cerr << "compare::structureSearch: orig structure " << xstr1 << endl;
              cerr << "compare::structureSearch: structure " << proto << endl;
            }
            deque<_atom> proto_atoms;
            for(uint i=0;i<species_str2.size();i++){
              for(uint j=0;j<proto.atoms.size();j++){
                if(species_str2[i]==proto.atoms[j].name){
                  proto_atoms.push_back(proto.atoms[j]);
                }
              }
            }
            vector<uint> im1, im2;
            vector<double> min_dists;
            if(findMatch(xstr1_atoms,proto_atoms,proto.lattice,im1,im2,min_dists,type_match)){;
              if(VERBOSE){
                for(uint m=0;m<im1.size();m++){
                  cerr << "compare::structureSearch: " << im1[m] << " == " << im2[m] << " : dist=" << min_dists[m] << endl;
                }
              }
              double cd, f;
              // Only calculate the NN for the proto if we found suitable matches.  
              // Only calculate once, nothing changes between shifts to origin (affine)
              if(!all_nn_calculated){
                all_nn_proto = computeNearestNeighbors(proto);
                if(VERBOSE){
                  cerr << "compare::structureSearch: Nearest neighbors:" << endl;
                  for(uint a=0;a<all_nn_proto.size();a++){
                    cerr << "compare::structureSearch: Nearest neighbor distance from " << a << " atom: " << all_nn_proto[a] << endl;
                  }
                }
                all_nn_calculated = true;
              }
              coordinateDeviation(xstr1,proto,all_nn1,all_nn_proto,im1,im2,min_dists,cd,f);
              if((xstr1.atoms[0].spin_is_given && proto.atoms[0].spin_is_given) || 
                 (xstr1.atoms[0].noncoll_spin_is_given && proto.atoms[0].noncoll_spin_is_given)){ 
                magneticDeviation(xstr1,proto,im1,im2,mag_dis,mag_fail);
                mis=computeMagneticMisfit(latt_devs[p],cd,f,mag_dis,mag_fail);
                if(LDEBUG){
                  cerr << "with spin: mis,latt_dev,cd,f,mag_dis: " << mis << ", " <<latt_devs[p] << ", " << cd << ", " << f << ", " << mag_dis << ", " << mag_fail <<  endl;
                  double tmp_mis=computeMisfit(latt_devs[p],cd,f);
                  cerr << "without spin: mis,latt_dev,cd,f: " << tmp_mis << ", " <<latt_devs[p] << ", " << cd << ", " << f <<  endl;
                }
              }
              else{
                mis=computeMisfit(latt_devs[p],cd,f);
                if(LDEBUG){
                  cerr << "mis,latt_dev,cd,f: " << mis << ", " <<latt_devs[p] << ", " << cd << ", " << f <<  endl;
                }
              }
              if(mis<min_misfit){
                //cerr << "storing: " << proto << endl;
                vprotos.clear();
                vprotos.push_back(proto);
                min_misfit=mis;
                min_latt_dev=latt_devs[p];
                min_coordinate_dis=cd;
                min_failure=f;
                min_magnetic_dis=mag_dis;
                min_magnetic_failure=mag_fail;
                index_match_1 = im1;
                index_match_2 = im2;
                min_distances = min_dists;
              }
              // If we want to simply find a match and not find the best match, we can exit early
	            if(mis<0.1 && !optimize_match) {
                return true;
              }
            }
          }
        }
      }// end of if protos.size()...
      else{
        if(LDEBUG){
          cerr << "compare::structureSearch: Atom counts do not match: orig=" << aurostd::joinWDelimiter(xstr1.num_each_type,",") << " vs test=" << aurostd::joinWDelimiter(proto.num_each_type,",") << endl;
        }
      }
    }
    return true;
  }  
}

// ***************************************************************************
// Quadruplet Search
// ***************************************************************************
namespace compare{
  void quadrupletSearch(const xmatrix<double>& q1, const xstructure& xstr_LFA_supercell, 
      const xstructure& xstr,
      vector<xvector<double> >& lattice_vecs, vector<vector<uint> >& ij_index){

    // This function scans through the possible quadruplets (sets of 4 LFA atoms) i
    // to find a lattice which is commensurate with the reference structure (xstr1). 
    // This function is parallelized since it is the time-limiting function.

    bool LDEBUG=(false || XHOST.DEBUG);
    string function_name = "compare::quadrupletSearch()";
    bool relative_tolerance=true;

    double min_q1_a = 0.0; double max_q1_a = 0.0;
    double min_q1_b = 0.0; double max_q1_b = 0.0;
    double min_q1_c = 0.0; double max_q1_c = 0.0;

    if(relative_tolerance){
      min_q1_a = aurostd::modulus(q1(1))-aurostd::modulus(q1(1))*0.1;
      max_q1_a = aurostd::modulus(q1(1))+aurostd::modulus(q1(1))*0.1;
      min_q1_b = aurostd::modulus(q1(2))-aurostd::modulus(q1(2))*0.1;
      max_q1_b = aurostd::modulus(q1(2))+aurostd::modulus(q1(2))*0.1;
      min_q1_c = aurostd::modulus(q1(3))-aurostd::modulus(q1(3))*0.1;
      max_q1_c = aurostd::modulus(q1(3))+aurostd::modulus(q1(3))*0.1;
    }
    else{
      min_q1_a = aurostd::modulus(q1(1))-1.0;
      max_q1_a = aurostd::modulus(q1(1))+1.0;
      min_q1_b = aurostd::modulus(q1(2))-1.0;
      max_q1_b = aurostd::modulus(q1(2))+1.0;
      min_q1_c = aurostd::modulus(q1(3))-1.0;
      max_q1_c = aurostd::modulus(q1(3))+1.0;
    }

    if(LDEBUG) {
      cerr << function_name << ": Lattice parameters: " << aurostd::modulus(q1(1)) << ", " << aurostd::modulus(q1(2)) << ", " << aurostd::modulus(q1(3)) << endl;
      cerr << function_name << ": Modulus search range for lattice vector a: " << min_q1_a << " - " << max_q1_a << endl;
      cerr << function_name << ": Modulus search range for lattice vector b: " << min_q1_b << " - " << max_q1_b << endl;
      cerr << function_name << ": Modulus search range for lattice vector c: " << min_q1_c << " - " << max_q1_c << endl;
    }

    xvector<double> tmp_vec;
    double tmp_mod = 0.0;
    
    // Search all possible vectors with modulus comparable to one of the lattice vectors 
    for(uint i=0; i<xstr_LFA_supercell.atoms.size(); i++){
      //DX TEST for(uint j=i+1; j<xstr_LFA_supercell.atoms.size(); j++){ //upper triangular
      //DX TEST   tmp_vec = xstr_LFA_supercell.atoms[j].cpos-xstr_LFA_supercell.atoms[i].cpos;
      uint j=0;
      tmp_vec = xstr_LFA_supercell.atoms[i].cpos;
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
      //DX TEST }
    }

    if(LDEBUG) {
      cerr << function_name << ": Number of potential lattice vectors: " << lattice_vecs.size() << endl;
    } 

    // ---------------------------------------------------------------------------
    // check if vectors preserve crystal periodicity
    // if only one LFA atom in unit cell -> lattice point -> vectors=lattice vectors (by definition)
    if(aurostd::min(xstr.num_each_type)!=1){
      // Removing non-periodic lattice vectors
      vector<xvector<double> > lattice_vecs_periodic;
      for(uint i=0;i<lattice_vecs.size();i++){
        if(vectorPeriodic(lattice_vecs[i],xstr,ij_index[i][0],ij_index[i][1])){ //DX 20190701 - xstr_LFA_supercell to xstr
          lattice_vecs_periodic.push_back(lattice_vecs[i]);
          //DX 20190318 [OBSOLETE] lattice_vecs_periodic.push_back(-lattice_vecs[i]);
        }
      }
      //vector<xvector<double> > final_lattice_vecs = lattice_vecs_periodic; //DX 20190320
      //DX NOT NEEDED ANYMORE, ALREADY ACCOUNTED FOR - START
      //DX 20190318 - only store negative if not a duplicate - START
      /*
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
      */
      //DX NOT NEEDED ANYMORE, ALREADY ACCOUNTED FOR - END
      lattice_vecs = lattice_vecs_periodic; //DX 20190320
    }
    if(LDEBUG) {
      cerr << function_name << ": Number of lattice vectors (preserves periodicity): " << lattice_vecs.size() << endl;
      for(uint i=0;i<lattice_vecs.size();i++){
        cerr << function_name << ": lattice vector " << i << ": " << lattice_vecs[i] << " (" << aurostd::modulus(lattice_vecs[i]) << ")" << endl; 
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
    bool VERBOSE=false;
    string function_name = "compare::buildSimilarLattices():";

    // ---------------------------------------------------------------------------
    // sort via smallest misfit for speed up
    bool sort_via_lattice_deviation = true; //DX 20190626 - speed increase
    bool relative_tolerance = false; //DX 20190703

    vector<double> D1,F1;
    cellDiagonal(q1,D1,F1,1);

    double tol_vol=0.1;
    double det_tol=tol_vol*abs_det_q1;

    // ---------------------------------------------------------------------------
    // tolerance for lattice vectors: relative or absolute 
    double tol_a=1e9; double tol_b=1e9; double tol_c=1e9;
    if(relative_tolerance){
      tol_a=abc_angles_q1[1]*0.3;
      tol_b=abc_angles_q1[2]*0.3;
      tol_c=abc_angles_q1[3]*0.3;
    }
    else{
      tol_a=1.0; tol_b=1.0; tol_c=1.0; // 1 Angstrom
    }

    // ---------------------------------------------------------------------------
    // LDEBUG: print lattice tolerances 
    if(LDEBUG) {
      if(relative_tolerance){
        cerr << function_name << " Tolerance for a (Angstroms): " << tol_a << endl;
        cerr << function_name << " Tolerance for b (Angstroms): " << tol_b << endl;
        cerr << function_name << " Tolerance for c (Angstroms): " << tol_c << endl;
        cerr << function_name << " Tolerance for alpha (degrees): " << abc_angles_q1[4]*0.3 << endl;
        cerr << function_name << " Tolerance for beta (degrees): " << abc_angles_q1[5]*0.3 << endl;
        cerr << function_name << " Tolerance for gamma (degrees): " << abc_angles_q1[6]*0.3 << endl;
      }
      else{
        cerr << function_name << " Tolerance for a (Angstroms): " << tol_a << endl;
        cerr << function_name << " Tolerance for b (Angstroms): " << tol_b << endl;
        cerr << function_name << " Tolerance for c (Angstroms): " << tol_c << endl;
        cerr << function_name << " Tolerance for alpha (degrees): " << 5 << endl;
        cerr << function_name << " Tolerance for beta (degrees): " << 5 << endl;
        cerr << function_name << " Tolerance for gamma (degrees): " << 5 << endl;
      }
    }

    xmatrix<double> tmp_lattice(3,3);
    xmatrix<double> tmp_clatt(3,3);

    // ---------------------------------------------------------------------------
    // compute lenths of possible lattices before-hand (faster than on-the-fly)
    int store=0;
    vector<double> translations_mod;
    for(uint i=0;i<translation_vectors.size();i++){
      translations_mod.push_back(aurostd::modulus(translation_vectors[i]));
    }
    if(LDEBUG) { cerr << function_name << " Number of lattice vectors: " << translation_vectors.size() << endl; }

    // ---------------------------------------------------------------------------
    // build all possible unit cells with combinations of lattice vectors 
    // (order matters, hence not upper triangular)
    for(uint i=0;i<translation_vectors.size();i++){
      // ---------------------------------------------------------------------------
      // check lattice vector length: a
      if(abs(translations_mod[i]-abc_angles_q1[1])<tol_a){ //check a
        for(uint j=0;j<translation_vectors.size();j++){
          if(j!=i){
            // ---------------------------------------------------------------------------
            // check lattice vector length: b
            if(abs(translations_mod[j]-abc_angles_q1[2])<tol_b){ // check b
              for(uint k=0;k<translation_vectors.size();k++){
                if(k!=i && k!=j){
                  // ---------------------------------------------------------------------------
                  // check lattice vector length: c
                  if(abs(translations_mod[k]-abc_angles_q1[3])<tol_c){ //check c
                    tmp_lattice = SYM::xvec2xmat(translation_vectors[i],translation_vectors[j],translation_vectors[k]);         
                    // ---------------------------------------------------------------------------
                    // check determinant 
                    if(abs(abs_det_q1-abs(det(tmp_lattice))) < det_tol){ //check determinant/volume
                      xvector<double> abc_angles_q2=Getabc_angles(tmp_lattice,DEGREES);
                      // ---------------------------------------------------------------------------
                      // check angles
                      if(checkTolerance(abc_angles_q1,abc_angles_q2)==false){
                        double tmp_latt_dev = checkLatticeDeviation(xstr1_vol,tmp_lattice,D1,F1);
                        // ---------------------------------------------------------------------------
                        // check lattice deviation (time-saver) 
                        // case 1: NOT optimize match: keep lattices with deviation smaller than Burzlaff's matching requirement)
                        //         otherwise, there is no possible way that it could match with anything 
                        // case 2: optimize match: keep lattices with deviation smaller than Burzlaff's same-family requirement)
                        // otherwise, there is no possible way that it could match with anything or be in the same-family
                        if((!optimize_match && tmp_latt_dev <= 0.1) || (optimize_match && tmp_latt_dev <= 0.2)) { //fast match doesn't care about finding same family information //DX 20190318 - removed unique since it doesn't exist yet
                          // ---------------------------------------------------------------------------
                          // now check uniqueness (this is more expensive than checking lattice deviation, hence why it is further in nesting) 
                          bool unique = true;
                          uint placement_index = lattices.size(); //DX 20190626 //default to the end
                          for(uint t=0;t<lattices.size();t++){
                            if(identical(tmp_lattice,lattices[t],1e-10)){
                              unique=false;
                              break;
                            }
                            // store placement/index based on lattice deviation 
                            if(sort_via_lattice_deviation && tmp_latt_dev < latt_devs[t] && placement_index == lattices.size()){ // third condition necessary; otherwise it could move down in order
                              placement_index = t;
                            }
                          }
                          if(unique){
                            // ---------------------------------------------------------------------------
                            // order/re-order by minimum lattice deviation (speed increase: more likely to find matches faster) 
                            // append to the end 
                            if(placement_index==lattices.size()){ // push_back: either it is the first or it goes to the end
                              lattices.push_back(tmp_lattice); // stores original original orientation
                              tmp_clatt=GetClat(abc_angles_q2[1],abc_angles_q2[2],abc_angles_q2[3],abc_angles_q2[4],abc_angles_q2[5],abc_angles_q2[6]);
                              clattices.push_back(tmp_clatt); // store Cartesian lattice, alignes with XYZ coordinates
                              latt_devs.push_back(tmp_latt_dev);
                            }
                            // insert to certain location via index 
                            else{
                              lattices.insert(lattices.begin()+placement_index, tmp_lattice); // stores original original orientation
                              tmp_clatt=GetClat(abc_angles_q2[1],abc_angles_q2[2],abc_angles_q2[3],abc_angles_q2[4],abc_angles_q2[5],abc_angles_q2[6]);
                              clattices.insert(clattices.begin()+placement_index, tmp_clatt); // store Cartesian lattice, alignes with XYZ coordinates
                              latt_devs.insert(latt_devs.begin()+placement_index, tmp_latt_dev);
                            }
                            store++;
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
    }

    // ---------------------------------------------------------------------------
    // print lattices and corresponding lattice deviation 
    if(VERBOSE){
      cerr << "q1: " << q1 << endl;
      cerr << "det(q1): " << det(q1) << endl;
      cerr << "abc angles q1: " << abc_angles_q1 << endl;
      for(uint i=0;i<lattices.size();i++){
        cerr << function_name << endl << " lattice: " << endl << lattices[i] << endl << " clattice: " << clattices[i] << endl << " volume: " << det(lattices[i]) << endl << " lattice deviation: " << endl << latt_devs[i] << endl;
        xvector<double> abc_angles_q2=Getabc_angles(lattices[i],DEGREES);
        cerr << "abc angles: " << abc_angles_q2 << endl;
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

// [OBSOLETE - DX 20190717] // ***************************************************************************
// [OBSOLETE - DX 20190717] // Internal structure
// [OBSOLETE - DX 20190717] // ***************************************************************************
// [OBSOLETE - DX 20190717] namespace compare{
// [OBSOLETE - DX 20190717]   bool structureSearch(const string& lfa, 
// [OBSOLETE - DX 20190717]                         const vector<double>& all_nn1, 
// [OBSOLETE - DX 20190717] 			const xstructure& xstr_supercell, //DX 20190530 - added "_supercell"; more descriptive 
// [OBSOLETE - DX 20190717] 			vector<xstructure>& vprotos, xstructure& xstr1, const xstructure& xstr2, 
// [OBSOLETE - DX 20190717] 			const int& type_match, double& possible_minMis,
// [OBSOLETE - DX 20190717] 			vector<xmatrix<double> >& lattices,
// [OBSOLETE - DX 20190717] 			vector<xmatrix<double> >& clattices, 
// [OBSOLETE - DX 20190717] 			vector<double>& latt_devs, 
// [OBSOLETE - DX 20190717] 			const bool& optimize_match){ 
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]     bool LDEBUG=(false || XHOST.DEBUG);
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]     double mis=1;  
// [OBSOLETE - DX 20190717]     xstructure proto;
// [OBSOLETE - DX 20190717]     int flag=0;
// [OBSOLETE - DX 20190717]     xstructure xstr2_tmp = xstr2;
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]     for(uint p=0;p<lattices.size();p++){
// [OBSOLETE - DX 20190717]       if(LDEBUG) {
// [OBSOLETE - DX 20190717]         cerr << "compare::structureSearch: Trying lattice " << p << endl;
// [OBSOLETE - DX 20190717]       }
// [OBSOLETE - DX 20190717]       proto=xstr_supercell; //DX 20190530 - added "_supercell"; more descriptive
// [OBSOLETE - DX 20190717]       proto.lattice=lattices[p];
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]       // Transform
// [OBSOLETE - DX 20190717]       for(uint iat=0;iat<proto.atoms.size();iat++){
// [OBSOLETE - DX 20190717] 	proto.atoms[iat].fpos=C2F(proto.lattice,proto.atoms[iat].cpos);
// [OBSOLETE - DX 20190717]       }
// [OBSOLETE - DX 20190717]       proto.lattice=clattices[p];
// [OBSOLETE - DX 20190717] 	for(uint iat=0;iat<proto.atoms.size();iat++){
// [OBSOLETE - DX 20190717] 	proto.atoms[iat].cpos=F2C(proto.lattice,proto.atoms[iat].fpos);
// [OBSOLETE - DX 20190717]       }
// [OBSOLETE - DX 20190717]       xstructure proto_new;
// [OBSOLETE - DX 20190717]       proto_new.title=proto.title;
// [OBSOLETE - DX 20190717]       proto_new.lattice=clattices[p];
// [OBSOLETE - DX 20190717] 
// [OBSOLETE - DX 20190717]       // DX NEW - START =======================
// [OBSOLETE - DX 20190717]       xmatrix<double> f2c = trasp(proto.lattice);
// [OBSOLETE - DX 20190717]       xmatrix<double> c2f = aurostd::inverse(trasp(proto.lattice));
// [OBSOLETE - DX 20190717]       bool skew = false;
// [OBSOLETE - DX 20190717]       double tol=0.01;
// [OBSOLETE - DX 20190717]       deque<_atom> new_basis;
// [OBSOLETE - DX 20190717]       for(uint j=0;j<proto.atoms.size();j++){
// [OBSOLETE - DX 20190717] 	if(new_basis.size()==0){
// [OBSOLETE - DX 20190717] 	  proto.atoms[j].fpos = BringInCell(proto.atoms[j].fpos,1e-10);
// [OBSOLETE - DX 20190717] 	  proto.atoms[j].cpos = f2c*proto.atoms[j].fpos;
// [OBSOLETE - DX 20190717] 	  new_basis.push_back(proto.atoms[j]);
// [OBSOLETE - DX 20190717] 	  //proto_new.AddAtom(proto.atoms[j]);
// [OBSOLETE - DX 20190717] 	}
// [OBSOLETE - DX 20190717] 	else {
// [OBSOLETE - DX 20190717] 	  bool duplicate_lattice_point=false;
// [OBSOLETE - DX 20190717] 	  for(uint a=0; a<new_basis.size(); a++){
// [OBSOLETE - DX 20190717] 	    xvector<double> tmp = BringInCell(proto.atoms[j].fpos,1e-10);
// [OBSOLETE - DX 20190717] 	    if(SYM::MapAtom(new_basis[a].fpos,tmp,c2f,f2c,skew,tol)){
// [OBSOLETE - DX 20190717] 	      duplicate_lattice_point=true;
// [OBSOLETE - DX 20190717] 	      break;
// [OBSOLETE - DX 20190717] 	    }
// [OBSOLETE - DX 20190717] 	  }
// [OBSOLETE - DX 20190717] 	  if(duplicate_lattice_point==false){
// [OBSOLETE - DX 20190717] 	    proto.atoms[j].fpos = BringInCell(proto.atoms[j].fpos,1e-10);
// [OBSOLETE - DX 20190717] 	    proto.atoms[j].cpos = f2c*proto.atoms[j].fpos;
// [OBSOLETE - DX 20190717] 	    new_basis.push_back(proto.atoms[j]);
// [OBSOLETE - DX 20190717] 	    //proto_new.AddAtom(proto.atoms[j]);
// [OBSOLETE - DX 20190717] 	  }
// [OBSOLETE - DX 20190717] 	}
// [OBSOLETE - DX 20190717]       }
// [OBSOLETE - DX 20190717]       proto_new.atoms = new_basis;
// [OBSOLETE - DX 20190717]       proto_new.BringInCell(1e-10); 
// [OBSOLETE - DX 20190717]       proto_new.FixLattices();
// [OBSOLETE - DX 20190717]       proto_new.SpeciesPutAlphabetic();
// [OBSOLETE - DX 20190717]       deque<int> sizes = SYM::arrange_atoms(new_basis);
// [OBSOLETE - DX 20190717]       proto_new = pflow::SetNumEachType(proto_new, sizes);
// [OBSOLETE - DX 20190717]       proto = proto_new;
// [OBSOLETE - DX 20190717]       if(sameSpecies(proto,xstr1,false)){
// [OBSOLETE - DX 20190717] 	vector<double> all_nn_proto;
// [OBSOLETE - DX 20190717] 	bool all_nn_calculated = false;
// [OBSOLETE - DX 20190717] 	for(uint iat=0; iat<proto.atoms.size();iat++){
// [OBSOLETE - DX 20190717] 	  if(proto.atoms[iat].name==lfa){
// [OBSOLETE - DX 20190717] 	    proto.ShifOriginToAtom(iat);
// [OBSOLETE - DX 20190717] 	    proto.BringInCell(1e-10);
// [OBSOLETE - DX 20190717]             if(LDEBUG){
// [OBSOLETE - DX 20190717]               cerr << "compare::structureSearch: orig structure " << xstr1 << endl;
// [OBSOLETE - DX 20190717]               cerr << "compare::structureSearch: structure " << proto << endl;
// [OBSOLETE - DX 20190717]             }
// [OBSOLETE - DX 20190717] 	    vector<uint> im1, im2;
// [OBSOLETE - DX 20190717] 	    vector<double> min_dists;
// [OBSOLETE - DX 20190717] 	    if(findMatch(xstr1,proto,im1,im2,min_dists,type_match)){;
// [OBSOLETE - DX 20190717]               for(uint m=0;m<im1.size();m++){
// [OBSOLETE - DX 20190717]                 cerr << im1[m] << " == " << im2[m] << " : dist=" << min_dists[m] << endl;
// [OBSOLETE - DX 20190717]               }
// [OBSOLETE - DX 20190717] 	      double cd, f;
// [OBSOLETE - DX 20190717] 	      // Only calculate the NN for the proto if we found suitable matches.  
// [OBSOLETE - DX 20190717] 	      // Only calculate once, nothing changes between shifts to origin (affine)
// [OBSOLETE - DX 20190717] 	      if(!all_nn_calculated){
// [OBSOLETE - DX 20190717]                 all_nn_proto = computeNearestNeighbors(proto);
// [OBSOLETE - DX 20190717] 		if(LDEBUG) {
// [OBSOLETE - DX 20190717] 		  cerr << "compare::structureSearch: Nearest neighbors:" << endl;
// [OBSOLETE - DX 20190717] 		  for(uint a=0;a<all_nn_proto.size();a++){
// [OBSOLETE - DX 20190717] 		    cerr << "compare::structureSearch: Nearest neighbor distance from " << a << " atom: " << all_nn_proto[a] << endl;
// [OBSOLETE - DX 20190717] 		  }
// [OBSOLETE - DX 20190717] 		}
// [OBSOLETE - DX 20190717] 		all_nn_calculated = true;
// [OBSOLETE - DX 20190717]               }
// [OBSOLETE - DX 20190717] 	      coordinateDeviation(xstr1,proto,all_nn1,all_nn_proto,im1,im2,min_dists,cd,f);
// [OBSOLETE - DX 20190717] 	      mis=computeMisfit(latt_devs[p],cd,f);
// [OBSOLETE - DX 20190717] 	      //if(LDEBUG) {
// [OBSOLETE - DX 20190717] 	      //  cerr << "mis,latt_dev,cd,f: " << mis << ", " <<latt_devs[p] << ", " << cd << ", " << f <<  endl;
// [OBSOLETE - DX 20190717] 	      //}
// [OBSOLETE - DX 20190717] 	      if(flag==0){
// [OBSOLETE - DX 20190717] 	        flag=1;
// [OBSOLETE - DX 20190717] 		//cerr << "storing: " << proto << endl;
// [OBSOLETE - DX 20190717] 		vprotos.push_back(proto);
// [OBSOLETE - DX 20190717] 		possible_minMis=mis;
// [OBSOLETE - DX 20190717] 	      }
// [OBSOLETE - DX 20190717] 	      else {
// [OBSOLETE - DX 20190717] 	        if(mis<possible_minMis){
// [OBSOLETE - DX 20190717] 		  vprotos.clear();
// [OBSOLETE - DX 20190717] 		  possible_minMis=mis;
// [OBSOLETE - DX 20190717] 		  //cerr << "storing: " << proto << endl;
// [OBSOLETE - DX 20190717] 		  vprotos.push_back(proto); //to here
// [OBSOLETE - DX 20190717] 		}
// [OBSOLETE - DX 20190717] 	      }
// [OBSOLETE - DX 20190717]               // If we want to simply find a match and not find the best match, we can exit early
// [OBSOLETE - DX 20190717] 	      if(mis<0.1 && !optimize_match) {
// [OBSOLETE - DX 20190717]                 return true;
// [OBSOLETE - DX 20190717] 	        //DEBUGGING
// [OBSOLETE - DX 20190717] 		//cerr <<"Winning combo: "<<i<<","<<j<<","<<k<<","<<w<<endl;
// [OBSOLETE - DX 20190717] 		//cerr << "proto.lattice: " << proto.lattice << endl;
// [OBSOLETE - DX 20190717] 		//cerr << "lattice(1): " << modulus(proto.lattice(1)) << endl;
// [OBSOLETE - DX 20190717] 		//cerr << "lattice(2): " << modulus(proto.lattice(2)) << endl;
// [OBSOLETE - DX 20190717] 		//cerr << "lattice(3): " << modulus(proto.lattice(3)) << endl;
// [OBSOLETE - DX 20190717] 	      }
// [OBSOLETE - DX 20190717] 	    }
// [OBSOLETE - DX 20190717] 	  }
// [OBSOLETE - DX 20190717] 	}
// [OBSOLETE - DX 20190717]       }// end of if protos.size()...
// [OBSOLETE - DX 20190717]     }
// [OBSOLETE - DX 20190717]     return true;
// [OBSOLETE - DX 20190717]   }
// [OBSOLETE - DX 20190717] }
//---------------------------------------------------------------

// ***************************************************************************
// get prototype designations 
// ***************************************************************************
namespace compare{
  void getPrototypeDesignations(vector<StructurePrototype>& prototypes){

    //DX: Perhaps create a version that doesn't recalulate symmetry 
    // i.e., take sym as input

    uint start_index=0; uint end_index=prototypes.size(); //DX 20191107 switching end index convention <= vs <
    getPrototypeDesignationsInRange(prototypes,start_index,end_index);
  }

  void getPrototypeDesignationsInRange(vector<StructurePrototype>& prototypes, uint start_index, uint end_index){
    for(uint i=start_index;i<end_index;i++){ //DX 20191107 switching end index convention <= vs <
      anrl::structure2anrl(prototypes[i].structure_representative,false); //DX 20190829 - false for do not recalulate symmetry, save time
    }
  }
}

// AFLOW-XTAL-MATCH (compare crystal structures) - Functions
// Written by David Hicks (david.hicks@duke.edu) 
// Contributors: Carlo De Santo
// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                                                                         *
// ***************************************************************************
