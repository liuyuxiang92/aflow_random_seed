// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *             Stefano Curtarolo - Duke University - 2003-2019             *
// *                 Kesong Yang - Duke University 2010-2011                 *
// *                  Corey Oses - Duke University 2013-2019                 *
// *                                                                         *
// ***************************************************************************
// aflow_pocc.h and aflow_pocc*.cpp*
//
// completely revised approach to KESONG YANG's original implementation
// no recursion needed --- too memory intensive
// issues with UFF (structure comparison), multiply occupied sites, and 
// vacancies are all addressed robustly here
// 2013-2019: corey.oses@duke.edu
// 2010-2011: kesong.yang@gmail.com (LEGACY)

#ifndef _AFLOW_POCC_H_
#define _AFLOW_POCC_H_

//precision defines tol and paddings
//const int _AFLOW_POCC_PRECISION_ = 8;    //moved to aflow.h
//tolerances
//const double _AFLOW_POCC_ZERO_TOL_ = pow(10,-_AFLOW_POCC_PRECISION_);  //moved to aflow.h

//uff param modes
const uint BOND_MODE    = 0;
const uint NONBOND_MODE = 1;

//files
const string AFLOW_POCC_TAG="[AFLOW_POCC]";

//elements
//to make POCC chemistry-independent, calculate UFF based on standard set of elements
//we should not need to have more than 30 species, but if you want to add to list
//I recommend adding IN ORDER so that elements become more different left to right (like periodic table)
const string STD_ELEMENTS_LIST="Sc Ti V Cr Mn Fe Co Ni Cu Zn "
                               "Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
                               "La Hf Ta W Re Os Ir Pt Au Hg ";

namespace pocc {
  bool structuresGenerated();
  xstructure extractPARTCAR(const string& AflowIn);
  
  void generateStructures(const string& AflowIn, ostream& oss=cout);
  void generateStructures(const string& AflowIn, ofstream& FileMESSAGE, ostream& oss=cout);
  void generateStructures(const _xvasp& xvasp,const string& AflowIn,const _aflags& aflags,const _kflags& kflags,const _vflags& vflags, ofstream& FileMESSAGE, ostream& oss=cout);

	void updateProgressBar(unsigned long long int current, unsigned long long int end, ostream& oss=cout);
  vector<string> getElementsList();

  class POccUnit;  //forward declaration
  class POccSiteConfiguration;    //forward declaration
  struct UFFParamAtom;  //forward declaration
  struct UFFParamBond;  //forward declaration
 
  //bool sortAtoms(const _atom& a1,const _atom& a2);
  //bool getNextSiteConfiguration(vector<vector<POccSiteConfiguration> >& vv_count_configs,vector<uint>& v_config_order,vector<int>& v_config_iterators,vector<vector<int> >& v_types_config);
  bool getNextSiteConfiguration(vector<vector<POccSiteConfiguration> >& vv_count_configs,vector<int>& v_config_iterators,vector<vector<int> >& v_types_config);
  bool getNextSiteConfiguration(vector<uint>& v_config_order,vector<vector<int> >& v_types_config);
  bool getNextSiteConfiguration(vector<int>& site_config);
  vector<uint> getConfigOrder(vector<vector<int> >& v_types_config);
  string getUFFParameterString(const string& element);
  vector<UFFParamAtom> getTypes2UFFParamsMap(const vector<string>& elements);
  xstructure createNonPOccStructure(const xstructure& xstr_pocc,const vector<POccUnit>& pocc_sites,bool use_automatic_volumes_in=true);             //convert pocc xstructure to non-pocc
  vector<POccUnit> getPOccSites(const xstructure& xstr_pocc,ostream& oss=cout);
  vector<POccUnit> getPOccSites(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss=cout);
  vector<POccUnit> getPOccSites(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss=cout);
  vector<POccUnit> getPOccSites(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
  vector<uint> getTypes2PCMap(const xstructure& xstr);

  bool sortPOccSites(const POccUnit& p1,const POccUnit& p2);
  bool sortPOccGroups(const POccUnit& p1,const POccUnit& p2);

} // namespace pocc

namespace pocc {
  //POccUnit could be
  //1) a bunch of atoms on the same site (site)
  //2) within a site, a bunch of atoms with the same partial_occupation_value (group)
  class POccUnit: public xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccUnit(ostream& oss=cout);
      POccUnit(ofstream& FileMESSAGE,ostream& oss=cout);
      POccUnit(const _aflags& aflags,ostream& oss=cout);
      POccUnit(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccUnit(const POccUnit& b);
      //constructors - STOP
      ~POccUnit();
      const POccUnit& operator=(const POccUnit& b);
      //bool operator<(const POccUnit& other) const;
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - END

      bool m_initialized;
      _aflags m_aflags;
      uint site;                                        //reflexive
      bool partial_occupation_flag;
      double partial_occupation_value;
      vector<uint> v_occupants;
      vector<uint> v_types;
      vector<POccUnit> m_pocc_groups;
      bool is_inequivalent;                             //from iatoms of non-pocc structure
      uint equivalent;                                  //from iatoms of non-pocc structure

      //general setters
      void create(ostream& oss=cout);
      void create(ofstream& FileMESSAGE,ostream& oss=cout);
      void create(const _aflags& aflags,ostream& oss=cout);
      void create(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      void setAFlags(const _aflags& aflags);
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccUnit& b);
      //NECESSARY END CLASS METHODS - END
  };
}  // namespace pocc

namespace pocc {  
  //simple structure for sorting by vacancy count
  struct SiteVacancyCount {
    uint site;
    uint vacancy_count;
    bool operator<(const SiteVacancyCount& other) const;
  };

  struct StructureConfiguration {
    vector<POccSiteConfiguration> site_configs;
    double max_stoich_error;
    double max_site_error;
  };
} // namespace pocc

namespace pocc {
  struct POccSuperCell {
    //no need (YET) to make a class, simple ints and double, no real methods
    POccSuperCell();
    POccSuperCell(const POccSuperCell& b);
    const POccSuperCell& operator=(const POccSuperCell& b);
    bool operator<(const POccSuperCell& other) const;
    void free();
    void copy(const POccSuperCell& b);
    
    unsigned long long int hnf_index;
    unsigned long long int site_config_index;
    unsigned long long int degeneracy;
    double m_uff_energy;
  };
  bool sortPSCsUFFEnergy(const POccSuperCell& a, const POccSuperCell& b);
} // namespace pocc

namespace pocc {
  struct POccSuperCellSet{
    //no need (YET) to make a class, simple ints and double, no real methods
    const POccSuperCellSet& operator=(const POccSuperCellSet& b);
    bool operator<(const POccSuperCellSet& other) const;
    unsigned long long int getDegeneracy() const;
    const POccSuperCell& getSuperCell() const;
    double getHNFIndex() const;
    double getSiteConfigIndex() const;
    double getUFFEnergy() const;
    
    vector<POccSuperCell> psc_set;
  };
  bool sortPSCSetsUFFEnergy(const POccSuperCellSet& a, const POccSuperCellSet& b);
} // namespace pocc

namespace pocc{
  struct UFFParamAtom{
    //no need (YET) to make a class, simple ints and double, no real methods
    const UFFParamAtom& operator=(const UFFParamAtom& b);
    
    string symbol;
    double r1;        //bond distance
    double theta0;
    double x1;        //nonbond distance
    double D1;        //nonbond energy
    double zeta;      //scale
    double Z1;        //effective charge
    double Vi;
    double Uj;
    double ChiI;       //electronegativity
    double hard;
    double radius;
  };

  struct UFFParamBond{
    //no need (YET) to make a class, simple ints and double, no real methods
    const UFFParamBond& operator=(const UFFParamBond& b);
    void calculate(UFFParamAtom& uffp1,UFFParamAtom& uffp2,double distij);
    
    double ren;
    double R0;
    double Kij;
    double Xij;
    double Dij;
    double delta;
    double X6;
    double X12;
  };
} // namespace pocc

namespace pocc {
  class POccSiteConfiguration {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccSiteConfiguration();
      POccSiteConfiguration(int _site,int _i_hnf,vector<POccUnit>& pocc_groups);
      POccSiteConfiguration(const POccSiteConfiguration& b);
      //constructors - STOP

      ~POccSiteConfiguration();
      const POccSiteConfiguration& operator=(const POccSiteConfiguration& b);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - END
      
      //debug
      //vector<int> types_configuration_debug;                  //atom types, vacancy is -1
      //debug
   
      int site;                                         //reflexive
      int i_hnf;                                        //reflexive
      bool partial_occupation_flag;
      //any vector or xvector is over pocc_groups
      vector<POccUnit> m_pocc_groups;                  //reflexive
      xvector<int> xv_occupation_count_input;      //pre multiplication, from xstr_pocc
      xvector<int> xv_occupation_multiple;              //increment with each pocc_group
      xvector<int> xv_occupation_count_supercell;           //post multiplication with multiple
      xvector<double> xv_partial_occupation_value;
      xvector<double> xv_site_error;
      //sum of occupation_count_total and vacancy_count yields i_hnf (total sites)
      int occupation_count_total;
      int vacancy_count;
      double max_site_error;
      //double error_total;
      
      void prepareNoPOccConfig();
      void preparePOccConfig();
      int getNextOccupationMultiple(int i_hnf,xvector<int>& xv_occupation_multiple);
      int calculateOccupationCountTotal(xvector<int>& xv_next_occupation_multiple);
      void updateOccupationCounts(int _i_hnf, xvector<int> & xv_next_occupation_multiple);
      void calculateError();
      bool isPartiallyOccupied() const;
      vector<int> getStartingTypesConfiguration() const;

      //const vector<int>& getTypesConfiguration() const;
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccSiteConfiguration& b);
      //NECESSARY END CLASS METHODS - END
  };
} // namespace pocc

namespace pocc {
  class POccCalculatorTemplate {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccCalculatorTemplate();
      //constructors - STOP
      ~POccCalculatorTemplate();
      //NECESSARY PUBLIC CLASS METHODS - END
      
      xstructure xstr_pocc;                   //input from PARTCAR
      aurostd::xoption m_p_flags;             //e.g., vpflow
      _aflags m_aflags;                       //standard aflow flags
      xvector<double> stoich_each_type;       //converting deque<double> to xvector<double>
      xstructure xstr_nopocc;                 //will contain symmetry objects (_sym_op, pgroups most important here)
      vector<uint> types2pc_map;              //list of atom indices where types2pc_map(0) is 1st type 0 atom, types2pc_map(1) is 1st type 1 atom, etc.
      vector<string> m_species_redecoration;  //species used to redecorate xstr's for consistent symmetry/uff energies
      
      void setPOccFlags(const aurostd::xoption& pocc_flags);
      void setAFlags(const _aflags& Aflags);                      //standard _aflags
      void setPOccStructure(const xstructure& xstr_pocc);
      void setNonPOccStructure(const xstructure& xstr_nonpocc);
      void setSpeciesRedecoration(const vector<string>& species_redecoration);

    protected:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccCalculatorTemplate& b);
      //NECESSARY PRIVATE CLASS METHODS - END
  };
} // namespace pocc

namespace pocc {
  vector<uint> getVacanciesSuperCell(const vector<int>& pc2sc_map,const vector<vector<int> >& v_types_config);
  void replaceRandomSitesSuperCell(const xstructure& xstr_pocc,const vector<uint>& types2pc_map,const vector<int>& pc2sc_map,const vector<vector<int> >& v_types_config,xstructure& supercell);
  void rebuildSuperCell(const xstructure& xstr_pocc,const vector<uint>& v_vacancies,xstructure& supercell);
} // namespace pocc

namespace pocc { 
  class POccUFFEnergyAnalyzer: public POccCalculatorTemplate,xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccUFFEnergyAnalyzer(ostream& oss=cout);
      POccUFFEnergyAnalyzer(const aurostd::xoption& pocc_flags,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const _aflags& aflags,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss=cout);
      POccUFFEnergyAnalyzer(ofstream& FileMESSAGE,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const _aflags& aflags,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,ofstream& FileMESSAGE,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccUFFEnergyAnalyzer(const POccUFFEnergyAnalyzer& b);
      //POccUFFEnergyAnalyzer(xmatrix<double>*& _hnf_mat,xstructure*& _xstr_nopocc);
      //constructors - STOP
      ~POccUFFEnergyAnalyzer();
      const POccUFFEnergyAnalyzer& operator=(const POccUFFEnergyAnalyzer& b);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - END
      
      //underlying data structures
      bool m_initialized;
      xmatrix<double> hnf_mat;
      vector<vector<int> > m_types_config;            //the config for which we determined bonding
      vector<UFFParamAtom> types2uffparams_map;
      vector<uint> m_vacancies;
      double m_exploration_radius;
      xmatrix<double> distance_matrix;                    //references xstr_nopocc
      vector<double> v_dist_nn;                           //references xstr_nopocc
      xstructure xstr_ss;       //superstructure
      vector<int> sc2pc_map;
      vector<int> pc2sc_map;
      
      void initialize(ostream& oss=cout);
      void initialize(const aurostd::xoption& pocc_flags,ostream& oss=cout); 
      void initialize(const _aflags& aflags,ostream& oss=cout);
      void initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss=cout);
      void initialize(ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const _aflags& aflags,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const xstructure& xstr_nopocc,const vector<string>& species_redecoration,const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);

      void setSpeciesRedecoration(const vector<string>& species_redecoration);
      void setExplorationRadius();
      void getCluster(xmatrix<double>& _hnf_mat);
      void setBonds(vector<vector<int> >& v_types_config);
      double getUFFEnergy();
      bool isVacancy(vector<uint>& v_vacancies,uint atom);
      double bondEnergyBond(const UFFParamBond& uffb);
      double bondEnergyNoBond(const UFFParamBond& uffb);
      double getUFFBondEnergy(xstructure& xstr,vector<vector<uint> >& v_bonded_atom_indices,uint MODE);
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccUFFEnergyAnalyzer& b);
      //NECESSARY PRIVATE CLASS METHODS - END

      //for bonding, we need to create a super-superstructure (cluster)
      xstructure xstr_cluster;                                //cluster of atoms within radius
      vector<vector<uint> > v_bonded_atom_indices;            //references xstr_cluster
      vector<vector<uint> > v_nonbonded_atom_indices;         //references xstr_cluster

      bool has_vacancies;                                     //are there vacancies present in m_types_config?
      bool bonding_set;                                       //have we already found bonding for this configuration?
      double m_uff_energy;

      uint NNDistancesMapPC(uint atom);
      uint NNDistancesMapSC(uint atom);
      void calculateNNDistances(xstructure& xstr,vector<uint>& v_vacancies);
  };
} // namespace pocc

namespace pocc {
  class POccCalculator : public POccCalculatorTemplate,xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      POccCalculator(ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss=cout);
      POccCalculator(ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const aurostd::xoption& pocc_flags,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ostream& oss=cout);
      POccCalculator(const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      POccCalculator(const POccCalculator& b);
      //constructors - STOP
      ~POccCalculator();
      const POccCalculator& operator=(const POccCalculator& b);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - END

      //inputs
      bool m_initialized;
      _kflags m_kflags;                         //standard aflow flags
      xstructure xstr_sym;                    //will contain symmetry objects (_sym_op, pgroups most important here)
      int n_hnf;
      vector<POccUnit> m_pocc_sites;                   //groupings of atoms that are on the same site, non-vacant sites only, relative to xstr_nopocc
      int pocc_atoms_total;
      vector<StructureConfiguration> v_str_configs;

      POccUFFEnergyAnalyzer energy_analyzer;
      unsigned long long int hnf_count;
      unsigned long long int types_config_permutations_count;
      unsigned long long int total_permutations_count;
      std::list<POccSuperCellSet> l_supercell_sets;
      //standard flags - ALL options will be handled via xoptions
      aurostd::xoption enumerator_mode;       //how do we determine duplicates - UFF, SNF, ...
      
      void initialize(ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const _aflags& aflags,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const _kflags& kflags,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ostream& oss=cout);
      void initialize(ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const aurostd::xoption& pocc_flags,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ostream& oss=cout);
      void initialize(const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);
      void initialize(const xstructure& xstr_pocc,const aurostd::xoption& pocc_flags,const _aflags& aflags,const _kflags& kflags,ofstream& FileMESSAGE,ostream& oss=cout);

      //external methods
      void setPOccFlags(const aurostd::xoption& pocc_flags);  //input flags, e.g., vpflow
      void setPOccStructure(const xstructure& xstr_pocc);
      void setKFlags(const _kflags& Kflags);                      //standard _kflags
      
      xstructure createXStructure(const POccSuperCell& psc,int n_hnf=0,unsigned long long int hnf_count=0,unsigned long long int types_config_permutations_count=0,bool clean_structure=false,bool primitivize=false);

      bool areEquivalentStructuresByUFF(std::list<POccSuperCellSet>::iterator it, const POccSuperCell& psc) const;
      void add2DerivativeStructuresList(const POccSuperCell& psc,std::list<POccSuperCellSet>::iterator i_start,std::list<POccSuperCellSet>::iterator i_end);
      void add2DerivativeStructuresList(const POccSuperCell& psc);
      void getHNFMatSiteConfig(const POccSuperCell& psc,xmatrix<double>& hnf_mat,vector<vector<int> >& v_types_config);
      bool areEquivalentStructures(const POccSuperCell& psc_a,const POccSuperCell& psc_b);
      bool areEquivalentStructures(const xstructure& a,const xstructure& b);
      unsigned long long int runRobustStructureComparison(std::list<POccSuperCellSet>::iterator it);
      void calculateHNF();
      void getTotalPermutationsCount();
      void calculate();
      string getARUNString(unsigned long long int i);
      xstructure getUniqueSuperCell(unsigned long long int i);
      vector<xstructure> getUniqueDerivativeStructures();
      unsigned long long int getUniqueSuperCellsCount() const;
      //bool printUniqueDerviativeStructures();
      //void resetMaxSLRadius();
      void resetHNFMatrices();
      void resetSiteConfigurations();
    private:
      //NECESSARY PRIVATE CLASS METHODS - START
      void free();
      void copy(const POccCalculator& b);
      //NECESSARY END CLASS METHODS - END
      
      //hnf matrices
      int a_start, b_start, c_start;
      int d_start, e_start, f_start;
      vector<xmatrix<double> > v_unique_superlattices;
      //double max_superlattice_radius;
      //configurations
      xmatrix<double> hnf_mat;
      vector<vector<int> > v_types_config;
      //vector<int> v_config_iterators;
      uint config_iterator;
      vector<uint> v_config_order;
      double m_uff_energy_tolerance;

      void initializePOccStructure();
      void cleanPOccStructure();
      void preparePOccStructure();  //do not write out sym stuff by default

      const xmatrix<double>& getLattice() const;
      const vector<_sym_op>& getPGroup() const;
      //const StructureConfiguration& getXStrCountConfiguration(uint i) const;
      
      //void calculateHNF();                   //get n_hnf
      
      //useful internal methods
      //vector<_sym_op> getPGroups();           //fetch pgroups of xstr_nopocc
      //uint getHNFCount();                     //get count of hnf matrices, refers to v_hnf
      //xmatrix<double> getHNFMatrix(uint i);   //fetch specific hnf matrix, refers to v_hnf
      //_atom getAtom(uint i);                  //grab specific atom, refers to xstr_pocc
      
      //table stuff
      //set some nice printing precisions and paddings, mostly definitions
      uint getHNFTabelPOCCPrecision() const;
      uint getHNFTableGeneralPrecision() const;
      uint getHNFTableIterationPadding() const;
      uint getHNFTableErrorPadding() const;
      uint getHNFTableColumnPadding() const;
      string getHeaderMaxSiteError() const;
      string getHeaderMaxStoichError() const;
      string getHeaderStoichiometry() const;
      
      string hnfTableHeader();
      void writeHNFTableOutput(int i_hnf,double& stoich_error,double& site_error);
      string hnfTableLineOutput(int i_hnf,int str_config);
      //void setHNFTablePadding(int _AFLOW_POCC_PRECISION_);
      
      void partitionPOccSites();              //get pocc_sites
      xvector<double> calculateStoichEachType(vector<vector<int> >& v_types_config);
      void calculateSymNonPOccStructure(bool verbose=true);          //calculate symmetry of non-pocc structure
      void propagateEquivalentAtoms2POccStructure();  //propagates equivalent atoms info from xstr_sym to xstr_pocc for sorting
      void redecorateXStructures();            //redecorate xstr_nopocc for standardization of UFF energies
      
      //void getSiteCountConfigurations(int i_hnf,double& stoich_error);
      void getSiteCountConfigurations(int i_hnf);
      void getOptimizedSiteCountConfigurations(int site,int i_hnf,vector<POccSiteConfiguration>& v_site_configs);



      //determines site occupancy and vacancy count, given n_hnf
      bool iterateHNFMatrix();                //calculate all unique hnf's
      void setConfigOrder();
      bool getNextSiteConfiguration();
      bool getNextSiteConfiguration(vector<vector<int> >& v_site_config);
  };

}

#endif  // _AFLOW_POCC_H_

// ***************************************************************************
// *                                                                         *
// *         aflow - Automatic FLOW for materials discovery project          *
// *             Stefano Curtarolo - Duke University - 2003-2019             *
// *                 Kesong Yang - Duke University 2010-2011                 *
// *                  Corey Oses - Duke University 2013-2019                 *
// *                                                                         *
// ***************************************************************************
