// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                                                                         *
// ***************************************************************************
// MAKEFILE FOR AFLOW_APL
// Written by Michal Jahnatek

#ifndef _AFLOW_APL_H_
#define _AFLOW_APL_H_

//#define USE_MKL 1

// Almost generally used precision in the whole apl, however, somewhere it is
// hard coded based on the tests and it work much better...
#define _AFLOW_APL_EPS_ 1E-6
extern bool _WITHIN_DUKE_;  //will define it immediately in kphonons

// Aflow core libraries
#include "../aflow.h"
#include "../AUROSTD/aurostd.h"
//The functions in the header file "aflow_qha_operations.h" include various vector and matrix operations.//
//These functions are used to calculate eigenvalues and eigenvectors of complex symmetric Hermitian matrices and .//
//Also perform other operations involved in calculating quasiharmonic properties.//
//These functions have been built based on vectorview concepts//
#include "aflow_qha_operations.h"
//#define _AFLOW_APL_REGISTER_ register   // register ?
#define _AFLOW_APL_REGISTER_

// Create the version of GCC, we will uset it for multithread parts of code,
// to check if the current compiling version of gcc is able to compile the
// std::thead features
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif

// ME20190219 - Define the checksum algorithm used for APL hibernate files
#define APL_CHECKSUM_ALGO  string("Fletcher32")

// Basic objects ...
// ***************************************************************************
// "aplexcept.h"
// in aurostd.h // [OBSOLETE]
// ME20200222 - APLStageBreak obsolete
//[OBSOLETE] #include <stdexcept>
//[OBSOLETE] namespace apl {
//[OBSOLETE]   //
//[OBSOLETE]   // OBSOLETE ME20191031 - use xerror
//[OBSOLETE]   //class APLRuntimeError : public std::runtime_error {
//[OBSOLETE]   // public:
//[OBSOLETE]   //  APLRuntimeError(const std::string& s) : std::runtime_error(s) {}
//[OBSOLETE]   //};
//[OBSOLETE]   //class APLLogicError : public std::logic_error {
//[OBSOLETE]   // public:
//[OBSOLETE]   //  APLLogicError(const std::string& s) : std::logic_error(s) {}
//[OBSOLETE]   //};
//[OBSOLETE]   //
//[OBSOLETE]   class APLStageBreak : public std::exception {
//[OBSOLETE]     public:
//[OBSOLETE]       APLStageBreak() {}
//[OBSOLETE]   };
//[OBSOLETE] }

// ***************************************************************************
// "logger.h"
namespace apl {
  class Logger;
  // Templates used for the definition of the one parameter manipulation
  // functions of Logger's streams
  template <class T>
    class LMANIP;
  template <class T>
    Logger& operator<<(Logger&, const LMANIP<T>&);
  template <class T>
    class LMANIP {
      public:
        LMANIP(Logger& (*a)(Logger&, T), T v) {
          _action = a;
          _value = v;
        }
        friend Logger& operator<<<>(Logger&, const LMANIP<T>&);

      private:
        Logger& (*_action)(Logger&, T);
        T _value;
    };
  template <class T>
    Logger& operator<<(Logger& s, const LMANIP<T>& m) {
      return (*m._action)(s, m._value);
    }
  // Function without parameter...
  typedef Logger& (*MyStreamManipulator)(Logger&);
  // this is the type of std::cout
  typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
  // this is the function signature of std::endl
  typedef CoutType& (*StandardEndLine)(CoutType&);
  //
  class Logger {
    private:
      ofstream* _os;
      std::string _barCode;
      std::string _typeofMessage;
      std::string _moduleName;
      _aflags _aflowFlags;
      stringstream _ss;
      int _progressBarLastPercent;
      double _progressBarPercent;  // ME20180831
      bool _isQuiet;

    public:
      Logger();
      Logger(std::ofstream&, const _aflags&);
      void initialize(std::ofstream&, const _aflags&);
      Logger(Logger&);
      ~Logger();
      ofstream& getOutputStream();
      void initProgressBar(const char*);
      void initProgressBar(const string&);
      void updateProgressBar(double);  // ME20180831
      void updateProgressBar(int, int);
      void finishProgressBar();
      void setTypeOfMessage(const string&);
      void setBarCode(const string&);
      void setBarCode(const char&);
      void setModuleName(const string&);
      void setQuietMode(bool);
      Logger& operator=(Logger&);
      Logger& operator<<(const string& s);
      Logger& operator<<(const int& s);
      Logger& operator<<(const uint& s);
      Logger& operator<<(const double& s);
      // Manipulation functions without parameter
      Logger& operator<<(StandardEndLine);
      Logger& operator<<(MyStreamManipulator);
      void endl();
      friend Logger& endl(Logger&);
      friend Logger& message(Logger&);
      friend Logger& notice(Logger&);
      friend Logger& warning(Logger&);
      friend Logger& error(Logger&);
      friend Logger& quiet(Logger&);
      // Manipulation functions with one parameter
      friend Logger& setwidth(Logger&, int);
      friend Logger& setformat(Logger&, const char*);
  };
  // Friend functions without parameters
  Logger& endl(Logger&);
  Logger& message(Logger&);
  Logger& notice(Logger&);
  Logger& warning(Logger&);
  Logger& error(Logger&);
  Logger& quiet(Logger&);
  // Friend functions with one parameter
  Logger& setwidth(Logger&, int);
  LMANIP<int> sw(int);
  Logger& setformat(Logger&, const char*);
  LMANIP<const char*> sf(const char*);
}  // namespace APL

// ***************************************************************************
// "hroutines.h"
// in aurostd.h // [OBSOLETE]
#include <typeinfo>

namespace apl {
  template <typename T>
    inline std::string stringify(const T& x) {
      std::ostringstream o;
      if (!(o << x)) {
        // ME20191031 - use xerror
        //throw APLRuntimeError(std::string("stringify(") + typeid(x).name() + ")");
        throw aurostd::xerror(_AFLOW_FILE_NAME_, "apl::stringify()", std::string("stringify(") + typeid(x).name() + ")");
      }
      return o.str();
    }
  // ME20190219 - BEGIN
  //[OBSOLETE] vector<vector<int> > getThreadDistribution(const int&, const int&);  // ME20180801
  void tokenize(const string& strin, vector<string>& tokens, const string& del);
  //[OBSOLETE] string getVASPVersionString(const string&);
  //[OBSOLETE] unsigned int getFileCheckSum(const string&);
  //[OBSOLETE] void printXVector(const xvector<double>&, bool = true);
  //[OBSOLETE] void printXVector(const xvector<xcomplex<double> >&);
  //[OBSOLETE] void printXMatrix(const xmatrix<double>&);
  //[OBSOLETE] void printXMatrix(const xmatrix<xcomplex<double> >&);
  //[OBSOLETE] void printXMatrix2(const xmatrix<xcomplex<double> >&);
  // ME20190219 - END
}  // namespace apl

// ***************************************************************************
// ***************************************************************************
// BEGIN ME: Anharmonic Force Constants (AAPL)
// ***************************************************************************

namespace apl {
  // Options for anharmonic IFC calculations
  // OBSOLETE ME20190501 - Replaced with xoption
  //struct _anharmonicIFCOptions {
  //  int max_iter;
  //  double sumrule_threshold;
  //  double mixing_coefficient;
  //};

  // _cluster holds a single cluster
  struct _cluster {
    vector<int> atoms;  // List of atoms inside the cluster
    int fgroup;  // Index pointing to the factor group that transforms the cluster into another cluster
    int permutation;  // Index pointing to the permutation that transforms the cluster into another cluster
  };

  // _ineq_distortions contains a list of inequivalent distortion and its equivalent
  // distortions for a given set of atoms
  struct _ineq_distortions {
    vector<int> atoms;  // A list of atoms involved in the distortions
    vector<int> clusters;  // A list of cluster sets that use these distortions for their force constant calculations
    vector<vector<vector<int> > > distortions; // Map of distortions. The distortions vectors need to be defined elsewhere. 
    vector<vector<int> > rotations;  // The factor group that holds the rotation to transform the distortions
    vector<vector<vector<int> > > transformation_maps;  // A map containing the transformation of the atoms for each distortion
  };

  // _linearCombinations is a structure to store linear combinations.
  struct _linearCombinations {
    vector<vector<int> > indices;  // Cartesian indices of each linear combination
    vector<vector<double> > coefficients;  // Coefficients of each linear combination
    vector<int> independent;  // The linearly independent values
    vector<int> dependent;  // The linearly dependent values
    // indep2depMap maps the independent coefficients to the coefficients that
    // depend on them. This is used for the IFC correction method.
    vector<vector<int> > indep2depMap;
  };

  class Supercell;  // Forward declaration
  class ClusterSet {
    // See aflow_aapl_cluster.cpp for detailed descriptions of the functions
    public:
      ClusterSet();
      ClusterSet(const Supercell&, int, double, ofstream&, _aflags&);  // Constructor
      ClusterSet(const string&, const Supercell&, int, double, int, ofstream&, _aflags&);  // From file
      ClusterSet(const ClusterSet&);  // Constructor from another ClusterSet instance
      ~ClusterSet();  // Destructor
      const ClusterSet& operator=(const ClusterSet&);  // Copy constructor
      void clear(ofstream&, _aflags&);
      void initialize(const Supercell&, int, double);

      vector<_cluster> clusters;
      vector<vector<int> > coordination_shells;  // Contains all coordinate shells. Central atoms is index 0.
      double cutoff;  // Cutoff radius in Angstroms
      vector<xvector<double> > distortion_vectors;  // List of distortion vectors
      vector<_ineq_distortions> higher_order_ineq_distortions;  // ME20190531 - for 3rd derivatives of higher order processes
      vector<vector<int> > ineq_clusters;  // Clusters rearranged into sets of equivalent clusters.  // ME20190520
      vector<_ineq_distortions> ineq_distortions; // List of inequivalent distortions
      vector<_linearCombinations> linear_combinations;  // List of linear combinations of the IFCs
      int nifcs;  // Number of force constants for each set of atoms.
      int order;  // Order of the cluster, i.e. the order of the force constant to be calculated.
      xstructure pcell;  // Structure of the primitive cell.
      vector<int> pc2scMap;  // Atom map from the primitive cell to the supercell.
      vector<vector<int> > permutations;  // List of possible permutations for the cluster
      xstructure scell;  // Structure of the supercell.
      vector<int> sc2pcMap;  // Atom map from the supercell to the primitive cell.
      xvector<int> sc_dim;  // Dimensions of the supercell.
      vector<vector<int> > symmetry_map;  // Symmetry atom map for the atoms in the clusters

      const _cluster& getCluster(const int& i) const;  // ME20190520
      void build(int);
      void buildDistortions();
      void writeClusterSetToFile(const string&);

    private:
      ofstream* messageFile;
      _aflags* aflags;

      void free();
      void copy(const ClusterSet&);

      double getMaxRad(const xstructure&, int);
      void buildShells();
      vector<_cluster> buildClusters();
      vector<vector<int> > getSymmetryMap();
      vector<vector<int> > getPermutations(int);

      // Clusters
      void getInequivalentClusters(vector<_cluster>&, vector<vector<int> >&);
      int getNumUniqueAtoms(const vector<int>&);
      vector<int> getComposition(const vector<int>&);
      bool sameComposition(const vector<int>&, const vector<int>&);
      int equivalenceCluster(const vector<int>&, const vector<int>&,
          const vector<vector<int> >&, const vector<vector<int> >&);
      vector<int> translateToPcell(const vector<int>&, int);
      int comparePermutations(const vector<int>&, const vector<int>&);
      bool atomsMatch(const vector<int>&, const vector<int>&, const vector<int>&, const int&);
      void getSymOp(_cluster&, const vector<int>&);

      // Distortions
      vector<xvector<double> > getCartesianDistortionVectors();
      vector<_ineq_distortions> initializeIneqDists();
      int sameDistortions(const _cluster&, const vector<_ineq_distortions>&);
      vector<vector<int> > getTestDistortions(const vector<int>&);
      void getInequivalentDistortions(const vector<vector<int> >&, _ineq_distortions&);
      void appendDistortion(_ineq_distortions&, vector<int>,
          const int& eq=-1, const int& fg=-1);
      bool allZeroDistortions(const vector<int>&, const vector<int>&);
      bool allAtomsMatch(const int&, const vector<int>&);
      int equivalenceDistortions(const xmatrix<double>&, const vector<int>&,
          const vector<vector<vector<int> > >&, const vector<int>&);
      vector<int> getTransformationMap(const int&, const int&);
      vector<_ineq_distortions> getHigherOrderDistortions();

      // Linear Combinations
      vector<_linearCombinations> getLinearCombinations();
      vector<vector<int> > getInvariantSymOps(const _cluster&);
      vector<vector<double> > buildCoefficientMatrix(const vector<vector<int> >&);
      vector<vector<double> > getRREF(vector<vector<double> >);

      // File I/O
      string writeParameters();
      string writeInequivalentClusters();
      string writeClusters(const vector<_cluster>&);
      string writeLinearCombinations(const _linearCombinations&);
      string writeInequivalentDistortions();
      string writeIneqDist(const _ineq_distortions&);
      string writeHigherOrderDistortions();

      void readClusterSetFromFile(const string&);
      bool checkCompatibility(uint&, const vector<string>&);
      void readInequivalentClusters(uint&, const vector<string>&);
      vector<_cluster> readClusters(uint&, const vector<string>&);
      _linearCombinations readLinearCombinations(uint&, const vector<string>&);
      void readInequivalentDistortions(uint&, const vector<string>&);
      _ineq_distortions readIneqDist(uint&, const vector<string>&);
      void readHigherOrderDistortions(uint&, const vector<string>&);
  };

  class AnharmonicIFCs {
    // See aflow_aapl_ifcs.cpp for detailed descriptions of the functions
    public:
      AnharmonicIFCs();
      AnharmonicIFCs(_xinput&, _aflags&, _kflags&, _xflags&, ClusterSet&, ofstream&);
      AnharmonicIFCs(const AnharmonicIFCs&);
      const AnharmonicIFCs& operator=(const AnharmonicIFCs&);
      ~AnharmonicIFCs();
      void clear(_xinput&,_aflags&, _kflags&, _xflags&, ClusterSet&, ofstream&);

      void setOptions(double, int, double, double, bool);
      int getOrder() const;

      bool runVASPCalculations(bool);
      bool calculateForceConstants();
      const vector<vector<double> >& getForceConstants() const;
      vector<vector<int> > getClusters() const;
      void writeIFCsToFile(const string&);

    private:
      _xinput* _xInput;
      _aflags* _aflowFlags;
      _kflags* _kbinFlags;
      _xflags* _xFlags;
      ClusterSet* clst;
      ofstream* messageFile;
      bool new_ofstream;

      vector<_xinput> xInputs;
      bool _useZeroStateForces;
      vector<vector<int> > cart_indices;  // A list of all Cartesian indices
      double distortion_magnitude;  // The magnitude of the distortions in Angstroms
      vector<vector<double> > force_constants;  // Symmetrized IFCs - ME20190520
      int max_iter;  // Number of iterations for the sum rules
      double mixing_coefficient;  // The mixing coefficient for the SCF procedure
      int order;  // The order of the IFCs
      double sumrule_threshold;  // Convergence threshold for the sum rules

      void free();
      void copy(const AnharmonicIFCs&);

      string buildRunName(const vector<int>&, const vector<int>&, int, int);
      void applyDistortions(_xinput&, const vector<xvector<double> >&, const vector<int>&, const vector<int>&, double=1.0);

      vector<vector<int> > getCartesianIndices();

      vector<vector<vector<xvector<double> > > > storeForces(vector<_xinput>&);
      vector<vector<xvector<double> > > getForces(int, int&, vector<_xinput>&);
      int getTransformedAtom(const vector<int>&, const int&);
      void addHigherOrderForces(vector<vector<vector<xvector<double> > > >&, int&, vector<_xinput>&);
      vector<vector<double> > calculateUnsymmetrizedIFCs(const vector<_ineq_distortions>&,
          const vector<vector<vector<xvector<double> > > >&);
      double finiteDifference(const vector<vector<xvector<double> > >&, int,
          const vector<int>&, const vector<int>&);

      // Symmetrization Functions
      vector<vector<double> > symmetrizeIFCs(vector<vector<double> >);
      typedef vector<std::pair<vector<int>, vector<double> > > tform;
      typedef vector<vector<vector<vector<int> > > > v4int;
      void getTensorTransformations(v4int&, vector<vector<tform> >&);
      vector<vector<int> > getReducedClusters();
      void applyLinCombs(vector<vector<double> >&);
      void transformIFCs(const vector<vector<tform> >&, vector<vector<double> >&);
      void applyPermutations(vector<int>, vector<vector<double> >&);
      void calcSums(const vector<vector<int> >&, const vector<vector<double> >&,
          vector<vector<double> >&, vector<vector<double> >&);
      void correctIFCs(vector<vector<double> >&, const vector<vector<double> >&,
          const vector<vector<double> >&,
          const vector<vector<int> >&, const v4int&);
      vector<double> getCorrectionTerms(const int&,
          const vector<vector<int> >&,
          const vector<vector<double> >&,
          const vector<vector<double> >&,
          const vector<vector<double> >&);
      uint findReducedCluster(const vector<vector<int> >&, const vector<int>&);

      // File I/O
      string writeParameters();
      string writeIFCs();
      bool checkCompatibility(uint&, const vector<string>&);
      vector<vector<double> > readIFCs(uint&, const vector<string>&);
  };

}  //namespace apl
// ***************************************************************************
// END ME: Anharmonic Force Constants (AAPL)
// ***************************************************************************

// ***************************************************************************
// Linear and nonlinear fitting functions, this functions compute fitting parameters correctly
//although chi-square matrix is not computing correctly [PINKU]
#define FIT_LIMIT 0.001
namespace apl {
  class aflowFITTING {
    public:
      aflowFITTING() {}
      ~aflowFITTING() { clear(); }
      void clear() {}

    private:
      double chisq_quadratic;
      double chisq_birch_murnaghan;
      uint Iteration_birch_murnaghan;
      double alamda_birch_murnaghan;
      vector<double> Uncertainties_birch_murnaghan;
      //user defined fitting functions
      void funcs(double x, xvector<double>& afunc);
      void birch_murnaghan_function(double x,
          const xvector<double> a,
          double* y,
          xvector<double>& dyda);

      //linear leatsquare fitting function
      template <class utype>
        bool lfit(xvector<utype>& x,
            xvector<utype>& y,
            xvector<utype>& sig,
            xvector<utype>& a,
            xvector<int>& ia,
            xmatrix<utype>& covar,
            utype& chisq,
            void (aflowFITTING::*funcs)(utype, xvector<utype>&));

      //nonlinear leatsquare fitting function
      bool mrqmin(xvector<double> x,
          xvector<double> y,
          xvector<double>& sig,
          xvector<double>& a,
          xvector<int>& ia,
          xmatrix<double>& covar,
          xmatrix<double>& alpha,
          double& chisq,
          void (aflowFITTING::*birch_murnaghan_function)(double, xvector<double>, double*, xvector<double>&),
          double* alamda);

      template <class utype>
        void covsrt(xmatrix<utype>& covar, xvector<int>& ia, int& mfit);
      template <class utype>
        bool gaussj(xmatrix<utype>& a, int& n, xmatrix<utype>& b, int m);

      void mrqcof(xvector<double> x,
          xvector<double> y,
          xvector<double>& sig,
          xvector<double>& a,
          xvector<int>& ia,
          xmatrix<double>& alpha,
          xvector<double>& beta,
          double& chisq,
          void (aflowFITTING::*birch_murnaghan_function)(double, xvector<double>, double*, xvector<double>&));

    public:
      bool
        quadraticfit(xvector<double> energy, xvector<double> volume, xvector<double>& params);
      bool
        birch_murnaghan_fitting(xvector<double> energy, xvector<double> volume, xvector<double> guess, xvector<double>& params);
      double get_chisq_quadratic() { return chisq_quadratic; }
      double get_chisq_birch_murnaghan() { return chisq_birch_murnaghan; }
      uint getIteration_birch_murnaghan() { return Iteration_birch_murnaghan; }
      double getalamda_birch_murnaghan() { return alamda_birch_murnaghan; }
      vector<double> getUncertainties_birch_murnaghan() { return Uncertainties_birch_murnaghan; }
  };
}
// ***************************************************************************
namespace apl {
#ifndef EIGEN_H
#define EIGEN_H
  //Functions in this class calculate eigenvalues and eigenvectors of complex symmetrix
  //nxn matrices and can sort them according to following options
  //eval and evec SORTING OPTIONS
  // [ 1. [APL_MV_EIGEN_SORT_VAL_ASC]  => ascending order   ]
  // [ 2. [APL_MV_EIGEN_SORT_VAL_DESC] => descending order ]
  // [ 3. [APL_MV_EIGEN_SORT_ABS_ASC]  => absolute ascending order ]
  // [ 4. [APL_MV_EIGEN_SORT_ABS_DESC] => absolute descending order]
  //Compute Eigenvalues and Eigenvectors for (nxn) Complex Hermitian matrices
  class aplEigensystems : public MVops {
    public:
      aplEigensystems() {}
      ~aplEigensystems() {}
      void clear() { } //this->clear(); //CO200106 - patching for auto-indenting
      void eigen_calculation(const aurostd::xmatrix<xcomplex<double> >& M, aurostd::xvector<double>& eval, aurostd::xmatrix<xcomplex<double> >& evec);
      void eigen_calculation(const xmatrix<xcomplex<double> >& M, xvector<double>& eval, xmatrix<xcomplex<double> >& evec, apl_eigen_sort_t t);
  };
#endif
}
// ***************************************************************************
namespace apl {
  //The functions in this class are used to perform linear and nonlinear fittings//
  //User defined functions can be used to perform fitting//
  //Both Levenberg-Marquardt and unscaled Levenberg-Marquardt algorithms can be performed
  //either with keyword "APL_multifit_fdfsolver_lmsder" or "APL_multifit_fdfsolver_lmder" //
#ifndef FIT_H
#define FIT_H
#define spread 1e-6               //accuracy of data fitting
#define allowed_fit_error 100.00  // if spread is small then allowed_fit_error is relatively take high value
#define Absolute_Error 1E-6
#define Relative_Error 1E-6
#define APL_multifit_fdfsolver_lmsder  //This is a robust and efficient version of the Levenberg-Marquardt algorithm
#undef APL_multifit_fdfsolver_lmder    //This is an unscaled version of the Levenberg-Marquardt algorithm

  class md_lsquares : public MVops {
    // multi-parameter linear regression
    // multidimensional nonlinear least-squares fitting
    public:
      md_lsquares() {}
      ~md_lsquares() {}

      //user define data sets
      vector<double> Xdata;  //usder defined input data
      vector<double> Ydata;  //usder defined input data
      //guess values calculated automatically
      vector<double> guess;  // guess values for nonlinear fitting, It will calculate automatically.
      //error managments
      bool data_read_error;   //return true if Xdata.size and Ydata.size are same
      int nl_success_status;  //return error message interms of int numbers
      string nl_err_msg;      //return type of error mesages

      //lfit output
      double luncertanity_V0;  //uncertainties linear fit(Quadratic function) in Volume
      double luncertanity_E0;  //uncertainties linear fit(Quadratic function) in Energies
      double luncertanity_B0;  //uncertainties linear fit(Quadratic function) in Bulk Modulus
      double lchisq;           //linear fit chisquare
      double leqmV0;           //Equilibrium Volume from linear fit (Quadratic function)
      double leqmE0;           //Equilibrium Energy from linear fit (Quadratic function)
      double leqmB0;           //Equilibrium Bulk Modulus from linear fit (Quadratic function)
      //nlfit output
      double uncertanity_V0;   //uncertainties non linear fit(Birch Murnaghan Function) in Volume
      double uncertanity_E0;   //uncertainties non linear fit(Birch Murnaghan Function) in Energy
      double uncertanity_B0;   //uncertainties non linear fit(Birch Murnaghan Function) in Bulk Modulus
      double uncertanity_Bp;   //uncertainties non linear fit(Birch Murnaghan Function) in Pressure derivatives Bulk Modulus
      double uncertanity_Bpp;  //uncertainties non linear fit(Birch Murnaghan Function) in Pressure derivatives Bp
      double chisq_dof;        //chi square per degrees of freedom
      double nleqmV0;          //Equilibrium Volume from nonlinear fit (Birch Murnaghan Function)
      double nleqmE0;          //Equilibrium Energy from nonlinear fit (Birch Murnaghan Function)
      double nleqmB0;          //Equilibrium Bulk Modulus from nonlinear fit (Birch Murnaghan Function)
      double nleqmBp;          //Equilibrium Pressure derivatives Bulk Modulus from nonlinear fit (Birch Murnaghan Function)
      double nleqmBpp;         //Equilibrium Pressure derivatives Bp from nonlinear fit (Birch Murnaghan Function)
      string fdfsolver_name;   //return nonlinear method name
      //fitting functions
      void cubic_polynomial_fit();                                            //cubic linear function fitting
      void birch_murnaghan_fit();                                             //Birch-Murnaghan nonlinear function fitting
      void birch_murnaghan_4th_order_fit(const xvector<double>& user_guess);  //Birch-Murnaghan 4th order nonlinear function fitting
      void birch_murnaghan_3rd_order_fit(const xvector<double>& user_guess);  //Birch-Murnaghan 3rd order nonlinear function fitting
      void clear();
  };
#endif
}
// ***************************************************************************
// ***************************************************************************
// Engine core for phonon calculation
// ***************************************************************************
// "supercell.h"
namespace apl {
  class Supercell {
    private:
      string _directory;  // for the logger
      ofstream* messageFile;
      xstructure _inStructure;
      xstructure _inStructure_original;  //CO
      xstructure _inStructure_light;     //CO, does not include HEAVY symmetry stuff
      //deque<_atom> _inStructure_atoms_original; //CO
      xstructure _pcStructure;           //CO20180406 - for the path
      xstructure _scStructure;
      xstructure _scStructure_original;  //CO
      xstructure _scStructure_light;     //CO, does not include HEAVY symmetry stuff
      //deque<_atom> _scStructure_atoms_original; //CO
      //CO - START
      bool _skew;                  //SYM::isLatticeSkewed(), same for pc and sc
      bool _derivative_structure;  //vs. simple expanded_lattice, derivative structure's lattice has LESS symmetry, so be careful ApplyAtom()'ing
      double _sym_eps;             //same for pc and sc
      //CO - END
      bool _isShellRestricted;
      int _maxShellID;
      vector<double> _maxShellRadius;
      bool _isConstructed;
      vector<vector<vector<xvector<double> > > > phase_vectors;  // ME20200116

    private:
      void calculateWholeSymmetry(xstructure&, bool=true);
      xstructure calculatePrimitiveStructure() const;
      bool getMaps(const xstructure&, const xstructure&, const xstructure&, vector<int>&, vector<int>&);  // ME20200117
      void free();
      void copy(const Supercell&);

    public:
      Supercell();
      Supercell(const xstructure&, ofstream&, string="./"); //CO20181226
      Supercell(const string&, ofstream&, string="./");  // ME20200112
      Supercell(const Supercell&);
      Supercell& operator=(const Supercell&);
      ~Supercell();
      void clear(ofstream&);
      void setDirectory(const string&);
      string getDirectory() const;
      void readFromStateFile(const string&);  // ME20200212
      void initialize(const xstructure&, bool=true);  // ME20191225
      void clearSupercell();
      //void LightCopy(const xstructure& a, xstructure& b);  // OBSOLETE ME20200220
      bool isConstructed();
      void reset();
      xvector<int> determineSupercellDimensions(const aurostd::xoption&);  // ME20191225
      void build(aurostd::xoption&, bool = true);  // ME20191225
      void build(const xvector<int>&, bool = true);  // ME20191225
      void build(int, int, int, bool = TRUE);
      void trimStructure(int, const xvector<double>&,
          const xvector<double>&, const xvector<double>&,
          bool = true);
      bool projectToPrimitive();  // ME20200117
      void projectToOriginal();  // ME20200117
      xvector<int> buildSuitableForShell(int, bool, bool VERBOSE);  // ME20200102
      void setupShellRestrictions(int);
      // ME20190715 BEGIN - added const to getter functions so they can be used with const Supercell &
      bool isShellRestricted() const;
      int getMaxShellID() const;
      int getNumberOfAtoms() const;
      int getNumberOfUniqueAtoms() const;
      int getNumberOfEquivalentAtomsOfType(int) const; //CO20190218
      int getUniqueAtomID(int) const;
      int getUniqueAtomID(int, int) const;
      const _atom& getUniqueAtom(int) const;
      string getUniqueAtomSymbol(int) const;
      double getUniqueAtomMass(int) const;
      double getAtomMass(int) const;
      int getAtomNumber(int) const;
      // ME20190715 END
      const xstructure& getSupercellStructure() const;
      const xstructure& getSupercellStructureLight() const;
      const xstructure& getPrimitiveStructure() const;
      const xstructure& getInputStructure() const;
      const xstructure& getInputStructureLight() const;
      const xstructure& getOriginalStructure() const;  // ME20200117
      int atomGoesTo(const _sym_op&, int, int, bool = TRUE); //CO20190218
      int atomComesFrom(const _sym_op&, int, int, bool = TRUE); //CO20190218
      const _sym_op& getSymOpWhichMatchAtoms(int, int, int);
      xvector<double> getFPositionItsNearestImage(const xvector<double>&,
          const xvector<double>&,
          const xmatrix<double>&);
      xvector<double> getFPositionItsNearestImage(int, int);
      xvector<double> getCPositionItsNearestImage(int, int);
      bool compareFPositions(xvector<double>&, xvector<double>&);          //CO
      bool compareFPositions(xvector<double>&, xvector<double>&, double);  //CO
      void calculatePhaseVectors();  // ME20200117
      bool calcShellPhaseFactor(int, int, const xvector<double>&, xcomplex<double>&);
      bool calcShellPhaseFactor(int, int, const xvector<double>&, xcomplex<double>&,
          int&, xvector<xcomplex<double> >&, bool);  // ME20180828
      int pc2scMap(int) const;
      int sc2pcMap(int) const;
      void center(int);
      //CO - START
      void center_original(void);
      //corey
      void getFullBasisAGROUP();  // ME20191218
      bool fullBasisCalculatedAGROUP();  // ME20191218
      const vector<vector<_sym_op> >& getAGROUP(void) const;
      const vector<_sym_op>& getFGROUP(void) const;
      const vector<_sym_op>& getAGROUP(int) const;
      bool isDerivativeStructure() const;
      double getEPS(void) const;
      // ME20190715 - END
      // CO - END
      // **** BEGIN JJPR *****
      xvector<int> scell;
      vector<int> _pc2scMap;
      vector<int> _sc2pcMap;
      // **** END  JJPR *****
  };
}  // namespace apl

// ***************************************************************************
// "ipc.h"

namespace apl {
  enum IPCFreqFlags {
    NONE = 0L,
    ALLOW_NEGATIVE = 1L << 1,
    OMEGA = 1L << 2,
    RAW = 1L << 3,  // eV/A/A/atomic_mass_unit
    HERTZ = 1L << 4,
    THZ = 1L << 5,
    RECIPROCAL_CM = 1L << 6,
    MEV = 1L << 7
  };
  inline IPCFreqFlags operator&(const IPCFreqFlags& __a, const IPCFreqFlags& __b) {
    return IPCFreqFlags(static_cast<int>(__a) & static_cast<int>(__b));
  }
  inline IPCFreqFlags operator|(const IPCFreqFlags& __a, const IPCFreqFlags& __b) {
    return IPCFreqFlags(static_cast<int>(__a) | static_cast<int>(__b));
  }
  inline IPCFreqFlags operator|=(IPCFreqFlags& __a, const IPCFreqFlags& __b) {
    return (__a = (__a | __b));
  }
  // //////////////////////////////////////////////////////////////////////////

}  // namespace apl

// ***************************************************************************
// "phoncalc.h"
#define _AFLOW_APL_BORN_EPSILON_RUNNAME_ string("LRBE")  // ME20190108
#define _AFLOW_APL_BORN_EPSILON_DIRECTORY_NAME_ string(ARUN_DIRECTORY_PREFIX + "APL_" + _AFLOW_APL_BORN_EPSILON_RUNNAME_) // ME20190108
//#define _AFLOW_APL_FORCEFIELDS_RUNNAME_ string("LRFF")  // ME20190108  // OBSOLETE ME20200213 - the calculation does not use force fields
//#define _AFLOW_APL_FORCEFIELDS_DIRECTORY_NAME_ string(ARUN_DIRECTORY_PREFIX + "APL_" + _AFLOW_APL_FORCEFIELDS_RUNNAME_) // ME20190108  // OBSOLETE ME20200213
#define _AFLOW_APL_DFPT_RUNNAME_ string("DFPT")  // ME20200213
#define _AFLOW_APL_DFPT_DIRECTORY_NAME_ string(ARUN_DIRECTORY_PREFIX + "APL_" + _AFLOW_APL_DFPT_RUNNAME_) // ME20200213

namespace apl {
  bool createAflowInPhonons(const _aflags&, const _kflags&, const _xflags&, _xinput&); // ME20190108
  void createAflowInPhononsAIMS(_aflags&, _kflags&, _xflags&, string&, _xinput&, ofstream&);
  bool filesExistPhonons(_xinput&);
  bool outfileFoundAnywherePhonons(vector<_xinput>&);
  bool outfileFoundEverywherePhonons(vector<_xinput>&, const string&, ofstream&, bool=false);  // ME20191029
  bool readForcesFromDirectory(_xinput&);  // ME20200219
  void subtractZeroStateForces(vector<_xinput>&, bool);
  void subtractZeroStateForces(vector<_xinput>&, _xinput&);  // ME20190114
}

namespace apl {
  class ForceConstantCalculator {
    protected:
      // Aflow's stuff required for running some routines
      Supercell* _supercell;
      _xinput* _xInput; //_xvasp& _vaspRun;
      _aflags* _aflowFlags;
      _kflags* _kbinFlags;
      _xflags* _xFlags; //_vflags& _vaspFlags;
      string* _AflowIn;
      ofstream* messageFile;

      vector<_xinput> xInputs;
      string _AflowInFileName;

      // Calculate forces at no distortion - since for some structure
      // (not well relaxed, or with other problems) these forces have to be
      // known and substracted from the calculated forces with distortion
      bool _calculateZeroStateForces;

      // For each atom of supercell, there is a full force field
      vector<vector<xmatrix<double> > > _forceConstantMatrices;
      // Stuff for polar materials
      bool _isPolarMaterial;
      // For each atom there is a matrix 3x3 of Born effective charge
      vector<xmatrix<double> > _bornEffectiveChargeTensor;
      // Dielectric tensor
      xmatrix<double> _dielectricTensor;

    private:
      void free();
      void copy(const ForceConstantCalculator&);

      virtual bool calculateForceConstants() {return false;} // ME20200211

      void symmetrizeForceConstantMatrices();
      void correctSumRules();

      void printForceConstantMatrices(ostream&);
      void printFCShellInfo(ostream&);

    public:
      ForceConstantCalculator();
      ForceConstantCalculator(Supercell&, _xinput&, _aflags&, _kflags&, _xflags&, string&, ofstream&);
      ForceConstantCalculator(const ForceConstantCalculator&);
      ForceConstantCalculator& operator=(const ForceConstantCalculator&);
      virtual ~ForceConstantCalculator() {};
      void clear(Supercell&, _xinput&, _aflags&, _kflags&, _xflags&, string&, ofstream&);

      virtual bool runVASPCalculations(bool) {return false;};  // ME20191029
      bool runVASPCalculationsBE(_xinput&, uint);
      void setPolarMaterial(bool b) { _isPolarMaterial = b; }  // ME20200218

      bool run();  // ME20191029
      virtual void hibernate(const string&) {};
      void writeHibernateHeader(stringstream&);
      void writeForceConstants(stringstream&);
      void writePolar(stringstream&);

      const vector<vector<xmatrix<double> > >& getForceConstants() const;
      const vector<xmatrix<double> >& getBornEffectiveChargeTensor() const;
      const xmatrix<double>& getDielectricTensor() const;
      bool isPolarMaterial() const;

      // Born charges + dielectric tensor
      bool calculateDielectricTensor(const _xinput&);  // ME20191029
      void readBornEffectiveChargesFromAIMSOUT(void);
      void readBornEffectiveChargesFromOUTCAR(const _xinput&);  // ME20190113
      void symmetrizeBornEffectiveChargeTensors(void);
      void readDielectricTensorFromAIMSOUT(void);
      void readDielectricTensorFromOUTCAR(const _xinput&);  // ME20190113
      virtual void saveState(const string&) {}  // ME20200112
      virtual void readFromStateFile(const string&) {};  // ME20200112
  };
}

// ***************************************************************************
// "dirphoncalc.h"
namespace apl {
  class DirectMethodPC : public ForceConstantCalculator {
    protected:
      //bool GENERATE_PLUS_MINUS;  //JAHNATEK ORIGINAL
      //CO - START
      bool AUTO_GENERATE_PLUS_MINUS;
      bool USER_GENERATE_PLUS_MINUS;
      //CO - END
      bool GENERATE_ONLY_XYZ;
      bool DISTORTION_SYMMETRIZE; //CO20190108
      double DISTORTION_MAGNITUDE;
      bool DISTORTION_INEQUIVONLY; //CO20190108
      string zerostate_dir;  // ME20191030
      // For each inequivalent atom, there is a set of unique distortions
      vector<vector<xvector<double> > > _uniqueDistortions;
      // For each inequivalent atom and unique distortion, there is a field
      // of forces (for each atom of the supercell)
      vector<vector<vector<xvector<double> > > > _uniqueForces;

    protected:
      void estimateUniqueDistortions(const xstructure&,
          vector<vector<xvector<double> > >&);
      void testDistortion(const xvector<double>&, const vector<_sym_op>&,
          vector<xvector<double> >&,
          vector<xvector<double> >&,
          bool integrate_equivalent_distortions=true);  //CO20190114
      bool needMinus(uint atom_index, uint distortion_index, bool inequiv_only=true);  //CO //CO20190218
      bool runVASPCalculations(bool);  // ME20190412

    public:
      DirectMethodPC();
      DirectMethodPC(Supercell&, _xinput&, _aflags&, _kflags&, _xflags&, string&, ofstream&);
      DirectMethodPC(const DirectMethodPC&);
      DirectMethodPC& operator=(const DirectMethodPC&);
      ~DirectMethodPC();
      void clear(Supercell&, _xinput&, _aflags&, _kflags&, _xflags&, string&, ofstream&);

      bool calculateForceFields();  // ME20190412  // ME20191029
      // Easy access to global parameters
      //void setGeneratePlusMinus(bool b) { GENERATE_PLUS_MINUS = b; } //JAHNATEK ORIGINAL
      void setDistortionMagnitude(double f) { DISTORTION_MAGNITUDE = f; }
      void setDistortionINEQUIVONLY(bool b) { DISTORTION_INEQUIVONLY = b; } //CO20190108
      void setCalculateZeroStateForces(bool b) { _calculateZeroStateForces = b; }
      void setGeneratePlusMinus(bool _auto_, bool _user_) {
        AUTO_GENERATE_PLUS_MINUS = _auto_;
        USER_GENERATE_PLUS_MINUS = _user_;
      }  //CO
      void setGenerateOnlyXYZ(bool b) { GENERATE_ONLY_XYZ = b; }
      void setDistortionSYMMETRIZE(bool b) { DISTORTION_SYMMETRIZE = b; } //CO20190108
      bool calculateForceConstants();  // ME20200211

      void hibernate(const string&);
      void writeFORCES();
      void writeDYNMAT();
      void writeXCrysDenForces();
      void saveState(const string&);  // ME200212
      void readFromStateFile(const string&);  // ME200212

    private:
      void copy(const DirectMethodPC&);
      void free();
      vector<vector<bool> > vvgenerate_plus_minus;  // ME20191029
      void completeForceFields();
      void projectToCartesianDirections();
      void buildForceConstantMatrices();

      void writeForceField(stringstream&);
  };
}  // namespace apl

//PINKU - START
// ***************************************************************************
//QHA-APL MACROS START PINKU
#define _AFLOW_QHA_PHONS_DIRECTORY_PREFIX_  string("ARUN.APL_PHONON_")
#define _EOS_AFLOW_DIRECTORY_PREFIX_ string("ARUN.APL_STATIC_")
#define MIN_FREQ_TRESHOLD -0.1//in AMU
#define RAW2Hz 15.6333046177
#define AMU2Kg 1.66053904 
//#define MIN_EIGEN_TRESHOLD -1.0e-2// eigenvalue treshold in AMU

#define _GP_AFLOW_DIRECTORY_PREFIX_ string("ARUN.APL_PHONON_")
#define _EOS_AFLOW_DIRECTORY_PREFIX_ string("ARUN.APL_STATIC_")
#define _PH_AFLOW_DIRECTORY_PREFIX_ string("ARUN.APL_PHONON_")
#define _aflowinpad_ 54
#define Set_QHA_Precision 15
//QHA-APL MACROS END PINKU

// ***************************************************************************
namespace apl {
  class QHA_AFLOWIN_CREATOR : public ForceConstantCalculator
  {
    private:
      ofstream _log;
      string _logfile;
      bool _is_gp_on, _is_gp_A_on, _is_gp_B_on, _is_gp_C_on;
      bool _is_sc_gp_on, _is_sc_gp_A_on,  _is_sc_gp_B_on,  _is_sc_gp_C_on;
      bool _is_eos, _is_eos_A, _is_eos_B, _is_eos_C;
      bool _is_edos_accurate_on;
      double _gp_vol_distortion;
      double _scqha_vol_distortion;

      vector<string> _scqha_dir_names;
      vector<string> _gp_dir_names;
      vector<string> _ph_dir_names;
      vector<string> _eos_dir_names;
      vector<double> _gp_volumes;
      vector<double> _scqha_volumes;
      vector<double> _eos_volumes;
      vector<double> _ph_volumes;

      int _zero_index, _lattice_index;
      string _pstress;

      double _EOS_VOL_START, _EOS_VOL_END, _EOS_VOL_INC;
      int _NEDOS;
      string _EOS_KSCHEME;
      string _EOS_STATIC_KSCHEME;
      int _EOS_STATIC_KPPRA;
      string _PSTRESS;  // ME20200220

    public:
      QHA_AFLOWIN_CREATOR(Supercell& sc, _xinput& xinput,
          _aflags& aflags, _kflags& kflags,
          _xflags& xflags,
          string&,
          ofstream& mf);
      ~QHA_AFLOWIN_CREATOR();
      void clear();
    public:
      //
      void run_qha();
      //
      void close_log();
      //functions to set user define keys 
      void setGP(bool, bool, bool, bool);
      void setSCGP(bool, bool, bool, bool);
      void setEOS(bool b);
      void setGP_VOL_DISTORTION(double d);
      void setSCGP_VOL_DISTORTION(double b);
      void setEOS_distortion_range(double a, double b, double c);
      void setEOS_NEDOS(int s);
      void setEOS_STATIC_KPPRA(int s);
      void setEOS_STATIC_KSCHEME(string s);
      void set_edos_accurate(bool b);
      //interface functions
      vector<string> get_scqha_dir_names();
      vector<string> get_gp_dir_names();
      vector<string> get_ph_dir_names();
      vector<string> get_eos_dir_names();
      vector<double> get_gp_volumes();
      vector<double> get_scqha_volumes();
      vector<double> get_eos_volumes();
      double get_scqha_vol_distortion();
      int get_zero_index();
    private:
      template <typename T> string NumToStr ( T Number );
      //create AFLOWIN for Gruneisen parameter
      //[phonon_option] 0->gp   || 1->sc-gp   || 2-> eos-phonon   || 3->eos-static 
      //[phonon_option] 4->gp_X || 5->sc-gp_X || 6-> eos-phonon_X || 7->eos-static-X
      void create_aflowin_phonon(const double distortion, const int phonon_option);
      //
      void create_aflowin_phonon_X(const double distortion, const int phonon_option);
      //
      void write_aflowin_phonon(const _xvasp& xvasp_input, const int phonon_option);
      void write_phonon_OUTPUT(const _xinput& xinput, const int phonon_option);
      //
      void write_static_AFLOWIN(const _xvasp& xvasp_input);
      void write_static_OUTPUT(const _xinput& xinput);

      void get_pstress();
      //
      void create_aflowin_zero_state(_xinput& xinput);
      void writeZEROOUTPUT(const _xinput& xinput);
      //
      void correcting_scqha_vol_distortion();
      //
      void create_aflowin_scqha_phonon(const _xvasp& xvasp_input);
      void writeSCQHAOUTPUT(const _xinput& xinput);
      //
      void create_aflowin_static_zero();
      void create_aflowin_static_zero_X();
      //
      string get_phonon_runname(const double i, const double distortion);
      string get_phonon_runname(const double i);
      string get_static_runname(const double i);
  };
}
// ***************************************************************************

// ***************************************************************************
// "lrphoncalc.h"
namespace apl {
  class LinearResponsePC : public ForceConstantCalculator {
    private:
      void free();
      void copy(const LinearResponsePC&);
      bool runVASPCalculationsDFPT(_xinput&);  // ME20190113  // ME20200213 - changed name
      bool readForceConstantsFromVasprun(_xinput&);  // ME20200211

    public:
      LinearResponsePC();
      LinearResponsePC(Supercell&, _xinput&, _aflags&, _kflags&, _xflags&, string&, ofstream&);
      LinearResponsePC(const LinearResponsePC&);
      LinearResponsePC& operator=(const LinearResponsePC&);
      ~LinearResponsePC();
      void clear(Supercell&, _xinput&, _aflags&, _kflags&, _xflags&, string&, ofstream&);

      bool runVASPCalculations(bool);  // ME20191029
      bool calculateForceConstants();  // ME20200211

      void hibernate(const string&);
      void saveState(const string&);  // ME200212
      void readFromStateFile(const string&); // ME200212
  };
}  // namespace apl

// ***************************************************************************
// "gsa.h"

//CO generally redirects to DM, the distinction between DM and GSA is obsolete
namespace apl {
  class GeneralizedSupercellApproach : public DirectMethodPC {
    public:
      GeneralizedSupercellApproach(Supercell&, _xinput&,
          _aflags&, _kflags&, _xflags&, //_vflags&, 
          string&, ofstream&);
      ~GeneralizedSupercellApproach();
      void clear();
  };
}  // namespace apl

namespace apl {
  class PhononCalculator {
    protected:
      // USER PARAMETERS
      string _system;  // ME20190614 - for VASP-style output files
      string _directory;  // for loggers
      int _ncpus;

      Supercell* _supercell;

      // harmonic IFCs
      vector<vector<xmatrix<double> > > _forceConstantMatrices;

      // Stuff for polar materials
      bool _isPolarMaterial;
      // For each atom there is a matrix 3x3 of Born effective charge
      vector<xmatrix<double> > _bornEffectiveChargeTensor;
      // Dielectric tensor
      xmatrix<double> _dielectricTensor;
      // Precomputed values used in non-analytical term (Gonze)
      xmatrix<double> _inverseDielectricTensor;
      double _recsqrtDielectricTensorDeterminant;
      // Precomputed Ewald sum at Gamma point
      bool _isGammaEwaldPrecomputed;
      vector<xmatrix<xcomplex<double> > > _gammaEwaldCorr;
      // Anharmonic IFCs
      vector<vector<vector<double> > > anharmonicIFCs;
      vector<vector<vector<int> > > clusters;
      ofstream* messageFile;


    private:
      void copy(const PhononCalculator&);  // ME20191228
      void free();
      xmatrix<xcomplex<double> > getNonanalyticalTermWang(const xvector<double>&);
      xmatrix<xcomplex<double> > getNonanalyticalTermWang(const xvector<double>&,
          vector<xmatrix<xcomplex<double> > >&, bool=true);  // ME20180829
      xmatrix<xcomplex<double> > getNonanalyticalTermGonze(const xvector<double>);
      xmatrix<xcomplex<double> > getEwaldSumDipoleDipoleContribution(const xvector<double>, bool = true);

    public:
      PhononCalculator();
      PhononCalculator(Supercell&, ofstream&);
      PhononCalculator(const PhononCalculator&);
      PhononCalculator& operator=(const PhononCalculator&);
      virtual ~PhononCalculator();
      void clear(Supercell&, ofstream&);

      // Getter functions
      const Supercell& getSupercell() const;
      const xstructure& getInputCellStructure() const;
      const xstructure& getSuperCellStructure() const;
      uint getNumberOfBranches() const;
      string getSystemName() const;  // ME20190614
      string getDirectory() const;
      int getNCPUs() const;
      ofstream& getOutputStream();
      bool isPolarMaterial() const;  // ME20200206
      const vector<vector<xmatrix<double> > >& getHarmonicForceConstants() const;
      const vector<vector<double> >& getAnharmonicForceConstants(int) const;
      const vector<vector<int> >& getClusters(int) const;

      // Set functions
      void setSystem(const string&);
      void setDirectory(const string&);
      void setNCPUs(const _kflags&);

      // IFCs
      void setHarmonicForceConstants(const ForceConstantCalculator&);
      void awake(const string&, bool=true);
      void setAnharmonicForceConstants(const AnharmonicIFCs&);
      void readAnharmonicIFCs(string, bool=true);

      // Dynamical Matrix/Frequencies
      xvector<double> getEigenvalues(const xvector<double>&, const xvector<double>&,
          xmatrix<xcomplex<double> >&, vector<xmatrix<xcomplex<double> > >&, bool=true);  // ME20180827
      xmatrix<xcomplex<double> > getDynamicalMatrix(const xvector<double>&);
      xmatrix<xcomplex<double> > getDynamicalMatrix(const xvector<double>&, const xvector<double>&);  // ME20200206
      xmatrix<xcomplex<double> > getDynamicalMatrix(const xvector<double>&, const xvector<double>&,
          vector<xmatrix<xcomplex<double> > >&, bool=true);  // ME20180827
      xvector<double> getFrequency(const xvector<double>&, const IPCFreqFlags&);  // ME20180827
      xvector<double> getFrequency(const xvector<double>&, const xvector<double>&, const IPCFreqFlags&);  // ME20200206
      xvector<double> getFrequency(const xvector<double>&, const IPCFreqFlags&, xmatrix<xcomplex<double> >&);  // ME20190624
      xvector<double> getFrequency(const xvector<double>&, const xvector<double>&, const IPCFreqFlags&, xmatrix<xcomplex<double> >&);  // ME20200206
      xvector<double> getFrequency(const xvector<double>&, const IPCFreqFlags&, xmatrix<xcomplex<double> >&,
          vector<xmatrix<xcomplex<double> > >&, bool=true);  // ME20180827
      xvector<double> getFrequency(const xvector<double>&, const xvector<double>&, const IPCFreqFlags&, xmatrix<xcomplex<double> >&,
          vector<xmatrix<xcomplex<double> > >&, bool=true);  // ME20200206
      double getFrequencyConversionFactor(IPCFreqFlags, IPCFreqFlags);
  };
}  // namespace apl


// ***************************************************************************

//PINKU - START
// OBSOLETE, moved to .aflow.rc - ME20181024
//#define AFLOW_APL_VASP_USE_LEPSILON
//#undef AFLOW_APL_VASP_USE_LCALCEPS  // HAS some problem [PINKU]
//PINKU - END

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// Supplementary classes for calculation of dispersion curves and density of states

// ***************************************************************************
// "pathbuilder.h"
namespace apl {
  class PathBuilder {
    public:
      enum StoreEnumType { RECIPROCAL_LATTICE,
        CARTESIAN_LATTICE };
      enum ModeEnumType { SINGLE_POINT_MODE,
        COUPLE_POINT_MODE };

    private:
      std::vector<aurostd::xvector<double> > _path;
      std::vector<aurostd::xvector<double> > _points;
      std::vector<std::string> _labels;
      aurostd::xmatrix<double> reciprocalLattice;
      aurostd::xmatrix<double> cartesianLattice;
      int _pointsVectorDimension;
      int _pointsVectorStartingIndex;
      uint _nPointsPerSubPath;
      ModeEnumType _mode;
      StoreEnumType _store;

    private:
      void buildPath();
      void free();
      void copy(const PathBuilder&);
      // [OBSOLETE ME20190219] void tokenize(const string&, vector<string>&, string);

    public:
      PathBuilder();
      PathBuilder(ModeEnumType);
      PathBuilder(const PathBuilder&);
      PathBuilder& operator=(const PathBuilder&);
      ~PathBuilder();
      void clear();
      void addPoint(const std::string& l, int dim, ...);
      void addPoint(const std::string&, const aurostd::xvector<double>&);
      void transform(const aurostd::xmatrix<double>&);
      void pointsAreDefinedFor(const xstructure&, StoreEnumType);
      void transformPointsFor(const xstructure&, StoreEnumType);
      void defineCustomPoints(const string&,const string&,const Supercell&,bool CARESTIAN_COORDS=false);
      void takeAflowElectronicPath(const string&,const Supercell&);//, const xstructure&, const xstructure&);
      void setMode(ModeEnumType);
      void setStore(StoreEnumType);
      const StoreEnumType& getStore() const;  // ME20190614
      void setDensity(int);
      int getDensity();
      uint getPathSize();
      uint getPointSize();
      aurostd::xvector<double> getPoint(uint);
      uint getPointIndexOnPath(uint);
      std::string getPointLabel(uint);
      std::vector<aurostd::xvector<double> > getPath();
      std::vector<aurostd::xvector<double> > getPath(ModeEnumType, const string&);
      double getPathLength();
      double getPathLength(uint);
      xKPOINTS createKPOINTS(const Supercell&);  // ME20190614
  };
}  // namespace apl

// ***************************************************************************
// "pdisc.h"
namespace apl {
  class PhononDispersionCalculator {
    private:
      PhononCalculator* _pc;
      PathBuilder _pb;
      void copy(const PhononDispersionCalculator&);
      void free();
      std::vector<xvector<double> > _qpoints;
      std::vector<xvector<double> > _freqs;
      IPCFreqFlags _frequencyFormat;
      double _temperature;  // ME20190614
      //[OBSOLETE PN180705]vector<double> path;       //[PINKU]
      //[OBSOLETE PN180705]vector<int> path_segment;  //[PINKU]
      void calculateInOneThread(int, int);

    public:
      PhononDispersionCalculator();
      PhononDispersionCalculator(PhononCalculator&);
      PhononDispersionCalculator(const PhononDispersionCalculator&);
      PhononDispersionCalculator& operator=(const PhononDispersionCalculator&);
      ~PhononDispersionCalculator();
      void clear(PhononCalculator&);
      void clear(PhononDispersionCalculator&);
      void initPathCoords(const string&,const string&,int,bool=false);  //CO20180406
      void initPathLattice(const string&, int);
      void setPath(const string&);
      void calc(const IPCFreqFlags);
      void writePDIS(const string&);
      bool isExactQPoint(const xvector<double>&, const xmatrix<double>&);
      std::vector<xvector<double> > get_qpoints() { return _qpoints; }  //[PINKU]
      // ME20190614 - START
      xEIGENVAL createEIGENVAL();
      void writePHEIGENVAL(const string&);
      void writePHKPOINTS(const string&);
      string _system;
      // ME2019614 - STOP
      //[OBSOLETE PN180705]std::vector<double> get_path() { return path; }                   //[PINKU]
      //[OBSOLETE PN180705]std::vector<int> get_path_segment() { return path_segment; }      //[PINKU]
  };
}  // namespace apl

// ***************************************************************************
namespace apl { //PN180705
  class PhononHSQpoints {
    private:
      Logger& _logger;
      vector< xvector<double> > _qpoints;
      vector<double> _path;
      vector<int> _path_segment;
      vector<xvector<double> >_hs_kpoints;
    private:
      template<typename T>
        std::vector<T> split(const std::string& line);
    public:
      PhononHSQpoints(Logger&);
      ~PhononHSQpoints();
      void clear();
      void read_qpointfile(const string&);
      //interface functions
      vector<xvector<double> > get_qpoints();
      vector<xvector<double> > get_hs_kpoints();
      vector<double> get_path();
      vector<int> get_path_segment();
  };
} // namespace apl
// ***************************************************************************
// ME20190417 - BEGIN
// ***************************************************************************
// OBSOLETE - all q-point grids are now described using the QMesh class
// "irpg.h"
//
//namespace apl {
//class IReciprocalPointGrid {
// public:
//  virtual ~IReciprocalPointGrid() {}
//  virtual int getN(int) = 0;
//  virtual aurostd::xmatrix<double> getReciprocalLattice() = 0;
//  virtual std::vector<aurostd::xvector<double> > getPoints() = 0;
//  virtual std::vector<double> getWeights() = 0;
//  virtual aurostd::xtensor<int> getAllToIrrPointMap() = 0; //ME20180705
//};
//}  // end namespace apl
//
// ***************************************************************************
// "mpmesh.h"
//namespace apl {
//class MonkhorstPackMesh : virtual public IReciprocalPointGrid {
// private:
//  Logger& _logger;
//  aurostd::xvector<int> _n;
//  std::vector<aurostd::xvector<double> > _kpoints;
//  std::vector<double> _weights;
//  aurostd::xvector<double> _shift;
//  aurostd::xmatrix<double> _rlattice;
//  aurostd::xmatrix<double> _klattice;
//  aurostd::xtensor<int> _allToIrrPointMap; //ME20180705
//  bool _isIrreducible;
//
// private:
//  void generateAllGridPoints();
//  void makeIrreducible(const std::vector<_sym_op>&);
//
// public:
//  MonkhorstPackMesh(int, int, int, const xstructure&, Logger&);
//  ~MonkhorstPackMesh();
//  void clear();
//  void setDensity(int);
//  void setDensity(int, int, int);
//  // Interface
//  int getN(int);
//  aurostd::xmatrix<double> getReciprocalLattice();
//  std::vector<aurostd::xvector<double> > getPoints();
//  std::vector<double> getWeights();
//  aurostd::xtensor<int> getAllToIrrPointMap(); //ME20180705
//};
//}  // namespace apl

namespace apl {

  struct _qpoint {
    xvector<double> cpos;  // Cartesian position of the q-point
    xvector<double> fpos;  // Fractional coordinates of the q-point
    int irredQpt;  // The irreducible q-point this q-point belongs to
    int ibzqpt;  // The index of the irreducible q-point in the _ibzqpt vector
    int symop;  // Symmetry operation to transform the irreducible q-point into this q-point
  };

  struct _kcell {
    xmatrix<double> lattice;  // The reciprocal lattice vectors
    xmatrix<double> rlattice;  // The real space lattice
    xmatrix<double> c2f;  // Conversion matrix from Cartesian to fractional
    xmatrix<double> f2c;  // Conversion matrix from fractional to Cartesian
    bool skewed;  // Is the lattice skewed?
    vector<_sym_op> pgroup;  // The point group operations of the reciprocal cell
  };

  class QMesh {
    public:
      QMesh();
      QMesh(ofstream&);
      QMesh(const xvector<int>&, const xstructure&, ofstream&, bool=true, bool=true, string="./");
      QMesh(const vector<int>&, const xstructure&, ofstream&, bool=true, bool=true, string="./");
      QMesh(const QMesh&);
      QMesh& operator=(const QMesh&);
      ~QMesh();
      void clear(ofstream&);
      void initialize(const vector<int>&, const xstructure& xs, bool=true, bool=true);
      void initialize(const xvector<int>&, const xstructure& xs, bool=true, bool=true);

      void setDirectory(const string& dir);
      void setModule(const string&);
      const string& getDirectory() const;
      const string& getModule() const;
      ofstream& getOutputStream();

      void makeIrreducible();
      void calculateLittleGroups();  // ME20200109
      void writeQpoints(string, bool=true);
      void writeIrredQpoints(string, bool=true);

      int getnIQPs() const;
      int getnQPs() const;
      int getGrid(int) const;
      const xvector<int>& getGrid() const;
      const _qpoint& getIrredQPoint(int) const;
      const _qpoint& getIrredQPoint(int, int, int) const;
      vector<xvector<double> > getIrredQPointsCPOS() const;
      vector<xvector<double> > getIrredQPointsFPOS() const;
      int getIrredQPointIndex(int) const;
      int getIrredQPointIndex(int, int, int) const;
      const _qpoint& getQPoint(int) const;
      const _qpoint& getQPoint(int, int, int) const;
      const _qpoint& getQPoint(const xvector<double>&) const;  // ME20190813
      int getQPointIndex(xvector<double>) const;  // ME20190813
      int getQPointIndex(int, int, int) const;
      vector<xvector<double> > getQPointsCPOS() const;
      vector<xvector<double> > getQPointsFPOS() const;
      int getIbzqpt(int) const;
      int getIbzqpt(int, int, int) const;
      const vector<int>& getIbzqpts() const;
      const vector<_qpoint>& getPoints() const;
      const _kcell& getReciprocalCell() const;
      bool isShifted() const;  // ME20190813
      const xvector<double>& getShift() const;
      const vector<int>& getWeights() const;
      bool isReduced() const;
      bool isGammaCentered() const;
      bool littleGroupsCalculated() const;  // ME20200109
      const vector<int>& getLittleGroup(int) const;  // ME20200109

    private:
      void free();
      void copy(const QMesh&);

      ofstream* messageFile;
      string _directory;

      vector<int> _ibzqpts;  // The indices of the irreducible q-points
      bool _isGammaCentered;  // Indicates whether the includes the Gamma point
      vector<vector<int> > _littleGroups;  // The little groups of the irreducible q-points
      bool _littleGroupsCalculated;  // Indicates whether the little groups have been calculated
      int _nIQPs;  // The number of irreducible q-points
      int _nQPs;  // The number of q-points
      xvector<int> _qptGrid;  // The dimensions of the q-point mesh
      vector<vector<vector<int> > > _qptMap;  // Maps a q-point triplet to a q-point index
      vector<_qpoint> _qpoints;  // The q-points of the mesh
      _kcell _recCell;  // The reciprocal cell
      bool _reduced;  // Indicates whether the q-point mesh has been reduced
      bool _shifted;  // Indicates whether the q-point mesh has been shifted
      xvector<double> _shift;  // The shift vector of the mesh
      vector<int> _weights;  // The weights of each irreducible q-point

      void setGrid(const xvector<int>&);
      void setupReciprocalCell(xstructure, bool);
      void generateGridPoints(bool);
      void shiftMesh(const xvector<double>&);
      void moveToBZ(xvector<double>&) const;
  };
}  // namespace apl

namespace apl {
  class LTMethod {
    public:
      LTMethod();
      LTMethod(QMesh&);
      LTMethod(const LTMethod&);
      LTMethod& operator=(const LTMethod&);
      ~LTMethod();
      void clear(QMesh&);

      void makeIrreducible();  // ME20190625

      const vector<vector<int> >& getTetrahedra() const;
      const vector<int>& getTetrahedron(int) const;
      const vector<int>& getIrredTetrahedron(int) const;
      int getCorner(int, int) const;
      // const int& getCornerIrred(const int&, const int&) const; OBSOLETE ME201906
      vector<vector<int> > getIrreducibleTetrahedra() const;  // ME20190625
      vector<vector<int> > getIrreducibleTetrahedraIbzqpt() const;  // ME20190625

      int getnTetrahedra() const;
      int getnIrredTetrahedra() const;
      double getVolumePerTetrahedron() const;
      const vector<int>& getWeights() const;
      int getWeight(int) const;
      bool isReduced() const;
    private:
      void free();
      void copy(const LTMethod&);

      QMesh* _qm;

      vector<vector<int> > _tetrahedra;  // The corners of the tetrahedra - ME20190625
      vector<int> _irredTetrahedra;  // List of irreducible tetrahedra - ME20190625
      bool _reduced;  // ME20190625
      int _nTetra;  // The number of tetrahedra
      int _nIrredTetra;  // The number of irreducible tetrahedra - ME20190625
      double _volumePerTetrahedron;  // The relative volume of each tetrahedron
      vector<int> _weights;  // The weights of each irreducible tetrahedron

      void generateTetrahedra();
      vector<vector<xvector<int> > > initializeTetrahedra();
      void findMostCompact(vector<vector<xvector<int> > >&);
      void generateAllTetrahedra(const vector<vector<xvector<int> > >&);  // ME20190625
  };
}  // namespace apl

// ***************************************************************************
// ME20190417 - END
// ***************************************************************************
// ***************************************************************************
// "doscalc.h"
namespace apl {
#define MIN_FREQ_TRESHOLD -0.1  //in AMU
#define RAW2Hz 15.6333046177
#define AMU2Kg 1.66053904
#define MIN_EIGEN_TRESHOLD -1e-2  // eigenvalue treshold in AMU
  class DOSCalculator  // ME20190424
  { //CO200106 - patching for auto-indenting
    protected:
      PhononCalculator* _pc;
      //  IReciprocalPointGrid& _rg;  OBSOLETE ME20190417
      QMesh* _rg;
      string _bzmethod;  // ME20190423
      std::vector<aurostd::xvector<double> > _qpoints;
      //std::vector<int> _qweights;  OBSOLETE ME20190423
      std::vector<aurostd::xvector<double> > _freqs;
      double _minFreq;
      double _maxFreq;
      double _stepDOS;
      double _halfStepDOS;
      std::vector<double> _bins;
      std::vector<double> _dos;
      std::vector<double> _idos;  // ME20190614
      std::vector<xmatrix<xcomplex<double> > > _eigen;  // ME20190624 - eigenvectors for projected DOS
      std::vector<vector<vector<double> > > _projectedDOS; // ME20190614 - projectedDOS.at(atom).at(direction).at(frequency)
      std::vector<xvector<double> > _projections;  // ME20190626 - the projection directions for the DOS in Cartesian coordinates
      double _temperature;  // ME20190614
      //CO - START
      //private:
      void copy(const DOSCalculator&);
      void calculateInOneThread(int, int);
      //CO - END
      void calculateFrequencies();
      void smearWithGaussian(vector<double>&, vector<double>&, double, double);  // ME20190614

    public:
      DOSCalculator();
      DOSCalculator(PhononCalculator&, QMesh&, string, const vector<xvector<double> >&);  // ME20190423 // ME20190624 - added projections
      DOSCalculator(const DOSCalculator&);
      DOSCalculator& operator=(const DOSCalculator&);
      ~DOSCalculator();
      void clear(PhononCalculator&, QMesh&);
      // ME20190423 - START
      //virtual void rawCalc(int) {} OBSOLETE ME20190419
      void calcDosRS();
      void calcDosLT();
      // ME20190423 - END
      void calc(int);
      void calc(int, double);
      void calc(int, double, double, double);  // ME20200203
      void writePDOS(const string&);
      void writePDOS(string, string);  //[PINKU]
      xDOSCAR createDOSCAR();  // ME20190614
      void writePHDOSCAR(const string&);  // ME20190614
      // Interface IDOSCalculator
      const std::vector<double>& getBins() const;  // ME20200108 - added const
      const std::vector<double>& getDOS() const;   // ME20200108 - added const
      const std::vector<double>& getIDOS() const;  // ME20200210
      bool hasNegativeFrequencies() const;  // ME20200108 - added const
      string _system;  // ME20190614
    private:
      void free();
  };
}  // namespace apl

// ***************************************************************************
// OBSOLETE ME20190423 - RootSamplingMethod and LinearTetrahedronMethod have
// been integrated into DOSCalculator
// "rsmdos.h"
//namespace apl {
//class RootSamplingMethod : public DOSCalculator {
// public:
//  RootSamplingMethod(IPhononCalculator&, IReciprocalPointGrid&, Logger&);
//  ~RootSamplingMethod();
//  void rawCalc(int);
//};
//}  // namespace apl
//
//// ***************************************************************************
//// "ltetdos.h"
//namespace apl {
//class LinearTetrahedronMethod : public DOSCalculator {
// private:
//  double _weightVolumeOfEachTetrahedron;
//  std::vector<std::vector<int> > _irrTetrahedraList;
//  std::vector<int> _irrTetrahedraWeightList;
//
// private:
//  void generateTetrahedras();
//
// public:
//  LinearTetrahedronMethod(IPhononCalculator&, IReciprocalPointGrid&, Logger&);
//  ~LinearTetrahedronMethod();
//  void clear();
//  void rawCalc(int);
//};
//}  // end namespace apl

// ***************************************************************************
// thermalpc.h
namespace apl {
  enum ThermalPropertiesUnits { eV,
    meV,
    ueV,
    eVK,
    meVK,
    ueVK,
    kB };

  class ThermalPropertiesCalculator {
    private:
      ofstream* messageFile;
      string _directory;
      std::vector<double> _freqs_0K;
      std::vector<double> _dos_0K;
      string system;

      void free();
      void copy(const ThermalPropertiesCalculator&);
 
      double getStepDOS(const vector<double>&);
      double getScalingFactor(const ThermalPropertiesUnits&);

    public:
      ThermalPropertiesCalculator();
      ThermalPropertiesCalculator(ofstream&);
      ThermalPropertiesCalculator(const DOSCalculator&, ofstream&, string="./");
      ThermalPropertiesCalculator(const xDOSCAR&, ofstream&, string="./");
      ThermalPropertiesCalculator(const ThermalPropertiesCalculator&);
      ThermalPropertiesCalculator& operator=(const ThermalPropertiesCalculator&);
      ~ThermalPropertiesCalculator();
      void clear(ofstream&);

      void setDirectory(const string&);
      string getDirectory() const;

      vector<double> temperatures;
      vector<double> Cv;
      vector<double> Fvib;
      vector<double> Svib;
      vector<double> U;
      double U0;

      void initialize(const vector<double>&, const vector<double>&, string="");
      void calculateThermalProperties(double, double, double);
      void addPoint(double, const xDOSCAR&);
      void addPoint(double, const vector<double>&, const vector<double>&);

      double getZeroPointEnergy();
      double getInternalEnergy(double, ThermalPropertiesUnits=apl::meV);
      double getInternalEnergy(double, const vector<double>&, const vector<double>&, ThermalPropertiesUnits=apl::meV);
      double getVibrationalFreeEnergy(double, ThermalPropertiesUnits=apl::meV);
      double getVibrationalFreeEnergy(double, const vector<double>&, const vector<double>&, ThermalPropertiesUnits=apl::meV);
      double getVibrationalEntropy(double, ThermalPropertiesUnits=apl::meV);
      double getVibrationalEntropy(double, const vector<double>&, const vector<double>&, ThermalPropertiesUnits=apl::kB);
      double getVibrationalEntropy(double, double, double, ThermalPropertiesUnits=apl::kB);
      double getIsochoricSpecificHeat(double, ThermalPropertiesUnits=apl::kB);
      double getIsochoricSpecificHeat(double, const vector<double>&, const vector<double>&, ThermalPropertiesUnits=apl::kB);
   
      void writePropertiesToFile(string);
  };
}  // namespace apl

// ***************************************************************************
// BEGIN ME: Lattice Thermal Conductivity (AAPL)
// ***************************************************************************

namespace apl {
  //[ME20190520 - MOVED UP]struct _qpoint {
  //[ME20190520 - MOVED UP]  xvector<double> cpos;  // Cartesian position of the q-point
  //[ME20190520 - MOVED UP]  xvector<double> fpos;  // Fractional coordinates of the q-point
  //[ME20190520 - MOVED UP]  xvector<int> indices;  // Indices of the q-point grid for this q-point
  //[ME20190520 - MOVED UP]  int symop;  // Symmetry operation to transform into an irreducible q-point
  //[ME20190520 - MOVED UP]};

  //[ME20190520 - MOVED UP]struct _kcell {
  //[ME20190520 - MOVED UP]  xmatrix<double> lattice;  // The reciprocal lattice vectors
  //[ME20190520 - MOVED UP]  xmatrix<double> c2f;  // Conversion matrix from Cartesian to fractional
  //[ME20190520 - MOVED UP]  xmatrix<double> f2c;  // Conversion matrix from fractional to Cartesian
  //[ME20190520 - MOVED UP]  bool skewed;  // Is the lattice skewed?
  //[ME20190520 - MOVED UP]};

  class TCONDCalculator {
    // See aflow_aapl_tcond.cpp for detailed descriptions of the functions
    public:
      TCONDCalculator();
      TCONDCalculator(PhononCalculator&, QMesh&, ofstream&, _aflags&);
      TCONDCalculator(const TCONDCalculator&);
      TCONDCalculator& operator=(const TCONDCalculator&);
      ~TCONDCalculator();
      void clear(PhononCalculator&, QMesh&, ofstream&, _aflags&);
      void initialize();

      aurostd::xoption calc_options; // Options for the the thermal conductivity calculation
      vector<xmatrix<xcomplex<double> > > eigenvectors;  // The eigenvectors at each q-point
      vector<vector<double> > freq;  // The frequencies at each q-point
      vector<vector<xvector<double> > > gvel;  // The group velocities
      int nBranches;  // The number of branches in the phonon spectrum
      int nIQPs;  // The total number of irreducible q-points in the grid
      int nQPs;  // The total number of q-points in the grid
      vector<vector<vector<int> > > processes;  // The sign, q-point and branch indices of the scattering processes
      vector<vector<double> > intr_trans_probs;  // The intrinsic transition probabilities
      vector<vector<vector<int> > > processes_iso;  // The q-point and branch indices of the isotope scattering processes
      vector<vector<double> > intr_trans_probs_iso;  // The intrinsic transition probabilities for isotope processes
      vector<vector<double> > rates_boundary;
      vector<vector<double> > rates_isotope;
      vector<double> temperatures;  // The mperatures for the thermal conductivity calculations
      vector<xmatrix<double> > thermal_conductivity;  // The thermal conductivity values

      void calculateThermalConductivity();

    private:
      PhononCalculator* _pc;  // Reference to the phonon calculator
      QMesh* _qm;  // Reference to the q-point mesh
      Logger _logger;  // The APL logger
      ofstream* messageFile;
      _aflags* aflags;

      void free();
      void copy(const TCONDCalculator&);

      vector<vector<double> > calculateModeGrueneisen(const vector<vector<vector<xcomplex<double> > > >& phases);
      double calculateAverageGrueneisen(double T, const vector<vector<double> >&);

      void calculateFrequenciesGroupVelocities();
      void calculateFreqGvel(int, int);
      void getWeightsLT(const LTMethod&, double, const vector<double>&, vector<double>&);
      void calculateTransitionProbabilities();
      vector<vector<vector<xcomplex<double> > > > calculatePhases(bool=false);
      void calculateTransitionProbabilitiesPhonon(int, int, const LTMethod&,
          vector<vector<vector<vector<double> > > >&,
          const vector<vector<vector<xcomplex<double> > > >&);
      void calculateTransitionProbabilitiesIsotope(int, int, const LTMethod&);
      vector<vector<double> > calculateTransitionProbabilitiesBoundary();
      void getProcess(const vector<int>&, vector<int>&, vector<int>&, int&);
      xmatrix<double> calculateThermalConductivityTensor(double,
          vector<vector<vector<double> > >&,
          vector<vector<vector<double> > >&);
      vector<vector<double> > getOccupationNumbers(double);
      vector<vector<double> > calculateAnharmonicRates(const vector<vector<double> >&);
      vector<vector<double> > calculateTotalRates(const vector<vector<double> >&, vector<vector<vector<double> > >&);
      double getOccupationTerm(const vector<vector<double> >&, int, const vector<int>&, const vector<int>&);
      void calcAnharmRates(int, int, const vector<vector<double> >&, vector<vector<double> >&);
      vector<vector<xvector<double> > > getMeanFreeDispRTA(const vector<vector<double> >&);
      xmatrix<double> calcTCOND(double, const vector<vector<double> >&,
          const vector<vector<xvector<double> > >&);
      void getMeanFreeDispFull(const vector<vector<double> >&,
          const vector<vector<double> >&, vector<vector<xvector<double> > >&);
      void calculateDelta(int, int, const vector<vector<double> >&,
          const vector<vector<xvector<double> > >&, vector<vector<xvector<double> > >&);
      void correctMFD(const vector<vector<double> >&, const vector<vector<xvector<double> > >&, vector<vector<xvector<double> > >&);

      void writeTempIndepOutput(const string&, string, const string&, const vector<vector<double> >&);
      void writeTempDepOutput(const string&, string, const string&, const vector<double>&, const vector<vector<vector<double> > >&);
      void writeDataBlock(stringstream&, const vector<vector<double> >&);
      void writeGroupVelocities(const string&);
      void writePhaseSpace(const string&, const vector<vector<vector<vector<double> > > >&);
      void writeGrueneisen(const string&, const vector<double>&, const vector<vector<double> >&);
      void writeThermalConductivity(const string&);
  };

}  // namespace apl

// ***************************************************************************
// END ME: Lattice Thermal Conductivity (AAPL)
// ***************************************************************************

// OBSOLETE - ME20190428 - all q-point grids are now described using the QMesh class
//namespace apl
//{
//  class UniformMesh
//  { //PN180705
//  private:
//    Logger& _logger;
//    vector<aurostd::xvector<double> > _kpoints;
//    vector<double> _weights;
//    vector<int> _sc_size;
//    xmatrix<double> _rlattice;
//    aurostd::xmatrix<double> _klattice;
//    int _k_index;//k index at gamma ponit
//    public:
//     UniformMesh(Logger&);
//    ~UniformMesh();
//    void clear();
//
//    public:
//        void create_uniform_mesh(int na, int nb, int nc, const xstructure& xs);
//        vector<int> get_sc_size();
//        vector<double> get_weights();
//        vector<aurostd::xvector<double> > get_kpoints();
//        xmatrix<double> get_rlattice();
//        xmatrix<double> get_klattice();
//        int get_k_index();
//  };
//}
// ***************************************************************************
//Functions in this class are used to calculate group velocities and related properties//
namespace apl
{
  class GroupVelocity
  { //PN180705
    private:
      PhononCalculator& _pc;
      //UniformMesh& _umesh;  OBSOLETE ME20190428
      QMesh& _umesh;  // ME20190428
      Logger& _logger;

      std::vector< aurostd::xvector<double> > _freq_kp;
      std::vector< aurostd::xvector<double> > _freq_km;
      std::vector< aurostd::xvector<double> > _freq;
      std::vector< aurostd::xmatrix<xcomplex<double> > > _eigenvectors;
      vector<bool> _freq_test;

      std::vector< aurostd::xmatrix<double> > _gv;//direction dependent group velocities
      std::vector< aurostd::xvector<double> > _phvel;//group velocities average over directions

      std::vector< aurostd::xvector<double> > _kpoints_kp;
      std::vector< aurostd::xvector<double> > _kpoints_km;
      vector<aurostd::xvector<double> > _kpoints;

      //std::vector< double > _weights;  OBSOLETE ME20190428 - not used
      uint  _nBranches;
      double _sound_speed;
      double _kshift;
      void populate_variables();
      void solve_eigenvalues_at_k(int startIndex, int endIndex, int cpuid, int ktype);
      bool eigen_solver(int ktype);
      bool eigen_solver();
      void sound_speed();
      int  indexofSmallestElement(const vector<double> &array);
      void clear_auxiliary_variables();

    public:
      GroupVelocity(PhononCalculator&, QMesh&, Logger&);  // ME20190428
      ~GroupVelocity();
      void clear();

    public:
      bool check_negative_frequencies();
      bool compute_group_velocities();
      void write();
      std::vector< aurostd::xvector<double> > get_freq(){return _freq;}
      std::vector< aurostd::xmatrix<double> > get_gv(){return _gv;}
      std::vector< aurostd::xvector<double> > get_phvel(){return _phvel;}
  };
}
// ***************************************************************************
namespace apl {

  class AtomicDisplacements {
    protected:
      PhononCalculator* _pc;

    private:
      void free();
      void copy(const AtomicDisplacements&);

      vector<vector<vector<xvector<xcomplex<double> > > > > _eigenvectors;
      vector<vector<double> > _frequencies;
      vector<vector<xmatrix<xcomplex<double> > > > _displacement_matrices;
      vector<vector<vector<xvector<xcomplex<double> > > > > _displacement_modes;
      vector<_qpoint> _qpoints;
      vector<double> _temperatures;

      void calculateEigenvectors();
      void calculateEigenvectorsInThread(int, int);
      void calculateMeanSquareDisplacementMatrices();
      void calculateNormalModeDisplacements();
      double getOccupationNumber(double, double);

    public:
      AtomicDisplacements();
      AtomicDisplacements(PhononCalculator&);
      AtomicDisplacements(const AtomicDisplacements&);
      AtomicDisplacements& operator=(const AtomicDisplacements&);
      ~AtomicDisplacements();
      void clear(PhononCalculator&);

      void calculateMeanSquareDisplacements(const QMesh&, double, double, double);
      void calculateNormalModeDisplacements(const vector<xvector<double> >& qpts, bool=true);

      const vector<double>& getTemperatures() const;
      const vector<vector<xmatrix<xcomplex<double> > > >& getDisplacementMatrices() const;
      vector<vector<xvector<double> > > getDisplacementVectors() const;
      const vector<vector<vector<xvector<xcomplex<double> > > > >& getModeDisplacements() const;

      vector<vector<vector<double> > > createDisplacementsXcrysden(const Supercell&, double, double, int, int, int);
      void getOrientedDisplacementsVsim(xstructure&, vector<vector<vector<xvector<xcomplex<double> > > > >&);

      void writeMeanSquareDisplacementsToFile(string);
      void writeSceneFileXcrysden(string, const vector<vector<vector<double> > >&, int);
      void writeSceneFileVsim(string, const xstructure&, const vector<vector<vector<xvector<xcomplex<double> > > > >&);
  };

  void createAtomicDisplacementSceneFile(const aurostd::xoption& vpflow);
}
// ***************************************************************************
// ***************************************************************************
namespace apl
{
  class QHA:public MVops
  { //PN180705
    protected:
      PhononCalculator& _pc;
      QHA_AFLOWIN_CREATOR&  _runeos;
      Logger& _logger;
    public:
      QHA(PhononCalculator&, QHA_AFLOWIN_CREATOR&, Logger&);
      ~QHA();
      void clear();

    private:
      double _cutoff_freq;
      bool _is_negative_freq;
      string _tmp_dir;
      //dynamical matrix at +ve volume
      vector<xmatrix<xcomplex<double> > > _DMp;
      //dynamical matrix at +ve volume
      vector<xmatrix<xcomplex<double> > > _DMm;
      //dynamical matrix at equilibrium volume
      vector <xmatrix<xcomplex<double> > > _DM0;
      //phonon branches 
      uint  _nBranches;
      //kpoints 
      std::vector< aurostd::xvector<double> > _kpoints;
      //kpoint weights
      std::vector<double> _weights;
      //check calculations at very highsymmetry kpoint
      vector<bool> _gp_path_test;
      //check calculations at kpoint in mesh
      vector<bool> _gp_mesh_test;
      //qha gruneisen in mesh
      vector<xvector<double> > _qha_gp_mesh;
      //qha gruneisen along path
      vector<xvector<double> > _qha_gp_path;
      //frequecies along high symmetry path in THz
      vector<xvector< double> >  _freqs_path;
      //frequecies along high symmetry path in THz for negative distortion
      vector<xvector< double> >  _freqs_pathM;
      //frequecies along high symmetry path in THz for positive distortion
      vector<xvector< double> >  _freqs_pathP;
      //fequecies in qmesh in THz
      vector<xvector< double> >  _freqs_mesh;

      //qhasiharmonic phonon directories
      vector<string> _qha_gpdir;
      //qhasiharmonic phonon volumes of all configurations
      vector<double> _qha_gpvol;
      //difference in volume
      double _delta_V;
      //equilibrium volume
      double _V0;
    public:
      //all QHA functions
      //
      //
      void get_tmp_dir_name(const string);
      //Calculate Geuneisen along high symmetry path
      bool calculation_gruneisen(const std::vector< aurostd::xvector<double> > &kpoints);
      //Calculate Gruneisen in mesh
      //bool calculation_gruneisen(apl::UniformMesh* umesh);  OBSOLETE ME20190428
      bool calculation_gruneisen(apl::QMesh* umesh);
      //get negative frequency index 
      bool get_is_negative_freq(){return _is_negative_freq;}
      //set cutoff frequency
      void set_cutoff_freq(double a){_cutoff_freq=a;}
      //Write Gruneisen parameter along path
      void write_gruneisen_parameter_path(const vector <double> &path, const vector <int> &path_seg);
      //Write Gruneisen parameter in mesh
      void write_gruneisen_parameter_mesh();
      //Write average Grineisen
      void Writeaverage_gp(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP);
      //Calculate average gruneisen
      double average_gruneisen_parameter(double temperature_in_kelvins);
      //
      bool set_imported_variables();
      //
      void print_freqs();

    private:
      //calculate gruneisen parameters along high symmetry path
      bool cal_gp_along_path();
      //
      bool cal_gp_in_mesh();
      //get Dynamical matrices along path
      bool get_dynamicalmatrices_along_path();
      //get Dynamical matrices in mesh
      bool get_dynamicalmatrices_in_mesh();
      //read dynamical matrices from different distortions
      bool read_matrix( vector<xmatrix<xcomplex<double> > >&A, const string file);
      //
      template<typename T> std::vector<T> split(const std::string& line);
      //calculate Gruneisen in path usong threads
      void calculate_gp_in_path(int startIndex, int endIndex, int cpuid);
      //calculate Gruneisen in mesh usong threads
      void calculate_gp_in_mesh(int startIndex, int endIndex, int cpuid);
      //read PDIS file
      bool read_PDIS(vector<string> &hash_lines);
      //
      void gruneisen_parameter_300K();
      //
      bool exists_test0 (const std::string& name);
  };
}
// ***************************************************************************
namespace apl
{
  class SCQHA_QHA3P:public MVops
  { //PN180705
    protected:
      PhononCalculator& _pc;
      QHA_AFLOWIN_CREATOR&  _runeos;
      Logger& _logger;
    public:
      SCQHA_QHA3P(PhononCalculator&, QHA_AFLOWIN_CREATOR&, Logger&);
      ~SCQHA_QHA3P();
      void clear();
    private:
      double _cutoff_freq;
      bool _is_negative_freq;
      bool _is_vol_err;
      string _tmp_dir;
      //dynamical matrix at +ve volume
      vector<xmatrix<xcomplex<double> > > _DMp;
      //dynamical matrix at +ve volume
      vector<xmatrix<xcomplex<double> > > _DMm;
      //dynamical matrix at equilibrium volume
      vector <xmatrix<xcomplex<double> > > _DM0;
      //phonon branches 
      uint  _nBranches;
      //kpoints 
      std::vector< aurostd::xvector<double> > _kpoints;
      //kpoint weights
      std::vector< double > _weights;
      //check calculations at very highsymmetry kpoint
      vector<bool> _gp_path_test;
      //check calculations at kpoint in mesh
      vector<bool> _gp_mesh_test;
      //qha gruneisen in mesh
      vector<xvector<double> > _qha_gp_mesh;
      //qha gruneisen along path
      vector<xvector<double> > _qha_gp_path;
      //frequecies along high symmetry path in THz
      vector<xvector< double> >  _freqs_path;
      vector<xvector< double> >  _freqs_pathM;
      vector<xvector< double> >  _freqs_pathP;
      //fequecies in qmesh in THz
      vector<xvector< double> >  _freqs_mesh;
      //fequecies in qmesh in THz
      vector<xvector< double> >  _freqs_meshM;
      //fequecies in qmesh in THz
      vector<xvector< double> >  _freqs_meshP;

      //qhasiharmonic phonon directories
      vector<string> _qha_gpdir;
      //qhasiharmonic phonon volumes of all configurations
      vector<double> _qha_gpvol;
      //difference in volume
      double _delta_V;
      //equilibrium volume
      double _V0;
    public:
      //all QHA functions
      //void setting_qha_type(bool qha, bool scqha, bool eos);
      //
      void get_tmp_dir_name(const string);
      //Calculate Gruneisen in given k-points 
      bool calculation_gruneisen(const std::vector< aurostd::xvector<double> > &kpoints);
      //Calculate Gruneisen in uniform k-mesh
      //bool calculation_gruneisen(apl::UniformMesh* umesh);  OBSOLETE ME20190428
      bool calculation_gruneisen(apl::QMesh* umesh);  // ME20190428
      //
      void set_cutoff_freq(double a){_cutoff_freq=a;}
      //Write Gruneisen
      void write_gruneisen_parameter_path(const vector <double> &path, const vector <int> &path_seg);
      //Write Gruneisen parameter in uniform mesh
      void write_gruneisen_parameter_mesh();
      //Write average Gruneisen parameter 
      void Writeaverage_gp(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP);
      //Calculate average Gruneisen parametars
      double average_gruneisen_parameter(double temperature_in_kelvins);
      //
      bool set_imported_variables();
      //interface functions 
      //get frequencies at 0-distortion
      vector<xvector< double> >  get_freqs_mesh();
      //get frequencies at negative-distortion
      vector<xvector< double> >  get_freqs_meshM();
      //get frequencies at positive-distortion
      vector<xvector< double> >  get_freqs_meshP();
      //get phonon volume
      vector<double> get_qha_gpvol();
      //get k-weights
      vector<double> get_weights();
      bool get_is_negative_freq();
      bool get_is_vol_err();
      void print_freqs();
      void print_path_freqs();
    private:
      //calculate gruneisen parameters along high symmetry path
      bool cal_gp_along_path();
      //Calculate Gruneisen parameter in mesh
      bool cal_gp_in_mesh();
      //get dynamical matrices along path
      bool get_dynamicalmatrices_along_path();
      //get dynamical matrices in mesh
      bool get_dynamicalmatrices_in_mesh();
      //read dynamical matrices from different distortions
      bool read_matrix( vector<xmatrix<xcomplex<double> > >&A, const string file);
      //
      template<typename T> std::vector<T> split(const std::string& line);
      //Calcualte Gruneisen parameter along high symmetry path in threads
      void calculate_gp_in_path(int startIndex, int endIndex, int cpuid);
      //Calcualte Gruneisen parameter in threads
      void calculate_gp_in_mesh(int startIndex, int endIndex, int cpuid);
      //Read PDIS
      bool read_PDIS(vector<string> &hash_lines);
      //
      void gruneisen_parameter_300K();
      //calculate Gruneisen parameters
      double calculate_gruneisen_with_freq_derivative(const double fp, const double fm, const double f0);
      //
      bool exists_test0 (const std::string& name);
  };
}
// ***************************************************************************
//It computes temperature dependent pdis
namespace apl
{
  class T_spectra_SCQHA_QHA3P:public MVops
  { //PN180705
    protected:
      PhononCalculator& _pc;
      QHA_AFLOWIN_CREATOR&  _runeos;
      Logger& _logger;
    public:
      T_spectra_SCQHA_QHA3P(PhononCalculator&, QHA_AFLOWIN_CREATOR&, Logger&);
      ~T_spectra_SCQHA_QHA3P();
      void clear();
    private:
      double _cutoff_freq;
      bool _is_negative_freq;
      bool _is_vol_err;
      string _tmp_dir;
      //user input data
      std::vector<std::vector<double> > _TV;
      //dynamical matrix at positive volume distortion
      vector<xmatrix<xcomplex<double> > > _DMp;
      //dynamical matrix at negative volume distortion
      vector<xmatrix<xcomplex<double> > > _DMm;
      //phonon branches 
      uint  _nBranches;
      //kpoints 
      std::vector< aurostd::xvector<double> > _kpoints;
      //frequecies along high symmetry path in THz
      vector<xvector< double> >  _freqs_path;
      //frequecies at negative volume distortion
      vector<xvector< double> >  _freqs_pathM;
      //frequecies at positive volume distortion
      vector<xvector< double> >  _freqs_pathP;
      //temperature dependent frequecies
      vector<xvector< double> >  _freqs_T;
      //The taylor 1st coefficients
      vector<xvector< double> >  _d1fdv1;
      //The taylor 2nd coefficients
      vector<xvector< double> >  _d2fdv2;

      //check Gruneisen calculation pass of fail
      vector<bool> _gp_path_test;
      //phonon directories
      vector<string> _qha_gpdir;
      //phonon volumes 
      vector<double> _qha_gpvol;
      //V-V0
      double _delta_V;
      //Relaxed volume, V0
      double _V0;
    public:
      //
      void get_tmp_dir_name(const string);
      //calculate frequencies  
      bool calculation_freqs(const std::vector< aurostd::xvector<double> > &kpoints);
      //
      bool set_imported_variables();
      //set frequency cutoff
      void set_cutoff_freq(double a){_cutoff_freq=a;}
      //Write PDIS
      void write_T_dispersion(double , double, const vector <double> &path, const vector <int> &path_seg);
      //
      void get_input_data(const std::vector<std::vector<double> >&m);
      //calculate temperature dependent pdis
      bool calculate_pdis_T(const vector <double> &path, const vector <int> &path_seg);
    private:
      //calculate frequencies
      bool calculation_freqs();
      //get dynamical matrices
      bool get_dynamicalmatrices_along_path();
      //
      //read dynamical matrices from different distortions
      bool read_matrix( vector<xmatrix<xcomplex<double> > >&A, const string file);
      //
      template<typename T> std::vector<T> split(const std::string& line);
      //calculate frequencies in threads
      void get_freqs_in_threads(int startIndex, int endIndex, int cpuid);
      //
      bool read_PDIS(vector<string> &hash_lines);
      //
      template <typename T>
        string NumToStr ( T Number );
      bool exists_test0 (const std::string& name);
  };
}
// ***************************************************************************
namespace apl
{
  class QH_ENERGIES
  { //PN180705
    private:
      //[CO20181202 - NOT USED]IPhononCalculator& _pc;
      QHA_AFLOWIN_CREATOR&  _runeos;
      Logger& _logger;
    public:
      //[CO20181202 - NOT USED]QH_ENERGIES(IPhononCalculator&, QHA_AFLOWIN_CREATOR&, Logger&);
      QH_ENERGIES(QHA_AFLOWIN_CREATOR&, Logger&);
      ~QH_ENERGIES();
      void clear();
    private:
      bool _is_magnetic;
      string _tmp_dir;
      //phonon configuration volumes
      vector<double> _ph_vols;
      //phonon configuration name
      vector<string> _ph_dirs;
      //electronic confogurationvolumes
      vector<double> _ele_vols;
      //electronic confoguration name
      vector<string> _ele_dirs;
      //index of equilibrium configuration
      int _eqm_ele_dir_index;
      //phonon dos
      vector<vector<vector<double> > > _pdos;
      //electronic dos
      vector<vector<vector<double> > > _edos;
      //corrected electronic dos
      vector<vector<vector<double> > > _cedos;
      //static energies
      vector<double> _eo;
      //Fermi energies
      vector<double> _fermi_energies;
      //pV energies
      vector<double> _pV;
      //magnetic cell
      vector<double> _mag_cell;
      //imaginary frequency index
      vector<uint> _index_imag_freq;
      //atomic species
      vector<string> _atomic_species;
    private:
      //get phonon dos
      bool get_pdos();
      //get electronic dos
      bool get_edos();
      //get phonon dos
      vector<vector<double> > get_pdos(const string file);
      //get electronic dos
      vector<vector<double> > get_edos(const string file, double &fermi);
      template<typename T>
        std::vector<T> split(const std::string& line);
      //get static energies 
      bool getE0K();
      //get static energies 
      double getE0K(const string file, double &pv, double &mag);
      void get_imaginary_freq_index();
      //remove imaginary frequencies from list
      bool remove_imaginary_freqs();
      //Write configurations having imaginary frequencies
      void print_imaginary_freq_msg();
      //chech sizes of vectors
      bool size_check();
      //chech sizes of SCQHA vectors
      bool size_scqha_check();
      //Write static energies
      bool write_energies();
      //Write configuration name after removing imaginary frequencies
      bool write_imag_freq_corrected_energies(const vector<double> &v);
      //bool eigenvalue_read(const string file, vector<vector< double> > &band_energies, vector<double> &weights);
      void get_electronic_corrected_energies();
      template<class T> vector<uint> sorted_order (const vector<T> & arr);
      bool exists_test0 (const std::string& name);
    public:
      //get QHA static energies
      bool get_qha_energies();
      //get SCQHA static energies
      bool get_scqha_energies();
      //get temporary directory name
      void get_tmp_dir_name(const string);
      //interface functions
      //get imaginary frequency index
      vector<uint>    get_index_imag_freq();
      //get magnetic cell
      vector<double>  get_mag_cell();
      //get pV energies
      vector<double>  get_pV();
      //get fremi energies
      vector<double>  get_fermi_energies();
      //get static energies
      vector<double>  get_eo();
      //get electronic dos
      vector<vector<vector<double> > >  get_edos_data();
      //get nopise free electronic dos
      vector<vector<vector<double> > >  get_cedos_data();
      //get phonon dos
      vector<vector<vector<double> > >  get_pdos_data();
      //get primitive cell volumes
      vector<double>  get_ele_vols();
      //magnetic check
      bool get_is_magnetic();
      //get relaxed volume index
      int  get_eqm_ele_dir_index();
      //get xtracture
      void get_xtracture(const xstructure& xs);
      //get atomic species
      vector<string> get_atomic_species();
  };
}
// ***************************************************************************
namespace apl
{
  class QHAEOS
  { //PN180705
    private:
      QHA&          _qha;
      QH_ENERGIES&  _qhen;
      Logger&       _logger;
    public:
      QHAEOS(QHA&, QH_ENERGIES&, Logger&);
      ~QHAEOS();
      void clear();

    private:
      int _eqm_ele_dir_index;
      bool _is_magnetic;
      bool _include_ele;
      vector<vector<double> > _TF; 
      //vector of static energies
      vector<double> _eo;
      //vector of primitive cell volumes
      vector<double> _ele_vols;
      //pdos data
      vector<vector<vector<double> > > _pdos;
      //edos data
      vector<vector<vector<double> > > _edos;
      //fermi energies
      vector<double> _fermi_energies;
      //pV energies when external pressure is applied
      vector<double> _pV;
      //zero point energies
      vector<double> _zpe;
      //list of atomic species
      vector<string> _atomic_species;
      //error managments
      bool _data_read_error;
      int _nl_success_status;
      string _nl_err_msg;
      //frequency cutoff
      double _cutoff_freq;
      //lfit output
      //uncertainties in V0 from linear least squares fit
      double _luncertanity_V0;
      //uncertainties in E0 from linear least squares fit
      double _luncertanity_E0;
      //uncertainties in B0 from linear least squares fit
      double _luncertanity_B0;
      //chisquare from linear least squares fit
      double _lchisq;
      //Equilibrium V0 from linear least squares fit
      double _leqmV0;
      //Equilibrium E0 from linear least squares fit
      double _leqmE0;
      //Equilibrium B0 from linear least squares fit
      double _leqmB0;
      //nlfit output
      //uncertainties in V0 from nonlinear least squares fit
      double _uncertanity_V0;
      //uncertainties in E0 from nonlinear least squares fit
      double _uncertanity_E0;
      //uncertainties in B0 from nonlinear least squares fit
      double _uncertanity_B0;
      //uncertainties in Bp from nonlinear least squares fit
      double _uncertanity_Bp;
      //uncertainties in Bpp from nonlinear least squares fit
      double _uncertanity_Bpp;
      //chi square per degrees of freedom
      double _chisq_dof;
      //Equilibrium V0 from linear least squares fit
      double _nleqmV0;
      //Equilibrium E0 from linear least squares fit
      double _nleqmE0;
      //Equilibrium B0 from linear least squares fit
      double _nleqmB0;
      //Equilibrium Bp from linear least squares fit
      double _nleqmBp;
      //Equilibrium Bpp from linear least squares fit
      double _nleqmBpp;
      //
      double chisq_quadratic;
      double chisq_birch_murnaghan;
      uint Iteration_birch_murnaghan;
      double alamda_birch_murnaghan;
      vector<double> Uncertainties_birch_murnaghan;
      //
      string _fdfsolver_name;
      //_fitting_type options
      // (1) BM1 => Murnaghan EOS
      // (2) BM2 => Birch-Murnaghan 3rd-order EOS
      // (3) BM3 => Birch-Murnaghan 4th-order EOS
      string _fitting_type;
      //equilibrium properties
      double _Feqm, _Beqm, _Veqm, _Bp, _Bpp;
      bool check_size();
      //calculate vibrational enengies
      double VibrationEnergy(double temperature_in_kelvins, uint dir_index);
      //calculate zeropoint energies
      void getZeroPointVibrationEnergy();
      void initialize_output_variables();
      //least square fitting
      void md_lsquares_call(const xvector<double> &V, const xvector<double> &E);
      //electronic energy calculation
      double ElectronicEnergy(double temperature_in_kelvins, uint dir_index);
      //electronic specific heat
      double Electronic_Cv(double temperature_in_kelvins, uint dir_index);
      //fermi dirac distribution
      double fermi_dirac_distribution(const double delE, const double t);
      bool more_refinement(const xvector<double> &E, const xvector<double> &V,
          xvector<double> &guess, xvector<double> &out);
      //calculate specific heat
      double getIsochoricSpecificHeat(double temperature_in_kelvins, uint dir_index);
      void calculate_zpe();
    public:
      bool setvariables();
      //calculate EOS with user defined temperatures
      void cal_qheos(double USER_TP_TSTART, double USER_TP_TEND,
          double USER_TP_TSTEP, ThermalPropertiesCalculator &e);
      void set_fitting_type(string s);
      void set_include_ele(bool b);
      //calculate total enthalypy w/o electronic contribution
      void total_enthalpy();
      //calculate total enthaly 
      void enthalpy_incuding_ele(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP);
      bool exists_test0 (const std::string& name);
  };
}
// ***************************************************************************
namespace apl
{
  class SCQHAEOS
  { //PN180705
    private:
      SCQHA_QHA3P&  _scqha;
      QH_ENERGIES& _qh_energies;
      Logger& _logger;
      //external pressure
      double _pext;
      //phonon branch
      uint _nBranches;
      //equilibrium properties
      double _Eeq, _Beq, _Veq, _Bp, _dE_dV, _d2E_dV2;
      //static energies
      vector<double> _eo;
      //static volumes
      vector<double> _ele_vols;
      //scqha phonon volumes
      vector<double> _scqha_volumes;
      //frequencies at undistorted volume
      vector<xvector<double> > _freq0;
      //frequencies at negative distorted volume
      vector<xvector<double> > _freqM;
      //frequencies at positive distorted volume
      vector<xvector<double> > _freqP;
      //Taylor's 1st order coefficient
      vector<xvector<double> > _d1fdv1;
      //Taylor's 2nd order coefficient
      vector<xvector<double> > _d2fdv2;
      //q-point weights
      vector<double> _weights;
      //temperature free energy data
      vector<vector<double> > _TF;
      //data to calculate temperature dependent PDIS
      vector<vector<double> > _TV;
      //user input temperature to calculate temperature dependent PDIS
      vector<double> _inpiut_T;
      //functions
      //calculate frequency derivative w.r.t volume
      void calculate_freq_derivative();
      //calculate frequency derivative w.r.t volume in threads
      void calculate_derivative(int startIndex, int endIndex);
      //error checking
      bool check_size();
      //polynomial fitting
      void fitting();
      //nonlinear fitting
      void md_lsquares_call(const xvector<double> &V, const xvector<double> &E);
      //
      bool more_refinement(const xvector<double> &E, const xvector<double> &V,
          xvector<double> &guess, xvector<double> &out);
      //energy at volume V
      double E_V(const double volume);
      //calculate heat capacity
      double heat_capacity(const double omeg, const double temp);
      //calculate entropy
      double entropy(const double omeg, const double temp);
      //calculate free energy
      double free_energy(const double omeg, const double temp);
      //calculate derivate of EOS w.r.t volume
      void derivatives(double vol);
      //calculate internal energy
      double internal_energy(const double omeg, const double temp);
      //print
      void print_freq_taylor_cofficients();
      //total enthalpy calculation
      void total_enthalpy();
    public:
      SCQHAEOS(SCQHA_QHA3P&, QH_ENERGIES&, Logger&);
      ~SCQHAEOS();
      void clear();
      bool import_variables();
      void sccycle(double Tmin, double Tmax, double delta_T);
      vector<vector<double> > get_TV_data();
      void  set_input_temperature(const vector<double> &a);
  };
}
// ***************************************************************************
namespace apl
{
  //It computes thermodynamic properties using QHA3P method
  class QHA3POINTS
  { //PN180705
    private:
      SCQHA_QHA3P&  _scqha;
      QH_ENERGIES& _qh_energies;
      Logger& _logger;
      //external pressure
      double _pext;
      bool _include_ele;
      bool _is_magnetic;
      uint _nBranches;
      //thermodynamic properties
      double _Eeq, _Beq, _Veq, _Bp; //, _dE_dV, _d2E_dV2; //CO20180817
      //temperature free energy data
      vector<vector<double> > _TF;
      //static energies
      vector<double> _eo;
      //pV energies
      vector<double> _pV;
      //primitive cell volues
      vector<double> _ele_vols;
      //extrapolated frequencies
      vector<vector<xvector<double> > > _ep_freqs;
      //volumes of Static calculations
      vector<double> _scqha_volumes;
      //frequencies at 0 distortion
      vector<xvector<double> > _freq0;
      //frequencies at negative distortion
      vector<xvector<double> > _freqM;
      //frequencies at positive distortion
      vector<xvector<double> > _freqP;
      //Taylor's 1st order coefficient
      vector<xvector<double> > _d1fdv1;
      //Taylor's 2nd order coefficient
      vector<xvector<double> > _d2fdv2;
      //k-points weights
      vector<double> _weights;
      //corrected electronic dos
      vector<vector<vector<double> > > _edos;
      //fermi energies
      vector<double> _fermi_energies;
      //edos at fermi energy
      //vector<double> _edosATfermi;
      //atomic species
      vector<string> _atomic_species;

      //functions
      //calculate frequency derivative
      void calculate_freq_derivative();
      //calculate frequency derivative in threads
      void calculate_derivative(int startIndex, int endIndex);
      //check errors
      bool check_size();
      //fitting function
      void fitting();
      //nonlinear fitting
      void md_lsquares_call(const xvector<double> &V, const xvector<double> &E);
      //fitting function
      bool more_refinement(const xvector<double> &E, const xvector<double> &V,
          xvector<double> &guess, xvector<double> &out);
      //calculate heat capacity
      double heat_capacity(const double omeg, const double temp);
      //calculate entropy
      double entropy(const double omeg, const double temp);
      //calculate free energy
      double free_energy(const double omeg, const double temp);
      //calculate internal energy
      double internal_energy(const double omeg, const double temp);
      //calculate electronic energy
      double ElectronicEnergy(double temperature_in_kelvins, uint dir_index);
      //double appElectronicEnergy(double temperature_in_kelvins, uint dir_index);
      double fermi_dirac_distribution(const double delE, const double t);
      //
      bool exists_test0 (const std::string& name);
      //
      void enthalpy_incuding_ele(double USER_TP_TSTART, double USER_TP_TEND, double USER_TP_TSTEP);
      //
      double Electronic_Cv(double temperature_in_kelvins, uint dir_index);
    public:
      QHA3POINTS(SCQHA_QHA3P&, QH_ENERGIES&, Logger&);
      ~QHA3POINTS();
      void clear();
      bool import_variables();
      void qha3pts_temperature_loop(double Tmin, double Tmax, double delta_T, ThermalPropertiesCalculator &e);
      void set_include_ele(bool b);
      void total_enthalpy();
  };
}
// ***************************************************************************
//Functions in this class are used to calculate Gruneisen Paramater related properties//
namespace apl
{
  //this class store dynamical matrices and phonon dispersion files from other directories
  class QHAsubdirectoryData
  { //PN180705
    private:
      PhononCalculator& _pc;
      Logger& _logger;
      string   _running_dir;
      string   _tmpdirfrefix;
      uint _nBranches;
      string _devnull;
      double _gp_vol_distortion;
      double _sc_vol_distortion;
      std::vector<aurostd::xvector<double> > _kpoints;
      vector <xmatrix<xcomplex<double> > > _dm;
      std::vector<double> _weights;
      aurostd::xmatrix<double> _rlattice;
      aurostd::xmatrix<double> _klattice;

    public:
      QHAsubdirectoryData( PhononCalculator&, Logger&);
      ~QHAsubdirectoryData();
      void clear();
      void createMPmesh(int, int, int, const xstructure&);
      void create_pdispath(std::vector<xvector<double> >& qpoints);
      void setdir_prefix(const string s);
      void set_gp_vol_distortion(const double d);
      void set_sc_vol_distortion(const double d);
      bool check_GP();
      bool check_SCQHA();
      void PDOSsave();
      void DMsave();
      void create_dm();
      string
        getdir_name(string path);
      std::vector<aurostd::xvector<double> >
        get_kpoints() { return _kpoints; }
      std::vector<double>
        get_weights() { return _weights; }
      vector<string>
        directory_list(const string path);

    private:
      template <typename T>
        string NumberToString(T Number);
      void
        get_dm(int startIndex, int endIndex);
      void write_dm(string file);
      void write_kpoints(string file);
      void write_weight(string file);
      template<class T>
        vector<T> splitWdelimiter(string s, string delimiter);
      template<typename T>
        std::vector<T> split(const std::string& line);
  };
}
// ***************************************************************************
// shellhandle.h

namespace apl {

  struct ShellData {
    ShellData();
    ShellData(const ShellData& b);

    int occupation;
    int occupationCapacity;
    bool isFull;
    double radius;
    double stdevRadius;
    vector<xvector<int> > index;
    vector<deque<_atom> > atoms;
    vector<deque<_atom> > ratoms;

    ~ShellData();
    ShellData& operator=(const ShellData&);

    void free();
    void copy(const ShellData& b);
  };

  class ShellHandle {
    private:
      int _idSafeGeneratedShell;
      int _idSafeMappedShell;

      int _centralAtomID;
      double _indexReductionConstant;
      xstructure _initStructure;
      xstructure _initStructure_original;  //corey, does not include HEAVY symmetry stuff
      //deque<_atom> _initStructure_atoms_original; //corey

      xvector<int> _safeDimension;
      vector<ShellData> _shells;

    private:
      xvector<double> getFPositionItsNearestImage(const xvector<double>&,
          const xvector<double>&,
          const xmatrix<double>&);

    public:
      ShellHandle();
      ShellHandle(const xstructure&, int, int);
      ~ShellHandle();

      void clear();
      void init(const xstructure&, int, int);

      double getShellRadius(int);
      int getShell(double);

      double getSafeShellRadius();
      void setSafeShell(int);
      int getSafeShell();

      void calcShells(const xstructure&, int, int);
      void splitBySymmetry();
      void removeSplitBySymmetry();
      void addAtomToShell(int, const _atom&, bool = true);
      void mapStructure(const xstructure&, int, bool = true);

      int getLastOccupiedShell();
      int getLastRegularShell();
      int getLastFullShell();

      int getNumberOfShells();
      int getNumberOfSubshells(int);
      std::deque<_atom> getAtomsAtSameShell(int, int = 0);
      const std::deque<_atom>& getReferenceAtomsAtSameShell(int, int = 0);

      double getIndexReductionConstant() { return _indexReductionConstant; }

      void printReport(ostream&);

      //COREY added here
      void center(int);
      void center_original(void);
      //COREY added here
  };

}  // end namespace apl

// ***************************************************************************
// aplmath.h

namespace aurostd {
  // Calculate vector projection b on a
  template <class utype>
    xvector<utype> getVectorProjection(const xvector<utype>& b, const xvector<utype>& a) {
      return (a * (utype)(scalar_product(a, b) / scalar_product(a, a)));
    }
  // Calculate vector projection c on a like b on a
  template <class utype>
    xvector<utype> getModeratedVectorProjection(const xvector<utype>& c, const xvector<utype>& b, const xvector<utype>& a) {
      return (c * (utype)(scalar_product(a, b) / scalar_product(a, a)));
    }
  // Calculate vector convolution
  template <class utype>
    xvector<utype> getVectorConvolution(const xvector<utype>& a, const xvector<utype>& b) {
      xvector<utype> v(a.rows);
      for (int i = v.lrows; i <= v.urows; i++) {
        v(i) = a(i) * b(i);
      }
      return (v);
    }
}


//OBSOLETE - ME20190815 - moved to aurostd::xmatrix
//[OBSOLETE] namespace apl {
//[OBSOLETE]   //void tred2(xmatrix<xcomplex<double> >&);
//[OBSOLETE]   void zheevByJacobiRotation(xmatrix<xcomplex<double> >&, xvector<double>&, xmatrix<xcomplex<double> >&);
//[OBSOLETE] #ifdef USE_MKL
//[OBSOLETE]   void zheevMKL(xmatrix<xcomplex<double> >&, xvector<double>&, xmatrix<xcomplex<double> >&);
//[OBSOLETE] #endif
//[OBSOLETE] }

// ***************************************************************************
// cursor.h

#define cursor_moveyx(y, x) printf("\033[%d;%dH", y, x) /*Move cursor to position y,x (rows, columns) with (1,1) as origin*/
#define cursor_moveup(y) printf("\033[%dA", y)          /*Move cursor up y*/
#define cursor_movedown(y) printf("\033[%dB", y)        /*Move cursor down y*/
#define cursor_moveright(x) printf("\033[%dC", x)       /*Move cursor right x*/
#define cursor_moveleft(x) printf("\033[%dD", x)        /*Move cursor left x*/
#define cursor_store() printf("\033[s")                 /*Store current cursor position and color*/
#define cursor_restore() printf("\033[u")               /*Restore cursor position and color from cursor_store()*/
#define cursor_clear() printf("\033[2J")                /*Clear screen and leave cursor where is*/
#define cursor_clearline() printf("\033[K")             /*Clear to end of line and leave cursor where is*/
#define cursor_fore_black() printf("\033[30m")          /*Change foreground color to black*/
#define cursor_fore_red() printf("\033[31m")            /*Change foreground color to red*/
#define cursor_fore_green() printf("\033[32m")          /*Change foreground color to green*/
#define cursor_fore_orange() printf("\033[33m")         /*Change foreground color to orange*/
#define cursor_fore_blue() printf("\033[34m")           /*Change foreground color to blue*/
#define cursor_fore_magenta() printf("\033[35m")        /*Change foreground color to magenta*/
#define cursor_fore_cyan() printf("\033[36m")           /*Change foreground color to cyan*/
#define cursor_fore_yellow() printf("\033[33m\033[1m")
#define cursor_fore_white() printf("\033[37m")    /*Change foreground color to white*/
#define cursor_back_black() printf("\033[40m")    /*Change background color to black*/
#define cursor_back_red() printf("\033[41m")      /*Change background color to red*/
#define cursor_back_green() printf("\033[42m")    /*Change background color to green*/
#define cursor_back_orange() printf("\033[43m")   /*Change background color to orange*/
#define cursor_back_blue() printf("\033[44m")     /*Change background color to blue*/
#define cursor_back_magenta() printf("\033[45m")  /*Change background color to magenta*/
#define cursor_back_cyan() printf("\033[46m")     /*Change background color to cyan*/
#define cursor_back_white() printf("\033[47m")    /*Change background color to white*/
#define cursor_attr_none() printf("\033[0m")      /*Turn off all cursor attributes*/
#define cursor_attr_bold() printf("\033[1m")      /*Make test bold*/
#define cursor_attr_underline() printf("\033[4m") /*Underline text*/
#define cursor_attr_blink() printf("\033[5m")     /*Supposed to make text blink, usually bolds it instead*/
#define cursor_attr_reverse() printf("\033[7m")   /*Swap background and foreground colors*/

// OBSOLETE ME20191031 - not used
// ***************************************************************************
// xtensor.hpp
//
//namespace apl {
//
////////////////////////////////////////////////////////////////////////////////
//template <typename T, unsigned int TENSOR_ORDER>
//class xtensor {
//  T* _data;
//  int _lindex;
//  int _hindex;
//  unsigned int _indexSize;
//  unsigned long long _arraySize;
//  std::vector<unsigned long long> _precomputedOffsets;
//
// private:
//  unsigned long long getOffset(const std::vector<unsigned int>&);
//
// public:
//  xtensor(int, int);
//  ~xtensor();
//  void zero();
//  void fill(const T&);
//  T& operator()(int, ...);
//};
//
////////////////////////////////////////////////////////////////////////////////
//template <typename T, unsigned int TENSOR_ORDER>
//xtensor<T, TENSOR_ORDER>::xtensor(int index1, int index2) {
//  _hindex = index1 > index2 ? index1 : index2;
//  _lindex = index1 < index2 ? index1 : index2;
//
//  _indexSize = _hindex - _lindex + 1;
//
//  // Allocate data
//  unsigned long long arraySizeBefore;
//  _arraySize = 1ULL;
//  _precomputedOffsets.push_back(_arraySize);
//  for (_AFLOW_APL_REGISTER_ unsigned int i = 0; i < TENSOR_ORDER; i++) {
//    arraySizeBefore = _arraySize;
//    _arraySize *= _indexSize;
//    // Detect overflow
//    if (arraySizeBefore > _arraySize)
//      throw APLRuntimeError("apl::xtensor<T>::xtensor(); The setting is producing an array which can not by handled by this implementation.");
//    _precomputedOffsets.push_back(_arraySize);
//  }
//
//  try {
//    // FIX: Problem; new(std::size_t), where size_t is hardware/platform specific,
//    // for 32 bit systems it is uint = 2^32, if our array is bigger than this,
//    // there is a overflow -> use more pages of the same block _data[PAGE][2^32]???,
//    // or go to the 64 bits systems -> 2^64 (unsigned long long)
//    if (_arraySize > (1ULL << (8 * sizeof(std::size_t) - 1)))  // dont go overboard  STEFANO (-1)
//      throw APLRuntimeError("apl::xtensor<T>::xtensor(); Problem to allocate a required array. Hardware specific problem.");
//    _data = new T[(std::size_t)_arraySize];
//  } catch (std::bad_alloc& e) {
//    throw APLRuntimeError("apl::xtensor<T>::xtensor(); Bad allocation.");
//  }
//
//  // Reverse precomputed for better manipulation and shift by 1
//  std::vector<unsigned long long> temp(_precomputedOffsets.rbegin() + 1, _precomputedOffsets.rend());
//  _precomputedOffsets = temp;
//  temp.clear();
//}
//
////////////////////////////////////////////////////////////////////////////////
//template <typename T, unsigned int TENSOR_ORDER>
//xtensor<T, TENSOR_ORDER>::~xtensor() {
//  _precomputedOffsets.clear();
//  delete[] _data;
//  _data = NULL;
//  _arraySize = 0ULL;
//}
//
////////////////////////////////////////////////////////////////////////////////
//template <typename T, unsigned int TENSOR_ORDER>
//void xtensor<T, TENSOR_ORDER>::zero() {
//  // Take care! The last argument is of std::size_t, hence platform specific
//  memset(_data, 0, _arraySize * sizeof(T));
//}
//
////////////////////////////////////////////////////////////////////////////////
//template <typename T, unsigned int TENSOR_ORDER>
//void xtensor<T, TENSOR_ORDER>::fill(const T& val) {
//  for (unsigned long long i = 0ULL; i < _arraySize; i++)
//    _data[i] = val;
//}
//
////////////////////////////////////////////////////////////////////////////////
//template <typename T, unsigned int TENSOR_ORDER>
//unsigned long long xtensor<T, TENSOR_ORDER>::getOffset(const std::vector<unsigned int>& idxs) {
//  unsigned long long offset = idxs.back();
//  for (_AFLOW_APL_REGISTER_ int i = idxs.size() - 2; i >= 0; i--)
//    offset += (unsigned long long)(idxs[i] * _precomputedOffsets[i]);
//  return offset;
//}
//
////////////////////////////////////////////////////////////////////////////////
//template <typename T, unsigned int TENSOR_ORDER>
//T& xtensor<T, TENSOR_ORDER>::operator()(int idx1, ...) {
//  // Get indices and transform to form 0...max
//  va_list arguments;
//  std::vector<unsigned int> idxs;
//
//  // The 1st index
//  idxs.push_back((unsigned int)(idx1 - _lindex));
//  //    if( idxs.back() < 0 || idxs.back() > _indexSize ) throw APLRuntimeError("apl::xtensor<T>::operator(); Index out of range.");
//  if (idxs.back() > _indexSize) throw APLRuntimeError("apl::xtensor<T>::operator(); Index out of range.");  // unsigned long long cant be negative
//
//  // The rest of indices
//  va_start(arguments, idx1);
//  for (_AFLOW_APL_REGISTER_ unsigned int i = 0; i < TENSOR_ORDER - 1; i++) {
//    idxs.push_back((unsigned int)(va_arg(arguments, int) - _lindex));
//    //   if( idxs.back() < 0 || idxs.back() > _indexSize ) throw APLRuntimeError("apl::xtensor<T>::operator(); Index out of range.");
//    if (idxs.back() > _indexSize) throw APLRuntimeError("apl::xtensor<T>::operator(); Index out of range.");  // unsigned long long cant be negative
//  }
//  va_end(arguments);
//
//  // Calculate offset
//  unsigned long long offset = getOffset(idxs);
//  idxs.clear();
//
//  // Return value
//  return _data[offset];
//}
////////////////////////////////////////////////////////////////////////////////
//}
// ***************************************************************************

#endif  // _AFLOW_APL_H_

// ***************************************************************************
