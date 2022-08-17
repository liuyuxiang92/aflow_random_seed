//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                  Simon Divilov - Duke University 2022                   *
// *                                                                         *
//****************************************************************************
// Written by Simon Divilov, 2022
// simon.divilov@duke.edu
// 
#ifndef _AFLOW_QCA_H_
#define _AFLOW_QCA_H_

#include "aflow.h"
#include "aflow_pocc.h"

#define CONC_SHIFT 0.01 // concentration shift away from 0
#define QCA_FILE_PREFIX string("aflow_qca_")
#define QCA_AFLOW_TAG string("[AFLOW_QCA]")

namespace qca {
  class QuasiChemApproxCalculator : public xStream {
    public:
      // constructors - START
      QuasiChemApproxCalculator(ostream& oss=cout);
      QuasiChemApproxCalculator(ofstream& FileMESSAGE, ostream& oss=cout);
      QuasiChemApproxCalculator(const aurostd::xoption& qca_flags, ostream& oss=cout);
      QuasiChemApproxCalculator(const aurostd::xoption& qca_flags, ofstream& FileMESSAGE, ostream& oss=cout);
      QuasiChemApproxCalculator(const QuasiChemApproxCalculator& b);
      // constructors - END
      ~QuasiChemApproxCalculator();
      const QuasiChemApproxCalculator& operator=(const QuasiChemApproxCalculator& b);
      void clear();

      bool initialized;
      _aflags m_aflags;
      bool image_only;
      bool calc_binodal;
      double cv_cluster; // UNIT: eV
      xvector<int> num_atom_cluster; // DIM: Nc
      xmatrix<double> conc_cluster; // UNIT: unitless | DIM: Nc, Ne
      xmatrix<double> num_elem_cluster; // UNIT: unitless | DIM: Nc, Ne
      xvector<double> excess_energy_cluster; // UNIT: eV | DIM: Nc
      xvector<long int> degeneracy_cluster; // DIM: Nc
      xmatrix<double> conc_macro; // UNIT: unitless | DIM: Nx, Ne
      xvector<double> temp; // UNIT: K | DIM: Nt
      xmatrix<double> prob_ideal_cluster; // DIM: Nx, Nc
      vector<xmatrix<double>> prob_cluster; // DIM: Nx, Nc, Nt
      std::pair<double, double> param_ec; // UNIT: unitless, K
      xmatrix<double> rel_s; // UNIT: unitless | DIM: Nx, Nt
      xvector<double> binodal_curve; // UNIT: K | DIM: Nx

      // initializers - START
      bool initialize(ostream& oss);
      bool initialize(ofstream& FileMESSAGE, ostream& oss);
      bool initialize();
      bool initialize(const aurostd::xoption& qca_flags, ostream& oss);
      bool initialize(const aurostd::xoption& qca_flags, ofstream& FileMESSAGE, ostream& oss);
      bool initialize(const aurostd::xoption& qca_flags);
      // initializers - END

      void printParams();
      void errorFix();
      void calculateBinodal();
      void writeData();
      void readData();
      void plotData();

    private:
      void free();
      void copy(const QuasiChemApproxCalculator& b);

      uint min_sleep;
      string print;
      bool screen_only;
      bool use_sg;
      string rootdirpath;
      string aflowlibpath;
      string plattice;
      vector<string> elements;
      int aflow_max_num_atoms;
      int max_num_atoms;
      int conc_npts;
      bool conc_curve;
      vector<double> conc_curve_range; // DIM: 2*Ne
      int temp_npts;
      vector<double> temp_range; // UNIT: K | DIM: 2
      double cv_cut; // UNIT: eV
      string alloyname;
      string rundirpath;
      vector<xstructure> vstr_aflow;
      string lat_atat;
      vector<xstructure> vstr_ce;
      vector<int> mapstr;
      xvector<int> skipstr;
      unsigned long int nelem;
      unsigned long int ncluster; // DIM: Nc
      unsigned long int nconc; // DIM: Nx
      xvector<double> beta; // UNIT: unitless | DIM: Nt
      xmatrix<double> soln0;

      void readQCAFlags(const aurostd::xoption& qca_flags);
      string getLatForATAT(bool scale=false);
      void readAFLOWXstructures();
      vector<xstructure> getAFLOWXstructuresCustom();
      vector<xstructure> getATATXstructures(const int max_num_atoms=0, bool fromfile=false);
      void calculateMapForXstructures(const vector<xstructure>& vstr1, const vector<xstructure>& vstr2);
      void generateFilesForATAT();
      void runATAT();
      void readCVCluster();
      void calculateNumAtomCluster();
      void calculateConcentrationCluster();
      void readExcessEnergyCluster();
      void setCongruentClusters();
      void calculateDegeneracyClusterSingle(vector<xstructure>::iterator& it, const vector<xstructure>::iterator& ib, const unsigned long int ic);
      void calculateDegeneracyCluster();
      void calculateConcentrationMacro();
      void calculateTemperatureRange();
      void calculateProbabilityIdealCluster();
      void checkProbabilityIdeal();
      void calculateProbabilityCluster1D(int iix, const int it);
      void calculateProbabilityClusterND(int iix, const int it);
      void calculateProbabilityCluster();
      double getProbabilityConstraint(const int it, const int ix, const int ie, const int ideq, const xvector<double>& xvar);
      void checkProbabilityEquilibrium();
      void calculateRelativeEntropyEC();
      void calculateRelativeEntropy();
      void calculateBinodalCurve();
  };
}

namespace qca {
  void quasiChemicalApprox(const aurostd::xoption& vpflow);
  void displayUsage();
}

#endif
