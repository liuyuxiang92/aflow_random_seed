// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu

#ifndef _AFLOW_ML_H_
#define _AFLOW_ML_H_

namespace aflowML {
  void insertElementalPropertiesCCE(const vector<string>& vproperties,const xelement::xelement& xel,vector<string>& vitems);
  void insertCrystalPropertiesCCE(const string& structure_path,const string& anion,const vector<string>& vheaders,vector<string>& vitems);
  double getStatistic(const xvector<double>& xvec,const string& stat);
  void insertElementalCombinationsCCE(const vector<string>& vproperties,vector<string>& vheaders);
  void insertElementalCombinationsCCE(const vector<string>& vproperties,const xelement::xelement& xel_cation,const xelement::xelement& xel_anion,const aflowlib::_aflowlib_entry& entry,double M_X_bonds,vector<string>& vheaders,vector<double>& vfeatures,bool vheaders_only=false);
  void writeCCECSV();
} // namespace aflowML
namespace aflowlib {
  void insertStoichStats(const vector<string> vstats,const xvector<double>& nspecies_xv,const xvector<double>& stoich_xv,vector<double>& vfeatures);
}

#endif  // _AFLOW_ML_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
