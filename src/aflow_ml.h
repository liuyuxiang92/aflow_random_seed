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

#define _VAR_THRESHOLD_STD_ 0.001
#define _Y_CORR_THRESHOLD_STD_ 0.0
#define _SELF_CORR_THRESHOLD_STD_ 0.95

namespace aflowML {
  void insertElementalProperties(const vector<string>& vproperties,const xelement::xelement& xel,vector<string>& vitems);
  void insertElementalPropertiesCCE(const vector<string>& vproperties,const xelement::xelement& xel,double M_X_bonds,double natoms_per_fu,vector<string>& vitems);
  void insertCrystalProperties(const string& structure_path,const string& anion,const vector<string>& vheaders,vector<string>& vitems,const string& e_props=_AFLOW_XELEMENT_PROPERTIES_ALL_);
  double getStatistic(const xvector<double>& xvec,const string& stat);
  void insertElementalCombinations(const vector<string>& vproperties,vector<string>& vheaders);
  void insertElementalCombinations(const vector<string>& vproperties,const xelement::xelement& xel_cation,const xelement::xelement& xel_anion,const aflowlib::_aflowlib_entry& entry,double M_X_bonds,double natoms_per_fu_cation,double natoms_per_fu_anion,vector<string>& vheaders,vector<double>& vfeatures,bool vheaders_only=false,uint count_vcols=AUROSTD_MAX_UINT);
  void getColumn(const vector<vector<string> >& table,uint icol,vector<string>& column,bool& isfloat,bool& isinteger,bool include_header=false);
  void delColumn(vector<vector<string> >& table,uint icol);
  void oneHotFeatures(vector<vector<string> >& table,const string& features_categories);
  void removeNaN(const xvector<double>& xvec,xvector<double>& xvec_new);
  void replaceNaN(xvector<double>& xvec,double val=0.0);
  void MinMaxScale(xvector<double>& xvec);
  void reduceFeatures(vector<vector<string> >& table,const string& yheader,double var_threshold=_VAR_THRESHOLD_STD_,double ycorr_threshold=_Y_CORR_THRESHOLD_STD_,double selfcorr_threshold=_SELF_CORR_THRESHOLD_STD_);
  void reduceFeatures(vector<vector<string> >& table,const string& yheader,const string& header2skip,double var_threshold=_VAR_THRESHOLD_STD_,double ycorr_threshold=_Y_CORR_THRESHOLD_STD_,double selfcorr_threshold=_SELF_CORR_THRESHOLD_STD_);
  void reduceFeatures(vector<vector<string> >& table,const string& yheader,const vector<string>& vheaders2skip,double var_threshold=_VAR_THRESHOLD_STD_,double ycorr_threshold=_Y_CORR_THRESHOLD_STD_,double selfcorr_threshold=_SELF_CORR_THRESHOLD_STD_);
  void reduceFeatures(vector<vector<string> >& table,const string& yheader,uint icol2skip,double var_threshold=_VAR_THRESHOLD_STD_,double ycorr_threshold=_Y_CORR_THRESHOLD_STD_,double selfcorr_threshold=_SELF_CORR_THRESHOLD_STD_);
  void reduceFeatures(vector<vector<string> >& table,const string& yheader,const vector<uint>& vicol2skip,double var_threshold=_VAR_THRESHOLD_STD_,double ycorr_threshold=_Y_CORR_THRESHOLD_STD_,double selfcorr_threshold=_SELF_CORR_THRESHOLD_STD_);
  string reduceEProperties(double var_threshold=_VAR_THRESHOLD_STD_,double selfcorr_threshold=_SELF_CORR_THRESHOLD_STD_);
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
