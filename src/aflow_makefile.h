// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *           Aflow COREY OSES - Duke University 2013-2019                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu

#ifndef _AFLOW_MAKEFILE_H_
#define _AFLOW_MAKEFILE_H_

namespace makefile {
  void getDependencies(const string& filename,vector<string>& files_already_explored,vector<string>& dfiles);
  void getDependencies(const string& filename,vector<string>& files_already_explored,vector<string>& dfiles,bool& mt_required);
  void replaceMakefileDefinitions(const vector<string>& vvariables,vector<string>& vdefinitions);
  void replaceMakefileDefinitions(const vector<string>& vvariables,vector<vector<string> >& vvdefinitions);
  void splitMakefileDefinitions(const string& definitions,vector<string>& vdefinitions);
  void splitMakefileDefinitions(const vector<string>& vdefinitions,vector<vector<string> >& vvdefinitions);
  void readMakefileVariables(const string& directory,vector<string>& vvariables,vector<vector<string> >& vvdefinitions);
  void readMakefileVariables(const vector<string>& vlines,vector<string>& vvariables,vector<vector<string> >& vvdefinitions);
  void updateDependenciesAUROSTD(vector<string>& vdependencies);
  void buildDependencies(const string& directory=".");
}

#endif  // _AFLOW_MAKEFILE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *           Aflow COREY OSES - Duke University 2013-2019                  *
// *                                                                         *
// ***************************************************************************
