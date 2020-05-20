// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu

#ifndef _AFLOW_MAKEFILE_H_
#define _AFLOW_MAKEFILE_H_

namespace makefile {
  void getDependencies(const string& filename,vector<string>& files_already_explored,vector<string>& dfiles);
  void getDependencies(const string& filename,vector<string>& files_already_explored,vector<string>& dfiles,bool& mt_required);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void replaceMakefileDefinitions(const vector<string>& vvariables,vector<string>& vdefinitions);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void replaceMakefileDefinitions(const vector<string>& vvariables,vector<vector<string> >& vvdefinitions);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void splitMakefileDefinitions(const string& definitions,vector<string>& vdefinitions);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void splitMakefileDefinitions(const vector<string>& vdefinitions,vector<vector<string> >& vvdefinitions);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void readMakefileVariables(const string& directory,vector<string>& vvariables,vector<vector<string> >& vvdefinitions);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void readMakefileVariables(const vector<string>& vlines,vector<string>& vvariables,vector<vector<string> >& vvdefinitions);
  void updateDependenciesAUROSTD(vector<string>& vdependencies);
  void createMakefileAFLOW(const string& directory=".");
}

#endif  // _AFLOW_MAKEFILE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
