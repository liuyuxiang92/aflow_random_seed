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
  void buildDependencies(const string& directory=".");
}

#endif  // _AFLOW_MAKEFILE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *           Aflow COREY OSES - Duke University 2013-2019                  *
// *                                                                         *
// ***************************************************************************
