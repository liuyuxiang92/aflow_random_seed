// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *           Aflow COREY OSES - Duke University 2013-2019                  *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses
// corey.oses@duke.edu

#ifndef _AFLOW_MAKEFILE_CPP_
#define _AFLOW_MAKEFILE_CPP_

#include "aflow.h"
#include "aflow_makefile.h"

namespace makefile {
  void trimPath(string& filename){
    if(filename.size()>1 && filename[0]=='.' && filename[1]=='/'){filename=filename.substr(2);} //removes ./ from beginning so j can start at 1 below
    if(filename.find("..")!=string::npos){ //might go up and down in tree structure
      filename=aurostd::CleanFileName(filename);  //get rid of /./
      vector<string> vtokens;
      aurostd::string2tokens(filename,vtokens,"/");

      uint j=0;
      for(j=1;j<vtokens.size();j++){  //start with 1
        if(vtokens[j]==".."){vtokens.erase(vtokens.begin()+j-1,vtokens.begin()+j+1);} //get rid of AUROSTD/../
      }
      filename=aurostd::joinWDelimiter(vtokens,"/");
    }
  }
  void getDependencies(const string& _filename,vector<string>& files_already_explored,vector<string>& dfiles){
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy="makefile::getDependencies():";
    
    string filename=aurostd::CleanFileName(_filename);
    trimPath(filename);

    if(LDEBUG){cerr << soliloquy << " fetching dependencies for " << filename << endl;};
    if(!aurostd::FileExist(filename)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"File \""+filename+"\" does not exist",_VALUE_ILLEGAL_);}

    if(aurostd::withinList(files_already_explored,filename)){return;}
    files_already_explored.push_back(filename);

    string dirname=".";
    if(filename.find("/")!=string::npos){ //this is a directory
      vector<string> vtokens;
      aurostd::string2tokens(filename,vtokens,"/");
      vtokens.pop_back(); //remove filename
      dirname=aurostd::joinWDelimiter(vtokens,"/");
    }
    if(LDEBUG){cerr << soliloquy << " dirname=" << dirname << endl;}

    vector<string> vlines;
    aurostd::file2vectorstring(filename,vlines);
    
    uint i=0;
    string::size_type loc;
    string line="",dfile="";
    //only check for line comments. does not work for block comments
    for(i=0;i<vlines.size();i++){
      line = vlines[i];
      loc = line.find("//");
      line = line.substr(0, loc);
      if(line.find("#include")!=string::npos && (line.find("<")==string::npos && line.find(">")==string::npos)){
        dfile=line;
        aurostd::StringSubst(dfile,"#include","");
        aurostd::StringSubst(dfile,"\"","");
        dfile=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(dfile);
        dfile=dirname+"/"+dfile;
        trimPath(dfile);
        if(LDEBUG){cerr << soliloquy << " dfile=" << dfile << endl;}
        dfiles.push_back(dfile);
        getDependencies(dfile,files_already_explored,dfiles);
      }
    }
    std::sort(dfiles.begin(),dfiles.end());dfiles.erase( std::unique( dfiles.begin(), dfiles.end() ), dfiles.end() );  //get unique set
    if(LDEBUG){cerr << soliloquy << " dfiles=" << aurostd::joinWDelimiter(dfiles,",") << endl;}
  }
  void getDependencies(const string& directory,vector<string>& dfiles){
  }
}

#endif  // _AFLOW_MAKEFILE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *           Aflow COREY OSES - Duke University 2013-2019                  *
// *                                                                         *
// ***************************************************************************
