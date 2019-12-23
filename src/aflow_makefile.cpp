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
      aurostd::CleanStringASCII_InPlace(filename);
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
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="makefile::getDependencies():";
    
    string filename=aurostd::CleanFileName(_filename);
    aurostd::CleanStringASCII_InPlace(filename);
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
      line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line);
      if(line.find("#include")==0 && (line.find("<")==string::npos && line.find(">")==string::npos)){ //needs to be at the beginning of the line, files with <> are standard libraries
        if(LDEBUG){cerr << soliloquy << " line with #include = \"" << line << "\"" << endl;}
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
  void fixDependencies(vector<string>& vdependencies){
  }
  void buildDependencies(const string& directory){
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="makefile::buildDependencies():";

    if(LDEBUG){cerr << soliloquy << " directory=" << directory << endl;}

    if(!aurostd::FileExist(directory+"/aflow.h")){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,directory+"/aflow.h not found",_VALUE_ILLEGAL_);}

    vector<string> vsubdirectories;
    vsubdirectories.push_back(directory);
    vsubdirectories.push_back(directory+"/APL");
    //[NOT NEEDED - only one command for AUROSTD cpp]vsubdirectories.push_back(directory+"/AUROSTD");
    
    uint i=0,j=0;
    string dir="",file="";
    vector<string> vfs,vfiles,files_already_explored;
    vector<vector<string> > vdependencies;
    std::sort(vsubdirectories.begin(),vsubdirectories.end());
    for(i=0;i<vsubdirectories.size();i++){
      dir=vsubdirectories[i];
      trimPath(dir);
      if(LDEBUG){cerr << soliloquy << " dir=" << dir << endl;}
      if(!aurostd::IsDirectory(dir)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,dir+" directory not found",_VALUE_ILLEGAL_);}
      aurostd::DirectoryLS(dir,vfs);
      std::sort(vfs.begin(),vfs.end());
      for(j=0;j<vfs.size();j++){
        file=dir+"/"+vfs[j];
        trimPath(file);
        if(LDEBUG){cerr << soliloquy << " file=" << file << endl;}
        if(file.find(".cpp")!=string::npos){
          if(LDEBUG){cerr << soliloquy << " building dependency for " << file << endl;}
          vfiles.push_back(file);
          vdependencies.push_back(vector<string>(0));files_already_explored.clear();
          getDependencies(vfiles.back(),files_already_explored,vdependencies.back());
          //get rid of any dependencies with AUROSTD and replace with $(AUROSTD_OBJ)
        }
      }
    }
    for(i=0;i<vfiles.size();i++){
      cout << vfiles[i] << ": " << aurostd::joinWDelimiter(vdependencies[i]," ") << endl;
      cout << "\t" << "$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\\\"\"$<\"\\\" $(INCLUDE) $(CCFLAGS) $(OPTS) $(ARCH) $< -c -o $@" << endl;
    }
  }
}

#endif  // _AFLOW_MAKEFILE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *           Aflow COREY OSES - Duke University 2013-2019                  *
// *                                                                         *
// ***************************************************************************
