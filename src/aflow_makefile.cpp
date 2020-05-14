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

#define _DEBUG_MAKEFILE_ true
#define EXCLUDE_DATA true

namespace makefile {
  void trimPath(string& filename){
    //get rid of ".", "..", etc.
    filename=aurostd::CleanFileName(filename);  //get rid of /./
    aurostd::CleanStringASCII_InPlace(filename);
    if(filename.size()>1 && filename[0]=='.' && filename[1]=='/'){filename=filename.substr(2,string::npos);} //removes ./ from beginning so j can start at 1 below
    if(filename.find("..")!=string::npos){ //might go up and down in tree structure
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
    bool mt_required=false;
    return getDependencies(_filename,files_already_explored,dfiles,mt_required);
  }
  void getDependencies(const string& _filename,vector<string>& files_already_explored,vector<string>& dfiles,bool& mt_required){
    //find #include and get files
    bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
    string soliloquy="makefile::getDependencies():";
    
    string filename=aurostd::CleanFileName(_filename);
    aurostd::CleanStringASCII_InPlace(filename);
    trimPath(filename);

    if(LDEBUG){cerr << soliloquy << " fetching dependencies for " << filename << endl;};
    if(!aurostd::FileExist(filename)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"File \""+filename+"\" does not exist",_VALUE_ILLEGAL_);}

    if(aurostd::WithinList(files_already_explored,filename)){return;}
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
    mt_required=false;
    //only check for line comments. does not work for block comments
    for(i=0;i<vlines.size();i++){
      line = vlines[i];
      loc=line.find("//");line=line.substr(0,loc);line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line); //remove comments and left/right whitespaces
      if(line.find("#include")==0 && (line.find("<")==string::npos && line.find(">")==string::npos)){ //needs to be at the beginning of the line, files with <> are standard libraries
        if(LDEBUG){cerr << soliloquy << " line with #include = \"" << line << "\"" << endl;}
        dfile=line;
        aurostd::StringSubst(dfile,"#include","");
        aurostd::StringSubst(dfile,"\"","");
        dfile=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(dfile);
        dfile=dirname+"/"+dfile;
        trimPath(dfile);  //resolve relative path as trim as possible
        if(LDEBUG){cerr << soliloquy << " dfile=" << dfile << endl;}
        dfiles.push_back(dfile);
        getDependencies(dfile,files_already_explored,dfiles); //do not propagate mt_required from sub-dependencies
      }
      if(line.find("#define")==0 && line.find("AFLOW_")!=string::npos && line.find("_MULTITHREADS_ENABLE")!=string::npos){mt_required=true;}
    }
    std::sort(dfiles.begin(),dfiles.end());dfiles.erase( std::unique( dfiles.begin(), dfiles.end() ), dfiles.end() );  //get unique set of dependent files
    if(LDEBUG){cerr << soliloquy << " dfiles=" << aurostd::joinWDelimiter(dfiles,",") << endl;}
  }
  
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void replaceMakefileDefinitions(const vector<string>& vvariables,vector<string>& vdefinitions){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string soliloquy="makefile::replaceMakefileDefinitions():";
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  bool replacement_made=false;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string::size_type loc_var_first,loc_var_last;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string variable="";
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  uint i=0,j=0;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  for(j=0;j<vdefinitions.size();j++){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    string& definitions=vdefinitions[j];
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    loc_var_first=definitions.find("$(");
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    if(loc_var_first!=string::npos){ //variable found in definition
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      if(LDEBUG){cerr << soliloquy << " vdefinitions[j=" << j << "]=\"" << definitions << "\"" << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      loc_var_last=definitions.find(')',loc_var_first);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      variable=definitions.substr(loc_var_first,loc_var_last-loc_var_first);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      aurostd::StringSubst(variable,"$(","");
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      aurostd::StringSubst(variable,")","");
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      //aurostd::StringSubst(variable,"\\\"","");  //weird functionality: $(TIME)\" needs to remove \"
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      if(LDEBUG){cerr << soliloquy << " variable inside vdefinitions[j=" << j << "]: \"" << variable << "\"" << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      if(variable.find("shell ")!=string::npos){continue;} //skip calls to shell
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      if(variable.find(":")!=string::npos && variable.find("=")!=string::npos){continue;} //skip subst_ref
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      replacement_made=false;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      for(i=0;i<vvariables.size();i++){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        if(variable==vvariables[i]){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          if(LDEBUG){cerr << soliloquy << " vdefinitions[j=" << j << "](pre )=\"" << definitions << "\"" << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          string& definitions_sub=vdefinitions[i];
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          if(LDEBUG){cerr << soliloquy << " vdefinitions[i=" << i << "]=\"" << definitions_sub << "\"" << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          if(definitions_sub.find("shell ")!=string::npos){break;} //skip replacing calls to shell
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          if(definitions_sub.find(":")!=string::npos && definitions_sub.find("=")!=string::npos){break;} //skip replacing subst_ref
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          aurostd::StringSubst(definitions,"$("+variable+")",definitions_sub);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          replacement_made=true;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          if(LDEBUG){cerr << soliloquy << " vdefinitions[j=" << j << "](post)=\"" << definitions << "\"" << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          break;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      if(replacement_made){j--;continue;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void replaceMakefileDefinitions(const vector<string>& vvariables,vector<vector<string> >& vvdefinitions){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string soliloquy="makefile::replaceMakefileDefinitions():";
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  vector<string> vdefinitions;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  for(uint j=0;j<vvdefinitions.size();j++){vdefinitions.push_back(aurostd::joinWDelimiter(vvdefinitions[j]," "));}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  replaceMakefileDefinitions(vvariables,vdefinitions);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  splitMakefileDefinitions(vdefinitions,vvdefinitions);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void splitMakefileDefinitions(const string& definitions,vector<string>& vdefinitions){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string soliloquy="makefile::splitMakefileDefinitions():";
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  //split by spaces, except for $(shell date ...)
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  vdefinitions.clear();
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  if(LDEBUG){cerr << soliloquy << " definitions=" << definitions << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string definition="";
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  uint count_start=0,count_end=0;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  //aurostd::string2tokens(definitions,vdefinitions," ");
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string::size_type loc_var_first,loc_var_middle,loc_var_last;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  loc_var_first=loc_var_last=0;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  loc_var_last=definitions.find(' ',loc_var_last+1);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  uint i=0;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  while(loc_var_last!=string::npos){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    if(LDEBUG){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      cerr << soliloquy << " loc_var_first=" << loc_var_first << endl;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      cerr << soliloquy << " loc_var_last=" << loc_var_last << endl;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    definition=definitions.substr(loc_var_first,loc_var_last-loc_var_first);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    loc_var_middle=loc_var_last+1;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    loc_var_last=definitions.find(' ',loc_var_middle);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    //check that you have enclosing $()
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    count_start=count_end=0;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    for(i=0;i<definition.size();i++){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      if(definition[i]=='$' && i+1<definition.size() && definition[i+1]=='('){count_start++;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      if(definition[i]==')'){count_end++;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    if(count_start!=count_end){continue;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    //if(definition.find("$(")!=string::npos && definition.find(")")==string::npos){continue;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    //if(definition.find("$(")==string::npos && definition.find(")")!=string::npos){continue;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    definition=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(definition); //repetita iuvant
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    if(definition.empty()){continue;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    if(LDEBUG){cerr << soliloquy << " splitting: " << definition << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    vdefinitions.push_back(definition);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    loc_var_first=loc_var_middle;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  //get end
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  definition=definitions.substr(loc_var_first,string::npos);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  definition=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(definition); //repetita iuvant
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  if(!definition.empty()){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    if(LDEBUG){cerr << soliloquy << " splitting: " << definition << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    vdefinitions.push_back(definition);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void splitMakefileDefinitions(const vector<string>& vdefinitions,vector<vector<string> >& vvdefinitions){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string soliloquy="makefile::splitMakefileDefinitions():";
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  //split by spaces, except for $(shell date ...)
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  for(uint j=0;j<vdefinitions.size();j++){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    if(LDEBUG){cerr << soliloquy << " vdefinitions[j=" << j << "]=\"" << vdefinitions[j] << "\"" << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    vvdefinitions.push_back(vector<string>(0));
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    splitMakefileDefinitions(vdefinitions[j],vvdefinitions.back());
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void readMakefileVariables(const string& directory,vector<string>& vvariables,vector<vector<string> >& vvdefinitions){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string soliloquy="makefile::readMakefileVariables():";
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  if(LDEBUG){cerr << soliloquy << " directory=" << directory << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string filename=aurostd::CleanFileName(directory+"/Makefile");
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  aurostd::CleanStringASCII_InPlace(filename);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  if(!aurostd::FileExist(filename)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"File \""+filename+"\" does not exist",_VALUE_ILLEGAL_);}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  vector<string> vlines;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  aurostd::file2vectorstring(filename,vlines);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  return readMakefileVariables(vlines,vvariables,vvdefinitions);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]void readMakefileVariables(const vector<string>& vlines,vector<string>& vvariables,vector<vector<string> >& vvdefinitions){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string soliloquy="makefile::readMakefileVariables():";
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  //find raw variables and definitions
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  vector<string> vtokens,vdefinitions;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string::size_type loc,loc_var_first,loc_var_last;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  string line="",variable="",definition="",subst_ref_string="";
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  bool found_subst_ref=false;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  uint i=0,j=0;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  for(i=0;i<vlines.size();i++){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    line = vlines[i];
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    loc=line.find("#");line=line.substr(0,loc);line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line);  //remove comments
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    if(LDEBUG){cerr << soliloquy << " line=\"" << line << "\"" << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    if(line.size()>0 && line[0]=='\t'){continue;} //a recipe
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    if(line.find("=")!=string::npos){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      //better to define definition first because you could have: TODAY=-DTODAY=\"$(TIME)\"
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      //definition
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      found_subst_ref=true;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      loc=string::npos;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      while(found_subst_ref){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        loc=line.find_last_of("=",loc);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        if(loc==string::npos){break;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        //check that it is not a substitution ref: $(AFLOW_CPP:.cpp=.o)
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        loc_var_first=line.find_last_of("$(",loc);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        if(loc_var_first==string::npos){break;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        loc_var_last=line.find(")",loc);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        if(LDEBUG){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          cerr << soliloquy << " loc=" << loc << endl;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          cerr << soliloquy << " loc_var_first=" << loc_var_first << endl;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          cerr << soliloquy << " loc_var_last=" << loc_var_last << endl;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        subst_ref_string=line.substr(loc_var_first-1,loc_var_last-loc_var_first+2);  //-1 because "$(" is length 2
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        if(LDEBUG){cerr << soliloquy << " subst_ref_string=" << subst_ref_string << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        if(subst_ref_string.find(":")==string::npos && subst_ref_string.find("=")==string::npos){break;} //not a subst_ref
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        loc-=1;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      if(loc==string::npos){continue;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      definition=line.substr(loc+1,string::npos); //+1 for equals sign
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      definition=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(definition);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      if(0&&LDEBUG){cerr << soliloquy << " definition=\"" << definition << "\"" << endl;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      //variable
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      variable=line.substr(0,loc);  //no +1 because of equals sign
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      variable=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(variable);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      if(variable.back()==':'){variable.pop_back();}  //in case the definition uses := instead of =
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      aurostd::string2tokens(variable,vtokens,"="); //double variable assignment: TODAY=-DTODAY=\"$(TIME)\"
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      for(j=0;j<vtokens.size();j++){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        variable=vtokens[j];
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        variable=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(variable);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        if(variable.empty()){continue;}
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        vvariables.push_back(variable);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        vdefinitions.push_back(definition);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        if(LDEBUG){
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          cerr << soliloquy << " variable=\"" << vvariables.back() << "\"" << endl;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]          cerr << soliloquy << " definition=\"" << vdefinitions.back() << "\"" << endl;
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]        }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]      }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  }
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  //perform replacements
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  //NB: some variables appear twice, e.g., ARCH
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  //so this replaces with the first instance (not necessarily the correct one)
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  //this does not affect current functionality, as we only care about prerequisites
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  //you are warned
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  replaceMakefileDefinitions(vvariables,vdefinitions);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  splitMakefileDefinitions(vdefinitions,vvdefinitions);
  //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]}
  
  void updateDependenciesAUROSTD(vector<string>& vdependencies){
    bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
    string soliloquy="makefile::updateDependenciesAUROSTD():";
    if(LDEBUG){cerr << soliloquy << " vdependencies(pre )=" << aurostd::joinWDelimiter(vdependencies,",") << endl;}
    //first look for AUROSTD dependencies and replace with single AUROSTD/aurostd.o
    uint i=0;
    bool found_aurostd=false;
    for(i=vdependencies.size()-1;i<vdependencies.size();i--){  //go backwards so you can erase
      if(vdependencies[i].find("AUROSTD/")!=string::npos){
        if(LDEBUG){cerr << soliloquy << " found AUROSTD dependencies" << endl;}
        vdependencies.erase(vdependencies.begin()+i);
        found_aurostd=true;
      }
    }
    if(found_aurostd){vdependencies.insert(vdependencies.begin(),"AUROSTD/aurostd.o");}
    if(LDEBUG){cerr << soliloquy << " vdependencies(post)=" << aurostd::joinWDelimiter(vdependencies,",") << endl;}
  }
  void createMakefileAFLOW(const string& directory){
    bool LDEBUG=(FALSE || _DEBUG_MAKEFILE_ || XHOST.DEBUG);
    string soliloquy="makefile::createMakefileAFLOW():";

    if(LDEBUG){cerr << soliloquy << " directory=" << directory << endl;}

    if(!aurostd::FileExist(directory+"/aflow.h")){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,directory+"/aflow.h not found",_VALUE_ILLEGAL_);} //check we are in an aflow directory
    if(!aurostd::IsDirectory(directory+"/AUROSTD")){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,directory+"/AUROSTD not found",_VALUE_ILLEGAL_);}  //check that AUROSTD exists (fast check that ANY object files exist)
    //[not a good idea, circular flow of information]if(!aurostd::FileExist(directory+"/aflow.o")){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,directory+"/aflow.o not found",_VALUE_ILLEGAL_);}  //check that aflow.o exists (fast check that ANY object files exist)

    //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]vector<string> vvariables;
    //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]vector<vector<string> > vvdefinitions;
    //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]readMakefileVariables(directory,vvariables,vvdefinitions);
    //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]if(LDEBUG){
    //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  for(uint i=0;i<vvariables.size();i++){
    //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    cerr << soliloquy << " variable  =\"" << vvariables[i] << "\"" << endl;
    //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]    cerr << soliloquy << " definition=\"" << aurostd::joinWDelimiter(vvdefinitions[i],",") << "\"" << endl;
    //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]  }
    //[CO20200508 - OBSOLETE, DON'T BUILD THE LATRINE UPSTREAM]}

    vector<string> vsubdirectories;
    vsubdirectories.push_back(directory);
    vsubdirectories.push_back(directory+"/APL");
    //DX20200512 vsubdirectories.push_back(directory+"/ANRL");
    vsubdirectories.push_back(directory+"/SQLITE");
    std::sort(vsubdirectories.begin(),vsubdirectories.end()); //sort before adding AUROSTD, which must come first
    
    //[do AUROSTD separately]vsubdirectories.insert(vsubdirectories.begin(),directory+"/AUROSTD");  //do first
    
    uint i=0,j=0;
    string dir="",file="",_file="";
    vector<string> vfs,vfiles,files_already_explored;
    vector<vector<string> > vvdependencies;
    vector<bool> vmt_required;
    bool mt_required=false;
    vector<string> vcpp_aurostd,vhpp_aurostd,vhpp_aflow; //SC variables - hack

    //do AUROSTD first - it gets compiled separately
    file=directory+"/AUROSTD/aurostd.cpp";
    trimPath(file);
    if(LDEBUG){cerr << soliloquy << " building dependency for " << file << endl;}
    vfiles.push_back(file);
    vvdependencies.push_back(vector<string>(0));files_already_explored.clear();
    mt_required=false;
    getDependencies(vfiles.back(),files_already_explored,vvdependencies.back(),mt_required); //[didn't compile for some reason]vmt_required.back());
    vmt_required.push_back(mt_required);

    //SC variables - hack
    dir=directory+"/AUROSTD";
    trimPath(dir);
    aurostd::DirectoryLS(dir,vfs);
    std::sort(vfs.begin(),vfs.end());
    for(j=0;j<vfs.size();j++){
      file=dir+"/"+vfs[j];
      if(aurostd::IsFile(file)){
        if(file.size()>4 && file.find(".cpp")==file.size()-4){vcpp_aurostd.push_back(file);}
        else if(file.size()>2 && file.find(".h")==file.size()-2){vhpp_aurostd.push_back(file);}
      }
    }

    vector<string> vaflow_deps; //populate me with missing skipped files

    bool skip_file=false; //these files do NOT get .o
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
        if(aurostd::IsFile(file)){
          //if(file.find(".cpp")!=string::npos) //NO - ignore .cpp.orig
          if(file.size()>2 && file.find(".h")==file.size()-2){vhpp_aflow.push_back(file);} //SC variables - hack
          else if(file.size()>4 && file.find(".cpp")==file.size()-4){
            skip_file=false;
            //BEGIN skipping
            //check that it's not a .cpp generated by another file
            _file=file;aurostd::StringSubst(_file,".cpp",".js");if(aurostd::FileExist(_file)){if(LDEBUG){cerr << soliloquy << " SKIPPING generated file=" << _file << endl;}skip_file=true;}
            _file=file;aurostd::StringSubst(_file,".cpp",".py");if(aurostd::FileExist(_file)){if(LDEBUG){cerr << soliloquy << " SKIPPING generated file=" << _file << endl;}skip_file=true;}
            _file=file;aurostd::StringSubst(_file,".cpp",".txt");if(aurostd::FileExist(_file)){if(LDEBUG){cerr << soliloquy << " SKIPPING generated file=" << _file << endl;}skip_file=true;}
            //check that it's not a git file
            if(file.find("_BASE_")!=string::npos){if(LDEBUG){cerr << soliloquy << " SKIPPING git file=" << file << endl;}continue;} //continue, we don't want these in vaflow_deps
            if(file.find("_REMOTE_")!=string::npos){if(LDEBUG){cerr << soliloquy << " SKIPPING git file=" << file << endl;}continue;} //continue, we don't want these in vaflow_deps
            if(file.find("_LOCAL_")!=string::npos){if(LDEBUG){cerr << soliloquy << " SKIPPING git file=" << file << endl;}continue;}  //continue, we don't want these in vaflow_deps
            if(file.find("_BACKUP_")!=string::npos){if(LDEBUG){cerr << soliloquy << " SKIPPING git file=" << file << endl;}continue;} //continue, we don't want these in vaflow_deps
            //remove anything with aflow_data
            if(EXCLUDE_DATA){if(file.find("aflow_data")!=string::npos){if(LDEBUG){cerr << soliloquy << " SKIPPING " << file << endl;}continue;}} //continue, we don't want these in vaflow_deps
            //skip misc cpp files
            if(file=="aflow_nomix.2014-01-15.cpp"){if(LDEBUG){cerr << soliloquy << " SKIPPING " << file << endl;}continue;} //continue, we don't want these in vaflow_deps  //no reason why we have this file in the build
            if(file=="aflow_xproto_gus_lib.cpp"){if(LDEBUG){cerr << soliloquy << " SKIPPING " << file << endl;}skip_file=true;}
            if(file=="aflow_test.cpp"){if(LDEBUG){cerr << soliloquy << " SKIPPING " << file << endl;}skip_file=true;}
            if(file=="aflow_matlab_funcs.cpp"){if(LDEBUG){cerr << soliloquy << " SKIPPING " << file << endl;}skip_file=true;}
            if(file=="aflow_gnuplot_funcs.cpp"){if(LDEBUG){cerr << soliloquy << " SKIPPING " << file << endl;}skip_file=true;}
            if(skip_file){
              vaflow_deps.push_back(file);  //add to aflow dependencies
              continue;
            }
            //END skipping
            if(LDEBUG){cerr << soliloquy << " building dependency for " << file << endl;}
            vfiles.push_back(file);
            vvdependencies.push_back(vector<string>(0));files_already_explored.clear();
            //[didn't compile for some reason]vmt_required.push_back(false);
            mt_required=false;
            getDependencies(vfiles.back(),files_already_explored,vvdependencies.back(),mt_required); //[didn't compile for some reason]vmt_required.back());
            //unfortunate hack that's needed
            if(file=="aflowlib_libraries.cpp"){
              if(!aurostd::WithinList(vvdependencies.back(),"aflow_matlab_funcs.cpp")){vvdependencies.back().push_back("aflow_matlab_funcs.cpp");}  //safety
              if(!aurostd::WithinList(vvdependencies.back(),"aflow_gnuplot_funcs.cpp")){vvdependencies.back().push_back("aflow_gnuplot_funcs.cpp");}  //safety
            }
            vmt_required.push_back(mt_required);
            //get rid of any dependencies with AUROSTD and replace with $(AUROSTD_OBJ)
            updateDependenciesAUROSTD(vvdependencies.back());
          }
        }
      }
    }
    //unfortunate hack that's needed
    if(!aurostd::WithinList(vaflow_deps,"aflow_matlab_funcs.cpp")){vaflow_deps.push_back("aflow_matlab_funcs.cpp");}  //safety
    if(!aurostd::WithinList(vaflow_deps,"aflow_gnuplot_funcs.cpp")){vaflow_deps.push_back("aflow_gnuplot_funcs.cpp");}  //safety
    stringstream makefile_rules_ss;
    vector<string> vobj_files;
    string obj_file="";
    for(i=0;i<vfiles.size();i++){
      if(vvdependencies[i].empty()){continue;} //we will define generic one at the end
      obj_file=vfiles[i];aurostd::StringSubst(obj_file,".cpp",".o");
      if(vfiles[i].find("aflow_data")==string::npos){ //CO20200508 - do not include anything for aflow_data in aflow dependencies
        vobj_files.push_back(obj_file);
        vaflow_deps.push_back(obj_file);
      }
      //[not a good idea, circular flow of information]if(!aurostd::FileExist(obj_file)){continue;}  //since we can only run this code once aflow is compiled, we can check that the obj_file is a real target or not
      makefile_rules_ss << obj_file << ": " << vfiles[i] << " " << aurostd::joinWDelimiter(vvdependencies[i]," ") << endl;
      makefile_rules_ss << "\t" << "$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\\\"\"$<\"\\\" $(INCLUDE) $(CCFLAGS" << (vmt_required[i] ? "_MT" : "") << ") $(OPTS" << (vmt_required[i] ? "_MT" : "") << ") $(ARCH) $< -c -o $@" << endl;  //(obj_file=="AUROSTD/aurostd.o"?"$^":"$<")
    }
    //unfortunate hack that's needed
    if(!aurostd::WithinList(vobj_files,"aflow_matlab_funcs.cpp")){vobj_files.push_back("aflow_matlab_funcs.cpp");}  //safety
    if(!aurostd::WithinList(vobj_files,"aflow_gnuplot_funcs.cpp")){vobj_files.push_back("aflow_gnuplot_funcs.cpp");}  //safety
    
    stringstream makefile_definitions_ss;
    if(vaflow_deps.size()){makefile_definitions_ss << "AFLOW_DEPS=" << aurostd::joinWDelimiter(vaflow_deps," ") << endl;}  //SC variables - hack
    if(vobj_files.size()){makefile_definitions_ss << "AFLOW_OBJS=" << aurostd::joinWDelimiter(vobj_files," ") << endl;}  //SC variables - hack
    if(vhpp_aflow.size()){makefile_definitions_ss << "AFLOW_HPPS=" << aurostd::joinWDelimiter(vhpp_aflow," ") << endl;} //SC variables - hack
    if(vcpp_aurostd.size()){makefile_definitions_ss << "AUROSTD_CPPS=" << aurostd::joinWDelimiter(vcpp_aurostd," ") << endl;} //SC variables - hack
    if(vhpp_aurostd.size()){makefile_definitions_ss << "AUROSTD_HPPS=" << aurostd::joinWDelimiter(vhpp_aurostd," ") << endl;} //SC variables - hack

    //join together
    stringstream makefile_aflow;
    makefile_aflow << makefile_definitions_ss.str() << endl;
    makefile_aflow << makefile_rules_ss.str() << endl;
    if(aurostd::FileExist(directory+"/"+"Makefile.aflow")){aurostd::file2file(directory+"/"+"Makefile.aflow",directory+"/"+"Makefile.aflow.OLD");}  //saving
    aurostd::stringstream2file(makefile_aflow,directory+"/"+"Makefile.aflow");
  }
}

#endif  // _AFLOW_MAKEFILE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *           Aflow COREY OSES - Duke University 2013-2019                  *
// *                                                                         *
// ***************************************************************************
