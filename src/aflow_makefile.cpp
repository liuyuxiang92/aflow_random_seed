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
    mt_required=false;
    //only check for line comments. does not work for block comments
    for(i=0;i<vlines.size();i++){
      line = vlines[i];
      loc=line.find("//");line=line.substr(0,loc);line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line); //remove comments
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
        getDependencies(dfile,files_already_explored,dfiles); //do not propagate mt_required from sub-dependencies
      }
      if(line.find("#define")==0 && line.find("AFLOW_")!=string::npos && line.find("_MULTITHREADS_ENABLE")!=string::npos){mt_required=true;}
    }
    std::sort(dfiles.begin(),dfiles.end());dfiles.erase( std::unique( dfiles.begin(), dfiles.end() ), dfiles.end() );  //get unique set
    if(LDEBUG){cerr << soliloquy << " dfiles=" << aurostd::joinWDelimiter(dfiles,",") << endl;}
  }
  void replaceMakefileDefinitions(const vector<string>& vvariables,vector<string>& vdefinitions){
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy="makefile::replaceMakefileDefinitions():";

    bool replacement_made=false;
    string::size_type loc_var_first,loc_var_last;
    string variable="";
    uint i=0,j=0;
    for(j=0;j<vdefinitions.size();j++){
      string& definitions=vdefinitions[j];
      loc_var_first=definitions.find("$(");
      if(loc_var_first!=string::npos){ //variable found in definition
        if(LDEBUG){cerr << soliloquy << " vdefinitions[j=" << j << "]=\"" << definitions << "\"" << endl;}
        loc_var_last=definitions.find(')',loc_var_first);
        variable=definitions.substr(loc_var_first,loc_var_last-loc_var_first);
        aurostd::StringSubst(variable,"$(","");
        aurostd::StringSubst(variable,")","");
        //aurostd::StringSubst(variable,"\\\"","");  //weird functionality: $(TIME)\" needs to remove \"
        if(LDEBUG){cerr << soliloquy << " variable inside vdefinitions[j=" << j << "]: \"" << variable << "\"" << endl;}
        if(variable.find("shell ")!=string::npos){continue;} //skip calls to shell
        if(variable.find(":")!=string::npos && variable.find("=")!=string::npos){continue;} //skip subst_ref
        replacement_made=false;
        for(i=0;i<vvariables.size();i++){
          if(variable==vvariables[i]){
            if(LDEBUG){cerr << soliloquy << " vdefinitions[j=" << j << "](pre )=\"" << definitions << "\"" << endl;}
            string& definitions_sub=vdefinitions[i];
            if(LDEBUG){cerr << soliloquy << " vdefinitions[i=" << i << "]=\"" << definitions_sub << "\"" << endl;}
            if(definitions_sub.find("shell ")!=string::npos){break;} //skip replacing calls to shell
            if(definitions_sub.find(":")!=string::npos && definitions_sub.find("=")!=string::npos){break;} //skip replacing subst_ref
            aurostd::StringSubst(definitions,"$("+variable+")",definitions_sub);
            replacement_made=true;
            if(LDEBUG){cerr << soliloquy << " vdefinitions[j=" << j << "](post)=\"" << definitions << "\"" << endl;}
            break;
          }
        }
        if(replacement_made){j--;continue;}
      }
    }
  }
  void replaceMakefileDefinitions(const vector<string>& vvariables,vector<vector<string> >& vvdefinitions){
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy="makefile::replaceMakefileDefinitions():";
    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    vector<string> vdefinitions;
    for(uint j=0;j<vvdefinitions.size();j++){vdefinitions.push_back(aurostd::joinWDelimiter(vvdefinitions[j]," "));}
    replaceMakefileDefinitions(vvariables,vdefinitions);
    splitMakefileDefinitions(vdefinitions,vvdefinitions);
  }
  void splitMakefileDefinitions(const string& definitions,vector<string>& vdefinitions){
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy="makefile::splitMakefileDefinitions():";

    //split by spaces, except for $(shell date ...)
    vdefinitions.clear();
    if(LDEBUG){cerr << soliloquy << " definitions=" << definitions << endl;}
    string definition="";
    uint count_start=0,count_end=0;
    //aurostd::string2tokens(definitions,vdefinitions," ");
    string::size_type loc_var_first,loc_var_middle,loc_var_last;
    loc_var_first=loc_var_last=0;
    loc_var_last=definitions.find(' ',loc_var_last+1);
    uint i=0;
    while(loc_var_last!=string::npos){
      if(LDEBUG){
        cerr << soliloquy << " loc_var_first=" << loc_var_first << endl;
        cerr << soliloquy << " loc_var_last=" << loc_var_last << endl;
      }
      definition=definitions.substr(loc_var_first,loc_var_last-loc_var_first);
      loc_var_middle=loc_var_last+1;
      loc_var_last=definitions.find(' ',loc_var_middle);
      //check that you have enclosing $()
      count_start=count_end=0;
      for(i=0;i<definition.size();i++){
        if(definition[i]=='$' && i+1<definition.size() && definition[i+1]=='('){count_start++;}
        if(definition[i]==')'){count_end++;}
      }
      if(count_start!=count_end){continue;}
      //if(definition.find("$(")!=string::npos && definition.find(")")==string::npos){continue;}
      //if(definition.find("$(")==string::npos && definition.find(")")!=string::npos){continue;}
      definition=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(definition); //repetita iuvant
      if(definition.empty()){continue;}
      if(LDEBUG){cerr << soliloquy << " splitting: " << definition << endl;}
      vdefinitions.push_back(definition);
      loc_var_first=loc_var_middle;
    }
    //get end
    definition=definitions.substr(loc_var_first,string::npos);
    definition=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(definition); //repetita iuvant
    if(!definition.empty()){
      if(LDEBUG){cerr << soliloquy << " splitting: " << definition << endl;}
      vdefinitions.push_back(definition);
    }
  }
  void splitMakefileDefinitions(const vector<string>& vdefinitions,vector<vector<string> >& vvdefinitions){
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy="makefile::splitMakefileDefinitions():";
    
    //split by spaces, except for $(shell date ...)
    for(uint j=0;j<vdefinitions.size();j++){
      if(LDEBUG){cerr << soliloquy << " vdefinitions[j=" << j << "]=\"" << vdefinitions[j] << "\"" << endl;}
      vvdefinitions.push_back(vector<string>(0));
      splitMakefileDefinitions(vdefinitions[j],vvdefinitions.back());
    }
  }
  void readMakefileVariables(const string& directory,vector<string>& vvariables,vector<vector<string> >& vvdefinitions){
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy="makefile::readMakefileVariables():";

    if(LDEBUG){cerr << soliloquy << " directory=" << directory << endl;}
    
    string filename=aurostd::CleanFileName(directory+"/Makefile");
    aurostd::CleanStringASCII_InPlace(filename);
    if(!aurostd::FileExist(filename)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"File \""+filename+"\" does not exist",_VALUE_ILLEGAL_);}
    vector<string> vlines;
    aurostd::file2vectorstring(filename,vlines);
    return readMakefileVariables(vlines,vvariables,vvdefinitions);
  }
  void readMakefileVariables(const vector<string>& vlines,vector<string>& vvariables,vector<vector<string> >& vvdefinitions){
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy="makefile::readMakefileVariables():";

    //find raw variables and definitions
    vector<string> vtokens,vdefinitions;
    string::size_type loc,loc_var_first,loc_var_last;
    string line="",variable="",definition="",subst_ref_string="";
    bool found_subst_ref=false;
    uint i=0,j=0;
    for(i=0;i<vlines.size();i++){
      line = vlines[i];
      loc=line.find("#");line=line.substr(0,loc);line=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(line);  //remove comments
      if(LDEBUG){cerr << soliloquy << " line=\"" << line << "\"" << endl;}
      if(line.size()>0 && line[0]=='\t'){continue;} //a recipe
      if(line.find("=")!=string::npos){
        //better to define definition first because you could have: TODAY=-DTODAY=\"$(TIME)\"
        //definition
        found_subst_ref=true;
        loc=string::npos;
        while(found_subst_ref){
          loc=line.find_last_of("=",loc);
          if(loc==string::npos){break;}
          //check that it is not a substitution ref: $(AFLOW_CPP:.cpp=.o)
          loc_var_first=line.find_last_of("$(",loc);
          if(loc_var_first==string::npos){break;}
          loc_var_last=line.find(")",loc);
          if(LDEBUG){
            cerr << soliloquy << " loc=" << loc << endl;
            cerr << soliloquy << " loc_var_first=" << loc_var_first << endl;
            cerr << soliloquy << " loc_var_last=" << loc_var_last << endl;
          }
          subst_ref_string=line.substr(loc_var_first-1,loc_var_last-loc_var_first+2);  //-1 because "$(" is length 2
          if(LDEBUG){cerr << soliloquy << " subst_ref_string=" << subst_ref_string << endl;}
          if(subst_ref_string.find(":")==string::npos && subst_ref_string.find("=")==string::npos){break;} //not a subst_ref
          loc-=1;
        }
        if(loc==string::npos){continue;}
        definition=line.substr(loc+1,string::npos); //+1 for equals sign
        definition=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(definition);
        if(0&&LDEBUG){cerr << soliloquy << " definition=\"" << definition << "\"" << endl;}
        //variable
        variable=line.substr(0,loc);  //no +1 because of equals sign
        variable=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(variable);
        if(variable.back()==':'){variable.pop_back();}  //in case the definition uses := instead of =
        aurostd::string2tokens(variable,vtokens,"="); //double variable assignment: TODAY=-DTODAY=\"$(TIME)\"
        for(j=0;j<vtokens.size();j++){
          variable=vtokens[j];
          variable=aurostd::RemoveWhiteSpacesFromTheFrontAndBack(variable);
          if(variable.empty()){continue;}
          vvariables.push_back(variable);
          vdefinitions.push_back(definition);
          if(LDEBUG){
            cerr << soliloquy << " variable=\"" << vvariables.back() << "\"" << endl;
            cerr << soliloquy << " definition=\"" << vdefinitions.back() << "\"" << endl;
          }
        }
      }
    }

    //perform replacements
    //NB: some variables appear twice, e.g., ARCH
    //so this replaces with the first instance (not necessarily the correct one)
    //this does not affect current functionality, as we only care about prerequisites
    //you are warned
    replaceMakefileDefinitions(vvariables,vdefinitions);

    splitMakefileDefinitions(vdefinitions,vvdefinitions);
  }
  void updateDependenciesAUROSTD(vector<string>& vdependencies){
    bool LDEBUG=(TRUE || XHOST.DEBUG);
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
  void buildDependencies(const string& directory){
    bool LDEBUG=(TRUE || XHOST.DEBUG);
    string soliloquy="makefile::buildDependencies():";

    if(LDEBUG){cerr << soliloquy << " directory=" << directory << endl;}

    if(!aurostd::FileExist(directory+"/aflow.h")){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,directory+"/aflow.h not found",_VALUE_ILLEGAL_);} //check we are in an aflow directory
    if(!aurostd::IsDirectory(directory+"/AUROSTD")){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,directory+"/AUROSTD not found",_VALUE_ILLEGAL_);}  //check that AUROSTD exists (fast check that ANY object files exist)
    //[not a good idea, circular flow of information]if(!aurostd::FileExist(directory+"/aflow.o")){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,directory+"/aflow.o not found",_VALUE_ILLEGAL_);}  //check that aflow.o exists (fast check that ANY object files exist)

    vector<string> vvariables;
    vector<vector<string> > vvdefinitions;
    readMakefileVariables(directory,vvariables,vvdefinitions);
    if(LDEBUG){
      for(uint i=0;i<vvariables.size();i++){
        cerr << soliloquy << " variable  =\"" << vvariables[i] << "\"" << endl;
        cerr << soliloquy << " definition=\"" << aurostd::joinWDelimiter(vvdefinitions[i],",") << "\"" << endl;
      }
    }

    vector<string> vsubdirectories;
    vsubdirectories.push_back(directory);
    vsubdirectories.push_back(directory+"/APL");
    std::sort(vsubdirectories.begin(),vsubdirectories.end()); //sort before adding AUROSTD, which must come first
    
    //[do AUROSTD separately]vsubdirectories.insert(vsubdirectories.begin(),directory+"/AUROSTD");  //do first
    
    uint i=0,j=0;
    string dir="",file="";
    vector<string> vfs,vfiles,files_already_explored;
    vector<vector<string> > vvdependencies;
    vector<bool> vmt_required;
    bool mt_required=false;

    //do AUROSTD first - it gets compiled separately
    file=directory+"/AUROSTD/aurostd.cpp";
    trimPath(file);
    if(LDEBUG){cerr << soliloquy << " building dependency for " << file << endl;}
    vfiles.push_back(file);
    vvdependencies.push_back(vector<string>(0));files_already_explored.clear();
    mt_required=false;
    getDependencies(vfiles.back(),files_already_explored,vvdependencies.back(),mt_required); //[didn't compile for some reason]vmt_required.back());
    vmt_required.push_back(mt_required);

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
          vvdependencies.push_back(vector<string>(0));files_already_explored.clear();
          //[didn't compile for some reason]vmt_required.push_back(false);
          mt_required=false;
          getDependencies(vfiles.back(),files_already_explored,vvdependencies.back(),mt_required); //[didn't compile for some reason]vmt_required.back());
          vmt_required.push_back(mt_required);
          //get rid of any dependencies with AUROSTD and replace with $(AUROSTD_OBJ)
          updateDependenciesAUROSTD(vvdependencies.back());
        }
      }
    }
    string obj_file="";
    for(i=0;i<vfiles.size();i++){
      if(vvdependencies[i].empty()){continue;} //we will define generic one at the end
      obj_file=vfiles[i];aurostd::StringSubst(obj_file,".cpp",".o");
      //[not a good idea, circular flow of information]if(!aurostd::FileExist(obj_file)){continue;}  //since we can only run this code once aflow is compiled, we can check that the obj_file is a real target or not
      cout << obj_file << ": " << vfiles[i] << " " << aurostd::joinWDelimiter(vvdependencies[i]," ") << endl;
      cout << "\t" << "$(CPP) $(VERS) -D_AFLOW_FILE_NAME_=\\\"\"$<\"\\\" $(INCLUDE) $(CCFLAGS) $(OPTS" << (vmt_required[i] ? "_MT" : "") << ") $(ARCH) $< -c -o $@" << endl;
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
