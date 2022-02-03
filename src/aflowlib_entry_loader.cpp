// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
// Corey Oses
// corey.oses@duke.edu

#ifndef _AFLOWLIB_ENTRY_LOADER_CPP_
#define _AFLOWLIB_ENTRY_LOADER_CPP_

#include <memory>

#include "aflow.h"
#include "aflowlib_entry_loader.h"


#define _DEBUG_ENTRY_LOADER_ true

// DEFINITIONS
const std::string AFLUX_DIRECTIVES=string("paging(0),format(aflow)");

namespace aflowlib {
  //--------------------------------------//
  // Helper functions for the EntryLoader //
  //--------------------------------------//


  string getAFLUX(const string& query){
    return "";
  }

  string getAFLUX(const vector<string>& matchbook){
    // Performs AFLUX call based on vector of matchbook entries
    string summons = aurostd::joinWDelimiter(matchbook,",");
    return getAFLUX(summons);
  }

  //--------------------//
  // Start Entry Loader //
  //--------------------//

  EntryLoader::EntryLoader(ostream& oss) : xStream(oss),m_initialized(false) {free();;}
  EntryLoader::EntryLoader(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {free();}
  EntryLoader::EntryLoader(const aurostd::xoption& flags,ostream& oss) : xStream(oss),m_initialized(false) {free();initialize(flags);}
  EntryLoader::EntryLoader(const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {free();initialize(flags);}
  EntryLoader::EntryLoader(const EntryLoader& b) : xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC

  EntryLoader::~EntryLoader() {xStream::free();free();}
  
  const EntryLoader& EntryLoader::operator=(const EntryLoader& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  void EntryLoader::clear() {free();}  //clear PUBLIC
  void EntryLoader::free() {
    m_initialized=false;
    m_input_processed=false;
    m_elflags.clear();
    m_sinput.clear();
    m_aflags.clear();
    m_entries_unique=true;
    m_auid_list.clear();
    uint i=0,j=0;
    for(i=0;i<m_ventries.size();i++){for(j=0;j<m_ventries[i].size();j++){m_ventries[i][j].clear();}m_ventries[i].clear();}m_ventries.clear(); //clear
  }

  void EntryLoader::copy(const EntryLoader& b) {  //copy PRIVATE
    m_initialized=b.m_initialized;
    m_input_processed=b.m_input_processed;
    m_elflags=b.m_elflags;
    m_sinput=b.m_sinput;
    m_aflags=b.m_aflags;
    m_auid_list=b.m_auid_list;
    m_entries_unique=b.m_entries_unique;
    uint i=0,j=0,k=0;
    for(i=0;i<m_ventries.size();i++){for(j=0;j<m_ventries[i].size();j++){m_ventries[i][j].clear();}m_ventries[i].clear();}m_ventries.clear(); //clear
    for(i=0;i<b.m_ventries.size();i++){for(j=0;j<b.m_ventries[i].size();j++){for(k=0;k<b.m_ventries[i][j].size();k++){m_ventries[i][j][k]=b.m_ventries[i][j][k];}}} //copy
  }
  
  bool EntryLoader::initialize(ostream& oss) {
    xStream::initialize(oss);
    return initialize();
  }

  bool EntryLoader::initialize(ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize();
  }
  
  bool EntryLoader::initialize() {
    free();
    loadAFlags();
    m_initialized=false;  //no entry
    return m_initialized;
  }
  
  bool EntryLoader::initialize(const aurostd::xoption& flags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(flags);
  }

  bool EntryLoader::initialize(const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(flags);
  }
  
  bool EntryLoader::initialize(const aurostd::xoption& flags) {
    free();
//    loadAFlags();
//    loadInput(flags);
    m_initialized=false;  //no point
    return m_initialized;
  }


  void EntryLoader::loadAFlags(){ //just a placeholder for necessary logger input
    if(XHOST.vflag_control.flag("DIRECTORY_CLEAN")){m_aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");} //CO20190402
    if(m_aflags.Directory.empty() || m_aflags.Directory=="./" || m_aflags.Directory=="."){m_aflags.Directory=aurostd::getPWD()+"/";} //".";  //CO20180220 //[CO20191112 - OBSOLETE]aurostd::execute2string(XHOST.command("pwd"))
  }

//  void EntryLoader::loadInput(const string& sinput){m_sinput=sinput;}
//  void EntryLoader::loadInput(const aurostd::xoption& flags){m_elflags=flags;}
//  void EntryLoader::loadInput(const string& sinput,const aurostd::xoption& flags){loadInput(sinput);loadInput(flags);}

  std::string EntryLoader::buildAFLUXQuery(const std::map<string, string> & matchbook){
    std::string query = "";
    for (std::map<string, string> part: {matchbook, m_aflux_directives}) {
      for (std::pair <std::string, std::string> it: part) {
        query += "," + it.first;
        if (!it.second.empty()) query += "(" + it.second + ")";
      }
    }
    query.erase(0,1);
    cout << query << endl;
    return "?" + aurostd::httpPercentEncodingFull(query);
  }

  void EntryLoader::selectSource(){
    if(aurostd::FileExist(m_sqlite_file)){
      m_current_source = Source::SQLITE;
      if(m_sqlite_db_ptr == nullptr){
        cout << "Init DB" << endl;
        m_sqlite_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_file);
      }
    } else if ("AFLUXtest"==aurostd::httpGet("http://aflowlib.duke.edu/test/?echo=AFLUXtest")){
      m_current_source = Source::AFLUX;
    } else if (false){
      m_current_source = Source::FILESYSTEM;
    } else if (false){
      m_current_source = Source::RESTAPI;
    } else {
      m_current_source = Source::FAILED;
    }
    cout << (int) m_current_source << endl;
  }

  void EntryLoader::loadSqliteWhere(const std::string & where) {
    stringstream message;
    string soliloquy=XPID+"EntryLoader::loadSqliteWhere():";
    vector<string> raw_lines = m_sqlite_db_ptr->getEntrySet(where, aflow_ft);
    cout << "Raw result size " << raw_lines.size() << " for "<< where << endl;
    aflowlib::_aflowlib_entry entry;
    size_t start_size=m_lib_entries.size();
    if(m_entries_unique) {
      for (std::string line: raw_lines) {
        entry.Load(line, *p_oss);
        if (!entry.auid.empty() && (m_auid_list.find(entry.auid)==m_auid_list.end())) {
          m_lib_entries.push_back(entry);
          m_auid_list.emplace(entry.auid);
        }
      }
    } else {
      for (std::string line: raw_lines) {
        entry.Load(line, *p_oss);
        if (!entry.auid.empty()) {
          m_lib_entries.push_back(entry);
          m_auid_list.emplace(entry.auid);
        }
      }
    }

    message << "Loaded " << m_lib_entries.size()-start_size << " new entries (overall "<< m_lib_entries.size() << " | " << m_auid_list.size() << " unique)";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);


  }
  void EntryLoader::loadAFLUXMatchbook(const std::map<string, string> & matchbook){
    loadAFLUXQuery(buildAFLUXQuery(matchbook));
  }

  void EntryLoader::loadAFLUXQuery(const std::string &query) {
    string soliloquy=XPID+"EntryLoader::loadAFLUXQuery():";
    std::string output = "";

    stringstream message;
    short status = aurostd::httpGetStatus(m_aflux_server, m_aflux_path, query, output);
    //TODO throw ERROR if web fails or try different source
//    cout << output << endl << endl;

    vector<string> raw_lines;
    aurostd::string2vectorstring(output,raw_lines);

    aflowlib::_aflowlib_entry entry;
    size_t start_size=m_lib_entries.size();
    if(m_entries_unique) {
      for (std::string line: raw_lines) {
        entry.Load(line, *p_oss);
        if (!entry.auid.empty() && (m_auid_list.find(entry.auid)==m_auid_list.end())) {
          m_lib_entries.push_back(entry);
          m_auid_list.emplace(entry.auid);
        }
      }
    } else {
      for (std::string line: raw_lines) {
        entry.Load(line, *p_oss);
        if (!entry.auid.empty()) {
          m_lib_entries.push_back(entry);
          m_auid_list.emplace(entry.auid);
        }
      }
    }

    message << "Loaded " << m_lib_entries.size()-start_size << " new entries (overall "<< m_lib_entries.size() << " | " << m_auid_list.size() << " unique)";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
  }

  void EntryLoader::loadAUID(string AUID) {
    if (m_current_source==Source::NONE) selectSource();
    if (AUID.substr(0, 5) == "auid:") AUID="aflow:" + AUID.substr(5);
    else if (AUID.substr(0, 6) != "aflow:") AUID="aflow:" + AUID;

    if (m_current_source==Source::SQLITE) {
      std::string where = "auid='\"" + AUID + "\"'";
      loadSqliteWhere(where);
    } else if (m_current_source==Source::AFLUX) {
       std::map <std::string, std::string> matchbook{{"*", ""},{"auid", "'" + AUID + "'"}};
       loadAFLUXMatchbook(matchbook);
    }
  }

  void EntryLoader::loadAUID(const vector<string> &AUID) {
    if (m_current_source==Source::NONE) selectSource();
    std::map<std::string, std::string> matchbook;
    std::string AUID_combined = "";
    for (std::string AUID_single: AUID){
      if (AUID_single.substr(0, 5) == "auid:") AUID_combined+=":'aflow:" + AUID_single.substr(5);
      else if (AUID_single.substr(0, 6) != "aflow:") AUID_combined+=":'aflow:" + AUID_single;
      else AUID_combined+=":'"+AUID_single;
      AUID_combined += "'";
    }
    AUID_combined.erase(0,1);
    loadAFLUXMatchbook({{"*", ""}, {"auid", AUID_combined}});
  }

  void EntryLoader::loadAlloy(const vector <std::string> & alloy, bool recursive) {
    if (m_current_source == Source::NONE) selectSource();
    vector <std::string> alloy_clean;

    // check that all parts of alloy are actually elements
    xelement::xelement xel;
    for (std::string element: alloy) {
      if (xel.isElement(element) == 0) {
        cerr << element << " is not an element" << endl;
      } else alloy_clean.emplace_back(element);
    }

    // build list of all needed sub-alloys
    vector<vector<string>> alloy_sub_list;
    if (recursive) {
      aurostd::xcombos xc;
      for (size_t i=1; i <= alloy_clean.size(); i++) {
        xc.reset(alloy_clean.size(), i);
        vector<int> choice;
        while (xc.increment()) {
          vector<string> alloy_tmp;
          choice = xc.getCombo();
          for (size_t idx=0; idx < choice.size(); ++idx) {
            if (choice[idx]) alloy_tmp.push_back(alloy_clean[idx]);
          }
          alloy_sub_list.push_back(alloy_tmp);
        }
      }

    } else {
      alloy_sub_list.push_back(alloy_clean);
    }

    // build sorted alloy strings
    std::vector<std::string> final_alloy_list;
    for (vector<string> a : alloy_sub_list){
      std::string alloy_tmp = "";
      std::sort(a.begin(), a.end());
      for (string e : a) alloy_tmp += e;
      final_alloy_list.push_back(alloy_tmp);
    }

    // load the alloys
    if (m_current_source==Source::SQLITE) {
      std::string where = aurostd::joinWDelimiter(final_alloy_list, "' OR alloy='");
      where = "alloy='" + where + "'";
      loadSqliteWhere(where);
    } else if (m_current_source==Source::AFLUX) {
      std::string alloy_match = aurostd::joinWDelimiter(final_alloy_list, ":");
      std::map <std::string, std::string> matchbook{{"*", ""},{"alloy", alloy_match}};
      loadAFLUXMatchbook(matchbook);
    }
  }


  void EntryLoader::loadAlloy(const std::string & alloy, bool recursive) {
    vector <string> alloy_elements;
    if (alloy.find(",") != string::npos) aurostd::string2tokens(alloy, alloy_elements, ",");
    else alloy_elements = aurostd::getElements(alloy);
    loadAlloy(alloy_elements, recursive);
  }

  void EntryLoader::retrieveOutput(string& soutput){
    string soliloquy=XPID+"EntryLoader::retrieveOutput():";
    soutput.clear();
    processInput();
    //for string output, only an AFLUX call makes sense
    if(!m_elflags.flag("ENTRY_LOADER::IS_SUMMONS")){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"AFLUX summons input not found",_INPUT_ILLEGAL_);}
    soutput=aflowlib::AFLUXCall(m_elflags.getattachedscheme("AFLUX::SUMMONS"));
  }



  void EntryLoader::retrieveOutput(vector<aflowlib::_aflowlib_entry>& ventries){
    ventries.clear();
    processInput();
    loadEntries();
    cout << "Testster" << endl;
    if(!aflowlib::mergeEntries(m_ventries,ventries)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"EntryLoader::retrieveOutput()[1]:","mergeEntries() failed",_RUNTIME_ERROR_);}
    cout << m_ventries.size() << " | " << ventries.size() << endl;
  }
  void EntryLoader::retrieveOutput(vector<vector<aflowlib::_aflowlib_entry> >& ventries){
    uint i=0;
    for(i=0;i<ventries.size();i++){ventries[i].clear();}ventries.clear();
    processInput();
    loadEntries();
    if(!aflowlib::mergeEntries(m_ventries,ventries)){throw aurostd::xerror(_AFLOW_FILE_NAME_,"EntryLoader::retrieveOutput()[2]:","mergeEntries() failed",_RUNTIME_ERROR_);}
  }
  void EntryLoader::retrieveOutput(vector<vector<vector<aflowlib::_aflowlib_entry> > >& ventries){
    uint i=0,j=0;
    for(i=0;i<ventries.size();i++){for(j=0;j<ventries[i].size();j++){ventries[i][j].clear();}ventries[i].clear();}ventries.clear();
    processInput();
    loadEntries();
    for(i=0;i<m_ventries.size();i++){ventries.push_back(vector<vector<aflowlib::_aflowlib_entry> >(0));for(j=0;j<m_ventries[i].size();j++){ventries.back().push_back(m_ventries[i][j]);}}
  }
  void EntryLoader::processInput(){
    bool LDEBUG=(FALSE || _DEBUG_ENTRY_LOADER_ || XHOST.DEBUG);
    string soliloquy=XPID+"EntryLoader::processInput():";
    stringstream message;

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    message << "Processing input=\"" << m_sinput << "\"";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    m_input_processed=false;
    
    //look for an aflux summons first
    if(m_input_processed==false){
      vector<string> vkeys_aflux=getDataNames();
      if(LDEBUG){
        for(uint i=0;i<vkeys_aflux.size();i++){
          cerr << soliloquy << " vkeys_aflux[" << i << "]=" << vkeys_aflux[i] << endl;
        }
      }
      for(uint i=0;i<vkeys_aflux.size()&&m_input_processed==false;i++){
        if(m_sinput.find(vkeys_aflux[i])!=string::npos){
          if(LDEBUG){cerr << soliloquy << " found aflux_key in m_sinput: " << vkeys_aflux[i] << endl;}
          //treating input as aflux summons
          m_input_processed=true;
          m_elflags.flag("ENTRY_LOADER::IS_SUMMONS",TRUE);
          m_elflags.push_attached("AFLUX::SUMMONS",m_sinput);
          message << "Input is AFLUX summons";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
        }
      }
    }

    //see if input is a bunch of elements
    if(m_input_processed==false){
      vector<string> velements;
      if(m_sinput.find(",")!=string::npos){aurostd::string2tokens(m_sinput,velements,",");}
      else{velements=aurostd::getElements(m_sinput);}
      xelement::xelement xel;
      m_input_processed=true; //unset later if false
      for(uint i=0;i<velements.size()&&m_input_processed==true;i++){
        if(xel.isElement(velements[i])==0){
          if(LDEBUG){cerr << soliloquy << " " << velements[i] << " is not an element" << endl;}
          m_input_processed=false;
        }
      }
      if(m_input_processed){
        if(LDEBUG){cerr << soliloquy << " found elements string in m_sinput: " << m_sinput << endl;}
        //treat as string of elements
        m_input_processed=true;
        m_elflags.flag("ENTRY_LOADER::IS_ELEMENTS_STRING",true);
        m_elflags.push_attached("ENTRY_LOADER::ELEMENTS_STRING",aurostd::joinWDelimiter(velements,""));  //input clean elements string
        message << "Input is elements string";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);
      }
    }

    if(m_input_processed==false){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Input cannot be processed",_INPUT_ILLEGAL_);}
  }
  void EntryLoader::sanitizeAFLUXSummons(){ //n from AFLUX paper
    bool LDEBUG=(FALSE || _DEBUG_ENTRY_LOADER_ || XHOST.DEBUG);
    string soliloquy=XPID+"EntryLoader::sanitizeAFLUXSummons():";
    stringstream message;

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    string summons=m_elflags.getattachedscheme("AFLUX::SUMMONS");
    if(summons.empty()){return;}

    uint i=0;
    vector<string> vtokens,vtokens_new;
    aurostd::string2tokens(summons,vtokens,",");
    //patch paging(n,k)
    for(i=0;i<vtokens.size();i++){
      if(vtokens[i].find("paging(")!=string::npos){
        if(vtokens[i].find(")")==string::npos){
          if(i>=vtokens.size()-1){
            throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"paging directive incomplete [1]",_INPUT_ILLEGAL_);
          }
          string paging=vtokens[i]+","+vtokens[i+1];
          if(paging.find(")")==string::npos){
            throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"paging directive incomplete [2]",_INPUT_ILLEGAL_);
          }
          vtokens[i]=paging;
          vtokens.erase(vtokens.begin()+i+1); //erase second piece
        }
      }
    }


    bool found_format=false,found_paging=false,found_all=false;
    for(i=0;i<vtokens.size();i++){
      //format
      if(vtokens[i].find("format(")!=string::npos){
        found_format=true;
        if(vtokens[i]!="format(aflow)"){
          vtokens_new.push_back("format(aflow)");
          message << "Converting \""+vtokens[i]+"\" to \"format(aflow)\"";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
          continue;  //skip
        }
      }
      //paging
      if(vtokens[i].find("paging(")!=string::npos){
        found_paging=true;
        if(vtokens[i]!="paging(0)"){  //paging(0) //"paging("+aurostd::utype2string(n)+",1000)"
          vtokens_new.push_back("paging(0)");
          message << "Converting \""+vtokens[i]+"\" to \"paging(0)\"";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
          continue;  //skip
        }
      }
      //all
      if(vtokens[i]=="*"){found_all=true;}
      vtokens_new.push_back(vtokens[i]);
    }
    if(!found_format){
      vtokens_new.push_back("format(aflow)");
      message << "Adding \"format(aflow)\" to summons";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
    }
    if(!found_paging){
      vtokens_new.push_back("paging(0)");
      message << "Adding \"paging(0)\" to summons";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
    }

    if(!found_all){
      //vtokens_new.push_back("*");
      vtokens_new.insert(vtokens_new.begin(),"*");  //insert at the beginning, better syntax for aflux
      message << "Adding \"*\" to summons";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
    }

    string summons_new=aurostd::joinWDelimiter(vtokens_new,",");
    message << "New summons=\""+summons_new+"\"";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    m_elflags.pop_attached("AFLUX::SUMMONS");
    m_elflags.push_attached("AFLUX::SUMMONS",summons_new);
  }
  uint EntryLoader::AFLUXSummons2Entries(){
    bool LDEBUG=(FALSE || _DEBUG_ENTRY_LOADER_ || XHOST.DEBUG);
    string soliloquy=XPID+"EntryLoader::AFLUXSummons2Entries():";
    stringstream message;
    
    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    string summons=m_elflags.getattachedscheme("AFLUX::SUMMONS");
    if(summons.empty()){return 0;}
    if(LDEBUG){cerr << soliloquy << " summons=\"" << summons << "\"" << endl;}
    string soutput=aflowlib::getAFLUX(m_elflags.getattachedscheme("AFLUX::SUMMONS"));
//    if(LDEBUG){
//      cerr << soliloquy << " response:" << endl;
//      cerr << soutput;
//      cerr << endl;
//    }

    vector<string> vlines;
    aurostd::string2vectorstring(soutput,vlines);
   
    uint i=0;
    aflowlib::_aflowlib_entry entry;

    vector<aflowlib::_aflowlib_entry> ventries; //tmp container
    for(i=0;i<vlines.size();i++){
      entry.Load(vlines[i],*p_oss);
      if(!entry.auid.empty()){ventries.push_back(entry);}
    }
    
    message << "Loaded " << ventries.size() << " entries";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);

    if(!aflowlib::mergeEntries(ventries,m_ventries)){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"mergeEntries() failed",_RUNTIME_ERROR_);}

    return ventries.size(); //count of new entries
  }
  void EntryLoader::setAFLUXSummons4ElementsString(){
    bool LDEBUG=(FALSE || _DEBUG_ENTRY_LOADER_ || XHOST.DEBUG);
    string soliloquy=XPID+"EntryLoader::setAFLUXSummons4ElementsString():";
    stringstream message;

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    string elements_string=m_elflags.getattachedscheme("ENTRY_LOADER::ELEMENTS_STRING");
    if(LDEBUG){cerr << soliloquy << " elements_string=" << elements_string << endl;}
    if(elements_string.empty()){return;}

    vector<string> velements=aurostd::getElements(elements_string);
    if(LDEBUG){cerr << soliloquy << " velements=" << aurostd::joinWDelimiter(velements,",") << endl;}
    uint nary=velements.size();

    //?((species(Mn:Pd:Pt),nspecies(1)):(species((Mn,Pd):(Mn,Pt):(Pd,Pt)),nspecies(2))),stoichiometry,enthalpy_formation_atom,paging(1,1000),format(json)
    vector<string> vtokens_naries;  //to be joined with ':' between naries
    vector<string> vtokens_elements;  //to be joined with ':' inside species()

    //do unary separately
    vtokens_naries.push_back( "nspecies(1),species("+aurostd::joinWDelimiter(velements,":")+")" );  //always put nspecies first, faster in sqlite

    uint i=0;
    aurostd::xcombos xc;
    for(i=1;i<nary;i++){  //start from binaries
      xc.reset(nary,i+1,'C');
      vtokens_elements.clear();
      while(xc.increment()){
        vtokens_elements.push_back( aurostd::joinWDelimiter(xc.applyCombo(velements),",") );
      }
      if(vtokens_elements.size()>1){vtokens_elements=aurostd::wrapVecEntries(vtokens_elements,"(",")");}
      vtokens_naries.push_back( "nspecies("+aurostd::utype2string(i+1)+"),species("+aurostd::joinWDelimiter(vtokens_elements,":")+")" );  //always put nspecies first, faster in sqlite
    }
    if(vtokens_naries.size()>1){vtokens_naries=aurostd::wrapVecEntries(vtokens_naries,"(",")");}

    string summons=aurostd::joinWDelimiter(vtokens_naries,":");
    if(LDEBUG){cerr << soliloquy << " summons=" << summons << endl;}
    
    m_elflags.pop_attached("AFLUX::SUMMONS");
    m_elflags.push_attached("AFLUX::SUMMONS",summons);
  }
  void EntryLoader::loadEntries(){
    bool LDEBUG=(FALSE || _DEBUG_ENTRY_LOADER_ || XHOST.DEBUG);
    string soliloquy=XPID+"EntryLoader::loadEntries():";
    stringstream message;
    
    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    if(m_input_processed==false){throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Input not processed",_INPUT_ILLEGAL_);}

    uint i=0,j=0;
    for(i=0;i<m_ventries.size();i++){for(j=0;j<m_ventries[i].size();j++){m_ventries[i][j].clear();}m_ventries[i].clear();}m_ventries.clear(); //clear
    
    if(m_elflags.flag("ENTRY_LOADER::IS_SUMMONS")){
      //uint n=1; //n from AFLUX paper
      //uint nentries=0;
      //while(true)
      //sanitizeAFLUXSummons((n++));
      //nentries=AFLUXSummons2Entries();
      //if(nentries==0){break;}
      sanitizeAFLUXSummons();
      AFLUXSummons2Entries();
      return;
    }
    else if(m_elflags.flag("ENTRY_LOADER::IS_ELEMENTS_STRING")){
      //assumes our input is chull-type: i.e., combinations of unary, binary, ternary, etc.
      setAFLUXSummons4ElementsString();
      sanitizeAFLUXSummons();
      AFLUXSummons2Entries();
      //uint n=1; //n from AFLUX paper
      //uint nentries=0;
      //while(true){
      //  sanitizeAFLUXSummons((n++));
      //  nentries=AFLUXSummons2Entries();
      //  if(nentries==0){break;}
      //}
    }
    else{throw aurostd::xerror(_AFLOW_FILE_NAME_,soliloquy,"Unknown input format",_INPUT_ILLEGAL_);}
  }



}

#endif  // _AFLOW_ENTRY_LOADER_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
