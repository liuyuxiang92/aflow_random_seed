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

#include "aflow.h"
#include "aflowlib_entry_loader.h"

#define _DEBUG_ENTRY_LOADER_ true

// DEFINITIONS
const std::string AFLUX_DIRECTIVES=string("paging(0),format(aflow)");

namespace aflowlib {
  EntryLoader::EntryLoader(ostream& oss) : xStream(oss),m_initialized(false) {free();;}
  EntryLoader::EntryLoader(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {free();}
  EntryLoader::EntryLoader(const aurostd::xoption& flags,ostream& oss) : xStream(oss),m_initialized(false) {free();initialize(flags);}
  EntryLoader::EntryLoader(const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {free();initialize(flags);}
  EntryLoader::EntryLoader(const string& sinput,ostream& oss) : xStream(oss),m_initialized(false) {free();initialize(sinput);}
  EntryLoader::EntryLoader(const string& sinput,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {free();initialize(sinput);}
  EntryLoader::EntryLoader(const string& sinput,const aurostd::xoption& flags,ostream& oss) : xStream(oss),m_initialized(false) {free();initialize(sinput,flags);}
  EntryLoader::EntryLoader(const string& sinput,const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss),m_initialized(false) {free();initialize(sinput,flags);}
  EntryLoader::EntryLoader(const EntryLoader& b) : xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC

  EntryLoader::~EntryLoader() {xStream::free();free();}
  
  const EntryLoader& EntryLoader::operator=(const EntryLoader& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  void EntryLoader::clear() {free();}  //clear PUBLIC
  void EntryLoader::free() {
    m_initialized=false;
    m_elflags.clear();
    m_sinput.clear();
    m_aflags.clear();
  }

  void EntryLoader::copy(const EntryLoader& b) {  //copy PRIVATE
    m_initialized=b.m_initialized;
    m_elflags=b.m_elflags;
    m_sinput=b.m_sinput;
    m_aflags=b.m_aflags;
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
    m_initialized=false;  //no point
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
    loadAFlags();
    loadInput(flags);
    m_initialized=false;  //no point
    return m_initialized;
  }
  
  bool EntryLoader::initialize(const string& sinput,ostream& oss) {
    xStream::initialize(oss);
    return initialize(sinput);
  }

  bool EntryLoader::initialize(const string& sinput,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(sinput);
  }

  bool EntryLoader::initialize(const string& sinput) {
    free();
    loadAFlags();
    loadInput(sinput);
    m_initialized=true;  //no point
    return m_initialized;
  }
  
  bool EntryLoader::initialize(const string& sinput,const aurostd::xoption& flags,ostream& oss) {
    xStream::initialize(oss);
    return initialize(sinput,flags);
  }

  bool EntryLoader::initialize(const string& sinput,const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize(sinput,flags);
  }

  bool EntryLoader::initialize(const string& sinput,const aurostd::xoption& flags) {
    free();
    loadAFlags();
    loadInput(sinput,flags);
    m_initialized=true;  //no point
    return m_initialized;
  }

  void EntryLoader::loadAFlags(){ //just a placeholder for necessary logger input
    if(XHOST.vflag_control.flag("DIRECTORY_CLEAN")){m_aflags.Directory=XHOST.vflag_control.getattachedscheme("DIRECTORY_CLEAN");} //CO20190402
    if(m_aflags.Directory.empty() || m_aflags.Directory=="./" || m_aflags.Directory=="."){m_aflags.Directory=aurostd::getPWD()+"/";} //".";  //CO20180220 //[CO20191112 - OBSOLETE]aurostd::execute2string(XHOST.command("pwd"))
  }

  void EntryLoader::loadInput(const string& sinput){m_sinput=sinput;}
  void EntryLoader::loadInput(const aurostd::xoption& flags){m_elflags=flags;}
  void EntryLoader::loadInput(const string& sinput,const aurostd::xoption& flags){loadInput(sinput);loadInput(flags);}

  void EntryLoader::retrieveOutput(string& soutput){
    soutput.clear();
    processInput();
    //string response = aflowlib::AFLUXCall(Summons);
  }
  void EntryLoader::retrieveOutput(vector<aflowlib::_aflowlib_entry>& entries){
    entries.clear();
  }
  void EntryLoader::retrieveOutput(vector<vector<aflowlib::_aflowlib_entry> >& entries){
    for(uint i=0;i<entries.size();i++){entries[i].clear();} entries.clear();
  }
  void EntryLoader::retrieveOutput(vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries){
    for(uint i=0;i<entries.size();i++){for(uint j=0;j<entries[i].size();j++){entries[i][j].clear();} entries[i].clear();} entries.clear();
  }

  void EntryLoader::processInput(){
    bool LDEBUG=(FALSE || _DEBUG_ENTRY_LOADER_ || XHOST.DEBUG);
    string soliloquy="EntryLoader()::processInput():";
    stringstream message;

    if(LDEBUG){cerr << soliloquy << " BEGIN" << endl;}

    message << "Processing input=\"" << m_sinput << "\"";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_MESSAGE_);

    vector<string> vkeys_aflux=getDataNames();

    if(LDEBUG){
      for(uint i=0;i<vkeys_aflux.size();i++){
        cerr << soliloquy << " vkeys_aflux[" << i << "]=" << vkeys_aflux[i] << endl;
      }
    }

    //look for an aflux summons first
    bool found_aflux=false;
    for(uint i=0;i<vkeys_aflux.size()&&found_aflux==false;i++){
      if(m_sinput.find(vkeys_aflux[i])!=string::npos){
        if(LDEBUG){cerr << soliloquy << " found aflux_key in m_sinput: " << vkeys_aflux[i] << endl;}
        //treating input as aflux summons
        found_aflux=true;
        m_elflags.push_attached("AFLUX::SUMMONS",m_sinput);
      }
    }

    //see if input is a bunch of elements
    if(found_aflux==false){
      vector<string> velements;
      aurostd::string2tokens(m_sinput,velements,",");
      for(uint i=0;i<velements.size();i++){

      }
    }

    if(found_aflux==false){}

  }

}

#endif  // _AFLOW_ENTRY_LOADER_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
