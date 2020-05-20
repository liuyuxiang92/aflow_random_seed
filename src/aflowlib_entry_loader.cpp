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
  }

  void EntryLoader::copy(const EntryLoader& b) {  //copy PRIVATE
    m_initialized=b.m_initialized;
    m_elflags=b.m_elflags;
    m_sinput=b.m_sinput;
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
    loadInput(sinput,flags);
    m_initialized=true;  //no point
    return m_initialized;
  }

  void EntryLoader::loadInput(const string& sinput){m_sinput=sinput;}
  void EntryLoader::loadInput(const aurostd::xoption& flags){m_elflags=flags;}
  void EntryLoader::loadInput(const string& sinput,const aurostd::xoption& flags){loadInput(sinput);loadInput(flags);}

  void EntryLoader::retrieveOutput(string& soutput){
    soutput.clear();
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

}

#endif  // _AFLOW_ENTRY_LOADER_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
