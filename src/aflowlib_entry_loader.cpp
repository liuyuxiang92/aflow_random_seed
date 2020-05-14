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
  EntryLoader::EntryLoader(ostream& oss) : xStream(),m_initialized(false) {initialize(oss);}
  EntryLoader::EntryLoader(ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(FileMESSAGE,oss);}
  EntryLoader::EntryLoader(const aurostd::xoption& flags,ostream& oss) : xStream(),m_initialized(false) {initialize(flags,oss);}
  EntryLoader::EntryLoader(const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(flags,FileMESSAGE,oss);}
  EntryLoader::EntryLoader(const string& sinput,ostream& oss) : xStream(),m_initialized(false) {initialize(sinput,oss);}
  EntryLoader::EntryLoader(const string& sinput,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(sinput,FileMESSAGE,oss);}
  EntryLoader::EntryLoader(const string& sinput,const aurostd::xoption& flags,ostream& oss) : xStream(),m_initialized(false) {initialize(sinput,flags,oss);}
  EntryLoader::EntryLoader(const string& sinput,const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss) : xStream(),m_initialized(false) {initialize(sinput,flags,FileMESSAGE,oss);}
  EntryLoader::EntryLoader(const EntryLoader& b) {copy(b);} // copy PUBLIC

  EntryLoader::~EntryLoader() {xStream::free();free();}
  
  const EntryLoader& EntryLoader::operator=(const EntryLoader& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  void EntryLoader::clear() {EntryLoader a;copy(a);}  //clear PUBLIC
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
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(*_p_FileMESSAGE,oss); 
    f_new_ofstream=true;  //override
    return m_initialized;
  }
  
  bool EntryLoader::initialize(ofstream& FileMESSAGE,ostream& oss) {
    free();
    setOFStream(FileMESSAGE); f_new_ofstream=false;
    setOSS(oss);
    m_initialized=false;  //no point
    return m_initialized;
  }
  
  bool EntryLoader::initialize(const aurostd::xoption& flags,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(flags,*_p_FileMESSAGE,oss); 
    f_new_ofstream=true;  //override
    return m_initialized;
  }
  
  bool EntryLoader::initialize(const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss) {
    free();
    setOFStream(FileMESSAGE); f_new_ofstream=false;
    setOSS(oss);
    m_elflags=flags;
    m_initialized=false;  //no point
    return m_initialized;
  }
  
  bool EntryLoader::initialize(const string& sinput,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(sinput,*_p_FileMESSAGE,oss); 
    f_new_ofstream=true;  //override
    return m_initialized;
  }
  
  bool EntryLoader::initialize(const string& sinput,ofstream& FileMESSAGE,ostream& oss) {
    free();
    setOFStream(FileMESSAGE); f_new_ofstream=false;
    setOSS(oss);
    m_sinput=sinput;
    m_initialized=true;  //no point
    return m_initialized;
  }
  
  bool EntryLoader::initialize(const string& sinput,const aurostd::xoption& flags,ostream& oss) {
    xStream::free();
    ofstream* _p_FileMESSAGE=new ofstream();f_new_ofstream=true;
    initialize(sinput,flags,*_p_FileMESSAGE,oss); 
    f_new_ofstream=true;  //override
    return m_initialized;
  }
  
  bool EntryLoader::initialize(const string& sinput,const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss) {
    free();
    setOFStream(FileMESSAGE); f_new_ofstream=false;
    setOSS(oss);
    m_sinput=sinput;
    m_elflags=flags;
    m_initialized=true;  //no point
    return m_initialized;
  }
}

#endif  // _AFLOW_ENTRY_LOADER_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
