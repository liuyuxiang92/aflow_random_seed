// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
// Corey Oses
// corey.oses@duke.edu

#ifndef _AFLOWLIB_ENTRY_LOADER_H_
#define _AFLOWLIB_ENTRY_LOADER_H_

namespace aflowlib {
  class EntryLoader : public xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      EntryLoader(ostream& oss=cout);
      EntryLoader(ofstream& FileMESSAGE,ostream& oss=cout);
      EntryLoader(const aurostd::xoption& flags,ostream& oss=cout);
      EntryLoader(const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss=cout);
      EntryLoader(const string& sinput,ostream& oss=cout);
      EntryLoader(const string& sinput,ofstream& FileMESSAGE,ostream& oss=cout);
      EntryLoader(const string& sinput,const aurostd::xoption& flags,ostream& oss=cout);
      EntryLoader(const string& sinput,const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss=cout);
      EntryLoader(const EntryLoader& b);
      //constructors - STOP
      ~EntryLoader();
      const EntryLoader& operator=(const EntryLoader& other);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - STOP
      
      //attributes
      bool m_initialized;
      aurostd::xoption m_elflags;
      string m_sinput;
      
      //initialization methods
      bool initialize(ostream& oss=cout);
      bool initialize(ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const aurostd::xoption& flags,ostream& oss=cout);
      bool initialize(const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const string& sinput,ostream& oss=cout);
      bool initialize(const string& sinput,ofstream& FileMESSAGE,ostream& oss=cout);
      bool initialize(const string& sinput,const aurostd::xoption& flags,ostream& oss=cout);
      bool initialize(const string& sinput,const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss=cout);

      //getters
      //
    private:
      //NECESSARY private CLASS METHODS - START
      void free();
      void copy(const EntryLoader& b);
      //NECESSARY END CLASS METHODS - END
  };
} // namespace aflowlib

#endif  // _AFLOW_ENTRY_LOADER_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
