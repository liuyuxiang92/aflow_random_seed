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
      
      //initializers
      bool initialize(ostream& oss);
      bool initialize(ofstream& FileMESSAGE,ostream& oss);
      bool initialize();
      bool initialize(const aurostd::xoption& flags,ostream& oss);
      bool initialize(const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss);
      bool initialize(const aurostd::xoption& flags);
      bool initialize(const string& sinput,ostream& oss);
      bool initialize(const string& sinput,ofstream& FileMESSAGE,ostream& oss);
      bool initialize(const string& sinput);
      bool initialize(const string& sinput,const aurostd::xoption& flags,ostream& oss);
      bool initialize(const string& sinput,const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss);
      bool initialize(const string& sinput,const aurostd::xoption& flags);

      //loaders
      void loadInput(const string& sinput);
      void loadInput(const aurostd::xoption& flags);
      void loadInput(const string& sinput,const aurostd::xoption& flags);

      //getters
      void retrieveOutput(string& soutput);
      void retrieveOutput(vector<aflowlib::_aflowlib_entry>& entries);
      void retrieveOutput(vector<vector<aflowlib::_aflowlib_entry> >& entries);
      void retrieveOutput(vector<vector<vector<aflowlib::_aflowlib_entry> > >& entries);
      
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
