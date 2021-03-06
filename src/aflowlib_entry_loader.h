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
      bool m_input_processed;
      aurostd::xoption m_elflags;
      string m_sinput;
      _aflags m_aflags; //NOT an input, it's not required, just for directory manipulation
      vector<vector<vector<aflowlib::_aflowlib_entry> > > m_ventries; //organize as 3-layer, always easier to go from 3-layer to 1-layer
      
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
      void loadAFlags();  //dummy
      void loadInput(const string& sinput);
      void loadInput(const aurostd::xoption& flags);
      void loadInput(const string& sinput,const aurostd::xoption& flags);

      //getters
      void retrieveOutput(string& soutput);
      void retrieveOutput(vector<aflowlib::_aflowlib_entry>& ventries);
      void retrieveOutput(vector<vector<aflowlib::_aflowlib_entry> >& ventries);
      void retrieveOutput(vector<vector<vector<aflowlib::_aflowlib_entry> > >& ventries);
    private:
      //NECESSARY private CLASS METHODS - START
      void free();
      void copy(const EntryLoader& b);
      //NECESSARY END CLASS METHODS - END

      void processInput();
      void sanitizeAFLUXSummons();  //uint n=1; return first page, n from AFLUX paper
      uint AFLUXSummons2Entries();
      void setAFLUXSummons4ElementsString();
      void loadEntries();
  };
} // namespace aflowlib

#endif  // _AFLOW_ENTRY_LOADER_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
