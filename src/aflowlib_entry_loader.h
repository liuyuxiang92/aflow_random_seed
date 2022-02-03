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

#include <unordered_set>

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

      //Defaults
      enum class Source {
        SQLITE,
        AFLUX,
        FILESYSTEM,
        RESTAPI,
        NONE,
        FAILED
      };
      vector<Source> m_source_order = {Source::SQLITE, Source::AFLUX,
                                       Source::FILESYSTEM, Source::RESTAPI};

      string m_sqlite_file = DEFAULT_AFLOW_DB_FILE;
      string m_aflux_server = "aflowlib.duke.edu";
      string m_aflux_path = "/API/aflux/v1.0/";
      std::map<std::string, std::string> m_aflux_directives {
          {"format", "aflow"},
          {"paging", "0"}};

      //Settings
      bool m_entries_unique = true;
      Source m_current_source = Source::NONE;

      //attributes
      bool m_initialized;
      bool m_input_processed;
      std::shared_ptr<aflowlib::AflowDB> m_sqlite_db_ptr;

      aurostd::xoption m_elflags;
      string m_sinput;
      _aflags m_aflags; //NOT an input, it's not required, just for directory manipulation
      vector<vector<vector<aflowlib::_aflowlib_entry> > > m_ventries; //organize as 3-layer, always easier to go from 3-layer to 1-layer
      vector<aflowlib::_aflowlib_entry> m_lib_entries;
      std::unordered_set<std::string> m_auid_list; // as sets are stored sorted enables faster find compared to a vector

      //initializers
      bool initialize(ostream& oss);
      bool initialize(ofstream& FileMESSAGE,ostream& oss);
      bool initialize();
      bool initialize(const aurostd::xoption& flags,ostream& oss);
      bool initialize(const aurostd::xoption& flags,ofstream& FileMESSAGE,ostream& oss);
      bool initialize(const aurostd::xoption& flags);

      //generic loader
      void loadAUID(string AUID);
      void loadAUID(const vector<std::string> &AUID);

      void loadAURL(const string& AURL);
      void loadAURL(const vector<std::string> &AURL);
      
      void loadAlloy(const string & alloy, bool recursive=true);
      void loadAlloy(const vector<std::string> &alloy, bool recursive=true);

      //specific loader
      void loadAFLUXQuery(const std::string & query);
      void loadAFLUXMatchbook(const std::map<string, string> & matchbook);

      void loadSqliteWhere(const std::string & where);

      void loadAFlags();

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

      string buildAFLUXQuery(const std::map<string, string> & matchbook);
      void selectSource();

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
