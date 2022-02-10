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
#include <unordered_map>
#include <regex>

namespace aflowlib {
  class EntryLoader : public xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      EntryLoader(ostream& oss=cout);
      EntryLoader(ofstream& FileMESSAGE,ostream& oss=cout);
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
        FILESYSTEM_RAW,
        RESTAPI,
        RESTAPI_RAW,
        NONE,
        FAILED
      };

      // Settings
      string m_sqlite_file = DEFAULT_AFLOW_DB_FILE;
      string m_sqlite_alloy_file = "../testing/aflowlib_lookup.db";
      string m_aflux_server = "aflowlib.duke.edu";
      string m_aflux_path = "/API/aflux/v1.0/";
      //TODO use server & path from aflux rc
      std::map<std::string, std::string> m_aflux_directives {
          {"format", "aflow"},
          {"paging", "0"}};
      string m_restapi_server = "aflowlib.duke.edu";
      string m_restapi_path = "/AFLOWDATA/";
      string m_restapi_directives = "?format=text";
      string m_restapi_listing = "?aflowlib_entries";
      string m_restapi_collection = "_WEB/";
      string m_filesystem_outfile = DEFAULT_FILE_AFLOWLIB_ENTRY_OUT;
      string m_filesystem_path = "/common/";
      string m_filesystem_collection = "RAW/";

    bool m_entries_unique = true;


      _aflags m_aflags; //NOT an input, it's not required, just for directory manipulation
      //TODO understand usage


      // Attributes
      std::shared_ptr<aflowlib::AflowDB> m_sqlite_db_ptr;
      std::shared_ptr<aflowlib::AflowDB> m_sqlite_alloy_db_ptr;
      std::unordered_set<std::string> m_auid_list; // as sets are stored sorted enables faster find compared to a vector

      // Data views
      std::shared_ptr<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>> m_entries_flat;
      std::shared_ptr<std::unordered_map<short, std::unordered_map<
                      std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>> m_entries_layered_map;



      //initializers
      void initialize(ostream& oss);
      void initialize(ofstream& FileMESSAGE,ostream& oss);
      void initialize();

      //generic loader
      void loadAUID(string AUID);
      void loadAUID(const vector<std::string> &AUID);

      //TODO add icsd:?

      void loadAURL(string AURL);
      void loadAURL(const vector<std::string> &AURL);
      
      void loadAlloy(const string & alloy, bool recursive=true);
      void loadAlloy(const vector<std::string> &alloy, bool recursive=true);

      //specific loader
      void loadAFLUXQuery(const std::string & query);
      void loadAFLUXMatchbook(const std::map<string, string> & matchbook);

      void loadRestAPIQueries(const std::vector<std::string> & queries, bool full_url=false);

      void loadFiles(const std::vector<std::string> & files);

      void loadSqliteWhere(const std::string & where);

      void loadText(const std::vector<std::string> & raw_data_lines);

      //setter
      bool setSource(EntryLoader::Source new_source);

      //getter for views
      void getEntriesViewFlat(std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>> & result);
      void getEntriesViewTwoLayer(std::vector<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>> & result);
      void getEntriesViewThreeLayer(std::vector<std::vector<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>> & result);
      void getEntriesViewMap(std::unordered_map<short, std::unordered_map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>> & result);

      //getter that copy data into new format
      void getEntriesFlat(std::vector<aflowlib::_aflowlib_entry> & result);
      void getEntriesTwoLayer(std::vector<std::vector<aflowlib::_aflowlib_entry>> & result);
      void getEntriesThreeLayer(std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>> & result);



    private:
      //NECESSARY private CLASS METHODS - START
      void free();
      void copy(const EntryLoader& b);
      //NECESSARY END CLASS METHODS - END

      Source m_current_source = Source::NONE;

      // helper for FILESYSTEM_RAW
      std::regex m_re_elements{"(A[cglmrstu]|B[aehikr]?|C[adeflmnorsu]?|D[bsy]|E[rsu]|F[elmr]?|G[ade]|H[efgos]?|I[nr]?|Kr?|L[airuv]|M[dgnot]|N[abdeiop]?|Os?|P[abdmortu]?|R[abefghnu]|S[bcegimnr]?|T[abcehilm]|U(u[opst])?|V|W|Xe|Yb?|Z[nr])"};
      std::regex m_re_ppclean{"(_GW|_AE|_200eV)"};
      std::regex m_re_aurl2file{"(LIB\\d{1,})_(?:(?:RAW)|(?:LIB)|(?:WEB))\\/"};

      std::vector<string> m_icsd_symmetries = {"BCC/","BCT/","CUB/","FCC/","HEX/","MCL/","MCLC/",
                                               "ORC/","ORCC/","ORCF/","ORCI/","RHL/","TET/","TRI/"};

      std::string buildAFLUXQuery(const std::map<string, string> & matchbook);
      bool cleanAUID(string & AUID);
      bool cleanAURL(string & AURL);
      void selectSource();
      void getAlloyAUIDList(const std::vector<string> & alloy, std::vector<string> & auid_list);
      void loadAlloySearchFSR(const std::vector<string> & alloy_list, uint lib_max, bool recursive); // FILESYSTEM_RAW
      void loadAlloySearchRR(const std::vector<string> & alloy_list, uint lib_max, bool recursive); // RESTAPI_RAW
      void listRestAPI(string url, vector<string> & result);
      std::string extractAlloy(std::string name, char lib_type);

  };
} // namespace aflowlib

#endif  // _AFLOW_ENTRY_LOADER_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
