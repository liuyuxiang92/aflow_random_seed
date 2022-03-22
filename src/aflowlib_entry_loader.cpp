// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
// Developed in fruitful collaboration by Corey Oses and Hagen Eckert.
// corey.oses@duke.edu | hagen.eckert@duke.edu


#ifndef _AFLOWLIB_ENTRY_LOADER_CPP_
#define _AFLOWLIB_ENTRY_LOADER_CPP_


#include "aflow.h"
#include "aflowlib_entry_loader.h"

#define _DEBUG_ENTRY_LOADER_ false


namespace aflowlib {

  /// @class EntryLoader
  /// @brief unified interface to gather AFLOW lib entries from different sources
  ///
  /// Basic usage
  /// @code
  /// std::vector <std::shared_ptr<aflowlib::_aflowlib_entry>> results;
  /// { // ensures that all resources are released as soon as possible
  ///   aflowlib::EntryLoader el;
  ///   el.m_out_debug = true; // change the behavior of the EntryLoader class
  ///   el.loadAlloy("NiCa"); // the load functions can be called multiple times
  ///   el.loadAUID(aflow:7dd846bc04c764e8); // no duplicates will be stored
  ///   el.getEntriesViewFlat(results);
  /// }
  /// for (auto & entry : results) std::cout << entry->auid << std::endl;
  /// @endcode

  // class constructor
  EntryLoader::EntryLoader(std::ostream &oss) : xStream(oss) { init(); }

  EntryLoader::EntryLoader(std::ofstream &FileMESSAGE, std::ostream &oss) : xStream(FileMESSAGE, oss) { init(); }

  /// @brief create a copy at initialization
  /// @note DB connections are reused and not recreated
  /// @note a new view on the data is created, but the same lib entries are reused
  EntryLoader::EntryLoader(const EntryLoader &b) : xStream(*b.getOFStream(), *b.getOSS()) { copy(b); } // copy PUBLIC

  // class de-constructor
  EntryLoader::~EntryLoader() { xStream::free(); }

  /// @brief create a copy through assigning
  /// @note DB connections are reused and not recreated
  /// @note a new view on the data is created, but the same lib entries are reused
  const EntryLoader &EntryLoader::operator=(const EntryLoader &other) {
    if (this != &other) { copy(other); }
    return *this;
  }

  /// @brief initialize the class (privat)
  /// create shared pointer used to store the data views and read default values from flags
  /// @TODO use flags to overwrite defaults
  void EntryLoader::init() {
    m_out_debug = (m_out_debug || XHOST.DEBUG || _DEBUG_ENTRY_LOADER_);
    m_entries_flat = std::make_shared<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>();
    m_entries_layered_map = std::make_shared<std::map<short,
                                             std::map<std::string,
                                             std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>>();
  }

  /// @brief reset a EntryLoader object
  void EntryLoader::clear() { *this = {}; }  // calling the constructor

  /// @brief create a copy (privat)
  void EntryLoader::copy(const EntryLoader &b) {  //copy PRIVATE
    m_out_debug = b.m_out_debug;
    m_out_silent = b. m_out_silent;
    m_xstructure_relaxed = b.m_xstructure_relaxed;
    m_xstructure_original = b.m_xstructure_original;
    m_xstructure_original_file_name = b.m_xstructure_original_file_name;
    m_xstructure_final_file_name = b.m_xstructure_final_file_name;
    m_sqlite_file = b.m_sqlite_file;
    m_sqlite_alloy_file = b.m_sqlite_alloy_file;
    m_sqlite_collection = b.m_sqlite_collection;
    m_aflux_server = b.m_aflux_server;
    m_aflux_path = b.m_aflux_path;
    m_aflux_collection = b.m_aflux_collection;
    m_aflux_directives = b.m_aflux_directives;
    m_restapi_server = b.m_restapi_server;
    m_restapi_path = b.m_restapi_path;
    m_restapi_directives = b.m_restapi_directives;
    m_restapi_listing_dirs = b.m_restapi_listing_dirs;
    m_restapi_listing_files = b.m_restapi_listing_files;
    m_restapi_collection = b.m_restapi_collection;
    m_filesystem_outfile = b.m_filesystem_outfile;
    m_filesystem_path = b.m_filesystem_path;
    m_filesystem_collection = b.m_filesystem_collection;
    m_current_source = b.m_current_source;
    m_filesystem_available = b.m_filesystem_available;
    m_auid_list = b.m_auid_list;
    aurostd::StringstreamClean(m_logger_message);

    m_sqlite_db_ptr = b.m_sqlite_db_ptr;
    m_sqlite_alloy_db_ptr = b.m_sqlite_db_ptr;

    m_entries_flat = std::make_shared < std::vector < std::shared_ptr < aflowlib::_aflowlib_entry>>>(*b.m_entries_flat);
    m_entries_layered_map = std::make_shared < std::map < short,
        std::map < std::string,
        std::vector < std::shared_ptr < aflowlib::_aflowlib_entry >>>>>(*b.m_entries_layered_map);
  }

  /// @brief load a single entry from a AUID
  /// @param AUID AFLOW unique ID
  /// @note following AUIDs would be equivalent:
  /// @note `aflow:d912e209c81aeb94`
  /// @note `auid:d912e209c81aeb94`
  /// @note `d912e209c81aeb94`
  void EntryLoader::loadAUID(std::string AUID) {
    selectSource();

    m_logger_message << "Try loading: " << AUID;
    outInfo(__func__);
    if (!cleanAUID(AUID)) {
      m_logger_message << "AUID cleaning failed";
      outError(__func__, __LINE__);
      return;
    }

    switch (m_current_source) {

      case Source::SQLITE: {
        std::string where = "auid='\"" + AUID + "\"'";
        loadSqliteWhere(where);
        break;
      }

      case Source::AFLUX: {
        std::map <std::string, std::string> matchbook{{"*",    ""},
                                                      {"auid", "'" + AUID + "'"}};
        loadAFLUXMatchbook(matchbook);
        break;
      }

      case Source::RESTAPI:
      case Source::RESTAPI_RAW: {
        loadAUID(std::vector<std::string>({AUID}));
        break;
      }

      case Source::FILESYSTEM:
      case Source::FILESYSTEM_RAW: {
        loadAUID(std::vector<std::string>({AUID}));
        break;
      }

      default:
        break;
    }
  }

  /// @brief load multiple entries from a vector of AUIDs
  /// @param AUID list of AURLs
  void EntryLoader::loadAUID(const std::vector <std::string> &AUID) {
    selectSource();
    std::vector <std::string> clean_AUID;

    m_logger_message << "Try loading " << AUID.size() <<" AUIDs";
    outInfo(__func__);

    for (std::string AUID_single: AUID) {
      if (cleanAUID(AUID_single)) clean_AUID.push_back(AUID_single);
    }

    if (clean_AUID.size()!=AUID.size()) {
      m_logger_message << "cleaning of " << AUID.size() - clean_AUID.size()  << " AUIDs failed";
      outError(__func__, __LINE__);
    }

    switch (m_current_source) {

      case Source::AFLUX: {
        std::string AUID_combined = aurostd::joinWDelimiter(clean_AUID, "':'");
        AUID_combined = "'" + AUID_combined + "'";
        loadAFLUXMatchbook({{"*",    ""},
                            {"auid", AUID_combined}});
        break;
      }

      case Source::SQLITE: {
        std::string where = aurostd::joinWDelimiter(clean_AUID, "\"','\"");
        where = "auid IN ('\"" + where + "\"')";
        loadSqliteWhere(where);
        break;
      }

      case Source::RESTAPI:
      case Source::RESTAPI_RAW: {
        std::vector <std::string> queries;
        for (std::string AUID_single: clean_AUID) {
          std::string rest_query = "AUID/aflow:";
          for (uint part = 0; part < 8; part++) rest_query += AUID_single.substr(6 + part * 2, 2) + "/";
          rest_query += "WEB/";
          queries.push_back(rest_query);
        }
        loadRestAPIQueries(queries);
        break;
      }

      case Source::FILESYSTEM:
      case Source::FILESYSTEM_RAW: {
        std::vector <std::string> files;
        for (std::string AUID_single: clean_AUID) {
          std::string file_path = m_filesystem_path + "AUID/aflow:";
          for (uint part = 0; part < 8; part++) file_path += AUID_single.substr(6 + part * 2, 2) + "/";
          file_path += m_filesystem_collection + "/" + m_filesystem_outfile;
          files.push_back(file_path);
        }
        loadFiles(files);
        break;
      }

      default:
        break;
    }
  }

  /// @brief load a single entry from a AURL
  /// @param AURL AFLOW unique URL
  /// @note following AURLs would be equivalent:
  /// @note `aflowlib.duke.edu:AFLOWDATA/LIB2_RAW/Ca_svCu_pv/138`
  /// @note `AFLOWDATA/LIB2_RAW/Ca_svCu_pv/138`
  /// @note `LIB2_RAW/Ca_svCu_pv/138`
  /// @note `WEB`, `RAW`, and `LIB` will be replaced by the selected collection
  ///       (#m_sqlite_collection, #m_aflux_collection, #m_restapi_collection, #m_filesystem_collection)
  void EntryLoader::loadAURL(std::string AURL) {
    selectSource();

    if (!cleanAURL(AURL)) return;

    m_logger_message << "Try loading " << AURL;
    outInfo(__func__);
    if (!cleanAURL(AURL)) {
      m_logger_message << "AURL cleaning failed";
      outError(__func__, __LINE__);
      return;
    }

    switch (m_current_source) {

      case Source::SQLITE: {
        std::string where = "aurl='\"" + AURL + "\"'";
        where = std::regex_replace(where, m_re_aurl2file, "$1_" + m_sqlite_collection + "/");
        loadSqliteWhere(where);
        break;
      }

      case Source::AFLUX: {
        AURL = std::regex_replace(AURL, m_re_aurl2file, "$1_" + m_aflux_collection + "/");
        std::map <std::string, std::string> matchbook{{"*",    ""},
                                                      {"aurl", "'" + AURL + "'"}};
        loadAFLUXMatchbook(matchbook);
        break;
      }

      case Source::RESTAPI:
      case Source::RESTAPI_RAW: {
        loadAURL(std::vector<std::string>({AURL}));
        break;
      }

      case Source::FILESYSTEM:
      case Source::FILESYSTEM_RAW: {
        loadAURL(std::vector<std::string>({AURL}));
        break;
      }

      default:
        break;
    }
  }

  /// @brief load multiple entries from a vector of AURLs
  /// @param AURL list of AURLs
  void EntryLoader::loadAURL(const std::vector <std::string> &AURL) {
    selectSource();

    m_logger_message << "Try loading " << AURL.size() <<" AURLs";
    outInfo(__func__);

    std::vector <std::string> clean_AURL;
    for (std::string AURL_single: AURL) {
      if (cleanAURL(AURL_single)) clean_AURL.push_back(AURL_single);
    }

    if (clean_AURL.size()!=AURL.size()) {
      m_logger_message << "cleaning of " << AURL.size() - clean_AURL.size()  << " AURLs failed";
      outError(__func__, __LINE__);
    }

    switch (m_current_source) {

      case Source::AFLUX: {
        std::string AURL_combined = aurostd::joinWDelimiter(clean_AURL, "':'");
        AURL_combined = "'" + AURL_combined + "'";
        AURL_combined = std::regex_replace(AURL_combined, m_re_aurl2file, "$1_" + m_aflux_collection + "/");
        loadAFLUXMatchbook({{"*",    ""},
                            {"aurl", AURL_combined}});
        break;
      }

      case Source::SQLITE: {
        std::string where = aurostd::joinWDelimiter(clean_AURL, "\"','\"");
        where = std::regex_replace(where, m_re_aurl2file, "$1_" + m_sqlite_collection + "/");
        where = "aurl IN ('\"" + where + "\"')";
        loadSqliteWhere(where);
        break;
      }

      case Source::RESTAPI:
      case Source::RESTAPI_RAW: {
        std::vector <std::string> queries;
        for (std::string AURL_single: clean_AURL) {
          std::string rest_query = std::regex_replace(AURL_single.substr(28), m_re_aurl2file, "$1_" + m_sqlite_collection + "/") + "/";
          queries.push_back(rest_query);
        }
        loadRestAPIQueries(queries);
        break;
      }

      case Source::FILESYSTEM:
      case Source::FILESYSTEM_RAW: {
        std::vector <std::string> files;
        for (std::string AURL_single: clean_AURL) {
          std::string file_path = m_filesystem_path + AURL_single.substr(28) + "/" + m_filesystem_outfile;
          file_path = std::regex_replace(file_path, m_re_aurl2file, "$1/" + m_filesystem_collection + "/");
          files.push_back(file_path);
        }
        loadFiles(files);
        break;
      }

      default:
        break;
    }
  }

  /// @brief load all entries that are available for a alloy set
  /// @param alloy alloy descriptor
  /// @param recursive include sub alloys
  /// @note an alloy can be given as comma separated list `"Cu, Au, Fe"`
  ///       or as combined string `"CuAuFe"` (elements don't need to be sorted)
  /// @note `recursive=true` (default) will generate sub-alloys so `CuAuFe`
  ///       will load `AuCuFe`, `AuCu`, `AuFe`, `CuFe`, `Au`, `Cu`, and `Fe`
  void EntryLoader::loadAlloy(const std::string &alloy, bool recursive) {
    std::vector <std::string> alloy_elements;
    if (alloy.find(", ") != std::string::npos) aurostd::string2tokens(alloy, alloy_elements, ", ");
    else if (alloy.find(",") != std::string::npos) aurostd::string2tokens(alloy, alloy_elements, ",");
    else alloy_elements = aurostd::getElements(alloy);
    loadAlloy(alloy_elements, recursive);
  }

  /// @brief load all entries that are available for a set alloy
  /// @param alloy split alloy elements
  /// @param recursive include sub alloys
  void EntryLoader::loadAlloy(const std::vector <std::string> &alloy, bool recursive) {
    selectSource();
    std::vector <std::string> alloy_clean;

    // check that all parts of alloy are actually elements
    for (std::string element: alloy) {
      if (xelement::xelement::isElement(element) == 0) {
        m_logger_message << element << " is not an element";
        outError(__func__,__LINE__);
      } else alloy_clean.emplace_back(element);
    }

    // build list of all needed sub-alloys
    std::vector <std::vector<std::string>> alloy_sub_list;
    if (recursive) {
      aurostd::xcombos xc;
      for (size_t i = 1; i <= alloy_clean.size(); i++) {
        xc.reset(alloy_clean.size(), i);
        std::vector<int> choice;
        while (xc.increment()) {
          std::vector <std::string> alloy_tmp;
          choice = xc.getCombo();
          for (size_t idx = 0; idx < choice.size(); ++idx) {
            if (choice[idx]) alloy_tmp.push_back(alloy_clean[idx]);
          }
          alloy_sub_list.push_back(alloy_tmp);
        }
      }
    } else {
      alloy_sub_list.push_back(alloy_clean);
    }

    // build sorted alloy strings
    std::vector <std::string> final_alloy_list;
    for (std::vector <std::string> &a: alloy_sub_list) {
      std::string alloy_tmp = "";
      std::sort(a.begin(), a.end());
      for (std::string e: a) alloy_tmp += e;
      final_alloy_list.push_back(alloy_tmp);
    }
    m_logger_message << "Try loading " << final_alloy_list.size() << " systems: ";
    for (const std::string & alloy : final_alloy_list){
      m_logger_message << "\"" << alloy << "\" ";
    }
    outInfo(__func__);

    // load the alloys
    switch (m_current_source) {
      case Source::SQLITE: {
        std::string where = aurostd::joinWDelimiter(final_alloy_list, "','");
        where = "alloy IN ('" + where + "')";
        loadSqliteWhere(where);
        break;
      }

      case Source::AFLUX: {
        std::string alloy_match = aurostd::joinWDelimiter(final_alloy_list, ":");
        std::map <std::string, std::string> matchbook{{"*",     ""},
                                                      {"alloy", alloy_match}};
        loadAFLUXMatchbook(matchbook);
        break;
      }

      case Source::RESTAPI: {
        std::vector <std::string> auid_list;
        getAlloyAUIDList(final_alloy_list, auid_list);
        loadAUID(auid_list);
        break;
      }

      case Source::RESTAPI_RAW: {
        loadAlloySearchRR(final_alloy_list, alloy_clean.size());
        break;
      }

      case Source::FILESYSTEM: {
        std::vector <std::string> auid_list;
        getAlloyAUIDList(final_alloy_list, auid_list);
        loadAUID(auid_list);
        break;
      }

      case Source::FILESYSTEM_RAW: {
        loadAlloySearchFSR(final_alloy_list, alloy_clean.size());
        break;
      }

      default:
        break;
    }

  }

  /// @brief load entries based on a custom query against the internal AFLUX SQLITE DB
  /// @param where part of the SQLITE query after a `WHERE` keyword
  /// @note use setSource() to change the current source to Source::SQLITE to ensure that the database is connected
  void EntryLoader::loadSqliteWhere(const std::string &where) {
    std::vector <std::string> raw_lines = getRawSqliteWhere(where);
    m_logger_message << raw_lines.size() << " entries found in DB for " << where;
    outDebug(__func__);

    size_t start_size = m_entries_flat->size();
    loadText(raw_lines);

    m_logger_message << "Loaded " << m_entries_flat->size() - start_size << " new entries";
    outInfo(__func__);

  }

  /// @brief load entries from a costum AFLUX query
  /// @param query AFLUX query
  /// @note #m_aflux_server and #m_aflux_path will be added
  void EntryLoader::loadAFLUXQuery(const std::string &query) {
    size_t start_size = m_entries_flat->size();
    loadText(getRawAFLUXQuery(query));
    m_logger_message << "Loaded " << m_entries_flat->size() - start_size << " new entries";
    outInfo(__func__);
  }

  /// @brief load entries from an AFLUX matchbook
  /// @param matchbook map with keyword and modifier
  /// @note a pair creates `<keyword>(<modifier>)`
  /// @note if modifier is empty just the keyword is used
  void EntryLoader::loadAFLUXMatchbook(const std::map <std::string, std::string> &matchbook) {
    loadAFLUXQuery(buildAFLUXQuery(matchbook));
  }

  /// @brief load entries from a list of REST API queries
  /// @param queries vector of queries
  /// @param full_url if `true` don't add #m_restapi_server, #m_restapi_path, and #m_restapi_directives
  void EntryLoader::loadRestAPIQueries(const std::vector <std::string> &queries, bool full_url) {
    size_t start_size = m_entries_flat->size();
    size_t done_downloads = 0;
    for (std::string query: queries) {
      loadText({getRawRestAPIQuery(query, full_url)});
      done_downloads++;
      if (done_downloads%100 == 0){
        m_logger_message << "Loaded " << done_downloads << " of " << queries.size();
        outInfo(__func__);
      }
    }
    m_logger_message << "Loaded " << m_entries_flat->size() - start_size << " new entries";
    outInfo(__func__);
  }

  /// @brief load entries from a list of file paths
  /// @param files list of file paths
  /// @note doesn't add #m_filesystem_path or #m_filesystem_collection
  void EntryLoader::loadFiles(const std::vector <std::string> &files) {
    size_t start_size = m_entries_flat->size();
    for (std::string file_path: files) {
       std::string file_content;
       if (aurostd::file2string(file_path, file_content) > 0) {
         loadText({file_content});
       }
    }
    m_logger_message << "Loaded " << m_entries_flat->size() - start_size << " new entries";
    outInfo(__func__);
  }

  /// @brief load entries from a list of strings
  /// @param raw_data_lines list of data strings
  /// @note each entry in raw_data_lines should correspond to one AFLOW lib entry
  void EntryLoader::loadText(const std::vector <std::string> &raw_data_lines) {
    for (std::string line: raw_data_lines) {
      std::shared_ptr <aflowlib::_aflowlib_entry> entry = std::make_shared<aflowlib::_aflowlib_entry>();
      entry->Load(line, *p_oss);
      if (!entry->auid.empty() && (std::find(m_auid_list.begin(),m_auid_list.end(), entry->auid) == m_auid_list.end())) {
        m_entries_flat->push_back(entry);
        (*m_entries_layered_map)[entry->nspecies][entry->species_pp].push_back(entry);
        m_auid_list.emplace_back(entry->auid);
        if (m_xstructure_original) addXstructure(*entry, true);
        if (m_xstructure_relaxed) addXstructure(*entry);
      }
    }
  }

  /// @brief returns the active data source
  EntryLoader::Source EntryLoader::getSource() const{
    return m_current_source;
  }

  /// @brief change the currently used source and prepares them
  /// @param new_source source to change to
  /// @param silent silent all logging when called from selectSource() (default false)
  /// @note if the public alloy SQLITE DB is not found Source::FILESYSTEM and Source::RESTAPI
  ///       are mapped to Source::FILESYSTEM_RAW and Source::RESTAPI_RAW respectively
  bool EntryLoader::setSource(EntryLoader::Source new_source) {
    if (new_source == m_current_source) return true;
    m_current_source = Source::FAILED;
    m_filesystem_available = aurostd::IsDirectory(m_filesystem_path + "AUID/");

    switch (new_source) {
      case Source::SQLITE: {
        if (aurostd::FileExist(m_sqlite_file)) {
          if (m_sqlite_db_ptr == nullptr) {
            m_sqlite_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_file);
          }
          m_current_source = Source::SQLITE;
          return true;
        } else {
          m_logger_message << "Internal AFLUX SQLITE DB not found at " << m_sqlite_file;
          outError(__func__,__LINE__);
        }
        break;
      }

      case Source::AFLUX: {
        if ("AFLUXtest" == aurostd::httpGet(m_aflux_server+"/test/?echo=AFLUXtest")) {
          m_current_source = Source::AFLUX;
          return true;
        } else {
          m_logger_message << "AFLUX API could not be reached at " << m_aflux_server;
          outError(__func__,__LINE__);
        }
        break;
      }

      case Source::FILESYSTEM: {
        if (m_filesystem_available) {
          if (aurostd::FileExist(m_sqlite_alloy_file)) { //
            if (m_sqlite_alloy_db_ptr == nullptr) {
              m_sqlite_alloy_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_alloy_file);
            }
            m_current_source = Source::FILESYSTEM;
            return true;
          } else {
            m_logger_message << "Could not find public alloy SQLITE DB at " << m_sqlite_alloy_file << "; switching to Source::FILESYSTEM_RAW";
            outInfo(__func__);
            return setSource(Source::FILESYSTEM_RAW);
          }
        } else {
          m_logger_message << "Could not find AFLOW entries in the filesystem at " << m_filesystem_path;
          outError(__func__,__LINE__);
        }
      }

      case Source::FILESYSTEM_RAW: {
        if (m_filesystem_available) { //
          if (aurostd::FileExist(m_sqlite_alloy_file)) { // use the alloy DB to optimize the RAW results
            if (m_sqlite_alloy_db_ptr == nullptr) {
              m_sqlite_alloy_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_alloy_file);
            }
          }
          m_current_source = Source::FILESYSTEM_RAW;
          return true;
        } else {
          m_logger_message << "Could not find AFLOW entries in the filesystem at " << m_filesystem_path;
          outError(__func__,__LINE__);
        }
        break;
      }

      case Source::RESTAPI: {
        if (200 == aurostd::httpGetStatus(m_restapi_server+m_restapi_path+"ICSD_WEB/")) {
          if (aurostd::FileExist(m_sqlite_alloy_file)) {
            if (m_sqlite_alloy_db_ptr == nullptr) {
              m_sqlite_alloy_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_alloy_file);
            }
            m_current_source = Source::RESTAPI;
            return true;
          } else {
            m_logger_message << "Could not find public alloy SQLITE DB at " << m_sqlite_alloy_file << "; switching to Source::RESTAPI_RAW";
            outInfo(__func__);
            return setSource(Source::RESTAPI_RAW);
          }
        } else {
          m_logger_message << "AFLOW REST API could not be reached at " << m_aflux_server+m_restapi_path;
          outError(__func__,__LINE__);
        }
        break;
      }

      case Source::RESTAPI_RAW: {
        if (200 == aurostd::httpGetStatus(m_restapi_server+m_restapi_path+"ICSD_WEB/")) {
          if (aurostd::FileExist(m_sqlite_alloy_file)) {
            if (m_sqlite_alloy_db_ptr == nullptr) {
              m_sqlite_alloy_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_alloy_file);
            }
          }
          m_current_source = Source::RESTAPI_RAW;
          return true;
        }
        else {
          m_logger_message << "AFLOW REST API could not be reached at " << m_aflux_server+m_restapi_path;
          outError(__func__,__LINE__);
        }
        break;
      }

      case Source::NONE: {
        m_current_source = Source::NONE;
        return true;
      }

      case Source::FAILED: {
        m_current_source = Source::FAILED;
        return true;
      }

    }
    return false;
  }

  /// @brief query the internal AFLUX SQLITE DB
  /// @param where part of the SQLITE query after a `WHERE` keyword
  /// @return list of result query strings
  std::vector <std::string> EntryLoader::getRawSqliteWhere(const std::string &where) const{
    if (m_sqlite_db_ptr != nullptr) return m_sqlite_db_ptr->getEntrySet(where, aflow_ft);
    else return {"",};
  }

  /// @brief query AFLUX using a matchbook
  /// @param matchbook map with keyword and modifier
  /// @return list of query result strings
  std::vector <std::string> EntryLoader::getRawAFLUXMatchbook(const std::map <std::string, std::string> &matchbook){
    return getRawAFLUXQuery(buildAFLUXQuery(matchbook));
  }

  /// @brief query AFLUX using a raw query
  /// @param matchbook map with keyword and modifier
  /// @return list of query result strings
  /// @note #m_aflux_server and #m_aflux_path will be added
  std::vector <std::string> EntryLoader::getRawAFLUXQuery(const std::string &query){
    std::string output = "";
    std::vector <std::string> raw_lines;
    if (200 == aurostd::httpGetStatus(m_aflux_server, m_aflux_path, query, output)){
      aurostd::string2vectorstring(output, raw_lines);
    } else {
      m_logger_message << "Failed to get AFLUX query ";
      m_logger_message << "(" << m_aflux_server << " | " << m_aflux_path << " | " << query << ")";
      outError(__func__, __LINE__);
    }
    return raw_lines;
  }

  /// @brief query the REST API
  /// @param query REST API query
  /// @param full_url if `true` don't add #m_restapi_server, #m_restapi_path, and #m_restapi_directives
  /// @return query result string
  std::string EntryLoader::getRawRestAPIQuery(const std::string &query, bool full_url) {
    std::string output = "";
    std::string url;
    if (full_url) url = query;
    else url = m_restapi_server + m_restapi_path + query + m_restapi_directives;
    if (200 == aurostd::httpGetStatus(url, output)){
      return output;
    } else {
      m_logger_message << "Failed to get REST API query " << "(" << url << ")";
      outError(__func__, __LINE__);
      return "";
    }
  }

  /// @brief add an xstructure to an AFLOW lib entry
  /// @param entry AFLOW lib entry
  /// @param orig if `true` load the original structure otherwise the final structure (default)
  /// @note use this function only if you want the xstructure for a subset of entries;
  ///       if you want xstructures for all entries set #m_xstructure_relaxed or #m_xstructure_original to `true`
  /// @note adds the structure to entry.vstr
  void EntryLoader::addXstructure(aflowlib::_aflowlib_entry &entry, bool orig) {
    if (orig) {
      if (entry.catalog == "ICSD") return; // no original structures available for ICSD entries
      xstructure new_structure;
      if (loadXstructureAflowIn(entry, new_structure)) { // load from aflow.in
        entry.vstr.push_back(new_structure);
      } else { // if no valid entry in aflow.in try to load files from m_xstructure_original_file_name list
        if (loadXstructureFile(entry, new_structure, true)) {
          entry.vstr.push_back(new_structure);
        } else {
          m_logger_message << "Failed to add original structure to " << entry.auid << " (" << entry.aurl << ")";
          outError(__func__, __LINE__);
        }
      }
    } else {
      xstructure new_structure;
      if (loadXstructureLibEntry(entry, new_structure)) { // load directly from the entry (AFLUX, SQLITE)
        entry.vstr.push_back(new_structure);
      } else { // if positions are not present in aflowlib entry try to load files from m_xstructure_final_file_name list
        if (loadXstructureFile(entry, new_structure, false)) {
          entry.vstr.push_back(new_structure);
        }
        else {
          m_logger_message << "Failed to add relaxed structure to " << entry.auid << " (" << entry.aurl << ")";
          outError(__func__, __LINE__);
        }
      }
    }
  }

  /// @brief load a structure from a structure file
  /// @param entry AFLOW lib entry
  /// @param new_structure save-to structure
  /// @param orig if `true` load the original structure otherwise the final structure (default)
  /// @return xstructure
  /// @note does not add the structure to entry.vstr
  bool EntryLoader::loadXstructureFile(const aflowlib::_aflowlib_entry &entry, xstructure &new_structure, bool orig) {
    std::string base_url = m_restapi_server + m_restapi_path + entry.aurl.substr(28) + "/";
    std::string base_folder = m_filesystem_path + entry.aurl.substr(28) + "/";
    base_folder = std::regex_replace(base_folder, m_re_aurl2file, "$1/" + m_filesystem_collection + "/");
    std::string poscar;
    if (entry.catalog =="LIB0" && !m_filesystem_available && entry.aurl.substr(entry.aurl.size() - 2)=="/0"){
      return false; // no entries in RESTAPI for WEB0 /0
    }

    std::vector <std::string> possible_files;
    if (orig) possible_files = m_xstructure_original_file_name;
    else possible_files = m_xstructure_final_file_name;

    std::vector <std::string> available_files;
    if (m_filesystem_available) {
      aurostd::DirectoryLS(base_folder, available_files);
    } else {
      listRestAPI(base_url, available_files, false);
    }

    std::string selected_file;
    for (const std::string &file_name: possible_files) {
      if (aurostd::EWithinList(available_files, file_name, selected_file)) {
        if (m_filesystem_available) aurostd::efile2string(base_folder + selected_file, poscar);
        else poscar = getRawRestAPIQuery(base_url + selected_file, true);
        if (!poscar.empty()) { // load from aflow.in
          new_structure = xstructure((std::stringstream) poscar, IOVASP_AUTO);
          return true;
        }
      }
    }
    return false;
  }

  /// @brief load a xstructure directly from a AFLOW lib entry
  /// @param entry AFLOW lib entry
  /// @param new_structure save-to structure
  /// @return xstructure
  /// @note this is always the final structure
  /// @note does not add the structure to entry.vstr
  bool EntryLoader::loadXstructureLibEntry(const aflowlib::_aflowlib_entry &entry, xstructure &new_structure) {
    if (entry.positions_cartesian.empty()) return false;
    new_structure.title = "Relaxed AFLUX: " + entry.aurl;
    new_structure.scale = 1.0;
    new_structure.neg_scale = false;
    new_structure.iomode = IOAFLOW_AUTO;
    std::vector<double> geometry;
    aurostd::string2tokens(entry.geometry, geometry, ",");
    new_structure.a = geometry[0];
    new_structure.b = geometry[1];
    new_structure.c = geometry[2];
    new_structure.alpha = geometry[3];
    new_structure.beta = geometry[4];
    new_structure.gamma = geometry[5];
    new_structure.lattice = GetClat(new_structure.a, new_structure.b, new_structure.c,
                                    new_structure.alpha, new_structure.beta, new_structure.gamma);
    new_structure.spacegroupnumber = entry.spacegroup_relax;
    new_structure.spacegrouplabel = GetSpaceGroupLabel(new_structure.spacegroupnumber);
    if (new_structure.spacegroupnumber > 0 && new_structure.spacegroupnumber <= 230) {
      new_structure.spacegroup = GetSpaceGroupName(new_structure.spacegroupnumber, new_structure.directory);
    } else new_structure.spacegroup = "";
    new_structure.coord_flag = _COORDS_FRACTIONAL_;
    strcpy(new_structure.coord_type, "D");
    new_structure.FixLattices();

    std::vector<int> composition;
    aurostd::string2tokens(entry.composition, composition, ",");
    uint num_atoms = 0;
    std::vector <std::vector<double>> positions_fractional;
    std::vector<double> positions_fractional_raw;
    std::vector <std::string> positions_fractional_raw_str;
    aurostd::string2tokens(entry.positions_fractional, positions_fractional_raw_str, ";");
    for (std::string pf_str: positions_fractional_raw_str) {
      positions_fractional_raw.clear();
      aurostd::string2tokens(pf_str, positions_fractional_raw, ",");
      positions_fractional.emplace_back(positions_fractional_raw);
    }

    std::vector <std::string> species;
    aurostd::string2tokens(entry.species, species, ",");
    xvector<double> v(3);
    for (uint type_idx = 0; type_idx < composition.size(); type_idx++) {
      for (int atom_idx = 0; atom_idx < composition[type_idx]; atom_idx++) {
        _atom atom;
        v.clear();
        v(1) = positions_fractional[num_atoms][0];
        v(2) = positions_fractional[num_atoms][1];
        v(3) = positions_fractional[num_atoms][2];
        atom.fpos = v;
        atom.cpos = new_structure.f2c * atom.fpos;
        atom.basis = num_atoms;
        atom.ijk(1) = 0;
        atom.ijk(2) = 0;
        atom.ijk(3) = 0;
        atom.corigin(1) = 0.0;
        atom.corigin(2) = 0.0;
        atom.corigin(3) = 0.0; // inside the zero cell
        atom.coord(1) = 0.0;
        atom.coord(2) = 0.0;
        atom.coord(3) = 0.0; // inside the zero cell
        atom.spin = 0.0;
        atom.noncoll_spin.clear();
        atom.type = type_idx;
        atom.order_parameter_value = 0;
        atom.order_parameter_atom = false;
        atom.partial_occupation_value = 1.0;
        atom.partial_occupation_flag = false;
        atom.name = species[type_idx];
        atom.cleanname = species[type_idx];
        atom.name_is_given = true;
        new_structure.AddAtom(atom);
        num_atoms++;
      }
    }
    return true;
  }

  /// @brief load the first structure in an `aflow.in` file
  /// @param entry AFLOW lib entry
  /// @param new_structure save-to structure
  /// @return xstructure
  /// @note this is always the original structure
  /// @note does not add the structure to entry.vstr
  bool EntryLoader::loadXstructureAflowIn(const aflowlib::_aflowlib_entry &entry, xstructure &new_structure) {
    std::string base_url = m_restapi_server + m_restapi_path + entry.aurl.substr(28) + "/";
    std::string base_folder = m_filesystem_path + entry.aurl.substr(28) + "/";
    base_folder = std::regex_replace(base_folder, m_re_aurl2file, "$1/" + m_filesystem_collection + "/");
    std::string aflowin_content;

    if (m_filesystem_available) {
      std::stringstream buffer;
      std::ifstream open_file(base_folder + "aflow.in");
      buffer << open_file.rdbuf();
      aflowin_content = buffer.str();
    } else { // load form REST API
      m_out_super_silent = true;
      aflowin_content = getRawRestAPIQuery(base_url + "aflow.in", true);
      m_out_super_silent = false;
    }
    std::string poscar = aflowin2poscar(aflowin_content);
    if (!poscar.empty()) { // load from aflow.in
      new_structure = xstructure((std::stringstream) poscar, IOVASP_AUTO);
      return true;
    }
    return false;
  }

  /// @brief save the shared pointers to the AFLOW lib entries into an external vector
  /// \param result save-to vector
  /// @note the underlying entries will not be copied and are likely not in a continuous chunk of memory
  /// @note creating a copy the smart pointer #m_entries_flat is a bit more efficient way to use
  ///       the loaded entries after the EntryLoader class goes out of scope
  void EntryLoader::getEntriesViewFlat(std::vector <std::shared_ptr<aflowlib::_aflowlib_entry>> &result) const {
    result = *m_entries_flat;
  }

  /// @brief save the shared pointers to the AFLOW lib entries into a two layer vector
  /// @param result save-to vector
  /// @note the underlying entries will not be copied and are likely not in a continuous chunk of memory
  /// @note the entries are grouped by number of elements (smallest to largest)
  void EntryLoader::getEntriesViewTwoLayer(vector<vector<std::shared_ptr < aflowlib::_aflowlib_entry>> > &result) const{
    for (auto layer1: *m_entries_layered_map) {
      std::vector <std::shared_ptr<aflowlib::_aflowlib_entry>> collected_entries;
      for (auto layer2: layer1.second) {
        for (auto entry: layer2.second) collected_entries.push_back(entry);
      }
      result.push_back(collected_entries);
    }
  }

  /// @brief save the shared pointers to the AFLOW lib entries into a three layer vector
  /// @param result save-to vector
  /// @note the underlying entries will not be copied and are likely not in a continuous chunk of memory
  /// @note the entries are grouped first by number of elements (smallest to largest) and then by alloy
  void EntryLoader::getEntriesViewThreeLayer(std::vector<std::vector<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>> &result) const{
    for (auto layer1: *m_entries_layered_map) {
      std::vector < std::vector < std::shared_ptr < aflowlib::_aflowlib_entry>>> collected_entries_l1;
      for (auto layer2: layer1.second) {
        std::vector <std::shared_ptr<aflowlib::_aflowlib_entry>> collected_entries_l2;
        for (auto entry: layer2.second) collected_entries_l2.push_back(entry);
        collected_entries_l1.push_back(collected_entries_l2);
      }
      result.push_back(collected_entries_l1);
    }
  }

  /// @brief save the shared pointers to the AFLOW lib entries into a three layer map
  /// @param result save-to map
  /// @note the underlying entries will not be copied and are likely not in a continuous chunk of memory
  /// @note the entries are grouped first by number of elements and then by alloy
  /// @note creating a copy of the smart pointer #m_entries_layered_map is a bit more efficient way to use
  ///       the loaded entries after the EntryLoader class goes out of scope
  void EntryLoader::getEntriesViewMap(std::map < short, std::map < std::string,
                                      std::vector < std::shared_ptr < aflowlib::_aflowlib_entry >> >> &result) const {
    result = *m_entries_layered_map;
  }

  /// @brief copy the AFLOW lib entries into a vector
  /// @param result save-to vector
  /// @note the underlying entries are copied into a continuous chunk of memory
  void EntryLoader::getEntriesFlat(std::vector <aflowlib::_aflowlib_entry> &result) const {
    for (auto entry: *m_entries_flat) {
      result.push_back(*entry);
    }
  }

  /// @brief copy the AFLOW lib entries into a two layer vector
  /// @param result save-to vector
  /// @note the underlying entries are copied into a continuous chunk of memory
  /// @note the entries are grouped by number of elements (smallest to largest)
  void EntryLoader::getEntriesTwoLayer(std::vector <std::vector<aflowlib::_aflowlib_entry>> &result) const {
    for (auto layer1: *m_entries_layered_map) {
      std::vector <aflowlib::_aflowlib_entry> collected_entries;
      for (auto layer2: layer1.second) {
        for (auto entry: layer2.second) collected_entries.push_back(*entry);
      }
      result.push_back(collected_entries);
    }
  }


  /// @brief copy the AFLOW lib entries into a three layer vector
  /// @param result save-to vector
  /// @note the underlying entries are copied into a continuous chunk of memory
  /// @note the entries are grouped first by number of elements (smallest to largest) and then by alloy
  void EntryLoader::getEntriesThreeLayer(std::vector<std::vector<vector<aflowlib::_aflowlib_entry>>> & result) const {
    for (auto layer1: *m_entries_layered_map) {
      std::vector <std::vector<aflowlib::_aflowlib_entry>> collected_entries_l1;
      for (auto layer2: layer1.second) {
        std::vector <aflowlib::_aflowlib_entry> collected_entries_l2;
        for (auto entry: layer2.second) collected_entries_l2.push_back(*entry);
        collected_entries_l1.push_back(collected_entries_l2);
      }
      result.push_back(collected_entries_l1);
    }
  }

  // private functions
  /// @brief find the best available source
  void EntryLoader::selectSource() {

    if (m_current_source == Source::NONE || m_current_source == Source::FAILED) {
      m_out_super_silent = true;
      if (EntryLoader::setSource(Source::SQLITE)) {
        m_out_super_silent = false;
        m_logger_message << "Automatically selected Source::SQLITE";
        outInfo(__func__);
        return;
      }
      if (EntryLoader::setSource(Source::AFLUX)) {
        m_out_super_silent = false;
        m_logger_message << "Automatically selected Source::AFLUX";
        outInfo(__func__);
        return;
      }
      if (EntryLoader::setSource(Source::FILESYSTEM)) {
        m_out_super_silent = false;
        m_logger_message << "Automatically selected";
        if (m_current_source==Source::FILESYSTEM) m_logger_message << " Source::FILESYSTEM";
        else if (m_current_source==Source::FILESYSTEM_RAW) m_logger_message << " Source::FILESYSTEM_RAW";
        outInfo(__func__);
        return;
      }
      if (EntryLoader::setSource(Source::RESTAPI)) {
        m_out_super_silent = false;
        m_logger_message << "Automatically selected";
        if (m_current_source==Source::RESTAPI) m_logger_message << " Source::RESTAPI";
        else if (m_current_source==Source::RESTAPI_RAW) m_logger_message << " Source::RESTAPI_RAW";
        outInfo(__func__);
        return;
      }
      EntryLoader::setSource(Source::FAILED);
      m_out_super_silent = false;
      m_logger_message << "Failed to find a working source!";
      outHardError(__func__, __LINE__, _GENERIC_ERROR_);
    }
  }

  /// @brief list directories or files for the REST API
  /// @param url entry url
  /// @param result save-to vector
  /// @param directories list directories if `true` (default)
  void EntryLoader::listRestAPI(std::string url, std::vector<std::string> &result, bool directories) {
    result.clear();
    if (directories) url += m_restapi_listing_dirs;
    else url += m_restapi_listing_files;
    std::string output;
    if (200 == aurostd::httpGetStatus(url, output)) {
      if (output[output.size()-1]=='\n') output.erase(output.size()-1);
      aurostd::string2tokens(output, result, ",");
    } else {
      m_logger_message << "Could not list content for " << url;
      outInfo(__func__);
    };
  }

  /// @brief generate a list of AUID for a given list of alloys by querying the public alloy SQLITE DB
  /// @param alloy_list list of cleaned and sorted alloys
  /// @param auid_list resulting AUIDs
  void EntryLoader::getAlloyAUIDList(const std::vector<std::string> & alloy_list, std::vector<std::string> &auid_list) {
    if (m_sqlite_alloy_db_ptr == nullptr) {
      m_logger_message << "Public alloy SQLITE DB is not ready";
      outHardError(__func__, __LINE__, _RUNTIME_ERROR_);
      return;
    }
    auid_list.clear();
    std::string where = aurostd::joinWDelimiter(alloy_list, "','");
    where = "alloy IN ('" + where + "')";
    std::vector <std::string> raw_lines = m_sqlite_alloy_db_ptr->getEntrySet(where, aflow_ft);
    for (std::string line: raw_lines) {
      auid_list.push_back(line.substr(5, 22));
    }
  }

  /// @brief load AFLOW lib entries for a given alloy directly from the filesystem (source FILESYSTEM_RAW)
  /// @param alloy_list list of cleaned and sorted alloys
  /// @param lib_max number of elements in the largest alloy
  /// @note if the public alloy SQLITE DB is available the results are expanded with missed entries
  void EntryLoader::loadAlloySearchFSR(const std::vector <std::string> &alloy_list, uint lib_max) {
    std::vector <std::string> found_entries;
    std::vector <std::string> search_path_list;
    search_path_list.emplace_back(m_filesystem_path + "ICSD/" + m_filesystem_collection + "/");
    if (alloy_list.size()>1) {
      for (uint lib_idx = 1; lib_idx <= lib_max; lib_idx++) {
        search_path_list.emplace_back(
            m_filesystem_path + "LIB" + std::to_string(lib_idx) + "/" + m_filesystem_collection + "/");
      }
    } else if (alloy_list.size()==1){
      search_path_list.emplace_back(
          m_filesystem_path + "LIB" + std::to_string(lib_max) + "/" + m_filesystem_collection + "/");
    } else {
      return;
    }

    char *paths[search_path_list.size() + 1];
    for (size_t list_idx = 0; list_idx < search_path_list.size(); list_idx++) {
      paths[list_idx] = (char *) search_path_list[list_idx].c_str();
    }
    paths[search_path_list.size() + 1] = nullptr;

    size_t scanned = 0;
    size_t found = 0;
    size_t check_idx = m_filesystem_path.size();
    struct stat file_stats{};

    // initialize a file tree scan
    // https://man7.org/linux/man-pages/man3/fts.3.html
    // FTS_PHYSICAL - don't follow symlinks
    // FTS_NOCHDIR - don't change the workdir of the program
    // FTS_XDEV - don't descend into folders that are on a different device
    FTS *tree = fts_open(paths, FTS_PHYSICAL | FTS_NOCHDIR | FTS_XDEV, nullptr);
    FTSENT *node;
    if (tree == nullptr) {
      m_logger_message << "Failed to initialized the file tree used to search for alloys!";
      outHardError(__func__,__LINE__,_FILE_ERROR_);
      return;
    }

    // Iterate over all entries found in the folders listed in paths
    while ((node = fts_read(tree))) {
      scanned += 1;
      if ((scanned-1) % 10000 == 0) {
        m_logger_message << (scanned-1) << " objects scanned; " << found << " entries found; next scan: " << node->fts_path;
        outDebug(__func__);
      }
      if ((node->fts_info & FTS_D)) {
        if (node->fts_path[check_idx] == 'I') { // ICSD
          if (node->fts_level == 2) {
            if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(node->fts_name, 'I')) == alloy_list.end()) {
              fts_set(tree, node, FTS_SKIP);
              continue;
            } else {
              std::string base_path = node->fts_path;
              std::string full_path = base_path + "/" + m_filesystem_outfile;
              if (stat(full_path.c_str(), &file_stats) == 0) {
                found += 1;
                found_entries.emplace_back(full_path);
                if (base_path.find("POCC") == std::string::npos) fts_set(tree, node, FTS_SKIP);
              }
            }
          }
        } else if (node->fts_path[check_idx] == 'L') { // LIBX
          if (node->fts_level == 1) {
            if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(node->fts_name, 'L')) == alloy_list.end()) {
              fts_set(tree, node, FTS_SKIP);
              continue;
            } else {
              std::string base_path = node->fts_path;
              std::string full_path = base_path + "/" + m_filesystem_outfile;
              if (stat(full_path.c_str(), &file_stats) == 0) {
                found += 1;
                found_entries.emplace_back(full_path);
                if (base_path.find("POCC") == std::string::npos) fts_set(tree, node, FTS_SKIP);
              }
            }
          } else if (node->fts_level > 1) {
            std::string base_path = node->fts_path;
            std::string full_path = base_path + "/" + m_filesystem_outfile;
            if (stat(full_path.c_str(), &file_stats) == 0) {
              found += 1;
              found_entries.emplace_back(full_path);
              if (base_path.find("POCC") == std::string::npos) fts_set(tree, node, FTS_SKIP);
            }
          }
        }
      }
    }
    m_logger_message << "Finishing search in the filesystem after scanning " << scanned << " objects";
    outInfo(__func__);
    loadFiles(found_entries);
    if (m_sqlite_alloy_db_ptr != nullptr) {
      std::vector <std::string> known_AUID_list;
      getAlloyAUIDList(alloy_list, known_AUID_list);
      std::vector <std::string> missing_AUID;
      for (const std::string &AUID: known_AUID_list) {
        if (std::find(m_auid_list.begin(), m_auid_list.end(), AUID) == m_auid_list.end()) missing_AUID.emplace_back(AUID);
      }
      m_logger_message << "Found " << missing_AUID.size() << " entries in the public alloy SQLITE DB that where not found in the filesystem search.";
      outDebug(__func__);
      loadAUID(missing_AUID);
    }
  }

  /// @brief load AFLOW lib entries for a given alloy set directly from the AFLOW REST API (source RESTAPI_RAW)
  /// @param alloy_list list of cleaned and sorted alloys
  /// @param lib_max number of elements in the largest alloy
  void EntryLoader::loadAlloySearchRR(const std::vector <std::string> & alloy_list, uint lib_max) {
    std::vector<std::string> rest_api_queries;
    std::vector<std::string> icsd_search_path_list;
    std::vector<std::string> libx_search_path_list;
    std::vector<std::string> libx_crawl_list;
    for (const std::string & sym : m_icsd_symmetries) {
      icsd_search_path_list.emplace_back(m_restapi_server + m_restapi_path + "ICSD" + "_" + m_restapi_collection + "/" + sym);
    }
    if (alloy_list.size()>1) {
      for (uint lib_idx = 1; lib_idx <= lib_max; lib_idx++) {
        libx_search_path_list.emplace_back(m_restapi_server + m_restapi_path + "LIB" + std::to_string(lib_idx) + "_" + m_restapi_collection + "/" );
      }
    } else if (alloy_list.size()==1){
      libx_search_path_list.emplace_back(m_restapi_server + m_restapi_path + "LIB" + std::to_string(lib_max) + "_" + m_restapi_collection + "/" );
    } else {
      return;
    }
    std::vector<std::string> listing;

    for (const std::string & url : icsd_search_path_list){
      listRestAPI(url, listing);
      for (const std::string & name: listing ){
        if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(name, 'I')) != alloy_list.end()){
          rest_api_queries.emplace_back(url + name + "/" + m_restapi_directives);
        }
      }
    }

    m_logger_message << "Found " << rest_api_queries.size() << " ICSD entries in REST API search";
    outDebug(__func__);

    for (const std::string & url : libx_search_path_list){
      listRestAPI(url, listing);
      for (const std::string & name: listing){
        if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(name, 'L')) != alloy_list.end()){
          libx_crawl_list.emplace_back(url + name + "/" );
        }
      }
    }

    m_logger_message << "Found " << libx_crawl_list.size() << " LIBX base folders in the REST API search";
    outDebug(__func__);

    uint scanned=0;
    uint found=0;
    while(!libx_crawl_list.empty()){
      scanned += 1;
      std::string url = libx_crawl_list.back();
      libx_crawl_list.pop_back();
      if ((scanned-1)%100 == 0) {
        m_logger_message << (scanned-1) << " folders scanned; " << found << " entries found; next scan: " << url;
        outDebug(__func__);
      }
      listRestAPI(url, listing);
      if (listing.empty()) {
        rest_api_queries.emplace_back(url + m_restapi_directives);
        found += 1;
      } else {
        for (const std::string & name: listing ) {
          libx_crawl_list.emplace_back(url + name + "/");
        }
      }
    }


    loadRestAPIQueries(rest_api_queries, true);
    if (m_sqlite_alloy_db_ptr != nullptr) {
      std::vector <std::string> known_AUID_list;
      getAlloyAUIDList(alloy_list, known_AUID_list);
      std::vector <std::string> missing_AUID;
      for (const std::string &AUID: known_AUID_list) {
        if (std::find(m_auid_list.begin(), m_auid_list.end(), AUID)== m_auid_list.end()) missing_AUID.emplace_back(AUID);
      }
      m_logger_message << "Found " << missing_AUID.size() << " entries in the public alloy SQLITE DB that where not found in the REST API search.";
      outDebug(__func__);
      loadAUID(missing_AUID);
    }

  }

  /// @brief map different AUID styles to a default
  /// @param AUID AUID to clean
  /// @return `false` if AUID is not valid
  /// @note AUID is not probably validated
  /// @note the AUID is cleaned in place
  bool EntryLoader::cleanAUID(std::string & AUID) {
    if (AUID.substr(0, 5) == "auid:") AUID="aflow:" + AUID.substr(5);
    else if (AUID.substr(0, 6) != "aflow:") AUID="aflow:" + AUID;
    if (aurostd::_ishex(AUID.substr(6)) && AUID.size()==22) return true;
    else return false;
  }

  /// @brief map different AURL styles to a default
  /// @param AUID AURL to clean
  /// @return `false` if AURL is not valid
  /// @note AURL is not validated
  /// @note the AURL is cleaned in place
  bool EntryLoader::cleanAURL(std::string & AURL) {
    if (AURL.substr(0, 5) == "aurl:")  AURL = AURL.substr(5);
    if (AURL.substr(0, 28)== "aflowlib.duke.edu:AFLOWDATA/") return true;
    if (AURL.substr(0, 9) == "AFLOWDATA") {
      AURL = "aflowlib.duke.edu:" + AURL;
      return true;
    }

    AURL = "aflowlib.duke.edu:AFLOWDATA/" + AURL;
    return true;
    //TODO how to check if AURL could be valid
  }

  /// @brief create a AFLUX query string from a matchbook
  /// @param matchbook map with keyword and modifier
  /// @returns AFLUX query
  /// @note does not include #m_aflux_server or #m_aflux_path
  /// @note the string is percent encoded
  std::string EntryLoader::buildAFLUXQuery(const std::map<std::string, std::string> & matchbook) const {
    std::string query = "";
    for (std::map<std::string, std::string> part: {matchbook, m_aflux_directives}) {
      for (std::pair <std::string, std::string> it: part) {
        query += "," + it.first;
        if (!it.second.empty()) query += "(" + it.second + ")";
      }
    }
    query.erase(0,1);
    return "?" + aurostd::httpPercentEncodingFull(query);
  }

  /// @brief extract alloy string from AFLOW run name
  /// @param name AFLOW run name (from path structure)
  /// @param lib_type indicate the LIB type | `I` for `ICSD` or `L` for 'LIBx'
  /// @return alloy string
  /// @note the alloy sting is sorted
  std::string EntryLoader::extractAlloy(std::string name, char lib_type) const {
    if (lib_type == 'I') { //ICSD
      name.erase(std::min(name.find("_ICSD_"), name.size()));
    } else if (lib_type == 'L') { //LIBX
      name.erase(std::min(name.find(':'), name.size()));
      name = std::regex_replace (name, m_re_ppclean,"");
    } else {
      return "";
    }
    std::vector <std::string> element_match(std::sregex_token_iterator(name.begin(), name.end(), m_re_elements),
                                            std::sregex_token_iterator());
    std::sort(element_match.begin(), element_match.end());
    return aurostd::joinWDelimiter(element_match, "");
  }

  /// @brief extract the first POSCAR entry from `aflow.in`
  /// \param raw_in content of `aflow.in`
  /// \return poscar content
  /// @note can only read `VASP_POSCAR_MODE_EXPLICIT` variants
  std::string EntryLoader::aflowin2poscar(std::string raw_in) const {
    size_t explicit_start = raw_in.find("[VASP_POSCAR_MODE_EXPLICIT]");
    std::string poscar;
    if (explicit_start != std::string::npos) {
      raw_in.erase(0, explicit_start);
      size_t explicit_stop = raw_in.find("[VASP_POSCAR_MODE_EXPLICIT]STOP");
      // START/STOP variant
      if (explicit_stop != std::string::npos) {
        poscar = raw_in.substr(38, explicit_stop - 38);
        // VASP_POSCAR_FILE variant
      } else {
        std::vector <std::string> raw_lines;
        aurostd::string2vectorstring(raw_in, raw_lines);
        for (std::string &line: raw_lines) {
          if (line.substr(0, 18) == "[VASP_POSCAR_FILE]") poscar += line.substr(16) + "\n";
          else break;
        }
      }
    }
    return poscar;
  }

  // TODO why _AFLOW_FILE_NAME_ and not __FILE__?
  // https://www.cprogramming.com/reference/preprocessor/__FILE__.html
  // also the use of __func__ and __LINE__ would be useful

  void EntryLoader::outInfo(const std::string & function_name) {
    if ((!m_out_silent || m_out_debug) && !m_out_super_silent){
      std::string soliloquy = XPID + "EntryLoader::" + function_name + "():";
      pflow::logger(_AFLOW_FILE_NAME_, soliloquy, m_logger_message, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);
    }
    aurostd::StringstreamClean(m_logger_message);
  }

  void EntryLoader::outDebug(const std::string & function_name) {
    if (m_out_debug && !m_out_super_silent){
      std::string soliloquy = XPID + "EntryLoader::" + function_name + "():";
      pflow::logger(_AFLOW_FILE_NAME_, soliloquy, m_logger_message, *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    }
    aurostd::StringstreamClean(m_logger_message);
  }

  void EntryLoader::outError(const std::string &function_name, int line_number) {
    if (!m_out_super_silent) {
      std::string soliloquy = XPID + "EntryLoader::" + function_name + "():";
      pflow::logger((std::string) _AFLOW_FILE_NAME_ + ":" + std::to_string(line_number),
                    soliloquy, m_logger_message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
    }
    aurostd::StringstreamClean(m_logger_message);
  }

  void EntryLoader::outHardError(const std::string &function_name, int line_number, int error_type) {
    std::string soliloquy = XPID + "EntryLoader::" + function_name + "():";
    throw aurostd::xerror((std::string) _AFLOW_FILE_NAME_ + ":" + std::to_string(line_number), soliloquy, m_logger_message, error_type);
    aurostd::StringstreamClean(m_logger_message);
  }


}

#endif  // _AFLOWLIB_ENTRY_LOADER_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
