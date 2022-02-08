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


namespace aflowlib {

  //--------------------//
  // Start Entry Loader //
  //--------------------//

  EntryLoader::EntryLoader(ostream& oss) : xStream(oss) {free();}
  EntryLoader::EntryLoader(ofstream& FileMESSAGE,ostream& oss) : xStream(FileMESSAGE,oss) {free();}
  EntryLoader::EntryLoader(const EntryLoader& b) : xStream(*b.getOFStream(),*b.getOSS()) {copy(b);} // copy PUBLIC

  EntryLoader::~EntryLoader() {xStream::free();free();}

  const EntryLoader& EntryLoader::operator=(const EntryLoader& other) {
    if(this!=&other) {copy(other);}
    return *this;
  }

  void EntryLoader::clear() {free();}  //clear PUBLIC
  void EntryLoader::free() {
    m_sqlite_file = DEFAULT_AFLOW_DB_FILE;
    m_aflux_server = "aflowlib.duke.edu";
    m_aflux_path = "/API/aflux/v1.0/";
    //TODO use server & path from aflux rc
    m_aflux_directives = {
        {"format", "aflow"},
        {"paging", "0"}};
    m_current_source = Source::NONE;
    m_entries_unique=true;
    m_auid_list.clear();
    m_sqlite_db_ptr = nullptr;
    m_entries_flat = std::make_shared<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>();
    m_entries_layered_map = std::make_shared <std::unordered_map<short, std::unordered_map<
                                              std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>>();
  }

  void EntryLoader::copy(const EntryLoader& b) {  //copy PRIVATE
    m_sqlite_file = b.m_sqlite_file;
    m_aflux_server = b.m_aflux_server;
    m_aflux_path = b.m_aflux_path;
    //TODO use server & path from aflux rc
    m_aflux_directives = b.m_aflux_directives;
    m_current_source = b.m_current_source;
    m_entries_unique=b.m_entries_unique;
    m_auid_list=b.m_auid_list;
    m_sqlite_db_ptr = b.m_sqlite_db_ptr;
    m_entries_flat = b.m_entries_flat;
    m_entries_layered_map = b.m_entries_layered_map;
  }

  void EntryLoader::initialize(ostream& oss) {
    xStream::initialize(oss);
    return initialize();
  }

  void EntryLoader::initialize(ofstream& FileMESSAGE,ostream& oss) {
    xStream::initialize(FileMESSAGE,oss);
    return initialize();
  }

  void EntryLoader::initialize() {
    free();
  }


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

  bool EntryLoader::cleanAUID(string & AUID) {
    if (AUID.substr(0, 5) == "auid:") AUID="aflow:" + AUID.substr(5);
    else if (AUID.substr(0, 6) != "aflow:") AUID="aflow:" + AUID;
    if (aurostd::_ishex(AUID.substr(6)) && AUID.size()==22) return true;
    else {
        cout << AUID << " has the wrong size" << endl; //TODO error logger
        return false;
    }
  }

  bool EntryLoader::setSource(EntryLoader::Source new_source) {
    // TODO error messages
    if (new_source == m_current_source) return true;
    m_current_source = Source::FAILED;

    switch (new_source) {

      case Source::SQLITE: {
        if (aurostd::FileExist(m_sqlite_file)) {
          if (m_sqlite_db_ptr == nullptr) {
            cout << "Init full DB" << endl;
            m_sqlite_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_file);
          }
          m_current_source = Source::SQLITE;
          return true;
        }
        break;
      }
      case Source::AFLUX: {
        if ("AFLUXtest" == aurostd::httpGet("http://aflowlib.duke.edu/test/?echo=AFLUXtest")) {
          m_current_source = Source::AFLUX;
          return true;
        }
        break;
      }
      case Source::FILESYSTEM: {
        //TODO if m_sqlite_alloy_file is missing AURL and AUID load should still work
        if (aurostd::IsDirectory(m_filesystem_path + "AUID/") && aurostd::FileExist(m_sqlite_alloy_file)) { //
          if (m_sqlite_alloy_db_ptr == nullptr) {
            cout << "Init alloy DB" << endl;
            m_sqlite_alloy_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_alloy_file);
          }
          m_current_source = Source::FILESYSTEM;
          return true;
        }
        break;
      }
      case Source::RESTAPI: {
        if (aurostd::FileExist(m_sqlite_alloy_file) &&
            200 == aurostd::httpGetStatus("http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/")) {
          if (m_sqlite_alloy_db_ptr == nullptr) {
            cout << "Init alloy DB" << endl;
            m_sqlite_alloy_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_alloy_file);
          }
          m_current_source = Source::RESTAPI;
          return true;
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

  void EntryLoader::selectSource(){
    if (m_current_source==Source::NONE || m_current_source==Source::FAILED) {
      if (EntryLoader::setSource(Source::SQLITE)) return;
      if (EntryLoader::setSource(Source::AFLUX)) return;
      if (EntryLoader::setSource(Source::FILESYSTEM)) return;
      if (EntryLoader::setSource(Source::RESTAPI)) return;
      EntryLoader::setSource(Source::FAILED);
    }
  }

  void EntryLoader::getAlloyAUIDList(const std::vector<string> & final_alloy_list, std::vector<string> & auid_list) {
    if (m_sqlite_alloy_db_ptr == nullptr) {
      cout << "Small DB not loaded" << endl; //TODO error handling
      return;
    }
    auid_list.clear();
    std::string where = aurostd::joinWDelimiter(final_alloy_list, "','");
    where = "alloy IN ('" + where + "')";
    vector<string> raw_lines = m_sqlite_alloy_db_ptr->getEntrySet(where, aflow_ft);
    for (string line: raw_lines){
      auid_list.push_back(line.substr(5,22));
    }
  }


  void EntryLoader::loadText(const std::vector<std::string> & raw_data_lines) {
    if(m_entries_unique) {
      for (std::string line: raw_data_lines) {
        std::shared_ptr<aflowlib::_aflowlib_entry> entry = std::make_shared<aflowlib::_aflowlib_entry>();
        entry->Load(line, *p_oss);
        if (!entry->auid.empty() && (m_auid_list.find(entry->auid)==m_auid_list.end())) {
          m_entries_flat->push_back(entry);
          (* m_entries_layered_map)[entry->nspecies][entry->species_pp].push_back(entry);
          m_auid_list.emplace(entry->auid);
        }
      }
    } else {
      for (std::string line: raw_data_lines) {
        std::shared_ptr<aflowlib::_aflowlib_entry> entry = std::make_shared<aflowlib::_aflowlib_entry>();
        entry->Load(line, *p_oss);
        if (!entry->auid.empty()) {
          m_entries_flat->push_back(entry);
          (* m_entries_layered_map)[entry->nspecies][entry->species_pp].push_back(entry);
          m_auid_list.emplace(entry->auid);
        }
      }
    }
  }

  void EntryLoader::loadSqliteWhere(const std::string & where) {
    stringstream message;
    string soliloquy=XPID+"EntryLoader::loadSqliteWhere():";
    vector<string> raw_lines = m_sqlite_db_ptr->getEntrySet(where, aflow_ft);
    cout << "Raw result size " << raw_lines.size() << " for "<< where << endl;

    size_t start_size=m_entries_flat->size();
    loadText(raw_lines);

    message << "Loaded " << m_entries_flat->size()-start_size << " new entries (overall "<< m_entries_flat->size() << " | " << m_auid_list.size() << " unique)";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
  }


  void EntryLoader::loadRestAPIQueries(const std::vector <std::string> &queries) {
    string soliloquy=XPID+"EntryLoader::loadRestAPIQueries():";
    stringstream message;
    std::string output = "";
    std::string url ="";
    cout << "Loading " << queries.size() << " entries:" << endl;
    uint n = 0;
    size_t start_size=m_entries_flat->size();
    for (string query: queries){
      n+=1;
      url = m_restapi_server+m_restapi_path + query + m_restapi_directives;
      cout << " Load #" << n << " from " << url << "\n";
      short status = aurostd::httpGetStatus(url, output);
      //TODO throw ERROR if web fails or try different source
      loadText({output});
    }
    cout << endl;

    message << "Loaded " << m_entries_flat->size()-start_size << " new entries (overall "<< m_entries_flat->size() << " | " << m_auid_list.size() << " unique)";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
  }

  void EntryLoader::loadFiles(const std::vector <std::string> &files) {
    string soliloquy=XPID+"EntryLoader::loadFiles():";
    stringstream message;
    std::string file_content;
    cout << "Loading " << files.size() << " entries:" << endl;
    size_t start_size=m_entries_flat->size();
    for (string file_path : files){
//      cout << "file path: " << file_path << endl;
      aurostd::efile2string(file_path, file_content);
      loadText({file_content});
    }
    message << "Loaded " << m_entries_flat->size()-start_size << " new entries (overall "<< m_entries_flat->size() << " | " << m_auid_list.size() << " unique)";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
  }

  void EntryLoader::loadAFLUXMatchbook(const std::map<string, string> & matchbook){
    loadAFLUXQuery(buildAFLUXQuery(matchbook));
  }

  void EntryLoader::loadAFLUXQuery(const std::string &query) {
    string soliloquy=XPID+"EntryLoader::loadAFLUXQuery():";
    stringstream message;

    std::string output = "";
    short status = aurostd::httpGetStatus(m_aflux_server, m_aflux_path, query, output);
    //TODO throw ERROR if web fails or try different source

    vector<string> raw_lines;
    aurostd::string2vectorstring(output,raw_lines);
    size_t start_size=m_entries_flat->size();
    loadText(raw_lines);

    message << "Loaded " << m_entries_flat->size()-start_size << " new entries (overall "<< m_entries_flat->size() << " | " << m_auid_list.size() << " unique)";
    pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,m_aflags,*p_FileMESSAGE,*p_oss,_LOGGER_NOTICE_);
  }

  void EntryLoader::loadAUID(string AUID) {
    selectSource();

    if (!cleanAUID(AUID)) return; //TODO error message

    switch (m_current_source) {

      case Source::SQLITE: {
        std::string where = "auid='\"" + AUID + "\"'";
        loadSqliteWhere(where);
        break;
        }

      case Source::AFLUX: {
        std::map <std::string, std::string> matchbook{{"*", ""},{"auid", "'" + AUID + "'"}};
        loadAFLUXMatchbook(matchbook);
        break;
        }

      case Source::RESTAPI:{
        loadAUID({AUID});
        break;
      }

      case Source::FILESYSTEM:{
        loadAUID({AUID});
        break;
      }

      default:
        break;
      }
  }

  void EntryLoader::loadAUID(const vector<string> &AUID) {
    selectSource();
    vector<std::string> clean_AUID;
    for (std::string AUID_single: AUID) {
      if (cleanAUID(AUID_single)) clean_AUID.push_back(AUID_single);
    }


    switch (m_current_source) {

      case Source::AFLUX:{
        std::string AUID_combined = aurostd::joinWDelimiter(clean_AUID, "':'");
        AUID_combined = "'"+ AUID_combined + "'";
        loadAFLUXMatchbook({{"*", ""}, {"auid", AUID_combined}});
      }

      case Source::SQLITE:{
        std::string where = aurostd::joinWDelimiter(clean_AUID, "\"','\"");
        where = "auid IN ('\"" + where + "\"')";
        loadSqliteWhere(where);
        break;
      }

      case Source::RESTAPI:{
        std::vector<std::string> queries;
        for (std::string AUID_single: clean_AUID) {
          string rest_query = "AUID/aflow:";
          for (uint part = 0; part < 8; part++) rest_query += AUID_single.substr(6 + part * 2, 2) + "/";
          rest_query += "WEB/";
          queries.push_back(rest_query);
        }
        loadRestAPIQueries(queries);
        break;
      }

      case Source::FILESYSTEM:{
        std::vector<std::string> files;
        for (std::string AUID_single: clean_AUID) {
          string file_path = m_filesystem_path + "AUID/aflow:";
          for (uint part = 0; part < 8; part++) file_path += AUID_single.substr(6 + part * 2, 2) + "/";
          file_path += "WEB/" + m_filesystem_outfile;
          files.push_back(file_path);
        }
        loadFiles(files);
        break;
      }

      default:
        break;
    }
  }

  void EntryLoader::loadAlloy(const vector <std::string> & alloy, bool recursive) {
    selectSource();
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
    switch(m_current_source){
      case Source::SQLITE: {
        std::string where = aurostd::joinWDelimiter(final_alloy_list, "','");
        where = "alloy IN ('" + where + "')";
        loadSqliteWhere(where);
        break;
      }

      case Source::AFLUX: {
        std::string alloy_match = aurostd::joinWDelimiter(final_alloy_list, ":");
        std::map <std::string, std::string> matchbook{{"*", ""},{"alloy", alloy_match}};
        loadAFLUXMatchbook(matchbook);
        break;
      }

      case Source::RESTAPI: {
        std::vector<string> auid_list;
        getAlloyAUIDList(final_alloy_list, auid_list);
        loadAUID(auid_list);
        break;
      }

      case Source::FILESYSTEM: {
        std::vector<string> auid_list;
        getAlloyAUIDList(final_alloy_list, auid_list);
        loadAUID(auid_list);
        break;
      }
      default:
        break;
    }

  }


  void EntryLoader::loadAlloy(const std::string & alloy, bool recursive) {
    vector <string> alloy_elements;
    if (alloy.find(",") != string::npos) aurostd::string2tokens(alloy, alloy_elements, ",");
    else alloy_elements = aurostd::getElements(alloy);
    loadAlloy(alloy_elements, recursive);
  }


  // View getter functions
  void EntryLoader::getEntriesViewFlat(std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>> & result) {
    result = *m_entries_flat;
  }

  void EntryLoader::getEntriesViewTwoLayer(std::vector<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>> & result) {
    for (auto layer1: *m_entries_layered_map){
      std::vector<std::shared_ptr < aflowlib::_aflowlib_entry>> collected_entries;
      for (auto layer2: layer1.second){
        for (auto entry: layer2.second) collected_entries.push_back(entry);
      }
      result.push_back(collected_entries);
    }
  }

  void EntryLoader::getEntriesViewThreeLayer(std::vector<std::vector<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>> & result) {
    for (auto layer1: *m_entries_layered_map){
      std::vector<std::vector<std::shared_ptr < aflowlib::_aflowlib_entry>>> collected_entries_l1;
      for (auto layer2: layer1.second){
        std::vector<std::shared_ptr < aflowlib::_aflowlib_entry>> collected_entries_l2;
        for (auto entry: layer2.second) collected_entries_l2.push_back(entry);
        collected_entries_l1.push_back(collected_entries_l2);
      }
      result.push_back(collected_entries_l1);
    }
  }

  void EntryLoader::getEntriesViewMap(std::unordered_map<short, std::unordered_map<std::string, std::vector < std::shared_ptr < aflowlib::_aflowlib_entry>>>> & result) {
    result = *m_entries_layered_map;
  }

  // Copy getter functions
   void EntryLoader::getEntriesFlat(std::vector <aflowlib::_aflowlib_entry> & result) {
    for (auto entry: *m_entries_flat){
      result.push_back(*entry);
    }
  }

  void EntryLoader::getEntriesTwoLayer(std::vector<std::vector<aflowlib::_aflowlib_entry>> & result) {
    for (auto layer1: *m_entries_layered_map){
      std::vector<aflowlib::_aflowlib_entry> collected_entries;
      for (auto layer2: layer1.second){
        for (auto entry: layer2.second) collected_entries.push_back(*entry);
      }
      result.push_back(collected_entries);
    }
  }

  void EntryLoader::getEntriesThreeLayer(std::vector<std::vector<std::vector<aflowlib::_aflowlib_entry>>> & result) {
    for (auto layer1: *m_entries_layered_map){
      std::vector<std::vector<aflowlib::_aflowlib_entry>> collected_entries_l1;
      for (auto layer2: layer1.second){
        std::vector<aflowlib::_aflowlib_entry> collected_entries_l2;
        for (auto entry: layer2.second) collected_entries_l2.push_back(*entry);
        collected_entries_l1.push_back(collected_entries_l2);
      }
      result.push_back(collected_entries_l1);
    }
  }
}

#endif  // _AFLOW_ENTRY_LOADER_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *           Aflow COREY OSES - Duke University 2013-2020                  *
// *                                                                         *
// ***************************************************************************
