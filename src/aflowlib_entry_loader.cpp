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
#include <fts.h>
#include <sys/stat.h>

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

  bool EntryLoader::cleanAURL(string & AURL) {
    if (AURL.substr(0, 28)== "aflowlib.duke.edu:AFLOWDATA/") return true;
    if (AURL.substr(0, 9) == "AFLOWDATA") {
      AURL = "aflowlib.duke.edu:" + AURL;
      return true;
    }

    AURL = "aflowlib.duke.edu:AFLOWDATA/" + AURL;
    return true;

    //TODO how to check if AURL could be valid
  }

  void EntryLoader::listRestAPI(string url, vector<string> & result) {
    result.clear();
    url+= m_restapi_listing;
    string output;
    short status = aurostd::httpGetStatus(url, output);
    if (status==200) aurostd::string2tokens(output, result, ",");
    else cout << status << " | " << url << endl;
  }

  string EntryLoader::extractAlloy(std::string name, char lib_type) {
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
        } else {
          return setSource(Source::FILESYSTEM_RAW);
        }
        break;
      }

      case Source::FILESYSTEM_RAW: {
        if (aurostd::IsDirectory(m_filesystem_path + "AUID/")) { //
          if (aurostd::FileExist(m_sqlite_alloy_file)) { // use the alloy DB to optemize the RAW results
            if (m_sqlite_alloy_db_ptr == nullptr) {
              cout << "Init alloy DB" << endl;
              m_sqlite_alloy_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_alloy_file);
            }
          }
          m_current_source = Source::FILESYSTEM_RAW;
          return true;
        }
        break;
      }

      case Source::RESTAPI: {
        if (200 == aurostd::httpGetStatus("http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/")) {
          if (aurostd::FileExist(m_sqlite_alloy_file)) {
            if (m_sqlite_alloy_db_ptr == nullptr) {
              cout << "Init alloy DB" << endl;
              m_sqlite_alloy_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_alloy_file);
            }
            m_current_source = Source::RESTAPI;
            return true;
          } else {
            return setSource(Source::RESTAPI_RAW);
          }
        }
        break;
      }

      case Source::RESTAPI_RAW: {
        if (200 == aurostd::httpGetStatus("http://aflowlib.duke.edu/AFLOWDATA/ICSD_WEB/")) {
          if (m_sqlite_alloy_db_ptr == nullptr) {
            cout << "Init alloy DB" << endl;
            m_sqlite_alloy_db_ptr = std::make_shared<aflowlib::AflowDB>(m_sqlite_alloy_file);
          }
          m_current_source = Source::RESTAPI_RAW;
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


  void EntryLoader::loadRestAPIQueries(const std::vector <std::string> &queries, bool full_url) {
    string soliloquy=XPID+"EntryLoader::loadRestAPIQueries():";
    stringstream message;
    std::string output = "";
    std::string url ="";
    cout << "Loading " << queries.size() << " entries:" << endl;
    uint n = 0;
    size_t start_size=m_entries_flat->size();
    for (string query: queries){
      n+=1;
      if (full_url)url = query;
      else url = m_restapi_server+m_restapi_path + query + m_restapi_directives;
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
//      aurostd::efile2string(file_path, file_content);
// TODO aurostd::efile2string is very slow (10 entries/s)
// this variant is ca 700 E/s
      std::ifstream open_file(file_path);
      std::stringstream buffer;
      buffer << open_file.rdbuf();
      loadText({buffer.str()});
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
    cout << query << endl;
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

  void EntryLoader::loadAURL(string AURL) {
    selectSource();

    if (!cleanAURL(AURL)) return; //TODO error message

    switch (m_current_source) {

      case Source::SQLITE: {
        std::string where = "aurl='\"" + AURL + "\"'";
        loadSqliteWhere(where);
        break;
      }

      case Source::AFLUX: {
        std::map <std::string, std::string> matchbook{{"*", ""},{"aurl", "'" + AURL + "'"}};
        loadAFLUXMatchbook(matchbook);
        break;
      }

      case Source::RESTAPI: case Source::RESTAPI_RAW: {
        loadAURL({AURL});
        break;
      }

      case Source::FILESYSTEM: case Source::FILESYSTEM_RAW: {
        loadAURL({AURL});
        break;
      }

      default:
        break;
    }
  }

  void EntryLoader::loadAURL(const vector <std::string> &AURL) {
    selectSource();
    vector<std::string> clean_AURL;
    for (std::string AURL_single: AURL) {
      if (cleanAURL(AURL_single)) clean_AURL.push_back(AURL_single);
    }
    switch (m_current_source) {

      case Source::AFLUX:{
        std::string AURL_combined = aurostd::joinWDelimiter(clean_AURL, "':'");
        AURL_combined = "'"+ AURL_combined + "'";
        loadAFLUXMatchbook({{"*", ""}, {"aurl", AURL_combined}});
        break;
      }

      case Source::SQLITE:{
        std::string where = aurostd::joinWDelimiter(clean_AURL, "\"','\"");
        where = "aurl IN ('\"" + where + "\"')";
        loadSqliteWhere(where);
        break;
      }

      case Source::RESTAPI: case Source::RESTAPI_RAW:{
        std::vector<std::string> queries;
        for (std::string AURL_single: clean_AURL) {
          string rest_query = AURL_single.erase(0,28);
          queries.push_back(rest_query);
        }
        loadRestAPIQueries(queries);
        break;
      }

      case Source::FILESYSTEM: case Source::FILESYSTEM_RAW:{
        std::vector<std::string> files;
        for (std::string AURL_single: clean_AURL) {
          string file_path = m_filesystem_path + AURL_single.erase(0,28) + "/" + m_filesystem_outfile;
          file_path = std::regex_replace (file_path, m_re_aurl2file,"$1/"+m_filesystem_collection);
          files.push_back(file_path);
          cout << file_path << endl;
        }
        loadFiles(files);
        break;

      }

      default:
        break;
    }
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

        case Source::RESTAPI: case Source::RESTAPI_RAW: {
        loadAUID({AUID});
        break;
      }

      case Source::FILESYSTEM: case Source::FILESYSTEM_RAW: {
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
        break;
      }

      case Source::SQLITE:{
        std::string where = aurostd::joinWDelimiter(clean_AUID, "\"','\"");
        where = "auid IN ('\"" + where + "\"')";
        loadSqliteWhere(where);
        break;
      }

      case Source::RESTAPI: case Source::RESTAPI_RAW:{
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

      case Source::FILESYSTEM: case Source::FILESYSTEM_RAW:{
        std::vector<std::string> files;
        for (std::string AUID_single: clean_AUID) {
          string file_path = m_filesystem_path + "AUID/aflow:";
          for (uint part = 0; part < 8; part++) file_path += AUID_single.substr(6 + part * 2, 2) + "/";
          file_path += m_filesystem_collection + m_filesystem_outfile;
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

      case Source::RESTAPI_RAW: {
        loadAlloySearchRR(final_alloy_list, alloy_clean.size(), recursive);
        break;
      }

      case Source::FILESYSTEM: {
        std::vector<string> auid_list;
        getAlloyAUIDList(final_alloy_list, auid_list);
        loadAUID(auid_list);
        break;
      }

      case Source::FILESYSTEM_RAW: {
        loadAlloySearchFSR(final_alloy_list, alloy_clean.size(), recursive);
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

  void EntryLoader::loadAlloySearchRR(const std::vector <string> & alloy_list, uint lib_max, bool recursive) {
    std::vector<std::string> rest_api_queries;
    std::vector<std::string> icsd_search_path_list;
    std::vector<std::string> libx_search_path_list;
    std::vector<std::string> libx_crawl_list;
    for (const string & sym : m_icsd_symmetries) {
      icsd_search_path_list.emplace_back(m_restapi_server + m_restapi_path + "ICSD" + m_restapi_collection + sym);
    }
    if (recursive) {
      for (uint lib_idx = 1; lib_idx <= lib_max; lib_idx++) {
        libx_search_path_list.emplace_back(m_restapi_server + m_restapi_path + "LIB" + std::to_string(lib_idx) + m_restapi_collection);
      }
    } else {
      libx_search_path_list.emplace_back(m_restapi_server + m_restapi_path + "LIB" + std::to_string(lib_max) + m_restapi_collection);
    }
    vector<string> listing;

    for (const string & url : icsd_search_path_list){
      listRestAPI(url, listing);
      cout << url << endl;
      for (const string & name: listing ){
        if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(name, 'I')) != alloy_list.end()){
          rest_api_queries.emplace_back(url + name + "/" + m_restapi_directives);
        }
      }
    }

    for (const string & url : libx_search_path_list){
      listRestAPI(url, listing);
      cout << url << endl;
      for (const string & name: listing ){
        if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(name, 'L')) != alloy_list.end()){
          cout << " " << name << " | " << extractAlloy(name, 'L') << endl;
          libx_crawl_list.emplace_back(url + name + "/" );
        }
      }
    }

    uint scanned=0;
    while(!libx_crawl_list.empty()){
      scanned += 1;
      string url = libx_crawl_list.back();
      libx_crawl_list.pop_back();
      if (scanned%100 == 0) cout << scanned << " | " << libx_crawl_list.size() << " | " << url << endl;
      listRestAPI(url, listing);
      if (listing.empty()) {
        rest_api_queries.emplace_back(url + m_restapi_directives);
        cout << " " << url << endl;
      } else {
        for (const string & name: listing ) {
          libx_crawl_list.emplace_back(url + name + "/");
        }
      }
    }


    loadRestAPIQueries(rest_api_queries, true);
    if (m_sqlite_alloy_db_ptr != nullptr) {
      vector <string> known_AUID_list;
      getAlloyAUIDList(alloy_list, known_AUID_list);
      vector <string> missing_AUID;
      for (const string &AUID: known_AUID_list) {
        if (m_auid_list.find(AUID) == m_auid_list.end()) missing_AUID.emplace_back(AUID);
      }
      std::cout << "Missed: " << missing_AUID.size() << std::endl;
      loadAUID(missing_AUID);
    }

  }


  void EntryLoader::loadAlloySearchFSR(const std::vector <string> & alloy_list, uint lib_max, bool recursive) {
    std::vector<std::string> found_entries;

    std::vector<std::string> search_path_list;
    search_path_list.emplace_back(m_filesystem_path+"ICSD/"+m_filesystem_collection);
    if (recursive) {
      for (uint lib_idx = 1; lib_idx <= lib_max; lib_idx++) {
        search_path_list.emplace_back(m_filesystem_path + "LIB" + std::to_string(lib_idx) + "/" + m_filesystem_collection);
      }
    } else {
      search_path_list.emplace_back(m_filesystem_path + "LIB" + std::to_string(lib_max) + "/" + m_filesystem_collection);
    }

    char * paths[search_path_list.size()+1];
    for (size_t list_idx=0; list_idx<search_path_list.size(); list_idx++){
      paths[list_idx] = (char*) search_path_list[list_idx].c_str();
    }
    paths[search_path_list.size()+1] = nullptr;

    size_t scanned = 0;
    size_t found = 0;
    size_t check_idx = m_filesystem_path.size();
    struct stat file_stats{};

    // initialize a file tree scan
    // https://man7.org/linux/man-pages/man3/fts.3.html
    // FTS_PHYSICAL - don't follow symlinks
    // FTS_NOCHDIR - don't change the workdir of the program
    // FTS_XDEV - don't descending into folders that are on a different device
    FTS *tree = fts_open(paths, FTS_PHYSICAL | FTS_NOCHDIR  | FTS_XDEV, nullptr);
    cout << "Loaded Tree" << endl;
    FTSENT *node;
    if(tree == nullptr)
    {
      //TODO error handling
      cout << "Could not create search tree";
      return;
    }

    // Iterate over all entries found in the folders listed in paths
    while((node = fts_read(tree))){
      scanned += 1;
//      cout << node->fts_level << " | " << node->fts_path  << endl;
      if ((node->fts_info & FTS_D)){
        if (node->fts_path[check_idx]=='I') {
          if (node->fts_level == 2) {
            if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(node->fts_name, 'I')) == alloy_list.end()) {
              fts_set(tree, node, FTS_SKIP);
              continue;
            } else {
              string base_path = node->fts_path;
              string full_path = base_path + "/aflowlib.out";
              if (stat(full_path.c_str(), &file_stats) == 0) {
                found += 1;
                found_entries.emplace_back(full_path);
                if (base_path.find("POCC") == std::string::npos) fts_set(tree, node, FTS_SKIP);
              }
            }
          }
        } else if (node->fts_path[check_idx]=='L'){
          if (node->fts_level == 1) {
            if (find(alloy_list.begin(), alloy_list.end(), extractAlloy(node->fts_name, 'L')) == alloy_list.end()) {
              fts_set(tree, node, FTS_SKIP);
              continue;
            } else {
              string base_path = node->fts_path;
              string full_path = base_path + "/aflowlib.out";
              if (stat(full_path.c_str(), &file_stats) == 0) {
                found += 1;
                found_entries.emplace_back(full_path);
                if (base_path.find("POCC") == std::string::npos) fts_set(tree, node, FTS_SKIP);
              }
            }
          }
          else if (node->fts_level > 1){
            string base_path = node->fts_path;
            string full_path = base_path + "/aflowlib.out";
            if (stat(full_path.c_str(), &file_stats) == 0) {
              found += 1;
              found_entries.emplace_back(full_path);
              if (base_path.find("POCC") == std::string::npos) fts_set(tree, node, FTS_SKIP);
            }
          }
        }
      }
      if (scanned%1000==0) {
        std::cout << scanned << " | " << found << " | " << node->fts_path << std::endl;
      }
    }
    std::cout << "Scanned: " << scanned << " | Found: " << found << std::endl;
    loadFiles(found_entries);
    if (m_sqlite_alloy_db_ptr != nullptr) {
      vector <string> known_AUID_list;
      getAlloyAUIDList(alloy_list, known_AUID_list);
      vector <string> missing_AUID;
      for (const string &AUID: known_AUID_list) {
        if (m_auid_list.find(AUID) == m_auid_list.end()) missing_AUID.emplace_back(AUID);
      }
      std::cout << "Missed: " << missing_AUID.size() << std::endl;
      loadAUID(missing_AUID);
    }
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
