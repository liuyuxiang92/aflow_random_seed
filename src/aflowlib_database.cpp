//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
//****************************************************************************

// Class for the AFLOW database, based on SQLite. It provides the framework
// for rebuilding and interacting with the database file using data in JSON
// format. This file is divided into the following sections:
//
// * Constructor
// * Temp file handling
// * Database rebuilder
// * Database analyzer
// * User-level functions to interact with the database
// * SQLite interface

#include "aflowlib.h"
#include "SQLITE/sqlite3.h"

// Some parts are written within the C++0x support in GCC, especially std::thread,
// which is implemented in gcc 4.4 and higher. For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400
#define AFLOW_DB_MULTITHREADS_ENABLE
#include <thread>
#include <mutex>
static std::mutex m;
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

#define _AFLOW_DB_DEBUG_     false
#define _SQL_COMMAND_DEBUG_  false  // debug SQL commands that are sent - verbose output
#define _SQL_CALLBACK_DEBUG_ false  // debug SQL callbacks - verbose output

using std::string;
using std::vector;

static const string _AFLOW_DB_ERR_PREFIX_ = "AflowDB::";
static const int _DEFAULT_SET_LIMIT_ = 16;
static const int _MAX_CPU_ = 16;

/************************** CONSTRUCTOR/DESTRUCTOR **************************/

namespace aflowlib {

// Open the database for read access
AflowDB::AflowDB(string db_file) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  free();
  database_file = db_file;
  if (LDEBUG) std::cerr << "AflowDB: reading database" << std::endl;
  open();
}

// Open the database for write access
AflowDB::AflowDB(string db_file, string dt_path, string schm_file, bool use_tmp_file) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  free();
  data_path = dt_path;
  database_file = db_file;
  schema_file = schm_file;
  use_tmp = use_tmp_file;
  if (LDEBUG) {
    std::cerr << "AflowDB: Database file: " << database_file << std::endl;
    std::cerr << "AflowDB: Data path: " << data_path << std::endl;
    std::cerr << "AflowDB: Schema file: " << schema_file << std::endl;
    std::cerr << "AflowDB: Use tmp file: " << use_tmp << std::endl;
  }
  open();
}

AflowDB::~AflowDB() {
  free();
  close();
}

void AflowDB::free() {
  data_path = "";
  database_file = "";
  schema_file = "";
  is_tmp = false;
}

// Opens the database file and creates a cursor
void AflowDB::open() {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "open(): Opening " << database_file << std::endl;
  int sql_code = sqlite3_open(database_file.c_str(), &db);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "open()";
    string message = "Could not open database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
}

// Closes the database file and removes the cursor
void AflowDB::close() {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  if (isTMP()) closeTmpFile(false);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "close() Closing " << database_file << std::endl;
  int sql_code = sqlite3_close(db);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "close()";
    string message = "Could not close database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
}

} // namespace aflowlib

/********************************* TMP FILE **********************************/

namespace aflowlib {

//openTmpFile/////////////////////////////////////////////////////////////////
// Opens a temporary database file 
void AflowDB::openTmpFile() {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "openTmpFile(): Opening " << database_file << ".tmp" << std::endl;
  int sql_code = sqlite3_close(db);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "openTmpFile()";
    string message = "Could not close main database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }

  string tmp_path = database_file + ".tmp";
  sql_code = sqlite3_open(tmp_path.c_str(), &db);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "openTmpFile()";
    string message = "Could not open tmp database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
  is_tmp = true;
}

//closeTmpFile////////////////////////////////////////////////////////////////
// Closes a database file and overwrites the original if nocopy. In any case,
// the temporary file is removed upon closing.
bool AflowDB::closeTmpFile(bool nocopy) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "closeTmpFile(): Closing " << database_file << ".tmp" << std::endl;
  int sql_code = sqlite3_close(db);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "openTmpFile()";
    string message = "Could not close tmp database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }

  bool copied;
  if (nocopy) {
    copied = false;
  } else {
    copied = aurostd::CopyFile(database_file + ".tmp", database_file);
  }
  aurostd::RemoveFile(database_file + ".tmp");

  sql_code = sqlite3_open(database_file.c_str(), &db);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "openTmpFile()";
    string message = "Could not open main database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
  is_tmp = false;
  return copied;
}

//isTMP///////////////////////////////////////////////////////////////////////
// Returns whether or not the cursor points to a temporary file or not
bool AflowDB::isTMP() {
  return is_tmp;
}

}  // namespace aflowlib

/***************************** REBUILD DATABASE *****************************/

namespace aflowlib {

//rebuildDatabase/////////////////////////////////////////////////////////////
// This function first checks if a rebuild is necessary and then initiates
// the rebuilding functions. Returns true if the database has been
// successfully rebuilt.
bool AflowDB::rebuildDatabase(bool force_rebuild) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);

  // Check whether to rebuild the database
  bool rebuild_db, file_empty;
  // Always rebuild when the user wants the rebuild.
  rebuild_db = force_rebuild;
  // Always rebuild if the database file doesn't exist. Note: When the
  // database is opened and the database file doesn't exist, SQLite creates
  // an empty file, so check if the file is empty and not if it exists.
  if (!rebuild_db) {
    file_empty = aurostd::FileEmpty(database_file); // file_empty needed for LDEBUG
    rebuild_db = file_empty;
  }

  // Check if any relevant files are newer than the database. No need to check
  // if the database has to be rebuilt for other reasons already.
  vector<vector<string> > entries = getJsonFiles(rebuild_db);

  // Rebuild the database if needed.
  if (rebuild_db) {
    if (LDEBUG) {
      string function = _AFLOW_DB_DEBUG_ + "rebuildDatabase(): ";
      if (force_rebuild) std::cerr << function << "Rebuilding database (force_rebuild)." << std::endl;
      if (file_empty) std::cerr << function << "Rebuilding database (file not found or empty)." << std::endl;
    }

    if (use_tmp) openTmpFile();
    // Drop all old tables before rebuilding
    vector<string> tables = getTables();
    for (uint t = 0; t < tables.size(); t++) {
      dropTable(tables[t]);
    }
    createSchemaTable(schema_file);
    createDataEntries(entries);
    if (use_tmp) return closeTmpFile();
    return true;
  } else {
    return false;
  }
}

//getJsonFiles////////////////////////////////////////////////////////////////
// Crawls through the AUID directory space and fetches all AUIDs that contain
// a JSON file.
vector<vector<string> > AflowDB::getJsonFiles(bool& rebuild_db) {
  bool LDEBUG = (TRUE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  vector<string> paths;
  aurostd::DirectoryLS(data_path, paths);
  if (paths.size() != 256) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "getJsonFiles()";
    string message = "Not a valid AUID directory.";
    throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
  }

  vector<vector<string> > entries(256);

#ifdef AFLOW_DB_MULTITHREADS_ENABLE
  int ncpus = init::GetCPUCores();
  if (ncpus < 1) ncpus = 1;
  if (ncpus > _MAX_CPU_) ncpus = _MAX_CPU_;
  vector<vector<int> > thread_dist = getThreadDistribution(32, ncpus);
  vector<std::thread*> threads;
  for (int t = 0; t < ncpus; t++) {
    threads.push_back(new std::thread(&AflowDB::getJsonFilesThread, this, thread_dist[t][0], thread_dist[t][1],
                                      std::ref(entries), std::ref(rebuild_db), std::ref(paths)));
  }
  for (int t = 0; t < ncpus; t++) {
    threads[t]->join();
    delete threads[t];
  }
#else
  getJsonFilesThread(0, 256, entries, rebuild_db, paths);
#endif
  if (LDEBUG) {
    std::cerr << _AFLOW_DB_ERR_PREFIX_ << "getJsonFile(): Number of entries:" << std::endl;
    for (int i = 0; i < 32; i++) {
      std::cerr << "    " << std::hex << std::setfill('0') << std::setw(2) << i
                << ": " << std::dec << entries[i].size() << std::endl;
    }
  }
  return entries;
}

void AflowDB::getJsonFilesThread(int startIndex, int endIndex, vector<vector<string> >& entries,
                                 bool& rebuild_db, const vector<string>& paths) {
  bool LDEBUG = (TRUE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  long int tm_db = aurostd::FileModificationTime(database_file);
  for (int i = startIndex; i < endIndex; i++) {
    // Fetch AUIDs
    crawl(entries[i], data_path + paths[i], 1, 8);

    // Check if there is a newer aflowlib.json (if necessary)
    if (!rebuild_db) {
      uint e = 0;
      uint nentries = entries[i].size();
      string filename;
      for (e = 0; !(rebuild_db) && (e < nentries); e++) {
        filename = data_path + "aflow:" + entries[i][e] + "/aflowlib.json";
        if (aurostd::FileModificationTime(filename) > tm_db) break;
      }
      if ((e < nentries) && (!rebuild_db)) {
      #ifdef AFLOW_DB_MULTITHREADS_ENABLE
        std::unique_lock<std::mutex> lk(m);
      #endif
        rebuild_db = true;
        if (LDEBUG) {
          string entry = aurostd::RemoveSubString(entries[i][e], "/");
          std::cerr << _AFLOW_DB_ERR_PREFIX_ << "getJsonFiles(): Found newer entry."
                    << " AUID = aflow:" << aurostd::toupper(entry) << std::endl;
        }
      }
    }
  }
}

//crawl///////////////////////////////////////////////////////////////////////
// Recursive function to crawl through the AUID directory space.
void AflowDB::crawl(vector<string>& entries, const string& parent_path, int depth, int maxdepth) {
  if (depth == maxdepth) { // Found the end of the directory tree
    string filename = parent_path + "/aflowlib.json";
    if (aurostd::FileExist(filename)) {
      entries.push_back(parent_path.substr(string(data_path + "aflow:").size(), string::npos));
    }
  } else {
    vector<string> child_paths;
    aurostd::DirectoryLS(parent_path, child_paths);
    for (uint p = 0, npaths = child_paths.size(); p < npaths; p++) {
      crawl(entries, parent_path + "/" + child_paths[p], depth + 1, maxdepth);
    }
  }
}

// Schema --------------------------------------------------------------------

//createSchemaTable///////////////////////////////////////////////////////////
// This function creates the table containing the necessary metadata from the
// schema. This table should only contain the metadata necessary to create
// the database entries and the data for the entry page.
void AflowDB::createSchemaTable(string schema_file) {
  if (aurostd::FileExist(schema_file)) {
    // Pull schema
    string schema = aurostd::file2string(schema_file);
    // To prevent SQL syntax errors, single quotes/apostrophes must be escaped
    // with an additional single quote.
    aurostd::StringSubst(schema, "'", "''");
    vector<string> cols;
    aurostd::string2tokens("name,function,status,type", cols, ",");
    vector<string> rows = getSchemaRows(schema);

    // Create table
    transaction(true);
    createTable(_DB_SCHEMA_TABLE_, cols, "TEXT");

    // Insert values
    for (uint r = 0; r < rows.size(); r++) {
      vector<string> values = getSchemaValues(rows[r], cols, schema);
      insertValues(_DB_SCHEMA_TABLE_, values);
    }
    createIndex("index_schema_name", _DB_SCHEMA_TABLE_, "name");
    transaction(false);
  } else {
    string function = _AFLOW_DB_ERR_PREFIX_ + "createSchemaTable()";
    string message = "Could not find file " + schema_file + ".";
    throw aurostd::xerror(function, message, _FILE_NOT_FOUND_);
  }
}

//getSchemaRows///////////////////////////////////////////////////////////////
// Gets the rows of the schema table, i.e. the keys of the schema file.
vector<string> AflowDB::getSchemaRows(const string& schema) {
  string command = "SELECT key FROM json_each('" + schema + "', '$.AAPI_schema')";
  command += " WHERE key IS NOT '__schema^2__'";
  return SQLexecuteCommandVECTOR(command);
}

//getSchemaValues/////////////////////////////////////////////////////////////
// Retrieves the desired metadata from the schema JSON file.
vector<string> AflowDB::getSchemaValues(string row, const vector<string>& cols, const string& schema) {
  string key, value;
  uint ncols = cols.size();
  vector<string> values(ncols, "NULL");
  for (uint c = 0; c < ncols; c++) {
    if (cols[c] == "name") {
      value = row;
    } else {
      key = "AAPI_schema." + row + "." + cols[c];
      value = extractJSONvalue(schema, key);
    }
    if (!value.empty()) values[c] = "'" + value + "'";
  }
  return values;
}

// Data ----------------------------------------------------------------------

//createDataEntries///////////////////////////////////////////////////////////
// This function inserts the data into the database.
void AflowDB::createDataEntries(const vector<vector<string> >& entries) {
  // Do not constantly synchronize the database with file on disk.
  // Increases performance significantly.
  SQLexecuteCommand("PRAGMA synchronous = OFF");

  // Get all relevant properties from the schema
  vector<string> cols, types;
  transaction(true);
  // Links and images are static, so storing them in the database would be wasteful
  string where = "function IS NOT 'link' AND type IS NOT 'image'";
  where += " AND name IS NOT 'icsd_number'";  // ICSD numbers are not stored in the database
  cols = getSet(_DB_SCHEMA_TABLE_, "name", false, where);
  types = getDataTypes(cols);

  // Insert data
#ifdef AFLOW_DB_ENABLE_MULTITHREADS
  int ncpus = init::GetCPUCores();
  if (ncpus < 1) ncpus = 1;
  if (ncpus > _MAX_CPU_) ncpus = _MAX_CPU_;
  vector<vector<int> > thread_dist = getThreadDistribution(256, ncpus);
  vector<std::thread*> threads;
  for (int t = 0; t < ncpus; t++) {
    threads.push_back(new std::thread(&AflowDB::populateTables, this, thread_dist[t][0], thread_dist[t][1],
                                      std::ref(entries), std::ref(cols), std::ref(types)));
  }
  for (int t = 0; t < ncpus; t++) {
    threads[t]->join();
    delete threads[t];
  }
#else
  populateTables(0, 256, entries, cols, types);
#endif

  // Create indices on the properties that will be used to retrieve database entries
  vector<string> tables, index_cols;
  tables = getTableSubset(_DB_SCHEMA_TABLE_);
  aurostd::string2tokens("auid,aurl,prototype,title,catalog", index_cols, ",");

  string index;
  for (uint t = 0; t < tables.size(); t++) {
    for (uint i = 0; i < index_cols.size(); i++) {
      index = "index_" + tables[t] + "_" + index_cols[i];
      createIndex(index, tables[t], index_cols[i]);
    }
  }
  transaction(false);
}

void AflowDB::populateTables(int startIndex, int endIndex, const vector<vector<string> >& entries,
                             const vector<string>& cols, const vector<string>& types) {
  string table, json, filename;
  vector<string> vals;
  int chunk_size = 1000;  // To insert entries in smaller chunks
  for (int i = startIndex; i < endIndex; i++) {
    uint nentries = entries[i].size();
    transaction(true);
    int counter = 0;
    stringstream t;
    t << std::setfill('0') << std::setw(2) << std::hex << i;
    table = t.str();
    createTable(table, cols, types);
    for (uint e = 0; e < nentries; e++) {
      json = aurostd::file2string(data_path + "aflow:" + entries[i][e] + "/aflowlib.json");
      vals = getDataValues(json, cols, types);
      insertValues(table, cols, vals);
      if (++counter % chunk_size == 0) {
        transaction(false);
        transaction(true);
      }
    }
    transaction(false);
  }
}

//getDataTypes////////////////////////////////////////////////////////////////
// Gets the data types of the schema keys and converts them into SQLite types.
// Note that SQLite does not recognize arrays, so they will be stored as text.
vector<string> AflowDB::getDataTypes(const vector<string>& cols) {
  string where;
  vector<string> types(cols.size());
  for (uint i = 0; i < cols.size(); i++) {
    where = "name='" + cols[i] + "'";
    types[i] = getValue(_DB_SCHEMA_TABLE_, "type", where);
    // AUID, AURL, and title have to be unique
    if ((cols[i] == "auid") || (cols[i] == "aurl") || (cols[i] == "title")) {
      types[i] = "TEXT UNIQUE NOT NULL";
    } else if (types[i] == "number") {
      types[i] = "REAL";
    } else {
      types[i] = "TEXT";
    }
  }
  return types;
}

//getDataValues///////////////////////////////////////////////////////////////
// Retrieves the values of each property from the aflowlib.json file.
vector<string> AflowDB::getDataValues(const string& entry, const vector<string>& cols, const vector<string>& types) {
  string value, id;
  uint ncols = cols.size();
  vector<string> values(ncols, "NULL");
  for (uint c = 0; c < ncols; c++) {
    value = extractJSONvalueAFLOW(entry, cols[c]);
    if (!value.empty()) {
      if (types[c] != "REAL") values[c] = "'" + value + "'";
      else values[c] = value;
    }
  }
  return values;
}

} // namespace aflowlib

/**************************** DATABASE ANALYTICS ****************************/

namespace aflowlib {

//analyzeDatabase/////////////////////////////////////////////////////////////
// Provides analytics for each database data table in JSON format.
void AflowDB::analyzeDatabase(string outfile) {
  std::stringstream json;
  string tab = "    ";
  json << "{" << std::endl << tab << "\"Aflow_DBs\": {" << std::endl;

  vector<string> tables = getTableSubset(_DB_SCHEMA_TABLE_);
  uint ntable = tables.size();
  for (uint t = 0; t < ntable; t++) {
    std::cout << "Starting analytics for table " << tables[t] << std::endl;
    DBstats db_stats = getColStats(tables[t]);
    if (db_stats.nentries > 0) {
      db_stats.loop_counts = getDBLoopCounts(tables[t]);
      vector<string> species = getSet(tables[t], "species", true);
      db_stats.nsystems = species.size();
      db_stats.species = getUniqueFromJSONArrays(species);
    }
    writeStatsToJSON(json, db_stats);
    if (t < ntable - 1) json << ",";
    json << std::endl;
  }
  json << tab << "}" << std::endl << "}" << std::endl;
  aurostd::stringstream2file(json, outfile);
}

//getColStats/////////////////////////////////////////////////////////////////
// Gets the number of non-null entries, max, min, and the set (if below a
// certain threshold) of each property in the schema.
DBstats AflowDB::getColStats(const string& table) {
  DBstats stats;
  stats.nentries = aurostd::string2utype<int>(getProperty("COUNT", table, "*"));
  stats.table = table;
  string where = "function IS NOT 'link' AND type IS NOT 'image'";
  where += " AND name IS NOT 'icsd_number'";  // ICSD numbers are not stored in the database
  stats.columns = getSet(_DB_SCHEMA_TABLE_, "name", false, where);
  uint ncols = stats.columns.size();
  stats.count.assign(ncols, 0);
  stats.max.assign(ncols, "");
  stats.min.assign(ncols, "");
  stats.set.resize(ncols);
  if (stats.nentries > 0) {
    for (uint c = 0; c < ncols; c++) {
      stats.count[c] = aurostd::string2utype<int>(getProperty("COUNT", table, stats.columns[c]));
      if (stats.count[c] > 0) {  // No need to determine properties if count is zero
        stats.max[c] = getProperty("MAX", table, stats.columns[c]);
        stats.min[c] = getProperty("MIN", table, stats.columns[c]);
        // Do not sort set yet - sorting is expensive
        string where = stats.columns[c] + " NOT NULL";
        stats.set[c] = getSet(table, stats.columns[c], true, where, _DEFAULT_SET_LIMIT_ + 1);
      }
    }
  }
  return stats;
}

//getDBLoopCounts/////////////////////////////////////////////////////////////
// Gets the number of times each loop type appears in the database. This is a
// proxy for the number of different property types that are stored in the
// database (e.g. elastic properties are stored as "agl" in the loop array).
vector<std::pair<string, int> > AflowDB::getDBLoopCounts(const string& table) {
  // First, get unique loop array elements
  vector<string> loop_entries, loops;
  loop_entries = getSet(table, "loop", true);
  loops = getUniqueFromJSONArrays(loop_entries);
  uint nloops = loops.size();

  vector<std::pair<string, int> > loop_counts(nloops);
  string where, count;
  string where_table = "catalog='" + table + "' AND ";
  for (uint l = 0; l < nloops; l++) {
    where = "loop LIKE '%" + loops[l] + "%'";
    count = getProperty("COUNT", table, "loop", where);
    loop_counts[l].first = loops[l];
    loop_counts[l].second = aurostd::string2utype<int>(count);
  }
  return loop_counts;
}

//getUniqueFromJSONArrays/////////////////////////////////////////////////////
// Determines the unique array elements in a set of 1D-array strings.
vector<string> AflowDB::getUniqueFromJSONArrays(const vector<string>& arrays) {
  vector<string> unique, tokens;
  string arr;
  for (uint a = 0; a < arrays.size(); a++) {
    vector<string> tokens;
    arr = aurostd::RemoveSubString(arr, "[");
    arr = aurostd::RemoveSubString(arr, "]");
    arr = aurostd::RemoveSubString(arr, "\"");
    aurostd::string2tokens(arr, tokens, ", ");
    for (uint t = 0; t < tokens.size(); t++) {
      if (unique.size() == 0) {
        unique.push_back(tokens[t]);
      } else {
        bool append = true;
        for (uint u = 0; u < unique.size() && append; u++) {
          if (unique[u] == tokens[t]) append = false;
        }
        if (append) unique.push_back(tokens[t]);
      }
    }
  }
  return unique;
}

//writeStatsToJSON////////////////////////////////////////////////////////////
// Writes the database statistics into a JSON-formatted string(stream).
void AflowDB::writeStatsToJSON(std::stringstream& json, const DBstats& db_stats) {
  string tab = "    ";
  string indent = tab + tab;
  json << indent << "\"" << db_stats.table << "\": {" << std::endl;
  json << indent << tab << "\"count\": " << db_stats.nentries << "," << std::endl;
  json << indent << tab << "\"systems\": " <<  db_stats.nsystems << "," << std::endl;
  for (uint l = 0; l < db_stats.loop_counts.size(); l++) {
    json << indent << tab << "\"" << db_stats.loop_counts[l].first << "\": ";
    json << db_stats.loop_counts[l].second << "," << std::endl;
  }
  json << indent << tab << "\"columns\": {" << std::endl;
  uint ncols = db_stats.columns.size();
  for (uint c = 0; c < ncols; c++) {
    json << indent << tab << tab << "\"" << db_stats.columns[c] << "\": {" << std::endl;
    json << indent << tab << tab << tab << "\"count\": ";
    json << aurostd::utype2string<int>(db_stats.count[c]) << "," << std::endl;
    json << indent << tab << tab << tab << "\"min\": ";
    json << (db_stats.min[c].empty()?"null":string("\"" + db_stats.min[c] + "\"")) << "," << std::endl;
    json << indent << tab << tab << tab << "\"max\": ";
    json << (db_stats.max[c].empty()?"null":string("\"" + db_stats.max[c] + "\"")) << "," << std::endl;

    // Write set
    uint nset = db_stats.set[c].size();
    json << indent << tab << tab << tab << "\"set\": ";
    if (nset > _DEFAULT_SET_LIMIT_) {
      json << "null" << std::endl;
    } else if (nset == 0) {
      json << "[]" << std::endl;
    } else {
      // Retrieve the set again, but sorted. Sorting via SQLite may not be
      // as efficient as std::sort, but std::sort does not sort the numbers
      // correctly without type hinting.
      string where = db_stats.columns[c] + " NOT NULL";
      vector<string> set = getSet(db_stats.table, db_stats.columns[c], true, where, _DEFAULT_SET_LIMIT_ + 1, db_stats.columns[c]);
      json << "[" << std::endl;
      for (uint s = 0; s < nset; s++) {
        json << indent << tab << tab << tab << tab << "\"" << set[s] << "\"";
        if (s < nset - 1) json << ",";
        json << std::endl;
      }
      json << indent << tab << tab << tab << "]" << std::endl;
    }

    json << indent << tab << tab << "}";
    if (c < ncols - 1) json << ",";
    json << std::endl;
  }
  json << indent << tab << "}," << std::endl;
  uint nspecies = db_stats.species.size();
  json << indent << tab << "\"species\": " << ((nspecies > 0)?"[":"null") << std::endl;
  for (uint s = 0; s < db_stats.species.size(); s++) {
    json << indent << tab << tab << "\"" << db_stats.species[s] << "\"";
    if (s < nspecies - 1) json << ",";
    json << std::endl;
  }
  if (nspecies > 0) json << indent << tab << "]" << std::endl;
  json << indent << "}";
}

}  // namespace aflowlib

/***************************** SQLite FUNCTIONS *****************************/

// These functions are higher level SQLite functions that call functions from
// the SQLite interface. They are essentially syntactic sugar to make the code
// easier to read and to facilitate the implementation of the AflowDB class
// into other parts of AFLOW.

namespace aflowlib {

// INDEX ---------------------------------------------------------------------

//createIndex/////////////////////////////////////////////////////////////////
// Creates an index.
void AflowDB::createIndex(const string& index, const string& table, const string& column) {
  string command = "CREATE INDEX " + index + " ON " + table + "(" + column +  ");";
  SQLexecuteCommand(command);
}

//dropIndex///////////////////////////////////////////////////////////////////
// Removes an index.
void AflowDB::dropIndex(const string& index) {
  string command = "DROP INDEX " + index + ";";
  SQLexecuteCommand(command);
}

// JSON ----------------------------------------------------------------------

//extractJSONvalue////////////////////////////////////////////////////////////
// Extracts a JSON value using SQLite's JSON extension.
string AflowDB::extractJSONvalue(const string& json, const string& key) {
  string command = "SELECT json_extract('" + json + "', '$." + key + "');";
  return SQLexecuteCommandSCALAR(command); 
}

//extractJSONvalueAFLOW///////////////////////////////////////////////////////
// This function extracts values from an aflowlib.json file. It is much faster
// than using SQLite's JSON extension, but has was designed to only work for
// the aflowlib.json. It cannot handle nested JSONs and there is no guarantee
// that it will work for any other flat JSON objects.
string AflowDB::extractJSONvalueAFLOW(const string& json, string key) {
  string value;
  string keyback = key;
  key = "\"" + key + "\":";
  string::size_type start, end;
  start = json.find(key);
  if (start != string::npos) {
    start += key.length();
    end = json.find("\":", start);
    if (end != string::npos) {
      value = json.substr(start, end - start);
      end = value.find_last_of(",\"") - 1;
      value = value.substr(0, end);
    } else {
      end = json.find("}", start);
      value = json.substr(start, end - start - 1);
    }
    if (value[0] == '"') {
      value = value.substr(1, end - 2);
    } else {
      value = value.substr(0, end);
    }
  } else {
    value = "";
  }
  return value;
}

// TRANSACTION ---------------------------------------------------------------

void AflowDB::transaction(bool begin) {
  string command = string(begin?"BEGIN":"END") + " TRANSACTION;";
  SQLexecuteCommand(command);
}

// TABLE ---------------------------------------------------------------------

//dropTable///////////////////////////////////////////////////////////////////
// Deletes a table from the database.
void AflowDB::dropTable(const string& table) {
  string command = "DROP TABLE IF EXISTS " + table + ";";
  SQLexecuteCommand(command);
}

//getTableSubset//////////////////////////////////////////////////////////////
// Gets all tables except for the table provided in exclude.
vector<string> AflowDB::getTableSubset(const string& exclude) {
  string where = "name IS NOT '" + exclude + "'";
  return getTables(where);
}

// Gets all tables except for the tables provided in exclude.
vector<string> AflowDB::getTableSubset(const vector<string>& exclude) {
  uint nexclude = exclude.size();
  string where;
  for (uint i = 0; i < nexclude; i++) {
    where += "name IS NOT '" + exclude[i] + "'";
    if (i < nexclude - 1) where += " AND ";
  }
  return getTables(where);
}

//getTables///////////////////////////////////////////////////////////////////
// Retrieves a set of tables. If where is empty, all tables in the database
// will be returned.
vector<string> AflowDB::getTables(string where) {
  string command = "SELECT name FROM sqlite_master WHERE type='table'";
  if (!where.empty()) command += " AND (" + where + ")";
  command += ";";
  return SQLexecuteCommandVECTOR(command);
}

//createTable/////////////////////////////////////////////////////////////////
// Creates a table where all columns have the same type.
void AflowDB::createTable(const string& table, const vector<string>& cols, const string& type) {
  vector<string> types(cols.size(), type);
  createTable(table, cols, types);
}

// Creates a table where each column is assigned its own type
void AflowDB::createTable(const string& table, const vector<string>& cols, const vector<string>& types, bool temp) {
  uint ncols = cols.size();
  if (ncols != types.size()) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "createTable()";
    string message = "Could not create table. ";
    message += "Number of columns and number of types do not match.";
    throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
  } else {
    string command = "CREATE " + string(temp?"TEMP ":"") + "TABLE IF NOT EXISTS " + table + "(";
    for (uint c = 0; c < ncols; c++) {
      command += cols[c] + " " + types[c];
      if (c < ncols - 1) command += ", ";
    }
    command += ");";
    SQLexecuteCommand(command);
  }
}

//createTableAs///////////////////////////////////////////////////////////////
// Creates a table using the AS option.
void AflowDB::createTableAs(const string& table, const string& as, bool temp) {
  string command = "CREATE " + string(temp?"TEMP ":"") + "TABLE IF NOT EXISTS " + table + " AS " + as + ";";
  SQLexecuteCommand(command);
}

//createTempTable/////////////////////////////////////////////////////////////
// Creates a temporary table that will be destroyed when the database cursor
// is closed.
void AflowDB::createTempTable(const string& table, const vector<string>& cols, const string& type) {
  vector<string> types(cols.size(), type);
  createTable(table, cols, types, true);
}

void AflowDB::createTempTable(const string& table, const vector<string>& cols, const vector<string>& types) {
  createTable(table, cols, types, true);
}

//createTableAs///////////////////////////////////////////////////////////////
// Creates a temporary table using the AS option.
void AflowDB::createTempTableAs(const string& table, const string& as) {
  createTableAs(table, as, true);
}

//dropTempTable///////////////////////////////////////////////////////////////
// Removes a temporary table.
void AflowDB::dropTempTable(const string& table) {
  // Make sure that the table that is dropped is a temporary table
  string command = "SELECT name FROM sqlite_temp_master WHERE type='table'";
  vector<string> temp_tables = SQLexecuteCommandVECTOR(command);
  if (aurostd::withinList(temp_tables, table)) {
    dropTable(table);
  } else {
    string function = "AflowDB::dropTempTable()";
    string message = "Could not find temporary table " + table + ".";
    throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
  }
}

// INSERT --------------------------------------------------------------------

//insertValues////////////////////////////////////////////////////////////////
// Inserts a set of values into a table.
void AflowDB::insertValues(const string& table, const vector<string>& vals) {
  vector<string> cols;
  insertValues(table, cols, vals);
}

// Inserts a set of values into a table (if cols is empty) or into specific
// columns of the table.
void AflowDB::insertValues(const string& table, const vector<string>& cols, const vector<string>& vals) {
  uint ncols = cols.size();
  uint nvals = vals.size();
  if ((ncols > 0) && (ncols != nvals)) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "insertValues()";
    string message = "Could not insert values. ";
    message += "Number of columns and number of values do not match.";
    throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
  } else {
    string command = "INSERT INTO " + table;
    if (ncols > 0) {
      command += " (";
      for (uint c = 0; c < ncols; c++) {
        command += cols[c];
        if (c < ncols - 1) command += ", ";
      }
      command += ")";
    }
    command += " VALUES(";
    for (uint v = 0; v < nvals; v++) {
      command += vals[v];
      if (v < nvals - 1) command += ", ";
    }
    command += ");";
    SQLexecuteCommand(command);
  }
}

// GET -----------------------------------------------------------------------

//getValue////////////////////////////////////////////////////////////////////
// Gets a value from a specific column. The row must be specified in the where
// condition, or else it just takes the first value.
string AflowDB::getValue(const string& table, const string& col, string where) {
  string command = prepareSELECT(table, "", col, where, 0, "");
  return SQLexecuteCommandSCALAR(command);
}

//getValueMultiCol////////////////////////////////////////////////////////////
// Gets a value from multiple columns. The row must be specified in the where
// condition, or else it just takes the first value.
vector<string> AflowDB::getValueMultiCol(const string& table, const vector<string>& cols, string where) {
  string command = prepareSELECT(table, "", cols, where, 0, "");
  vector<vector<string> > values = SQLexecuteCommand2DVECTOR(command);
  return values[0];
}

//getProperty/////////////////////////////////////////////////////////////////
// Gets a database property for a specific column.
string AflowDB::getProperty(const string& property, const string& table,
                            const string& col, string where) {
  string command = prepareSELECT(table, property, col, where, 0, "");
  return SQLexecuteCommandSCALAR(command);
}

//getColumn///////////////////////////////////////////////////////////////////
// Returns an entire column. The output can be sorted by another column using
// the "order_by" string.
vector<string> AflowDB::getColumn(const string& table, const string& col, string order_by) {
  return getSet(table, col, false, "", 0, order_by);
}

//getSet//////////////////////////////////////////////////////////////////////
// Retrieves a (distinct) set from a single column.
vector<string> AflowDB::getSet(const string& table, const string& col, bool distinct,
                               string where, int limit, string order_by) {
  string property = string((distinct?"DISTINCT":""));
  string command = prepareSELECT(table, property, col, where, limit, order_by);
  return SQLexecuteCommandVECTOR(command);
}

//getSetMultiCol//////////////////////////////////////////////////////////////
// Retrieves a (distinct) set from a multiple columns.
vector<vector<string> > AflowDB::getSetMultiCol(const string& table, const vector<string>& cols, bool distinct,
                                                string where, int limit, string order_by) {
  string property = string(distinct?"DISTINCT":"");
  string command = prepareSELECT(table, property, aurostd::joinWDelimiter(cols, ", "), where, limit, order_by);
  return SQLexecuteCommand2DVECTOR(command);
}

//prepateSELECT///////////////////////////////////////////////////////////////
// Lower level function to prepare a SELECT statement for all GET functions.
string AflowDB::prepareSELECT(const string& table, const string& property, const string& cols,
                              string where, int limit, string order_by) {
  stringstream command;
  command << "SELECT ";
  if (!property.empty()) command << property << ((property == "DISTINCT")?" ":"(");
  command << cols;
  if (!property.empty()) command << ((property == "DISTINCT")?"":")");
  command << " FROM " << table;
  if (!where.empty()) command << " WHERE (" << where << ")";
  if (!order_by.empty()) command << " ORDER BY " << order_by;
  if (limit > 0) command << " LIMIT " << limit;
  command << ";";
  return command.str();
}

// Wrapper function for a vector representation of the columns.
string AflowDB::prepareSELECT(const string& table, const string& property, const vector<string>& cols,
                              string where, int limit, string order_by) {
  return prepareSELECT(table, property, aurostd::joinWDelimiter(cols, ", "), where, limit, order_by);
}

}  // namespace aflowlib

/***************************** SQLite INTERFACE *****************************/

// These functions provide a direct interface to SQLite. They execute commands
// and handle callbacks.

namespace aflowlib {

// Execute command functions -------------------------------------------------
// These functions provide a framework to execute SQLite commands (including
// exception handling). The return types should all be void or string-typed to
// keep the number of functions to a minimum. Conversion to other data types
// should be handled outside these functions.

void AflowDB::SQLexecuteCommand(const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << "AflowDB::SQLexecuteCommand(): command = " << command << std::endl;
  char* sqlErrMsg = 0;
  int sql_code = sqlite3_exec(db, command.c_str(), SQLcallback, 0, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "SQLexecuteCommand()";
    string message = string(sqlErrMsg) + " in command " + command;
    throw aurostd::xerror(function, message, _RUNTIME_SQL_);
  }
}

string AflowDB::SQLexecuteCommandSCALAR(const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << "AflowDB::SQLexecuteCommandSCALAR(): command = " << command << std::endl;
  char* sqlErrMsg = 0;
  string returnstring;
  int sql_code = sqlite3_exec(db, command.c_str(), SQLcallbackSCALAR, &returnstring, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "SQLexecuteCommandSCALAR()";
    string message = string(sqlErrMsg) + " in command " + command;
    throw aurostd::xerror(function, message, _RUNTIME_SQL_);
  } else {
    return returnstring;
  }
}

vector<string> AflowDB::SQLexecuteCommandVECTOR(const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "SQLexecuteCommandVECTOR(): command = " << command << std::endl;
  char *sqlErrMsg = 0;
  vector<string> returnvector;
  int sql_code = sqlite3_exec(db, command.c_str(), SQLcallbackVECTOR, &returnvector, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "SQLexecuteCommandVECTOR()";
    string message = string(sqlErrMsg) + " in command " + command;
    throw aurostd::xerror(function, message, _RUNTIME_SQL_);
  } else {
    return returnvector;
  }
}

vector<vector<string> > AflowDB::SQLexecuteCommand2DVECTOR(const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "SQLexecuteCommand2DVECTOR(): command = " << command << std::endl;
  char *sqlErrMsg = 0;
  vector<vector<string> > returnvector;
  int sql_code = sqlite3_exec(db, command.c_str(), SQLcallback2DVECTOR, &returnvector, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "SQLexecuteCommand2DVECTOR()";
    string message = string(sqlErrMsg) + " in command " + command;
    throw aurostd::xerror(function, message, _RUNTIME_SQL_);
  } else {
    return returnvector;
  }
}

// Callback functions --------------------------------------------------------
// The following functions are the callback functions passed into sqlite_exec.
// Each executeCommand function should have its own callback function.

int AflowDB::SQLcallback(void* data, int argc, char** argv, char** col) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_CALLBACK_DEBUG_);
  (void) data;  // To suppress compiler warnings
  if (LDEBUG) {
    for (int i = 0; i < argc; i++) {
      std::cerr << _AFLOW_DB_ERR_PREFIX_ << "SQLcallback()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
    }
    std::cerr << std::endl;
  }
  return 0;
}

int AflowDB::SQLcallbackSCALAR(void* data, int argc, char** argv, char** col) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_CALLBACK_DEBUG_);
  if (LDEBUG) {
    for (int i = 0; i < argc; i++) {
      std::cerr << _AFLOW_DB_ERR_PREFIX_ << "SQLcallbackSCALAR()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
    }
    std::cerr << std::endl;
  }
  string* val = static_cast<string*>(data);
  if (argc == 1) {
    if (argv[0] == NULL) {
      *val = "";
    } else {
      *val = std::string(argv[0]);
    }
    return 0;
  } else {
    return 1;
  }
}

int AflowDB::SQLcallbackVECTOR(void* data, int argc, char** argv, char** col) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_CALLBACK_DEBUG_);
  if (LDEBUG) {
    for (int i = 0; i < argc; i++) {
      std::cerr << _AFLOW_DB_ERR_PREFIX_ << "SQLcallbackVECTOR()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
    }
    std::cerr << std::endl;
  }
  vector<string>* vec = static_cast<vector<string>*>(data);
  if (argc == 1) {
    if (argv[0] != NULL) vec->push_back(string(argv[0]));
    else vec->push_back("");
    return 0;
  } else {
    return 1;
  }
}

int AflowDB::SQLcallback2DVECTOR(void* data, int argc, char** argv, char** col) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_CALLBACK_DEBUG_);
  if (LDEBUG) {
    for (int i = 0; i < argc; i++) {
      std::cerr << _AFLOW_DB_ERR_PREFIX_ << "SQLcallback2DVECTOR()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
    }
    std::cerr << std::endl;
  }
  vector<vector<string> >* vec2d = static_cast<vector<vector<string> >*>(data);
  vector<string> vec;
  for (int i = 0; i < argc; i++) {
    if (argv[i] != NULL) vec.push_back(string(argv[i]));
    else vec.push_back("");
  }
  vec2d->push_back(vec);
  return 0;
}

}  // namespace aflowlib

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
//****************************************************************************
