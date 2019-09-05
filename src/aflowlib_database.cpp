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
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

#define _AFLOW_DB_DEBUG_     true
#define _SQL_COMMAND_DEBUG_  false  // debug SQL commands that are sent - verbose output
#define _SQL_CALLBACK_DEBUG_ false  // debug SQL callbacks - verbose output

using std::string;
using std::vector;

static const string _AFLOW_DB_ERR_PREFIX_ = "AflowDB::";
static const int _DEFAULT_SET_LIMIT_ = 16;
static const int _MAX_CPU_ = 16;
static const int _N_AUID_TABLES_ = 256;
static const int ntables_test = _N_AUID_TABLES_;
//static const int ntables_test = 2;

/************************** CONSTRUCTOR/DESTRUCTOR **************************/

namespace aflowlib {

// Open the database for read access
AflowDB::AflowDB(string db_file) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  free();
  database_file = db_file;
  if (LDEBUG) std::cerr << "AflowDB: reading database" << std::endl;
  open(SQLITE_OPEN_READONLY);
}

// Open the database for write access
AflowDB::AflowDB(string db_file, string dt_path, string schm_file) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  free();
  data_path = dt_path;
  database_file = db_file;
  schema_file = schm_file;
  if (LDEBUG) {
    std::cerr << "AflowDB: Database file: " << database_file << std::endl;
    std::cerr << "AflowDB: Data path: " << data_path << std::endl;
    std::cerr << "AflowDB: Schema file: " << schema_file << std::endl;
  }
  open();
}

AflowDB::~AflowDB() {
  close();
  free();
}

void AflowDB::free() {
  data_path = "";
  database_file = "";
  schema_file = "";
  is_tmp = false;
}

// Opens the database file and creates a cursor
void AflowDB::open(int open_flags) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "open(): Opening " << database_file << std::endl;
  int sql_code = sqlite3_open_v2(database_file.c_str(), &db, open_flags, nullptr);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "open()";
    string message = "Could not open database file " + database_file;
    message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
}

// Closes the database file and removes the cursor
void AflowDB::close() {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  if (isTMP()) closeTmpFile(true, true);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "close() Closing " << database_file << std::endl;
  int sql_code = sqlite3_close(db);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "close()";
    string message = "Could not close database file " + database_file;
    message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
}

} // namespace aflowlib

/********************************* TMP FILE **********************************/

namespace aflowlib {

//openTmpFile/////////////////////////////////////////////////////////////////
// Opens a temporary database file 
void AflowDB::openTmpFile(int open_flags) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  string tmp_path = database_file + ".tmp";
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "openTmpFile(): Opening " << tmp_path << std::endl;
  if (aurostd::FileExist(tmp_path)) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "openTmpFile()";
    string message = "Could not create temporary database file. File already exists.";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
  int sql_code = sqlite3_close(db);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "openTmpFile()";
    string message = "Could not close main database file " + database_file;
    message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }

  sql_code = sqlite3_open_v2(tmp_path.c_str(), &db, open_flags, nullptr);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "openTmpFile()";
    string message = "Could not open tmp database file " + database_file + ".tmp";
    message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
  is_tmp = true;
}

//closeTmpFile////////////////////////////////////////////////////////////////
// Closes a database file and overwrites the original unless nocopy is true.
bool AflowDB::closeTmpFile(bool nocopy, bool keep) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "closeTmpFile(): Closing " << database_file << ".tmp" << std::endl;
  int sql_code = sqlite3_close(db);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "openTmpFile()";
    string message = "Could not close tmp database file " + database_file;
    message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }

  bool copied;
  if (nocopy) {
    copied = false;
  } else {// if (aurostd::FileSize(database_file) < aurostd::FileSize(database_file + ".tmp")) {
    copied = aurostd::CopyFile(database_file + ".tmp", database_file);
  }
  if (!keep) aurostd::RemoveFile(database_file + ".tmp");

  open();
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
  string function = _AFLOW_DB_ERR_PREFIX_ + "rebuildDatabase()";  // for LDEBUG and xerror
  bool rebuild_db = false, update_db = false;

  // Always rebuild when the user wants the rebuild.
  rebuild_db = force_rebuild;
  if (LDEBUG && rebuild_db) std::cerr << function << ": Rebuilding database (user rebuild)." << std::endl;

  // Always rebuild if the database file doesn't exist. Note: When the
  // database is opened and the database file doesn't exist, SQLite creates
  // an empty file, so check if the file is empty and not if it exists.
  if (!rebuild_db) {
    rebuild_db = aurostd::FileEmpty(database_file); // file_empty needed for LDEBUG
    if (LDEBUG && rebuild_db) std::cerr << function << ": Rebuilding database (file not found or empty)." << std::endl;
  }

  // Rebuild when the schema file has new entries
  if (!rebuild_db) {
    vector<string> columns, keys;
    string table = getTables()[0];  // All tables have the same columns, so any table is good
    columns = getColumnNames(table);

    string schema = aurostd::file2string(schema_file);
    // To prevent SQL syntax errors, single quotes/apostrophes must be escaped
    // with an additional single quote.
    aurostd::StringSubst(schema, "'", "''");
    keys = getSchemaKeys(schema);

    uint nkeys = keys.size();
    rebuild_db = (nkeys != columns.size());
    if (!rebuild_db) {
      vector<string> types_db = getColumnTypes(table);
      vector<string> types_schema(nkeys);
      uint k;
      for (k = 0; k < nkeys; k++) types_schema[k] = getSchemaValues(keys[k], "type", schema)[0];
      for (k = 0; k < nkeys; k++) {
        if (!aurostd::withinList(columns, keys[k])) break;
        if (!aurostd::withinList(types_db, types_schema[k])) break;
      }
      rebuild_db = (k != nkeys);
    }
    if (LDEBUG && rebuild_db) std::cerr << function << ": Rebuilding database (schema updated)." << std::endl;
  }

  // Check if any relevant files are newer than the database. No need to check
  // if the database has to be rebuilt for other reasons already.
  if (!rebuild_db) {
    vector<string> paths;
    aurostd::DirectoryLS(data_path, paths);
    if (paths.size() != _N_AUID_TABLES_) {
      string message = "Directory " + data_path + " is not a valid AUID directory.";
      throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
    }

    long int tm_db = aurostd::FileModificationTime(database_file);
    string dir;
    int i;
    for (i = 0; i < _N_AUID_TABLES_; i++) {
      dir = data_path + "/" + paths[i];
      if (aurostd::FileModificationTime(dir) > tm_db) break;
    }
    update_db = (i != _N_AUID_TABLES_);

    if (LDEBUG) {
      if (update_db) std::cerr << function << ": Updating database (new files found)." << std::endl;
      else std::cerr << function << ": No new files found. Database will not be updated." << std::endl;
    }
  }
  std::cout << aurostd::get_time() << std::endl;

  if (rebuild_db || update_db) {
    if (update_db) aurostd::CopyFile(database_file, database_file + ".tmp");
    openTmpFile();
    performUpdate(rebuild_db);
    std::cout << aurostd::get_time() << std::endl;
    return closeTmpFile();
  } else {
    return false;
  }
}

// Update --------------------------------------------------------------------

void AflowDB::performUpdate(bool rebuild_db) {
  // Do not constantly synchronize the database with file on disk.
  // Increases performance significantly.
  SQLexecuteCommand(db, "PRAGMA synchronous = OFF");
  vector<string> columns, types, tables, paths_to_update;
  long int tm_db = 0;
  
  if (rebuild_db) { // Need schema for table columns
    string schema = aurostd::file2string(schema_file);
    aurostd::StringSubst(schema, "'", "''");
    columns = getSchemaKeys(schema);
    uint ncols = columns.size();
    types.resize(ncols);
    for (uint i = 0; i < ncols; i++) types[i] = getSchemaValues(columns[i], "type", schema)[0];

    tables.resize(ntables_test);
    paths_to_update.resize(ntables_test);
    for (int i = 0; i < ntables_test; i++) {
      stringstream t;
      t << "auid_" << std::setfill('0') << std::setw(2) << std::hex << i;
      string table = t.str();
      createTable(table, columns, types);
      tables[i] = table;
      table = aurostd::StringSubst(table, "auid_", "aflow:");
      paths_to_update[i] = aurostd::CleanFileName(data_path + "/" + table);
    }
  } else { // Tables already exist, so no need to get the schema involved
    vector<string> paths;
    aurostd::DirectoryLS(data_path, paths);
    for (uint i = 0; i < paths.size(); i++) paths[i] = aurostd::CleanFileName(data_path + "/" + paths[i]);
    tm_db = aurostd::FileModificationTime(database_file);
    for (int i = 0; i < _N_AUID_TABLES_; i++) {
      if (aurostd::FileModificationTime(paths[i]) > tm_db) {
        paths_to_update.push_back(paths[i]);
      }
    }
    columns = getColumnNames(tables[0]);
    types = getColumnTypes(tables[0]);
  }

  int chunk_size = 1000;
  for (uint i = 0; i < paths_to_update.size(); i++) {
    vector<vector<string> > values = getUpdatedJsonData(paths_to_update[i], columns, types, rebuild_db, tm_db);
    int count = 0;
    uint nvals = values.size();
    transaction(true);
    for (uint v = 0; v < nvals; v++) {
      insertValues(tables[i], columns, values[i], !rebuild_db);
      if (++count % chunk_size == 0) {
        transaction(false);
        transaction(true);
      }
    }
    transaction(false);
  }

  if (rebuild_db) {
    // Create indexes on important database properties
    vector<string> index_cols;
    aurostd::string2tokens("auid,catalog", index_cols, ",");
    string index, index_expression;
    for (int t = 0; t < ntables_test; t++) {
      for (uint i = 0; i < index_cols.size(); i++) {
        index = "index_" + tables[t] + "_" + index_cols[i];
        createIndex(index, tables[t], index_cols[i]);
      }
      // Create special indexes on the species strings to accelerate database queries
      for (int e = 1; e < NUM_ELEMENTS; e++) {
        index = "index_" + tables[t] + "_" + vatom_symbol[e];
        index_expression = "INSTR(species, '\"" + vatom_symbol[e] + "\"')";
        createIndex(index, tables[t], index_expression);
      }
    }
  }
}

//getUpdatedJsonData//////////////////////////////////////////////////////////
// Fetches and processes all data from the JSON files in the AUID file system.
vector<vector<string> > AflowDB::getUpdatedJsonData(const string& parent_path,
                                                    const vector<string>& columns,
                                                    const vector<string>& types,
                                                    bool rebuild_db, long int tm_db) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  if (LDEBUG) {
    std::cerr << "AflowDB::getUpdateJsonData():"
              << " Fetching aflowlib.json files from path " << parent_path << "." << std::endl;
  }
  vector<string> paths, lsdir;
  string path;
  aurostd::DirectoryLS(parent_path, lsdir);
  for (uint i = 0; i < lsdir.size(); i++) {
    path = parent_path + "/" + lsdir[i];
    if (rebuild_db || aurostd::FileModificationTime(path) > tm_db) {
      paths.push_back(path);
    }
  }
  uint npaths = paths.size();
  vector<vector<vector<string> > > data_threads(npaths);

#ifdef AFLOW_DB_MULTITHREADS_ENABLE
  int ncpus = init::GetCPUCores();
  if (ncpus < 1) ncpus = 1;
  if (ncpus > _MAX_CPU_) ncpus = _MAX_CPU_;
  vector<vector<int> > thread_dist = getThreadDistribution(npaths, ncpus);
  vector<std::thread*> threads;
  for (int i = 0; i < ncpus; i++) {
    threads.push_back(new std::thread(&AflowDB::getUpdatedJsonDataThread, this,
                                      thread_dist[i][0], thread_dist[i][1],
                                      std::ref(paths), std::ref(columns), std::ref(types),
                                      rebuild_db, tm_db, std::ref(data_threads)));
  }
  for (int i = 0; i < ncpus; i++) {
    threads[i]->join();
    delete threads[i];
  }
#else
  getUpdatedJsonDataThread(0, npaths, paths, columns, types, rebuild_db, tm_db, data_threads);
#endif
  vector<vector<string> > data;
  for (uint i = 0; i < npaths; i++) {
    for (uint j = 0; j < data_threads[i].size(); j++) {
      data.push_back(data_threads[i][j]);
    }
  }
  if (LDEBUG) {
    vector<string> tokens;
    aurostd::string2tokens(parent_path, tokens, "/");
    std::cerr << "Found " << data.size() << " new entries for AUID " << tokens.back() << std::endl;
  }
  return data;
}

void AflowDB::getUpdatedJsonDataThread(int startIndex, int endIndex,
                                       const vector<string>& paths,
                                       const vector<string>& columns,
                                       const vector<string>& types,
                                       bool rebuild_db, long int tm_db,
                                       vector<vector<vector<string> > >& data) {
  for (int i = startIndex; i < endIndex; i++) {
    vector<string> json_files;
    crawl(json_files, paths[i], 2, 8, rebuild_db, tm_db);
    uint njson = json_files.size();
    data[i].resize(njson);
    for (uint j = 0; j < njson; j++) {
      data[i][j] = getDataValues(json_files[j], columns, types);
    }
  }
}

//crawl///////////////////////////////////////////////////////////////////////
// Walks the file system and retrieves the contents of the JSON files that
// need to be added to the database.
void AflowDB::crawl(vector<string>& json_files, const string& parent_path,
                    int depth, int maxdepth, bool rebuild_db, long int tm_db) {
  if (depth == maxdepth) {
    string filename = parent_path + "/aflowlib.json";
    if (aurostd::FileExist(filename)) {
      json_files.push_back(aurostd::file2string(filename));
    }
  } else {
    vector<string> child_paths;
    string child_path;
    aurostd::DirectoryLS(parent_path, child_paths);
    for (uint p = 0, npaths = child_paths.size(); p < npaths; p++) {
      child_path = parent_path + "/" + child_paths[p];
      if (rebuild_db || aurostd::FileModificationTime(child_path) > tm_db) {
        crawl(json_files, child_path, depth + 1, maxdepth, rebuild_db, tm_db);
      }
    }
  }
}

// Schema --------------------------------------------------------------------

//getSchemaKeys///////////////////////////////////////////////////////////////
// Gets the keys of the schema file.
vector<string> AflowDB::getSchemaKeys(const string& schema) {
  string command = "SELECT key FROM json_each('" + schema + "', '$.AAPI_schema')";
  command += " WHERE key IS NOT '__schema^2__'";
  return SQLexecuteCommandVECTOR(db, command);
}

//getSchemaValues/////////////////////////////////////////////////////////////
// Retrieves the desired metadata from the schema JSON file.
vector<string> AflowDB::getSchemaValues(const string& row, const string& col, const string& schema) {
  vector<string> cols(1, col);
  return getSchemaValues(row, cols, schema);
}

vector<string> AflowDB::getSchemaValues(const string& row, const vector<string>& cols, const string& schema) {
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

//getDataTypes////////////////////////////////////////////////////////////////
// Gets the data types of the schema keys and converts them into SQLite types.
// Note that SQLite does not recognize arrays, so they will be stored as text.
vector<string> AflowDB::getDataTypes(const vector<string>& cols, const string& schema) {
  uint ncols = cols.size();
  vector<string> types(ncols);
  for (uint i = 0; i < ncols; i++) {
    types[i] = getSchemaValues(cols[i], "type", schema)[0];
    // AUID has to be unique
    if (cols[i] == "auid") {
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

/**************************** DATABASE ANALYSIS *****************************/

namespace aflowlib {

//analyzeDatabase/////////////////////////////////////////////////////////////
// Provides analytics for each database data table in JSON format.
void AflowDB::analyzeDatabase(string outfile) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);

  // Get properties and tables for which statistics need to be collected
  // Since all tables have the same columns, only one table needs to be searched
  vector<string> tables = getTables();
  vector<string> columns = getColumnNames(tables[0]);

  string tab = "    ";
  std::stringstream json;
  json << "{" << std::endl << tab << "\"Aflow_DBs\": {" << std::endl;

  vector<string> catalogs = getSetMultiTables(tables, "catalog", true);
  uint ncatalogs = catalogs.size();
  for (uint c = 0; c < ncatalogs; c++) {
    std::cout << aurostd::get_time() << std::endl;
    if (LDEBUG) std::cerr << "Starting analysis for catalog " << catalogs[c] << std::endl;
    DBStats db_stats = getCatalogStats(catalogs[c], tables, columns);
    writeStatsToJSON(json, db_stats);
    std::cout << aurostd::get_time() << std::endl;
    if (c < ncatalogs - 1) json << ",";
    json << std::endl;
  }

  json << tab << "}" << std::endl << "}" << std::endl;
  aurostd::stringstream2file(json, outfile);
}

DBStats AflowDB::getCatalogStats(const string& catalog, const vector<string>& tables, const vector<string>& cols) {
  bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
  uint ncols = cols.size();

  DBStats stats;
  stats.catalog = catalog;
  stats.columns = cols;
  stats.count.assign(ncols, 0);
  stats.max.assign(ncols, "");
  stats.min.assign(ncols, "");
  stats.set.resize(ncols);

  string where = "catalog='" + catalog + "'";
  vector<string> entries = getPropertyMultiTables("COUNT", tables, "*", where);
  for (int t = 0; t < _N_AUID_TABLES_; t++) stats.nentries += aurostd::string2utype<int>(entries[t]);
  if (LDEBUG) std::cerr << "catalog = " << catalog << ", nentries = " << stats.nentries << std::endl;

  if (stats.nentries > 0) {
#ifdef AFLOW_DB_MULTITHREADS_ENABLE
    int ncpus = init::GetCPUCores();
    if (ncpus < 1) ncpus = 1;
    //if (ncpus > _MAX_CPU_) ncpus = _MAX_CPU_;
    vector<vector<int> > thread_dist = getThreadDistribution(ncols, ncpus);
    vector<std::thread*> threads;
    for (int t = 0; t < ncpus; t++) {
      threads.push_back(new std::thread(&AflowDB::getColStats, this,
                                        thread_dist[t][0], thread_dist[t][1],
                                        std::ref(stats), std::ref(tables)));
    }
    for (int t = 0; t < ncpus; t++) {
      threads[t]->join();
      delete threads[t];
    }
#else
    getColStats(0, ncols, stats, tables);
#endif
    stats.loop_counts = getDBLoopCounts(catalog, tables);
    vector<string> species = getSetMultiTables(tables, "species", true, where);
    stats.nsystems = species.size();
    stats.species = getUniqueFromJSONArrays(species);
  }
  return stats;
}

void AflowDB::getColStats(int startIndex, int endIndex,
                          DBStats& stats, const vector<string>& tables) {
  vector<int> counts(_N_AUID_TABLES_);

  sqlite3* cursor;
  int sql_code = sqlite3_open_v2(database_file.c_str(), &cursor, SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX, nullptr);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "open()";
    string message = "Could not open cursor on database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }

  vector<string> types = getColumnTypes(cursor, tables[0]);

  string where;
  for (int i = startIndex; i < endIndex; i++) {
    where = "catalog='" + stats.catalog + "'";
    where += " AND " + stats.columns[i]  + " NOT NULL";

    vector<string> tbls;
    vector<string> results_counts = getPropertyMultiTables(cursor, "COUNT", tables, stats.columns[i], where);
    for (int t = 0; t < _N_AUID_TABLES_; t++) {
      counts[t] = aurostd::string2utype<int>(results_counts[t]);
      if (counts[t] > 0) {
        stats.count[i] += counts[t];
        tbls.push_back(tables[t]);
      }
    }
    uint ntbls = tbls.size();

    vector<string> results(ntbls);
    vector<double> results_dbl(ntbls);
    if (stats.count[i] > 0) {  // No need to determine properties if count is zero
      results = getPropertyMultiTables(cursor, "MAX", tbls, stats.columns[i], where);
      if (types[i] == "number")  {
        for (uint t = 0; t < ntbls; t++) results_dbl[t] = aurostd::string2utype<double>(results[t]);
        // Converting the max of results_dbl into string leads to ugly numbers
        // in the final json, so sort and take the string in results. This does
        // not increase the runtime significantly.
        aurostd::sort(results_dbl, results);
        stats.max[i] = results.back();
      } else {
        stats.max[i] = aurostd::max(results);
      }
      results = getPropertyMultiTables(cursor, "MIN", tbls, stats.columns[i], where);
      if (types[i] == "number") {
        for (uint t = 0; t < ntbls; t++) results_dbl[t] = aurostd::string2utype<double>(results[t]);
        aurostd::sort(results_dbl, results);
        stats.min[i] = results[0];
      } else {
        stats.min[i] = aurostd::min(results);
      }
      stats.set[i] = getSetMultiTables(cursor, tbls, stats.columns[i], true, where, _DEFAULT_SET_LIMIT_ + 1);
    }
  }

  sql_code = sqlite3_close(cursor);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "open()";
    string message = "Could not close cursor on database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
}

/*
void AflowDB::getColStats(int startIndex, int endIndex,
                          DBStats& stats, const vector<string>& tables) {
  vector<int> counts(_N_AUID_TABLES_);

  sqlite3* cursor;
  int sql_code = sqlite3_open_v2(database_file.c_str(), &cursor, SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX, nullptr);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "open()";
    string message = "Could not open cursor on database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }

  vector<string> types = getColumnTypes(cursor, tables[0]);

  uint t;
  string where;
  for (int i = startIndex; i < endIndex; i++) {
    where = "catalog='" + stats.catalog + "'";
    where += " AND " + stats.columns[i]  + " NOT NULL";

    vector<string> tbls;
    vector<string> results_counts = getPropertyMultiTables(cursor, "COUNT", tables, stats.columns[i], where);
    for (t = 0; t < _N_AUID_TABLES_; t++) {
      counts[t] = aurostd::string2utype<int>(results_counts[t]);
      if (counts[t] > 0) {
        tbls.push_back(tables[t]);
        stats.count[i] += counts[t];
      }
    }
    uint ntbls = tbls.size();

    if (stats.count[i] > 0) {  // No need to determine properties if count is zero
      vector<double> results_dbl(ntbls);
      vector<vector<string> > v2_results_set(ntbls);
      vector<string> results(ntbls), results_set;

      for (t = 0; t < ntbls; t++) {
        results[t] = getProperty("MAX", tbls[t], stats.columns[i], where);
      }
      if (types[i] == "number")  {
        for (t = 0; t < ntbls; t++) results_dbl[t] = aurostd::string2utype<double>(results[t]);
        // Converting the max of results_dbl into string leads to ugly numbers
        // in the final json, so sort and take the string in results. This does
        // not increase the runtime significantly.
        aurostd::sort(results_dbl, results);
        stats.max[i] = results.back();
      } else {
        stats.max[i] = aurostd::max(results);
      }
      for (t = 0; t < ntbls; t++) {
        results[t] = getProperty("MIN", tbls[t], stats.columns[i], where);
      }
      if (types[i] == "number") {
        for (t = 0; t < ntbls; t++) results_dbl[t] = aurostd::string2utype<double>(results[t]);
        aurostd::sort(results_dbl, results);
        stats.min[i] = results[0];
      } else {
        stats.min[i] = aurostd::min(results);
      }
      for (t = 0; t < ntbls; t++) {
        v2_results_set[t] = getSet(tbls[t], stats.columns[i], true, where, _DEFAULT_SET_LIMIT_ + 1);
        uint nset = v2_results_set[t].size();
        for (uint j = 0; j < nset; j++) results_set.push_back(v2_results_set[t][j]);
        if (nset > _DEFAULT_SET_LIMIT_) break;  // Limit exceeded, so no need to process further
      }
      if (t == ntbls) {
        aurostd::sort_remove_duplicates(results_set);
        if (type == "number") {
          uint nset = results_set.size();
          vector<double> results_set_dbl(nset);
          for (uint j = 0; j < nset; j++) {
            results_set_dbl[j] = aurostd::string2utype<double>(results_set[j]);
          }
          aurostd::sort(results_set_dbl, results_set);
        }
      }
      stats.set[i] = results_set;
    }
  }

  sql_code = sqlite3_close(cursor);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "open()";
    string message = "Could not close cursor on database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
}
*/

//getDBLoopCounts/////////////////////////////////////////////////////////////
// Gets the number of times each loop type appears in the database. This is a
// proxy for the number of different property types that are stored in the
// database (e.g. elastic properties are stored as "agl" in the loop array).
vector<std::pair<string, int> > AflowDB::getDBLoopCounts(const string& catalog, const vector<string>& tables) {
  // First, get unique loop array elements
  vector<string> loop_entries, loops;
  loop_entries = getSetMultiTables(tables, "loop", true);
  loops = getUniqueFromJSONArrays(loop_entries);
  uint nloops = loops.size();

  vector<std::pair<string, int> > loop_counts(nloops);
  string where;
  string where_catalog = "catalog='" + catalog + "' AND ";
  vector<string> counts;
  for (uint l = 0; l < nloops; l++) {
    where = where_catalog + "loop LIKE '%" + loops[l] + "%'";
    counts = getPropertyMultiTables("COUNT", tables, "loop", where);
    loop_counts[l].first = loops[l];
    loop_counts[l].second = 0;
    for (int t = 0; t < _N_AUID_TABLES_; t++) loop_counts[l].second += aurostd::string2utype<int>(counts[t]);
  }
  return loop_counts;
}

//getUniqueFromJSONArrays/////////////////////////////////////////////////////
// Determines the unique array elements in a set of 1D-array strings.
vector<string> AflowDB::getUniqueFromJSONArrays(const vector<string>& arrays) {
  vector<string> unique, tokens;
  string arr;
  int nunique = 0, u = 0;
  for (uint a = 0; a < arrays.size(); a++) {
    vector<string> tokens;
    arr = arrays[a];
    arr = aurostd::RemoveSubString(arr, "[");
    arr = aurostd::RemoveSubString(arr, "]");
    arr = aurostd::RemoveSubString(arr, "\"");
    aurostd::string2tokens(arr, tokens, ", ");
    for (uint t = 0; t < tokens.size(); t++) {
      if (nunique == 0) {
        unique.push_back(tokens[t]);
        nunique = 1;
      } else {
        for (u = 0; u < nunique; u++) {
          if (unique[u] == tokens[t]) break;
        }
        if (u == nunique) {
          unique.push_back(tokens[t]);
          nunique++;
        }
      }
    }
  }
  return unique;
}

//writeStatsToJSON////////////////////////////////////////////////////////////
// Writes the database statistics into a JSON-formatted string(stream).
void AflowDB::writeStatsToJSON(std::stringstream& json, const DBStats& db_stats) {
  string tab = "    ";
  string indent = tab + tab;
  json << indent << "\"" << db_stats.catalog << "\": {" << std::endl;
  json << indent << tab << "\"count\": " << db_stats.nentries << "," << std::endl;
  json << indent << tab << "\"systems\": " <<  db_stats.nsystems << "," << std::endl;
  for (uint l = 0; l < db_stats.loop_counts.size(); l++) {
    json << indent << tab << "\"" << db_stats.loop_counts[l].first << "\": ";
    json << db_stats.loop_counts[l].second << "," << std::endl;
  }
  json << indent << tab << "\"columns\": {" << std::endl;
  uint ncols = db_stats.columns.size();
  string str_formatted;
  for (uint c = 0; c < ncols; c++) {
    json << indent << tab << tab << "\"" << db_stats.columns[c] << "\": {" << std::endl;
    json << indent << tab << tab << tab << "\"count\": ";
    json << aurostd::utype2string<int>(db_stats.count[c]) << "," << std::endl;
    str_formatted = db_stats.min[c];
    str_formatted = aurostd::StringSubst(str_formatted, "\"", "\\\"");
    json << indent << tab << tab << tab << "\"min\": "
         << (db_stats.min[c].empty()?"null":string("\"" + str_formatted + "\"")) << "," << std::endl;
    str_formatted = db_stats.max[c];
    str_formatted = aurostd::StringSubst(str_formatted, "\"", "\\\"");
    json << indent << tab << tab << tab << "\"max\": "
         << (db_stats.max[c].empty()?"null":string("\"" + str_formatted + "\"")) << "," << std::endl;

    // Write set
    uint nset = db_stats.set[c].size();
    json << indent << tab << tab << tab << "\"set\": ";
    if (nset > _DEFAULT_SET_LIMIT_) {
      json << "null" << std::endl;
    } else if (nset == 0) {
      json << "[]" << std::endl;
    } else {
      json << "[" << std::endl;
      for (uint s = 0; s < nset; s++) {
        str_formatted = db_stats.set[c][s];
        str_formatted = aurostd::StringSubst(str_formatted, "\"", "\\\"");
        json << indent << tab << tab << tab << tab << "\"" << str_formatted << "\"";
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

// These functions are higher level SQLite functions that call functions of
// the SQLite interface. They are essentially syntactic sugar to make the code
// easier to read and to facilitate the implementation of the AflowDB class
// into other parts of AFLOW.
//
// Getter functions are overloaded to take a cursor other than the default
// cursor. Since SQLite does not allow concurrent writing, the writer
// functions should not be overloaded.

namespace aflowlib {

// INDEX ---------------------------------------------------------------------

//createIndex/////////////////////////////////////////////////////////////////
// Creates an index.
void AflowDB::createIndex(const string& index, const string& table, const string& column) {
  string command = "CREATE INDEX " + index + " ON " + table + "(" + column +  ");";
  SQLexecuteCommand(db, command);
}

//dropIndex///////////////////////////////////////////////////////////////////
// Removes an index.
void AflowDB::dropIndex(const string& index) {
  string command = "DROP INDEX " + index + ";";
  SQLexecuteCommand(db, command);
}

// JSON ----------------------------------------------------------------------

//extractJSONvalue////////////////////////////////////////////////////////////
// Extracts a JSON value using SQLite's JSON extension.
string AflowDB::extractJSONvalue(const string& json, const string& key) {
  string command = "SELECT json_extract('" + json + "', '$." + key + "');";
  return SQLexecuteCommandSCALAR(db, command); 
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
  SQLexecuteCommand(db, command);
}

// TABLE ---------------------------------------------------------------------

//dropTable///////////////////////////////////////////////////////////////////
// Deletes a table from the database.
void AflowDB::dropTable(const string& table) {
  string command = "DROP TABLE IF EXISTS " + table + ";";
  SQLexecuteCommand(db, command);
}

//getTables///////////////////////////////////////////////////////////////////
// Retrieves a set of tables. If where is empty, all tables in the database
// will be returned.
vector<string> AflowDB::getTables(string where) {
  return getTables(db, where);
}

vector<string> AflowDB::getTables(sqlite3* cursor, string where) {
  string command = "SELECT name FROM sqlite_master WHERE type='table'";
  if (!where.empty()) command += " AND (" + where + ")";
  command += ";";
  return SQLexecuteCommandVECTOR(cursor, command);
}

//createTable/////////////////////////////////////////////////////////////////
// Creates a table where all columns have the same type.
void AflowDB::createTable(const string& table, const vector<string>& cols, const string& type) {
  vector<string> types(cols.size(), type);
  createTable(table, cols, types);
}

// Creates a table where each column is assigned its own type
void AflowDB::createTable(const string& table, const vector<string>& cols, const vector<string>& types) {
  uint ncols = cols.size();
  if (ncols != types.size()) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "createTable()";
    string message = "Could not create table. ";
    message += "Number of columns and number of types do not match.";
    throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
  } else {
    string command = "CREATE TABLE IF NOT EXISTS " + table + " (";
    for (uint c = 0; c < ncols; c++) {
      command += cols[c] + " " + types[c];
      if (c < ncols - 1) command += ", ";
    }
    command += ");";
    SQLexecuteCommand(db, command);
  }
}

// INSERT --------------------------------------------------------------------

//insertValues////////////////////////////////////////////////////////////////
// Inserts a set of values into a table.
void AflowDB::insertValues(const string& table, const vector<string>& vals, bool update) {
  vector<string> cols;
  insertValues(table, cols, vals, update);
}

// Inserts a set of values into a table (if cols is empty) or into specific
// columns of the table.
void AflowDB::insertValues(const string& table, const vector<string>& cols,
                           const vector<string>& vals, bool update) {
  uint ncols = cols.size();
  uint nvals = vals.size();
  if ((ncols > 0) && (ncols != nvals)) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "insertValues()";
    string message = "Could not insert values. ";
    message += "Number of columns and number of values do not match.";
    throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
  } else {
    string command;
    if (update) command = "REPLACE";
    else command = "INSERT";
    command += " INTO " + table;
    if (ncols > 0) command += " (" + aurostd::joinWDelimiter(cols, ", ") + ")";
    command += " VALUES(" + aurostd::joinWDelimiter(vals, ", ") + ");";
    SQLexecuteCommand(db, command);
  }
}

// GET -----------------------------------------------------------------------

//getColumnNames//////////////////////////////////////////////////////////////
// Returns the names of all columns in a specific table.
vector<string> AflowDB::getColumnNames(const string& table) {
  return getColumnNames(db, table);
}

vector<string> AflowDB::getColumnNames(sqlite3* cursor, const string& table) {
  string command = "PRAGMA table_info(" + table + ")";
  vector<vector<string> > pragma_results = SQLexecuteCommand2DVECTOR(cursor, command);
  uint nresults = pragma_results.size();
  vector<string> columns(nresults);
  for (uint r = 0; r < nresults; r++) {
    columns[r] = pragma_results[r][1];
  }
  return columns;
}

//getColumnTypes//////////////////////////////////////////////////////////////
// Returns the data types of all columns in a specific table.
vector<string> AflowDB::getColumnTypes(const string& table) {
  return getColumnTypes(db, table);
}

vector<string> AflowDB::getColumnTypes(sqlite3* cursor, const string& table) {
  string command = "PRAGMA table_info(" + table + ")";
  vector<vector<string> > pragma_results = SQLexecuteCommand2DVECTOR(cursor, command);
  uint nresults = pragma_results.size();
  vector<string> columns(nresults);
  for (uint r = 0; r < nresults; r++) {
    columns[r] = pragma_results[r][2];
  }
  return columns;
}

//getValue////////////////////////////////////////////////////////////////////
// Gets a value from a specific column. The row must be specified in the where
// condition, or else it just takes the first value.
string AflowDB::getValue(const string& table, const string& col, string where) {
  return getValue(db, table, col, where);
}

string AflowDB::getValue(sqlite3* cursor, const string& table, const string& col, string where) {
  string command = prepareSELECT(table, "", col, where, 0, "");
  return SQLexecuteCommandSCALAR(cursor, command);
}

//getValueMultiCol////////////////////////////////////////////////////////////
// Gets a value from multiple columns. The row must be specified in the where
// condition, or else it just takes the first value.
vector<string> AflowDB::getValueMultiCol(const string& table, const vector<string>& cols, string where) {
  return getValueMultiCol(db, table, cols, where);
}

vector<string> AflowDB::getValueMultiCol(sqlite3* cursor, const string& table, const vector<string>& cols, string where) {
  string command = prepareSELECT(table, "", cols, where, 0, "");
  vector<vector<string> > values = SQLexecuteCommand2DVECTOR(cursor, command);
  return values[0];
}

//getProperty/////////////////////////////////////////////////////////////////
// Gets a database property for a specific column.
string AflowDB::getProperty(const string& property, const string& table,
                            const string& col, string where) {
  return getProperty(db, property, table, col, where);
}
string AflowDB::getProperty(sqlite3* cursor, const string& property, const string& table,
                            const string& col, string where) {
  string command = prepareSELECT(table, property, col, where, 0, "");
  return SQLexecuteCommandSCALAR(cursor, command);
}

//getPropertyMultiTables//////////////////////////////////////////////////////
// Gets a database property for a specific column across multiple tables.
vector<string> AflowDB::getPropertyMultiTables(const string& property, const vector<string>& tables,
                                               const string& col, string where) {
  return getPropertyMultiTables(db, property, tables, col, where);
}

vector<string> AflowDB::getPropertyMultiTables(sqlite3* cursor, const string& property, const vector<string>& tables,
                                               const string& col, string where) {
  uint ntables = tables.size();
  vector<string> commands(ntables);
  for (uint t = 0; t < ntables; t++) {
    commands[t] = prepareSELECT(tables[t], property, col, where);
    commands[t].pop_back();  // Remove semicolon at end
  }
  return SQLexecuteCommandVECTOR(cursor, aurostd::joinWDelimiter(commands, " UNION ALL ") + ";");
}

//getColumn///////////////////////////////////////////////////////////////////
// Returns an entire column. The output can be sorted by another column using
// the "order_by" string.
vector<string> AflowDB::getColumn(const string& table, const string& col, string order_by) {
  return getColumn(db, table, col, order_by);
}

vector<string> AflowDB::getColumn(sqlite3* cursor, const string& table, const string& col, string order_by) {
  return getSet(cursor, table, col, false, "", 0, order_by);
}

//getSet//////////////////////////////////////////////////////////////////////
// Retrieves a (distinct) set from a single column.
vector<string> AflowDB::getSet(const string& table, const string& col, bool distinct,
                               string where, int limit, string order_by) {
  return getSet(db, table, col, distinct, where, limit, order_by);
}

vector<string> AflowDB::getSet(sqlite3* cursor, const string& table, const string& col, bool distinct,
                               string where, int limit, string order_by) {
  string property = string((distinct?"DISTINCT":""));
  string command = prepareSELECT(table, property, col, where, limit, order_by);
  return SQLexecuteCommandVECTOR(cursor, command);
}

//getSetMulitTables///////////////////////////////////////////////////////////
// Retrieves a (distinct) set from a single column across multiple tables.
// The result is sorted already, so there is not need for order_by.
vector<string> AflowDB::getSetMultiTables(const vector<string>& tables, const string& col,
                                          bool distinct, string where, int limit) {
  return getSetMultiTables(db, tables, col, distinct, where, limit);
}
vector<string> AflowDB::getSetMultiTables(sqlite3* cursor, const vector<string>& tables, const string& col,
                                          bool distinct, string where, int limit) {
  string property = string((distinct?"DISTINCT":""));
  uint ntables = tables.size();
  vector<string> commands(ntables);
  for (uint t = 0; t < ntables; t++) {
    commands[t] = prepareSELECT(tables[t], property, col, where, 0);
    commands[t].pop_back();  // Remove semicolon at end
  }
  string union_string = " UNION ";
  if (!distinct) union_string += "ALL ";
  string command = aurostd::joinWDelimiter(commands, union_string);
  if (limit > 0) command += " LIMIT " + aurostd::utype2string<int>(limit);
  command += ";";
  return SQLexecuteCommandVECTOR(cursor, command);
}

//getSetMultiCol//////////////////////////////////////////////////////////////
// Retrieves a (distinct) set from a multiple columns.
vector<vector<string> > AflowDB::getSetMultiCol(const string& table, const vector<string>& cols, bool distinct,
                                                string where, int limit, string order_by) {
  return getSetMultiCol(db, table, cols, distinct, where, limit, order_by);
}

vector<vector<string> > AflowDB::getSetMultiCol(sqlite3* cursor, const string& table, const vector<string>& cols, bool distinct,
                                                string where, int limit, string order_by) {
  string property = string(distinct?"DISTINCT":"");
  string command = prepareSELECT(table, property, aurostd::joinWDelimiter(cols, ", "), where, limit, order_by);
  return SQLexecuteCommand2DVECTOR(cursor, command);
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

void AflowDB::SQLexecuteCommand(sqlite3* cursor, const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << "AflowDB::SQLexecuteCommand(): command = " << command << std::endl;
  char* sqlErrMsg = 0;
  int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallback, 0, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "SQLexecuteCommand()";
    string message = string(sqlErrMsg) + " in command " + command;
    throw aurostd::xerror(function, message, _RUNTIME_SQL_);
  }
}

string AflowDB::SQLexecuteCommandSCALAR(sqlite3* cursor, const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << "AflowDB::SQLexecuteCommandSCALAR(): command = " << command << std::endl;
  char* sqlErrMsg = 0;
  string returnstring;
  int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallbackSCALAR, &returnstring, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "SQLexecuteCommandSCALAR()";
    string message = string(sqlErrMsg) + " in command " + command;
    message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
    throw aurostd::xerror(function, message, _RUNTIME_SQL_);
  } else {
    return returnstring;
  }
}

vector<string> AflowDB::SQLexecuteCommandVECTOR(sqlite3* cursor, const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "SQLexecuteCommandVECTOR(): command = " << command << std::endl;
  char *sqlErrMsg = 0;
  vector<string> returnvector;
  int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallbackVECTOR, &returnvector, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "SQLexecuteCommandVECTOR()";
    string message = string(sqlErrMsg) + " in command " + command;
    message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
    throw aurostd::xerror(function, message, _RUNTIME_SQL_);
  } else {
    return returnvector;
  }
}

vector<vector<string> > AflowDB::SQLexecuteCommand2DVECTOR(sqlite3* cursor, const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << _AFLOW_DB_ERR_PREFIX_ << "SQLexecuteCommand2DVECTOR(): command = " << command << std::endl;
  char *sqlErrMsg = 0;
  vector<vector<string> > returnvector;
  int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallback2DVECTOR, &returnvector, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "SQLexecuteCommand2DVECTOR()";
    string message = string(sqlErrMsg) + " in command " + command;
    message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
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
