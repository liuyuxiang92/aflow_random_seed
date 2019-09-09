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

#define _AFLOW_DB_DEBUG_     false

using std::string;
using std::vector;

static const string _AFLOW_DB_ERR_PREFIX_ = "AflowDB::";
static const int _DEFAULT_SET_LIMIT_ = 16;
static const int _N_AUID_TABLES_ = 256;

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

// Opens the database file and creates the main a cursor
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

// Closes the database file and removes the main cursor
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
  bool rebuild_db = false;

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
      for (k = 0; k < nkeys; k++) types_schema[k] = extractJsonValue(schema, "AAPI_schema." + keys[k] + ".type");
      for (k = 0; k < nkeys; k++) {
        if (!aurostd::withinList(columns, keys[k])) break;
        if (!aurostd::withinList(types_db, types_schema[k])) break;
      }
      rebuild_db = (k != nkeys);
    }
    if (LDEBUG && rebuild_db) std::cerr << function << ": Rebuilding database (schema updated)." << std::endl;
  }

  // Check if any relevant files are newer than the database.
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
    rebuild_db = (i != _N_AUID_TABLES_);

    if (LDEBUG) {
      if (rebuild_db) std::cerr << function << ": Rebuilding database (new files found)." << std::endl;
      else std::cerr << function << ": No new files found. Database will not be rebuilt." << std::endl;
    }
  }
  std::cerr << aurostd::get_time() << std::endl;

  if (rebuild_db) {
    openTmpFile();
    updateDatabaseJsonFiles(data_path);
    std::cerr << aurostd::get_time() << std::endl;
    rebuildDB();
    std::cerr << aurostd::get_time() << std::endl;
    return closeTmpFile();
  } else {
    return false;
  }
}

// Rebuild -------------------------------------------------------------------

void AflowDB::rebuildDB() {
  // Do not constantly synchronize the database with file on disk.
  // Increases performance significantly.
  SQLexecuteCommand(db, "PRAGMA synchronous = OFF");

  // Get columns from schema
  string schema = aurostd::file2string(schema_file);
  aurostd::StringSubst(schema, "'", "''");
  vector<string> columns = getSchemaKeys(schema);
  vector<string> types = getDataTypes(schema, columns);

#ifdef AFLOW_DB_MULTITHREADS_ENABLE
  int ncpus = init::GetCPUCores();
  if (ncpus < 1) ncpus = 1;
  vector<vector<int> > thread_dist = getThreadDistribution(_N_AUID_TABLES_, ncpus);
  vector<std::thread*> threads;
  for (int i = 0; i < ncpus; i++) {
    threads.push_back(new std::thread(&AflowDB::buildTables, this,
                                      thread_dist[i][0], thread_dist[i][1],
                                      std::ref(columns), std::ref(types)));
  }
  for (int i = 0; i < ncpus; i++) {
    threads[i]->join();
    delete threads[i];
  }
#else
  buildTables(0, _N_AUID_TABLES_, columns, types);
#endif
}

void AflowDB::buildTables(int startIndex, int endIndex, const vector<string>& columns, const vector<string>& types) {
  int chunk_size = 1000;

  for (int i = startIndex; i < endIndex; i++) {
    stringstream t;
    t << std::setfill('0') << std::setw(2) << std::hex << i;
    string table = "auid_" + t.str();
    createTable(table, columns, types);

    string jsonfile = aurostd::CleanFileName(data_path + "/aflow:" + t.str() + ".json");
    vector<string> data;
    aurostd::file2vectorstring(jsonfile, data);

    int count = 0;
    vector<vector<string> > values(chunk_size);
    for (uint d = 0, ndata = data.size(); d < ndata; d++) {
      values[count] = getDataValues(data[d], columns, types);
      if (++count % chunk_size == 0) {
        transaction(true);
        for (int v = 0; v < chunk_size; v++) insertValues(table, columns, values[v]);
        transaction(false);
        count = 0;
      }
    }

    transaction(true);
    for (int v = 0; v < count; v++) insertValues(table, columns, values[v]);
    // Create indexes on important database properties
    vector<string> index_cols;
    aurostd::string2tokens("auid,catalog,nspecies", index_cols, ",");
    string index, index_expression;
    for (int t = 0; t < _N_AUID_TABLES_; t++) {
      for (uint i = 0; i < index_cols.size(); i++) {
        index = "index_" + table + "_" + index_cols[i];
        createIndex(index, table, index_cols[i]);
      }
      // Create special indexes on the species strings to accelerate database queries
      for (int e = 1; e < NUM_ELEMENTS; e++) {
        index = "index_" + table + "_" + vatom_symbol[e];
        index_expression = "INSTR(species, '\"" + vatom_symbol[e] + "\"')";
        createIndex(index, table, index_expression);
      }
    }
    transaction(false);
  }
}

// Schema --------------------------------------------------------------------

vector<string> AflowDB::getSchemaKeys(const string& schema) {
  vector<string> keys_unfiltered = getJsonKeys(schema, "AAPI_schema");

  vector<string> keys;
  string function;
  for (uint k = 0; k < keys_unfiltered.size(); k++) {
    if (keys_unfiltered[k] != "__schema^2__") {
      function = extractJsonValue(schema, "AAPI_schema." + keys_unfiltered[k] + ".function");
      if ((function != "link") && (function != "image")) {
        keys.push_back(keys_unfiltered[k]);
      }
    }
  }
  return keys;
}

// Data ----------------------------------------------------------------------

//getDataTypes////////////////////////////////////////////////////////////////
// Gets the data types of the schema keys and converts them into SQLite types.
// Note that SQLite does not recognize arrays, so they will be stored as text.
vector<string> AflowDB::getDataTypes(const string& schema, const vector<string>& cols) {
  uint ncols = cols.size();
  vector<string> types(ncols);
  for (uint i = 0; i < ncols; i++) {
    types[i] = extractJsonValue(schema, "AAPI_schema." + cols[i] + ".type");
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
    value = extractJsonValueAflow(entry, cols[c]);
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
// Provides analytics for the database in JSON format.
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
    std::cerr << aurostd::get_time() << std::endl;
    if (LDEBUG) std::cerr << "Starting analysis for catalog " << catalogs[c] << std::endl;
    DBStats db_stats = getCatalogStats(catalogs[c], tables, columns);
    writeStatsToJSON(json, db_stats);
    std::cerr << aurostd::get_time() << std::endl;
    if (c < ncatalogs - 1) json << ",";
    json << std::endl;
  }

  json << tab << "}" << std::endl << "}" << std::endl;
  aurostd::stringstream2file(json, outfile);
}

//getCatalogStats/////////////////////////////////////////////////////////////
// Gets the statistics for all properties in the catalog.
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
    vector<vector<vector<string> >  > maxmin(_N_AUID_TABLES_, vector<vector<string> >(ncols, vector<string>(2))), sets(_N_AUID_TABLES_, vector<vector<string> >(ncols));
    vector<vector<int> > counts(_N_AUID_TABLES_, vector<int>(ncols));
#ifdef AFLOW_DB_MULTITHREADS_ENABLE
    int ncpus = init::GetCPUCores();
    if (ncpus < 1) ncpus = 1;
    // The maximum number of CPUs are empirically found values that appear
    // to result in the lowest run times. Further testing may be necessary.
    int cpu_max = 32;
    if (stats.nentries < 100000) cpu_max = 16;
    if (stats.nentries < 50000) cpu_max = 8;
    if (stats.nentries < 10000) cpu_max = 4;
    if (ncpus > cpu_max) ncpus = cpu_max;
    vector<vector<int> > thread_dist = getThreadDistribution(_N_AUID_TABLES_, ncpus);
    vector<std::thread*> threads;
    for (int i = 0; i < ncpus; i++) {
      threads.push_back(new std::thread(&AflowDB::getColStats, this,
                                        thread_dist[i][0], thread_dist[i][1],
                                        std::ref(catalog), std::ref(tables),
                                        std::ref(cols), std::ref(counts),
                                        std::ref(maxmin), std::ref(sets)));
    }
    for (int i = 0; i < ncpus; i++) {
      threads[i]->join();
      delete threads[i];
    }
#else
    getColStats(0, _N_AUID_TABLES_, catalog, tables, cols, counts, maxmin, sets);
#endif

    // Post-processing
    vector<string> types = getColumnTypes(tables[0]);
    for (uint c = 0; c < ncols; c++) {
      for (int t = 0; t < _N_AUID_TABLES_; t++) stats.count[c] += counts[t][c];
      if (stats.count[c] > 0) {
        string max, min;
        vector<string> set;
        uint nset, n;
        for (int t = 0; t < _N_AUID_TABLES_; t++) {
          if (counts[t][c] > 0) {
            if (max.empty()) {
              max = maxmin[t][c][0];
            } else {
              if (types[c] == "number") {
                if (aurostd::string2utype<double>(maxmin[t][c][0]) > aurostd::string2utype<double>(max)) max = maxmin[t][c][0];
              } else {
                if (maxmin[t][c][0] > max) max = maxmin[t][c][0];
              }
            }
            if (min.empty()) {
              min = maxmin[t][c][1];
            } else {
              if (types[c] == "number") {
                if (aurostd::string2utype<double>(maxmin[t][c][1]) < aurostd::string2utype<double>(min)) min = maxmin[t][c][1];
              } else {
                if (maxmin[t][c][1] > min) min = maxmin[t][c][1];
              }
            }
            if (nset <= _DEFAULT_SET_LIMIT_) {
              n = sets[t][c].size();
              if (n > _DEFAULT_SET_LIMIT_) {
                set = sets[t][c];
                nset = n;
              } else {
                for (uint i = 0; i < n; i++) {
                  if (!aurostd::withinList(set, sets[t][c][i])) {
                    set.push_back(sets[t][c][i]);
                    nset++;
                  }
                  if (nset > _DEFAULT_SET_LIMIT_) break;
                }
              }
            }
          }
        }
        stats.max[c] = max;
        stats.min[c] = min;
        if (nset <= _DEFAULT_SET_LIMIT_) {
          if (types[c] == "number") {
            vector<double> set_dbl(nset);
            for (uint i = 0; i < nset; i++) set_dbl[i] = aurostd::string2utype<double>(set[i]);
            aurostd::sort(set_dbl, set);
          } else {
            std::sort(set.begin(), set.end());
          }
        }
        stats.set[c] = set;
      }
    }
  }
  return stats;
}

//getColStats/////////////////////////////////////////////////////////////////
// Retrieves the statistics for each database property.
void AflowDB::getColStats(int startIndex, int endIndex,
                          const string& catalog, const vector<string>& tables,
                          const vector<string>& cols, vector<vector<int> >& counts,
                          vector<vector<vector<string> > >& maxmin, vector<vector<vector<string> > >& sets) {
  sqlite3* cursor;
  int sql_code = sqlite3_open_v2(database_file.c_str(), &cursor, SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX, nullptr);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "open()";
    string message = "Could not open cursor on database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }

  uint ncols = cols.size();
  string where;
  for (int i = startIndex; i < endIndex; i++) {
    for (uint c = 0; c < ncols; c++) {
      where = "catalog='" + catalog + "' AND " + cols[c] + " NOT NULL";
      counts[i][c] = aurostd::string2utype<int>(getProperty(cursor, "COUNT", tables[i], cols[c], where));
      if (counts[i][c] > 0) {
        maxmin[i][c][0] = getProperty(cursor, "MAX", tables[i], cols[c], where);
        maxmin[i][c][1] = getProperty(cursor, "MIN", tables[i], cols[c], where);
        sets[i][c] = getSet(cursor, tables[i], cols[c], true, where, _DEFAULT_SET_LIMIT_ + 1);
      }
    }
  }

  sql_code = sqlite3_close(cursor);
  if (sql_code != SQLITE_OK) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "open()";
    string message = "Could not close cursor on database file " + database_file + ".";
    throw aurostd::xerror(function, message, _FILE_ERROR_);
  }
}

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
  string command = "CREATE INDEX " + index + " ON " + table + "(" + column +  ")";
  SQLexecuteCommand(db, command);
}

//dropIndex///////////////////////////////////////////////////////////////////
// Removes an index.
void AflowDB::dropIndex(const string& index) {
  string command = "DROP INDEX " + index;
  SQLexecuteCommand(db, command);
}

// TRANSACTION ---------------------------------------------------------------

void AflowDB::transaction(bool begin) {
  string command = string(begin?"BEGIN":"END") + " TRANSACTION";
  SQLexecuteCommand(db, command);
}

// TABLE ---------------------------------------------------------------------

//dropTable///////////////////////////////////////////////////////////////////
// Deletes a table from the database.
void AflowDB::dropTable(const string& table) {
  string command = "DROP TABLE IF EXISTS " + table;
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
    command += ")";
    SQLexecuteCommand(db, command);
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
void AflowDB::insertValues(const string& table, const vector<string>& cols,
                           const vector<string>& vals) {
  uint ncols = cols.size();
  uint nvals = vals.size();
  if ((ncols > 0) && (ncols != nvals)) {
    string function = _AFLOW_DB_ERR_PREFIX_ + "insertValues()";
    string message = "Could not insert values. ";
    message += "Number of columns and number of values do not match.";
    throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
  } else {
    string command = "INSERT INTO " + table;
    if (ncols > 0) command += " (" + aurostd::joinWDelimiter(cols, ", ") + ")";
    command += " VALUES(" + aurostd::joinWDelimiter(vals, ", ") + ")";
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
  return SQLexecuteCommandVECTOR(cursor, aurostd::joinWDelimiter(commands, " UNION ALL "));
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
  return command.str();
}

// Wrapper function for a vector representation of the columns.
string AflowDB::prepareSELECT(const string& table, const string& property, const vector<string>& cols,
                              string where, int limit, string order_by) {
  return prepareSELECT(table, property, aurostd::joinWDelimiter(cols, ", "), where, limit, order_by);
}

}  // namespace aflowlib

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
//****************************************************************************
