//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
//****************************************************************************

#include "aflowlib.h"
#include "SQLITE/sqlite3.h"

#define _SQL_COMMAND_DEBUG_  false  // debug SQL commands that are sent - verbose output
#define _SQL_CALLBACK_DEBUG_ false  // debug SQL callbacks - verbose output

// Some parts are written within the C++0x support in GCC, especially std::thread,
// which is implemented in gcc 4.4 and higher. For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400
#define AFLOW_DB_MULTITHREADS_ENABLE
#include <thread>
static const int _MAX_CPU_ = 16;
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif

static const int _N_AUID_TABLES_ = 256;

/*********************************** JSON ***********************************/

namespace aflowlib {

void updateDatabaseJsonFiles(const string& data_path) {
  string auid_path = vAFLOW_PROJECTS_DIRECTORIES.at(XHOST_LIBRARY_AUID) + "/RAW/";

  // Fetch the files to update first
  long int tm;
  string json, path;
  vector<long int> mod_times;
  vector<string> json_files, entry_paths;
  for (int i = 0; i < _N_AUID_TABLES_; i++) {
    stringstream t;
    t << std::setfill('0') << std::setw(2) << std::hex << i;
    json = aurostd::CleanFileName(data_path + "/aflow:" + t.str() + ".json");
    if (aurostd::FileExist(json)) tm = aurostd::FileModificationTime(json);
    else tm = 0;
    path = auid_path + "aflow:" + t.str();
    if (aurostd::FileModificationTime(path) > tm) {
      json_files.push_back(json);
      entry_paths.push_back(path);
      mod_times.push_back(tm);
    }
  }

  int njson = (int) json_files.size();
#ifdef AFLOW_DB_MULTITHREADS_ENABLE
  int ncpu = init::GetCPUCores();
  if (ncpu < 1) ncpu = 1;
  if (ncpu > _MAX_CPU_) ncpu = _MAX_CPU_;
  if (ncpu > njson) ncpu = njson;
  vector<vector<int> > thread_dist = getThreadDistribution(njson, ncpu);
  vector<std::thread*> threads;
  for (int i = 0; i < ncpu; i++) {
    threads.push_back(new std::thread(&updateDatabaseJsonFilesThread,
                                      thread_dist[i][0], thread_dist[i][1],
                                      std::ref(json_files), std::ref(entry_paths), std::ref(mod_times)));
  }
  for (int i = 0; i < ncpu; i++) {
    threads[i]->join();
    delete threads[i];
  }
#else
  updateDatabaseJsonFilesThread(0, njson, json_files, entry_paths, mod_times);
#endif
}

void updateDatabaseJsonFilesThread(int startIndex, int endIndex, const vector<string>& json_files,
                                   const vector<string>& entry_paths, const vector<long int>& mod_times) {
  for (int i = startIndex; i < endIndex; i++) {
    vector<string> entries;
    fetchUpdatedJsonFiles(entries, entry_paths[i], 1, 8, mod_times[i]);
    if (mod_times[i] > 0) {
      vector<string> json;
      aurostd::file2vectorstring(json_files[i], json);
      uint njson = json.size();
      uint nentry = entries.size();
      vector<string> auids(njson);
      uint j;
      for (j = 0; j < njson; j++) auids[j] = extractJsonValueAflow(json[j], "auid");

      string auid;
      for (uint e = 0; e < nentry; e++) {
        for (j = 0; j < njson; j++) {
          auid = extractJsonValueAflow(entries[e], "auid");
          if (auid == auids[j]) break;
        }
        if (j == njson) json.push_back(entries[e]);
        else json[j] = entries[e];
      }
      aurostd::string2file(aurostd::joinWDelimiter(json, ","), json_files[i]);
    } else {
      aurostd::string2file(aurostd::joinWDelimiter(entries, ","), json_files[i]);
    }
  }
}

//walkAuidPathDB//////////////////////////////////////////////////////////////
// Walks the file system and retrieves the contents of the JSON files that
// need to be added to the database.
void fetchUpdatedJsonFiles(vector<string>& json_files, const string& parent_path,
                           int depth, int maxdepth, long int mod_time) {
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
      if (mod_time == 0 || aurostd::FileModificationTime(child_path) > mod_time) {
        fetchUpdatedJsonFiles(json_files, child_path, depth + 1, maxdepth, mod_time);
      }
    }
  }
}

//getJsonKeys/////////////////////////////////////////////////////////////////
// Gets the keys of a JSON object using SQLite's JSON extension.
vector<string> getJsonKeys(const string& json, string address) {
  sqlite3* cursor;
  sqlite3_open(":memory:", &cursor);
  string command = "SELECT key FROM json_each('" + json + "'";
  if (!address.empty()) command += ", '$." + address + "'";
  command += ")";
  return SQLexecuteCommandVECTOR(cursor, command);
}

//extractJsonValue////////////////////////////////////////////////////////////
// Extracts a JSON value using SQLite's JSON extension.
string extractJsonValue(const string& json, const string& key) {
  sqlite3* cursor;
  sqlite3_open(":memory:", &cursor);
  string command = "SELECT json_extract('" + json + "', '$." + key + "');";
  return SQLexecuteCommandSCALAR(cursor, command);
}

//extractJsonValueAflow///////////////////////////////////////////////////////
// This function extracts values from an aflowlib.json file. It is much faster
// than using SQLite's JSON extension, but has was designed to only work for
// the aflowlib.json. It cannot handle nested JSONs!
string extractJsonValueAflow(const string& json, string key) {
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

}  // namespace aflowlib

/********************************** SQLITE **********************************/

// These functions provide a direct interface to SQLite. They execute commands
// and handle callbacks.

namespace aflowlib {

// Execute command functions -------------------------------------------------
// These functions provide a framework to execute SQLite commands (including
// exception handling). The return types should all be void or string-typed to
// keep the number of functions to a minimum. Conversion to other data types
// should be handled outside these functions.

void SQLexecuteCommand(sqlite3* cursor, const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << "aflowlib::SQLexecuteCommand(): command = " << command << std::endl;
  char* sqlErrMsg = 0;
  int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallback, 0, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = "aflowlib::SQLexecuteCommand()";
    string message = string(sqlErrMsg) + " in command " + command;
    throw aurostd::xerror(function, message, _RUNTIME_SQL_);
  }
}

string SQLexecuteCommandSCALAR(sqlite3* cursor, const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << "aflowlib::SQLexecuteCommandSCALAR(): command = " << command << std::endl;
  char* sqlErrMsg = 0;
  string returnstring;
  int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallbackSCALAR, &returnstring, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = "aflowlib::SQLexecuteCommandSCALAR()";
    string message = string(sqlErrMsg) + " in command " + command;
    message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
    throw aurostd::xerror(function, message, _RUNTIME_SQL_);
  } else {
    return returnstring;
  }
}

vector<string> SQLexecuteCommandVECTOR(sqlite3* cursor, const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << "aflowlib::SQLexecuteCommandVECTOR(): command = " << command << std::endl;
  char *sqlErrMsg = 0;
  vector<string> returnvector;
  int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallbackVECTOR, &returnvector, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = "aflowlib::SQLexecuteCommandVECTOR()";
    string message = string(sqlErrMsg) + " in command " + command;
    message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
    throw aurostd::xerror(function, message, _RUNTIME_SQL_);
  } else {
    return returnvector;
  }
}

vector<vector<string> > SQLexecuteCommand2DVECTOR(sqlite3* cursor, const string& command) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG) && _SQL_COMMAND_DEBUG_);
  if (LDEBUG) std::cerr << "aflowlib::SQLexecuteCommand2DVECTOR(): command = " << command << std::endl;
  char *sqlErrMsg = 0;
  vector<vector<string> > returnvector;
  int sql_code = sqlite3_exec(cursor, command.c_str(), SQLcallback2DVECTOR, &returnvector, &sqlErrMsg);
  if (sql_code != SQLITE_OK) {
    string function = "aflowlib::SQLexecuteCommand2DVECTOR()";
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

int SQLcallback(void* data, int argc, char** argv, char** col) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG) && _SQL_CALLBACK_DEBUG_);
  (void) data;  // To suppress compiler warnings
  if (LDEBUG) {
    for (int i = 0; i < argc; i++) {
      std::cerr << "aflowlib::SQLcallback()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
    }
    std::cerr << std::endl;
  }
  return 0;
}

int SQLcallbackSCALAR(void* data, int argc, char** argv, char** col) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG) && _SQL_CALLBACK_DEBUG_);
  if (LDEBUG) {
    for (int i = 0; i < argc; i++) {
      std::cerr << "aflowlib::SQLcallbackSCALAR()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
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

int SQLcallbackVECTOR(void* data, int argc, char** argv, char** col) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG) && _SQL_CALLBACK_DEBUG_);
  if (LDEBUG) {
    for (int i = 0; i < argc; i++) {
      std::cerr << "aflowlib::SQLcallbackVECTOR()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
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

int SQLcallback2DVECTOR(void* data, int argc, char** argv, char** col) {
  bool LDEBUG = ((FALSE || XHOST.DEBUG) && _SQL_CALLBACK_DEBUG_);
  if (LDEBUG) {
    for (int i = 0; i < argc; i++) {
      std::cerr << "aflowlib::SQLcallback2DVECTOR()[" << i << "]: " << col[i] << " = " << (argv[i]?argv[i]:"NULL") << std::endl;
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
