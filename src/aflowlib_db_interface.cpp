//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
//****************************************************************************

// A collection of functions that interface directly with AFLOW JSON files and
// sqlite3 databases.

#include "aflowlib.h"
#include "SQLITE/sqlite3.h"

#define _SQL_COMMAND_DEBUG_  false  // debug SQL commands that are sent - verbose output
#define _SQL_CALLBACK_DEBUG_ false  // debug SQL callbacks - extremely verbose output
#define _DB_INTERFACE_DEBUG_ true   // switched on for logging

using std::string;
using std::vector;

static const int _N_AUID_TABLES_ = 256;

/*********************************** JSON ***********************************/

namespace aflowlib {

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
  string value = "";
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
      value = json.substr(start, end - start);
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
  string returnstring = "";
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
