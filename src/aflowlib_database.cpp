//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *            Aflow MARCO ESTERS - Duke University 2019-2020               *
// *                                                                         *
//****************************************************************************

// Class for the AFLOW database, based on SQLite. It provides the framework
// for rebuilding and interacting with the database file using data in JSON
// format. This file is divided into the following sections:
//
// * Constructor
// * Temp file handling
// * Database builder
// * Database analyzer
// * User-level functions to interact with the database
//
// ---------------------------------- USAGE ----------------------------------
//
// To update the database, i.e. to rebuild the database only when new entries
// are available, use: aflow --update_database
//
// To force a database rebuild, use: aflow --rebuild_database
//
// Database statistics will be collected in each case. To only perform the
// analysis, use: aflow --analyze_database
//
// ----------------------------- REBUILD PROCESS -----------------------------
//
// The database rebuild process has the following steps:
//   * Check if the database needs to be rebuild. This can be triggered either
//     by the user, by new entries in the schema, or by new/updated dat files.
//     The dat files are collections of aflowlib.json files for a set of
//     AUIDs (i.e. aflow:00.jsonl, aflow:01.jsonl, etc.).
//   * Create a temporary database file and populate with data from these JSONs.
//     Rebuilding from scratch instead of incrementally adding into the existing
//     database protects the database from corruption and injection attacks.
//     Using a temporary database file prevents that the active database file
//     get corrupted due to an error in the build process.
//   * Compare temporary and current database file and copy over if necessary.
//     The database file only gets copied after successful completion of the
//     rebuild process. A rebuild is considered successful when the number of
//     entries and properties of the new database file is greater than or equal
//     to that of the old file.
//   * Analyze the database and output the database statistics.
//
// ------------------------------ RETURN CODES -------------------------------
//
// The database rebuild uses return codes to convey to outside processes why
// it terminated. The codes are listed below.

// Return code definitions

// Success
#define _AFLOW_DB_SUCCESS_            0  // Database successfully updated

// Codes that do not trigger a database analysis
#define _AFLOW_DB_NOT_UPDATED_       100  // Nothing to update or patch
#define _AFLOW_DB_NOT_COPIED_        101  // Temp database file was not copied
#define _AFLOW_DB_NOT_PATCHED_       102  // No patches applied even though files were given

// Codes that trigger a database analysis
#define _AFLOW_DB_PATCH_SKIPPED_     200  // Returned when rebuild was successful but no patch was applied
#define _AFLOW_DB_PATCH_INCOMPLETE_  201  // At least one patch was incomplete

#include "aflowlib.h"
#include "SQLITE/aflow_sqlite.h"

// Some parts are written within the C++0x support in GCC, especially std::thread,
// which is implemented in gcc 4.4 and higher. For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400
#define AFLOW_DB_MULTITHREADS_ENABLE
#include <thread>
#include <mutex>
std::mutex m;
#else
#warning "The multithread parts of AflowDB will not be included, since they need gcc 4.4 and higher (C++0x support)."
#endif

#define _AFLOW_DB_DEBUG_ false

using std::string;
using std::vector;

static const string _AFLOW_DB_ERR_PREFIX_ = "AflowDB::";
static const uint _DEFAULT_SET_LIMIT_ = 16; //DX20200317 - int -> uint
static const int _N_AUID_TABLES_ = 256;

/************************** CONSTRUCTOR/DESTRUCTOR **************************/

namespace aflowlib {

  // Open the database for read access
  AflowDB::AflowDB(const string& db_file, ostream& oss) : xStream(oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    free();
    database_file = db_file;
    if (LDEBUG) std::cerr << "AflowDB: reading database" << std::endl;
    open(SQLITE_OPEN_READONLY);
  }

  // Open the database for write access
  AflowDB::AflowDB(const string& db_file, const string& dt_path, const string& lck_file, ostream& oss) : xStream (oss) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    free();
    data_path = dt_path;
    database_file = db_file;
    lock_file = lck_file;
    if (LDEBUG) {
      std::cerr << "AflowDB: Database file: " << database_file << std::endl;
      std::cerr << "AflowDB: Data path: " << data_path << std::endl;
      std::cerr << "AflowDB:: Lock file: " << lock_file << std::endl;
    }
    open();
  }

  // Copy constructors
  AflowDB::AflowDB(const AflowDB& that) : xStream(*that.getOFStream(), *that.getOSS()) {
    copy(that);
  }

  AflowDB& AflowDB::operator=(const AflowDB& that) {
    copy(that);
    return *this;
  }


  void AflowDB::copy(const AflowDB& that) {
    if (this == &that) return;
    xStream::copy(that);
    data_path = that.data_path;
    database_file = that.database_file;
    lock_file = that.lock_file;
    // Two databases should never operate on a temporary
    // file or have write access at the same time.
    is_tmp = false;
    open(SQLITE_OPEN_READONLY);
  }

  // Destructor
  AflowDB::~AflowDB() {
    close();
    free();
  }

  void AflowDB::free() {
    data_path = "";
    database_file = "";
    lock_file = "";
    is_tmp = false;
  }

  void AflowDB::clear() {
    close();
    free();
  }

  // Opens the database file and creates the main cursor
  void AflowDB::open(int open_flags) {
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "open():";
    string message = "Opening " + database_file + ".";
    pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
    int sql_code = sqlite3_open_v2(database_file.c_str(), &db, open_flags, NULL); //DX20200319 - nullptr -> NULL
    if (sql_code != SQLITE_OK) {
      message = "Could not open database file " + database_file;
      message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
  }

  // Closes the database file and removes the main cursor
  void AflowDB::close() {
    if (isTMP()) closeTmpFile(false, true, true);
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "close():";
    string message = "Closing " + database_file + ".";
    pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
    int sql_code = sqlite3_close(db);
    if (sql_code != SQLITE_OK) {
      message = "Could not close database file " + database_file;
      message += " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
  }

}  // namespace aflowlib

/********************************* TMP FILE **********************************/

namespace aflowlib {

  //openTmpFile///////////////////////////////////////////////////////////////
  // Opens a temporary database file 
  void AflowDB::openTmpFile(int open_flags, bool copy_original) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "openTmpFile():";
    stringstream message;

    string tmp_file = database_file + ".tmp";
    message << "Opening " << tmp_file;
    pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);

    // Create a symbolic link to lock the tmp file creation process. Since creating
    // symbolic links is an atomic process, it can be used to prevent race conditions.
    string lock_link = lock_file + ".lnk";
    if (symlink(lock_file.c_str(), lock_link.c_str())) {
      message << "Could not create symbolic link to lock file " + lock_file + ": ";
      if (errno == EEXIST) {
        message << "Another process has created it already.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_BUSY_);
      } else {
        message << "Process exited with errno " << errno << ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
      }
    }

    // If there is a temporary database file, a rebuild process is either in
    // progress or failed.
    if (aurostd::FileExist(tmp_file)) {
      long int tm_tmp = aurostd::FileModificationTime(tmp_file);
      time_t t = std::time(NULL); //DX20200319 - nullptr -> NULL
      long int tm_curr = (long int) t;
      int pid = -1;
      bool del_tmp = false;
      if (LDEBUG) {
        std::cerr << function << " Found temporary database file with time stamp " << tm_tmp
          << " (current time: " << tm_curr << ")." << std::endl;
      }
      // Check if the rebuild process that is saved in the lock file
      // is still active.
      if (aurostd::FileExist(lock_file) && !aurostd::FileEmpty(lock_file)) {
        pid = aurostd::string2utype<int>(aurostd::file2string(lock_file));
        if (kill(pid, 0)) pid = -1;
      }
      if (tm_curr - tm_tmp < DEFAULT_AFLOW_DB_STALE_THRESHOLD) {
        if (pid > -1) del_tmp = false;
        else del_tmp = true;
      } else {
        del_tmp = true;
      }

      if (del_tmp) {
        if (LDEBUG) {
          message << "Temporary database file already exists, but process has become stale."
            << " Removing temporary files and killing outstanding process.";
          pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
        }
        if (pid > -1) {
          int killreturn = kill(pid, 9);
          if (killreturn) {
            aurostd::RemoveFile(lock_link);
            message << "Could not kill active process " << pid << "(errno =  " << errno << ").";
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
          }
        }
        aurostd::RemoveFile(tmp_file);
        aurostd::RemoveFile(tmp_file + "-journal");
      } else {
        aurostd::RemoveFile(lock_link);
        message << "Could not create temporary database file. File already exists and rebuild process is active.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_BUSY_);
      }
    }

    aurostd::string2file(aurostd::utype2string<int>(getpid()), lock_file);
    int sql_code = sqlite3_close(db);
    if (sql_code != SQLITE_OK) {
      aurostd::RemoveFile(lock_link);
      message << "Could not close main database file " + database_file << " (SQL code " << sql_code << ").";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }

    if (copy_original) {
      bool copied = aurostd::CopyFile(database_file, tmp_file);
      // Make sure the file can be modified
      if (copied) {
        aurostd::ChmodFile("664", tmp_file);
      } else {
        aurostd::RemoveFile(lock_link);
        message << "Unable to copy original database file.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
      }
    }

    sql_code = sqlite3_open_v2(tmp_file.c_str(), &db, open_flags, NULL); //DX20200319 - nullptr -> NULL
    if (sql_code != SQLITE_OK) {
      message << "Could not open tmp database file " + tmp_file << " (SQL code " << sql_code << ").";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }

    is_tmp = true;
    aurostd::RemoveFile(lock_link);
  }

  //closeTmpFile//////////////////////////////////////////////////////////////
  // Closes a temporary database file and overwrites the original if the
  // build was successful. Unless --rebuild_database was selected by the user,
  // this function tests whether the new database has less entries than the
  // old one. To avoid data loss, the original file should not be overwritten.
  // nocopy: just close the tmp file without copying
  bool AflowDB::closeTmpFile(bool force_copy, bool keep, bool nocopy) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    string tmp_file = database_file + ".tmp";
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "closeTmpFile():";
    stringstream message;
    message << "Closing " << tmp_file;
    pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);

    // To throw
    bool found_error = false;
    stringstream message_error;

    // If the tmp file is empty, the build failed, so all other tests will fail as well.
    if (aurostd::FileEmpty(tmp_file)) {
      message_error << "Temporary database file empty. Database will not be overwritten.";
      nocopy = true;
      found_error = true;
    }
    if (!found_error && !nocopy) {
      if ((int) getTables().size() != _N_AUID_TABLES_) {
        message_error << "Temporary database file has the wrong number of tables."
          << " Database file will not be copied.";
        nocopy = true;
        found_error = true;
      }
    }

    // If the user does not force a rebuild, test that the new database does not
    // result in less data than the old one. To build a database with less
    // properties, aflow needs to be run with the --rebuild_database option.
    uint nentries_tmp = 0, ncols_tmp = 0, nentries_old = 0, ncols_old = 0;
    if (!found_error && !force_copy && !nocopy) {
      vector<string> props, tables;
      tables = getTables();
      // Get number of properties
      props = getColumnNames(tables[0]);  // All tables have the same columns, so any table is good
      ncols_tmp = props.size();

      // Get number of entries
      props = getPropertyMultiTables("COUNT", tables, "*");
      for (int i = 0; i < _N_AUID_TABLES_; i++) nentries_tmp += aurostd::string2utype<uint>(props[i]);
    }

    int sql_code = sqlite3_close(db);
    if (sql_code != SQLITE_OK) {
      message_error << "Could not close tmp database file " + tmp_file
        << " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message_error, _FILE_ERROR_);
    }
    is_tmp = false;

    bool copied = false;
    if (found_error || nocopy) {
      copied = false;
    } else if (force_copy) {
      message << "Force copy selected. Database will be overwritten.";
      pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
      copied = aurostd::CopyFile(tmp_file, database_file);
    } else {
      if (!aurostd::FileEmpty(database_file)) {
        message << "Old database file found. Determining number of entries and properties to prevent data loss.";
        pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);

        sqlite3* main_db;
        sql_code = sqlite3_open_v2(database_file.c_str(), &main_db, SQLITE_OPEN_READONLY, NULL);
        if (sql_code != SQLITE_OK) {
          message_error << "Could not close main database file " + database_file
            << " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message_error, _FILE_ERROR_);
        }
        vector<string> props, tables;
        tables = getTables(main_db);

        // Get number of properties
        props = getColumnNames(main_db, tables[0]);
        ncols_old = props.size();

        // Get number of entries
        props = getPropertyMultiTables(main_db, "COUNT", tables, "*");
        for (int i = 0; i < _N_AUID_TABLES_; i++) nentries_old += aurostd::string2utype<uint>(props[i]);

        sql_code = sqlite3_close(main_db);
        if (sql_code != SQLITE_OK) {
          message_error << "Could not close main database file " + database_file
            << " (SQL code " + aurostd::utype2string<int>(sql_code) + ").";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message_error, _FILE_ERROR_);
        }
      } else {
        if (LDEBUG) std::cerr << function << "No old datbase found." << std::endl;
        ncols_old = 0;
        nentries_old = 0;
      }

      if (LDEBUG) {
        std::cerr << function << "Number of entries: " << nentries_tmp << " (current database: " << nentries_old << ")." << std::endl;
        std::cerr << function << "Number of properties: " << ncols_tmp << " (current database: " << ncols_old << ")." << std::endl;
      }

      if (nentries_tmp < nentries_old) {
        message_error << "The rebuild process resulted in less entries"
          << " than in the current database. To prevent accidental"
          << " data loss, the temporary database will not be copied."
          << " Rerun as aflow --rebuild_database if the database should"
          << " be updated regardless. The temporary database file will be"
          << " kept to allow debugging.";
        keep = true;
        copied = false;
        found_error = true;
      } else if (ncols_tmp < ncols_old) {
        message_error << "The rebuild process resulted in less properties"
          << " than in the current database. To prevent accidental"
          << " data loss, the temporary database will not be copied."
          << " Rerun as aflow --rebuild_database if the database should"
          << " be updated regardless. The temporary database file will be"
          << " kept to allow debugging.";
        keep = true;
        copied = false;
        found_error = true;
      } else if (!found_error) {
        // No problems found, so replace old database file with new file
        copied = aurostd::CopyFile(tmp_file, database_file);
      }
      message << "Database file " << (copied?"successfully":"not") << " copied.";
      pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
    }
    if (!found_error && !keep) {
      aurostd::RemoveFile(tmp_file);
      message << "Temporary database file removed.";
      pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
    } else if (LDEBUG) {
      message << "Temporary database file " << tmp_file << " will not be deleted.";
      pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
    }

    if (found_error) {
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }

    open();

    aurostd::string2file("", lock_file);
    if (copied) aurostd::ChmodFile("444", database_file);
    return copied;
  }

  //isTMP/////////////////////////////////////////////////////////////////////
  // Returns whether or not the cursor points to a temporary file or not
  bool AflowDB::isTMP() {
    return is_tmp;
  }

}  // namespace aflowlib

/***************************** REBUILD DATABASE *****************************/

namespace aflowlib {

  //rebuildDatabase///////////////////////////////////////////////////////////
  // This function first checks if a rebuild is necessary and then initiates
  // the rebuilding functions. Returns true if the database has been
  // successfully rebuilt.
  int AflowDB::rebuildDatabase(bool force_rebuild) {
    vector<string> patch_placeholder;
    return rebuildDatabase(patch_placeholder, force_rebuild);
  }

  int AflowDB::rebuildDatabase(const string& patch_files, bool force_rebuild) {
    vector<string> patch_files_input;
    if (!patch_files.empty()) aurostd::string2tokens(patch_files, patch_files_input, ",");
    return rebuildDatabase(patch_files_input, force_rebuild);
  }

  int AflowDB::rebuildDatabase(const vector<string>& patch_files_input, bool force_rebuild) {
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "rebuildDatabase():";
    stringstream message;
    bool rebuild_db = false;

    // Always rebuild when the user wants the rebuild.
    rebuild_db = force_rebuild;
    if (rebuild_db) {
      message << "Rebuilding database (user rebuild).";
      pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
    }

    // Always rebuild if the database file doesn't exist. Note: When the
    // database is opened and the database file doesn't exist, SQLite creates
    // an empty file, so check if the file is empty and not if it exists.
    if (!rebuild_db) {
      rebuild_db = aurostd::FileEmpty(database_file); // file_empty needed for LDEBUG
      if (rebuild_db) {
        message << "Rebuilding database (file not found or empty).";
        pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
      }
    }

    // Rebuild when the schema has been updated
    if (!rebuild_db) {
      vector<string> columns;
      string table = getTables()[0];  // All tables have the same columns, so any table is good
      columns = getColumnNames(table);

      // Properties
      vector<string> keys_schema = getSchemaKeys();
      uint nkeys = keys_schema.size();
      if (nkeys > columns.size()) {
        rebuild_db = true;
      } else if (nkeys < columns.size()) {
        message << "The number of properties (" << nkeys << ") is smaller than in the"
          << " current database (" << columns.size() << "). To prevent accidental"
          << " data loss, AFLOW will abort the update process. Rerun as"
          << " aflow --rebuild_database if the database should be updated regardless.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
      }

      // Data types
      if (!rebuild_db) {
        vector<string> types_db = getColumnTypes(table);
        vector<string> types_schema = getDataTypes(keys_schema, false);
        int index = -1;
        string key = "";
        uint k = 0;
        for (k = 0; k < nkeys; k++) {
          key = XHOST.vschema.getattachedscheme("SCHEMA::NAME:" + keys_schema[k]);
          if (!aurostd::WithinList(columns, key, index)) break;
          if (types_db[index] != types_schema[k]) break;
        }
        rebuild_db = (k != nkeys);
      }
      if (rebuild_db) {
        message << " Rebuilding database (schema updated).";
        pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
      }
    }

    // Check if any relevant files are newer than the database.
    if (!rebuild_db) {
      vector<string> json_files(_N_AUID_TABLES_);
      for (int i = 0; i < _N_AUID_TABLES_; i++) {
        stringstream t;
        t << std::setfill('0') << std::setw(2) << std::hex << i;
        json_files[i] = aurostd::CleanFileName(data_path + "/aflow:" + t.str() + ".jsonl");
        if (!aurostd::EFileExist(json_files[i]) && !aurostd::FileExist(json_files[i])) {
          message << data_path << " is not a valid data path. Missing file for aflow:" << t.str() << ".";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_NOT_FOUND_);
        }
      }

      long int tm_db = aurostd::FileModificationTime(database_file);
      int i = 0;
      for (i = 0; i < _N_AUID_TABLES_; i++) {
        if (aurostd::FileModificationTime(json_files[i]) > tm_db) break;
      }
      rebuild_db = (i != _N_AUID_TABLES_);

      if (rebuild_db) message << "Rebuilding database (updated files found).";
      else message << "No updated files found. Database will not be rebuilt.";
      pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
    }

    // Check if the database may need to be patched
    uint npatch_input = patch_files_input.size();

    if (rebuild_db) {
      openTmpFile();
      rebuildDB();
      bool close_tmp = closeTmpFile(force_rebuild);
      // If the tmp database was not copied, something went wrong - return 1
      if (!close_tmp) return _AFLOW_DB_NOT_COPIED_;

      // No patches needed, so return
      if (npatch_input == 0) return (close_tmp?_AFLOW_DB_SUCCESS_:_AFLOW_DB_NOT_COPIED_);
    }

    if (npatch_input > 0) {
      // Do not check timestamps after a rebuild because patches must be applied
      int patch_code = patchDatabase(patch_files_input, !rebuild_db);
      if (rebuild_db) {
        if (patch_code == _AFLOW_DB_NOT_PATCHED_) return _AFLOW_DB_PATCH_SKIPPED_; // Database was rebuilt but not patched
        else return patch_code;
      } else if (patch_code == _AFLOW_DB_NOT_PATCHED_) {
        return _AFLOW_DB_NOT_UPDATED_;  // Database not rebuilt and not patched
      } else {
        return patch_code;
      }
    }

    return _AFLOW_DB_NOT_UPDATED_;
  }

  //patchDatabase/////////////////////////////////////////////////////////////
  // Patches the database using jsonl files. This can be used for emergency
  // patches before a lib2raw run or to add data that lib2raw does not cover.
  int AflowDB::patchDatabase(const string& patch_files, bool check_timestamps) {
    vector<string> patch_files_input;
    if (!patch_files.empty()) aurostd::string2tokens(patch_files, patch_files_input, ",");
    return patchDatabase(patch_files_input, check_timestamps);
  }

  int AflowDB::patchDatabase(const vector<string>& patch_files_input, bool check_timestamps) {
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "patchDatabase():";
    stringstream message;
    int patch_code = 0;

    vector<string> patch_files;
    uint npatch_input = patch_files_input.size();
    if (npatch_input == 0) return 0;  // No files, so nothing to patch

    long int tm_db = 0;
    if (check_timestamps) tm_db = aurostd::FileModificationTime(database_file);

    // Check files that need to be patched
    for (uint i = 0; i < npatch_input; i++) {
      string filename = patch_files_input[i];
      if (!aurostd::EFileExist(filename)) {
        // Check if the file is in the data path
        if (aurostd::EFileExist(data_path + "/" + filename)) {
          filename = aurostd::CleanFileName(data_path + "/" + filename);
        } else {
          message << "File " << filename << " not found. Skipping.";
          pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);

          // If a file is missing, then the patch is incomplete
          patch_code = _AFLOW_DB_PATCH_INCOMPLETE_;
          continue;
        }
      }

      if (!check_timestamps || (aurostd::FileModificationTime(filename) > tm_db)) {
        patch_files.push_back(filename);
        message << "Adding file " << filename << " to patch list.";
      } else {
        message << "Skipping file " << filename << ". File is older than the database.";
      }
      pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
    }

    uint npatch = patch_files.size();
    if (npatch == 0) {  // Nothing to patch
      if (check_timestamps) return _AFLOW_DB_NOT_PATCHED_;  // No new files
      else return _AFLOW_DB_PATCH_INCOMPLETE_;  // All (new) files missing
    }

    // Start patching
    uint patches_applied = 0;
    uint entries_patched = 0;
    for (uint i = 0; i < npatch; i++) {
      message << "Patching database using " << patch_files[i] << ".";
      pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
      openTmpFile(SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, true);  // Copy original

      // Do not constantly synchronize the database with file on disk (increases
      // performance significantly).
      sql::SQLexecuteCommand(db, "PRAGMA synchronous = OFF");

      vector<string> data;
      aurostd::efile2vectorstring(patch_files[i], data);
      entries_patched = applyPatchFromJsonl(data);
      if (entries_patched > 0) {
        if (closeTmpFile(false)) {  // Do not force copy or a bad patch may break the database
          message << "Patched " << entries_patched << "/" << data.size() << " entries.";
          pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
          // Patch incomplete, so adjust code (unless an error was found before)
          if ((patch_code < 200) && (entries_patched < data.size())) patch_code = _AFLOW_DB_PATCH_INCOMPLETE_;
          patches_applied++;
        } else {
          message << "Could not close temporary database file.";
          pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        }
      } else {
        closeTmpFile(false, false, true);  // Do not copy tmp file
        message << "No patches applied.";
        // Patch incomplete, so adjust code (unless an error was found before)
        if ((patch_code < 200) && (entries_patched < data.size())) patch_code = _AFLOW_DB_PATCH_INCOMPLETE_;
        pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      }
    }

    if (patches_applied == 0) patch_code = _AFLOW_DB_NOT_PATCHED_;

    return patch_code;
  }

  // Rebuild -----------------------------------------------------------------

  //rebuildDB/////////////////////////////////////////////////////////////////
  // Rebuilds the database from scratch.
  void AflowDB::rebuildDB() {
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "rebuildDB():";
    stringstream message;

    message << "Starting rebuild.";
    pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);

    // Do not constantly synchronize the database with file on disk (increases
    // performance significantly).
    sql::SQLexecuteCommand(db, "PRAGMA synchronous = OFF");

    // Get columns and types from schema
    vector<string> keys = getSchemaKeys();
    uint nkeys = keys.size();
    vector<string> columns(nkeys);
    for (uint k = 0; k < nkeys; k++) columns[k] = XHOST.vschema.getattachedscheme("SCHEMA::NAME:" + keys[k]);
    vector<string> types = getDataTypes(columns, true);

    // Rebuild
#ifdef AFLOW_DB_MULTITHREADS_ENABLE
    int ncpus = init::GetCPUCores();
    int max_cpu = 16;
    if (ncpus < 1) ncpus = 1;
    if (ncpus > max_cpu) ncpus = max_cpu;
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
    message << function << "Finished rebuild.";
    pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);
  }

  //buildTables///////////////////////////////////////////////////////////////
  // Reads the .jsonl files and processes the JSONs for the database writer.
  void AflowDB::buildTables(int startIndex, int endIndex, const vector<string>& columns, const vector<string>& types) {
    for (int i = startIndex; i < endIndex; i++) {
      stringstream t;
      t << std::setfill('0') << std::setw(2) << std::hex << i;
      string table = "auid_" + t.str();
      createTable(table, columns, types);

      string jsonfile = aurostd::CleanFileName(data_path + "/aflow:" + t.str() + ".jsonl");
      vector<string> data;
      aurostd::efile2vectorstring(jsonfile, data);
      uint ndata = data.size();
      vector<vector<string> > values(ndata);
      for (uint d = 0; d < ndata; d++) {
        values[d] = getDataValues(data[d], columns, types);
      }
      populateTable(table, columns, values);
    }
  }

  //populateTable/////////////////////////////////////////////////////////////
  // Populates the database tables and creates indexes. This function uses a
  // mutex to make the writing thread-safe. While SQLITE does not allow
  // concurrent writing anyway, multiple threads may open a transaction, which
  // causes the build to fail (hence the mutex). The mutex is also the reason
  // why this function should never do any processing.
  void AflowDB::populateTable(const string& table, const vector<string>& columns, const vector<vector<string> >& values) {
    bool LDEBUG = (FALSE || XHOST.DEBUG || _AFLOW_DB_DEBUG_);
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "populateTable():";
#ifdef AFLOW_DB_MULTITHREADS_ENABLE
    std::unique_lock<std::mutex> lk(m);
#endif
    int chunk_size = 1000, count = 0;
    transaction(true);
    for (uint v = 0, nvals = values.size(); v < nvals; v++) {
      insertValues(table, columns, values[v]);
      if (++count % chunk_size == 0 || v == nvals - 1) {
        transaction(false);
        transaction(true);
        count = 0;
      }
    }

    // Create indexes on important database properties
    vector<string> index_cols;
    string indexes = "auid,aurl,catalog,compound,nspecies,Pearson_symbol_relax,spacegroup_relax";
    aurostd::string2tokens(indexes, index_cols, ",");
    string index, index_expression;
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
    transaction(false);
    if (LDEBUG) std::cerr << function <<  " Finished building table " << table << std::endl;
  }

  // Patch -------------------------------------------------------------------

  //applyPatchFromJsonl///////////////////////////////////////////////////////
  // Applies database patches from data in jsonl format. Returns the number
  // of entries that were patched.
  uint AflowDB::applyPatchFromJsonl(const vector<string>& data) {
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "applyPatchFromJsonl():";
    stringstream message;

    uint nlines = data.size();
    if (nlines == 0) return 0;  // Nothing to do

    uint chunk_size = 1000;

    vector<string> schema_keys = getSchemaKeys();

    // Patch
    vector<string> keys, vals, cols, types;
    string auid = "";
    uint npatches = 0, nkeys = 0;
    transaction(true);
    for (uint l = 0; l < nlines; l++) {
      cols.clear();
      keys = extractJsonKeysAflow(data[l]);
      nkeys = keys.size();

      // Check
      for (uint k = 0; k < nkeys; k++) {
        if (keys[k] == "auid") {
          auid = extractJsonValueAflow(data[l], keys[k]);
          auid = auid.substr(1, auid.size() - 2);  // Remove quotes
        } else if (aurostd::WithinList(schema_keys, keys[k])) {
          cols.push_back(keys[k]);
        }
      }

      if (auid.empty()) {
        message << "Skipping line " << l << ". No AUID found in JSON.";
        pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        continue;
      }

      if (!auidInDatabase(auid)) {
        message << "Skipping line " << l << ". AUID " << auid << " not in database.";
        pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        continue;
      }

      if (cols.size() == 0) {
        message << "Skipping line " << l << ". No keys match properties in schema.";
        pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        continue;
      }

      types = getDataTypes(cols, false);
      vals = getDataValues(data[l], cols, types);

      // Update
      try {
        updateEntry(auid, cols, vals);
        npatches++;
      } catch (aurostd::xerror& e) {
        message << "Failed to patch using line " << l << " (SQL error): " << e.error_message << ".";
        pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        continue;
      }

      if ((npatches + 1) % chunk_size == 0) {
        transaction(false);
        transaction(true);
      }
    }

    transaction(false);

    return npatches;
  }

  //auidInDatabase////////////////////////////////////////////////////////////
  // Test if an AUID is in the database by retrieving a value of the row (here
  // AUID since it is indexed).
  bool AflowDB::auidInDatabase(const string& auid) {
    string table = "auid_" + auid.substr(6, 2);
    string where = "auid='\"" + auid + "\"'";
    string result = getValue(table, "auid", where);
    return (!result.empty());
  }

  //updateEntry///////////////////////////////////////////////////////////////
  // Updates columns of an entry with a specific AUID.
  void AflowDB::updateEntry(const string& auid, const vector<string>& cols, const vector<string>& vals) {
    string table = "auid_" + auid.substr(6, 2);
    string where = "AUID='\"" + auid + "\"'";
    updateRow(table, cols, vals, where);
  }

  // Schema ------------------------------------------------------------------

  //getSchemaKeys/////////////////////////////////////////////////////////////
  // Returns the keys from the AFLOW schema.
  vector<string> AflowDB::getSchemaKeys() {
    vector<string> keys;
    string key = "";
    for (uint i = 0, n = XHOST.vschema.vxsghost.size(); i < n; i += 2) {
      if(XHOST.vschema.vxsghost[i].find("::NAME:") != string::npos) {
        key=aurostd::RemoveSubString(XHOST.vschema.vxsghost[i], "SCHEMA::NAME:");
        // json keys are lower case
        keys.push_back(aurostd::tolower(key));
      }
    }
    return keys;
  }

  // Data --------------------------------------------------------------------

  //getDataTypes//////////////////////////////////////////////////////////////
  // Gets the data types of the schema keys and converts them into SQLite
  // types. Note that SQLite does not recognize arrays or Booleans, so they
  // will be stored as text.
  vector<string> AflowDB::getDataTypes(const vector<string>& keys, bool unique) {
    uint nkeys = keys.size();
    vector<string> types(nkeys);
    string type = "";
    for (uint k = 0; k < nkeys; k++) {
      // AUID has to be unique
      type = XHOST.vschema.getattachedscheme("SCHEMA::TYPE:" + aurostd::toupper(keys[k]));
      if (unique && (keys[k] == "AUID")) {
        types[k] = "TEXT UNIQUE NOT NULL";
      } else if (type == "number") {
        types[k] = "REAL";
      } else {
        types[k] = "TEXT";
      }
    }
    return types;
  }

  //getDataValues/////////////////////////////////////////////////////////////
  // Retrieves the values of each property from the aflowlib.json file.
  vector<string> AflowDB::getDataValues(const string& entry, const vector<string>& cols, const vector<string>& types) {
    string value = "", id = "";
    uint ncols = cols.size();
    vector<string> values(ncols, "NULL");
    for (uint c = 0; c < ncols; c++) {
      value = extractJsonValueAflow(entry, cols[c]);
      if (!value.empty() && (aurostd::toupper(value) != "NULL")) {
        if (types[c] != "REAL") values[c] = "'" + value + "'";
        else values[c] = value;
      }
    }
    return values;
  }

  //extractJsonKeysAflow//////////////////////////////////////////////////////
  // This function extracts keys from an aflowlib.json file. It is much
  // faster than using SQLite's JSON extension, but has was designed to only
  // work for the aflowlib.json. It cannot handle nested JSONs!
  vector<string> AflowDB::extractJsonKeysAflow(const string& json) {
    vector<string> keys;
    string substring = "";
    string::size_type pos = 0, lastPos = 0, dpos = 0, quote1  = 0, quote2 = 0, colon = 0;

    // Find the first comma - this is either the end of the key-value pair
    // or part of an array. Either way, the key is inside.
    pos = json.find(",");
    lastPos = 1;  // First character is a curly brace, so skip
    dpos = pos - lastPos;
    while ((pos != string::npos) || (lastPos != string::npos)) {
      // A comma could be separating a key-value pair an array
      // or numbers or strings
      substring = json.substr(lastPos, dpos);

      // Find the colon - if there is no colon, it cannot be a key-value pair
      colon = substring.find(":");
      if (colon != string::npos) {
        // A key is enclosed in quotes, so there must be at least two of them
        quote1 = substring.find("\"");
        if (quote1 != string::npos) {
          quote2 = substring.find("\"", quote1 + 1);
          // Most non-keys are filtered out by now. There could still be array
          // elements left. In that case, however, the colon is between the quotes,
          // so make sure that the first two quotes appear before the colon and
          // take everything in-between as the key. This breaks if quotes, colons,
          // and commas are inside a string in the right sequence, but should not
          // be the case in AFLOW's JSON files.
          if ((quote2 != string::npos) && (quote1 < colon) && (quote2 < colon)) {
            substring = substring.substr(quote1 + 1, quote2 - quote1 - 1);
            if (!substring.empty()) keys.push_back(substring);
          }
        }
      }
      // Move on to the next comma
      lastPos = json.find_first_not_of(",", pos);
      pos = json.find(",", lastPos);
      dpos = pos - lastPos;
    }
    return keys;
  }

  //extractJsonValueAflow/////////////////////////////////////////////////////
  // This function extracts values from an aflowlib.json file. It is much
  // faster than using SQLite's JSON extension, but has was designed to only
  // work for the aflowlib.json. It cannot handle nested JSONs!
  string AflowDB::extractJsonValueAflow(const string& json, string key) {
    string value = "";
    key = "\"" + key + "\":";
    string::size_type start = 0, end = 0;
    start = json.find(key);
    if (start != string::npos) {
      start += key.length();
      end = json.find("\":", start);
      if (end != string::npos) {
        // In case there is any white space between key and value
        value = aurostd::RemoveWhiteSpacesFromTheFront(json.substr(start, end - start));
        // If we have a nested object, "value" should only be '{' + white space by now.
        if (value[0] == '{') {
          string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "extractJsonValueAflow():";
          string message = "JSON parser cannot read nested objects.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
        }
        end = value.find_last_of(",");
        value = value.substr(0, end);
      } else {
        end = json.find("}", start);
        // In case there is any white space between key and value
        value = aurostd::RemoveWhiteSpacesFromTheFront(json.substr(start, end - start));
        // If we have a nested object, it should start with '{'
        if (value[0] == '{') {
          string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "extractJsonValueAflow():";
          string message = "JSON parser cannot read nested objects.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
        }
      }
    } else {
      value = "";
    }
    return value;
  }

}  // namespace aflowlib

/**************************** DATABASE ANALYSIS *****************************/

namespace aflowlib {

  //analyzeDatabase///////////////////////////////////////////////////////////
  // Provides analytics for the database in JSON format.
  void AflowDB::analyzeDatabase(const string& outfile) {

    if (aurostd::FileEmpty(database_file)) {
      string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "analyzeDatabase():";
      string message = "Cannot analyze database. File empty.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_CORRUPT_);
    }

    // Get properties and tables for which statistics need to be collected
    // Since all tables have the same columns, only one table needs to be searched
    vector<string> tables = getTables();
    vector<string> columns = getColumnNames(tables[0]);

    vector<string> catalogs = getSetMultiTables(tables, "catalog", true);
    uint ncatalogs = catalogs.size();

    vector<string> loop_entries, loops;
    loop_entries = getSetMultiTables(tables, "loop", true);
    loops = getUniqueFromJsonArrays(loop_entries);

    vector<DBStats> db_stats(ncatalogs + 1);
    DBStats& total_stats = db_stats[0];  // Declare to make code more legible
    total_stats = initDBStats("\"total\"", columns, loops);

    for (uint c = 0; c < ncatalogs; c++) {
      DBStats& catalog_stats = db_stats[c+1];  // Declare to make code more legible
      catalog_stats = getCatalogStats(catalogs[c], tables, columns, loops);

      // Add to totals
      total_stats.nentries += catalog_stats.nentries;
      // Since we cannot determine which systems are in the ICSD but not
      // in other catalogs, skip it for the system count. The error should
      // be small though.
      if (catalog_stats.catalog != "ICSD") total_stats.nsystems += catalog_stats.nsystems;

      // Columns
      for (uint i = 0; i < total_stats.columns.size(); i++) {
        total_stats.count[i][0] += catalog_stats.count[i][0];
        total_stats.count[i][1] += catalog_stats.count[i][1];
        if ((catalog_stats.types[i] != "bool") && (catalog_stats.count[i][0] + catalog_stats.count[i][1] > 0)) {
          if (total_stats.types[i] == "number") {
            if (total_stats.max[i].empty()) total_stats.max[i] = catalog_stats.max[i];
            else if (aurostd::string2utype<double>(catalog_stats.max[i]) > aurostd::string2utype<double>(total_stats.max[i])) total_stats.max[i] = catalog_stats.max[i];

            if (total_stats.min[i].empty()) total_stats.min[i] = catalog_stats.min[i];
            else if (aurostd::string2utype<double>(total_stats.min[i]) > aurostd::string2utype<double>(catalog_stats.min[i])) total_stats.min[i] = catalog_stats.min[i];
          }

          // Just append to the sets for now and sort later
          uint nset = total_stats.set[i].size();
          for (uint s = 0; s < catalog_stats.set[i].size(); s++) {
            if ((nset <= _DEFAULT_SET_LIMIT_) && !aurostd::WithinList(total_stats.set[i], catalog_stats.set[i][s])) {
              total_stats.set[i].push_back(catalog_stats.set[i][s]);
              nset++;
            }
          }
        }
      }

      // Loops
      for (uint i = 0; i < total_stats.loop_counts.size(); i++){
        total_stats.loop_counts[i].second += catalog_stats.loop_counts[i].second;
      }


      // Species
      for (uint s = 0; s < catalog_stats.species.size(); s++) {
        if (!aurostd::WithinList(total_stats.species, catalog_stats.species[s])) {
          total_stats.species.push_back(catalog_stats.species[s]);
        }
      }
    }

    // Get the sets for the totals
    for (uint i = 0; i < total_stats.columns.size(); i++) {
      uint nset = total_stats.set[i].size();
      if (nset <= _DEFAULT_SET_LIMIT_) {
        if (total_stats.types[i] == "bool") {
          total_stats.set[i].clear();
          total_stats.set[i].push_back("true");
          total_stats.set[i].push_back("false");
        } else if (total_stats.types[i] == "number") {
          vector<double> set_dbl(nset);
          for (uint s = 0; s < nset; s++) set_dbl[s] = aurostd::string2utype<double>(total_stats.set[i][s]);
          aurostd::sort(set_dbl, total_stats.set[i]);
        } else {
          std::sort(total_stats.set[i].begin(), total_stats.set[i].end());
        }
      }
    }
    std::sort(total_stats.species.begin(), total_stats.species.end());

    // Write everything now
    ncatalogs = ncatalogs + 1; // Include totals now

    string tab = "    ";
    std::stringstream json;
    json << "{" << std::endl << tab << "\"Aflow_DBs\": {" << std::endl;

    for (uint c = 0; c < ncatalogs; c++) {
      writeStatsToJson(json, db_stats[c]);
      if (c < ncatalogs - 1) json << ",";
      json << std::endl;
    }
    json << tab << "}" << std::endl << "}" << std::endl;
    aurostd::stringstream2file(json, outfile);
  }

  DBStats AflowDB::initDBStats(const string& catalog, const vector<string>& cols, const vector<string>& loops) {
    uint ncols = cols.size();
    uint nloops = loops.size();

    DBStats stats;
    stats.catalog = catalog;
    stats.columns = cols;
    stats.count.assign(ncols, vector<int>(2, 0));
    stats.max.assign(ncols, "");
    stats.min.assign(ncols, "");
    stats.nentries = 0;
    stats.nsystems = 0;
    stats.set.resize(ncols);
    stats.loop_counts.resize(nloops);
    for (uint l = 0; l < nloops; l++) {
      stats.loop_counts[l].first = loops[l];
      stats.loop_counts[l].second = 0;
    }

    // Get types for post-processing
    vector<string> types(ncols);
    for (uint i = 0; i < ncols; i++) {
      types[i] = XHOST.vschema.getattachedscheme("SCHEMA::TYPE:" + aurostd::toupper(cols[i]));
    }
    stats.types = types;

    return stats;
  }

  //getCatalogStats///////////////////////////////////////////////////////////
  // Gets the statistics for all properties in the catalog.
  DBStats AflowDB::getCatalogStats(const string& catalog, const vector<string>& tables,
      const vector<string>& cols, const vector<string>& loops) {
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "getCatalogStats():";
    stringstream message;

    uint ncols = cols.size();
    uint nloops = loops.size();

    DBStats stats = initDBStats(catalog, cols, loops);
    const vector<string>& types = stats.types;

    string where = "catalog='" + catalog + "'";
    vector<string> entries = getPropertyMultiTables("COUNT", tables, "*", where);
    for (int t = 0; t < _N_AUID_TABLES_; t++) stats.nentries += aurostd::string2utype<int>(entries[t]);
    message << "Starting analysis for catalog " << catalog << " (" << stats.nentries << " entries).";
    pflow::logger(_AFLOW_FILE_NAME_, function, message, *p_FileMESSAGE, *p_oss);

    if (stats.nentries > 0) {
      vector<vector<vector<string> >  > maxmin(_N_AUID_TABLES_, vector<vector<string> >(ncols, vector<string>(2))), sets(_N_AUID_TABLES_, vector<vector<string> >(ncols));
      vector<vector<vector<int> > > counts(_N_AUID_TABLES_, vector<vector<int> >(ncols, vector<int>(2, 0)));
      vector<vector<int> > loop_counts(_N_AUID_TABLES_, vector<int>(nloops));
#ifdef AFLOW_DB_MULTITHREADS_ENABLE
      int ncpus = init::GetCPUCores();
      if (ncpus < 1) ncpus = 1;
      // The maximum number of CPUs are empirically found values that appear
      // to result in the shortest run times. Further testing may be necessary.
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
              std::ref(cols), std::ref(types), std::ref(loops),
              std::ref(counts), std::ref(loop_counts),
              std::ref(maxmin), std::ref(sets)));
      }
      for (int i = 0; i < ncpus; i++) {
        threads[i]->join();
        delete threads[i];
      }
#else
      getColStats(0, _N_AUID_TABLES_, catalog, tables, cols, types,
          loops, counts, loop_counts, maxmin, sets);
#endif

      // Properties: count, max, min, set
      vector<string> set;
      string max = "", min = "";
      uint nset = 0, n = 0;
      for (uint c = 0; c < ncols; c++) {
        for (int t = 0; t < _N_AUID_TABLES_; t++) {
          stats.count[c][0] += counts[t][c][0];
          stats.count[c][1] += counts[t][c][1];
        }
        if (stats.count[c][0] + stats.count[c][1] > 0) {
          set.clear();
          max = ""; min = "";
          nset = 0; n = 0;
          if (types[c] != "bool") {  // No max, min, or set for bool
            for (int t = 0; t < _N_AUID_TABLES_; t++) {
              if (counts[t][c][0] > 0) {
                if (types[c] == "number") {
                  if (max.empty()) max = maxmin[t][c][0];
                  else if (aurostd::string2utype<double>(maxmin[t][c][0]) > aurostd::string2utype<double>(max)) max = maxmin[t][c][0];

                  if (min.empty()) min = maxmin[t][c][1];
                  else if (aurostd::string2utype<double>(maxmin[t][c][1]) < aurostd::string2utype<double>(min)) min = maxmin[t][c][1];
                }
                if (nset <= _DEFAULT_SET_LIMIT_) {
                  n = sets[t][c].size();
                  if (n > _DEFAULT_SET_LIMIT_) {
                    set = sets[t][c];
                    nset = n;
                  } else {
                    for (uint i = 0; i < n; i++) {
                      if (!aurostd::WithinList(set, sets[t][c][i])) {
                        set.push_back(sets[t][c][i]);
                        nset++;
                      }
                      if (nset > _DEFAULT_SET_LIMIT_) break;
                    }
                  }
                }
              }
            }
          }
          stats.max[c] = max;
          stats.min[c] = min;
          if (types[c] == "bool") {
            set.clear();
            set.push_back("true");
            set.push_back("false");
          } else if (nset <= _DEFAULT_SET_LIMIT_) {
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

      // Species, systems
      vector<string> species = getSetMultiTables(tables, "species", true, where);
      stats.nsystems = species.size();
      stats.species = getUniqueFromJsonArrays(species);

      // Loop counts
      for (uint l = 0; l < nloops; l++) {
        for (int t = 0; t < _N_AUID_TABLES_; t++) stats.loop_counts[l].second += loop_counts[t][l];
      }
    }

    return stats;
  }

  //getColStats///////////////////////////////////////////////////////////////
  // Retrieves the statistics for each database property and the loops.
  void AflowDB::getColStats(int startIndex, int endIndex, const string& catalog,
      const vector<string>& tables, const vector<string>& cols, const vector<string>& types,
      const vector<string>& loops, vector<vector<vector<int> > >& counts, vector<vector<int> >& loop_counts,
      vector<vector<vector<string> > >& maxmin, vector<vector<vector<string> > >& sets) {
    sqlite3* cursor;
    string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "getColStats():";
    string message = "";
    int sql_code = sqlite3_open_v2(database_file.c_str(), &cursor, SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX, NULL); //DX20200319 - nullptr -> NULL
    if (sql_code != SQLITE_OK) {
      message = "Could not open cursor on database file " + database_file + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }

    uint ncols = cols.size();
    string where = "", where_loop = "";
    for (int i = startIndex; i < endIndex; i++) {
      for (uint c = 0; c < ncols; c++) {
        if (types[c] == "bool") {
          where = "catalog='" + catalog + "' AND " + cols[c] + "=";
          counts[i][c][0] = aurostd::string2utype<int>(getProperty(cursor, "COUNT", tables[i], cols[c], where + "'true'"));
          counts[i][c][1] = aurostd::string2utype<int>(getProperty(cursor, "COUNT", tables[i], cols[c], where + "'false'"));
        } else {
          where = "catalog='" + catalog + "' AND " + cols[c] + " NOT NULL";
          counts[i][c][0] = aurostd::string2utype<int>(getProperty(cursor, "COUNT", tables[i], cols[c], where));
        }
        // No need to determine max, min, or set for bool
        if ((types[c] != "bool") && (counts[i][c][0] > 0)) {
          // Max and  min only make sense for numbers
          if (types[c] == "number") {
            maxmin[i][c][0] = getProperty(cursor, "MAX", tables[i], cols[c], where);
            maxmin[i][c][1] = getProperty(cursor, "MIN", tables[i], cols[c], where);
          }
          sets[i][c] = getSet(cursor, tables[i], cols[c], true, where, _DEFAULT_SET_LIMIT_ + 1);
        }
      }
      for (uint l = 0, nloops = loops.size(); l < nloops; l++) {
        where_loop = where + " AND loop LIKE '%" + loops[l] + "%'";
        loop_counts[i][l] = aurostd::string2utype<int>(getProperty("COUNT", tables[i], "loop", where_loop));
      }
    }

    sql_code = sqlite3_close(cursor);
    if (sql_code != SQLITE_OK) {
      message = "Could not close cursor on database file " + database_file + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
  }

  //getUniqueFromJsonArrays///////////////////////////////////////////////////
  // Determines the unique array elements in a set of 1D-array strings.
  vector<string> AflowDB::getUniqueFromJsonArrays(const vector<string>& arrays) {
    vector<string> unique, tokens;
    string arr = "";
    int nunique = 0, u = 0;
    for (uint a = 0; a < arrays.size(); a++) {
      vector<string> tokens;
      arr = arrays[a];
      arr=aurostd::RemoveSubString(arr, "[");
      arr=aurostd::RemoveSubString(arr, "]");
      arr=aurostd::RemoveSubString(arr, "\"");
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

  //writeStatsToJson//////////////////////////////////////////////////////////
  // Writes the database statistics into a JSON-formatted string(stream).
  void AflowDB::writeStatsToJson(std::stringstream& json, const DBStats& db_stats) {
    string tab = "    ";
    string indent = tab + tab;
    json << indent << db_stats.catalog << ": {" << std::endl;
    json << indent << tab << "\"count\": " << db_stats.nentries << "," << std::endl;
    json << indent << tab << "\"systems\": " <<  db_stats.nsystems << "," << std::endl;
    for (uint l = 0; l < db_stats.loop_counts.size(); l++) {
      json << indent << tab << "\"" << db_stats.loop_counts[l].first << "\": ";
      json << db_stats.loop_counts[l].second << "," << std::endl;
    }
    json << indent << tab << "\"columns\": {" << std::endl;
    uint ncols = db_stats.columns.size();
    string str_formatted = "";
    for (uint c = 0; c < ncols; c++) {
      json << indent << tab << tab << "\"" << db_stats.columns[c] << "\": {" << std::endl;
      json << indent << tab << tab << tab << "\"count\": ";
      if (db_stats.types[c] == "bool") {
        json << "{" << std::endl;
        json << indent << tab << tab << tab << tab << "\"true\": " << db_stats.count[c][0] << "," << std::endl;
        json << indent << tab << tab << tab << tab << "\"false\": " << db_stats.count[c][1] << std::endl;
        json << indent << tab << tab << tab << "}," << std::endl;
      } else {
        json << aurostd::utype2string<int>(db_stats.count[c][0]) << "," << std::endl;
      }
      str_formatted = db_stats.min[c];
      if ((str_formatted[0] == '\"') && (str_formatted[str_formatted.size()-1] == '\"')) {  // Don't escape enclosing strings //DX20200317 - back() doesn't work for old GCC versions
        str_formatted = str_formatted.substr(1, str_formatted.size() - 2);
        str_formatted = aurostd::StringSubst(str_formatted, "\"", "\\\"");
        str_formatted = "\"" + str_formatted + "\"";
      // ME20201024 - Breaks arrays of strings. Need to find a way to account for quotes
      // in strings for arrays of strings. This corner case has not occurred yet though.
      //} else {
      //  str_formatted = aurostd::StringSubst(str_formatted, "\"", "\\\"");
      }
      json << indent << tab << tab << tab << "\"min\": "
        << (db_stats.min[c].empty()?"null":str_formatted) << "," << std::endl;
      str_formatted = db_stats.max[c];
      if ((str_formatted[0] == '\"') && (str_formatted[str_formatted.size()-1] == '\"')) {  // Don't escape enclosing strings //DX20200317 - back() doesn't work for old GCC versions
        str_formatted = str_formatted.substr(1, str_formatted.size() - 2);
        str_formatted = aurostd::StringSubst(str_formatted, "\"", "\\\"");
        str_formatted = "\"" + str_formatted + "\"";
      // ME20201024 - Breaks arrays of strings. Need to find a way to account for quotes
      // in strings for arrays of strings. This corner case has not occurred yet though.
      //} else {
      //  str_formatted = aurostd::StringSubst(str_formatted, "\"", "\\\"");
      }
      json << indent << tab << tab << tab << "\"max\": "
        << (db_stats.max[c].empty()?"null":str_formatted) << "," << std::endl;

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
          if ((str_formatted[0] == '\"') && (str_formatted[str_formatted.size()-1] == '\"')) {  // Don't escape enclosing strings //DX20200317 - back() doesn't work for old GCC versions
            str_formatted = str_formatted.substr(1, str_formatted.size() - 2);
            str_formatted = aurostd::StringSubst(str_formatted, "\"", "\\\"");
            str_formatted = "\"" + str_formatted + "\"";
          // ME20201024 - Breaks arrays of strings. Need to find a way to account for quotes
          // in strings for arrays of strings. This corner case has not occurred yet though.
          //} else {
          //  str_formatted = aurostd::StringSubst(str_formatted, "\"", "\\\"");
          }
          json << indent << tab << tab << tab << tab << str_formatted;
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

  // INDEX -------------------------------------------------------------------

  //createIndex///////////////////////////////////////////////////////////////
  // Creates an index.
  void AflowDB::createIndex(const string& index, const string& table, const string& column) {
    string command = "CREATE INDEX " + index + " ON " + table + "(" + column +  ")";
    sql::SQLexecuteCommand(db, command);
  }

  //dropIndex/////////////////////////////////////////////////////////////////
  // Removes an index.
  void AflowDB::dropIndex(const string& index) {
    string command = "DROP INDEX " + index;
    sql::SQLexecuteCommand(db, command);
  }

  // TRANSACTION -------------------------------------------------------------

  //transaction///////////////////////////////////////////////////////////////
  // Begings (begin == true) or ends a database transaction.
  void AflowDB::transaction(bool begin) {
    string command = string(begin?"BEGIN":"END") + " TRANSACTION";
    sql::SQLexecuteCommand(db, command);
  }

  // TABLE -------------------------------------------------------------------

  //dropTable/////////////////////////////////////////////////////////////////
  // Deletes a table from the database.
  void AflowDB::dropTable(const string& table) {
    string command = "DROP TABLE IF EXISTS " + table;
    sql::SQLexecuteCommand(db, command);
  }

  //getTables/////////////////////////////////////////////////////////////////
  // Retrieves a set of tables. If where is empty, all tables in the database
  // will be returned.
  vector<string> AflowDB::getTables(const string& where) {
    return getTables(db, where);
  }

  vector<string> AflowDB::getTables(sqlite3* cursor, const string& where) {
    string command = "SELECT name FROM sqlite_master WHERE type='table'";
    if (!where.empty()) command += " AND (" + where + ")";
    return sql::SQLexecuteCommandVECTOR(cursor, command);
  }

  //createTable///////////////////////////////////////////////////////////////
  // Creates a table where all columns have the same type.
  void AflowDB::createTable(const string& table, const vector<string>& cols, const string& type) {
    vector<string> types(cols.size(), type);
    createTable(table, cols, types);
  }

  // Creates a table where each column is assigned its own type
  void AflowDB::createTable(const string& table, const vector<string>& cols, const vector<string>& types) {
    uint ncols = cols.size();
    if (ncols != types.size()) {
      string function = XPID + _AFLOW_DB_ERR_PREFIX_ + "createTable():";
      string message = "Could not create table. ";
      message += "Number of columns and number of types do not match.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    } else {
      string command = "CREATE TABLE IF NOT EXISTS " + table + " (";
      for (uint c = 0; c < ncols; c++) {
        command += cols[c] + " " + types[c];
        if (c < ncols - 1) command += ", ";
      }
      command += ")";
      sql::SQLexecuteCommand(db, command);
    }
  }

  // INSERT ------------------------------------------------------------------

  //insertValues//////////////////////////////////////////////////////////////
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
      string function = XPID +  _AFLOW_DB_ERR_PREFIX_ + "insertValues():";
      string message = "Could not insert values. ";
      message += "Number of columns and number of values do not match.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    } else {
      string command = "INSERT INTO " + table;
      if (ncols > 0) command += " (" + aurostd::joinWDelimiter(cols, ", ") + ")";
      command += " VALUES(" + aurostd::joinWDelimiter(vals, ", ") + ")";
      sql::SQLexecuteCommand(db, command);
    }
  }

  // UPDATE ------------------------------------------------------------------
  // Performs an UPDATE command.
  // WARNING: if WHERE is empty, all rows will updated! Never default to an
  // empty string!
  void AflowDB::updateRow(const string& table, const vector<string>& cols,
      const vector<string>& vals, const string& where) {
    uint ncols = cols.size();
    uint nvals = vals.size();
    if (ncols == 0) {
      string function = XPID +  _AFLOW_DB_ERR_PREFIX_ + "updateRow():";
      string message = "Could not update row. ";
      message += "No columns selected.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    } else if (ncols != nvals) {
      string function = XPID +  _AFLOW_DB_ERR_PREFIX_ + "updateRow():";
      string message = "Could not update row. ";
      message += "Number of columns and number of values do not match.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    } else {
      string command = "UPDATE " + table + " SET ";
      for (uint i = 0; i < ncols; i++) {
        command += cols[i] + "=" + vals[i];
        if (i < ncols - 1) command += ",";
      }
      if (!where.empty()) command += " WHERE " + where;
      sql::SQLexecuteCommand(db, command);
    }
  }


  // GET ---------------------------------------------------------------------

  //getColumnNames////////////////////////////////////////////////////////////
  // Returns the names of all columns in a specific table.
  vector<string> AflowDB::getColumnNames(const string& table) {
    return getColumnNames(db, table);
  }

  vector<string> AflowDB::getColumnNames(sqlite3* cursor, const string& table) {
    string command = "PRAGMA table_info(" + table + ")";
    vector<vector<string> > pragma_results = sql::SQLexecuteCommand2DVECTOR(cursor, command);
    uint nresults = pragma_results.size();
    vector<string> columns(nresults);
    for (uint r = 0; r < nresults; r++) {
      columns[r] = pragma_results[r][1];
    }
    return columns;
  }

  //getColumnTypes////////////////////////////////////////////////////////////
  // Returns the data types of all columns in a specific table.
  vector<string> AflowDB::getColumnTypes(const string& table) {
    return getColumnTypes(db, table);
  }

  vector<string> AflowDB::getColumnTypes(sqlite3* cursor, const string& table) {
    string command = "PRAGMA table_info(" + table + ")";
    vector<vector<string> > pragma_results = sql::SQLexecuteCommand2DVECTOR(cursor, command);
    uint nresults = pragma_results.size();
    vector<string> columns(nresults);
    for (uint r = 0; r < nresults; r++) {
      columns[r] = pragma_results[r][2];
    }
    return columns;
  }

  //getValue//////////////////////////////////////////////////////////////////
  // Gets a value from a specific column. The row must be specified in the
  // where condition, or else it just takes the first value.
  string AflowDB::getValue(const string& table, const string& col, const string& where) {
    return getValue(db, table, col, where);
  }

  string AflowDB::getValue(sqlite3* cursor, const string& table, const string& col, const string& where) {
    string command = prepareSELECT(table, "", col, where, 0, "");
    return sql::SQLexecuteCommandSCALAR(cursor, command);
  }

  //getProperty///////////////////////////////////////////////////////////////
  // Gets a database property for a specific column.
  string AflowDB::getProperty(const string& property, const string& table,
      const string& col, const string& where) {
    return getProperty(db, property, table, col, where);
  }
  string AflowDB::getProperty(sqlite3* cursor, const string& property, const string& table,
      const string& col, const string& where) {
    string command = prepareSELECT(table, property, col, where, 0, "");
    return sql::SQLexecuteCommandSCALAR(cursor, command);
  }

  //getPropertyMultiTables////////////////////////////////////////////////////
  // Gets a database property for a specific column across multiple tables.
  vector<string> AflowDB::getPropertyMultiTables(const string& property, const vector<string>& tables,
      const string& col, const string& where) {
    return getPropertyMultiTables(db, property, tables, col, where);
  }

  vector<string> AflowDB::getPropertyMultiTables(sqlite3* cursor, const string& property, const vector<string>& tables,
      const string& col, const string& where) {
    uint ntables = tables.size();
    vector<string> commands(ntables);
    for (uint t = 0; t < ntables; t++) commands[t] = prepareSELECT(tables[t], property, col, where);
    return sql::SQLexecuteCommandVECTOR(cursor, aurostd::joinWDelimiter(commands, " UNION ALL "));
  }

  //getSet////////////////////////////////////////////////////////////////////
  // Retrieves a (distinct) set from a single column.
  vector<string> AflowDB::getSet(const string& table, const string& col, bool distinct,
      const string& where, int limit, const string& order_by) {
    return getSet(db, table, col, distinct, where, limit, order_by);
  }

  vector<string> AflowDB::getSet(sqlite3* cursor, const string& table, const string& col, bool distinct,
      const string& where, int limit, const string& order_by) {
    string property = string((distinct?"DISTINCT":""));
    string command = prepareSELECT(table, property, col, where, limit, order_by);
    return sql::SQLexecuteCommandVECTOR(cursor, command);
  }

  //getSetMulitTables/////////////////////////////////////////////////////////
  // Retrieves a (distinct) set from a single column across multiple tables.
  // The result is sorted already, so there is not need for order_by.
  vector<string> AflowDB::getSetMultiTables(const vector<string>& tables, const string& col,
      bool distinct, const string& where, int limit) {
    return getSetMultiTables(db, tables, col, distinct, where, limit);
  }
  vector<string> AflowDB::getSetMultiTables(sqlite3* cursor, const vector<string>& tables, const string& col,
      bool distinct, const string& where, int limit) {
    string property = string((distinct?"DISTINCT":""));
    uint ntables = tables.size();
    vector<string> commands(ntables);
    for (uint t = 0; t < ntables; t++) commands[t] = prepareSELECT(tables[t], property, col, where, 0);
    string union_string = " UNION ";
    if (!distinct) union_string += "ALL ";
    string command = aurostd::joinWDelimiter(commands, union_string);
    if (limit > 0) command += " LIMIT " + aurostd::utype2string<int>(limit);
    return sql::SQLexecuteCommandVECTOR(cursor, command);
  }

  // SELECT ------------------------------------------------------------------

  //prepateSELECT/////////////////////////////////////////////////////////////
  // Lower level function to prepare a SELECT statement for all GET functions.
  string AflowDB::prepareSELECT(const string& table, const string& property, const string& cols,
      const string& where, int limit, const string& order_by) {
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
      const string& where, int limit, const string& order_by) {
    return prepareSELECT(table, property, aurostd::joinWDelimiter(cols, ", "), where, limit, order_by);
  }

}  // namespace aflowlib

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *            Aflow MARCO ESTERS - Duke University 2019-2020               *
// *                                                                         *
//****************************************************************************
