// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
// Written by Marco Esters, 2022
// Based on prior work by Hagen Eckert
//
// Unit test class to conduct parallelized unit tests with unified output.
//
// Tests are grouped into categories so that entire suites can be tested.
// For example, the test "aurostd" calls xscalar, xmatrix, etc. The individual
// tests can be called as well.

#include "aflow.h"
#include "aflow_anrl.h"  //DX20201104
#include "aflow_compare_structure.h"

using namespace std::placeholders;

namespace unittest {

  UnitTest::UnitTest(ostream& oss) : xStream(oss) {
    initialize();
  }

  UnitTest::UnitTest(ofstream& mf, ostream& oss) : xStream(mf, oss) {
    initialize();
  }

  UnitTest::UnitTest(const UnitTest& ut) : xStream(*ut.getOFStream(), *ut.getOSS()) {
    copy(ut);
  }

  const UnitTest& UnitTest::operator=(const UnitTest& ut) {
    copy(ut);
    return *this;
  }

  UnitTest::~UnitTest() {
    xStream::free();
    free();
  }

  void UnitTest::clear() {
    free();
  }

  void UnitTest::free() {
    aflags.clear();
    test_functions.clear();
    test_groups.clear();
  }

  void UnitTest::copy(const UnitTest& ut) {
    if (this == &ut) return;
    aflags = ut.aflags;
    test_functions = ut.test_functions;
    test_groups = ut.test_groups;
  }

  void UnitTest::initialize() {
    free();
    string dir = XHOST.vflag_control.getattachedscheme("DIRECTORY");
    aflags.Directory = (!dir.empty()?dir:aurostd::getPWD());

    initializeTestFunctions();
    initializeTestGroups();
  }

  /// @brief Initialize unit tests and add them to map of test functions.
  void UnitTest::initializeTestFunctions() {
    xcheck xchk;

    // aurostd
    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::xscalarTest, this, _1, _2, _3);
    xchk.function_name = "xscalarTest():";
    xchk.task_description = "xscalar functions";
    test_functions["xscalar"] = xchk;

    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::xvectorTest, this, _1, _2, _3);
    xchk.function_name = "xvectorTest():";
    xchk.task_description = "xvector functions";
    test_functions["xvector"] = xchk;

    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::xmatrixTest, this, _1, _2, _3);
    xchk.function_name = "xmatrixTest():";
    xchk.task_description = "xmatrix functions";
    test_functions["xmatrix"] = xchk;

    // database
    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::schemaTest, this, _1, _2, _3);
    xchk.function_name = "schemaTest():";
    xchk.task_description = "AFLOW schema";
    test_functions["schema"] = xchk;

    // xstructure
    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::atomicEnvironmentTest, this, _1, _2, _3);
    xchk.function_name = "atomicEnvironmentTest():";
    xchk.task_description = "Creating atomic environments";
    test_functions["atomic_environment"] = xchk;

    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::xstructureParserTest, this, _1, _2, _3);
    xchk.function_name = "xstructureParserTest():";
    xchk.task_description = "xstructure parsers";
    test_functions["xstructure_parser"] = xchk;

    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::xstructureTest, this, _1, _2, _3);
    xchk.function_name= "xstructureTest():";
    xchk.task_description = "xstructure functions";
    test_functions["xstructure"] = xchk;

    // structure generation
    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::ceramgenTest, this, _1, _2, _3);
    xchk.function_name= "ceramgenTest():";
    xchk.task_description = "pflow::GENERATE_CERAMICS()";
    test_functions["ceramgen"] = xchk;

    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::prototypeGeneratorTest, this, _1, _2, _3);
    xchk.function_name= "prototypeGeneratorTest():";
    xchk.task_description = "Generate all prototypes and test symmetry";
    test_functions["proto"] = xchk;

    // ovasp
    //Not working yet because we cannot load OUTCARs via the RestAPI
    //xchk = initializeXCheck();
    //xchk.func = std::bind(&UnitTest::xoutcarTest, this, _1, _2, _3);
    //xchk.function_name= "xoutcarTest():";
    //xchk.task_description = "xOUTCAR class";
    //test_functions["outcar"] = xchk;
  }

  /// @brief Initialize xcheck struct to an empty object.
  xcheck UnitTest::initializeXCheck() {
    xcheck xt;
    resetUnitTest(xt);
    xt.func = nullptr;
    xt.function_name = "";
    xt.task_description = "";
    xt.test_group = "";
    return xt;
  }

  /// @brief Reset a unit test to its pre-run state based on test name.
  void UnitTest::resetUnitTest(const string& test_name) {
    if (test_functions.count(test_name)) {
      resetUnitTest(test_functions[test_name]);
    }
  }

  /// @brief Reset an xcheck object to its pre-run state.
  void UnitTest::resetUnitTest(xcheck& test) {
    test.errors.clear();
    test.finished = false;
    test.passed_checks = 0;
    test.results.clear();
  }

  /// @brief Create unit test groups.
  void UnitTest::initializeTestGroups() {
    test_groups.clear();

    test_groups["aurostd"] = {"xscalar", "xvector", "xmatrix"};
    test_groups["database"] = {"schema"};
    test_groups["structure"] = {"atomic_environment", "xstructure", "xstructure_parser"};
    test_groups["structure_gen"] = {"ceramgen", "proto"};
    //test_groups["ovasp"] = {"outcar"};

    for (std::map<string, vector<string> >::iterator it = test_groups.begin(); it != test_groups.end(); ++it) {
      const vector<string>& members = (*it).second;
      for (size_t i = 0; i < members.size(); i++) {
        test_functions[members[i]].test_group = (*it).first;
      }
    }
  }

}

// Run functions
namespace unittest {

  /// @brief Run single unit test.
  ///
  /// @param unit_test Unit test name
  ///
  /// @return Whether all tests were successful.
  bool UnitTest::runTestSuites(const string& unit_test) {
    vector<string> vunit_test(1, unit_test);
    return runTestSuites(vunit_test);
  }

  /// @brief Run a set of unit tests.
  ///
  /// @param unit_tests_in Set of unit test names.
  ///
  /// @return Whether all tests were successful.
  ///
  /// The function consists of three steps:
  ///   1) Collect the set of tasks and expand test groups into individual tests.
  ///   2) Run all requested test functions, outputting results as test groups finish.
  ///   3) Print final summary.
  bool UnitTest::runTestSuites(const vector<string>& unit_tests_in) {
    stringstream  message;
    // Create task lists (groups or individual tests).
    // unit_test is the individual small tests over which to parallelize
    vector<string> unit_tests, tasks;
    for (size_t t = 0; t < unit_tests_in.size(); t++) {
      const string& test = unit_tests_in[t];
      bool isgroup = (test_groups.find(test) != test_groups.end());
      if (test == "all") {
        tasks.clear();
        for (std::map<string, vector<string> >::iterator it = test_groups.begin(); it != test_groups.end(); ++it) {
          tasks.push_back((*it).first);
          for (size_t m = 0; m < (*it).second.size(); m++) {
            unit_tests.push_back((*it).second[m]);
          }
        }
        break;
      } else if (!isgroup && (test_functions.find(test) == test_functions.end())) {
        message << "Skipping unrecognized test name " << test << ".";
        pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      } else if (isgroup && !aurostd::WithinList(tasks, test)) {
        tasks.push_back(test);
        const vector<string>& members = test_groups[test];
        for (size_t m = 0; m < members.size(); m++) {
          unit_tests.push_back(members[m]);
        }
      } else if (!aurostd::WithinList(tasks, test_functions[test].test_group)) {
        unit_tests.push_back(test);
        tasks.push_back(test);
      }
    }
    uint ntasks = tasks.size();
    if (ntasks == 0) {
      message << "No unit tests to run.";
      pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);
      return true;
    }

    // Run tests
    // Many AFLOW functions produce output to screen without the opportunity
    // to silence it, which makes unit test output harder to read and can
    // lead to garbled output when run in parallel. To get clean output,
    // silence output globally except for the functions that produce desired
    // unit test output, unless --quiet is requested or --debug is run.
    bool quiet_copy = XHOST.QUIET;
    vector<string> whitelist;
    if (!XHOST.QUIET && !XHOST.DEBUG) {
      whitelist.push_back("unittest::UnitTest::runUnitTest()");
      // Add function names to whitelist for displayResults
      for (uint i = 0; i < unit_tests.size(); i++) whitelist.push_back(test_functions[unit_tests[i]].function_name);
      XHOST.QUIET = true;
      for (size_t i = 0; i < whitelist.size(); i++) XHOST.LOGGER_WHITELIST.push_back(whitelist[i]);
    }
#ifdef AFLOW_MULTITHREADS_ENABLE
    std::mutex mtx;
    xthread::xThread xt(KBIN::get_NCPUS());
    std::function<void(vector<string>::iterator&, const vector<string>&)> fn = std::bind(&UnitTest::runUnitTest, this, _1, _2);
    xt.run(unit_tests, fn, tasks);
#else
    for (vector<string>::iterator it = unit_tests.begin(); it != unit_tests.end(); ++it) runUnitTest(it, tasks);
#endif
    XHOST.QUIET = quiet_copy;
    if (!XHOST.QUIET && !XHOST.DEBUG) {
      XHOST.QUIET = quiet_copy;
      for (size_t i = 0; i < whitelist.size(); i++) XHOST.LOGGER_WHITELIST.pop_back();
    }

    // Print final summary
    uint nsuccess = 0;
    vector<vector<string> > summary(tasks.size(), vector<string>(2));
    for (size_t t = 0; t < tasks.size(); t++) {
      bool success = taskSuccessful(tasks[t]);
      if (success) nsuccess++;
      summary[t][0] = tasks[t];
      summary[t][1] = (success?"pass":"fail");
    }

    if (nsuccess == ntasks) {
      message << "Unit tests passed successfully (passing " << ntasks << " task" + string((ntasks == 1)?"":"s") + ").";
      pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    } else {
      message << "Some unit tests failed (" << (ntasks - nsuccess) << " of " << ntasks << " failed).";
      pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
    }
    pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, formatResultsTable(summary), aflags, *p_FileMESSAGE, *p_oss, _LOGGER_RAW_);
    return (nsuccess == ntasks);
  }

  /// @brief Run unit test function inside a thread.
  ///
  /// @param it    Iterator pointing to the unit test name.
  /// @param tasks Set of task names.
  void UnitTest::runUnitTest(vector<string>::iterator& it, const vector<string>& tasks) {
    const string& test_name = (*it);
    xcheck& test = test_functions[test_name];
    resetUnitTest(test);
    test.func(test.passed_checks, test.results, test.errors);
    // Output results
    std::lock_guard<std::mutex> lk(mtx);
    test.finished = true;
    if (aurostd::WithinList(tasks, test_name)) {
      // If the test name is in the task list, it is not part
      // of a group, so no need to check if other members are done
      displayResult(test);
    } else {
      // Test if part of a test group, so check if all members
      // of the group have finished before producing output
      const string& group = test_functions[test_name].test_group;
      const vector<string>& vtests_group = test_groups[group];
      uint ntests_group = vtests_group.size();
      uint i = 0;
      for (; i < ntests_group; i++) {
        const string& test_name_group = vtests_group[i];
        if (!test_functions[test_name_group].finished) break;
      }
      if (i == ntests_group) {
        // All finished - print output
        uint nsuccess = 0;
        for (size_t t = 0; t < vtests_group.size(); t++) {
          if (taskSuccessful(vtests_group[t])) nsuccess++;
        }
        stringstream message;
        if (nsuccess == ntests_group) {
          message << "Unit tests of group " << group << " passed successfully (passing " << ntests_group << " test" << ((ntests_group == 1)?"":"s") << ").";
          pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
        } else {
          message << "Some unit tests of group " << group << " failed (" << (ntests_group - nsuccess) << " of " << ntests_group << " failed).";
          pflow::logger(_AFLOW_FILE_NAME_, __AFLOW_FUNC__, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        }
        for (size_t t = 0; t < vtests_group.size(); t++) {
          displayResult(test_functions[vtests_group[t]]);
        }
      }
    }
  }

  /// @brief Check if a task has finished successfully.
  ///
  /// @param task The task to check.
  ///
  /// @return Whether a tasks has finished without errors.
  ///
  /// Criteria for returning true:
  ///   1) Task is finished.
  ///   2) All unit tests passed.
  ///   3) There are no additional errors.
  bool UnitTest::taskSuccessful(const string& task) {
    std::map<string, vector<string> >::iterator it = test_groups.find(task);
    if (it != test_groups.end()) {
      const vector<string> members = (*it).second;
      for (size_t m = 0; m < members.size(); m++) {
        const string& member = members[m];
        const xcheck& xchk = test_functions[member];
        if (!xchk.finished || (xchk.errors.size() > 0) || (xchk.passed_checks != xchk.results.size())) return false;
      }
      return true;
    } else {
      const xcheck& xchk = test_functions[task];
      return (xchk.finished && (xchk.errors.size() == 0) && (xchk.passed_checks == xchk.results.size()));
    }
  }
}

// Output formatters
namespace unittest {
  /// @brief Convert results vector into a formatted table.
  ///
  /// @param table Structured table data.
  ///
  /// @return Formatted table string.
  ///
  /// This function makes no assumption about the table dimensions, so tables
  /// rows can have different sizes.
  /// Empty columns will be skipped, but not empty rows.
  /// The first column will have spaces prepended to indent the table.
  string UnitTest::formatResultsTable(const vector<vector<string> >& table) {
    size_t nrows = table.size();
    if (nrows == 0) return "";  // Empty table

    // Determine dimensions of the table
    vector<size_t> col_sizes;
    size_t str_length = 0;
    size_t maxcol = 0;
    for (size_t r = 0; r < nrows; r++) {
      maxcol = std::max(table[r].size(), maxcol);
      for (size_t c = 0; c < table[r].size(); c++) {
        str_length = table[r][c].length();
        if (c == col_sizes.size()) {
          col_sizes.push_back(str_length);
        } else if (str_length > col_sizes[c]) {
          col_sizes[c] = str_length;
        }
      }
    }
    if (maxcol == 0) return "";  // Empty rows

    vector<string> output(nrows), row;
    string col = "";
    for (size_t r = 0; r < nrows; r++) {
      row.clear();
      // Last column does not need padding
      for (size_t c = 0; c < maxcol - 1; c++) {
        if (col_sizes[c] > 0) {
          col = (c == 0?"  ":"") + aurostd::PaddedPOST((c < table[r].size()?table[r][c]:""), col_sizes[c]);
          row.push_back(col);
        }
      }
      if (table[r].size() == maxcol) row.push_back(table[r].back());
      output[r] = aurostd::joinWDelimiter(row, " | ");
    }
    return aurostd::joinWDelimiter(output, "\n");
  }

  /// @brief Display results of a unit test.
  ///
  /// @param xchk Unit test object containing all results.
  void UnitTest::displayResult(const xcheck& xchk) {
    stringstream message;
    size_t check_num = xchk.results.size();
    if (xchk.passed_checks == check_num) {
      if (xchk.errors.size() > 0) {
        // All attempted checks passed, but there were errors.
        // This happens when a prerequisite for a test fails (e.g. file loading).
        message << "FAIL " << xchk.task_description << " due to runtime errors" << std::endl;
        pflow::logger(_AFLOW_FILE_NAME_,xchk.function_name,message,aflags,*p_FileMESSAGE,*p_oss,_LOGGER_ERROR_);
      } else {
        message << "SUCCESS " << xchk.task_description << " (passing " << check_num << " check" << ((check_num == 1)?"":"s") << ")" << std::endl;
        pflow::logger(_AFLOW_FILE_NAME_,xchk.function_name,message,aflags,*p_FileMESSAGE,*p_oss,_LOGGER_COMPLETE_);
      }
    } else {
      message << "FAIL " << xchk.task_description << " (" << (check_num - xchk.passed_checks) << " of " << check_num << " checks failed)" << std::endl;
      pflow::logger(_AFLOW_FILE_NAME_,xchk.function_name,message,aflags,*p_FileMESSAGE,*p_oss,_LOGGER_ERROR_);
    }
    pflow::logger(_AFLOW_FILE_NAME_,xchk.function_name,formatResultsTable(xchk.results),aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    if (xchk.errors.size() > 0) {
      message << "\nAdditional error messages:\n" << aurostd::joinWDelimiter(xchk.errors, "\n");
      pflow::logger(_AFLOW_FILE_NAME_,xchk.function_name,message,aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    }
  }
}

// Collection of generic check functions, to streamline testing.
namespace unittest {
  template <typename utype>
  void UnitTest::checkEqual(const vector<utype>& calculated, const vector<utype>& expected, const string& check_function,
      const string& check_description, uint& passed_checks, vector<vector<string> >& results) {
    bool passed = (calculated.size() == expected.size());
    for (size_t i = 0; i < calculated.size() && passed; i++) {
      passed = aurostd::isequal(calculated[i], expected[i]);
    }
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }
  void UnitTest::checkEqual(const vector<string>& calculated, const vector<string>& expected, const string& check_function,
      const string& check_description, uint& passed_checks, vector<vector<string> >& results) {
    bool passed = (calculated.size() == expected.size());
    for (size_t i = 0; i < calculated.size() && passed; i++) {
      passed = (calculated[i] == expected[i]);
    }
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }

  template <typename utype>
  void UnitTest::checkEqual(const utype& calculated, const utype& expected, const string& check_function,
      const string& check_description, uint& passed_checks, vector<vector<string> >& results) {
    bool passed = (aurostd::isequal(calculated, expected));
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }
  void UnitTest::checkEqual(const string& calculated, const string& expected, const string& check_function,
      const string& check_description, uint & passed_checks, vector<vector<string> >& results) {
    bool passed = (calculated == expected);
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }
  void UnitTest::checkEqual(const bool calculated, const bool expected, const string& check_function,
      const string& check_description, uint& passed_checks, vector<vector<string> >& results) {
    bool passed = (calculated == expected);
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }

  template <typename utype>
  void UnitTest::check(const bool passed, const vector<utype>& calculated, const vector<utype>& expected, const string& check_function,
    const string& check_description, uint& passed_checks, vector<vector<string> >& results) {
    check(passed, aurostd::joinWDelimiter(calculated, ","), aurostd::joinWDelimiter(expected, ","), check_function, check_description, passed_checks, results);
  }
  void UnitTest::check(const bool passed, const vector<double>& calculated, const vector<double>& expected, const string& check_function,
    const string& check_description, uint& passed_checks, vector<vector<string> >& results) {
    check(passed, aurostd::joinWDelimiter(aurostd::vecDouble2vecString(calculated), ","), aurostd::joinWDelimiter(aurostd::vecDouble2vecString(expected), ","), check_function, check_description, passed_checks, results);
  }

  template <typename utype>
  void UnitTest::check(const bool passed, const xmatrix<utype>& calculated, const xmatrix<utype>& expected, const string& check_function,
    const string& check_description, uint& passed_checks, vector<vector<string> >& results) {
    check(passed, aurostd::xmat2String(calculated), aurostd::xmat2String(expected), check_function, check_description, passed_checks, results);
  }
  void UnitTest::check(const bool passed, const xmatrix<double>& calculated, const xmatrix<double>& expected, const string& check_function,
    const string& check_description, uint& passed_checks, vector<vector<string> >& results) {
    check(passed, aurostd::xmatDouble2String(calculated), aurostd::xmatDouble2String(expected), check_function, check_description, passed_checks, results);
  }

  /// @brief Base function to check results and update results.
  ///
  /// @param passed            Whether the test has passed.
  /// @param calculated        Calculated value.
  /// @param expected          Expected value.
  /// @param check_function    Function called for the test.
  /// @param check_description Description of the performed test.
  /// @param passed_checks     Number of passed checks.
  /// @param results           Results data - doubles as number of performed checks.
  template <typename utype>
  void UnitTest::check(const bool passed, const utype& calculated, const utype& expected, const string& check_function,
      const string& check_description, uint& passed_checks, vector<vector<string> >& results) {
    vector<string> result;
    uint check_num = results.size() + 1;
    result.push_back(aurostd::utype2string<uint>(check_num));
    if (passed) {
      passed_checks++;
      result.push_back("pass");
    } else {
      result.push_back("FAIL");
    }
    result.push_back(check_function);
    result.push_back(check_description);
    if (!passed) {
      stringstream failstring;
      failstring << " (result: " << calculated << " | expected: " << expected << ")";
      result.back() += failstring.str();
    }
    results.push_back(result);
  }

}

// aurostd
namespace unittest {

  void UnitTest::xscalarTest(uint& passed_checks, vector<vector<string> >& results, vector<string>& errors) {
    if (errors.size()) {}  // Suppress compiler warnings
    // setup test environment
    string check_function = "", check_description = "";
    double calculated_dbl = 0.0, expected_dbl = 0.0;
    int calculated_int = 0, expected_int = 0;
    vector<int> calculated_vint, expected_vint;

    // ---------------------------------------------------------------------------
    // Check | double2fraction conversion //DX20210908
    // ---------------------------------------------------------------------------
    check_function = "aurostd::double2fraction()";
    check_description = "convert a double to a fraction.";

    double test_double = 1.625;
    int numerator = 1, denominator = 1;
    string answer = "13/8";
    aurostd::double2fraction(test_double,numerator,denominator);
    stringstream result_ss; result_ss << numerator << "/" << denominator;

    checkEqual(result_ss.str(), answer, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | mod_floored (int) //SD20220124
    // ---------------------------------------------------------------------------
    check_function = "aurostd::mod_floored()";
    check_description = "floored mod; numbers as int";
    expected_int = -1;

    calculated_int = aurostd::mod_floored(5, -3);
    checkEqual(calculated_int, expected_int, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | mod_floored (double) //SD20220124
    // ---------------------------------------------------------------------------
    check_function = "aurostd::mod_floored()";
    check_description = "floored mod; numbers as double";
    expected_dbl = 1.4;

    calculated_dbl = aurostd::mod_floored(-5.2, 3.3);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | mod_floored (divisor 0) //SD20220124
    // ---------------------------------------------------------------------------
    check_function = "aurostd::mod_floored()";
    check_description = "floored mod; divisor is 0";
    expected_dbl = 11.11;

    calculated_dbl = aurostd::mod_floored(11.11, 0.0);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | mod_floored (divisor inf) //SD20220124
    // ---------------------------------------------------------------------------
    check_function = "aurostd::mod_floored()";
    check_description = "floored mod; divisor is inf";
    expected_dbl = 11.11;

    calculated_dbl = aurostd::mod_floored(11.11, (double)INFINITY);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | gcd //CO20190520
    // ---------------------------------------------------------------------------
    check_function = "aurostd::GCD()";
    int a=0,b=0,x1=0,y1=0,gcd=0;

    check_description = "gcd(25,15)";
    a=25;b=15;
    expected_vint = {5, -1, 2};
    aurostd::GCD(a,b,gcd,x1,y1);
    calculated_vint = {gcd, x1, y1};
    checkEqual(calculated_vint, expected_vint, check_function, check_description, passed_checks, results);

    check_description = "gcd(25,0)";
    a=25;b=0;
    expected_vint = {25, 1, 0};
    aurostd::GCD(a,b,gcd,x1,y1);
    calculated_vint = {gcd, x1, y1};
    checkEqual(calculated_vint, expected_vint, check_function, check_description, passed_checks, results);

    check_description = "gcd(0,15)";
    a=0;b=15;
    expected_vint = {15, 0, 1};
    aurostd::GCD(a,b,gcd,x1,y1);
    calculated_vint = {gcd, x1, y1};
    checkEqual(calculated_vint, expected_vint, check_function, check_description, passed_checks, results);

    check_description = "gcd(-5100,30450)";
    a=-5100;b=30450;
    expected_vint = {150, -6, -1};
    aurostd::GCD(a,b,gcd,x1,y1);
    calculated_vint = {gcd, x1, y1};
    checkEqual(calculated_vint, expected_vint, check_function, check_description, passed_checks, results);
  }

  void UnitTest::xvectorTest(uint& passed_checks, vector<vector<string> >& results, vector<string>& errors) {
    if (errors.size()) {}  // Suppress compiler warnings
    // setup test environment
    string check_function = "", check_description = "";

    //HE20210511
    double expected_dbl = 0.0, calculated_dbl = 0.0;
    int expected_int = 0;
    string expected_str = "";
    vector<xvector<double> > points;
    vector<xvector<int> > ipoints;
    vector<vector<uint> > facets;
    vector<uint> facet;

    // variables to store examples as doubles (p#) and int (p#i) variants
    xvector<double> p0(3,1); xvector<int> p0i(3,1);
    xvector<double> p1(3,1); xvector<int> p1i(3,1);
    xvector<double> p2(3,1); xvector<int> p2i(3,1);
    xvector<double> p3(3,1); xvector<int> p3i(3,1);
    xvector<double> p4(3,1); xvector<int> p4i(3,1);
    xvector<double> p5(3,1); xvector<int> p5i(3,1);
    xvector<double> p6(3,1); xvector<int> p6i(3,1);
    xvector<double> p7(3,1); xvector<int> p7i(3,1);
    xvector<double> p8(3,1); xvector<int> p8i(3,1);
    xvector<double> p9(3,1); xvector<int> p9i(3,1);
    xvector<double> p10(3,1); xvector<int> p10i(3,1);
    xvector<double> p11(3,1); xvector<int> p11i(3,1);

    // define convex solid
    p0i(1) = p0(1) = 0.0; p0i(2) = p0(2) = 0.0; p0i(3) = p0(3) = 0.0;
    p1i(1) = p1(1) = 1.0; p1i(2) = p1(2) = 0.0; p1i(3) = p1(3) = 0.0;
    p2i(1) = p2(1) = 1.0; p2i(2) = p2(2) = 1.0; p2i(3) = p2(3) = 0.0;
    p3i(1) = p3(1) = 0.0; p3i(2) = p3(2) = 1.0; p3i(3) = p3(3) = 0.0;
    p4i(1) = p4(1) = 0.0; p4i(2) = p4(2) = 0.0; p4i(3) = p4(3) = 2.0;
    p5i(1) = p5(1) = 1.0; p5i(2) = p5(2) = 0.0; p5i(3) = p5(3) = 2.0;
    p6i(1) = p6(1) = 1.0; p6i(2) = p6(2) = 1.0; p6i(3) = p6(3) = 3.0;
    p7i(1) = p7(1) = 0.0; p7i(2) = p7(2) = 1.0; p7i(3) = p7(3) = 3.0;

    // transfer data into vectors
    points.clear(); ipoints.clear(); facets.clear();
    points.push_back(p0); points.push_back(p1); points.push_back(p2); points.push_back(p3); points.push_back(p4);
    points.push_back(p5); points.push_back(p6); points.push_back(p7);
    ipoints.push_back(p0i); ipoints.push_back(p1i); ipoints.push_back(p2i); ipoints.push_back(p3i); ipoints.push_back(p4i);
    ipoints.push_back(p5i); ipoints.push_back(p6i); ipoints.push_back(p7i);
    facet.clear(); facet.resize(4); facet[0]=0; facet[1]=1; facet[2]=2; facet[3]=3; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=4; facet[1]=5; facet[2]=6; facet[3]=7; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=1; facet[1]=2; facet[2]=6; facet[3]=5; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=0; facet[1]=3; facet[2]=7; facet[3]=7; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=0; facet[1]=1; facet[2]=5; facet[3]=4; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=3; facet[1]=2; facet[2]=6; facet[3]=7; facets.push_back(facet);

    // ---------------------------------------------------------------------------
    // Check | convex solid volume (double)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "convex solid, points as doubles";
    expected_dbl = 2.5;

    calculated_dbl = aurostd::volume(points, facets, true);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | convex solid volume (int)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "convex solid, points as int";

    calculated_dbl = aurostd::volume(ipoints, facets, true);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);


    // define non convex solid
    p0i(1) = p0(1)   = 0.0; p0i(2) = p0(2)   = 0.0; p0i(3) = p0(3)   = 0.0;
    p1i(1) = p1(1)   = 0.0; p1i(2) = p1(2)   = 4.0; p1i(3) = p1(3)   = 0.0;
    p2i(1) = p2(1)   = 2.0; p2i(2) = p2(2)   = 4.0; p2i(3) = p2(3)   = 0.0;
    p3i(1) = p3(1)   = 1.0; p3i(2) = p3(2)   = 1.0; p3i(3) = p3(3)   = 0.0;
    p4i(1) = p4(1)   = 4.0; p4i(2) = p4(2)   = 2.0; p4i(3) = p4(3)   = 0.0;
    p5i(1) = p5(1)   = 4.0; p5i(2) = p5(2)   = 0.0; p5i(3) = p5(3)   = 0.0;
    p6i(1) = p6(1)   = 0.0; p6i(2) = p6(2)   = 0.0; p6i(3) = p6(3)   = 4.0;
    p7i(1) = p7(1)   = 0.0; p7i(2) = p7(2)   = 4.0; p7i(3) = p7(3)   = 4.0;
    p8i(1) = p8(1)   = 2.0; p8i(2) = p8(2)   = 4.0; p8i(3) = p8(3)   = 4.0;
    p9i(1) = p9(1)   = 1.0; p9i(2) = p9(2)   = 1.0; p9i(3) = p9(3)   = 4.0;
    p10i(1) = p10(1) = 4.0; p10i(2) = p10(2) = 2.0; p10i(3) = p10(3) = 4.0;
    p11i(1) = p11(1) = 4.0; p11i(2) = p11(2) = 0.0; p11i(3) = p11(3) = 4.0;

    // transfer data into vectors
    points.clear(); ipoints.clear(); facets.clear();
    points.push_back(p0); points.push_back(p1); points.push_back(p2); points.push_back(p3); points.push_back(p4);
    points.push_back(p5); points.push_back(p6); points.push_back(p7); points.push_back(p8); points.push_back(p9);
    points.push_back(p10); points.push_back(p11);
    ipoints.push_back(p0i); ipoints.push_back(p1i); ipoints.push_back(p2i); ipoints.push_back(p3i); ipoints.push_back(p4i);
    ipoints.push_back(p5i); ipoints.push_back(p6i); ipoints.push_back(p7i); ipoints.push_back(p8i); ipoints.push_back(p9i);
    ipoints.push_back(p10i); ipoints.push_back(p11i);

    facet.clear(); facet.resize(6); facet[0]=5; facet[1]=4; facet[2]=3; facet[3]=2; facet[4]=1; facet[5]=0; facets.push_back(facet);
    facet.clear(); facet.resize(6); facet[0]=6; facet[1]=7; facet[2]=8; facet[3]=9; facet[4]=10; facet[5]=11; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=0; facet[1]=6; facet[2]=11; facet[3]=5; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=4; facet[1]=5; facet[2]=11; facet[3]=10; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=3; facet[1]=4; facet[2]=10; facet[3]=9; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=3; facet[1]=9; facet[2]=8; facet[3]=2; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=1; facet[1]=2; facet[2]=8; facet[3]=7; facets.push_back(facet);
    facet.clear(); facet.resize(4); facet[0]=0; facet[1]=1; facet[2]=7; facet[3]=6; facets.push_back(facet);

    // ---------------------------------------------------------------------------
    // Check | non convex solid volume (double)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "non convex solid, points as doubles";
    expected_dbl = 40.0;

    calculated_dbl = aurostd::volume(points, facets);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | error facet/normals mismatch
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "error: facet/normals mismatch";
    vector<xvector<double> > normals;
    expected_str = "xerror code 30 (VALUE_ERROR)";
    expected_int = _VALUE_ERROR_;

    try {
      calculated_dbl = aurostd::volume(points, facets, normals);
      check(false, std::string("no error"), expected_str, check_function, check_description, passed_checks, results);
    }
    catch (aurostd::xerror e)
    {
      if (e.whatCode() == expected_int) check(true, "", "", check_function, check_description, passed_checks, results);
      else check(false, aurostd::utype2string(e.whatCode()), expected_str, check_function, check_description, passed_checks, results);
    }
    catch (...) {
      check(false, std::string("not an xerror"), expected_str, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | non convex solid volume (int)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "non convex solid, points as int";
    expected_dbl = 40.0;

    calculated_dbl = aurostd::volume(ipoints, facets);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | error facet size
    // ---------------------------------------------------------------------------
    check_function = "aurostd::volume()";
    check_description = "error: wrong facet size";
    expected_str = "xerror code 30 (VALUE_ERROR)";
    expected_int = _VALUE_ERROR_;

    facet.clear(); facet.resize(2); facet[0]=1; facet[1]=2; facets.push_back(facet);
    try {
      calculated_dbl = aurostd::volume(points, facets);
      check(false, std::string("no error"), expected_str, check_function, check_description, passed_checks, results);
    }
    catch (aurostd::xerror e)
    {
      if (e.whatCode() == expected_int) check(true, "", "", check_function, check_description, passed_checks, results);
      else check(false, aurostd::utype2string(e.whatCode()), expected_str, check_function, check_description, passed_checks, results);
    }
    catch (...) {
      check(false, std::string("not an xerror"), expected_str, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | non convex area (double)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::areaPointsOnPlane()";
    check_description = "non convex area; points as double";
    expected_dbl = 10.0;

    //fill vectors with data
    points.clear(); ipoints.clear(); facets.clear();
    points.push_back(p0); points.push_back(p1); points.push_back(p2); points.push_back(p3); points.push_back(p4);
    points.push_back(p5);
    ipoints.push_back(p0i); ipoints.push_back(p1i); ipoints.push_back(p2i); ipoints.push_back(p3i); ipoints.push_back(p4i);
    ipoints.push_back(p5i);

    calculated_dbl = aurostd::areaPointsOnPlane(points);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | non convex area (int)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::areaPointsOnPlane()";
    check_description = "non convex area; points as int";
    expected_dbl = 10.0;

    calculated_dbl = aurostd::areaPointsOnPlane(ipoints);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);


    // define triangle in 3D to better test int handling
    p0i(1) = p0(1) = 0.0; p0i(2) = p0(2) = 0.0; p0i(3) = p0(3) = 0.0;
    p1i(1) = p1(1) = 1.0; p1i(2) = p1(2) = 1.0; p1i(3) = p1(3) = 1.0;
    p2i(1) = p2(1) = 5.0; p2i(2) = p2(2) = 0.0; p2i(3) = p2(3) = 5.0;

    points.clear(); ipoints.clear(); facets.clear();
    points.push_back(p0); points.push_back(p1); points.push_back(p2);
    ipoints.push_back(p0i); ipoints.push_back(p1i); ipoints.push_back(p2i);

    // ---------------------------------------------------------------------------
    // Check | 3d triangle area (double)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::areaPointsOnPlane()";
    check_description = "3d triangle; points as double";
    expected_dbl = 3.5355339059;

    calculated_dbl = aurostd::areaPointsOnPlane(points);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | 3d triangle area (int)
    // ---------------------------------------------------------------------------
    check_function = "aurostd::areaPointsOnPlane()";
    check_description = "3d triangle; points as int";
    expected_dbl = 3.5355339059;

    calculated_dbl = aurostd::areaPointsOnPlane(ipoints);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);
  }

  void UnitTest::xmatrixTest(uint& passed_checks, vector<vector<string> >& results, vector<string>& errors) {
    if (errors.size()) {}  // Suppress compiler warnings
    // setup test environment
    string check_function = "", check_description = "";

    xmatrix<int> calculated_xmatint, expected_xmatint;

    // ---------------------------------------------------------------------------
    // Check | reshape //SD20220319
    // ---------------------------------------------------------------------------
    check_function = "aurostd::reshape()";
    check_description = "reshape a rectangular matrix";
    expected_xmatint = xmatrix<int>(3,4);
    expected_xmatint(1,1) = 1; expected_xmatint(1,2) =  2; expected_xmatint(1,3) =  3; expected_xmatint(1,4) = 4;
    expected_xmatint(2,1) = 5; expected_xmatint(2,2) =  6; expected_xmatint(2,3) =  7; expected_xmatint(2,4) = 8;
    expected_xmatint(3,1) = 9; expected_xmatint(3,2) = 10; expected_xmatint(3,3) = 11; expected_xmatint(3,4) = 12;
    calculated_xmatint = aurostd::reshape(aurostd::reshape(expected_xmatint,4,3),3,4);
    checkEqual(calculated_xmatint, expected_xmatint, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | ehermite //CO20190520
    // ---------------------------------------------------------------------------
    check_function = "aurostd::getEHermite()";
    check_description = "calculate elementary Hermite transformation";
    expected_xmatint = xmatrix<int>(2, 2);
    expected_xmatint[1][1] =   5; expected_xmatint[1][2] = -2;
    expected_xmatint[2][1] = -12; expected_xmatint[2][2] =  5;
    calculated_xmatint = xmatrix<int>(2, 2);
    aurostd::getEHermite(5, 12, calculated_xmatint);
    checkEqual(calculated_xmatint, expected_xmatint, check_function, check_description, passed_checks, results);
  }

}

// database
namespace unittest {

  void UnitTest::schemaTest(uint& passed_checks, vector<vector<string> >& results, vector<string>& errors) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if (errors.size()) {}  // Suppress compiler warnings

    // Set up test environment
    string check_function = "";
    string check_description = "";

    // ---------------------------------------------------------------------------
    // Check | internal consistency
    // ---------------------------------------------------------------------------
    check_function = "XHOST.vschema";
    check_description = "Internal consistency of vschema";
    vector<string> vschema_keys;
    vector<string> vschema_types = {"UNIT", "TYPE"};
    string key = "";
    uint ninconsistent = 0;
    for (size_t i = 0; i < XHOST.vschema.vxsghost.size(); i+= 2) {
      if(XHOST.vschema.vxsghost[i].find("::NAME:") != string::npos) {
        key=aurostd::RemoveSubString(XHOST.vschema.vxsghost[i], "SCHEMA::NAME:");
        vschema_keys.push_back(XHOST.vschema.getattachedscheme("SCHEMA::NAME:" + key));
        for (size_t j = 0; j < vschema_types.size(); j++) {
          if (!XHOST.vschema.isdefined("SCHEMA::" + vschema_types[j] + ":" + key)) {
            ninconsistent++;
            if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " SCHEMA::" << vschema_types[j] << ":" << key << " not found." << std::endl;
          }
        }
      }
    }
    checkEqual(ninconsistent, 0, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | consistency between _aflowlib_entry and schema
    // ---------------------------------------------------------------------------
    check_function = "_aflowlib_entry";
    check_description = "Consistency between _aflowlib_entry json and schema";
    aflowlib::_aflowlib_entry aentry;
    string aflowlib_json = aentry.aflowlib2string("JSON", true);
    vector<string> json_keys = aurostd::extractJsonKeysAflow(aflowlib_json);

    vector<string> vkeys_ignore = {"data_language", "error_status", "natoms_orig",
                                   "density_orig", "volume_cell_orig", "volume_atom_orig",
                                   "spinD_magmom_orig"};
    for (size_t k = 0; k < json_keys.size(); k++) {
      const string& key = json_keys[k];
      if (!aurostd::WithinList(vkeys_ignore, key) && !aurostd::WithinList(vschema_keys, key)) {
        ninconsistent++;
        if (LDEBUG) std::cerr << __AFLOW_FUNC__ << " " << key << " not found in schema." << std::endl;
      }
    }
    checkEqual(ninconsistent, 0, check_function, check_description, passed_checks, results);
  }

}

// structure tests
namespace unittest {

  //HE20210511
  void UnitTest::atomicEnvironmentTest(uint& passed_checks, vector<vector<string> >& results, vector<string>& errors) {
    if (errors.size()) {}  // Suppress compiler warnings

    // setup test environment
    string check_function = "";
    string check_description = "";

    // ---------------------------------------------------------------------------
    // Test 1: create AE - mode 1
    // ---------------------------------------------------------------------------

    // load test system
    string auid = "aflow:d912e209c81aeb94";
    string aurl = "aflowlib.duke.edu:AFLOWDATA/LIB2_RAW/Ca_svCu_pv/138";
    aflowlib::_aflowlib_entry entry;
    xstructure str;

    entry.aurl=aurl;
    entry.auid=auid;
    pflow::loadXstructures(entry);
    str = entry.vstr.back();
    vector<AtomEnvironment> AE = getAtomEnvironments(str, 1);
    check_function = "getAtomEnvironments()";

    // ---------------------------------------------------------------------------
    // Check | number of created AEs
    check_description = "number of created AEs";
    checkEqual(uint(AE.size()), uint(6), check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | point index mapping
    check_description = "point index mapping";
    checkEqual(AE[1].index2Point(10), AE[1].coordinates_neighbor[1][2], check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | coordinate matching
    check_description = "coordinate matching";
    xvector<double> compare_point(3,1);
    compare_point(1) = -1.9551925593108925e0;
    compare_point(2) = -2.2642136090979212e0;
    compare_point(3) = 2.4896636484942385e0;

    checkEqual(AE[1].index2Point(2), compare_point, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | center element
    check_description = "center element";
    checkEqual(std::string(AE[2].element_center), std::string("Ca"), check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | center element id
    check_description = "center element id";
    checkEqual(AE[4].type_center, uint(1), check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Test 2: create AE convex hull
    // ---------------------------------------------------------------------------

    // setup test environment
    const uint test_AE = 4;

    // create hull
    AE[test_AE].constructAtomEnvironmentHull();
    check_function = "contructAtomEnvironmentHull()";

    // ---------------------------------------------------------------------------
    // Check | hull bit set
    check_description = "hull bit set";
    checkEqual(AE[test_AE].has_hull, true, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | hull volume
    check_description = "hull volume";
    checkEqual(AE[test_AE].volume, 31.4622167689, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | hull area
    check_description = "hull area";
    checkEqual(AE[test_AE].area, 60.4979100628, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | triangle count
    check_description = "triangle count";
    checkEqual(AE[test_AE].facet_order[0], uint(6), check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | tetragon count
    check_description = "tetragon count";
    checkEqual(AE[test_AE].facet_order[1], uint(2), check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | pentagon count
    check_description = "pentagon count";
    checkEqual(AE[test_AE].facet_order[2], uint(0), check_function, check_description, passed_checks, results);
  }

  void UnitTest::xstructureParserTest(uint& passed_checks, vector<vector<string> >& results, vector<string>& errors) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if (errors.size()) {}  // Suppress compiler warnings

    if (LDEBUG) std::cerr << "Running " << __AFLOW_FUNC__ << std::endl;

    // Set up test environment
    string check_function = "", check_description = "";
    bool calculated_bool = false, expected_bool = false;

    // Set up structure variables
    _aflags aflags; aflags.Directory = aurostd::getPWD();
    bool check_passed = true;
    stringstream xstrss;
    xstructure xstr_cif, xstr_poscar;
    string str_cif = "", str_poscar = "";

    bool same_species = true;
    bool scale_volume = false;
    bool optimize_match = false;
    double misfit = 0.0;

    // ---------------------------------------------------------------------------
    // Check | CIF parser
    // ---------------------------------------------------------------------------

    // ---------------------------------------------------------------------------
    // Test 1: parse structure with known settings
    // ---------------------------------------------------------------------------

    // CrO3 was a problematic structure in the past
    str_cif =
      "data_Cr4O12\n"
      "_pd_phase_name Cr4O12\n"
      "_cell_length_a 5.743\n"
      "_cell_length_b 8.557\n"
      "_cell_length_c 4.789\n"
      "_cell_angle_alpha 90.\n"
      "_cell_angle_beta 90.\n"
      "_cell_angle_gamma 90.\n"
      "_cell_volume 235.35\n"
      "_cell_formula_units_Z 4\n"
      "_space_group_name_H-M_alt 'A m a 2'\n"
      "_space_group_IT_number 40\n"
      "loop_\n"
      "_space_group_symop_id\n"
      "_space_group_symop_operation_xyz\n"
      "1 'x+1/2,-y,z'\n"
      "2 '-x+1/2,y,z'\n"
      "3 '-x,-y,z'\n"
      "4 'x,y,z'\n"
      "5 'x+1/2,-y+1/2,z+1/2'\n"
      "6 '-x+1/2,y+1/2,z+1/2'\n"
      "7 '-x,-y+1/2,z+1/2'\n"
      "8 'x,y+1/2,z+1/2'\n"
      "loop_\n"
      "_atom_type_symbol\n"
      "_atom_type_oxidation_number\n"
      "Cr6+ 6\n"
      "O2- -2\n"
      "loop_\n"
      "_atom_site_label\n"
      "_atom_site_type_symbol\n"
      "_atom_site_symmetry_multiplicity\n"
      "_atom_site_Wyckoff_symbol\n"
      "_atom_site_fract_x\n"
      "_atom_site_fract_y\n"
      "_atom_site_fract_z\n"
      "_atom_site_B_iso_or_equiv\n"
      "_atom_site_occupancy\n"
      "Cr1 Cr6+ 4 b 0.25 0.09676 0.5 . 1.\n"
      "O1 O2- 4 a 0. 0. 0.3841 . 1.\n"
      "O2 O2- 4 b 0.25 0.2677 0.3755 . 1.\n"
      "O3 O2- 4 b 0.25 0.6078 0.3284 . 1.\n";

    str_poscar =
      "Cr4O12\n"
      "1.000000\n"
      "   5.74300000000000   0.00000000000000   0.00000000000000\n"
      "   0.00000000000000   8.55700000000000   0.00000000000000\n"
      "   0.00000000000000   0.00000000000000   4.78900000000000\n"
      "4 12\n"
      "Direct(16) [A4B12]\n"
      "   0.25000000000000   0.09676000000000   0.50000000000000  Cr\n"
      "   0.75000000000000   0.90324000000000   0.50000000000000  Cr\n"
      "   0.25000000000000   0.59676000000000   0.00000000000000  Cr\n"
      "   0.75000000000000   0.40324000000000   0.00000000000000  Cr\n"
      "   0.00000000000000   0.00000000000000   0.38410000000000  O \n"
      "   0.50000000000000   0.00000000000000   0.38410000000000  O \n"
      "   0.00000000000000   0.50000000000000   0.88410000000000  O \n"
      "   0.50000000000000   0.50000000000000   0.88410000000000  O \n"
      "   0.25000000000000   0.26770000000000   0.37550000000000  O \n"
      "   0.75000000000000   0.73230000000000   0.37550000000000  O \n"
      "   0.25000000000000   0.76770000000000   0.87550000000000  O \n"
      "   0.75000000000000   0.23230000000000   0.87550000000000  O \n"
      "   0.25000000000000   0.60780000000000   0.32840000000000  O \n"
      "   0.75000000000000   0.39220000000000   0.32840000000000  O \n"
      "   0.25000000000000   0.10780000000000   0.82840000000000  O \n"
      "   0.75000000000000   0.89220000000000   0.82840000000000  O \n";

    check_function = "xstructure::operator<<";
    check_description = "Parsing CIF file with recognized setting (CrO3)";

    aurostd::StringstreamClean(xstrss);
    xstrss << str_cif;
    xstr_cif = xstructure(xstrss);

    aurostd::StringstreamClean(xstrss);
    xstrss << str_poscar;
    xstr_poscar = xstructure(xstrss);

    // ---------------------------------------------------------------------------
    // test: parse structure
    expected_bool = true;
    calculated_bool = compare::aflowCompareStructure(xstr_cif, xstr_poscar, same_species, scale_volume, optimize_match, misfit);
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // test: compare Wyckoff positions
    check_description = "Compare parsed Wyckoff positions of CrO3";
    vector<wyckoffsite_ITC> vwyckoff(4);
    xvector<double> coords;
    vwyckoff[0].type = "Cr"; vwyckoff[0].letter = "b"; vwyckoff[0].site_symmetry = "m.."; vwyckoff[0].multiplicity = 4; vwyckoff[0].coord[1] = 0.25; vwyckoff[0].coord[2] = 0.09676; vwyckoff[0].coord[3] = 0.5;
    vwyckoff[1].type = "O"; vwyckoff[1].letter = "a"; vwyckoff[1].site_symmetry = "..2"; vwyckoff[1].multiplicity = 4; vwyckoff[1].coord[1] = 0.0; vwyckoff[1].coord[2] = 0.0; vwyckoff[1].coord[3] = 0.3841;
    vwyckoff[2].type = "O"; vwyckoff[2].letter = "b"; vwyckoff[2].site_symmetry = "m.."; vwyckoff[2].multiplicity = 4; vwyckoff[2].coord[1] = 0.25; vwyckoff[2].coord[2] = 0.2677; vwyckoff[2].coord[3] = 0.3755;
    vwyckoff[3].type = "O"; vwyckoff[3].letter = "b"; vwyckoff[3].site_symmetry = "m.."; vwyckoff[3].multiplicity = 4; vwyckoff[3].coord[1] = 0.25; vwyckoff[3].coord[2] = 0.6078; vwyckoff[3].coord[3] = 0.3284;
    check_passed = (vwyckoff.size() == xstr_cif.wyckoff_sites_ITC.size());
    for (size_t i = 0; i < vwyckoff.size() && check_passed; i++) {
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].type == vwyckoff[i].type));
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].letter == vwyckoff[i].letter));
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].multiplicity == vwyckoff[i].multiplicity));
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].coord == vwyckoff[i].coord));
      if (LDEBUG && !check_passed) {
        std::cerr << "Failed site:" << std::endl;
        std::cerr << "type = " << xstr_cif.wyckoff_sites_ITC[i].type << ", ";
        std::cerr << "letter = " << xstr_cif.wyckoff_sites_ITC[i].letter << ", ";
        std::cerr << "multiplicity = " << xstr_cif.wyckoff_sites_ITC[i].multiplicity << ", ";
        std::cerr << "coord = " << xstr_cif.wyckoff_sites_ITC[i].coord << std::endl << std::endl;
        std::cerr << "Should be:" << std::endl;
        std::cerr << "type = " << vwyckoff[i].type << ", ";
        std::cerr << "letter = " << vwyckoff[i].letter << ", ";
        std::cerr << "multiplicity = " << vwyckoff[i].multiplicity << ", ";
        std::cerr << "coord = " << vwyckoff[i].coord << std::endl;
      }
    }

    expected_bool = true;
    calculated_bool = check_passed;
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);

    // May need a better test case where the labels actually change
    check_description = "Compare calculated Wyckoff positions of CrO3";
    xstr_cif.SpaceGroup_ITC();
    vwyckoff[0].coord[1] = 0.25; vwyckoff[0].coord[2] = 0.59676; vwyckoff[0].coord[3] = 0.0000;
    vwyckoff[1].coord[1] = 0.00; vwyckoff[1].coord[2] = 0.00000; vwyckoff[1].coord[3] = 0.3841;
    vwyckoff[2].coord[1] = 0.25; vwyckoff[2].coord[2] = 0.76770; vwyckoff[2].coord[3] = 0.8755;
    vwyckoff[3].coord[1] = 0.25; vwyckoff[3].coord[2] = 0.60780; vwyckoff[3].coord[3] = 0.3284;
    check_passed = (vwyckoff.size() == xstr_cif.wyckoff_sites_ITC.size());
    for (size_t i = 0; i < vwyckoff.size() && check_passed; i++) {
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].type == vwyckoff[i].type));
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].letter == vwyckoff[i].letter));
      check_passed = (check_passed && (xstr_cif.wyckoff_sites_ITC[i].multiplicity == vwyckoff[i].multiplicity));
      check_passed = (check_passed && aurostd::isequal(xstr_cif.wyckoff_sites_ITC[i].coord, vwyckoff[i].coord));
      if (LDEBUG && !check_passed) {
        std::cerr << "Failed site:" << std::endl;
        std::cerr << "type = " << xstr_cif.wyckoff_sites_ITC[i].type << ", ";
        std::cerr << "letter = " << xstr_cif.wyckoff_sites_ITC[i].letter << ", ";
        std::cerr << "multiplicity = " << xstr_cif.wyckoff_sites_ITC[i].multiplicity << ", ";
        std::cerr << "coord = " << xstr_cif.wyckoff_sites_ITC[i].coord << std::endl << std::endl;
        std::cerr << "Should be:" << std::endl;
        std::cerr << "type = " << vwyckoff[i].type << ", ";
        std::cerr << "letter = " << vwyckoff[i].letter << ", ";
        std::cerr << "multiplicity = " << vwyckoff[i].multiplicity << ", ";
        std::cerr << "coord = " << vwyckoff[i].coord << std::endl;
      }
    }

    expected_bool = true;
    calculated_bool = check_passed;
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Test 2: parse structure with unrecognized (old) settings
    // ---------------------------------------------------------------------------
    check_description = "Parsing CIF file with unrecognized setting (GePt3)";
    aurostd::StringstreamClean(xstrss);
    str_cif =
      "data_Ge8Pt24\n"
      "_pd_phase_name Ge8Pt24\n"
      "_cell_length_a 7.93\n"
      "_cell_length_b 7.767\n"
      "_cell_length_c 7.767\n"
      "_cell_angle_alpha 90.\n"
      "_cell_angle_beta 90.06\n"
      "_cell_angle_gamma 90.\n"
      "_cell_volume 478.39\n"
      "_cell_formula_units_Z 8\n"
      "_space_group_name_H-M_alt 'F 1 2/m 1'\n"
      "_space_group_IT_number 12\n"
      "loop_\n"
      "_space_group_symop_id\n"
      "_space_group_symop_operation_xyz\n"
      "1 '-x, -y, -z'\n"
      "2 'x, -y, z'\n"
      "3 '-x, y, -z'\n"
      "4 'x, y, z'\n"
      "5 '-x, -y+1/2, -z+1/2'\n"
      "6 'x, -y+1/2, z+1/2'\n"
      "7 '-x, y+1/2, -z+1/2'\n"
      "8 'x, y+1/2, z+1/2'\n"
      "9 '-x+1/2, -y, -z+1/2'\n"
      "10 'x+1/2, -y, z+1/2'\n"
      "11 '-x+1/2, y, -z+1/2'\n"
      "12 'x+1/2, y, z+1/2'\n"
      "13 '-x+1/2, -y+1/2, -z'\n"
      "14 'x+1/2, -y+1/2, z'\n"
      "15 '-x+1/2, y+1/2, -z'\n"
      "16 'x+1/2, y+1/2, z'\n"
      "loop_\n"
      "_atom_type_symbol\n"
      "_atom_type_oxidation_number\n"
      "Pt0+ 0\n"
      "Ge0+ 0\n"
      "loop_\n"
      "_atom_site_label\n"
      "_atom_site_type_symbol\n"
      "_atom_site_symmetry_multiplicity\n"
      "_atom_site_Wyckoff_symbol\n"
      "_atom_site_fract_x\n"
      "_atom_site_fract_y\n"
      "_atom_site_fract_z\n"
      "_atom_site_B_iso_or_equiv\n"
      "_atom_site_occupancy\n"
      "Pt1 Pt0+ 8 g 0. 0.2 0. . 1.\n"
      "Pt2 Pt0+ 8 h 0.25 0.25 0.25 . 1.\n"
      "Pt3 Pt0+ 8 i 0. 0. 0.3 . 1.\n"
      "Ge1 Ge0+ 8 i 0.25 0. 0. . 1.\n";

    str_poscar =
      "Ge8Pt24\n"
      "1.0000000\n"
      "7.93000000000000   0.00000000000000   0.00000000000000\n"
      "0.00000000000000   7.76700000000000   0.00000000000000\n"
      "-0.00813358189357   0.00000000000000   7.76699574126609\n"
      "Ge Pt\n"
      "8 24\n"
      "Direct(32) [A8B24]\n"
      "0.25000000000000   0.00000000000000   0.00000000000000  Ge\n"
      "0.75000000000000   0.00000000000000   0.00000000000000  Ge\n"
      "0.75000000000000   0.50000000000000   0.50000000000000  Ge\n"
      "0.25000000000000   0.50000000000000   0.50000000000000  Ge\n"
      "0.25000000000000   0.00000000000000   0.50000000000000  Ge\n"
      "0.75000000000000   0.00000000000000   0.50000000000000  Ge\n"
      "0.25000000000000   0.50000000000000   0.00000000000000  Ge\n"
      "0.75000000000000   0.50000000000000   0.00000000000000  Ge\n"
      "0.00000000000000   0.70000000000000   0.50000000000000  Pt\n"
      "0.50000000000000   0.80000000000000   0.50000000000000  Pt\n"
      "0.50000000000000   0.20000000000000   0.50000000000000  Pt\n"
      "0.50000000000000   0.30000000000000   0.00000000000000  Pt\n"
      "0.50000000000000   0.70000000000000   0.00000000000000  Pt\n"
      "0.25000000000000   0.25000000000000   0.25000000000000  Pt\n"
      "0.00000000000000   0.00000000000000   0.30000000000000  Pt\n"
      "0.00000000000000   0.20000000000000   0.00000000000000  Pt\n"
      "0.00000000000000   0.00000000000000   0.70000000000000  Pt\n"
      "0.00000000000000   0.50000000000000   0.20000000000000  Pt\n"
      "0.00000000000000   0.50000000000000   0.80000000000000  Pt\n"
      "0.50000000000000   0.00000000000000   0.20000000000000  Pt\n"
      "0.50000000000000   0.00000000000000   0.80000000000000  Pt\n"
      "0.50000000000000   0.50000000000000   0.70000000000000  Pt\n"
      "0.50000000000000   0.50000000000000   0.30000000000000  Pt\n"
      "0.25000000000000   0.75000000000000   0.25000000000000  Pt\n"
      "0.75000000000000   0.75000000000000   0.75000000000000  Pt\n"
      "0.75000000000000   0.25000000000000   0.75000000000000  Pt\n"
      "0.25000000000000   0.25000000000000   0.75000000000000  Pt\n"
      "0.75000000000000   0.25000000000000   0.25000000000000  Pt\n"
      "0.25000000000000   0.75000000000000   0.75000000000000  Pt\n"
      "0.75000000000000   0.75000000000000   0.25000000000000  Pt\n"
      "0.00000000000000   0.80000000000000   0.00000000000000  Pt\n"
      "0.00000000000000   0.30000000000000   0.50000000000000  Pt\n";

    bool quiet_tmp = XHOST.QUIET;
    XHOST.QUIET = !LDEBUG;  // Suppress warnings
    xstrss << str_cif;
    xstr_cif = xstructure(xstrss);
    XHOST.QUIET = quiet_tmp;

    // ---------------------------------------------------------------------------
    // test: parse structure
    aurostd::StringstreamClean(xstrss);
    xstrss << str_poscar;
    xstr_poscar = xstructure(xstrss);
    expected_bool = true;
    calculated_bool = compare::aflowCompareStructure(xstr_cif, xstr_poscar, same_species, scale_volume, optimize_match, misfit);
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);
  }

  void UnitTest::xstructureTest(uint& passed_checks, vector<vector<string> >& results, vector<string>& errors) {
    if (errors.size()) {} // Suppress compiler warnings
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    // Set up test environment
    string check_function = "", check_description = "";
    bool expected_bool = false, calculated_bool = false;
    uint expected_uint = 0, calculated_uint = 0;

    // ---------------------------------------------------------------------------
    // Check | getCoordinations() //CO20190520
    // ---------------------------------------------------------------------------
    xstructure xstr("aflowlib.duke.edu:AFLOWDATA/ICSD_WEB/FCC/Cl1Na1_ICSD_240599","CONTCAR.relax.vasp",IOAFLOW_AUTO);
    deque<deque<uint> > coordinations;
    xstr.GetCoordinations(coordinations);

    check_function = "xstructure::getCoordinations()";
    check_description = "Number of iatoms";
    calculated_uint = coordinations.size();
    expected_uint = 2;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    check_description = "Number of coordination environments atom 1 > 2";
    expected_bool = true;
    calculated_bool = (coordinations[0].size() > 2);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    check_description = "Number of coordination environments atom 2 > 2";
    expected_bool = true;
    calculated_bool = (coordinations[1].size() > 2);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    //first iatom
    //first shell
    check_description = "First shell atom 1";
    calculated_uint = coordinations[0][0];
    expected_uint = 6;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    //second shell
    check_description = "Second shell atom 1";
    calculated_uint = coordinations[0][1];
    expected_uint = 12;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    //second iatom
    //first shell
    check_description = "First shell atom 2";
    calculated_uint = coordinations[1][0];
    expected_uint = 6;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    //second shell
    check_description = "Second shell atom 2";
    calculated_uint = coordinations[1][1];
    expected_uint = 12;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | foldAtomdInCell //DX20210129
    // ---------------------------------------------------------------------------

    // set fold atoms in cell variables
    bool skew = false;
    double tol = 0.01;
    bool check_min_dists = false;

    // ---------------------------------------------------------------------------
    // test 1: expand cell
    // create 3x1x1 supercell expansion matrix
    check_function = "xstructure::foldAtomsInCell()";
    check_description = "expand cell (3x1x1)";
    xmatrix<double> supercell_matrix = aurostd::eye<double>(3,3);
    supercell_matrix(1,1)=3.0;
    xmatrix<double> lattice_new = supercell_matrix*xstr.lattice; // transform lattice

    xstructure xstr_supercell = xstr;
    xstr_supercell.foldAtomsInCell(lattice_new, skew, tol, check_min_dists);

    bool same_species = true;
    expected_bool = true;
    calculated_bool = compare::structuresMatch(xstr,xstr_supercell,same_species);
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // test 2: reduce cell
    // convert supercell back to original lattice
    check_function = "xstructure::foldAtomsInCell()";
    check_description = "reduce cell into primitive";
    xstructure xstr_reduced = xstr_supercell;
    xstr_reduced.foldAtomsInCell(xstr.lattice, skew, tol, check_min_dists);

    expected_bool = true;
    calculated_bool = compare::structuresMatch(xstr,xstr_reduced,same_species);
    checkEqual(expected_bool, calculated_bool, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | slab test //CO20190520
    // ---------------------------------------------------------------------------
    //See W. Sun and G. Ceder, Surface Science 617 (2013) 53-59
    string xstr_str = "";
    stringstream xstrss;
    double min_dist = 0.0, min_dist_orig = 0.0;
    xvector<int> hkl(3);
    hkl[1] = 1; hkl[2] = 0; hkl[3] = 4;

    //create input structure
    xstr_str =
      "FeO\n"
      "1.0\n"
      "  0.000000   4.736235   0.000000\n"
      "  4.101698  -2.368118   0.000000\n"
      "  0.000000   0.000000  13.492372\n"
      "12 18\n"
      "Direct\n"
      "   0.33333333333332   0.66666666666665   0.31943934568918  Fe\n"
      "  -0.00000000000002  -0.00000000000001   0.65277267902252  Fe\n"
      "   0.66666666666665   0.33333333333332   0.98610601235585  Fe\n"
      "   0.66666666666668   0.33333333333335   0.18056065431082  Fe\n"
      "   0.33333333333335   0.66666666666668   0.51389398764415  Fe\n"
      "   0.00000000000002   0.00000000000001   0.84722732097748  Fe\n"
      "  -0.00000000000001  -0.00000000000001   0.15277363902251  Fe\n"
      "   0.66666666666665   0.33333333333332   0.48610697235585  Fe\n"
      "   0.33333333333332   0.66666666666665   0.81944030568918  Fe\n"
      "   0.33333333333335   0.66666666666668   0.01389302764415  Fe\n"
      "   0.00000000000001   0.00000000000001   0.34722636097749  Fe\n"
      "   0.66666666666668   0.33333333333335   0.68055969431082  Fe\n"
      "   0.00000000000000   0.31486811660227   0.25000000000000  O \n"
      "   0.66666666666667   0.64820144993561   0.58333333333333  O \n"
      "   0.33333333333333   0.98153478326894   0.91666666666667  O \n"
      "   0.68513188339780   0.68513188339776   0.24999999999999  O \n"
      "   0.35179855006446   0.01846521673110   0.58333333333332  O \n"
      "   0.01846521673113   0.35179855006443   0.91666666666665  O \n"
      "   0.31486811660220  -0.00000000000004   0.25000000000001  O \n"
      "   0.98153478326887   0.33333333333330   0.58333333333335  O \n"
      "   0.64820144993554   0.66666666666663   0.91666666666668  O \n"
      "   0.35179771006447   0.33333333333337   0.08333333333332  O \n"
      "   0.01846437673113   0.66666666666670   0.41666666666665  O \n"
      "   0.68513104339780   0.00000000000004   0.74999999999999  O \n"
      "   0.98153562326887   0.64820228993557   0.08333333333335  O \n"
      "   0.64820228993553   0.98153562326890   0.41666666666668  O \n"
      "   0.31486895660220   0.31486895660224   0.75000000000001  O \n"
      "   0.66666666666667   0.01846437673106   0.08333333333333  O \n"
      "   0.33333333333333   0.35179771006440   0.41666666666667  O \n"
      "   0.00000000000000   0.68513104339773   0.75000000000000  O \n";
    aurostd::StringstreamClean(xstrss);
    xstrss << xstr_str;
    xstructure xstr_in;
    try {
      xstrss >> xstr_in;          //CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
    } catch (aurostd::xerror& excpt) {} //CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
    min_dist = min_dist_orig = xstr_in.MinDist();
    if(LDEBUG){
      std::cerr << __AFLOW_FUNC__ << " xstr_in=\n" << xstr_in << std::endl;
      std::cerr << __AFLOW_FUNC__ << " xstr_in.MinDist()=" << min_dist << std::endl;
    }

    //create xstr_slab (correct answer)
    xstr_str =
    "FeO\n"
    "1.0\n"
    " -4.73623366665202   0.00000000000000   0.00000000000000\n"
    " -9.47247466669728  21.24210623923067   0.00000000000000\n"
    " -2.36811866667432   3.16803536641816   2.60528076661565\n"
    "12 18\n"
    "Direct\n"
    "   0.66666667      0.31943935     0.38890928   Fe\n"
    "  -0.00000000      0.65277268     0.38890928   Fe\n"
    "   0.33333333      0.98610601     0.38890928   Fe\n"
    "   0.33333333      0.18056065     0.61109072   Fe\n"
    "   0.66666667      0.51389399     0.61109072   Fe\n"
    "   0.00000000      0.84722732     0.61109072   Fe\n"
    "  -0.00000000      0.15277364     0.38890544   Fe\n"
    "   0.33333333      0.48610697     0.38890544   Fe\n"
    "   0.66666667      0.81944031     0.38890544   Fe\n"
    "   0.66666667      0.01389303     0.61109456   Fe\n"
    "   0.00000000      0.34722636     0.61109456   Fe\n"
    "   0.33333333      0.68055969     0.61109456   Fe\n"
    "   0.31486812      0.25000000     0.00000000   O \n"
    "   0.64820145      0.58333333     0.00000000   O \n"
    "   0.98153478      0.91666667     0.00000000   O \n"
    "   0.68513188      0.25000000     0.31486812   O \n"
    "   0.01846522      0.58333333     0.31486812   O \n"
    "   0.35179855      0.91666667     0.31486812   O \n"
    "  -0.00000000      0.25000000     0.68513188   O \n"
    "   0.33333333      0.58333333     0.68513188   O \n"
    "   0.66666667      0.91666667     0.68513188   O \n"
    "   0.33333333      0.08333333     0.31486896   O \n"
    "   0.66666667      0.41666667     0.31486896   O \n"
    "   0.00000000      0.75000000     0.31486896   O \n"
    "   0.64820229      0.08333333     0.68513104   O \n"
    "   0.98153562      0.41666667     0.68513104   O \n"
    "   0.31486896      0.75000000     0.68513104   O \n"
    "   0.01846438      0.08333333     0.00000000   O \n"
    "   0.35179771      0.41666667     0.00000000   O \n"
    "   0.68513104      0.75000000     0.00000000   O \n";
    xstructure xstr_slab_correct;
    aurostd::StringstreamClean(xstrss);
    xstrss << xstr_str;
    try {
      xstrss >> xstr_slab_correct; //CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
    } catch (aurostd::xerror& excpt) {} //CO20200404 - this WILL throw an error because det(lattice)<0.0, leave alone
    // ---------------------------------------------------------------------------
    // test 1: compare min distance of bulk and correct slab
    min_dist = xstr_slab_correct.MinDist();
    if(LDEBUG) {
      std::cerr << __AFLOW_FUNC__ << " xstr_slab_correct=\n" << xstr_slab_correct << std::endl;
      min_dist=xstr_slab_correct.MinDist();
      std::cerr << __AFLOW_FUNC__ << " xstr_slab_correct.MinDist()=" << min_dist << std::endl;
    }
    check_function = "xstructure::CreateSlab_SurfaceLattice()";
    check_description = "Compare minimum distance bulk vs. correct slab";
    checkEqual(min_dist, min_dist_orig, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // test 2: compare min distance of generated slab and correct slab
    xmatrix<double> lattice_slab_origbasis;
    xstructure xstr_slab_test = slab::CreateSlab_SurfaceLattice(xstr_in,hkl,1,0,5.0); //the v3len_max_strict is very important here, as the test from Sun et al. takes a shortcut here
    min_dist = xstr_slab_test.MinDist();
    check_function = "xstructure::CreateSlab_SurfaceLattice(hkl = 104)";
    check_description = "Compare minimum distance slab vs. correct slab";
    if(LDEBUG){
      std::cerr << __AFLOW_FUNC__ << " xstr_slab_test=\n" << xstr_slab_test << std::endl;
      min_dist=xstr_slab_test.MinDist();
      std::cerr << __AFLOW_FUNC__ << " xstr_slab_test.MinDist()=" << min_dist << std::endl;
    }
    checkEqual(min_dist, min_dist_orig, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // test 3: compare structures of generated slab and correct slab
    check_description = "Match structures of generated slab and correct slab";
    expected_bool = true;
    calculated_bool = compare::structuresMatch(xstr_slab_correct,xstr_slab_test,true,false,false);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

  }

}

// structure generation
namespace unittest {

  //CO20190520
  void UnitTest::ceramgenTest(uint& passed_checks, vector<vector<string> >& results, vector<string>& errors) {
    if (errors.size()) {} // Suppress compiler warnings

    // setup test environment
    string check_function = "", check_description = "";
    bool calculated_bool = false, expected_bool = false;
    uint calculated_uint = 0, expected_uint = 0;

    //./aflow --generate_ceramics --nm=N,C --m=Co,Mo,Fe,Ru,Ni,Rh,Pt,Cu,Cr,V --N=5
    vector<string> vnonmetals,vmetals;
    aurostd::string2tokens("N,C",vnonmetals,",");
    aurostd::string2tokens("Co,Mo,Fe,Ru,Ni,Rh,Pt,Cu,Cr,V",vmetals,",");

    check_function = "pflow::GENERATE_CERAMICS()";
    vector<string> commands=pflow::GENERATE_CERAMICS(vnonmetals,vmetals,5);

    check_description = "number of commands";
    calculated_uint = commands.size();
    expected_uint = 6;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N
    //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V
    //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
    //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
    //C:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
    //C:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V

    string search = "";

    search = "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N";
    check_description = "find " + search;
    expected_bool = true;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    search = "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V";
    check_description = "find " + search;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    search = "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V";
    check_description = "find " + search;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    search = "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V";
    check_description = "find " + search;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    search = "C:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V";
    check_description = "find " + search;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    search = "C:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V";
    check_description = "find " + search;
    calculated_bool = aurostd::WithinList(commands, search);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);
  }

  //DX20200925
  //ME20220324 - refactored to run in parallel
  void _testPrototype(uint i, const vector<string>& prototype_labels, vector<uint>& nprotos, vector<string>& errors
#ifdef AFLOW_MULTITHREADS_ENABLE
    , std::mutex& m
#endif
  ) {
    double tolerance_sym = 0.0;
    string label_input = "";
    bool generated = false, sym = false, unique = false;
    xstructure xstr;
    ofstream ofs("/dev/null");

    vector<string> parameter_sets = anrl::getANRLParameters(prototype_labels[i], "all");
    for (size_t j = 0; j < parameter_sets.size(); j++) {
      string error = "";
      try {
        xstr = aflowlib::PrototypeLibraries(ofs,prototype_labels[i],parameter_sets[j],1);
        generated = true;
      } catch(aurostd::xerror& excpt) {
        error = "Could not generate prototype=" + prototype_labels[i] + ", params=" + parameter_sets[j];
      }

      if (error.empty()) {
        stringstream label_input_ss;
        label_input_ss << prototype_labels[i] << "-" << std::setw(3) << std::setfill('0') << j+1;
        label_input = label_input_ss.str();
        tolerance_sym = anrl::specialCaseSymmetryTolerances(label_input);
        if(tolerance_sym != AUROSTD_MAX_DOUBLE) {
          xstr.sym_eps = tolerance_sym;
          xstr.sym_eps_calculated = true;
        }
      }

      if (error.empty()) {
        string updated_label_and_params = "";
        if(!anrl::structureAndLabelConsistent(xstr, prototype_labels[i], updated_label_and_params, tolerance_sym)){ //DX20201105 - added symmetry tolerance
          error = "The structure has a higher symmetry than indicated by the label (orig: proto="
              + prototype_labels[i] + ", params=" + parameter_sets[j] + ")."
              + " The correct label and parameters for this structure are:\n" + updated_label_and_params; 
        } else {
          sym = true;
        }
      }
      if (error.empty()) {
        aurostd::xoption vpflow;
        // check if the prototype matches to more than one prototype
        // (i.e., a prototype should match with itself, but no others)
        string catalog = "anrl";
        vector<string> protos_matching = compare::getMatchingPrototypes(xstr, vpflow, catalog);
        // if it matches to more than one
        if(protos_matching.size() > 1 && !anrl::isSpecialCaseEquivalentPrototypes(protos_matching)) {
          error = prototype_labels[i] + ", params=" + parameter_sets[j]
              + " matches multiple prototypes (and not a documented special case): "
              + aurostd::joinWDelimiter(protos_matching,",") + "."
              + " If the prototype was newly added, ONLY include it in the encyclopedia"
              + " for a valid reason (e.g., historical, special designation, etc.)"
              + " and document this in anrl::isSpecialCaseEquivalentPrototypes().";
        // if it doesn't match with ITSELF
        } else if (protos_matching.size() == 0) {
          error = prototype_labels[i] + ", params=" + parameter_sets[j]
              + " does not match to any prototypes"
              + " (requires special symmetry tolerance or there is a bug with XtalFinder).";
        } else {
          unique = true;
        }
      }
#ifdef AFLOW_MULTITHREADS_ENABLE
      std::lock_guard<std::mutex> lk(m);
#endif
      if (generated) nprotos[0]++;
      if (sym) nprotos[1]++;
      if (unique) nprotos[2]++;
      if (!error.empty()) errors.push_back(error);
    }
  }

  void UnitTest::prototypeGeneratorTest(uint& passed_checks, vector<vector<string> >& results, vector<string>& errors) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    // Set up test environment
    string check_function = "", check_description = "";

    vector<string> prototype_labels, compositions;
    vector<uint> space_group_numbers;
    vector<vector<vector<string> > > grouped_Wyckoff_letters;
    string library = "anrl";

    uint num_protos = aflowlib::GetAllPrototypeLabels(prototype_labels,
        compositions,
        space_group_numbers,
        grouped_Wyckoff_letters,
        library);
    if (LDEBUG) std::cerr << __AFLOW_FUNC__ << "Number of prototype labels = " << num_protos << " (each may have multiple parameter sets)";

    // Test
    // 1: if the prototype can be generated,
    // 2: if symmetry and label are consistent
    // 3: if it is a unique prototye
    // Keep results in vector to simplify function input
    vector<uint> nprotos(3, 0);
#ifdef AFLOW_MULTITHREADS_ENABLE
    xthread::xThread xt(KBIN::get_NCPUS()); // Okay to be greedy - xThread will manage number of threads
    std::mutex m;
    xt.run(num_protos, _testPrototype, prototype_labels, nprotos, errors, m);
#else
    for (uint i = 0; i < num_protos; i++) _testPrototype(i, prototype_labels, nprotos, errors);
#endif

    // Get number of all protoypes + parameter sets
    uint expected_uint = 0;
    for (uint i = 0; i < num_protos; i++) {
      expected_uint += anrl::getANRLParameters(prototype_labels[i], "all").size();
    }

    check_function = "aflowlib::PrototypeLibraries()";
    check_description = "generate prototypes";
    checkEqual(nprotos[0], expected_uint, check_function, check_description, passed_checks, results);

    check_function = "anrl::structureAndLabelConsistent()";
    check_description = "symmetry consistent with prototype label";
    checkEqual(nprotos[1], expected_uint, check_function, check_description, passed_checks, results);

    check_function = "compare::getMatchingPrototypes()";
    check_description = "protoypes are unique";
    checkEqual(nprotos[2], expected_uint, check_function, check_description, passed_checks, results);
  }

}

// ovasp
namespace unittest {
  void UnitTest::xoutcarTest(uint& passed_checks, vector<vector<string> >& results, vector<string>& errors) {
    bool LDEBUG = (FALSE || XHOST.DEBUG);

    // setup test environment
    string check_function = "", check_description = "";
    string calculated_str = "", expected_str = "";
    bool calculated_bool = false, expected_bool = false;
    double calculated_dbl = 0.0, expected_dbl = 0.0;

    string system="",path="",query="",file="",efile="",ext="",tfile="",Egap_type="";
    vector<string> files;
    xOUTCAR xout;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //FCC/Si1_ICSD_150530
    system="ICSD_WEB/FCC/Si1_ICSD_150530";

    // Fetch and parse required files first - abort test if unsuccessful
    path=AFLOWLIB_SERVER_DEFAULT+"/AFLOWDATA/"+system;
    query=path+"/?files";
    aurostd::url2tokens(query,files,",");
    if(files.size()==0){
      errors.push_back("Could not fetch query: " + query);
      return;
    }

    //OUTCAR.static
    file="OUTCAR.static";
    if(!aurostd::EWithinList(files,file,efile)){
      errors.push_back("No " + file + " found within " + query);
      return;
    }

    ext=aurostd::GetCompressionExtension(efile);
    tfile=aurostd::TmpFileCreate("Egap_file1")+ext;
    query=path+"/"+efile;
    if(!(aurostd::url2file(query,tfile,LDEBUG) && aurostd::FileExist(tfile))){
      errors.push_back("Could not fetch query: " + query);
      return;
    }

    check_function = "xOUTCAR::GetProperties()";
    check_description = "load OUTCAR.static";
    if (!xout.GetPropertiesFile(tfile,!LDEBUG)) {
      errors.push_back("Could not parse " + file);
      return;
    }
#ifndef _AFLOW_TEMP_PRESERVE_
    aurostd::RemoveFile(tfile);
#endif
    double EFERMI=xout.Efermi;

    //OUTCAR.bands
    file="OUTCAR.bands";
    if(!aurostd::EWithinList(files,file,efile)){
      query=path+"/?files"; //reload query for error message
      errors.push_back("No " + file + " found within " + query);
      return;
    }
    ext=aurostd::GetCompressionExtension(efile);
    tfile=aurostd::TmpFileCreate("Egap_file1")+ext;
    query=path+"/"+efile;
    if(!(aurostd::url2file(query,tfile,LDEBUG) && aurostd::FileExist(tfile))){
      errors.push_back("Could not fetch query: " + query);
      return;
    }
    if(!xout.GetPropertiesFile(tfile,!LDEBUG)){
      errors.push_back("Could not parse " + file);
      return;
    }
#ifndef _AFLOW_TEMP_PRESERVE_
    aurostd::RemoveFile(tfile);
#endif

    //GetBandGap
    check_function = "xOUTCAR::GetBandGap()";

    check_description = "Calculate Egap successfully";
    expected_bool = true;
    bool parsed = calculated_bool = xout.GetBandGap(EFERMI);
    checkEqual(calculated_bool, expected_bool, check_function, check_description, passed_checks, results);

    check_description = "Egap value";
    expected_dbl = 6.1000e-01;
    calculated_dbl = (parsed?xout.Egap[0]:AUROSTD_MAX_DOUBLE);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    check_description = "Egap type";
    expected_str = "insulator-indirect";
    calculated_str = (parsed?xout.Egap_type[0]:"N/A");
    checkEqual(calculated_str, expected_str, check_function, check_description, passed_checks, results);
  }

}

bool smithTest(ostream& oss){ofstream FileMESSAGE;return smithTest(FileMESSAGE,oss);}  //CO20190520
bool smithTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  xmatrix<int> A1(3,3),U1,V1,S1;
  A1[1][1]=3;A1[1][2]=2;A1[1][3]=1;
  A1[2][1]=5;A1[2][2]=3;A1[2][3]=1;
  A1[3][1]=6;A1[3][2]=8;A1[3][3]=9;

  aurostd::getSmithNormalForm(A1,U1,V1,S1);

  if(LDEBUG){
    cerr << __AFLOW_FUNC__ << " A=" << endl;cerr << A1 << endl;
    cerr << __AFLOW_FUNC__ << " U=" << endl;cerr << U1 << endl;
    cerr << __AFLOW_FUNC__ << " V=" << endl;cerr << V1 << endl;
    cerr << __AFLOW_FUNC__ << " S=" << endl;cerr << S1 << endl;
  }

  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]if(!(
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      U1[1][1]==24 && U1[1][2]==-13 && U1[1][3]==-1 &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      U1[2][1]==13 && U1[2][2]==-7  && U1[2][3]==-1 &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      U1[3][1]==2  && U1[3][2]==-1  && U1[3][3]==0  &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      TRUE
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]    )
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  ){
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  if(LDEBUG){cerr << __AFLOW_FUNC__ << " U1(1) failed of getSmithNormalForm()" << endl;}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  return FALSE;
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]if(!(
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      V1[1][1]==0  && V1[1][2]==1  && V1[1][3]==3  &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      V1[2][1]==-1 && V1[2][2]==-1 && V1[2][3]==-1 &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      V1[3][1]==1  && V1[3][2]==0  && V1[3][3]==-1 &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      TRUE
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]    )
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  ){
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  if(LDEBUG){cerr << __AFLOW_FUNC__ << " V1(1) failed of getSmithNormalForm()" << endl;}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  return FALSE;
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]if(!(
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      S1[1][1]==1 && S1[1][2]==0 && S1[1][3]==0 &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      S1[2][1]==0 && S1[2][2]==1 && S1[2][3]==0 &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      S1[3][1]==0 && S1[3][2]==0 && S1[3][3]==1 &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      TRUE
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]    )
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  ){
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  if(LDEBUG){cerr << __AFLOW_FUNC__ << " S1(1) failed of getSmithNormalForm()" << endl;}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  return FALSE;
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]}

  xmatrix<long long int> A2(5,5),U2,V2,S2;  //long long int is CRUCIAL, Matlab actually gets this wrong because it uses long int by default
  A2[1][1]=25;    A2[1][2]=-300;   A2[1][3]=1050;    A2[1][4]=-1400;   A2[1][5]=630;
  A2[2][1]=-300;  A2[2][2]=4800;   A2[2][3]=-18900;  A2[2][4]=26880;   A2[2][5]=-12600;
  A2[3][1]=1050;  A2[3][2]=-18900; A2[3][3]=79380;   A2[3][4]=-117600; A2[3][5]=56700;
  A2[4][1]=-1400; A2[4][2]=26880;  A2[4][3]=-117600; A2[4][4]=179200;  A2[4][5]=-88200;
  A2[5][1]=630;   A2[5][2]=-12600; A2[5][3]=56700;   A2[5][4]=-88200;  A2[5][5]=44100;

  aurostd::getSmithNormalForm(A2,U2,V2,S2);

  if(LDEBUG){ //COME BACK AND PATCH FOR ANSWERS
    cerr << __AFLOW_FUNC__ << " A=" << endl;cerr << A2 << endl;
    cerr << __AFLOW_FUNC__ << " U=" << endl;cerr << U2 << endl;
    cerr << __AFLOW_FUNC__ << " V=" << endl;cerr << V2 << endl;
    cerr << __AFLOW_FUNC__ << " S=" << endl;cerr << S2 << endl;
  }

  message << "smith test successful";pflow::logger(_AFLOW_FILE_NAME_,__AFLOW_FUNC__,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
  return TRUE; //CO20180419
}

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
