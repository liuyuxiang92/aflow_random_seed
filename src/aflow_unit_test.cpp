// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************

#include "aflow.h"
#include "aflow_anrl.h"  //DX20201104
#include "aflow_compare_structure.h"  //ME20220125

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
  }

  void UnitTest::copy(const UnitTest& ut) {
    if (this == &ut) return;
    aflags = ut.aflags;
    test_functions = ut.test_functions;
  }

  void UnitTest::initialize() {
    free();
    string dir = XHOST.vflag_control.getattachedscheme("DIRECTORY");
    aflags.Directory = (!dir.empty()?dir:aurostd::getPWD());

    initializeTestFunctions();
    initializeTestGroups();
  }

  void UnitTest::initializeTestFunctions() {
    // Initialize unit tests
    xcheck xchk;

    // aurostd
    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::xscalarTest, this, _1, _2, _3);
    xchk.function_name = XPID + "xscalarTest():";
    xchk.task_description = "Testing xscalar";
    test_functions["schema"] = xchk;

    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::xvectorTest, this, _1, _2, _3);
    xchk.function_name = XPID + "xvectorTest():";
    xchk.task_description = "Testing xvector";
    test_functions["schema"] = xchk;

    // database
    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::schemaTest, this, _1, _2, _3);
    xchk.function_name = XPID + "SchemaTest():";
    xchk.task_description = "Testing schema";
    test_functions["schema"] = xchk;

    // xstructure
    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::atomicEnvironmentTest, this, _1, _2, _3);
    xchk.function_name = XPID + "AtomicEnvironmentTest():";
    xchk.task_description = "Creating atomic environments";
    test_functions["atomic_environment"] = xchk;

    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::atomicEnvironmentTest, this, _1, _2, _3);
    xchk.function_name = XPID + "coordinationTest():";
    xchk.task_description = "Testing coordination numbers";
    test_functions["coordination"] = xchk;

    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::cifParserTest, this, _1, _2, _3);
    xchk.function_name = XPID + "cifParserTest():";
    xchk.task_description = "Testing CIF file parser";
    test_functions["cif_parser"] = xchk;

    xchk = initializeXCheck();
    xchk.func = std::bind(&UnitTest::foldAtomsInCellTest, this, _1, _2, _3);
    xchk.function_name=XPID+"foldAtomsInCellTest():";
    xchk.task_description = "Folding atoms into cell";
    test_functions["fold_atoms"] = xchk;
  }

  xcheck UnitTest::initializeXCheck() {
    xcheck xt;
    resetUnitTest(xt);
    xt.func = nullptr;
    xt.function_name = "";
    xt.task_description = "";
    return xt;
  }

  void UnitTest::resetUnitTest(const string& test_name) {
    if (test_functions.count(test_name)) {
      resetUnitTest(test_functions[test_name]);
    }
  }

  void UnitTest::resetUnitTest(xcheck& test) {
    test.errors.clear();
    test.finished = false;
    test.passed_checks = 0;
    test.results.clear();
  }

  void UnitTest::initializeTestGroups() {
    test_groups.clear();
    test2group.clear();

    test_groups["aurostd"] = {"xscalar", "xvector"};
    test_groups["database"] = {"schema"};
    test_groups["xstructure"] = {"atomic_environment", "coordination", "cif_parser", "fold_atoms"};

    for (const auto& group : test_groups) {
      for (const string& member : group.second) {
        test2group[member] = group.first;
      }
    }
  }

}

namespace unittest {

  bool UnitTest::runTestSuites(const string& unit_test) {
    vector<string> vunit_test = {unit_test};
    return runTestSuites(vunit_test);
  }

  bool UnitTest::runTestSuites(const vector<string>& unit_tests_in) {
    string function_name = XPID + "runUnitTests():";
    stringstream  message;
    // Create task lists (groups or individual tests)
    // unit_test is the individual small tests over
    // which to parallelize
    vector<string> unit_tests, tasks;
    for (const string& test : unit_tests_in) {
      bool isgroup = (test_groups.find(test) != test_groups.end());
      if (test == "all") {
        tasks.clear();
        for (const auto& group : test_groups) {
          tasks.push_back(group.first);
          for (const string& member : group.second ) {
            unit_tests.push_back(member);
          }
        }
        break;
      } else if (!isgroup && (test_functions.find(test) != test_functions.end())) {
        message << "Skipping unrecognized test name " << test << ".";
        pflow::logger(_AFLOW_FILE_NAME_, function_name, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      } else if (isgroup && !aurostd::WithinList(tasks, test)) {
        tasks.push_back(test);
        for (const string& member : test_groups[test]) {
          unit_tests.push_back(member);
        }
      } else if (!aurostd::WithinList(tasks, test2group[test])) {
        unit_tests.push_back(test);
        tasks.push_back(test);
      }
    }
    uint ntasks = tasks.size();
    if (ntasks == 0) {
      message << "No unit tests to run.";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_NOTICE_);
      return true;
    }

    // Run
#ifdef AFLOW_MULTITHREADS_ENABLE
    std::mutex mtx;
    xthread::xThread xt(KBIN::get_NCPUS());
    std::function<void(vector<string>::iterator&, const vector<string>&)> fn = std::bind(&UnitTest::runUnitTest, this, _1, _2);
    xt.run(unit_tests, fn, tasks);
#else
    for (vector<string>::iterator it = unit_tests.begin(); it != unit_tests.end(); ++it) runUnitTest(it, tasks);
#endif

    // Print final summary
    uint nsuccess = 0;
    stringstream summary;
    for (const string& task : tasks) {
      bool success = taskSuccessful(task);
      if (success) nsuccess++;
      summary << "\t" << task << " | " << (success?"pass":"fail") << "\n";
    }

    if (nsuccess == ntasks) {
      message << "Unit tests passed successfully (passsing " << ntasks << " tests).";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
    } else {
      message << "Some unit tests failed (" << (ntasks - nsuccess) << " of " << ntasks << " failed).";
      pflow::logger(_AFLOW_FILE_NAME_, function_name, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
    }
    pflow::logger(_AFLOW_FILE_NAME_, function_name, summary, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_RAW_);
    return (nsuccess == ntasks);
  }

  void UnitTest::runUnitTest(vector<string>::iterator& it, const vector<string>& tasks) {
    string function_name = XPID + "UnitTest::runUnitTest():";
    const string& test_name = (*it);
    xcheck& test = test_functions[test_name];
    resetUnitTest(test);
    test.func(test.passed_checks, test.results, test.errors);
    std::lock_guard<std::mutex> lk(mtx);
    test.finished = true;
    // Output results
    if (aurostd::WithinList(tasks, test_name)) {
      // If the test name is in the task list, it is not part
      // of a group, so no need to check if other members are done
      displayResult(test);
    } else {
      // Test if part of a test group, so check if all members
      // of the group have finished before producing output
      const string& group = test2group[test_name];
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
        for (const string& test_name_group : vtests_group) {
          if (taskSuccessful(test_name_group)) nsuccess++;
        }
        stringstream message;
        if (nsuccess == ntests_group) {
          message << "Unit tests of group " << group << " passed successfully (passsing " << ntests_group << " tests).";
          pflow::logger(_AFLOW_FILE_NAME_, function_name, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_COMPLETE_);
        } else {
          message << "Some unit tests of group " << group << " failed (" << (ntests_group - nsuccess) << " of " << ntests_group << " failed).";
          pflow::logger(_AFLOW_FILE_NAME_, function_name, message, aflags, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        }
        for (const string& test_name_group : vtests_group) {
          displayResult(test_functions[test_name_group]);
        }
      }
    }
  }

  bool UnitTest::taskSuccessful(const string& task) {
    auto it = test_groups.find(task);
    if (it != test_groups.end()) {
      const vector<string> members = (*it).second;
      for (const string& member : members ) {
        const xcheck& xchk = test_functions[member];
        if (!xchk.finished || (test_functions[member].passed_checks != test_functions[member].results.size())) return false;
      }
      return true;
    } else {
      const xcheck& xchk = test_functions[task];
      return (xchk.finished && (xchk.passed_checks == xchk.results.size()));
    }
  }

}

namespace unittest {
  // Collection of generic check functions, to streamline testing.
  // HE20210616
  void UnitTest::displayResult(const xcheck& xchk) {
    stringstream message;
    uint check_num = xchk.results.size();
    if (xchk.passed_checks == check_num) {
      message << "SUCCESS " << xchk.task_description << " (passing " << check_num << " checks)" << endl;
      pflow::logger(_AFLOW_FILE_NAME_,xchk.function_name,message,aflags,*p_FileMESSAGE,*p_oss,_LOGGER_COMPLETE_);
    } else {
      message << "FAIL " << xchk.task_description << " (" << (check_num - xchk.passed_checks) << " of " << check_num << " checks failed)" << endl;
      pflow::logger(_AFLOW_FILE_NAME_,xchk.function_name,message,aflags,*p_FileMESSAGE,*p_oss,_LOGGER_ERROR_);
    }
    message << "\t" << aurostd::joinWDelimiter(xchk.results, "\n\t");
    pflow::logger(_AFLOW_FILE_NAME_,xchk.function_name,message,aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    if (xchk.errors.size() > 0) {
      message << "\nAdditional error messages:\n" << aurostd::joinWDelimiter(xchk.errors, "\n");
      pflow::logger(_AFLOW_FILE_NAME_,xchk.function_name,message,aflags,*p_FileMESSAGE,*p_oss,_LOGGER_RAW_);
    }
  }

  template <typename utype>
  void UnitTest::check(const bool passed, const utype &calculated, const utype &expected, const string &check_function,
      const string check_description, uint &passed_checks, vector<string> &results){
    stringstream result;
    uint check_num = results.size() + 1;
    if (passed) {
      passed_checks++;
      if (check_function.empty()){
        result << std::setw(3) << check_num << " | pass | " << check_description;
      } else {
        result << std::setw(3) << check_num << " | pass | " << check_function << " | " << check_description;
      }
    }
    else {
      if (check_function.empty()) {
        result << std::setw(3) << check_num << " | FAIL | " << check_description
          << " (result: " << calculated << " | expected: " << expected << ")";
      } else {
        result << std::setw(3) << check_num << " | FAIL | " << check_function << " | " << check_description
          << " (result: " << calculated << " | expected: " << expected << ")";
      }
    }
    results.push_back(result.str());
  }

  template <typename utype>
  void UnitTest::checkEqual(const utype &calculated, const utype &expected, const string &check_function,
      const string check_description, uint &passed_checks, vector<string> &results){
    bool passed = (aurostd::isequal(calculated, expected));
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }
  void UnitTest::checkEqual(const string &calculated, const string &expected, const string &check_function,
      const string check_description, uint &passed_checks, vector<string> &results){
    bool passed = (calculated == expected);
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }
  void UnitTest::checkEqual(const bool calculated, const bool expected, const string &check_function,
      const string check_description, uint &passed_checks, vector<string> &results){
    bool passed = (calculated == expected);
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }

  void UnitTest::checkSimilar(const double calculated, const double expected, const string &check_function,
      const string &check_description, uint &passed_checks, vector<string> &results, const double relative){
    bool passed = (std::abs(expected - calculated) <= expected * relative);
    check(passed, calculated, expected, check_function, check_description, passed_checks, results);
  }

}

// aurostd
namespace unittest {

  void UnitTest::xscalarTest(uint& passed_checks, vector<string>& results, vector<string>& errors) {
    if (errors.size()) {}  // Suppress compiler warnings
    // setup test environment
    string check_function = "";
    string check_description = "";

    // ---------------------------------------------------------------------------
    // Check | double2fraction conversion //DX20210908
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------
    check_function = "aurostd::double2fraction()";
    check_description = "convert a double to a fraction.";

    double test_double = 1.625;
    int numerator=1, denominator=1;
    string answer = "13/8";
    aurostd::double2fraction(test_double,numerator,denominator);
    stringstream result_ss; result_ss << numerator << "/" << denominator;

    checkEqual(result_ss.str(), answer, check_function, check_description, passed_checks, results);
  }

  void UnitTest::xvectorTest(uint& passed_checks, vector<string>& results, vector<string>& errors) {
    if (errors.size()) {}  // Suppress compiler warnings
    // setup test environment
    string check_function = "";
    string check_description = "";

    //HE20210511
    double expected_dbl = 0.0;
    double calculated_dbl = 0.0;

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
    check_function = "aurostd::volume()";
    check_description = "convex solid, points as doubles";
    expected_dbl = 2.5;

    calculated_dbl = aurostd::volume(points, facets, true);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | convex solid volume (int)
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
    check_function = "aurostd::volume()";
    check_description = "non convex solid, points as doubles";
    expected_dbl = 40.0;

    calculated_dbl = aurostd::volume(points, facets);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | error facet/normals mismatch
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
      if (e.error_code == expected_int) check(true, "", "", check_function, check_description, passed_checks, results);
      else check(false, aurostd::utype2string(e.error_code), expected_str, check_function, check_description, passed_checks, results);
    }
    catch (...) {
      check(false, std::string("not an xerror"), expected_str, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | non convex solid volume (int)
    check_function = "aurostd::volume()";
    check_description = "non convex solid, points as int";
    expected_dbl = 40.0;

    calculated_dbl = aurostd::volume(ipoints, facets);
    checkEqual(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | error facet size
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
      if (e.error_code == expected_int) check(true, "", "", check_function, check_description, passed_checks, results);
      else check(false, aurostd::utype2string(e.error_code), expected_str, check_function, check_description, passed_checks, results);
    }
    catch (...) {
      check(false, std::string("not an xerror"), expected_str, check_function, check_description, passed_checks, results);
    }

    // ---------------------------------------------------------------------------
    // Check | non convex area (double)
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
    check_function = "aurostd::areaPointsOnPlane()";
    check_description = "3d triangle; points as double";
    expected_dbl = 3.5355339059;

    calculated_dbl = aurostd::areaPointsOnPlane(points);
    checkSimilar(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | 3d triangle area (int)
    check_function = "aurostd::areaPointsOnPlane()";
    check_description = "3d triangle; points as int";
    expected_dbl = 3.5355339059;

    calculated_dbl = aurostd::areaPointsOnPlane(ipoints);
    checkSimilar(calculated_dbl, expected_dbl, check_function, check_description, passed_checks, results);
  }

}

// database
namespace unittest {

  void UnitTest::schemaTest(uint& passed_checks, vector<string>& results, vector<string>& errors) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "SchemaTest():";
    if (errors.size()) {}  // Suppress compiler warnings

    // Set up test environment
    string check_function = "";
    string check_description = "";

    check_function = "XHOST.vschema";
    check_description = "Internal consistency of vschema";
    vector<string> vschema_keys;
    vector<string> vschema_types = {"UNIT", "TYPE"};
    string key = "";
    uint ninconsistent = 0;
    for (uint i = 0; i < XHOST.vschema.vxsghost.size(); i+= 2) {
      if(XHOST.vschema.vxsghost[i].find("::NAME:") != string::npos) {
        key=aurostd::RemoveSubString(XHOST.vschema.vxsghost[i], "SCHEMA::NAME:");
        vschema_keys.push_back(XHOST.vschema.getattachedscheme("SCHEMA::NAME:" + key));
        for (uint j = 0; j < vschema_types.size(); j++) {
          if (!XHOST.vschema.isdefined("SCHEMA::" + vschema_types[j] + ":" + key)) {
            ninconsistent++;
            if (LDEBUG) std::cerr << function_name << " SCHEMA::" << vschema_types[j] << ":" << key << " not found." << std::endl;
          }
        }
      }
    }
    checkEqual(ninconsistent, 0, check_function, check_description, passed_checks, results);

    check_function = "_aflowlib_entry";
    check_description = "Consistency between _aflowlib_entry json and schema";
    aflowlib::_aflowlib_entry aentry;
    string aflowlib_json = aentry.aflowlib2string("JSON", true);
    vector<string> json_keys = aurostd::extractJsonKeysAflow(aflowlib_json);

    vector<string> vkeys_ignore = {"data_language", "error_status", "natoms_orig",
                                   "density_orig", "volume_cell_orig", "volume_atom_orig",
                                   "spinD_magmom_orig"};
    for (const string& key : json_keys) {
      if (!aurostd::WithinList(vkeys_ignore, key) && !aurostd::WithinList(vschema_keys, key)) {
        ninconsistent++;
        if (LDEBUG) std::cerr << function_name << " " << key << " not found in schema." << std::endl;
      }
    }
    checkEqual(ninconsistent, 0, check_function, check_description, passed_checks, results);
  }

}

// xstructure tests
namespace unittest {

  //HE20210511
  void UnitTest::atomicEnvironmentTest(uint& passed_checks, vector<string> & results, vector<string>& errors) {
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
    results.clear();
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
    checkSimilar(AE[test_AE].volume, 31.4622167689, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // Check | hull area
    check_description = "hull area";
    checkSimilar(AE[test_AE].area, 60.4979100628, check_function, check_description, passed_checks, results);

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

  void UnitTest::cifParserTest(uint& passed_checks, vector<string>& results, vector<string>& errors) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string function_name = XPID + "cifParserTest():";
    if (errors.size()) {}  // Suppress compiler warnings

    if (LDEBUG) std::cerr << "Running " << function_name << std::endl;

    // Set up test environment
    string check_function = "", check_description = "", expected_str = "", calculated_str = "";

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
    bool match = false;

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
    match = compare::aflowCompareStructure(xstr_cif, xstr_poscar, same_species, scale_volume, optimize_match, misfit);
    expected_str = "Match";
    calculated_str = string(match?"":"No ") + "Match";
    checkEqual(expected_str, calculated_str, check_function, check_description, passed_checks, results);

    check_description = "Checking parsed Wyckoff positions of CrO3";
    vector<wyckoffsite_ITC> vwyckoff(4);
    xvector<double> coords;
    vwyckoff[0].type = "Cr"; vwyckoff[0].letter = "b"; vwyckoff[0].site_symmetry = "m.."; vwyckoff[0].multiplicity = 4; vwyckoff[0].coord[1] = 0.25; vwyckoff[0].coord[2] = 0.09676; vwyckoff[0].coord[3] = 0.5;
    vwyckoff[1].type = "O"; vwyckoff[1].letter = "a"; vwyckoff[1].site_symmetry = "..2"; vwyckoff[1].multiplicity = 4; vwyckoff[1].coord[1] = 0.0; vwyckoff[1].coord[2] = 0.0; vwyckoff[1].coord[3] = 0.3841;
    vwyckoff[2].type = "O"; vwyckoff[2].letter = "b"; vwyckoff[2].site_symmetry = "m.."; vwyckoff[2].multiplicity = 4; vwyckoff[2].coord[1] = 0.25; vwyckoff[2].coord[2] = 0.2677; vwyckoff[2].coord[3] = 0.3755;
    vwyckoff[3].type = "O"; vwyckoff[3].letter = "b"; vwyckoff[3].site_symmetry = "m.."; vwyckoff[3].multiplicity = 4; vwyckoff[3].coord[1] = 0.25; vwyckoff[3].coord[2] = 0.6078; vwyckoff[3].coord[3] = 0.3284;
    check_passed = (vwyckoff.size() == xstr_cif.wyckoff_sites_ITC.size());
    for (uint i = 0; i < vwyckoff.size() && check_passed; i++) {
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

    expected_str = "Wyckoff positions match";
    calculated_str = "Wyckoff positions " + string(check_passed?"match":"do not match");
    checkEqual(expected_str, calculated_str, check_function, check_description, passed_checks, results);

    // May need a better test case where the labels actually change
    check_description = "Checking calculated Wyckoff positions of CrO3";
    xstr_cif.SpaceGroup_ITC();
    vwyckoff[0].coord[1] = 0.25; vwyckoff[0].coord[2] = 0.59676; vwyckoff[0].coord[3] = 0.0000;
    vwyckoff[1].coord[1] = 0.00; vwyckoff[1].coord[2] = 0.00000; vwyckoff[1].coord[3] = 0.3841;
    vwyckoff[2].coord[1] = 0.25; vwyckoff[2].coord[2] = 0.76770; vwyckoff[2].coord[3] = 0.8755;
    vwyckoff[3].coord[1] = 0.25; vwyckoff[3].coord[2] = 0.60780; vwyckoff[3].coord[3] = 0.3284;
    check_passed = (vwyckoff.size() == xstr_cif.wyckoff_sites_ITC.size());
    for (uint i = 0; i < vwyckoff.size() && check_passed; i++) {
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

    expected_str = "Wyckoff positions match";
    calculated_str = "Wyckoff positions " + string(check_passed?"match":"do not match");
    checkEqual(expected_str, calculated_str, check_function, check_description, passed_checks, results);

    // Test that the CIF parser works for structures with old settings
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

    aurostd::StringstreamClean(xstrss);
    xstrss << str_poscar;
    xstr_poscar = xstructure(xstrss);
    match = compare::aflowCompareStructure(xstr_cif, xstr_poscar, same_species, scale_volume, optimize_match, misfit);
    expected_str = "Match";
    calculated_str = string(match?"":"No ") + "Match";
    checkEqual(expected_str, calculated_str, check_function, check_description, passed_checks, results);
  }

  //CO20190520
  void UnitTest::coordinationTest(uint& passed_checks, vector<string>& results, vector<string>& errors) {
    if (errors.size()) {} // Suppress compiler warnings
    // Set up test environment
    string check_function = "", check_description = "";
    uint expected_uint = 0, calculated_uint = 0;

    xstructure str("aflowlib.duke.edu:AFLOWDATA/ICSD_WEB/FCC/Cl1Na1_ICSD_240599","CONTCAR.relax.vasp",IOAFLOW_AUTO);
    deque<deque<uint> > coordinations;
    str.GetCoordinations(coordinations);

    check_function = "coordinations.size()";
    check_description = "Number of iatoms";
    calculated_uint = coordinations.size();
    expected_uint = 2;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    check_function = "coordinations[0].size()";
    check_description = "Number of coordination environments atom 1";
    calculated_uint = coordinations[0].size();
    expected_uint = 2;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    check_function = "coordinations[1].size()";
    check_description = "Number of coordination environments atom 2";
    calculated_uint = coordinations[1].size();
    expected_uint = 2;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    //first iatom
    //first shell
    check_function = "coordinations[0][0]";
    check_description = "First shell atom 1";
    calculated_uint = coordinations[0][0];
    expected_uint = 6;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    //second shell
    check_function = "coordinations[0][1]";
    check_description = "Second shell atom 1";
    calculated_uint = coordinations[0][1];
    expected_uint = 12;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    //second iatom
    //first shell
    check_function = "coordinations[1][0]";
    check_description = "First shell atom 2";
    calculated_uint = coordinations[1][0];
    expected_uint = 6;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);

    //second shell
    check_function = "coordinations[1][1]";
    check_description = "Second shell atom 2";
    calculated_uint = coordinations[1][1];
    expected_uint = 12;
    checkEqual(calculated_uint, expected_uint, check_function, check_description, passed_checks, results);
  }

  //DX20210129
  void UnitTest::foldAtomsInCellTest(uint& passed_checks, vector<string>& results, vector<string>& errors) {
    // setup test environment
    string check_function = "", check_description = "", expected_str = "", calculated_str = "";
    stringstream message;

    // ---------------------------------------------------------------------------
    // generate rocksalt structure
    check_function = "anrl::getANRLParameters()";
    check_description = "prototype parameters for AB_cF8_225_a_b";
    string prototype_label = "AB_cF8_225_a_b"; 
    vector<string> parameter_sets = anrl::getANRLParameters(prototype_label,"all");
    checkEqual(parameter_sets.size(), 1, check_function, check_description, passed_checks, results);

    check_function = "aflowlib::PrototypeLibraries()";
    check_description = "structure from prototype label and params";
    xstructure xstr;
    bool generated = false;
    try{
      xstr = aflowlib::PrototypeLibraries(*p_oss,prototype_label,parameter_sets[0],1);
      generated = true;
    }
    catch(aurostd::xerror& excpt){
      message << "Could not generate prototype=" << prototype_label << " given parameters=" << parameter_sets[0] << "; check inputs or the symbolic generator.";
      errors.push_back(message.str());
      generated = false;
    }
    expected_str = "Prototype generated";
    calculated_str = "Prototype " + string(generated?"":"not ") + "generated";

    checkEqual(expected_str, calculated_str, check_function, check_description, passed_checks, results);

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
    bool match = compare::structuresMatch(xstr,xstr_supercell,same_species);
    expected_str = "Match";
    calculated_str = "Match", string(match?"":"No ") + "Match";
    checkEqual(expected_str, calculated_str, check_function, check_description, passed_checks, results);

    // ---------------------------------------------------------------------------
    // test 2: reduce cell
    // convert supercell back to original lattice
    check_function = "xstructure::foldAtomsInCell()";
    check_description = "reduce cell into primitive";
    xstructure xstr_reduced = xstr_supercell;
    xstr_reduced.foldAtomsInCell(xstr.lattice, skew, tol, check_min_dists);

    match = compare::structuresMatch(xstr,xstr_reduced,same_species);
    expected_str = "Match";
    calculated_str = "Match", string(match?"":"No ") + "Match";
    checkEqual(expected_str, calculated_str, check_function, check_description, passed_checks, results);
  }

}

bool CeramGenTest(ostream& oss){ofstream FileMESSAGE;return CeramGenTest(FileMESSAGE,oss);}  //CO20190520
bool CeramGenTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy=XPID+"CeramGenTest():";
  //bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  message << "Performing ceramics generation test";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  //./aflow --generate_ceramics --nm=N,C --m=Co,Mo,Fe,Ru,Ni,Rh,Pt,Cu,Cr,V --N=5
  vector<string> vnonmetals,vmetals;
  aurostd::string2tokens("N,C",vnonmetals,",");aurostd::string2tokens("Co,Mo,Fe,Ru,Ni,Rh,Pt,Cu,Cr,V",vmetals,",");

  vector<string> commands=pflow::GENERATE_CERAMICS(vnonmetals,vmetals,5);
  if(commands.size()!=6){
    message << "commands.size()!=6";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N
  //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V
  //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
  //C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
  //C:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V
  //C:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V

  if(!aurostd::WithinList(commands,"C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N")){
    message << "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  if(!aurostd::WithinList(commands,"C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V")){
    message << "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  if(!aurostd::WithinList(commands,"C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V")){
    message << "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  if(!aurostd::WithinList(commands,"C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V")){
    message << "C:Co,Cr,Cu,Fe,Mo:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  if(!aurostd::WithinList(commands,"C:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V")){
    message << "C:Co,Cr,Cu,Fe,Mo:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  if(!aurostd::WithinList(commands,"C:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V")){
    message << "C:N:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V:Ni,Pt,Rh,Ru,V not found";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }

  message << "Ceramics generation test successful";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
  return true;
}
bool EgapTest(ostream& oss){ofstream FileMESSAGE;return EgapTest(FileMESSAGE,oss);}  //CO20190520
bool EgapTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy=XPID+"EgapTest():";
  bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  message << "Performing Egap test";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  string system="",path="",query="",file="",efile="",ext="",tfile="",Egap_type="";
  vector<string> files;
  xOUTCAR xout(FileMESSAGE,oss);
  double EFERMI=AUROSTD_MAX_DOUBLE,Egap=0.0;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  //FCC/Si1_ICSD_150530
  system="ICSD_WEB/FCC/Si1_ICSD_150530";

  path=AFLOWLIB_SERVER_DEFAULT+"/AFLOWDATA/"+system;
  query=path+"/?files";
  message << "Fetching: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  aurostd::url2tokens(query,files,",");
  if(files.size()==0){
    message << "Could not fetch query: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  //OUTCAR.static
  file="OUTCAR.static";
  if(!aurostd::EWithinList(files,file,efile)){
    message << "No " << file << " found within " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  ext=aurostd::GetCompressionExtension(efile);
  tfile=aurostd::TmpFileCreate("Egap_file1")+ext;
  query=path+"/"+efile;
  message << "Fetching: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  if(!(aurostd::url2file(query,tfile,LDEBUG) && aurostd::FileExist(tfile))){
    message << "Could not fetch query: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  message << "Loaded file to: " << tfile;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  if(!xout.GetPropertiesFile(tfile,!LDEBUG)){
    message << "xOUTCAR::GetProperties() failed";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
#ifndef _AFLOW_TEMP_PRESERVE_
  aurostd::RemoveFile(tfile);
#endif
  EFERMI=xout.Efermi;
  //OUTCAR.bands
  file="OUTCAR.bands";
  if(!aurostd::EWithinList(files,file,efile)){
    query=path+"/?files"; //reload query for error message
    message << "No " << file << " found within " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  ext=aurostd::GetCompressionExtension(efile);
  tfile=aurostd::TmpFileCreate("Egap_file1")+ext;
  query=path+"/"+efile;
  message << "Fetching: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  if(!(aurostd::url2file(query,tfile,LDEBUG) && aurostd::FileExist(tfile))){
    message << "Could not fetch query: " << query;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  message << "Loaded file to: " << tfile;pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  if(!xout.GetPropertiesFile(tfile,!LDEBUG)){
    message << "xOUTCAR::GetProperties() failed";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
#ifndef _AFLOW_TEMP_PRESERVE_
  aurostd::RemoveFile(tfile);
#endif
  //GetBandGap
  message << "Running bandgap code";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  if(!xout.GetBandGap(EFERMI)){
    message << "xOUTCAR::GetBandGap() failed";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  Egap=+6.1000e-01;
  if(!aurostd::isequal(Egap,xout.Egap[0])){
    message << "xOUTCAR::GetBandGap() did not find Egap==" << Egap << ", found instead xOUTCAR.Egap[0]==" << xout.Egap[0];pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  Egap_type="insulator-indirect";
  if(xout.Egap_type[0]!="insulator-indirect"){
    message << "xOUTCAR::GetBandGap() did not find type==" << Egap_type << ", found instead xOUTCAR.Egap_type[0]==" << xout.Egap_type[0];pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////

  message << "Egap test successful";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
  return true;
}

bool gcdTest(ostream& oss){ofstream FileMESSAGE;return gcdTest(FileMESSAGE,oss);}  //CO20190520
bool gcdTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy = XPID + "gcdTest():";
  bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  message << "Performing gcd test";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  int a=0,b=0,x1=0,y1=0,gcd=0;

  a=25;b=15;
  aurostd::GCD(a,b,gcd,x1,y1);
  if(!(gcd==5 && x1==-1 && y1==2)){
    if(LDEBUG){
      cerr << soliloquy << " gcd(25,15) failed" << endl;
      cerr << soliloquy << " gcd=" << gcd << endl;
      cerr << soliloquy << " x=" << x1 << endl;
      cerr << soliloquy << " y=" << y1 << endl;
    }
    return FALSE;
  }

  a=25;b=0;
  aurostd::GCD(a,b,gcd,x1,y1);
  if(!(gcd==25 && x1==1 && y1==0)){
    if(LDEBUG){
      cerr << soliloquy << " gcd(25,0) failed" << endl;
      cerr << soliloquy << " gcd=" << gcd << endl;
      cerr << soliloquy << " x=" << x1 << endl;
      cerr << soliloquy << " y=" << y1 << endl;
    }
    return FALSE;
  }

  a=0;b=15;
  aurostd::GCD(a,b,gcd,x1,y1);
  if(!(gcd==15 && x1==0 && y1==1)){
    if(LDEBUG){
      cerr << soliloquy << " gcd(0,15) failed" << endl;
      cerr << soliloquy << " gcd=" << gcd << endl;
      cerr << soliloquy << " x=" << x1 << endl;
      cerr << soliloquy << " y=" << y1 << endl;
    }
    return FALSE;
  }

  a=-5100;b=30450;
  aurostd::GCD(a,b,gcd,x1,y1);
  if(!(gcd==150 && x1==-6 && y1==-1)){
    if(LDEBUG){
      cerr << soliloquy << " gcd(-5100,30450) failed" << endl;
      cerr << soliloquy << " gcd=" << gcd << endl;
      cerr << soliloquy << " x=" << x1 << endl;
      cerr << soliloquy << " y=" << y1 << endl;
    }
    return FALSE;
  }

  message << "gcd test successful";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
  return TRUE; //CO20180419
}

bool smithTest(ostream& oss){ofstream FileMESSAGE;return smithTest(FileMESSAGE,oss);}  //CO20190520
bool smithTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy = XPID + "smithTest():";
  bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  //test ehermite
  xmatrix<int> ehermite(2,2);
  aurostd::getEHermite(5,12,ehermite);
  if(!(
        ehermite[1][1]==5 &&
        ehermite[1][2]==-2 &&
        ehermite[2][1]==-12 &&
        ehermite[2][2]==5 &&
        TRUE
      )
    ){
    if(LDEBUG){cerr << soliloquy << " getEHermite(5,12) failed" << endl;}
    return FALSE;
  }

  xmatrix<int> A1(3,3),U1,V1,S1;
  A1[1][1]=3;A1[1][2]=2;A1[1][3]=1;
  A1[2][1]=5;A1[2][2]=3;A1[2][3]=1;
  A1[3][1]=6;A1[3][2]=8;A1[3][3]=9;

  aurostd::getSmithNormalForm(A1,U1,V1,S1);

  if(LDEBUG){
    cerr << soliloquy << " A=" << endl;cerr << A1 << endl;
    cerr << soliloquy << " U=" << endl;cerr << U1 << endl;
    cerr << soliloquy << " V=" << endl;cerr << V1 << endl;
    cerr << soliloquy << " S=" << endl;cerr << S1 << endl;
  }

  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]if(!(
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      U1[1][1]==24 && U1[1][2]==-13 && U1[1][3]==-1 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      U1[2][1]==13 && U1[2][2]==-7  && U1[2][3]==-1 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      U1[3][1]==2  && U1[3][2]==-1  && U1[3][3]==0  && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      TRUE
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]    )
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  ){
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  if(LDEBUG){cerr << soliloquy << " U1(1) failed of getSmithNormalForm()" << endl;}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  return FALSE;
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]if(!(
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      V1[1][1]==0  && V1[1][2]==1  && V1[1][3]==3  && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      V1[2][1]==-1 && V1[2][2]==-1 && V1[2][3]==-1 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      V1[3][1]==1  && V1[3][2]==0  && V1[3][3]==-1 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      TRUE
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]    )
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  ){
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  if(LDEBUG){cerr << soliloquy << " V1(1) failed of getSmithNormalForm()" << endl;}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  return FALSE;
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]}
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]if(!(
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      S1[1][1]==1 && S1[1][2]==0 && S1[1][3]==0 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      S1[2][1]==0 && S1[2][2]==1 && S1[2][3]==0 && 
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      S1[3][1]==0 && S1[3][2]==0 && S1[3][3]==1 &&
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]      TRUE
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]    )
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  ){
  //[CO20191201 - OBSOLETE: robust check inside getSmithNormalForm()]  if(LDEBUG){cerr << soliloquy << " S1(1) failed of getSmithNormalForm()" << endl;}
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
    cerr << soliloquy << " A=" << endl;cerr << A2 << endl;
    cerr << soliloquy << " U=" << endl;cerr << U2 << endl;
    cerr << soliloquy << " V=" << endl;cerr << V2 << endl;
    cerr << soliloquy << " S=" << endl;cerr << S2 << endl;
  }

  message << "smith test successful";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
  return TRUE; //CO20180419
}

bool PrototypeGeneratorTest(ostream& oss, bool check_symmetry, bool check_uniqueness){ofstream FileMESSAGE;return PrototypeGeneratorTest(FileMESSAGE,oss,check_symmetry,check_uniqueness);} //DX20200925
bool PrototypeGeneratorTest(ofstream& FileMESSAGE,ostream& oss,bool check_symmetry, bool check_uniqueness){  //DX20200925
  string function_name=XPID+"PrototypeGeneratorTest():";
  bool LDEBUG=FALSE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=aurostd::getPWD();

  message << "Testing generation of all AFLOW prototypes" << (check_symmetry?" AND checking symmetry of all generated AFLOW prototypes":check_uniqueness?" AND checking all AFLOW prototypes are unique":"");
  pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  vector<string> prototype_labels, compositions;
  vector<uint> space_group_numbers;
  vector<vector<vector<string> > > grouped_Wyckoff_letters;
  string library = "anrl";

  uint num_protos = aflowlib::GetAllPrototypeLabels(prototype_labels,
      compositions,
      space_group_numbers,
      grouped_Wyckoff_letters,
      library);

  message << "Number of prototype labels = " << num_protos << " (each may have multiple parameter sets)";
  pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);

  string catalog="anrl";
  for(uint i=0;i<num_protos;i++){
    // get parameters
    vector<string> parameter_sets = anrl::getANRLParameters(prototype_labels[i],"all");
    if(LDEBUG){ cerr << "Number of parameters for label=" << prototype_labels[i] << ": " << parameter_sets.size() << endl; }

    for(uint j=0;j<parameter_sets.size();j++){
      xstructure xstr;
      try{
        xstr = aflowlib::PrototypeLibraries(oss,prototype_labels[i],parameter_sets[j],1);
      }
      catch(aurostd::xerror& excpt){
        message << "Could not generate prototype=" << prototype_labels[i] << " given parameters=" << parameter_sets[j] << "; check inputs or the symbolic generator.";
        pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
        return false;
      }

      // check symmetry
      if(check_symmetry){
        if(LDEBUG){ cerr << "Check that the generated structure is consistent with the label=" << prototype_labels[i] << ": " << parameter_sets.size() << endl; }

        // symmetry tolerances
        // some prototype require special tolerance values
        stringstream label_input_ss; label_input_ss << prototype_labels[i] << "-" << std::setw(3) << std::setfill('0') << j+1;
        string label_input = label_input_ss.str();
        double tolerance_sym = anrl::specialCaseSymmetryTolerances(label_input);

        string updated_label_and_params = "";
        if(!anrl::structureAndLabelConsistent(xstr, prototype_labels[i], updated_label_and_params, tolerance_sym)){ //DX20201105 - added symmetry tolerance
          // if changes symmetry, give the appropriate label
          message << "The structure has a higher symmetry than indicated by the label ";
          message << "(orig: proto=" << prototype_labels[i] << " and " << parameter_sets[j] << "). ";
          message << "The correct label and parameters for this structure are:" << endl;
          message << updated_label_and_params << endl;
          message << "Please feed this label and set of parameters into the prototype generator.";
          pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
          return false;
        }
      }
      // check uniqueness
      if(check_uniqueness){
        aurostd::xoption vpflow;
        stringstream label_input_ss; label_input_ss << prototype_labels[i] << "-" << std::setw(3) << std::setfill('0') << j+1;
        string label_input = label_input_ss.str();
        // use special symmetry tolerance if necessary (otherwise, we won't check the prototypes with the correct symmetry)
        double sym_eps = anrl::specialCaseSymmetryTolerances(label_input);
        if(sym_eps!=AUROSTD_MAX_DOUBLE){;
          xstr.sym_eps = sym_eps;
          xstr.sym_eps_calculated = true;
        }
        // check if the prototype matches to more than one prototype
        // (i.e., a prototype should match with itself, but no others)
        vector<string> protos_matching = compare::getMatchingPrototypes(xstr, vpflow, catalog);
        // if it matches to more than one
        if(protos_matching.size()>1 && !anrl::isSpecialCaseEquivalentPrototypes(protos_matching)){
          message << "ERROR: " << prototype_labels[i] << " given parameters=" << parameter_sets[j] << " matches to more than one prototype (and not a documented special case): ";
          message << aurostd::joinWDelimiter(protos_matching,",") << ". ";
          message << "If the prototype was newly added, ONLY include it in the encyclopedia for a valid reason (e.g., historical, special designation, etc.)";
          message << " and document this in anrl::isSpecialCaseEquivalentPrototypes().";
          pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
          return false;
        }
        // if it doesn't match with ITSELF
        if(protos_matching.size()==0){
          message << "ERROR: " << prototype_labels[i] << " given parameters=" << parameter_sets[j] << " does NOT match to any prototypes ";
          message << "(either this system requires a special symmetry tolerance or there is a bug with XtalFinder)." << endl;
          pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
        }
      }
    }
  }
  message << "Successfully generated all prototypes!";
  pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);

  return true;
}

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
