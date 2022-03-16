// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************

#include "aflow.h"
#include "aflow_anrl.h"  //DX20201104
#include "aflow_compare_structure.h"  //ME20220125
#include "aflow_apdc.h"  //SD20220202


// Collection of generic check functions, to streamline testing.
// HE20210616
template <typename utype>
void check(const bool &passed, const utype &calculated, const utype &expected, const string &check_function,
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
void check_equal(const utype &calculated, const utype &expected, const string &check_function,
    const string check_description, uint &passed_checks, vector<string> &results){
  bool passed = false;
  if (aurostd::isequal(calculated, expected)) passed = true;
  check(passed, calculated, expected, check_function, check_description, passed_checks, results);
}
void check_equal(const string &calculated, const string &expected, const string &check_function, 
    const string check_description, uint &passed_checks, vector<string> &results){
  bool passed = false;
  if (calculated == expected) passed = true;
  check(passed, calculated, expected, check_function, check_description, passed_checks, results);
}
void check_equal(const bool &calculated, const bool &expected, const string &check_function,
    const string check_description, uint &passed_checks, vector<string> &results){
  bool passed = false;
  if (calculated == expected) passed = true;
  check(passed, calculated, expected, check_function, check_description, passed_checks, results);
}

void check_similar(const double &calculated, const double &expected, const string &check_function,
    const string check_description, uint &passed_checks, vector<string> &results, const double &relative=1E-10){
  bool passed = false;
  if (std::abs(expected - calculated) <= expected * relative) passed = true;
  check(passed, calculated, expected, check_function, check_description, passed_checks, results);
}

bool display_result(const uint passed_checks, const string & task_description, const vector<string> & results,
    const string & function_name, ofstream & FileMESSAGE, ostream & oss){
  stringstream message;
  _aflags aflags;
  aflags.Directory=aurostd::getPWD();
  uint check_num = results.size();
  if (passed_checks == check_num) {
    message << "SUCCESS " << task_description << " (passing " << check_num << " checks)" << endl;
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
    message << "\t" << aurostd::joinWDelimiter(results, "\n\t");
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);
  }
  else {
    message << "FAIL " << task_description << " (" << check_num-passed_checks << " of " << check_num << " checks failed)" << endl;
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    message << "\t" << aurostd::joinWDelimiter(results, "\n\t");
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_RAW_);
    return false;
  }
  return true;
}

// This should become a collection of tests regarding aurostd.
// At the moment, just the functions aurostd::volume and aurostd::area are tested here.
bool aurostdTest(ostream& oss){ofstream FileMESSAGE; return aurostdTest(FileMESSAGE,oss);} //HE20210511
bool aurostdTest(ofstream& FileMESSAGE, ostream& oss) { //HE20210511

  string function_name = XPID + "aurostdTest():";
  stringstream message;


  // setup test environment
  string task_description = "Testing aurostd";
  vector<string> results;
  stringstream result;
  uint passed_checks = 0;
  string check_function = "";
  string check_description = "";

  double expected = 0.0;
  double calculated = 0.0;

  int expected_int = 0;
  string expected_error = "";
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
  expected = 2.5;

  calculated = aurostd::volume(points, facets, true);
  check_equal(calculated, expected, check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | convex solid volume (int)
  check_function = "aurostd::volume()";
  check_description = "convex solid, points as int";

  calculated = aurostd::volume(ipoints, facets, true);
  check_equal(calculated, expected, check_function, check_description, passed_checks, results);


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
  expected = 40.0;

  calculated = aurostd::volume(points, facets);
  check_equal(calculated, expected, check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | error facet/normals mismatch
  check_function = "aurostd::volume()";
  check_description = "error: facet/normals mismatch";
  vector<xvector<double> > normals;
  expected_error = "xerror code 30 (VALUE_ERROR)";
  expected_int = _VALUE_ERROR_;

  try {
    calculated = aurostd::volume(points, facets, normals);
    check(false, std::string("no error"), expected_error, check_function, check_description, passed_checks, results);
  }
  catch (aurostd::xerror e)
  {
    if (e.error_code == expected_int) check(true, "", "", check_function, check_description, passed_checks, results);
    else check(false, aurostd::utype2string(e.error_code), expected_error, check_function, check_description, passed_checks, results);
  }
  catch (...) {
    check(false, std::string("not an xerror"), expected_error, check_function, check_description, passed_checks, results);
  }

  // ---------------------------------------------------------------------------
  // Check | non convex solid volume (int)
  check_function = "aurostd::volume()";
  check_description = "non convex solid, points as int";
  expected = 40.0;

  calculated = aurostd::volume(ipoints, facets);
  check_equal(calculated, expected, check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | error facet size
  check_function = "aurostd::volume()";
  check_description = "error: wrong facet size";
  expected_error = "xerror code 30 (VALUE_ERROR)";
  expected_int = _VALUE_ERROR_;

  facet.clear(); facet.resize(2); facet[0]=1; facet[1]=2; facets.push_back(facet);
  try {
    calculated = aurostd::volume(points, facets);
    check(false, std::string("no error"), expected_error, check_function, check_description, passed_checks, results);
  }
  catch (aurostd::xerror e)
  {
    if (e.error_code == expected_int) check(true, "", "", check_function, check_description, passed_checks, results);
    else check(false, aurostd::utype2string(e.error_code), expected_error, check_function, check_description, passed_checks, results);
  }
  catch (...) {
    check(false, std::string("not an xerror"), expected_error, check_function, check_description, passed_checks, results);
  }

  // ---------------------------------------------------------------------------
  // Check | non convex area (double)
  check_function = "aurostd::areaPointsOnPlane()";
  check_description = "non convex area; points as double";
  expected = 10.0;

  //fill vectors with data
  points.clear(); ipoints.clear(); facets.clear();
  points.push_back(p0); points.push_back(p1); points.push_back(p2); points.push_back(p3); points.push_back(p4);
  points.push_back(p5);
  ipoints.push_back(p0i); ipoints.push_back(p1i); ipoints.push_back(p2i); ipoints.push_back(p3i); ipoints.push_back(p4i);
  ipoints.push_back(p5i);

  calculated = aurostd::areaPointsOnPlane(points);
  check_equal(calculated, expected, check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | non convex area (int)
  check_function = "aurostd::areaPointsOnPlane()";
  check_description = "non convex area; points as int";
  expected = 10.0;

  calculated = aurostd::areaPointsOnPlane(ipoints);
  check_equal(calculated, expected, check_function, check_description, passed_checks, results);


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
  expected = 3.5355339059;

  calculated = aurostd::areaPointsOnPlane(points);
  check_similar(calculated, expected, check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | 3d triangle area (int)
  check_function = "aurostd::areaPointsOnPlane()";
  check_description = "3d triangle; points as int";
  expected = 3.5355339059;

  calculated = aurostd::areaPointsOnPlane(ipoints);
  check_similar(calculated, expected, check_function, check_description, passed_checks, results);

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

  check_equal(result_ss.str(), answer, check_function, check_description, passed_checks, results);

  // present overall result
  return display_result(passed_checks, task_description, results, function_name, FileMESSAGE, oss);
}

bool AtomicEnvironmentTest(ostream& oss){ofstream FileMESSAGE;return AtomicEnvironmentTest(FileMESSAGE,oss);} //HE20210511
bool AtomicEnvironmentTest(ofstream& FileMESSAGE, ostream& oss){ //HE20210511

  string function_name = XPID + "AtomicEnvironmentTest():";

  // setup test environment
  string task_description = "Creating a atomic environment [aflow:d912e209c81aeb94]";
  vector<string> results;
  stringstream result;
  uint passed_checks = 0;
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

  // ---------------------------------------------------------------------------
  // Check | number of created AEs
  check_description = "number of created AEs";
  check_equal(uint(AE.size()), uint(6), check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | point index mapping
  check_description = "point index mapping";
  check_equal(AE[1].index2Point(10), AE[1].coordinates_neighbor[1][2], check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | coordinate matching
  check_description = "coordinate matching";
  xvector<double> compare_point(3,1);
  compare_point(1) = -1.9551925593108925e0;
  compare_point(2) = -2.2642136090979212e0;
  compare_point(3) = 2.4896636484942385e0;

  check_equal(AE[1].index2Point(2), compare_point, check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | center element
  check_description = "center element";
  check_equal(std::string(AE[2].element_center), std::string("Ca"), check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | center element id
  check_description = "center element id";
  check_equal(AE[4].type_center, uint(1), check_function, check_description, passed_checks, results);

  // present overall result
  bool isPassed = display_result(passed_checks, task_description, results, function_name, FileMESSAGE, oss);
  if (!isPassed) return isPassed;
  // ---------------------------------------------------------------------------
  // Test 2: create AE convex hull
  // ---------------------------------------------------------------------------

  // setup test environment
  results.clear();
  result.str("");
  result.clear();
  task_description = "Creating convex hull with constructAtomEnvironmentHull() [aflow:d912e209c81aeb94]";
  passed_checks = 0;
  const uint test_AE = 4;

  // create hull
  AE[test_AE].constructAtomEnvironmentHull();

  // ---------------------------------------------------------------------------
  // Check | hull bit set
  check_description = "hull bit set";
  check_equal(AE[test_AE].has_hull, true, check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | hull volume
  check_description = "hull volume";
  check_similar(AE[test_AE].volume, 31.4622167689, check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | hull area
  check_description = "hull area";
  check_similar(AE[test_AE].area, 60.4979100628, check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | triangle count
  check_description = "triangle count";
  check_equal(AE[test_AE].facet_order[0], uint(6), check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | tetragon count
  check_description = "tetragon count";
  check_equal(AE[test_AE].facet_order[1], uint(2), check_function, check_description, passed_checks, results);

  // ---------------------------------------------------------------------------
  // Check | pentagon count
  check_description = "pentagon count";
  check_equal(AE[test_AE].facet_order[2], uint(0), check_function, check_description, passed_checks, results);

  // present overall result
  return display_result(passed_checks, task_description, results, function_name, FileMESSAGE, oss);
}

bool SchemaTest(ostream& oss){ofstream FileMESSAGE;return SchemaTest(FileMESSAGE,oss);}
bool SchemaTest(ofstream& FileMESSAGE,ostream& oss) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string function_name = XPID + "SchemaTest():";
  _aflags aflags; aflags.Directory = aurostd::getPWD();

  // Set up test environment
  string task_description = "Testing schema";
  vector<string> results;
  uint passed_checks = 0;
  string check_function = "";
  string check_description = "";

  check_function = "XHOST.vschema";
  check_description = "Internal consistency of vschema";
  vector<string> vschema_keys;
  string schema_types = "UNIT,TYPE";
  vector<string> vschema_types;
  aurostd::string2tokens(schema_types, vschema_types, ",");
  string key = "";
  uint ninconsistent = 0;
  for (uint i = 0; i < XHOST.vschema.vxsghost.size(); i+= 2) {
    if(XHOST.vschema.vxsghost[i].find("::NAME:") != string::npos) {
      key=aurostd::RemoveSubString(XHOST.vschema.vxsghost[i], "SCHEMA::NAME:");
      vschema_keys.push_back(XHOST.vschema.getattachedscheme("SCHEMA::NAME:" + key));
      for (uint j = 0; j < vschema_types.size(); j++) {
        if (!XHOST.vschema.isdefined("SCHEMA::" + vschema_types[j] + ":" + key)) {
          ninconsistent++;
          if (LDEBUG) std::cerr << "SCHEMA::" << vschema_types[j] << ":" << key << " not found." << std::endl;
        }
      }
    }
  }
  check_equal(ninconsistent, 0, check_function, check_description, passed_checks, results);

  check_function = "_aflowlib_entry";
  check_description = "Consistency between _aflowlib_entry json and schema";
  aflowlib::_aflowlib_entry aentry;
  string aflowlib_json = aentry.aflowlib2string("JSON", true);
  vector<string> json_keys = aurostd::extractJsonKeysAflow(aflowlib_json);

  string keys_ignore = "data_language,error_status,natoms_orig,density_orig,volume_cell_orig,volume_atom_orig,spinD_magmom_orig";
  vector<string> vkeys_ignore;
  aurostd::string2tokens(keys_ignore, vkeys_ignore, ",");
  for (uint i = 0; i < json_keys.size(); i++) {
    const string& key = json_keys[i];
    if (!aurostd::WithinList(vkeys_ignore, key) && !aurostd::WithinList(vschema_keys, key)) {
      ninconsistent++;
      if (LDEBUG) std::cerr << key << " not found in schema." << std::endl;
    }
  }
  check_equal(ninconsistent, 0, check_function, check_description, passed_checks, results);

  return display_result(passed_checks, task_description, results, function_name, FileMESSAGE, oss);
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

bool coordinationTest(ostream& oss){ofstream FileMESSAGE;return coordinationTest(FileMESSAGE,oss);}  //CO20190520
bool coordinationTest(ofstream& FileMESSAGE,ostream& oss){  //CO20190520
  string soliloquy=XPID+"coordinationTest():";
  bool LDEBUG=TRUE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=".";

  xstructure str("aflowlib.duke.edu:AFLOWDATA/ICSD_WEB/FCC/Cl1Na1_ICSD_240599","CONTCAR.relax.vasp",IOAFLOW_AUTO);
  deque<deque<uint> > coordinations;
  str.GetCoordinations(coordinations);
  if(coordinations.size()<2){
    if(LDEBUG){cerr << soliloquy << " coordinations not found" << endl;}
    return FALSE;
  }
  if(coordinations[0].size()<2){
    if(LDEBUG){cerr << soliloquy << " coordinations[0] not found" << endl;}
    return FALSE;
  }
  if(coordinations[1].size()<2){
    if(LDEBUG){cerr << soliloquy << " coordinations[1] not found" << endl;}
    return FALSE;
  }
  //first iatom
  //first shell
  if(coordinations[0][0]!=6){
    if(LDEBUG){cerr << soliloquy << " coordinations[0][0]!=6 (==" << coordinations[0][0] << ")" << endl;}
    return FALSE;
  }
  //second shell
  if(coordinations[0][1]!=12){
    if(LDEBUG){cerr << soliloquy << " coordinations[0][1]!=12 (==" << coordinations[0][1] << ")" << endl;}
    return FALSE;
  }
  //second iatom
  //first shell
  if(coordinations[1][0]!=6){
    if(LDEBUG){cerr << soliloquy << " coordinations[1][0]!=6 (==" << coordinations[1][0] << ")" << endl;}
    return FALSE;
  }
  //second shell
  if(coordinations[1][1]!=12){
    if(LDEBUG){cerr << soliloquy << " coordinations[1][1]!=12 (==" << coordinations[1][1] << ")" << endl;}
    return FALSE;
  }

  message << "coordination test successful";pflow::logger(_AFLOW_FILE_NAME_,soliloquy,message,aflags,FileMESSAGE,oss,_LOGGER_COMPLETE_);
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

bool FoldAtomsInCellTest(ostream& oss){ofstream FileMESSAGE;return FoldAtomsInCellTest(FileMESSAGE,oss);} //DX20210129
bool FoldAtomsInCellTest(ofstream& FileMESSAGE,ostream& oss){ //DX20210129
  string function_name=XPID+"FoldAtomsInCellTest():";
  //bool LDEBUG=FALSE; // TRUE;
  stringstream message;
  _aflags aflags;aflags.Directory=aurostd::getPWD();

  // ---------------------------------------------------------------------------
  // generate rocksalt structure
  string prototype_label = "AB_cF8_225_a_b"; 
  vector<string> parameter_sets = anrl::getANRLParameters(prototype_label,"all");

  if(parameter_sets.size() != 1){
    message << "Expected only one parameter set for the rocksalt structure (" << prototype_label << ") # different parameter sets=" << parameter_sets.size() << ".";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }

  xstructure xstr;
  try{
    xstr = aflowlib::PrototypeLibraries(oss,prototype_label,parameter_sets[0],1);
  }
  catch(aurostd::xerror& excpt){
    message << "Could not generate prototype=" << prototype_label << " given parameters=" << parameter_sets[0] << "; check inputs or the symbolic generator.";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_ERROR_);
    return false;
  }

  // set fold atoms in cell variables
  bool skew = false;
  double tol = 0.01;
  bool check_min_dists = false;

  // ---------------------------------------------------------------------------
  // test 1: expand cell
  // create 3x1x1 supercell expansion matrix
  xmatrix<double> supercell_matrix = aurostd::eye<double>(3,3);
  supercell_matrix(1,1)=3.0;
  xmatrix<double> lattice_new = supercell_matrix*xstr.lattice; // transform lattice

  xstructure xstr_supercell = xstr;
  xstr_supercell.foldAtomsInCell(lattice_new, skew, tol, check_min_dists);

  bool same_species = true;
  if(compare::structuresMatch(xstr,xstr_supercell,same_species)){
    message << "Successfully expanded rocksalt structure into a 3x1x1 supercell via foldAtomsInCell().";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  }
  else{
    message << "Expanded rocksalt structure is not equivalent to the original structure; bad expansion.";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);
    return false;
  }

  // ---------------------------------------------------------------------------
  // test 2: reduce cell
  // convert supercell back to original lattice
  xstructure xstr_reduced = xstr_supercell;
  xstr_reduced.foldAtomsInCell(xstr.lattice, skew, tol, check_min_dists);

  if(compare::structuresMatch(xstr,xstr_reduced,same_species)){
    message << "Successfully reduced 3x1x1 rocksalt structure into a primitive form via foldAtomsInCell().";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  }
  else{
    message << "Reduced 3x1x1 rocksalt structure is not equivalent to the original structure; bad reduction.";
    pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_WARNING_);
    return false;
  }

  return true;
}

//ME20220125
bool cifParserTest(ostream& oss){ofstream FileMESSAGE;return cifParserTest(FileMESSAGE,oss);}
bool cifParserTest(ofstream& FileMESSAGE, ostream& oss) {
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  string function_name = XPID + "cifParserTest():";

  // Set up test environment
  string task_description = "Testing CIF parser";
  vector<string> results;
  uint passed_checks = 0;
  string check_function = "";
  string check_description = "";
  string calculated = "";
  string expected = "";

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
  expected = "Match";
  calculated = string((match?"Match":"No Match")) +  " (misfit = " + aurostd::utype2string<double>(misfit) + ")";
  check(match && (misfit < _ZERO_TOL_), calculated, expected, check_function, check_description, passed_checks, results);

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

  expected = "Wyckoff positions match";
  calculated = "Wyckoff positions " + string(check_passed?"matched":"did not match") + ".";
  check(check_passed, calculated, expected, check_function, check_description, passed_checks, results);

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

  expected = "Wyckoff positions match";
  calculated = "Wyckoff positions " + string(check_passed?"matched":"did not match") + ".";
  check(check_passed, calculated, expected, check_function, check_description, passed_checks, results);

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
  expected = "Match";
  calculated = string(match?"Match":"No Match") +  " (misfit = " + aurostd::utype2string<double>(misfit) + ")";
  check(match && (misfit < _ZERO_TOL_), calculated, expected, check_function, check_description, passed_checks, results);

  return display_result(passed_checks, task_description, results, function_name, FileMESSAGE, oss);
}

//SD20220202
bool APDCTest(ostream& oss){ofstream FileMESSAGE;return APDCTest(FileMESSAGE,oss);}
bool APDCTest(ofstream& FileMESSAGE, ostream& oss) {
  _apdc_data apdc_data;
  // logger
  string function_name = XPID + "APDCTest():";
  stringstream message;
  _aflags aflags;
  message << "Testing APDC";
  // initalize data
  apdc_data.rootdirpath = "/home/sd453/tmp/testing2/";
  apdc_data.plattice = "FCC";
  apdc_data.elements = vector<string>(2);
  apdc_data.elements[0] = "Pt";
  apdc_data.elements[1] = "Au";
  // run functions
  apdc::GetPhaseDiagram(apdc_data);

  // return
  pflow::logger(_AFLOW_FILE_NAME_,function_name,message,aflags,FileMESSAGE,oss,_LOGGER_MESSAGE_);
  return TRUE;
}
