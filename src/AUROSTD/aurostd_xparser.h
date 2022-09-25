// ***************************************************************************
// *                                                                         *
// *              Aflow COREY OSES - Duke University 2003-2021               *
// *                                                                         *
// ***************************************************************************
// Written by Corey Oses 2020

#ifndef _AUROSTD_XPARSER_H_
#define _AUROSTD_XPARSER_H_

#include "aurostd.h"

//compound specification is how a compound is specified
//composition (Mn2Pt3) is ORTHOGONAL to pseudopotential string (Mn_pvPt)
//for instance, H1.25 can be a pseudopotential and NOT a composition
enum elements_string_type {
  composition_string,
  pp_string,
};

//CO20190712 - see VASP_PseudoPotential_CleanName_InPlace() in aflow_ivasp.cpp
const string CAPITAL_LETTERS_PP_LIST="_GW2"    //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Li_AE_GW2
",_GW"    //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/As_GW
",_ZORA"  //CO20190712 - potpaw_PBE/potpaw_PBE.20100506/Pt_ZORA
",_LDApU" //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Zn_sv_LDApU
",_AE"    //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/Li_AE_GW2
",_NC2"   //CO20190712 - potpaw_LDA/potpaw_LDA.20100505/As_NC2
",_200eV"
"";

namespace aurostd {
  void VASP_PseudoPotential_CleanName_InPlace(string& species,bool capital_letters_only=false,bool remove_floats=true); //CO20190712  //CO20210623 - added remove_floats
  ////////////////////////////////////////////////////////////////////////////////
  void elementsFromCompositionString(const string& input);  //CO20190712
  template<class utype> void elementsFromCompositionString(const string& input,vector<string>& velements,vector<utype>& vcomposition); //CO20190712
  void elementsFromPPString(const string& input,vector<string>& velements,bool keep_pp=false); //CO20190712
  ////////////////////////////////////////////////////////////////////////////////
  // returns UNSORTED vector<string> from string
  vector<string> getElements(const string& input); //CO20190712
  vector<string> getElements(const string& input,elements_string_type e_str_type,bool clean=true,bool sort_elements=false,bool keep_pp=false,ostream& oss=cout);
  template<class utype> vector<string> getElements(const string& input,vector<utype>& vcomposition,bool clean=true,bool sort_elements=false,bool keep_pp=false,ostream& oss=cout);
  vector<string> getElements(const string& input,elements_string_type e_str_type,ofstream& FileMESSAGE,bool clean=true,bool sort_elements=false,bool keep_pp=false,ostream& oss=cout);
  template<class utype> vector<string> getElements(const string& input,vector<utype>& vcomposition,elements_string_type e_str_type,bool clean=true,bool sort_elements=false,bool keep_pp=false,ostream& oss=cout);
  template<class utype> vector<string> getElements(const string& input,vector<utype>& vcomposition,elements_string_type e_str_type,ofstream& FileMESSAGE,bool clean=true,bool sort_elements=false,bool keep_pp=false,ostream& oss=cout);
} // namespace aurostd

//AS20201214 BEGIN JSONwriter
namespace aurostd {
  /// Class-container to output data into JSON format.
  class JSONwriter{
    private:
      vector<string> content;
      void free();
      void copy(const JSONwriter &jw);
    public:
      JSONwriter();
      JSONwriter(const JSONwriter &jw);
      ~JSONwriter();
      const JSONwriter& operator=(const JSONwriter &jw);
      void clear();
      template <typename utype> void addNumber(const string &key, const utype value);
      template <typename utype> void addVector(const string &key, const utype &value);
      void addVector(const string &key, const vector<double> &value,
          int precision = AUROSTD_DEFAULT_PRECISION, bool roundoff = false,
          double tol = AUROSTD_ROUNDOFF_TOL);
      void addVector(const string &key, const deque<double> &value,
          int precision = AUROSTD_DEFAULT_PRECISION, bool roundoff = false,
          double tol = AUROSTD_ROUNDOFF_TOL);
      void addVector(const string &key, const xvector<double> &value,
          int precision = AUROSTD_DEFAULT_PRECISION, bool roundoff = false,
          double tol = AUROSTD_ROUNDOFF_TOL); //DX20210308 - added xvector variant
      void addVector(const string &key, const vector<string> &value, bool wrap=true); //DX20210301 - added wrap
      void addVector(const string &key, const deque<string> &value, bool wrap=true); //DX20210301 - added wrap
      void addVector(const string &key, vector<JSONwriter> &value); //AS20210309
      template <typename utype> void addMatrix(const string &key, const utype &value);
      void addMatrix(const string &key, const vector<vector<double> > &value,
          int precision = AUROSTD_DEFAULT_PRECISION, bool roundoff = false,
          double tol = AUROSTD_ROUNDOFF_TOL);
      void addMatrix(const string &key, const deque<deque<double> > &value,
          int precision = AUROSTD_DEFAULT_PRECISION, bool roundoff = false,
          double tol = AUROSTD_ROUNDOFF_TOL);
      void addMatrix(const string &key, const xmatrix<double> &value,
          int precision = AUROSTD_DEFAULT_PRECISION, bool roundoff = false,
          double tol = AUROSTD_ROUNDOFF_TOL); //DX20210308 - added xmatrix variant
      void addMatrix(const string &key, const vector<vector<string> > &value); //DX20210211
      void addMatrix(const string &key, const deque<deque<string> > &value); //DX20210211
      void addString(const string &key, const string &value);
      void addBool(const string &key, bool value);
      void mergeRawJSON(const string &value); //DX20210304 - changed from addRaw to mergeRawJSON
      void addRawJSON(const string &key, const string& value); //DX20210301
      void addNull(const string &key); //DX20210301
      void addJSON(const string &key, JSONwriter &value);
      void mergeJSON(JSONwriter &value);
      string toString(bool wrap=true);
  };
}
//AS20201214 END

//AS20210309 - added JSON benchmarks (keeping to easily test later, but commented out for now)
//AS20210309 [BENCHMARK] namespace aurostd {
//AS20210309 [BENCHMARK] 	 void test_vector_string(int niterations);
//AS20210309 [BENCHMARK]	 void test_vector_json(int niterations);
//AS20210309 [BENCHMARK]   void run_vector_string_vs_json_test(void);
//AS20210309 [BENCHMARK] }

//ME2020408 - JSON reader
//Moved from the AflowDB class

namespace aurostd {
  class JSONReader{
  public:
    enum class object_types{ // order is important!
      DICTIONARY,
      LIST,
      STRING,
      FLOAT,
      INTEGER,
      T,
      F,
      NONE
    };
    struct storage_object {
      object_types type = object_types::NONE;
      std::shared_ptr<void> obj;
      JSONReader::storage_object &operator[](const size_t index);
      JSONReader::storage_object &operator[](const std::string key);
      JSONReader::storage_object &operator[](const char* key);

      explicit operator bool() const;
      explicit operator std::string() const;
      explicit operator double() const;
      explicit operator float() const;
      operator long long() const;
      explicit operator unsigned long long () const;
      explicit operator unsigned long () const;
      template<class utype> operator std::vector<utype>() const;
      template<class utype> operator std::map<std::string, utype>() const;
      template<class utype> operator aurostd::xvector<utype>() const;
      template<class utype> operator aurostd::xmatrix<utype>() const;
      std::string toString(const bool json_format=true) const;
    };
  private:

    void copy(const JSONReader &jsonReader);

    typedef std::map<std::string, storage_object> Dictionary;
    typedef std::vector<storage_object> List;

    storage_object root;

    std::pair<size_t, size_t> find_string(const std::string & raw_content, std::pair<size_t, size_t> border={0,0}) const;
    std::pair<size_t, size_t> find_bracket(const std::string & raw_content, char kind_open, std::pair<size_t, size_t> border={0,0}) const;
    std::pair<size_t, size_t> find_strip(const std::string & raw_content, std::pair<size_t, size_t> border={0,0}) const;

    storage_object parse(const std::string &raw_content, std::pair<size_t, size_t> border={0,0}) const;
    std::string parse_string(const std::string & raw_content, std::pair<size_t, size_t> border={0,0}) const;
    std::string utf8(const std::string & cp) const;

    friend struct storage_object;
  public:
    JSONReader::storage_object &operator[](const size_t index);
    JSONReader::storage_object &operator[](const std::string key);
    JSONReader::storage_object &operator[](const char *key);
    void loadFile(const std::string & file_path);
    void loadString(const std::string & content);
    std::string toString() const;
    JSONReader();
    JSONReader(const JSONReader &jr);
    ~JSONReader();
    explicit operator storage_object() const;
    const JSONReader& operator=(const JSONReader &jr);
    void clear();

  };
  ostream& operator<<(ostream& os, const JSONReader::storage_object& so);
  ostream& operator<<(ostream& os, const JSONReader& jr);
}

namespace aurostd {
  vector<string> extractJsonKeysAflow(const string& json);
  string extractJsonValueAflow(const string& json, string key);
}

#endif // _AUROSTD_XPARSER_H_

// **************************************************************************
// *                                                                        *
// *              Aflow COREY OSES - Duke University 2003-2021              *
// *                                                                        *
// **************************************************************************
