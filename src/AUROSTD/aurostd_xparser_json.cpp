// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
// This JSON class is the evolution of different prior solutions to integrate JSON with the AFLOW source base.
// hagen.eckert@duke.edu


#ifndef _AUROSTD_XPARSER_JSON_CPP_
#define _AUROSTD_XPARSER_JSON_CPP_

#ifndef _AUROSTD_XPARSER_H_

#include "aurostd_xparser.h"

#endif


namespace aurostd {

  /// @struct JSON::object
  /// @brief storge container for a JSON object
  ///
  /// @authors
  /// @mod{HE,20220924,created struct}

  /// @brief direct index access to JSON::object_types::LIST objects
  /// @param index list index
  /// @return JSON::object stored at index
  /// @note throws an error if used on something other than JSON::object_types::LIST
  /// @note the returned JSON::object can be directly be saved into an appropriate variable
  /// like `string content = json_obj[3];`
  ///
  /// @authors
  /// @mod{HE,20220924,created}
  JSON::object &JSON::object::operator[](const size_t index) const{
    if (this->type == JSON::object_types::LIST) {
      std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(this->obj);
      return content->operator[](index);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Not a list", _INDEX_ILLEGAL_);
    }
  }

  /// @brief direct key access to JSON::object_types::DICTIONARY objects
  /// @param key dictionary key
  /// @return JSON::object stored at key
  /// @note throws an error if used on something other than JSON::object_types::DICTIONARY
  /// @note the returned JSON::object can be directly be saved into an appropriate variable
  /// like `string content = json_obj['my_key'];`
  ///
  /// @authors
  /// @mod{HE,20220924,created}
  JSON::object &JSON::object::operator[](const std::string key) const{
    if (this->type == JSON::object_types::DICTIONARY) {
      std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(this->obj);
      return content->operator[](key);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Not a dictionary", _INDEX_ILLEGAL_);
    }
  }

  /// @brief direct key access to JSON::object_types::DICTIONARY objects
  JSON::object &JSON::object::operator[](const char *key) const {
    if (this->type == JSON::object_types::DICTIONARY) {
      std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(this->obj);
      return content->operator[](key);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Not a dictionary", _INDEX_ILLEGAL_);
    }
  }

  /// @brief enables the direct use of JSON::object with `cout`
  ostream &operator<<(ostream &os, const JSON::object &so) {
    os << so.toString();
    return os;
  }

  /// @brief converts a JSON::object into a string
  /// @param json_format if `true` add encapsulate strings in `"` (default: `true`)
  /// @param escape if `true` escape unicode in strings - `αβγ` becomes `\u03b1\u03b2\u03b3` (default: `true`)
  /// @return JSON formatted string
  ///
  /// @authors
  /// @mod{HE,20220924,created}
  std::string JSON::object::toString(const bool json_format, const bool escape) const {
    bool first = true;
    stringstream result;
    switch (type) {
      case object_types::DICTIONARY: {
        result << "{";
        std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
        for (const auto &entry: *content) {
          if (first) first = false;
          else result << ",";
          result << "\"" << entry.first << "\":" << entry.second; //
        }
        result << "}";
      }
        break;
      case object_types::LIST: {
        result << "[";
        std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
        for (const auto &entry: *content) {
          if (first) first = false;
          else result << ",";
          result << entry; //
        }
        result << "]";
      }
        break;
      case object_types::STRING: {
        std::shared_ptr <std::string> content = std::static_pointer_cast<std::string>(obj);
        if (json_format && escape) result << "\"" << JSON::escape_unicode(*content) << "\"";
        else if (json_format) result << "\"" << *content << "\"";
        else return *content;
      }
        break;
      case object_types::INTEGER: {
        std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
        result << *content;
      }
        break;
      case object_types::FLOAT: {
        std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
        result << *content;
      }
        break;
      case object_types::T: {
        result << "true";
      }
        break;
      case object_types::F: {
        result << "false";
      }
        break;
      case object_types::NONE: {
        result << "null";
      }
        break;
      default:
        break;
    }
    return result.str();
  }


  /// @brief change this JSON::object to a JSON::object_types::STRING
  /// @authors
  /// @mod{HE,20220926,created}
  void JSON::object::fromString(const std::string & content){
    std::shared_ptr <std::string> new_string = std::make_shared<std::string>();
    *new_string = content;
    this->obj = new_string;
    this->type = object_types::STRING;
  }

  /// @brief change this JSON::object to a JSON::object_types::INTEGER
  /// @authors
  /// @mod{HE,20220926,created}
  template<class utype>  void JSON::object::fromNumber(const utype content){
    if (std::is_integral<utype>::value){
      this->type = object_types::INTEGER;
      std::shared_ptr <long long> new_content = std::make_shared<long long>();
      *new_content = (long long) content;
      this->obj = new_content;
      return;
    } else if (std::is_floating_point<utype>::value){
      if (std::isnan(content)) {
        this->type = object_types::NONE;
        this->obj = nullptr;
        return;
      }
      this->type = object_types::FLOAT;
      std::shared_ptr <double> new_content = std::make_shared<double>();
      *new_content = (double) content;
      this->obj = new_content;
      return;
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on a xvector
  /// @authors
  /// @mod{HE,20220926,created}
  template<class utype> void JSON::object::fromXvector(const xvector<utype> &content) {
    std::shared_ptr <JSON::List> new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = object_types::LIST;
    for (int idx=content.lrows; idx<=content.urows; idx++) {
      new_list->push_back(content[idx]);
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on a xmatrix
  /// @authors
  /// @mod{HE,20220926,created}
  template<class utype> void JSON::object::fromXmatrix(const xmatrix<utype> &content) {
    std::shared_ptr <JSON::List> new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = object_types::LIST;
    for (int idx=content.lrows; idx<=content.urows; idx++) {
      new_list->push_back(content(idx));
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on a vector
  /// @authors
  /// @mod{HE,20220926,created}
  template<class utype> void JSON::object::fromVector(const vector<utype> & content) {
    std::shared_ptr <JSON::List> new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = object_types::LIST;
    for (utype entry: content) {
      new_list->push_back(entry);
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::DICTIONARY based on a map
  /// @authors
  /// @mod{HE,20220926,created}
  template<class utype> void JSON::object::fromMap(const std::map<std::string, utype> & content) {
    std::shared_ptr <JSON::Dictionary> new_dictionary = std::make_shared<JSON::Dictionary>();
    this->obj = new_dictionary;
    this->type = object_types::DICTIONARY;
    for (std::pair<std::string, utype> entry: content) {
      new_dictionary->insert({entry.first, entry.second});
    }
  }

  /// @brief converting constructor: set the content of this JSON::object based on a char
  JSON::object::object(const char* content) {
    this->fromString(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a string
  JSON::object::object(const std::string &content) {
    this->fromString(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a bool
  JSON::object::object(bool content) {
    if (content) this->type = object_types::T;
    else this->type = object_types::F;
  }

  /// @brief converting constructor: set the content of this JSON::object to null
  /// @note pass in a nullptr
  JSON::object::object(std::nullptr_t content) {
    this->type = object_types::NONE;
    this->obj = content;
  }

  /// @brief converting constructor: set the content of this JSON::object based on a number
  template<class utype> JSON::object::object(const utype content) {
    this->fromNumber(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a xvector
  template<class utype> JSON::object::object(const xvector<utype> & content) {
    this->fromXvector(content);
  }
  // template instantiation for xvector types
  template JSON::object::object(const xvector<char> &);
  template JSON::object::object(const xvector<int> &);
  template JSON::object::object(const xvector<unsigned int> &);
  template JSON::object::object(const xvector<long long> &);
  template JSON::object::object(const xvector<unsigned long long> &);
  template JSON::object::object(const xvector<float> &);
  template JSON::object::object(const xvector<double> &);

  /// @brief converting constructor: set the content of this JSON::object based on a xmatrix
  template<class utype> JSON::object::object(const xmatrix<utype> & content) {
    this->template fromXmatrix(content);
  }
  // template instantiation for xmatrix types
  template JSON::object::object(const xmatrix<int> &);
  template JSON::object::object(const xmatrix<unsigned int> &);
  template JSON::object::object(const xmatrix<long long> &);
  template JSON::object::object(const xmatrix<unsigned long long> &);
  template JSON::object::object(const xmatrix<float> &);
  template JSON::object::object(const xmatrix<double> &);

  /// @brief converting constructor: set the content of this JSON::object based on a vector
  template<class utype> JSON::object::object(const vector<utype> & content) {
    this->fromVector(content);
  }
  // template instantiation for vector types
  template JSON::object::object(const std::vector<int> &);
  template JSON::object::object(const std::vector<unsigned int> &);
  template JSON::object::object(const std::vector<long long> &);
  template JSON::object::object(const std::vector<unsigned long long> &);
  template JSON::object::object(const std::vector<float> &);
  template JSON::object::object(const std::vector<double> &);
  template JSON::object::object(const std::vector<std::string> &);
  template JSON::object::object(const std::vector<std::vector<int>> &);
  template JSON::object::object(const std::vector<std::vector<unsigned int>> &);
  template JSON::object::object(const std::vector<std::vector<long long>> &);
  template JSON::object::object(const std::vector<std::vector<unsigned long long>> &);
  template JSON::object::object(const std::vector<std::vector<float>> &);
  template JSON::object::object(const std::vector<std::vector<double>> &);
  template JSON::object::object(const std::vector<std::vector<std::string>> &);

  /// @brief converting constructor: set the content of this JSON::object based on a map
  template<class utype> JSON::object::object(const std::map<std::string, utype> & content) {
    this->template fromMap(content);
  }
  // template instantiation for map types (key is always a sting in JSON)
  template JSON::object::object(const std::map<std::string, int> &);
  template JSON::object::object(const std::map<std::string, unsigned int> &);
  template JSON::object::object(const std::map<std::string, long long> &);
  template JSON::object::object(const std::map<std::string, unsigned long long> &);
  template JSON::object::object(const std::map<std::string, float> &);
  template JSON::object::object(const std::map<std::string, double> &);
  template JSON::object::object(const std::map<std::string, std::string> &);
  template JSON::object::object(const std::map<std::string, std::vector<float>> &);

  /// @brief create an empty JSON::object of a given JSON::object_types
  /// @authors
  /// @mod{HE,20221031,created}
  JSON::object::object(object_types create_type) {
    this->type = create_type;
    switch (create_type) {
      {
        case object_types::DICTIONARY: {
          std::shared_ptr <JSON::Dictionary> content = std::make_shared<JSON::Dictionary>();
          this->obj = content;
          cout << "define content" << endl;
          break;
        }
        case object_types::LIST: {
          std::shared_ptr <JSON::List> content = std::make_shared<JSON::List>();
          this->obj = content;
          break;
        }
        case object_types::STRING: {
          std::shared_ptr <std::string> content = std::make_shared<std::string>();
          this->obj = content;
          break;
        }
        case object_types::INTEGER: {
          std::shared_ptr<long long int> content = std::make_shared<long long int>();
          this->obj = content;
          break;
        }
        case object_types::FLOAT: {
          std::shared_ptr<double> content = std::make_shared<double>();
          this->obj = content;
          break;
        }
        default: {
          this->obj = nullptr;
        }
      }
    }
  }


  ///@brief assignment operator for char
  void JSON::object::operator=(const char* content) {
    this->fromString(content);
  }

  ///@brief assignment operator for string
  void JSON::object::operator=(const std::string & content) {
    this->fromString(content);
  }

  ///@brief assignment operator for bool
  void JSON::object::operator=(bool content) {
      if (content) this->type = object_types::T;
      else this->type = object_types::F;
  }

  ///@brief assignment operator for nullptr
  void JSON::object::operator=(std::nullptr_t content) {
    this->type = object_types::NONE;
    this->obj = content;
  }

  ///@brief assignment operator for numbers
  template<class utype> void JSON::object::operator=(const utype content){
    this->template fromNumber(content);
  }
  // template instantiation for number types
  template void JSON::object::operator=(const char);
  template void JSON::object::operator=(const int);
  template void JSON::object::operator=(const unsigned int);
  template void JSON::object::operator=(const long long);
  template void JSON::object::operator=(const unsigned long long);
  template void JSON::object::operator=(const double);
  template void JSON::object::operator=(const float);

  template<class utype> void JSON::object::operator=(const xvector<utype> & content){
    fromXvector(content);
  }

  template void JSON::object::operator=(const xvector<int> &);
  template void JSON::object::operator=(const xvector<unsigned int> &);
  template void JSON::object::operator=(const xvector<long long> &);
  template void JSON::object::operator=(const xvector<unsigned long long> &);
  template void JSON::object::operator=(const xvector<float> &);
  template void JSON::object::operator=(const xvector<double> &);


  template<class utype> void JSON::object::operator=(const xmatrix<utype> & content){
    fromXmatrix(content);
  }

  template void JSON::object::operator=(const xmatrix<int> &);
  template void JSON::object::operator=(const xmatrix<unsigned int> &);
  template void JSON::object::operator=(const xmatrix<long long> &);
  template void JSON::object::operator=(const xmatrix<unsigned long long> &);
  template void JSON::object::operator=(const xmatrix<float> &);
  template void JSON::object::operator=(const xmatrix<double> &);

  template<class utype> void JSON::object::operator=(const std::vector<utype> & content){
    fromVector(content);
  }

  template void JSON::object::operator=(const std::vector<int> &);
  template void JSON::object::operator=(const std::vector<unsigned int> &);
  template void JSON::object::operator=(const std::vector<long long> &);
  template void JSON::object::operator=(const std::vector<unsigned long long> &);
  template void JSON::object::operator=(const std::vector<float> &);
  template void JSON::object::operator=(const std::vector<double> &);
  template void JSON::object::operator=(const std::vector<std::string> &);
  template void JSON::object::operator=(const std::vector<vector<int>> &);
  template void JSON::object::operator=(const std::vector<vector<unsigned int>> &);
  template void JSON::object::operator=(const std::vector<vector<long long>> &);
  template void JSON::object::operator=(const std::vector<vector<unsigned long long>> &);
  template void JSON::object::operator=(const std::vector<vector<float>> &);
  template void JSON::object::operator=(const std::vector<vector<double>> &);
  template void JSON::object::operator=(const std::vector<vector<std::string>> &);

  template<class utype> void JSON::object::operator=(const std::map<std::string, utype> & content){
    fromMap(content);
  }

  template void JSON::object::operator=(const std::map<std::string, int> &);
  template void JSON::object::operator=(const std::map<std::string, unsigned int> &);
  template void JSON::object::operator=(const std::map<std::string, long long> &);
  template void JSON::object::operator=(const std::map<std::string, unsigned long long> &);
  template void JSON::object::operator=(const std::map<std::string, float> &);
  template void JSON::object::operator=(const std::map<std::string, double> &);
  template void JSON::object::operator=(const std::map<std::string, std::string> &);
  template void JSON::object::operator=(const std::map<std::string, std::vector<float>> &);


  //conversion functions
  JSON::object::operator bool() const {
    switch (type) {
      {
        case object_types::DICTIONARY: {
          std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
          if (content->empty()) return false;
          else return true;
        }
        case object_types::LIST: {
          std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
          if (content->empty()) return false;
          else return true;
        }
        case object_types::STRING: {
          std::shared_ptr <std::string> content = std::static_pointer_cast<std::string>(obj);
          if (content->empty()) return false;
          else return true;
        }
        case object_types::INTEGER: {
          std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
          if (*content) return true;
          else return false;
        }
        case object_types::FLOAT: {
          std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
          if (*content) return true;
          else return false;
        }
        case object_types::T: {
          return true;
        }
        case object_types::F: {
          return false;
        }
        case object_types::NONE: {
          return false;
        }
        default: {
          return false;
        }
      }
    }
  }

  JSON::object::operator double() const {
    switch (type) {
      {
        case object_types::INTEGER: {
          std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
          return (double) *content;
        }
        case object_types::FLOAT: {
          std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
          return *content;
        }
        case object_types::T:
          return 1.0;
        case object_types::F:
          return 0.0;
        case object_types::NONE:
          return NAN;
        default:
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON double conversion failed: is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
      }
    }
  }

  JSON::object::operator float() const {
    return (double) *this;
  }

  JSON::object::operator long long() const {
    switch (type) {
      {
        case object_types::INTEGER: {
          std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
          return *content;
        }
        case object_types::FLOAT: {
          std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
          return (long long) *content;
        }
        case object_types::T:
          return 1;
        case object_types::F:
          return 0;
        default:
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON long long conversion failed: element is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
      }
    }
  }

  JSON::object::operator unsigned long long() const {
    long long result = (long long) *this;
    if (result >= 0) return result;
    else throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON unsigned long long conversion failed: element is not a positive number: " + (std::string) * this, _VALUE_ILLEGAL_);
  }

  JSON::object::operator unsigned long() const {
    return (unsigned long long) *this;
  }

  JSON::object::operator unsigned int() const {
    return (unsigned long long) *this;
  }

  JSON::object::operator long() const {
    return (long long) *this;
  }

  JSON::object::operator int() const {
    return (long long) *this;
  }

  JSON::object::operator std::string() const {
    return this->toString(false, false);
  }

  template<class utype> JSON::object::operator std::vector<utype>() const {
    if (type != object_types::LIST)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a LIST: " + (std::string) * this, _VALUE_ILLEGAL_);
    std::vector <utype> result;
    std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    if (content->empty()) return result;
    if (typeid(utype) == typeid(std::string) || typeid(utype) == typeid(bool)) {
      for (const JSON::object &so: *content)
        result.push_back((utype) so); // should be emplace_back - change in cpp14 (error for bool)
    } else if (std::is_floating_point<utype>::value) {
      for (const JSON::object &so: *content) {
        if (so.type > object_types::STRING) result.push_back((utype) so);
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (const JSON::object &so: *content) {
        if (so.type > object_types::STRING && so.type < object_types::NONE) result.push_back((utype) so);
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
      }
    }
    return result;
  }

  template JSON::object::operator std::vector<std::string>() const;
  template JSON::object::operator std::vector<double>() const;
  template JSON::object::operator std::vector<float>() const;
  template JSON::object::operator std::vector<long long>() const;
  template JSON::object::operator std::vector<long>() const;
  template JSON::object::operator std::vector<int>() const;
  template JSON::object::operator std::vector<uint>() const;
  template JSON::object::operator std::vector<bool>() const;

  JSON::object::operator std::map<std::string, JSON::object> () const {
    if (type != object_types::DICTIONARY)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a DICTIONARY: " + (std::string) *this, _VALUE_ILLEGAL_);
    std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
    return *content;
  }

  template<class utype> JSON::object::operator std::map<std::string, utype>() const {
    if (type != object_types::DICTIONARY)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a DICTIONARY: " + (std::string) *this, _VALUE_ILLEGAL_);
    std::map <std::string, utype> result;
    std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
    if (content->empty()) return result;
    if (typeid(utype) == typeid(std::string) || typeid(utype) == typeid(bool)) {
      for (std::pair<const std::string, object> &entry: *content)
        result.insert({entry.first, (utype) entry.second});
    } else if (std::is_floating_point<utype>::value) {
      for (std::pair<const std::string, object> &entry: *content) {
        if (entry.second.type > object_types::STRING) result.insert({entry.first, (utype) entry.second});
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a number: " + (std::string) entry.second, _VALUE_ILLEGAL_);
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (std::pair<const std::string, object> &entry: *content) {
        if (entry.second.type > object_types::STRING && entry.second.type < object_types::NONE)
          result.insert({entry.first, (utype) entry.second});
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a number: " + (std::string) entry.second, _VALUE_ILLEGAL_);
      }
    }
    return result;
  }

  template JSON::object::operator std::map<std::string, std::string>() const;
  template JSON::object::operator std::map<std::string, double>() const;
  template JSON::object::operator std::map<std::string, float>() const;
  template JSON::object::operator std::map<std::string, long long>() const;
  template JSON::object::operator std::map<std::string, long>() const;
  template JSON::object::operator std::map<std::string, int>() const;
  template JSON::object::operator std::map<std::string, uint>() const;
  template JSON::object::operator std::map<std::string, bool>() const;

  template<class utype> JSON::object::operator aurostd::xvector<utype>() const {
    if (type != object_types::LIST)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xvector conversion failed: element is not a LIST: " + (std::string) * this, _VALUE_ILLEGAL_);
    std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    aurostd::xvector<utype> result(content->size(), 1);
    if (content->empty()) return result;
    size_t idx = 1;
    if (typeid(utype) == typeid(bool)) {
      for (const JSON::object &so: *content) {
        result[idx] = (utype) so;
        idx++;
      }
    } else if (std::is_floating_point<utype>::value) {
      for (const JSON::object &so: *content) {
        if (so.type > object_types::STRING) result[idx] = (utype) so;
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
        idx++;
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (const JSON::object &so: *content) {
        if (so.type > object_types::STRING && so.type < object_types::NONE) result[idx] = (utype) so;
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
        idx++;
      }
    }
    return result;
  }

  template JSON::object::operator aurostd::xvector<double>() const;
  template JSON::object::operator aurostd::xvector<float>() const;
  template JSON::object::operator aurostd::xvector<long long>() const;
  template JSON::object::operator aurostd::xvector<long>() const;
  template JSON::object::operator aurostd::xvector<int>() const;
  template JSON::object::operator aurostd::xvector<uint>() const;
  template JSON::object::operator aurostd::xvector<bool>() const;

  template<class utype> JSON::object::operator aurostd::xmatrix<utype>() const {
    if (type != object_types::LIST)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: element is not a LIST: " + (std::string) *this, _VALUE_ILLEGAL_);
    std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    if (content->empty()) return aurostd::xmatrix<utype>();

    // size scan
    size_t rows = 0;
    size_t cols = 0;
    std::vector <size_t> cols_v;
    for (const JSON::object &so: *content) {
      if (so.type != object_types::LIST)
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: sub element is not a LIST: " + (std::string) *this, _VALUE_ILLEGAL_);
      std::shared_ptr <JSON::List> row = std::static_pointer_cast<JSON::List>(so.obj);
      rows++;
      cols_v.template emplace_back(row->size());
    }
    if (std::adjacent_find(cols_v.begin(), cols_v.end(), std::not_equal_to<size_t>()) == cols_v.end()) cols = cols_v[0];
    else
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: sub elements have different length", _VALUE_ILLEGAL_);

    aurostd::xmatrix<utype> result(rows, cols);
    if (typeid(utype) == typeid(bool)) {
      for (int r = result.lrows; r <= result.urows; r++) {
        std::shared_ptr <JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          result[r][c] = (utype) row->operator[](c - 1);
        }
      }
    } else if (std::is_floating_point<utype>::value) {
      for (int r = result.lrows; r <= result.urows; r++) {
        std::shared_ptr <JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          if (row->operator[](c - 1).type > object_types::STRING) result[r][c] = (utype) row->operator[](c - 1);
          else
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: elements is not a number: " + (string) row->operator[](c - 1), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (int r = result.lrows; r <= result.urows; r++) {
        std::shared_ptr <JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          if (row->operator[](c - 1).type > object_types::STRING && row->operator[](c - 1).type < object_types::NONE)
            result[r][c] = (utype) row->operator[](c - 1);
          else
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: elements is not a number: " + (string) row->operator[](c - 1), _VALUE_ILLEGAL_);
        }
      }
    }
    return result;
  }

  template JSON::object::operator aurostd::xmatrix<double>() const;
  template JSON::object::operator aurostd::xmatrix<float>() const;
  template JSON::object::operator aurostd::xmatrix<long long>() const;
  template JSON::object::operator aurostd::xmatrix<long>() const;
  template JSON::object::operator aurostd::xmatrix<int>() const;
  template JSON::object::operator aurostd::xmatrix<uint>() const;
  template JSON::object::operator aurostd::xmatrix<bool>() const;


  template<class utype>
    void JSON::object::push_back(const utype content){
    if (type == JSON::object_types::LIST) {
      std::shared_ptr <JSON::List> list_obj = std::static_pointer_cast<JSON::List>(obj);
      list_obj->push_back(content);
    }
    else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "push_back is just allowed for a JSON LIST ", _VALUE_ILLEGAL_);

    }

  }
  template void JSON::object::push_back(const JSON::object content);


}
namespace aurostd {
  /// @class JSON
  /// @brief unified class to read and write JSON
  ///
  /// @authors
  /// @mod{AS,2020,created JSONWriter}
  /// @mod{HE,20220924,rewrite to enable parsing}
  ///
  /// Basic usage
  /// @code
  /// aurostd::JSON jr;
  /// jr.loadFile("testcases.json");
  /// // output complete file as JSON
  /// cout << jr << endl;
  /// // output element
  /// cout << jr["xvector"][3] << endl;
  /// // save element
  /// xvector<double> xvd = jr["xvector"];
  /// cout << xvd << endl;
  /// // save sub element
  /// uint el_uint = jr["xvector"][3];
  /// cout << el_uint << endl;
  /// // cast to type
  /// cout << (float) jr["xvector"][3] << endl;
  /// @endcode

  void JSON::clear() { *this = {}; }  // calling the constructor
  void JSON::copy(const JSON &jr) {
    root = jr.root;
  }

  JSON::JSON() {}

  JSON::JSON(const JSON &jr) {
    copy(jr);
  }

  JSON::~JSON() {}

  const JSON &JSON::operator=(const JSON &jr) {
    copy(jr);
    return *this;
  }

  ostream &operator<<(ostream &os, const JSON &jr) {
    os << jr.toString(false);
    return os;
  }

  JSON::object &JSON::operator[](const size_t index) {
    return root[index];
  }

  JSON::object &JSON::operator[](const std::string key) {
    return root[key];
  }

  JSON::object &JSON::operator[](const char *key) {
    return root[key];
  }

  static inline size_t range_find(const char *content_ptr, std::pair <size_t, size_t> border, char to_find) {
    return (char *) memchr(content_ptr + border.first, to_find, border.second - border.first + 1) - content_ptr;
  }


  std::string JSON::parse_string(const std::string &raw_content, std::pair <size_t, size_t> border) const {
    if (border.second == 0) {
      border.second = raw_content.size() - 1;
    }
    size_t last_escape_pos = border.first;
    size_t current_escape_pos = range_find(raw_content.c_str(), border, '\\');

    if (current_escape_pos > border.second) return raw_content.substr(border.first, border.second - border.first + 1);

    std::string result = "";
    while (current_escape_pos <= border.second) {
      result += raw_content.substr(last_escape_pos, current_escape_pos - last_escape_pos);
      switch (raw_content[current_escape_pos + 1]) {
        case '"': {
          result += '"';
        }
          break;
        case '\\': {
          result += '\\';
        }
          break;
        case '/': {
          result += '/';
        }
          break;
        case 'b': {
          result += '\b';
        }
          break;
        case 'f': {
          result += '\f';
        }
          break;
        case 'n': {
          result += '\n';
        }
          break;
        case 'r': {
          result += '\r';
        }
          break;
        case 't': {
          result += '\t';
        }
          break;
        case 'u': {
          if (current_escape_pos + 6 > raw_content.size())
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
          result += unescape_unicode(raw_content, current_escape_pos);
          current_escape_pos += 2;
        }
          break;
        default: {
          stringstream message;
          message << "JSON parsing failed: string contains undefined escape '\\" << raw_content[current_escape_pos + 1] << "'";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message.str(), _FILE_WRONG_FORMAT_);
        }
      }
      last_escape_pos = current_escape_pos + 2;
      current_escape_pos = range_find(raw_content.c_str(), {last_escape_pos, border.second}, '\\');
    }
    result += raw_content.substr(last_escape_pos, border.second - last_escape_pos + 1);
    return result;
  }

  std::string JSON::escape_unicode(const std::string & raw) {
    std::stringstream out;
    std::u16string utf16 = std::wstring_convert<std::codecvt_utf8_utf16<char16_t>, char16_t>{}.from_bytes(raw.data());
    for (char16_t c : utf16){
      if (c<0x80) out << (char)c;
      else out << "\\u"  << std::hex << std::setw(4)  << std::setfill('0') <<  c;
    }
    return out.str();
  }

  std::string JSON::char32_to_string(const char32_t cp)  {
    string out;
    if (cp < 0x80) {
      out += static_cast<char>(cp);
      return out;
    }

    if (cp < 0x800) {
      out += static_cast<char>((cp >> 6) | 0xc0);
      out += static_cast<char>((cp & 0x3f) | 0x80);
      return out;
    }

    if (cp < 0x10000) {
      out += static_cast<char>((cp >> 12) | 0xe0);
      out += static_cast<char>(((cp >> 6) & 0x3f) | 0x80);
      out += static_cast<char>((cp & 0x3f) | 0x80);
      return out;
    }

    {
      out += static_cast<char>((cp >> 18) | 0xf0);
      out += static_cast<char>(((cp >> 12) & 0x3f) | 0x80);
      out += static_cast<char>(((cp >> 6) & 0x3f) | 0x80);
      out += static_cast<char>((cp & 0x3f) | 0x80);
      return out;
    }
  }

  std::string JSON::unescape_unicode(const std::string & raw, size_t & pos) {
    pos += 2;
    char32_t cp = (char32_t) aurostd::string2utype<uint>(raw.substr(pos, 4), 16);

    if (cp < 0xd800 || cp > 0xdfff) return char32_to_string(cp);
    else if (cp > 0xdbff) throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
    else {
      if (raw[pos+4] == '\\' and raw[pos+5] == 'u'){
        pos += 6;
        char32_t trailing_cp = (char32_t) aurostd::string2utype<uint>(raw.substr(pos, 4), 16);
        if (trailing_cp < 0xdc00 || trailing_cp > 0xdfff) throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
        char32_t combo_cp = ((cp - 0xd800) << 10) + (trailing_cp - 0xdc00) + 0x10000;
        return char32_to_string(combo_cp);
      }
      else throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
    }
  }

  std::pair <size_t, size_t> JSON::find_strip(const std::string &raw_content, std::pair <size_t, size_t> border) const {
    if (border.second == 0) border.second = raw_content.size() - 1;
    size_t start = raw_content.find_first_not_of(" \n\t\r\v\f", border.first);
    size_t end = raw_content.find_last_not_of(" \n\t\r\v\f", border.second);
    return {start, min(end, border.second)};
  }


  std::pair <size_t, size_t> JSON::find_string(const std::string &raw_content, std::pair <size_t, size_t> border) const {
    if (border.second == 0) border.second = raw_content.size() - 1;
    size_t start = range_find(raw_content.c_str(), border, '"');
    if (start >= border.second) return {std::string::npos, std::string::npos};
    size_t end = start;
    if (start == std::string::npos) return {start, end};
    do {
      end = range_find(raw_content.c_str(), {end + 1, border.second}, '"');
    } while (raw_content[end - 1] == '\\' && (end <= border.second));
    if (end > border.second)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: string not closed", _FILE_WRONG_FORMAT_);
    return {start, end};
  }

  std::pair <size_t, size_t> JSON::find_bracket(const std::string &raw_content, char kind_open, std::pair <size_t, size_t> border) const {
    char kind_close;
    switch (kind_open) {
      case '[':
        kind_close = ']';
        break;
      case '{':
        kind_close = '}';
        break;
      case '(':
        kind_close = ')';
        break;
      default:
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: bracket kind not allowed", _FILE_WRONG_FORMAT_);
    }
    size_t start = raw_content.find(kind_open, border.first);
    size_t end = start;
    if (start > border.second) return {std::string::npos, std::string::npos};
    size_t next_open;
    size_t next_close;
    std::pair <size_t, size_t> string_section;
    size_t open_count = 1;
    do {
//       cout <<"end: " << end << " | "  << "open: " << open_count << " | current selection: " << raw_content.substr(start, end-start+1) << endl;
      string_section = find_string(raw_content, {end + 1, border.second});
      next_close = range_find(raw_content.c_str(), {end + 1, border.second}, kind_close);
      next_open = range_find(raw_content.c_str(), {end + 1, border.second}, kind_open);
      if (next_close > string_section.first && next_open > string_section.first) {
        end = string_section.second;
        continue;
      }
      if (next_close < next_open) {
        if (next_close < string_section.first) {
          end = next_close;
          open_count--;
        }
      } else if (next_close > next_open) {
        if (next_open < string_section.first) {
          end = next_open;
          open_count++;
        }
      } else { // when == std::string::npos
        break;
      }

    } while (open_count > 0 && (end < border.second));

    return {start, end};

  }

  JSON::object JSON::parse(const std::string &raw_content, std::pair <size_t, size_t> border) const {
    if (border.second == 0) border.second = raw_content.size() - 1;
    object result;
    result.type = object_types::NONE;
    result.obj = nullptr;
    border = find_strip(raw_content, border);
    std::pair <size_t, size_t> section({0, 0});
    switch (raw_content[border.first]) {
      case '{': {
        if (raw_content[border.second] != '}')
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: object not closed by '}'", _FILE_WRONG_FORMAT_);
        std::shared_ptr <JSON::Dictionary> new_dictionary = std::make_shared<JSON::Dictionary>();
        border = {border.first + 1, border.second - 1};
        section = {0, 0};
        do {
          border = find_strip(raw_content, border);
//          cout  << border.first << " | "<< border.second << " - current: " << raw_content.substr(border.first, border.second-border.first+1) << endl;
          if (border.first >= border.second) break;
          section = find_string(raw_content, border);
          if (section.first == std::string::npos) break;
          std::string key_string = parse_string(raw_content, {section.first + 1, section.second - 1}); //cut " already
//          cout << "key: " << key_string << endl;
          border.first = raw_content.find(':', section.second) + 1;
          border = find_strip(raw_content, border);
//          cout  << border.first << " | "<< border.second << " - switch: " << raw_content.substr(border.first, border.second-border.first+1) << endl;

          switch (raw_content[border.first]) {
            case '"': {
              section = find_string(raw_content, border);
            }
              break;
            case '[': {
              section = find_bracket(raw_content, '[', border);
            }
              break;
            case '{': {
              section = find_bracket(raw_content, '{', border);
            }
              break;
            default: {
              section.first = border.first;
              section.second = min(border.second, range_find(raw_content.c_str(), border, ',') - 1);
            }
          }
          new_dictionary->insert({key_string, parse(raw_content, section)});
          border.first = range_find(raw_content.c_str(), {section.second, border.second}, ',');
          if (border.first == std::string::npos) border.first = border.second;
          else border.first++;
        } while (border.first < border.second);
        result.obj = new_dictionary;
        result.type = object_types::DICTIONARY;
      }
        break;
      case '[': {
        if (raw_content[border.second] != ']')
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: list not closed by ']'", _FILE_WRONG_FORMAT_);
        std::shared_ptr <JSON::List> new_list = std::make_shared<JSON::List>();
        border = {border.first + 1, border.second - 1};
        section = {0, 0};
        do {
          border = find_strip(raw_content, border);
          switch (raw_content[border.first]) {
            case '"': {
              section = find_string(raw_content, border);
            }
              break;
            case '[': {
              section = find_bracket(raw_content, '[', border);
            }
              break;
            case '{': {
              section = find_bracket(raw_content, '{', border);
            }
              break;
            default: {
              section.first = border.first;
              section.second = min(border.second, raw_content.find(',', border.first) - 1);
            }
          }
          new_list->emplace_back(parse(raw_content, section));
          border.first = raw_content.find(',', section.second);
          if (border.first == std::string::npos) border.first = border.second;
          else border.first++;
        } while (section.second < border.second);
        result.obj = new_list;
        result.type = object_types::LIST;
      }
        break;
      case '"': {
        if (raw_content[border.second] != '"')
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: string not enclosed by '\\\"'\"", _FILE_WRONG_FORMAT_);
        std::shared_ptr <std::string> new_string = std::make_shared<std::string>();
        *new_string = parse_string(raw_content, {border.first + 1, border.second - 1});
        result.obj = new_string;
        result.type = object_types::STRING;
      }
        break;
      case 'n': {
        result.obj = nullptr;
        result.type = object_types::NONE;
      }
        break;
      case 't': {
        result.obj = nullptr;
        result.type = object_types::T;
      }
        break;
      case 'f': {
        result.obj = nullptr;
        result.type = object_types::F;
      }
        break;
      default: { //number
        if (raw_content.find('.', border.first) != std::string::npos || raw_content.find('E', border.first) != std::string::npos || raw_content.find('e', border.first) != std::string::npos) {
          std::shared_ptr<double> new_number = std::make_shared<double>();
          *new_number = std::strtod(raw_content.c_str() + border.first, nullptr);
          result.obj = new_number;
          result.type = object_types::FLOAT;
        } else {
          std::shared_ptr<long long> new_number = std::make_shared<long long>();
          *new_number = std::strtoll(raw_content.c_str() + border.first, nullptr, 10);
          result.obj = new_number;
          result.type = object_types::INTEGER;
        }
      }
    }
//    cout << "result: " << obj_to_string(result) << endl;
    return result;
  }

  void JSON::loadFile(const std::string &file_path) {
    std::string raw_content = aurostd::file2string(file_path);
    root = parse(raw_content);
  }

  void JSON::loadString(const std::string &content) {
    root = parse(content);
  }

  std::string JSON::toString(const bool escape) const {
    return root.toString(true, escape);
  }

  JSON::operator object() const{
    return this->root;
  }

}

#endif // _AUROSTD_XPARSER_JSON_CPP_