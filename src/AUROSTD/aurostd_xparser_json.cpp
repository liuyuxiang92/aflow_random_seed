// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************
// This JSON class is the evolution of different prior solutions to integrate with the AFLOW source base.
// hagen.eckert@duke.edu


#ifndef _AUROSTD_XPARSER_JSON_CPP_
#define _AUROSTD_XPARSER_JSON_CPP_

#ifndef _AUROSTD_XPARSER_H_

#include "aurostd_xparser.h"

#endif


namespace aurostd {

  /// @struct JSONReader::storage_object
  /// @brief storge container for JSONReader
  ///
  /// @authors
  /// @mod{HE,20220924,created struct}

  JSONReader::storage_object &JSONReader::storage_object::operator[](const size_t index) {
    if (this->type == JSONReader::object_types::LIST) {
      std::shared_ptr <JSONReader::List> content = std::static_pointer_cast<JSONReader::List>(this->obj);
      return content->operator[](index);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Not a list", _INDEX_ILLEGAL_);
    }
  }

  JSONReader::storage_object &JSONReader::storage_object::operator[](const std::string key) {
    if (this->type == JSONReader::object_types::DICTIONARY) {
      std::shared_ptr <JSONReader::Dictionary> content = std::static_pointer_cast<JSONReader::Dictionary>(this->obj);
      return content->operator[](key);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Not a dictionary", _INDEX_ILLEGAL_);
    }
  }

  JSONReader::storage_object &JSONReader::storage_object::operator[](const char *key) {
    if (this->type == JSONReader::object_types::DICTIONARY) {
      std::shared_ptr <JSONReader::Dictionary> content = std::static_pointer_cast<JSONReader::Dictionary>(this->obj);
      return content->operator[](key);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Not a dictionary", _INDEX_ILLEGAL_);
    }
  }

  ostream &operator<<(ostream &os, const JSONReader::storage_object &so) {
    os << so.toString();
    return os;
  }

  std::string JSONReader::storage_object::toString(const bool json_format) const {
    bool first = true;
    stringstream result;
    switch (type) {
      case object_types::DICTIONARY: {
        result << "{";
        std::shared_ptr <JSONReader::Dictionary> content = std::static_pointer_cast<JSONReader::Dictionary>(obj);
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
        std::shared_ptr <JSONReader::List> content = std::static_pointer_cast<JSONReader::List>(obj);
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
        if (json_format) result << "\"" << *content << "\""; //TODO escapes
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

  JSONReader::storage_object::operator bool() const {
    switch (type) {
      {
        case object_types::DICTIONARY: {
          std::shared_ptr <JSONReader::Dictionary> content = std::static_pointer_cast<JSONReader::Dictionary>(obj);
          if (content->empty()) return false;
          else return true;
        }
        case object_types::LIST: {
          std::shared_ptr <JSONReader::List> content = std::static_pointer_cast<JSONReader::List>(obj);
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

  JSONReader::storage_object::operator double() const {
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

  JSONReader::storage_object::operator float() const {
    return (double) *this;
  }

  JSONReader::storage_object::operator long long() const {
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

  JSONReader::storage_object::operator unsigned long long() const {
    long long result = (long long) *this;
    if (result >= 0) return result;
    else throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON unsigned long long conversion failed: element is not a positive number: " + (std::string) * this, _VALUE_ILLEGAL_);
  }

  JSONReader::storage_object::operator unsigned long() const {
    return (unsigned long long) *this;
  }

  JSONReader::storage_object::operator std::string() const {
    return this->toString(false);
  }

  template<class utype> JSONReader::storage_object::operator std::vector<utype>() const {
    if (type != object_types::LIST)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a LIST: " + (std::string) * this, _VALUE_ILLEGAL_);
    std::vector <utype> result;
    std::shared_ptr <JSONReader::List> content = std::static_pointer_cast<JSONReader::List>(obj);
    if (content->empty()) return result;
    if (typeid(utype) == typeid(std::string) || typeid(utype) == typeid(bool)) {
      for (const JSONReader::storage_object &so: *content)
        result.push_back((utype) so); // should be emplace_back - change in cpp14 (error for bool)
    } else if (std::is_floating_point<utype>::value) {
      for (const JSONReader::storage_object &so: *content) {
        if (so.type > object_types::STRING) result.push_back((utype) so);
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (const JSONReader::storage_object &so: *content) {
        if (so.type > object_types::STRING && so.type < object_types::NONE) result.push_back((utype) so);
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
      }
    }
    return result;
  }

  template JSONReader::storage_object::operator std::vector<std::string>() const;
  template JSONReader::storage_object::operator std::vector<double>() const;
  template JSONReader::storage_object::operator std::vector<float>() const;
  template JSONReader::storage_object::operator std::vector<long long>() const;
  template JSONReader::storage_object::operator std::vector<long>() const;
  template JSONReader::storage_object::operator std::vector<int>() const;
  template JSONReader::storage_object::operator std::vector<uint>() const;
  template JSONReader::storage_object::operator std::vector<bool>() const;

  template<class utype> JSONReader::storage_object::operator std::map<std::string, utype>() const {
    if (type != object_types::DICTIONARY)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a DICTIONARY: " + (std::string) *this, _VALUE_ILLEGAL_);
    std::map <std::string, utype> result;
    std::shared_ptr <JSONReader::Dictionary> content = std::static_pointer_cast<JSONReader::Dictionary>(obj);
    if (content->empty()) return result;
    if (typeid(utype) == typeid(std::string) || typeid(utype) == typeid(bool)) {
      for (std::pair<const std::string, storage_object> &entry: *content)
        result.insert({entry.first, (utype) entry.second});
    } else if (std::is_floating_point<utype>::value) {
      for (std::pair<const std::string, storage_object> &entry: *content) {
        if (entry.second.type > object_types::STRING) result.insert({entry.first, (utype) entry.second});
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a number: " + (std::string) entry.second, _VALUE_ILLEGAL_);
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (std::pair<const std::string, storage_object> &entry: *content) {
        if (entry.second.type > object_types::STRING && entry.second.type < object_types::NONE)
          result.insert({entry.first, (utype) entry.second});
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a number: " + (std::string) entry.second, _VALUE_ILLEGAL_);
      }
    }
    return result;
  }

  template JSONReader::storage_object::operator std::map<std::string, std::string>() const;
  template JSONReader::storage_object::operator std::map<std::string, double>() const;
  template JSONReader::storage_object::operator std::map<std::string, float>() const;
  template JSONReader::storage_object::operator std::map<std::string, long long>() const;
  template JSONReader::storage_object::operator std::map<std::string, long>() const;
  template JSONReader::storage_object::operator std::map<std::string, int>() const;
  template JSONReader::storage_object::operator std::map<std::string, uint>() const;
  template JSONReader::storage_object::operator std::map<std::string, bool>() const;

  template<class utype> JSONReader::storage_object::operator aurostd::xvector<utype>() const {
    if (type != object_types::LIST)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xvector conversion failed: element is not a LIST: " + (std::string) * this, _VALUE_ILLEGAL_);
    std::shared_ptr <JSONReader::List> content = std::static_pointer_cast<JSONReader::List>(obj);
    aurostd::xvector<utype> result(content->size(), 1);
    if (content->empty()) return result;
    size_t idx = 1;
    if (typeid(utype) == typeid(bool)) {
      for (const JSONReader::storage_object &so: *content) {
        result[idx] = (utype) so;
        idx++;
      }
    } else if (std::is_floating_point<utype>::value) {
      for (const JSONReader::storage_object &so: *content) {
        if (so.type > object_types::STRING) result[idx] = (utype) so;
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
        idx++;
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (const JSONReader::storage_object &so: *content) {
        if (so.type > object_types::STRING && so.type < object_types::NONE) result[idx] = (utype) so;
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
        idx++;
      }
    }
    return result;
  }

  template JSONReader::storage_object::operator aurostd::xvector<double>() const;
  template JSONReader::storage_object::operator aurostd::xvector<float>() const;
  template JSONReader::storage_object::operator aurostd::xvector<long long>() const;
  template JSONReader::storage_object::operator aurostd::xvector<long>() const;
  template JSONReader::storage_object::operator aurostd::xvector<int>() const;
  template JSONReader::storage_object::operator aurostd::xvector<uint>() const;
  template JSONReader::storage_object::operator aurostd::xvector<bool>() const;

  template<class utype> JSONReader::storage_object::operator aurostd::xmatrix<utype>() const {
    if (type != object_types::LIST)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: element is not a LIST: " + (std::string) *this, _VALUE_ILLEGAL_);
    std::shared_ptr <JSONReader::List> content = std::static_pointer_cast<JSONReader::List>(obj);
    if (content->empty()) return aurostd::xmatrix<utype>();

    // size scan
    size_t rows = 0;
    size_t cols = 0;
    std::vector <size_t> cols_v;
    for (const JSONReader::storage_object &so: *content) {
      if (so.type != object_types::LIST)
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: sub element is not a LIST: " + (std::string) *this, _VALUE_ILLEGAL_);
      std::shared_ptr <JSONReader::List> row = std::static_pointer_cast<JSONReader::List>(so.obj);
      rows++;
      cols_v.template emplace_back(row->size());
    }
    if (std::adjacent_find(cols_v.begin(), cols_v.end(), std::not_equal_to<size_t>()) == cols_v.end()) cols = cols_v[0];
    else
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: sub elements have different length", _VALUE_ILLEGAL_);

    aurostd::xmatrix<utype> result(rows, cols);
    if (typeid(utype) == typeid(bool)) {
      for (int r = result.lrows; r <= result.urows; r++) {
        std::shared_ptr <JSONReader::List> row = std::static_pointer_cast<JSONReader::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          result[r][c] = (utype) row->operator[](c - 1);
        }
      }
    } else if (std::is_floating_point<utype>::value) {
      for (int r = result.lrows; r <= result.urows; r++) {
        std::shared_ptr <JSONReader::List> row = std::static_pointer_cast<JSONReader::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          if (row->operator[](c - 1).type > object_types::STRING) result[r][c] = (utype) row->operator[](c - 1);
          else
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: elements is not a number: " + (string) row->operator[](c - 1), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (int r = result.lrows; r <= result.urows; r++) {
        std::shared_ptr <JSONReader::List> row = std::static_pointer_cast<JSONReader::List>(content->operator[](r - 1).obj);
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

  template JSONReader::storage_object::operator aurostd::xmatrix<double>() const;
  template JSONReader::storage_object::operator aurostd::xmatrix<float>() const;
  template JSONReader::storage_object::operator aurostd::xmatrix<long long>() const;
  template JSONReader::storage_object::operator aurostd::xmatrix<long>() const;
  template JSONReader::storage_object::operator aurostd::xmatrix<int>() const;
  template JSONReader::storage_object::operator aurostd::xmatrix<uint>() const;
  template JSONReader::storage_object::operator aurostd::xmatrix<bool>() const;

}
namespace aurostd {
  /// @class JSONReader
  /// @brief unified class to read and write JSON
  ///
  /// @authors
  /// @mod{AS,2020,created JSONWriter}
  /// @mod{HE,20220924,rewrite to enable parsing}
  ///
  /// Basic usage
  /// @code
  /// aurostd::JSONReader jr;
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

  void JSONReader::clear() { *this = {}; }  // calling the constructor
  void JSONReader::copy(const JSONReader &jr) {
    root = jr.root;
  }

  JSONReader::JSONReader() {}

  JSONReader::JSONReader(const JSONReader &jr) {
    copy(jr);
  }

  JSONReader::~JSONReader() {}

  const JSONReader &JSONReader::operator=(const JSONReader &jr) {
    copy(jr);
    return *this;
  }

  ostream &operator<<(ostream &os, const JSONReader &jr) {
    os << jr.toString();
    return os;
  }

  JSONReader::storage_object &JSONReader::operator[](const size_t index) {
    return root[index];
  }

  JSONReader::storage_object &JSONReader::operator[](const std::string key) {
    return root[key];
  }

  JSONReader::storage_object &JSONReader::operator[](const char *key) {
    return root[key];
  }

  static inline size_t range_find(const char *content_ptr, std::pair <size_t, size_t> border, char to_find) {
    return (char *) memchr(content_ptr + border.first, to_find, border.second - border.first + 1) - content_ptr;
  }


  std::string JSONReader::parse_string(const std::string &raw_content, std::pair <size_t, size_t> border) const {
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
          result += utf8(raw_content.substr(current_escape_pos + 2, 4));
          current_escape_pos += 4;
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

  std::string JSONReader::utf8(const std::string &cp) const {
    char32_t unicode_character = (char32_t) aurostd::string2utype<uint>(cp, 16);
    char32_t const *from = &unicode_character;
    char utf8_chars[4];
    char *end_of_utf8;

    std::mbstate_t mbs;
    std::codecvt_utf8 <char32_t> ccv;
//    std::codecvt_base::result
    ccv.out(mbs, from, from + 1, from, utf8_chars, utf8_chars + 4, end_of_utf8);
    // For debug
//    if (result!=std::codecvt_base::result::ok) cerr << "\\u" << cp << " is skipped" << endl;
    return {utf8_chars, end_of_utf8};

  }

  std::pair <size_t, size_t> JSONReader::find_strip(const std::string &raw_content, std::pair <size_t, size_t> border) const {
    if (border.second == 0) border.second = raw_content.size() - 1;
    size_t start = raw_content.find_first_not_of(" \n\t\r\v\f", border.first);
    size_t end = raw_content.find_last_not_of(" \n\t\r\v\f", border.second);
    return {start, min(end, border.second)};
  }


  std::pair <size_t, size_t> JSONReader::find_string(const std::string &raw_content, std::pair <size_t, size_t> border) const {
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

  std::pair <size_t, size_t> JSONReader::find_bracket(const std::string &raw_content, char kind_open, std::pair <size_t, size_t> border) const {
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

  JSONReader::storage_object JSONReader::parse(const std::string &raw_content, std::pair <size_t, size_t> border) const {
    if (border.second == 0) border.second = raw_content.size() - 1;
    storage_object result;
    result.type = object_types::NONE;
    result.obj = nullptr;
    border = find_strip(raw_content, border);
    std::pair <size_t, size_t> section({0, 0});
    switch (raw_content[border.first]) {
      case '{': {
        if (raw_content[border.second] != '}')
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: object not closed by '}'", _FILE_WRONG_FORMAT_);
        std::shared_ptr <JSONReader::Dictionary> new_dictionary = std::make_shared<JSONReader::Dictionary>();
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
        std::shared_ptr <JSONReader::List> new_list = std::make_shared<JSONReader::List>();
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

  void JSONReader::loadFile(const std::string &file_path) {
    std::string raw_content = aurostd::file2string(file_path);
    root = parse(raw_content);
  }

  void JSONReader::loadString(const std::string &content) {
    root = parse(content);
  }

  std::string JSONReader::toString() const {
    return root.toString();
  }

  JSONReader::operator storage_object() const{
    return this->root;
  }

}

#endif // _AUROSTD_XPARSER_JSON_CPP_