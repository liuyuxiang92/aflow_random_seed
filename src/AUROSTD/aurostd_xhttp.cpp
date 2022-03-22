// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *           Aflow HAGEN ECKERT - Duke University 2021-2022                *
// *                                                                         *
// ***************************************************************************
// Written by Hagen Eckert
// hagen.eckert@duke.edu

#ifndef _AUROSTD_XHTTP_CPP_
#define _AUROSTD_XHTTP_CPP_

#include "aurostd_xhttp.h"

#define _DEBUG_XHTTP_ false

namespace aurostd {

  URL httpConstructURL(const std::string &host, const std::string &path="/", const std::string &query="", const unsigned int &port=80){
    URL url;
    url.scheme = "http";
    url.port = port;
    url.host = host;
    url.path = path;
    url.query = query;
    return url;
  }

  /// @brief split an URL string into its parts
  /// @param url url string
  /// @param strict if set, a full URL needs a schema (http://, https://), else it is interpreted as path
  /// @return URL struct
  ///
  /// While a `user` can be extracted by this parser `user:password` is not supported, and should not be added!
  /// https://datatracker.ietf.org/doc/html/rfc3986#section-3.2.1
  URL httpParseURL(const std::string &url, const bool strict) {

    bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XHTTP_);
    std::string soliloquy = XPID + "aurostd::httpParseResponse():";

    if (LDEBUG) {
      if (strict) cerr << soliloquy << " Parse '" << url << "' strict" << endl;
      else cerr << soliloquy << " Parse '" << url << "' lenient" << endl;
    }

    URL result;
    std::string delimiter = "";
    size_t start = 0;
    size_t location = 0;

    // locate schema
    delimiter = "://";
    location = url.find(delimiter);
    if (location != std::string::npos) {
      result.scheme = aurostd::tolower(url.substr(start, location - start));
      start = location + delimiter.length();
    }
    if (!strict || !result.scheme.empty()) {
      // locate username
      location = url.find('@', start);
      if (location != std::string::npos) {
        result.user = url.substr(start, location - start);
        start = location + 1;
      }

      // locate host
      location = url.find('/', start);
      if (location != std::string::npos) {
        result.host = url.substr(start, location - start);
        start = location + 1;
      } else {
        result.host = url.substr(start, url.length() - start);
        start = url.length();
      }

      // locate port (default to 80 or 443)
      location = result.host.find(':');
      if (location != std::string::npos) {
        result.port = aurostd::string2utype(result.host.substr(location + 1));
        result.host.erase(location);
      } else {
        if (result.scheme == "https") result.port = 443;
        else if (result.scheme == "ftp") result.port = 21;
        else if (result.scheme == "sftp") result.port = 22;
        else result.port = 80;
      }

      // locate path
      result.path = "/" + url.substr(start);
    } else {
      result.path = url.substr(start);
      result.port = 0;
    }

    // split query from raw path
    location = result.path.find('?');
    if (location != std::string::npos) {
      result.query = result.path.substr(location);
      result.path = result.path.substr(0, location);
    } else {
      result.query = "";
    }

    if (LDEBUG) {
      cerr << soliloquy << " split up URL: " << endl;
      cerr << "    " << "scheme: " << result.scheme << std::endl;
      cerr << "    " << "user: " << result.user << std::endl;
      cerr << "    " << "host: " << result.host << std::endl;
      cerr << "    " << "port: " << result.port << std::endl;
      cerr << "    " << "path: " << result.path << std::endl;
      cerr << "    " << "query: " << result.query << std::endl;
    }

    return result;
  }

  /// @brief build a new URL struct from a redirect location
  /// @param base_url url that initiated the redirection
  /// @param new_location new location - url or path (absolute or relative)
  /// @return constructed URL struct
  ///
  /// @note https://en.wikipedia.org/wiki/HTTP_location
  URL httpParseRedirect(const URL &base_url, const std::string &new_location) {
    URL new_url = httpParseURL(new_location, true);

    // redirect contains just a path
    // add information from the base_url
    if (new_url.port == 0) {
      new_url.port = base_url.port;
      new_url.scheme = base_url.scheme;
      new_url.host = base_url.host;
      new_url.user = base_url.user;
      // if relative path add new content behind the last '/'
      if (new_url.path.at(0) != '/') {
        int base_split = base_url.path.find_last_of('/');
        new_url.path = base_url.path.substr(0, base_split) + "/" + new_url.path;
      }
    }
    // some full redirects miss the original port, add it back in
    else if ((new_url.port != base_url.port) &&
             (base_url.port != 80 && base_url.port != 443) &&
             (new_url.host == base_url.host) &&
             (new_url.scheme == base_url.scheme)) {
      new_url.port = base_url.port;
    }
    return new_url;
  }

  /// @brief Parsing the output of a raw http response
  /// @param response raw http response
  /// @param output cleaned message body of the http response
  /// @param header map of header information
  /// @return http status code
  ///
  /// @note Message Header definition: https://datatracker.ietf.org/doc/html/rfc7230#section-3.2
  /// @note Chunked Transfer Coding definition: https://www.w3.org/Protocols/rfc2616/rfc2616-sec3.html#sec3.6.1
  int httpParseResponse(std::string response,
                        std::string &output, std::map <std::string, std::string> &header) {

    string soliloquy = XPID + "aurostd::httpParseResponse():";

    std::string header_raw = "";
    int status_code = -1;
    std::string status_line = "";
    std::string const delimiter = "\r\n";
    std::string chunk_length_octet = "";
    size_t chuck_length = 0;
    size_t border = 0;
    size_t split = 0;
    size_t full_request_length = 0;
    size_t offset = delimiter.length();

    // separate headers from message body
    border = response.find(delimiter + delimiter);
    header_raw = response.substr(0, border);
    response.erase(0, border + offset + offset);

    // check if Transfer-Encoding is chunked
    bool isChunked = ( header_raw.find("chunked") != std::string::npos );

    // extract first line that contains the status
    border = header_raw.find(delimiter);
    status_line = header_raw.substr(0, border);
    header_raw.erase(0, border + offset);

    // save status code
    split = status_line.find(' ') + 1;
    border = status_line.find(' ', split + 1);
    status_code = aurostd::string2utype(status_line.substr(split, border - split));

    // extract the header data
    std::string key = "";
    std::string item = "";
    while (!header_raw.empty()) {
      split = header_raw.find(':');
      border = header_raw.find(delimiter);
      if (border != std::string::npos) {
        header_raw.substr(0, split);
        // Header field name are case-insensitive
        key = aurostd::tolower(header_raw.substr(0, split));
        item = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(header_raw.substr(split + 1, border - split - 1));
        header_raw.erase(0, border + offset);

        // Headers that appear multiple times are equal to a list separated by a comma
        if (header.find(key) == header.end()) {
          header.emplace(key, item);
        } else {
          header[key] += "," + item;
        }
      } else { header_raw.clear(); }
    }

    // in HTTP1.1 data is mostly send in chunks
    if (isChunked) {
      while (!response.empty() and chunk_length_octet != "0") {
        border = response.find(delimiter);
        chunk_length_octet = response.substr(0, border);
        response.erase(0, border + offset);
        if (chunk_length_octet != "0") {
          chuck_length = std::stoi(chunk_length_octet, 0, 16);
          full_request_length += chuck_length;
          // sanity check one - every chunk should end with \r\n
          if (response.substr(chuck_length, offset) != delimiter) {
            cerr << soliloquy << " Chunk did not end in `\\r\\n`." << endl;
            return -1;
          }
          output += response.substr(0, chuck_length);
          response.erase(0, chuck_length + offset);
        }
      }
      // sanity check two - the constructed string should have the length reported by the server
      if (full_request_length != output.length()) {
        cerr << soliloquy << " Final string size is different from the chunk size sum." << endl;
        return -1;
      }
      return status_code;
    } else {
      output = response;
      return status_code;
    }
  }


  /// @brief convert characters into their percent representation
  /// @param raw_str sting to escape
  /// @param characters replace just the given characters
  /// @return escaped string
  /// @note https://www.rfc-editor.org/rfc/rfc3986#section-2.1
  /// @note characters to be escaped need to be a single byte long (for utf-8 strings use httpPercentEncodingFull)
  string httpPercentEncodingSelected(const string &raw_str, const string &characters){

    bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XHTTP_);
    string soliloquy = XPID + "aurostd::httpPercentEncodingSelected():";

    const char *reserved=characters.c_str();

    size_t start=0;
    size_t border=0;
    int to_replace=0;

    char * str_position;
    char str[raw_str.length()];
    strcpy(str, raw_str.c_str());
    str_position=std::strpbrk(str, reserved);

    std::stringstream output;
    if (LDEBUG) cerr << soliloquy << " Escaping '" << raw_str << "'" << endl;

    while (str_position!=NULL){
      border = str_position-str;
      to_replace = str[border];
      if (to_replace<0) to_replace+=256;
      output << raw_str.substr(start, border-start) << "%" << std::uppercase << std::hex << std::setfill('0') << std::setw(2) << to_replace;
      if (LDEBUG) cerr << " Match '" << str[border] << "' (%" << std::uppercase << std::hex << std::setfill('0') << to_replace << std::dec << ") at " << border << endl;
      start = border+1;
      str_position=std::strpbrk(str_position+1,reserved);
    }
    output << raw_str.substr(start);
    return output.str();

  }


  /// @brief Fully percent encode a string
  /// @param work_str sting to escape
  /// @return escaped string
  /// @note https://www.rfc-editor.org/rfc/rfc3986#section-2.1
  /// @note just leave unreserved characters (ALPHA / DIGIT / "-" / "." / "_" / "~")
  string httpPercentEncodingFull(string work_str){

    bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XHTTP_);
    string soliloquy = XPID + "aurostd::httpPercentEncoding():";

    const char *allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                          "abcdefghijklmnopqrstuvwxyz"
                          "0123456789"
                          "-_.~";

    size_t border=0;
    int to_replace=0;
    std::stringstream output;
    if (LDEBUG) cerr << soliloquy << " Escaping '" << work_str << "'" << std::endl;

    while (!work_str.empty()){
      border = std::strspn(work_str.c_str(), allowed);
      to_replace = work_str[border];
      if (to_replace<0) to_replace+=256;
      output << work_str.substr(0, border) << "%" << std::uppercase << std::hex << std::setfill('0') << std::setw(2) << to_replace;
      if (LDEBUG) cerr << " Match '" << work_str[border] << "' (%" << std::uppercase << std::hex << std::setfill('0') << to_replace << std::dec << ")" << std::endl; // at " << border << endl;
      work_str.erase(0,border+1);
    }
    return output.str();
  }


  /// @brief Get a raw http response
  /// @param hostname hostname to establish a connection to
  /// @param query query to run on the server (like "/API/aflux/?nspecies(4),paging(1,300)")
  /// @param response raw http response
  /// @param success was the request successful
  ///
  /// @todo implementing timeout
  /// @todo write response into a file
  ///
  /// This function is primitive and is aimed to communicate with AFLOW's API servers.
  bool httpGetResponse(URL url, std::string &response) {

    bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XHTTP_);
    string soliloquy = XPID + "aurostd::httpGetResponse():";

    response.clear();

    char buffer[BUFSIZ];
    char request_template[] = "GET %s HTTP/1.0\r\nHost: %s\r\nUser-Agent: aflow/%s (https://aflow.org)\r\n\r\n";

    struct protoent *protocol_entry;
    struct hostent *host_entry;

    in_addr_t ip_address;
    struct sockaddr_in socket_entry{};
    char ip_address_str[20];
    char *request;
    int request_len;
    int socket_file_descriptor;
    ssize_t nbytes_total, nbytes_last;
    unsigned long loaded_bytes = 0;

    // build request and save its length
    request_len = asprintf(&request, request_template, (url.path + url.query).c_str(), url.host.c_str(), AFLOW_VERSION);


    // get the TCP protocol entry
    protocol_entry = getprotobyname("tcp");
    if (protocol_entry == nullptr) {
      cerr << soliloquy << " Failed to get TCP protocol entry!" << endl;
      return false;
    }

    // open the socket
    socket_file_descriptor = socket(AF_INET, SOCK_STREAM, protocol_entry->p_proto);
    if (socket_file_descriptor == -1) {
      cerr << soliloquy << " Failed to open socket!" << endl;
      return false;
    }

    // find the IP of the host
    host_entry = gethostbyname(url.host.c_str());
    if (host_entry == nullptr) {
      cerr << soliloquy << " Failed to find information on '" << url.host << "'!" << endl;
      return false;
    }

    ip_address = inet_addr(inet_ntoa(*(struct in_addr *) *(host_entry->h_addr_list)));
    if (ip_address == (in_addr_t) - 1) {
      cerr << soliloquy << " Failed to generate IP address for '" << url.host << "'!" << endl;
      return false;
    }
    inet_ntop(AF_INET, &ip_address, ip_address_str, sizeof(ip_address_str));

    socket_entry.sin_addr.s_addr = ip_address;
    socket_entry.sin_family = AF_INET;
    if (url.port > 65535) {
      cerr << soliloquy << " Failed to connect to " << url.host << " as port is out of range (" << url.port << ">65535)" << endl;
      return false;
    }
    socket_entry.sin_port = htons((short) url.port);

    // start connection
    if (connect(socket_file_descriptor, (struct sockaddr *) &socket_entry, sizeof(socket_entry)) == -1) {
      cerr << soliloquy << " Failed to connect to " << url.host << " (" << ip_address_str << ":" << url.port << ")" << endl;
      return false;
    }
    if (LDEBUG)
      cerr << soliloquy << " Connected to " << url.host << " (" << ip_address_str << ":" << url.port << ")" << endl;

    // send the request
    nbytes_total = 0;
    while (nbytes_total < request_len) {
      nbytes_last = write(socket_file_descriptor, request + nbytes_total, request_len - nbytes_total);
      if (nbytes_last == -1) {
        cerr << soliloquy << " Failed to send the request to " << url.host << "( " << ip_address_str << ":" << url.port
             << ")" << endl;
        return false;
      }
      nbytes_total += nbytes_last;
    }
    if (LDEBUG) cerr << soliloquy << " Request written (" << nbytes_total << " bytes)" << endl;

    // read the response
    if (LDEBUG) cerr << soliloquy << " Start reading response:";
    while ((nbytes_total = read(socket_file_descriptor, buffer, BUFSIZ)) > 0) {
      if (LDEBUG) {
        loaded_bytes += nbytes_total;
        cerr << "    read " << loaded_bytes << " bytes \r";
      }
      response += std::string(buffer, nbytes_total);
    }
    if (LDEBUG) cerr << soliloquy << " Finished reading response (" << loaded_bytes << " bytes)" << endl;

    if (nbytes_total == -1) {
      cerr << soliloquy << " Failed to read response from " << url.host << "( " << ip_address_str << ":" << url.port
           << ")" << endl;
      response.clear();
      return false;
    }

    // close the socket
    close(socket_file_descriptor);
    return true;
  }

  void httpGet(URL url,
               std::string &output, int &status_code, std::map <std::string, std::string> &header) {

    bool LDEBUG = (false || XHOST.DEBUG || _DEBUG_XHTTP_);
    string soliloquy = XPID + "aurostd::httpGet():";

    if (LDEBUG) cerr << soliloquy << " GET '" << url.path << "' from " << url.host << endl;

    output.clear();
    header.clear();
    int redirect_codes[] = {301, 302, 307, 308};
    std::string response = "";
    const int max_number_of_tries = 3;
    int number_of_tries = 0;
    bool success = false;

    while (status_code != 200 and number_of_tries < max_number_of_tries) {
      // read server response
      success = httpGetResponse(url, response);
      number_of_tries += 1;
      if (success) {
        // parse the response to extract header and status code
        status_code = httpParseResponse(response, output, header);
        if (LDEBUG) cerr << soliloquy << " Received status code " << status_code << endl;
        if (LDEBUG && status_code!=-1) {
          cerr << soliloquy << " Header content: " << endl;
          for (auto it = header.begin(); it != header.end(); it++)
            cerr << "    " << it->first << ": " << it->second << std::endl;
        }
        if (std::find(std::begin(redirect_codes), std::end(redirect_codes), status_code) != std::end(redirect_codes)) {
          if (header.find("location")!=header.end()) {
            url = httpParseRedirect(url, header["location"]);
            if (LDEBUG) cerr << soliloquy << " Redirected to '" << header["location"] << "'" << endl;
          }
          else {
            cerr << soliloquy << " GET request failed due to a redirect without location." << endl;
            return;
          }
        }
      } else if (number_of_tries < max_number_of_tries) {
        cerr << soliloquy << " Retry failed GET in 5s (" << number_of_tries << " of " << max_number_of_tries << ")" << endl;
        aurostd::sleep(5);
        continue;
      } else {
        cerr << soliloquy << " GET request failed after " << number_of_tries << " tries." << endl;
        return;
      }
    }
    if (LDEBUG) {
      if (number_of_tries > 1) cerr << soliloquy << " GET request done after " << number_of_tries << " tries." << endl;
      else cerr << soliloquy << " GET request done" << endl;
    }
  }

  // Start of external callable functions


  /// @brief Retrieve data from an url string
  /// @param url_str content url
  /// @return HTTP status code (-1 on failure)
  int httpGetStatus(const std::string &url_str) {
    URL url = httpParseURL(url_str);
    std::string output = "";
    int status_code = -1;
    std::map <std::string, std::string> header;
    httpGet(url, output, status_code, header);
    return status_code;
  }

  /// @brief Retrieve data from an url string
  /// @param url_str content url
  /// @param output message body
  /// @return HTTP status code (-1 on failure)
  int httpGetStatus(const std::string &url_str, std::string &output) {
    URL url = httpParseURL(url_str);
    int status_code = -1;
    std::map <std::string, std::string> header;
    httpGet(url, output, status_code, header);
    return status_code;
  }

  /// @brief Retrieve data from an url string
  /// @param url_str content url
  /// @param output message body
  /// @param header response header
  /// @return HTTP status code (-1 on failure)
  int httpGetStatus(const std::string &url_str, std::string &output, std::map <std::string, std::string> &header) {
    URL url = httpParseURL(url_str);
    int status_code = -1;
    httpGet(url, output, status_code, header);
    return status_code;
  }

  /// @brief Retrieve data from an url string
  /// @param host server name or IP to contact
  /// @param query GET query
  /// @param output message body
  /// @return HTTP status code (-1 on failure)
  int httpGetStatus(const std::string &host, const std::string &path, const std::string &query, std::string &output) {
    URL url = httpConstructURL(host, path, query);
    int status_code = -1;
    std::map <std::string, std::string> header;
    httpGet(url, output, status_code, header);
    return status_code;
  }

  /// @brief Retrieve data from an url string
  /// @param host server name or IP to contact
  /// @param query GET query
  /// @param output message body
  /// @param header response header
  /// @return HTTP status code (-1 on failure)
  int httpGetStatus(const std::string &host, const std::string &path, const std::string &query, std::string &output, std::map <std::string, std::string> &header) {
    URL url = httpConstructURL(host, path, query);
    int status_code = -1;
    httpGet(url, output, status_code, header);
    return status_code;
  }

  /// @brief Retrieve data from an url string
  /// @param url_str content url
  /// @return message body
  std::string httpGet(const std::string &url_str) {
    URL url = httpParseURL(url_str);
    std::string output = "";
    int status_code = -1;
    std::map <std::string, std::string> header;

    httpGet(url, output, status_code, header);
    return output;
  }

  /// @brief Retrieve data from an url string
  /// @param url_str content url
  /// @param status_code HTTP status code (-1 on failure)
  /// @return message body
  std::string httpGet(const std::string &url_str, int &status_code) {
    URL url = httpParseURL(url_str);
    std::string output = "";
    std::map <std::string, std::string> header;

    status_code = -1;
    httpGet(url, output, status_code, header);
    return output;
  }

  /// @brief Retrieve data from an url string
  /// @param url_str content url
  /// @param status_code HTTP status code (-1 on failure)
  /// @param header response header
  /// @return message body
  std::string httpGet(const std::string &url_str, int &status_code, std::map <std::string, std::string> &header) {
    URL url = httpParseURL(url_str);

    std::string output="";
    status_code = -1;
    httpGet(url, output, status_code, header);
    return output;
  }

  /// @brief Retrieve data for host + query
  /// @param host server name or IP to contact
  /// @param query GET query
  /// @return message body
  std::string httpGet(const std::string &host, const std::string &path, const std::string &query) {
    URL url = httpConstructURL(host, path, query);
    std::string output="";
    int status_code = -1;
    std::map <std::string, std::string> header;

    httpGet(url, output, status_code, header);
    return output;
  }

  /// @brief Retrieve data for host + query
  /// @param host server name or IP to contact
  /// @param query GET query
  /// @param status_code HTTP status code (-1 on failure)
  /// @return message body
  std::string httpGet(const std::string &host, const std::string &path, const std::string &query, int &status_code) {
    URL url = httpConstructURL(host, path, query);
    std::string output="";
    std::map <std::string, std::string> header;

    status_code = -1;
    httpGet(url, output, status_code, header);
    return output;
  }

  /// @brief Retrieve data for host + query
  /// @param host server name or IP to contact
  /// @param query GET query
  /// @param status_code HTTP status code (-1 on failure)
  /// @param header response header
  /// @return message body
  std::string httpGet(const std::string &host, const std::string &path, const std::string &query, int &status_code, std::map <std::string, std::string> &header) {
    URL url = httpConstructURL(host, path, query);
    std::string output="";

    status_code = -1;
    httpGet(url, output, status_code, header);
    return output;
  }

}

#endif  // _AUROSTD_XHTTP_CPP_
