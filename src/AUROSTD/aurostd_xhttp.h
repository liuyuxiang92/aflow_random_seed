// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *           Aflow HAGEN ECKERT - Duke University 2021-2022                *
// *                                                                         *
// ***************************************************************************
// Written by Hagen Eckert
// hagen.eckert@duke.edu

#ifndef AFLOW_SRC_AUROSTD_XHTTP_H
#define AFLOW_SRC_AUROSTD_XHTTP_H

namespace aurostd{
  struct URL {
    std::string scheme;
    std::string user;
    std::string host;
    unsigned short port;
    std::string path;
    std::string query;
  };

  int httpGetStatus(const std::string &url, std::string &output);
  int httpGetStatus(const std::string &url, std::string &output, std::map<std::string, std::string> &header);

  int httpGetStatus(const std::string &host, const std::string &path, const std::string &query, std::string &output);
  int httpGetStatus(const std::string &host, const std::string &path, const std::string &query, std::string &output, std::map<std::string, std::string> &header);

  std::string httpGet(const std::string &url);
  std::string httpGet(const std::string &url, int &status_code);
  std::string httpGet(const std::string &url, int &status_code, std::map<std::string, std::string> &header);

  std::string httpGet(const std::string &host, const std::string &path, const std::string &query);
  std::string httpGet(const std::string &host, const std::string &path, const std::string &query, int &status_code);
  std::string httpGet(const std::string &host, const std::string &path, const std::string &query, int &status_code, std::map<std::string, std::string> &header);

  std::string httpPercentEncodingSelected(const std::string &raw_str, const std::string &characters);
  std::string httpPercentEncodingFull(std::string work_str);

  URL httpParseURL(const std::string &url, const bool strict = false);
}

#endif //AFLOW_SRC_AUROSTD_XHTTP_H
