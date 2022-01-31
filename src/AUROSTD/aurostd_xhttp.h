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
  int httpGet(const std::string &url, std::string &output);
  int httpGet(const std::string &url, std::string &output, std::map<std::string, std::string> &header);

  int httpGet(const std::string &host, const std::string &query, std::string &output);
  int httpGet(const std::string &host, const std::string &query, std::string &output, std::map<std::string, std::string> &header);

  std::string httpGet(const std::string &url);
  std::string httpGet(const std::string &url, int &status_code);
  std::string httpGet(const std::string &url, int &status_code, std::map<std::string, std::string> &header);

  std::string httpGet(const std::string &host, const std::string &query);
  std::string httpGet(const std::string &host, const std::string &query, int &status_code);
  std::string httpGet(const std::string &host, const std::string &query, int &status_code, std::map<std::string, std::string> &header);

  std::string httpPercentEncoding(const std::string &raw_str, const std::string &characters="");
}

#endif //AFLOW_SRC_AUROSTD_XHTTP_H
