//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *           Aflow Andriy Smolyanyuk - Duke University 2021-2021           *
// *                                                                         *
//****************************************************************************
// Written by Andriy Smolyanyuk, 2021.
//
// This file provides a framework to calculate thermal properties for
// disordered materials modeled using the POCC + QHA methodology.

#ifndef _AFLOW_POCC_QHA_CPP_
#define _AFLOW_POCC_QHA_CPP_

#include <limits.h>
#include "aflow.h"
#include "aflow_pocc.h"

#define _DEBUG_POCC_QHA_ false

#define SW 5  // width of columns with blank space separator
#define TW 15 // width of columns containing label/number

namespace pocc {
  void writeQHAdatablock(stringstream &output, vector<vector<double> > data,
      string block_name)
  {
    string blockname = "[" + block_name + "]";

    output.precision(10);
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    //file << blockname + "SYSTEM=" << system_title << std::endl;
    output << blockname + "START" << std::endl;

    // print header
    output << setw(5)  << "#T[K]"                          << setw(SW) << ' ' <<
      setw(TW) << "V[A^3/atom]"                            << setw(SW) << ' ' <<
      setw(TW) << "F(V)[eV/atom]"                          << setw(SW) << ' ' <<
      setw(TW) << "B[GPa]"                                 << setw(SW) << ' ' <<
      setw(TW) << "beta[10^-5/K]"                          << setw(SW) << ' ' <<
      setw(TW) << "Cv(V)[kB/atom]"                         << setw(SW) << ' ' <<
      setw(TW) << "Cp(V)[kB/atom]"                         << setw(SW) << ' ' <<
      setw(TW) << "gamma(beta,B,Cv(V))"                    << setw(SW) << ' ' <<
      setw(TW) << "Bprime"                                 << setw(SW) << ' ' <<
      setw(TW) << "gamma(V,mesh)"                          << setw(SW) << ' ' <<
      setw(TW) << "beta(gamma(V,mesh),B,Cv(V))[10^-5/K]"   << setw(SW) << ' ' <<
      setw(TW) << "Cv(V,mesh)[kB/atom]"                    << setw(SW) << ' ' <<
      setw(TW) << "Cp(V,mesh)[kB/atom]"                    << setw(SW) << ' ' <<
      setw(TW) << "gamma(V0,mesh)"                         << setw(SW) << ' ' <<
      setw(TW) << "beta(gamma(V0,mesh),B,Cv(V0))[10^-5/K]" << setw(SW) << ' ' <<
      setw(TW) << "Cv(V0,mesh)[kB/atom]"                   << setw(SW) << ' ' <<
      setw(TW) << "Cp(V0,mesh)[kB/atom]";
    output << std::endl;

    // print data
    for (uint row=0; row<data.size(); row++){
      for (uint col=0; col<data[row].size(); col++){
        output << setw(TW) << data[row][col] << setw(SW) << ' ';
      }
      output << std::endl;
    }

    // print footer
    output << blockname + "STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
  }

  /// Calcualates POCC-average of QHA-related properties.
  void POccCalculator::calculateQHAProperties(const vector<double>& v_temperatures) {
    string function = XPID +  "POccCalculator::calculateQHAProperties():";
    string msg = "", filename = "";

    bool LDEBUG = false || _DEBUG_POCC_QHA_ || XHOST.DEBUG;

    msg = "Performing POCC+QHA post-processing step.";
    pflow::logger(_AFLOW_FILE_NAME_, function, msg, m_aflags.Directory,
      *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    vector<vector<vector<double> > > pocc_qha_thermo_properties;
    vector<double> line;
    uint ncols = 0;
    string data;
    vector<string> aruns;
    vector<bool> arun_contains_data(m_ARUN_directories.size());
    for (uint i=0; i<m_ARUN_directories.size(); i++){
      data = "";

      filename = m_aflags.Directory + "/" + m_ARUN_directories[i] + "/";
      filename += DEFAULT_QHA_FILE_PREFIX + DEFAULT_QHA_THERMO_FILE;
      if (LDEBUG) cerr << function << " file: " << filename << endl;

      if (aurostd::EFileExist(filename)){
         aurostd::ExtractToStringEXPLICIT(aurostd::efile2string(filename), data,
          "[QHA_SJ_THERMO]START", "[QHA_SJ_THERMO]STOP");

        vector<string> lines = aurostd::string2vectorstring(data);
        if (LDEBUG) cerr << function << " datasize: " <<  lines.size() << endl;
        if (lines.size() > 1){
          lines.erase(lines.begin()); // remove the first line containing properties labels

          vector<vector<double> > thermo_properties;
          if (i==0){
            aurostd::string2tokens(lines[0], line);
            ncols = line.size();
          }

          for (uint j=0; j<lines.size(); j++){
            aurostd::string2tokens(lines[j], line);

            if (ncols==line.size()){
              thermo_properties.push_back(line);
              ncols = line.size();
            }
            else{
              msg="Inconsistent amount of properties (columns) among different";
              msg+="POCC-QHA calculations.";
              throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg,
                  _INDEX_MISMATCH_);
            }
          }
          arun_contains_data[i] = true;
          pocc_qha_thermo_properties.push_back(thermo_properties);
          aruns.push_back(m_ARUN_directories[i]);
        }
        else{
          arun_contains_data[i] = false;

          msg = "The file " + filename + " is empty. The possible reason is ";
          msg += "that the structure is dynamically unstable and QHA refused ";
          msg += "to proceed with the calculation. ";
          msg += "Please check if that's the case.";
          pflow::logger(_AFLOW_FILE_NAME_, function, msg, m_aflags.Directory,
            *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
        }
      }
      else{
        arun_contains_data[i] = false;

        msg = "The file " + filename + " is missing. The possible reason is ";
        msg += "that the structure is dynamically unstable and QHA refused ";
        msg += "to proceed with the calculation. ";
        msg += "Please check if that's the case.";
        pflow::logger(_AFLOW_FILE_NAME_, function, msg, m_aflags.Directory,
            *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      }
    }

    if (!ncols) {
      msg = "POCC+QHA was not able to extract any data.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_CORRUPT_);
    }

    if (LDEBUG){
      cerr << function << " ncols: " << ncols << endl;
      cerr << function << " pocc_qha size: " << pocc_qha_thermo_properties.size();
      cerr << endl;
    }

    msg = "The data from the following ARUN directories will be used:\n";
    for (uint i=0; i<aruns.size(); i++) msg += aruns[i] + '\n';
    pflow::logger(_AFLOW_FILE_NAME_, function, msg, m_aflags.Directory,
        *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    uint n = pocc_qha_thermo_properties.size();
    uint minsize = pocc_qha_thermo_properties[0].size();
    for (uint i=1; i<n; i++){
      if (pocc_qha_thermo_properties[i].size() < minsize)
        minsize = pocc_qha_thermo_properties[i].size();

      if (minsize == 0) minsize = pocc_qha_thermo_properties[i].size();
    }

    for (uint i=0; i<n; i++){
      if (pocc_qha_thermo_properties[i].size()){
        pocc_qha_thermo_properties[i].resize(minsize);
      }
      else{
        pocc_qha_thermo_properties[i] = vector<vector<double> >(minsize,
          vector<double> (ncols));
      }
    }

    // check that temperatures are the same for all data blocks
    for (uint i=0; i<n-1; i++){
      for (uint row=0; row<minsize; row++){
        if (!aurostd::isequal(pocc_qha_thermo_properties[i][row][0],
                              pocc_qha_thermo_properties[i+1][row][0])){
          msg="Inconsistent list of temperatures among different";
          msg+="POCC-QHA calculations.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg,
              _VALUE_ILLEGAL_);
        }
      }
    }

    // collect degeneracies only for stable structures
    vector<int> degeneracies;
    int n_total = 0;
    unsigned long long int isupercell=0;
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
         it != l_supercell_sets.end(); ++it)
    {
      isupercell=std::distance(l_supercell_sets.begin(),it);
      if (arun_contains_data[isupercell]){
        degeneracies.push_back((*it).getDegeneracy());
        n_total += degeneracies.back();
      }
    }

    // average for T_pocc=inf
    vector<vector<double> > averaged_data(minsize, vector<double> (ncols));
    for (uint i=0; i<n; i++){
      for (uint row=0; row<minsize; row++){
        for (uint col=0; col<ncols; col++){
          averaged_data[row][col] += degeneracies[i]*pocc_qha_thermo_properties[i][row][col]/n_total;
        }
      }
    }

    // output averaged data
    stringstream file;
    writeQHAdatablock(file, averaged_data, "POCC_QHA_SJ_THERMO_T=INFTY");

    // average for const Tpocc, where Tpocc is a given (sintering) temperature
    double T = 0;
    for (uint i=0; i<v_temperatures.size(); i++){
      T = v_temperatures[i];
      setPOccStructureProbabilities(T);

      vector<double> probabilities;
      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
           it != l_supercell_sets.end(); ++it)
      {
        isupercell=std::distance(l_supercell_sets.begin(),it);
        if (arun_contains_data[isupercell]){
          probabilities.push_back((*it).m_probability);
        }
      }

      double total_prob = 0;
      for (uint j=0; j<probabilities.size(); j++) total_prob += probabilities[j];
      for (uint j=0; j<probabilities.size(); j++) probabilities[j]/=total_prob;

      cout << ": " << T << aurostd::vector2xvector(probabilities) << endl;

      vector<vector<double> > averaged_qha_data(minsize, vector<double> (ncols));
      for (uint i=0; i<n; i++){
        for (uint row=0; row<minsize; row++){
          for (uint col=0; col<ncols; col++){
            averaged_qha_data[row][col] += probabilities[i]*
              pocc_qha_thermo_properties[i][row][col];
          }
        }
      }

      writeQHAdatablock(file, averaged_qha_data, "POCC_QHA_SJ_THERMO_T=" +
          aurostd::utype2string(T));
    }

    // average for Tpocc = T
    vector<vector<double> > averaged_qha_data_T(minsize, vector<double> (ncols));
    for (uint row=0; row<minsize; row++){
      T = pocc_qha_thermo_properties[0][row][0];
      setPOccStructureProbabilities(T);

      vector<double> probabilities;

      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
           it != l_supercell_sets.end(); ++it)
      {
        isupercell=std::distance(l_supercell_sets.begin(),it);
        if (arun_contains_data[isupercell]){
          probabilities.push_back((*it).m_probability);
        }
      }

      double total_prob = 0;
      for (uint j=0; j<probabilities.size(); j++) total_prob += probabilities[j];
      for (uint j=0; j<probabilities.size(); j++) probabilities[j]/=total_prob;

      cout << ":: " << T << aurostd::vector2xvector(probabilities) << endl;

      for (uint i=0; i<n; i++){
        for (uint col=0; col<ncols; col++){
          averaged_qha_data_T[row][col] += probabilities[i]*
                  pocc_qha_thermo_properties[i][row][col];
        }
      }
    }

    writeQHAdatablock(file, averaged_qha_data_T, "POCC_QHA_SJ_THERMO_T=T");

    string output_file = "aflow.pocc.qha.thermo.out";
    msg = "Writing the averaged POCC+QHA data to " + output_file+" file.";
    pflow::logger(_AFLOW_FILE_NAME_, function, msg, m_aflags.Directory,
        *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);
    output_file = m_aflags.Directory + "/" + output_file;

    if (!aurostd::stringstream2file(file, output_file)){
      msg = "Error writing to " + output_file + " file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
    }
  }
}

#endif
