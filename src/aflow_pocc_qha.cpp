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

namespace pocc {
  void POccCalculator::calculateQHAProperties(const vector<double>& v_temperatures) {
    cout << "QHA and POCC" << endl;
    double T = 0;
    for (uint i=0; i<v_temperatures.size(); i++){
      T = v_temperatures[i];
      setPOccStructureProbabilities(T);
      cout << "T= " << T << " ";

      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); 
        it != l_supercell_sets.end(); ++it) {
        cout << (*it).m_probability << " ";
      }
      cout << endl;
    }

    double LDEBUG = true;

    string filename;
    vector<vector<vector<double> > > pocc_qha_thermo_properties;
    vector<double> line;
    uint ncols = 0;
    string data;
    for (uint i=0; i<m_ARUN_directories.size(); i++){
      data = "";
      cout << "=== " << i << " ===" << endl;

      filename = m_aflags.Directory + "/" + m_ARUN_directories[i] + "/";
      filename += DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_THERMO_FILE;
      if (LDEBUG) cerr << filename << endl;

      if (aurostd::EFileExist(filename)){
         aurostd::ExtractToStringEXPLICIT(aurostd::efile2string(filename), data, 
          "[QHA_SJ_THERMO]START", "[QHA_SJ_THERMO]STOP");

        vector<string> lines = aurostd::string2vectorstring(data);
        if (LDEBUG) cerr << "datasize: " <<  lines.size() << endl;
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
              cerr << "Inconsistent amount of properties" << endl;
              // TODO: throw
            }
          }
          pocc_qha_thermo_properties.push_back(thermo_properties);
        }
        else{
          pocc_qha_thermo_properties.push_back(vector<vector<double> >());
          cout << "Empty file. Dynamically unstable structure?" << endl;
        }
      }
      else{
        pocc_qha_thermo_properties.push_back(vector<vector<double> >());
        cout << "Missing file. Dynamically unstable structure?" << endl;
      }
    }

    if (!ncols) {
      cerr << "No properties!" << endl;
      //TODO: throw; 
    }

    cout << "ncols: " << ncols << endl;
    cout << "size: " << pocc_qha_thermo_properties.size() << endl;
    for (uint i=0; i<pocc_qha_thermo_properties.size(); i++){
      cout << pocc_qha_thermo_properties[i].size() << " ";
//      cout << pocc_qha_thermo_properties[i][0].size() << endl;
    }

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

//    for (uint i=0; i<n; i++){
//      cout << "================================================================================" << endl;
//      for (uint j=0; j<minsize; j++){
//        for (uint k=0; k<ncols; k++){
//          cout << pocc_qha_thermo_properties[i][j][k] << " ";
//        }
//        cout << endl;
//      }
//      cout << endl;
//    }

    if (LDEBUG) cerr << "minsize = " << minsize << " " << "ncols = " << ncols << endl;

    cout << endl << endl;

    // complete random
    vector<vector<double> > averaged_qha_data(minsize, vector<double> (ncols));
    unsigned long long int isupercell=0;
    for (uint i=0; i<minsize; i++){
      setPOccStructureProbabilities(100000);

      averaged_qha_data[i][0] = pocc_qha_thermo_properties[0][i][0];
      for (uint j=1; j<ncols; j++){
        for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); 
          it != l_supercell_sets.end(); ++it) 
        {
          isupercell=std::distance(l_supercell_sets.begin(),it);

//          averaged_qha_data[i][j] += (*it).m_probability *
          averaged_qha_data[i][j] += 1.0/6.0 *
            pocc_qha_thermo_properties[isupercell][i][j];
        }
      }
    }

    for (uint i=0; i<minsize; i++){
      for (uint j=0; j<ncols; j++){
        cout << averaged_qha_data[i][j] << " ";
      }
      cout << endl;
    }

    // plot probabilites
    setPOccStructureProbabilities(100000);
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin(); 
          it != l_supercell_sets.end(); ++it) 
    {
      cout << (*it).m_probability  << " ";
    }
    cout << endl;
    cout << "1/7 = " << 1.0/7.0 << endl;

//////////////////////////////////////////////////////////////////////////////////////
  }
}

#endif
