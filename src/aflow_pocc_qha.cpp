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
#define PRECISION 10

namespace pocc {
  void writeQHAdatablock(stringstream &output, vector<vector<double> > data,
      string block_name)
  {
    string blockname = "[" + block_name + "]";

    output << AFLOWIN_SEPARATION_LINE << std::endl;
    //file << blockname + "SYSTEM=" << system_title << std::endl;
    output << blockname + "START" << std::endl;

    // print header
    output << setw(5)  << "#T[K]"  << setw(SW) << ' ' <<
      setw(TW) << "V[A^3/atom]"          << setw(SW) << ' ' <<
      setw(TW) << "F(V)[eV/atom]"        << setw(SW) << ' ' <<
      setw(TW) << "B[GPa]"               << setw(SW) << ' ' <<
      setw(TW) << "beta[10^-5/K]"        << setw(SW) << ' ' <<
      setw(TW) << "Cv(V)[kB/atom]"       << setw(SW) << ' ' <<
      setw(TW) << "Cp(V)[kB/atom]"       << setw(SW) << ' ' <<
      setw(TW) << "gamma"                << setw(SW) << ' ' <<
      setw(TW) << "Bprime";
    output << std::endl;

    // print data
    for (uint row=0; row<data.size(); row++){
      // special formatting for temperature
      output.unsetf(ios_base::floatfield);
      output << setw(5) << data[row][0] << setw(SW) << ' '
             << std::fixed << std::setprecision(PRECISION);
      // output other properties
      for (uint col=1; col<data[row].size(); col++){
        output << setw(TW) << data[row][col] << setw(SW) << ' ';
      }
      output << std::endl;
    }

    // print footer
    output << blockname + "STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;
  }

  /// Calculates the thermal expansion coefficient as a logarithmic derivative
  /// of equilibrium volume employing the Savitzky-Golay filter for the
  /// differentiation.
  xvector<double> calcThermalExpansionSGloc(const xvector<double> &volumes, double dT)
  {
    // Convolution weights for the Savitzky-Golay 5pt cubic filter as reported in:
    // "General least-squares smoothing and differentiation by the convolution (Savitzky-Golay) method"
    // Peter A. Gorry Analytical Chemistry 1990 62 (6), 570-573
    // https://doi.org/10.1021/ac00205a007
    const static xmatrix<double> SGmat(5,5);
    SGmat[1][1]=-125.0/84.0; SGmat[1][2]=-19.0/42.0; SGmat[1][3]= 1.0/12.0; SGmat[1][4]=  5.0/42.0; SGmat[1][5]= -29.0/84.0;
    SGmat[2][1]= 136.0/84.0; SGmat[2][2]= -1.0/42.0; SGmat[2][3]=-8.0/12.0; SGmat[2][4]=-13.0/42.0; SGmat[2][5]=  88.0/84.0;
    SGmat[3][1]=  48.0/84.0; SGmat[3][2]= 12.0/42.0; SGmat[3][3]= 0.0/12.0; SGmat[3][4]=-12.0/42.0; SGmat[3][5]= -48.0/84.0;
    SGmat[4][1]= -88.0/84.0; SGmat[4][2]= 13.0/42.0; SGmat[4][3]= 8.0/12.0; SGmat[4][4]=  1.0/42.0; SGmat[4][5]=-136.0/84.0;
    SGmat[5][1]=  29.0/84.0; SGmat[5][2]= -5.0/42.0; SGmat[5][3]=-1.0/12.0; SGmat[5][4]= 19.0/42.0; SGmat[5][5]= 125.0/84.0;
    const static xvector<double> SGvec(5);
    SGvec[1]= 1.0/12.0; SGvec[2]=-8.0/12.0; SGvec[3]= 0.0/12.0; SGvec[4]= 8.0/12.0; SGvec[5]=-1.0/12.0;
    ////////////////////////////////////////////////////////////////////////////////

    string function = "calcThermalExpansionSG():", msg = "";
    int npoints = volumes.rows;
    if (npoints<5){
      msg = "Savitzky-Golay filter requires at least 5 points: only ";
      msg += aurostd::utype2string(npoints) + " were provided.";
      msg += " Make sure that the minimum range among the set of POCC structures,";
      msg += "defined by the [AFLOW_APL]TPT parameter for each, is reasonable.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg, _INDEX_ILLEGAL_);
    }

    xvector<double> endpoints(5), dummy(5);
    xvector<double> beta(npoints);

    // calculate derivatives for the first 2 points
    for (int i=1; i<=5; i++) dummy[i] = volumes[i];
    endpoints = dummy*SGmat;
    for (int i=1; i<=2; i++) beta[i] = endpoints[i];

    // calculate derivatives for the [3:end-3] points
    int id = 0;
    for (int i=3; i<=npoints-2; i++){
      beta[i] = 0.0;
      for (int j=1; j<=5; j++){
        id = i - 3 + j;
        beta[i] += SGvec[j]*volumes[id];
      }
    }

    // calculate derivatives for the last 2 points
    for (int i=1; i<=5; i++) dummy[i] = volumes[npoints-5+i];
    endpoints = dummy*SGmat;
    for (int i=4; i<=5; i++) beta[npoints-5+i] = endpoints[i];

    // calculate the coefficient of thermal expansion
    for (int i=1; i<=npoints; i++) beta[i] /= (volumes[i]*dT);

    return beta;
  }

  /// Calcualates POCC-average of QHA-related properties.
  void POccCalculator::calculateQHAPropertiesAVG(const vector<double>& v_temperatures) {
    string function = XPID +  "POccCalculator::calculateQHAProperties():";
    string msg = "", filename = "";

    bool LDEBUG = false || _DEBUG_POCC_QHA_ || XHOST.DEBUG;

    msg = "Performing POCC+QHA post-processing step.";
    pflow::logger(_AFLOW_FILE_NAME_, function, msg, m_aflags.Directory,
      *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

////////////////////////////////////////////////////////////////////////////////
    apl::QHAmethod qha_method = apl::QHA_CALC;
    vector<string> dirs2ignore;

    // check if there are unstable structures we want to ignore
    for (uint i=0; i<m_ARUN_directories.size(); i++){
      filename = m_aflags.Directory + "/" + m_ARUN_directories[i] + "/";
      filename += DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_IMAG_FILE;

      if (apl::hasImaginary(filename, apl::QHAmethod2label(qha_method))){
        dirs2ignore.push_back(m_ARUN_directories[i]);
      }
    }

    // if there are: ignore and reload POCC data
    if (dirs2ignore.size()){
      XHOST.vflag_control.flag("ARUNS2SKIP", true);
      XHOST.vflag_control.addattachedscheme("ARUNS2SKIP",
           aurostd::joinWDelimiter(dirs2ignore, ","), true);

      loadDataIntoCalculator();
      setDFTEnergies();
    }
////////////////////////////////////////////////////////////////////////////////

    vector<vector<vector<double> > > pocc_qha_thermo_properties;
    vector<double> line;
    uint ncols = 0;
    string data;
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
          pocc_qha_thermo_properties.push_back(thermo_properties);
        }
        else{
          msg = "The " + filename + " file is empty.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg, _FILE_ERROR_);
        }
      }
      else{
        msg = "The " + filename + " file is missing.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg, _FILE_NOT_FOUND_);
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
      degeneracies.push_back((*it).getDegeneracy());
      n_total += degeneracies.back();
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

      vector<vector<double> > averaged_qha_data(minsize, vector<double> (ncols));
      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
           it != l_supercell_sets.end(); ++it)
      {
        isupercell=std::distance(l_supercell_sets.begin(),it);
        for (uint row=0; row<minsize; row++){
          for (uint col=0; col<ncols; col++){
            averaged_qha_data[row][col] += (*it).m_probability *
              pocc_qha_thermo_properties[isupercell][row][col];
          }
        }
      }
      writeQHAdatablock(file, averaged_qha_data, "POCC_QHA_SJ_THERMO_T=" +
          aurostd::utype2string(T));
    }

    // average for Tpocc = T using the DFT energies to calculate the probabilities
    vector<vector<double> > averaged_qha_data_T(minsize, vector<double> (ncols));
    for (uint row=0; row<minsize; row++){
      T = pocc_qha_thermo_properties[0][row][0];
      setPOccStructureProbabilities(T);

      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
           it != l_supercell_sets.end(); ++it)
      {
        isupercell=std::distance(l_supercell_sets.begin(),it);
        for (uint col=0; col<ncols; col++){
          averaged_qha_data_T[row][col] += (*it).m_probability *
                  pocc_qha_thermo_properties[isupercell][row][col];
        }
      }
    }
    writeQHAdatablock(file, averaged_qha_data_T, "POCC_QHA_SJ_THERMO_E_T=T");

    // average for Tpocc = T using the free energies to calculate the probabilities
    vector<vector<double> > averaged_qha_data_F_T(minsize, vector<double> (ncols));
    for (uint row=0; row<minsize; row++){
      m_energy_dft_ground = AUROSTD_MAX_DOUBLE;
      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
           it != l_supercell_sets.end(); ++it)
      {
        isupercell=std::distance(l_supercell_sets.begin(),it);
        (*it).m_energy_dft = pocc_qha_thermo_properties[isupercell][row][2];
        m_energy_dft_ground = std::min(m_energy_dft_ground, (*it).m_energy_dft);
      }

      T = pocc_qha_thermo_properties[0][row][0];
      setPOccStructureProbabilities(T);

      for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
           it != l_supercell_sets.end(); ++it)
      {
        isupercell=std::distance(l_supercell_sets.begin(),it);
        for (uint col=0; col<ncols; col++){
          averaged_qha_data_F_T[row][col] += (*it).m_probability *
                  pocc_qha_thermo_properties[isupercell][row][col];
        }
      }
    }
    writeQHAdatablock(file, averaged_qha_data_F_T, "POCC_QHA_SJ_THERMO_F_T=T");

    // calculate the thermal expansion as a logarithmic derivative of volume
    double dT = (pocc_qha_thermo_properties[0][minsize-1][0] -
                 pocc_qha_thermo_properties[0][0][0])/(minsize - 1);

    xvector<double> vols(minsize);
    for (uint row=0; row<minsize; row++){
      vols[row+1] = averaged_qha_data_T[row][1];
    }
    xvector<double> beta = calcThermalExpansionSGloc(vols, dT);
    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[POCC_QHA_SJ_BETA_E]START" << std::endl;
    file << setw(5)  << "#T[K]"         << setw(SW) << ' ' <<
            setw(TW) << "beta[10^-5/K]" << std::endl;
    for (uint row=0; row<minsize; row++){
      file.unsetf(ios_base::floatfield);
      file << setw(5)    << averaged_qha_data_T[row][0] << setw(SW) << ' '
           << std::fixed << std::setprecision(PRECISION)
           << setw(TW)   << beta[row+1]*1e5
           << std::endl;
    }
    file << "[POCC_QHA_SJ_BETA_E]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    for (uint row=0; row<minsize; row++){
      vols[row+1] = averaged_qha_data_F_T[row][1];
    }
    beta = calcThermalExpansionSGloc(vols, dT);
    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << "[POCC_QHA_SJ_BETA_F]START" << std::endl;
    file << setw(5)  << "#T[K]"         << setw(SW) << ' ' <<
            setw(TW) << "beta[10^-5/K]" << std::endl;
    for (uint row=0; row<minsize; row++){
      file.unsetf(ios_base::floatfield);
      file << setw(5)    << averaged_qha_data_F_T[row][0] << setw(SW) << ' '
           << std::fixed << std::setprecision(PRECISION)
           << setw(TW)   << beta[row+1]*1e5
           << std::endl;
    }
    file << "[POCC_QHA_SJ_BETA_F]STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    string output_file = "aflow.pocc.qha.avgthermo.out";
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
