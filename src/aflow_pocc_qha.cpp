//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *           Aflow Andriy Smolyanyuk - Duke University 2021-2021           *
// *                                                                         *
//****************************************************************************
// Written by Andriy Smolyanyuk, 2021.
//
// This file provides a framework to calculate thermodynamic properties for
// disordered materials modeled using the POCC + QHA methodology.

#ifndef _AFLOW_POCC_QHA_CPP_
#define _AFLOW_POCC_QHA_CPP_

#include <cfloat>
#include "aflow.h"
#include "aflow_pocc.h"

#define _DEBUG_POCC_QHA_ false

#define SW 5  // width of columns with blank space separator
#define TW 15 // width of columns containing label/number

namespace pocc {
  /// Reads fitting coefficients from the filename file.
  /// QHA and EOS methods are defined by blockname.
  bool readCoeffData(const string& filename, const string& blockname,
      xvector<double> &T, xmatrix<double> &coeffs)
  {
    string function = "readCoeffData():", msg = "";
    string block = "[" + blockname + "_COEFF]";
    string data;

    bool LDEBUG = false || _DEBUG_POCC_QHA_ || XHOST.DEBUG;

    // extract data
    aurostd::ExtractToStringEXPLICIT(aurostd::efile2string(filename), data,
        block + "START", block + "STOP");
    vector<string> v_data = aurostd::string2vectorstring(data);

    // determine number of columns in the file and do some consistency check
    uint ncols = 0;
    uint nrows = v_data.size();
    vector<double> tokens;
    if (nrows){
      aurostd::string2tokens(v_data[0], tokens);
      ncols = tokens.size();
      if (ncols <= 1){
        msg = "Data block " + blockname + " contains only one column or no data";
        msg += " at all: the file might be corrupt.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg, _FILE_CORRUPT_);
      }
    }
    else{
      // there is no data in the block, but this might be a nominal operation
      return false;
    }

    // extract coefficients and list of temperatures
    coeffs = xmatrix<double>(nrows, ncols-1);
    T = xvector<double>(nrows);

    for (uint i=0; i<nrows; i++){
      aurostd::string2tokens(v_data[i], tokens);
      if (tokens.size() != ncols){
        msg = "Data block " + blockname + " does not have the same amount of";
        msg += " columns in each line: the file " + filename + " might be corrupt.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg, _FILE_CORRUPT_);
      }

      T[i+1] = tokens[0];
      for (uint j=1; j<ncols; j++) coeffs[i+1][j] = tokens[j];
    }

    if (LDEBUG){
      cerr << function << "nrows="  << nrows << " ncols=" << ncols << std::endl;
      cerr << function << "coeffs=" << coeffs << std::endl;
    }

    return true;
  }

  /// Checks if the given POCC structure contains imaginary frequencies by
  /// reading the corresponding flag from filename.
  /// The flag is IMAG and is set to YES if it contains imaginary frequencies.
  bool hasImaginary(const string& filename, const string &QHA_method)
  {
    string function = "hasImaginary():", msg = "";

    vector<string> vlines;
    bool has_imaginary = false;
    if (!aurostd::efile2vectorstring(filename, vlines)){
      msg = "File " + filename + " does not exist.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg, _FILE_NOT_FOUND_);
    }

    vector<string> tokens;
    for (uint i=0; i<vlines.size(); i++){
      if (vlines[i].find("["+QHA_method+"]") != std::string::npos){
        if (vlines[i].find("IMAG") != std::string::npos){
          aurostd::string2tokens(vlines[i],tokens,"=");
          if (tokens.size() != 2){
            msg = "Incorrect number of tokens: should be 2 instead of ";
            msg += tokens.size();
            throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg, _FILE_CORRUPT_);
          }

          has_imaginary = tokens[1].find("YES") != std::string::npos;
          break;
        }
      }
    }
    return has_imaginary;
  }

  /// Reads POCC structure calculation parameters (minimum and maximum volumes)
  /// from the file filename.
  void readCoeffParameters(const string& filename, double &Vmin, double &Vmax)
  {
    string function = "readCoeffParameters():", msg = "";

    vector<string> vlines;
    Vmin = 0; Vmax = 0;
    if (!aurostd::efile2vectorstring(filename, vlines)) return;

    bool Vmin_found = false, Vmax_found = false;

    vector<double> tokens;
    vector<string> stokens;
    for (uint i=0; i<vlines.size(); i++){
      if (vlines[i].find("VMIN") != std::string::npos){
        aurostd::string2tokens(vlines[i],tokens,"=");
        if (tokens.size() != 2){
          msg = "Number of tokens for VMIN parameter is ";
          msg += aurostd::utype2string(tokens.size()) + " instead of 2:";
          msg += "file " + filename + " might be corrupt.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg, _FILE_CORRUPT_);
        }
        Vmin = tokens[1];
        Vmin_found = true;
      }

      if (vlines[i].find("VMAX") != std::string::npos){
        aurostd::string2tokens(vlines[i],tokens,"=");
        if (tokens.size() != 2){
          msg = "Number of tokens for VMAX parameter is ";
          msg += aurostd::utype2string(tokens.size()) + " instead of 2:";
          msg += "file " + filename + " might be corrupt.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg, _FILE_CORRUPT_);
        }
        Vmax = tokens[1];
        Vmax_found = true;
      }

      if (Vmin_found && Vmax_found) break;
    }
  }

  /// Calculates the logarithm of the partition function.
  double logZ(const xvector<double> &E, const vector<int> &degeneracies, double T)
  {
    xvector<double> khi = E;
    for (int i=khi.lrows; i<=khi.urows; i++) khi[i] /= -(KBOLTZEV*T);

    double shift = max(khi);
    double sum = 0.0;

    for (int i=khi.lrows; i<=khi.urows; i++){
      sum += degeneracies[i-khi.lrows] * exp(khi[i]-shift);
    }

    return shift + std::log(sum);
  }

  /// Calculates the thermal expansion coefficient as a logarithmic derivative
  /// of equilibrium volume employing Savitzky-Golay filter for the differentiation.
  xvector<double> calcThermalExpansionSG(const xvector<double> &volumes, double dT)
  {
    // Convolution weights for Savitzky-Golay 5pt cubic filter as reported in:
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

  /// Calculates the isobaric specif heat as a second derivative of the free
  /// energy multiplied by (-1)*temperature employing Savitzky-Golay filter for 
  /// the differentiation.
  xvector<double> calcIsobaricSpecificHeatSG(const xvector<double> &free_energies, double dT)
  {
    // Convolution weights for Savitzky-Golay 5pt cubic filter as reported in:
    // "General least-squares smoothing and differentiation by the convolution (Savitzky-Golay) method"
    // Peter A. Gorry Analytical Chemistry 1990 62 (6), 570-573
    // https://doi.org/10.1021/ac00205a007
    const static xmatrix<double> SGmat(5,5);
    SGmat[1][1]= 1.2857142857142856E+000; SGmat[1][2]= 7.8571428571428559E-001; SGmat[1][3]= 2.8571428571428570E-001;
    SGmat[2][1]=-2.1428571428571428E+000; SGmat[2][2]=-1.1428571428571428E+000; SGmat[2][3]=-1.4285714285714285E-001;
    SGmat[3][1]=-2.8571428571428570E-001; SGmat[3][2]=-2.8571428571428570E-001; SGmat[3][3]=-2.8571428571428570E-001;
    SGmat[4][1]= 1.8571428571428570E+000; SGmat[4][2]= 8.5714285714285698E-001; SGmat[4][3]=-1.4285714285714285E-001;
    SGmat[5][1]=-7.1428571428571419E-001; SGmat[5][2]=-2.1428571428571425E-001; SGmat[5][3]= 2.8571428571428570E-001;

    SGmat[1][4]=-2.1428571428571425E-001; SGmat[1][5]=-7.1428571428571419E-001;
    SGmat[2][4]= 8.5714285714285698E-001; SGmat[2][5]= 1.8571428571428570E+000;
    SGmat[3][4]=-2.8571428571428570E-001; SGmat[3][5]=-2.8571428571428570E-001;
    SGmat[4][4]=-1.1428571428571428E+000; SGmat[4][5]=-2.1428571428571428E+000;
    SGmat[5][4]= 7.8571428571428559E-001; SGmat[5][5]= 1.2857142857142856E+000;

    const static xvector<double> SGvec(5);
    SGvec[1]= 2.8571428571428570E-001;
    SGvec[2]=-1.4285714285714285E-001;
    SGvec[3]=-2.8571428571428570E-001;
    SGvec[4]=-1.4285714285714285E-001;
    SGvec[5]= 2.8571428571428570E-001;
    ////////////////////////////////////////////////////////////////////////////////

    string function = "calcThermalExpansionSG():", msg = "";
    int npoints = free_energies.rows;
    if (npoints<5){
      msg = "Savitzky-Golay filter requires at least 5 points: only ";
      msg += aurostd::utype2string(npoints) + " were provided.";
      msg += " Make sure that the minimum range among the set of POCC structures,";
      msg += "defined by the [AFLOW_APL]TPT parameter for each, is reasonable.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg, _INDEX_ILLEGAL_);
    }

    xvector<double> endpoints(5), dummy(5);
    xvector<double> Cp(npoints);

    // calculate derivatives for the first 2 points
    for (int i=1; i<=5; i++) dummy[i] = free_energies[i];
    endpoints = dummy*SGmat;
    for (int i=1; i<=2; i++) Cp[i] = endpoints[i];

    // calculate derivatives for the [3:end-3] points
    int id = 0;
    for (int i=3; i<=npoints-2; i++){
      Cp[i] = 0.0;
      for (int j=1; j<=5; j++){
        id = i - 3 + j;
        Cp[i] += SGvec[j]*free_energies[id];
      }
    }

    // calculate derivatives for the last 2 points
    for (int i=1; i<=5; i++) dummy[i] = free_energies[npoints-5+i];
    endpoints = dummy*SGmat;
    for (int i=4; i<=5; i++) Cp[npoints-5+i] = endpoints[i];

    // calculate the coefficient of thermal expansion
    for (int i=1; i<=npoints; i++){
      Cp[i] /= std::pow(dT, 2);
    }

    return Cp;
  }


  /// Calculates POCC-average of QHA-related properties.
  void POccCalculator::calculateQHAProperties(apl::QHAmethod qha_method,
      apl::EOSmethod eos_method)
  {
    string function = XPID +  "POccCalculator::calculateQHAProperties():";
    string msg = "", filename = "";

    bool LDEBUG = false || _DEBUG_POCC_QHA_ || XHOST.DEBUG;

    msg = "Performing POCC+QHA post-processing step for ";
    msg += apl::QHAmethod2label(qha_method) + " with ";
    msg += apl::EOSmethod2label(eos_method) + " EOS.";
    pflow::logger(_AFLOW_FILE_NAME_, function, msg, m_aflags.Directory,
      *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    double POCC_Vmin = DBL_MAX, POCC_Vmax = 0;
    double Vmin = 0, Vmax = 0;
    vector<bool> v_is_unstable(m_ARUN_directories.size());

    xvector<double> T;
    xmatrix<double> coeffs;

    vector<xvector<double> > T_list;
    vector<xmatrix<double> > coeffs_list;

    for (uint i=0; i<m_ARUN_directories.size(); i++){
      filename = m_aflags.Directory + "/" + m_ARUN_directories[i] + "/";
      filename += DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_IMAG_FILE;

      if (!aurostd::EFileExist(filename)){
        msg = "File " + filename + " does not exists:";
        msg += " the calculation will be stopped.";
        pflow::logger(_AFLOW_FILE_NAME_, function, msg, m_aflags.Directory,
         *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
        return;
      }

      v_is_unstable[i] = hasImaginary(filename, apl::QHAmethod2label(qha_method));

      if (!v_is_unstable[i]){
        filename = m_aflags.Directory + "/" + m_ARUN_directories[i] + "/";
        filename += DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_COEFF_FILE;

        readCoeffParameters(filename, Vmin, Vmax);

        POCC_Vmin = min(Vmin, POCC_Vmin);
        POCC_Vmax = max(Vmax, POCC_Vmax);

        if (!readCoeffData(filename, apl::QHAmethod2label(qha_method) + "_" + 
              apl::EOSmethod2label(eos_method), T, coeffs))
        {
          msg = "No data was extracted for " + apl::QHAmethod2label(qha_method);
          msg += " with " + apl::EOSmethod2label(eos_method) + " EOS.";
          msg += " The calculation will be stopped.";
          pflow::logger(_AFLOW_FILE_NAME_, function, msg, m_aflags.Directory,
           *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
          return;
        }
        T_list.push_back(T);
        coeffs_list.push_back(coeffs);
      }
    }

    // to be sure, make the region wider: +-5%
    POCC_Vmin *= 0.95;
    POCC_Vmax *= 1.05;

    if (LDEBUG) cerr << function << " POCC_Vmin = " << POCC_Vmin << " POCC_Vmax = " << POCC_Vmax << std::endl;

    // collect degeneracies only for stable structures
    vector<int> degeneracies;
    unsigned long long int isupercell=0;
    for (std::list<POccSuperCellSet>::iterator it = l_supercell_sets.begin();
         it != l_supercell_sets.end(); ++it)
    {
      isupercell=std::distance(l_supercell_sets.begin(),it);
      if (!v_is_unstable[isupercell]) degeneracies.push_back((*it).getDegeneracy());
    }

    if (LDEBUG){
      cerr << function << " Degeneracies:" << std::endl;
      for (uint i=0; i<degeneracies.size(); i++){
        cerr << function << " id = " << i << " deg = " << degeneracies[i] << std::endl;
      }
    }

    // check that there is data to work with
    uint Nstructures = T_list.size();
    if (!Nstructures){
      msg="No data was extracted: check that all QHA calculations are complete.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
    }

    // the temperature range, where QHA is able to calculate the thermodynamic
    // properties, might differ for different POCC structures:
    // the temperature region for the averaged POCC "material" is at least the
    // lowest range among the structures
    int nrows = T_list[0].rows;
    for (uint i=1; i<Nstructures; i++){
      nrows = std::min(nrows, T_list[i].rows);
    }

    // check that calculations for POCC structures are consistent and the same
    // set of temperatures was used for QHA calculation for each of them
    for (uint i=0; i<Nstructures-1; i++){
      for (int row=1; row<=nrows; row++){
        if (!aurostd::isequal(T_list[i][row], T_list[i+1][row])){
          msg="Inconsistent list of temperatures among different ";
          msg+="POCC::QHA calculations.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg,
              _VALUE_ILLEGAL_);
        }
      }
    }
    T = T_list[0];

    // sanity/corruption check: the number of fitting coefficients should
    // be the same for all POCC structures
    int ncols = coeffs_list[0].cols;
    for (uint i=1; i<Nstructures; i++){
      if (ncols != coeffs_list[i].cols){
        msg="Inconsistent number of fitting coefficients among different ";
        msg+="POCC::QHA calculations.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, msg,
              _VALUE_ILLEGAL_);
      }
    }

    // generate a set of volumes to determine the EOS for the POCC "material"
    static int Nvolumes = 5;
    xvector<double> volumes(Nvolumes);
    for (int i=1; i<=Nvolumes; i++){
      volumes[i] = POCC_Vmin + (POCC_Vmax-POCC_Vmin)*(i-1)/(Nvolumes-1);
    }

    // use QHA to calculate the thermodynamic properties of the POCC "material"
    apl::QHA qha;

    xvector<double> F(Nvolumes), E(Nstructures);
    xvector<double> Feq(nrows), Veq(nrows), B(nrows), Bprime(nrows);
    try{
      for (int i=1; i<=nrows; i++){
        if (aurostd::isequal(T[i], 0.0)) T[i] = _ZERO_TOL_;

        for (int v=volumes.lrows; v<=volumes.urows; v++){
          for (uint s=0; s<Nstructures; s++){
            // EOS for each POCC structure
            E[s+1] = qha.evalEOSmodel(volumes[v], coeffs_list[s](i), eos_method);
          }
          if (LDEBUG) cerr << function << " E: " << E << std::endl;
          // the set of free energies (includes the contribution from the configurational
          // entropy) to be used in the fitting procedure
          F[v] = -KBOLTZEV * T[i] * logZ(E, degeneracies, T[i]);
        }

        // fit EOS for the POCC "material"
        qha.fitToEOSmodel(volumes, F, eos_method);
        Feq[i]    = qha.EOS_energy_at_equilibrium;
        Veq[i]    = qha.EOS_volume_at_equilibrium;
        B[i]      = qha.EOS_bulk_modulus_at_equilibrium;
        Bprime[i] = qha.EOS_Bprime_at_equilibrium;
        if (LDEBUG){
          cerr << function << " V: " << volumes << std::endl;
          cerr << function << " F: " << F << std::endl;
        }
      }
    } catch (aurostd::xerror e){
      // QHA throws _VALUE_RANGE_ exception only when there is no minimum in
      // the energy-volume relation: at this point the calculation of 
      // thermodynamic properties should be stopped and a warning should be
      // printed, and all calculated data should be saved to the file
      if (e.error_code == _VALUE_RANGE_){
        pflow::logger(e.whereFileName(), e.whereFunction(), e.error_message, 
          m_aflags.Directory, *p_FileMESSAGE, *p_oss, _LOGGER_WARNING_);
      }
      else{
        throw;
      }
    }

    // calculation of the thermal expansion and heat capacity involves derivatives:
    // dT is a temperature step for the derivative
    double dT = (max(T)-min(T))/(nrows-1);
    if (LDEBUG) cerr << function << " dT = " << dT << std::endl;

    // calculate the thermal expansion coefficients
    xvector<double> beta = calcThermalExpansionSG(Veq, dT);

    // calculate the heat capacity
    xvector<double> Cp = calcIsobaricSpecificHeatSG(Feq, dT);
    for (int i=1; i<=nrows; i++) Cp[i] *= -T[i]/KBOLTZEV;

    xvector<double> Cv(nrows);
    for (int i=1; i<=nrows; i++){
      Cv[i] = Cp[i] - Veq[i]*T[i]*B[i]*pow(beta[i],2)/eV2GPa/KBOLTZEV;
    }

    // calculate the average Grueneisen parameters
    xvector<double> gamma(nrows);
    for (int i=1; i<=nrows; i++){
      gamma[i] = (beta[i]/Cv[i])*B[i]*Veq[i]/eV2GPa/KBOLTZEV;
    }

    // prepare the data for the output
    stringstream file;
    file.precision(10);

    // write header
    string block = "[" + apl::QHAmethod2label(qha_method) + "_";
    block += apl::EOSmethod2label(eos_method) + "_THERMO]";
    file << AFLOWIN_SEPARATION_LINE << std::endl;
    file << block + "START" << std::endl;
    file << setw(5)  << "#T[K]"          << setw(SW) << ' ' <<
      setw(TW) << "V[A^3/atom]"          << setw(SW) << ' ' <<
      setw(TW) << "F(V)[eV/atom]"        << setw(SW) << ' ' <<
      setw(TW) << "B[GPa]"               << setw(SW) << ' ' <<
      setw(TW) << "beta[10^-5/K]"        << setw(SW) << ' ' <<
      setw(TW) << "Cv(V)[kB/atom]"       << setw(SW) << ' ' <<
      setw(TW) << "Cp(V)[kB/atom]"       << setw(SW) << ' ' <<
      setw(TW) << "gamma(beta,B,Cv(V))"  << setw(SW) << ' ' <<
      setw(TW) << "Bprime"
      << std::endl;

    for (int i=1; i<=nrows; i++){
      file << setw(5) << T[i]         << setw(SW) << ' ' <<
        setw(TW) << Veq[i]            << setw(SW) << ' ' <<
        setw(TW) << Feq[i]            << setw(SW) << ' ' <<
        setw(TW) << B[i]              << setw(SW) << ' ' <<
        setw(TW) << beta[i] * 1e5     << setw(SW) << ' ' <<
        setw(TW) << Cv[i]             << setw(SW) << ' ' <<
        setw(TW) << Cp[i]             << setw(SW) << ' ' <<
        setw(TW) << gamma[i]          << setw(SW) << ' ' <<
        setw(TW) << Bprime[i]
        << std::endl;
    }
    file << block + "STOP" << std::endl;
    file << AFLOWIN_SEPARATION_LINE << std::endl;

    // save data into the file
    filename = POCC_FILE_PREFIX + "qha." + DEFAULT_QHA_THERMO_FILE;
    if (aurostd::FileExist(filename)){
      if (!aurostd::stringstream2file(file, filename, "APPEND")){
        msg = "Error writing to " + filename + " file.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
      }
    }
    else{
      if (!aurostd::stringstream2file(file, filename)){
        msg = "Error writing to " + filename + " file.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_,function,msg,_FILE_ERROR_);
      }
    }
  }

  /// Calculates POCC-average of QHA-related properties for various QHA methods
  /// with various EOS models.
  void POccCalculator::calculateQHAProperties()
  {
    static const int N_QHA_methods = 3;
    static const apl::QHAmethod QHA_methods[N_QHA_methods] = {apl::QHA_CALC,
      apl::QHA3P_CALC, apl::QHANP_CALC};
    static const int N_EOS_methods = 5;
    static const apl::EOSmethod EOS_methods[N_EOS_methods] = {apl::EOS_SJ,
      apl::EOS_BIRCH_MURNAGHAN2, apl::EOS_BIRCH_MURNAGHAN3, apl::EOS_BIRCH_MURNAGHAN4,
      apl::EOS_MURNAGHAN};

    string function = XPID +  "POccCalculator::calculateQHAProperties():";
    string msg = "", filename = "";

    msg = "Performing POCC+QHA post-processing step.";
    pflow::logger(_AFLOW_FILE_NAME_, function, msg, m_aflags.Directory,
      *p_FileMESSAGE, *p_oss, _LOGGER_MESSAGE_);

    filename = POCC_FILE_PREFIX + "qha." + DEFAULT_QHA_THERMO_FILE;
    if (aurostd::FileExist(filename)) aurostd::RemoveFile(filename);

    try{
      for (int i=0; i<N_QHA_methods; i++){
        for (int j=0; j<N_EOS_methods; j++){
          POccCalculator::calculateQHAProperties(QHA_methods[i], EOS_methods[j]);
        }
      }
    } catch (aurostd::xerror e) {
      pflow::logger(e.whereFileName(), e.whereFunction(), e.error_message,
          m_aflags.Directory, *p_FileMESSAGE, *p_oss, _LOGGER_ERROR_);
    }
  }
}

#endif
