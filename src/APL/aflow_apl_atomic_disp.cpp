//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *               Pinku Nath - Duke University 2014 - 2016                  *
// *                  Marco Esters - Duke University 2020                    *
// *                                                                         *
//****************************************************************************
// Written by Marco Esters based on work by Pinku Nath

#include "aflow_apl.h"

// Some parts are written within the C++0x support in GCC, especially std::thread,
// which is implemented in gcc 4.4 and higher. For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION>= 40400
#define AFLOW_APL_MULTITHREADS_ENABLE
#include <thread>
#else
#warning "The multithread parts of APL will not be included, since they need gcc 4.4 and higher (C++0x support)."
#endif

static const string _APL_ADISP_MODULE_ = "APL";  // for the logger
static const xcomplex<double> iONE(0.0, 1.0);

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         CONSTRUCTORS/DESTRUCTORS                         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  AtomicDisplacements::AtomicDisplacements() {
    free();
  }

  AtomicDisplacements::AtomicDisplacements(PhononCalculator& pc) {
    free();
    _pc = &pc;
  }

  AtomicDisplacements::AtomicDisplacements(const AtomicDisplacements& that) {
    free();
    copy(that);
  }

  AtomicDisplacements& AtomicDisplacements::operator=(const AtomicDisplacements& that) {
    if (this != &that) {
      free();
      copy(that);
    }
    return *this;
  }

  AtomicDisplacements::~AtomicDisplacements() {
    free();
  }

  void AtomicDisplacements::copy(const AtomicDisplacements& that) {
    _eigenvectors = that._eigenvectors;
    _frequencies = that._frequencies;
    _displacement_matrices = that._displacement_matrices;
    _displacement_modes = that._displacement_modes;
    _pc = that._pc;
    _qpoints = that._qpoints;
    _temperatures = that._temperatures;
  }

  void AtomicDisplacements::free() {
    _eigenvectors.clear();
    _frequencies.clear();
    _displacement_matrices.clear();
    _displacement_modes.clear();
    _qpoints.clear();
    _temperatures.clear();
  }

  void AtomicDisplacements::clear(PhononCalculator& pc) {
    free();
    _pc = &pc;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                             EIGENVECTORS                                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  void AtomicDisplacements::calculateEigenvectors() {
    _eigenvectors.clear();
    _frequencies.clear();
    int nq = (int) _qpoints.size();
    if (nq == 0) return;
    uint natoms = _pc->getInputCellStructure().atoms.size();
    uint nbranches = _pc->getNumberOfBranches();
    _eigenvectors.resize(nq, vector<vector<xvector<xcomplex<double> > > >(nbranches, vector<xvector<xcomplex<double> > >(natoms, xvector<xcomplex<double> >(3))));
    _frequencies.resize(nq, vector<double> (nbranches, 0.0));
#ifdef AFLOW_APL_MULTITHREADS_ENABLE
    int ncpus = _pc->getNCPUs();
    if (ncpus > nq) ncpus = nq;
    if (ncpus < 1) ncpus = 1;
    if (ncpus > 1) {
      vector<vector<int> > thread_dist = getThreadDistribution(nq, ncpus);
      vector<std::thread*> threads;
      for (int i = 0; i < ncpus; i++) {
        threads.push_back(new std::thread(&AtomicDisplacements::calculateEigenvectorsInThread, this, thread_dist[i][0], thread_dist[i][1]));
      }
      for (int i = 0; i < ncpus; i++) {
        threads[i]->join();
        delete threads[i];
      }
    } else {
      calculateEigenvectorsInThread(0, nq);
    }
#else
    calculateEigenvectorsInThread(0, nq);
#endif
  }

  void AtomicDisplacements::calculateEigenvectorsInThread(int startIndex, int endIndex) {
    uint nbranches = _pc->getNumberOfBranches();
    uint natoms = _pc->getInputCellStructure().atoms.size();
    xvector<double> freq(nbranches);
    xmatrix<xcomplex<double> > eig(nbranches, nbranches, 1, 1);
    for (int i = startIndex; i < endIndex; i++) {
      freq = _pc->getFrequency(_qpoints[i].cpos, apl::THZ, eig);
      for (uint br = 0; br < nbranches; br++) {
        _frequencies[i][br] = freq[br + 1];
        for (uint at = 0; at < natoms; at++) {
          for (int j = 1; j < 4; j++) {
            _eigenvectors[i][br][at][j] = eig[3 * at + j][br + 1];
          }
        }
      }
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                             DISPLACEMENTS                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  void AtomicDisplacements::calculateMeanSquareDisplacements(const QMesh& qmesh, double Tstart, double Tend, double Tstep) {
    _qpoints.clear();
    _temperatures.clear();

    if (Tstart > Tend) {
      string function = "AtomicDisplacements::calculateDisplacements()";
      string message = "Tstart cannot be higher than Tend.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
    }
    for (double T = Tstart; T <= Tend; T += Tstep) _temperatures.push_back(T);

    _qpoints = qmesh.getPoints();
    calculateMeanSquareDisplacementMatrices();
  }

  void AtomicDisplacements::calculateMeanSquareDisplacementMatrices() {
    string message = "Calculating mean square displacement matrices.";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_ADISP_MODULE_, message, _pc->getDirectory(), _pc->getOutputStream(), std::cout);
    _displacement_matrices.clear();
    _displacement_modes.clear();
    calculateEigenvectors();

    uint ntemps = _temperatures.size();
    uint natoms = _pc->getInputCellStructure().atoms.size();
    uint nq = _qpoints.size();
    uint nbranches = _pc->getNumberOfBranches();

    _displacement_matrices.resize(ntemps, vector<xmatrix<xcomplex<double> > >(natoms, xmatrix<xcomplex<double> >(3, 3)));
    vector<xmatrix<xcomplex<double> > > outer(natoms, xmatrix<xcomplex<double> >(3, 3));
    vector<double> masses(natoms);
    for (uint at = 0; at < natoms; at++) masses[at] = AMU2KILOGRAM * _pc->getSupercell().getAtomMass(_pc->getSupercell().pc2scMap(at));

    // Factor 2pi necessary because frequencies are raw
    double prefactor = (PLANCKSCONSTANT_hbar * Hz2THz * std::pow(1e10, 2)/(2 * PI * (double) _qpoints.size()));
    for (uint q = 0; q < nq; q++) {
      for (uint br = 0; br < nbranches; br++) {
        if (_frequencies[q][br] > _AFLOW_APL_EPS_) {
          for (uint at = 0; at < natoms; at++) {
            outer[at] = aurostd::outer_product(_eigenvectors[q][br][at], conj(_eigenvectors[q][br][at]));
          }
          for (uint t = 0; t < ntemps; t++) {
            double prefactor_T = prefactor * ((0.5 + getOccupationNumber(_temperatures[t], _frequencies[q][br]))/_frequencies[q][br]);
            for (uint at = 0; at < natoms; at++) {
              // Add element-wise (much faster)
              for (int i = 1; i < 4; i++) {
                for (int j = 1; j < 4; j++) {
                  _displacement_matrices[t][at][i][j].re += (prefactor_T/masses[at]) * outer[at][i][j].re;
                  _displacement_matrices[t][at][i][j].im += (prefactor_T/masses[at]) * outer[at][i][j].im;
                }
              }
            }
          }
        }
      }
    }
  }

  void AtomicDisplacements::calculateNormalModeDisplacements(const vector<xvector<double> >& qpts, bool coords_are_fractional) {
    _qpoints.clear();
    uint nq = qpts.size();
    _qpoints.resize(nq);
    if (coords_are_fractional) {
      xmatrix<double> f2c = trasp(ReciprocalLattice(_pc->getInputCellStructure()));
      for (uint q = 0; q < nq; q++) {
        _qpoints[q].fpos = qpts[q];
        _qpoints[q].cpos = f2c * qpts[q];
      }
    } else {
      xmatrix<double> c2f = inverse(trasp(ReciprocalLattice(_pc->getInputCellStructure())));
      for (uint q = 0; q < nq; q++) {
        _qpoints[q].cpos = qpts[q];
        _qpoints[q].fpos = c2f * qpts[q];
      }
    }
    calculateNormalModeDisplacements();
  }

  void AtomicDisplacements::calculateNormalModeDisplacements() {
    _displacement_matrices.clear();
    _displacement_modes.clear();
    _temperatures.clear();
    calculateEigenvectors();

    uint nq = _qpoints.size();
    uint nbranches = _pc->getNumberOfBranches();
    uint natoms = _pc->getInputCellStructure().atoms.size();
    _displacement_modes.resize(nq, vector<vector<xvector<xcomplex<double> > > >(nbranches, vector<xvector<xcomplex<double> > >(natoms)));

    vector<double> masses(natoms);
    for (uint at = 0; at < natoms; at++) masses[at] = _pc->getSupercell().getAtomMass(_pc->getSupercell().pc2scMap(at));

    for (uint q = 0; q < nq; q++) {
      for (uint br = 0; br < nbranches; br++) {
        if (_frequencies[q][br] > _AFLOW_APL_EPS_) {
          for (uint at = 0; at < natoms; at++) {
            _displacement_modes[q][br][at] = _eigenvectors[q][br][at]/sqrt(masses[at]);
          }
        }
      }
    }
  }

  double AtomicDisplacements::getOccupationNumber(double T, double f) {
    if (T < _AFLOW_APL_EPS_) return 0.0;
    else return (1.0/(exp(BEfactor_h_THz * f/T) - 1));
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            GETTER FUNCTIONS                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  const vector<double>& AtomicDisplacements::getTemperatures() const {
    return _temperatures;
  }

  const vector<vector<xmatrix<xcomplex<double> > > >& AtomicDisplacements::getDisplacementMatrices() const {
    return _displacement_matrices;
  }

  vector<vector<xvector<double> > > AtomicDisplacements::getDisplacementVectors() const {
    vector<vector<xvector<double> > >  disp_vec;
    uint ntemps = _displacement_matrices.size();
    if (ntemps > 0) {
      uint natoms = _displacement_matrices[0].size();
      disp_vec.resize(ntemps, vector<xvector<double> >(natoms, xvector<double>(3)));
      for (uint t = 0; t < ntemps; t++) {
        for (uint at = 0; at < natoms; at++) {
          for (int i = 1; i < 4; i++) {
            disp_vec[t][at][i] = _displacement_matrices[t][at][i][i].re;
          }
        }
      }
    }
    return disp_vec;
  }

  const vector<vector<vector<xvector<xcomplex<double> > > > >& AtomicDisplacements::getModeDisplacements() const {
    return _displacement_modes;
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                FILE I/O                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  void AtomicDisplacements::writeMeanSquareDisplacementsToFile(string filename) {
    filename = aurostd::CleanFileName(filename);
    string message = "Writing mean square displacements into file " + filename + ".";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_ADISP_MODULE_, message, _pc->getDirectory(), _pc->getOutputStream(), std::cout);
    vector<vector<xvector<double> > > disp_vec = getDisplacementVectors();
    stringstream output;
    string tag = "[APL_DISPLACEMENTS]";

    output << AFLOWIN_SEPARATION_LINE << std::endl;
    output << tag << "SYSTEM=" << _pc->getSystemName() << std::endl;
    output << tag << "START" << std::endl;
    output << "#" << std::setw(9) << "T (K)" << setw(15) << "Species"
      << std::setw(15) << "x (A^2)" << std::setw(15) << "y (A^2)" << std::setw(15) << "z (A^2)" << std::endl;
    output << std::fixed << std::showpoint;
    for (uint t = 0; t < _temperatures.size(); t++) {
      output << std::setw(10) << std::setprecision(2) << _temperatures[t];
      for (uint at = 0; at < disp_vec[t].size(); at++) {
        if (at == 0) {
          output << std::setw(15) << _pc->getInputCellStructure().atoms[at].cleanname;
        } else {
          output << std::setw(25) << _pc->getInputCellStructure().atoms[at].cleanname;
        }
        for (int i = 1; i < 4; i++) output << std::setw(15) << std::setprecision(8) << disp_vec[t][at][i];
        output << std::endl;
      }
    }
    output << tag << "STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;

    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      string function = "AtomicDisplacements::writeMeanSquareDisplacementsToFile()";
      message = "Could not write to file " + filename + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
  }

  void AtomicDisplacements::writeSceneFileXcrysden(string filename, const vector<vector<vector<double> > >& disp, int nperiods) {
    filename = aurostd::CleanFileName(filename);
    string message = "Writing atomic displacements in XCRYSDEN format into file " + filename + ".";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_ADISP_MODULE_, message, _pc->getDirectory(), _pc->getOutputStream(), std::cout);

    uint natoms = disp.size();
    uint nsteps = disp[0].size();

    stringstream output;
    output << "ANIMSTEPS " << (nperiods * nsteps) << std::endl;
    output << "CRYSTAL" << std::endl;
    output << "PRIMVEC" << std::endl;
    output << _pc->getInputCellStructure().lattice;  // no std::endl necessary

    int step = 1;
    for (int i = 0; i < nperiods; i++) {
      for (uint j = 0; j < nsteps; j++) {
        output << "PRIMCOORD " << step << std::endl;
        output << std::setw(4) << natoms << " 1" << std::endl;
        for (uint at = 0; at < natoms; at++) {
          output << std::setw(5) << GetAtomName((uint) disp[j][at][0]);
          for (int k = 1; k < 7; k++) {
            output << std::setw(15) << std::fixed << std::setprecision(8) << disp[j][at][k];
          }
          output << std::endl;
        }
        step++;
      }
    }

    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      string function = "AtomicDisplacements::writeMeanSquareDisplacementsToFile()";
      message = "Could not write to file " + filename + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
  }

  void AtomicDisplacements::writeSceneFileVsim(string filename, const xstructure& xstr_projected,
      const vector<vector<vector<xvector<xcomplex<double> > > > >& displacements) {
    filename = aurostd::CleanFileName(filename);
    string message = "Writing atomic displacements in V_SIM format into file " + filename + ".";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_ADISP_MODULE_, message, _pc->getDirectory(), _pc->getOutputStream(), std::cout);

    stringstream output;
    // Lattice
    output << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[1][1]
      << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[2][1]
      << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[2][2] << std::endl;
    output << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[3][1]
      << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[3][2]
      << std::setw(15) << std::setprecision(8) << std::fixed << xstr_projected.lattice[3][3] << std::endl;
    // Atoms
    uint natoms = xstr_projected.atoms.size();
    for (uint at = 0; at < natoms; at++) {
      for (int i = 1; i < 4; i++) output << std::setw(15) << std::setprecision(8) << xstr_projected.atoms[at].cpos[i];
      output << std::setw(4) << xstr_projected.atoms[at].cleanname << std::endl;
    }

    for (uint q = 0; q < displacements.size(); q++) {
      for (uint br = 0; br < displacements[q].size(); br++) {
        output << "#metaData: qpt=[";
        for (int i = 1; i < 4; i++) output << _qpoints[q].fpos[i] << ";";
        output << _frequencies[q][br] << "\\" << std::endl;
        for (uint at = 0; at < natoms; at++) {
          output << "#;";
          for (int i = 1; i < 4; i++) output << displacements[q][br][at][i].re << ";";
          for (int i = 1; i < 4; i++) {
            output << displacements[q][br][at][i].im;
            if (i < 3) output << ";";
            else output << "\\" << std::endl;
          }
          output << std::endl;
        }
        output << "# ]" << std::endl;
      }
    }

    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      string function = "AtomicDisplacements::writeMeanSquareDisplacementsToFile()";
      message = "Could not write to file " + filename + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _FILE_ERROR_);
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                               INTERFACE                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  void createAtomicDisplacementSceneFile(const aurostd::xoption& vpflow) {
    string function = "apl::createAtomicDisplacementSceneFile()";
    string message = "";
    ofstream mf("/dev/null");

    // Parse command line options
    string directory = vpflow.getattachedscheme("ADISP::DIRECTORY");
    if (directory.empty()) directory = "./";
    else directory = aurostd::CleanFileName(directory + "/");

    // Format
    string format = aurostd::toupper(vpflow.getattachedscheme("ADISP::FORMAT"));
    if (format.empty()) format = aurostd::toupper(DEFAULT_APL_ADISP_SCENE_FORMAT);
    string allowed_formats_str = "XCRYSDEN,V_SIM";
    vector<string> allowed_formats;
    aurostd::string2tokens(allowed_formats_str, allowed_formats, ",");
    if (!aurostd::withinList(allowed_formats, format)) {
      message = "Unregonized format " + format + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ILLEGAL_);
    }
    
    // Amplitude
    string amplitude_str = vpflow.getattachedscheme("ADISP::AMPLITUDE");
    double amplitude = 0.0;
    if (amplitude_str.empty()) amplitude = DEFAULT_APL_ADISP_AMPLITUDE;
    else amplitude = aurostd::string2utype<double>(amplitude_str);
    if (amplitude < _AFLOW_APL_EPS_) {
      message = "Amplitude must be positive.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ILLEGAL_);
    }

    // Number of steps per period
    string nsteps_str = vpflow.getattachedscheme("ADISP::STEPS");
    int nsteps = 0;
    if (format != "V_SIM") {
      if (nsteps_str.empty()) nsteps = DEFAULT_APL_ADISP_NSTEPS;
      else nsteps = aurostd::string2utype<int>(nsteps_str);
      if (nsteps < 1) {
        message = "Number of steps must be a positive integer";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ILLEGAL_);
      }
    }

    // Number of periods
    string nperiods_str = vpflow.getattachedscheme("ADISP::PERIODS");
    int nperiods = 0;
    if (format != "V_SIM") {
      if (nperiods_str.empty()) nperiods = DEFAULT_APL_ADISP_NPERIODS;
      else nperiods = aurostd::string2utype<int>(nperiods_str);
      if (nsteps < 1) {
        message = "Number of periods must be a positive integer";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ILLEGAL_);
      }
    }

    // Range/supercell
    string supercell_str = vpflow.getattachedscheme("ADISP::SUPERCELL");
    string range_str = vpflow.getattachedscheme("ADISP::RANGE");
    if (format != "V_SIM") {
      if (supercell_str.empty() && range_str.empty()) {
        supercell_str = "1x1x1";
      } else if (!supercell_str.empty() && !range_str.empty()) {
        message = "Cannot specify supercell and range at the same time.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_AMBIGUOUS_);
      } else if (!supercell_str.empty()) {
        vector<string> tokens;
        aurostd::string2tokens(supercell_str, tokens, "xX");
        if (tokens.size() != 3) {
          message = "Broken supercell format.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ILLEGAL_);
        }
      }
    }

    // Branches
    string branches_str = vpflow.getattachedscheme("ADISP::BRANCHES");
    vector<int> branches;
    if (format != "V_SIM") {
      if (branches_str.empty()) {
        message = "No branches selected. Displacements will be calculated for all.";
        pflow::logger(_AFLOW_FILE_NAME_, _APL_ADISP_MODULE_, directory, mf, std::cout);
      } else {
        aurostd::string2tokens(branches_str, branches, ",");
      }
    }

    // q-points
    string qpoints_str = vpflow.getattachedscheme("ADISP::QPOINTS");
    vector<xvector<double> > qpoints;
    if (qpoints_str.empty()) {
      message = "No q-points given.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_MISSING_);
    } else {
      vector<string> tokens, tokens_qpt;
      aurostd::string2tokens(qpoints_str, tokens, ";");
      for (uint i = 0; i < tokens.size(); i++) {
        aurostd::string2tokens(tokens[i], tokens_qpt, ",");
        if (tokens_qpt.size() == 3) {
          xvector<double> q(3);
          for (int j = 0; j < 3; j++) q[j + 1] = aurostd::string2utype<double>(tokens_qpt[j]);
          BringInCellInPlace(q, _ZERO_TOL_, 0.5, -0.5);
          qpoints.push_back(q);
        } else {
          message = "Broken q-points format.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ILLEGAL_);
        }
      }
    }

    // Initialize phonon calculator
    string statefile = directory + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_STATE_FILE;
    Supercell sc_pcalc(statefile, mf);
    PhononCalculator pc(sc_pcalc, mf);
    pc.setDirectory(directory);
    string hibfile = directory + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_HARMIFC_FILE;
    pc.awake(hibfile, false);
    // Must project to primitive or the vibrations will be incorrect
    if (!sc_pcalc.projectToPrimitive()) {
      message = "Could not project to primitive structure.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }
    // Check branches
    int nbr = (int) branches.size();
    if (nbr == 0) {
      nbr = (int) pc.getNumberOfBranches();
      for (int br = 0; br < nbr; br++) branches.push_back(br);
    } else {
      int nbranches = pc.getNumberOfBranches();
      for (int br = 0; br < nbr; br++) {
        if (branches[br] >= nbranches) {
          message = "Index " + aurostd::utype2string<int>(branches[br]) + " out of range.";
          throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_BOUNDS_);
        }
      }
    }

    // Done with setup - calculate displacements
    AtomicDisplacements ad(pc);
    ad.calculateNormalModeDisplacements(qpoints);

    if (format == "XCRYSDEN") {
      // Create supercell for the scene file
      Supercell scell(pc.getInputCellStructure(), mf);
      double range = AUROSTD_MAX_DOUBLE;
      if (supercell_str.empty()) {
        range = aurostd::string2utype<double>(range_str);
        xvector<int> dims = LatticeDimensionSphere(scell.getInputStructure().lattice, range);
        // Sphere is inside, so increase dimension to be safe
        for (int i = 1; i < 4; i++) dims[i] += 1;
        scell.build(dims, false);
      } else {
        vector<int> tokens;
        aurostd::string2tokens(supercell_str, tokens, "xX");
        scell.build(aurostd::vector2xvector(tokens), false);
      }
      if (!scell.projectToPrimitive()) {
        message = "Could not project to primitive structure.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
      }

      vector<vector<vector<double> > > disp;
      for (uint q = 0; q < qpoints.size(); q++) {
        for (uint br = 0; br < branches.size(); br++) {
          disp = ad.createDisplacementsXcrysden(scell, range, amplitude, q, br, nsteps);
          stringstream filename;
          filename << directory <<  DEFAULT_APL_FILE_PREFIX;
          for (int i = 1; i < 4; i++) filename << "_" << qpoints[q][i];
          filename << "_" << br << ".axsf";
          ad.writeSceneFileXcrysden(filename.str(), disp, nperiods);
        }
      }
    } else if (format == "V_SIM") {
      xstructure xstr_oriented;
      vector<vector<vector<xvector<xcomplex<double> > > > > displacements_oriented;
      ad.getOrientedDisplacementsVsim(xstr_oriented, displacements_oriented);
      string filename = directory + DEFAULT_APL_FILE_PREFIX + "displacements.ascii";
      ad.writeSceneFileVsim(filename, xstr_oriented, displacements_oriented);
    }
  }

  vector<vector<vector<double> > > AtomicDisplacements::createDisplacementsXcrysden(const Supercell& scell,
      double range, double amplitude, int q, int br, int nsteps) {
    const xvector<double>& qpt = _qpoints[q].cpos;
    const vector<xvector<xcomplex<double> > >& disp = _displacement_modes[q][br];

    // Generate a list of atoms that are inside the desired sphere
    const xstructure& scell_str = scell.getInputStructure();
    uint natoms = scell_str.atoms.size();
    vector<int> atoms;
    if (range == AUROSTD_MAX_DOUBLE) {
      for (uint iat = 0; iat < natoms; iat++) atoms.push_back(iat);
    } else {
      for (uint iat = 0; iat < natoms; iat++) {
        if (aurostd::modulus(scell_str.atoms[iat].cpos) <= range) atoms.push_back(iat);
      }
      natoms = atoms.size();
    }

    // Calculate original displacements
    const vector<int>& sc2pcMap = scell._sc2pcMap;
    const vector<int>& pc2scMap = scell._pc2scMap;

    // Calculate original displacements
    vector<vector<vector<double> > > displacements(nsteps, vector<vector<double> >(natoms, vector<double>(7, 0.0)));
    vector<xvector<xcomplex<double> > > displacements_orig(natoms, xvector<xcomplex<double> >(3));
    int at = -1, at_pc = -1, at_eq_sc = -1;
    xvector<double> dist_scell(3);
    for (uint iat = 0; iat < natoms; iat++) {
      at = atoms[iat];
      at_pc = sc2pcMap[at];
      at_eq_sc = pc2scMap[at_pc];
      dist_scell = scell_str.atoms[at].cpos - scell_str.atoms[at_eq_sc].cpos;
      displacements_orig[iat] = amplitude * exp(iONE * scalar_product(qpt, dist_scell)) * disp[at_pc];
      displacements[0][iat][0] = (double) scell_str.atoms[at].atomic_number;
      for (int i = 1; i < 4; i++) {
        displacements[0][iat][i] = scell_str.atoms[at].cpos[i];
        displacements[0][iat][i + 3] = displacements_orig[iat][i].re;
      }
    }

    // Calculate displacements for the rest of the period
    xcomplex<double> phase;
    xvector<xcomplex<double> > disp_step;
    for (int s = 1; s < nsteps; s++) {
      phase = exp(-iONE * ((double) s/(double) nsteps));
      for (uint iat = 0; iat < natoms; iat++) {
        disp_step = phase * displacements_orig[iat];
        at = atoms[iat];
        displacements[s][at][0] = (double) scell_str.atoms[at].atomic_number;
        for (int i = 1; i < 4; i++) {
          displacements[s][at][i] += displacements[s - 1][at][i + 3];
          displacements[s][at][i + 3] = disp_step[i].re;
        }
      }
    }

    return displacements;
  }

  void AtomicDisplacements::getOrientedDisplacementsVsim(xstructure& xstr_oriented,
      vector<vector<vector<xvector<xcomplex<double> > > > >& displacements) {
    // Project the structure such that a points along x and b along
    // the x-y plane as required by V_sim
    const xstructure& xstr_orig = _pc->getInputCellStructure();
    xstr_oriented.clear();

    // Project lattice
    xvector<double> params = Getabc_angles(xstr_orig.lattice, RADIANS);
    xmatrix<double> lattice(3, 3);
    double cosalpha = cos(params[4]);
    double cosbeta = cos(params[5]);
    double cosgamma = cos(params[6]);
    double singamma = sin(params[6]);
    double l32 = (2* cosalpha - 2 * cosbeta * cosgamma)/(2 * singamma);
    lattice[1][1] = params[1];
    lattice[2][1] = params[2] * cosgamma;
    lattice[2][2] = params[2] * singamma;
    lattice[3][1] = params[3] * cosbeta;
    lattice[3][2] = params[3] * l32;
    lattice[3][3] = params[3] * sqrt(1 - std::pow(cosbeta, 2) - std::pow(l32, 2));
    xstr_oriented.lattice = lattice;
    xstr_oriented.f2c = trasp(lattice);
    xstr_oriented.c2f = inverse(xstr_oriented.f2c);

    // Project positions
    uint natoms = xstr_orig.atoms.size();
    for (uint i = 0; i < natoms; i++) {
      _atom at = xstr_orig.atoms[i];
      at.cpos = xstr_oriented.f2c * at.fpos;
      xstr_oriented.atoms.push_back(at);
    }

    // Project displacements
    uint nq = _qpoints.size();
    uint nbr = _pc->getNumberOfBranches();
    displacements.clear();
    displacements.resize(nq, vector<vector<xvector<xcomplex<double> > > >(nbr, vector<xvector<xcomplex<double> > >(natoms)));

    xmatrix<double> transf = xstr_oriented.f2c * xstr_orig.c2f;
    for (uint q = 0; q < nq; q++) {
      for (uint br = 0; br < nbr; br++) {
        for (uint at = 0; at < natoms; at++) {
          displacements[q][br][at] = transf * _displacement_modes[q][br][at];
        }
      }
    }
  }

}  // namespace apl


//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2020           *
// *               Pinku Nath - Duke University 2014 - 2016                  *
// *                  Marco Esters - Duke University 2020                    *
// *                                                                         *
//****************************************************************************
