// [OBSOLETE] #include <iostream>
// [OBSOLETE] #include <sstream>
// [OBSOLETE] #include <string>
// [OBSOLETE] #include <limits>

#include "aflow_apl.h"

//CO - START
// Some parts are written within the C++0x support in GCC, especially the std::thread,
// which is implemented in gcc 4.4 and higher.... For multithreads with std::thread see:
// http://www.justsoftwaresolutions.co.uk/threading/multithreading-in-c++0x-part-1-starting-threads.html
#if GCC_VERSION >= 40400  // added two zeros
#define AFLOW_APL_MULTITHREADS_ENABLE 1
#include <thread>
#else
#warning "The multithread parts of APL will be not included, since they need gcc 4.4 and higher (C++0x support)."
#endif
//CO - END

#define MIN_FREQ_TRESHOLD -0.1

using namespace std;

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  //DOSCalculator::DOSCalculator(IPhononCalculator& pc, IReciprocalPointGrid& rg, Logger& l)  OBSOLETE ME190423
  DOSCalculator::DOSCalculator(IPhononCalculator& pc, QMesh& rg, Logger& l, string method,
      const vector<xvector<double> >& projections)  // ME190624
    : _pc(pc), _rg(rg), _logger(l) {
      clear();
      _bzmethod = method;
      _system = _pc.getSystemName();  // ME190614
      _projections = projections;  // ME190626
      calculateFrequencies();
    }

  // ///////////////////////////////////////////////////////////////////////////

  DOSCalculator::~DOSCalculator() {
    clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void DOSCalculator::clear() {
    _qpoints.clear();
    //_qweights.clear();  OBSOLETE ME190423
    _freqs.clear();
    _bins.clear();
    _dos.clear();
    _idos.clear();  // ME190614
    _eigen.clear();  // ME190624
    _projectedDOS.clear(); // ME190614
    _projections.clear();  // ME190624
    _bzmethod = "";
    _temperature = 0.0;  // ME190614
  }

  // ///////////////////////////////////////////////////////////////////////////

  //CO - START
  void DOSCalculator::calculateInOneThread(int startIndex, int endIndex) {
    //cout << "Thread: from " << startIndex << " to " <<  endIndex << std::endl;
    for (int iqp = startIndex; iqp < endIndex; iqp++) {
      _logger.updateProgressBar(iqp, _qpoints.size());
      _freqs[iqp] = _pc.getFrequency(_qpoints[iqp], apl::THZ | apl::ALLOW_NEGATIVE, _eigen[iqp]);  // ME190624
      //std::this_thread::yield();
    }
  }
  //CO - END

  //////////////////////////////////////////////////////////////////////////////

  void DOSCalculator::calculateFrequencies() {
    // Get q-points for which to calculate the frequencies
    // ME190419 - BEGIN
    //_qpoints = _rg.getPoints();
    //_qweights = _rg.getWeights();
    _qpoints = _rg.getIrredQPointsCPOS();
    // ME190419 - END

    // 07-08-2010 We do not need it anymore, since mpmesh has all points in catesian coords now
    // Transform points to the form as it is expected by CPC:
    // Transform from the reciprocal lattice of the primitive cell to cartesian coordinates
    // for(uint t = 0; t < _qpoints.size(); t++)
    //    _qpoints[t] = trasp(ReciprocalLattice(_pc.getPrimitiveCellStructure().lattice)) * _qpoints[t];

    //CO - START
#ifdef AFLOW_APL_MULTITHREADS_ENABLE

    // Get the number of CPUS
    int ncpus; //= sysconf(_SC_NPROCESSORS_ONLN);  // AFLOW_MachineNCPUs;  //CO 180214
    _pc.get_NCPUS(ncpus);  //CO 180214
    if (ncpus < 1) ncpus = 1;
    //  int qpointsPerCPU = _qpoints.size() / ncpus;  OBSOLETE ME 180801

    // Show info
    if (ncpus == 1)
      _logger.initProgressBar("Calculating frequencies for DOS");
    else
      _logger.initProgressBar("Calculating frequencies for DOS (" + stringify(ncpus) + " threads)");

    // Prepare storage
    _freqs.clear();
    xvector<double> zero(_pc.getNumberOfBranches());
    for (uint i = 0; i < _qpoints.size(); i++)
      _freqs.push_back(zero);
    _eigen.resize(_qpoints.size(), xmatrix<xcomplex<double> >(_pc.getNumberOfBranches(), _pc.getNumberOfBranches()));  // ME190624

    // Distribute the calculation
    int startIndex, endIndex;
    std::vector<std::thread*> threads;
    vector<vector<int> > thread_dist = getThreadDistribution((int) _qpoints.size(), ncpus);
    for (int icpu = 0; icpu < ncpus; icpu++) {
      startIndex = thread_dist[icpu][0];
      endIndex = thread_dist[icpu][1];
      threads.push_back(new std::thread(&DOSCalculator::calculateInOneThread, this, startIndex, endIndex));
    }

    // OBSOLETE ME 180801
    //   for (int icpu = 0; icpu < ncpus; icpu++) {
    //   startIndex = icpu * qpointsPerCPU;
    //   endIndex = startIndex + qpointsPerCPU;
    //   if (((uint)endIndex > _qpoints.size()) ||
    //   ((icpu == ncpus - 1) && ((uint)endIndex < _qpoints.size())))
    //   endIndex = _qpoints.size();
    //   threads.push_back(new std::thread(&DOSCalculator::calculateInOneThread, this, startIndex, endIndex));
    //   }

    // Wait to finish all threads here!
    for (uint i = 0; i < threads.size(); i++) {
      threads[i]->join();
      delete threads[i];
    }
    threads.clear();

    // Done
    _logger.finishProgressBar();

#else

    // Calculate frequencies
    // ME190624 - added eigenvectors for projected DOS
    xmatrix<xcomplex<double> > xmtrx(_pc.getNumberOfBranches(), _pc.getNumberOfBranches());
    _logger.initProgressBar("Calculating frequencies for DOS");
    for (uint iqp = 0; iqp < _qpoints.size(); iqp++) {
      _logger.updateProgressBar(iqp, _qpoints.size());
      _freqs.push_back(_pc.getFrequency(_qpoints[iqp], apl::THZ | apl::ALLOW_NEGATIVE, xmtrx));
      _eigen.push_back(xmtrx);
    }
    _logger.finishProgressBar();

#endif
    //CO - END

    //if freq > MIN_FREQ_TRESHOLD considerd as +ve freq [PINKU]
    for (uint i = 0; i < _freqs.size(); i++) {
      for (int j = _freqs[i].lrows; j <= _freqs[i].urows; j++) {
        if ((_freqs[i][j] < 0.00) && (_freqs[i][j] > MIN_FREQ_TRESHOLD)) _freqs[i][j] = 0.00;
      }
    }
    //PINKU END

    // Get min and max values
    _maxFreq = -1.0;
    _minFreq = 0.0;
    for (uint i = 0; i < _freqs.size(); i++) {
      for (int j = _freqs[i].lrows; j <= _freqs[i].urows; j++) {
        if (_freqs[i](j) > _maxFreq) _maxFreq = _freqs[i](j);
        if (_freqs[i](j) < _minFreq) _minFreq = _freqs[i](j);
      }
    }
    _maxFreq += 1.0;
    if (_minFreq < MIN_FREQ_TRESHOLD) _minFreq -= 1.0;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME190614 - added integrated DOS
  void DOSCalculator::smearWithGaussian(vector<double>& dos, vector<double>& idos, double h, double sigma) {
    // Construct table for gaussian function
    int ng = (int)(6.0 * sigma / h + 1.0);
    double fact = 1.0 / (sqrt(2.0 * M_PI) * sigma);
    vector<double> gauss, igauss;
    double gnorm = 0.0;
    for (int ig = -ng; ig <= ng; ig++) {
      double eg = ig * h;
      double arg = eg * eg / (sigma * sigma) / 2.0;
      gauss.push_back(fact * exp(-arg));
      igauss.push_back(fact * erf(-arg));
      gnorm += gauss.back();
    }

    // Norm gauss table to one
    gnorm *= h;
    for (int ig = -ng; ig <= ng; ig++) {
      gauss[ig + ng] /= gnorm;
      igauss[ig + ng] /= gnorm;
    }

    // Prepare new dos
    vector<double> newdos;
    for (uint i = 0; i < dos.size(); i++) {
      newdos.push_back(0.0);
    }
    vector<double> newidos(newdos.size(), 0.0);

    // Convolute...
    for (int ie = 0; ie < (int)dos.size(); ie++) {
      double wt = dos[ie] * h;
      double wti = idos[ie] * h;
      for (int jg = -ng; jg <= ng; jg++) {
        int je = ie + jg;
        if (je < 0) continue;
        if (je >= (int)dos.size()) continue;
        newdos[je] += gauss[jg + ng] * wt;
        newidos[je] += igauss[jg + ng] * wti;
      }
    }

    //
    dos.clear();
    for (uint i = 0; i < newdos.size(); i++) {
      dos.push_back(newdos[i]);
    }
    newdos.clear();
    gauss.clear();
  }

  // ME200203 - DOS can now be calculated within any frequency range
  // ///////////////////////////////////////////////////////////////////////////

  void DOSCalculator::calc(int USER_DOS_NPOINTS) {
    calc(USER_DOS_NPOINTS, 0.0, _minFreq, _maxFreq);
  }

  // ///////////////////////////////////////////////////////////////////////////

  void DOSCalculator::calc(int USER_DOS_NPOINTS, double USER_DOS_SMEAR) {
    calc(USER_DOS_NPOINTS, USER_DOS_SMEAR, _minFreq, _maxFreq);
  }

  // ///////////////////////////////////////////////////////////////////////////
  void DOSCalculator::calc(int USER_DOS_NPOINTS, double USER_DOS_SMEAR,
      double fmin, double fmax) {
    // Check parameters
    if (aurostd::isequal(fmax, fmin, _AFLOW_APL_EPS_)) {
      string function = "apl::DOSCalculator::calc()";
      string message = "Frequency range of phonon DOS is nearly zero.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
    } else if (fmin > fmax) {
      double tmp = fmax;
      fmax = fmin;
      fmin = tmp;
    }
    // Calculate steps
    _stepDOS = (fmax - fmin) / (double)USER_DOS_NPOINTS;
    _halfStepDOS = 0.5 * _stepDOS;

    // Clear old stuff
    _dos.clear();
    _idos.clear();  // ME190614
    _bins.clear();

    // Prepare storagearrays
    for (int k = 0; k < USER_DOS_NPOINTS; k++) {
      _dos.push_back(0);
      _idos.push_back(0);  // ME190614
      _bins.push_back(fmin + k * _stepDOS + _halfStepDOS);
    }
    // ME190624
    if (_projections.size() > 0)
      _projectedDOS.resize(_pc.getInputCellStructure().atoms.size(),
          vector<vector<double> >(_projections.size(), vector<double>(USER_DOS_NPOINTS)));

    // Perform the raw specific calculation by method
    // ME190423 - START
    //rawCalc(USER_DOS_NPOINTS);  OBSOLETE
    if (_bzmethod == "LT") calcDosLT();
    else if (_bzmethod == "RS") calcDosRS();
    // ME190423 - END

    // Smooth DOS by gaussians
    if (USER_DOS_SMEAR > 1E-6)
      smearWithGaussian(_dos, _idos, _stepDOS, USER_DOS_SMEAR);  // ME190614

    // Normalize to number of branches
    double sum = 0.0;
    for (int k = 0; k < USER_DOS_NPOINTS; k++)
      sum += _dos[k];
    sum /= _pc.getNumberOfBranches();

    for (int k = 0; k < USER_DOS_NPOINTS; k++) {
      _dos[k] /= (sum * _stepDOS);
      _idos[k] /= (sum * _stepDOS);  // ME190614
    }
  }

  // ME190423 - START

  // ///////////////////////////////////////////////////////////////////////////

  void DOSCalculator::calcDosRS() {
    _logger << "Calculating phonon DOS using the root sampling method." << apl::endl;
    for (uint k = 0; k < _bins.size(); k++) {
      for (uint i = 0; i < _freqs.size(); i++) {
        for (int j = _freqs[i].lrows; j <= _freqs[i].urows; j++) {
          if ((_freqs[i][j] > (_bins[k] - _halfStepDOS)) &&
              (_freqs[i][j] < (_bins[k] + _halfStepDOS))) {
            _dos[k] += (double) _rg.getWeights()[i];
          }
        }
      }
    }
    // ME190614 - calculate integrated DOS
    _idos[0] = _dos[0];
    for (uint k = 1; k < _dos.size(); k++) {
      _idos[k] = _idos[k-1] + _dos[k];
    }
  }

// ///////////////////////////////////////////////////////////////////////////

  // ME190614 - added integrated DOS
  // ME190625 - rearranged and added projected DOS
  void DOSCalculator::calcDosLT() {
    _logger << "Calculating phonon DOS using the linear tetrahedron method." << apl::endl;
    // Procompute projections for each q-point and branch to save time
    uint nproj = _projections.size();
    vector<xvector<double> > proj_norm(nproj, xvector<double>(3));
    for (uint p = 0; p < nproj; p++) {
      proj_norm[p] = _projections[p]/aurostd::modulus(_projections[p]);
    }
    const xstructure& xstr = _pc.getInputCellStructure();
    uint natoms = xstr.atoms.size();
    vector<vector<vector<vector<double> > > > parts;
    if (nproj > 0) {
      // Precompute eigenvector projections
      xcomplex<double> eig;
      parts.assign(_rg.getnQPs(), vector<vector<vector<double> > >(_pc.getNumberOfBranches(), vector<vector<double> >(nproj, vector<double>(natoms, 0))));
      for (int q = 0; q < _rg.getnQPs(); q++) {
        for (uint br = 0; br < _pc.getNumberOfBranches(); br++) {
          int ibranch = br + _freqs[0].lrows;
          for (uint p = 0; p < nproj; p++) {
            for (uint at = 0; at < natoms; at++) {
              eig.re = 0.0;
              eig.im = 0.0;
              for (int i = 1; i < 4; i++) eig += proj_norm[p][i] * _eigen[q][3*at + i][ibranch];
              parts[q][br][p][at] = aurostd::magsqr(eig);
            }
          }
        }
      }
    }

    vector<double> f(4);
    double max_freq = _bins.back() + _halfStepDOS;
    double min_freq = _bins.front() - _halfStepDOS;
    vector<vector<int> > tet_corners;
    LTMethod _lt(_rg, _logger);
    if (nproj == 0) {
      _lt.makeIrreducible();
      tet_corners = _lt.getIrreducibleTetrahedraIbzqpt();
    } else {
      tet_corners = _lt.getTetrahedra();
    }
    for (int itet = 0; itet < _lt.getnIrredTetrahedra(); itet++) {
      double weightVolumeTetrahedron = _lt.getWeight(itet) * _lt.getVolumePerTetrahedron();
      const vector<int>& corners = tet_corners[itet];
      for (int ibranch = _freqs[0].lrows; ibranch <= _freqs[0].urows; ibranch++) {
        for (int icorner = 0; icorner < 4; icorner++) {
          f[icorner] = _freqs[corners[icorner]][ibranch];
        }
        std::sort(f.begin(), f.end());
        double fmin = f[0];
        double fmax = f[3];
        if (fmax > max_freq) continue;
        if (fmin < min_freq) continue;

        int kstart = (int) ((fmin - min_freq)/_stepDOS) + 1;
        if (kstart < 0) kstart = 0;
        if (kstart > (int) _bins.size()) kstart = _bins.size() - 1;
        int kstop = (int) ((fmax - min_freq)/_stepDOS) + 1;
        if (kstop < 0) kstop = 0;
        if (kstop < (int) _bins.size()) kstop = _bins.size() - 1;

        double f21 = f[1]  - f[0];
        double f31 = f[2]  - f[0];
        double f32 = f[2]  - f[1];
        double f41 = f[3]  - f[0];
        double f42 = f[3]  - f[1];
        double f43 = f[3]  - f[2];
        double cc12 = 3.0 * weightVolumeTetrahedron/(f21 * f31 * f41);
        double cc23 = weightVolumeTetrahedron/(f31 * f41);
        double cc23a = 3.0 * f21 * cc23;
        double cc23b = 6.0 * cc23;
        double cc23c = -3.0 * cc23 * (f31 + f42)/(f32 * f42);
        double cc34 = 3.0 * weightVolumeTetrahedron/(f41 * f42 * f43);

        double fbin, dos, part;
        int br = ibranch - _freqs[0].lrows;
        for (int k = kstart; k <= kstop; k++) {
          // ME200203 - Use bins to accommodate different frequency range
          fbin = _bins[k]; // _minFreq + k * _stepDOS + _halfStepDOS;
          dos = 0.0;
          if ((f[0] <= fbin) && (fbin <= f[1])) {
            double df = fbin - f[0];
            dos = cc12 * df * df;
            _idos[k] += cc12 * (df * df * df)/3.0;
          } else if ((f[1] < fbin) && (fbin <= f[2])) {
            double df = fbin - f[1];
            dos = cc23a + cc23b * df + cc23c * df * df;
            _idos[k] += cc23a * f21/3.0 + 3.0 * cc23 * f21 * df + 3.0 * cc23 * df * df + cc23c * (df * df * df)/3.0;
          } else if ((f[2] < fbin) && (fbin <= f[3])) {
            double df = f[3] - fbin;
            dos = cc34 * df * df;
            _idos[k] += weightVolumeTetrahedron + cc34 * (df * df * df)/3.0;
          } else if (f[3] < fbin) {
            _idos[k] += weightVolumeTetrahedron;
          }
          _dos[k] += dos;
          for (uint p = 0; p < nproj; p++) {
            for (uint at = 0; at < natoms; at++) {
              part = 0.0;
              for (int icorner = 0; icorner < 4; icorner++) {
                part += parts[corners[icorner]][br][p][at];
              }
              _projectedDOS[at][p][k] += 0.25 * dos * part;
            }
          }
        }
      }
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME190423 - END

  void DOSCalculator::writePDOS(const string& directory) {
    // Write PHDOS file
    //CO - START
    //ofstream outfile("PDOS",ios_base::out);
    stringstream outfile;
    //if( !outfile.is_open() )
    //{
    //    throw apl::APLRuntimeError("DOSCalculator::writePDOS(); Cannot open output PDOS file.");
    //}
    //CO - END

    string filename = directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_PDOS_FILE;
    double factorTHz2Raw = _pc.getFrequencyConversionFactor(apl::THZ, apl::RAW);
    double factorRaw2rcm = _pc.getFrequencyConversionFactor(apl::RAW, apl::RECIPROCAL_CM);
    double factorRaw2meV = _pc.getFrequencyConversionFactor(apl::RAW, apl::MEV);

    _logger << "Writing phonon density of states into file " << aurostd::CleanFileName(filename) << "." << apl::endl; //ME181226
    //outfile << "############### ############### ############### ###############" << std::endl;
    outfile << "#    f(THz)      1/lambda(cm-1)      E(meV)          pDOS      " << std::endl;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(8);
    for (uint k = 0; k < _dos.size(); k++) {
      outfile << setw(15) << _bins[k] << " "
        << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2rcm << " "
        << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2meV << " "
        << setw(15) << _dos[k] << std::endl;
    }

    //CO - START
    aurostd::stringstream2file(outfile, filename); //ME181226
    if (!aurostd::FileExist(filename)) { //ME181226
      string function = "DOSCalculator::writePDOS()";
      string message = "Cannot open output file " + filename + "."; //ME181226
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
      //    throw apl::APLRuntimeError("DOSCalculator::writePDOS(); Cannot open output PDOS file.");
    }
    //outfile.clear();
    //outfile.close();
    //CO - END
  }

  // ME190614 - writes phonon DOS in DOSCAR format
  void DOSCalculator::writePHDOSCAR(const string& directory) {
    string filename = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_PHDOSCAR_FILE);
    _logger << "Writing phonon density of states into file " << filename << "." << apl::endl;
    stringstream doscar;
    xDOSCAR xdos = createDOSCAR();
    doscar << xdos;
    aurostd::stringstream2file(doscar, filename);
    if (!aurostd::FileExist(filename)) {
      string function = "PhononDispersionCalculator::writePHDOSCAR()";
      string message = "Cannot open output file " + filename + ".";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
    }
    // OBSOLETE ME191219 - PHPOSCAR is already written in KBIN::RunPhonons_APL
    // if (xdos.partial) {  // Write PHPOSCAR if there are projected DOS
    //   filename = aurostd::CleanFileName(directory + "/" + DEFAULT_APL_PHPOSCAR_FILE);
    //   xstructure xstr = _pc.getInputCellStructure();
    //   xstr.is_vasp5_poscar_format = true;
    //   stringstream poscar;
    //   poscar << xstr;
    //   aurostd::stringstream2file(poscar, filename);
    //   if (!aurostd::FileExist(filename)) {
    //     string function = "PhononDispersionCalculator::writePHPOSCAR()";
    //     string message = "Cannot open output file " + filename + ".";
    //     throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
    //   }
  }

  xDOSCAR DOSCalculator::createDOSCAR() {
    xDOSCAR xdos;
    xdos.spin = 0;
    // Header values
    xdos.number_atoms = _pc.getInputCellStructure().atoms.size();
    xdos.partial = (_projectedDOS.size() > 0);
    xdos.Vol = GetVolume(_pc.getInputCellStructure())/xdos.number_atoms;
    xvector<double> lattice;
    lattice[1] = _pc.getInputCellStructure().a * 1E-10;
    lattice[2] = _pc.getInputCellStructure().b * 1E-10;
    lattice[3] = _pc.getInputCellStructure().c * 1E-10;
    xdos.lattice = lattice;
    xdos.POTIM = 0.5E-15;
    xdos.temperature = _temperature;
    xdos.carstring = "PHON";
    xdos.title = _system;

    // Data
    double factorTHz2Raw = _pc.getFrequencyConversionFactor(apl::THZ, apl::RAW);
    double factorRaw2meV = _pc.getFrequencyConversionFactor(apl::RAW, apl::MEV);
    double conv = factorTHz2Raw * factorRaw2meV/1000;
    // ME200203 - use _bins instead of _minFreq in case the DOS was calculated
    // using different frequency ranges
    xdos.energy_max = _bins.back() * conv;
    xdos.energy_min = _bins[0] * conv;
    xdos.number_energies = _dos.size();
    xdos.Efermi = 0.0;  // phonon DOS have no Fermi energy
    xdos.venergy = aurostd::vector2deque(_bins);
    for (uint i = 0; i < xdos.number_energies; i++) xdos.venergy[i] *= conv;
    xdos.viDOS.resize(1);
    xdos.viDOS[0] = aurostd::vector2deque(_idos);
    deque<deque<deque<deque<double> > > > vDOS;
    // ME190625
    if (_projections.size() > 0) {
      vDOS.resize(xdos.number_atoms + 1, deque<deque<deque<double> > >(_projections.size() + 1, deque<deque<double> >(1)));
    } else {
      vDOS.resize(1, deque<deque<deque<double> > >(1, deque<deque<double> >(1)));
    }

    vDOS[0][0][0] = aurostd::vector2deque<double>(_dos);

    // ME190624 - projected DOS
    if (_projections.size() > 0) {
      for (uint at = 0; at < xdos.number_atoms; at++) {
        for (uint p = 0; p < _projections.size(); p++) {
          vDOS[at + 1][p + 1][0] = aurostd::vector2deque(_projectedDOS[at][p]);
        }
      }
    }
    xdos.vDOS = vDOS;

    return xdos;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME200108 - added const
  vector<double> DOSCalculator::getBins() const {
    return _bins;
  }

  vector<double> DOSCalculator::getDOS() const {
    return _dos;
  }

  bool DOSCalculator::hasNegativeFrequencies() const {
    return (_minFreq < MIN_FREQ_TRESHOLD ? true : false);
  }

  // ///////////////////////////////////////////////////////////////////////////
  //PINKU - START
  void DOSCalculator::writePDOS(string path, string ex)  //[PINKU]
  {
    //CO - START
    // Write PHDOS file
    //ofstream outfile(file.c_str(),ios_base::out);
    stringstream outfile;
    //if( !outfile.is_open() )
    //{
    //    throw apl::APLRuntimeError("DOSCalculator::writePDOS(); Cannot open output PDOS file.");
    //}
    //CO - END

    double factorTHz2Raw = _pc.getFrequencyConversionFactor(apl::THZ, apl::RAW);
    double factorRaw2rcm = _pc.getFrequencyConversionFactor(apl::RAW, apl::RECIPROCAL_CM);
    double factorRaw2meV = _pc.getFrequencyConversionFactor(apl::RAW, apl::MEV);

    string filename = DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_PDOS_FILE; //ME181226
    _logger << "Writing phonon density of states into file " << filename << "." << apl::endl; //ME181226
    //outfile << "############### ############### ############### ###############" << std::endl;
    outfile << "#    f(THz)      1/lambda(cm-1)      E(meV)          pDOS      " << std::endl;
    outfile << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    outfile << setprecision(8);
    for (uint k = 0; k < _dos.size(); k++) {
      outfile << setw(15) << _bins[k] << " "
        << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2rcm << " "
        << setw(15) << _bins[k] * factorTHz2Raw * factorRaw2meV << " "
        << setw(15) << _dos[k] << std::endl;
    }

    //CO - START
    string file = path + "/" + filename + "." + ex; //ME181226
    aurostd::stringstream2file(outfile, file);
    if (!aurostd::FileExist(file)) {
      string function = "DOSCalculator::writePDOS()";
      string message = "Cannot open output file " + filename + "."; //ME181226
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
      //    throw apl::APLRuntimeError("DOSCalculator::writePDOS(); Cannot open output PDOS file.");
    }
    //outfile.clear();
    //outfile.close();
    //CO - END
  }
  //PINKU - END
  // ///////////////////////////////////////////////////////////////////////////

}  // namespace apl
