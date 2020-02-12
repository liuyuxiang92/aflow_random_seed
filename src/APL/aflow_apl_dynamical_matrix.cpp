#include "aflow_apl.h"

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  // ME180827 - Overloaded to calculate derivative and eigenvectors for AAPL
  // ME200206 - Added variants for the case near the Gamma point where the
  // non-analytical correction also needs a direction.
  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const IPCFreqFlags& flags) {
    return getFrequency(kpoint, kpoint, flags);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const xvector<double>& kpoint_nac, const IPCFreqFlags& flags) {
    const xstructure& pc = _supercell.getInputStructureLight();  //CO
    uint nBranches = 3 * pc.atoms.size();
    xmatrix<xcomplex<double> > placeholder_eigen(nBranches, nBranches, 1, 1);
    return getFrequency(kpoint, kpoint_nac, flags, placeholder_eigen);
  }

  // ME190624 - get eigenvectors and frequencies
  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const IPCFreqFlags& flags,
      xmatrix<xcomplex<double> >& eigenvectors) {
    return getFrequency(kpoint, kpoint, flags, eigenvectors);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const xvector<double>& kpoint_nac,
      const IPCFreqFlags& flags, xmatrix<xcomplex<double> >& eigenvectors) {
    vector<xmatrix<xcomplex<double> > > placeholder_mat;
    return getFrequency(kpoint, kpoint_nac, flags, eigenvectors, placeholder_mat, false);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint,
      const IPCFreqFlags& flags, xmatrix<xcomplex<double> >& eigenvectors,
      vector<xmatrix<xcomplex<double> > >& dDynMat, bool calc_derivative) {
    return getFrequency(kpoint, kpoint, flags, eigenvectors, dDynMat, calc_derivative);
  }

  xvector<double> PhononCalculator::getFrequency(const xvector<double>& kpoint, const xvector<double>& kpoint_nac,
      const IPCFreqFlags& flags, xmatrix<xcomplex<double> >& eigenvectors,
      vector<xmatrix<xcomplex<double> > >& dDynMat, bool calc_derivative) {
    // Compute frequency(omega) from eigenvalues [in eV/A/A/atomic_mass_unit]
    xvector<double> omega = getEigenvalues(kpoint, kpoint_nac, eigenvectors, dDynMat, calc_derivative);

    // Get value of conversion factor
    double conversionFactor = getFrequencyConversionFactor(apl::RAW | apl::OMEGA, flags);

    // Transform values to desired format
    for (_AFLOW_APL_REGISTER_ int i = omega.lrows; i <= omega.urows; i++) {
      if (omega(i) < 0) {
        if (flags & ALLOW_NEGATIVE)
          omega(i) = -sqrt(-omega(i));
        else
          omega(i) = 0.0;
      } else {
        omega(i) = sqrt(omega(i));
      }

      // Convert to desired units
      omega(i) *= conversionFactor;
    }

    // Return
    return (omega);
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME200108 - replaced with constants in xscalar
  double PhononCalculator::getFrequencyConversionFactor(IPCFreqFlags inFlags, IPCFreqFlags outFlags) {
    double conversionFactor = 1.0;

    // Conversion from eV/A/A/atomic_mass_unit -> something
    if (inFlags & apl::RAW) {
      if (outFlags & apl::RAW) {
        // Transform to eV/A/A/atomic_mass_unit
        conversionFactor = 1.0;
      } else if (outFlags & apl::HERTZ) {
        // Transform to s-1; sqrt(EV_TO_JOULE / (ANGSTROM_TO_METER*ANGSTROM_TO_METER) / AMU_TO_KG);
        conversionFactor = au2Hz;
      } else if (outFlags & apl::THZ) {
        // Transform to THz; (in Hertz) / 1E12;
        conversionFactor = au2Hz * Hz2THz;
      } else if (outFlags & apl::RECIPROCAL_CM) {
        // Transform to cm-1; 1/lambda(m) = freq.(s-1) / light_speed(m/s);
        conversionFactor = au2rcm;
      } else if (outFlags & apl::MEV) {
        // Transform to meV; E(eV) = h(eV.s) * freq(s-1); h[(from J.s) -> (eV.s)] = 4.1356673310E-15
        conversionFactor = 1000 * au2eV;
      } else {
        // ME191031 - use xerror
        //throw APLRuntimeError("apl::PhononCalculator:convertFrequencyUnit(); Not implemented conversion.");
        string function = "apl::PhononCalculator:convertFrequencyUnit()";
        string message = "Not implemented conversion.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
      }
    }

    // Conversion from THz -> something
    else if (inFlags & apl::THZ) {
      if (outFlags & apl::RAW) {
        // Transform to eV/A/A/atomic_mass_unit
        conversionFactor = 1.0 / (au2Hz * Hz2THz);
      } else if (outFlags & apl::THZ) {
        conversionFactor = 1.0;
      } else if (outFlags & apl::MEV) {
        conversionFactor = 1000 * PLANCKSCONSTANTEV_h * THz2Hz;
      } else {
        // ME191031 - use xerror
        //throw APLRuntimeError("apl::PhononCalculator:convertFrequencyUnit(); Not implemented conversion.");
        string function = "apl::PhononCalculator:convertFrequencyUnit()";
        string message = "Not implemented conversion.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
      }
    }

    // Nothing suits?
    else {
      // ME191031 - use xerror
      //throw APLRuntimeError("apl::PhononCalculator:convertFrequencyUnit(); Not implemented conversion.");
      string function = "apl::PhononCalculator:convertFrequencyUnit()";
      string message = "Not implemented conversion.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
    }

    //
    if ((outFlags & OMEGA) && !(inFlags & OMEGA))
      conversionFactor *= 2.0 * M_PI;
    if (!(outFlags & OMEGA) && (inFlags & OMEGA))
      conversionFactor /= 2.0 * M_PI;

    //
    return (conversionFactor);
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME180827 - Overloaded to calculate derivative and eigenvectors for AAPL
  // OBSOLETE ME200206 - not used anywhere and notuseful for debugging (use get Frequency)
  //[OBSOLETE] xvector<double> PhononCalculator::getEigenvalues(const xvector<double>& kpoint) {
  //[OBSOLETE]   const xstructure& pc = _supercell.getInputStructureLight();  //CO
  //[OBSOLETE]   uint nBranches = 3 * pc.atoms.size();
  //[OBSOLETE]   xmatrix<xcomplex<double> > placeholder_eigen(nBranches, nBranches, 1, 1);
  //[OBSOLETE]   vector<xmatrix<xcomplex<double> > > placeholder_mat;
  //[OBSOLETE]   return getEigenvalues(kpoint, kpoint, placeholder_eigen, placeholder_mat, false);
  //[OBSOLETE] }

  xvector<double> PhononCalculator::getEigenvalues(const xvector<double>& kpoint,
      const xvector<double>& kpoint_nac,
      xmatrix<xcomplex<double> >& eigenvectors,
      vector<xmatrix<xcomplex<double> > >& dDynMat,
      bool calc_derivative) {
    // Get dynamical matrix
    xmatrix<xcomplex<double> > dynamicalMatrix = getDynamicalMatrix(kpoint, kpoint_nac, dDynMat, calc_derivative);

    // Diagonalize
    xvector<double> eigenvalues(dynamicalMatrix.rows, 1);
    //    xmatrix<xcomplex<double> > unitaryMatrix;  OBSOLETE ME 180827

    // OBSOLETE ME190815 - moved to aurostd::xmatrix
    //#ifdef USE_MKL
    //    zheevMKL(dynamicalMatrix, eigenvalues, eigenvectors);
    //#else
    //    //tred2(dynamicalMatrix);
    //    zheevByJacobiRotation(dynamicalMatrix2, eigenvalues2, eigenvectors2);
    //    eigenvectors2 = trasp(eigenvectors2);
    //#endif

    // ME 180828; OBSOLETE ME190815 - use Jacobi algorithm in aurostd::xmatrix, which
    // is much, much faster than aplEigensystems for large systems
    //    apl::aplEigensystems e;
    //    e.eigen_calculation(dynamicalMatrix, eigenvalues, eigenvectors, APL_MV_EIGEN_SORT_VAL_ASC);

    eigenvalues = jacobiHermitian(dynamicalMatrix, eigenvectors);  // ME190815

    return eigenvalues;
  }

  //  // ///////////////////////////////////////////////////////////////////////////
  // ME180827 - Overloaded to calculate derivative for AAPL
  // ME200206 - Added variants for the case near the Gamma point where the
  // non-analytical correction also needs a direction. While dynamical matrices
  // are not used directly, these functions are helpful debugging tools.
  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint) {
    return getDynamicalMatrix(kpoint, kpoint);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint, const xvector<double>& kpoint_nac) {
    vector<xmatrix<xcomplex<double> > > placeholder;
    return getDynamicalMatrix(kpoint, kpoint_nac, placeholder, false);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getDynamicalMatrix(const xvector<double>& kpoint,
      const xvector<double>& kpoint_nac,
      vector<xmatrix<xcomplex<double> > >& dDynMat,
      bool calc_derivative) {
    uint scAtomsSize = _supercell.getSupercellStructure().atoms.size();
    uint pcAtomsSize = _supercell.getInputStructure().atoms.size();

    uint nBranches = 3 * pcAtomsSize;
    xmatrix<xcomplex<double> > dynamicalMatrix(nBranches, nBranches, 1, 1);
    xmatrix<xcomplex<double> > dynamicalMatrix0(nBranches, nBranches, 1, 1);

    xcomplex<double> phase;
    double value;
    // ME 180828 - Prepare derivative calculation
    xvector<xcomplex<double> > derivative(3);
    vector<xmatrix<xcomplex<double> > > dDynMat_NAC;
    if (calc_derivative) {  // reset dDynMat
      dDynMat.clear();
      xmatrix<xcomplex<double> > mat(nBranches, nBranches, 1, 1);
      dDynMat.assign(3, mat);
    }

    // Calculate nonanalytical contribution
    xmatrix<xcomplex<double> > dynamicalMatrixNA(nBranches, nBranches, 1, 1);
    if (_isPolarMaterial)
      dynamicalMatrixNA = getNonanalyticalTermWang(kpoint_nac, dDynMat_NAC, calc_derivative);

    // Loop over primitive cell
    xcomplex<double> nac;
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      uint isc1 = _supercell.pc2scMap(ipc1);

      for (uint isc2 = 0; isc2 < scAtomsSize; isc2++) {
        uint ipc2 = _supercell.sc2pcMap(isc2);
        int neq;  // Important for NAC derivative
        if (_supercell.calcShellPhaseFactor(isc2, isc1, kpoint, phase, neq, derivative, calc_derivative)) {  // ME180827
          for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
            for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
              value = 0.5 * (_forceConstantMatrices[isc1][isc2](ix, iy) + _forceConstantMatrices[isc2][isc1](iy, ix));
              dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += value * phase;
              if (_isPolarMaterial)
                dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += dynamicalMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy) * phase;
              dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy) += value;
              if (_isPolarMaterial)
                dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy) += dynamicalMatrixNA(3 * ipc1 + ix, 3 * ipc2 + iy);
              if (calc_derivative) {
                for (int d = 0; d < 3; d++) {
                  dDynMat[d](3 * ipc1 + ix, 3 * ipc2 + iy) += value * derivative[d+1];
                  if (_isPolarMaterial && (aurostd::modulus(kpoint) > _AFLOW_APL_EPS_)) {
                    nac = ((double) neq) * phase * dDynMat_NAC[d](3 * ipc1 + ix, 3 * ipc2 + iy);
                    dDynMat[d](3 * ipc1 + ix, 3 * ipc2 + iy) += nac;
                  }
                }
              }
            }
          }
        }
      }
    }
    //printXMatrix2(dynamicalMatrix);

    // Subtract the sum of all "forces" from the central atom, this is like an automatic sum rule...
    for (uint i = 0; i < pcAtomsSize; i++) {
      for (uint j = 0; j < pcAtomsSize; j++) {
        for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
          for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
            dynamicalMatrix(3 * i + ix, 3 * i + iy) = dynamicalMatrix(3 * i + ix, 3 * i + iy) - dynamicalMatrix0(3 * i + ix, 3 * j + iy);
          }
        }
      }
    }

    // Get correction for polar materials
    //if( _isPolarMaterial )
    // dynamicMatrix += getNonanalyticalTermGonze(kpoint);

    // Make it hermitian
    for (uint i = 0; i <= pcAtomsSize - 1; i++) {
      for (uint j = 0; j <= i; j++) {
        for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
          for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
            dynamicalMatrix(3 * i + ix, 3 * j + iy) += conj(dynamicalMatrix(3 * j + iy, 3 * i + ix));
            dynamicalMatrix(3 * i + ix, 3 * j + iy) *= 0.5;
            dynamicalMatrix(3 * j + iy, 3 * i + ix) = conj(dynamicalMatrix(3 * i + ix, 3 * j + iy));
            if (calc_derivative) {
              for (int d = 0; d < 3; d++) {
                dDynMat[d](3 * i + ix, 3 * j + iy) += conj(dDynMat[d](3 * j + iy, 3 * i + ix));
                dDynMat[d](3 * i + ix, 3 * j + iy) *= 0.5;
                dDynMat[d](3 * j + iy, 3 * i + ix) = conj(dDynMat[d](3 * i + ix, 3 * j + iy));
              }
            }
          }
        }
      }
    }

    // Divide by masses
    for (uint i = 0; i < pcAtomsSize; i++) {
      double mass_i = _supercell.getAtomMass(_supercell.pc2scMap(i));
      for (uint j = 0; j < pcAtomsSize; j++) {
        double mass_j = _supercell.getAtomMass(_supercell.pc2scMap(j));
        for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
          for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
            dynamicalMatrix(3 * i + ix, 3 * j + iy) *= 1.0 / sqrt(mass_i * mass_j);
            if (calc_derivative) {
              for (int d = 0; d < 3; d++) {
                dDynMat[d](3 * i + ix, 3 * j + iy) *= 1.0/sqrt(mass_i * mass_j);
              }
            }
          }
        }
      }
    }

    return dynamicalMatrix;
  }

 ///////////////////////////////////////////////////////////////////////////

  // Y. Wang et.al, J. Phys.:Condens. Matter 22, 202201 (2010)
  // DOI: 10.1088/0953-8984/22/20/202201

  // ME180827 - Overloaded to calculate derivative for AAPL
  // ME200207 - This function assummed that Born charges were stored for each type,
  // but it is actually stored for each iatom.
  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermWang(const xvector<double>& _q) {
    vector<xmatrix<xcomplex<double> > > placeholder;
    return getNonanalyticalTermWang(_q, placeholder, false);
  }

  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermWang(const xvector<double>& _q,
      vector<xmatrix<xcomplex<double> > >& derivative,
      bool calc_derivative) {
    const xstructure& sc = _supercell.getSupercellStructureLight();           //CO
    const xstructure& pc = _supercell.getInputStructure();  //CO  // ME200207 - grab input structure (need iatoms)

    // to correct the q=\Gamma as a limit
    xvector<double> q(_q);
    if (aurostd::modulus(q) < _AFLOW_APL_EPS_) {
      q(1) = _AFLOW_APL_EPS_ * 1.001;
    }

    uint pcAtomsSize = pc.atoms.size();
    uint pcIAtomsSize = pc.iatoms.size();
    uint nBranches = 3 * pcAtomsSize;

    xmatrix<xcomplex<double> > dynamicalMatrix(nBranches, nBranches);

    if (calc_derivative) {  // reset derivative
      derivative.clear();
      xmatrix<xcomplex<double> > mat(nBranches, nBranches, 1, 1);
      derivative.assign(3, mat);
    }

    // Calculation
    double fac0 = hartree2eV * bohr2angst;  // from a.u. to eV/A  // ME200206 - replaced with xscalar constants
    double volume = det(pc.lattice);
    double fac1 = 4.0 * PI / volume;
    double nbCells = det(sc.lattice) / volume;


    if (aurostd::modulus(q) > _AFLOW_APL_EPS_) {
      // Precompute product of q-point with charge tensor
      vector<xvector<double> > qZ(pcIAtomsSize);
      for (uint at = 0; at < pcIAtomsSize; at++) qZ[at] = q * _bornEffectiveChargeTensor[at];

      double dotprod = scalar_product(q, _dielectricTensor * q);
      double prefactor = fac0 * fac1/(dotprod * nbCells);
      for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
        int iat1 = pc.atoms[ipc1].index_iatoms;
        for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
          int iat2 = pc.atoms[ipc2].index_iatoms;
          for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
            for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
              //int typei = pc.atoms[ipc1].type;
              //int typej = pc.atoms[ipc2].type;
              //double borni = (q * _bornEffectiveChargeTensor[typei])(ix);
              //double bornj = (q * _bornEffectiveChargeTensor[typej])(iy);
              double borni = qZ[iat1][ix];
              double bornj = qZ[iat2][iy];
              dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) = prefactor * borni * bornj;
              if (calc_derivative) {
                for (int d = 0; d < 3; d++) {
                  xcomplex<double> coeff(0, 0);
                  //coeff += borni * _bornEffectiveChargeTensor[ipc1](iy, d + 1);
                  //coeff += bornj * _bornEffectiveChargeTensor[ipc2](ix, d + 1);
                  coeff += borni * _bornEffectiveChargeTensor[iat1](iy, d + 1);
                  coeff += bornj * _bornEffectiveChargeTensor[iat2](ix, d + 1);
                  coeff -= 2 * borni * bornj * scalar_product(_dielectricTensor(d + 1), q)/dotprod;
                  derivative[d](3 * ipc1 + ix, 3 * ipc2 + iy) = prefactor * coeff;
                }
              }
            }
          }
        }
      }
    }

    //
    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // X. Gonze et al., Phys. Rev. B 50, 13035 (1994)
  // X. Gonze and Ch. Lee, Phys. Rev. B 55, 10355 (1997)

  xmatrix<xcomplex<double> > PhononCalculator::getNonanalyticalTermGonze(const xvector<double> kpoint) {
    uint pcAtomsSize = _supercell.getInputStructure().atoms.size();
    //    uint nBranches = 3 * pcAtomsSize; // not needed

    if (!_isGammaEwaldPrecomputed) {
      xvector<double> zero(3);
      xmatrix<xcomplex<double> > dynamicalMatrix0(getEwaldSumDipolDipolContribution(zero, false));

      _gammaEwaldCorr.clear();
      for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
        xmatrix<xcomplex<double> > sum(3, 3);
        for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
          for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
            for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
              sum(ix, iy) += dynamicalMatrix0(3 * ipc1 + ix, 3 * ipc2 + iy);
        }
        _gammaEwaldCorr.push_back(sum);
      }

      _isGammaEwaldPrecomputed = true;
    }

    //
    xmatrix<xcomplex<double> > dynamicalMatrix(getEwaldSumDipolDipolContribution(kpoint));

    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
        for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
          dynamicalMatrix(3 * ipc1 + ix, 3 * ipc1 + iy) -= _gammaEwaldCorr[ipc1](ix, iy);
    }

    //
    return dynamicalMatrix;
  }

  // ///////////////////////////////////////////////////////////////////////////

  // ME200207 - This function assummed that Born charges were stored for each type,
  // but it is actually stored for each iatom.
  xmatrix<xcomplex<double> > PhononCalculator::getEwaldSumDipolDipolContribution(const xvector<double> qpoint, bool includeTerm1) {
    // Definitions
    const xstructure& sc = _supercell.getSupercellStructureLight();           //CO
    const xstructure& pc = _supercell.getInputStructure();  //CO  // ME200207 - grab input structure (need iatoms)

    uint pcAtomsSize = pc.atoms.size();
    uint nBranches = 3 * pcAtomsSize;

    xmatrix<xcomplex<double> > dynamicalMatrix(nBranches, nBranches);

    double gmax = 14.0;
    double lambda = 1.0;
    double lambda2 = lambda * lambda;
    double lambda3 = lambda2 * lambda;
    double geg = gmax * lambda2 * 4.0;

    // Reciprocal Space
    xmatrix<double> klattice = trasp(ReciprocalLattice(pc.lattice));

    // Grid
    int n1 = (int)(sqrt(geg) / aurostd::modulus(klattice(1))) + 1;
    int n2 = (int)(sqrt(geg) / aurostd::modulus(klattice(2))) + 1;
    int n3 = (int)(sqrt(geg) / aurostd::modulus(klattice(3))) + 1;

    // Calculation
    double fac0 = hartree2eV * bohr2angst;  // from a.u. to eV/A  // ME200207 - replaced with xscalar constants
    double SQRTPI = sqrt(PI);
    double volume = det(pc.lattice);
    double fac = 4.0 * PI / volume;
    xcomplex<double> iONE(0.0, 1.0);

    // Term 1 - Reciprocal space sum

    if (includeTerm1) {
      for (_AFLOW_APL_REGISTER_ int m1 = -n1; m1 <= n1; m1++) {
        for (_AFLOW_APL_REGISTER_ int m2 = -n2; m2 <= n2; m2++) {
          for (_AFLOW_APL_REGISTER_ int m3 = -n3; m3 <= n3; m3++) {
            xvector<double> g = m1 * klattice(1) + m2 * klattice(2) + m3 * klattice(3) + qpoint;

            geg = scalar_product(g, _dielectricTensor * g);

            if (aurostd::abs(geg) > _AFLOW_APL_EPS_ && geg / lambda2 / 4.0 < gmax) {
              double fac2 = fac * exp(-geg / lambda2 / 4.0) / geg;

              for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
                //xvector<double> zag = g * _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
                int iat1 = pc.atoms[ipc1].index_iatoms;
                xvector<double> zag = g * _bornEffectiveChargeTensor[iat1];

                for (uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
                  //xvector<double> zbg = g * _bornEffectiveChargeTensor[pc.atoms[ipc2].type];
                  int iat2 = pc.atoms[ipc2].index_iatoms;
                  xvector<double> zbg = g * _bornEffectiveChargeTensor[iat2];

                  //xcomplex<double> e;
                  //(void)_supercell.calcShellPhaseFactor(ipc2,ipc1,g,e);
                  //xcomplex<double> facg = fac2 * e;
                  xcomplex<double> facg = fac2 * exp(iONE * scalar_product(g, sc.atoms[ipc2].cpos - sc.atoms[ipc1].cpos));

                  for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++) {
                    for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++) {
                      dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) += fac0 * facg * zag(ix) * zbg(iy);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // Term 2 - Real space sum
    //for(_AFLOW_APL_REGISTER_ int m1 = -n1; m1 <= n1; m1++)
    //  for(_AFLOW_APL_REGISTER_ int m2 = -n2; m2 <= n2; m2++)
    //    for(_AFLOW_APL_REGISTER_ int m3 = -n2; m3 <= n3; m3++) {
    //      xvector<double> rc = m1 * pc.lattice(1) + m2 * pc.lattice(2)
    //        + m3 * pc.lattice(3);

    //      //xvector<double> zero(3);
    //      //xvector<double> rf = _supercell.getFPositionItsNearestImage(rc,zero,pc.lattice);
    //      //rc = F2C(pc.lattice,rf);

    //      if( aurostd::modulus(rc) < _AFLOW_APL_EPS_ ) continue;

    //      //
    //      xvector<double> delta = _inverseDielectricTensor * rc;
    //      double D = sqrt( scalar_product(delta,rc) );

    //      //
    //      xmatrix<double> H(3,3);
    //      xvector<double> x = lambda * delta;
    //      double y = lambda * D;
    //      double y2 = y * y;
    //      double ym2 = 1.0 / y2;
    //      double emy2dpi = 2.0 * exp( -y2 ) / SQRTPI;
    //      double erfcdy = erfc(y) / y;
    //      double c1 = ym2 * ( 3.0 * erfcdy * ym2 + ( emy2dpi * ( 3.0 * ym2 + 2.0 ) ) );
    //      double c2 = ym2 * ( erfcdy + emy2dpi );
    //      for(_AFLOW_APL_REGISTER_ int a = 1; a <= 3; a++)
    //        for(_AFLOW_APL_REGISTER_ int b = 1; b <= 3; b++) {
    //          H(a,b) = x(a) * x(b) * c1 - _inverseDielectricTensor(a,b) * c2;
    //        }

    //      //
    //      xcomplex<double> e = exp( iONE * scalar_product(qpoint,rc) );
    //      xcomplex<double> fac = fac0 * lambda3 * _recsqrtDielectricTensorDeterminant * e;

    //      //
    //      for(uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
    //        xmatrix<double> zh = _bornEffectiveChargeTensor[pc.atoms[ipc1].type] * H;

    //        for(uint ipc2 = 0; ipc2 < pcAtomsSize; ipc2++) {
    //          xmatrix<double> zhz = zh * _bornEffectiveChargeTensor[pc.atoms[ipc2].type];

    //          for(_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
    //            for(_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
    //              dynamicalMatrix(3*ipc1+ix,3*ipc2+iy) -= fac * zhz(ix,iy);
    //        }
    //      }
    //    }

    // Term 2
    uint scAtomsSize = sc.atoms.size();
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      uint isc1 = _supercell.pc2scMap(ipc1);

      for (uint isc2 = 0; isc2 < scAtomsSize; isc2++) {
        uint ipc2 = _supercell.sc2pcMap(isc2);

        xvector<double> rf = _supercell.getFPositionItsNearestImage(isc2, isc1);
        xvector<double> rc = F2C(sc.lattice, rf);

        if (aurostd::modulus(rc) < _AFLOW_APL_EPS_) continue;

        //
        xvector<double> delta = _inverseDielectricTensor * rc;
        double D = sqrt(scalar_product(delta, rc));

        //
        xmatrix<double> H(3, 3);
        xvector<double> x = lambda * delta;
        double y = lambda * D;
        double y2 = y * y;
        double ym2 = 1.0 / y2;
        double emy2dpi = 2.0 * exp(-y2) / SQRTPI;
        double erfcdy = erfc(y) / y;
        double c1 = ym2 * (3.0 * erfcdy * ym2 + (emy2dpi * (3.0 * ym2 + 2.0)));
        double c2 = ym2 * (erfcdy + emy2dpi);
        for (_AFLOW_APL_REGISTER_ int a = 1; a <= 3; a++) {
          for (_AFLOW_APL_REGISTER_ int b = 1; b <= 3; b++) {
            H(a, b) = x(a) * x(b) * c1 - _inverseDielectricTensor(a, b) * c2;
          }
        }

        //xmatrix<double> za = _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
        //xmatrix<double> zb = _bornEffectiveChargeTensor[pc.atoms[ipc2].type];
        int iat1 = pc.atoms[ipc1].index_iatoms;
        int iat2 = pc.atoms[ipc2].index_iatoms;
        xmatrix<double> za = _bornEffectiveChargeTensor[iat1];
        xmatrix<double> zb = _bornEffectiveChargeTensor[iat2];
        xmatrix<double> zhz = za * H * zb;

        //
        xcomplex<double> e;  // = exp( iONE * scalar_product(qpoint,rc) );
        (void)_supercell.calcShellPhaseFactor(isc2, isc1, qpoint, e);

        //
        xcomplex<double> fac = fac0 * lambda3 * _recsqrtDielectricTensorDeterminant * e;
        for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
          for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
            dynamicalMatrix(3 * ipc1 + ix, 3 * ipc2 + iy) -= fac * zhz(ix, iy);
      }
    }

    // Term 3 - Limiting contribution

    double facterm3 = fac0 * 4.0 * lambda3 * _recsqrtDielectricTensorDeterminant / (3.0 * SQRTPI);
    for (uint ipc1 = 0; ipc1 < pcAtomsSize; ipc1++) {
      //xmatrix<double> z = _bornEffectiveChargeTensor[pc.atoms[ipc1].type];
      int iat1 = pc.atoms[ipc1].index_iatoms;
      xmatrix<double> z = _bornEffectiveChargeTensor[iat1];
      xmatrix<double> zez = z * _inverseDielectricTensor * z;

      for (_AFLOW_APL_REGISTER_ int ix = 1; ix <= 3; ix++)
        for (_AFLOW_APL_REGISTER_ int iy = 1; iy <= 3; iy++)
          dynamicalMatrix(3 * ipc1 + ix, 3 * ipc1 + iy) -= facterm3 * zez(ix, iy);
    }

    //
    return dynamicalMatrix;
  }

}  // namespace apl
