//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2017                    *
// *                                                                         *
//****************************************************************************
// Written by Marco Esters, 2018. Based on work by Jose J. Plata (AFLOW AAPL,
// DOI: 10.1038/s41524-017-0046-7) and Jesus Carrete (ShengBTE, 
// DOI: 10.1016/j.cpc.2014.02.015).
//
// This class calculates the anharmonic interatomic force constants (IFCs). It
// reads the forces of the inequivalent distortions first, and then calculates
// the forces of the equivalent distortions. These forces are used to
// calculate the IFCs of the inequivalent clusters that are then symmetrized
// in an iterative procedure using sum rules.
//
// See aflow_apl.h for descriptions of the classes and their members, and for
// the struct _linearCombinations.

#include "aflow_apl.h"

#define _DEBUG_AAPL_IFCS_ false

using std::vector;
using aurostd::xcombos;
using aurostd::xerror;

static const string _AAPL_IFCS_ERR_PREFIX_ = "apl::AnharmonicIFCs::";

// tform represents a tensor transformation containing the index and the
// coefficients. vector<vector<int> > holds the indices, vector<double>
// the coefficients.
typedef vector<std::pair<vector<int>, vector<double> > > tform;
// v4int defined for brevity
typedef vector<vector<vector<vector<int> > > > v4int;

// Coefficients for the finite difference method
static const double C12[2][3] = {{0.5,  0.0, -0.5},
                                 {1.0, -2.0,  1.0}};
static const double C3[5] = {-0.5, 1.0, 0.0, -1.0, 0.5};

/************************ CONSTRUCTORS/DESTRUCTOR ***************************/

namespace apl {

//Constructors////////////////////////////////////////////////////////////////
// Default constructor
AnharmonicIFCs::AnharmonicIFCs(vector<_xinput>& xInp,
                               ClusterSet& _clst,
                               const double& dist_mag, 
                               const aurostd::xoption& options,  // ME190501
                               Logger& l, _aflags& a) : clst(_clst), _logger(l), aflags(a) {
  free();
  distortion_magnitude = dist_mag;
  max_iter = aurostd::string2utype<int>(options.getattachedscheme("MAX_ITER"));  // ME190501
  order = clst.order;
  mixing_coefficient = aurostd::string2utype<double>(options.getattachedscheme("MIXING_COEFFICIENT"));  // ME190501
  sumrule_threshold = aurostd::string2utype<double>(options.getattachedscheme("THRESHOLD"));  // ME190501
  cart_indices = getCartesianIndices();

  vector<vector<vector<xvector<double> > > > force_tensors;
  force_tensors = storeForces(xInp);

  _logger << "Caulating anharmonic IFCs." << apl::endl;
  vector<vector<double> > ifcs_unsym = calculateUnsymmetrizedIFCs(clst.ineq_distortions, force_tensors);
  force_tensors.clear();

  _logger << "Symmetrizing IFCs." << apl::endl;
  force_constants = symmetrizeIFCs(ifcs_unsym);
}

// From file
AnharmonicIFCs::AnharmonicIFCs(const string& filename,
                               ClusterSet& _clst,
                               const double& dist_mag,
                               const aurostd::xoption& options,  // ME190501
                               Logger& l, _aflags& a) : clst(_clst), _logger(l), aflags(a) {
  free();
  distortion_magnitude = dist_mag;
  max_iter = aurostd::string2utype<int>(options.getattachedscheme("MAX_ITER"));  // ME190501
  order = clst.order;
  mixing_coefficient = aurostd::string2utype<double>(options.getattachedscheme("MIXING_COEFFICIENT"));  // ME190501
  sumrule_threshold = aurostd::string2utype<double>(options.getattachedscheme("THRESHOLD"));  // ME190501
  readIFCsFromFile(filename);
}

//Copy constructor
const AnharmonicIFCs& AnharmonicIFCs::operator=(const AnharmonicIFCs& that) {
  if (this != &that) {
    _logger = that._logger;
    aflags = that.aflags;
    cart_indices = that.cart_indices;
    clst = that.clst;
    distortion_magnitude = that.distortion_magnitude;
    force_constants = that.force_constants;
    max_iter = that.max_iter;
    mixing_coefficient = that.mixing_coefficient;
    order = that.order;
    sumrule_threshold = that.sumrule_threshold;
  }
  return *this;
}

//Destructor//////////////////////////////////////////////////////////////////
AnharmonicIFCs::~AnharmonicIFCs() {
  free();
}

//free////////////////////////////////////////////////////////////////////////
// Clears all vectors and resets all values.
void AnharmonicIFCs::free() {
  cart_indices.clear();
  distortion_magnitude = 0.0;
  force_constants.clear();
  max_iter = 0;
  mixing_coefficient = 0.0;
  order = 0;
  sumrule_threshold = 0.0;
}

}  // namespace apl

/*************************** INITIAL CALCULATIONS ***************************/

namespace apl {

//getCartesianIndices/////////////////////////////////////////////////////////
// Returns a list of Cartesian indices. Since they are looped over frequently,
// it is quicker to calculate them once at the beginning.
vector<vector<int> > AnharmonicIFCs::getCartesianIndices() {
  vector<vector<int> > indices;
  xcombos crt_ind(3, order, 'E', true);
  while (crt_ind.increment()) {
    indices.push_back(crt_ind.getCombo());
  }
  return indices;
}

// BEGIN Forces
//storeForces/////////////////////////////////////////////////////////////////
// Stores the forces from the VASP calculations. Each item in the vector holds
// the force tensor for a set of distorted atoms.
vector<vector<vector<xvector<double> > > > AnharmonicIFCs::storeForces(vector<_xinput>& xInp) {
  vector<vector<vector<xvector<double> > > > force_tensors(clst.ineq_distortions.size());
  int idxRun = 0;
  for (uint id = 0; id < clst.ineq_distortions.size(); id++) {
    force_tensors[id] = getForces(id, idxRun, xInp);
  }
  if (clst.order == 4) addHigherOrderForces(force_tensors, idxRun, xInp);
  return force_tensors;
}

//getForces///////////////////////////////////////////////////////////////////
// Retrieves all forces from the calculations. Also transforms the forces
// into the forces of the equivalent distortions.
vector<vector<xvector<double> > > AnharmonicIFCs::getForces(int id, int& idxRun,
                                                            vector<_xinput>& xInp) {
  const _ineq_distortions& ineq_dists = clst.ineq_distortions[id];
  int natoms = (int) clst.scell.atoms.size();
  vector<int> powers(order - 1, 1);
  int ndist = 1;
  for (int i = 0; i < order - 1; i++) {
    ndist *= 6;
    for (int j = 0; j < order - 2 - i; j++) {
      powers[i] *= 6;
    }
  }
  vector<vector<xvector<double> > > force_tensor(ndist, vector<xvector<double> >(natoms, xvector<double>(3)));

  int attrans = 0, fg = 0, index = 0;
  for (uint ineq = 0; ineq < ineq_dists.distortions.size(); ineq++) {
    // For the inequivalent distortion, just read the forces from VASP 
    const vector<xvector<double> >& qmforces = xInp[idxRun].getXStr().qm_forces;
    index = 0;
    for (int ind = 0; ind < order - 1; ind++) {
      index += powers[ind] * ineq_dists.distortions[ineq][0][ind];
    }
    for (int at = 0; at < natoms; at++) { 
      force_tensor[index][at] = qmforces[at];
    }
    for (uint i = 1; i < ineq_dists.distortions[ineq].size(); i++) {
      // For each equivalent distortion, transform the forces using symmetry
      vector<xvector<double> > qmforces_trans(natoms, xvector<double>(3));
      index = 0;
      for (int ind = 0; ind < order - 1; ind++) {
        index += powers[ind] * ineq_dists.distortions[ineq][i][ind];
      }
      fg = ineq_dists.rotations[ineq][i];
      const vector<int>& transformation_map = ineq_dists.transformation_maps[ineq][i];
      for (int at = 0; at < natoms; at++) {
        attrans = getTransformedAtom(transformation_map, at);
        force_tensor[index][at] = clst.pcell.fgroup[fg].Uc * qmforces[attrans];
      }
    }
    idxRun++;
  }
  return force_tensor;
}

void AnharmonicIFCs::addHigherOrderForces(vector<vector<vector<xvector<double> > > >& force_tensor,
                                          int& idxRun, vector<_xinput>& xInp) {
  const vector<_ineq_distortions>& ineq_dists = clst.higher_order_ineq_distortions;
  uint ndist = force_tensor[0].size();
  uint natoms = clst.scell.atoms.size();
  vector<xvector<double> > forces(natoms, xvector<double>(3));
  for (uint ineq = 0; ineq < ineq_dists.size(); ineq++) {
    uint idist;
    int at = ineq_dists[ineq].atoms[0];
    for (idist = 0; idist < clst.ineq_distortions.size(); idist++) {
      int a;
      for (a = 0; a < clst.order - 1; a++) {
        if (clst.ineq_distortions[idist].atoms[a] != at) break;
      }
      if (a == clst.order - 1) break;
    }
    for (int i = 0; i < 6; i++) force_tensor[idist].push_back(forces);

    for (uint dist = 0; dist < ineq_dists[ineq].distortions.size(); dist++) {
      const vector<xvector<double> >& qmforces = xInp[idxRun].getXStr().qm_forces;
      int d = ineq_dists[ineq].distortions[dist][0][0];
      force_tensor[idist][ndist + d][at] = qmforces[at];
      int fg = 0;
      for (uint i = 1; i < ineq_dists[ineq].distortions[dist].size(); i++) {
        d = ineq_dists[ineq].distortions[dist][i][0];
        fg = ineq_dists[ineq].rotations[dist][i];
        force_tensor[idist][ndist + d][at] = clst.pcell.fgroup[fg].Uc * qmforces[at];
      }
      xInp[idxRun].clear();
      idxRun++;
    }
  }
}

//getTransformedAtom//////////////////////////////////////////////////////////
// Retrieves the atom that is obtained by the transformation in the given
// symmetry map.
int AnharmonicIFCs::getTransformedAtom(const vector<int>& symmap, const int& at) {
  for (uint i = 0; i < symmap.size(); i++) {
    if (symmap[i] == at) {
      return i;
    }
  }
  // If the for-loop runs until the end, the atom was not found
  string function = _AAPL_IFCS_ERR_PREFIX_ + "getTransformedAtom";
  stringstream message;
  message << "Could not transform atom " << at;
  throw xerror(_AFLOW_FILE_NAME_,function, message, _RUNTIME_ERROR_);
}
//END Forces

// BEGIN Calculate unsymmetrized force constants
//calculateUnsymmetrizedIFCs//////////////////////////////////////////////////
// Calculates the IFCs of the inequivalent clusters from the forces.
vector<vector<double> >
    AnharmonicIFCs::calculateUnsymmetrizedIFCs(const vector<_ineq_distortions>& idist,
                                               const vector<vector<vector<xvector<double> > > >& forces) {
  vector<vector<double> > ifcs(clst.clusters.size(), vector<double>(cart_indices.size()));
  int at = 0, cl = 0, ic = 0;
  double denom = std::pow(distortion_magnitude, order - 1);
  for (uint f = 0; f < forces.size(); f++) {
    for (uint c = 0; c < idist[f].clusters.size(); c++) {
      ic = idist[f].clusters[c];
      cl = clst.ineq_clusters[ic][0];
      at = clst.getCluster(cl).atoms[order - 1];
      for (int cart = 0; cart < clst.nifcs; cart++) {
        //ifcs[cl][cart] = calculateIFC(forces[f], at, cart_indices[cart], idist[f].atoms);
        ifcs[cl][cart] = finiteDifference(forces[f], at, cart_indices[cart], idist[f].atoms)/denom;
      }
    }
  }
  return ifcs;
}

//calculateIFC////////////////////////////////////////////////////////////////
// Calculates a specific IFC from the forces using the central difference
// method. The numerator and denominator are calculated separately.
//
// The numerator is a sum of the forces for the set of atoms distorted along
// the Cartesian indices into positive and negative directions. The sign of
// the force depends on the number of negative distortions in the set.
//
// The denominator is the product of the distortion lengths times two (the
// factor of two comes from the central difference method). If the same atom
// gets distorted into the same direction multiple times, the distortion
// length needs to be adjusted so that the total length is the same as the
// distortion magnitude chosen by the user.
double AnharmonicIFCs::finiteDifference(const vector<vector<xvector<double> > >& forces, int at,
                                        const vector<int>& cart_ind, const vector<int>& atoms) {
  vector<int> powers(order - 1, 1);
  for (int i = 0; i < order - 2; i++) {
    for (int j = 0; j < order - 2 - i; j++) {
      powers[i] *= 6;
    }
  }

  vector<int> derivatives;
  int count = 0;
  for (uint a = 0; a < atoms.size(); a++) {
    count = 1;
    while (((int) a + count < order - 1) &&
           (atoms[a + count] == atoms[a]) &&
           (cart_ind[a + count] == cart_ind[a])) {
      count++;
    }
    derivatives.push_back(count);
    a += count - 1;
  }
  
  double diff = 0.0;
  vector<int> atm(order - 1);
  if (derivatives[0] == 3) {
    vector<int> dists(5);
    int pwr = 0;
    for (int i = 0; i < order - 1; i++) pwr += powers[i];
    dists[0] = clst.nifcs + cart_ind[0];
    dists[1] = cart_ind[0] * pwr;
    // No need to occupy dists[2] because C3[2] is zero
    dists[3] = (cart_ind[0] + 3) * pwr;
    dists[4] = clst.nifcs + cart_ind[0] + 3;
    for (int i = 0; i < 5; i++) {
      if (C3[i] != 0.0) diff -= C3[i] * forces[dists[i]][at][cart_ind[order - 1] + 1];
    }
  } else {
    uint nder = derivatives.size();
    vector<int> index(nder);
    xcombos ind(3, nder, 'E', true);
    double coeff = 0.0;
    int d = 0, dist = 0;
    while (ind.increment()) {
      coeff = 1.0;
      index = ind.getCombo();
      for (uint i = 0; i < nder; i++) {
        coeff *= C12[derivatives[i] - 1][index[i]];
      }
      if (coeff != 0.0) {
        dist = 0;
        d = 0;
        for (uint i = 0; i < nder; i++) {
          for (int j = 0; j < derivatives[i]; j++) {
            if (index[i] == 1) dist += powers[d] * (cart_ind[d] + 3 * j);
            else dist += powers[d] * (cart_ind[d] + 3 * index[i]/2);
            d++;
          }
        }
        diff -= coeff * forces[dist][at][cart_ind[order - 1] + 1];
      }
    }
  }
  return diff;
}
// END Calculate unsymmetrized force constants

}  // namespace apl

/****************************** SYMMETRIZATION ******************************/

namespace apl {

//symmetrizeIFCs//////////////////////////////////////////////////////////////
// Symmetrizes the IFCs using an iterative procedure to ensure 
// that the force constants sum to zero.
// 1. The IFCs of the inequivalent clusters will be symmetrized according to
//    the linear combinations found while determining ClusterSets.
// 2. The force constants will be transformed to the other clusters using the
//    symmetry of the crystal.
// 3. Determine the deviations of the IFC sums from zero.
// 4. If at least one deviation is above the chosen threshold, correct the
//    linearly dependent IFCs. If none are above the threshold or too many
//    iterations have been performed, stop the iteration procedure.
//
// Check the typedefs at the beginning of the file for tform and v4int
vector<vector<double> > AnharmonicIFCs::symmetrizeIFCs(vector<vector<double> > ifcs) {
  // Initialize tensors
  vector<vector<int> > reduced_clusters = getReducedClusters();
  vector<vector<double> > dev_from_zero(reduced_clusters.size(), vector<double>(cart_indices.size()));
  vector<vector<double> > abssum = dev_from_zero;

  // Tensor transformations
  vector<vector<tform> > transformations(clst.ineq_clusters.size());
  v4int eq_ifcs(clst.ineq_clusters.size());
  getTensorTransformations(eq_ifcs, transformations);

  // Do iterations
  int num_iter = 0;
  double max_err = 0.0;
  _logger << "Begin SCF for anharmonic force constants." << apl::endl;
  std::cout << std::setiosflags(std::ios::fixed | std::ios::right);
  std::cout << std::setw(15) << "Iteration";
  std::cout << std::setiosflags(std::ios::fixed | std::ios::right);
  std::cout << std::setw(20) << "Abs. max. error" << std::endl;
  do {
    // 1. Symmetrize using linear combinations
    applyLinCombs(ifcs);  

    // 2. Transform IFCs using symmetry and permutations
    transformIFCs(transformations, ifcs);

    // 3. Determine deviations from zero
    calcSums(reduced_clusters, ifcs, dev_from_zero, abssum);

    max_err = -1E30;
    for (uint i = 0; i < dev_from_zero.size(); i++) {
      for (uint j = 0; j < dev_from_zero[i].size(); j++) {
        if (std::abs(dev_from_zero[i][j]) > max_err) max_err = std::abs(dev_from_zero[i][j]);
      }
    }

    std::cout << std::setiosflags(std::ios::fixed | std::ios::right);
    std::cout << std::setw(15) << num_iter;
    std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
    std::cout << std::setw(20) << max_err << std::endl;

    // 4. Correct IFCs
    if (max_err > sumrule_threshold) {
      correctIFCs(ifcs, dev_from_zero, abssum, reduced_clusters, eq_ifcs);
    }
    num_iter++;
  } while ((num_iter <= max_iter) && (max_err > sumrule_threshold));
  _logger << "End SCF for anharmonic force constants." << apl::endl;
  if (num_iter > max_iter) {
    string function = _AAPL_IFCS_ERR_PREFIX_ + "symmetrizeIFCs";
    stringstream message;
    message << "Anharmonic force constants did not converge within " << max_iter << " iterations.";
    throw xerror(_AFLOW_FILE_NAME_,function, message, _RUNTIME_ERROR_);
  } else {
    return ifcs;
  }
}

// BEGIN Initializers
//getReducedClusters//////////////////////////////////////////////////////////
// Determines, for each set of inequivalent clusters, a uinque set of clusters
// that do not contain the last atom of the clusters. This set is important
// for the sum rules as they frequently require summations over the last atom
// within a set of inequivalent clusters.
vector<vector<int> > AnharmonicIFCs::getReducedClusters() {
  vector<vector<int> > reduced_clusters;
  vector<int> cluster(order - 1);
  for (uint c = 0; c < clst.clusters.size(); c++) {
    const vector<int>& atoms = clst.clusters[c].atoms;
    bool append = true;
    for (uint r = 0; r < reduced_clusters.size(); r++) {
      append = false;
      for (int i = 0; i < order - 1; i++) {
        if (atoms[i] != reduced_clusters[r][i]) {
          append = true;
          i = order;
        }
      }
      if (!append) {  // If append stays false, the reduced cluster is not new
        r = reduced_clusters.size();
      }
    }
    if (append) {
      for (int i = 0; i < order - 1; i++) cluster[i] = atoms[i];
      reduced_clusters.push_back(cluster);
    }
  }
  return reduced_clusters;
}

//getTensorTransformations////////////////////////////////////////////////////
// This algorithm does two things: it generates the tensor transformations
// for each cluster to transform the IFCs of the inequivalent clusters; and
// it generates a list of equivalent IFCs for each inequivalent cluster to
// calculate the corrections.
//
// Check the typedefs at the beginning of the file for tform and v5int
void AnharmonicIFCs::getTensorTransformations(v4int& eq_ifcs,
                                              vector<vector<tform> >& transformations) {
  for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
    vector<tform> transform(clst.ineq_clusters[ineq].size() - 1);
    vector<vector<vector<int> > > eq(clst.nifcs, vector<vector<int> >(clst.ineq_clusters[ineq].size()));
    int ind = 0;
    for (int crt = 0; crt < clst.nifcs; crt++) {
      eq[ind][0].push_back(crt);
      ind++;
    }
    for (uint c = 1; c < clst.ineq_clusters[ineq].size(); c++) {
      tform trf;
      _cluster cluster_trans = clst.getCluster(clst.ineq_clusters[ineq][c]);
      int fg = 0, perm = 0, rw = 0, cl = 0, p = 0;
      vector<int> atoms_trans = cluster_trans.atoms;
      atoms_trans[0] = clst.sc2pcMap[atoms_trans[0]];  // transfer to pcell
      fg = cluster_trans.fgroup;
      perm = cluster_trans.permutation;
      for (int itrans = 0; itrans < clst.nifcs; itrans++) {
        std::pair<vector<int>, vector<double> > t;
        int ind_orig = 0;
        for (int iorig = 0; iorig < clst.nifcs; iorig++) {
          double coeff = 1.0;
          for (int o = 0; o < order; o++) {
            rw = cart_indices[itrans][o] + 1;
            p = clst.permutations[perm][o];
            cl = cart_indices[iorig][p] + 1;
            coeff *= clst.pcell.fgroup[fg].Uc[rw][cl];
            if (abs(coeff) < _ZERO_TOL_) {
              coeff = 0.0;
              o = order;
            }
          }
          if (abs(coeff) > _ZERO_TOL_) {
            t.first.push_back(iorig);
            t.second.push_back(coeff);
            eq[ind_orig][c].push_back(itrans);
          }
          ind_orig++;
        }
        trf.push_back(t);
      }
      transform[c-1] = trf;
    }
    eq_ifcs[ineq] = eq;
    transformations[ineq] =  transform;
  }
}

// END Initializers

// BEGIN Iterations
//applyLinCombs///////////////////////////////////////////////////////////////
// Sets the lineraly dependent IFCs according to the obtained linear
// combinations.
void AnharmonicIFCs::applyLinCombs(vector<vector<double> >& ifcs) {
  int c = 0, cart_ind = 0, cart_ind_indep = 0;
  for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
    c = clst.ineq_clusters[ineq][0];
    _linearCombinations lcomb = clst.linear_combinations[ineq];
    for (uint d = 0; d < lcomb.dependent.size(); d++) {
      cart_ind = lcomb.dependent[d];
      ifcs[c][cart_ind] = 0.0;  // reset
      for (uint lc = 0; lc < lcomb.indices[d].size(); lc++) {
        cart_ind_indep = lcomb.indices[d][lc];
        ifcs[c][cart_ind] += lcomb.coefficients[d][lc] * ifcs[c][cart_ind_indep];
      }
    }
  }
}

//transformIFCs///////////////////////////////////////////////////////////////
// Transforms all IFCs using symmetry. The first two loops go over all
// equivalent clusters and transform the inequivalent clusters using tensor
// transformations.
//
// See the top of this file for the typedef of tform.
void AnharmonicIFCs::transformIFCs(const vector<vector<tform> >& transformations,
                                   vector<vector<double> >& ifcs) {
  int clst_orig = 0, clst_trans = 0, cart_indices_orig = 0;
  double coeff = 0.0;
  for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
    clst_orig = clst.ineq_clusters[ineq][0];
    for (uint c = 1; c < clst.ineq_clusters[ineq].size(); c++) {
      clst_trans = clst.ineq_clusters[ineq][c];
      for (int itrans = 0; itrans < clst.nifcs; itrans++) {
        ifcs[clst_trans][itrans] = 0.0;  // reset
        const std::pair<vector<int>, vector<double> >& transf = transformations[ineq][c-1][itrans];
        for (uint t = 0; t < transf.first.size(); t++) {
          cart_indices_orig = transf.first[t];
          coeff = transf.second[t];
          ifcs[clst_trans][itrans] += coeff * ifcs[clst_orig][cart_indices_orig];
        }
      }
    }
  }
}

//calcSums////////////////////////////////////////////////////////////////////
// Determines the deviations from zero and the sum of the absolute values
// for each reduced cluster. Both are used for the correction of the IFCs
// while the deviations from zero are also used as errors for the sum rules.
void AnharmonicIFCs::calcSums(const vector<vector<int> >& reduced_clusters,
                              const vector<vector<double> >& ifcs,
                              vector<vector<double> >& dev_from_zero,
                              vector<vector<double> >& abssum) {
  dev_from_zero.assign(reduced_clusters.size(), vector<double>(clst.nifcs, 0));
  abssum.assign(reduced_clusters.size(), vector<double>(clst.nifcs, 0));
  for (uint i = 0; i < clst.ineq_clusters.size(); i++) {
    for (uint j = 0; j < clst.ineq_clusters[i].size(); j++) {
      uint r = findReducedCluster(reduced_clusters, clst.getCluster(clst.ineq_clusters[i][j]).atoms);
      for (int c = 0; c < clst.nifcs; c++) {
        dev_from_zero[r][c] += ifcs[clst.ineq_clusters[i][j]][c];
        abssum[r][c] += abs(ifcs[clst.ineq_clusters[i][j]][c]);
      }
    }
  }
}

//correctIFCs/////////////////////////////////////////////////////////////////
// Corrects the IFCs the comply with the sum rule using weighted averages.
//
// Check the typedef at the beginning of the file for v5int.
void AnharmonicIFCs::correctIFCs(vector<vector<double> >& ifcs,
                                 const vector<vector<double> >& dev_from_zero, 
                                 const vector<vector<double> >& abssum,
                                 const vector<vector<int> >& reduced_clusters,
                                 const v4int& eq_ifcs) {
  vector<int> eq;
  for (uint ineq = 0; ineq < clst.ineq_clusters.size(); ineq++) {
    uint nclusters = clst.ineq_clusters[ineq].size();
    int ic = clst.ineq_clusters[ineq][0];
    // Calculate correction terms
    vector<vector<double> > correction_terms(nclusters);
    for (uint c = 0; c < nclusters; c++) {
      correction_terms[c] = getCorrectionTerms(clst.ineq_clusters[ineq][c],
                                               reduced_clusters, ifcs,
                                               dev_from_zero, abssum);
    }
    
    // Correct the linearly independent IFCs
    const _linearCombinations& lcomb = clst.linear_combinations[ineq];
    const vector<int>& indep = lcomb.independent;
    for (uint i = 0; i < indep.size(); i++) {
      int neq = 0, cl = 0;
      double corrected_ifc = 0.0;
      for (uint c = 0; c < nclusters; c++) {
        cl = clst.ineq_clusters[ineq][c];
        eq = eq_ifcs[ineq][indep[i]][c];
        uint eqsize = eq.size();
        for (uint e = 0; e < eqsize; e++) {
          if (ifcs[cl][eq[e]] != 0.0) {
            neq++;
            corrected_ifc += correction_terms[c][eq[e]] * ifcs[ic][indep[i]] / ifcs[cl][eq[e]];
          }
        }
        for (uint dep = 0; dep < lcomb.indep2depMap[i].size(); dep++) {
          eq = eq_ifcs[ineq][lcomb.indep2depMap[i][dep]][c];
          for (uint e = 0; e < eqsize; e++) {
            if (ifcs[cl][eq[e]] != 0.0) {
              neq++;
              corrected_ifc += correction_terms[c][eq[e]] * ifcs[ic][indep[i]] / ifcs[cl][eq[e]];
            }
          }
        }
      }
      if (neq > 0) {
        ifcs[ic][indep[i]] = mixing_coefficient * ifcs[ic][indep[i]] + (1 - mixing_coefficient) * corrected_ifc/((double) neq);
      }
    }
  }
}

//getCorrectionTerms//////////////////////////////////////////////////////////
// Calculates the correction term for each IFC.
vector<double>
    AnharmonicIFCs::getCorrectionTerms(const int& c,
                                       const vector<vector<int> >& reduced_clusters,
                                       const vector<vector<double> >& ifcs,
                                       const vector<vector<double> >& dev_from_zero,
                                       const vector<vector<double> >& abssum) {
  vector<double> correction_terms = ifcs[c];
  vector<double> correction(clst.nifcs);
  for (int i = 0; i < clst.nifcs; i++) correction[i] = std::abs(correction_terms[i]);
  uint r = findReducedCluster(reduced_clusters, clst.getCluster(c).atoms);
  for (int crt = 0; crt < clst.nifcs; crt++) {
    if (abssum[r][crt] != 0.0) correction[crt] *= dev_from_zero[r][crt]/abssum[r][crt];
    else correction[crt] = 0.0;
    correction_terms[crt] -= correction[crt];
  }
  return correction_terms;
}

uint AnharmonicIFCs::findReducedCluster(const vector<vector<int> >& reduced_clusters,
                                        const vector<int>& atoms) {
  for (uint r = 0; r < reduced_clusters.size(); r++) {
    int at = 0;
    for (at = 0; at < order - 1; at++) {
      if (reduced_clusters[r][at] != atoms[at]) break;
    }
    if (at == order - 1) return r;
  }
  return -1;
}

// END Iterations

}  // namespace apl

/********************************* FILE I/O *********************************/

namespace apl {

// BEGIN Write files
//writeIFCsToFile/////////////////////////////////////////////////////////////
// Writes the AnharmonicIFCs object and minimal structure information to an
// XML file.
void AnharmonicIFCs::writeIFCsToFile(const string& filename) {
  stringstream output;
  // Header
  output << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << std::endl;
  output << "<anharmonicifcs>" << std::endl;

  output << writeParameters();
  output << writeIFCs();
  output << "</anharmonicifcs>" << std::endl;
  aurostd::stringstream2file(output, filename);
  if (!aurostd::FileExist(filename)) {
    string function = _AAPL_IFCS_ERR_PREFIX_ + "writeIFCsToFile";
    string message = "Could not write tensor to file.";
    throw xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
  }
}

//writeParameters/////////////////////////////////////////////////////////////
// Writes the calculation parameters and minimal structure information to the
// XML file.
string AnharmonicIFCs::writeParameters() {
  stringstream parameters;
  string tab = " ";

  // Info about calculation run
  parameters << tab << "<generator>" << std::endl;
  string time = aflow_get_time_string();
  if (time[time.size() - 1] == '\n') time.erase(time.size() - 1);
  parameters << tab << tab << "<i name=\"date\" type=\"string\">" << time << "</i>" << std::endl;
  parameters << tab << tab << "<i name=\"checksum\" file=\"" << _AFLOWIN_;
  parameters << "\" type=\"" << APL_CHECKSUM_ALGO << "\">" << std::hex << aurostd::getFileCheckSum(aflags.Directory + "/" + _AFLOWIN_ + "", APL_CHECKSUM_ALGO);  // ME190219
  parameters.unsetf(std::ios::hex);  // ME190125 - Remove hexadecimal formatting
  parameters  << "</i>" << std::endl;
  parameters << tab << "</generator>" << std::endl;

  // Distortion magnitude
  parameters << tab << "<distortion_magnitude units=\"Angstrom\" cs=\"cartesian\">";
  parameters << distortion_magnitude << "</distortion_magnitude>" << std::endl;

  //Order
  parameters << tab << "<order>" << order << "</order>" << std::endl;

  // Iteration parameters
  parameters << tab << "<iteration>" << std::endl;
  // std::dec prevents hexadecimal output
  parameters << tab << tab << "<max_iter>" << std::dec << max_iter << "</max_iter>" << std::endl;
  parameters << tab << tab << "<mixing_coefficient>";
  parameters << std::setprecision(8) << mixing_coefficient;
  parameters << "</mixing_coefficient>" << std::endl;
  parameters << tab << tab << "<sumrule_threshold>";
  parameters << std::setprecision(15) << sumrule_threshold;
  parameters << "</sumrule_threshold>" << std::endl;
  parameters << tab << "</iteration>" << std::endl;

  // Structure
  parameters << tab << "<structure units=\"Angstrom\" cs=\"fractional\">" << std::endl;
  parameters << tab << tab << "<varray name=\"pcell lattice\">" << std::endl;
  for (int i = 1; i < 4; i++) {
    parameters << tab << tab << tab << "<v>";
    for (int j = 1; j < 4; j++) {
      parameters << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      parameters << std::setprecision(8);
      parameters << std::setw(15) << clst.pcell.lattice[i][j];
    }
    parameters << "</v>" << std::endl;
  }
  parameters << tab << tab << "</varray>" << std::endl;
  parameters << tab << tab << "<varray name=\"positions\">" << std::endl;
  for (uint i = 0; i < clst.pcell.atoms.size(); i++) {
    int t = clst.pcell.atoms[i].type;
    parameters << tab << tab << tab << "<v species=\"" << clst.pcell.species[t] << "\">";
    for (int j = 1; j < 4; j++) {
      parameters << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
      parameters << std::setprecision(8);
      parameters << std::setw(15) << clst.pcell.atoms[i].fpos[j];
    }
    parameters << "</v>" << std::endl;
  }
  parameters << tab << tab << "</varray>" << std::endl;
  parameters << tab << tab << "<varray name=\"supercell\">" << std::endl;
  parameters << tab << tab << tab << "<v>";
  for (int i = 1; i < 4; i++) {
    parameters << tab << clst.sc_dim[i];
  }
  parameters << "</v>" << std::endl;
  parameters << tab << tab << "</varray>" << std::endl;
  parameters << tab << "</structure>" << std::endl;
  return parameters.str();
}

//writeIFCs///////////////////////////////////////////////////////////////////
// Writes the force constants part of the XML file.
string AnharmonicIFCs::writeIFCs() {
  stringstream ifcs;
  string tab = " ";
  int precision = 15;
  int extra = 5; // first digit, decimal point, minus sign + 2 spaces
  double max_ifc = -1E30;
  for (uint i = 0; i < force_constants.size(); i++) {
    for (uint j = 0; j < force_constants[i].size(); j++) {
      if (std::abs(force_constants[i][j]) > max_ifc) max_ifc = std::abs(force_constants[i][j]);
    }
  }
  int width = (int) log10(max_ifc);
  if (width < 0) {
    width = precision + extra;
  } else {
    width += precision + extra;
  }

  ifcs << tab << "<force_constants>" << std::endl;
  for (uint c = 0; c < clst.clusters.size(); c++) {
    ifcs << tab << tab << "<varray atoms=\""
         << aurostd::joinWDelimiter(clst.clusters[c].atoms, " ") << "\">" << std::endl;
    int crt = 0;
    for (int i = 0; i < clst.nifcs/9; i++) {
      ifcs << tab << tab << tab << "<varray slice=\"" << cart_indices[crt][0];
      for (int j = 1; j < order - 2; j++) {
        ifcs << " " << cart_indices[crt][j];
      }
      ifcs << "\">" << std::endl;
      for (int j = 0; j < 3; j++) {
        ifcs << tab << tab << tab << tab << "<v>";
        for (int k = 0; k < 3; k++) {
          ifcs << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
          ifcs << std::setprecision(precision);
          ifcs << std::setw(width) << force_constants[c][crt];
          crt++;
        }
        ifcs << "</v>" << std::endl;
      }
      ifcs << tab << tab << tab << "</varray>" << std::endl;
    }
    ifcs << tab << tab << "</varray>" << std::endl;
  }
  ifcs << tab << "</force_constants>" << std::endl;
  return ifcs.str();
}

// END Write files

// BEGIN Read files
//readIFCsFromFile////////////////////////////////////////////////////////////
// Reads an AnharmonicIFCs object from an XML file.
void AnharmonicIFCs::readIFCsFromFile(const string& filename) {
  // Open file and handle exceptions
  string function = _AAPL_IFCS_ERR_PREFIX_ + "readIFCsFromFile";
  stringstream message;

  if (!aurostd::EFileExist(filename) && !aurostd::FileExist(filename)) {
    message << "Could not open file " << filename << ". File not found.";
    throw xerror(_AFLOW_FILE_NAME_,function, message, _FILE_NOT_FOUND_);
  }
  vector<string> vlines;
  aurostd::efile2vectorstring(filename, vlines);
  if (vlines.size() == 0) {
    message << "Cannot open file " << filename << ". File empty or corrupt.";
    throw xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
  }

  // Start reading
  uint line_count = 0;
  string line = vlines[line_count++];

  // Check that this is a valid xml file
  if (line.find("xml") == string::npos) {
    message << "File is not a valid xml file.";
    throw xerror(_AFLOW_FILE_NAME_,function, message, _FILE_WRONG_FORMAT_);
  }

  // Check if xml file can be used to read anharmonic IFCs
  if (checkCompatibility(line_count, vlines)) {
    force_constants = readIFCs(line_count, vlines);
  } else {
    message << "The settings in the hibernate file and the aflow.in file are incompatible.";
    throw xerror(_AFLOW_FILE_NAME_,function, message, _RUNTIME_ERROR_);
  }
}

//checkCompatibility//////////////////////////////////////////////////////////
// Checks if the hibernate XML file is compatible with the aflow.in file. If
// the checksum in the XML file is the same as the checksum of the aflow.in
// file, then the parameters are the same. If not, the function checks if the
// parameters relevant for the IFC calculation (supercell, order, calculation
// parameters) are the same. This prevents the anharmonic IFCs from being
// recalculated when only post-processing parameters are changed.
bool AnharmonicIFCs::checkCompatibility(uint& line_count, 
                                        const vector<string>& vlines) {
  string function = _AAPL_IFCS_ERR_PREFIX_ + "checkCompatibility";
  string line = "";
  stringstream message;
  bool compatible = true;
  int t = 0;
  vector<string> tokens;
  uint vsize = vlines.size();

  // Compare checksum
  while (true) {
    if (line_count == vsize) {
      message << "Checksum not found in hibernate file.";
      throw xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
    }
    line = vlines[line_count++];
    if (line.find("checksum") != string::npos) {
      break;
    }
  }

  t = line.find_first_of(">") + 1;
  tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
  if (strtoul(tokens[0].c_str(), NULL, 16) != aurostd::getFileCheckSum(aflags.Directory + "/" + _AFLOWIN_, APL_CHECKSUM_ALGO)) {  // ME190219
    message << "The " << _AFLOWIN_ << " file has been changed from the hibernated state. ";

    tokens.clear();
    // Compare calculation parameters
    //// distortion_magnitude
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find distortion_magnitude tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("distortion_magnitude") != string::npos) {
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        double dist_mag = aurostd::string2utype<double>(tokens[0]);
        tokens.clear();
        if (dist_mag != distortion_magnitude) {
          message << "Hibernate file and aflow.in have different distortion magnitudes. ";
          compatible = false;
        }
        break;
      }
    }

    //// order
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find order tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("order") != string::npos) {
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        int ord = aurostd::string2utype<int>(tokens[0]);
        tokens.clear();
        if (ord != order) {
          message << "Hibernate file and aflow.in have different IFC order. ";
          compatible = false;
        }
        break;
      }
    }

    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find iteration tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("iteration") != string::npos) {
        break;
      }
    }
    //// max_iter
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find max_iter tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("max_iter") != string::npos) {
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        int iter = aurostd::string2utype<int>(tokens[0]);
        tokens.clear();
        if (iter != max_iter) {
          message << "Hibernate file and aflow.in have different max. number of iterations. ";
          compatible = false;
        }
        break;
      }
    }

    //// mixing_coefficient
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find mixing_coefficient tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("mixing_coefficient") != string::npos) {
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        double mix = aurostd::string2utype<bool>(tokens[0]);
        tokens.clear();
        if (mix != mixing_coefficient) {
          message << "Hibernate file and aflow.in have different mixing coefficients. ";
          compatible = false;
        }
        break;
      }
    }

    //// sumrule_threshold
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find sumrule_threshold tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("sumrule_threshold") != string::npos) {
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        double thresh = aurostd::string2utype<double>(tokens[0]);
        tokens.clear();
        if (thresh != sumrule_threshold) {
          message << "Hibernate file and aflow.in have different convergence criteria. ";
          compatible = false;
        }
        break;
      }
    }

    // Compare supercells
    //// Lattice
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find structure tag. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("structure") != string::npos) {
        break;
      }
    }

    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find primitive lattice vectors. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("varray name=\"pcell lattice\"") != string::npos) {
        xmatrix<double> latt(3, 3);
        for (int i = 1; i < 4; i++) {
          line = vlines[line_count++];
          if (line_count == vsize) {
            message << "pcell lattice tag is corrupt. ";
          }
          t = line.find_first_of(">") + 1;
          tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
          for (int j = 1; j < 3; j++) {
            latt(i, j) = aurostd::string2utype<double>(tokens[j - 1]);
          }
          tokens.clear();
        }
        line = vlines[line_count++];
        if (line.find("</varray>") == string::npos) {
          message << "pcell lattice tag is corrupt. ";
          compatible = false;
        } else {
          break;
        }
        if (latt != clst.pcell.lattice) {
          message << "Hibernate file and aflow.in do not have the same lattice. ";
          compatible = false;
        }
      }
    }

    //// Atomic positions
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find atomic positions. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("varray name=\"positions\"") != string::npos) {
        // First check for tag corruption and extract everything
        vector<string> species;
        vector<xvector<double> > positions;
        while (compatible) {
          if (line_count == vsize) {
            message << "positions tag is corrupt. ";
            compatible = false;
          }
          line = vlines[line_count++];
          if (line.find("species=\"") != string::npos) {
            // Extract species
            t = line.find_first_of("\"") + 1;
            tokenize(line.substr(t, line.find_last_of("\"") - t), tokens, string(" "));
            species.push_back(tokens[0]);
            tokens.clear();
            // Extract positions
            t = line.find_first_of(">") + 1;
            tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
            xvector<double> fpos(3);
            for (int i = 1; i < 4; i++) {
              fpos(i) = aurostd::string2utype<double>(tokens[i-1]);
            }
            positions.push_back(fpos);
            tokens.clear();
          } else if (line.find("/varray") != string::npos) {
            break;
          } else {
            message << "positions tag is corrupt. ";
            compatible = false;
          }
        }
        // Now compare with the primitive cell from the aflow.in file
        uint nspecies = species.size();
        if (nspecies == clst.pcell.atoms.size()) {
          for (uint sp = 0; sp < nspecies; sp++) {
            int type = clst.pcell.atoms[sp].type;
            string spec = clst.pcell.species[type];
            xvector<double> pos = clst.pcell.atoms[sp].fpos;
            if (species[sp] != spec) {
              message << "The structures in the hibernate file and aflow.in ";
              message << "have different species. ";
              compatible = false;
              sp = nspecies;
            } else if (positions[sp] != pos) {
              message << "The structures in the hibernate file and aflow.in ";
              message << "have different atomic positions.";
              compatible = false;
              sp = nspecies;
            }
          }
        } else {
          message << "The structures in the hibernate file and in aflow.in ";
          message << "do not have the same number of atoms. ";
          compatible = false;
        }
        break;
      }
    }

    //// Supercell dimensions
    while (compatible) {
      if (line_count == vsize) {
        message << "Could not find supercell dimensions. ";
        compatible = false;
      }
      line = vlines[line_count++];
      if (line.find("varray name=\"supercell\"") != string::npos) {
        line = vlines[line_count++];
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        xvector<int> sc(3);
        for (int i = 1; i < 4; i++) {
          sc(i) = aurostd::string2utype<double>(tokens[i-1]);
        }
        if (sc == clst.sc_dim) {
          break;
        } else {
          message << "The supercells in the hibernate file and in aflow.in ";
          message << "have different dimensions. ";
          compatible = false;
        }
      }
    }

    if (compatible) {
      message << "The relevant settings appear to be the same, ";
      message << "so the anharmonic IFCs will not be recalculated. ";
      message << "Make sure that the changes do not impact the IFCs.";
    } else {
      message << "Anharmonic IFCs need to be recalculated.";
    }
    _logger << apl::warning << "apl::AnharmonicIFCs::readIFCsfromFile(); ";
    _logger << message.str() << apl::endl;
  }
  return compatible;
}

//readIFCs////////////////////////////////////////////////////////////////////
// Reads the IFCs from the hibernate file.
vector<vector<double> > AnharmonicIFCs::readIFCs(uint& line_count,
                                                 const vector<string>& vlines) {
  string function = _AAPL_IFCS_ERR_PREFIX_ + "readIFCs";
  string line = "", message = "";
  vector<int> atoms(order), slice(order - 2);
  vector<string> tokens;
  int t = 0;
  uint vsize = vlines.size();
  vector<vector<double> > ifcs;

  // Find force_constants tag
  while (true) {
    if (line_count == vsize) {
      message = "force_constants tag not found.";
      throw xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
    }
    line = vlines[line_count++];
    if (line.find("force_constants") != string::npos) {
      break;
    }
  }

  // Read IFCs
  vector<double> ifc;
  while (line.find("/force_constants") == string::npos) {
    if (line_count == vsize) {
      message = "force_constants tag incomplete.";
      throw xerror(_AFLOW_FILE_NAME_,function, message, _FILE_CORRUPT_);
    }
    line = vlines[line_count++];
    if (line.find("atoms") != string::npos) {
      ifc.clear();
    } else if (line.find("slice") != string::npos) {
      // If a slice is found, populate tensor
      for (int i = 0; i < 3; i++) {
        line = vlines[line_count++];
        t = line.find_first_of(">") + 1;
        tokenize(line.substr(t, line.find_last_of("<") - t), tokens, string(" "));
        for (int j = 0; j < 3; j++) {
          double val = aurostd::string2utype<double>(tokens[j]);
          ifc.push_back(val);
        }
        tokens.clear();
      }
      line_count++;
    } else if (line.find("</varray>") != string::npos) {
      ifcs.push_back(ifc);
    }
  }

  return ifcs;
}

// END Read files

} // namespace apl

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                Aflow Marco Esters - Duke University 2018                *
// *                                                                         *
// ***************************************************************************
