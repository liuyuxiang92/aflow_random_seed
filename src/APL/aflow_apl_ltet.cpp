#include "aflow_apl.h"

#define _DEBUG_APL_LTET_ false

using aurostd::xvector;
using std::vector;
using std::string;

static const string _APL_LTET_ERR_PREFIX_ = "apl::LTMethod::";

namespace apl {

LTMethod::LTMethod(QMesh& qm, Logger& l) : _qm(qm), _logger(l) {
  if (_qm.getnQPs() < 4) {
    string function = _APL_LTET_ERR_PREFIX_ + "LTMethod()";
    string message = "At least four q-points are required for the linear tetrahedron method.";
    throw aurostd::xerror(function, message, _RUNTIME_ERROR_);
  }
  free();
  _volumePerTetrahedron = 1.0/(double)(6 * _qm.getGrid(1) * _qm.getGrid(2) * _qm.getGrid(3));
  generateTetrahedra();
  _logger << "Found " << _nTetra << " irreducible tetrahedra." << apl::endl;
}

LTMethod::LTMethod(const LTMethod& that) : _qm(that._qm), _logger(that._logger) {
  *this = that;
}

LTMethod& LTMethod::operator=(const LTMethod& that) {
  if (this != &that) {
    _irredTetrahedra = that._irredTetrahedra;
    _logger = that._logger;
    _nTetra = that._nTetra;
    _qm = that._qm;
    _volumePerTetrahedron = that._volumePerTetrahedron;
    _weights = that._weights;
  }
  return *this;
}

LTMethod::~LTMethod() {
  free();
}

void LTMethod::free() {
  _irredTetrahedra.clear();
  _nTetra = 0;
  _volumePerTetrahedron = 0.0;
  _weights.clear();
}

void LTMethod::generateTetrahedra() {
  vector<vector<xvector<int> > > tetra_init = initializeTetrahedra();
  findMostCompact(tetra_init);
  getIrreducibleTetrahedra(tetra_init);
}

vector<vector<xvector<int> > > LTMethod::initializeTetrahedra() {
  vector<vector<xvector<int> > > tetra_init(6, vector<xvector<int> >(4, xvector<int>(3)));
  //Tetrahedron 3126
  tetra_init[0][0][1] = 0; tetra_init[0][0][2] = 0; tetra_init[0][0][3] = 0;
  tetra_init[0][1][1] = 0; tetra_init[0][1][2] = 1; tetra_init[0][1][3] = 0;
  tetra_init[0][2][1] = 1; tetra_init[0][2][2] = 1; tetra_init[0][2][3] = 0;
  tetra_init[0][3][1] = 1; tetra_init[0][3][2] = 1; tetra_init[0][3][3] = 1;
  //Tetrahedron 3426
  tetra_init[1][0][1] = 0; tetra_init[1][0][2] = 0; tetra_init[1][0][3] = 0;
  tetra_init[1][1][1] = 1; tetra_init[1][1][2] = 0; tetra_init[1][1][3] = 0;
  tetra_init[1][2][1] = 1; tetra_init[1][2][2] = 1; tetra_init[1][2][3] = 0;
  tetra_init[1][3][1] = 1; tetra_init[1][3][2] = 1; tetra_init[1][3][3] = 1;
  //Tetrahedron 3486
  tetra_init[2][0][1] = 0; tetra_init[2][0][2] = 0; tetra_init[2][0][3] = 0;
  tetra_init[2][1][1] = 1; tetra_init[2][1][2] = 0; tetra_init[2][1][3] = 0;
  tetra_init[2][2][1] = 1; tetra_init[2][2][2] = 0; tetra_init[2][2][3] = 1;
  tetra_init[2][3][1] = 1; tetra_init[2][3][2] = 1; tetra_init[2][3][3] = 1;
  //Tetrahedron 3156
  tetra_init[3][0][1] = 0; tetra_init[3][0][2] = 0; tetra_init[3][0][3] = 0;
  tetra_init[3][1][1] = 0; tetra_init[3][1][2] = 1; tetra_init[3][1][3] = 0;
  tetra_init[3][2][1] = 0; tetra_init[3][2][2] = 1; tetra_init[3][2][3] = 1;
  tetra_init[3][3][1] = 1; tetra_init[3][3][2] = 1; tetra_init[3][3][3] = 1;
  //Tetrahedron 3756
  tetra_init[4][0][1] = 0; tetra_init[4][0][2] = 0; tetra_init[4][0][3] = 0;
  tetra_init[4][1][1] = 0; tetra_init[4][1][2] = 0; tetra_init[4][1][3] = 1;
  tetra_init[4][2][1] = 0; tetra_init[4][2][2] = 1; tetra_init[4][2][3] = 1;
  tetra_init[4][3][1] = 1; tetra_init[4][3][2] = 1; tetra_init[4][3][3] = 1;
  //Tetrahedron 3786
  tetra_init[5][0][1] = 0; tetra_init[5][0][2] = 0; tetra_init[5][0][3] = 0;
  tetra_init[5][1][1] = 0; tetra_init[5][1][2] = 0; tetra_init[5][1][3] = 1;
  tetra_init[5][2][1] = 1; tetra_init[5][2][2] = 0; tetra_init[5][2][3] = 1;
  tetra_init[5][3][1] = 1; tetra_init[5][3][2] = 1; tetra_init[5][3][3] = 1;

  return tetra_init;
}

void LTMethod::findMostCompact(vector<vector<xvector<int> > >& tetrahedra) {
  // Determine the configuration that yields the most compact tetrahedra
  int lxx = 0;
  int lyy = 0;
  double gmax = 1e30;
  vector<xvector<double> > tet(4, xvector<double>(3));
  for (int lx = 0; lx <= 1; lx++) {
    for (int ly = 0; ly <= 1; ly++) {
      double d, lmax = 0.0;
      for (int i = 0; i < 6; i++) {
        // Transform tetrahedra 
        for (int j = 0; j < 4; j++) {
          for (int k = 1; k < 4; k++) tet[j][k] = (double) tetrahedra[i][j][k];
          if (lx == 1) tet[j][1] = 1 - tet[j][1];
          if (ly == 1) tet[j][2] = 1 - tet[j][2];
          tet[j] = F2C(_qm.getReciprocalCell().lattice, tet[j]);
        }
        // Measure tetrahedra sides and determine minimum
        for (int j = 0; j < 3; j++) {
          for (int k = j + 1; k < 4; k++) {
            d = aurostd::modulus(tet[j] - tet[k]);
            if (d > lmax) lmax = d;
          }
        }
      }
      if (lmax < gmax) {
        lxx = lx;
        lyy = ly;
        gmax = lmax;
      }
    }
  }

  // Apply the most compact configuration to existing tetrahedra
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 4; j++) {
      if (lxx == 1) tetrahedra[i][j][1] = 1 - tetrahedra[i][j][1];
      if (lyy == 1) tetrahedra[i][j][2] = 1 - tetrahedra[i][j][2];
    }
  }
}

void LTMethod::getIrreducibleTetrahedra(const vector<vector<xvector<int> > >& tetrahedra) {
  vector<vector<vector<int> > > cornerMap(2, vector<vector<int> >(2, vector<int>(2)));
  vector<int> tet(4);
  int j1, j2, j3, m, t;
  for (int q3 = 0; q3 < _qm.getGrid(3); q3++) {
    for (int q2 = 0; q2 < _qm.getGrid(2); q2++) {
      for (int q1 = 0; q1 < _qm.getGrid(1); q1++) {
        // Get the corners of the microcell and map them to the irreducible
        // q-points
        for (int k3 = 0; k3 <= 1; k3++) {
          j3 = (q3 + k3) % _qm.getGrid(3);
          for (int k2 = 0; k2 <= 1; k2++) {
            j2 = (q2 + k2) % _qm.getGrid(2);
            for (int k1 = 0; k1 <= 1; k1++) {
              j1 = (q1 + k1) % _qm.getGrid(1);
              cornerMap[k1][k2][k3] = _qm.getIrredQPointIndex(j1, j2, j3);
            }
          }
        }
        // Get the tetrahedra associated with the corners
        for (int i = 0; i < 6; i++) {
          for (int j = 0; j < 4; j++) {
            tet[j] = cornerMap[tetrahedra[i][j][1]][tetrahedra[i][j][2]][tetrahedra[i][j][3]];
          }
          std::sort(tet.begin(), tet.end());
          for (t = 0; t < _nTetra; t++) {
            for (m = 0; m < 4; m++) {
              if (_irredTetrahedra[t][m] != tet[m]) break;
            }
            if (m == 4) break;
          }
          if (t == _nTetra) {
            _irredTetrahedra.push_back(tet);
            _weights.push_back(1);
            _nTetra++;
          } else {
            _weights[t]++;
          }
        }
      }
    }
  }
  _nTetra = (int) _irredTetrahedra.size();
}

const vector<vector<int> >& LTMethod::getTetrahedra() const {
  return _irredTetrahedra;
}

const vector<int>& LTMethod::getTetrahedron(const int& tetrahedron) const {
  return _irredTetrahedra[tetrahedron];
}

const int& LTMethod::getCorner(const int& tetrahedron, const int& corner) const {
  return _irredTetrahedra[tetrahedron][corner];
}

const int& LTMethod::getCornerIrred(const int& tetrahedron, const int& corner) const {
  return _qm.getIrredQPointIndex(_irredTetrahedra[tetrahedron][corner]);
}

const int& LTMethod::getnTetrahedra() const {
  return _nTetra;
}

const double& LTMethod::getVolumePerTetrahedron() const {
  return _volumePerTetrahedron;
}

const vector<int>& LTMethod::getWeights() const {
  return _weights;
}

const int& LTMethod::getWeight(const int& i) const {
  return _weights[i];
}

} // namespace apl
