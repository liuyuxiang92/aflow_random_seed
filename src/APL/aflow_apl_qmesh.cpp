//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
//****************************************************************************

// This class describes a mesh of q-points.

#include "aflow_apl.h"

#define _DEBUG_APL_QMESH_ false

using aurostd::xvector;
using std::vector;
using std::string;

static const string _APL_QMESH_ERR_PREFIX_ = "apl::QMesh::";
static const string _APL_QMESH_MODULE_ = "QMESH";  // for the logger

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         CONSTRUCTORS/DESTRUCTORS                         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  // Default Constructor
  QMesh::QMesh(ofstream& mf) {
    free();
    messageFile = &mf;
    _directory = "./";
  }

  QMesh::QMesh(const xvector<int>& grid, const xstructure& xs, ofstream& mf,
      bool include_inversions, bool gamma_centered) {
    free();
    messageFile = &mf;
    initialize(grid, xs, include_inversions, gamma_centered);
  }

  QMesh::QMesh(const vector<int>& vgrid, const xstructure& xs, ofstream& mf,
      bool include_inversions, bool gamma_centered) {
    free();
    messageFile = &mf;
    initialize(aurostd::vector2xvector(vgrid), xs, include_inversions, gamma_centered);
  }

  // Copy constructors
  QMesh::QMesh(const QMesh& that) {
    copy(that);
  }

  QMesh& QMesh::operator=(const QMesh& that) {
    if (this != &that) copy(that);
    return *this;
  }

  void QMesh::copy(const QMesh& that) {
    _ibzqpts = that._ibzqpts;
    _isGammaCentered = that._isGammaCentered;
    _littleGroups = that._littleGroups;
    _littleGroupsCalculated = that._littleGroupsCalculated;
    messageFile = that.messageFile;
    _nIQPs = that._nIQPs;
    _nQPs = that._nQPs;
    _qptGrid = that._qptGrid;
    _qptMap = that._qptMap;
    _qpoints = that._qpoints;
    _recCell = that._recCell;
    _reduced = that._reduced;
    _shifted = that._shifted;  // ME190813
    _shift = that._shift;
    _weights = that._weights;
  }

  // Destructor
  QMesh::~QMesh() {
    free();
  }

  void QMesh::free() {
    xvector<int> zeroint(3);
    xvector<double> zerodbl(3);
    xmatrix<double> zeroMatrix(3, 3);
    _directory = "";
    _ibzqpts.clear();
    _isGammaCentered = false;
    _littleGroups.clear();
    _littleGroupsCalculated = false;
    _reduced = false;
    _nIQPs = 0;
    _nQPs = 0;
    _qpoints.clear();
    _qptGrid = zeroint;
    _qptMap.clear();
    _recCell.lattice = zeroMatrix;
    _recCell.rlattice = zeroMatrix;
    _recCell.c2f = zeroMatrix;
    _recCell.f2c = zeroMatrix;
    _recCell.skewed = false;
    _recCell.pgroup.clear();
    _shifted = false;  // ME190701
    _shift = zerodbl;
    _weights.clear();
  }

  void QMesh::clear(ofstream& mf) {
    QMesh that(mf);
    copy(that);
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                              INTERFACE                                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// For the logger
namespace apl {

  void QMesh::setDirectory(const string& dir) {
    _directory = dir;
  }

  const string& QMesh::getDirectory() const {
    return _directory;
  }

  ofstream& QMesh::getOutputStream() {
    return *messageFile;
  }

}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                          Q-POINT FUNCTIONS                               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //initialize////////////////////////////////////////////////////////////////
  // Initializes the q-point grid
  void QMesh::initialize(const xvector<int>& grid, const xstructure& xs,
      bool include_inversions, bool gamma_centered) {
    setGrid(grid);
    setupReciprocalCell(xs, include_inversions);
    generateGridPoints(gamma_centered);
  }

  //setGrid///////////////////////////////////////////////////////////////////
  // Sets up the grid size
  void QMesh::setGrid(const xvector<int>& grid) {
    _qptGrid = grid;
    _nQPs = _qptGrid[1] * _qptGrid[2] * _qptGrid[3];
    _nIQPs = _nQPs;
    _qptMap.assign(_qptGrid[1], vector<vector<int> >(_qptGrid[2], vector<int>(_qptGrid[3])));
  }

  //setupReciprocalCell///////////////////////////////////////////////////////
  // Sets up the reciprocal cell that belongs to the q-mesh.
  void QMesh::setupReciprocalCell(xstructure xs, bool include_inversions) {
    _recCell.rlattice = xs.lattice;
    _recCell.lattice = ReciprocalLattice(_recCell.rlattice);
    _recCell.f2c = trasp(_recCell.lattice);
    _recCell.c2f = inverse(_recCell.f2c);

    // Determine skewedness
    xvector<double> min_distances(3);
    for (int i = 1; i < 4; i++) {
      min_distances[i] = aurostd::modulus(_recCell.lattice(i))/((double) _qptGrid[i]);
    }
    double min_dist = aurostd::min(min_distances);
    double tol = _AFLOW_APL_EPS_;
    _recCell.skewed = SYM::isLatticeSkewed(_recCell.lattice, min_dist, tol);

    // Calculate the point group of the reciprocal cell. Literature and phonon
    // codes are not consistent about whether to use pgroupk or pgroupk_xtal.
    // To get all symmetry properties (see DOI: 10.1103/RevModPhys.40.1),
    // pgroupk_xtal must be used or else the transformation properties of the
    // dynamical matrix cannot be captured since they require symmetry
    // operations that map atoms in real space. Thus, pgroupk_xtal needs to be
    // used to get the irreducible wedge - using pgroupk is not correct.
    // However, observables such as the phonon frequencies, eigenvectors, or
    // phonon-phonon scattering matrices also have inversion symmetry, which
    // is not always present in pgroupk_xtal. So, unless the dynamical matrix
    // itself is needed, the Patterson symmetry (pgroupk_Patterson) can be
    // used to create the irreducible wedge.
    if (include_inversions && !xs.pgroupk_Patterson_calculated) {
      xs.CalculateSymmetryPointGroupKPatterson(false);
    } else if (!xs.pgroupk_xtal_calculated) {
      xs.CalculateSymmetryPointGroupKCrystal(false);
    }

    if ((include_inversions && !xs.pgroupk_Patterson_calculated) ||
        (!include_inversions && !xs.pgroupk_xtal_calculated)) {
      string function = _APL_QMESH_ERR_PREFIX_ + "setupReciprocalCell()";
      string message = "Calculation of the point group of the reciprocal cell unsuccessful.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _RUNTIME_ERROR_);
    }

    if (include_inversions) _recCell.pgroup=xs.pgroupk_Patterson;
    else _recCell.pgroup = xs.pgroupk_xtal;
  }

  //generateGridPoints////////////////////////////////////////////////////////
  // Generates all the grid points. No reductions is performed yet since not
  // every purpose requires the irreducible q-points.
  void QMesh::generateGridPoints(bool force_gamma) {
    stringstream message;
    message << "Generating a " << _qptGrid[1] << "x" << _qptGrid[2] << "x" << _qptGrid[3] << " q-point mesh.";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_QMESH_MODULE_, message, _directory, *messageFile, std::cout);
    _qpoints.resize(_nQPs);
    _ibzqpts.resize(_nQPs);  // Before making the mesh irreducible, treat all q-points as irreducible q-points
    _weights.assign(_nQPs, 1);

    // Use Monkhorst-Pack formula to generate a mesh - do not center yet. Cartesian
    // coordinates will be calculated after all shifts have been performed.
    double q1 = 0.0, q2 = 0.0, q3 = 0.0;
    _qpoint qpt;
    int q = 0;
    for (int s = 1; s <= _qptGrid[3]; s++) {
      q3 = (2.0 * s - _qptGrid[3] - 1)/(2.0 * _qptGrid[3]);
      for (int r = 1; r <= _qptGrid[2]; r++) {
        q2 = (2.0 * r - _qptGrid[2] - 1)/(2.0 * _qptGrid[2]);
        for (int p = 1; p <= _qptGrid[1]; p++) {
          q1 = (2.0 * p - _qptGrid[1] - 1)/(2.0 * _qptGrid[1]);
          qpt.fpos[1] = q1;
          qpt.fpos[2] = q2;
          qpt.fpos[3] = q3;
          qpt.symop = 0;
          qpt.ibzqpt = q;
          qpt.irredQpt = q;
          _qptMap[p-1][r-1][s-1] = q;
          _qpoints[q] = qpt;
          _ibzqpts[q] = q;
          q++;
        }
      }
    }

    // Determine if the grid is gamma-centered and which dimensions are not
    bool gamma = true;
    xvector<double> shift(3);
    for (int i = 1; i < 4; i++) {
      if (_qptGrid[i] % 2 == 0) {
        gamma = false;
        shift[i] = _qpoints[0].fpos[i];
      }
    }
    // Center if necessary
    if (!gamma && force_gamma) {
      shiftMesh(shift);
      _shift = shift;
      gamma = true;
    }
    _isGammaCentered = gamma;
    _shifted = !aurostd::iszero(shift);  // ME190813

    // Obtain Cartesian coordinates
    if (!_shifted) {
      for (int q = 0; q < _nQPs; q++) {
        _qpoints[q].cpos = _recCell.f2c * _qpoints[q].fpos;
      }
    }
  }

  //shiftMesh/////////////////////////////////////////////////////////////////
  // Shifts the entire q-point mesh along a specific vector. This is useful
  // to center the q-point mesh around the Gamma point.
  void QMesh::shiftMesh(const xvector<double>& shift) {
    for (int q = 0; q < _nQPs; q++) {
      _qpoints[q].fpos -= shift;
      moveToBZ(_qpoints[q].fpos);
      _qpoints[q].cpos = _recCell.f2c * _qpoints[q].fpos;
    }
  }

  //moveToBZ//////////////////////////////////////////////////////////////////
  // Moves a q-point into the first Brillouin zone.
  // ME190702 - made more robust
  void QMesh::moveToBZ(xvector<double>& qpt) const {
    BringInCellInPlace(qpt, _ZERO_TOL_, 0.5, -0.5); //DX 20190905 - removed SYM namespace
  }

  //makeIrreducible///////////////////////////////////////////////////////////
  // Makes the q-point mesh irreducible
  // ME190813 - Changed algorithm to be much faster
  void QMesh::makeIrreducible() {
    if (_reduced) return;  // ME190701 - don't reduce if it's already reduced
    stringstream message;

    _ibzqpts.clear();
    _weights.clear();
    _nIQPs = 0;
    int nsym = (int) _recCell.pgroup.size();
    vector<vector<int> > irred_trans;
    vector<int> trans(nsym, -1);
    for (int q = 0; q < _nQPs; q++) {
      bool append = true;
      for (int sym = 0; sym < nsym; sym++) {
        for (int iq = 0; iq < _nIQPs; iq++) {
          if (irred_trans[iq][sym] == q) {
            append = false;
            _weights[iq]++;
            _qpoints[q].symop = sym;
            _qpoints[q].ibzqpt = iq;
            _qpoints[q].irredQpt = _ibzqpts[iq];
            sym = nsym;
            iq = _nIQPs;
          }
        }
      }
      if (append) {
        _qpoints[q].ibzqpt = _nIQPs;
        _ibzqpts.push_back(q);
        _weights.push_back(1);
        _nIQPs++;
        // Calculate the transformed irreducible q-point once to avoid repeated
        // matrix multiplications
        for (int sym = 0; sym < nsym; sym++) {
          trans[sym] = getQPointIndex(_recCell.pgroup[sym].Uf * _qpoints[q].fpos);
        }
        irred_trans.push_back(trans);
      }
    }
    _reduced = true;
    message << "Found " << _nIQPs << " irreducible qpoints.";
    pflow::logger(_AFLOW_FILE_NAME_, _APL_QMESH_MODULE_, message, _directory, *messageFile, std::cout);
  }

  // ME200109
  //calculateLittleGroups/////////////////////////////////////////////////////
  // Calculates little/small groups for each irreducible q-point.
  void QMesh::calculateLittleGroups() {
    _littleGroups.resize(_nIQPs, vector<int>(1, 0));  // Identity is always invariant
    uint nsymops = _recCell.pgroup.size();
    int q = -1;
    for (int iq = 0; iq < _nIQPs; iq++) {
      q = _ibzqpts[iq];
      const xvector<double>& fpos = _qpoints[q].fpos;
      for (uint isym = 1; isym < nsymops; isym++) {
        if (getQPointIndex(_recCell.pgroup[isym].Uf * fpos) == q) _littleGroups[iq].push_back(isym);
      }
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                            GETTER FUNCTIONS                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  int QMesh::getnIQPs() const {
    return _nIQPs;
  }

  int QMesh::getnQPs() const {
    return _nQPs;
  }

  int QMesh::getGrid(int i) const {
    return _qptGrid[i];
  }

  const xvector<int>& QMesh::getGrid() const {
    return _qptGrid;
  }

  const _qpoint& QMesh::getIrredQPoint(int i) const {
    return _qpoints[_qpoints[i].irredQpt];
  }

  const _qpoint& QMesh::getIrredQPoint(int i, int j, int k) const {
    return _qpoints[_qpoints[_qptMap[i][j][k]].irredQpt];
  }

  int QMesh::getIrredQPointIndex(int i) const {
    return _qpoints[i].irredQpt;
  }

  int QMesh::getIrredQPointIndex(int i, int j, int k) const {
    return _qpoints[_qptMap[i][j][k]].irredQpt;
  }

  vector<xvector<double> > QMesh::getIrredQPointsCPOS() const {
    vector<xvector<double> > cpos(_nIQPs, xvector<double>(3));
    for (int q = 0; q < _nIQPs; q++) {
      cpos[q] = _qpoints[_ibzqpts[q]].cpos;
    }
    return cpos;
  }

  vector<xvector<double> > QMesh::getIrredQPointsFPOS() const {
    vector<xvector<double> > fpos(_nIQPs, xvector<double>(3));
    for (int q = 0; q < _nIQPs; q++) {
      fpos[q] = _qpoints[_ibzqpts[q]].fpos;
    }
    return fpos;
  }

  const _qpoint& QMesh::getQPoint(int i) const {
    return _qpoints[i];
  }

  const _qpoint& QMesh::getQPoint(int i, int j, int k) const {
    return _qpoints[_qptMap[i][j][k]];
  }

  const _qpoint& QMesh::getQPoint(const xvector<double>& fpos) const {
    return _qpoints[getQPointIndex(fpos)];
  }

  // ME190813
  // Returns the index of the qpoint based on the fractional
  // position. It assumes that the point is already on the grid.
  int QMesh::getQPointIndex(xvector<double> fpos) const {
    // Shift back to original Monkhorst-Pack positions
    if (_shifted) fpos += _shift;
    moveToBZ(fpos);
    // invert Monkhorst-Pack formula;
    int p = (int) aurostd::nint((fpos[1] * 2 * _qptGrid[1] + _qptGrid[1] + 1)/2);
    int r = (int) aurostd::nint((fpos[2] * 2 * _qptGrid[2] + _qptGrid[2] + 1)/2);
    int s = (int) aurostd::nint((fpos[3] * 2 * _qptGrid[3] + _qptGrid[3] + 1)/2);
    return _qptMap[p - 1][r - 1][s - 1];
  }

  int QMesh::getQPointIndex(int i, int j, int k) const {
    return _qptMap[i][j][k];
  }

  vector<xvector<double> > QMesh::getQPointsCPOS() const {
    vector<xvector<double> > cpos(_nQPs, xvector<double>(3));
    for (int q = 0; q < _nQPs; q++) {
      cpos[q] = _qpoints[q].cpos;
    }
    return cpos;
  }

  vector<xvector<double> > QMesh::getQPointsFPOS() const {
    vector<xvector<double> > fpos(_nQPs, xvector<double>(3));
    for (int q = 0; q < _nQPs; q++) {
      fpos[q] = _qpoints[q].fpos;
    }
    return fpos;
  }

  int QMesh::getIbzqpt(int i) const {
    return _qpoints[i].ibzqpt;
  }

  int QMesh::getIbzqpt(int i, int j, int k) const {
    return _qpoints[_qptMap[i][j][k]].ibzqpt;
  }

  const vector<int>& QMesh::getIbzqpts() const {
    return _ibzqpts;
  }

  const vector<_qpoint>& QMesh::getPoints() const {
    return _qpoints;
  }

  const _kcell& QMesh::getReciprocalCell() const {
    return _recCell;
  }

  // ME190813
  bool QMesh::isShifted() const {
    return _shifted;
  }

  const xvector<double>& QMesh::getShift() const {
    return _shift;
  }

  const vector<int>& QMesh::getWeights() const {
    return _weights;
  }

  bool QMesh::isReduced() const {
    return _reduced;
  }

  bool QMesh::isGammaCentered() const {
    return _isGammaCentered;
  }

  // ME200109
  bool QMesh::littleGroupsCalculated() const {
    return _littleGroupsCalculated;
  }

  const vector<int>& QMesh::getLittleGroup(int iq) const {
    if (iq < _nIQPs) {
      return _littleGroups[iq];
    } else {
      string function = _APL_QMESH_ERR_PREFIX_ + "getLittleGroup()";
      string message = "Little groups are only calculated for irreducible q-points.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_RANGE_);
    }
  }

}  // namespace apl

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                FILE I/O                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace apl {

  //writeQpoints//////////////////////////////////////////////////////////////
  // Writes the Cartesian coordinates of each q-point into a file.
  void QMesh::writeQpoints(string filename, bool cartesian) {
    stringstream output;

    // Header
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    output << "[APL_QPOINTS]START" << std::endl;
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << std::setw(10) << "# Index";
    output << std::setw(20) << " ";
    output << "Q-points " << (cartesian?"(1/Angstrom)":"fractional") << std::endl;

    // Body
    for (int q = 0; q < _nQPs; q++) {
      output << std::setiosflags(std::ios::fixed | std::ios::right);
      output << std::setw(10) << q;
      for (int i = 1; i < 4; i++) {
        output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        output << std::setw(20) << std::setprecision(10) << std::scientific << (cartesian?_qpoints[q].cpos[i]:_qpoints[q].fpos[i]);
      }
      output << std::endl;
    }

    output << "[APL_QPOINTS]STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;

    // Write to file
    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      string function = _APL_QMESH_ERR_PREFIX_ + "writeQpoints";
      string message = "Could not write q-points to file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
    }
  }

  //writeIrredQpoints/////////////////////////////////////////////////////////
  // Writes the Cartesian coordinates and the multiplicity of the irreducible
  // q-points into a file.
  void QMesh::writeIrredQpoints(string filename, bool cartesian) {
    stringstream output;

    // Header
    output << AFLOWIN_SEPARATION_LINE << std::endl;
    output << "[APL_IRREDUCIBLE_QPOINTS]START" << std::endl;
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << std::setw(10) << "# Index";
    output << std::setiosflags(std::ios::fixed | std::ios::right);
    output << std::setw(15) << "Multiplicity";
    output << std::setw(20) << " ";
    output << "Q-points (1/Angstrom)" << std::endl;

    // Body
    for (int iq = 0; iq < _nIQPs; iq++) {
      int q = _ibzqpts[iq];
      output << std::setiosflags(std::ios::fixed | std::ios::right);
      output << std::setw(10) << q;
      output << std::setiosflags(std::ios::fixed | std::ios::right);
      output << std::setw(15) << _weights[iq];
      for (int i = 1; i < 4; i++) {
        output << std::setiosflags(std::ios::fixed | std::ios::showpoint | std::ios::right);
        output << std::setw(20) << std::setprecision(10) << std::scientific << (cartesian?_qpoints[q].cpos[i]:_qpoints[q].fpos[i]);
      }
      output << std::endl;
    }

    output << "[APL_IRREDUCIBLE_QPOINTS]STOP" << std::endl;
    output << AFLOWIN_SEPARATION_LINE << std::endl;

    // Write to file
    aurostd::stringstream2file(output, filename);
    if (!aurostd::FileExist(filename)) {
      string function = _APL_QMESH_ERR_PREFIX_ + "writeIrredQpoints";
      string message = "Could not write irreducible q-points to file.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_,function, message, _FILE_ERROR_);
    }
  }

} // namespace apl

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
//****************************************************************************
