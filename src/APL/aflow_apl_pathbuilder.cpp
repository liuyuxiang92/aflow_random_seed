#include "aflow_apl.h"

using namespace std;

namespace apl {

  // ///////////////////////////////////////////////////////////////////////////

  PathBuilder::PathBuilder() {
    this->clear();
  }

  PathBuilder::PathBuilder(ModeEnumType mode) {
    this->clear();
    setMode(mode);
  }

  PathBuilder::~PathBuilder() {
    this->clear();
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::clear() {
    _mode = SINGLE_POINT_MODE;
    _store = CARTESIAN_LATTICE;
    _path.clear();
    _points.clear();
    _labels.clear();
    _pointsVectorDimension = 0;
    _pointsVectorStartingIndex = 0;
    _nPointsPerSubPath = 0;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::setMode(ModeEnumType mode) {
    _mode = mode;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::setStore(StoreEnumType store) {
    _store = store;
  }

  // ///////////////////////////////////////////////////////////////////////////
  //ME20190614
  const apl::PathBuilder::StoreEnumType& PathBuilder::getStore() const {
    return _store;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::addPoint(const string& l, int dim, ...) {
    va_list arguments;
    xvector<double> point(dim,1);

    va_start(arguments, dim);
    for(int i = 0; i < dim; i++) {
      point(i+1) = va_arg(arguments, double);
    }
    va_end(arguments);
    //
    addPoint(l, point);
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::addPoint(const string& l, const xvector<double>& p) {
    if( _points.empty() ) {
      _pointsVectorDimension = p.rows;
      _pointsVectorStartingIndex = p.lrows;
    }

    if( p.rows != _pointsVectorDimension ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::addPoint(); Wrong dimension of the point.");
      string function = "apl::PathBuilder::addPoint()";
      string message = "Wrong dimension of the point.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }

    if( p.lrows != _pointsVectorStartingIndex ) {
      xvector<double> pp(_pointsVectorDimension+_pointsVectorStartingIndex,_pointsVectorStartingIndex);
      for(int i = 0; i < _pointsVectorDimension; i++){
        pp(_pointsVectorStartingIndex+i) = p(p.lrows+i);
      }
      _points.push_back(pp);
    } else {
      _points.push_back(p);
    }

    _labels.push_back(l);
  }

  // ///////////////////////////////////////////////////////////////////////////

  int PathBuilder::getDensity() {
    return _nPointsPerSubPath;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::setDensity(int n) {
    if( n < 0 ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::setDensity(); The density should be >= 0.");
      string function = "apl::PathBuilder::setDensity()";
      string message = "The density should be >= 0.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }
    _nPointsPerSubPath = n;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::buildPath() {
    // Remove the old path
    _path.clear();

    // Quick solution...
    if( _points.empty() ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::buildPath; There are no points.");
      string function = "apl::PathBuilder::buildPath";
      string message ="There are no points.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    };

    // Create path in the SINGLE_POINT_MODE
    if( _mode == SINGLE_POINT_MODE ) {
      for(uint i = 1; i < _points.size(); i++) {
        xvector<double> startPoint = _points[i-1];
        _path.push_back(startPoint);
        if( _nPointsPerSubPath != 0 ) {
          xvector<double> dPoint = ( _points[i] - _points[i-1] ) * ( 1.0 / (_nPointsPerSubPath) );
          for(int j = 1; j <= (int)_nPointsPerSubPath; j++)
          {
            xvector<double> p = startPoint + j*dPoint;
            _path.push_back(p);
          }
        }
      }
      _path.push_back(_points[_points.size()-1]);
    }

    // Create path in the COUPLE_POINT_MODE
    if( _mode == COUPLE_POINT_MODE ) {
      // If the number of points is odd -> escape
      if( _points.size() % 2 != 0 ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::buildPath(); The number of points is odd.");
        string function = "apl::PathBuilder::buildPath()";
        string message = "The number of points is odd.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
      }

      xvector<double> startPoint(3), endPoint(3);
      endPoint = _points[0];
      for(uint i = 0; i < _points.size(); i+=2) {
        startPoint = _points[i];
        //if( aurostd::modulus( startPoint - endPoint ) > _AFLOW_APL_EPS_ )
        //    _path.push_back(endPoint);
        _path.push_back(startPoint);
        endPoint = _points[i+1];

        if( _nPointsPerSubPath != 0 ) {
          xvector<double> dPoint = ( endPoint - startPoint ) * ( 1.0 / (_nPointsPerSubPath) );
          for(int j = 1; j <= (int)_nPointsPerSubPath; j++) {
            xvector<double> p = startPoint + j*dPoint;
            _path.push_back(p);
          }
        }
      }
      //_path.push_back(endPoint);
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  uint PathBuilder::getPathSize() {
    return _path.size();
  }

  // ///////////////////////////////////////////////////////////////////////////

  uint PathBuilder::getPointSize() {
    if( _mode == SINGLE_POINT_MODE ) {
      return( _points.size() );
    }

    if( _mode == COUPLE_POINT_MODE ) {
      return( ( _points.size() / 2 ) + 1 );
    }

    //ME20191031 - use xerror
    //throw APLRuntimeError("apl::PathBuilder::getPointSize(); Unknown mode.");
    string function = "apl::PathBuilder::getPointSize()";
    string message = "Unknown mode.";
    throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
  }

  // ///////////////////////////////////////////////////////////////////////////

  double PathBuilder::getPathLength() {
    // Quick solution...
    if( _points.empty() ) return 0;

    //
    uint npaths;
    if( _mode == SINGLE_POINT_MODE )
      npaths = _points.size()-1;
    else if( _mode == COUPLE_POINT_MODE )
      npaths = _points.size()/2;
    else {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::getPointLength(); Unknown mode.");
      string function = "apl::PathBuilder::getPointLength()";
      string message = "Unknown mode.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
    }

    //
    double length = 0.0;

    for(uint i = 1; i <= npaths; i++) {
      length += getPathLength(i);
    }

    return length;
  }

  // ///////////////////////////////////////////////////////////////////////////

  double PathBuilder::getPathLength(uint i) {
    if( i <= 0 ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::getPathLength(); Wrong index. The index has to start from 1.");
      string function = "apl::PathBuilder::getPathLength()";
      string message = "Wrong index. The index has to start from 1.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_BOUNDS_);
    }

    // Quick solution 1...
    if( _points.empty() ) return 0;

    double length;
    if( _mode == SINGLE_POINT_MODE ) {
      // Quick solution 2...
      if( i > _points.size() ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPathLength(); Wrong index.");
        string function = "apl::PathBuilder::getPathLength()";
        string message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_BOUNDS_);
      }
      if( _store == RECIPROCAL_LATTICE ) {
        length = aurostd::modulus( F2C(trasp(reciprocalLattice),_points[i]) - F2C(trasp(reciprocalLattice),_points[i-1]) );
      } else {
        length = aurostd::modulus( _points[i] - _points[i-1] );
      }
      return length;
    }
    else if( _mode == COUPLE_POINT_MODE ) {
      // Quick solution 2...
      if( i > _points.size()/2 ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPathLength(); Wrong index.");
        string function = "apl::PathBuilder::getPathLength()";
        string message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_BOUNDS_);
      }
      if( _store == RECIPROCAL_LATTICE ) {
        length = aurostd::modulus( F2C(trasp(reciprocalLattice),_points[(i*2)-1]) - F2C(trasp(reciprocalLattice),_points[(i-1)*2]) );
      } else {
        length = aurostd::modulus( _points[(i*2)-1] - _points[(i-1)*2] );
      }
      return length;
    }

    //ME20191031 - use xerror
    //throw APLRuntimeError("apl::PathBuilder::getPathLength(); Unknown mode.");
    string function = "apl::PathBuilder::getPathLength()";
    string message = "Unknown mode.";
    throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
  }

  // ///////////////////////////////////////////////////////////////////////////

  xvector<double> PathBuilder::getPoint(uint i) {
    if( i <= 0 ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::getPoint(); Wrong index. The index has to start from 1.");
      string function = "apl::PathBuilder::getPoint()";
      string message = "Wrong index. The index has to start from 1.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_BOUNDS_);
    }
    if( _mode == SINGLE_POINT_MODE ) {
      if( i > _points.size() ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPoint(); Wrong index.");
        string function = "apl::PathBuilder::getPoint()";
        string message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_BOUNDS_);
      }
      return _points[i-1];
    } else if( _mode == COUPLE_POINT_MODE ) {
      if( i > (_points.size()/2)+1 ) {
        //ME20191031 - use xerro
        //throw APLRuntimeError("apl::PathBuilder::getPoint(); Wrong index.");
        string function = "apl::PathBuilder::getPoint()";
        string message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_BOUNDS_);
      }
      if( i == 1 ) return _points[0];
      if( i == (_points.size()/2)+1 ) return _points[_points.size()-1];
      return _points[2*i-3];
    }

    //ME20191031 - use xerror
    //throw APLRuntimeError("apl::PathBuilder::getPoint(); Unknown mode.");
    string function = "apl::PathBuilder::getPoint()";
    string message = "Unknown mode.";
    throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
  }

  // ///////////////////////////////////////////////////////////////////////////

  string PathBuilder::getPointLabel(uint i) {
    if( i <= 0 ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Wrong index. The index has to start from 1.");
      string function = "apl::PathBuilder::getPointLabel()";
      string message = "Wrong index. The index has to start from 1.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_BOUNDS_);
    }

    if( _mode == SINGLE_POINT_MODE ) {
      if( i > _labels.size() ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Wrong index.");
        string function = "apl::PathBuilder::getPointLabel()";
        string message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_BOUNDS_);
      }
      return _labels[i-1];
    } else if( _mode == COUPLE_POINT_MODE ) {
      if( i > (_labels.size()/2)+1 ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Wrong index.");
        string function = "apl::PathBuilder::getPointLabel()";
        string message = "Wrong index " + aurostd::utype2string<int>(i) + ".";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INDEX_BOUNDS_);
      }
      if( i == 1 ) return _labels[0];
      if( i == (_labels.size()/2)+1 ) return _labels[_labels.size()-1];

      int h = 2*i-2;
      int l = 2*i-3;

      if( _labels[l] == _labels[h] )
        return _labels[l];
      else
        return(_labels[l]+"|"+_labels[h]);
    }

    //ME20191031 - use xerror
    //throw APLRuntimeError("apl::PathBuilder::getPointLabel(); Unknown mode.");
    string function = "apl::PathBuilder::getPointLabel()";
    string message = "Unknown mode.";
    throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _VALUE_ILLEGAL_);
  }

  // ///////////////////////////////////////////////////////////////////////////

  vector< aurostd::xvector<double> > PathBuilder::getPath(ModeEnumType mode, const string& userPath) {
    vector< xvector<double> > new_points;
    vector<string> new_labels;
    vector<string> tokens;

    // Generate users path by points
    tokenize(userPath,tokens,"|-");
    for(uint i = 0; i < tokens.size(); i++) {
      uint j = 0;
      for( ; j < _labels.size(); j++) {
        // This works fine, but notorious changes of the Gamma point label
        // in aflow_kpoints.cpp lead me to do it like: G, Gamma, of \Gamma or
        // or anything like GAMMA \G or GaMmA or I do not know what else is
        // still a gamma point labeled as G!
        //if( _labels[j] == tokens.at(i) ) break;
        if( _labels[j].find(tokens.at(i)) != string::npos ) break;
      }
      if( j == _labels.size() ) {
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::getPath(); Undefined label of the point.");
        string function = "apl::PathBuilder::getPath()";
        string message = "Undefined label of the point.";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
      }

      new_points.push_back(_points[j]);      
      //new_labels.push_back(_labels[j]);
      // I will use my labels...
      new_labels.push_back(tokens.at(i));
    }

    // update
    _points.clear();
    _labels.clear();
    for(uint i = 0; i < new_points.size(); i++) {
      _points.push_back(new_points[i]);
      _labels.push_back(new_labels[i]);
    }
    new_points.clear();
    new_labels.clear();

    // Build full path
    _mode = mode;
    buildPath();

    //
    return _path;
  }

  // ///////////////////////////////////////////////////////////////////////////

  std::vector< aurostd::xvector<double> > PathBuilder::getPath() {
    buildPath();
    return _path;
  }

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::transform(const aurostd::xmatrix<double>& m) {
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="PathBuilder::transform():";
    if(LDEBUG){
      cerr << "PathBuilder::transform():" << " matrix:" << std::endl;
      cerr << m;
    }
    for(uint i = 0; i < _points.size(); i++) {
      if(LDEBUG){cerr << soliloquy << " IN: " << _labels[i] << ": " << _points[i] << std::endl;}
      _points[i] = m * _points[i];
      if(LDEBUG){cerr << soliloquy << " OUT: " << _labels[i] << ": " << _points[i] << std::endl;}
    }
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190219 - OBSOLETE duplicate
  //void PathBuilder::tokenize(const string& str,vector<string>& tokens, string del) {
  //  string delimiters = del;

  // Skip delimiters at beginning.
  //  string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // Find first "non-delimiter".
  //  string::size_type pos = str.find_first_of(delimiters, lastPos);

  //  while (string::npos != pos || string::npos != lastPos) {
  // Found a token, add it to the vector.
  //    tokens.push_back(str.substr(lastPos, pos - lastPos));

  // Skip delimiters. Note the "not_of"
  //    lastPos = str.find_first_not_of(delimiters, pos);

  // Find next "non-delimiter"
  //    pos = str.find_first_of(delimiters, lastPos);
  //  }
  //}

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::pointsAreDefinedFor(const xstructure& primitiveStructure, StoreEnumType store) {
    // Transform from the reciprocal lattice of the primitive cell to cartesian coordinates
    if( store == RECIPROCAL_LATTICE ) {
      transform( trasp(ReciprocalLattice(primitiveStructure.lattice)) );
    }

    //
    _store = store;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //void PathBuilder::transformPointsFor(const xstructure& supercellStructure, StoreEnumType store) {
  //  if( store == CARTESIAN_LATTICE ) {
  //    if( _store == RECIPROCAL_LATTICE )
  //      transform( trasp(ReciprocalLattice(supercellStructure.lattice)) );
  //    else
  //  }
  //  else if( store == RECIPROCAL_LATTICE ) {
  //    transform( inverse(trasp(ReciprocalLattice(supercellStructure.lattice))) );
  //  }

  //  _store = store;
  //  cartesianLattice = supercellStructure.lattice;
  //  reciprocalLattice = ReciprocalLattice(supercellStructure.lattice);

  //}

  // ///////////////////////////////////////////////////////////////////////////

  void PathBuilder::defineCustomPoints(const string& coords,const string& labels,const Supercell& sc,bool CARTESIAN_COORDS){ //CO20180409
    bool LDEBUG=(FALSE || XHOST.DEBUG);
    string soliloquy="PathBuilder::defineCustomPoints():";
    vector<string> vcoords,vlabels;
    vector<double> coordinate;
    aurostd::string2tokens(coords,vcoords,";");
    aurostd::string2tokens(labels,vlabels,",;");

    if(vcoords.empty()){
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::defineCustomPoints(); No input coordinates found");
      string function = "apl::PathBuilder::defineCustomPoints()";
      string message = "No input coordinates found";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ERROR_);
    }
    if(vlabels.empty()){
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::defineCustomPoints(); No input labels found");
      string function = "apl::PathBuilder::defineCustomPoints()";
      string message = "No input labels found";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ERROR_);
    }
    if(vcoords.size()!=vlabels.size()){
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::defineCustomPoints(); Size of input coordinates does not match size of input labels");
      string function = "apl::PathBuilder::defineCustomPoints()";
      string message = "Size of input coordinates does not match size of input labels";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ERROR_);
    }
    if(vcoords.size()<2){
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::defineCustomPoints(); Path size should include at least two points");
      string function = "apl::PathBuilder::defineCustomPoints()";
      string message = "Path size should include at least two points";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ERROR_);
    }

    for(uint i=0;i<vcoords.size();i++){
      aurostd::string2tokens<double>(vcoords[i],coordinate,",");
      if(coordinate.size()!=3){
        //ME20191031 - use xerror
        //throw APLRuntimeError("apl::PathBuilder::defineCustomPoints(); Input coordinate["+aurostd::utype2string(i)+"] is not dimension==3");
        string function = "apl::PathBuilder::defineCustomPoints()";
        string message = "Input coordinate["+aurostd::utype2string(i)+"] is not dimension==3";
        throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _INPUT_ERROR_);
      }
      if(LDEBUG){cerr << soliloquy << " found point " << vlabels[i] << ": (" << (CARTESIAN_COORDS?"cartesian":"fractional") << ") " << coordinate[0] << "," << coordinate[1] << "," << coordinate[2] << std::endl;}
      addPoint(vlabels[i],3,coordinate[0],coordinate[1],coordinate[2]);
    }

    if(!CARTESIAN_COORDS){
      // ME20200203 - these are custom points, i.e. they are user generated based
      // on a user-defined input structure. AFLOW should not switch the basis
      // behind the scenes. Assume that the user knows what they are doing.
      //transform( trasp(ReciprocalLattice(sc.getPrimitiveStructure().lattice)) ); //must be reciprocal
      transform( trasp(ReciprocalLattice(sc.getInputStructure().lattice)) );
    }

    //set lattices
    cartesianLattice = sc.getSupercellStructure().lattice;
    reciprocalLattice = ReciprocalLattice(sc.getSupercellStructure().lattice);

    setMode(SINGLE_POINT_MODE); //unless modified later
    _store = CARTESIAN_LATTICE;
  }

  void PathBuilder::takeAflowElectronicPath(const string& latticeName,const Supercell& sc) //CO20180409
    //const xstructure& inStructure,
    //const xstructure& scStructure) 
  {
    // Get path from electronic structure...
    stringstream fileKPOINTS;
    // input is not the reciprocal lattice!
    bool foundBZ;
    fileKPOINTS << LATTICE::KPOINTS_Directions(latticeName,sc.getInputStructure().lattice,10,sc.getInputStructure().iomode,foundBZ);
    fileKPOINTS.flush();

    if( !foundBZ ) {
      //ME20191031 - use xerror
      //throw APLRuntimeError("apl::PathBuilder::takeAflowElectronicPath(); The BZ not found for this lattice.");
      string function = "apl::PathBuilder::takeAflowElectronicPath()";
      string message = "The Brillouin zone was not found for this lattice.";
      throw aurostd::xerror(_AFLOW_FILE_NAME_, function, message, _RUNTIME_ERROR_);
    }

    //   cerr << fileKPOINTS.str() << std::endl;
    string line;
    vector<string> tokens;
    int nLine = 0;
    while( getline( fileKPOINTS, line ) ) {
      if( nLine++ < 4 ) continue;
      if( line.empty() ) continue;
      if( line.size() == 1 ) continue;
      //cout << line << std::endl;
      tokenize(line,tokens,string(" !"));
      addPoint(tokens[3],3,atof(tokens[0].c_str()),atof(tokens[1].c_str()),atof(tokens[2].c_str()));
      tokens.clear();
      line.clear();
    }
    fileKPOINTS.clear();

    //CO20180406 - with help from Mike Mehl
    //these coordinates are wrt PRIMITIVE lattice, so if input is not PRIMITIVE, the conversion will not work
    // Transform from the reciprocal lattice of the primitive cell to
    // cartesian coordinates (note; in klattice vectors are rows! )
    transform( trasp(ReciprocalLattice(sc.getPrimitiveStructure().lattice)) );  //inStructure.lattice)) );  //CO20180406

    //set lattices
    cartesianLattice = sc.getSupercellStructure().lattice;
    reciprocalLattice = ReciprocalLattice(sc.getSupercellStructure().lattice);

    // set mode
    setMode(COUPLE_POINT_MODE);
    _store = CARTESIAN_LATTICE;
  }

  // ///////////////////////////////////////////////////////////////////////////

  //ME20190614
  // Writes the k-points path into a VASP KPOINTS-formatted file
  xKPOINTS PathBuilder::createKPOINTS(const Supercell& sc) {
    vector<xvector<double> > points(2 * (_points.size() - 1), xvector<double>(3));
    vector<string> labels(2 * (_points.size() - 1));
    if (_mode == SINGLE_POINT_MODE) {
      points[0] = _points[0];
      labels[0] = _labels[0];
      for (uint i = 1; i < _points.size() - 1; i++) {
        points[2 * i - 1] = _points[i];
        points[2 * i] = _points[i];
        labels[2 * i - 1] = _labels[i];
        labels[2 * i] = _labels[i];
      }
      points.back() = _points.back();
      labels.back() = _labels.back();
    } else {
      points = _points;
      labels = _labels;
    }
    // ME20200117 - Convert to reciprocal coordinates of the
    // original structure or the distances will be wrong
    if (_store == CARTESIAN_LATTICE) {
      xmatrix<double> c2f = inverse(trasp(ReciprocalLattice(sc.getOriginalStructure().lattice)));
      for (uint i = 0; i < points.size(); i++) points[i] = c2f * points[i];
    } else {
      xmatrix<double> f2c = trasp(ReciprocalLattice(sc.getInputStructure().lattice));
      xmatrix<double> c2f_orig = inverse(trasp(ReciprocalLattice(sc.getOriginalStructure().lattice)));
      for (uint i = 0; i < points.size(); i++) points[i] = c2f_orig * f2c * points[i];
    }

    xKPOINTS xkpts;
    xkpts.is_KPOINTS_PATH = true;
    xkpts.is_KPOINTS_NNN = false;
    xkpts.vpath = labels;
    xkpts.vkpoints = points;
    xkpts.title = xkpts.createStandardTitlePath(sc.getPrimitiveStructure());
    xkpts.path_grid = _path.size()/(_points.size()/2);
    xkpts.grid_type = "Line-mode";
    xkpts.path_mode = "reciprocal";

    return xkpts;
  }

} // namespace apl
