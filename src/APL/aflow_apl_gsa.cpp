#include "aflow_apl.h"

using namespace std;

namespace apl {

// ///////////////////////////////////////////////////////////////////////////

GeneralizedSupercellApproach::GeneralizedSupercellApproach(
    Supercell& sc, _xinput& xinput,
    _aflags& aflags, _kflags& kflags,
    _xflags& xflags, //_vflags& vflags, 
    string& AflowIn, ofstream& mf)
    : DirectMethodPC(sc, xinput, aflags, kflags, xflags, AflowIn, mf) {
}

// ///////////////////////////////////////////////////////////////////////////

GeneralizedSupercellApproach::~GeneralizedSupercellApproach() {
  clear();
}

// ///////////////////////////////////////////////////////////////////////////

void GeneralizedSupercellApproach::clear() {
}

//////////////////////////////////////////////////////////////////////////////

}  // namespace apl
