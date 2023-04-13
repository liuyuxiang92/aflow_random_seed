// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by hagen.eckert@duke.edu

#ifndef _AUROSTD_DATA_CPP_
#define _AUROSTD_DATA_CPP_

#ifndef _AUROSTD_INCBIN_H_
#define INCBIN_PREFIX afdata
#define INCBIN_SILENCE_BITCODE_WARNING
#include "aurostd_incbin.h"

#endif // _AUROSTD_INCBIN_H_

// Include files either as text (INCTXT) or binary (INCBIN) into the aflow binary
// see AUROSTD/aurostd_incbin.h
extern "C" {
  // strings used for unittests
  INCTXT(UnitTest, "DATA/unit_test.json");
}

// Helper functions to retrieve embedded data //HE20230413
namespace aurostd{
  /// @brief return data for unittests
  /// @return content of `DATA/unit_test.json`
  JSON::object get_aflow_data_unit_test(){
    return JSON::loadString(afdataUnitTestData);
  }
}

#endif  //_AUROSTD_DATA_CPP_