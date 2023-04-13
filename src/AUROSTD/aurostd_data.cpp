

#ifndef _AUROSTD_DATA_CPP_
#define _AUROSTD_DATA_CPP_

#ifndef _AUROSTD_INCBIN_H_
#define INCBIN_PREFIX afdata
#define INCBIN_SILENCE_BITCODE_WARNING
#include "aurostd_incbin.h"

#endif // _AUROSTD_INCBIN_H_

extern "C" {
INCTXT(UnitTest, "DATA/unit_test.json");
}

namespace aurostd{
  JSON::object get_unit_test_data(){
    return JSON::loadString(afdataUnitTestData);
  }
}

#endif  //_AUROSTD_DATA_CPP_