

#ifndef _AUROSTD_DATA_CPP_
#define _AUROSTD_DATA_CPP_

#ifndef _AUROSTD_INCBIN_H_
#define INCBIN_PREFIX afdata
#define INCBIN_SILENCE_BITCODE_WARNING
#include "aurostd_incbin.h"

#endif // _AUROSTD_INCBIN_H_

// https://stackoverflow.com/a/13059195
struct membuf: std::streambuf {
  membuf(char const* base, size_t size) {char* p(const_cast<char*>(base));
    this->setg(p, p, p + size);
  }
};

struct imemstream: virtual membuf, std::istream {
  imemstream(char const* base, size_t size) : membuf(base, size), std::istream(static_cast<std::streambuf*>(this)) {
  }
};


extern "C" {
INCTXT(Prototest, "/Users/nathan/Projects/ecp/data/prototypes/Part1/AB_cP2_221_b_a/info.json");
INCBIN(char, Vlibs, "/Users/nathan/Projects/AFLOW-orig/src/Library_HTQC.txt.bz2");
}

namespace aurostd{
  string get_testdata(){
    return afdataPrototestData;
  }
  string get_vLIBS(){
//    std::streambuf user(afdatavLIBSData, afdatavLIBSSize);
    imemstream sbuf(&afdataVlibsData[0], afdataVlibsSize);
    //std::istream in(&sbuf);
    bxz::istream here(sbuf);
    string data;
    stream2string(here, data);
    cout << "hello" << endl;
    return data;
  }
}

#endif  //_AUROSTD_DATA_CPP_