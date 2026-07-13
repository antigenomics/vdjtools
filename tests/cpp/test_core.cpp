#include "vdjtools/core.hpp"

#include <cassert>
#include <cstdio>
#include <string>

int main() {
    using namespace vdjtools;
    assert(hamming("CASSL", "CASSL") == 0);
    assert(hamming("CASSL", "CASSF") == 1);
    assert(hamming("CASSLAP", "CFSSLAP") == 1);
    assert(hamming("CASS", "CASSL") == -1);   // length mismatch
    assert(std::string(version()) == "2.2.1");
    std::puts("ok");
    return 0;
}
