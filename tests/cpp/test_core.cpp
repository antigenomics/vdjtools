#include "vdjtools/core.hpp"

#include <cassert>
#include <cstdio>
#include <string>

int main() {
    using namespace vdjtools;
    assert(std::string(version()) == "2.7.0");
    std::puts("ok");
    return 0;
}
