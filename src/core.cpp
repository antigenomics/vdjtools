#include "vdjtools/core.hpp"

namespace vdjtools {

int hamming(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return -1;
    int d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) ++d;
    }
    return d;
}

const char* version() { return "2.0.0"; }

}  // namespace vdjtools
