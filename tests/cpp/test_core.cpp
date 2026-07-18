#include "vdjtools/core.hpp"

#include <cassert>
#include <cstdio>
#include <string>

// VDJTOOLS_VERSION is the single source of truth (parsed from pyproject.toml by CMake and
// injected into vdjtools_core as a PUBLIC compile definition). Asserting against the macro
// rather than a hand-copied literal verifies that version() actually returns the compiled-in
// version, and can never drift out of sync on a release bump — the failure this replaces.
int main() {
    using namespace vdjtools;
    assert(std::string(version()) == VDJTOOLS_VERSION);
    assert(std::string(version()).find('.') != std::string::npos);  // looks like X.Y.Z
    std::puts("ok");
    return 0;
}
