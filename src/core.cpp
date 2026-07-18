#include "vdjtools/core.hpp"

#ifndef VDJTOOLS_VERSION
// CMake injects the real version (parsed from pyproject.toml) as a PUBLIC compile definition;
// this guard only covers an unusual build that compiles core.cpp without CMakeLists.txt.
#define VDJTOOLS_VERSION "0.0.0+unknown"
#endif

namespace vdjtools {

const char* version() { return VDJTOOLS_VERSION; }

}  // namespace vdjtools
