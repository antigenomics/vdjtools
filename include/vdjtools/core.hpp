#pragma once
#include <string>

namespace vdjtools {

// Hamming distance between two equal-length strings; -1 if their lengths differ.
// Placeholder hot-path primitive; the Pgen DP / generation sampler / EM E-step
// (Phase 1) land alongside it in this native core.
int hamming(const std::string& a, const std::string& b);

// Native core version string (kept in sync with the Python package version).
const char* version();

}  // namespace vdjtools
