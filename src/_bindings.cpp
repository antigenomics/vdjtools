#include "vdjtools/core.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = "vdjtools native core (C++). Pgen DP, generation sampler, and EM E-step "
              "inner loop land here (Phase 1); currently exposes a testable primitive.";
    m.def("hamming", &vdjtools::hamming,
          py::arg("a"), py::arg("b"),
          "Hamming distance between two equal-length strings; -1 if lengths differ.");
    m.def("version", &vdjtools::version, "Native core version string.");
}
