#include "vdjtools/core.hpp"
#include "vdjtools/inext.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = "vdjtools native core (C++). Pgen DP, generation sampler, and EM E-step "
              "inner loops land here (Phase 1); currently exposes the hamming primitive "
              "and the iNEXT size-based diversity kernel (curve + bootstrap + batch).";
    m.def("hamming", &vdjtools::hamming,
          py::arg("a"), py::arg("b"),
          "Hamming distance between two equal-length strings; -1 if lengths differ.");
    m.def("version", &vdjtools::version, "Native core version string.");

    // --- iNEXT size-based diversity kernel ---
    py::class_<vdjtools::InextCurve>(m, "InextCurve")
        .def_readonly("qD", &vdjtools::InextCurve::qD)
        .def_readonly("coverage", &vdjtools::InextCurve::coverage);

    py::class_<vdjtools::InextSample>(m, "InextSample")
        .def_readonly("qD", &vdjtools::InextSample::qD)
        .def_readonly("coverage", &vdjtools::InextSample::coverage)
        .def_readonly("se", &vdjtools::InextSample::se);

    m.def("inext_digamma", &vdjtools::digamma, py::arg("x"),
          "Digamma (psi) function; matches scipy.special.digamma.");

    m.def("inext_curve", &vdjtools::inext_curve,
          py::arg("counts"), py::arg("q_list"), py::arg("sizes"),
          py::call_guard<py::gil_scoped_release>(),
          "Deterministic size-based R/E point curve + sample coverage. Returns an "
          "InextCurve with .qD ([n_orders][n_sizes]) and .coverage ([n_sizes]).");

    m.def("inext_bootstrap", &vdjtools::inext_bootstrap,
          py::arg("counts"), py::arg("q_list"), py::arg("sizes"),
          py::arg("nboot"), py::arg("seed"),
          py::call_guard<py::gil_scoped_release>(),
          "Bootstrap standard errors of qD(m) via the augmented assemblage; returns "
          "an [n_orders][n_sizes] SE matrix.");

    m.def("inext_batch", &vdjtools::inext_batch,
          py::arg("samples"), py::arg("sample_sizes"), py::arg("q_list"),
          py::arg("nboot"), py::arg("seed"), py::arg("threads"),
          py::call_guard<py::gil_scoped_release>(),
          "Point curve + bootstrap SE for many samples, parallelized across samples. "
          "Returns a list of InextSample (one per input sample).");
}
