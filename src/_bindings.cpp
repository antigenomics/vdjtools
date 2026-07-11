#include "vdjtools/core.hpp"
#include "vdjtools/model.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace vdjtools;

PYBIND11_MODULE(_core, m) {
    m.doc() = "vdjtools native core (C++): packed V(D)J model + Pgen hot loops.";

    m.def("hamming", &vdjtools::hamming, py::arg("a"), py::arg("b"),
          "Hamming distance between two equal-length strings; -1 if lengths differ.");
    m.def("version", &vdjtools::version, "Native core version string.");

    // PackedModel is built field-by-field from Python (see model/native.py) and passed to the
    // hot loops. Fields are read/write so the Python packer can populate them directly.
    py::class_<PackedModel>(m, "PackedModel")
        .def(py::init<>())
        .def_readwrite("vdj", &PackedModel::vdj)
        .def_readwrite("maxpal_v3", &PackedModel::maxpal_v3)
        .def_readwrite("maxpal_j5", &PackedModel::maxpal_j5)
        .def_readwrite("maxpal_d5", &PackedModel::maxpal_d5)
        .def_readwrite("maxpal_d3", &PackedModel::maxpal_d3)
        .def_readwrite("cut_v", &PackedModel::cut_v)
        .def_readwrite("cut_j", &PackedModel::cut_j)
        .def_readwrite("cut_d", &PackedModel::cut_d)
        .def_readwrite("func_v", &PackedModel::func_v)
        .def_readwrite("func_j", &PackedModel::func_j)
        .def_readwrite("func_d", &PackedModel::func_d)
        .def_readwrite("pv", &PackedModel::pv)
        .def_readwrite("pj", &PackedModel::pj)
        .def_readwrite("pjv", &PackedModel::pjv)
        .def_readwrite("pd_given_j", &PackedModel::pd_given_j)
        .def_readwrite("nbins_v", &PackedModel::nbins_v)
        .def_readwrite("nbins_j", &PackedModel::nbins_j)
        .def_readwrite("nbins_d5", &PackedModel::nbins_d5)
        .def_readwrite("nbins_d3", &PackedModel::nbins_d3)
        .def_readwrite("del_v", &PackedModel::del_v)
        .def_readwrite("del_j", &PackedModel::del_j)
        .def_readwrite("del_d", &PackedModel::del_d)
        .def_readwrite("ins_vd", &PackedModel::ins_vd)
        .def_readwrite("ins_dj", &PackedModel::ins_dj)
        .def_readwrite("ins_vj", &PackedModel::ins_vj)
        .def_readwrite("R_vd", &PackedModel::R_vd)
        .def_readwrite("R_dj", &PackedModel::R_dj)
        .def_readwrite("R_vj", &PackedModel::R_vj)
        .def_readwrite("bias_vd", &PackedModel::bias_vd)
        .def_readwrite("bias_dj", &PackedModel::bias_dj)
        .def_readwrite("bias_vj", &PackedModel::bias_vj);

    m.def("pgen_nt", &vdjtools::pgen_nt, py::arg("model"), py::arg("cdr3"),
          py::arg("v_idx") = -1, py::arg("j_idx") = -1,
          "Generation probability of an int-coded nt CDR3; v_idx/j_idx = -1 sums over all genes.");
    m.def("pgen_aa", &vdjtools::pgen_aa, py::arg("model"), py::arg("aa"),
          py::arg("v_idx") = -1, py::arg("j_idx") = -1,
          "Generation probability of an amino-acid CDR3; v_idx/j_idx = -1 sums over all genes.");
    m.def("pgen_aa_hamming1", &vdjtools::pgen_aa_hamming1, py::arg("model"), py::arg("aa"),
          py::arg("v_idx") = -1, py::arg("j_idx") = -1,
          "Total Pgen of the amino-acid CDR3 and all its Hamming-1 neighbours (one substitution).");

    py::class_<Counts>(m, "Counts")
        .def_readonly("v_choice", &Counts::v_choice)
        .def_readonly("j_choice", &Counts::j_choice)
        .def_readonly("d_gene", &Counts::d_gene)
        .def_readonly("v_3_del", &Counts::v_3_del)
        .def_readonly("j_5_del", &Counts::j_5_del)
        .def_readonly("d_del", &Counts::d_del)
        .def_readonly("ins_vd", &Counts::ins_vd)
        .def_readonly("ins_dj", &Counts::ins_dj)
        .def_readonly("ins_vj", &Counts::ins_vj)
        .def_readonly("dinucl_vd", &Counts::dinucl_vd)
        .def_readonly("dinucl_dj", &Counts::dinucl_dj)
        .def_readonly("dinucl_vj", &Counts::dinucl_vj);

    m.def("make_counts", &vdjtools::make_counts, py::arg("model"));
    m.def("estep_batch", &vdjtools::estep_batch, py::arg("model"), py::arg("seqs"),
          py::arg("vmasks"), py::arg("jmasks"), py::arg("dmasks"), py::arg("counts"),
          "One EM E-step: accumulate soft counts, return summed log-Pgen.");
}
