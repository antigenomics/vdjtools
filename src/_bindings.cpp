#include "vdjtools/core.hpp"
#include "vdjtools/inext.hpp"
#include "vdjtools/kmer.hpp"
#include "vdjtools/model.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace vdjtools;

PYBIND11_MODULE(_core, m) {
    m.doc() = "vdjtools native core (C++): packed V(D)J model + Pgen / EM hot loops, and the "
              "iNEXT size-based diversity kernel (curve + bootstrap + batch).";

    m.def("version", &vdjtools::version, "Native core version string.");

    // --- Packed V(D)J recombination model + Pgen / EM (Phase 1) ---
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
        .def_readwrite("bias_vj", &PackedModel::bias_vj)
        .def_readwrite("dd", &PackedModel::dd)
        .def_readwrite("p_nd1", &PackedModel::p_nd1)
        .def_readwrite("p_nd2", &PackedModel::p_nd2)
        .def_readwrite("pd2_given_d1", &PackedModel::pd2_given_d1)
        .def_readwrite("del_d2", &PackedModel::del_d2)
        .def_readwrite("ins_dd", &PackedModel::ins_dd)
        .def_readwrite("R_dd", &PackedModel::R_dd)
        .def_readwrite("bias_dd", &PackedModel::bias_dd);

    m.def("pgen_nt", &vdjtools::pgen_nt, py::arg("model"), py::arg("cdr3"),
          py::arg("v_idx") = -1, py::arg("j_idx") = -1,
          py::call_guard<py::gil_scoped_release>(),
          "Generation probability of an int-coded nt CDR3; v_idx/j_idx = -1 sums over all genes.");
    m.def("pgen_aa", &vdjtools::pgen_aa, py::arg("model"), py::arg("aa"),
          py::arg("v_idx") = -1, py::arg("j_idx") = -1,
          py::call_guard<py::gil_scoped_release>(),
          "Generation probability of an amino-acid CDR3; v_idx/j_idx = -1 sums over all genes.");
    m.def("pgen_aa_hamming1", &vdjtools::pgen_aa_hamming1, py::arg("model"), py::arg("aa"),
          py::arg("v_idx") = -1, py::arg("j_idx") = -1,
          py::call_guard<py::gil_scoped_release>(),
          "Total Pgen of the amino-acid CDR3 and all its Hamming-1 neighbours (one substitution).");
    m.def("pgen_aa_batch", &vdjtools::pgen_aa_batch, py::arg("model"), py::arg("seqs"),
          py::arg("v_idxs") = std::vector<int>{}, py::arg("j_idxs") = std::vector<int>{},
          py::arg("mismatches") = 0, py::arg("threads") = 0,
          py::call_guard<py::gil_scoped_release>(),
          "Batch aa Pgen over many CDR3s, parallelized across sequences (mismatches=1 -> Hamming-1 "
          "ball). Bitwise-identical to per-sequence pgen_aa/pgen_aa_hamming1; threads=0 -> auto.");

    py::class_<vdjtools::AaScenario>(m, "AaScenario")
        .def_readonly("w", &vdjtools::AaScenario::w)
        .def_readonly("v", &vdjtools::AaScenario::v)
        .def_readonly("len_v", &vdjtools::AaScenario::len_v)
        .def_readonly("j", &vdjtools::AaScenario::j)
        .def_readonly("len_j", &vdjtools::AaScenario::len_j)
        .def_readonly("d", &vdjtools::AaScenario::d)
        .def_readonly("idx5", &vdjtools::AaScenario::idx5)
        .def_readonly("idx3", &vdjtools::AaScenario::idx3)
        .def_readonly("pos", &vdjtools::AaScenario::pos);
    m.def("best_aa_scenarios", &vdjtools::best_aa_scenarios, py::arg("model"), py::arg("aa"),
          py::arg("v_idx") = -1, py::arg("j_idx") = -1, py::arg("k") = 8,
          py::call_guard<py::gil_scoped_release>(),
          "Top-k recombination scenarios for an amino-acid CDR3 by joint max-product weight — the "
          "argmax counterpart of pgen_aa, over the same Pi_L*Pi_R transfer matrix.");

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
        .def_readonly("dinucl_vj", &Counts::dinucl_vj)
        .def_readonly("n_d", &Counts::n_d)
        .def_readonly("d2_gene", &Counts::d2_gene)
        .def_readonly("d2_del", &Counts::d2_del)
        .def_readonly("ins_dd", &Counts::ins_dd)
        .def_readonly("dinucl_dd", &Counts::dinucl_dd);

    m.def("make_counts", &vdjtools::make_counts, py::arg("model"));
    m.def("estep_batch", &vdjtools::estep_batch, py::arg("model"), py::arg("seqs"),
          py::arg("vmasks"), py::arg("jmasks"), py::arg("dmasks"), py::arg("counts"),
          py::arg("threads") = 0, py::arg("dd_allowed") = std::vector<int>{},
          py::call_guard<py::gil_scoped_release>(),
          "One EM E-step: accumulate soft counts, return summed log-Pgen. threads=0 -> auto; "
          "dd_allowed (per-read 0/1, empty=all) gates the n_D=2 tandem E-step.");

    // --- iNEXT size-based diversity kernel (Phase 2) ---
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

    // --- (V gene x k-mer) sparse profile kernel ---
    py::class_<KmerRow>(m, "KmerRow")
        .def_readonly("codes", &vdjtools::KmerRow::codes)
        .def_readonly("weights", &vdjtools::KmerRow::weights);

    m.def("kmer_code_space", &vdjtools::kmer_code_space,
          py::arg("pattern"), py::arg("n_alphabet"), py::arg("n_v"),
          "Size of the (V x k-mer) code space for a pattern/alphabet/V configuration. Raises "
          "rather than wrapping: an overflow would alias distinct k-mers onto one code.");

    m.def("kmer_row", &vdjtools::kmer_row,
          py::arg("junctions"), py::arg("v_codes"), py::arg("weights"), py::arg("pattern"),
          py::arg("alphabet"), py::arg("n_alphabet"), py::arg("n_v"), py::arg("flank"),
          py::call_guard<py::gil_scoped_release>(),
          "Aggregated sparse (V, k-mer) profile of one repertoire, without materialising the "
          "k-mer explosion. Returns a KmerRow with sorted unique .codes and aligned .weights.");

    m.def("kmer_rows", &vdjtools::kmer_rows,
          py::arg("junctions"), py::arg("v_codes"), py::arg("weights"), py::arg("pattern"),
          py::arg("alphabet"), py::arg("n_alphabet"), py::arg("n_v"), py::arg("flank"),
          py::arg("threads"),
          py::call_guard<py::gil_scoped_release>(),
          "kmer_row over many repertoires, parallelized across samples; input order preserved.");

    // `lookup`, `v_codes` and `weights` come in as numpy buffers and are read in place. Bound as
    // std::vector pybind11 would copy them on every call; for an 8.16M-entry lookup that alone
    // was 5x the kernel's own cost.
    m.def("kmer_gather",
          [](const std::vector<std::string>& junctions,
             py::array_t<int32_t, py::array::c_style | py::array::forcecast> v_codes,
             py::array_t<double, py::array::c_style | py::array::forcecast> weights,
             const std::vector<uint8_t>& pattern, const std::vector<int8_t>& alphabet,
             py::array_t<int32_t, py::array::c_style | py::array::forcecast> lookup,
             int32_t n_alphabet, int32_t n_v, int32_t flank, int32_t n_columns) {
              if (v_codes.ndim() != 1 || weights.ndim() != 1 || lookup.ndim() != 1)
                  throw std::invalid_argument("v_codes, weights and lookup must be 1-D");
              if (v_codes.shape(0) != weights.shape(0))
                  throw std::invalid_argument("v_codes and weights must be the same length");
              const int32_t* vptr = v_codes.data();
              const double* wptr = weights.data();
              const int32_t* lptr = lookup.data();
              const auto n = static_cast<std::size_t>(v_codes.shape(0));
              const auto lsize = static_cast<int64_t>(lookup.shape(0));
              py::gil_scoped_release release;
              return vdjtools::kmer_gather(junctions, vptr, wptr, n, pattern, alphabet, lptr,
                                           lsize, n_alphabet, n_v, flank, n_columns);
          },
          py::arg("junctions"), py::arg("v_codes"), py::arg("weights"), py::arg("pattern"),
          py::arg("alphabet"), py::arg("lookup"), py::arg("n_alphabet"), py::arg("n_v"),
          py::arg("flank"), py::arg("n_columns"),
          "Accumulate one repertoire directly onto a frozen vocabulary: one array lookup and "
          "one add per k-mer occurrence, no sort. The scoring hot path.");

    m.def("kmer_document_frequency", &vdjtools::kmer_document_frequency,
          py::arg("rows"), py::arg("n_codes"),
          py::call_guard<py::gil_scoped_release>(),
          "Per-code count of how many rows contain it -- the corpus statistic a vocabulary cut "
          "and an IDF are built from.");
}
