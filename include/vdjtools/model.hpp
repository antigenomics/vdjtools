#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace vdjtools {

// A V(D)J recombination model packed into contiguous arrays for the native Pgen / generation /
// EM hot loops. Built once from the Python polars model (see python/vdjtools/model/native.py);
// the field layout mirrors OLGA's processed arrays so the port matches to numerical tolerance.
//
// Conventions (identical to the Python reference):
//   - nucleotides int-coded A,C,G,T = 0..3
//   - deletion arrays are dense, indexed by  array_idx = ndel + max_palindrome
//   - dinucleotide R is stored row-major [to*4 + from] = P(next=to | prev=from) (column-stochastic)
struct PackedModel {
    bool vdj = false;
    int maxpal_v3 = 0, maxpal_j5 = 0, maxpal_d5 = 0, maxpal_d3 = 0;

    // Germline: palindrome-extended CDR3-region cut segments, pre-encoded to 0..3.
    std::vector<std::vector<int8_t>> cut_v, cut_j, cut_d;
    std::vector<int> func_v, func_j, func_d;  // indices of functional (usable) genes

    // Gene choice.
    std::vector<double> pv;          // [nV]
    std::vector<double> pj;          // VDJ: P(J) [nJ]
    std::vector<double> pjv;         // VJ:  P(J|V) [nV*nJ] row-major
    std::vector<double> pd_given_j;  // VDJ: P(D|J) [nJ*nD] row-major

    // Deletions (dense, array_idx = ndel + max_palindrome).
    int nbins_v = 0, nbins_j = 0, nbins_d5 = 0, nbins_d3 = 0;
    std::vector<double> del_v;   // [nV*nbins_v]
    std::vector<double> del_j;   // [nJ*nbins_j]
    std::vector<double> del_d;   // [nD*nbins_d5*nbins_d3]

    // Insertions + dinucleotide Markov (R row-major [to*4+from], bias = steady state).
    std::vector<double> ins_vd, ins_dj, ins_vj;
    std::vector<double> R_vd, R_dj, R_vj;        // each length 16
    std::vector<double> bias_vd, bias_dj, bias_vj;  // each length 4

    // D-D (tandem, n_D=2) extension — populated only when `dd` is true (P(n_D=2)>0). D2 draws
    // from the same D germline (cut_d) as D1, so it reuses nbins_d5/nbins_d3 and maxpal_d5/d3.
    bool dd = false;
    double p_nd1 = 1.0, p_nd2 = 0.0;        // P(n_D=1) (0-D folds in via a fully-trimmed D), P(n_D=2)
    std::vector<double> pd2_given_d1;       // [nD*nD] row-major, P(D2|D1)
    std::vector<double> del_d2;             // [nD*nbins_d5*nbins_d3], second-D joint 5'/3' deletion
    std::vector<double> ins_dd, R_dd, bias_dd;  // DD insertion junction (same layout as vd/dj)

    int nV() const { return static_cast<int>(cut_v.size()); }
    int nJ() const { return static_cast<int>(cut_j.size()); }
    int nD() const { return static_cast<int>(cut_d.size()); }
};

// Generation probability of a nucleotide CDR3, optionally restricted to a V and/or J (index into
// the gene lists; -1 = sum over all functional genes of that segment). Matches OLGA / the Python
// reference exactly.
double pgen_nt(const PackedModel& m, const std::vector<int8_t>& cdr3, int v_idx, int j_idx);

// Generation probability of an amino-acid CDR3 (codon-marginalizing transfer matrix). ``aa`` is the
// CDR3 amino-acid string; v_idx/j_idx as for ``pgen_nt`` (-1 = agnostic / sum over all genes).
double pgen_aa(const PackedModel& m, const std::string& aa, int v_idx, int j_idx);

// Total Pgen of ``aa`` and all its Hamming-distance-1 amino-acid neighbours (one substitution),
// via the inclusion-exclusion identity over per-position wildcards. v_idx/j_idx as above.
double pgen_aa_hamming1(const PackedModel& m, const std::string& aa, int v_idx, int j_idx);

// EM soft counts — one accumulator per event realization, laid out like the PackedModel prob
// arrays so the Python M-step can renormalize them directly.
struct Counts {
    std::vector<double> v_choice, j_choice, d_gene, v_3_del, j_5_del, d_del;
    std::vector<double> ins_vd, ins_dj, ins_vj, dinucl_vd, dinucl_dj, dinucl_vj;
};

// A zeroed :class:`Counts` sized to ``m``.
Counts make_counts(const PackedModel& m);

// One EM E-step over a batch of CDR3s: accumulate soft counts into ``counts`` and return the
// summed log-Pgen (over scoreable reads). Per-read ``vmasks[i]`` / ``jmasks[i]`` / ``dmasks[i]``
// are gene-index lists restricting enumeration (empty => all functional genes of that segment).
double estep_batch(const PackedModel& m,
                   const std::vector<std::vector<int8_t>>& seqs,
                   const std::vector<std::vector<int>>& vmasks,
                   const std::vector<std::vector<int>>& jmasks,
                   const std::vector<std::vector<int>>& dmasks,
                   Counts& counts);

}  // namespace vdjtools
