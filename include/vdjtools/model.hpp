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

// Translate junction nucleotides the way the legacy converters do, and collapse the non-coding
// runs, for a whole column at once: the ``to_unified_cdr3aa(translate(nt))`` every format reader
// applies. An in-frame sequence is a plain codon walk with ``*`` for a stop and ``X`` for a codon
// that is not clean ACGT. An OUT-OF-FRAME one is padded in the middle with ``?`` and translated
// inward from both ends, leaving the untranslatable middle lower-cased -- then every run of
// lower-case nucleotides and of ``# ~ _ ?`` collapses to a single ``_``.
//
// It lives here because it needs the genetic code table, and it is batched because it was the last
// per-row Python left in the readers: a 42,877-row immunoSEQ export walks ~40,000 codon strings.
std::vector<std::string> translate_junctions(const std::vector<std::string>& seqs, int nthreads);

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

// The degenerate-sequence DP the two above are special cases of. ``allowed[c]`` is a 64-bit mask
// over codon indices (a*16+b*4+c, nt A,C,G,T = 0..3) accepted at amino-acid position c, so a
// position may permit any residue subset at no extra cost — the transfer matrix contracts over the
// allowed codons either way. ``L`` is the length in amino acids; v_idx/j_idx as above.
double pgen_aa_masked(const PackedModel& m, const uint64_t* allowed, int L, int v_idx, int j_idx);

// Total Pgen of every amino-acid CDR3 matching a motif given as per-position residue sets:
// ``allowed[c]`` is the string of residues permitted at position c, where an empty string or one
// containing 'X' means a wildcard (any of the 20 amino acids). One masked DP pass whatever the sets
// contain — the sequences the motif matches are never enumerated. Throws std::invalid_argument on a
// character the genetic code does not name (an empty mask would score a silent 0). v/j as above.
double pgen_aa_degenerate(const PackedModel& m, const std::vector<std::string>& allowed,
                          int v_idx, int j_idx);

// Batch ``pgen_aa_degenerate``, parallelized across queries like ``pgen_aa_batch`` and
// bitwise-identical to the per-query calls. Empty v_idxs/j_idxs → all -1 (gene-agnostic).
std::vector<double> pgen_aa_degenerate_batch(const PackedModel& m,
                                             const std::vector<std::vector<std::string>>& allowed,
                                             const std::vector<int>& v_idxs,
                                             const std::vector<int>& j_idxs, int nthreads);

// Batch nucleotide Pgen over many sequences, parallelized across sequences exactly as
// ``pgen_aa_batch`` is, and bitwise-identical to the per-sequence ``pgen_nt`` calls for any
// ``nthreads``. Sequences are validated and encoded serially up front: a worker may not throw, and
// a non-ACGT base has to be an error rather than a silent 0. Per-sequence v_idxs/j_idxs; empty →
// all -1 (gene-agnostic). nthreads=0 → auto (hw-2); batches under 64 stay single-threaded.
std::vector<double> pgen_nt_batch(const PackedModel& m, const std::vector<std::string>& seqs,
                                  const std::vector<int>& v_idxs, const std::vector<int>& j_idxs,
                                  int nthreads);

// Batch aa Pgen over many sequences, parallelized across sequences (mismatches: 0 = exact,
// 1 = Hamming-1 ball). Bitwise-identical to the per-sequence calls. Per-sequence v/j indices;
// empty v_idxs/j_idxs → all -1 (gene-agnostic). nthreads=0 → auto (hw-2).
std::vector<double> pgen_aa_batch(const PackedModel& m, const std::vector<std::string>& seqs,
                                  const std::vector<int>& v_idxs, const std::vector<int>& j_idxs,
                                  int mismatches, int nthreads);

// One recombination scenario for an amino-acid CDR3: the argmax counterpart of what ``pgen_aa``
// sums over. Coordinates are 0-based, half-open, in CDR3 nucleotide space. ``d < 0`` for a VJ
// chain (then idx5/idx3/pos are unused). ``w`` is the joint max weight P(nt, scenario) EXCLUDING
// nothing — it is directly comparable across scenarios of the same sequence.
struct AaScenario {
    double w = 0.0;
    int v = -1, len_v = 0, j = -1, len_j = 0;
    int d = -1, idx5 = 0, idx3 = 0, pos = 0;
};

// Top-``k`` scenarios by joint max-product weight for an amino-acid CDR3 — the same Pi_L*Pi_R
// transfer-matrix DP ``pgen_aa`` uses, with ``max`` in place of the sums and the winning
// (V, delV) / (J, delJ) carried through the state.
//
// Returns SCENARIOS, not nucleotides: recovering the nt string is one cheap per-scenario DP over
// the free positions (see vdjtools.model.viterbi), so the expensive search need not carry paths.
//
// Tandem-D (n_D=2) is deliberately not enumerated: a single D trimmed to zero length already
// reaches every middle, so D-D can only reorder candidates, never add one.
std::vector<AaScenario> best_aa_scenarios(const PackedModel& m, const std::string& aa,
                                          int v_idx, int j_idx, int k);

// Batch ``best_aa_scenarios`` over many CDR3s, parallelized across sequences exactly as
// ``pgen_aa_batch`` is: disjoint work, no reduction, so the result is identical to the per-sequence
// calls for any ``nthreads``. Returned as parallel COLUMNS with one entry per scenario (not per
// query), so the Python side builds a frame without materializing k objects per row. ``row`` indexes
// ``seqs``, ``rank`` is the 0-based position within that row's top-k; a query the DP cannot explain
// (unknown residue, or no reachable rearrangement under the pinned V/J) contributes no entries.
struct AaScenarioBatch {
    std::vector<int> row, rank;
    std::vector<double> w;
    std::vector<int> v, len_v, j, len_j, d, idx5, idx3, pos;
};

AaScenarioBatch best_aa_scenarios_batch(const PackedModel& m,
                                        const std::vector<std::string>& seqs,
                                        const std::vector<int>& v_idxs,
                                        const std::vector<int>& j_idxs, int k, int nthreads);

// One row of :func:`infer_nt` -- the most likely nucleotide CDR3 for an amino-acid one, with the
// scenario behind it. ``ok`` is false for a row no rearrangement explains (an unknown residue, or
// nothing reachable under the pinned V/J); every other field is then unset.
struct InferNt {
    bool ok = false;
    std::string nt;
    int v = -1, len_v = 0, j = -1, len_j = 0, d = -1, d_start = 0, d_end = 0;
    int n_candidates = 0;
    double pgen = 0.0, scenario_p = 0.0, runner_up_pgen = 0.0;
};

// ``infer_nt`` over a whole table: scenario search, codon reconstruction and the exact marginal
// re-score, all inside one call and parallelized across rows.
//
// This is the entry point a caller with more than one junction wants. Everything the Python side
// used to do per candidate -- laying a scenario out and picking its codons -- is the DP below, so
// the GIL is released for the whole batch rather than reacquired once per row.
//
// ``v_alts[i]`` / ``j_alts[i]`` are the gene indices to search for row i: one entry pins that
// allele, several search all of them and let the model choose, and ``-1`` marginalizes over the
// segment. An EMPTY list is not ``-1``: it names no allele, so it explains nothing and the row is
// declined -- as does an empty ``aas[i]``, which is how a caller passes through a row it has
// already refused. ``k`` is how many distinct candidates are re-scored with the exact ``pgen_nt``.
std::vector<InferNt> infer_nt_batch(const PackedModel& m, const std::vector<std::string>& aas,
                                    const std::vector<std::vector<int>>& v_alts,
                                    const std::vector<std::vector<int>>& j_alts,
                                    int k, int nthreads);

// One batch of ancestral draws from the model -- the native counterpart of
// ``vdjtools.model.generate``. Columns, one entry per returned sequence.
struct GenBatch {
    std::vector<std::string> nt, aa;          // the CDR3/junction nucleotides and their translation
    std::vector<int> v, d, d2, j;             // gene indices; d/d2 = -1 when none contributed
    std::vector<int8_t> productive;           // in-frame and stop-free
};

// Sample ``n`` recombinations: pick genes, deletions, insertion lengths and non-templated
// nucleotides (through the dinucleotide Markov chain), assemble the CDR3.
//
// **Every draw is seeded from ``(seed, slot)``**, so the result depends on the seed and the slot
// index alone -- never on ``nthreads``. That is a correctness property rather than a nicety: a
// corpus whose contents depend on the builder's core count cannot be compared with one built
// anywhere else. ``productive_only`` re-draws a slot (advancing that slot's own stream) until the
// draw is in-frame and stop-free, and throws once a slot exhausts ``kDrawAttempts``.
//
// NOTE: this is a DIFFERENT random stream from the Python reference sampler, which draws from
// numpy's PCG64. The two agree in distribution, not sequence by sequence.
GenBatch generate_batch(const PackedModel& m, int n, uint64_t seed, bool productive_only,
                        int nthreads);

// Alignment-vote seeding for EM: each read votes its longest-matching V and J germline, the vote
// split across every gene tied for the longest match. Returns ``(v_votes, j_votes, vj_votes)``,
// the last a flat [nV*nJ] joint filled only when ``joint`` is set (a VJ chain needs P(J|V)).
//
// Deliberately single-threaded: the votes accumulate in read order, which makes the result
// bitwise-identical to the Python reference rather than merely equal to tolerance. It is ~50-70
// microseconds per read there (nV + nJ germline comparisons, each a Python character loop) and it
// runs before EVERY fit, native path included.
//
// Splitting a tie matters and is not a detail: germline-identical paralogs (TRBV6-2/6-5/6-6,
// IGKV2-28/2D-28) tie exactly, and handing the whole family to one representative seeds the rest
// at P(V)=0 -- which the E-step's zero-probability skip then makes absorbing.
struct AlignVotes {
    std::vector<double> v, j, vj;
};

AlignVotes align_init_votes(const PackedModel& m, const std::vector<std::string>& seqs, bool joint);

// EM soft counts — one accumulator per event realization, laid out like the PackedModel prob
// arrays so the Python M-step can renormalize them directly.
struct Counts {
    std::vector<double> v_choice, j_choice, d_gene, v_3_del, j_5_del, d_del;
    std::vector<double> ins_vd, ins_dj, ins_vj, dinucl_vd, dinucl_dj, dinucl_vj;
    // D-D (tandem) soft counts — sized/populated only for a `dd` model. n_d[k] = P(n_D=k) mass
    // (k in {1,2}); d2_gene = P(D2|D1) [nD*nD]; d2_del [nD*nbins_d5*nbins_d3]; ins_dd/dinucl_dd
    // like the vd junction. Single-D E-step fills n_d[1] only (renormalizes to delta(1), a no-op).
    std::vector<double> n_d, d2_gene, d2_del, ins_dd, dinucl_dd;
};

// A zeroed :class:`Counts` sized to ``m``.
Counts make_counts(const PackedModel& m);

// One EM E-step over a batch of CDR3s: accumulate soft counts into ``counts`` and return the
// summed log-Pgen (over scoreable reads). Per-read ``vmasks[i]`` / ``jmasks[i]`` / ``dmasks[i]``
// are gene-index lists restricting enumeration (empty => all functional genes of that segment).
// ``nthreads`` partitions the reads across worker threads (each with a private accumulator, reduced
// in a fixed order); 0 = auto (hardware_concurrency - 2). Small batches run single-threaded so their
// result stays bitwise-identical. The soft counts are exact regardless of thread count (up to the
// float summation order, which the fixed reduction keeps deterministic per ``nthreads``).
// ``dd_allowed`` optionally gates the tandem (n_D=2) E-step per read (1 = may be tandem, 0 = single-D
// only); empty = every read may be tandem. Used to anchor D-D learning to arda's tandem calls.
double estep_batch(const PackedModel& m,
                   const std::vector<std::vector<int8_t>>& seqs,
                   const std::vector<std::vector<int>>& vmasks,
                   const std::vector<std::vector<int>>& jmasks,
                   const std::vector<std::vector<int>>& dmasks,
                   Counts& counts,
                   int nthreads = 0,
                   const std::vector<int>& dd_allowed = {});

}  // namespace vdjtools
