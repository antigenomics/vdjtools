#pragma once
#include <cstdint>
#include <string>
#include <vector>

// Native (V gene x k-mer) profile kernel.
//
// The Python path (features/kmer.py::_explode_kmers) materialises one row per k-mer occurrence
// and then groups it. That is the right shape for a tidy frame and the wrong shape for a
// corpus: a 100k-clonotype repertoire at k=4, flank=4 explodes to ~600k rows, and a 23k-sample
// reference corpus to ~1.4e10. The frame is never wanted here -- only the aggregated sparse row
// is -- so this builds that row directly.
//
// Two things make it cheap. K-mers are encoded as integers in a positional base-A numbering
// rather than hashed as strings, so the alphabet reduction (below) shrinks the code space
// instead of merely renaming keys; and the V gene is folded into the same integer, so the
// (V, k-mer) pair -- the feature a V+k-mer search is actually named for -- costs nothing extra.
//
// Gapped patterns are a per-position keep-mask over a window, so {1,1,0,1} is a gapped 3-mer
// spanning 4 residues. Ungapped k is just an all-ones mask; there is one code path.

namespace vdjtools {

// One sample's aggregated profile. `codes` is sorted ascending and unique; `weights` aligns.
struct KmerRow {
    std::vector<int64_t> codes;
    std::vector<double> weights;
};

// Number of distinct codes a (pattern, alphabet, V) configuration can produce. The caller needs
// this to size a vocabulary, and kmer_row needs it to fold V in; computing it in one place stops
// the two from disagreeing. Throws if the space would overflow int64.
int64_t kmer_code_space(const std::vector<uint8_t>& pattern, int32_t n_alphabet, int32_t n_v);

// Aggregate one repertoire.
//
//   junctions  per-clonotype junction_aa
//   v_codes    per-clonotype V index in [0, n_v); NEGATIVE means the V call did not resolve,
//              which is folded into its own bucket n_v rather than dropped -- an unresolvable
//              call is a fact about the sample, not a reason to silently lose its k-mers
//   weights    per-clonotype clone weight (already the caller's chosen ladder)
//   pattern    per-window-position keep mask; window width is pattern.size()
//   alphabet   256-entry char -> group id in [0, n_alphabet), or -1 for a residue we do not
//              model. A window containing one voids the whole window: a k-mer with an X in it
//              is not a k-mer, and substituting a wildcard would invent a count
//   flank      residues trimmed from EACH end before windowing, matching features/kmer.py.
//              A junction whose core is shorter than the window contributes nothing
KmerRow kmer_row(const std::vector<std::string>& junctions,
                 const std::vector<int32_t>& v_codes,
                 const std::vector<double>& weights,
                 const std::vector<uint8_t>& pattern,
                 const std::vector<int8_t>& alphabet,
                 int32_t n_alphabet, int32_t n_v, int32_t flank);

// Aggregate many repertoires, one thread per sample (0 = hardware concurrency). Sample order is
// preserved regardless of completion order.
std::vector<KmerRow> kmer_rows(const std::vector<std::vector<std::string>>& junctions,
                               const std::vector<std::vector<int32_t>>& v_codes,
                               const std::vector<std::vector<double>>& weights,
                               const std::vector<uint8_t>& pattern,
                               const std::vector<int8_t>& alphabet,
                               int32_t n_alphabet, int32_t n_v, int32_t flank,
                               int32_t threads);

// Document frequency: for each code, how many samples contain it at all. This is the corpus
// statistic a vocabulary cut and an IDF are built from, and it is a reduction over the same rows
// -- doing it here avoids handing Python 23k sparse rows just to count their supports.
std::vector<int64_t> kmer_document_frequency(const std::vector<KmerRow>& rows, int64_t n_codes);

// Accumulate one repertoire DIRECTLY onto a frozen vocabulary.
//
// This is the scoring hot path, and it is a different problem from kmer_row. Once the vocabulary
// is frozen, the full aggregation is waste: of a code space in the millions only the ~1e4-1e5
// selected columns are ever read, so sorting every occurrence to produce codes that are then
// discarded is most of the cost. Here each occurrence does one array lookup and one add -- no
// sort, no allocation proportional to the k-mer count.
//
//   lookup  code -> column index, or -1 for "not in the vocabulary". Length must be the code
//           space. Taken as a RAW POINTER, not a vector, and that is load-bearing: bound as
//           std::vector<int32_t> pybind11 copies the whole thing on every call, and a code space
//           of 8.16M turned a 26 ms kernel into a 124 ms one -- 5x slower than the polars path
//           it was written to beat. Read-only, so one numpy buffer is shared across every thread
//           and every sample and its memory is paid once for a whole corpus.
//
// Returns a dense vector of length n_columns. Dense because the vocabulary is already the
// sparsification: a selected column absent from this sample is a real zero, and the caller needs
// it to stay aligned with every other sample's vector.
std::vector<double> kmer_gather(const std::vector<std::string>& junctions,
                                const int32_t* v_codes,
                                const double* weights,
                                std::size_t n_clonotypes,
                                const std::vector<uint8_t>& pattern,
                                const std::vector<int8_t>& alphabet,
                                const int32_t* lookup, int64_t lookup_size,
                                int32_t n_alphabet, int32_t n_v, int32_t flank,
                                int32_t n_columns);

}  // namespace vdjtools
