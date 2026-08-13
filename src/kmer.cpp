#include "vdjtools/kmer.hpp"

#include <algorithm>
#include <atomic>
#include <stdexcept>
#include <mutex>
#include <thread>

namespace vdjtools {

namespace {

int popcount_mask(const std::vector<uint8_t>& pattern) {
    int k = 0;
    for (auto p : pattern) {
        if (p) ++k;
    }
    return k;
}

}  // namespace

int64_t kmer_code_space(const std::vector<uint8_t>& pattern, int32_t n_alphabet, int32_t n_v) {
    if (pattern.empty()) throw std::invalid_argument("pattern must have at least one position");
    if (n_alphabet < 2) throw std::invalid_argument("n_alphabet must be >= 2");
    if (n_v < 0) throw std::invalid_argument("n_v must be >= 0");
    const int k = popcount_mask(pattern);
    if (k == 0) throw std::invalid_argument("pattern must keep at least one position");

    // span = n_alphabet^k, then one V bucket beyond n_v for unresolved calls. Checked rather
    // than wrapped: a silent overflow here would alias distinct k-mers onto one code, which
    // reads downstream as a real (and very confident) shared feature.
    const int64_t kLimit = static_cast<int64_t>(1) << 62;
    int64_t span = 1;
    for (int i = 0; i < k; ++i) {
        if (span > kLimit / n_alphabet) throw std::overflow_error("k-mer code space overflows");
        span *= n_alphabet;
    }
    const int64_t buckets = static_cast<int64_t>(n_v) + 1;
    if (span > kLimit / buckets) throw std::overflow_error("k-mer code space overflows");
    return span * buckets;
}

KmerRow kmer_row(const std::vector<std::string>& junctions,
                 const std::vector<int32_t>& v_codes,
                 const std::vector<double>& weights,
                 const std::vector<uint8_t>& pattern,
                 const std::vector<int8_t>& alphabet,
                 int32_t n_alphabet, int32_t n_v, int32_t flank) {
    if (junctions.size() != v_codes.size() || junctions.size() != weights.size())
        throw std::invalid_argument("junctions, v_codes and weights must be the same length");
    if (alphabet.size() != 256) throw std::invalid_argument("alphabet must have 256 entries");
    if (flank < 0) throw std::invalid_argument("flank must be >= 0");

    const int width = static_cast<int>(pattern.size());
    const int k = popcount_mask(pattern);
    int64_t span = 1;
    for (int i = 0; i < k; ++i) span *= n_alphabet;   // validated by kmer_code_space
    kmer_code_space(pattern, n_alphabet, n_v);

    // Precompute the (offset, place-value) of each kept position so the inner loop skips gaps
    // instead of testing them.
    std::vector<std::pair<int, int64_t>> kept;
    kept.reserve(k);
    {
        int64_t mul = 1;
        for (int j = 0; j < width; ++j) {
            if (!pattern[j]) continue;
            kept.emplace_back(j, mul);
            mul *= n_alphabet;
        }
    }

    std::vector<std::pair<int64_t, double>> hits;
    hits.reserve(junctions.size() * 4);

    for (size_t i = 0; i < junctions.size(); ++i) {
        const std::string& s = junctions[i];
        // Signed throughout. The Python path had to cast before subtracting because a UInt32
        // length underflows to ~4e9 on a junction shorter than its flanks and then silently
        // yields a germline tail; the same trap exists here with size_t.
        const int64_t len = static_cast<int64_t>(s.size());
        const int64_t core_len = len - 2 * static_cast<int64_t>(flank);
        if (core_len < width) continue;

        const int64_t v = v_codes[i] < 0 ? static_cast<int64_t>(n_v)
                                         : static_cast<int64_t>(v_codes[i]);
        if (v > n_v) throw std::out_of_range("v_code >= n_v");
        const int64_t base = v * span;
        const char* core = s.data() + flank;
        const double w = weights[i];

        for (int64_t pos = 0; pos + width <= core_len; ++pos) {
            int64_t code = 0;
            bool ok = true;
            for (const auto& kp : kept) {
                const int8_t a = alphabet[static_cast<unsigned char>(core[pos + kp.first])];
                if (a < 0) { ok = false; break; }
                code += static_cast<int64_t>(a) * kp.second;
            }
            if (!ok) continue;
            hits.emplace_back(base + code, w);
        }
    }

    std::sort(hits.begin(), hits.end(),
              [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
                  return a.first < b.first;
              });

    KmerRow out;
    for (const auto& h : hits) {
        if (!out.codes.empty() && out.codes.back() == h.first) {
            out.weights.back() += h.second;
        } else {
            out.codes.push_back(h.first);
            out.weights.push_back(h.second);
        }
    }
    return out;
}

std::vector<KmerRow> kmer_rows(const std::vector<std::vector<std::string>>& junctions,
                               const std::vector<std::vector<int32_t>>& v_codes,
                               const std::vector<std::vector<double>>& weights,
                               const std::vector<uint8_t>& pattern,
                               const std::vector<int8_t>& alphabet,
                               int32_t n_alphabet, int32_t n_v, int32_t flank,
                               int32_t threads) {
    const size_t n = junctions.size();
    if (v_codes.size() != n || weights.size() != n)
        throw std::invalid_argument("junctions, v_codes and weights must be the same length");

    std::vector<KmerRow> out(n);
    if (n == 0) return out;

    unsigned nthreads = threads > 0 ? static_cast<unsigned>(threads)
                                    : std::thread::hardware_concurrency();
    if (nthreads == 0) nthreads = 1;
    nthreads = std::min<unsigned>(nthreads, static_cast<unsigned>(n));

    std::atomic<size_t> next{0};
    std::atomic<bool> failed{false};
    std::string error;
    std::mutex error_mutex;

    auto worker = [&]() {
        for (;;) {
            const size_t i = next.fetch_add(1);
            if (i >= n || failed.load()) return;
            try {
                out[i] = kmer_row(junctions[i], v_codes[i], weights[i], pattern, alphabet,
                                  n_alphabet, n_v, flank);
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!failed.exchange(true)) error = e.what();
                return;
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    for (unsigned t = 0; t < nthreads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    // Rethrow on the calling thread: a worker that died silently would leave an empty row, which
    // is indistinguishable from a sample that genuinely had no k-mers.
    if (failed.load()) throw std::runtime_error(error);
    return out;
}

std::vector<double> kmer_gather(const std::vector<std::string>& junctions,
                                const int32_t* v_codes,
                                const double* weights,
                                std::size_t n_clonotypes,
                                const std::vector<uint8_t>& pattern,
                                const std::vector<int8_t>& alphabet,
                                const int32_t* lookup, int64_t lookup_size,
                                int32_t n_alphabet, int32_t n_v, int32_t flank,
                                int32_t n_columns) {
    if (junctions.size() != n_clonotypes)
        throw std::invalid_argument("junctions, v_codes and weights must be the same length");
    if (alphabet.size() != 256) throw std::invalid_argument("alphabet must have 256 entries");
    if (flank < 0) throw std::invalid_argument("flank must be >= 0");
    if (n_columns < 0) throw std::invalid_argument("n_columns must be >= 0");

    const int64_t space = kmer_code_space(pattern, n_alphabet, n_v);
    if (lookup_size != space)
        throw std::invalid_argument("lookup must have one entry per code in the space");

    const int width = static_cast<int>(pattern.size());
    const int k = popcount_mask(pattern);
    int64_t span = 1;
    for (int i = 0; i < k; ++i) span *= n_alphabet;

    std::vector<std::pair<int, int64_t>> kept;
    kept.reserve(k);
    {
        int64_t mul = 1;
        for (int j = 0; j < width; ++j) {
            if (!pattern[j]) continue;
            kept.emplace_back(j, mul);
            mul *= n_alphabet;
        }
    }

    std::vector<double> out(static_cast<size_t>(n_columns), 0.0);
    for (size_t i = 0; i < junctions.size(); ++i) {
        const std::string& s = junctions[i];
        const int64_t len = static_cast<int64_t>(s.size());
        const int64_t core_len = len - 2 * static_cast<int64_t>(flank);
        if (core_len < width) continue;

        const int64_t v = v_codes[i] < 0 ? static_cast<int64_t>(n_v)
                                         : static_cast<int64_t>(v_codes[i]);
        if (v > n_v) throw std::out_of_range("v_code >= n_v");
        const int64_t base = v * span;
        const char* core = s.data() + flank;
        const double w = weights[i];

        for (int64_t pos = 0; pos + width <= core_len; ++pos) {
            int64_t code = 0;
            bool ok = true;
            for (const auto& kp : kept) {
                const int8_t a = alphabet[static_cast<unsigned char>(core[pos + kp.first])];
                if (a < 0) { ok = false; break; }
                code += static_cast<int64_t>(a) * kp.second;
            }
            if (!ok) continue;
            const int32_t col = lookup[static_cast<size_t>(base + code)];
            if (col >= 0) {
                if (col >= n_columns) throw std::out_of_range("lookup points past n_columns");
                out[static_cast<size_t>(col)] += w;
            }
        }
    }
    return out;
}

std::vector<int64_t> kmer_document_frequency(const std::vector<KmerRow>& rows, int64_t n_codes) {
    if (n_codes < 0) throw std::invalid_argument("n_codes must be >= 0");
    std::vector<int64_t> df(static_cast<size_t>(n_codes), 0);
    for (const auto& row : rows) {
        // Codes within a row are unique by construction, so one pass is the document count.
        for (int64_t c : row.codes) {
            if (c < 0 || c >= n_codes) throw std::out_of_range("code outside the declared space");
            ++df[static_cast<size_t>(c)];
        }
    }
    return df;
}

}  // namespace vdjtools
