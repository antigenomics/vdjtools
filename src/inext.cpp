#include "vdjtools/inext.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <random>
#include <thread>

// Native iNEXT size-based engine. See include/vdjtools/inext.hpp and the Python
// reference python/vdjtools/stats/inext.py — every function below mirrors a
// helper there (name in parentheses).

namespace vdjtools {
namespace {

constexpr double NEG_INF = -std::numeric_limits<double>::infinity();

// log C(a, b) via lgamma; -inf when b < 0 or a < b (so exp() -> 0). (``_logbinom``)
inline double logbinom(double a, double b) {
    if (b < -1e-12 || a - b < -1e-12) return NEG_INF;
    return std::lgamma(a + 1.0) - std::lgamma(b + 1.0) - std::lgamma(a - b + 1.0);
}

// Precomputed abundance spectrum + summary statistics for one assemblage.
struct Sample {
    Vec x;        // positive counts
    Vec c, f;     // distinct counts (ascending) and their multiplicities f_k
    double n = 0; // total individuals
    double f1 = 0, f2 = 0;
    int S = 0;    // observed richness
};

Sample make_sample(const Vec& raw) {
    Sample s;
    s.x.reserve(raw.size());
    for (double v : raw)
        if (v > 0) s.x.push_back(v);
    Vec sorted = s.x;
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t i = 0; i < sorted.size();) {
        std::size_t j = i;
        while (j < sorted.size() && sorted[j] == sorted[i]) ++j;
        s.c.push_back(sorted[i]);
        s.f.push_back(static_cast<double>(j - i));
        i = j;
    }
    double n = 0;
    for (double v : s.x) n += v;
    s.n = n;
    s.S = static_cast<int>(s.x.size());
    for (std::size_t i = 0; i < s.c.size(); ++i) {
        if (s.c[i] == 1.0) s.f1 = s.f[i];
        else if (s.c[i] == 2.0) s.f2 = s.f[i];
    }
    return s;
}

// Empirical (MLE plug-in) Hill number of order q on x/n. (``_plugin``)
double plugin(const Sample& s, int q) {
    const double n = s.n;
    if (q == 0) return static_cast<double>(s.S);
    if (q == 1) {
        double h = 0.0;
        for (double xi : s.x) {
            double p = xi / n;
            h -= p * std::log(p);
        }
        return std::exp(h);
    }
    double sum = 0.0;
    for (double xi : s.x) sum += std::pow(xi / n, q);
    return std::exp(std::log(sum) / (1 - q));
}

// Asymptotic Hill-number estimator of order q. (``_asymptotic``)
double asymptotic(const Sample& s, int q) {
    const double n = s.n, f1 = s.f1, f2 = s.f2;
    if (q == 0) {
        double f0 = (f2 > 0) ? f1 * f1 / 2.0 / f2 : f1 * (f1 - 1.0) / 2.0;
        return s.S + (n - 1) / n * f0;
    }
    if (q == 1) {
        double first = 0.0;
        double dig_n = digamma(n);
        for (std::size_t i = 0; i < s.c.size(); ++i)
            first += s.f[i] * s.c[i] / n * (dig_n - digamma(s.c[i]));
        double a_cwj;
        if (f2 > 0) a_cwj = 2 * f2 / ((n - 1) * f1 + 2 * f2);
        else if (f1 > 0) a_cwj = 2 / ((n - 1) * (f1 - 1) + 2);
        else a_cwj = 1.0;
        double second = 0.0;
        if (!(f1 == 0.0 || a_cwj == 1.0)) {
            double base = 1 - a_cwj;
            double acc = 0.0, powr = 1.0;
            int rmax = static_cast<int>(n) - 1;
            for (int r = 1; r <= rmax; ++r) {
                powr *= base;
                acc += powr / r;
            }
            second = f1 / n * std::pow(base, -n + 1) * (-std::log(a_cwj) - acc);
        }
        return std::exp(first + second);
    }
    double moment = 0.0;
    double lb_nq = logbinom(n, q);
    for (std::size_t i = 0; i < s.c.size(); ++i)
        if (s.c[i] >= q) moment += s.f[i] * std::exp(logbinom(s.c[i], q) - lb_nq);
    return std::pow(moment, 1.0 / (1 - q));
}

// Sigma_k g(k) * fhat_k(m) over the spectrum, integer m < n. (``_rtd_moment``)
double rtd_moment(const Sample& s, int m, int q) {
    const double n = s.n;
    const double log_cnm = logbinom(n, static_cast<double>(m));
    double total = 0.0;
    for (std::size_t i = 0; i < s.c.size(); ++i) {
        const double cval = s.c[i], fcount = s.f[i];
        int kmax = static_cast<int>(std::min(cval, static_cast<double>(m)));
        for (int k = 1; k <= kmax; ++k) {
            double lt = logbinom(cval, k) + logbinom(n - cval, m - k) - log_cnm;
            double g;
            if (q == 0) {
                g = 1.0;
            } else if (q == 1) {
                double r = static_cast<double>(k) / m;
                g = -r * std::log(r);
            } else {
                g = std::pow(static_cast<double>(k) / m, q);
            }
            total += fcount * g * std::exp(lt);
        }
    }
    return total;
}

// Closed-form order-2 Hill number at size m. (``_d2``)
double d2(const Sample& s, double m) {
    const double n = s.n;
    double sum = 0.0;
    for (double xi : s.x) sum += xi * (xi - 1.0);
    double ss = sum / (n * (n - 1));
    return 1.0 / (1.0 / m + (m - 1.0) / m * ss);
}

// Rarefied (interpolated) Hill number of order q at integer size m. (``_rtd``)
double rtd(const Sample& s, int m, int q) {
    if (q == 0) return rtd_moment(s, m, 0);
    if (q == 1) return std::exp(rtd_moment(s, m, 1));
    if (q == 2) return d2(s, static_cast<double>(m));
    return std::pow(rtd_moment(s, m, q), 1.0 / (1 - q));
}

// Interpolated rtd with a linear floor/ceil blend for non-integer m.
double rtd_at(const Sample& s, double m, int q) {
    if (m == std::round(m)) return rtd(s, static_cast<int>(std::llround(m)), q);
    int lo = static_cast<int>(std::floor(m)), hi = static_cast<int>(std::ceil(m));
    return (hi - m) * rtd(s, lo, q) + (m - lo) * rtd(s, hi, q);
}

// Sample coverage Chat(m). (``_chat``)
double chat(const Sample& s, double m) {
    const double n = s.n, f1 = s.f1, f2 = s.f2;
    double f0 = (f2 == 0) ? (n - 1) / n * f1 * (f1 - 1) / 2.0
                          : (n - 1) / n * f1 * f1 / 2.0 / f2;
    double a_ext = (f1 > 0) ? n * f0 / (n * f0 + f1) : 1.0;
    auto at_int = [&](double k) -> double {
        if (k == n) return 1 - f1 / n * a_ext;
        double acc = 0.0;
        for (std::size_t i = 0; i < s.c.size(); ++i) {
            double xx = s.c[i];
            if (n - xx >= k)
                acc += s.f[i] * (xx / n)
                       * std::exp(std::lgamma(n - xx + 1) - std::lgamma(n - xx - k + 1)
                                  - std::lgamma(n) + std::lgamma(n - k));
        }
        return 1 - acc;
    };
    if (m == n) return 1 - f1 / n * a_ext;
    if (m > n) return 1 - f1 / n * std::pow(a_ext, m - n + 1);
    if (m == std::round(m)) return at_int(std::round(m));
    double lo = std::floor(m), hi = std::ceil(m);
    return (hi - m) * at_int(lo) + (m - lo) * at_int(hi);
}

inline double clamp01(double v) { return std::min(std::max(v, 0.0), 1.0); }

// Fill qD[i][*] for one order across all sizes; extrapolation constants (which
// are m-independent) are computed once per order. Reproduces ``_diversity_at``.
void curve_for_order(const Sample& s, int q, const Vec& sizes, std::vector<double>& out) {
    const double n = s.n;
    out.resize(sizes.size());
    if (q == 2) {
        for (std::size_t j = 0; j < sizes.size(); ++j) {
            double m = sizes[j];
            out[j] = (m == n) ? plugin(s, 2) : d2(s, m);
        }
        return;
    }
    double obs = plugin(s, q);
    bool need_extrap = false;
    for (double m : sizes)
        if (m > n) need_extrap = true;
    double asy = 0.0, rfd = 0.0, beta = 0.0;
    if (need_extrap) {
        asy = asymptotic(s, q);
        rfd = rtd(s, static_cast<int>(n - 1), q);
        beta = (asy == rfd) ? 0.0 : (obs - rfd) / (asy - rfd);
    }
    for (std::size_t j = 0; j < sizes.size(); ++j) {
        double m = sizes[j];
        if (m < n) out[j] = rtd_at(s, m, q);
        else if (m == n) out[j] = obs;
        else out[j] = obs + (asy - obs) * (1 - std::pow(1 - beta, m - n));
    }
}

void fill_curve(const Sample& s, const std::vector<int>& qs, const Vec& sizes,
                Mat& qD, Vec& coverage) {
    qD.assign(qs.size(), Vec(sizes.size()));
    coverage.assign(sizes.size(), 0.0);
    for (std::size_t j = 0; j < sizes.size(); ++j)
        coverage[j] = clamp01(chat(s, sizes[j]));
    for (std::size_t i = 0; i < qs.size(); ++i)
        curve_for_order(s, qs[i], sizes, qD[i]);
}

// Augmented-assemblage detection probabilities (unnormalized). (``_bootstrap_probs``)
Vec bootstrap_probs(const Sample& s) {
    const double n = s.n, f1 = s.f1, f2 = s.f2;
    double f0 = (f2 == 0) ? (n - 1) / n * f1 * (f1 - 1) / 2.0
                          : (n - 1) / n * f1 * f1 / 2.0 / f2;
    double a_ext = (f1 > 0) ? n * f0 / (n * f0 + f1) : 1.0;
    double a = f1 / n * a_ext;
    double b = 0.0;
    for (double xi : s.x) b += xi / n * std::pow(1 - xi / n, n);
    double w = (f0 == 0 || b == 0) ? 0.0 : a / b;
    int k = static_cast<int>(std::ceil(f0));
    Vec p;
    p.reserve(s.x.size() + (k > 0 ? k : 0));
    for (double xi : s.x) p.push_back(xi / n * (1 - w * std::pow(1 - xi / n, n)));
    for (int i = 0; i < k; ++i) p.push_back(a / k);
    return p;
}

// Sequential-conditional multinomial(n, p) draw. (numpy rng.multinomial analog)
void multinomial(std::mt19937_64& rng, int n, const Vec& p, std::vector<int>& out) {
    out.assign(p.size(), 0);
    int remaining = n;
    double remaining_p = 1.0;
    for (std::size_t i = 0; i + 1 < p.size(); ++i) {
        if (remaining <= 0) { out[i] = 0; continue; }
        double pi = (remaining_p > 0) ? p[i] / remaining_p : 0.0;
        pi = std::min(std::max(pi, 0.0), 1.0);
        std::binomial_distribution<int> binom(remaining, pi);
        int cnt = binom(rng);
        out[i] = cnt;
        remaining -= cnt;
        remaining_p -= p[i];
        if (remaining_p < 0) remaining_p = 0;
    }
    if (!p.empty()) out[p.size() - 1] = remaining;
}

// Bootstrap SE of qD(m); shape [n_orders][n_sizes]. (``_bootstrap_se``)
Mat bootstrap_se(const Sample& s, const std::vector<int>& qs, const Vec& sizes,
                 int nboot, std::uint64_t seed) {
    const std::size_t nq = qs.size(), nm = sizes.size();
    Vec p = bootstrap_probs(s);
    double psum = 0.0;
    for (double v : p) psum += v;
    for (double& v : p) v /= psum;
    int n = static_cast<int>(std::llround(s.n));

    // running mean + M2 (Welford) per (q, m)
    Mat mean(nq, Vec(nm, 0.0));
    Mat m2(nq, Vec(nm, 0.0));

    std::mt19937_64 rng(seed);
    std::vector<int> reps;
    Vec xb;
    Mat qD;
    Vec cov;  // discarded
    for (int r = 0; r < nboot; ++r) {
        multinomial(rng, n, p, reps);
        xb.clear();
        for (int v : reps)
            if (v > 0) xb.push_back(static_cast<double>(v));
        Sample sb = make_sample(xb);
        qD.assign(nq, Vec(nm));
        for (std::size_t i = 0; i < nq; ++i)
            curve_for_order(sb, qs[i], sizes, qD[i]);
        double count = r + 1;
        for (std::size_t i = 0; i < nq; ++i)
            for (std::size_t j = 0; j < nm; ++j) {
                double val = qD[i][j];
                double delta = val - mean[i][j];
                mean[i][j] += delta / count;
                m2[i][j] += delta * (val - mean[i][j]);
            }
    }
    Mat se(nq, Vec(nm, 0.0));
    if (nboot > 1) {
        for (std::size_t i = 0; i < nq; ++i)
            for (std::size_t j = 0; j < nm; ++j)
                se[i][j] = std::sqrt(m2[i][j] / (nboot - 1));
    }
    return se;
}

}  // namespace

// --------------------------------------------------------------------------- //
// public API
// --------------------------------------------------------------------------- //
double digamma(double x) {
    double result = 0.0;
    while (x < 6.0) {
        result -= 1.0 / x;
        x += 1.0;
    }
    double r = 1.0 / x;
    result += std::log(x) - 0.5 * r;
    double rr = r * r;
    result -= rr * (1.0 / 12.0
                    - rr * (1.0 / 120.0
                            - rr * (1.0 / 252.0
                                    - rr * (1.0 / 240.0
                                            - rr * (1.0 / 132.0
                                                    - rr * (691.0 / 32760.0))))));
    return result;
}

InextCurve inext_curve(const Vec& counts, const std::vector<int>& qs, const Vec& sizes) {
    Sample s = make_sample(counts);
    InextCurve out;
    fill_curve(s, qs, sizes, out.qD, out.coverage);
    return out;
}

Mat inext_bootstrap(const Vec& counts, const std::vector<int>& qs, const Vec& sizes,
                    int nboot, std::uint64_t seed) {
    Sample s = make_sample(counts);
    return bootstrap_se(s, qs, sizes, nboot, seed);
}

std::vector<InextSample> inext_batch(const std::vector<Vec>& samples,
                                     const std::vector<Vec>& sample_sizes,
                                     const std::vector<int>& qs,
                                     int nboot, std::uint64_t seed, int threads) {
    const int N = static_cast<int>(samples.size());
    std::vector<InextSample> results(N);

    auto worker = [&](int idx) {
        Sample s = make_sample(samples[idx]);
        const Vec& sizes = sample_sizes[idx];
        InextSample& res = results[idx];
        fill_curve(s, qs, sizes, res.qD, res.coverage);
        if (nboot > 0)
            res.se = bootstrap_se(s, qs, sizes, nboot,
                                  seed + static_cast<std::uint64_t>(idx));
    };

    int nthreads = (threads > 0)
                       ? threads
                       : static_cast<int>(std::thread::hardware_concurrency());
    if (nthreads < 1) nthreads = 1;
    if (nthreads > N) nthreads = N;

    if (nthreads <= 1 || N <= 1) {
        for (int i = 0; i < N; ++i) worker(i);
        return results;
    }
    std::atomic<int> next{0};
    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    for (int t = 0; t < nthreads; ++t)
        pool.emplace_back([&]() {
            int i;
            while ((i = next.fetch_add(1)) < N) worker(i);
        });
    for (auto& th : pool) th.join();
    return results;
}

}  // namespace vdjtools
