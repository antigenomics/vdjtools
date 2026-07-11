#pragma once
#include <cstdint>
#include <vector>

// Native iNEXT size-based rarefaction/extrapolation kernel.
//
// A faithful C++ port of the validated numpy engine in
// ``python/vdjtools/stats/inext.py`` (the ``TD.m.est`` / ``Chat.Ind`` /
// ``Diversity_profile`` / ``EstiBootComm.Ind`` internals of the iNEXT R package
// v3.0.2). Only the SIZE-BASED machinery is ported: the MVUE interpolation
// kernel, the q=2 closed form, the beta-method extrapolation, the asymptotic
// estimators (Chao1 q=0; Chao-Wang-Jost Shannon q=1; MVUE Simpson q=2), and the
// sample-coverage Chat(m). Coverage-based inversion stays in Python.
//
// Point curves are deterministic and reproduce the Python/oracle numbers
// exactly; the bootstrap draws replicates with a seeded std::mt19937_64 and
// returns per-(q,m) standard deviations.

namespace vdjtools {

using Vec = std::vector<double>;
using Mat = std::vector<std::vector<double>>;  // [n_orders][n_sizes]

// Deterministic point curve for a single assemblage.
struct InextCurve {
    Mat qD;         // qD[i][j] = Hill number of order qs[i] at size sizes[j]
    Vec coverage;   // coverage[j] = clamped Chat(sizes[j])
};

// Point curve + bootstrap standard errors for a single assemblage (batch item).
struct InextSample {
    Mat qD;         // point estimates
    Vec coverage;   // clamped Chat(m)
    Mat se;         // bootstrap SE per (order, size); empty when nboot == 0
};

// Digamma (psi) via asymptotic expansion + recurrence; matches scipy.special.digamma.
double digamma(double x);

// Deterministic size-based R/E point curve (interp + q2 closed form + beta
// extrapolation) plus the sample-coverage curve. ``counts`` are clonotype
// abundances (non-positive entries ignored); ``sizes`` are sampling depths.
InextCurve inext_curve(const Vec& counts, const std::vector<int>& qs, const Vec& sizes);

// Bootstrap standard errors of qD(m): build the iNEXT augmented assemblage,
// draw ``nboot`` seeded multinomial replicates, recompute the curve per
// replicate, return the per-(q,m) standard deviation (ddof=1).
Mat inext_bootstrap(const Vec& counts, const std::vector<int>& qs, const Vec& sizes,
                    int nboot, std::uint64_t seed);

// Point curve + bootstrap SE for many samples, parallelized across samples.
// Each sample carries its own sizes grid; sample ``i`` is seeded ``seed + i``.
// ``threads <= 0`` uses hardware_concurrency(); it is capped at the sample count.
std::vector<InextSample> inext_batch(const std::vector<Vec>& samples,
                                     const std::vector<Vec>& sample_sizes,
                                     const std::vector<int>& qs,
                                     int nboot, std::uint64_t seed, int threads);

}  // namespace vdjtools
