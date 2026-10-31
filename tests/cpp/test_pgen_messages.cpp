// Exercise private message algebra without exposing a testing API.
#include "../../src/pgen.cpp"
#include <iostream>

namespace {
using namespace vdjtools;
void extend_ins_reference(const std::vector<double>& dp0, int q0, int N, const uint64_t* allowed,
                     const std::vector<double>& pins, const std::vector<double>& R,
                     const std::vector<double>& bias, std::vector<double>& Mf) {
    // At an insertion END the next block is fixed germline. Completed
    // codons need no past bases; a one-base partial codon needs only that
    // base. Keep the full state in cur for the insertion Markov recurrence.
    auto accumulate = [&](int q, const double* state, double weight) {
        if (q % 3 == 0) {
            double mass = 0.0;
            for (int s = 0; s < 25; ++s) mass += state[s];
            Mf[q * 25] += weight * mass;
        } else if (q % 3 == 1) {
            for (int s = 0; s < 25; ++s) {
                int a = s / 5 - 1;
                if (a >= 0) Mf[q * 25 + sidx(a, -1)] += weight * state[s];
            }
        } else {
            for (int s = 0; s < 25; ++s) Mf[q * 25 + s] += weight * state[s];
        }
    };
    if (!pins.empty() && pins[0] != 0.0) accumulate(q0, dp0.data(), pins[0]);
    std::array<double, 25> first{}, second{};
    double* cur = first.data();
    double* ncur = second.data();
    std::copy_n(dp0.begin(), 25, cur);
    for (int ell = 1; q0 + ell <= N && ell < (int)pins.size(); ++ell) {
        int mm = q0 + ell - 1;
        bool ce = (mm % 3 == 2);
        uint64_t am = allowed[mm / 3];
        std::fill_n(ncur, 25, 0.0);
        bool any = false;
        for (int s = 0; s < 25; ++s) {
            double w = cur[s];
            if (w == 0.0) continue;
            int p1 = s / 5 - 1, p2 = s % 5 - 1;
            for (int nt = 0; nt < 4; ++nt) {
                double mult = (ell == 1) ? bias[nt] : R[nt * 4 + p1];
                if (mult == 0.0) continue;
                if (ce && !ok_codon(am, p2 * 16 + p1 * 4 + nt)) continue;
                ncur[sidx(nt, p1)] += w * mult;
                any = true;
            }
        }
        std::swap(cur, ncur);
        if (!any) break;
        if (pins[ell] != 0.0) accumulate(q0 + ell, cur, pins[ell]);
    }
}

}

int main() {
    using namespace vdjtools;
    constexpr int N = 12;
    const uint64_t masks[] = {mask_for_aa('A'), mask_for_aa('R'), mask_wildcard(), mask_for_aa('G')};
    std::vector<double> right((N + 1) * 25), left(25);
    for (int q = 1; q < N; ++q)
        for (int s = 0; s < 25; ++s)
            if ((q + s) % 3 == 0) right[q * 25 + s] = (s + 1) / 31.0;
    const auto adjoint = right_coefficients(right, N, masks);
    for (int q = 1; q < N; ++q) {
        for (int s = 0; s < 25; ++s)
            left[s] = (s + q) % 4 == 0 ? (s + 2) / 29.0 : 0.0;
        const double forward = combine_tm(left.data(), right.data() + q * 25, q, masks);
        const double reverse = contract_right(left.data(), adjoint.data() + q * 25, -1);
        if (std::abs(forward - reverse) > 1e-13 * std::max(forward, reverse)) {
            std::cerr << "boundary adjoint failed at " << q << '\n';
            return 1;
        }
    }
    // Every endpoint/state must agree with the full 25-state recurrence.
    const std::vector<double> bias = {.1, .2, .3, .4};
    std::vector<double> transitions(16);
    for (int a = 0; a < 4; ++a)
        for (int c = 0; c < 4; ++c) transitions[c * 4 + a] = (c + a) % 4 == 0 ? 0.0 : (1 + (c + a) % 4) / 9.0;
    for (int q0 = 1; q0 <= 6; ++q0) for (int sparse = 0; sparse < 2; ++sparse)
        for (int support = 0; support < 4; ++support) {
            std::vector<double> seed(25), pins(10);
            for (int a = -1; a < 4; ++a) for (int b = -1; b < 4; ++b) {
                if ((q0 % 3 >= 1 && a < 0) || (q0 % 3 == 2 && b < 0)) continue;
                int s = sidx(a, b);
                if (!sparse || s % 3 == 0) seed[s] = (s + 1) / 31.0;
            }
            for (int n = 0; n < 10; ++n)
                if (support == 0 || n == support - 1 || n == 8) pins[n] = (n + 1) / 55.0;
            std::vector<double> actual((N + 1) * 25), expected(actual.size());
            extend_ins_reference(seed, q0, N, masks, pins, transitions, bias, expected);
            extend_ins_into(seed, q0, N, masks, pins, transitions, bias, actual);
            for (size_t i = 0; i < actual.size(); ++i)
                if (std::abs(actual[i] - expected[i]) > 1e-13 * std::max(actual[i], expected[i])) {
                    std::cerr << "insertion endpoint failed: " << q0 << ' ' << i << '\n';
                    return 3;
                }
        }
    PackedModel history;
    history.vdj = history.dd = history.joint_history = true;
    history.p_nd2 = 1.; history.cut_d.resize(3); history.cut_j.resize(1);
    history.func_d = {0, 1, 2}; history.func_j = {0}; history.pj = {1.};
    history.q_tandem = {0., .2, .3, .1, 0., .05, .25, .1, 0.};
    const auto original_history = history.q_tandem;
    const auto factors = rank_one_history(history);
    double approximate_mass = 0.0;
    for (int d = 0; d < 3; ++d) for (int e = 0; e < 3; ++e)
        if (d != e) approximate_mass += factors.u[d] * factors.v[e];
    if (std::abs(approximate_mass - 1.) > 1e-14 || history.q_tandem != original_history) return 4;
    const std::vector<double> cancellation_weights = {1., 1., 1.};
    std::vector<std::vector<double>> middle = {{1e300}, {1e-200}, {0.}};
    if (rank_one_residual(1e300, 1e300, cancellation_weights.data(), 0, 0, history.func_d, middle) != 1e-200) return 5;
    middle[1][0] = 0.;
    if (rank_one_residual(1e300, 1e300, cancellation_weights.data(), 0, 0, history.func_d, middle) != 0.) return 6;
    // Reusing preparation after a shorter/zero input must clear all active cells.
    PrefixPreparation prepared;
    std::vector<double> seeds((N + 1) * 25, 0.25);
    prepared.prepare(seeds, N - 1, masks);
    std::fill(seeds.begin(), seeds.end(), 0.0);
    prepared.prepare(seeds, 4, masks);
    if (std::any_of(prepared.first.begin(), prepared.first.end(), [](double x) { return x != 0.0; }) ||
        std::any_of(prepared.second.begin(), prepared.second.end(), [](double x) { return x != 0.0; }) ||
        std::any_of(prepared.first_any.begin(), prepared.first_any.end(), [](char x) { return x != 0; }))
        return 2;
    return 0;
}
