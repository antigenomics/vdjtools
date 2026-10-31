// Exercise private message algebra without exposing a testing API.
#include "../../src/pgen.cpp"
#include <iostream>

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
