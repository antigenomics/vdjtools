#include "vdjtools/model.hpp"

#include <algorithm>

namespace vdjtools {
namespace {

struct Opt {
    int len;
    double p;
};

int common_prefix(const std::vector<int8_t>& cut, const int8_t* s, int slen) {
    int n = std::min<int>(cut.size(), slen), i = 0;
    while (i < n && cut[i] == s[i]) ++i;
    return i;
}

int common_suffix(const std::vector<int8_t>& cut, const int8_t* s, int slen) {
    int n = std::min<int>(cut.size(), slen), i = 0;
    int lc = cut.size();
    while (i < n && cut[lc - 1 - i] == s[slen - 1 - i]) ++i;
    return i;
}

// (len_v, P(delV|V)) for 3' trims of V whose germline prefixes s; len_v >= 1 (never fully deleted).
void v_options(const PackedModel& m, int v, const int8_t* s, int slen, std::vector<Opt>& out) {
    out.clear();
    const auto& cut = m.cut_v[v];
    int L = cut.size();
    for (int len_v = common_prefix(cut, s, slen); len_v >= 1; --len_v) {
        int di = L - len_v;  // deletion array index = ndel + maxpal
        if (di < 0 || di >= m.nbins_v) continue;
        double p = m.del_v[v * m.nbins_v + di];
        if (p > 0.0) out.push_back({len_v, p});
    }
}

void j_options(const PackedModel& m, int j, const int8_t* s, int slen, std::vector<Opt>& out) {
    out.clear();
    const auto& cut = m.cut_j[j];
    int L = cut.size();
    for (int len_j = common_suffix(cut, s, slen); len_j >= 1; --len_j) {
        int di = L - len_j;
        if (di < 0 || di >= m.nbins_j) continue;
        double p = m.del_j[j * m.nbins_j + di];
        if (p > 0.0) out.push_back({len_j, p});
    }
}

// Probability of the N-region s[lo, hi): P(len) * first-nt bias * dinucleotide Markov chain.
double p_insert(const std::vector<double>& p_len, const std::vector<double>& R,
                const std::vector<double>& bias, const int8_t* s, int lo, int hi, bool from_right) {
    int n = hi - lo;
    if (n >= static_cast<int>(p_len.size())) return 0.0;
    double p = p_len[n];
    if (p == 0.0 || n == 0) return p;
    if (from_right) {
        p *= bias[s[hi - 1]];
        for (int k = n - 2; k >= 0; --k) p *= R[s[lo + k] * 4 + s[lo + k + 1]];
    } else {
        p *= bias[s[lo]];
        for (int k = 1; k < n; ++k) p *= R[s[lo + k] * 4 + s[lo + k - 1]];
    }
    return p;
}

// Sum over D, its 5'/3' trims and its position of P(D|J) * P(delD|D) * Pins(VD) * Pins(DJ).
double d_middle(const PackedModel& m, int j, const int8_t* mid, int mlen) {
    double total = 0.0;
    int nD = m.nD();
    for (int d : m.func_d) {
        double pdj = m.pd_given_j[j * nD + d];
        if (pdj == 0.0) continue;
        const auto& cut = m.cut_d[d];
        int L = cut.size();
        double acc = 0.0;
        for (int idx5 = 0; idx5 <= L && idx5 < m.nbins_d5; ++idx5) {
            for (int idx3 = 0; idx3 <= L - idx5 && idx3 < m.nbins_d3; ++idx3) {
                double pdel = m.del_d[(d * m.nbins_d5 + idx5) * m.nbins_d3 + idx3];
                if (pdel == 0.0) continue;
                int ld = L - idx5 - idx3;
                for (int pos = 0; pos <= mlen - ld; ++pos) {
                    bool ok = true;
                    for (int k = 0; k < ld; ++k) {
                        if (cut[idx5 + k] != mid[pos + k]) { ok = false; break; }
                    }
                    if (!ok) continue;
                    double w = p_insert(m.ins_vd, m.R_vd, m.bias_vd, mid, 0, pos, false);
                    if (w == 0.0) continue;
                    w *= p_insert(m.ins_dj, m.R_dj, m.bias_dj, mid, pos + ld, mlen, true);
                    acc += pdel * w;
                }
            }
        }
        total += pdj * acc;
    }
    return total;
}

}  // namespace

double pgen_nt(const PackedModel& m, const std::vector<int8_t>& cdr3, int v_idx, int j_idx) {
    int N = cdr3.size();
    const int8_t* s = cdr3.data();

    std::vector<int> vmask = (v_idx >= 0) ? std::vector<int>{v_idx} : m.func_v;
    std::vector<int> jmask = (j_idx >= 0) ? std::vector<int>{j_idx} : m.func_j;

    struct JCand {
        int j;
        std::vector<Opt> opts;
    };
    std::vector<JCand> jcands;
    std::vector<Opt> tmp;
    for (int j : jmask) {
        j_options(m, j, s, N, tmp);
        if (!tmp.empty()) jcands.push_back({j, tmp});
    }

    double total = 0.0;
    std::vector<Opt> vopt;
    for (int v : vmask) {
        double pv = m.pv[v];
        if (pv == 0.0) continue;
        v_options(m, v, s, N, vopt);
        if (vopt.empty()) continue;
        for (const auto& jc : jcands) {
            double pj = m.vdj ? m.pj[jc.j] : m.pjv[v * m.nJ() + jc.j];
            if (pj == 0.0) continue;
            for (const auto& vo : vopt) {
                for (const auto& jo : jc.opts) {
                    if (vo.len + jo.len > N) continue;
                    int midlen = N - vo.len - jo.len;
                    const int8_t* mid = s + vo.len;
                    double inner = m.vdj
                        ? d_middle(m, jc.j, mid, midlen)
                        : p_insert(m.ins_vj, m.R_vj, m.bias_vj, mid, 0, midlen, false);
                    total += pv * pj * vo.p * jo.p * inner;
                }
            }
        }
    }
    return total;
}

}  // namespace vdjtools
