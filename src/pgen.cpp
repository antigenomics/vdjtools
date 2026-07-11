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

namespace {

// Standard genetic code, indexed [a*16 + b*4 + c] with a,b,c in A,C,G,T = 0..3.
static const char CODON[65] =
    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

// Per-position insertion spec for the aa DP: kind -1 germline, 0 = VD/VJ (5'->3'), 1 = DJ (3'->5').
struct Spec {
    int kind = -1;
    const double* R = nullptr;
    const double* bias = nullptr;
    bool first = false;
    bool last = false;
};

// Codon-constrained left-to-right sum over N-region nt consistent with the aa string, Markov-
// weighted; state = (nt[i-1], nt[i-2]) with -1 sentinel encoded as index (p1+1)*5 + (p2+1).
double aa_dp(const char* aa, const std::vector<int8_t>& tmpl, const std::vector<Spec>& specs) {
    int N = tmpl.size();
    double dp[25] = {0.0};
    dp[0] = 1.0;  // (p1=-1, p2=-1)
    for (int i = 0; i < N; ++i) {
        int fixed = tmpl[i];
        const Spec& sp = specs[i];
        bool ce = (i % 3 == 2);
        char aai = aa[i / 3];
        double ndp[25] = {0.0};
        bool any = false;
        for (int ps = 0; ps < 25; ++ps) {
            double w = dp[ps];
            if (w == 0.0) continue;
            int p1 = ps / 5 - 1, p2 = ps % 5 - 1;
            int lo = (fixed >= 0) ? fixed : 0, hi = (fixed >= 0) ? fixed : 3;
            for (int nt = lo; nt <= hi; ++nt) {
                double ww = w;
                if (sp.kind == 0) {
                    ww *= sp.first ? sp.bias[nt] : sp.R[nt * 4 + p1];
                } else if (sp.kind == 1) {
                    if (!sp.first) ww *= sp.R[p1 * 4 + nt];
                    if (sp.last) ww *= sp.bias[nt];
                }
                if (ce && CODON[p2 * 16 + p1 * 4 + nt] != aai) continue;
                ndp[(nt + 1) * 5 + (p1 + 1)] += ww;
                any = true;
            }
        }
        if (!any) return 0.0;
        for (int k = 0; k < 25; ++k) dp[k] = ndp[k];
    }
    double s = 0.0;
    for (int k = 0; k < 25; ++k) s += dp[k];
    return s;
}

std::vector<int> mask_or_all(int idx, const std::vector<int>& all) {
    if (idx >= 0) return {idx};
    return all;
}

double pgen_aa_vj(const PackedModel& m, const char* aa, int alen, int v_idx, int j_idx) {
    int N = 3 * alen;
    double total = 0.0;
    for (int v : mask_or_all(v_idx, m.func_v)) {
        double pv = m.pv[v];
        if (pv == 0.0) continue;
        const auto& cutv = m.cut_v[v];
        int Lv = cutv.size();
        for (int j : mask_or_all(j_idx, m.func_j)) {
            double pj = m.pjv[v * m.nJ() + j];
            if (pj == 0.0) continue;
            const auto& cutj = m.cut_j[j];
            int Lj = cutj.size();
            for (int len_v = 1; len_v <= std::min(Lv, N); ++len_v) {
                int div = Lv - len_v;
                if (div < 0 || div >= m.nbins_v) continue;
                double pdv = m.del_v[v * m.nbins_v + div];
                if (pdv == 0.0) continue;
                for (int len_j = 1; len_j <= std::min(Lj, N - len_v); ++len_j) {
                    int dij = Lj - len_j;
                    if (dij < 0 || dij >= m.nbins_j) continue;
                    double pdj = m.del_j[j * m.nbins_j + dij];
                    if (pdj == 0.0) continue;
                    int ins_len = N - len_v - len_j;
                    if (ins_len < 0 || ins_len >= (int)m.ins_vj.size() || m.ins_vj[ins_len] == 0.0) continue;
                    std::vector<int8_t> tmpl(N);
                    std::vector<Spec> specs(N);
                    for (int k = 0; k < len_v; ++k) tmpl[k] = cutv[k];
                    for (int k = 0; k < ins_len; ++k) {
                        int p = len_v + k;
                        tmpl[p] = -1;
                        specs[p] = {0, m.R_vj.data(), m.bias_vj.data(), k == 0, k == ins_len - 1};
                    }
                    for (int k = 0; k < len_j; ++k) tmpl[len_v + ins_len + k] = cutj[(Lj - len_j) + k];
                    double w = aa_dp(aa, tmpl, specs);
                    if (w > 0.0) total += pv * pj * pdv * pdj * m.ins_vj[ins_len] * w;
                }
            }
        }
    }
    return total;
}

double d_aa_middle(const PackedModel& m, const char* aa, const std::vector<int8_t>& gv,
                   const std::vector<int8_t>& gj, int len_v, int len_j, int N, int j) {
    int right = N - len_j;
    int nD = m.nD();
    double out = 0.0;
    for (int d : m.func_d) {
        double pdg = m.pd_given_j[j * nD + d];
        if (pdg == 0.0) continue;
        const auto& cutd = m.cut_d[d];
        int L = cutd.size();
        for (int idx5 = 0; idx5 <= L && idx5 < m.nbins_d5; ++idx5) {
            for (int idx3 = 0; idx3 <= L - idx5 && idx3 < m.nbins_d3; ++idx3) {
                double pdel = m.del_d[(d * m.nbins_d5 + idx5) * m.nbins_d3 + idx3];
                if (pdel == 0.0) continue;
                int ld = L - idx5 - idx3;
                for (int pos = len_v; pos <= right - ld; ++pos) {
                    int lvd = pos - len_v, ldj = right - pos - ld;
                    if (lvd >= (int)m.ins_vd.size() || m.ins_vd[lvd] == 0.0) continue;
                    if (ldj >= (int)m.ins_dj.size() || m.ins_dj[ldj] == 0.0) continue;
                    bool dok = true;  // D-germline full codons must translate (cheap pre-prune)
                    for (int c = (pos + 2) / 3; 3 * c + 2 < pos + ld; ++c) {
                        int o = 3 * c - pos;
                        if (CODON[cutd[idx5 + o] * 16 + cutd[idx5 + o + 1] * 4 + cutd[idx5 + o + 2]] != aa[c]) {
                            dok = false;
                            break;
                        }
                    }
                    if (!dok) continue;
                    std::vector<int8_t> tmpl(N);
                    std::vector<Spec> specs(N);
                    for (int k = 0; k < len_v; ++k) tmpl[k] = gv[k];
                    for (int k = 0; k < lvd; ++k) {
                        int p = len_v + k;
                        tmpl[p] = -1;
                        specs[p] = {0, m.R_vd.data(), m.bias_vd.data(), k == 0, k == lvd - 1};
                    }
                    for (int k = 0; k < ld; ++k) tmpl[pos + k] = cutd[idx5 + k];
                    for (int k = 0; k < ldj; ++k) {
                        int p = pos + ld + k;
                        tmpl[p] = -1;
                        specs[p] = {1, m.R_dj.data(), m.bias_dj.data(), k == 0, k == ldj - 1};
                    }
                    for (int k = 0; k < len_j; ++k) tmpl[right + k] = gj[k];
                    double w = aa_dp(aa, tmpl, specs);
                    if (w > 0.0) out += pdg * pdel * m.ins_vd[lvd] * m.ins_dj[ldj] * w;
                }
            }
        }
    }
    return out;
}

double pgen_aa_vdj(const PackedModel& m, const char* aa, int alen, int v_idx, int j_idx) {
    int N = 3 * alen;
    double total = 0.0;
    for (int v : mask_or_all(v_idx, m.func_v)) {
        double pv = m.pv[v];
        if (pv == 0.0) continue;
        const auto& cutv = m.cut_v[v];
        int Lv = cutv.size();
        for (int j : mask_or_all(j_idx, m.func_j)) {
            double pj = m.pj[j];
            if (pj == 0.0) continue;
            const auto& cutj = m.cut_j[j];
            int Lj = cutj.size();
            for (int len_v = 1; len_v <= std::min(Lv, N); ++len_v) {
                int div = Lv - len_v;
                if (div < 0 || div >= m.nbins_v) continue;
                double pdv = m.del_v[v * m.nbins_v + div];
                if (pdv == 0.0) continue;
                bool vok = true;  // V-germline prefix must translate
                for (int c = 0; c < len_v / 3; ++c) {
                    if (CODON[cutv[3 * c] * 16 + cutv[3 * c + 1] * 4 + cutv[3 * c + 2]] != aa[c]) {
                        vok = false;
                        break;
                    }
                }
                if (!vok) continue;
                std::vector<int8_t> gv(cutv.begin(), cutv.begin() + len_v);
                for (int len_j = 1; len_j <= std::min(Lj, N - len_v); ++len_j) {
                    int dij = Lj - len_j;
                    if (dij < 0 || dij >= m.nbins_j) continue;
                    double pdj = m.del_j[j * m.nbins_j + dij];
                    if (pdj == 0.0) continue;
                    std::vector<int8_t> gj(cutj.begin() + (Lj - len_j), cutj.end());
                    total += pv * pj * pdv * pdj * d_aa_middle(m, aa, gv, gj, len_v, len_j, N, j);
                }
            }
        }
    }
    return total;
}

}  // namespace

double pgen_aa(const PackedModel& m, const std::string& aa, int v_idx, int j_idx) {
    return m.vdj ? pgen_aa_vdj(m, aa.c_str(), aa.size(), v_idx, j_idx)
                 : pgen_aa_vj(m, aa.c_str(), aa.size(), v_idx, j_idx);
}

}  // namespace vdjtools
