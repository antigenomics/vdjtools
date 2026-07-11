#include "vdjtools/model.hpp"

#include <algorithm>
#include <cmath>
#include <thread>

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

// Sum over the two-D (tandem) scenarios producing `mid` = [insVD] D1 [insDD] D2 [insDJ]. Each D
// contributes >=1 nt (disjoint n_D=1/n_D=2 partition; mirrors the Python reference _dd_middle).
// Factorized per D1 into left/right partial sums so it is O(nD^2 L^2 N + nD N^2), not the naive
// O(nD^2 L^4 N^2): A[e1] = weight of [insVD] D1 ending at position e1; B[s2] = weight of D2
// starting at s2 then [insDJ] to the end (with P(D2|D1)); the DD insertion couples them in one
// O(N^2) sweep that accumulates the DD Markov product incrementally.
double dd_middle(const PackedModel& m, int j, const int8_t* mid, int mlen) {
    double total = 0.0;
    int nD = m.nD();
    std::vector<double> A(mlen + 1), B(mlen + 1);
    for (int d1 : m.func_d) {
        double pd1 = m.pd_given_j[j * nD + d1];
        if (pd1 == 0.0) continue;
        const auto& cut1 = m.cut_d[d1];
        int L1 = cut1.size();
        // A[e1]: [insVD] D1 with D1 (>=1 nt) ending at e1.
        A.assign(mlen + 1, 0.0);
        for (int i5 = 0; i5 <= L1 && i5 < m.nbins_d5; ++i5)
            for (int i3 = 0; i3 <= L1 - i5 && i3 < m.nbins_d3; ++i3) {
                int ld1 = L1 - i5 - i3;
                if (ld1 < 1) continue;
                double pdel1 = m.del_d[(d1 * m.nbins_d5 + i5) * m.nbins_d3 + i3];
                if (pdel1 == 0.0) continue;
                for (int pos1 = 0; pos1 + ld1 <= mlen; ++pos1) {
                    bool ok = true;
                    for (int k = 0; k < ld1; ++k)
                        if (cut1[i5 + k] != mid[pos1 + k]) { ok = false; break; }
                    if (!ok) continue;
                    double left = p_insert(m.ins_vd, m.R_vd, m.bias_vd, mid, 0, pos1, false);
                    if (left != 0.0) A[pos1 + ld1] += pd1 * pdel1 * left;
                }
            }
        // B[s2]: D2 (>=1 nt) starting at s2 then [insDJ] to the end, weighted by P(D2|D1).
        B.assign(mlen + 1, 0.0);
        for (int d2 : m.func_d) {
            double pd2 = m.pd2_given_d1[d1 * nD + d2];
            if (pd2 == 0.0) continue;
            const auto& cut2 = m.cut_d[d2];
            int L2 = cut2.size();
            for (int k5 = 0; k5 <= L2 && k5 < m.nbins_d5; ++k5)
                for (int k3 = 0; k3 <= L2 - k5 && k3 < m.nbins_d3; ++k3) {
                    int ld2 = L2 - k5 - k3;
                    if (ld2 < 1) continue;
                    double pdel2 = m.del_d2[(d2 * m.nbins_d5 + k5) * m.nbins_d3 + k3];
                    if (pdel2 == 0.0) continue;
                    for (int pos2 = 0; pos2 + ld2 <= mlen; ++pos2) {
                        bool ok = true;
                        for (int k = 0; k < ld2; ++k)
                            if (cut2[k5 + k] != mid[pos2 + k]) { ok = false; break; }
                        if (!ok) continue;
                        double right = p_insert(m.ins_dj, m.R_dj, m.bias_dj, mid, pos2 + ld2, mlen, true);
                        if (right != 0.0) B[pos2] += pd2 * pdel2 * right;
                    }
                }
        }
        // Combine: sum_{e1 <= s2} A[e1] * P_insDD(mid[e1:s2]) * B[s2].
        int ddlen = static_cast<int>(m.ins_dd.size());
        for (int e1 = 0; e1 <= mlen; ++e1) {
            if (A[e1] == 0.0) continue;
            double markov = 1.0;  // product of DD transitions over mid[e1+1..s2-1]
            for (int s2 = e1; s2 <= mlen; ++s2) {
                int len = s2 - e1;
                if (len >= ddlen) break;
                double pins = m.ins_dd[len];
                if (pins != 0.0 && B[s2] != 0.0) {
                    double w = (len == 0) ? pins : pins * m.bias_dd[mid[e1]] * markov;
                    total += A[e1] * w * B[s2];
                }
                if (s2 >= e1 + 1) markov *= m.R_dd[mid[s2] * 4 + mid[s2 - 1]];
            }
        }
    }
    return total;
}

// P(mid) mixed over the D-count prior: P(n_D=1)*single-D + P(n_D=2)*tandem.
double vdj_middle(const PackedModel& m, int j, const int8_t* mid, int mlen) {
    double t = m.p_nd1 * d_middle(m, j, mid, mlen);
    if (m.dd && m.p_nd2 > 0.0) t += m.p_nd2 * dd_middle(m, j, mid, mlen);
    return t;
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
                        ? vdj_middle(m, jc.j, mid, midlen)
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

// Per-codon allowed-codon set: allowed[c] is a 64-bit mask over codon indices (a*16+b*4+c) that
// pass at amino-acid position c. A singleton mask = one aa (exact query); a wildcard mask (all
// non-stop codons) = "any amino acid" here, which drives motif/regex and Hamming-ball Pgens
// without enumerating the set. ok_codon tests membership; the mask replaces the CODON[...]==aa char
// compare and costs the same.
inline bool ok_codon(uint64_t allowed_c, int codon_idx) {
    return (allowed_c >> codon_idx) & 1ULL;
}

uint64_t mask_for_aa(char x) {
    uint64_t m = 0;
    for (int t = 0; t < 64; ++t)
        if (CODON[t] == x) m |= (1ULL << t);
    return m;
}

uint64_t mask_wildcard() {  // any of the 20 amino acids (excludes stop) — OLGA's degenerate 'X'
    uint64_t m = 0;
    for (int t = 0; t < 64; ++t)
        if (CODON[t] != '*') m |= (1ULL << t);
    return m;
}

std::vector<int> mask_or_all(int idx, const std::vector<int>& all) {
    if (idx >= 0) return {idx};
    return all;
}

// ---- aa Pgen: Murugan/OLGA transfer-matrix (Pi_L * Pi_R split) -----------------------------
// We factorize: the *only* cross-boundary coupling between the V/VD-insertion side and the
// D-weighted J/DJ-insertion side is codon translation (never the insertion Markov chain), so the
// left DP (Lf) and right DP (Rb) are each built once and stitched at the D placement. State is the
// trailing two nt (nt[i-1], nt[i-2]), encoded s = (a+1)*5 + (b+1) in [0,25). Matches OLGA exactly.
inline int sidx(int a, int b) { return (a + 1) * 5 + (b + 1); }

// Lf[p*25 + s] = weight of V germline (>=1 nt) + N1 insertion filling CDR3 nt positions [0,p),
// trailing state s = (nt[p-1], nt[p-2]); Pins + first-nt bias baked in; complete codons in [0,p)
// constrained. p is where the next segment (D for VDJ, J for VJ) starts. The insertion arrays are
// passed in (VD for VDJ, VJ for VJ). ``jbind`` >= 0 folds the VJ gene weight P(J|V) into the V seed
// (so Lf sums over V at fixed J); jbind < 0 uses P(V) alone (VDJ, where J decouples from V).
void mk_left_tm(const PackedModel& m, const uint64_t* allowed, int N, int v_idx,
                const std::vector<double>& pins, const std::vector<double>& R,
                const std::vector<double>& bias, int jbind,
                std::vector<double>& Lf, std::vector<char>& lf_any) {
    Lf.assign((N + 1) * 25, 0.0);
    std::vector<double> seed((N + 1) * 25, 0.0);
    std::vector<char> has_seed(N + 1, 0);
    for (int v : mask_or_all(v_idx, m.func_v)) {
        double pv = m.pv[v];
        if (pv == 0.0) continue;
        if (jbind >= 0) {  // VJ: weight the V seed by P(J|V)
            double pjv = m.pjv[v * m.nJ() + jbind];
            if (pjv == 0.0) continue;
            pv *= pjv;
        }
        const auto& gv = m.cut_v[v];
        int Lv = gv.size();
        for (int len_v = 1; len_v <= std::min(Lv, N); ++len_v) {
            int div = Lv - len_v;
            if (div < 0 || div >= m.nbins_v) continue;
            double pdv = m.del_v[v * m.nbins_v + div];
            if (pdv == 0.0) continue;
            bool ok = true;  // complete codons within [0,len_v) must translate
            for (int c = 0; c < len_v / 3; ++c) {
                if (!ok_codon(allowed[c], gv[3 * c] * 16 + gv[3 * c + 1] * 4 + gv[3 * c + 2])) { ok = false; break; }
            }
            if (!ok) continue;
            int a = gv[len_v - 1], b = (len_v >= 2) ? gv[len_v - 2] : -1;
            seed[len_v * 25 + sidx(a, b)] += pv * pdv;
            has_seed[len_v] = 1;
        }
    }
    std::vector<double> cur(25), ncur(25);
    for (int q = 0; q <= N; ++q) {
        if (!has_seed[q]) continue;
        for (int s = 0; s < 25; ++s) cur[s] = seed[q * 25 + s];
        if (!pins.empty() && pins[0] != 0.0)
            for (int s = 0; s < 25; ++s) Lf[q * 25 + s] += pins[0] * cur[s];
        for (int ell = 1; q + ell <= N && ell < (int)pins.size(); ++ell) {
            int mm = q + ell - 1;
            bool ce = (mm % 3 == 2);
            uint64_t am = allowed[mm / 3];
            std::fill(ncur.begin(), ncur.end(), 0.0);
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
            cur.swap(ncur);
            if (!any) break;
            if (pins[ell] != 0.0)
                for (int s = 0; s < 25; ++s) Lf[(q + ell) * 25 + s] += pins[ell] * cur[s];
        }
    }
    lf_any.assign(N + 1, 0);
    for (int p = 0; p <= N; ++p)
        for (int s = 0; s < 25; ++s)
            if (Lf[p * 25 + s] != 0.0) { lf_any[p] = 1; break; }
}

// Rb[p*25 + s] = sum_J P(D|J) P(J) P(delJ) over DJ insertion + J germline filling [p,N), trailing
// state s = (nt[p], nt[p+1]); Pins(DJ) + bias baked in; complete codons in [p,N) constrained. The
// DJ Markov reads 3'->5', so the backward pass extends leftward from the J boundary.
void mk_right_tm(const PackedModel& m, const uint64_t* allowed, int N, int D, int j_idx,
                 std::vector<double>& Rb, std::vector<char>& rb_any) {
    Rb.assign((N + 1) * 25, 0.0);
    const auto& pins = m.ins_dj;
    const auto& R = m.R_dj;
    const auto& bias = m.bias_dj;
    int nD = m.nD();
    std::vector<double> seed((N + 1) * 25, 0.0);
    std::vector<char> has_seed(N + 1, 0);
    for (int j : mask_or_all(j_idx, m.func_j)) {
        double pj = m.pj[j];
        if (pj == 0.0) continue;
        double pdg = m.pd_given_j[j * nD + D];
        if (pdg == 0.0) continue;
        const auto& gj = m.cut_j[j];
        int Lj = gj.size();
        for (int len_j = 1; len_j <= std::min(Lj, N); ++len_j) {
            int right = N - len_j, idxj = Lj - len_j;
            if (idxj < 0 || idxj >= m.nbins_j) continue;
            double pdj = m.del_j[j * m.nbins_j + idxj];
            if (pdj == 0.0) continue;
            bool ok = true;  // complete codons within [right,N) must translate
            for (int c = (right + 2) / 3; c < N / 3; ++c) {
                if (3 * c >= right) {
                    int o = 3 * c - right;
                    if (!ok_codon(allowed[c], gj[idxj + o] * 16 + gj[idxj + o + 1] * 4 + gj[idxj + o + 2])) { ok = false; break; }
                }
            }
            if (!ok) continue;
            int cc = gj[idxj], dd = (len_j >= 2) ? gj[idxj + 1] : -1;
            seed[right * 25 + sidx(cc, dd)] += pdg * pj * pdj;
            has_seed[right] = 1;
        }
    }
    std::vector<double> cur(25), ncur(25);
    for (int right = N; right >= 0; --right) {
        if (!has_seed[right]) continue;
        for (int s = 0; s < 25; ++s) cur[s] = seed[right * 25 + s];
        if (!pins.empty() && pins[0] != 0.0)
            for (int s = 0; s < 25; ++s) Rb[right * 25 + s] += pins[0] * cur[s];
        for (int ell = 1; right - ell >= 0 && ell < (int)pins.size(); ++ell) {
            int mm = right - ell;
            bool ce = (mm % 3 == 0);  // codon [mm,mm+1,mm+2] completes now
            uint64_t am = allowed[mm / 3];
            std::fill(ncur.begin(), ncur.end(), 0.0);
            bool any = false;
            for (int s = 0; s < 25; ++s) {
                double w = cur[s];
                if (w == 0.0) continue;
                int c = s / 5 - 1, d = s % 5 - 1;  // nt[mm+1], nt[mm+2]
                for (int nt = 0; nt < 4; ++nt) {
                    double mult = (ell == 1) ? bias[nt] : R[nt * 4 + c];
                    if (mult == 0.0) continue;
                    if (ce && !ok_codon(am, nt * 16 + c * 4 + d)) continue;
                    ncur[sidx(nt, c)] += w * mult;
                    any = true;
                }
            }
            cur.swap(ncur);
            if (!any) break;
            if (pins[ell] != 0.0)
                for (int s = 0; s < 25; ++s) Rb[(right - ell) * 25 + s] += pins[ell] * cur[s];
        }
    }
    rb_any.assign(N + 1, 0);
    for (int p = 0; p <= N; ++p)
        for (int s = 0; s < 25; ++s)
            if (Rb[p * 25 + s] != 0.0) { rb_any[p] = 1; break; }
}

// Stitch left dp (state nt[p-1],nt[p-2]) with right rb (state nt[p],nt[p+1]) at boundary p,
// checking the single codon (if any) that straddles p against its allowed-codon mask.
double combine_tm(const std::vector<double>& dp, const double* rb, int p, const uint64_t* allowed) {
    if (p % 3 == 0) {  // clean cut — no straddling codon
        double a = 0.0, b = 0.0;
        for (int s = 0; s < 25; ++s) { a += dp[s]; b += rb[s]; }
        return a * b;
    }
    if (p % 3 == 1) {  // codon [p-1,p,p+1] = (a,c,d), completes at p+1
        uint64_t am = allowed[(p - 1) / 3];
        double dpA[4] = {0, 0, 0, 0};
        for (int s = 0; s < 25; ++s)
            if (dp[s] != 0.0) { int a = s / 5 - 1; if (a >= 0) dpA[a] += dp[s]; }
        double tot = 0.0;
        for (int s = 0; s < 25; ++s) {
            double w = rb[s];
            if (w == 0.0) continue;
            int c = s / 5 - 1, d = s % 5 - 1;
            if (c < 0 || d < 0) continue;
            for (int a = 0; a < 4; ++a)
                if (dpA[a] != 0.0 && ok_codon(am, a * 16 + c * 4 + d)) tot += dpA[a] * w;
        }
        return tot;
    }
    // p % 3 == 2: codon [p-2,p-1,p] = (b,a,c), completes at p
    uint64_t am = allowed[p / 3];
    double rbC[4] = {0, 0, 0, 0};
    for (int s = 0; s < 25; ++s)
        if (rb[s] != 0.0) { int c = s / 5 - 1; if (c >= 0) rbC[c] += rb[s]; }
    double tot = 0.0;
    for (int s = 0; s < 25; ++s) {
        double w = dp[s];
        if (w == 0.0) continue;
        int a = s / 5 - 1, b = s % 5 - 1;
        if (a < 0 || b < 0) continue;
        for (int c = 0; c < 4; ++c)
            if (rbC[c] != 0.0 && ok_codon(am, b * 16 + a * 4 + c)) tot += w * rbC[c];
    }
    return tot;
}

double pgen_aa_vdj(const PackedModel& m, const uint64_t* allowed, int alen, int v_idx, int j_idx) {
    int N = 3 * alen;
    std::vector<double> Lf, Rb;
    std::vector<char> lf_any, rb_any;
    mk_left_tm(m, allowed, N, v_idx, m.ins_vd, m.R_vd, m.bias_vd, -1, Lf, lf_any);
    double total = 0.0;
    std::vector<double> dp(25), ndp(25);
    for (int D : m.func_d) {
        mk_right_tm(m, allowed, N, D, j_idx, Rb, rb_any);
        const auto& cutd = m.cut_d[D];
        int Ld = cutd.size();
        for (int idx5 = 0; idx5 <= Ld && idx5 < m.nbins_d5; ++idx5) {
            for (int idx3 = 0; idx3 <= Ld - idx5 && idx3 < m.nbins_d3; ++idx3) {
                double pdel = m.del_d[(D * m.nbins_d5 + idx5) * m.nbins_d3 + idx3];
                if (pdel == 0.0) continue;
                int ld = Ld - idx5 - idx3;
                for (int pos = 1; pos <= N - ld; ++pos) {
                    if (!lf_any[pos]) continue;
                    int p = pos + ld;
                    if (p > N || !rb_any[p]) continue;
                    for (int s = 0; s < 25; ++s) dp[s] = Lf[pos * 25 + s];
                    bool ok = true;
                    for (int k = 0; k < ld; ++k) {  // thread the fixed D germline nt
                        int mm = pos + k, nt = cutd[idx5 + k];
                        bool ce = (mm % 3 == 2);
                        uint64_t am = allowed[mm / 3];
                        std::fill(ndp.begin(), ndp.end(), 0.0);
                        bool any = false;
                        for (int s = 0; s < 25; ++s) {
                            double w = dp[s];
                            if (w == 0.0) continue;
                            int p1 = s / 5 - 1, p2 = s % 5 - 1;
                            if (ce && !ok_codon(am, p2 * 16 + p1 * 4 + nt)) continue;
                            ndp[sidx(nt, p1)] += w;
                            any = true;
                        }
                        dp.swap(ndp);
                        if (!any) { ok = false; break; }
                    }
                    if (ok) total += pdel * combine_tm(dp, &Rb[p * 25], p, allowed);
                }
            }
        }
    }
    return total;
}

// VJ aa Pgen (TRA/TRG/IGK/IGL): no D, one VJ insertion. J plays the role D does in the VDJ split.
// For each J we build the left DP (V germline + VJ insertion, weighted by P(V)P(J|V)), then thread
// the fixed J germline suffix over [N-len_j, N) and sum. O(nJ * N * 25) — replaces the old
// per-scenario enumeration, which was ~1000x slower on VJ loci with many V/J alleles.
double pgen_aa_vj(const PackedModel& m, const uint64_t* allowed, int alen, int v_idx, int j_idx) {
    int N = 3 * alen;
    double total = 0.0;
    std::vector<double> Lf, dp(25), ndp(25);
    std::vector<char> lf_any;
    for (int j : mask_or_all(j_idx, m.func_j)) {
        mk_left_tm(m, allowed, N, v_idx, m.ins_vj, m.R_vj, m.bias_vj, j, Lf, lf_any);
        const auto& gj = m.cut_j[j];
        int Lj = gj.size();
        for (int len_j = 1; len_j <= std::min(Lj, N); ++len_j) {
            int p = N - len_j, idxj = Lj - len_j;
            if (p < 1 || !lf_any[p]) continue;  // >=1 V/insertion nt before J
            if (idxj < 0 || idxj >= m.nbins_j) continue;
            double pdj = m.del_j[j * m.nbins_j + idxj];
            if (pdj == 0.0) continue;
            for (int s = 0; s < 25; ++s) dp[s] = Lf[p * 25 + s];
            bool ok = true;
            for (int k = 0; k < len_j; ++k) {  // thread the fixed J germline suffix
                int mm = p + k, nt = gj[idxj + k];
                bool ce = (mm % 3 == 2);
                uint64_t am = allowed[mm / 3];
                std::fill(ndp.begin(), ndp.end(), 0.0);
                bool any = false;
                for (int s = 0; s < 25; ++s) {
                    double w = dp[s];
                    if (w == 0.0) continue;
                    int p1 = s / 5 - 1, p2 = s % 5 - 1;
                    if (ce && !ok_codon(am, p2 * 16 + p1 * 4 + nt)) continue;
                    ndp[sidx(nt, p1)] += w;
                    any = true;
                }
                dp.swap(ndp);
                if (!any) { ok = false; break; }
            }
            if (!ok) continue;
            double sum = 0.0;
            for (int s = 0; s < 25; ++s) sum += dp[s];
            total += pdj * sum;
        }
    }
    return total;
}

}  // namespace

namespace {
double pgen_aa_masked(const PackedModel& m, const uint64_t* allowed, int L, int v_idx, int j_idx) {
    return m.vdj ? pgen_aa_vdj(m, allowed, L, v_idx, j_idx)
                 : pgen_aa_vj(m, allowed, L, v_idx, j_idx);
}
}  // namespace

double pgen_aa(const PackedModel& m, const std::string& aa, int v_idx, int j_idx) {
    int L = aa.size();
    std::vector<uint64_t> allowed(L);
    for (int c = 0; c < L; ++c) allowed[c] = mask_for_aa(aa[c]);
    return pgen_aa_masked(m, allowed.data(), L, v_idx, j_idx);
}

// Total Pgen of the amino-acid CDR3 and every sequence within Hamming distance 1 of it (one aa
// substitution). By inclusion-exclusion this is  sum_k Pgen(a with position k wildcarded)
//   - (L-1) Pgen(a)  (OLGA's identity), but each term here is one fast transfer-matrix pass and the
// wildcard is a single mask (no 19x enumeration). Done entirely in C++ so the packed model is reused.
double pgen_aa_hamming1(const PackedModel& m, const std::string& aa, int v_idx, int j_idx) {
    int L = aa.size();
    if (L == 0) return 0.0;
    std::vector<uint64_t> allowed(L);
    for (int c = 0; c < L; ++c) allowed[c] = mask_for_aa(aa[c]);
    double base = pgen_aa_masked(m, allowed.data(), L, v_idx, j_idx);
    uint64_t wild = mask_wildcard();
    double total = 0.0;
    for (int k = 0; k < L; ++k) {
        uint64_t save = allowed[k];
        allowed[k] = wild;
        total += pgen_aa_masked(m, allowed.data(), L, v_idx, j_idx);
        allowed[k] = save;
    }
    return total - (L - 1) * base;
}

// ---- EM E-step ----------------------------------------------------------------------------
namespace {

void accum_dinucl(std::vector<double>& dn, const int8_t* s, int lo, int hi, bool from_right, double w) {
    int n = hi - lo;
    if (from_right) {
        for (int k = n - 2; k >= 0; --k) dn[s[lo + k] * 4 + s[lo + k + 1]] += w;
    } else {
        for (int k = 1; k < n; ++k) dn[s[lo + k] * 4 + s[lo + k - 1]] += w;
    }
}

// One VJ scenario: weight w, accumulate its realizations into local counts; return w.
double accum_vj(const PackedModel& m, int v, int j, int div, int dij,
                const int8_t* mid, int L, double base, Counts& c) {
    double pins = p_insert(m.ins_vj, m.R_vj, m.bias_vj, mid, 0, L, false);
    double w = base * pins;
    if (w <= 0.0) return 0.0;
    c.v_choice[v] += w;
    c.j_choice[v * m.nJ() + j] += w;
    c.v_3_del[v * m.nbins_v + div] += w;
    c.j_5_del[j * m.nbins_j + dij] += w;
    c.ins_vj[L] += w;
    accum_dinucl(c.dinucl_vj, mid, 0, L, false, w);
    return w;
}

double accum_vdj(const PackedModel& m, int j, int v, int div, int dij,
                 const int8_t* mid, int mlen, double base, const std::vector<int>& dm, Counts& c) {
    int nD = m.nD();
    double seq_total = 0.0;
    for (int d : dm) {
        double pdj = m.pd_given_j[j * nD + d];
        if (pdj == 0.0) continue;
        const auto& cut = m.cut_d[d];
        int L = cut.size();
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
                    int lvd = pos, ldj = mlen - pos - ld;
                    double wvd = p_insert(m.ins_vd, m.R_vd, m.bias_vd, mid, 0, pos, false);
                    if (wvd == 0.0) continue;
                    double wdj = p_insert(m.ins_dj, m.R_dj, m.bias_dj, mid, pos + ld, mlen, true);
                    double w = base * pdj * pdel * wvd * wdj;
                    if (w <= 0.0) continue;
                    seq_total += w;
                    c.v_choice[v] += w;
                    c.j_choice[j] += w;
                    c.d_gene[j * nD + d] += w;
                    c.v_3_del[v * m.nbins_v + div] += w;
                    c.j_5_del[j * m.nbins_j + dij] += w;
                    c.d_del[(d * m.nbins_d5 + idx5) * m.nbins_d3 + idx3] += w;
                    c.ins_vd[lvd] += w;
                    c.ins_dj[ldj] += w;
                    accum_dinucl(c.dinucl_vd, mid, 0, pos, false, w);
                    accum_dinucl(c.dinucl_dj, mid, pos + ld, mlen, true, w);
                }
            }
        }
    }
    return seq_total;
}

// One tandem (n_D=2) block: soft counts for [insVD] D1 [insDD] D2 [insDJ], each D >=1 nt. Factorized
// like `dd_middle` (per D1: left partial sums A[e1], right partial sums B[s2] carrying P(D2|D1)),
// so seq_total is byte-identical to the Pgen. The combine sweep additionally builds the backward
// message C[e1]=sum_{s2} insDD(e1->s2)*B[s2] and forward message Dmsg[s2]=sum_{e1} A[e1]*insDD(e1->s2);
// re-enumerating each block once then attributes every per-realization soft count (matches the pure-
// Python reference _accum_dd exactly). V/J and their trims are constant per (V,J) => the whole DD mass.
// Returns base*seq_total (the n_D=2 Pgen contribution for this V,J,vop,jop); caller folds it into n_d[2].
double accum_dd(const PackedModel& m, int j, int v, int div, int dij,
                const int8_t* mid, int mlen, double base, const std::vector<int>& dm, Counts& c) {
    int nD = m.nD();
    int ddlen = static_cast<int>(m.ins_dd.size());
    std::vector<double> A(mlen + 1), B(mlen + 1), C(mlen + 1), Dmsg(mlen + 1);
    double seq_total = 0.0;
    for (int d1 : dm) {
        double pd1 = m.pd_given_j[j * nD + d1];
        if (pd1 == 0.0) continue;
        const auto& cut1 = m.cut_d[d1];
        int L1 = cut1.size();
        // A[e1]: [insVD] D1 (>=1 nt) ending at e1.
        A.assign(mlen + 1, 0.0);
        for (int i5 = 0; i5 <= L1 && i5 < m.nbins_d5; ++i5)
            for (int i3 = 0; i3 <= L1 - i5 && i3 < m.nbins_d3; ++i3) {
                int ld1 = L1 - i5 - i3;
                if (ld1 < 1) continue;
                double pdel1 = m.del_d[(d1 * m.nbins_d5 + i5) * m.nbins_d3 + i3];
                if (pdel1 == 0.0) continue;
                for (int pos1 = 0; pos1 + ld1 <= mlen; ++pos1) {
                    bool ok = true;
                    for (int k = 0; k < ld1; ++k)
                        if (cut1[i5 + k] != mid[pos1 + k]) { ok = false; break; }
                    if (!ok) continue;
                    double wvd = p_insert(m.ins_vd, m.R_vd, m.bias_vd, mid, 0, pos1, false);
                    if (wvd != 0.0) A[pos1 + ld1] += pd1 * pdel1 * wvd;
                }
            }
        // B[s2]: D2 (>=1 nt) starting at s2 then [insDJ] to the end, weighted by P(D2|D1).
        B.assign(mlen + 1, 0.0);
        for (int d2 : dm) {
            double pd2 = m.pd2_given_d1[d1 * nD + d2];
            if (pd2 == 0.0) continue;
            const auto& cut2 = m.cut_d[d2];
            int L2 = cut2.size();
            for (int k5 = 0; k5 <= L2 && k5 < m.nbins_d5; ++k5)
                for (int k3 = 0; k3 <= L2 - k5 && k3 < m.nbins_d3; ++k3) {
                    int ld2 = L2 - k5 - k3;
                    if (ld2 < 1) continue;
                    double pdel2 = m.del_d2[(d2 * m.nbins_d5 + k5) * m.nbins_d3 + k3];
                    if (pdel2 == 0.0) continue;
                    for (int pos2 = 0; pos2 + ld2 <= mlen; ++pos2) {
                        bool ok = true;
                        for (int k = 0; k < ld2; ++k)
                            if (cut2[k5 + k] != mid[pos2 + k]) { ok = false; break; }
                        if (!ok) continue;
                        double wdj = p_insert(m.ins_dj, m.R_dj, m.bias_dj, mid, pos2 + ld2, mlen, true);
                        if (wdj != 0.0) B[pos2] += pd2 * pdel2 * wdj;
                    }
                }
        }
        // Combine: total_D1, C[e1] (backward), Dmsg[s2] (forward); attribute dd_ins/dd_dinucl here.
        C.assign(mlen + 1, 0.0);
        Dmsg.assign(mlen + 1, 0.0);
        double total_D1 = 0.0;
        for (int e1 = 0; e1 <= mlen; ++e1) {
            if (A[e1] == 0.0) continue;
            double markov = 1.0;  // product of DD transitions over mid[e1+1..s2-1]
            for (int s2 = e1; s2 <= mlen; ++s2) {
                int len = s2 - e1;
                if (len >= ddlen) break;
                double pins = m.ins_dd[len];
                if (pins != 0.0 && B[s2] != 0.0) {
                    double wdd = (len == 0) ? pins : pins * m.bias_dd[mid[e1]] * markov;
                    double comb = A[e1] * wdd * B[s2];
                    total_D1 += comb;
                    C[e1] += wdd * B[s2];
                    Dmsg[s2] += A[e1] * wdd;
                    double bw = base * comb;
                    c.ins_dd[len] += bw;
                    accum_dinucl(c.dinucl_dd, mid, e1, s2, false, bw);
                }
                if (s2 >= e1 + 1) markov *= m.R_dd[mid[s2] * 4 + mid[s2 - 1]];
            }
        }
        if (total_D1 == 0.0) continue;
        c.d_gene[j * nD + d1] += base * total_D1;
        seq_total += total_D1;
        // Re-enumerate LEFT: attribute d_del(D1), vd_ins, vd_dinucl with weight (left realization)*C[e1].
        for (int i5 = 0; i5 <= L1 && i5 < m.nbins_d5; ++i5)
            for (int i3 = 0; i3 <= L1 - i5 && i3 < m.nbins_d3; ++i3) {
                int ld1 = L1 - i5 - i3;
                if (ld1 < 1) continue;
                double pdel1 = m.del_d[(d1 * m.nbins_d5 + i5) * m.nbins_d3 + i3];
                if (pdel1 == 0.0) continue;
                for (int pos1 = 0; pos1 + ld1 <= mlen; ++pos1) {
                    int e1 = pos1 + ld1;
                    if (C[e1] == 0.0) continue;
                    bool ok = true;
                    for (int k = 0; k < ld1; ++k)
                        if (cut1[i5 + k] != mid[pos1 + k]) { ok = false; break; }
                    if (!ok) continue;
                    double wvd = p_insert(m.ins_vd, m.R_vd, m.bias_vd, mid, 0, pos1, false);
                    if (wvd == 0.0) continue;
                    double bw = base * pd1 * pdel1 * wvd * C[e1];
                    c.d_del[(d1 * m.nbins_d5 + i5) * m.nbins_d3 + i3] += bw;
                    c.ins_vd[pos1] += bw;
                    accum_dinucl(c.dinucl_vd, mid, 0, pos1, false, bw);
                }
            }
        // Re-enumerate RIGHT: attribute d2_gene, d2_del(D2), dj_ins, dj_dinucl with weight (right)*Dmsg[s2].
        for (int d2 : dm) {
            double pd2 = m.pd2_given_d1[d1 * nD + d2];
            if (pd2 == 0.0) continue;
            const auto& cut2 = m.cut_d[d2];
            int L2 = cut2.size();
            for (int k5 = 0; k5 <= L2 && k5 < m.nbins_d5; ++k5)
                for (int k3 = 0; k3 <= L2 - k5 && k3 < m.nbins_d3; ++k3) {
                    int ld2 = L2 - k5 - k3;
                    if (ld2 < 1) continue;
                    double pdel2 = m.del_d2[(d2 * m.nbins_d5 + k5) * m.nbins_d3 + k3];
                    if (pdel2 == 0.0) continue;
                    for (int pos2 = 0; pos2 + ld2 <= mlen; ++pos2) {
                        if (Dmsg[pos2] == 0.0) continue;
                        bool ok = true;
                        for (int k = 0; k < ld2; ++k)
                            if (cut2[k5 + k] != mid[pos2 + k]) { ok = false; break; }
                        if (!ok) continue;
                        double wdj = p_insert(m.ins_dj, m.R_dj, m.bias_dj, mid, pos2 + ld2, mlen, true);
                        if (wdj == 0.0) continue;
                        double bw = base * pd2 * pdel2 * wdj * Dmsg[pos2];
                        c.d2_gene[d1 * nD + d2] += bw;
                        c.d2_del[(d2 * m.nbins_d5 + k5) * m.nbins_d3 + k3] += bw;
                        c.ins_dj[mlen - pos2 - ld2] += bw;
                        accum_dinucl(c.dinucl_dj, mid, pos2 + ld2, mlen, true, bw);
                    }
                }
        }
    }
    if (seq_total > 0.0) {
        double bw = base * seq_total;
        c.v_choice[v] += bw;
        c.j_choice[j] += bw;
        c.v_3_del[v * m.nbins_v + div] += bw;
        c.j_5_del[j * m.nbins_j + dij] += bw;
    }
    return base * seq_total;
}

double estep_one(const PackedModel& m, const std::vector<int8_t>& cdr3,
                 const std::vector<int>& vmask, const std::vector<int>& jmask,
                 const std::vector<int>& dmask, Counts& local) {
    int N = cdr3.size();
    const int8_t* s = cdr3.data();
    const std::vector<int>& vm = vmask.empty() ? m.func_v : vmask;
    const std::vector<int>& jm = jmask.empty() ? m.func_j : jmask;
    const std::vector<int>& dm = (m.vdj && !dmask.empty()) ? dmask : m.func_d;

    struct JC {
        int j;
        std::vector<Opt> opts;
    };
    std::vector<JC> jc;
    std::vector<Opt> tmp;
    for (int j : jm) {
        j_options(m, j, s, N, tmp);
        if (!tmp.empty()) jc.push_back({j, tmp});
    }
    double total = 0.0;
    std::vector<Opt> vo;
    for (int v : vm) {
        double pv = m.pv[v];
        if (pv == 0.0) continue;
        v_options(m, v, s, N, vo);
        if (vo.empty()) continue;
        int Lv = m.cut_v[v].size();
        for (const auto& J : jc) {
            double pj = m.vdj ? m.pj[J.j] : m.pjv[v * m.nJ() + J.j];
            if (pj == 0.0) continue;
            int Lj = m.cut_j[J.j].size();
            for (const auto& vop : vo) {
                for (const auto& jop : J.opts) {
                    if (vop.len + jop.len > N) continue;
                    int midlen = N - vop.len - jop.len;
                    const int8_t* mid = s + vop.len;
                    double base = pv * pj * vop.p * jop.p;
                    int div = Lv - vop.len, dij = Lj - jop.len;
                    if (!m.vdj) {
                        total += accum_vj(m, v, J.j, div, dij, mid, midlen, base, local);
                    } else {
                        // n_D=1 (0-D folds in via a fully-trimmed D); weighted by P(n_D<=1).
                        double c1 = accum_vdj(m, J.j, v, div, dij, mid, midlen, base * m.p_nd1, dm, local);
                        local.n_d[1] += c1;
                        total += c1;
                        if (m.dd && m.p_nd2 > 0.0) {  // n_D=2 tandem, weighted by P(n_D=2)
                            double c2 = accum_dd(m, J.j, v, div, dij, mid, midlen, base * m.p_nd2, dm, local);
                            local.n_d[2] += c2;
                            total += c2;
                        }
                    }
                }
            }
        }
    }
    return total;
}

void zero(Counts& c) {
    for (auto* v : {&c.v_choice, &c.j_choice, &c.d_gene, &c.v_3_del, &c.j_5_del, &c.d_del,
                    &c.ins_vd, &c.ins_dj, &c.ins_vj, &c.dinucl_vd, &c.dinucl_dj, &c.dinucl_vj,
                    &c.n_d, &c.d2_gene, &c.d2_del, &c.ins_dd, &c.dinucl_dd}) {
        std::fill(v->begin(), v->end(), 0.0);
    }
}

void add_scaled(Counts& dst, const Counts& src, double s) {
    auto go = [s](std::vector<double>& a, const std::vector<double>& b) {
        for (size_t i = 0; i < a.size(); ++i) a[i] += b[i] * s;
    };
    go(dst.v_choice, src.v_choice); go(dst.j_choice, src.j_choice); go(dst.d_gene, src.d_gene);
    go(dst.v_3_del, src.v_3_del); go(dst.j_5_del, src.j_5_del); go(dst.d_del, src.d_del);
    go(dst.ins_vd, src.ins_vd); go(dst.ins_dj, src.ins_dj); go(dst.ins_vj, src.ins_vj);
    go(dst.dinucl_vd, src.dinucl_vd); go(dst.dinucl_dj, src.dinucl_dj); go(dst.dinucl_vj, src.dinucl_vj);
    go(dst.n_d, src.n_d); go(dst.d2_gene, src.d2_gene); go(dst.d2_del, src.d2_del);
    go(dst.ins_dd, src.ins_dd); go(dst.dinucl_dd, src.dinucl_dd);
}

}  // namespace

Counts make_counts(const PackedModel& m) {
    Counts c;
    int nV = m.nV(), nJ = m.nJ(), nD = m.nD();
    c.v_choice.assign(nV, 0.0);
    c.v_3_del.assign(nV * m.nbins_v, 0.0);
    c.j_5_del.assign(nJ * m.nbins_j, 0.0);
    if (m.vdj) {
        c.j_choice.assign(nJ, 0.0);
        c.d_gene.assign(nJ * nD, 0.0);
        c.d_del.assign(nD * m.nbins_d5 * m.nbins_d3, 0.0);
        c.ins_vd.assign(m.ins_vd.size(), 0.0);
        c.ins_dj.assign(m.ins_dj.size(), 0.0);
        c.dinucl_vd.assign(16, 0.0);
        c.dinucl_dj.assign(16, 0.0);
        c.n_d.assign(3, 0.0);  // buckets for n_D in {0,1,2}; only 1 and 2 ever accumulate
        if (m.dd) {
            c.d2_gene.assign(nD * nD, 0.0);
            c.d2_del.assign(nD * m.nbins_d5 * m.nbins_d3, 0.0);
            c.ins_dd.assign(m.ins_dd.size(), 0.0);
            c.dinucl_dd.assign(16, 0.0);
        }
    } else {
        c.j_choice.assign(nV * nJ, 0.0);
        c.ins_vj.assign(m.ins_vj.size(), 0.0);
        c.dinucl_vj.assign(16, 0.0);
    }
    return c;
}

namespace {
// One thread's slice of the read batch → its private accumulator ``acc`` and summed log-Pgen.
double estep_range(const PackedModel& m,
                   const std::vector<std::vector<int8_t>>& seqs,
                   const std::vector<std::vector<int>>& vmasks,
                   const std::vector<std::vector<int>>& jmasks,
                   const std::vector<std::vector<int>>& dmasks,
                   size_t lo, size_t hi, Counts& acc) {
    Counts local = make_counts(m);
    bool have_masks = !vmasks.empty();
    std::vector<int> empty;
    double ll = 0.0;
    for (size_t i = lo; i < hi; ++i) {
        zero(local);
        const std::vector<int>& vm = have_masks ? vmasks[i] : empty;
        const std::vector<int>& jm = have_masks ? jmasks[i] : empty;
        const std::vector<int>& dm = (have_masks && !dmasks.empty()) ? dmasks[i] : empty;
        double total = estep_one(m, seqs[i], vm, jm, dm, local);
        if (total > 0.0) {
            ll += std::log(total);
            add_scaled(acc, local, 1.0 / total);
        }
    }
    return ll;
}
constexpr size_t kEstepThreadMin = 64;  // batches smaller than this run single-threaded (stay bitwise-exact)
}  // namespace

double estep_batch(const PackedModel& m,
                   const std::vector<std::vector<int8_t>>& seqs,
                   const std::vector<std::vector<int>>& vmasks,
                   const std::vector<std::vector<int>>& jmasks,
                   const std::vector<std::vector<int>>& dmasks,
                   Counts& counts,
                   int nthreads) {
    size_t n = seqs.size();
    int T = nthreads;
    if (T <= 0) {
        unsigned hw = std::thread::hardware_concurrency();
        T = hw > 3 ? static_cast<int>(hw) - 2 : 1;
    }
    if (n < kEstepThreadMin || T <= 1) return estep_range(m, seqs, vmasks, jmasks, dmasks, 0, n, counts);
    if (static_cast<size_t>(T) > n) T = static_cast<int>(n);

    std::vector<Counts> acc(T);
    for (int t = 0; t < T; ++t) acc[t] = make_counts(m);
    std::vector<double> lls(T, 0.0);
    std::vector<std::thread> pool;
    size_t chunk = (n + T - 1) / T;
    for (int t = 0; t < T; ++t) {
        size_t lo = static_cast<size_t>(t) * chunk, hi = std::min(n, lo + chunk);
        if (lo >= hi) break;
        pool.emplace_back([&, t, lo, hi] {
            lls[t] = estep_range(m, seqs, vmasks, jmasks, dmasks, lo, hi, acc[t]);
        });
    }
    for (auto& th : pool) th.join();
    double ll = 0.0;  // reduce in fixed thread order → deterministic for a given nthreads
    for (int t = 0; t < T; ++t) { ll += lls[t]; add_scaled(counts, acc[t], 1.0); }
    return ll;
}

}  // namespace vdjtools
