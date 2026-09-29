#include "vdjtools/model.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cctype>
#include <stdexcept>
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

// NOTE: the aa DP's prefix-sharing trick does not apply here and was measured not to be the cost.
// This inner match loop exits on the first mismatching nucleotide, so a 3' trim that fails costs
// one comparison, not len(D). What makes a V/J-marginalized nt Pgen expensive (36.9 ms/junction on
// TRB, 40.7 on TRD, measured 2026-09-25) is the V x J enumeration wrapped around this call, not
// the D state space inside it.
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

    // Fast path: an in-frame VDJ nt CDR3 reduces to an aa query with singleton codon masks, so the
    // Pi_L*Pi_R transfer matrix gives the identical value ~50x faster than the enumeration below.
    // `pgen_aa_masked` mixes the D-count prior itself (p_nd1*single-D + p_nd2*tandem), so this covers
    // single-D *and* D-D. The enumeration below only runs for non-in-frame nt (N%3!=0) and VJ loci
    // (the latter are already fast via the direct insertion product).
    if (m.vdj && N % 3 == 0) {
        int Lc = N / 3;
        std::vector<uint64_t> allowed(Lc);
        for (int c = 0; c < Lc; ++c) {
            int a = s[3 * c], b = s[3 * c + 1], cc = s[3 * c + 2];
            allowed[c] = (a >= 0 && a < 4 && b >= 0 && b < 4 && cc >= 0 && cc < 4)
                             ? (1ULL << (a * 16 + b * 4 + cc))
                             : 0ULL;  // non-ACGT codon: unsatisfiable, contributes 0 (matches enum)
        }
        return pgen_aa_masked(m, allowed.data(), Lc, v_idx, j_idx);
    }

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
        double pdg = (D < 0) ? 1.0 : m.pd_given_j[j * nD + D];  // D<0: no D weight (tandem right DP)
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

// One fixed germline nucleotide ``nt`` at CDR3 position ``mm``, through the 25-state DP. False
// when the codon masks leave no mass: every longer prefix is then zero too, so callers break.
inline bool thread_one(std::vector<double>& dp, std::vector<double>& ndp, int nt, int mm,
                       const uint64_t* allowed) {
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
    return any;
}

double pgen_aa_vdj(const PackedModel& m, const uint64_t* allowed, int alen, int v_idx, int j_idx) {
    int N = 3 * alen;
    std::vector<double> Lf, Rb, res;
    std::vector<char> lf_any, rb_any, live;
    mk_left_tm(m, allowed, N, v_idx, m.ins_vd, m.R_vd, m.bias_vd, -1, Lf, lf_any);
    double total = 0.0;
    std::vector<double> dp(25), ndp(25);
    for (int D : m.func_d) {
        mk_right_tm(m, allowed, N, D, j_idx, Rb, rb_any);
        const auto& cutd = m.cut_d[D];
        int Ld = cutd.size();
        const double* del_row = &m.del_d[D * m.nbins_d5 * m.nbins_d3];
        for (int idx5 = 0; idx5 <= Ld && idx5 < m.nbins_d5; ++idx5) {
            // The D emitted by trim (idx5, idx3) is a PREFIX of the one emitted by (idx5, idx3-1),
            // so one forward pass covers every 3' trim of this (D, idx5) and the len(D) factor
            // leaves the cost. The smallest 3' trim carrying mass fixes how far to walk.
            int hi3 = std::min(Ld - idx5, m.nbins_d3 - 1), ld_hi = -1;
            for (int idx3 = 0; idx3 <= hi3; ++idx3)
                if (del_row[idx5 * m.nbins_d3 + idx3] != 0.0) { ld_hi = Ld - idx5 - idx3; break; }
            if (ld_hi < 0) continue;  // no live trim for this 5' cut
            live.assign(ld_hi + 1, 0);
            for (int idx3 = 0; idx3 <= hi3; ++idx3)
                if (del_row[idx5 * m.nbins_d3 + idx3] != 0.0) live[Ld - idx5 - idx3] = 1;
            res.assign((ld_hi + 1) * (N + 1), 0.0);
            for (int pos = 1; pos <= N; ++pos) {
                if (!lf_any[pos]) continue;
                for (int s = 0; s < 25; ++s) dp[s] = Lf[pos * 25 + s];
                for (int ld = 0; ld <= ld_hi && pos + ld <= N; ++ld) {
                    if (ld > 0 && !thread_one(dp, ndp, cutd[idx5 + ld - 1], pos + ld - 1, allowed))
                        break;  // mass gone; every longer prefix stays 0 and res is pre-zeroed
                    int p = pos + ld;
                    if (live[ld] && rb_any[p])
                        res[ld * (N + 1) + pos] = combine_tm(dp, &Rb[p * 25], p, allowed);
                }
            }
            // Second pass, in the ORIGINAL (idx3 ascending, pos ascending) order. That is what
            // keeps `total` bitwise identical rather than merely equal to 1e-16 -- summing
            // straight out of the walk reassociates, and vsig:pgen:*:frac_atypical is compared
            // against a frozen pgen_q05 reference where a last-bit move is a data bug.
            for (int idx3 = 0; idx3 <= hi3; ++idx3) {
                double pdel = del_row[idx5 * m.nbins_d3 + idx3];
                if (pdel == 0.0) continue;
                int ld = Ld - idx5 - idx3;
                const double* r = &res[ld * (N + 1)];
                for (int pos = 1; pos <= N - ld; ++pos) total += pdel * r[pos];
            }
        }
    }
    return total;
}

// NOTE: `combine_tm` is the bottleneck now that the D walk shares its prefixes -- unchanged at
// ~437k calls per 20-aa IGH junction against 588k threading steps. The available win is hoisting
// its per-boundary reduction (rb -> 4 or 16 numbers, loop-invariant in dp) out to one pass per
// (D, p), worth roughly another 2x. It is NOT free: that regroups the sum, so Pgen moves in the
// last ~1e-16 and the bitwise gate on `tests/python/fixtures/pgen_golden.json` has to be relaxed
// to a tolerance. Its own piece of work, with its own concordance run.
// Thread a fixed germline block g[0..ld) starting at CDR3 position `start` forward through the
// codon-constrained state DP (25 states). Returns false if the codon masks kill all mass.
inline bool thread_fixed(std::vector<double>& dp, std::vector<double>& ndp,
                         const int8_t* g, int ld, int start, const uint64_t* allowed) {
    for (int k = 0; k < ld; ++k)
        if (!thread_one(dp, ndp, g[k], start + k, allowed)) return false;
    return true;
}

// From a seed state `dp0` at boundary q0, extend forward through an insertion block (pins/R/bias,
// read 5'->3' like VD), accumulating Pins·bias·markov (codon-constrained) into Mf[q*25+s] for
// q in [q0, N]. Same recurrence as mk_left_tm's insertion loop but seeded from an arbitrary state.
void extend_ins_into(const std::vector<double>& dp0, int q0, int N, const uint64_t* allowed,
                     const std::vector<double>& pins, const std::vector<double>& R,
                     const std::vector<double>& bias, std::vector<double>& Mf) {
    if (!pins.empty() && pins[0] != 0.0)
        for (int s = 0; s < 25; ++s) Mf[q0 * 25 + s] += pins[0] * dp0[s];
    std::vector<double> cur(dp0), ncur(25);
    for (int ell = 1; q0 + ell <= N && ell < (int)pins.size(); ++ell) {
        int mm = q0 + ell - 1;
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
            for (int s = 0; s < 25; ++s) Mf[(q0 + ell) * 25 + s] += pins[ell] * cur[s];
    }
}

// NOTE: both D loops below carry the same nested-prefix redundancy `pgen_aa_vdj` just shed, and
// the D2 loop would take the same two-pass treatment. Left alone on purpose: every bundled model
// reports P(n_D=2) = 0, so nothing shipped reaches this function. The D1 loop is also harder --
// it accumulates into the shared `Mf` through `extend_ins_into`, so preserving the summation
// order (and with it bitwise equality) is a different problem from the single-D case.
// Tandem (n_D=2) aa Pgen middle: V+insVD (Lf, reused) → D1 → insDD (Mf) → D2 → insDJ+J (Rb per J).
// P(D1|J) couples D1 to J, so J is looped explicitly (few J on the D-bearing loci) and its P(J),
// P(delJ), insDJ+J go into a D-less right DP; P(D1|J)·P(D2|D1) apply at the D placements. Matches
// the pure-Python `_dd_aa_middle` reference (which matches Σ nt over synonymous codons) exactly.
double pgen_aa_vdj_dd(const PackedModel& m, const uint64_t* allowed, int alen, int v_idx, int j_idx) {
    int N = 3 * alen, nD = m.nD();
    std::vector<double> Lf, RbJ, Mf, dp(25), ndp(25);
    std::vector<char> lf_any, rb_any;
    mk_left_tm(m, allowed, N, v_idx, m.ins_vd, m.R_vd, m.bias_vd, -1, Lf, lf_any);  // V + insVD, once
    int8_t g[128];
    double total = 0.0;
    for (int j : mask_or_all(j_idx, m.func_j)) {
        if (m.pj[j] == 0.0) continue;
        mk_right_tm(m, allowed, N, /*D=*/-1, /*j_idx=*/j, RbJ, rb_any);  // insDJ + J, P(J)·P(delJ), no D
        for (int d1 : m.func_d) {
            double pd1 = m.pd_given_j[j * nD + d1];
            if (pd1 == 0.0) continue;
            Mf.assign((N + 1) * 25, 0.0);
            const auto& cut1 = m.cut_d[d1];
            int L1 = cut1.size();
            for (int i5 = 0; i5 <= L1 && i5 < m.nbins_d5; ++i5)
                for (int i3 = 0; i3 <= L1 - i5 && i3 < m.nbins_d3; ++i3) {
                    int ld1 = L1 - i5 - i3;
                    if (ld1 < 1) continue;
                    double pdel1 = m.del_d[(d1 * m.nbins_d5 + i5) * m.nbins_d3 + i3];
                    if (pdel1 == 0.0) continue;
                    for (int k = 0; k < ld1; ++k) g[k] = cut1[i5 + k];
                    for (int pos1 = 1; pos1 + ld1 < N; ++pos1) {  // leave ≥1 nt for D2 after
                        if (!lf_any[pos1]) continue;
                        double wgt = pd1 * pdel1;
                        for (int s = 0; s < 25; ++s) dp[s] = Lf[pos1 * 25 + s] * wgt;
                        if (!thread_fixed(dp, ndp, g, ld1, pos1, allowed)) continue;
                        extend_ins_into(dp, pos1 + ld1, N, allowed, m.ins_dd, m.R_dd, m.bias_dd, Mf);
                    }
                }
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
                        for (int k = 0; k < ld2; ++k) g[k] = cut2[k5 + k];
                        for (int pos2 = 1; pos2 + ld2 <= N; ++pos2) {
                            int p2 = pos2 + ld2;
                            if (!rb_any[p2]) continue;
                            bool anyq = false;
                            for (int s = 0; s < 25; ++s) if (Mf[pos2 * 25 + s] != 0.0) { anyq = true; break; }
                            if (!anyq) continue;
                            double wgt = pd2 * pdel2;
                            for (int s = 0; s < 25; ++s) dp[s] = Mf[pos2 * 25 + s] * wgt;
                            if (!thread_fixed(dp, ndp, g, ld2, pos2, allowed)) continue;
                            total += combine_tm(dp, &RbJ[p2 * 25], p2, allowed);
                        }
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

double pgen_aa_masked(const PackedModel& m, const uint64_t* allowed, int L, int v_idx, int j_idx) {
    if (!m.vdj) return pgen_aa_vj(m, allowed, L, v_idx, j_idx);
    // n_D mixture: P(n_D≤1)·single-D + P(n_D=2)·tandem. Non-D-D models have p_nd1=1, p_nd2=0.
    double r = m.p_nd1 * pgen_aa_vdj(m, allowed, L, v_idx, j_idx);
    if (m.dd && m.p_nd2 > 0.0) r += m.p_nd2 * pgen_aa_vdj_dd(m, allowed, L, v_idx, j_idx);
    return r;
}

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

namespace {

// out[i] = one(i) for i in [0,n), split across `nthreads` workers (<=0 → auto = hw-2). The work is
// embarrassingly disjoint (no shared state, no reduction), so the result is bitwise-identical to
// the serial loop regardless of thread count; `one` must not throw (an escaping exception in a
// worker terminates the process — validate before calling).
template <class F>
auto run_batch(size_t n, int nthreads, F&& one) -> std::vector<decltype(one(size_t{0}))> {
    std::vector<decltype(one(size_t{0}))> out(n);
    int T = nthreads;
    if (T <= 0) {
        unsigned hw = std::thread::hardware_concurrency();
        T = hw > 3 ? static_cast<int>(hw) - 2 : 1;
    }
    constexpr size_t kBatchThreadMin = 64;  // small batches stay single-threaded (bitwise-exact)
    if (n < kBatchThreadMin || T <= 1) {
        for (size_t i = 0; i < n; ++i) out[i] = one(i);
        return out;
    }
    if (static_cast<size_t>(T) > n) T = static_cast<int>(n);
    std::vector<std::thread> pool;
    size_t chunk = (n + T - 1) / T;
    for (int t = 0; t < T; ++t) {
        size_t lo = static_cast<size_t>(t) * chunk, hi = std::min(n, lo + chunk);
        if (lo >= hi) break;
        pool.emplace_back([&, lo, hi] { for (size_t i = lo; i < hi; ++i) out[i] = one(i); });
    }
    for (auto& th : pool) th.join();
    return out;
}

}  // namespace

// Batch aa Pgen over many sequences, parallelized across sequences. Each output is the exact
// per-sequence pgen_aa (mismatches=0) or pgen_aa_hamming1 (mismatches=1) — the parallelism is
// embarrassingly disjoint (no shared state, no reduction), so the result is bitwise-identical to
// the serial computation regardless of thread count. This is the clean exact speedup for the real
// workload (Pgen/1mm over many clonotypes, e.g. TCRnet / biomarker matching); the packed model is
// shared read-only across threads. Per-sequence v_idxs/j_idxs (empty → all -1 = gene-agnostic).
std::vector<double> pgen_aa_batch(const PackedModel& m, const std::vector<std::string>& seqs,
                                  const std::vector<int>& v_idxs, const std::vector<int>& j_idxs,
                                  int mismatches, int nthreads) {
    auto vi = [&](size_t i) { return v_idxs.empty() ? -1 : v_idxs[i]; };
    auto ji = [&](size_t i) { return j_idxs.empty() ? -1 : j_idxs[i]; };
    return run_batch(seqs.size(), nthreads, [&](size_t i) {
        return mismatches == 1 ? pgen_aa_hamming1(m, seqs[i], vi(i), ji(i))
                               : pgen_aa(m, seqs[i], vi(i), ji(i));
    });
}

// Batch nucleotide Pgen, parallelized across sequences exactly as ``pgen_aa_batch`` is. The
// encoding pass is deliberately serial and ahead of the threads: `run_batch`'s callable must not
// throw (an escaping exception in a worker terminates the process), and a non-ACGT base must raise
// rather than score 0 -- a silent zero here is indistinguishable from a real answer.
std::vector<double> pgen_nt_batch(const PackedModel& m, const std::vector<std::string>& seqs,
                                  const std::vector<int>& v_idxs, const std::vector<int>& j_idxs,
                                  int nthreads) {
    std::vector<std::vector<int8_t>> code(seqs.size());
    for (size_t i = 0; i < seqs.size(); ++i) {
        code[i].resize(seqs[i].size());
        for (size_t p = 0; p < seqs[i].size(); ++p) {
            switch (seqs[i][p]) {
                case 'A': code[i][p] = 0; break;
                case 'C': code[i][p] = 1; break;
                case 'G': code[i][p] = 2; break;
                case 'T': code[i][p] = 3; break;
                default:
                    throw std::invalid_argument("pgen_nt_batch: sequence " + std::to_string(i) +
                                                " has a non-ACGT base");
            }
        }
    }
    auto vi = [&](size_t i) { return v_idxs.empty() ? -1 : v_idxs[i]; };
    auto ji = [&](size_t i) { return j_idxs.empty() ? -1 : j_idxs[i]; };
    return run_batch(seqs.size(), nthreads,
                     [&](size_t i) { return pgen_nt(m, code[i], vi(i), ji(i)); });
}

namespace {

// One codon -> its residue, or 'X' when it is not clean ACGT (the legacy converters' contract).
inline char codon_aa(char a, char b, char c) {
    auto code = [](char x) {
        switch (x) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    };
    const int i = code(a), j = code(b), k = code(c);
    return (i < 0 || j < 0 || k < 0) ? 'X' : CODON[i * 16 + j * 4 + k];
}

inline bool non_coding_marker(char c) {
    return c == '#' || c == '~' || c == '_' || c == '?' ||
           c == 'a' || c == 'c' || c == 'g' || c == 't';
}

// `translate` then `to_unified_cdr3aa`, in one pass over one sequence.
std::string translate_junction(const std::string& in) {
    if (in.empty()) return "";
    const int oof = static_cast<int>(in.size() % 3);
    std::string seq = in;
    if (oof) {                                   // pad the middle so both ends are in frame
        const size_t mid = in.size() / 2;
        seq = in.substr(0, mid) + std::string(static_cast<size_t>(3 - oof), '?') + in.substr(mid);
    }
    const int n = static_cast<int>(seq.size());
    std::string left;
    int left_end = -1;
    for (int i = 0; i + 2 < n; i += 3) {
        if (seq[i] == '?' || seq[i + 1] == '?' || seq[i + 2] == '?') { left_end = i; break; }
        left.push_back(codon_aa(seq[i], seq[i + 1], seq[i + 2]));
    }
    std::string out;
    if (oof == 0) {
        out = left;
    } else {
        std::string right;
        int right_end = -1;
        for (int i = n; i > 2; i -= 3) {
            if (seq[i - 3] == '?' || seq[i - 2] == '?' || seq[i - 1] == '?') { right_end = i; break; }
            right.push_back(codon_aa(seq[i - 3], seq[i - 2], seq[i - 1]));
        }
        std::reverse(right.begin(), right.end());
        std::string mid;
        if (left_end >= 0 && right_end > left_end)
            for (int i = left_end; i < right_end; ++i)
                mid.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(seq[i]))));
        out = left + mid + right;
    }
    // Collapse each run of non-coding markers to a single '_'.
    std::string unified;
    unified.reserve(out.size());
    bool in_run = false;
    for (char c : out) {
        if (non_coding_marker(c)) {
            if (!in_run) { unified.push_back('_'); in_run = true; }
        } else {
            unified.push_back(c);
            in_run = false;
        }
    }
    return unified;
}

}  // namespace

std::vector<std::string> translate_junctions(const std::vector<std::string>& seqs, int nthreads) {
    return run_batch(seqs.size(), nthreads, [&](size_t i) { return translate_junction(seqs[i]); });
}

namespace {

// Per-position residue sets -> per-position allowed-codon masks. ``allowed[c]`` is the string of
// residues permitted at position c; an empty string or one containing 'X' is a wildcard (any of the
// 20 amino acids, OLGA's degenerate X). Throws on a character the genetic code does not name: the
// mask would otherwise be empty and the whole query would score a silent 0.
std::vector<uint64_t> masks_for_sets(const std::vector<std::string>& allowed) {
    std::vector<uint64_t> out(allowed.size(), 0ULL);
    for (size_t c = 0; c < allowed.size(); ++c) {
        const std::string& set = allowed[c];
        if (set.empty() || set.find('X') != std::string::npos) {
            out[c] = mask_wildcard();
            continue;
        }
        uint64_t mk = 0;
        for (char x : set) {
            uint64_t one = mask_for_aa(x);
            if (one == 0ULL)
                throw std::invalid_argument(
                    "allowed[" + std::to_string(c) + "]: '" + std::string(1, x) +
                    "' is not an amino acid; use '' or 'X' for a wildcard position");
            mk |= one;
        }
        out[c] = mk;
    }
    return out;
}

}  // namespace

// Total Pgen of every amino-acid CDR3 matching a per-position residue set — the general masked DP
// (``pgen_aa`` is the all-singletons case, ``pgen_aa_hamming1`` the one-wildcard-at-a-time one).
// Cost is one transfer-matrix pass whatever the sets contain, so a motif is scored without
// enumerating the sequences it matches. See ``masks_for_sets`` for the '' / 'X' wildcard rule.
double pgen_aa_degenerate(const PackedModel& m, const std::vector<std::string>& allowed,
                          int v_idx, int j_idx) {
    std::vector<uint64_t> masks = masks_for_sets(allowed);
    return pgen_aa_masked(m, masks.data(), static_cast<int>(masks.size()), v_idx, j_idx);
}

// Batch degenerate aa Pgen, parallelized across queries exactly as ``pgen_aa_batch`` is (and
// bitwise-identical to the per-query calls). Masks are built serially up front so a bad residue
// raises here rather than inside a worker thread.
std::vector<double> pgen_aa_degenerate_batch(const PackedModel& m,
                                             const std::vector<std::vector<std::string>>& allowed,
                                             const std::vector<int>& v_idxs,
                                             const std::vector<int>& j_idxs, int nthreads) {
    std::vector<std::vector<uint64_t>> masks;
    masks.reserve(allowed.size());
    for (const auto& a : allowed) masks.push_back(masks_for_sets(a));
    auto vi = [&](size_t i) { return v_idxs.empty() ? -1 : v_idxs[i]; };
    auto ji = [&](size_t i) { return j_idxs.empty() ? -1 : j_idxs[i]; };
    return run_batch(masks.size(), nthreads, [&](size_t i) {
        return pgen_aa_masked(m, masks[i].data(), static_cast<int>(masks[i].size()), vi(i), ji(i));
    });
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

// NOTE: same nested-prefix structure as `pgen_aa_vdj`, deliberately not restructured. The output
// here is soft counts spread across a dozen accumulator arrays rather than one scalar, so the
// order-preserving second pass that keeps Pgen bitwise identical has no cheap analogue; the tests
// pin these at atol=1e-12, not bitwise. It runs during model regeneration on Aldan-3, not in the
// signature hot path.
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
                 const std::vector<int>& dmask, Counts& local, bool allow_dd = true) {
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
                        if (m.dd && m.p_nd2 > 0.0 && allow_dd) {  // n_D=2 tandem, weighted by P(n_D=2)
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
                   const std::vector<int>& dd_allowed,
                   size_t lo, size_t hi, Counts& acc) {
    Counts local = make_counts(m);
    bool have_masks = !vmasks.empty();
    bool have_dd_gate = !dd_allowed.empty();
    std::vector<int> empty;
    double ll = 0.0;
    for (size_t i = lo; i < hi; ++i) {
        zero(local);
        const std::vector<int>& vm = have_masks ? vmasks[i] : empty;
        const std::vector<int>& jm = have_masks ? jmasks[i] : empty;
        const std::vector<int>& dm = (have_masks && !dmasks.empty()) ? dmasks[i] : empty;
        bool allow_dd = have_dd_gate ? dd_allowed[i] != 0 : true;
        double total = estep_one(m, seqs[i], vm, jm, dm, local, allow_dd);
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
                   int nthreads,
                   const std::vector<int>& dd_allowed) {
    size_t n = seqs.size();
    int T = nthreads;
    if (T <= 0) {
        unsigned hw = std::thread::hardware_concurrency();
        T = hw > 3 ? static_cast<int>(hw) - 2 : 1;
    }
    if (n < kEstepThreadMin || T <= 1) return estep_range(m, seqs, vmasks, jmasks, dmasks, dd_allowed, 0, n, counts);
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
            lls[t] = estep_range(m, seqs, vmasks, jmasks, dmasks, dd_allowed, lo, hi, acc[t]);
        });
    }
    for (auto& th : pool) th.join();
    double ll = 0.0;  // reduce in fixed thread order → deterministic for a given nthreads
    for (int t = 0; t < T; ++t) { ll += lls[t]; add_scaled(counts, acc[t], 1.0); }
    return ll;
}


// ---- aa argmax: the same Pi_L*Pi_R transfer matrix with `max` instead of the sums -------------
//
// `pgen_aa` marginalizes over every recombination AND every synonymous codon assignment. The most
// likely nucleotide CDR3 is the maximum of the same product, so this is a faithful mirror of
// mk_left_tm / mk_right_tm / combine_tm / pgen_aa_vdj with `+=` replaced by `max` — same states,
// same factors, same codon masks, same pruning. Nothing here re-derives the model.
//
// It returns SCENARIOS rather than nucleotide paths. Once the scenario is known the free positions
// are two short insertion blocks, and picking their best codons is a cheap per-scenario DP that
// Python already has; carrying full back-pointers through the sweep would cost far more than the
// handful of reconstructions the caller actually needs.
namespace {

// Which (V, delV) / (J, delJ) achieved the max at each (position, state). Parallel to Lf/Rb.
struct ArgTM {
    std::vector<double> w;        // [(N+1)*25]
    std::vector<int> gene, len;   // the winning gene index and germline length
    std::vector<char> any;        // [N+1] — is any state live at this boundary
};

inline void relax(ArgTM& t, int idx, double w, int gene, int len) {
    if (w > t.w[idx]) { t.w[idx] = w; t.gene[idx] = gene; t.len[idx] = len; }
}

// Max counterpart of mk_left_tm: best V germline (>=1 nt) + N1 insertion filling [0,p).
void mk_left_max(const PackedModel& m, const uint64_t* allowed, int N, int v_idx,
                 const std::vector<double>& pins, const std::vector<double>& R,
                 const std::vector<double>& bias, int jbind, ArgTM& L) {
    int W = (N + 1) * 25;
    L.w.assign(W, 0.0); L.gene.assign(W, -1); L.len.assign(W, 0);
    ArgTM seed;
    seed.w.assign(W, 0.0); seed.gene.assign(W, -1); seed.len.assign(W, 0);
    std::vector<char> has_seed(N + 1, 0);
    for (int v : mask_or_all(v_idx, m.func_v)) {
        double pv = m.pv[v];
        if (pv == 0.0) continue;
        if (jbind >= 0) {                       // VJ: fold P(J|V) into the V seed
            double pjv = m.pjv[v * m.nJ() + jbind];
            if (pjv == 0.0) continue;
            pv *= pjv;
        }
        const auto& gv = m.cut_v[v];
        int Lv = static_cast<int>(gv.size());
        for (int len_v = 1; len_v <= std::min(Lv, N); ++len_v) {
            int div = Lv - len_v;
            if (div < 0 || div >= m.nbins_v) continue;
            double pdv = m.del_v[v * m.nbins_v + div];
            if (pdv == 0.0) continue;
            bool ok = true;                     // complete codons within [0,len_v) must translate
            for (int c = 0; c < len_v / 3; ++c)
                if (!ok_codon(allowed[c], gv[3 * c] * 16 + gv[3 * c + 1] * 4 + gv[3 * c + 2])) { ok = false; break; }
            if (!ok) continue;
            int a = gv[len_v - 1], b = (len_v >= 2) ? gv[len_v - 2] : -1;
            relax(seed, len_v * 25 + sidx(a, b), pv * pdv, v, len_v);
            has_seed[len_v] = 1;
        }
    }
    std::vector<double> cur(25), ncur(25);
    std::vector<int> cg(25), cl(25), ng(25), nl(25);
    for (int q = 0; q <= N; ++q) {
        if (!has_seed[q]) continue;
        for (int s = 0; s < 25; ++s) {
            cur[s] = seed.w[q * 25 + s]; cg[s] = seed.gene[q * 25 + s]; cl[s] = seed.len[q * 25 + s];
        }
        if (!pins.empty() && pins[0] != 0.0)
            for (int s = 0; s < 25; ++s)
                if (cur[s] > 0.0) relax(L, q * 25 + s, pins[0] * cur[s], cg[s], cl[s]);
        for (int ell = 1; q + ell <= N && ell < static_cast<int>(pins.size()); ++ell) {
            int mm = q + ell - 1;
            bool ce = (mm % 3 == 2);
            uint64_t am = allowed[mm / 3];
            std::fill(ncur.begin(), ncur.end(), 0.0);
            std::fill(ng.begin(), ng.end(), -1);
            std::fill(nl.begin(), nl.end(), 0);
            bool any = false;
            for (int s = 0; s < 25; ++s) {
                double w = cur[s];
                if (w == 0.0) continue;
                int p1 = s / 5 - 1, p2 = s % 5 - 1;
                for (int nt = 0; nt < 4; ++nt) {
                    double mult = (ell == 1) ? bias[nt] : R[nt * 4 + p1];
                    if (mult == 0.0) continue;
                    if (ce && !ok_codon(am, p2 * 16 + p1 * 4 + nt)) continue;
                    int t = sidx(nt, p1);
                    double nw = w * mult;
                    if (nw > ncur[t]) { ncur[t] = nw; ng[t] = cg[s]; nl[t] = cl[s]; }
                    any = true;
                }
            }
            cur.swap(ncur); cg.swap(ng); cl.swap(nl);
            if (!any) break;
            if (pins[ell] != 0.0)
                for (int s = 0; s < 25; ++s)
                    if (cur[s] > 0.0) relax(L, (q + ell) * 25 + s, pins[ell] * cur[s], cg[s], cl[s]);
        }
    }
    L.any.assign(N + 1, 0);
    for (int p = 0; p <= N; ++p)
        for (int s = 0; s < 25; ++s)
            if (L.w[p * 25 + s] != 0.0) { L.any[p] = 1; break; }
}

// Max counterpart of mk_right_tm: best DJ insertion + J germline filling [p,N), weighted by
// P(D|J)P(J)P(delJ). D<0 leaves out the P(D|J) factor (the VJ chain has none).
void mk_right_max(const PackedModel& m, const uint64_t* allowed, int N, int D, int j_idx, ArgTM& Rt) {
    int W = (N + 1) * 25;
    Rt.w.assign(W, 0.0); Rt.gene.assign(W, -1); Rt.len.assign(W, 0);
    const auto& pins = m.ins_dj; const auto& R = m.R_dj; const auto& bias = m.bias_dj;
    int nD = m.nD();
    ArgTM seed;
    seed.w.assign(W, 0.0); seed.gene.assign(W, -1); seed.len.assign(W, 0);
    std::vector<char> has_seed(N + 1, 0);
    for (int j : mask_or_all(j_idx, m.func_j)) {
        double pj = m.pj[j];
        if (pj == 0.0) continue;
        double pdg = (D < 0) ? 1.0 : m.pd_given_j[j * nD + D];
        if (pdg == 0.0) continue;               // genomically impossible D-J pair: prune
        const auto& gj = m.cut_j[j];
        int Lj = static_cast<int>(gj.size());
        for (int len_j = 1; len_j <= std::min(Lj, N); ++len_j) {
            int right = N - len_j, idxj = Lj - len_j;
            if (idxj < 0 || idxj >= m.nbins_j) continue;
            double pdj = m.del_j[j * m.nbins_j + idxj];
            if (pdj == 0.0) continue;
            bool ok = true;                     // complete codons within [right,N) must translate
            for (int c = (right + 2) / 3; c < N / 3; ++c) {
                if (3 * c >= right) {
                    int o = 3 * c - right;
                    if (!ok_codon(allowed[c], gj[idxj + o] * 16 + gj[idxj + o + 1] * 4 + gj[idxj + o + 2])) { ok = false; break; }
                }
            }
            if (!ok) continue;
            int cc = gj[idxj], dd = (len_j >= 2) ? gj[idxj + 1] : -1;
            relax(seed, right * 25 + sidx(cc, dd), pdg * pj * pdj, j, len_j);
            has_seed[right] = 1;
        }
    }
    std::vector<double> cur(25), ncur(25);
    std::vector<int> cg(25), cl(25), ng(25), nl(25);
    for (int right = N; right >= 0; --right) {
        if (!has_seed[right]) continue;
        for (int s = 0; s < 25; ++s) {
            cur[s] = seed.w[right * 25 + s]; cg[s] = seed.gene[right * 25 + s]; cl[s] = seed.len[right * 25 + s];
        }
        if (!pins.empty() && pins[0] != 0.0)
            for (int s = 0; s < 25; ++s)
                if (cur[s] > 0.0) relax(Rt, right * 25 + s, pins[0] * cur[s], cg[s], cl[s]);
        for (int ell = 1; right - ell >= 0 && ell < static_cast<int>(pins.size()); ++ell) {
            int mm = right - ell;
            bool ce = (mm % 3 == 0);
            uint64_t am = allowed[mm / 3];
            std::fill(ncur.begin(), ncur.end(), 0.0);
            std::fill(ng.begin(), ng.end(), -1);
            std::fill(nl.begin(), nl.end(), 0);
            bool any = false;
            for (int s = 0; s < 25; ++s) {
                double w = cur[s];
                if (w == 0.0) continue;
                int c = s / 5 - 1, d = s % 5 - 1;
                for (int nt = 0; nt < 4; ++nt) {
                    double mult = (ell == 1) ? bias[nt] : R[nt * 4 + c];
                    if (mult == 0.0) continue;
                    if (ce && !ok_codon(am, nt * 16 + c * 4 + d)) continue;
                    int t = sidx(nt, c);
                    double nw = w * mult;
                    if (nw > ncur[t]) { ncur[t] = nw; ng[t] = cg[s]; nl[t] = cl[s]; }
                    any = true;
                }
            }
            cur.swap(ncur); cg.swap(ng); cl.swap(nl);
            if (!any) break;
            if (pins[ell] != 0.0)
                for (int s = 0; s < 25; ++s)
                    if (cur[s] > 0.0) relax(Rt, (right - ell) * 25 + s, pins[ell] * cur[s], cg[s], cl[s]);
        }
    }
    Rt.any.assign(N + 1, 0);
    for (int p = 0; p <= N; ++p)
        for (int s = 0; s < 25; ++s)
            if (Rt.w[p * 25 + s] != 0.0) { Rt.any[p] = 1; break; }
}

// Max counterpart of combine_tm: best join of a left state (nt[p-1],nt[p-2]) with a right state
// (nt[p],nt[p+1]) at boundary p, honouring the single codon that straddles p. Writes back which
// left/right states won so the caller can read off (V,delV) and (J,delJ).
double combine_max(const double* lw, const double* rw, int p, const uint64_t* allowed,
                   int& bs_l, int& bs_r) {
    double best = 0.0; bs_l = bs_r = -1;
    int phase = p % 3;
    uint64_t am = 0;
    if (phase == 1) am = allowed[(p - 1) / 3];
    else if (phase == 2) am = allowed[p / 3];
    for (int sl = 0; sl < 25; ++sl) {
        double a = lw[sl];
        if (a == 0.0) continue;
        int la = sl / 5 - 1, lb = sl % 5 - 1;
        for (int sr = 0; sr < 25; ++sr) {
            double b = rw[sr];
            if (b == 0.0) continue;
            double w = a * b;
            if (w <= best) continue;
            int rc = sr / 5 - 1, rd = sr % 5 - 1;
            if (phase == 1) {                   // codon [p-1,p,p+1] = (la, rc, rd)
                if (la < 0 || rc < 0 || rd < 0) continue;
                if (!ok_codon(am, la * 16 + rc * 4 + rd)) continue;
            } else if (phase == 2) {            // codon [p-2,p-1,p] = (lb, la, rc)
                if (la < 0 || lb < 0 || rc < 0) continue;
                if (!ok_codon(am, lb * 16 + la * 4 + rc)) continue;
            }
            best = w; bs_l = sl; bs_r = sr;
        }
    }
    return best;
}

// Keep the k best scenarios seen, cheapest possible: k is single digits.
void offer(std::vector<AaScenario>& top, size_t k, const AaScenario& s) {
    if (s.w <= 0.0) return;
    if (top.size() < k) { top.push_back(s); }
    else {
        size_t worst = 0;
        for (size_t i = 1; i < top.size(); ++i) if (top[i].w < top[worst].w) worst = i;
        if (s.w <= top[worst].w) return;
        top[worst] = s;
    }
}

void best_vdj(const PackedModel& m, const uint64_t* allowed, int alen, int v_idx, int j_idx,
              size_t k, std::vector<AaScenario>& top) {
    int N = 3 * alen;
    ArgTM L, Rt;
    mk_left_max(m, allowed, N, v_idx, m.ins_vd, m.R_vd, m.bias_vd, -1, L);
    std::vector<double> dp(25), ndp(25);
    for (int D : m.func_d) {
        mk_right_max(m, allowed, N, D, j_idx, Rt);
        const auto& cutd = m.cut_d[D];
        int Ld = static_cast<int>(cutd.size());
        for (int idx5 = 0; idx5 <= Ld && idx5 < m.nbins_d5; ++idx5) {
            std::vector<int> vg(25), vl(25), ng(25), nl(25);
            for (int pos = 1; pos <= N; ++pos) {
                if (!L.any[pos]) continue;
                for (int s = 0; s < 25; ++s) {
                    dp[s] = L.w[pos * 25 + s];
                    vg[s] = L.gene[pos * 25 + s];
                    vl[s] = L.len[pos * 25 + s];
                }
                // Thread the D germline one nt at a time and offer a join after each prefix: the
                // prefix of length `ld` IS the trim (idx5, Ld-idx5-ld), so every 3' trim of this
                // (D, idx5) is covered by one forward pass instead of one pass per trim.
                // (V, delV) rides along, because the join's winning state is only known after the
                // D has been threaded and the identity would otherwise be lost.
                for (int ld = 0; ld <= Ld - idx5 && pos + ld <= N; ++ld) {
                    if (ld > 0) {
                        int mm = pos + ld - 1, nt = cutd[idx5 + ld - 1];
                        bool ce = (mm % 3 == 2);
                        uint64_t am = allowed[mm / 3];
                        std::fill(ndp.begin(), ndp.end(), 0.0);
                        std::fill(ng.begin(), ng.end(), -1);
                        std::fill(nl.begin(), nl.end(), 0);
                        bool any = false;
                        for (int s = 0; s < 25; ++s) {
                            double w = dp[s];
                            if (w == 0.0) continue;
                            int p1 = s / 5 - 1, p2 = s % 5 - 1;
                            if (ce && !ok_codon(am, p2 * 16 + p1 * 4 + nt)) continue;
                            int t = sidx(nt, p1);
                            if (w > ndp[t]) { ndp[t] = w; ng[t] = vg[s]; nl[t] = vl[s]; }
                            any = true;
                        }
                        dp.swap(ndp); vg.swap(ng); vl.swap(nl);
                        if (!any) break;
                    }
                    int idx3 = Ld - idx5 - ld;
                    if (idx3 < 0 || idx3 >= m.nbins_d3) continue;
                    double pdel = m.del_d[(D * m.nbins_d5 + idx5) * m.nbins_d3 + idx3];
                    if (pdel == 0.0) continue;
                    int p = pos + ld;
                    if (!Rt.any[p]) continue;
                    int sl = -1, sr = -1;
                    double w = combine_max(&dp[0], &Rt.w[p * 25], p, allowed, sl, sr);
                    if (w <= 0.0 || sl < 0 || sr < 0) continue;
                    AaScenario sc;
                    sc.w = m.p_nd1 * pdel * w;
                    sc.d = D; sc.idx5 = idx5; sc.idx3 = idx3; sc.pos = pos;
                    sc.v = vg[sl]; sc.len_v = vl[sl];
                    sc.j = Rt.gene[p * 25 + sr]; sc.len_j = Rt.len[p * 25 + sr];
                    offer(top, k, sc);
                }
            }
        }
    }
}

void best_vj(const PackedModel& m, const uint64_t* allowed, int alen, int v_idx, int j_idx,
             size_t k, std::vector<AaScenario>& top) {
    int N = 3 * alen;
    ArgTM L;
    std::vector<double> dp(25), ndp(25);
    for (int j : mask_or_all(j_idx, m.func_j)) {
        mk_left_max(m, allowed, N, v_idx, m.ins_vj, m.R_vj, m.bias_vj, j, L);
        const auto& gj = m.cut_j[j];
        int Lj = static_cast<int>(gj.size());
        for (int len_j = 1; len_j <= std::min(Lj, N); ++len_j) {
            int p = N - len_j, idxj = Lj - len_j;
            if (p < 1 || !L.any[p]) continue;   // >=1 V/insertion nt before J
            if (idxj < 0 || idxj >= m.nbins_j) continue;
            double pdj = m.del_j[j * m.nbins_j + idxj];
            if (pdj == 0.0) continue;
            for (int s = 0; s < 25; ++s) dp[s] = L.w[p * 25 + s];
            std::vector<int> vg(25), vl(25), ng(25), nl(25);
            for (int s = 0; s < 25; ++s) { vg[s] = L.gene[p * 25 + s]; vl[s] = L.len[p * 25 + s]; }
            bool ok = true;
            for (int kk = 0; kk < len_j; ++kk) {   // thread the fixed J germline suffix
                int mm = p + kk, nt = gj[idxj + kk];
                bool ce = (mm % 3 == 2);
                uint64_t am = allowed[mm / 3];
                std::fill(ndp.begin(), ndp.end(), 0.0);
                std::fill(ng.begin(), ng.end(), -1);
                std::fill(nl.begin(), nl.end(), 0);
                bool any = false;
                for (int s = 0; s < 25; ++s) {
                    double w = dp[s];
                    if (w == 0.0) continue;
                    int p1 = s / 5 - 1, p2 = s % 5 - 1;
                    if (ce && !ok_codon(am, p2 * 16 + p1 * 4 + nt)) continue;
                    int t = sidx(nt, p1);
                    if (w > ndp[t]) { ndp[t] = w; ng[t] = vg[s]; nl[t] = vl[s]; }
                    any = true;
                }
                dp.swap(ndp); vg.swap(ng); vl.swap(nl);
                if (!any) { ok = false; break; }
            }
            if (!ok) continue;
            int bs = -1;
            double best = 0.0;
            for (int s = 0; s < 25; ++s) if (dp[s] > best) { best = dp[s]; bs = s; }
            if (bs < 0) continue;
            AaScenario sc;
            sc.w = pdj * best;
            sc.v = vg[bs]; sc.len_v = vl[bs];
            sc.j = j; sc.len_j = len_j;
            offer(top, k, sc);
        }
    }
}

}  // namespace

std::vector<AaScenario> best_aa_scenarios(const PackedModel& m, const std::string& aa,
                                          int v_idx, int j_idx, int k) {
    int L = static_cast<int>(aa.size());
    std::vector<AaScenario> top;
    if (L == 0 || k <= 0) return top;
    std::vector<uint64_t> allowed(L);
    for (int c = 0; c < L; ++c) {
        allowed[c] = mask_for_aa(aa[c]);
        if (allowed[c] == 0ULL) return top;      // an unknown residue has no codons
    }
    if (m.vdj) best_vdj(m, allowed.data(), L, v_idx, j_idx, static_cast<size_t>(k), top);
    else best_vj(m, allowed.data(), L, v_idx, j_idx, static_cast<size_t>(k), top);
    std::sort(top.begin(), top.end(),
              [](const AaScenario& a, const AaScenario& b) { return a.w > b.w; });
    return top;
}

AaScenarioBatch best_aa_scenarios_batch(const PackedModel& m, const std::vector<std::string>& seqs,
                                        const std::vector<int>& v_idxs,
                                        const std::vector<int>& j_idxs, int k, int nthreads) {
    auto vi = [&](size_t i) { return v_idxs.empty() ? -1 : v_idxs[i]; };
    auto ji = [&](size_t i) { return j_idxs.empty() ? -1 : j_idxs[i]; };
    auto per = run_batch(seqs.size(), nthreads, [&](size_t i) {
        return best_aa_scenarios(m, seqs[i], vi(i), ji(i), k);
    });
    size_t total = 0;
    for (const auto& one : per) total += one.size();
    AaScenarioBatch out;
    for (auto* c : {&out.row, &out.rank, &out.v, &out.len_v, &out.j, &out.len_j,
                    &out.d, &out.idx5, &out.idx3, &out.pos}) c->reserve(total);
    out.w.reserve(total);
    for (size_t i = 0; i < per.size(); ++i) {
        for (size_t r = 0; r < per[i].size(); ++r) {
            const AaScenario& sc = per[i][r];
            out.row.push_back(static_cast<int>(i));
            out.rank.push_back(static_cast<int>(r));
            out.w.push_back(sc.w);
            out.v.push_back(sc.v); out.len_v.push_back(sc.len_v);
            out.j.push_back(sc.j); out.len_j.push_back(sc.len_j);
            out.d.push_back(sc.d);
            out.idx5.push_back(sc.idx5); out.idx3.push_back(sc.idx3); out.pos.push_back(sc.pos);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------------------------
// infer_nt: the scenario search, the codon reconstruction and the marginal re-score, in one call.
//
// `best_aa_scenarios` above returns SCENARIOS; turning one into nucleotides is a 25-state
// max-product DP over the CDR3, and until #181 that DP lived in Python (`viterbi._aa_dp_max`).
// Measured there at 35% of a human TRB row and 87% of a TRA one -- and it held the GIL, which is
// what capped the batch entry point at 1.33x however many threads it was handed. It is the same DP
// here, factor for factor, including its tie-break: two paths of exactly equal weight resolve to
// the lexicographically larger nucleotide string, so a run stays reproducible.
// ---------------------------------------------------------------------------------------------

namespace {

constexpr int kMaxNt = 3 * 128;  // CDR3 nucleotides; a longer junction is not a junction

enum : int8_t { kGerm = 0, kLeft = 1, kRight = 2 };

// One scenario laid out per CDR3 position: the germline nucleotide, or -1 plus which insertion
// junction scores the position when it is free. `kLeft` reads 5'->3' (VD, or VJ on a VJ chain),
// `kRight` reads 3'->5' (DJ) -- the same two spec kinds `viterbi._d_template` builds.
struct Layout {
    int N = 0;
    int8_t tpl[kMaxNt], kind[kMaxNt], first[kMaxNt], last[kMaxNt];
    const double *Rl = nullptr, *bl = nullptr, *Rr = nullptr, *br = nullptr;
};

bool layout_of(const PackedModel& m, int N, const AaScenario& sc, Layout& o) {
    if (N <= 0 || N > kMaxNt) return false;
    if (sc.v < 0 || sc.v >= m.nV() || sc.j < 0 || sc.j >= m.nJ()) return false;
    const auto& cv = m.cut_v[sc.v];
    const auto& cj = m.cut_j[sc.j];
    int lcv = static_cast<int>(cv.size()), lcj = static_cast<int>(cj.size());
    if (sc.len_v < 0 || sc.len_j < 0) return false;
    if (sc.len_v > lcv || sc.len_j > lcj || sc.len_v + sc.len_j > N) return false;

    o.N = N;
    std::fill(o.tpl, o.tpl + N, static_cast<int8_t>(-1));
    std::fill(o.kind, o.kind + N, static_cast<int8_t>(kGerm));
    std::fill(o.first, o.first + N, static_cast<int8_t>(0));
    std::fill(o.last, o.last + N, static_cast<int8_t>(0));
    for (int p = 0; p < sc.len_v; ++p) o.tpl[p] = cv[p];
    int right = N - sc.len_j;
    for (int p = 0; p < sc.len_j; ++p) o.tpl[right + p] = cj[lcj - sc.len_j + p];

    auto block = [&](int lo, int hi, int8_t kind) {
        for (int p = lo; p < hi; ++p) {
            o.kind[p] = kind;
            o.first[p] = static_cast<int8_t>(p == lo);
            o.last[p] = static_cast<int8_t>(p == hi - 1);
        }
    };
    if (sc.d < 0) {  // VJ chain: the whole gap is one VJ insertion
        int ins = N - sc.len_v - sc.len_j;
        if (ins >= static_cast<int>(m.ins_vj.size()) || m.ins_vj[ins] == 0.0) return false;
        block(sc.len_v, sc.len_v + ins, kLeft);
        o.Rl = m.R_vj.data();
        o.bl = m.bias_vj.data();
        return true;
    }
    if (sc.d >= m.nD()) return false;
    const auto& cd = m.cut_d[sc.d];
    int ld = std::max(0, static_cast<int>(cd.size()) - sc.idx5 - sc.idx3);
    if (sc.pos < sc.len_v || sc.pos + ld > right) return false;   // lvd >= 0 and ldj >= 0
    for (int p = 0; p < ld; ++p) o.tpl[sc.pos + p] = cd[sc.idx5 + p];
    block(sc.len_v, sc.pos, kLeft);
    block(sc.pos + ld, right, kRight);
    o.Rl = m.R_vd.data();
    o.bl = m.bias_vd.data();
    o.Rr = m.R_dj.data();
    o.br = m.bias_dj.data();
    return true;
}

// The max-product codon DP for one laid-out scenario: germline positions are pinned to their
// segment and each free position takes the nucleotide maximising the insertion model covering it,
// subject to every complete codon translating to the requested residue.
//
// State s = (nt[i-1], nt[i-2]) packed as (p1+1)*5 + (p2+1); placing `nt` moves to (nt+1)*5+(p1+1).
// The states carry NO path -- a backpointer table does, and a position's nucleotide is read
// straight off its state -- so the hot loop copies no strings. A tie (two paths at exactly equal
// weight) is resolved by walking both backpointer chains, which costs nothing on the paths that do
// not tie, and is the same lexicographic-max rule the Python reference used.
struct CodonDP {
    std::vector<int8_t> bp;  // [i*25 + s] = the state at position i-1 this path came from
    double w[25], nw[25];
    char pa[kMaxNt], pb[kMaxNt];

    void trace(int i, int s, char* buf) const {
        for (int p = i; p >= 0; --p) {
            buf[p] = "ACGT"[s / 5 - 1];
            s = bp[static_cast<size_t>(p) * 25 + s];
        }
    }

    double run(const Layout& L, const std::string& aa, std::string& out) {
        const int N = L.N;
        bp.assign(static_cast<size_t>(N) * 25, 0);
        std::fill(w, w + 25, 0.0);
        w[0] = 1.0;  // (nt[-1], nt[-2]) = (-1, -1)
        for (int i = 0; i < N; ++i) {
            std::fill(nw, nw + 25, 0.0);
            int8_t* bpi = &bp[static_cast<size_t>(i) * 25];
            const bool codon_end = (i % 3 == 2);
            const char aa_i = aa[i / 3];
            const int lo = L.tpl[i] < 0 ? 0 : L.tpl[i];
            const int hi = L.tpl[i] < 0 ? 3 : L.tpl[i];
            for (int s = 0; s < 25; ++s) {
                if (w[s] == 0.0) continue;
                const int p1 = s / 5 - 1, p2 = s % 5 - 1;
                for (int nt = lo; nt <= hi; ++nt) {
                    double ww = w[s];
                    if (L.kind[i] == kLeft) {
                        ww *= L.first[i] ? L.bl[nt] : L.Rl[nt * 4 + p1];
                    } else if (L.kind[i] == kRight) {
                        if (!L.first[i]) ww *= L.Rr[p1 * 4 + nt];
                        if (L.last[i]) ww *= L.br[nt];
                    }
                    if (ww <= 0.0) continue;
                    if (codon_end && CODON[p2 * 16 + p1 * 4 + nt] != aa_i) continue;
                    const int ns = (nt + 1) * 5 + (p1 + 1);
                    if (ww > nw[ns]) {
                        nw[ns] = ww;
                        bpi[ns] = static_cast<int8_t>(s);
                    } else if (ww == nw[ns] && i > 0) {
                        trace(i - 1, s, pa);
                        trace(i - 1, bpi[ns], pb);
                        if (std::memcmp(pa, pb, static_cast<size_t>(i)) > 0)
                            bpi[ns] = static_cast<int8_t>(s);
                    }
                }
            }
            std::copy(nw, nw + 25, w);
            bool any = false;
            for (int s = 0; s < 25 && !any; ++s) any = w[s] != 0.0;
            if (!any) return 0.0;
        }
        int bs = -1;
        for (int s = 0; s < 25; ++s) {
            if (w[s] == 0.0) continue;
            if (bs < 0 || w[s] > w[bs]) {
                bs = s;
            } else if (w[s] == w[bs]) {
                trace(N - 1, s, pa);
                trace(N - 1, bs, pb);
                if (std::memcmp(pa, pb, static_cast<size_t>(N)) > 0) bs = s;
            }
        }
        if (bs < 0) return 0.0;
        trace(N - 1, bs, pa);
        out.assign(pa, pa + N);
        return w[bs];
    }
};

// -2 = "whatever allele this scenario chose"; anything else is the caller's own pin (-1 included,
// which marginalizes over the segment). A call naming exactly one allele is a pin; a call offering
// several is not, so the re-score follows the winner.
constexpr int kFromScenario = -2;

int pin_of(const std::vector<int>& alts) {
    return alts.size() == 1 ? alts[0] : kFromScenario;
}

bool alts_ok(const std::vector<int>& alts, int n) {
    for (int a : alts) if (a < -1 || a >= n) return false;
    return true;
}

InferNt infer_nt_one(const PackedModel& m, const std::string& aa, const std::vector<int>& valts,
                     const std::vector<int>& jalts, int k, CodonDP& dp) {
    InferNt r;
    const int N = 3 * static_cast<int>(aa.size());
    if (aa.empty() || k <= 0 || N > kMaxNt) return r;
    if (!alts_ok(valts, m.nV()) || !alts_ok(jalts, m.nJ())) return r;

    // Every (V alternative, J alternative) the call offers, pooled. Ragged on purpose: a multi-
    // allele call is several searches whose candidates compete on the model's own weight.
    std::vector<AaScenario> pool;
    for (int vv : valts) {
        for (int jj : jalts) {
            auto got = best_aa_scenarios(m, aa, vv, jj, k);
            pool.insert(pool.end(), got.begin(), got.end());
        }
    }
    std::stable_sort(pool.begin(), pool.end(),
                     [](const AaScenario& a, const AaScenario& b) { return a.w > b.w; });

    // Distinct nucleotide sequence -> the best-weight scenario that produced it. A linear scan,
    // because the pool is |V alts| x |J alts| x k entries -- single digits in every real call.
    std::vector<std::pair<std::string, AaScenario>> best;
    Layout L;
    std::string nt;
    for (const AaScenario& sc : pool) {
        if (!layout_of(m, N, sc, L)) continue;
        if (dp.run(L, aa, nt) <= 0.0) continue;
        auto it = std::find_if(best.begin(), best.end(),
                               [&](const std::pair<std::string, AaScenario>& kv) {
                                   return kv.first == nt;
                               });
        if (it == best.end()) best.emplace_back(nt, sc);
        else if (sc.w > it->second.w) it->second = sc;
    }
    if (best.empty()) return r;
    // Distinct candidates reconstructed, BEFORE the top-k truncation: that is the count
    // `Scenario.n_candidates` reports, and a multi-allele call pools more of them than k.
    r.n_candidates = static_cast<int>(best.size());
    std::sort(best.begin(), best.end(), [](const std::pair<std::string, AaScenario>& a,
                                           const std::pair<std::string, AaScenario>& b) {
        return a.second.w != b.second.w ? a.second.w > b.second.w : a.first < b.first;
    });
    if (static_cast<int>(best.size()) > k) best.resize(static_cast<size_t>(k));

    // Stage 2: the exact marginal. Stage 1 maximises the JOINT P(nt, scenario); this function's
    // contract is the marginal P(nt), and the re-score is what turns one into the other.
    const int pv = pin_of(valts), pj = pin_of(jalts);
    std::vector<int8_t> code(static_cast<size_t>(N));
    std::vector<std::pair<double, size_t>> scored(best.size());
    for (size_t c = 0; c < best.size(); ++c) {
        const std::string& s = best[c].first;
        for (int p = 0; p < N; ++p)
            code[p] = static_cast<int8_t>(s[p] == 'A' ? 0 : s[p] == 'C' ? 1 : s[p] == 'G' ? 2 : 3);
        scored[c] = {pgen_nt(m, code, pv == kFromScenario ? best[c].second.v : pv,
                             pj == kFromScenario ? best[c].second.j : pj), c};
    }
    std::sort(scored.begin(), scored.end(),
              [&](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) {
                  return a.first != b.first ? a.first > b.first
                                            : best[a.second].first < best[b.second].first;
              });

    const AaScenario& sc = best[scored[0].second].second;
    const int ld = sc.d >= 0 && sc.d < m.nD()
                       ? std::max(0, static_cast<int>(m.cut_d[sc.d].size()) - sc.idx5 - sc.idx3)
                       : 0;
    r.ok = true;
    r.nt = best[scored[0].second].first;
    r.v = sc.v;
    r.len_v = sc.len_v;
    r.j = sc.j;
    r.len_j = sc.len_j;
    r.d = ld ? sc.d : -1;
    r.d_start = sc.pos;
    r.d_end = sc.pos + ld;
    r.pgen = scored[0].first;
    r.scenario_p = sc.w;
    r.runner_up_pgen = scored.size() > 1 ? scored[1].first : 0.0;
    return r;
}

}  // namespace

std::vector<InferNt> infer_nt_batch(const PackedModel& m, const std::vector<std::string>& aas,
                                    const std::vector<std::vector<int>>& v_alts,
                                    const std::vector<std::vector<int>>& j_alts,
                                    int k, int nthreads) {
    if (v_alts.size() != aas.size() || j_alts.size() != aas.size())
        throw std::invalid_argument("v_alts and j_alts must have one entry per sequence");
    return run_batch(aas.size(), nthreads, [&](size_t i) {
        CodonDP dp;  // per row, so nothing is shared between workers
        return infer_nt_one(m, aas[i], v_alts[i], j_alts[i], k, dp);
    });
}

// ---------------------------------------------------------------------------------------------
// Ancestral sampling, natively. The Python reference is a per-sequence loop of ~20 numpy calls --
// 4,200 seq/s measured, and pool generation is the dominant serial cost of a synthetic corpus
// build (77 s of single-core work per 5.2M sequences). Same Bayes net, same factors, in C++.
// ---------------------------------------------------------------------------------------------

namespace {

constexpr int kDrawAttempts = 10000;  // per slot, under `productive_only`

// splitmix64: seeds a slot's generator from (seed, slot) so a draw is a pure function of both.
inline uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// xoshiro256++ -- small state, so seeding it per slot is free (mt19937_64's 2.5 KB is not).
struct Rng {
    uint64_t s[4];

    explicit Rng(uint64_t seed, uint64_t slot) {
        uint64_t x = seed ^ (slot * 0xD1B54A32D192ED03ULL);
        for (uint64_t& w : s) w = splitmix64(x);
    }
    static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    inline uint64_t next() {
        const uint64_t r = rotl(s[0] + s[3], 23) + s[0];
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return r;
    }
    // Uniform in [0,1) from the top 53 bits, as numpy does -- same construction, different stream.
    inline double uniform() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }
};

// A categorical as (value, cumulative) with the zero-probability atoms dropped and the rest
// renormalized -- `generate._cum`'s contract.
struct Cat {
    std::vector<int> value;
    std::vector<double> cum;

    bool empty() const { return value.empty(); }
    int pick(Rng& r) const {
        double u = r.uniform();
        size_t i = static_cast<size_t>(
            std::lower_bound(cum.begin(), cum.end(), u) - cum.begin());
        return value[std::min(i, value.size() - 1)];
    }
};

Cat make_cat(const double* p, int n, const std::vector<int>* only = nullptr) {
    Cat c;
    double tot = 0.0;
    auto add = [&](int i) { if (p[i] > 0.0) { c.value.push_back(i); tot += p[i]; } };
    if (only) { for (int i : *only) add(i); } else { for (int i = 0; i < n; ++i) add(i); }
    double run = 0.0;
    for (int i : c.value) { run += p[i] / tot; c.cum.push_back(run); }
    return c;
}

// Everything the sampler needs, built once per call and shared read-only across the workers.
struct GenPrep {
    bool vdj = false;
    Cat v, j_marg;
    std::vector<Cat> j_given_v, d_given_j, del_v, del_j, del_d, del_d2, d2_given_d1;
    Cat ins_vd, ins_dj, ins_vj, ins_dd, n_d;
    double Rcum_vd[16], Rcum_dj[16], Rcum_vj[16], Rcum_dd[16];      // [prev*4 + to], cumulative
    double bcum_vd[4], bcum_dj[4], bcum_vj[4], bcum_dd[4];
};

void fill_chain(const std::vector<double>& R, const std::vector<double>& bias,
                double* Rcum, double* bcum) {
    if (R.size() < 16 || bias.size() < 4) {
        std::fill(Rcum, Rcum + 16, 0.0);
        std::fill(bcum, bcum + 4, 0.0);
        return;
    }
    for (int prev = 0; prev < 4; ++prev) {          // R is [to*4 + from]; normalize the column
        double tot = 0.0;
        for (int to = 0; to < 4; ++to) tot += R[to * 4 + prev];
        double run = 0.0;
        for (int to = 0; to < 4; ++to) {
            run += tot > 0.0 ? R[to * 4 + prev] / tot : 0.0;
            Rcum[prev * 4 + to] = run;
        }
    }
    double tot = 0.0, run = 0.0;
    for (int i = 0; i < 4; ++i) tot += bias[i];
    for (int i = 0; i < 4; ++i) { run += tot > 0.0 ? bias[i] / tot : 0.0; bcum[i] = run; }
}

int pick_chain(Rng& r, const double* cum4) {
    double u = r.uniform();
    for (int i = 0; i < 3; ++i) if (u < cum4[i]) return i;
    return 3;
}

GenPrep build_gen_prep(const PackedModel& m) {
    GenPrep g;
    g.vdj = m.vdj;
    g.v = make_cat(m.pv.data(), m.nV(), &m.func_v);
    if (m.vdj) {
        g.j_marg = make_cat(m.pj.data(), m.nJ(), &m.func_j);
        g.d_given_j.resize(m.nJ());
        for (int j : m.func_j)
            g.d_given_j[j] = make_cat(m.pd_given_j.data() + static_cast<size_t>(j) * m.nD(),
                                      m.nD(), &m.func_d);
    } else {
        g.j_given_v.resize(m.nV());
        for (int v : m.func_v)
            g.j_given_v[v] = make_cat(m.pjv.data() + static_cast<size_t>(v) * m.nJ(),
                                      m.nJ(), &m.func_j);
    }
    g.del_v.resize(m.nV());
    for (int v : m.func_v)
        g.del_v[v] = make_cat(m.del_v.data() + static_cast<size_t>(v) * m.nbins_v, m.nbins_v);
    g.del_j.resize(m.nJ());
    for (int j : m.func_j)
        g.del_j[j] = make_cat(m.del_j.data() + static_cast<size_t>(j) * m.nbins_j, m.nbins_j);
    if (m.vdj) {
        const int nd = m.nbins_d5 * m.nbins_d3;     // the 5'/3' trim is drawn jointly
        g.del_d.resize(m.nD());
        for (int d : m.func_d)
            g.del_d[d] = make_cat(m.del_d.data() + static_cast<size_t>(d) * nd, nd);
        if (m.dd) {
            g.del_d2.resize(m.nD());
            g.d2_given_d1.resize(m.nD());
            for (int d : m.func_d) {
                g.del_d2[d] = make_cat(m.del_d2.data() + static_cast<size_t>(d) * nd, nd);
                g.d2_given_d1[d] = make_cat(
                    m.pd2_given_d1.data() + static_cast<size_t>(d) * m.nD(), m.nD(), &m.func_d);
            }
            const double nd_p[2] = {m.p_nd1, m.p_nd2};
            g.n_d = make_cat(nd_p, 2);              // value 0 -> n_D=1, value 1 -> n_D=2
        }
        g.ins_vd = make_cat(m.ins_vd.data(), static_cast<int>(m.ins_vd.size()));
        g.ins_dj = make_cat(m.ins_dj.data(), static_cast<int>(m.ins_dj.size()));
        if (m.dd) g.ins_dd = make_cat(m.ins_dd.data(), static_cast<int>(m.ins_dd.size()));
        fill_chain(m.R_vd, m.bias_vd, g.Rcum_vd, g.bcum_vd);
        fill_chain(m.R_dj, m.bias_dj, g.Rcum_dj, g.bcum_dj);
        fill_chain(m.R_dd, m.bias_dd, g.Rcum_dd, g.bcum_dd);
    } else {
        g.ins_vj = make_cat(m.ins_vj.data(), static_cast<int>(m.ins_vj.size()));
        fill_chain(m.R_vj, m.bias_vj, g.Rcum_vj, g.bcum_vj);
    }
    return g;
}

// An N-region: first nucleotide from the steady-state bias, the rest from the dinucleotide chain.
// `from_right` walks 3'->5', which is how the DJ junction is parameterized.
void draw_insert(Rng& r, int len, const double* Rcum, const double* bcum, bool from_right,
                 std::string& out) {
    if (len <= 0) return;
    size_t base = out.size();
    out.resize(base + static_cast<size_t>(len));
    if (from_right) {
        int prev = pick_chain(r, bcum);
        out[base + len - 1] = "ACGT"[prev];
        for (int k = len - 2; k >= 0; --k) {
            prev = pick_chain(r, Rcum + prev * 4);
            out[base + k] = "ACGT"[prev];
        }
    } else {
        int prev = pick_chain(r, bcum);
        out[base] = "ACGT"[prev];
        for (int k = 1; k < len; ++k) {
            prev = pick_chain(r, Rcum + prev * 4);
            out[base + k] = "ACGT"[prev];
        }
    }
}

void append_seg(const std::vector<int8_t>& cut, int from, int to, std::string& out) {
    for (int p = from; p < to; ++p) out.push_back("ACGT"[cut[p]]);
}

// D's contribution under a jointly drawn 5'/3' trim. `min1` re-draws until it is non-empty, which
// is what a tandem D needs to stay in the Pgen partition; it gives up rather than looping forever.
int draw_d_seg(Rng& r, const PackedModel& m, const Cat& del, int d, bool min1, std::string& out) {
    const auto& cut = m.cut_d[d];
    const int L = static_cast<int>(cut.size());
    for (int attempt = 0; attempt < 64; ++attempt) {
        if (del.empty()) return 0;
        const int flat = del.pick(r);
        const int i5 = flat / m.nbins_d3, i3 = flat % m.nbins_d3;
        const int from = i5, to = L - i3;
        if (to > from) {
            append_seg(cut, from, to, out);
            return to - from;
        }
        if (!min1) return 0;
    }
    return 0;
}

// One draw into `nt`; returns the genes used (d/d2 = -1 when they contributed nothing).
struct Drawn { int v, d, d2, j; };

Drawn draw_one(Rng& r, const PackedModel& m, const GenPrep& g, std::string& nt) {
    nt.clear();
    Drawn out{-1, -1, -1, -1};
    out.v = g.v.pick(r);
    out.j = g.vdj ? g.j_marg.pick(r) : g.j_given_v[out.v].pick(r);

    // V and J each contribute at least one nucleotide -- the Pgen model's own invariant.
    const auto& cutv = m.cut_v[out.v];
    const int lv = std::max(1, static_cast<int>(cutv.size()) - g.del_v[out.v].pick(r));
    append_seg(cutv, 0, std::min(lv, static_cast<int>(cutv.size())), nt);

    const auto& cutj = m.cut_j[out.j];
    const int Lj = static_cast<int>(cutj.size());
    const int dj = std::min(g.del_j[out.j].pick(r), Lj - 1);

    if (!g.vdj) {
        draw_insert(r, g.ins_vj.empty() ? 0 : g.ins_vj.pick(r), g.Rcum_vj, g.bcum_vj, false, nt);
        append_seg(cutj, dj, Lj, nt);
        return out;
    }

    const int d = g.d_given_j[out.j].pick(r);
    const bool tandem = m.dd && !g.n_d.empty() && g.n_d.pick(r) == 1 &&
                        d >= 0 && !g.d2_given_d1[d].empty();
    std::string mid;
    if (tandem) {
        const int d2 = g.d2_given_d1[d].pick(r);
        std::string d1c, d2c;
        const int l1 = draw_d_seg(r, m, g.del_d[d], d, true, d1c);
        const int l2 = draw_d_seg(r, m, g.del_d2[d2], d2, true, d2c);
        draw_insert(r, g.ins_vd.pick(r), g.Rcum_vd, g.bcum_vd, false, mid);
        mid += d1c;
        draw_insert(r, g.ins_dd.empty() ? 0 : g.ins_dd.pick(r), g.Rcum_dd, g.bcum_dd, false, mid);
        mid += d2c;
        draw_insert(r, g.ins_dj.pick(r), g.Rcum_dj, g.bcum_dj, true, mid);
        out.d = l1 ? d : -1;
        out.d2 = l2 ? d2 : -1;
    } else {
        std::string d1c;
        const int l1 = draw_d_seg(r, m, g.del_d[d], d, false, d1c);
        draw_insert(r, g.ins_vd.pick(r), g.Rcum_vd, g.bcum_vd, false, mid);
        mid += d1c;
        draw_insert(r, g.ins_dj.pick(r), g.Rcum_dj, g.bcum_dj, true, mid);
        out.d = l1 ? d : -1;
    }
    nt += mid;
    append_seg(cutj, dj, Lj, nt);
    return out;
}

// Translate whole codons, dropping any trailing partial one -- `reference.translate`'s contract.
std::string translate_nt(const std::string& nt) {
    auto code = [](char c) { return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3; };
    std::string aa;
    aa.reserve(nt.size() / 3);
    for (size_t p = 0; p + 2 < nt.size(); p += 3)
        aa.push_back(CODON[code(nt[p]) * 16 + code(nt[p + 1]) * 4 + code(nt[p + 2])]);
    return aa;
}

bool is_productive(const std::string& nt, const std::string& aa) {
    return !nt.empty() && nt.size() % 3 == 0 && aa.find('*') == std::string::npos;
}

}  // namespace

GenBatch generate_batch(const PackedModel& m, int n, uint64_t seed, bool productive_only,
                        int nthreads) {
    if (n < 0) throw std::invalid_argument("generate_batch: n must be >= 0");
    const GenPrep g = build_gen_prep(m);
    if (g.v.empty() || (!m.vdj && n > 0 && g.j_given_v.empty()))
        throw std::invalid_argument("generate_batch: the model has no usable V genes");

    struct One { std::string nt, aa; Drawn genes; bool productive; };
    auto per = run_batch(static_cast<size_t>(n), nthreads, [&](size_t slot) {
        Rng r(seed, slot + 1);                 // the slot index alone fixes this draw
        One o;
        for (int attempt = 0; attempt < kDrawAttempts; ++attempt) {
            o.genes = draw_one(r, m, g, o.nt);
            o.aa = translate_nt(o.nt);
            o.productive = is_productive(o.nt, o.aa);
            if (!productive_only || o.productive) return o;
        }
        o.nt.clear();                          // signals the budget ran out; reported below
        return o;
    });

    GenBatch out;
    out.nt.reserve(per.size());
    out.aa.reserve(per.size());
    for (auto* c : {&out.v, &out.d, &out.d2, &out.j}) c->reserve(per.size());
    out.productive.reserve(per.size());
    for (const One& o : per) {
        if (productive_only && o.nt.empty())
            throw std::runtime_error("generate_batch: a slot exhausted its attempt budget "
                                     "(are productive draws this rare?)");
        out.nt.push_back(o.nt);
        out.aa.push_back(o.aa);
        out.v.push_back(o.genes.v);
        out.d.push_back(o.genes.d);
        out.d2.push_back(o.genes.d2);
        out.j.push_back(o.genes.j);
        out.productive.push_back(static_cast<int8_t>(o.productive));
    }
    return out;
}

AlignVotes align_init_votes(const PackedModel& m, const std::vector<std::string>& seqs,
                            bool joint) {
    AlignVotes out;
    out.v.assign(static_cast<size_t>(m.nV()), 0.0);
    out.j.assign(static_cast<size_t>(m.nJ()), 0.0);
    if (joint) out.vj.assign(static_cast<size_t>(m.nV()) * m.nJ(), 0.0);

    std::vector<int8_t> read;
    std::vector<int> best_v, best_j;
    for (const std::string& raw : seqs) {
        read.resize(raw.size());
        for (size_t p = 0; p < raw.size(); ++p) {
            switch (raw[p]) {
                case 'A': case 'a': read[p] = 0; break;
                case 'C': case 'c': read[p] = 1; break;
                case 'G': case 'g': read[p] = 2; break;
                case 'T': case 't': read[p] = 3; break;
                default: read[p] = -1; break;   // matches no germline base, as in the reference
            }
        }
        const int slen = static_cast<int>(read.size());
        best_v.clear();
        int top = -1;
        for (int v : m.func_v) {
            const int sc = common_prefix(m.cut_v[v], read.data(), slen);
            if (sc > top) { top = sc; best_v.clear(); best_v.push_back(v); }
            else if (sc == top) best_v.push_back(v);
        }
        best_j.clear();
        top = -1;
        for (int j : m.func_j) {
            const int sc = common_suffix(m.cut_j[j], read.data(), slen);
            if (sc > top) { top = sc; best_j.clear(); best_j.push_back(j); }
            else if (sc == top) best_j.push_back(j);
        }
        if (best_v.empty() || best_j.empty()) continue;
        const double wv = 1.0 / static_cast<double>(best_v.size());
        const double wj = 1.0 / static_cast<double>(best_j.size());
        for (int v : best_v) out.v[v] += wv;
        for (int j : best_j) out.j[j] += wj;
        if (joint)
            for (int v : best_v)
                for (int j : best_j) out.vj[static_cast<size_t>(v) * m.nJ() + j] += wv * wj;
    }
    return out;
}

}  // namespace vdjtools
