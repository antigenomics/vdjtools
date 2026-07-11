"""Generation probability (Pgen) of a CDR3 from a :class:`~vdjtools.model.model.Model`.

This is the reference **nucleotide** implementation: a transparent sum over every
recombination scenario (V/J/D choice × deletions × insertion partition) that could produce
the observed CDR3 nt, reading straight from the model's polars tables (no OLGA at runtime).
It reproduces OLGA's ``compute_nt_CDR3_pgen`` to numerical tolerance (see the tests).

It is correctness-first, not speed-first — the fast amino-acid transfer-matrix DP is the job
of the native ``_core`` port (Phase 1f). The nt scenario sum here is exactly the quantity the
EM E-step needs (the bootstrap training data is out-of-frame *nucleotide* reads).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .model import Model
from .reference import _CODON_TABLE

_NT2NUM = {"A": 0, "C": 1, "G": 2, "T": 3}
_NUM2NT = "ACGT"
# aa of a codon given three nucleotide integers (A,C,G,T = 0..3).
_CODON_AA = {
    (a, b, c): _CODON_TABLE[_NUM2NT[a] + _NUM2NT[b] + _NUM2NT[c]]
    for a in range(4)
    for b in range(4)
    for c in range(4)
}


def _steady_state(R: np.ndarray) -> np.ndarray:
    """Stationary distribution of a column-stochastic Markov matrix ``R[next, prev]``."""
    w, v = np.linalg.eig(R)
    x = np.real(v[:, np.argmin(np.abs(w - 1.0))])
    return x / x.sum()


@dataclass(slots=True)
class _Prepared:
    """Model unpacked into plain dicts/arrays for the Pgen sum (built once, reused)."""

    chain_type: str
    functional_v: list[str]
    functional_j: list[str]
    functional_d: list[str]
    cut: dict[str, dict[str, str]]          # "v"/"j"/"d" -> {allele: palindrome-extended segment}
    p_v: dict[str, float]                   # P(V)
    p_j: dict[tuple, float] | dict[str, float]  # VDJ: {j: P(J)}; VJ: {(v, j): P(J|V)}
    p_d_given_j: dict[tuple, float]         # {(j, d): P(D|J)}  (VDJ only)
    p_del: dict[str, dict]                  # "v"/"j" -> {(allele, ndel): p}; "d" -> {(allele, n5, n3): p}
    maxpal: dict[str, int]                  # "v_3"/"j_5"/"d_5"/"d_3" -> palindrome max
    p_ins: dict[str, np.ndarray]            # "vd"/"dj"/"vj"/"dd" -> P(length)
    R: dict[str, np.ndarray]                # "vd"/"dj"/"vj"/"dd" -> 4×4 Markov
    bias: dict[str, np.ndarray]             # "vd"/"dj"/"vj"/"dd" -> steady-state first-nt bias
    # D-D (tandem) extension — populated only when the model declares an ``n_d`` mass at 2.
    p_nd: dict[int, float] = field(default_factory=lambda: {1: 1.0})   # P(n_D)
    p_d2_given_d1: dict[tuple, float] = field(default_factory=dict)     # {(d1, d2): P(D2|D1)}
    p_del_d2: dict = field(default_factory=dict)                        # {(d2, n5, n3): p}


def _dinucl_matrix(df) -> np.ndarray:
    R = np.zeros((4, 4))
    for frm, to, p in df.select(["from_nt", "to_nt", "p"]).iter_rows():
        R[to, frm] = p  # R[next, prev]
    return R


def prepare(model: Model) -> _Prepared:
    """Unpack a model's polars tables into the lookup structures the Pgen sum uses."""
    t = model.tables
    vdj = model.chain_type == "VDJ"

    def cut_map(seg: str) -> dict[str, str]:
        g = model.genomic[f"genes_{seg}"]
        return dict(zip(g[f"{seg}_allele"], g["cut_segment"]))

    def functional(seg: str) -> list[str]:
        g = model.genomic[f"genes_{seg}"]
        return g.filter(g["functional"])[f"{seg}_allele"].to_list()

    p_v = dict(zip(t["v_choice"]["v_allele"], t["v_choice"]["p"]))
    if vdj:
        p_j = dict(zip(t["j_choice"]["j_allele"], t["j_choice"]["p"]))
        p_d_given_j = {(j, d): p for j, d, p in t["d_gene"].select(["j_allele", "d_allele", "p"]).iter_rows()}
    else:
        p_j = {(v, j): p for v, j, p in t["j_choice"].select(["v_allele", "j_allele", "p"]).iter_rows()}
        p_d_given_j = {}

    p_del: dict[str, dict] = {
        "v": {(a, n): p for a, n, p in t["v_3_del"].select(["v_allele", "ndel", "p"]).iter_rows()},
        "j": {(a, n): p for a, n, p in t["j_5_del"].select(["j_allele", "ndel", "p"]).iter_rows()},
    }
    if vdj:
        p_del["d"] = {
            (a, n5, n3): p
            for a, n5, n3, p in t["d_del"].select(["d_allele", "ndel5", "ndel3", "p"]).iter_rows()
        }

    juncs = ["vd", "dj"] if vdj else ["vj"]
    if vdj and "dd_ins" in t:  # tandem-D middle junction
        juncs.append("dd")
    p_ins, R, bias = {}, {}, {}
    for junc in juncs:
        p_ins[junc] = t[f"{junc}_ins"].sort("length")["p"].to_numpy()
        R[junc] = _dinucl_matrix(t[f"{junc}_dinucl"])
        bias[junc] = _steady_state(R[junc])

    # D-D tables — present only when the model declares the tandem extension.
    p_nd = {int(n): float(p) for n, p in t["n_d"].select(["n_d", "p"]).iter_rows()} if "n_d" in t else {1: 1.0}
    p_d2_given_d1: dict[tuple, float] = {}
    p_del_d2: dict = {}
    if vdj and "d2_gene" in t:
        p_d2_given_d1 = {(d1, d2): p for d1, d2, p in t["d2_gene"].select(["d_allele", "d2_allele", "p"]).iter_rows()}
        p_del_d2 = {(a, n5, n3): p for a, n5, n3, p in t["d2_del"].select(["d2_allele", "ndel5", "ndel3", "p"]).iter_rows()}
    # A model that places mass on n_D=2 must ship the tandem tables, else _vdj_middle would silently
    # drop that mass (Pgen under-counts). Reject the malformed state at prepare time.
    if p_nd.get(2, 0.0) > 0.0 and not p_d2_given_d1:
        raise ValueError("model has P(n_D=2)>0 but no d2_gene/d2_del/dd tables (malformed tandem model)")

    return _Prepared(
        chain_type=model.chain_type,
        functional_v=functional("v"),
        functional_j=functional("j"),
        functional_d=functional("d") if vdj else [],
        cut={"v": cut_map("v"), "j": cut_map("j"), **({"d": cut_map("d")} if vdj else {})},
        p_v=p_v,
        p_j=p_j,
        p_d_given_j=p_d_given_j,
        p_del=p_del,
        maxpal=model.manifest.palindrome_max,
        p_ins=p_ins,
        R=R,
        bias=bias,
        p_nd=p_nd,
        p_d2_given_d1=p_d2_given_d1,
        p_del_d2=p_del_d2,
    )


def _p_insert(seq: str, p_len: np.ndarray, R: np.ndarray, bias: np.ndarray, *, from_right: bool) -> float:
    """Probability of an N-region ``seq``: P(len) × first-nt bias × Markov chain product.

    ``from_right`` walks the chain from the 3' (J-adjacent) end — the DJ convention.
    """
    n = len(seq)
    if n >= len(p_len):
        return 0.0
    p = p_len[n]
    if p == 0.0 or n == 0:
        return p
    nums = [_NT2NUM[c] for c in seq]
    if from_right:
        p *= bias[nums[-1]]
        for k in range(n - 2, -1, -1):
            p *= R[nums[k], nums[k + 1]]
    else:
        p *= bias[nums[0]]
        for k in range(1, n):
            p *= R[nums[k], nums[k - 1]]
    return p


def _common_prefix(a: str, b: str) -> int:
    n = min(len(a), len(b))
    i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def _common_suffix(a: str, b: str) -> int:
    n = min(len(a), len(b))
    i = 0
    while i < n and a[-1 - i] == b[-1 - i]:
        i += 1
    return i


def _v_options(prep: _Prepared, v: str, s: str) -> list[tuple[int, int, float]]:
    """(len_v, ndel, P(delV|V)) for every 3' trim of V whose germline prefixes ``s``.

    ``cut[:len_v] == s[:len_v]`` holds exactly for ``len_v <= `` the common-prefix length.
    """
    cut = prep.cut["v"][v]
    maxp = prep.maxpal["v_3"]
    out = []
    # len_v >= 1: the V gene must contribute at least one nt to the CDR3 (a valid V(D)J
    # rearrangement is never fully V-deleted). Matches OLGA, which excludes len_v == 0.
    for len_v in range(_common_prefix(cut, s), 0, -1):
        p = prep.p_del["v"].get((v, len(cut) - len_v - maxp), 0.0)
        if p:
            out.append((len_v, len(cut) - len_v - maxp, p))
    return out


def _j_options(prep: _Prepared, j: str, s: str) -> list[tuple[int, int, float]]:
    """(len_j, ndel, P(delJ|J)) for every 5' trim of J whose germline suffixes ``s``."""
    cut = prep.cut["j"][j]
    maxp = prep.maxpal["j_5"]
    out = []
    # len_j >= 1: the J gene must contribute at least one nt to the CDR3 (matches OLGA).
    for len_j in range(_common_suffix(cut, s), 0, -1):
        p = prep.p_del["j"].get((j, len(cut) - len_j - maxp), 0.0)
        if p:
            out.append((len_j, len(cut) - len_j - maxp, p))
    return out


def _v_candidates(prep: _Prepared, s: str, v_mask):
    """[(V, P(V), v_opts)] restricted to genes whose germline can prefix ``s`` (the pruning step)."""
    out = []
    for v in v_mask:
        pv = prep.p_v.get(v, 0.0)
        if pv == 0.0:
            continue
        opts = _v_options(prep, v, s)
        if opts:
            out.append((v, pv, opts))
    return out


def _j_candidates(prep: _Prepared, s: str, j_mask):
    """[(J, j_opts)] restricted to genes whose germline can suffix ``s``."""
    out = []
    for j in j_mask:
        opts = _j_options(prep, j, s)
        if opts:
            out.append((j, opts))
    return out


def _d_middle(prep: _Prepared, j: str, middle: str) -> float:
    """Σ over D, its 5'/3' trims and its position of P(D|J)·P(delD|D)·Pins(VD)·Pins(DJ)."""
    total = 0.0
    m = len(middle)
    maxdl, maxdr = prep.maxpal["d_5"], prep.maxpal["d_3"]
    pins_vd, R_vd, b_vd = prep.p_ins["vd"], prep.R["vd"], prep.bias["vd"]
    pins_dj, R_dj, b_dj = prep.p_ins["dj"], prep.R["dj"], prep.bias["dj"]
    for d in prep.functional_d:
        pdj = prep.p_d_given_j.get((j, d), 0.0)
        if pdj == 0.0:
            continue
        cut = prep.cut["d"][d]
        acc = 0.0
        for idx5 in range(len(cut) + 1):
            for idx3 in range(len(cut) - idx5 + 1):
                d_contrib = cut[idx5:len(cut) - idx3]
                pdel = prep.p_del["d"].get((d, idx5 - maxdl, idx3 - maxdr), 0.0)
                if pdel == 0.0:
                    continue
                ld = len(d_contrib)
                for pos in range(0, m - ld + 1):
                    if middle[pos:pos + ld] != d_contrib:
                        continue
                    w = _p_insert(middle[:pos], pins_vd, R_vd, b_vd, from_right=False)
                    if w == 0.0:
                        continue
                    w *= _p_insert(middle[pos + ld:], pins_dj, R_dj, b_dj, from_right=True)
                    acc += pdel * w
        total += pdj * acc
    return total


def _dd_middle(prep: _Prepared, j: str, middle: str) -> float:
    """Σ over the *tandem* (two-D) scenarios producing ``middle`` = [insVD] D1 [insDD] D2 [insDJ] J.

    Both D segments must contribute ≥1 nt (a genuine tandem — a fully-trimmed second D collapses
    to the single-D case, already summed by :func:`_d_middle`, so excluding it keeps the n_D=1 and
    n_D=2 scenario sets disjoint and identifiable). ``P(D2|D1)`` carries the genomic-order mask.

    The DD-junction insertion is read 5'→3' (``from_right=False``), like the VD junction; the future
    native port and the tandem generator must use the same orientation for its dinucleotide Markov.
    """
    m = len(middle)
    maxdl, maxdr = prep.maxpal["d_5"], prep.maxpal["d_3"]
    pins_vd, R_vd, b_vd = prep.p_ins["vd"], prep.R["vd"], prep.bias["vd"]
    pins_dd, R_dd, b_dd = prep.p_ins["dd"], prep.R["dd"], prep.bias["dd"]
    pins_dj, R_dj, b_dj = prep.p_ins["dj"], prep.R["dj"], prep.bias["dj"]

    def d_placements(d, pdel_tbl):
        """[(contrib_str, pdel)] for each 5'/3' trim of D that leaves ≥1 nt."""
        cut = prep.cut["d"][d]
        out = []
        for idx5 in range(len(cut) + 1):
            for idx3 in range(len(cut) - idx5 + 1):
                contrib = cut[idx5:len(cut) - idx3]
                if not contrib:  # tandem requires each D to contribute ≥1 nt
                    continue
                pdel = pdel_tbl.get((d, idx5 - maxdl, idx3 - maxdr), 0.0)
                if pdel:
                    out.append((contrib, pdel))
        return out

    total = 0.0
    for d1 in prep.functional_d:
        pd1 = prep.p_d_given_j.get((j, d1), 0.0)
        if pd1 == 0.0:
            continue
        pl1 = d_placements(d1, prep.p_del["d"])
        if not pl1:
            continue
        for d2 in prep.functional_d:
            pd2 = prep.p_d2_given_d1.get((d1, d2), 0.0)
            if pd2 == 0.0:
                continue
            pl2 = d_placements(d2, prep.p_del_d2)
            if not pl2:
                continue
            w_gene = pd1 * pd2
            for contrib1, pdel1 in pl1:
                ld1 = len(contrib1)
                for pos1 in range(0, m - ld1):  # leave room for D2 (≥1 nt) after
                    if middle[pos1:pos1 + ld1] != contrib1:
                        continue
                    left = _p_insert(middle[:pos1], pins_vd, R_vd, b_vd, from_right=False)
                    if left == 0.0:
                        continue
                    for contrib2, pdel2 in pl2:
                        ld2 = len(contrib2)
                        for pos2 in range(pos1 + ld1, m - ld2 + 1):
                            if middle[pos2:pos2 + ld2] != contrib2:
                                continue
                            mid = _p_insert(middle[pos1 + ld1:pos2], pins_dd, R_dd, b_dd, from_right=False)
                            if mid == 0.0:
                                continue
                            right = _p_insert(middle[pos2 + ld2:], pins_dj, R_dj, b_dj, from_right=True)
                            total += w_gene * pdel1 * pdel2 * left * mid * right
    return total


def _vdj_middle(prep: _Prepared, j: str, middle: str) -> float:
    """P(middle) mixed over the D-count prior: P(n_D=1)·single-D + P(n_D=2)·tandem-D."""
    p1 = prep.p_nd.get(1, 0.0) + prep.p_nd.get(0, 0.0)  # 0-D folds into 1-D (fully-trimmed D)
    total = p1 * _d_middle(prep, j, middle) if p1 else 0.0
    p2 = prep.p_nd.get(2, 0.0)
    if p2 and prep.p_d2_given_d1:
        total += p2 * _dd_middle(prep, j, middle)
    return total


def _aa_dp(aa: str, template: list[int], specs: list) -> float:
    """Left-to-right codon-constrained sum over N-region nt, Markov-weighted.

    ``template`` is the CDR3 as nt integers with ``-1`` at insertion positions. ``specs[i]`` is
    ``None`` for a germline (fixed) position, else ``(kind, R, bias, is_first, is_last)`` where
    ``kind`` is ``"L"`` (VD/VJ, read 5'->3': first nt uses bias, then ``R[cur, prev]``) or ``"R"``
    (DJ, read 3'->5': ``R[prev, cur]`` transitions and the block's *last* nt uses bias). Codons are
    checked as they complete; mistranslations prune. Excludes the P(insLen) factors.
    """
    dp = {(-1, -1): 1.0}  # state = (nt[i-1], nt[i-2])
    for i, fixed in enumerate(template):
        spec = specs[i]
        codon_end = i % 3 == 2
        aa_i = aa[i // 3]
        ndp: dict = {}
        for (p1, p2), w in dp.items():
            for nt in (range(4) if fixed < 0 else (fixed,)):
                ww = w
                if spec is not None:
                    kind, R, bias, is_first, is_last = spec
                    if kind == "L":
                        ww *= bias[nt] if is_first else R[nt, p1]
                    else:  # "R" — DJ, read from the 3' end
                        if not is_first:
                            ww *= R[p1, nt]
                        if is_last:
                            ww *= bias[nt]
                if codon_end and _CODON_AA[(p2, p1, nt)] != aa_i:
                    continue
                key = (nt, p1)
                ndp[key] = ndp.get(key, 0.0) + ww
        if not ndp:
            return 0.0
        dp = ndp
    return sum(dp.values())


def _pgen_aa_vj(prep: _Prepared, aa: str, v: str | None, j: str | None) -> float:
    N = 3 * len(aa)
    pins, R, bias = prep.p_ins["vj"], prep.R["vj"], prep.bias["vj"]
    total = 0.0
    for V in ([v] if v else prep.functional_v):
        pv = prep.p_v.get(V, 0.0)
        if pv == 0.0:
            continue
        cutv = prep.cut["v"][V]
        maxpv = prep.maxpal["v_3"]
        for J in ([j] if j else prep.functional_j):
            pj = prep.p_j.get((V, J), 0.0)
            if pj == 0.0:
                continue
            cutj = prep.cut["j"][J]
            maxpj = prep.maxpal["j_5"]
            for len_v in range(1, min(len(cutv), N) + 1):  # V contributes >= 1 nt (matches OLGA)
                pdv = prep.p_del["v"].get((V, len(cutv) - len_v - maxpv), 0.0)
                if pdv == 0.0:
                    continue
                gv = [_NT2NUM[c] for c in cutv[:len_v]]
                for len_j in range(1, min(len(cutj), N - len_v) + 1):  # J contributes >= 1 nt
                    idxj = len(cutj) - len_j
                    pdj = prep.p_del["j"].get((J, idxj - maxpj), 0.0)
                    if pdj == 0.0:
                        continue
                    ins_len = N - len_v - len_j
                    if ins_len < 0 or ins_len >= len(pins) or pins[ins_len] == 0.0:
                        continue
                    template = gv + [-1] * ins_len + [_NT2NUM[c] for c in cutj[idxj:]]
                    specs = [None] * N
                    for p in range(len_v, len_v + ins_len):
                        specs[p] = ("L", R, bias, p == len_v, p == len_v + ins_len - 1)
                    w = _aa_dp(aa, template, specs)
                    if w > 0.0:
                        total += pv * pj * pdv * pdj * pins[ins_len] * w
    return total


def _aa_germline_prefix_ok(cut: str, aa: str, len_v: int) -> bool:
    """The V germline's *full* codons (within [0, len_v)) must translate to the aa prefix."""
    full = len_v // 3
    return all(_CODON_TABLE[cut[3 * c:3 * c + 3]] == aa[c] for c in range(full))


def _pgen_aa_vdj(prep: _Prepared, aa: str, v: str | None, j: str | None) -> float:
    N = 3 * len(aa)
    pins_vd, R_vd, b_vd = prep.p_ins["vd"], prep.R["vd"], prep.bias["vd"]
    pins_dj, R_dj, b_dj = prep.p_ins["dj"], prep.R["dj"], prep.bias["dj"]
    maxdl, maxdr = prep.maxpal["d_5"], prep.maxpal["d_3"]
    total = 0.0
    for V in ([v] if v else prep.functional_v):
        pv = prep.p_v.get(V, 0.0)
        if pv == 0.0:
            continue
        cutv, maxpv = prep.cut["v"][V], prep.maxpal["v_3"]
        for J in ([j] if j else prep.functional_j):
            pj = prep.p_j.get(J, 0.0)
            if pj == 0.0:
                continue
            cutj, maxpj = prep.cut["j"][J], prep.maxpal["j_5"]
            for len_v in range(1, min(len(cutv), N) + 1):
                pdv = prep.p_del["v"].get((V, len(cutv) - len_v - maxpv), 0.0)
                if pdv == 0.0 or not _aa_germline_prefix_ok(cutv, aa, len_v):
                    continue
                gv = [_NT2NUM[c] for c in cutv[:len_v]]
                for len_j in range(1, min(len(cutj), N - len_v) + 1):
                    idxj = len(cutj) - len_j
                    pdj_del = prep.p_del["j"].get((J, idxj - maxpj), 0.0)
                    if pdj_del == 0.0:
                        continue
                    gj_seq = cutj[idxj:]
                    # J full codons (aligned to the CDR3 frame) must translate to the aa suffix.
                    if not all(
                        _CODON_TABLE[gj_seq[3 * c - (N - len_j):][:3]] == aa[c]
                        for c in range((N - len_j + 2) // 3, N // 3)
                        if 3 * c >= N - len_j
                    ):
                        continue
                    gj = [_NT2NUM[c] for c in gj_seq]
                    args = (aa, gv, gj, len_v, len_j, N, J, prep,
                            pins_vd, R_vd, b_vd, pins_dj, R_dj, b_dj, maxdl, maxdr)
                    p1 = prep.p_nd.get(1, 0.0) + prep.p_nd.get(0, 0.0)  # 0-D folds into 1-D
                    mid_aa = p1 * _d_aa_middle(*args) if p1 else 0.0
                    p2 = prep.p_nd.get(2, 0.0)
                    if p2 and prep.p_d2_given_d1:  # tandem-D contribution
                        mid_aa += p2 * _dd_aa_middle(*args)
                    total += pv * pj * pdv * pdj_del * mid_aa
    return total


def _d_codons_ok(dc_nt: list[int], pos: int, aa: str) -> bool:
    """D-germline full codons (aligned to the CDR3 frame) translate to aa — cheap pre-prune."""
    ld = len(dc_nt)
    c = -(-pos // 3)  # first codon starting at/after pos
    while 3 * c + 2 < pos + ld:
        if _CODON_AA[(dc_nt[3 * c - pos], dc_nt[3 * c + 1 - pos], dc_nt[3 * c + 2 - pos])] != aa[c]:
            return False
        c += 1
    return True


def _d_aa_middle(aa, gv, gj, len_v, len_j, N, J, prep, pins_vd, R_vd, b_vd,
                 pins_dj, R_dj, b_dj, maxdl, maxdr) -> float:
    """Σ over D, its 5'/3' trims and its position of P(D|J)·P(delD|D) × the N-region codon DP."""
    out = 0.0
    right = N - len_j  # DJ insertion + J start here
    for D in prep.functional_d:
        pdgene = prep.p_d_given_j.get((J, D), 0.0)
        if pdgene == 0.0:
            continue
        cutd = prep.cut["d"][D]
        for idx5 in range(len(cutd) + 1):
            for idx3 in range(len(cutd) - idx5 + 1):
                pdel = prep.p_del["d"].get((D, idx5 - maxdl, idx3 - maxdr), 0.0)
                if pdel == 0.0:
                    continue
                dc = [_NT2NUM[c] for c in cutd[idx5:len(cutd) - idx3]]
                ld = len(dc)
                for pos in range(len_v, right - ld + 1):
                    lvd, ldj = pos - len_v, right - pos - ld
                    if lvd >= len(pins_vd) or pins_vd[lvd] == 0.0:
                        continue
                    if ldj >= len(pins_dj) or pins_dj[ldj] == 0.0:
                        continue
                    if not _d_codons_ok(dc, pos, aa):
                        continue
                    template = gv + [-1] * lvd + dc + [-1] * ldj + gj
                    specs = [None] * N
                    for p in range(len_v, pos):
                        specs[p] = ("L", R_vd, b_vd, p == len_v, p == pos - 1)
                    for p in range(pos + ld, right):
                        specs[p] = ("R", R_dj, b_dj, p == pos + ld, p == right - 1)
                    w = _aa_dp(aa, template, specs)
                    if w > 0.0:
                        out += pdgene * pdel * pins_vd[lvd] * pins_dj[ldj] * w
    return out


def _dd_aa_middle(aa, gv, gj, len_v, len_j, N, J, prep, pins_vd, R_vd, b_vd,
                  pins_dj, R_dj, b_dj, maxdl, maxdr) -> float:
    """Σ over the *tandem* (two-D) scenarios × the N-region codon DP — the aa analogue of
    :func:`_dd_middle`. Middle = [insVD] D1 [insDD] D2 [insDJ]; each D contributes ≥1 nt. The
    DD-junction insertion is read 5'→3' (kind ``"L"``, like VD), matching the nt reference."""
    pins_dd, R_dd, b_dd = prep.p_ins["dd"], prep.R["dd"], prep.bias["dd"]
    right = N - len_j  # DJ insertion + J start here
    out = 0.0
    for d1 in prep.functional_d:
        pd1 = prep.p_d_given_j.get((J, d1), 0.0)
        if pd1 == 0.0:
            continue
        cut1 = prep.cut["d"][d1]
        for i5 in range(len(cut1) + 1):
            for i3 in range(len(cut1) - i5 + 1):
                pdel1 = prep.p_del["d"].get((d1, i5 - maxdl, i3 - maxdr), 0.0)
                if pdel1 == 0.0:
                    continue
                dc1 = [_NT2NUM[c] for c in cut1[i5:len(cut1) - i3]]
                ld1 = len(dc1)
                if ld1 < 1:
                    continue
                for pos1 in range(len_v, right - ld1):  # leave ≥1 nt (and a slot) for D2 after
                    lvd = pos1 - len_v
                    if lvd >= len(pins_vd) or pins_vd[lvd] == 0.0:
                        continue
                    if not _d_codons_ok(dc1, pos1, aa):
                        continue
                    for d2 in prep.functional_d:
                        pd2g = prep.p_d2_given_d1.get((d1, d2), 0.0)
                        if pd2g == 0.0:
                            continue
                        cut2 = prep.cut["d"][d2]
                        for k5 in range(len(cut2) + 1):
                            for k3 in range(len(cut2) - k5 + 1):
                                pdel2 = prep.p_del_d2.get((d2, k5 - maxdl, k3 - maxdr), 0.0)
                                if pdel2 == 0.0:
                                    continue
                                dc2 = [_NT2NUM[c] for c in cut2[k5:len(cut2) - k3]]
                                ld2 = len(dc2)
                                if ld2 < 1:
                                    continue
                                for pos2 in range(pos1 + ld1, right - ld2 + 1):
                                    ldd, ldj = pos2 - (pos1 + ld1), right - (pos2 + ld2)
                                    if ldd >= len(pins_dd) or pins_dd[ldd] == 0.0:
                                        continue
                                    if ldj >= len(pins_dj) or pins_dj[ldj] == 0.0:
                                        continue
                                    if not _d_codons_ok(dc2, pos2, aa):
                                        continue
                                    template = gv + [-1] * lvd + dc1 + [-1] * ldd + dc2 + [-1] * ldj + gj
                                    specs = [None] * N
                                    for p in range(len_v, pos1):
                                        specs[p] = ("L", R_vd, b_vd, p == len_v, p == pos1 - 1)
                                    for p in range(pos1 + ld1, pos2):
                                        specs[p] = ("L", R_dd, b_dd, p == pos1 + ld1, p == pos2 - 1)
                                    for p in range(pos2 + ld2, right):
                                        specs[p] = ("R", R_dj, b_dj, p == pos2 + ld2, p == right - 1)
                                    w = _aa_dp(aa, template, specs)
                                    if w > 0.0:
                                        out += (pd1 * pd2g * pdel1 * pdel2
                                                * pins_vd[lvd] * pins_dd[ldd] * pins_dj[ldj] * w)
    return out


def pgen_aa(prep: _Prepared, cdr3_aa: str, v: str | None = None, j: str | None = None) -> float:
    """Generation probability of an amino-acid CDR3, marginalizing over synonymous codons.

    Args:
        prep: A :func:`prepare` -d model.
        cdr3_aa: The CDR3 amino-acid sequence (conserved Cys → conserved Phe/Trp inclusive).
        v, j: Optional gene names to restrict the sum to.

    Returns:
        Pgen as a float. Sums tandem-D (``n_D=2``) scenarios for a D-D model, mixed with the
        single-D term by ``P(n_D)`` — consistent with :func:`pgen_nt`.
    """
    aa = cdr3_aa.upper()
    if prep.chain_type == "VJ":
        return _pgen_aa_vj(prep, aa, v, j)
    return _pgen_aa_vdj(prep, aa, v, j)


def pgen_nt(prep: _Prepared, cdr3_nt: str, v: str | None = None, j: str | None = None) -> float:
    """Generation probability of a nucleotide CDR3, optionally restricted to a V and/or J.

    Args:
        prep: A :func:`prepare` -d model.
        cdr3_nt: The CDR3 nucleotide sequence (conserved-Cys → conserved-Phe/Trp inclusive).
        v, j: Optional gene names to restrict the sum to (an OLGA-style usage mask of size 1).

    Returns:
        Pgen as a float.
    """
    s = cdr3_nt.upper()
    vdj = prep.chain_type == "VDJ"
    v_cands = _v_candidates(prep, s, [v] if v else prep.functional_v)
    j_cands = _j_candidates(prep, s, [j] if j else prep.functional_j)
    total = 0.0
    for V, pv, v_opts in v_cands:
        for J, j_opts in j_cands:
            pj = prep.p_j.get(J if vdj else (V, J), 0.0)
            if pj == 0.0:
                continue
            for len_v, _nv, p_dv in v_opts:
                for len_j, _nj, p_dj in j_opts:
                    if len_v + len_j > len(s):
                        continue
                    middle = s[len_v:len(s) - len_j]
                    if vdj:
                        inner = _vdj_middle(prep, J, middle)
                    else:
                        inner = _p_insert(middle, prep.p_ins["vj"], prep.R["vj"], prep.bias["vj"], from_right=False)
                    total += pv * pj * p_dv * p_dj * inner
    return total
