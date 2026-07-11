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

from dataclasses import dataclass

import numpy as np

from .model import Model

_NT2NUM = {"A": 0, "C": 1, "G": 2, "T": 3}


def _steady_state(R: np.ndarray) -> np.ndarray:
    """Stationary distribution of a column-stochastic Markov matrix ``R[next, prev]``."""
    w, v = np.linalg.eig(R)
    x = np.real(v[:, np.argmin(np.abs(w - 1.0))])
    return x / x.sum()


def _markov(R: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(R, steady-state bias) from a dinucleotide table's implied 4×4 matrix."""
    return R, _steady_state(R)


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
    p_ins: dict[str, np.ndarray]            # "vd"/"dj"/"vj" -> P(length)
    R: dict[str, np.ndarray]                # "vd"/"dj"/"vj" -> 4×4 Markov
    bias: dict[str, np.ndarray]             # "vd"/"dj"/"vj" -> steady-state first-nt bias


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

    p_ins, R, bias = {}, {}, {}
    for junc in (("vd", "dj") if vdj else ("vj",)):
        p_ins[junc] = t[f"{junc}_ins"].sort("length")["p"].to_numpy()
        R[junc] = _dinucl_matrix(t[f"{junc}_dinucl"])
        bias[junc] = _steady_state(R[junc])

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
    for len_v in range(_common_prefix(cut, s), -1, -1):
        p = prep.p_del["v"].get((v, len(cut) - len_v - maxp), 0.0)
        if p:
            out.append((len_v, len(cut) - len_v - maxp, p))
    return out


def _j_options(prep: _Prepared, j: str, s: str) -> list[tuple[int, int, float]]:
    """(len_j, ndel, P(delJ|J)) for every 5' trim of J whose germline suffixes ``s``."""
    cut = prep.cut["j"][j]
    maxp = prep.maxpal["j_5"]
    out = []
    for len_j in range(_common_suffix(cut, s), -1, -1):
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
                        inner = _d_middle(prep, J, middle)
                    else:
                        inner = _p_insert(middle, prep.p_ins["vj"], prep.R["vj"], prep.bias["vj"], from_right=False)
                    total += pv * pj * p_dv * p_dj * inner
    return total
