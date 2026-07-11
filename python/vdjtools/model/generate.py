"""Sample recombined sequences from a :class:`~vdjtools.model.model.Model` (OLGA-style generation).

Ancestral sampling over the model's Bayes net: pick genes, deletions, insertion lengths and
non-templated nt (via the dinucleotide Markov chain), assemble the CDR3, translate. Returns a
polars DataFrame. This is the reference sampler; the native port is Phase 1f.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import polars as pl

from .model import Model
from .reference import translate

_NT = "ACGT"


def _cum(values: np.ndarray, probs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(values, cumulative) for a categorical, dropping zero-prob atoms; renormalized."""
    m = probs > 0
    v, p = values[m], probs[m]
    return v, np.cumsum(p / p.sum())


@dataclass(slots=True)
class _GenPrep:
    chain_type: str
    v: tuple[np.ndarray, np.ndarray]                 # (alleles, cum) P(V)
    j_marg: tuple[np.ndarray, np.ndarray] | None     # VDJ: (alleles, cum) P(J)
    j_given_v: dict[str, tuple]                       # VJ: v -> (alleles, cum) P(J|V)
    d_given_j: dict[str, tuple]                        # VDJ: j -> (alleles, cum) P(D|J)
    cut: dict[str, dict[str, str]]                    # "v"/"j"/"d" -> {allele: cut segment}
    delv: dict[str, tuple]                             # v -> (ndel, cum)
    delj: dict[str, tuple]                             # j -> (ndel, cum)
    deld: dict[str, tuple]                             # d -> (idx into pairs, cum), pairs[i]=(n5,n3)
    deld_pairs: dict[str, np.ndarray]
    maxpal: dict[str, int]
    ins: dict[str, tuple]                             # junction -> (lengths, cum)
    R: dict[str, np.ndarray]
    bias: dict[str, np.ndarray]


def _dinucl_matrix(df) -> np.ndarray:
    R = np.zeros((4, 4))
    for frm, to, p in df.select(["from_nt", "to_nt", "p"]).iter_rows():
        R[to, frm] = p
    return R


def _steady_state(R: np.ndarray) -> np.ndarray:
    w, v = np.linalg.eig(R)
    x = np.real(v[:, np.argmin(np.abs(w - 1.0))])
    return x / x.sum()


def prepare_generation(model: Model) -> _GenPrep:
    """Precompute cumulative distributions for fast ancestral sampling."""
    t = model.tables
    vdj = model.chain_type == "VDJ"

    def cut_map(seg):
        g = model.genomic[f"genes_{seg}"]
        return dict(zip(g[f"{seg}_allele"], g["cut_segment"]))

    def usable(seg):
        # Only sample functional alleles (non-empty germline + a defined deletion distribution).
        # OLGA keeps non-functional genes in the marginal with residual mass but can't recombine
        # them; drawing them would emit sequences with zero generation probability.
        g = model.genomic[f"genes_{seg}"]
        return set(g.filter(g["functional"])[f"{seg}_allele"].to_list())

    uv, uj = usable("v"), usable("j")
    ud = usable("d") if vdj else set()

    def _choice(df, allele_col, keep):
        sub = df.filter(df[allele_col].is_in(list(keep)))
        return _cum(sub[allele_col].to_numpy(), sub["p"].to_numpy())

    v = _choice(t["v_choice"], "v_allele", uv)

    j_marg, j_given_v, d_given_j = None, {}, {}
    if vdj:
        j_marg = _choice(t["j_choice"], "j_allele", uj)
        for j, sub in t["d_gene"].filter(t["d_gene"]["d_allele"].is_in(list(ud))).group_by("j_allele"):
            key = j[0] if isinstance(j, tuple) else j
            if key in uj:
                d_given_j[key] = _cum(sub["d_allele"].to_numpy(), sub["p"].to_numpy())
    else:
        jgv = t["j_choice"].filter(t["j_choice"]["j_allele"].is_in(list(uj)))
        for vv, sub in jgv.group_by("v_allele"):
            key = vv[0] if isinstance(vv, tuple) else vv
            if key in uv:
                j_given_v[key] = _cum(sub["j_allele"].to_numpy(), sub["p"].to_numpy())

    def del_map(tab, allele_col):
        out = {}
        for a, sub in tab.group_by(allele_col):
            key = a[0] if isinstance(a, tuple) else a
            out[key] = _cum(sub["ndel"].to_numpy(), sub["p"].to_numpy())
        return out

    delv = del_map(t["v_3_del"], "v_allele")
    delj = del_map(t["j_5_del"], "j_allele")
    deld, deld_pairs = {}, {}
    if vdj:
        for d, sub in t["d_del"].group_by("d_allele"):
            key = d[0] if isinstance(d, tuple) else d
            pairs = np.stack([sub["ndel5"].to_numpy(), sub["ndel3"].to_numpy()], axis=1)
            idx, cum = _cum(np.arange(len(pairs)), sub["p"].to_numpy())
            deld[key] = (idx, cum)
            deld_pairs[key] = pairs

    ins, R, bias = {}, {}, {}
    for junc in (("vd", "dj") if vdj else ("vj",)):
        it = t[f"{junc}_ins"].sort("length")
        ins[junc] = _cum(it["length"].to_numpy(), it["p"].to_numpy())
        R[junc] = _dinucl_matrix(t[f"{junc}_dinucl"])
        bias[junc] = _steady_state(R[junc])

    return _GenPrep(
        chain_type=model.chain_type, v=v, j_marg=j_marg, j_given_v=j_given_v, d_given_j=d_given_j,
        cut={"v": cut_map("v"), "j": cut_map("j"), **({"d": cut_map("d")} if vdj else {})},
        delv=delv, delj=delj, deld=deld, deld_pairs=deld_pairs,
        maxpal=model.manifest.palindrome_max, ins=ins, R=R, bias=bias,
    )


def _pick(rng, prep_pair):
    values, cum = prep_pair
    return values[np.searchsorted(cum, rng.random())]


def _insert(rng, length: int, R: np.ndarray, bias: np.ndarray, *, from_right: bool) -> str:
    if length == 0:
        return ""
    out = [0] * length
    if from_right:
        out[-1] = np.searchsorted(np.cumsum(bias), rng.random())
        for k in range(length - 2, -1, -1):
            out[k] = np.searchsorted(np.cumsum(R[:, out[k + 1]]), rng.random())
    else:
        out[0] = np.searchsorted(np.cumsum(bias), rng.random())
        for k in range(1, length):
            out[k] = np.searchsorted(np.cumsum(R[:, out[k - 1]]), rng.random())
    return "".join(_NT[i] for i in out)


def _draw(prep: _GenPrep, rng) -> tuple[str, str, str, str]:
    """One recombination draw -> (cdr3_nt, v, d, j)."""
    vdj = prep.chain_type == "VDJ"
    v = _pick(rng, prep.v)
    if vdj:
        j = _pick(rng, prep.j_marg)
        d = _pick(rng, prep.d_given_j[j])
    else:
        j = _pick(rng, prep.j_given_v[v])
        d = ""

    cutv = prep.cut["v"][v]
    len_v = max(0, len(cutv) - (_pick(rng, prep.delv[v]) + prep.maxpal["v_3"]))
    v_contrib = cutv[:len_v]
    cutj = prep.cut["j"][j]
    dj = _pick(rng, prep.delj[j]) + prep.maxpal["j_5"]
    j_contrib = cutj[dj:]

    if vdj:
        idx, _cumd = prep.deld[d]
        n5, n3 = prep.deld_pairs[d][_pick(rng, prep.deld[d])]
        cutd = prep.cut["d"][d]
        start = n5 + prep.maxpal["d_5"]
        end = len(cutd) - (n3 + prep.maxpal["d_3"])
        d_contrib = cutd[start:end] if start < end else ""
        ins_vd = _insert(rng, int(_pick(rng, prep.ins["vd"])), prep.R["vd"], prep.bias["vd"], from_right=False)
        ins_dj = _insert(rng, int(_pick(rng, prep.ins["dj"])), prep.R["dj"], prep.bias["dj"], from_right=True)
        cdr3 = v_contrib + ins_vd + d_contrib + ins_dj + j_contrib
    else:
        ins_vj = _insert(rng, int(_pick(rng, prep.ins["vj"])), prep.R["vj"], prep.bias["vj"], from_right=False)
        cdr3 = v_contrib + ins_vj + j_contrib
    return cdr3, v, d, j


def generate(model: Model, n: int, *, seed: int | None = None, productive_only: bool = False) -> pl.DataFrame:
    """Sample ``n`` recombined CDR3s from the model.

    Args:
        model: The generative model.
        n: Number of sequences to return.
        seed: RNG seed for reproducibility.
        productive_only: If True, reject out-of-frame / stop-codon draws and keep sampling.

    Returns:
        DataFrame with ``cdr3_nt, cdr3_aa, v_call, d_call, j_call, productive``.
    """
    prep = prepare_generation(model)
    rng = np.random.default_rng(seed)
    rows = {k: [] for k in ("cdr3_nt", "cdr3_aa", "v_call", "d_call", "j_call", "productive")}
    got = 0
    guard = 0
    max_guard = n * 10000 + 1000
    while got < n:
        guard += 1
        if guard > max_guard:
            raise RuntimeError("generation exceeded attempt budget (productive draws too rare?)")
        cdr3, v, d, j = _draw(prep, rng)
        aa = translate(cdr3)
        productive = len(cdr3) % 3 == 0 and "*" not in aa and len(cdr3) > 0
        if productive_only and not productive:
            continue
        rows["cdr3_nt"].append(cdr3)
        rows["cdr3_aa"].append(aa)
        rows["v_call"].append(v)
        rows["d_call"].append(d if d else None)
        rows["j_call"].append(j)
        rows["productive"].append(productive)
        got += 1
    return pl.DataFrame(rows)
