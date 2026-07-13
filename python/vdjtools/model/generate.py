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
from .pgen import _dinucl_matrix, _steady_state
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
    # D-D (tandem) generation — populated only for a model with P(n_D=2)>0, else p_nd is None.
    p_nd: tuple | None                                # (n_d values, cum) P(n_D)
    d2_given_d1: dict                                 # d1 -> (d2 alleles, cum) P(D2|D1)
    deld2: dict                                       # d2 -> (idx into pairs, cum)
    deld2_pairs: dict                                 # d2 -> pairs[i]=(n5,n3)


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
    juncs = ["vd", "dj"] if vdj else ["vj"]
    if vdj and "dd_ins" in t:
        juncs.append("dd")
    for junc in juncs:
        it = t[f"{junc}_ins"].sort("length")
        ins[junc] = _cum(it["length"].to_numpy(), it["p"].to_numpy())
        R[junc] = _dinucl_matrix(t[f"{junc}_dinucl"])
        bias[junc] = _steady_state(R[junc])

    # D-D (tandem) generation structures — only when the model declares them.
    p_nd, d2_given_d1, deld2, deld2_pairs = None, {}, {}, {}
    if vdj and "d2_gene" in t:
        ndf = t["n_d"].sort("n_d")
        p_nd = _cum(ndf["n_d"].to_numpy(), ndf["p"].to_numpy())
        d2t = t["d2_gene"].filter(t["d2_gene"]["d2_allele"].is_in(list(ud)))
        for d1, sub in d2t.group_by("d_allele"):
            key = d1[0] if isinstance(d1, tuple) else d1
            if key in ud:
                d2_given_d1[key] = _cum(sub["d2_allele"].to_numpy(), sub["p"].to_numpy())
        for d, sub in t["d2_del"].group_by("d2_allele"):
            key = d[0] if isinstance(d, tuple) else d
            pairs = np.stack([sub["ndel5"].to_numpy(), sub["ndel3"].to_numpy()], axis=1)
            idx, cum = _cum(np.arange(len(pairs)), sub["p"].to_numpy())
            deld2[key] = (idx, cum)
            deld2_pairs[key] = pairs

    return _GenPrep(
        chain_type=model.chain_type, v=v, j_marg=j_marg, j_given_v=j_given_v, d_given_j=d_given_j,
        cut={"v": cut_map("v"), "j": cut_map("j"), **({"d": cut_map("d")} if vdj else {})},
        delv=delv, delj=delj, deld=deld, deld_pairs=deld_pairs,
        maxpal=model.manifest.palindrome_max, ins=ins, R=R, bias=bias,
        p_nd=p_nd, d2_given_d1=d2_given_d1, deld2=deld2, deld2_pairs=deld2_pairs,
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


def _d_seg(prep: _GenPrep, d: str, rng, deld: dict, deld_pairs: dict, *, min1: bool) -> str:
    """Sample D's 5'/3' deletions and return its CDR3 contribution.

    If ``min1``, resample until the D contributes ≥1 nt (a tandem D must be non-empty to match the
    Pgen partition); gives up after a few tries if the deletion distribution rarely leaves any.
    """
    cutd = prep.cut["d"][d]
    for _ in range(64):
        n5, n3 = deld_pairs[d][_pick(rng, deld[d])]
        start = n5 + prep.maxpal["d_5"]
        end = len(cutd) - (n3 + prep.maxpal["d_3"])
        contrib = cutd[start:end] if start < end else ""
        if contrib or not min1:
            return contrib
    return ""


def _draw(prep: _GenPrep, rng) -> tuple[str, str, str, str, str]:
    """One recombination draw -> (junction_nt, v, d, j, d2). ``d2`` is the second D for a tandem, else ""."""
    vdj = prep.chain_type == "VDJ"
    v = _pick(rng, prep.v)
    d, d2 = "", ""
    if vdj:
        j = _pick(rng, prep.j_marg)
        d = _pick(rng, prep.d_given_j[j])
    else:
        j = _pick(rng, prep.j_given_v[v])

    # V and J each contribute >= 1 nt to the CDR3 (matches the Pgen model — never fully deleted).
    cutv = prep.cut["v"][v]
    len_v = max(1, len(cutv) - (_pick(rng, prep.delv[v]) + prep.maxpal["v_3"]))
    v_contrib = cutv[:len_v]
    cutj = prep.cut["j"][j]
    dj = min(_pick(rng, prep.delj[j]) + prep.maxpal["j_5"], len(cutj) - 1)
    j_contrib = cutj[dj:]

    if not vdj:
        ins_vj = _insert(rng, int(_pick(rng, prep.ins["vj"])), prep.R["vj"], prep.bias["vj"], from_right=False)
        return v_contrib + ins_vj + j_contrib, v, d, j, d2

    n_d = int(_pick(rng, prep.p_nd)) if prep.p_nd is not None else 1
    d2_dist = prep.d2_given_d1.get(d) if (n_d == 2 and prep.d2_given_d1) else None
    if d2_dist is not None and len(d2_dist[0]):  # tandem D, and D1 has >=1 possible D2 partner
        d2 = _pick(rng, d2_dist)  # P(D2 | D1)
        d1c = _d_seg(prep, d, rng, prep.deld, prep.deld_pairs, min1=True)
        d2c = _d_seg(prep, d2, rng, prep.deld2, prep.deld2_pairs, min1=True)
        ins_vd = _insert(rng, int(_pick(rng, prep.ins["vd"])), prep.R["vd"], prep.bias["vd"], from_right=False)
        ins_dd = _insert(rng, int(_pick(rng, prep.ins["dd"])), prep.R["dd"], prep.bias["dd"], from_right=False)
        ins_dj = _insert(rng, int(_pick(rng, prep.ins["dj"])), prep.R["dj"], prep.bias["dj"], from_right=True)
        cdr3 = v_contrib + ins_vd + d1c + ins_dd + d2c + ins_dj + j_contrib
    else:  # single D (n_D=1; n_D=0 folds in as a fully-trimmed D)
        d1c = _d_seg(prep, d, rng, prep.deld, prep.deld_pairs, min1=False)
        ins_vd = _insert(rng, int(_pick(rng, prep.ins["vd"])), prep.R["vd"], prep.bias["vd"], from_right=False)
        ins_dj = _insert(rng, int(_pick(rng, prep.ins["dj"])), prep.R["dj"], prep.bias["dj"], from_right=True)
        cdr3 = v_contrib + ins_vd + d1c + ins_dj + j_contrib
    return cdr3, v, d, j, d2


def generate(model: Model, n: int, *, seed: int | None = None, productive_only: bool = False) -> pl.DataFrame:
    """Sample ``n`` recombined CDR3s from the model.

    Args:
        model: The generative model.
        n: Number of sequences to return.
        seed: RNG seed for reproducibility.
        productive_only: If True, reject out-of-frame / stop-codon draws and keep sampling.

    Returns:
        DataFrame with ``junction_nt, junction_aa, v_call, d_call, d2_call, j_call, productive``. ``d2_call``
        is the second D of a tandem (``n_D=2``) draw, else null; for a single-D model it is all-null.
    """
    prep = prepare_generation(model)
    rng = np.random.default_rng(seed)
    cols = ("junction_nt", "junction_aa", "v_call", "d_call", "d2_call", "j_call", "productive")
    rows = {k: [] for k in cols}
    got = 0
    guard = 0
    max_guard = n * 10000 + 1000
    while got < n:
        guard += 1
        if guard > max_guard:
            raise RuntimeError("generation exceeded attempt budget (productive draws too rare?)")
        cdr3, v, d, j, d2 = _draw(prep, rng)
        aa = translate(cdr3)
        productive = len(cdr3) % 3 == 0 and "*" not in aa and len(cdr3) > 0
        if productive_only and not productive:
            continue
        rows["junction_nt"].append(cdr3)
        rows["junction_aa"].append(aa)
        rows["v_call"].append(v)
        rows["d_call"].append(d if d else None)
        rows["d2_call"].append(d2 if d2 else None)
        rows["j_call"].append(j)
        rows["productive"].append(productive)
        got += 1
    return pl.DataFrame(rows)
