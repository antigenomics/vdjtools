"""Bridge from the polars :class:`Model` to the native ``_core`` hot loops.

``pack`` reconstructs the model's dense arrays (gene choice, deletion, insertion, dinucleotide)
and germline cut segments into a C++ :class:`PackedModel`; ``pgen_nt`` calls the native Pgen.
The result matches the pure-Python reference (and OLGA) exactly — the native path is just faster.
"""
from __future__ import annotations

import numpy as np

from .model import Model
from .pgen import prepare

_NT2NUM = {"A": 0, "C": 1, "G": 2, "T": 3}
_pack_cache: dict[int, tuple] = {}


def _encode(seq: str) -> list[int]:
    return [_NT2NUM[c] for c in seq]


def _del_dense(pdel: dict, idx_of: dict, maxpal: int) -> tuple[list[float], int]:
    """{(allele, ndel): p} -> (flat [n_allele * nbins], nbins), index = ndel + maxpal."""
    nbins = max((n + maxpal for (_a, n) in pdel), default=0) + 1
    arr = [0.0] * (len(idx_of) * nbins)
    for (a, n), p in pdel.items():
        arr[idx_of[a] * nbins + (n + maxpal)] = p
    return arr, nbins


def _del_dense_d(pdel: dict, idx_of: dict, maxdl: int, maxdr: int) -> tuple[list[float], int, int]:
    n5 = max((a + maxdl for (_d, a, _b) in pdel), default=0) + 1
    n3 = max((b + maxdr for (_d, _a, b) in pdel), default=0) + 1
    arr = [0.0] * (len(idx_of) * n5 * n3)
    for (d, a, b), p in pdel.items():
        arr[(idx_of[d] * n5 + (a + maxdl)) * n3 + (b + maxdr)] = p
    return arr, n5, n3


def pack(model: Model):
    """Build (and cache) the native :class:`PackedModel` for ``model``."""
    key = id(model)
    if key in _pack_cache:
        return _pack_cache[key]

    from .._core import PackedModel

    prep = prepare(model)
    vdj = model.chain_type == "VDJ"
    v_alleles = model.genomic["genes_v"]["v_allele"].to_list()
    j_alleles = model.genomic["genes_j"]["j_allele"].to_list()
    d_alleles = model.genomic["genes_d"]["d_allele"].to_list() if vdj else []
    vi = {a: i for i, a in enumerate(v_alleles)}
    ji = {a: i for i, a in enumerate(j_alleles)}
    di = {a: i for i, a in enumerate(d_alleles)}

    pm = PackedModel()
    pm.vdj = vdj
    pm.maxpal_v3 = prep.maxpal["v_3"]
    pm.maxpal_j5 = prep.maxpal["j_5"]
    pm.cut_v = [_encode(prep.cut["v"][a]) for a in v_alleles]
    pm.cut_j = [_encode(prep.cut["j"][a]) for a in j_alleles]
    pm.func_v = [vi[a] for a in prep.functional_v]
    pm.func_j = [ji[a] for a in prep.functional_j]
    pm.pv = [float(prep.p_v.get(a, 0.0)) for a in v_alleles]
    pm.del_v, pm.nbins_v = _del_dense(prep.p_del["v"], vi, prep.maxpal["v_3"])
    pm.del_j, pm.nbins_j = _del_dense(prep.p_del["j"], ji, prep.maxpal["j_5"])

    if vdj:
        pm.maxpal_d5 = prep.maxpal["d_5"]
        pm.maxpal_d3 = prep.maxpal["d_3"]
        pm.cut_d = [_encode(prep.cut["d"][a]) for a in d_alleles]
        pm.func_d = [di[a] for a in prep.functional_d]
        pm.pj = [float(prep.p_j.get(a, 0.0)) for a in j_alleles]
        pm.pd_given_j = [float(prep.p_d_given_j.get((j, d), 0.0)) for j in j_alleles for d in d_alleles]
        pm.del_d, pm.nbins_d5, pm.nbins_d3 = _del_dense_d(prep.p_del["d"], di, prep.maxpal["d_5"], prep.maxpal["d_3"])
        pm.ins_vd = prep.p_ins["vd"].tolist()
        pm.ins_dj = prep.p_ins["dj"].tolist()
        pm.R_vd = prep.R["vd"].reshape(-1).tolist()
        pm.R_dj = prep.R["dj"].reshape(-1).tolist()
        pm.bias_vd = prep.bias["vd"].tolist()
        pm.bias_dj = prep.bias["dj"].tolist()
    else:
        pm.pjv = [float(prep.p_j.get((v, j), 0.0)) for v in v_alleles for j in j_alleles]
        pm.ins_vj = prep.p_ins["vj"].tolist()
        pm.R_vj = prep.R["vj"].reshape(-1).tolist()
        pm.bias_vj = prep.bias["vj"].tolist()

    _pack_cache[key] = (pm, vi, ji)
    return _pack_cache[key]


def pgen_nt(model: Model, cdr3_nt: str, v: str | None = None, j: str | None = None) -> float:
    """Native nucleotide Pgen — same result as :func:`vdjtools.model.pgen.pgen_nt`, faster."""
    from .._core import pgen_nt as _pgen_nt

    pm, vi, ji = pack(model)
    return _pgen_nt(pm, _encode(cdr3_nt.upper()),
                    vi.get(v, -1) if v else -1, ji.get(j, -1) if j else -1)


def pgen_aa(model: Model, cdr3_aa: str, v: str | None = None, j: str | None = None) -> float:
    """Native amino-acid Pgen — same result as :func:`vdjtools.model.pgen.pgen_aa`, much faster."""
    from .._core import pgen_aa as _pgen_aa

    pm, vi, ji = pack(model)
    return _pgen_aa(pm, cdr3_aa.upper(), vi.get(v, -1) if v else -1, ji.get(j, -1) if j else -1)
