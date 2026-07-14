"""Paired-chain generation probability for single-cell repertoires.

Under chain independence the paired generation probability of a cell is
``Pgen(α) · Pgen(β)`` — the product of each chain's junction Pgen under the native
:mod:`vdjtools.model` engine (bundled per-locus models). This is the single-cell
paired-Pgen residual from Phase 7; it is computed entirely from the native model
(no ``vdjmatch`` dependency).

The paired frame is the :func:`vdjtools.sc.pair.pair_chains` layout — ``alpha_v_call``,
``alpha_j_call``, ``alpha_junction_aa`` and the ``beta_*`` counterparts (α/light and
β/heavy). Each chain's locus is inferred from its V-call prefix (``TRA``/``TRB``, or
``IGK``/``IGL`` + ``IGH`` for BCR) unless given explicitly.

Conditioning on V/J requires the call to match a model **allele** (e.g. ``TRBV20-1*01``);
a gene-level or unmatched call marginalises over all V/J for that chain (still a valid,
if less specific, Pgen). Pass ``condition_vj=False`` to marginalise unconditionally.
"""
from __future__ import annotations

import polars as pl

from ..model import load_bundled, native

ALPHA_V, ALPHA_J, ALPHA_AA = "alpha_v_call", "alpha_j_call", "alpha_junction_aa"
BETA_V, BETA_J, BETA_AA = "beta_v_call", "beta_j_call", "beta_junction_aa"


def _infer_locus(vcalls: pl.Series) -> str | None:
    """Most common three-letter locus prefix among the non-null V calls."""
    pref = vcalls.drop_nulls().str.slice(0, 3)
    if pref.len() == 0:
        return None
    m = pref.mode()
    return m[0] if m.len() else None


def _chain_pgen(model, aa, v, j, condition_vj: bool) -> float | None:
    if not aa:
        return None
    try:
        return native.pgen_aa(model, aa, v if condition_vj else None,
                              j if condition_vj else None)
    except Exception:
        return None


def paired_pgen(
    paired: pl.DataFrame,
    *,
    source: str = "olga",
    condition_vj: bool = True,
    alpha_locus: str | None = None,
    beta_locus: str | None = None,
) -> pl.DataFrame:
    """Add ``pgen_alpha``, ``pgen_beta`` and ``pgen_paired`` to a paired single-cell frame.

    Args:
        paired: A paired-chain frame (:func:`vdjtools.sc.pair.pair_chains` layout).
        source: Bundled model set — ``"olga"`` (OLGA-derived) or ``"learned"`` (native EM).
        condition_vj: Condition each chain's Pgen on its V/J call (must match a model
            allele; otherwise it marginalises). ``False`` marginalises unconditionally.
        alpha_locus: Locus of the α/light chain (e.g. ``"TRA"``, ``"IGK"``); inferred from
            the ``alpha_v_call`` prefix if ``None``.
        beta_locus: Locus of the β/heavy chain (e.g. ``"TRB"``, ``"IGH"``); inferred from
            the ``beta_v_call`` prefix if ``None``.

    Returns:
        ``paired`` with three added Float64 columns. ``pgen_alpha`` / ``pgen_beta`` are null
        for a cell missing that chain's junction; ``pgen_paired`` is null unless both are set.
    """
    a_loc = alpha_locus or (_infer_locus(paired[ALPHA_V]) if ALPHA_V in paired.columns else None)
    b_loc = beta_locus or (_infer_locus(paired[BETA_V]) if BETA_V in paired.columns else None)
    ma = load_bundled(a_loc, source) if a_loc else None
    mb = load_bundled(b_loc, source) if b_loc else None

    pa: list[float | None] = []
    pb: list[float | None] = []
    pp: list[float | None] = []
    # Memoize each chain's Pgen over its distinct clonotype key — cells sharing a clonotype
    # (expanded clones) otherwise recompute the identical native Pgen. Exact: native Pgen is
    # deterministic in (junction, v, j), so a cached value equals the per-row call.
    ca: dict = {}
    cb: dict = {}

    def _memo(cache, model, aa, v, j):
        if model is None or not aa:
            return None
        k = (aa, v, j) if condition_vj else (aa,)
        if k not in cache:
            cache[k] = _chain_pgen(model, aa, v, j, condition_vj)
        return cache[k]

    for r in paired.iter_rows(named=True):
        a = _memo(ca, ma, r.get(ALPHA_AA), r.get(ALPHA_V), r.get(ALPHA_J))
        b = _memo(cb, mb, r.get(BETA_AA), r.get(BETA_V), r.get(BETA_J))
        pa.append(a)
        pb.append(b)
        pp.append(a * b if (a is not None and b is not None) else None)

    return paired.with_columns(
        pl.Series("pgen_alpha", pa, dtype=pl.Float64),
        pl.Series("pgen_beta", pb, dtype=pl.Float64),
        pl.Series("pgen_paired", pp, dtype=pl.Float64),
    )
