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


def _chain_scores(model, aas: list, vs: list, js: list, condition_vj: bool) -> list:
    """Pgen per row for one chain -- one native batch over the **distinct** clonotype keys.

    Three things this has to get right, and a per-row loop over ``native.pgen_aa`` only got the
    first two:

    * **Deduplicate, do not cache.** Cells of an expanded clone share a clonotype key and the
      native Pgen is deterministic in ``(junction, v, j)``, so scoring the distinct keys within one
      call is exact deduplication -- nothing survives the call.
    * **An unnameable V/J is a null for its own rows only.** The model is keyed by allele and
      raises on anything else, deliberately: mapping an unknown call to "marginalise" returned a
      value 2.38x too high with no error. ``pgen_aa_batch`` resolves names itself and would raise
      for the **whole batch**, so names are resolved here first and the unresolvable keys never
      enter it.
    * **Batch it.** ``pgen_aa_batch`` releases the GIL and threads across the input; the per-row
      loop was serial, and a whole barcoded dataset is tens of thousands of receptors.
    """
    n = len(aas)
    if model is None:
        return [None] * n
    _pm, vi, ji = native.pack(model)
    row_key: list = [None] * n
    nameable: dict = {}
    for i in range(n):
        aa = aas[i]
        if not isinstance(aa, str) or not aa:
            continue
        k = (aa, vs[i], js[i]) if condition_vj else (aa, None, None)
        row_key[i] = k
        if k not in nameable:
            try:
                native._gene_idx(vi, k[1], "V")
                native._gene_idx(ji, k[2], "J")
                nameable[k] = True
            except KeyError:
                nameable[k] = False
    keys = [k for k, good in nameable.items() if good]
    scored = dict(zip(keys, native.pgen_aa_batch(
        model, [k[0] for k in keys], [k[1] for k in keys], [k[2] for k in keys]))) if keys else {}
    return [None if k is None else scored.get(k) for k in row_key]


def _warn_if_all_null(values, model, locus, col) -> None:
    """A whole column of nulls is almost always a naming mismatch -- say so, don't ship it."""
    if model is None or not values or any(v is not None for v in values):
        return
    import warnings

    warnings.warn(
        f"paired_pgen: every {locus} chain scored null ({len(values)} rows). The usual cause "
        f"is a {col} naming the model does not carry; check a value against the model's "
        "alleles, or pass condition_vj=False to marginalise over V/J deliberately.",
        UserWarning, stacklevel=3,
    )


def paired_pgen(
    paired: pl.DataFrame,
    *,
    source: str = "olga",
    condition_vj: bool = True,
    resolve_genes: bool = True,
    alpha_locus: str | None = None,
    beta_locus: str | None = None,
) -> pl.DataFrame:
    """Add ``pgen_alpha``, ``pgen_beta`` and ``pgen_paired`` to a paired single-cell frame.

    Args:
        paired: A paired-chain frame (:func:`vdjtools.sc.pair.pair_chains` layout).
        source: Bundled model set — ``"olga"`` (OLGA-derived) or ``"learned"`` (native EM).
        condition_vj: Condition each chain's Pgen on its V/J call. ``False`` marginalises
            over all V/J unconditionally.
        resolve_genes: Resolve a **gene**-level call (``TRBV10-3``) to a representative
            model allele (``TRBV10-3*01``) before scoring -- see
            :func:`vdjtools.model.native.gene_to_allele`.
            Default ``True``, because CellRanger reports genes and without this every 10x
            row scores ``None``. Set ``False`` to score only exact allele matches.
        alpha_locus: Locus of the α/light chain (e.g. ``"TRA"``, ``"IGK"``); inferred from
            the ``alpha_v_call`` prefix if ``None``.
        beta_locus: Locus of the β/heavy chain (e.g. ``"TRB"``, ``"IGH"``); inferred from
            the ``beta_v_call`` prefix if ``None``.

    Returns:
        ``paired`` with three added Float64 columns. ``pgen_alpha`` / ``pgen_beta`` are null
        for a cell missing that chain's junction, or carrying a V/J call the model does not
        know; ``pgen_paired`` is null unless both are set.

    Warns:
        UserWarning: If every chain of a locus scored null -- the usual cause is a V/J
            naming the model does not recognise, which would otherwise be an entire column
            of silent nulls.
    """
    a_loc = alpha_locus or (_infer_locus(paired[ALPHA_V]) if ALPHA_V in paired.columns else None)
    b_loc = beta_locus or (_infer_locus(paired[BETA_V]) if BETA_V in paired.columns else None)
    ma = load_bundled(a_loc, source) if a_loc else None
    mb = load_bundled(b_loc, source) if b_loc else None

    # gene -> representative allele, per model (the two loci have disjoint gene names).
    aliases: dict[str, str] = {}
    if resolve_genes and condition_vj:
        for m in (ma, mb):
            if m is not None:
                aliases.update(native.gene_to_allele(m))

    def _call(name):
        return aliases.get(name, name) if name else name

    def col(name):
        return (paired[name].to_list() if name in paired.columns
                else [None] * paired.height)

    pa = _chain_scores(ma, col(ALPHA_AA), [_call(x) for x in col(ALPHA_V)],
                       [_call(x) for x in col(ALPHA_J)], condition_vj)
    pb = _chain_scores(mb, col(BETA_AA), [_call(x) for x in col(BETA_V)],
                       [_call(x) for x in col(BETA_J)], condition_vj)
    pp = [a * b if (a is not None and b is not None) else None for a, b in zip(pa, pb)]

    _warn_if_all_null(pa, ma, a_loc, ALPHA_V)
    _warn_if_all_null(pb, mb, b_loc, BETA_V)

    return paired.with_columns(
        pl.Series("pgen_alpha", pa, dtype=pl.Float64),
        pl.Series("pgen_beta", pb, dtype=pl.Float64),
        pl.Series("pgen_paired", pp, dtype=pl.Float64),
    )
