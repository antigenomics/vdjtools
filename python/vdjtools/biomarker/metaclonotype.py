"""Metaclonotype grouping — collapse edit-scope-neighbour clonotypes into one feature.

Emerson-style incidence association tests each *exact* public TCRβ. A **metaclonotype**
instead groups clonotypes whose CDR3s are within an edit scope (default one substitution)
— optionally requiring the same V and/or J — so a family of near-variants counts as one
biomarker feature. The fuzzy search is **not reimplemented**: it is delegated to
:func:`vdjmatch.cluster.overlap` (native ``seqtree`` engine), exactly as
:mod:`vdjtools.overlap.fuzzy` does. Connected components of the within-scope neighbour
graph (single-linkage, union-find) become the metaclonotypes.

Scale: the unique clonotype keys are clustered **once** (not per subject). Partitioning by
V and/or J means the all-pairs search only ever runs within one gene group, and the native
call releases the GIL — so ~1M unique CDR3s are grouped in a handful of multi-threaded
passes.
"""
from __future__ import annotations

import polars as pl

from ..io.schema import JUNCTION_AA, J_CALL, V_CALL

_VDJMATCH_HINT = (
    "vdjmatch is required for vdjtools.biomarker.metaclonotype; install the extra with "
    "`pip install 'vdjtools[overlap]'` (or `pip install vdjmatch>=0.0.1`)."
)


def _require_vdjmatch():
    """Import and return ``vdjmatch.cluster``; raise a helpful error if missing."""
    try:
        import vdjmatch.cluster as cluster  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without vdjmatch
        raise ImportError(_VDJMATCH_HINT) from exc
    return cluster


def metaclonotypes(clonotypes: pl.DataFrame, *, scope: str = "1,0,0,1",
                   match_v: bool = True, match_j: bool = True,
                   threads: int = 0) -> pl.DataFrame:
    """Group unique clonotype keys into metaclonotypes by CDR3 edit-scope neighbourhood.

    Two keys share a ``meta_id`` iff their ``junction_aa`` are within ``scope`` **and** they
    share a V call (when ``match_v``) and a J call (when ``match_j``). Grouping is
    single-linkage (connected components of the neighbour graph).

    Args:
        clonotypes: A clonotype frame; must carry ``junction_aa`` and, when the corresponding
            ``match_*`` flag is set, ``v_call`` / ``j_call``.
        scope: vdjmatch edit-distance scope ``"subs,ins,dels,total"`` (default one
            substitution, length-preserving). ``"0,0,0,0"`` reduces to exact grouping.
        match_v: Require the same ``v_call`` for two keys to be grouped.
        match_j: Require the same ``j_call`` for two keys to be grouped.
        threads: Worker threads for the native search (``0`` = all cores).

    Returns:
        The distinct grouping keys (``junction_aa`` plus ``v_call``/``j_call`` as applicable)
        with an added ``meta_id`` column (compact 0-based integer, singletons included).

    Raises:
        ImportError: If vdjmatch is not installed (see the ``overlap`` extra).
    """
    cluster = _require_vdjmatch()

    group_cols = []
    if match_v and V_CALL in clonotypes.columns:
        group_cols.append(V_CALL)
    if match_j and J_CALL in clonotypes.columns:
        group_cols.append(J_CALL)

    uniq = (clonotypes.select([JUNCTION_AA, *group_cols])
            .drop_nulls(JUNCTION_AA)
            .unique(maintain_order=True)
            .with_row_index("_gid"))
    n = uniq.height
    if n == 0:
        return uniq.drop("_gid").with_columns(pl.Series("meta_id", [], dtype=pl.Int64))

    parent = list(range(n))

    def find(x: int) -> int:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:  # path compression
            parent[x], x = root, parent[x]
        return root

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    # One native fuzzy-overlap per V/J partition (unique CDR3s within the partition), then
    # union-find over the returned positional pairs. No partition (CDR3-only) => one pass.
    groups = uniq.group_by(group_cols, maintain_order=True) if group_cols else [((), uniq)]
    for _, sub in groups:
        cdr3s = sub[JUNCTION_AA].to_list()
        if len(cdr3s) < 2:
            continue
        gids = sub["_gid"].to_list()
        pairs = cluster.overlap(cdr3s, scope=scope, threads=threads)
        for a, b in zip(pairs["a_idx"].to_list(), pairs["b_idx"].to_list()):
            union(gids[a], gids[b])

    # Compact-relabel component roots to 0..K-1 in first-seen order.
    remap: dict[int, int] = {}
    meta = []
    for i in range(n):
        r = find(i)
        if r not in remap:
            remap[r] = len(remap)
        meta.append(remap[r])

    return uniq.drop("_gid").with_columns(pl.Series("meta_id", meta, dtype=pl.Int64))
