"""Chain resolution, paired-receptor assembly, and doublet / mispairing QC.

A 10x cell's contigs are noisy: barcodes collide (doublets), ambient mRNA leaks a
spurious light chain, and a genuine α/β cell can legitimately carry two productive α.
These functions clean that up on the flat single-cell frame from :mod:`vdjtools.sc.read`.

Locus roles (heavy / light) follow the standard receptor families::

    heavy   = {TRB, TRD, IGH}        # one per cell
    light   = {TRA, TRG}             # one, sometimes two (dual-α)
    B-light = {IGK, IGL}             # one, sometimes two

The per-cell chain ranking key is ``(-duplicate_count, -umi_count, sequence_id)``
everywhere: most reads first, then most UMIs, then a stable id tie-break. The
thresholds encode the mirpy-derived rule *one heavy but possibly two light*, and are
reimplemented here on polars (no mirpy dependency).
"""
from __future__ import annotations

import polars as pl

from ..io.schema import CDR3_AA, J_CALL, LOCUS, V_CALL
from .read import CELL_ID, COUNT, SEQUENCE_ID, UMI_COUNT

HEAVY_LOCI: tuple[str, ...] = ("TRB", "TRD", "IGH")
LIGHT_LOCI: tuple[str, ...] = ("TRA", "TRG")
B_LIGHT_LOCI: tuple[str, ...] = ("IGK", "IGL")

#: locus-pair family -> (light/chain1 locus, heavy/chain2 locus).
LOCUS_PAIR_TO_LOCI: dict[str, tuple[str, str]] = {
    "TRA_TRB": ("TRA", "TRB"),
    "TRG_TRD": ("TRG", "TRD"),
    "IGH_IGK": ("IGK", "IGH"),
    "IGH_IGL": ("IGL", "IGH"),
}

#: Rank key: most reads, then most UMIs, then stable sequence_id.
_SORT = [COUNT, UMI_COUNT, SEQUENCE_ID]
_DESC = [True, True, False]


def _role(locus_expr: pl.Expr) -> pl.Expr:
    """Map a locus to its role bucket (``heavy`` / ``light`` / ``b_light`` / null)."""
    return (
        pl.when(locus_expr.is_in(list(HEAVY_LOCI))).then(pl.lit("heavy"))
        .when(locus_expr.is_in(list(LIGHT_LOCI))).then(pl.lit("light"))
        .when(locus_expr.is_in(list(B_LIGHT_LOCI))).then(pl.lit("b_light"))
        .otherwise(pl.lit(None, dtype=pl.Utf8))
    )


def _keep_light(
    group: pl.DataFrame,
    *,
    secondary_ratio: float,
    secondary_min_umi: int,
    secondary_min_dup: int,
) -> pl.DataFrame:
    """Keep the top light chain, plus a second only if it clears every threshold."""
    if group.height <= 1:
        return group
    ranked = group.sort(_SORT, descending=_DESC, nulls_last=True)
    first, second = ranked.row(0, named=True), ranked.row(1, named=True)
    first_dup = max(1, int(first[COUNT] or 0))
    first_umi = max(1, int(first[UMI_COUNT] or 0))
    second_dup = int(second[COUNT] or 0)
    second_umi = int(second[UMI_COUNT] or 0)
    keep_two = (
        second_dup / first_dup > secondary_ratio
        and second_umi / first_umi > secondary_ratio
        and second_umi >= secondary_min_umi
        and second_dup >= secondary_min_dup
    )
    return ranked.head(2) if keep_two else ranked.head(1)


def resolve_chains(
    rearr: pl.DataFrame,
    *,
    secondary_ratio: float = 0.1,
    secondary_min_umi: int = 2,
    secondary_min_dup: int = 5,
) -> pl.DataFrame:
    """Reduce over-expanded per-cell chains to one heavy and one (or two) light.

    Per ``cell_id``:

    - keep **exactly the top-1** heavy chain (``TRB`` / ``TRD`` / ``IGH``);
    - keep the **top-1** light chain, and a **second** light chain only when all of
      ``second_dup/first_dup > secondary_ratio``, ``second_umi/first_umi >
      secondary_ratio``, ``second_umi >= secondary_min_umi`` and ``second_dup >=
      secondary_min_dup`` hold (the dual-α allowance);
    - the same secondary rule applies jointly across ``IGK`` + ``IGL`` for B-cells.

    Args:
        rearr: Single-cell long frame (:data:`vdjtools.sc.read.SC_COLUMNS`).
        secondary_ratio: Minimum second/first ratio (on both reads and UMIs) to admit
            a second light chain.
        secondary_min_umi: Minimum absolute UMI count for a second light chain.
        secondary_min_dup: Minimum absolute read count for a second light chain.

    Returns:
        The cleaned per-cell contigs (same columns as the input), ordered by cell then
        rank. Contigs on loci outside the receptor roles are dropped.
    """
    withrole = rearr.with_columns(_role(pl.col(LOCUS)).alias("_role"))
    withrole = withrole.filter(pl.col("_role").is_not_null())
    kept: list[pl.DataFrame] = []
    for _, cell in withrole.group_by(CELL_ID, maintain_order=True):
        heavies = cell.filter(pl.col("_role") == "heavy")
        if heavies.height:
            kept.append(heavies.sort(_SORT, descending=_DESC, nulls_last=True).head(1))
        for role in ("light", "b_light"):
            grp = cell.filter(pl.col("_role") == role)
            if grp.height:
                kept.append(_keep_light(
                    grp, secondary_ratio=secondary_ratio,
                    secondary_min_umi=secondary_min_umi,
                    secondary_min_dup=secondary_min_dup,
                ))
    if not kept:
        return rearr.head(0)
    return pl.concat(kept).drop("_role")


def pair_chains(
    rearr: pl.DataFrame,
    *,
    locus_pair: str = "TRA_TRB",
    resolve: bool = True,
) -> pl.DataFrame:
    """Assemble paired receptors as the Cartesian product of a cell's light × heavy.

    After (optionally) :func:`resolve_chains`, each cell forms one paired receptor per
    (light, heavy) combination of its chains in the requested family — so a cell with
    two α and one β yields **two** pairs (``<cell>_1``, ``<cell>_2``). Cells missing
    either side of the family are **counted but not emitted** (see
    :func:`chain_multiplicity`).

    Args:
        rearr: Single-cell long frame.
        locus_pair: Family to pair — one of ``"TRA_TRB"``, ``"TRG_TRD"``, ``"IGH_IGK"``,
            ``"IGH_IGL"``. The first locus is the α/light side (``alpha_*`` columns),
            the second the β/heavy side (``beta_*`` columns).
        resolve: Run :func:`resolve_chains` first (default ``True``).

    Returns:
        One row per paired receptor with ``cell_id, pair_id, alpha_v_call,
        alpha_j_call, alpha_cdr3_aa, alpha_umi_count, alpha_duplicate_count`` and the
        matching ``beta_*`` columns.

    Raises:
        ValueError: If ``locus_pair`` is not a recognised family.
    """
    if locus_pair not in LOCUS_PAIR_TO_LOCI:
        raise ValueError(
            f"locus_pair must be one of {sorted(LOCUS_PAIR_TO_LOCI)}; got {locus_pair!r}"
        )
    light_locus, heavy_locus = LOCUS_PAIR_TO_LOCI[locus_pair]
    if resolve:
        rearr = resolve_chains(rearr)

    rows: list[dict] = []
    for cell_id, cell in rearr.group_by(CELL_ID, maintain_order=True):
        cid = cell_id[0] if isinstance(cell_id, tuple) else cell_id
        alphas = cell.filter(pl.col(LOCUS) == light_locus).sort(
            _SORT, descending=_DESC, nulls_last=True)
        betas = cell.filter(pl.col(LOCUS) == heavy_locus).sort(
            _SORT, descending=_DESC, nulls_last=True)
        if alphas.height == 0 or betas.height == 0:
            continue  # incomplete cell: counted in chain_multiplicity, not emitted
        pairs = [(a, b) for b in betas.to_dicts() for a in alphas.to_dicts()]
        multi = len(pairs) > 1
        for idx, (a, b) in enumerate(pairs, start=1):
            rows.append({
                "cell_id": str(cid),
                "pair_id": f"{cid}_{idx}" if multi else str(cid),
                "alpha_v_call": a.get(V_CALL), "alpha_j_call": a.get(J_CALL),
                "alpha_cdr3_aa": a.get(CDR3_AA),
                "alpha_umi_count": a.get(UMI_COUNT), "alpha_duplicate_count": a.get(COUNT),
                "beta_v_call": b.get(V_CALL), "beta_j_call": b.get(J_CALL),
                "beta_cdr3_aa": b.get(CDR3_AA),
                "beta_umi_count": b.get(UMI_COUNT), "beta_duplicate_count": b.get(COUNT),
            })

    schema = {
        "cell_id": pl.Utf8, "pair_id": pl.Utf8,
        "alpha_v_call": pl.Utf8, "alpha_j_call": pl.Utf8, "alpha_cdr3_aa": pl.Utf8,
        "alpha_umi_count": pl.Int64, "alpha_duplicate_count": pl.Int64,
        "beta_v_call": pl.Utf8, "beta_j_call": pl.Utf8, "beta_cdr3_aa": pl.Utf8,
        "beta_umi_count": pl.Int64, "beta_duplicate_count": pl.Int64,
    }
    if not rows:
        return pl.DataFrame(schema=schema)
    return pl.from_dicts(rows, schema=schema)


def chain_multiplicity(rearr: pl.DataFrame, *, locus_pair: str = "TRA_TRB") -> pl.DataFrame:
    """Presence-quadrant histogram ``(n_light, n_heavy) -> cell_count`` for a family.

    Counts, over cells, how many carry each ``(n_light, n_heavy)`` combination of chain
    multiplicities in ``locus_pair`` — the α/β quadrant table used to diagnose
    doublets and dropout. Cells with neither chain in the family contribute a
    ``(0, 0)`` row.

    Args:
        rearr: Single-cell long frame.
        locus_pair: Family to tabulate (see :func:`pair_chains`).

    Returns:
        A ``pl.DataFrame`` with columns ``n_light, n_heavy, cell_count``, sorted by
        ``n_light`` then ``n_heavy``.

    Raises:
        ValueError: If ``locus_pair`` is not a recognised family.
    """
    if locus_pair not in LOCUS_PAIR_TO_LOCI:
        raise ValueError(
            f"locus_pair must be one of {sorted(LOCUS_PAIR_TO_LOCI)}; got {locus_pair!r}"
        )
    light_locus, heavy_locus = LOCUS_PAIR_TO_LOCI[locus_pair]
    per_cell = rearr.group_by(CELL_ID).agg(
        (pl.col(LOCUS) == light_locus).sum().cast(pl.Int64).alias("n_light"),
        (pl.col(LOCUS) == heavy_locus).sum().cast(pl.Int64).alias("n_heavy"),
    )
    return (per_cell.group_by("n_light", "n_heavy").len()
            .rename({"len": "cell_count"})
            .with_columns(pl.col("cell_count").cast(pl.Int64))
            .sort("n_light", "n_heavy"))


def flag_mispairing(
    paired: pl.DataFrame,
    *,
    max_slaves_per_master: int | None = None,
    drop: bool = False,
) -> pl.DataFrame:
    """Flag suspected mispaired / ambient α chains against a master(β) → slave(α) graph.

    Builds, across all cells, how often each master (β) heavy chain co-occurs with each
    slave (α) light chain. For every master its **canonical** slave is the one with the
    most co-occurrences (ties broken by summed read+UMI support). Any pairing whose α is
    **not** the master's canonical slave is flagged as suspected mispairing /
    contamination. If a master pairs with more than ``max_slaves_per_master`` distinct α
    across the dataset, the master itself is flagged as **ambient** (a β smeared across
    too many barcodes).

    Chains are keyed on ``(v_call, j_call, cdr3_aa)`` per side, so identical clonotypes
    across cells are recognised as the same master / slave.

    Args:
        paired: Output of :func:`pair_chains` (``alpha_*`` / ``beta_*`` columns).
        max_slaves_per_master: Distinct-α ceiling above which a master is called
            ambient; ``None`` disables the ambient check.
        drop: If ``True``, remove flagged rows instead of annotating them.

    Returns:
        The paired frame plus ``mispairing_flag`` (bool) and ``mispairing_reason``
        (``"ok"`` / ``"noncanonical_alpha"`` / ``"ambient_master"``). When ``drop`` is
        set, flagged rows are removed and the two columns omitted.

    Raises:
        ValueError: If ``max_slaves_per_master`` is not a positive integer.
    """
    if max_slaves_per_master is not None and int(max_slaves_per_master) <= 0:
        raise ValueError("max_slaves_per_master must be a positive integer when provided")
    if paired.height == 0:
        out = paired.with_columns(
            pl.lit(False).alias("mispairing_flag"),
            pl.lit("ok").alias("mispairing_reason"),
        )
        return out.drop("mispairing_flag", "mispairing_reason") if drop else out

    beta_key = pl.concat_str("beta_v_call", "beta_j_call", "beta_cdr3_aa",
                             separator="|", ignore_nulls=False)
    alpha_key = pl.concat_str("alpha_v_call", "alpha_j_call", "alpha_cdr3_aa",
                              separator="|", ignore_nulls=False)
    work = paired.with_columns(
        beta_key.alias("_mkey"),
        alpha_key.alias("_skey"),
        (pl.col("alpha_duplicate_count").fill_null(0)
         + pl.col("alpha_umi_count").fill_null(0)).alias("_support"),
    )

    # master -> slave edge co-occurrence (one contribution per row) + support.
    edges = (work.group_by("_mkey", "_skey")
             .agg(pl.len().alias("_edge_count"), pl.col("_support").sum().alias("_edge_support")))
    # canonical slave per master = argmax (edge_count, support, skey).
    canonical = (edges.sort(["_edge_count", "_edge_support", "_skey"],
                            descending=[True, True, False])
                 .group_by("_mkey", maintain_order=True)
                 .agg(pl.col("_skey").first().alias("_canonical_skey")))
    # distinct slaves per master (for the ambient check).
    degree = edges.group_by("_mkey").agg(pl.col("_skey").n_unique().alias("_master_degree"))

    work = work.join(canonical, on="_mkey", how="left").join(degree, on="_mkey", how="left")

    # Within-cell dual-α is legitimate biology (10-30% of T cells carry two productive
    # TRA), not contamination: if a master's canonical α is ALSO present in this cell,
    # the cell's other α for that master is a real second chain, not a mispairing. Only
    # flag a non-canonical α whose cell LACKS the canonical α (a cross-barcode smear).
    canon_in_cell = (work.group_by("cell_id", "_mkey")
                     .agg((pl.col("_skey") == pl.col("_canonical_skey")).any()
                          .alias("_canon_in_cell")))
    work = work.join(canon_in_cell, on=["cell_id", "_mkey"], how="left")

    ambient = (
        (pl.col("_master_degree") > max_slaves_per_master)
        if max_slaves_per_master is not None
        else pl.lit(False)
    )
    noncanon = (pl.col("_skey") != pl.col("_canonical_skey")) & ~pl.col("_canon_in_cell")
    work = work.with_columns(
        pl.when(ambient).then(pl.lit("ambient_master"))
        .when(noncanon).then(pl.lit("noncanonical_alpha"))
        .otherwise(pl.lit("ok")).alias("mispairing_reason"),
    ).with_columns(
        (pl.col("mispairing_reason") != "ok").alias("mispairing_flag"),
    )

    work = work.drop("_mkey", "_skey", "_support", "_canonical_skey", "_master_degree",
                     "_canon_in_cell")
    if drop:
        return work.filter(~pl.col("mispairing_flag")).drop("mispairing_flag", "mispairing_reason")
    return work
