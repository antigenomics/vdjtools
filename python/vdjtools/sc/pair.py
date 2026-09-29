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

from ..io.schema import JUNCTION_AA, J_CALL, LOCUS, V_CALL
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
    # One pass with window functions, not a loop over cells. The per-cell loop concatenated a
    # frame per (cell, role) -- on 10^4-10^6 cells that is that many slices and that many
    # allocations, for arithmetic that is two `over(cell, role)` expressions.
    withrole = (rearr.with_row_index("_i")
                .with_columns(_role(pl.col(LOCUS)).alias("_role"))
                .filter(pl.col("_role").is_not_null()))
    if withrole.height == 0:
        return rearr.head(0)

    grp = [CELL_ID, "_role"]
    ranked = (withrole.sort(_SORT, descending=_DESC, nulls_last=True)
              .with_columns(pl.int_range(pl.len()).over(grp).alias("_rank")))
    # A light chain's runner-up is kept only if it clears EVERY threshold against its own group's
    # best -- which `first().over(...)` is, after that sort. Denominators are floored at 1, exactly
    # as the per-group scalar form's `max(1, ...)` was, so a zero-count best cannot divide by zero.
    # A heavy chain never keeps a second: `keep` below gates the rank-1 branch on the role.
    first_dup = pl.max_horizontal(pl.col(COUNT).fill_null(0).first().over(grp), pl.lit(1))
    first_umi = pl.max_horizontal(pl.col(UMI_COUNT).fill_null(0).first().over(grp), pl.lit(1))
    keep_two = ((pl.col(COUNT).fill_null(0) / first_dup > secondary_ratio)
                & (pl.col(UMI_COUNT).fill_null(0) / first_umi > secondary_ratio)
                & (pl.col(UMI_COUNT).fill_null(0) >= secondary_min_umi)
                & (pl.col(COUNT).fill_null(0) >= secondary_min_dup))
    keep = ((pl.col("_rank") == 0)
            | ((pl.col("_role") != "heavy") & (pl.col("_rank") == 1) & keep_two))

    # Restore the loop's output order: cells in first-appearance order, and within a cell
    # heavy, then light, then b_light, each by rank. Sorting by cell id instead would silently
    # reorder every caller's rows.
    return (ranked.filter(keep)
            .with_columns(pl.col("_i").min().over(CELL_ID).alias("_cell"),
                          pl.col("_role").replace_strict({"heavy": 0, "light": 1, "b_light": 2},
                                                         return_dtype=pl.Int8).alias("_ro"))
            .sort(["_cell", "_ro", "_rank"])
            .drop("_i", "_role", "_rank", "_cell", "_ro"))


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
        alpha_j_call, alpha_junction_aa, alpha_umi_count, alpha_duplicate_count`` and the
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

    # The Cartesian product is a join on cell_id, not a loop that calls `.to_dicts()` per cell --
    # which built a Python dict per contig per cell and then one per emitted pair.
    schema = {
        "cell_id": pl.Utf8, "pair_id": pl.Utf8,
        "alpha_v_call": pl.Utf8, "alpha_j_call": pl.Utf8, "alpha_junction_aa": pl.Utf8,
        "alpha_umi_count": pl.Int64, "alpha_duplicate_count": pl.Int64,
        "beta_v_call": pl.Utf8, "beta_j_call": pl.Utf8, "beta_junction_aa": pl.Utf8,
        "beta_umi_count": pl.Int64, "beta_duplicate_count": pl.Int64,
    }
    keep = [CELL_ID, V_CALL, J_CALL, JUNCTION_AA, UMI_COUNT, COUNT]

    def side(locus: str, tag: str) -> pl.DataFrame:
        """One side's contigs, ranked within the cell exactly as the loop's sort ranked them."""
        return (rearr.with_row_index("_i")
                .filter(pl.col(LOCUS) == locus)
                .sort(_SORT, descending=_DESC, nulls_last=True)
                .select(*keep, "_i")
                .with_columns(pl.int_range(pl.len()).over(CELL_ID).alias(f"_{tag}r"),
                              pl.len().over(CELL_ID).alias(f"_{tag}n"))
                .rename({c: f"{tag}_{c}" for c in keep if c != CELL_ID} | {"_i": f"_{tag}i"}))

    a, b = side(light_locus, "alpha"), side(heavy_locus, "beta")
    if a.height == 0 or b.height == 0:      # incomplete cells are counted, never emitted
        return pl.DataFrame(schema=schema)

    # `pair_id` numbering must match the loop's `for b in betas for a in alphas`: beta-major, from
    # 1, and bare when the cell yields exactly one pair.
    idx = pl.col("_betar") * pl.col("_alphan") + pl.col("_alphar") + 1
    cid = pl.col(CELL_ID).cast(pl.Utf8)
    pairs = (a.join(b, on=CELL_ID, how="inner")
             .with_columns(pl.when(pl.col("_alphan") * pl.col("_betan") > 1)
                           .then(cid + pl.lit("_") + idx.cast(pl.Utf8))
                           .otherwise(cid).alias("pair_id"))
             # cells in first-appearance order, then the loop's (beta, alpha) pair order
             .sort([pl.col("_alphai").min().over(CELL_ID), "_betar", "_alphar"]))
    return pairs.select(
        cid.alias("cell_id"), "pair_id",
        *[pl.col(f"{tag}_{c}").cast(dt).alias(f"{tag}_{name}")
          for tag in ("alpha", "beta")
          for c, name, dt in ((V_CALL, "v_call", pl.Utf8), (J_CALL, "j_call", pl.Utf8),
                              (JUNCTION_AA, "junction_aa", pl.Utf8),
                              (UMI_COUNT, "umi_count", pl.Int64), (COUNT, "duplicate_count", pl.Int64))],
    )


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

    Chains are keyed on ``(v_call, j_call, junction_aa)`` per side, so identical clonotypes
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

    beta_key = pl.concat_str("beta_v_call", "beta_j_call", "beta_junction_aa",
                             separator="|", ignore_nulls=False)
    alpha_key = pl.concat_str("alpha_v_call", "alpha_j_call", "alpha_junction_aa",
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
