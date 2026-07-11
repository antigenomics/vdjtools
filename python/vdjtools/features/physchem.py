"""CDR3 physicochemical-property profiles from the legacy amino-acid property table.

For each clonotype the mean of a property over the residues of a chosen CDR3 region
is computed, then averaged (weighted by reads/frequency, or unweighted) within a
group (e.g. per V-J pairing, or per locus).

Region definitions (over ``cdr3_aa``, length ``L``):
    ``all``     — the entire CDR3.
    ``trimmed`` — ``cdr3_aa[3:-3]`` (conserved-anchor-trimmed core); clonotypes with
        ``L <= 6`` have an empty core and are skipped for this region.
    ``center``  — the middle five residues ``cdr3_aa[L//2-2 : L//2+3]``; clonotypes
        with ``L < 5`` are skipped for this region.
"""
from __future__ import annotations

import functools
import io
from importlib import resources

import polars as pl

from ..io.schema import CDR3_AA, LOCUS, add_locus, weight_expr

#: Default property subset: Kidera-factor-free physicochemistry + the 10 Kidera factors.
DEFAULT_PROPERTIES = (
    "hydropathy", "charge", "polarity", "volume", "strength",
    *(f"kf{i}" for i in range(1, 11)),
)


@functools.lru_cache(maxsize=1)
def load_property_table() -> pl.DataFrame:
    """Load the legacy amino-acid property table (cached).

    The shipped ``resources/aa_property_table.txt`` uses classic-Mac ``\\r`` line
    endings and a leading ``##`` reference comment; both are handled here.

    Returns:
        A ``pl.DataFrame`` with an ``amino_acid`` column and one column per property
        (all property columns cast to ``Float64``).
    """
    raw = resources.files("vdjtools.resources").joinpath("aa_property_table.txt").read_text()
    text = raw.replace("\r\n", "\n").replace("\r", "\n")
    lines = [ln for ln in text.split("\n") if ln.strip() and not ln.startswith("##")]
    tbl = pl.read_csv(io.StringIO("\n".join(lines)), separator="\t")
    prop_cols = [c for c in tbl.columns if c != "amino_acid"]
    return tbl.with_columns(pl.col(c).cast(pl.Float64, strict=False) for c in prop_cols)


def _region_expr(region: str) -> pl.Expr:
    """Return the CDR3-region substring expression (null where the region is empty)."""
    length = pl.col(CDR3_AA).str.len_chars()
    if region == "all":
        return pl.col(CDR3_AA)
    if region == "trimmed":
        return (pl.when(length > 6)
                  .then(pl.col(CDR3_AA).str.slice(3, length - 6))
                  .otherwise(None))
    if region == "center":
        return (pl.when(length >= 5)
                  .then(pl.col(CDR3_AA).str.slice(length // 2 - 2, 5))
                  .otherwise(None))
    raise ValueError(f"region must be 'all', 'trimmed' or 'center'; got {region!r}")


def physchem_profile(df: pl.DataFrame, group_by=("v_call", "j_call"),
                     region: str = "all", weight: str = "reads",
                     properties: "tuple[str, ...] | None" = None) -> pl.DataFrame:
    """Group-wise weighted mean of CDR3 physicochemical properties.

    For each clonotype the region residues are looked up in the property table and
    averaged per property (the per-clonotype property mean). These are then combined
    per group by a weighted mean ``Σ w_c m_c / Σ w_c`` (``w_c`` = clonotype weight).

    Args:
        df: A clonotype frame.
        group_by: Grouping column(s). Either the string ``"locus"`` (derived if
            absent) or an iterable of column names present in ``df`` (default
            ``("v_call", "j_call")``).
        region: ``"all"``, ``"trimmed"``, or ``"center"`` (see module docstring).
        weight: ``"reads"``, ``"unique"``, or ``"freq"``.
        properties: Property names to compute. Defaults to
            :data:`DEFAULT_PROPERTIES`.

    Returns:
        Tidy ``pl.DataFrame`` with the group columns, ``region``, ``property`` and
        ``mean_value``, sorted by group then property.

    Raises:
        ValueError: If ``region`` is unknown or a requested property is missing.
    """
    props = list(properties) if properties is not None else list(DEFAULT_PROPERTIES)
    group_cols = [group_by] if isinstance(group_by, str) else list(group_by)

    df = df if LOCUS in df.columns else add_locus(df)
    tbl = load_property_table()
    missing = [p for p in props if p not in tbl.columns]
    if missing:
        raise ValueError(f"unknown properties {missing}; available: {tbl.columns}")

    work = df.with_row_index("_cid").with_columns(
        _region_expr(region).alias("_region"),
        weight_expr(weight).cast(pl.Float64).alias("_w"),
    )
    work = work.filter(pl.col("_region").is_not_null()
                       & (pl.col("_region").str.len_chars() > 0))
    if work.height == 0:
        return pl.DataFrame(schema={**{c: pl.Utf8 for c in group_cols},
                                    "region": pl.Utf8, "property": pl.Utf8,
                                    "mean_value": pl.Float64})

    # explode the region into one row per residue
    work = work.with_columns(pl.col("_region").str.len_chars().alias("_L"))
    work = work.with_columns(
        pl.int_ranges(0, pl.col("_L")).alias("_pos")).explode("_pos", empty_as_null=True)
    work = work.with_columns(pl.col("_region").str.slice(pl.col("_pos"), 1).alias("amino_acid"))
    work = work.join(tbl.select(["amino_acid", *props]), on="amino_acid", how="inner")

    # per-clonotype mean of each property over its residues
    per_clone = work.group_by(["_cid", *group_cols], maintain_order=True).agg(
        pl.col("_w").first(),
        *[pl.col(p).mean().alias(p) for p in props],
    )
    # weighted mean per group, then melt to tidy long form
    grouped = per_clone.group_by(group_cols, maintain_order=True).agg(
        [((pl.col(p) * pl.col("_w")).sum() / pl.col("_w").sum()).alias(p) for p in props]
    )
    tidy = grouped.unpivot(index=group_cols, on=props,
                           variable_name="property", value_name="mean_value")
    return tidy.with_columns(pl.lit(region).alias("region")).select(
        [*group_cols, "region", "property", "mean_value"]
    ).sort([*group_cols, "property"])
