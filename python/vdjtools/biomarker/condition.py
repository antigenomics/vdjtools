"""Turn subject metadata into the **design frame** the association engine consumes.

A design frame has one row per subject (or per subject × level) with the reserved columns

- ``_pos`` — boolean, the condition-positive indicator for that (subject, level);
- ``_level`` — optional, the category level (HLA allele, homo/hetero, …); absent ⇒ a single
  binary test;
- ``_stratum`` — optional, the Cochran–Mantel–Haenszel stratum; absent ⇒ unstratified.

:func:`vdjtools.biomarker.association` reads exactly these columns, so any of the builders
below (or a hand-built frame with the same columns) can drive it. All keep the subject key
column ``sample_id``.
"""
from __future__ import annotations

import polars as pl

from ..io.cohort import SAMPLE_ID

_TRUE = {"1", "true", "yes", "y", "+", "pos", "positive"}
_FALSE = {"0", "false", "no", "n", "-", "neg", "negative"}


def _to_bool(col: str) -> pl.Expr:
    """Coerce a free-form metadata column to boolean (unknown → null)."""
    s = pl.col(col).cast(pl.String).str.strip_chars().str.to_lowercase()
    return (pl.when(s.is_in(list(_TRUE))).then(True)
            .when(s.is_in(list(_FALSE))).then(False)
            .otherwise(None))


def binary(meta: pl.DataFrame, col: str, *, sample_col: str = SAMPLE_ID) -> pl.DataFrame:
    """Design frame for a binary condition (e.g. CMV ``+``/``−``); unknown labels dropped."""
    return (meta.lazy()
            .select(pl.col(sample_col).alias(SAMPLE_ID), _to_bool(col).alias("_pos"))
            .drop_nulls("_pos").unique(subset=SAMPLE_ID).collect())


def categorical(meta: pl.DataFrame, col: str, *, min_level_size: int = 1,
                sample_col: str = SAMPLE_ID) -> pl.DataFrame:
    """One-vs-rest design for a single-label category (each subject in exactly one level).

    Emits the full subject × level cross with ``_pos = (col == level)``. Levels carried by
    fewer than ``min_level_size`` subjects are dropped. Null category values are excluded.
    """
    m = (meta.select(pl.col(sample_col).alias(SAMPLE_ID), pl.col(col).cast(pl.String).alias("_level"))
         .drop_nulls("_level").unique(subset=SAMPLE_ID))
    levels = (m.group_by("_level").len()
              .filter(pl.col("len") >= min_level_size)["_level"].to_list())
    subjects = m.select(SAMPLE_ID)
    grid = subjects.join(pl.DataFrame({"_level": levels}), how="cross")
    return (grid.join(m, on=SAMPLE_ID, how="left", suffix="_true")
            .with_columns((pl.col("_level") == pl.col("_level_true")).alias("_pos"))
            .select(SAMPLE_ID, "_level", "_pos"))


def hla_alleles(meta: pl.DataFrame, cols: "list[str]", *, resolution: int | None = None,
                min_level_size: int = 1, sample_col: str = SAMPLE_ID) -> pl.DataFrame:
    """Multi-label design over HLA alleles: per allele, carriers vs non-carriers.

    ``cols`` are the allele columns for a locus (e.g. ``["HLA-A.1", "HLA-A.2"]`` or the
    4-digit ``["sample.HLA-A.1", "sample.HLA-A.2"]``). A subject is ``_pos`` for every allele
    it carries. ``resolution`` trims the allele to that many colon-separated fields
    (``resolution=1`` → ``A*02``); ``None`` keeps the field as written. Alleles carried by
    fewer than ``min_level_size`` subjects are dropped.

    **HLA-untyped subjects are dropped**, not counted as non-carriers — consistent with
    :func:`binary` / :func:`zygosity`, which drop unknown phenotypes. Treating an untyped
    subject as a non-carrier of every allele silently inflates the negative arm and biases
    the odds ratio anticonservatively.
    """
    def norm(c: str) -> pl.Expr:
        e = pl.col(c).cast(pl.String).str.strip_chars()
        e = pl.when(e.is_in(["", "NA", "nan", "NaN"])).then(None).otherwise(e)
        if resolution is not None:
            e = e.str.split(":").list.slice(0, resolution).list.join(":")
        return e

    long = (meta.select(pl.col(sample_col).alias(SAMPLE_ID),
                        *[norm(c).alias(f"_a{i}") for i, c in enumerate(cols)])
            .unpivot(index=SAMPLE_ID, on=[f"_a{i}" for i in range(len(cols))],
                     value_name="_level").drop("variable")
            .drop_nulls("_level").unique())
    levels = (long.group_by("_level").agg(pl.col(SAMPLE_ID).n_unique().alias("n"))
              .filter(pl.col("n") >= min_level_size)["_level"].to_list())
    # Typed subjects only — `long` has already dropped null/NA alleles, so a subject with no
    # typing at all is absent here and is excluded rather than becoming a phantom non-carrier.
    subjects = long.select(SAMPLE_ID).unique()
    carried = long.filter(pl.col("_level").is_in(levels)).with_columns(pl.lit(True).alias("_c"))
    grid = subjects.join(pl.DataFrame({"_level": levels},
                                      schema={"_level": pl.String}), how="cross")
    return (grid.join(carried, on=[SAMPLE_ID, "_level"], how="left")
            .with_columns(pl.col("_c").fill_null(False).alias("_pos"))
            .select(SAMPLE_ID, "_level", "_pos"))


def zygosity(meta: pl.DataFrame, locus_cols: "tuple[str, str]", *,
             sample_col: str = SAMPLE_ID) -> pl.DataFrame:
    """Binary homozygous(``_pos=True``)/heterozygous design for a locus's two allele columns."""
    c1, c2 = locus_cols
    return (meta.lazy()
            .select(pl.col(sample_col).alias(SAMPLE_ID),
                    pl.col(c1).cast(pl.String).str.strip_chars().alias("_x"),
                    pl.col(c2).cast(pl.String).str.strip_chars().alias("_y"))
            .filter(pl.col("_x").is_not_null() & pl.col("_y").is_not_null()
                    & (pl.col("_x") != "") & (pl.col("_y") != ""))
            .select(SAMPLE_ID, (pl.col("_x") == pl.col("_y")).alias("_pos"))
            .unique(subset=SAMPLE_ID).collect())


def stratified(meta: pl.DataFrame, pheno_col: str, stratum_col: str, *,
               sample_col: str = SAMPLE_ID) -> pl.DataFrame:
    """Design for a paired condition: binary ``pheno_col`` stratified by ``stratum_col`` (CMH).

    E.g. CMV association conditioned on an HLA group. Subjects with an unknown phenotype or
    stratum are dropped.
    """
    return (meta.lazy()
            .select(pl.col(sample_col).alias(SAMPLE_ID), _to_bool(pheno_col).alias("_pos"),
                    pl.col(stratum_col).cast(pl.String).alias("_stratum"))
            .drop_nulls(["_pos", "_stratum"]).unique(subset=SAMPLE_ID).collect())
