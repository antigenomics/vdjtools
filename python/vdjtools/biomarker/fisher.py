"""Incidence-based biomarker association (Fisher's exact) — the Emerson 2017 method.

For each clonotype **feature**, a 2×2 subject-incidence table

===============  ==========  ==========
                 phenotype+  phenotype−
===============  ==========  ==========
feature present  ``a``       ``b``
feature absent   ``c``       ``d``
===============  ==========  ==========

is tested with Fisher's exact test, where ``a`` = number of phenotype-positive subjects
whose repertoire *contains* the feature, etc. (Emerson et al., *Nat Genet* 2017,
doi:10.1038/ng.3822). A feature is public if it is present in ``>= min_incidence`` subjects.

Two axes are exposed as first-class options:

- **V/J match requirement** — the ``key`` tuple. ``(junction_aa,)`` matches on CDR3 alone;
  ``(junction_aa, v_call)`` also requires the V gene; ``(junction_aa, v_call, j_call)`` (the
  default, Emerson's definition) requires both.
- **exact vs 1-mismatch CDR3 matching** — ``match``. ``"exact"`` groups on the literal
  key (pure polars, no extra deps). ``"1mm"`` first groups near-variant keys into
  metaclonotypes (:func:`vdjtools.biomarker.metaclonotypes`, needs the ``[overlap]`` extra)
  and tests one row per metaclonotype.

Scale: the incidence table is one streamed ``group_by`` over a
:func:`vdjtools.io.scan_cohort` LazyFrame (a cohort far larger than RAM never
materialises), and the Fisher p-values are computed **vectorised** over the whole feature
table via the hypergeometric tail — no per-feature Python loop — so millions of features
score in a single ``scipy`` call.
"""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control, hypergeom

from ..io import schema as S
from ..io.cohort import SAMPLE_ID
from ..overlap.metrics import DEFAULT_KEY
from .metaclonotype import metaclonotypes

_ALTERNATIVES = ("greater", "less", "two-sided")

#: Result columns appended to the feature key (in order).
_STAT_COLS = ["incidence", "n_pos_present", "n_neg_present", "n_pos", "n_neg",
              "odds_ratio", "log2_or", "p_value", "q_value", "direction"]


def fisher_association(
    cohort: pl.LazyFrame | pl.DataFrame,
    phenotype: pl.DataFrame | pl.LazyFrame,
    *,
    pheno_col: str,
    key: tuple[str, ...] = DEFAULT_KEY,
    match: str = "exact",
    min_incidence: int = 2,
    alternative: str = "greater",
    productive_only: bool = True,
    strip_allele: bool = True,
    scope: str = "1,0,0,1",
    threads: int = 0,
) -> pl.DataFrame:
    """Test each clonotype feature's subject incidence against a binary phenotype.

    Args:
        cohort: A clonotype cohort — a :class:`polars.LazyFrame` from
            :func:`vdjtools.io.scan_cohort` (recommended, streamed) or an in-memory
            :class:`polars.DataFrame`. Must carry a ``sample_id`` column plus the ``key``
            columns.
        phenotype: One row per subject with ``sample_id`` and the binary ``pheno_col``
            (bool / 0-1); subjects with a null label are dropped from both classes.
        pheno_col: Name of the binary phenotype column in ``phenotype``.
        key: Feature key — a subset of ``(junction_aa, v_call, j_call)`` (must include
            ``junction_aa``). This is the **V/J match requirement**.
        match: ``"exact"`` (literal key) or ``"1mm"`` (group near-variants into
            metaclonotypes first; needs the ``[overlap]`` extra).
        min_incidence: Minimum number of subjects a feature must appear in to be tested.
        alternative: ``"greater"`` (enrichment in phenotype+, one-tailed — the CMV setting),
            ``"less"`` (depletion), or ``"two-sided"`` (the HLA setting, by doubling).
        productive_only: Drop CDR3s containing a stop (``*``) or frameshift (``_``) marker.
        strip_allele: Strip the IMGT allele suffix (``*01``) from V/J calls in ``key``.
        scope: Edit scope for ``match="1mm"`` (see :func:`metaclonotypes`).
        threads: Worker threads for the native metaclonotype search (``0`` = all cores).

    Returns:
        One row per tested feature: the key columns (or ``meta_id`` + a representative key
        and ``n_members`` when ``match="1mm"``) followed by ``incidence, n_pos_present,
        n_neg_present, n_pos, n_neg, odds_ratio, log2_or, p_value, q_value`` (Benjamini-
        Hochberg) and ``direction`` (``"enriched"``/``"depleted"``), sorted by ``p_value``.

    Raises:
        ValueError: If ``alternative``/``match`` is unrecognised, ``key`` omits ``junction_aa``,
            or the phenotype has only one class after dropping unknowns.
    """
    if alternative not in _ALTERNATIVES:
        raise ValueError(f"alternative must be one of {_ALTERNATIVES}; got {alternative!r}")
    if match not in ("exact", "1mm"):
        raise ValueError(f"match must be 'exact' or '1mm'; got {match!r}")
    key = tuple(key)
    if S.JUNCTION_AA not in key:
        raise ValueError(f"key must include {S.JUNCTION_AA!r}; got {key}")

    lf = cohort.lazy()

    # Phenotype -> boolean; subjects with an unknown label are excluded entirely.
    ph = (phenotype.lazy()
          .select(SAMPLE_ID, pl.col(pheno_col).cast(pl.Boolean).alias("_ph"))
          .drop_nulls("_ph"))
    ph_df = ph.collect()
    n_pos = int(ph_df["_ph"].sum())
    n_tot = ph_df.height
    n_neg = n_tot - n_pos
    if n_pos == 0 or n_neg == 0:
        raise ValueError("phenotype has only one class after dropping unknown labels")

    # Feature engineering: productive filter + allele stripping on the key's V/J columns.
    lf = lf.filter(pl.col(S.JUNCTION_AA).is_not_null())
    if productive_only:
        lf = lf.filter(~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
    if strip_allele:
        for c in (S.V_CALL, S.J_CALL):
            if c in key:
                lf = lf.with_columns(S.strip_allele(pl.col(c)).alias(c))

    rep = None
    if match == "exact":
        idcols = list(key)
        feat_lf = lf.select([*idcols, SAMPLE_ID])
    else:
        uniq_keys = lf.select(list(key)).unique().collect()
        mm = metaclonotypes(uniq_keys, scope=scope,
                            match_v=(S.V_CALL in key), match_j=(S.J_CALL in key),
                            threads=threads)
        # Representative per metaclonotype = the highest-count member; carry the size too.
        kc = lf.group_by(list(key)).agg(pl.col(S.COUNT).sum().alias("_kc")).collect()
        on = [c for c in key if c in mm.columns]
        rep = (mm.join(kc, on=on, how="left")
               .with_columns(pl.col("_kc").fill_null(0))
               .sort("_kc", descending=True)
               .group_by("meta_id", maintain_order=True)
               .agg([pl.col(c).first().alias(c) for c in key] + [pl.len().alias("n_members")]))
        feat_lf = lf.join(mm.lazy(), on=on, how="inner").select(["meta_id", SAMPLE_ID])
        idcols = ["meta_id"]

    inc_lf = feat_lf.join(ph, on=SAMPLE_ID, how="inner")
    if min_incidence > 1:
        # Cheap superset prefilter: a feature's distinct-subject incidence is <= its raw row
        # count, so features with fewer than min_incidence rows cannot be public. Dropping them
        # with a light pl.len() count before the (memory-heavier) n_unique pass keeps peak RAM
        # bounded on huge cohorts — the millions of private clonotypes never build a hash-set.
        # Correctness-preserving: the kept set is a superset of the true public features.
        public = (inc_lf.group_by(idcols).agg(pl.len().alias("_rows"))
                  .filter(pl.col("_rows") >= min_incidence).select(idcols))
        inc_lf = inc_lf.join(public, on=idcols, how="semi")

    # Subject-incidence table in one streamed aggregation: distinct subjects carrying the
    # feature (incidence) and distinct phenotype-positive ones (a). n_unique folds the
    # per-subject dedup into the group_by, so no global .unique() materialises.
    inc = (inc_lf.group_by(idcols)
           .agg(pl.col(SAMPLE_ID).n_unique().alias("incidence"),
                pl.col(SAMPLE_ID).filter(pl.col("_ph")).n_unique().alias("_a"))
           .filter(pl.col("incidence") >= min_incidence)
           .collect(engine="streaming"))

    # A zero-row `inc` (nothing reaches min_incidence) flows through unchanged: the scipy
    # tails and numpy ops are empty-safe, so the result keeps the same column layout as a
    # non-empty one (important for pl.concat across per-partition calls).
    a = inc["_a"].to_numpy().astype(np.int64)
    present = inc["incidence"].to_numpy().astype(np.int64)
    b = present - a                     # phenotype- subjects with the feature
    c = n_pos - a                       # phenotype+ subjects without it
    d = n_neg - b                       # phenotype- subjects without it

    # Fisher one-tailed tails == hypergeometric survival/cdf, vectorised over all features.
    p_greater = hypergeom.sf(a - 1, n_tot, n_pos, present)   # P(X >= a): enrichment in +
    p_less = hypergeom.cdf(a, n_tot, n_pos, present)         # P(X <= a): depletion
    if alternative == "greater":
        p = p_greater
    elif alternative == "less":
        p = p_less
    else:
        p = np.minimum(1.0, 2.0 * np.minimum(p_greater, p_less))

    orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))  # Haldane-Anscombe corrected
    q = false_discovery_control(p, method="bh") if p.size else p
    direction = np.where(a * n_neg > b * n_pos, "enriched", "depleted")

    res = inc.drop("_a").with_columns(
        pl.Series("n_pos_present", a),
        pl.Series("n_neg_present", b),
        pl.lit(n_pos, dtype=pl.Int64).alias("n_pos"),
        pl.lit(n_neg, dtype=pl.Int64).alias("n_neg"),
        pl.Series("odds_ratio", orr),
        pl.Series("log2_or", np.log2(orr)),
        pl.Series("p_value", p),
        pl.Series("q_value", q),
        pl.Series("direction", direction),
    )
    if rep is not None:
        res = res.join(rep, on="meta_id", how="left")
    return res.sort("p_value")
