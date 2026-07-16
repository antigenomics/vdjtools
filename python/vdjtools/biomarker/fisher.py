"""Incidence-based Fisher association — the Emerson-2017 method (back-compat entry point).

This is now a thin wrapper over :func:`vdjtools.biomarker.association` with ``test="fisher"``;
it exists so the original column schema and call signature keep working. New code that wants
other tests (χ², Bayesian, permutation), category / stratified conditions, candidate
selection, or co-occurrence should call :func:`vdjtools.biomarker.association` /
:func:`vdjtools.biomarker.cooccurrence` directly.

For each clonotype **feature**, a 2×2 subject-incidence table (present/absent × phenotype±) is
tested with Fisher's exact test (Emerson et al., *Nat Genet* 2017, doi:10.1038/ng.3822).
"""
from __future__ import annotations

import polars as pl

from ..overlap.metrics import DEFAULT_KEY
from .association import association

#: Legacy result columns, in order, appended after the feature key.
_LEGACY = ["incidence", "n_pos_present", "n_neg_present", "n_pos", "n_neg",
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
    """Test each clonotype feature's subject incidence against a binary phenotype (Fisher).

    See :func:`vdjtools.biomarker.association` for the arguments (this fixes ``test="fisher"``).

    Returns:
        One row per tested feature — the key columns (or ``meta_id`` + a representative key and
        ``n_members`` for ``match="1mm"``) followed by ``incidence, n_pos_present, n_neg_present,
        n_pos, n_neg, odds_ratio, log2_or, p_value, q_value`` (Benjamini-Hochberg) and
        ``direction``, sorted by ``p_value``.
    """
    res = association(cohort, phenotype, pheno_col=pheno_col, test="fisher", key=key,
                      match=match, min_incidence=min_incidence, alternative=alternative,
                      productive_only=productive_only, strip_allele=strip_allele,
                      scope=scope, threads=threads)
    key = tuple(key)
    lead = list(key) if match == "exact" else ["meta_id"]
    tail = [] if match == "exact" else [*key, "n_members"]
    order = [*lead, *_LEGACY, *tail]
    return res.select([c for c in order if c in res.columns])
