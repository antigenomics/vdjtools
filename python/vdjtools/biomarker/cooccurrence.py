"""Feature-vs-feature co-occurrence across a cohort — in-silico α-β pairing & co-specificity.

Two clonotype features that co-occur across subjects more than expected under independence are
candidate **pairs**: a TRA and a TRB from the same clone (Howie 2015 pairSEQ, with subjects
playing the role of wells; Vlasova 2026 in-silico α-β pairing), or two same-chain clonotypes
that recognise the same antigen (De Witt 2018 co-occurrence patterns). This is the same
subject-incidence machinery as :func:`vdjtools.biomarker.association`, applied to a pair of
features: per pair a 2×2 incidence table over the ``n`` subjects profiled for both chains,

    ===============  =========  =========
                     has B      no B
    ===============  =========  =========
    has A            ``n_AB``   …
    no A             …          …
    ===============  =========  =========

is scored by the lift ``θ = n·n_AB/(n_A·n_B)`` (Vlasova; observed / expected co-occurrence),
Fisher's exact / χ², and Benjamini-Hochberg FDR. ``evalue=True`` adds the expected count and a
Poisson upper-tail E-value (the classic control-calibrated co-occurrence significance). The
candidate features per chain are bounded by an incidence threshold (and ``max_features``); the
returned pairs are those with at least ``min_cooccurrence`` co-occurrences.
"""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy.stats import poisson

from ..io import schema as S
from ..io.cohort import SAMPLE_ID
from ..overlap.metrics import DEFAULT_KEY
from . import stats
from .association import _feature_frame, select_candidates


def _incidence_matrix(cohort, chain, key, match, cand, min_incidence, min_incidence_frac,
                      productive_only, strip_allele, scope, threads, max_features):
    """(sample→row, feature-frame, boolean [n_samples_seen × n_features]) for one chain."""
    sub = cohort.lazy().filter(pl.col(S.V_CALL).str.slice(0, 3) == chain)
    feats = select_candidates(sub, key=key, match=match, min_incidence=min_incidence,
                              min_incidence_frac=min_incidence_frac, productive_only=productive_only,
                              strip_allele=strip_allele, scope=scope, threads=threads)
    if cand is not None:
        feats = feats.join(cand.select(list(key)).unique(), on=list(key), how="semi")
    truncated = feats.height > max_features
    feats = feats.head(max_features)                       # top by incidence (select_candidates sorts)
    feat_lf, idcols, _ = _feature_frame(sub, key, match, productive_only=productive_only,
                                        strip_allele=strip_allele, scope=scope, threads=threads,
                                        candidates=(feats if match == "exact" else None))
    fi = feats.select(idcols).with_row_index("_fi")
    pairs = (feat_lf.select([*idcols, SAMPLE_ID]).unique().collect()
             .join(fi, on=idcols, how="inner"))
    return feats, idcols, fi, pairs, truncated


def cooccurrence(
    cohort: pl.LazyFrame | pl.DataFrame,
    *,
    chain_a: str = "TRA",
    chain_b: str | None = "TRB",
    key: tuple[str, ...] = DEFAULT_KEY,
    match: str = "exact",
    test: str = "fisher",
    min_incidence: int = 2,
    min_incidence_frac: float | None = None,
    min_cooccurrence: int = 2,
    candidates_a: pl.DataFrame | None = None,
    candidates_b: pl.DataFrame | None = None,
    evalue: bool = False,
    alternative: str = "greater",
    max_features: int = 2000,
    productive_only: bool = True,
    strip_allele: bool = True,
    scope: str = "1,0,0,1",
    threads: int = 0,
) -> pl.DataFrame:
    """Score co-occurrence of feature pairs across the subjects profiled for both chains.

    Args:
        cohort: Clonotype cohort with a ``sample_id`` column (both chains present per subject
            for α-β pairing).
        chain_a, chain_b: Loci to pair (e.g. ``"TRA"``/``"TRB"``). ``chain_b=None`` (or equal to
            ``chain_a``) does same-chain pairs (upper triangle, self-pairs excluded).
        key, match: Feature definition and exact/1mm scope (as in :func:`association`).
        test: ``"fisher"`` or ``"chi2"`` for the per-pair 2×2 p-value.
        min_incidence, min_incidence_frac: Candidate-feature incidence threshold per chain.
        min_cooccurrence: Keep only pairs co-occurring in ≥ this many subjects.
        candidates_a, candidates_b: Restrict each chain's features to these keys.
        evalue: Also report ``expected`` (= ``n_A·n_B/n``) and a Poisson upper-tail ``e_value``.
        max_features: Cap candidate features per chain (top by incidence); logs if it truncates.

    Returns:
        One row per surviving pair: ``a_<key>`` / ``b_<key>`` columns, ``n`` (subjects with both
        chains), ``n_a``, ``n_b``, ``n_ab``, ``theta`` (lift), ``odds_ratio``, ``log2_or``,
        ``p_value``, ``q_value`` (BH), and — with ``evalue`` — ``expected``, ``e_value``. Sorted
        by ``p_value``.
    """
    if test not in ("fisher", "chi2"):
        raise ValueError(f"test must be 'fisher' or 'chi2'; got {test!r}")
    chain_b = chain_b or chain_a
    same = chain_a == chain_b

    fa, ida, fia, pa, ta = _incidence_matrix(cohort, chain_a, key, match, candidates_a,
                                             min_incidence, min_incidence_frac, productive_only,
                                             strip_allele, scope, threads, max_features)
    if same:
        fb, idb, fib, pb, tb = fa, ida, fia, pa, ta
    else:
        fb, idb, fib, pb, tb = _incidence_matrix(cohort, chain_b, key, match, candidates_b,
                                                 min_incidence, min_incidence_frac, productive_only,
                                                 strip_allele, scope, threads, max_features)
    if ta or tb:
        import warnings
        warnings.warn(f"candidate features capped at max_features={max_features}; some low-"
                      "incidence features were dropped", stacklevel=2)
    if not fa.height or not fb.height:
        return _empty(key)

    # Universe = subjects profiled for BOTH chains (else missing-chain confounds co-occurrence).
    sa = set(pa[SAMPLE_ID].to_list())
    sb = set(pb[SAMPLE_ID].to_list())
    universe = sorted(sa & sb) if not same else sorted(sa)
    n = len(universe)
    if n == 0:
        return _empty(key)
    si = {s: i for i, s in enumerate(universe)}

    M_a = _dense(pa, si, fa.height)
    M_b = M_a if same else _dense(pb, si, fb.height)
    cooc = (M_a.T.astype(np.int64) @ M_b.astype(np.int64))     # [n_a_feat, n_b_feat]
    na = M_a.sum(0).astype(np.int64)
    nb = M_b.sum(0).astype(np.int64)

    ia, ib = np.nonzero(cooc >= min_cooccurrence)
    if same:
        keep = ia < ib                                          # upper triangle, drop self-pairs
        ia, ib = ia[keep], ib[keep]
    if ia.size == 0:
        return _empty(key)
    n_ab = cooc[ia, ib]
    n_a, n_b = na[ia], nb[ib]
    a = n_ab
    b = n_a - n_ab
    c = n_b - n_ab
    d = n - n_a - n_b + n_ab
    theta = n * n_ab / (n_a * n_b)
    p = (stats.fisher_p(a, b, c, d, alternative=alternative) if test == "fisher"
         else stats.chi2_p(a, b, c, d))

    left = fa.select([pl.col(col).alias(f"a_{col}") for col in ida])[ia.tolist()]
    right = fb.select([pl.col(col).alias(f"b_{col}") for col in idb])[ib.tolist()]
    out = left.hstack(right).with_columns(
        pl.lit(n, dtype=pl.Int64).alias("n"), pl.Series("n_a", n_a), pl.Series("n_b", n_b),
        pl.Series("n_ab", n_ab), pl.Series("theta", theta),
        pl.Series("odds_ratio", stats.odds_ratio(a, b, c, d)),
        pl.Series("log2_or", np.log2(stats.odds_ratio(a, b, c, d))),
        pl.Series("p_value", p), pl.Series("q_value", stats.fdr_bh(p)))
    if evalue:
        expected = n_a * n_b / n
        out = out.with_columns(pl.Series("expected", expected),
                               pl.Series("e_value", poisson.sf(n_ab - 1, expected)))
    return out.sort("p_value")


def _dense(pairs: pl.DataFrame, si: dict, n_feat: int) -> np.ndarray:
    """Boolean [n_universe × n_feat] matrix from a (feature-idx, sample) long frame."""
    rows = np.array([si.get(s, -1) for s in pairs[SAMPLE_ID].to_list()], dtype=np.int64)
    keep = rows >= 0                                   # drop samples outside the shared universe
    m = np.zeros((len(si), n_feat), dtype=bool)
    m[rows[keep], pairs["_fi"].to_numpy()[keep]] = True
    return m


def _empty(key: tuple[str, ...]) -> pl.DataFrame:
    cols = {f"{s}_{c}": pl.Series([], dtype=pl.String) for s in ("a", "b") for c in key}
    cols.update(n=pl.Series([], dtype=pl.Int64), n_a=pl.Series([], dtype=pl.Int64),
                n_b=pl.Series([], dtype=pl.Int64), n_ab=pl.Series([], dtype=pl.Int64),
                theta=pl.Series([], dtype=pl.Float64), odds_ratio=pl.Series([], dtype=pl.Float64),
                log2_or=pl.Series([], dtype=pl.Float64), p_value=pl.Series([], dtype=pl.Float64),
                q_value=pl.Series([], dtype=pl.Float64))
    return pl.DataFrame(cols)
