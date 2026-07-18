"""General incidence-based biomarker association across a cohort of repertoires.

Generalises the Emerson-2017 Fisher test (:func:`vdjtools.biomarker.fisher_association`) along
four axes, all sharing the same streamed subject-incidence table:

- **statistical test** — Fisher, χ², a Bayesian log-odds posterior, a Beta-Binomial Bayes
  factor, or a label permutation (:mod:`vdjtools.biomarker.stats`);
- **condition type** — binary, a category expanded one-vs-rest per level (HLA allele, zygosity),
  or a paired condition combined by Cochran–Mantel–Haenszel (built with
  :mod:`vdjtools.biomarker.condition`);
- **match scope** — exact ``cdr3aa`` / ``+v`` / ``+v+j`` (the ``key``), a ``fuzzy`` 1-mismatch
  *search* (each candidate keeps its identity and gains incidence), or ``1mm`` metaclonotypes
  (candidates are *merged* into groups);
- **candidate set** — all public features (``min_incidence`` count or ``min_incidence_frac``
  fraction of subjects), or an explicit ``candidates`` list.

The heavy step (the incidence table) is one streamed ``group_by`` over a
:func:`vdjtools.io.scan_cohort` LazyFrame; every test is then vectorised numpy over the whole
feature table.
"""
from __future__ import annotations

import math

import numpy as np
import polars as pl

from ..io import schema as S
from ..io.cohort import SAMPLE_ID
from ..overlap.metrics import DEFAULT_KEY
from . import stats
from .metaclonotype import metaclonotypes

_TESTS = ("fisher", "chi2", "bayes_logodds", "bayes_bf", "permutation")


def _fuzzy_feature_frame(lf: pl.LazyFrame, key: tuple[str, ...], *, scope: str, threads: int,
                         candidates: pl.DataFrame | None):
    """Fuzzy-**search** incidence: each candidate keeps its identity and gains incidence.

    ``incidence(c) = # subjects carrying ANY feature within `scope` of c`` (itself included).
    The 1-mismatch is an incidence-*estimation* trick, not a grouping: the output is a list of
    individual biomarker clonotypes. Contrast :func:`metaclonotypes` (``match="1mm"``), which
    MERGES candidates into groups and tests the group — a legitimate but different operation,
    and the one to use *downstream* of a biomarker list, not to build it.

    Non-``junction_aa`` key columns (V/J) must match EXACTLY between the candidate and the
    neighbour, so ``key=(junction_aa, v_call)`` pins the germline half of the contact surface
    while the CDR3 is allowed to vary by one residue.

    The search delegates to :func:`vdjmatch.cluster.overlap` (seqtree-backed) — never a
    hand-rolled Hamming scan.
    """
    aa = S.JUNCTION_AA
    other = [c for c in key if c != aa]
    rows = lf.select([*key, SAMPLE_ID]).unique().collect()
    universe = rows[aa].unique().sort().to_list()
    cand = (candidates.select(list(key)).unique() if candidates is not None
            else rows.select(list(key)).unique()).sort(list(key))
    cq = cand[aa].unique().sort().to_list()
    if not cq or not universe:
        return rows.lazy().select([*key, SAMPLE_ID]), list(key), None

    import vdjmatch.cluster as vc
    pairs = vc.overlap(cq, universe, scope=scope, threads=threads)
    qmap = pl.DataFrame({"a_idx": np.arange(len(cq), dtype=np.int64), "_cand": cq})
    umap = pl.DataFrame({"b_idx": np.arange(len(universe), dtype=np.int64), aa: universe})
    nb = (pairs.join(qmap, on="a_idx", how="inner").join(umap, on="b_idx", how="inner")
          .select("_cand", aa))
    feat = (nb.join(cand.rename({aa: "_cand"}), on="_cand", how="inner")
            .join(rows, on=[aa, *other], how="inner")
            .select([pl.col("_cand").alias(aa), *other, SAMPLE_ID]).unique())
    return feat.lazy(), list(key), None


def _feature_frame(cohort: pl.LazyFrame | pl.DataFrame, key: tuple[str, ...], match: str, *,
                   productive_only: bool, strip_allele: bool, scope: str, threads: int,
                   candidates: pl.DataFrame | None = None):
    """Feature-engineer the cohort and resolve the feature id columns.

    Returns ``(feat_lf, idcols, rep)`` where ``feat_lf`` has ``[*idcols, sample_id]``, ``idcols``
    is ``list(key)`` (exact/fuzzy) or ``["meta_id"]`` (1mm), and ``rep`` is the per-metaclonotype
    representative frame (or ``None``).
    """
    key = tuple(key)
    if match not in ("exact", "fuzzy", "1mm"):
        raise ValueError(f"match must be 'exact', 'fuzzy' or '1mm'; got {match!r}")
    # `fuzzy`/`1mm` SEARCH on the CDR3, so they need it in the key. `exact` only groups by the
    # key, so any column works there — which is what lets a derived feature (e.g. the `kmer`
    # column from `features.kmer.kmer_cohort`) be tested by the same machinery.
    if match in ("fuzzy", "1mm") and S.JUNCTION_AA not in key:
        raise ValueError(f"match={match!r} searches on {S.JUNCTION_AA!r}, so it must be in key; "
                         f"got {key}. Use match='exact' to test a non-CDR3 feature key.")
    # ...but a key of germline calls ALONE is segment usage wearing a biomarker costume: it asks
    # "is TRBV9 enriched in cases", which is `stats.segment_usage`, not an incidence test over
    # clonotypes. Keep that door shut; the derived-feature door (junction_aa, or `kmer`, or any
    # non-germline column) stays open.
    if S.JUNCTION_AA not in key and set(key) <= {S.V_CALL, S.D_CALL, S.J_CALL, S.C_CALL}:
        raise ValueError(
            f"key {key} is germline calls only — that is segment usage, not a clonotype "
            f"biomarker; use vdjtools.stats.segment_usage. Add {S.JUNCTION_AA!r} or a derived "
            f"feature column (e.g. 'kmer' from features.kmer.kmer_cohort).")
    have = set(cohort.lazy().collect_schema().names())
    missing = [c for c in key if c not in have]
    if missing:
        raise ValueError(f"key columns absent from the cohort: {missing}")

    lf = cohort.lazy()
    if S.JUNCTION_AA in have:
        lf = lf.filter(pl.col(S.JUNCTION_AA).is_not_null())
        if productive_only:
            lf = lf.filter(~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
    if strip_allele:
        for c in (S.V_CALL, S.J_CALL):
            if c in key:
                lf = lf.with_columns(S.strip_allele(pl.col(c)).alias(c))
    if candidates is not None and match != "fuzzy":
        # NOT for fuzzy: `candidates` is the QUERY set there, while the search universe must
        # stay the whole cohort — a candidate's 1mm neighbours are usually not candidates
        # themselves, and dropping them would silently under-count fuzzy incidence.
        cand = candidates.lazy().select(list(key)).unique()
        lf = lf.join(cand, on=list(key), how="semi")

    if match == "exact":
        return lf.select([*key, SAMPLE_ID]), list(key), None
    if match == "fuzzy":
        return _fuzzy_feature_frame(lf, key, scope=scope, threads=threads,
                                    candidates=candidates)

    uniq_keys = lf.select(list(key)).unique().collect()
    mm = metaclonotypes(uniq_keys, scope=scope, match_v=(S.V_CALL in key),
                        match_j=(S.J_CALL in key), threads=threads)
    kc = lf.group_by(list(key)).agg(pl.col(S.COUNT).sum().alias("_kc")).collect()
    on = [c for c in key if c in mm.columns]
    rep = (mm.join(kc, on=on, how="left").with_columns(pl.col("_kc").fill_null(0))
           .sort("_kc", descending=True).group_by("meta_id", maintain_order=True)
           .agg([pl.col(c).first().alias(c) for c in key] + [pl.len().alias("n_members")]))
    return lf.join(mm.lazy(), on=on, how="inner").select(["meta_id", SAMPLE_ID]), ["meta_id"], rep


def select_candidates(cohort: pl.LazyFrame | pl.DataFrame, *, key: tuple[str, ...] = DEFAULT_KEY,
                      match: str = "exact", min_incidence: int = 2,
                      min_incidence_frac: float | None = None, productive_only: bool = True,
                      strip_allele: bool = True, scope: str = "1,0,0,1",
                      threads: int = 0) -> pl.DataFrame:
    """Public features whose subject incidence clears a count and/or fraction threshold.

    ``min_incidence_frac`` (e.g. ``0.05`` = 5% of subjects) is resolved against the cohort's
    distinct ``sample_id`` count and combined with ``min_incidence`` by ``max``. Returns the
    feature key columns (or ``meta_id`` + representative key for ``match="1mm"``) plus
    ``incidence``, sorted by descending incidence.
    """
    feat_lf, idcols, rep = _feature_frame(cohort, key, match, productive_only=productive_only,
                                          strip_allele=strip_allele, scope=scope, threads=threads)
    thr = min_incidence
    if min_incidence_frac is not None:
        n = feat_lf.select(pl.col(SAMPLE_ID).n_unique()).collect().item()
        thr = max(thr, math.ceil(min_incidence_frac * n))
    inc = (feat_lf.group_by(idcols).agg(pl.col(SAMPLE_ID).n_unique().alias("incidence"))
           .filter(pl.col("incidence") >= thr).collect(engine="streaming"))
    if rep is not None:
        inc = inc.join(rep, on="meta_id", how="left")
    # Tie-break on the feature key. `incidence` is a small integer, so a `max_features` cut almost
    # always lands *inside* a tie band; sorting on it alone (polars sort is multithreaded and
    # `maintain_order=False`) would pick a different subset on every call, making cooccurrence()
    # irreproducible run-to-run.
    return inc.sort(["incidence", *idcols], descending=[True] + [False] * len(idcols))


def _normalize_design(phenotype: pl.DataFrame | pl.LazyFrame, pheno_col: str | None,
                      level_col: str | None, stratum_col: str | None) -> pl.DataFrame:
    """Coerce a phenotype/design frame to the reserved ``_pos``/``_level``/``_stratum`` columns."""
    d = phenotype.lazy()
    cols = d.collect_schema().names()
    pos = "_pos" if "_pos" in cols else pheno_col
    if pos is None:
        raise ValueError("pass pheno_col (or a design frame with a '_pos' column)")
    sel = [pl.col(SAMPLE_ID), pl.col(pos).cast(pl.Boolean).alias("_pos")]
    sel.append((pl.col(level_col) if level_col else pl.col("_level") if "_level" in cols
                else pl.lit("_all")).cast(pl.String).alias("_level"))
    sel.append((pl.col(stratum_col) if stratum_col else pl.col("_stratum") if "_stratum" in cols
                else pl.lit("_all")).cast(pl.String).alias("_stratum"))
    return d.select(sel).drop_nulls("_pos").collect()


def _assemble(base: pl.DataFrame, a, b, c, d, present, n_pos, n_neg, tests, alternative,
              n_perm, seed, perm_present, perm_labels):
    """Long-format result: one row per (feature[/level], test) with shared + per-test columns."""
    a, b, c, d, present = (np.asarray(x, dtype=np.int64) for x in (a, b, c, d, present))
    shared = base.with_columns(
        pl.Series("incidence", present), pl.Series("n_pos_present", a),
        pl.Series("n_neg_present", b), pl.lit(n_pos, dtype=pl.Int64).alias("n_pos"),
        pl.lit(n_neg, dtype=pl.Int64).alias("n_neg"),
        pl.Series("odds_ratio", stats.odds_ratio(a, b, c, d)),
        pl.Series("log2_or", np.log2(stats.odds_ratio(a, b, c, d))),
        pl.Series("direction", stats.direction(a, b, c, d)))
    nul = pl.lit(None, dtype=pl.Float64)
    out = []
    for t in tests:
        cols = {"p_value": nul, "q_value": nul, "logor": nul, "logor_ci_lo": nul,
                "logor_ci_hi": nul, "p_or_gt1": nul, "log_bf": nul}
        if t in ("fisher", "chi2", "permutation"):
            if t == "fisher":
                p = stats.fisher_p(a, b, c, d, alternative=alternative)
            elif t == "chi2":
                p = stats.chi2_p(a, b, c, d)
            elif perm_present is not None and a.size:
                if alternative == "two-sided":
                    # A real two-sided permutation p (doubling convention), not a silent
                    # substitution of the upper tail. Same seed -> same permutation draws for
                    # both tails; 2*min(greater, less) clamped to 1.
                    pg = stats.permutation_p(perm_present, perm_labels, n_perm=n_perm, seed=seed,
                                             alternative="greater")
                    plt = stats.permutation_p(perm_present, perm_labels, n_perm=n_perm, seed=seed,
                                              alternative="less")
                    p = np.minimum(1.0, 2.0 * np.minimum(pg, plt))
                else:
                    p = stats.permutation_p(perm_present, perm_labels, n_perm=n_perm, seed=seed,
                                            alternative=alternative)
            else:
                p = np.zeros(a.size)
            cols["p_value"] = pl.Series(p, dtype=pl.Float64)
            cols["q_value"] = pl.Series(stats.fdr_bh(p), dtype=pl.Float64)
        elif t == "bayes_logodds":
            r = stats.bayes_logodds(a, b, c, d)
            cols.update(logor=pl.Series(r["logor"]), logor_ci_lo=pl.Series(r["ci_lo"]),
                        logor_ci_hi=pl.Series(r["ci_hi"]), p_or_gt1=pl.Series(r["p_or_gt1"]))
        elif t == "bayes_bf":
            cols["log_bf"] = pl.Series(stats.bayes_bf(a, b, c, d))
        out.append(shared.with_columns(pl.lit(t).alias("test"), **cols))
    return pl.concat(out)


def association(
    cohort: pl.LazyFrame | pl.DataFrame,
    phenotype: pl.DataFrame | pl.LazyFrame,
    *,
    pheno_col: str | None = None,
    level_col: str | None = None,
    stratum_col: str | None = None,
    test: str | list[str] = "fisher",
    key: tuple[str, ...] = DEFAULT_KEY,
    match: str = "exact",
    min_incidence: int = 2,
    min_incidence_frac: float | None = None,
    candidates: pl.DataFrame | None = None,
    alternative: str = "greater",
    productive_only: bool = True,
    strip_allele: bool = True,
    scope: str = "1,0,0,1",
    n_perm: int = 1000,
    seed: int = 0,
    threads: int = 0,
) -> pl.DataFrame:
    """Test each clonotype feature's subject incidence against a condition.

    Args:
        cohort: Clonotype cohort — a streamed :func:`vdjtools.io.scan_cohort` LazyFrame or a
            :class:`polars.DataFrame` with ``sample_id`` + the ``key`` columns.
        phenotype: A design frame — one row per subject (or subject × level) with ``sample_id``
            and either the reserved ``_pos``/``_level``/``_stratum`` columns (from
            :mod:`vdjtools.biomarker.condition`) or a plain binary ``pheno_col``.
        pheno_col: Binary phenotype column (if ``phenotype`` has no ``_pos``).
        level_col: Category-level column → one test per level (adds a ``level`` column).
        stratum_col: Stratum column → the tests are combined by Cochran–Mantel–Haenszel.
        test: One test or a list of ``{"fisher", "chi2", "bayes_logodds", "bayes_bf",
            "permutation"}`` (long output, one row per feature×level×test). Ignored when
            ``stratum_col`` is set (CMH is used).
        key: Feature key (subset of ``(junction_aa, v_call, j_call)``) — the V/J match scope.
        match: ``"exact"``, ``"fuzzy"`` or ``"1mm"``.

            * ``"exact"`` — the ``key`` itself.
            * ``"fuzzy"`` — a 1-mismatch **search**: a candidate's incidence counts every subject
              carrying anything within ``scope`` of it, but the candidate keeps its own identity.
              Non-``junction_aa`` key columns still have to match exactly, so
              ``key=("junction_aa", "v_call")`` pins the germline half while the CDR3 varies. This
              is the discovery mode.
            * ``"1mm"`` — metaclonotypes: candidates within ``scope`` are **merged** and the group
              is tested (``meta_id`` replaces the key). A downstream step, not discovery — merging
              dilutes a strong member into its neighbourhood.
        min_incidence: Minimum subjects a feature must appear in.
        min_incidence_frac: Alternative/added fraction-of-subjects threshold (e.g. ``0.05``).
        candidates: Restrict testing to these feature keys (a frame with the ``key`` columns).
        alternative: ``"greater"`` / ``"less"`` / ``"two-sided"`` (Fisher/permutation).
        n_perm, seed: Permutation settings.

    Returns:
        Long frame: feature key (or ``meta_id`` + representative key + ``n_members``), optional
        ``level``, then ``incidence, n_pos_present, n_neg_present, n_pos, n_neg, odds_ratio,
        log2_or, direction, test`` and the per-test statistics ``p_value, q_value`` (Fisher/χ²/
        permutation), ``logor, logor_ci_lo, logor_ci_hi, p_or_gt1`` (Bayesian log-odds), or
        ``log_bf`` (Beta-Binomial). CMH adds ``or_mh``. Sorted by ``p_value`` where present.
    """
    tests = [test] if isinstance(test, str) else list(test)
    for t in tests:
        if t not in _TESTS:
            raise ValueError(f"test must be a subset of {_TESTS}; got {t!r}")

    feat_lf, idcols, rep = _feature_frame(cohort, key, match, productive_only=productive_only,
                                          strip_allele=strip_allele, scope=scope,
                                          threads=threads, candidates=candidates)
    design = _normalize_design(phenotype, pheno_col, level_col, stratum_col)
    # Restrict the design to subjects we actually OBSERVED. n_pos/n_neg below are the arms of the
    # 2x2 (c = n_pos - a is "condition-positive subjects without the feature"), so a labelled
    # subject with no repertoire in `cohort` would silently vote "feature absent" for every
    # feature rather than being unobserved. On a metadata sheet that outruns the sequenced cohort
    # this is severe and one-sided: covid19 ships 1,212 labelled subjects for 572 repertoires, so
    # 438 of the 472 "negatives" (93%) were phantoms and nearly every feature looked enriched.
    seen = cohort.lazy().select(pl.col(SAMPLE_ID)).unique().collect()
    design = design.join(seen, on=SAMPLE_ID, how="semi")
    if not design.height:
        raise ValueError("no labelled subject appears in the cohort — check the sample_id join")
    n_labeled = design[SAMPLE_ID].n_unique()
    thr = min_incidence
    if min_incidence_frac is not None:
        thr = max(thr, math.ceil(min_incidence_frac * n_labeled))

    joined = feat_lf.join(design.lazy(), on=SAMPLE_ID, how="inner")
    if thr > 1:
        # Cheap superset prefilter (carried over from the v2.6.0 fisher_association it replaced):
        # a feature's distinct-subject incidence is <= its raw row count, so features with fewer
        # than `thr` rows cannot be public. Dropping them with a light pl.len() count before the
        # memory-heavier n_unique pass keeps peak RAM bounded on huge cohorts — the millions of
        # private clonotypes never build a hash-set. Correctness-preserving: the kept set is a
        # superset of the true public features, and the real threshold is applied below.
        public = (joined.group_by(idcols).agg(pl.len().alias("_rows"))
                  .filter(pl.col("_rows") >= thr).select(idcols))
        joined = joined.join(public, on=idcols, how="semi")
    stratified = design["_stratum"].n_unique() > 1
    gcols = [*idcols, "_level", "_stratum"]
    per = (joined.group_by(gcols)
           .agg(pl.col(SAMPLE_ID).n_unique().alias("present"),
                pl.col(SAMPLE_ID).filter(pl.col("_pos")).n_unique().alias("a"))
           .collect(engine="streaming"))
    tot = (design.group_by(["_level", "_stratum"])
           .agg(pl.col("_pos").sum().alias("n_pos_s"), pl.len().alias("n_s")))

    if not stratified:
        per = per.join(tot, on=["_level", "_stratum"], how="left")
        per = per.filter(pl.col("present") >= thr)
        res = _run_unstratified(per, idcols, rep, tests, alternative, n_perm, seed,
                                joined, design, thr)
    else:
        res = _run_cmh(per, tot, idcols, rep, thr)

    if res.height and "_level" in res.columns and res["_level"].n_unique() == 1 \
            and res["_level"][0] == "_all":
        res = res.drop("_level")
    elif "_level" in res.columns:
        res = res.rename({"_level": "level"})
    sort_c = "p_value" if "p_value" in res.columns else ("chi2" if "chi2" in res.columns else None)
    return res.sort(sort_c) if sort_c and res.height else res


def _run_unstratified(per, idcols, rep, tests, alternative, n_perm, seed, joined, design, thr):
    """Per (feature, level) 2×2 → the requested tests. Permutation gets a per-level matrix."""
    levels = design["_level"].unique().to_list()
    single_binary = levels == ["_all"]
    frames = []
    for (level,), sub in per.group_by(["_level"], maintain_order=True):
        row = design.filter(pl.col("_level") == level)
        n_pos = int(row["_pos"].sum())
        n_neg = row.height - n_pos
        if n_pos == 0 or n_neg == 0:
            if single_binary:
                raise ValueError("phenotype has only one class after dropping unknown labels")
            continue                                 # a level with no contrast can't be tested
        a = sub["a"].to_numpy().astype(np.int64)
        present = sub["present"].to_numpy().astype(np.int64)
        b, c, d = present - a, n_pos - a, n_neg - (present - a)
        base = sub.select([*idcols, "_level"])
        pp = pl_ = None
        if "permutation" in tests:
            pp, pl_ = _perm_matrix(joined, row, level, base, idcols)
        frames.append(_assemble(base, a, b, c, d, present, n_pos, n_neg, tests, alternative,
                                n_perm, seed, pp, pl_))
    if not frames:                                   # keep the full schema when nothing passes
        z = np.array([], dtype=np.int64)
        frames = [_assemble(per.head(0).select([*idcols, "_level"]), z, z, z, z, z,
                            0, 0, tests, alternative, n_perm, seed, None, None)]
    res = pl.concat(frames)
    if rep is not None:
        res = res.join(rep, on="meta_id", how="left")
    return res


def _perm_matrix(joined, row, level, base, idcols):
    """Boolean (n_subjects × n_features) incidence + label vector, feature order = ``base``."""
    feats = base.select(idcols).with_row_index("_fi")
    subs = row.select(SAMPLE_ID, "_pos").with_row_index("_si")
    inc = (joined.filter(pl.col("_level") == level).select([*idcols, SAMPLE_ID]).unique()
           .collect().join(feats, on=idcols, how="inner").join(subs, on=SAMPLE_ID, how="inner"))
    mat = np.zeros((subs.height, feats.height), dtype=bool)
    mat[inc["_si"].to_numpy(), inc["_fi"].to_numpy()] = True
    return mat, subs["_pos"].to_numpy().astype(bool)


def _run_cmh(per, tot, idcols, rep, thr):
    """Cochran–Mantel–Haenszel over strata, per (feature, level)."""
    strata = sorted(tot["_stratum"].unique().to_list())
    tot_m = {(r["_level"], r["_stratum"]): (r["n_pos_s"], r["n_s"] - r["n_pos_s"])
             for r in tot.iter_rows(named=True)}
    frames = []
    for (level,), sub in per.group_by(["_level"], maintain_order=True):
        feat = sub.select(idcols).unique()
        inc_tot = sub.group_by(idcols).agg(pl.col("present").sum().alias("incidence"))
        feat = feat.join(inc_tot, on=idcols).filter(pl.col("incidence") >= thr).select(idcols)
        if not feat.height:
            continue
        wide_a = _pivot(sub, feat, idcols, "a", strata)
        wide_p = _pivot(sub, feat, idcols, "present", strata)
        npos = np.array([tot_m.get((level, s), (0, 0))[0] for s in strata], dtype=np.float64)
        nneg = np.array([tot_m.get((level, s), (0, 0))[1] for s in strata], dtype=np.float64)
        A = wide_a
        B = wide_p - wide_a
        C = npos[None, :] - A
        D = nneg[None, :] - B
        r = stats.cmh(A, B, C, D)
        base = feat.with_columns(
            pl.lit(level).alias("_level"),
            pl.Series("incidence", wide_p.sum(1).astype(np.int64)),
            pl.Series("n_pos_present", A.sum(1).astype(np.int64)),
            pl.Series("n_neg_present", B.sum(1).astype(np.int64)),
            pl.Series("or_mh", r["or_mh"]), pl.Series("log2_or", np.log2(r["or_mh"])),
            pl.lit("cmh").alias("test"), pl.Series("chi2", r["chi2"]),
            pl.Series("p_value", r["p_value"]),
            pl.Series("q_value", stats.fdr_bh(r["p_value"])))
        frames.append(base)
    res = pl.concat(frames) if frames else per.head(0)
    if rep is not None:
        res = res.join(rep, on="meta_id", how="left")
    return res


def _pivot(sub, feat, idcols, value, strata):
    """(n_features × n_strata) dense array of ``value`` in ``feat`` order, missing → 0."""
    w = (feat.join(sub, on=idcols, how="left")
         .pivot(values=value, index=idcols, on="_stratum", aggregate_function="first"))
    for s in strata:
        if s not in w.columns:
            w = w.with_columns(pl.lit(0).alias(s))
    return w.select(strata).fill_null(0).to_numpy().astype(np.float64)
