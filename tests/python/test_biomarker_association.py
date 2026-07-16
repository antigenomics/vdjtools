"""Tests for vdjtools.biomarker.association — tests, category expansion, CMH, candidates."""
import numpy as np
import polars as pl
import pytest

from vdjtools.biomarker import association, condition, select_candidates


def _clono(sid, v, j, aa):
    return dict(sample_id=sid, v_call=v, j_call=j, junction_aa=aa,
                junction_nt="ACG", duplicate_count=10, frequency=0.0)


def _cohort(seed=0):
    """40 subjects (20 CMV+, 20 CMV−); CASSXF enriched, CASSZF depleted, CASSBG background."""
    rng = np.random.default_rng(seed)
    rows, meta = [], []
    for i in range(40):
        sid, pos = f"S{i:02d}", i < 20
        meta.append(dict(sample_id=sid, cmv="+" if pos else "-",
                         hla="A*02" if i % 2 == 0 else "A*01"))
        rows.append(_clono(sid, "TRBV1", "TRBJ1", "CASSBG"))                 # ubiquitous
        if (pos and rng.random() < 0.85) or (not pos and rng.random() < 0.1):
            rows.append(_clono(sid, "TRBV1", "TRBJ1", "CASSXF"))            # enriched in +
        if (not pos and rng.random() < 0.85) or (pos and rng.random() < 0.1):
            rows.append(_clono(sid, "TRBV4", "TRBJ4", "CASSZF"))           # depleted in +
    return pl.DataFrame(rows), pl.DataFrame(meta)


def test_all_tests_agree_on_enriched_feature():
    cohort, meta = _cohort()
    ph = condition.binary(meta, "cmv")
    r = association(cohort, ph, test=["fisher", "chi2", "bayes_logodds", "bayes_bf", "permutation"],
                    n_perm=1000, seed=1)
    x = r.filter(pl.col("junction_aa") == "CASSXF")
    assert (x.filter(pl.col("test") == "fisher")["direction"][0]) == "enriched"
    for t in ("fisher", "chi2", "permutation"):
        assert x.filter(pl.col("test") == t)["p_value"][0] < 0.01
    assert x.filter(pl.col("test") == "bayes_bf")["log_bf"][0] > 3        # strong evidence
    assert x.filter(pl.col("test") == "bayes_logodds")["p_or_gt1"][0] > 0.97
    # background clonotype: not significant on any frequentist test
    bg = r.filter((pl.col("junction_aa") == "CASSBG") & (pl.col("test") == "fisher"))
    assert bg.height == 0 or bg["p_value"][0] > 0.2


def test_depleted_two_sided():
    cohort, meta = _cohort()
    ph = condition.binary(meta, "cmv")
    r = association(cohort, ph, test="fisher", alternative="two-sided")
    z = r.filter(pl.col("junction_aa") == "CASSZF")
    assert z["direction"][0] == "depleted" and z["log2_or"][0] < 0 and z["p_value"][0] < 0.05


def test_categorical_hla_alleles_expands_per_allele():
    cohort, meta = _cohort()
    des = condition.hla_alleles(meta, ["hla"], min_level_size=5)
    r = association(cohort, des, level_col="_level", pheno_col="_pos", test="fisher")
    assert set(r["level"].unique().to_list()) == {"A*01", "A*02"}
    assert "level" in r.columns

    zyg = condition.zygosity(meta, ("hla", "hla"))                          # all homozygous here
    assert set(zyg["_pos"].to_list()) == {True}


def test_hla_alleles_drops_untyped_subjects_rather_than_calling_them_non_carriers():
    """An HLA-untyped subject must be excluded, not counted as a non-carrier of every allele.

    Counting untyped subjects as non-carriers inflates the negative arm and biases the odds
    ratio anticonservatively (regression: it silently added 106 phantom non-carriers per
    allele on the 572-subject covid19 cohort).
    """
    _, meta = _cohort()
    typed = set(meta.filter(pl.col("hla").is_not_null())["sample_id"].to_list())
    untyped = meta.head(3).with_columns(
        (pl.col("sample_id") + "_untyped").alias("sample_id"),
        pl.lit(None, dtype=pl.String).alias("hla"))          # no typing at all
    des = condition.hla_alleles(pl.concat([meta, untyped]), ["hla"], min_level_size=5)

    got = set(des["sample_id"].unique().to_list())
    assert got == typed, "untyped subjects must not appear in the design at any level"
    assert not any(s.endswith("_untyped") for s in got)
    # Empty-string / NA spellings are normalised to null and dropped the same way.
    blank = meta.head(2).with_columns((pl.col("sample_id") + "_blank").alias("sample_id"),
                                      pl.lit("NA").alias("hla"))
    des2 = condition.hla_alleles(pl.concat([meta, blank]), ["hla"], min_level_size=5)
    assert not any(s.endswith("_blank") for s in des2["sample_id"].unique().to_list())


def test_cmh_controls_a_confounded_association():
    # HLA is a confounder: A*02+ subjects are both more often CMV+ AND carry CASSHLA more.
    # The marginal Fisher is significant; after stratifying by HLA (CMH) it vanishes.
    rows, meta = [], []
    sid = 0
    for si, (hla, pos_rate, carry_rate, n) in enumerate(
            [("A*02", 0.85, 0.85, 40), ("A*01", 0.15, 0.15, 40)]):
        rng = np.random.default_rng(si)                 # deterministic seed (not hash())
        for _ in range(n):
            s = f"S{sid:03d}"
            sid += 1
            cmv = rng.random() < pos_rate
            meta.append(dict(sample_id=s, cmv="+" if cmv else "-", hla=hla))
            rows.append(_clono(s, "TRBV1", "TRBJ1", "CASSBG"))
            if rng.random() < carry_rate:                                   # tied to HLA, not CMV
                rows.append(_clono(s, "TRBV1", "TRBJ1", "CASSHLA"))
    cohort, meta = pl.DataFrame(rows), pl.DataFrame(meta)
    marg = association(cohort, condition.binary(meta, "cmv"), test="fisher", min_incidence=5)
    strat = association(cohort, condition.stratified(meta, "cmv", "hla"),
                        stratum_col="_stratum", min_incidence=5)
    p_marg = marg.filter(pl.col("junction_aa") == "CASSHLA")["p_value"][0]
    p_cmh = strat.filter(pl.col("junction_aa") == "CASSHLA")["p_value"][0]
    assert p_marg < 0.05                                                    # confounded signal
    assert p_cmh > 0.2                                                      # gone after stratifying
    assert "or_mh" in strat.columns


def test_min_incidence_frac_and_candidates():
    cohort, meta = _cohort()
    ph = condition.binary(meta, "cmv")
    cand = select_candidates(cohort, min_incidence_frac=0.5)                # >=20/40 subjects
    assert cand["junction_aa"].to_list() == ["CASSBG"]                      # only the ubiquitous one
    r = association(cohort, ph, candidates=cand.select("junction_aa", "v_call", "j_call"))
    assert set(r["junction_aa"].to_list()) <= {"CASSBG"}


def test_single_class_raises_and_unknown_labels_dropped():
    cohort, meta = _cohort()
    meta1 = meta.with_columns(pl.lit("+").alias("cmv"))                     # only one class
    with pytest.raises(ValueError, match="only one class"):
        association(cohort, condition.binary(meta1, "cmv"))
    # unknown labels are dropped, shrinking the tested cohort
    meta2 = meta.with_columns(pl.when(pl.col("sample_id") == "S00").then(pl.lit("?"))
                              .otherwise(pl.col("cmv")).alias("cmv"))
    r = association(cohort, condition.binary(meta2, "cmv"))
    assert int(r["n_pos"][0] + r["n_neg"][0]) == 39
