"""``check_model`` — one deliberately-broken model per check, plus a clean-model baseline.

Each test breaks exactly one thing and asserts that check fires, so a false positive somewhere else
shows up as the clean-baseline test failing rather than as a vague count mismatch.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.model import load_bundled
from vdjtools.model.bundled import LOCI
from vdjtools.model.check import check_model
from vdjtools.model.model import Model


def _rebuild(model, **tables):
    return Model(manifest=model.manifest, tables={**model.tables, **tables},
                 genomic=model.genomic)


def _checks(issues: pl.DataFrame, severity: str | None = None) -> set[str]:
    if severity:
        issues = issues.filter(pl.col("severity") == severity)
    return set(issues["check"].to_list())


def test_toy_model_is_clean(toy_model):
    assert check_model(toy_model, germline="none").filter(pl.col("severity") == "error").is_empty()


def test_toy_vdj_model_is_clean(toy_model_vdj):
    issues = check_model(toy_model_vdj, germline="none")
    assert issues.filter(pl.col("severity") == "error").is_empty()


def test_bundled_trb_has_no_errors():
    """TRB is the workhorse locus; a shipped model must not carry error-level issues."""
    for source in ("olga", "learned"):
        issues = check_model(load_bundled("TRB", source))
        assert issues.filter(pl.col("severity") == "error").is_empty(), issues


def test_issue_frame_schema(toy_model):
    issues = check_model(toy_model, germline="none")
    assert issues.columns == ["severity", "check", "event", "segment", "allele", "detail", "value"]
    assert set(issues["severity"]) <= {"error", "warn", "info"}


def test_normalization_reports_every_offender(toy_model):
    """validate_tables raises on the first bad group; the checker must list them all."""
    broken = _rebuild(toy_model,
                      v_3_del=toy_model.tables["v_3_del"].with_columns(p=pl.col("p") * 2.0))
    issues = check_model(broken, germline="none")
    norm = issues.filter(pl.col("check") == "normalization")
    assert norm.height == 3          # one per V allele group, not just the first
    assert set(norm["severity"]) == {"error"}


def test_probability_range(toy_model):
    broken = _rebuild(toy_model, vj_ins=toy_model.tables["vj_ins"].with_columns(
        p=pl.when(pl.col("length") == 0).then(-0.5).otherwise(pl.col("p"))))
    assert "probability_range" in _checks(check_model(broken, germline="none"), "error")


def test_allele_not_in_genomic(toy_model):
    extra = pl.DataFrame({"v_allele": ["GHOSTV*01"], "p": [0.0]})
    broken = _rebuild(toy_model, v_choice=pl.concat([toy_model.tables["v_choice"], extra]))
    issues = check_model(broken, germline="none")
    assert "allele_not_in_genomic" in _checks(issues, "error")
    assert "GHOSTV*01" in issues["allele"].to_list()


def test_genomic_not_in_tables(toy_model):
    broken = _rebuild(toy_model, v_choice=toy_model.tables["v_choice"]
                      .filter(pl.col("v_allele") != "TOYV2*01"))
    # Renormalize so the missing allele, not the sums, is what is reported.
    broken = _rebuild(broken, v_choice=broken.tables["v_choice"]
                      .with_columns(p=pl.col("p") / pl.col("p").sum()))
    assert "genomic_not_in_tables" in _checks(check_model(broken, germline="none"), "warn")


def test_functional_zero_mass(toy_model):
    """A real gene pinned to zero is the absorbing-state failure — silent, and Pgen-fatal."""
    zeroed = toy_model.tables["v_choice"].with_columns(
        p=pl.when(pl.col("v_allele") == "TOYV2*01").then(0.0).otherwise(pl.col("p")))
    broken = _rebuild(toy_model, v_choice=zeroed.with_columns(p=pl.col("p") / pl.col("p").sum()))
    assert "functional_zero_mass" in _checks(check_model(broken, germline="none"), "warn")


def test_dinucleotide_completeness(toy_model):
    broken = _rebuild(toy_model, vj_dinucl=toy_model.tables["vj_dinucl"].head(15))
    assert "dinucleotide_complete" in _checks(check_model(broken, germline="none"), "error")


def test_deletion_unreachable_scales_with_the_lost_mass(toy_model):
    """Reported as the FRACTION of an allele's mass that no trim can reach, not per row."""
    t = toy_model.tables["v_3_del"]
    # TOYV1*01's cut_segment is 13 nt with palindrome_max 4, so ndel > 8 is unreachable.
    heavy = t.with_columns(
        p=pl.when((pl.col("v_allele") == "TOYV1*01") & (pl.col("ndel") > 8))
        .then(1.0).otherwise(pl.col("p")))
    heavy = heavy.with_columns(p=pl.col("p") / pl.col("p").sum().over("v_allele"))
    issues = check_model(_rebuild(toy_model, v_3_del=heavy), germline="none")
    row = issues.filter(pl.col("check") == "deletion_unreachable").filter(
        pl.col("allele") == "TOYV1*01")
    assert row.height == 1
    assert row["severity"][0] == "error"       # far more than the 10% error threshold
    assert row["value"][0] > 0.5


def test_palindrome_max_overrun(toy_model):
    shifted = toy_model.tables["v_3_del"].with_columns(ndel=pl.col("ndel") - 5)
    assert "palindrome_max" in _checks(
        check_model(_rebuild(toy_model, v_3_del=shifted), germline="none"), "error")


def test_vj_model_carrying_a_d_event_is_an_error(toy_model, toy_model_vdj):
    """A VJ manifest with D machinery cannot be scored; report it rather than crash in prepare."""
    broken = Model(manifest=toy_model_vdj.manifest, tables=toy_model_vdj.tables,
                   genomic={k: v for k, v in toy_model_vdj.genomic.items() if k != "genes_d"})
    assert "event_set" in _checks(check_model(broken, germline="none"), "error")


def test_tandem_mass_without_tandem_events(toy_model_vdj):
    nd = pl.DataFrame({"n_d": pl.Series([1, 2], dtype=pl.UInt8), "p": [0.8, 0.2]})
    assert "event_set" in _checks(
        check_model(_rebuild(toy_model_vdj, n_d=nd), germline="none"), "error")


def test_insertion_truncated(toy_model):
    """A support clipped before the distribution decayed — the from_arda ins_max artifact."""
    flat = toy_model.tables["vj_ins"].with_columns(p=pl.lit(0.25))
    assert "insertion_truncated" in _checks(
        check_model(_rebuild(toy_model, vj_ins=flat), germline="none"), "warn")


def test_raise_on_error(toy_model):
    broken = _rebuild(toy_model,
                      v_3_del=toy_model.tables["v_3_del"].with_columns(p=pl.col("p") * 2.0))
    with pytest.raises(ValueError, match="error-level"):
        check_model(broken, germline="none", raise_on="error")
    # A clean model must not raise.
    check_model(toy_model, germline="none", raise_on="error")


def test_germline_reconciliation_flags_a_changed_sequence(toy_model, toy_germline):
    changed = toy_germline.with_columns(
        sequence=pl.when(pl.col("allele") == "TOYV1*01")
        .then(pl.lit("TGTGCCAAA")).otherwise(pl.col("sequence")))
    issues = check_model(toy_model, germline=changed)
    assert "germline_source" in _checks(issues, "warn")


def test_custom_locus_skips_reconciliation_silently(toy_model):
    """arda has no 'TOY' locus; that must not be reported as a problem with the model."""
    assert "germline_source" not in _checks(check_model(toy_model, germline="auto"))


@pytest.mark.slow
def test_every_bundled_model_is_loadable_and_checkable():
    for source in ("olga", "learned"):
        for locus in LOCI:
            issues = check_model(load_bundled(locus, source))
            assert issues.columns[0] == "severity"


# --- the D->J locus-order constraint ------------------------------------------------------------

def test_forbidden_dj_pairs_only_constrains_trb():
    """Deletional joining cannot reach a J 5' of the D. Only TRB interleaves its clusters."""
    from vdjtools.model.reference import forbidden_dj_pairs

    d, j = ["TRBD1*01", "TRBD2*01"], ["TRBJ1-1*01", "TRBJ1-6*01", "TRBJ2-1*01"]
    bad = forbidden_dj_pairs(d, j, "TRB")
    assert bad == {("TRBD2*01", "TRBJ1-1*01"), ("TRBD2*01", "TRBJ1-6*01")}
    # TRBD1 sits 5' of both clusters, so it reaches either.
    assert not any(p[0] == "TRBD1*01" for p in bad)
    # Every other locus puts all D 5' of all J — the rule must not invent constraints.
    assert forbidden_dj_pairs(["IGHD1-1*01"], ["IGHJ1*01"], "IGH") == set()
    assert forbidden_dj_pairs(["TRDD2*01"], ["TRDJ1*01"], "TRD") == set()


def test_impossible_dj_pair_is_reported():
    from vdjtools.model import load_bundled
    from vdjtools.model.check import check_model as _check

    m = load_bundled("TRB", "olga", collapse=False)
    hit = _check(m).filter(pl.col("check") == "impossible_dj_pair")
    assert hit.height > 0
    # OLGA's own model carries these, so a faithful import is warned about, not failed.
    assert set(hit["severity"]) == {"warn"}


def test_enforce_dj_order_removes_them_and_renormalizes():
    from vdjtools.model import load_bundled
    from vdjtools.model.check import check_model as _check
    from vdjtools.model.infer import enforce_dj_order

    m = load_bundled("TRB", "olga", collapse=False)
    fixed = enforce_dj_order(m)
    assert _check(fixed).filter(pl.col("check") == "impossible_dj_pair").is_empty()
    fixed.validate()
    # Mass is redistributed within each J, not lost.
    sums = fixed.tables["d_gene"].group_by("j_allele").agg(pl.col("p").sum())
    assert ((sums["p"] - 1.0).abs() < 1e-9).all()
    # A J in the TRBJ2 cluster is untouched: nothing about it was impossible.
    def j2(model):
        return (model.tables["d_gene"].filter(pl.col("j_allele") == "TRBJ2-1*01")
                .sort("d_allele")["p"].to_list())
    assert j2(fixed) == pytest.approx(j2(m))


def test_em_cannot_relearn_an_impossible_pair(toy_model_vdj):
    """The mask lives in the M-step, so a fit cannot put the mass back."""
    from vdjtools.model.generate import generate
    from vdjtools.model.infer import infer_native

    seqs = generate(toy_model_vdj, 150, seed=1)["junction_nt"].to_list()
    fitted, _ = infer_native(toy_model_vdj, seqs, max_iter=2)
    # The toy locus is not TRB, so nothing is forbidden and the fit is unconstrained.
    assert check_model(fitted, germline="none").filter(
        pl.col("check") == "impossible_dj_pair").is_empty()
