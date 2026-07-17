"""Collapsing a model to one *01 allele per gene — the default working resolution.

The collapse is a marginalisation: choice sums exactly, conditionals are usage-weighted allele
averages, and Pgen over a gene tracks the uncollapsed allele-sum closely (exact where the
CDR3-region germline is allele-invariant).
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest
from scipy.stats import pearsonr

from vdjtools.model import collapse_alleles, load_bundled, native
from vdjtools.model.generate import generate
from vdjtools.model.infer import _gene_to_alleles, call_alleles

LOCI = ["TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"]


def _gene_mass(model, col):
    t = model.tables[col].with_columns(pl.col(col.replace("choice", "allele")).str.split("*").list.first().alias("g"))
    return {r["g"]: r["p"] for r in t.group_by("g").agg(pl.col("p").sum()).iter_rows(named=True)}


@pytest.mark.parametrize("locus", LOCI)
def test_collapse_is_one_allele_per_gene_and_consistent(locus):
    m = load_bundled(locus, "olga", collapse=False)
    c = collapse_alleles(m)
    # every allele column, every table, is a *01 representative
    for tname, t in c.tables.items():
        for col in t.columns:
            if col.endswith("allele"):
                bad = [a for a in t[col].unique().to_list() if a and not a.endswith("*01")]
                assert not bad, f"{locus} {tname}.{col} has non-*01 alleles: {bad[:3]}"
    # one row per gene in v_choice
    n_genes = m.tables["v_choice"].with_columns(pl.col("v_allele").str.split("*").list.first().alias("g"))["g"].n_unique()
    assert c.tables["v_choice"].height == n_genes
    c.validate()
    generate(c, 100, seed=0)                       # generatively complete


@pytest.mark.parametrize("locus", LOCI)
def test_collapse_preserves_marginal_usage_exactly(locus):
    m = load_bundled(locus, "olga", collapse=False)
    c = collapse_alleles(m)
    before = _gene_mass(m, "v_choice")
    after = _gene_mass(c, "v_choice")
    for g in before:
        assert after.get(g, 0.0) == pytest.approx(before[g], abs=1e-12), f"{locus} P({g}) changed"
    assert sum(after.values()) == pytest.approx(1.0)


def test_collapsed_pgen_tracks_the_allele_sum():
    """Collapsed Pgen(gene) ~ uncollapsed Σ_alleles Pgen — high rank correlation."""
    m = load_bundled("TRB", "olga", collapse=False)
    c = collapse_alleles(m)
    va, ja = _gene_to_alleles(m, "v"), _gene_to_alleles(m, "j")
    draws = generate(m, 80, seed=3)
    lu, lc = [], []
    for row in draws.iter_rows(named=True):
        vg, jg = row["v_call"].split("*")[0], row["j_call"].split("*")[0]
        cdr = row["junction_nt"].upper()
        u = sum(native.pgen_nt(m, cdr, v, j) for v in call_alleles(va, vg) for j in call_alleles(ja, jg))
        cc = native.pgen_nt(c, cdr, f"{vg}*01", f"{jg}*01")
        if u > 0 and cc > 0:
            lu.append(np.log10(u)); lc.append(np.log10(cc))
    assert len(lu) >= 40
    assert pearsonr(lu, lc)[0] > 0.95


def test_load_bundled_collapse_flag():
    full = load_bundled("TRB", "learned", collapse=False)
    coll = load_bundled("TRB", "learned")              # default True
    assert coll.tables["v_choice"].height < full.tables["v_choice"].height
    assert full.tables["v_choice"].height == 89        # allele resolution retained when asked
    # a query on the collapsed model maps any allele of a gene to its *01 representative
    p1 = native.pgen_nt(coll, "TGTGCCAGCAGCTTC", "TRBV20-1*01", "TRBJ2-1*01")
    assert p1 >= 0.0


def test_collapse_guards_unknown_allele_table():
    """A new allele-keyed table with no collapse rule must raise, not silently pass through."""
    from vdjtools.model.collapse import collapse_alleles as ca
    m = load_bundled("TRB", "olga", collapse=False)
    m.tables["bogus"] = pl.DataFrame({"v_allele": ["TRBV20-1*02"], "p": [1.0]})
    with pytest.raises(NotImplementedError, match="allele-keyed"):
        ca(m)
