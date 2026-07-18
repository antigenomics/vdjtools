"""V/J usage rescaling — swap a model's protocol-dependent usage, keep the junction model.

Usage is protocol-dependent (5'RACE vs DNA multiplex amplify different V's); the junction model
is not. These pin that rescale_usage moves ONLY the usage.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.model import load_bundled, rescale_usage


def _sample(v_calls, j_calls=None) -> pl.DataFrame:
    j_calls = j_calls or ["TRBJ2-7*01"] * len(v_calls)
    return pl.DataFrame({
        "junction_aa": [f"CASS{i}EQYF" for i in range(len(v_calls))],
        "junction_nt": ["ACG"] * len(v_calls),
        "v_call": v_calls, "j_call": j_calls,
        "duplicate_count": [1] * len(v_calls), "frequency": [1.0 / len(v_calls)] * len(v_calls),
    })


def _gene_mass(model, table, col):
    t = model.tables[table].with_columns(pl.col(col).str.split("*").list.first().alias("g"))
    return {r["g"]: r["p"] for r in t.group_by("g").agg(pl.col("p").sum()).iter_rows(named=True)}


def test_rescale_sets_v_usage_to_the_sample():
    m = load_bundled("TRB", "olga")
    s = _sample(["TRBV19*01"] * 75 + ["TRBV20-1*01"] * 25)
    r = rescale_usage(m, s)
    g = _gene_mass(r, "v_choice", "v_allele")
    assert g["TRBV19"] == pytest.approx(0.75)
    assert g["TRBV20-1"] == pytest.approx(0.25)
    assert sum(g.values()) == pytest.approx(1.0)
    assert g.get("TRBV5-1", 0.0) == 0.0          # absent from the sample -> no mass


def test_rescale_leaves_the_junction_model_untouched():
    """The whole point: usage moves, the recombination machinery does not."""
    m = load_bundled("TRB", "olga")
    r = rescale_usage(m, _sample(["TRBV19*01"] * 10 + ["TRBV20-1*01"] * 10))
    for ev in ("v_3_del", "j_5_del", "d_del", "vd_ins", "dj_ins", "vd_dinucl", "dj_dinucl",
               "d_gene", "n_d"):
        assert r.tables[ev].equals(m.tables[ev]), f"{ev} changed — rescale touched the junction model"


def test_rescale_preserves_within_gene_allele_split():
    """A gene's new mass is split across its alleles in the model's existing proportions.

    The sample cannot resolve alleles (mismapping on short reads), so the model's own split is
    the best available and must survive.
    """
    m = load_bundled("TRB", "olga", collapse=False)      # this test is specifically about alleles
    v = m.tables["v_choice"].with_columns(pl.col("v_allele").str.split("*").list.first().alias("g"))
    multi = [r["g"] for r in v.group_by("g").len().iter_rows(named=True) if r["len"] > 1]
    gene = next(g for g in multi if v.filter(pl.col("g") == g)["p"].sum() > 0)
    before = v.filter(pl.col("g") == gene)
    ratio_before = (before["p"] / before["p"].sum()).to_list()

    r = rescale_usage(m, _sample([f"{gene}*01"] * 10))
    after = r.tables["v_choice"].with_columns(pl.col("v_allele").str.split("*").list.first().alias("g")).filter(pl.col("g") == gene)
    assert after["p"].sum() == pytest.approx(1.0)
    assert (after["p"] / after["p"].sum()).to_list() == pytest.approx(ratio_before)


def test_rescale_j_and_selective_flags():
    m = load_bundled("TRB", "olga")
    s = _sample(["TRBV19*01"] * 10, ["TRBJ2-7*01"] * 6 + ["TRBJ1-1*01"] * 4)
    r = rescale_usage(m, s)
    gj = _gene_mass(r, "j_choice", "j_allele")
    assert gj["TRBJ2-7"] == pytest.approx(0.6)
    assert gj["TRBJ1-1"] == pytest.approx(0.4)
    # v=False leaves v_choice alone
    assert rescale_usage(m, s, v=False).tables["v_choice"].equals(m.tables["v_choice"])
    assert rescale_usage(m, s, j=False).tables["j_choice"].equals(m.tables["j_choice"])


def test_rescale_rejects_an_unusable_sample():
    m = load_bundled("TRB", "olga")
    # An all-null V column is genuinely unusable.
    s = _sample(["TRBV6-2*01"] * 3).with_columns(pl.lit(None, dtype=pl.String).alias("v_call"))
    with pytest.raises(ValueError, match="no usable"):
        rescale_usage(m, s)


def test_ambiguous_only_sample_is_usable_via_fractional_split():
    """An all-ambiguous sample is NOT unusable — the ties still carry gene information."""
    m = load_bundled("TRB", "olga")
    r = rescale_usage(m, _sample(["TRBV6-2*01,TRBV6-3*01"] * 5))    # was: raised; now: 0.5/0.5
    g = (r.tables["v_choice"].with_columns(pl.col("v_allele").str.split("*").list.first().alias("g"))
         .group_by("g").agg(pl.col("p").sum()))
    mass = {row["g"]: row["p"] for row in g.iter_rows(named=True)}
    assert mass.get("TRBV6-2", 0) == pytest.approx(0.5)
    assert mass.get("TRBV6-3", 0) == pytest.approx(0.5)


def test_ambiguous_calls_are_split_fractionally_not_dropped():
    """A comma tie is allocated 1/k per named gene, not discarded (duplicated loci tie constantly)."""
    from vdjtools.model.rescale import _empirical

    # 4 clonotypes: 2 unambiguous IGHV3-23, 2 ties naming IGHV3-23 + IGHV3-23D.
    s = pl.DataFrame({
        "junction_aa": [f"CAS{i}F" for i in range(4)],
        "junction_nt": ["ACG"] * 4,
        "v_call": ["IGHV3-23*01", "IGHV3-23*01",
                   "IGHV3-23*01,IGHV3-23D*01", "IGHV3-23*01,IGHV3-23D*01"],
        "j_call": ["IGHJ4*02"] * 4,
        "duplicate_count": [1] * 4, "frequency": [0.25] * 4,
    })
    p = _empirical(s, "v_call")
    # votes: IGHV3-23 = 2 + 0.5 + 0.5 = 3 ; IGHV3-23D = 0.5 + 0.5 = 1 ; total 4
    assert p["IGHV3-23"] == pytest.approx(0.75)
    assert p["IGHV3-23D"] == pytest.approx(0.25)
    assert sum(p.values()) == pytest.approx(1.0)
    # dropping the ties (the old behaviour) would have given IGHV3-23=1.0 and hidden IGHV3-23D
    assert "IGHV3-23D" in p, "the duplicated-locus paralog was dropped"


def test_rescale_transfers_usage_across_protocols():
    """Cross-protocol transfer: DNA-multiplex usage onto the 5'RACE junction model.

    Our learned models are 5'RACE (isalgo/airr_model_read); OLGA is DNA-multiplex. V/J usage is
    protocol-dependent, the junction (recombination) model is not — so ``rescale_usage`` must take a
    DNA-protocol repertoire's usage and keep the 5'RACE deletions/insertions/dinucleotides untouched.
    Shown on loci both protocols cover; gamma/delta (TRG/TRD) are ours-only (no DNA-based OLGA model),
    so they are excluded from the cross-protocol comparison.
    """
    import numpy as np

    from vdjtools.model.generate import generate

    def gene_usage(m):
        return {r["g"]: r["p"] for r in
                m.tables["v_choice"].with_columns(pl.col("v_allele").str.split("*").list.first().alias("g"))
                .group_by("g").agg(pl.col("p").sum()).iter_rows(named=True)}

    for locus, junction_events in (
        ("TRB", ["v_3_del", "j_5_del", "d_del", "vd_ins", "dj_ins", "vd_dinucl", "dj_dinucl", "d_gene", "n_d"]),
        ("IGK", ["v_3_del", "j_5_del", "vj_ins", "vj_dinucl"]),
    ):
        race = load_bundled(locus, "learned")            # 5'RACE junction model (our lab)
        dna = load_bundled(locus, "olga")                # DNA-multiplex usage
        dna_lib = generate(dna, 12000, seed=1)           # a DNA-protocol repertoire
        out = rescale_usage(race, dna_lib)

        # the recombination machinery stays the 5'RACE model's, byte for byte
        for ev in junction_events:
            assert out.tables[ev].equals(race.tables[ev]), f"{locus} {ev} changed — rescale touched the junction model"

        # ...while V usage has moved to the DNA protocol: closer to DNA than to the 5'RACE original
        u_out, u_dna, u_race = gene_usage(out), gene_usage(dna), gene_usage(race)
        genes = sorted(set(u_dna) | set(u_race))
        a = np.array([u_out.get(g, 0.0) for g in genes])
        d = np.array([u_dna.get(g, 0.0) for g in genes])
        r = np.array([u_race.get(g, 0.0) for g in genes])
        assert np.corrcoef(a, d)[0, 1] > 0.9, f"{locus}: rescaled usage does not match the DNA sample"
        assert np.corrcoef(a, d)[0, 1] > np.corrcoef(a, r)[0, 1], f"{locus}: usage did not move toward DNA"


def test_rescale_to_a_dataset_pool_vs_mean():
    """Rescaling to a whole dataset: 'pool' weights by depth, 'mean' weights samples equally."""
    m = load_bundled("TRB", "olga")
    deep = _sample(["TRBV19*01"] * 100)          # 100 clonotypes, all TRBV19
    shallow = _sample(["TRBV20-1*01"] * 10)      # 10 clonotypes, all TRBV20-1

    pooled = _gene_mass(rescale_usage(m, [deep, shallow], aggregate="pool"), "v_choice", "v_allele")
    meaned = _gene_mass(rescale_usage(m, [deep, shallow], aggregate="mean"), "v_choice", "v_allele")

    # pool: 100/110 vs 10/110 — the deep sample dominates
    assert pooled["TRBV19"] == pytest.approx(100 / 110, abs=1e-6)
    assert pooled["TRBV20-1"] == pytest.approx(10 / 110, abs=1e-6)
    # mean: each sample is one distribution averaged equally -> 50/50 regardless of depth
    assert meaned["TRBV19"] == pytest.approx(0.5)
    assert meaned["TRBV20-1"] == pytest.approx(0.5)
    # a one-element dataset equals passing that frame directly
    one = _gene_mass(rescale_usage(m, [deep]), "v_choice", "v_allele")
    direct = _gene_mass(rescale_usage(m, deep), "v_choice", "v_allele")
    assert one == pytest.approx(direct)


def test_rescale_streams_a_dataset_without_materializing_it():
    """The dataset is consumed one sample at a time — a one-shot generator and LazyFrames both work.

    Rescaling to a whole dataset must not require holding it in memory: V and J usage are accumulated
    in a single streaming pass, so a generator (usable once) is not exhausted by computing V then J,
    and per-sample ``LazyFrame``s (``scan_parquet``) stream each sample's group-by straight from disk.
    """
    m = load_bundled("TRB", "olga")
    samples = [_sample(["TRBV19*01"] * 30, ["TRBJ2-7*01"] * 30),
               _sample(["TRBV20-1*01"] * 30, ["TRBJ1-1*01"] * 30)]
    ref = rescale_usage(m, samples)                              # list baseline

    # a one-shot generator: rescaling BOTH v and j must not exhaust it mid-way
    gen = (s for s in samples)
    out_gen = rescale_usage(m, gen, v=True, j=True)
    assert _gene_mass(out_gen, "v_choice", "v_allele") == pytest.approx(_gene_mass(ref, "v_choice", "v_allele"))
    assert _gene_mass(out_gen, "j_choice", "j_allele") == pytest.approx(_gene_mass(ref, "j_choice", "j_allele"))

    # per-sample LazyFrames (as scan_parquet would yield) — never materialized as one frame
    out_lazy = rescale_usage(m, (s.lazy() for s in samples))
    assert _gene_mass(out_lazy, "v_choice", "v_allele") == pytest.approx(_gene_mass(ref, "v_choice", "v_allele"))
