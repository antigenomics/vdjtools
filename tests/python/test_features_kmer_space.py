"""The frozen (V x k-mer) feature space.

What is protected here is mostly *frozenness*: a vocabulary, an IDF and a rotation that move
between two runs are worse than useless, because the numbers still look fine. The native kernel
is checked against the polars explode path it replaces rather than against itself.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjtools.features.kmer import _explode_kmers
from vdjtools.features.kmer_space import (
    AMINO_ACIDS,
    alphabet_table,
    fit_kmer_space,
    pattern_of,
    reduced_alphabet,
)

AA = list(AMINO_ACIDS)


def sample(n: int, motif: str | None = None, seed: int = 0) -> pl.DataFrame:
    r = np.random.default_rng(seed)
    js = []
    for i in range(n):
        s = "".join(r.choice(AA, r.integers(10, 18)))
        if motif and i % 5 == 0:
            p = int(r.integers(4, max(5, len(s) - len(motif) - 4)))
            s = s[:p] + motif + s[p + len(motif):]
        js.append("C" + s + "F")
    return pl.DataFrame({
        "junction_aa": js,
        "v_call": [f"TRBV{1 + int(x) % 20}" for x in r.integers(0, 20, n)],
        "j_call": ["TRBJ2-2"] * n,
        "duplicate_count": r.integers(1, 50, n).tolist(),
    }).with_columns((pl.col("duplicate_count") / pl.col("duplicate_count").sum()).alias("frequency"))


@pytest.fixture(scope="module")
def cohort():
    return [sample(400, "WGGD", seed=100 + i) for i in range(20)] + \
           [sample(400, None, seed=200 + i) for i in range(20)]


class TestReducedAlphabet:
    def test_groups_are_chemically_coherent(self):
        """BLOSUM62 should isolate the residues that substitute for nothing."""
        g = reduced_alphabet(8)
        assert g["I"] == g["L"] == g["V"], "the aliphatics must group"
        assert g["F"] == g["Y"], "the aromatics must group"
        assert g["K"] == g["R"], "the positives must group"
        assert len({g["C"], g["P"], g["W"]}) == 3, "C, P and W each substitute for little"
        assert g["C"] not in {g[a] for a in AA if a != "C"}

    def test_the_partition_is_complete_and_contiguous(self):
        for n in (2, 5, 8, 12, 19):
            g = reduced_alphabet(n)
            assert set(g) == set(AA)
            assert sorted(set(g.values())) == list(range(n))

    def test_twenty_groups_is_the_identity(self):
        g = reduced_alphabet(20)
        assert len(set(g.values())) == 20

    def test_it_is_deterministic(self):
        assert reduced_alphabet(8) == reduced_alphabet(8)

    def test_out_of_range_is_refused(self):
        with pytest.raises(ValueError):
            reduced_alphabet(1)
        with pytest.raises(ValueError):
            reduced_alphabet(21)

    def test_unmodelled_residues_void_a_window(self):
        table, n = alphabet_table(reduced_alphabet(8))
        assert n == 8
        assert table[ord("X")] == -1 and table[ord("*")] == -1 and table[ord("a")] == -1
        assert table[ord("A")] >= 0


class TestPatternSpec:
    def test_ungapped_and_gapped_parse(self):
        assert pattern_of("xxx") == [1, 1, 1]
        assert pattern_of("xx.x") == [1, 1, 0, 1]
        assert pattern_of("x.x.x") == [1, 0, 1, 0, 1]

    def test_a_leading_or_trailing_gap_is_refused(self):
        # It is the shorter pattern shifted, so admitting it puts two names on one column.
        for bad in (".xx", "xx.", ".x."):
            with pytest.raises(ValueError, match="start and end"):
                pattern_of(bad)

    def test_junk_is_refused(self):
        with pytest.raises(ValueError):
            pattern_of("xxAx")


class TestNativeKernelMatchesThePolarsPath:
    """The kernel exists to avoid materialising the explosion; it must not change the answer."""

    @pytest.mark.parametrize("spec,flank", [("xxx", 0), ("xxxx", 4), ("xx.x", 4), ("x.x.x", 4)])
    def test_codes_and_weights_agree_exactly(self, spec, flank):
        from vdjtools import _core

        df = sample(600, seed=11)
        # An X, and a junction shorter than its own flanks: the first must void only its window,
        # the second must contribute nothing at all rather than a germline tail.
        j = df["junction_aa"].to_list()
        j[0], j[1] = "CASXSLKF", "CF"
        df = df.with_columns(pl.Series("junction_aa", j))
        pat = pattern_of(spec)
        groups = reduced_alphabet(20)
        table, n_alpha = alphabet_table(groups)
        vlist = sorted(set(df["v_call"]))
        vidx = {v: i for i, v in enumerate(vlist)}
        vcodes = [vidx[v] for v in df["v_call"]]
        w = df["frequency"].to_numpy().astype(float)

        row = _core.kmer_row(df["junction_aa"].to_list(), vcodes, w.tolist(), pat,
                             table.tolist(), n_alpha, len(vlist), flank)
        got = dict(zip(row.codes, row.weights))

        span = n_alpha ** sum(pat)
        ref: dict[int, float] = {}
        ex = _explode_kmers(df.with_columns(pl.Series("_w", w)), len(pat), flank)
        for kmer, vc, ww in zip(ex["kmer"], ex["v_call"], ex["_w"]):
            sel = [kmer[j] for j, p in enumerate(pat) if p]
            if any(c not in AMINO_ACIDS for c in sel):
                continue
            code = vidx[vc] * span + sum(AMINO_ACIDS.index(c) * n_alpha ** j
                                         for j, c in enumerate(sel))
            ref[code] = ref.get(code, 0.0) + ww

        assert set(got) == set(ref)
        for c, v in ref.items():
            assert got[c] == pytest.approx(v)

    def test_the_code_space_refuses_to_overflow(self):
        from vdjtools import _core

        # Aliasing two k-mers onto one code would read downstream as a real shared feature.
        with pytest.raises(OverflowError):
            _core.kmer_code_space([1] * 20, 20, 30)


class TestFittedSpaceIsFrozen:
    def test_refit_reproduces_the_basis(self, cohort):
        a = fit_kmer_space(cohort, pattern="xxxx", n_groups=8, n_components=8, threads=4)
        b = fit_kmer_space(cohort, pattern="xxxx", n_groups=8, n_components=8, threads=1)
        assert np.array_equal(a.codes, b.codes)
        assert np.allclose(a.idf, b.idf)
        # Thread count must not reach the numbers, and ARPACK must not start from a random vector.
        assert np.allclose(a.components, b.components, atol=1e-10)
        assert np.allclose(a.transform(cohort[0]), b.transform(cohort[0]))

    def test_transform_width_follows_the_basis(self, cohort):
        proj = fit_kmer_space(cohort, n_groups=8, n_components=12, threads=2)
        assert proj.transform(cohort[0]).shape == (12,)
        raw = fit_kmer_space(cohort, n_groups=8, n_components=0, threads=2)
        assert raw.transform(cohort[0]).shape == (raw.n_columns,)

    def test_rows_are_depth_normalised(self, cohort):
        """A deeper sample must not simply have a longer vector."""
        sp = fit_kmer_space(cohort, n_groups=8, n_components=0, threads=2)
        x = sp.transform(cohort[0])
        assert np.linalg.norm(x) == pytest.approx(1.0)

    def test_an_unseen_v_gene_is_bucketed_not_dropped(self, cohort):
        """A collaborator on other nomenclature must not silently lose their k-mers."""
        sp = fit_kmer_space(cohort, n_groups=8, n_components=0, threads=2)
        odd = cohort[0].with_columns(pl.lit("TCRBV09-01").alias("v_call"))
        assert np.isfinite(sp.transform(odd)).all()

    def test_the_document_frequency_window_binds(self, cohort):
        wide = fit_kmer_space(cohort, n_groups=8, min_df=0.0, max_df=1.0, n_components=0,
                              threads=2)
        narrow = fit_kmer_space(cohort, n_groups=8, min_df=0.5, max_df=0.9, n_components=0,
                                threads=2)
        assert narrow.n_columns < wide.n_columns

    def test_the_column_cap_is_reported_when_it_binds(self, cohort):
        sp = fit_kmer_space(cohort, n_groups=8, max_columns=500, n_components=0, threads=2)
        assert sp.n_columns == 500
        assert sp.meta["dropped_by_cap"] > 0, "a silent truncation reads as full coverage"

    def test_an_empty_window_raises_rather_than_shipping_nothing(self, cohort):
        with pytest.raises(ValueError, match="no k-mer survives"):
            fit_kmer_space(cohort, n_groups=8, min_df=0.99, max_df=0.995, n_components=0,
                           threads=2)


class TestItFindsAPlantedMotif:
    @pytest.mark.parametrize("spec", ["xxxx", "xx.x", "x.x.x"])
    def test_a_single_component_separates_carriers(self, cohort, spec):
        y = np.array([1] * 20 + [0] * 20)
        sp = fit_kmer_space(cohort, pattern=spec, n_groups=8, min_df=0.05, max_df=0.95,
                            n_components=16, threads=4)
        Z = np.array([sp.transform(f) for f in cohort])
        d = np.abs(Z[y == 1].mean(0) - Z[y == 0].mean(0)) / (Z.std(0) + 1e-12)
        assert d.max() > 1.0, f"planted motif not recovered by any component (best {d.max():.2f})"

    def test_the_reduced_alphabet_makes_more_columns_estimable(self, cohort):
        """The point of the reduction: at 20 letters almost every 4-mer is too rare to keep."""
        full = fit_kmer_space(cohort, pattern="xxxx", n_groups=20, min_df=0.05, max_df=0.95,
                              n_components=0, threads=4)
        red = fit_kmer_space(cohort, pattern="xxxx", n_groups=8, min_df=0.05, max_df=0.95,
                             n_components=0, threads=4)
        assert red.meta["code_space"] < full.meta["code_space"]
        assert red.n_columns > full.n_columns


class TestTheBlockIsHeldOutUntilFitted:
    """Same rule as `pub`: no dead columns in the contract."""

    def test_the_layout_does_not_declare_kmer_before_a_space_exists(self):
        from vdjtools.signature import layout as L

        assert not any(b.name == "kmer" for b in L.registry()), \
            "a block whose artifact does not exist must not be in the contract"

    def test_a_fitted_set_of_spaces_produces_a_well_formed_block(self, cohort):
        from vdjtools.signature.kmer import kmer_block, kmer_spec

        sp = fit_kmer_space(cohort, n_groups=8, n_components=6, threads=2)
        spec = kmer_spec({"TRA": sp, "TRB": sp}, tier="full")
        assert spec.sig == "vsig" and spec.name == "kmer" and spec.attributable
        assert list(spec.features) == [f"PC0{i}" for i in range(1, 7)] + ["resid"]
        out = kmer_block(cohort[0], "TRB", sp)
        assert set(out) == set(spec.features)
        assert np.isfinite(list(out.values())).all()

    def test_ragged_ranks_are_refused(self, cohort):
        from vdjtools.signature.kmer import kmer_spec

        a = fit_kmer_space(cohort, n_groups=8, n_components=6, threads=2)
        b = fit_kmer_space(cohort, n_groups=8, n_components=4, threads=2)
        with pytest.raises(ValueError, match="same rank"):
            kmer_spec({"TRA": a, "TRB": b})

    def test_a_missing_space_is_a_hole_not_a_zero(self, cohort):
        from vdjtools.signature.kmer import kmer_block

        sp = fit_kmer_space(cohort, n_groups=8, n_components=6, threads=2)
        empty = cohort[0].clear()
        assert all(np.isnan(v) for v in kmer_block(empty, "TRB", sp).values())


class TestTheSpaceCanBeFrozen:
    def test_round_trip_reproduces_the_vector_exactly(self, cohort, tmp_path):
        from vdjtools.features.kmer_space import load_kmer_spaces, save_kmer_spaces

        a = fit_kmer_space(cohort, pattern="xx.x", n_groups=8, n_components=6, threads=2)
        b = fit_kmer_space(cohort, pattern="xxxx", n_groups=12, n_components=6, threads=2)
        p = tmp_path / "spaces.npz"
        save_kmer_spaces({"TRA": a, "TRB": b}, p)
        back = load_kmer_spaces(p)

        assert set(back) == {"TRA", "TRB"}
        for name, orig in (("TRA", a), ("TRB", b)):
            got = back[name]
            assert got.pattern == orig.pattern
            assert got.groups == orig.groups
            assert got.v_genes == orig.v_genes
            assert got.flank == orig.flank
            assert np.array_equal(got.codes, orig.codes)
            assert np.allclose(got.idf, orig.idf)
            # The point of freezing: the same repertoire must give the same numbers.
            assert np.allclose(got.transform(cohort[0]), orig.transform(cohort[0]))

    def test_the_alphabet_travels_as_a_partition_not_a_group_count(self, cohort, tmp_path):
        """A later change to the clustering must not re-partition an already-frozen space."""
        from vdjtools.features.kmer_space import save_kmer_spaces

        sp = fit_kmer_space(cohort, n_groups=8, n_components=4, threads=2)
        p = tmp_path / "s.npz"
        save_kmer_spaces({"TRB": sp}, p)
        z = np.load(p, allow_pickle=True)
        assert z["TRB/groups"].shape == (20,), "the partition itself must be stored"

    def test_an_unprojected_space_round_trips_too(self, cohort, tmp_path):
        from vdjtools.features.kmer_space import load_kmer_spaces, save_kmer_spaces

        sp = fit_kmer_space(cohort, n_groups=8, n_components=0, threads=2)
        p = tmp_path / "s.npz"
        save_kmer_spaces({"TRB": sp}, p)
        assert load_kmer_spaces(p)["TRB"].components is None


class TestTheAssemblerEmitsTheBlockOnceRegistered:
    def test_registered_columns_appear_and_are_finite(self, cohort):
        import importlib

        from vdjtools.signature import assemble, layout
        from vdjtools.signature.kmer import register_kmer

        # register_kmer mutates the process-wide registry, so reload it afterwards rather than
        # leaking a `kmer` block into every other test in the session.
        try:
            sp = fit_kmer_space(cohort, n_groups=8, n_components=4, threads=2)
            register_kmer({"TRB": sp}, tier="full")
            cols = layout.columns("full", "vsig")
            assert "vsig:kmer:TRB:PC01" in cols
            out = assemble.vsig({"TRB": cohort[0]}, tier="full", weight="log2p1",
                                kmer_spaces={"TRB": sp})
            vals = [out[f"vsig:kmer:TRB:PC{i:02d}"] for i in range(1, 5)]
            assert np.isfinite(vals).all(), vals
        finally:
            importlib.reload(layout)

    def test_without_registration_passing_spaces_is_a_no_op(self, cohort):
        from vdjtools.signature import assemble, layout

        sp = fit_kmer_space(cohort, n_groups=8, n_components=4, threads=2)
        assert not any(c.startswith("vsig:kmer") for c in layout.columns("full", "vsig"))
        out = assemble.vsig({"TRB": cohort[0]}, tier="full", kmer_spaces={"TRB": sp})
        assert not any(k.startswith("vsig:kmer") for k in out), "width must not change silently"


class TestAmbiguousVCallsResolveToOneGene:
    """A comma-ambiguous call must not become its own V bucket."""

    def test_resolve_gene_takes_the_first_after_stripping(self):
        from vdjtools.io.schema import resolve_gene, strip_allele

        s = pl.DataFrame({"v": ["TRBV9*01", "TRBV10-1,TRBV10-2", "TRBV1*01,TRBV23-1*02",
                                "TRBV1-2*02,TRBV1-2*04", None]})
        got = s.select(resolve_gene(pl.col("v")).alias("g"))["g"].to_list()
        assert got == ["TRBV9", "TRBV10-1", "TRBV1", "TRBV1-2", None]

        # strip_allele must still keep the tie -- the two serve different purposes
        kept = s.select(strip_allele(pl.col("v")).alias("g"))["g"].to_list()
        assert kept[2] == "TRBV1,TRBV23-1"

    def test_first_listed_means_first_as_written_not_alphabetically_first(self):
        """strip_allele sorts, so composing it with .list.first() picks the wrong gene."""
        from vdjtools.io.schema import resolve_gene

        s = pl.DataFrame({"v": ["TRBV5*01,TRBV19*03", "TRBV19*03,TRBV5*01"]})
        assert s.select(resolve_gene(pl.col("v")).alias("g"))["g"].to_list() == ["TRBV5", "TRBV19"]

    def test_the_fitted_v_vocabulary_has_no_comma_strings(self, cohort):
        amb = cohort[0].with_columns(
            pl.when(pl.int_range(pl.len()) % 3 == 0)
            .then(pl.col("v_call") + pl.lit(",TRBV30") )
            .otherwise(pl.col("v_call")).alias("v_call"))
        sp = fit_kmer_space([amb] + cohort[1:6], n_groups=8, n_components=0, threads=2)
        assert not any("," in v for v in sp.v_genes), \
            "each ambiguity string would otherwise be its own V bucket, draining the real one"

    def test_an_ambiguous_call_scores_as_its_first_gene(self, cohort):
        sp = fit_kmer_space(cohort, n_groups=8, n_components=0, threads=2)
        assert "TRBV5" in sp.v_genes
        clean = cohort[0].with_columns(pl.lit("TRBV5").alias("v_call"))
        amb = cohort[0].with_columns(pl.lit("TRBV5*01,TRBV19*03").alias("v_call"))
        assert np.allclose(sp.transform(clean), sp.transform(amb))


class TestTheResidualCarriesWhatTheBasisDiscards:
    def test_pythagoras_holds(self, cohort):
        """||x||^2 = ||Vx||^2 + resid^2 for an orthonormal basis of an L2-normalised vector."""
        sp = fit_kmer_space(cohort, n_groups=8, n_components=8, threads=2)
        z = sp.transform(cohort[0], residual=True)
        comp, resid = z[:-1], z[-1]
        assert np.isclose(np.dot(comp, comp) + resid ** 2, 1.0, atol=1e-6)

    def test_a_full_rank_basis_leaves_no_residual(self, cohort):
        sp = fit_kmer_space(cohort, n_groups=8, n_components=8, threads=2)
        sp.components = np.eye(sp.n_columns)
        assert sp.transform(cohort[0], residual=True)[-1] == pytest.approx(0.0, abs=1e-6)

    def test_the_residual_moves_when_a_rare_motif_is_planted(self, cohort):
        """The point of keeping it: a rare feature is orthogonal to the leading components."""
        sp = fit_kmer_space(cohort, n_groups=20, n_components=8, threads=2)
        plain = cohort[-1]
        spiked = plain.with_columns(
            pl.when(pl.int_range(pl.len()) < 40)
            .then(pl.lit("CWWWHWWWHWWWHWWF"))
            .otherwise(pl.col("junction_aa")).alias("junction_aa"))
        assert (sp.transform(spiked, residual=True)[-1]
                > sp.transform(plain, residual=True)[-1])

    def test_the_block_ships_it_alongside_the_components(self, cohort):
        from vdjtools.signature.kmer import kmer_block, kmer_spec

        sp = fit_kmer_space(cohort, n_groups=8, n_components=4, threads=2)
        assert list(kmer_spec({"TRB": sp}).features) == ["PC01", "PC02", "PC03", "PC04", "resid"]
        out = kmer_block(cohort[0], "TRB", sp)
        assert set(out) == {"PC01", "PC02", "PC03", "PC04", "resid"}
        assert np.isfinite(list(out.values())).all()
