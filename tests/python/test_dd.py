"""Tandem-D (D-D) extension: enumeration correctness + single-D backward compatibility.

The reference ``_dd_middle`` enumeration is validated on a hand-computable tiny model (fast and
exact); on the real TRD model we only assert that ``p_nd2 = 0`` reproduces the single-D Pgen
(the ``n_D=2`` enumeration is deliberately not exercised on full TRD — it is the native port's
job to make it fast).
"""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from vdjtools.model import Event, EventKind, Manifest, Model, from_olga, generate
from vdjtools.model import pgen
from vdjtools.model.dd import to_dd
from vdjtools.model.pgen import _dd_middle

OLGA = Path("/Users/mikesh/vcs/code/mirpy/mir/resources/olga/default_models")


def _tiny_dd_model(*, p_nd2: float = 0.5, pdel_full_only: bool = True, dd_ins1: float = 0.0) -> Model:
    """A minimal VDJ D-D model: 1 V, 1 J, 1 D (germline "AT"), zero-length insertions.

    With ``pdel_full_only`` the D can only contribute its whole "AT" (no partial trims), which
    makes the tandem enumeration exactly hand-checkable. ``dd_ins1`` sets ``P(dd insertion len=1)``
    so the DD-junction Markov path can be exercised.
    """
    def genes(seg, allele, cut):
        return pl.DataFrame({f"{seg}_allele": [allele], "cut_segment": [cut], "functional": [True]})

    ins0 = pl.DataFrame({"length": pl.Series([0, 1], dtype=pl.Int16), "p": [1.0, 0.0]})
    dd_ins = pl.DataFrame({"length": pl.Series([0, 1], dtype=pl.Int16), "p": [1.0 - dd_ins1, dd_ins1]})
    dinucl = pl.DataFrame({
        "from_nt": pl.Series([f for f in range(4) for _ in range(4)], dtype=pl.UInt8),
        "to_nt": pl.Series([t for _ in range(4) for t in range(4)], dtype=pl.UInt8),
        "p": [0.25] * 16,
    })
    ddel = pl.DataFrame({
        "d_allele": ["D", "D", "D"], "ndel5": pl.Series([0, 0, 1], dtype=pl.Int16),
        "ndel3": pl.Series([0, 1, 0], dtype=pl.Int16),
        "p": [1.0, 0.0, 0.0] if pdel_full_only else [0.5, 0.25, 0.25],
    })
    tables = {
        "v_choice": pl.DataFrame({"v_allele": ["V"], "p": [1.0]}),
        "j_choice": pl.DataFrame({"j_allele": ["J"], "p": [1.0]}),
        "d_gene": pl.DataFrame({"j_allele": ["J"], "d_allele": ["D"], "p": [1.0]}),
        "n_d": pl.DataFrame({"n_d": pl.Series([1, 2], dtype=pl.UInt8), "p": [1.0 - p_nd2, p_nd2]}),
        "v_3_del": pl.DataFrame({"v_allele": ["V"], "ndel": pl.Series([0], dtype=pl.Int16), "p": [1.0]}),
        "j_5_del": pl.DataFrame({"j_allele": ["J"], "ndel": pl.Series([0], dtype=pl.Int16), "p": [1.0]}),
        "d_del": ddel,
        "vd_ins": ins0, "dj_ins": ins0, "dd_ins": dd_ins,
        "vd_dinucl": dinucl, "dj_dinucl": dinucl, "dd_dinucl": dinucl,
        "d2_gene": pl.DataFrame({"d_allele": ["D"], "d2_allele": ["D"], "p": [1.0]}),
        "d2_del": ddel.rename({"d_allele": "d2_allele"}),
    }
    events = {
        "v_choice": Event("v_choice", EventKind.GENE_CHOICE),
        "j_choice": Event("j_choice", EventKind.GENE_CHOICE),
        "d_gene": Event("d_gene", EventKind.GENE_CHOICE, ("j_choice",)),
        "n_d": Event("n_d", EventKind.N_D),
        "v_3_del": Event("v_3_del", EventKind.DELETION, ("v_choice",)),
        "j_5_del": Event("j_5_del", EventKind.DELETION, ("j_choice",)),
        "d_del": Event("d_del", EventKind.DELETION_2D, ("d_gene",)),
        "vd_ins": Event("vd_ins", EventKind.INS_LENGTH),
        "dj_ins": Event("dj_ins", EventKind.INS_LENGTH),
        "vd_dinucl": Event("vd_dinucl", EventKind.DINUCLEOTIDE),
        "dj_dinucl": Event("dj_dinucl", EventKind.DINUCLEOTIDE),
        "d2_gene": Event("d2_gene", EventKind.GENE_CHOICE, ("d_gene",)),
        "d2_del": Event("d2_del", EventKind.DELETION_2D, ("d2_gene",)),
        "dd_ins": Event("dd_ins", EventKind.INS_LENGTH),
        "dd_dinucl": Event("dd_dinucl", EventKind.DINUCLEOTIDE),
    }
    manifest = Manifest(locus="TEST", organism="synthetic", chain_type="VDJ", events=events,
                        palindrome_max={"v_3": 0, "j_5": 0, "d_5": 0, "d_3": 0, "d2_5": 0, "d2_3": 0},
                        source="tiny")
    genomic = {"genes_v": genes("v", "V", "GG"), "genes_j": genes("j", "J", "CC"),
               "genes_d": genes("d", "D", "AT")}
    return Model(manifest=manifest, tables=tables, genomic=genomic)


def test_dd_middle_exact_tandem():
    """middle = "ATAT" = D1("AT") + D2("AT"), no insertions → P = P(D1)·P(D2|D1) = 1."""
    prep = pgen.prepare(_tiny_dd_model())
    assert _dd_middle(prep, "J", "ATAT") == pytest.approx(1.0)


def test_dd_middle_requires_two_ds():
    """A single D fits "AT" but a tandem needs two ≥1-nt Ds → the n_D=2 middle is 0."""
    prep = pgen.prepare(_tiny_dd_model())
    assert _dd_middle(prep, "J", "AT") == 0.0


def test_dd_middle_trims_each_d_at_5p_and_3p():
    """Each tandem D is trimmed independently at 5' and 3'.

    With partial deletions allowed — P(del)= "AT":0.5, 3'-trim→"A":0.25, 5'-trim→"T":0.25 — the
    middle "AT" now IS a tandem: D1 3'-trimmed to "A" + D2 5'-trimmed to "T", zero insertions, so
    P = P(del_{D1}=(0,1))·P(del_{D2}=(1,0)) = 0.25·0.25 = 0.0625. (0 under full-D-only, above.)
    """
    prep = pgen.prepare(_tiny_dd_model(pdel_full_only=False))
    assert _dd_middle(prep, "J", "AT") == pytest.approx(0.0625)


def test_dd_middle_with_dd_insertion():
    """"ATGAT" = D1("AT") + insDD("G") + D2("AT"); with P(insDD len1)=0 the tandem is impossible."""
    prep = pgen.prepare(_tiny_dd_model())  # ins0 gives P(len=1)=0
    assert _dd_middle(prep, "J", "ATGAT") == 0.0


def test_pgen_partitions_over_nd():
    """Full pgen_nt = P(n_D=1)·single-D + P(n_D=2)·tandem, on the tiny model."""
    m = _tiny_dd_model(p_nd2=0.5)
    prep = pgen.prepare(m)
    # sequence: V("G") + D1("AT") + D2("AT") + J("C")  → "GATATC"
    s = "GATATC"
    p_full = pgen.pgen_nt(prep, s)
    p_single = pgen.pgen_nt(pgen.prepare(_tiny_dd_model(p_nd2=0.0)), s)
    # tandem-only piece via a p2=1 model
    p_tandem = pgen.pgen_nt(pgen.prepare(_tiny_dd_model(p_nd2=1.0)), s)
    assert p_full == pytest.approx(0.5 * p_single + 0.5 * p_tandem)


@pytest.mark.slow
@pytest.mark.skipif(not OLGA.exists(), reason="OLGA models not available")
def test_dd_backward_compatible_on_trd():
    """to_dd(p_nd2=0) reproduces the single-D TRD Pgen exactly (no n_D=2 enumeration).

    Slow tier: the pure-Python reference ``pgen_nt`` on real TRD junctions is seconds/seq.
    """
    m = from_olga(OLGA / "human_T_delta", locus="TRD")
    prep0 = pgen.prepare(m)
    prep_dd0 = pgen.prepare(to_dd(m, p_nd2=0.0))
    seqs = [s.upper() for s in generate.generate(m, 20, seed=3)["cdr3_nt"].to_list()]
    checked = 0
    for s in seqs:
        p1 = pgen.pgen_nt(prep0, s)
        if p1 == 0.0:
            continue
        assert pgen.pgen_nt(prep_dd0, s) == pytest.approx(p1, rel=1e-12)
        checked += 1
    assert checked >= 3


def test_dd_middle_nonzero_dd_insertion():
    """Exercise the DD-junction Markov path: "ATGAT" = D1("AT") + insDD("G") + D2("AT").

    With P(insDD len=1)=0.5 and a uniform dinucleotide bias (0.25), the single tandem tiling gives
    P = pdel(AT)·pdel(AT)·P(insDD len1)·bias[G] = 1·1·0.5·0.25 = 0.125.
    """
    prep = pgen.prepare(_tiny_dd_model(dd_ins1=0.5))
    assert _dd_middle(prep, "J", "ATGAT") == pytest.approx(0.125)


def test_vdj_middle_is_nd_mixture():
    """_vdj_middle == P(n_D=1)·single-D + P(n_D=2)·tandem, directly (not just end-to-end)."""
    from vdjtools.model.pgen import _d_middle, _vdj_middle
    prep = pgen.prepare(_tiny_dd_model(p_nd2=0.3))
    mid = "GATATC"[1:-1]  # "ATAT": producible as a tandem (and not as single-D here)
    assert _vdj_middle(prep, "J", mid) == pytest.approx(
        0.7 * _d_middle(prep, "J", mid) + 0.3 * _dd_middle(prep, "J", mid)
    )


def test_prepare_rejects_mass_at_2_without_d2_tables():
    """A model with P(n_D=2)>0 but no d2_gene table is malformed → prepare raises (no silent drop)."""
    m = _tiny_dd_model(p_nd2=0.5)
    bad = Model(manifest=m.manifest,
                tables={k: v for k, v in m.tables.items() if k not in ("d2_gene", "d2_del")},
                genomic=m.genomic)
    with pytest.raises(ValueError, match="malformed tandem"):
        pgen.prepare(bad)


def test_to_dd_rejects_bad_inputs():
    m = _tiny_dd_model()  # already-DD
    with pytest.raises(ValueError, match="already"):
        to_dd(m)
    # a fresh single-D VDJ model for the p_nd2-bound checks
    if OLGA.exists():
        sd = from_olga(OLGA / "human_T_delta", locus="TRD")
        for bad in (1.0, -0.1):
            with pytest.raises(ValueError):
                to_dd(sd, p_nd2=bad)
        with pytest.raises(ValueError, match="VDJ"):
            to_dd(from_olga(OLGA / "human_T_alpha", locus="TRA"))  # VJ locus


@pytest.mark.skipif(not OLGA.exists(), reason="OLGA models not available")
def test_to_dd_tables_wellformed():
    m2 = to_dd(from_olga(OLGA / "human_T_delta", locus="TRD"), p_nd2=0.05)
    m2.validate()  # columns + normalization
    assert set(m2.tables["d2_gene"].columns) == {"d_allele", "d2_allele", "p"}
    assert set(m2.tables["d2_del"].columns) == {"d2_allele", "ndel5", "ndel3", "p"}
    nd = dict(zip(m2.tables["n_d"]["n_d"].to_list(), m2.tables["n_d"]["p"].to_list()))
    assert nd[2] == pytest.approx(0.05) and nd[1] == pytest.approx(0.95)


@pytest.mark.skipif(not OLGA.exists(), reason="OLGA models not available")
def test_dd_guards_raise_and_single_d_untouched():
    """Not-yet-tandem paths refuse a D-D model; single-D (incl. to_dd p_nd2=0) flows through."""
    from vdjtools.model import generate, infer, native
    m = from_olga(OLGA / "human_T_delta", locus="TRD")
    m_dd = to_dd(m, p_nd2=0.1)
    s = generate.generate(m, 1, seed=1)["cdr3_nt"][0].upper()
    # aa Pgen (transfer matrix) and EM do not yet support tandems → must raise.
    for fn in (
        lambda: pgen.pgen_aa(pgen.prepare(m_dd), "CAAAF"),
        lambda: native.pgen_aa(m_dd, "CAAAF"),
        lambda: infer.infer(m_dd, [s], max_iter=1),
        lambda: infer.infer_native(m_dd, [s], max_iter=1),
    ):
        with pytest.raises(NotImplementedError):
            fn()
    # native *nt* Pgen and generation DO support tandems now: must not raise.
    assert native.pgen_nt(m_dd, s) >= 0.0
    assert generate.generate(m_dd, 3, seed=1).height == 3
    # to_dd(p_nd2=0) is generatively single-D → native stays byte-identical to the single-D model.
    m_dd0 = to_dd(m, p_nd2=0.0)
    assert native.pgen_nt(m_dd0, s) == pytest.approx(native.pgen_nt(m, s), rel=1e-12)


def test_generate_dd_fraction_and_scoreable():
    """Tandem generation: the n_D=2 draw fraction matches P(n_D=2), and every tandem draw is
    scoreable by the native D-D Pgen (closed-loop consistency needed for D-D EM)."""
    from vdjtools.model import generate, native
    m = _tiny_dd_model(p_nd2=0.4, pdel_full_only=False)
    df = generate.generate(m, 4000, seed=1)
    assert df["d2_call"].is_not_null().mean() == pytest.approx(0.4, abs=0.03)  # ~ P(n_D=2)
    tand = df.filter(df["d2_call"].is_not_null()).head(20)
    assert all(native.pgen_nt(m, r["cdr3_nt"].upper(), r["v_call"], r["j_call"]) > 0
               for r in tand.to_dicts())


def test_native_dd_matches_python_reference():
    """Native tandem-D nt Pgen == the pure-Python reference, exactly (tiny hand-checkable model)."""
    from vdjtools.model import native
    for kw in (dict(p_nd2=0.3), dict(p_nd2=0.3, pdel_full_only=False), dict(p_nd2=0.5, dd_ins1=0.5)):
        m = _tiny_dd_model(**kw)
        prep = pgen.prepare(m)
        for s in ("GATATC", "GATC", "GATGATC", "GATGGATC"):
            assert native.pgen_nt(m, s) == pytest.approx(pgen.pgen_nt(prep, s), rel=1e-12, abs=1e-18)


@pytest.mark.skipif(not OLGA.exists(), reason="OLGA models not available")
def test_dd_diagnostics_render():
    """A D-D model exposes the n_d / d2_gene / dd_* nodes to the diagnostics."""
    from vdjtools.model import analyze
    m_dd = to_dd(from_olga(OLGA / "human_T_delta", locus="TRD"), p_nd2=0.05)
    ent = analyze.entropy_table(m_dd)
    assert set(ent["event"]) >= {"n_d", "d2_gene", "d2_del", "dd_ins", "dd_dinucl"}
    assert ent.filter(pl.col("event") == "n_d")["H_bits"][0] > 0.0  # non-degenerate D-count
    dot = analyze.bayes_net_dot(m_dd)
    assert '"d2_gene"' in dot and '"d_gene" -> "d2_gene"' in dot
