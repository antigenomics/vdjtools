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
    # V("GG") + D1("AT") + D2("AT") + J("CC") → "GGATATCC" (tandem-only; single-D can't span "ATAT").
    s = "GGATATCC"
    p_full = pgen.pgen_nt(pgen.prepare(_tiny_dd_model(p_nd2=0.5)), s)
    p_single = pgen.pgen_nt(pgen.prepare(_tiny_dd_model(p_nd2=0.0)), s)
    p_tandem = pgen.pgen_nt(pgen.prepare(_tiny_dd_model(p_nd2=1.0)), s)
    assert p_tandem > 0.0 and p_single == 0.0  # non-vacuous: this read is genuinely tandem-only
    assert p_full == pytest.approx(0.5 * p_single + 0.5 * p_tandem)


def test_dd_estep_total_equals_pgen():
    """The D-D E-step's soft-count normalizer equals the independently-computed Pgen (all scenarios)."""
    from collections import defaultdict
    from vdjtools.model.infer import _estep_seq, _fit_events
    m = _tiny_dd_model(p_nd2=0.5, pdel_full_only=False)
    prep = pgen.prepare(m)
    for s in ("GGATATCC", "GGATCC", "GGATGCC"):
        counts = {n: defaultdict(float) for n in _fit_events(m.manifest)}
        assert _estep_seq(prep, s, counts) == pytest.approx(pgen.pgen_nt(prep, s), rel=1e-12, abs=1e-18)


def test_em_recovers_p_nd2():
    """Closed-loop: generate tandem data from a known model → EM recovers P(n_D=2) from a wrong init."""
    from vdjtools.model import generate
    from vdjtools.model.infer import infer
    reads = [r.upper() for r in generate.generate(_tiny_dd_model(p_nd2=0.35, pdel_full_only=False),
                                                  1500, seed=11)["cdr3_nt"].to_list()]
    fitted, _ = infer(_tiny_dd_model(p_nd2=0.15, pdel_full_only=False), reads, max_iter=25, init="template")
    nd = dict(zip(fitted.tables["n_d"]["n_d"].to_list(), fitted.tables["n_d"]["p"].to_list()))
    assert nd[2] == pytest.approx(0.35, abs=0.05)  # recovered from a deliberately wrong 0.15 start


def test_native_dd_estep_matches_python_reference():
    """The native C++ D-D E-step accumulates the *same* soft counts as the pure-Python ``_accum_dd``
    reference — exactly (factorized forward/backward attribution vs the naive two-D enumeration)."""
    from collections import defaultdict

    from vdjtools.model import native
    from vdjtools.model.infer import _estep_seq, _fit_events
    _core = pytest.importorskip("vdjtools._core")

    m = _tiny_dd_model(p_nd2=0.5, pdel_full_only=False, dd_ins1=0.3)
    seqs = ["GGATATCC", "GGATCC", "GGATAATCC", "GGATATATCC"]
    prep = pgen.prepare(m)
    pyc = {n: defaultdict(float) for n in _fit_events(m.manifest)}
    for s in seqs:  # _estep_seq already stores the posterior-weighted counts (w / Pgen)
        _estep_seq(prep, s, pyc, None)

    pm, _vi, _ji = native.pack(m)
    cn = _core.make_counts(pm)
    _core.estep_batch(pm, [native._encode(s) for s in seqs], [], [], [], cn)

    d = m.genomic["genes_d"]["d_allele"].to_list()
    nD, n5, n3 = len(d), pm.nbins_d5, pm.nbins_d3
    m5, m3 = prep.maxpal["d_5"], prep.maxpal["d_3"]

    def diff(name, arr, key_of):
        nd_ = {key_of(i): arr[i] for i in range(len(arr))}
        return max((abs(nd_.get(k, 0.0) - pyc[name].get(k, 0.0)) for k in set(nd_) | set(pyc[name])), default=0.0)

    checks = {
        "v_choice": (cn.v_choice, lambda i: ("V",)),
        "j_choice": (cn.j_choice, lambda i: ("J",)),
        "d_gene": (cn.d_gene, lambda i: ("J", d[i % nD])),
        "d2_gene": (cn.d2_gene, lambda i: (d[i // nD], d[i % nD])),
        "d_del": (cn.d_del, lambda i: (d[i // (n5 * n3)], (i // n3) % n5 - m5, i % n3 - m3)),
        "d2_del": (cn.d2_del, lambda i: (d[i // (n5 * n3)], (i // n3) % n5 - m5, i % n3 - m3)),
        "vd_ins": (cn.ins_vd, lambda i: (i,)),
        "dd_ins": (cn.ins_dd, lambda i: (i,)),
        "dj_ins": (cn.ins_dj, lambda i: (i,)),
        "n_d": (cn.n_d, lambda i: (i,)),
        "vd_dinucl": (cn.dinucl_vd, lambda i: (i % 4, i // 4)),
        "dd_dinucl": (cn.dinucl_dd, lambda i: (i % 4, i // 4)),
        "dj_dinucl": (cn.dinucl_dj, lambda i: (i % 4, i // 4)),
    }
    for name, (arr, key_of) in checks.items():
        assert diff(name, arr, key_of) < 1e-12, f"{name} native vs python soft counts differ"


def test_native_dd_em_recovers_p_nd2():
    """Closed-loop with the native E-step: native EM recovers P(n_D=2) from a wrong init (fast path)."""
    from vdjtools.model import generate
    from vdjtools.model.infer import infer_native
    pytest.importorskip("vdjtools._core")
    reads = [r.upper() for r in generate.generate(_tiny_dd_model(p_nd2=0.35, pdel_full_only=False),
                                                  1500, seed=11)["cdr3_nt"].to_list()]
    fitted, _ = infer_native(_tiny_dd_model(p_nd2=0.15, pdel_full_only=False), reads, max_iter=25, init="template")
    nd = dict(zip(fitted.tables["n_d"]["n_d"].to_list(), fitted.tables["n_d"]["p"].to_list()))
    assert nd[2] == pytest.approx(0.35, abs=0.05)


def test_dd_anchor_and_prior_regularize_and_match():
    """The per-read D-D gate (``dd_allowed``) and the ``nd_prior`` Dirichlet pseudocount both pull the
    learned ``P(n_D=2)`` below the unregularized value, and native == Python for each."""
    from vdjtools.model import generate
    from vdjtools.model.infer import infer, infer_native
    pytest.importorskip("vdjtools._core")
    reads = [r.upper() for r in generate.generate(_tiny_dd_model(p_nd2=0.35, pdel_full_only=False),
                                                  1200, seed=11)["cdr3_nt"].to_list()]
    gate = [i % 2 == 0 for i in range(len(reads))]  # allow D-D on half the reads

    def p2(model_and_report):
        m = model_and_report[0]
        return dict(zip(m.tables["n_d"]["n_d"].to_list(), m.tables["n_d"]["p"].to_list()))[2]

    def seed():
        return _tiny_dd_model(p_nd2=0.15, pdel_full_only=False)

    base = p2(infer_native(seed(), reads, max_iter=20, init="template"))
    for kw in ({"dd_allowed": gate}, {"nd_prior": 200.0}):
        py = p2(infer(seed(), reads, max_iter=20, init="template", **kw))
        nat = p2(infer_native(seed(), reads, max_iter=20, init="template", **kw))
        assert nat == pytest.approx(py, abs=1e-4)  # native E-step == Python reference
        assert nat < base - 0.02                    # regularizer lowered the tandem estimate


def _tiny_aligned_dd_model(p_nd2=0.6) -> Model:
    """Codon-aligned tiny D-D model (V=TGT=Cys, D=CA, J=TTT=Phe) with nonzero VD/DD/DJ insertions
    (len 0/1/2) and a non-uniform DD dinucleotide, so an aa CDR3 can require a real insDD block."""
    def genes(seg, a, cut):
        return pl.DataFrame({f"{seg}_allele": [a], "cut_segment": [cut], "functional": [True]})
    ins = pl.DataFrame({"length": pl.Series([0, 1, 2], dtype=pl.Int16), "p": [0.2, 0.3, 0.5]})
    unif = pl.DataFrame({"from_nt": pl.Series([f for f in range(4) for _ in range(4)], dtype=pl.UInt8),
                         "to_nt": pl.Series([t for _ in range(4) for t in range(4)], dtype=pl.UInt8), "p": [0.25] * 16})
    ddn = pl.DataFrame({"from_nt": pl.Series([f for f in range(4) for _ in range(4)], dtype=pl.UInt8),
                        "to_nt": pl.Series([t for _ in range(4) for t in range(4)], dtype=pl.UInt8), "p": [0.4, 0.3, 0.2, 0.1] * 4})
    ddel = pl.DataFrame({"d_allele": ["D"], "ndel5": pl.Series([0], dtype=pl.Int16),
                         "ndel3": pl.Series([0], dtype=pl.Int16), "p": [1.0]})
    tables = {
        "v_choice": pl.DataFrame({"v_allele": ["V"], "p": [1.0]}),
        "j_choice": pl.DataFrame({"j_allele": ["J"], "p": [1.0]}),
        "d_gene": pl.DataFrame({"j_allele": ["J"], "d_allele": ["D"], "p": [1.0]}),
        "n_d": pl.DataFrame({"n_d": pl.Series([1, 2], dtype=pl.UInt8), "p": [1 - p_nd2, p_nd2]}),
        "v_3_del": pl.DataFrame({"v_allele": ["V"], "ndel": pl.Series([0], dtype=pl.Int16), "p": [1.0]}),
        "j_5_del": pl.DataFrame({"j_allele": ["J"], "ndel": pl.Series([0], dtype=pl.Int16), "p": [1.0]}),
        "d_del": ddel, "vd_ins": ins, "dj_ins": ins, "dd_ins": ins,
        "vd_dinucl": unif, "dj_dinucl": unif, "dd_dinucl": ddn,
        "d2_gene": pl.DataFrame({"d_allele": ["D"], "d2_allele": ["D"], "p": [1.0]}),
        "d2_del": ddel.rename({"d_allele": "d2_allele"})}
    ev = {"v_choice": Event("v_choice", EventKind.GENE_CHOICE), "j_choice": Event("j_choice", EventKind.GENE_CHOICE),
          "d_gene": Event("d_gene", EventKind.GENE_CHOICE, ("j_choice",)), "n_d": Event("n_d", EventKind.N_D),
          "v_3_del": Event("v_3_del", EventKind.DELETION, ("v_choice",)),
          "j_5_del": Event("j_5_del", EventKind.DELETION, ("j_choice",)),
          "d_del": Event("d_del", EventKind.DELETION_2D, ("d_gene",)), "vd_ins": Event("vd_ins", EventKind.INS_LENGTH),
          "dj_ins": Event("dj_ins", EventKind.INS_LENGTH), "vd_dinucl": Event("vd_dinucl", EventKind.DINUCLEOTIDE),
          "dj_dinucl": Event("dj_dinucl", EventKind.DINUCLEOTIDE),
          "d2_gene": Event("d2_gene", EventKind.GENE_CHOICE, ("d_gene",)),
          "d2_del": Event("d2_del", EventKind.DELETION_2D, ("d2_gene",)),
          "dd_ins": Event("dd_ins", EventKind.INS_LENGTH), "dd_dinucl": Event("dd_dinucl", EventKind.DINUCLEOTIDE)}
    man = Manifest(locus="TEST", organism="synthetic", chain_type="VDJ", events=ev,
                   palindrome_max={"v_3": 0, "j_5": 0, "d_5": 0, "d_3": 0, "d2_5": 0, "d2_3": 0}, source="tiny")
    gen = {"genes_v": genes("v", "V", "TGT"), "genes_j": genes("j", "J", "TTT"), "genes_d": genes("d", "D", "CA")}
    return Model(manifest=man, tables=tables, genomic=gen)


def test_aa_dd_equals_nt_sum():
    """aa D-D Pgen (Python and native) == Σ nt-Pgen over synonymous codons, on a codon-aligned model.

    Also pins native nt D-D Pgen: ``native.pgen_nt`` routes an in-frame CDR3 through the aa transfer
    matrix (singleton codon masks, incl. the ``p_nd2`` tandem term) and must equal the Python
    ``_dd_middle`` enumeration per nt. ``CHHF`` is tandem-only (single D can't tile it); ``CHIF``
    requires a real DD insertion.
    """
    import itertools
    from collections import defaultdict

    from vdjtools.model import native
    from vdjtools.model.reference import _CODON_TABLE
    pytest.importorskip("vdjtools._core")
    m = _tiny_aligned_dd_model(p_nd2=0.6)
    prep = pgen.prepare(m)
    syn = defaultdict(list)
    for cod, a in _CODON_TABLE.items():
        syn[a].append(cod)
    tandem_seen = False
    for aa in ("CHF", "CHHF", "CHIF"):
        brute = 0.0
        for c in itertools.product(*[syn[a] for a in aa]):
            nt = "".join(c)
            p_ref = pgen.pgen_nt(prep, nt, "V", "J")  # Python enumeration (single-D + _dd_middle)
            assert native.pgen_nt(m, nt, "V", "J") == pytest.approx(p_ref, rel=1e-9, abs=1e-300)
            brute += p_ref
        assert pgen.pgen_aa(prep, aa, "V", "J") == pytest.approx(brute, rel=1e-9, abs=1e-300)
        assert native.pgen_aa(m, aa, "V", "J") == pytest.approx(brute, rel=1e-8, abs=1e-300)
        if aa == "CHHF":
            tandem_seen = brute > 0.0
    assert tandem_seen  # the tandem-only CDR3 must have nonzero Pgen


def test_native_aa_dd_hamming1_and_agnostic():
    """Native aa D-D: the Hamming-1 ball equals the brute-force neighbour sum, and the v/j-agnostic
    query equals the Python reference."""
    from vdjtools.model import native
    pytest.importorskip("vdjtools._core")
    aas = "ACDEFGHIKLMNPQRSTVWY"
    m = _tiny_aligned_dd_model(p_nd2=0.6)
    prep = pgen.prepare(m)
    for aa in ("CHF", "CHHF", "CHIF"):
        assert native.pgen_aa(m, aa) == pytest.approx(pgen.pgen_aa(prep, aa), rel=1e-8, abs=1e-300)  # v/j-agnostic
        ball = native.pgen_aa(m, aa, "V", "J") + sum(
            native.pgen_aa(m, aa[:k] + x + aa[k + 1:], "V", "J")
            for k in range(len(aa)) for x in aas if x != aa[k])
        assert native.pgen_aa(m, aa, "V", "J", mismatches=1) == pytest.approx(ball, rel=1e-9, abs=1e-300)


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
def test_dd_full_native_support_and_single_d_untouched():
    """The whole pipeline supports tandem-D natively (nothing raises); single-D stays byte-identical."""
    from vdjtools.model import generate, infer, native
    m = from_olga(OLGA / "human_T_delta", locus="TRD")
    m_dd = to_dd(m, p_nd2=0.1)
    s = generate.generate(m, 1, seed=1)["cdr3_nt"][0].upper()
    # nt Pgen, aa Pgen, Hamming-1, v/j-agnostic, generation, native EM — all native, none raise.
    # (aa D-D numerical correctness vs the Python reference / Σnt is test_aa_dd_equals_nt_sum; the
    # pure-Python aa D-D is intractable on real-length TRD, so it is not exercised here.)
    assert native.pgen_aa(m_dd, "CAAAF") >= 0.0                     # native aa D-D transfer matrix
    assert native.pgen_aa(m_dd, "CAAAF", mismatches=1) >= 0.0       # native Hamming-1 D-D
    assert native.pgen_aa(m_dd, "CAAAF", v=None, j=None) >= 0.0     # native v/j-agnostic D-D
    assert native.pgen_nt(m_dd, s) >= 0.0
    assert generate.generate(m_dd, 3, seed=1).height == 3
    assert infer.infer_native(m_dd, [s], max_iter=1)[0].tables["n_d"].height == 2
    # to_dd(p_nd2=0) is generatively single-D → native stays byte-identical to the single-D model.
    m_dd0 = to_dd(m, p_nd2=0.0)
    assert native.pgen_nt(m_dd0, s) == pytest.approx(native.pgen_nt(m, s), rel=1e-12)
    assert native.pgen_aa(m_dd0, "CAAAF") == pytest.approx(native.pgen_aa(m, "CAAAF"), rel=1e-12)


@pytest.mark.skipif(not OLGA.exists(), reason="OLGA models not available")
def test_dd_default_for_d_loci():
    """EM defaults to a tandem-D model on the D-bearing loci (IGH/TRD/TRB); single_d opts out; VJ
    loci are unaffected."""
    from vdjtools.model import generate
    from vdjtools.model.infer import DD_DEFAULT_LOCI, gene_masks, infer_native
    pytest.importorskip("vdjtools._core")
    assert DD_DEFAULT_LOCI == {"TRB", "TRD", "IGH"}
    m = from_olga(OLGA / "human_T_delta", locus="TRD")
    df = generate.generate(m, 30, seed=3)
    seqs = [s.upper() for s in df["cdr3_nt"]]
    mk = gene_masks(m, df["v_call"].to_list(), df["j_call"].to_list())
    dd, _ = infer_native(m, seqs, masks=mk, max_iter=1)                 # default → D-D
    sd, _ = infer_native(m, seqs, masks=mk, max_iter=1, single_d=True)  # opt out → single-D
    assert dd.tables["n_d"].height == 2 and "d2_gene" in dd.tables
    assert sd.tables["n_d"].height == 1 and "d2_gene" not in sd.tables


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
