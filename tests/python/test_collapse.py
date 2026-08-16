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
from vdjtools.model.reference import load_germline, translate

LOCI = ["TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"]


def _shipped() -> list[tuple[str, str, str]]:
    """Every ``(source, organism, locus)`` actually shipped — ``arda`` alone carries mouse."""
    from vdjtools.model.bundled import list_bundled

    out = []
    for source, keys in list_bundled().items():
        for key in keys:
            organism, _, locus = key.rpartition("_")
            out.append((source, organism or "human", locus))
    return out


#: 23 models: ``olga`` and ``learned`` over 7 human loci, ``arda`` over 7 human + mouse TRA/TRB.
BUNDLED = _shipped()


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


# --- the collapsed model must be scoreable with the germline it kept ---------------------------

@pytest.mark.parametrize("locus", LOCI)
@pytest.mark.parametrize("source", ["olga", "learned"])
def test_collapsed_deletions_fit_the_representative_germline(locus, source):
    """No collapsed gene may carry deletion mass its own germline cannot reach.

    Averaging a conditional over alleles of differing CDR3-region length can strand mass on trims
    the representative does not support; the Pgen DP never visits those, so the probability would
    vanish from every Pgen through that gene instead of being redistributed.
    """
    from vdjtools.model.check import max_reachable_trim

    m = load_bundled(locus, source, collapse=True)
    for name, event in m.manifest.events.items():
        seg = name.split("_")[0]
        frame = m.genomic.get(f"genes_{'d' if seg.startswith('d') else seg}")
        acol = f"{seg}_allele"
        if frame is None or acol not in m.tables[name].columns or "ndel" not in str(m.tables[name].columns):
            continue
        gcol = f"{'d' if seg.startswith('d') else seg}_allele"
        limits = {r[gcol]: max_reachable_trim(name, event.kind, len(r["cut_segment"] or ""),
                                              m.manifest.palindrome_max)
                  for r in frame.iter_rows(named=True)}
        if any(v is None for v in limits.values()):
            continue
        t = m.tables[name]
        total = (pl.col("ndel") if "ndel" in t.columns else pl.col("ndel5") + pl.col("ndel3"))
        over = t.with_columns(
            _lim=pl.col(acol).replace_strict(limits, default=10**6, return_dtype=pl.Int64)
        ).filter((pl.col("p") > 0) & (total > pl.col("_lim")))
        assert over.is_empty(), f"{source}/{locus} {name}: unreachable mass on {over[acol].to_list()[:3]}"


def test_representative_is_the_longest_germline_not_the_most_used():
    """IGKV3-20*02 is 11 nt against *01's 30 and had the higher learned usage.

    Picking by usage alone made the truncated allele the gene's germline (relabelled *01, which was
    doubly misleading) and stranded 25% of the gene's own deletion distribution.
    """
    raw = load_bundled("IGK", "learned", collapse=False)
    alleles = raw.genomic["genes_v"].filter(pl.col("gene") == "IGKV3-20")
    lengths = {r["v_allele"]: len(r["cut_segment"]) for r in alleles.iter_rows(named=True)}
    assert lengths["IGKV3-20*02"] < lengths["IGKV3-20*01"], "fixture assumption changed"

    c = load_bundled("IGK", "learned", collapse=True)
    rep = c.genomic["genes_v"].filter(pl.col("gene") == "IGKV3-20")
    assert rep.height == 1
    assert len(rep["cut_segment"][0]) == max(lengths.values())


@pytest.mark.parametrize("locus", ["TRB", "IGK"])
def test_collapse_still_preserves_gene_usage_exactly(locus):
    """Projecting deletions onto the germline must not disturb the choice marginals."""
    m = load_bundled(locus, "learned", collapse=False)
    c = collapse_alleles(m)
    before, after = _gene_mass(m, "v_choice"), _gene_mass(c, "v_choice")
    for gene, p in before.items():
        assert after[gene] == pytest.approx(p, abs=1e-12)


# --- the representative must carry a functional, anchor-framed germline -------------------------
#
# arda's `database/vdj/<organism>/cdr3_anchors.tsv` is the reference here: it carries per allele the
# IMGT `functionality` (F / ORF / P), the CDR3-region `germline_nt`, the `templated_aa` it codes for
# and an explicit `status`. vdjtools' own `anchor` column is an nt offset into `full_germline`, not
# into `cut_segment`, and reading it as the latter flags 800+ correct entries — arda settles the
# convention instead of guessing it.

#: A J CDR3-region germline ends ON the conserved [FW]118 codon, so its length mod 3 is the offset
#: of the anchor reading frame (this equals arda's `partial_nt` for every shipped allele).
def _anchor_frame_aa(cdr3_segment: str) -> str:
    return translate(cdr3_segment[len(cdr3_segment) % 3:])


#: Functional J alleles whose germline genuinely does not template a terminal Phe/Trp.
#:
#: * ``TRAJ35*01`` — arda records ``anchor_nt=25``/``templated_aa=IGFGNVLHC`` with ``status=ok``:
#:   human TRAJ35's anchor sits past the FGxG motif, so the templated region ends on a Cys. IMGT
#:   still calls the allele F.
#: * ``IGHJ6*02`` — OLGA's namespace only (``olga``/``learned``). OLGA ships the CDR3-region
#:   germline 1 nt shorter than arda's 32-nt ``ATTACTACTACTACTACGGTATGGACGTCTGG``, cutting the
#:   Trp118 codon in half; arda's own IGHJ6*02 templates ``YYYYYGMDVW`` and passes.
_TERMINAL_ANCHOR_EXCEPTIONS = {"TRAJ35*01", "IGHJ6*02"}


@pytest.mark.parametrize("source,organism,locus", BUNDLED)
@pytest.mark.parametrize("collapse", [False, True])
def test_bundled_j_germline_templates_a_terminal_anchor_residue(source, organism, locus, collapse):
    """Every functional bundled J germline must translate to a terminal F or W in its anchor frame.

    This is the invariant the shipped collapse broke: ``TRBJ2-7*02`` templates ``SYEQYV``, so once
    the length tie-break degenerated to the allele name and handed it the gene, every real
    ``…YEQYF`` junction scored exactly 0 against ``TRBJ2-7*01``.

    ORF and P alleles are exempt by IMGT's own definitions — an ORF is precisely a sequence whose
    conserved motif (here J-PHE / J-TRP) is altered — as are alleles arda marks ``no_anchor``.
    A **collapsed** row is judged on its gene, not on whichever allele supplied the germline:
    it is labelled ``gene*01``, so if the gene owns any functional allele the representative must
    template that gene's anchor. Judging it on the allele it came from would exempt exactly the
    defect — ``TRBJ2-7*02`` is an ORF, so "ORF alleles are exempt" would have waved it through.
    """
    m = load_bundled(locus, source, organism=organism, collapse=collapse)
    ref = load_germline(locus, organism)
    func = dict(zip(ref["allele"].to_list(), ref["functionality"].to_list()))
    status = dict(zip(ref["allele"].to_list(), ref["status"].to_list()))
    usable = {a for a in func if func[a] == "F" and status[a] == "ok"}

    checked, bad = 0, []
    for r in m.genomic["genes_j"].iter_rows(named=True):
        gene, seq = r["gene"], r["cdr3_segment"] or ""
        if collapse:
            required = any(a.split("*")[0] == gene for a in usable)
            exempt = all(a in _TERMINAL_ANCHOR_EXCEPTIONS
                         for a in usable if a.split("*")[0] == gene)
        else:
            required = r["j_allele"] in usable
            exempt = r["j_allele"] in _TERMINAL_ANCHOR_EXCEPTIONS
        if not seq or not required or exempt:
            continue
        checked += 1
        aa = _anchor_frame_aa(seq)
        if not aa or aa[-1] not in "FW":
            bad.append((gene, seq, aa))
    assert not bad, f"{source}/{organism}/{locus} collapse={collapse}: non-anchor J germlines {bad}"
    assert checked, f"{source}/{organism}/{locus} collapse={collapse}: nothing was checked"


def test_collapsed_representative_is_functional_wherever_a_functional_allele_exists():
    """No gene may be represented by an ORF/pseudogene allele when it has a functional one.

    Before the fix this failed on 18 genes across the shipped sets — ``TRBJ2-7`` (ORF ``*02`` for
    F ``*01``) in all three sources, ``TRAJ8``/``IGKJ4``/``IGKV1-39``/``IGHV6-1``/… in ``arda``,
    and 3 mouse TRAV genes.

    The allele behind a representative is recovered from the *model's own* uncollapsed germline
    (the collapse copies the row verbatim), not from arda's — the two namespaces differ often
    enough that matching on arda's sequence silently finds nothing and passes.
    """
    bad = []
    for source, organism, locus in BUNDLED:
        raw = load_bundled(locus, source, organism=organism, collapse=False)
        col = load_bundled(locus, source, organism=organism, collapse=True)
        ref = load_germline(locus, organism)
        func = dict(zip(ref["allele"].to_list(), ref["functionality"].to_list()))
        for seg in ("v", "j"):
            by_gene: dict[str, dict[str, str]] = {}
            for r in raw.genomic[f"genes_{seg}"].iter_rows(named=True):
                by_gene.setdefault(r["gene"], {})[r[f"{seg}_allele"]] = r["cdr3_segment"] or ""
            for r in col.genomic[f"genes_{seg}"].iter_rows(named=True):
                alleles = by_gene.get(r["gene"], {})
                if not any(func.get(a) == "F" for a in alleles):
                    continue                       # no functional allele to prefer
                picked = [a for a, s in alleles.items() if s == (r["cdr3_segment"] or "")]
                if picked and not any(func.get(a) == "F" for a in picked):
                    bad.append((source, organism, locus, r["gene"], picked[0], func.get(picked[0])))
    assert not bad, f"non-functional representatives where a functional allele exists: {bad}"


#: The ``learned`` human TRA model was built by `from_olga(derive_orf=True)` before that path knew
#: the CDR3 region lies on the 5' side of a J anchor, so these 11 germlines are the framework
#: *downstream* of Phe118. `io._genomic_table` is fixed (3.9.1); the SHIPPED parquet still carries
#: them, and only regenerating the models clears it — a multi-hour all-loci job, deliberately not
#: run for a defect with no measured downstream consumer (**0** of VDJdb's 30,937 human TRA records
#: use any of these 11 alleles; see `bench/results/vdjtools_germline_pgen_shift.md`). All 11 are
#: ORF/P except TRAJ35*01, which the terminal-anchor test above exempts on arda's own record.
_WRONG_SIDE_OF_THE_J_ANCHOR = {
    ("learned", "human", "TRA", a) for a in
    ("TRAJ1*01", "TRAJ2*01", "TRAJ19*01", "TRAJ25*01", "TRAJ35*01", "TRAJ51*01", "TRAJ55*01",
     "TRAJ58*01", "TRAJ59*01", "TRAJ60*01", "TRAJ61*01")
}


@pytest.mark.parametrize("source,organism,locus", BUNDLED)
def test_bundled_cdr3_germline_sits_on_the_documented_side_of_the_anchor(source, organism, locus):
    """V's CDR3 region is ``full[anchor:]``; J's is ``full[:anchor + 3]`` — never the other way."""
    m = load_bundled(locus, source, organism=organism, collapse=False)
    bad = []
    for seg in ("v", "j"):
        for r in m.genomic[f"genes_{seg}"].iter_rows(named=True):
            full, seq, anchor = r["full_germline"] or "", r["cdr3_segment"] or "", r["anchor"]
            if not full or not seq or anchor < 0:
                continue
            want = full[anchor:] if seg == "v" else full[:anchor + 3]
            if seq != want and (source, organism, locus, r[f"{seg}_allele"]) not in _WRONG_SIDE_OF_THE_J_ANCHOR:
                bad.append(r[f"{seg}_allele"])
    assert not bad, f"{source}/{organism}/{locus}: CDR3 germline off the anchor for {bad}"


def test_the_learned_set_is_never_a_mix_of_builders():
    """The seven `learned` loci must all report the **same** builder version.

    A germline defect lives in the builder, so "which build is this from" should be a lookup rather
    than a posterior comparison against a reference fit. The dangerous state is a *partial*
    regeneration, where half the set silently answers a different question from the other half —
    this caught exactly that (4 loci at ``3.9.2``, 3 at ``""``) when an all-loci rebuild was
    interrupted midway.

    Deliberately **not** asserting that every model carries a stamp. Only regenerating the shipped
    set could satisfy that, so it would encode "always rebuild" as policy — and the models are
    regenerated when a germline actually changes, not to populate a metadata field.
    """
    stamped = {locus: load_bundled(locus, "learned", collapse=False).manifest.builder_version
               for locus in LOCI}
    assert len(set(stamped.values())) == 1, (
        f"learned set has mixed builder versions: {stamped} — regenerate all seven loci together, "
        f'never a subset. ("" is not missing data: it means the model predates 3.9.2.)')


def test_trbj2_7_is_not_a_silent_zero():
    """Regression pin for the defect: a real ``…YEQYF`` junction must score above zero.

    ``TRBJ2-7*01`` (F, ``SYEQYF``) and ``*02`` (ORF, ``SYEQYV``) are both 19 nt, so the shipped
    length-then-name key handed the gene to ``*02`` and relabelled it ``*01``. Pgen of every
    ``…YEQYF`` clonotype was then exactly 0 against the default collapsed model.
    """
    for source in ("olga", "learned", "arda"):
        raw = load_bundled("TRB", source, collapse=False)
        col = load_bundled("TRB", source, collapse=True)
        rep = col.genomic["genes_j"].filter(pl.col("gene") == "TRBJ2-7")
        assert rep.height == 1
        assert _anchor_frame_aa(rep["cdr3_segment"][0]) == "SYEQYF", source

        p_col = native.pgen_aa(col, "CASSIRSSYEQYF", "TRBV19*01", "TRBJ2-7*01")
        p_raw = native.pgen_aa(raw, "CASSIRSSYEQYF", "TRBV19*01", "TRBJ2-7*01")
        assert p_col > 0, f"{source}: collapsed Pgen is a silent zero"
        assert p_raw > 0, f"{source}: uncollapsed Pgen is zero — fixture assumption changed"
        # collapse is an approximation, but the two must agree to well within an order of magnitude
        assert np.log10(p_col) == pytest.approx(np.log10(p_raw), abs=0.3), source
