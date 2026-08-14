"""Named feature presets — which columns to use, and which not to.

A signature is 1,403 columns. Almost nobody wants all of them, and the right subset depends on the
question: a model that must transfer to another lab's cohort wants different columns from one
scoring samples inside a single study. This module names those choices, documents each one, and
**ranks it**, so a collaborator can pick by intent rather than by reading a column dictionary.

Every preset resolves to a column list from the frozen layout alone -- block names, loci and tier.
No corpus, no fitted artifact, no data of ours is needed to reproduce one, and two people selecting
the same preset get the same columns in the same order.

The ranking is the useful part:

``recommended``
    Use this unless you have a reason not to.
``specific``
    Correct for a stated purpose and wrong outside it.
``avoid``
    Present because it is a control, a baseline, or a measured dead end. Named so that choosing it
    is deliberate rather than accidental.

The rankings come from a benchmark over a public multi-study AIRR corpus: **182 study groups /
198 accessions, 14,553 samples**, scored with **study-disjoint folds** — fit on some studies,
predict on studies the fit never saw — so a column that merely encodes sequencing protocol scores
at chance. Anyone with a comparable SRA/AIRR corpus can reproduce the ranking; nothing here depends
on a private dataset. (An earlier revision said "several hundred study groups"; the counted figure
is 182, and the accession list is published in the analysis repo's ``heldout/``.)

Example:
    >>> from vdjtools.signature import presets
    >>> presets.get("transfer").rank
    'recommended'
    >>> cols = presets.columns("compact")          # a concrete column list
    >>> presets.table()                            # every preset, as a DataFrame
"""
from __future__ import annotations

from dataclasses import dataclass, field

from . import layout as L

#: Blocks that carry sequencing protocol rather than immunology: depth, presence masks and
#: call-quality fractions. Useful as covariates or as a control; not a feature set.
NUISANCE_BLOCKS = ("depth", "mask", "qc")

#: Blocks measured to carry the most study-to-study (batch) variance. Excluded from `transfer`.
#: These are the count-and-diversity statistics, whose absolute level moves with library prep.
HIGH_BATCH_BLOCKS = ("pair", "div", "clon", "iso")

#: The immunoglobulin loci. B-cell work lives here.
IG_LOCI = ("IGH", "IGK", "IGL")


@dataclass(frozen=True)
class Preset:
    """One named feature selection, with everything a user needs to judge it."""

    name: str
    rank: str                       # recommended | specific | avoid
    summary: str
    features: str                   # which columns, in words
    how: str                        # how they are computed
    use_cases: str
    tier: str = "standard"
    sig: tuple[str, ...] = ("vsig", "rsig")
    drop_blocks: tuple[str, ...] = ()
    keep_blocks: tuple[str, ...] = ()
    loci: tuple[str, ...] = ()
    scaling: str = "robust"         # the representation that measured best for this preset
    notes: str = ""
    _cache: dict = field(default_factory=dict, repr=False, compare=False)

    def columns(self) -> list[str]:
        """Resolve to a concrete, ordered column list."""
        out = []
        for c in L.columns(self.tier):
            sig, block, locus = L.parse(c)[0], L.parse(c)[1], (L.parse(c) + ("-",))[2]
            if sig not in self.sig:
                continue
            if self.keep_blocks and block not in self.keep_blocks:
                continue
            if block in self.drop_blocks:
                continue
            if self.loci and locus not in self.loci and locus != "-":
                continue
            out.append(c)
        return out

    @property
    def n_columns(self) -> int:
        return len(self.columns())


PRESETS: dict[str, Preset] = {}


def _add(p: Preset) -> Preset:
    PRESETS[p.name] = p
    return p


_add(Preset(
    name="compact",
    rank="recommended",
    summary="The smallest vector that still describes a repertoire. Start here.",
    features="Core tier, both halves: depth, diversity, clonality, CDR3 length, chain pairing, "
             "and the leading embedding coordinates per locus.",
    how="Counts and Hill numbers at a frozen coverage level, plus linear functionals of the "
        "prototype-sum embedding. Each column is transformed by its own variance-stabiliser and "
        "rescaled against the frozen reference.",
    use_cases="A first look; cohorts of tens to a few hundred samples; any model where the number "
              "of features must stay well under the number of samples.",
    tier="core",
    drop_blocks=NUISANCE_BLOCKS,
    notes="Width is a small multiple of the locus count, so it stays usable at n = 50.",
))

_add(Preset(
    name="classify",
    rank="recommended",
    summary="The general-purpose set. Best measured task performance when train and test come "
            "from comparable cohorts.",
    features="Standard tier, both halves, nuisance blocks removed.",
    how="Everything in `compact` plus amino-acid composition, physicochemistry, generation "
        "probability, isotype and the full per-locus embedding coordinates.",
    use_cases="Supervised classification and regression on a cohort of a few hundred samples or "
              "more; the default input to a random forest or gradient boosting.",
    drop_blocks=NUISANCE_BLOCKS,
    scaling="asinh",
    notes="Non-linear learners measured best on this set. Plain scaling beat PCA projection at "
          "every rank tested, so do not project it first.",
))

_add(Preset(
    name="transfer",
    rank="recommended",
    summary="For models that must work on another lab's samples. Drops the columns whose level "
            "moves most between studies.",
    features="Standard tier, both halves, minus the nuisance blocks and minus the "
             "count-and-diversity blocks (chain pairing, Hill diversity, clonality, isotype).",
    how="As `classify`, then excluding the blocks measured to carry the most study-to-study "
        "variance. What remains is dominated by composition and embedding geometry, which vary "
        "far more between donors than between studies.",
    use_cases="Cross-cohort prediction; meta-analysis over public data; any setting where the "
              "training and application cohorts were sequenced by different people.",
    drop_blocks=NUISANCE_BLOCKS + HIGH_BATCH_BLOCKS,
    notes="Buys robustness by giving up the diversity statistics, which are informative WITHIN a "
          "study. If train and test share a protocol, `classify` is the better choice.",
))

_add(Preset(
    name="geometry",
    rank="specific",
    summary="Embedding coordinates only — no count statistics at all.",
    features="The `rsig` half: per-locus embedding coordinates, their norms, band composition and "
             "the contrast against an unselected repertoire.",
    how="A weighted mean of fixed per-clonotype embedding vectors, rotated onto a per-locus basis "
        "derived from bundled reference receptors. No corpus is involved in the rotation.",
    use_cases="When batch is the main adversary; when comparing samples across tissue and blood; "
              "as the second view in a two-view model alongside `statistics`.",
    sig=("rsig",),
    drop_blocks=NUISANCE_BLOCKS,
    notes="Measured to carry the least study variance and the most donor variance of any block "
          "family, and to be almost unaffected by whether the sample is blood or tissue. It wins "
          "fewer supervised tasks outright than the statistics half — prefer it for robustness, "
          "not for raw accuracy.",
))

_add(Preset(
    name="statistics",
    rank="specific",
    summary="Classical repertoire statistics only. Needs no embedding, so vdjtools alone suffices.",
    features="The `vsig` half: diversity, clonality, CDR3 length, amino-acid composition, "
             "physicochemistry, generation probability, isotype, chain pairing.",
    how="Defined statistics of the clone-size vector and the germline vocabulary.",
    use_cases="When mirpy is not installed or the embedding is too costly; when every feature must "
              "have a textbook definition; as the first view in a two-view model.",
    sig=("vsig",),
    drop_blocks=NUISANCE_BLOCKS,
    notes="Wins the most supervised tasks of any single family, but carries several times more "
          "study-to-study variance than the geometry half. Check a batch label before trusting a "
          "between-cohort difference.",
))

_add(Preset(
    name="bcell",
    rank="specific",
    summary="B-cell receptor work: the immunoglobulin loci with somatic hypermutation and isotype.",
    features="Standard tier restricted to IGH/IGK/IGL, both halves, including the SHM and isotype "
             "blocks.",
    how="As `classify`, restricted to the Ig loci. The SHM column is a clone-weighted mean of "
        "V-region identity to germline and requires `v_identity` in the input (see `keep=` on the "
        "readers).",
    use_cases="Germinal-centre activity, class switching, ectopic lymphoid structure in tissue, "
              "any BCR-centric question.",
    loci=IG_LOCI,
    drop_blocks=NUISANCE_BLOCKS,
    notes="Pair with `vdjtools.stats.shm` for the mutation-level distribution, which a single mean "
          "cannot express. SHM is strongly depth-dependent — condition on depth or you measure "
          "library size.",
))

_add(Preset(
    name="full",
    rank="specific",
    summary="Every contract column. For feature selection, not for fitting.",
    features="The full tier, both halves, nothing removed.",
    how="The complete frozen layout.",
    use_cases="Feature-importance studies; building a new preset; regularised models that can "
              "handle far more columns than samples.",
    tier="full",
    notes="The ONLY preset that contains the nuisance blocks, because it is defined as every "
          "contract column. Every other preset drops them, so that a model trained on a preset "
          "cannot be reading sequencing depth. Measured to score WORSE than a few dozen "
          "well-chosen columns on real tasks: width is not free.",
))

_add(Preset(
    name="nuisance",
    rank="avoid",
    summary="Sequencing protocol only. A control, not a feature set.",
    features="Depth, presence masks and call-quality fractions.",
    how="Read and clonotype counts, per-locus presence flags, and the fraction of clones whose "
        "V/J call missed the germline index.",
    use_cases="As covariates to regress out; as the floor a real feature set must beat. If a model "
              "on this set matches your real model, your real model is reading library prep.",
    keep_blocks=NUISANCE_BLOCKS,
    tier="full",
    notes="Named and shipped precisely so it can be used as a control. Do not train a classifier "
          "on it and report the AUC as a finding.",
))


def get(name: str) -> Preset:
    """Look up a preset, with the valid names in the error."""
    if name not in PRESETS:
        raise KeyError(f"unknown preset {name!r}; available: {', '.join(sorted(PRESETS))}")
    return PRESETS[name]


def columns(name: str) -> list[str]:
    """The concrete column list for a preset."""
    return get(name).columns()


def table():
    """Every preset as a polars DataFrame — the documentation, machine-readable."""
    import polars as pl

    order = {"recommended": 0, "specific": 1, "avoid": 2}
    rows = [{"preset": p.name, "rank": p.rank, "columns": p.n_columns, "tier": p.tier,
             "halves": "+".join(p.sig), "scaling": p.scaling, "summary": p.summary,
             "features": p.features, "how": p.how, "use_cases": p.use_cases, "notes": p.notes}
            for p in PRESETS.values()]
    return pl.DataFrame(rows).sort(
        [pl.col("rank").replace_strict(order, default=3), pl.col("columns")])


def _demo() -> None:
    """Presets must resolve, be distinct, and never silently return nothing."""
    seen = {}
    for name, p in PRESETS.items():
        cols = p.columns()
        assert cols, f"{name} resolves to no columns"
        assert p.rank in ("recommended", "specific", "avoid"), (name, p.rank)
        assert all(c in set(L.columns("full")) for c in cols), f"{name} invented a column"
        seen[name] = tuple(cols)
    assert len(set(seen.values())) == len(seen), "two presets resolve identically"
    assert len(seen["compact"]) < len(seen["classify"]) < len(seen["full"])
    # `transfer` is a strict subset of `classify`: it only removes blocks.
    assert set(seen["transfer"]) < set(seen["classify"])
    # the halves partition the non-nuisance columns
    assert not set(seen["geometry"]) & set(seen["statistics"])
    # nuisance is disjoint from every recommended set
    assert not set(seen["nuisance"]) & set(seen["classify"])
    # `bcell` must carry no TCR-specific column, but the locus-agnostic ones stay: `pair` is a set
    # of cross-locus ratios with no single locus, and log(IGH/TRB) -- how B-dominated a sample is --
    # is exactly a B-cell readout, not a TCR one.
    tcr = [c for c in seen["bcell"] if any(f":{t}:" in c for t in ("TRA", "TRB", "TRG", "TRD"))]
    assert not tcr, f"bcell leaked TCR-specific columns: {tcr[:3]}"
    print("presets: " + ", ".join(f"{n}={len(c)}" for n, c in sorted(seen.items())))


if __name__ == "__main__":
    _demo()
