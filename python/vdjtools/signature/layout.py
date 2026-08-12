"""The repertoire-signature column contract — names, tiers, transforms.

A *signature* is a fixed-width, fixed-order, name-addressed feature vector for one sample,
emitted on a scale a downstream model can consume without fitting a scaler. This module owns
the **contract** only: which columns exist, in which order, at which tier, and which
variance-stabilising transform each carries. It computes nothing.

Two signatures share the contract and concatenate on ``sample_id``:

* ``vsig`` — **statistics**, owned here: diversity, clonality, junction length, segment usage,
  isotype fractions, residue composition, physicochemistry, Pgen, SHM, locus balance, and the
  public-clonotype burden. Every column is a defined statistic of the clone-size vector or the
  germline vocabulary.
* ``rsig`` — **geometry**, owned by ``mir.signature`` and registered into this same registry.
  Every column is a linear functional, a norm, or a mixture coefficient of the prototype-sum
  measure. ``mir`` registers its blocks with :func:`register`; nothing here imports ``mir``.

``depth`` and ``div`` exist under **both** signatures on purpose — the count-native and the
embedding-native readings of the same idea are different objects, and their head-to-head is a
result rather than a redundancy. The ``sig`` prefix keeps them from colliding.

Column names are four colon-separated parts, always::

    <sig>:<block>:<locus>:<feature>        vsig:div:TRB:log1D_c
    <sig>:<block>:-:<feature>              vsig:pair:-:log_TRA_TRB     (not per-locus)

Tiers are nested: ``core`` ⊂ ``standard`` ⊂ ``full``, and each tier is an exact **index subset**
of the one frozen full-width layout, so a tier can be sliced out of a wider matrix without
recomputing anything.
"""
from __future__ import annotations

from dataclasses import dataclass

from ..features.physchem import DEFAULT_PROPERTIES as _PCHEM_PROPERTIES

#: The seven human receptor loci, in canonical signature order.
LOCI: tuple[str, ...] = ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")

#: Tiers, narrowest first. Membership is cumulative.
TIERS: tuple[str, ...] = ("core", "standard", "full")

#: Placeholder used in the locus slot of a column that is not per-locus.
NO_LOCUS = "-"

_TIER_RANK = {t: i for i, t in enumerate(TIERS)}


#: Variance-stabilising transforms a feature may declare. Applied where the feature is
#: computed, before any reference rescaling. See :mod:`vdjtools.signature.transform`.
TRANSFORMS: tuple[str, ...] = ("none", "log10", "log1p", "logit", "clr", "arcsine")


def feats(tier: str, transform: str, *names: str) -> dict[str, tuple[str, str]]:
    """``{name: (tier, transform)}`` for a run of features that share both.

    Blocks are frequently heterogeneous — a clonality block carries CLR-transformed count
    fractions beside a logit-transformed top-clone share — so the transform is a property of
    the feature, not of the block. This keeps the homogeneous runs terse anyway.
    """
    if transform not in TRANSFORMS:
        raise ValueError(f"unknown transform {transform!r}; known: {TRANSFORMS}")
    return {n: (tier, transform) for n in names}


@dataclass(frozen=True)
class Block:
    """One named feature family.

    Args:
        sig: Owning signature — ``"vsig"`` (statistics) or ``"rsig"`` (geometry).
        name: Block name; may be declared by several ``Block`` entries as long as their
            features are disjoint (e.g. the all-loci and IGH-only halves of ``mask``).
        features: ``{feature_name: (minimum_tier, transform)}``, most easily built with
            :func:`feats`. A feature appears in every tier at or above its own; insertion
            order is the emitted column order.
        loci: Loci this block is emitted for. ``None`` means all of :data:`LOCI`; an empty
            tuple means the block is not per-locus and uses :data:`NO_LOCUS`.
        attributable: Whether a column has a clonotype pre-image, i.e. whether asking "which
            clones drive this" is a well-posed question. Declared here at build time, never
            inferred from the name — a Hill number and a read fraction are summaries and have
            no pre-image, so asking is a category error rather than an unanswered question.
        exempt: Skip the frozen reference rescaling entirely (masks, which are already 0/1).
        magnitude: Rescale by one frozen scalar for the whole block, with no centring, because
            the block's *magnitude* is its signal. Per-column standardisation would force every
            coordinate to unit variance and make a near-zero sample indistinguishable from a
            typical one — it deletes exactly the deficiency the block exists to carry.
    """

    sig: str
    name: str
    features: dict[str, tuple[str, str]]
    loci: tuple[str, ...] | None = None
    attributable: bool = False
    exempt: bool = False
    magnitude: bool = False

    def __post_init__(self) -> None:
        if self.sig not in ("vsig", "rsig"):
            raise ValueError(f"sig must be 'vsig' or 'rsig'; got {self.sig!r}")
        bad = {f: v for f, v in self.features.items()
               if v[0] not in _TIER_RANK or v[1] not in TRANSFORMS}
        if bad:
            raise ValueError(f"{self.sig}:{self.name} features have an unknown tier or "
                             f"transform: {bad}; tiers {TIERS}, transforms {TRANSFORMS}")
        if self.exempt and self.magnitude:
            raise ValueError(f"{self.sig}:{self.name} cannot be both exempt and magnitude-scaled")

    def transform(self, feature: str) -> str:
        """The transform declared for one of this block's features."""
        return self.features[feature][1]

    @property
    def emitted_loci(self) -> tuple[str, ...]:
        """Loci this block emits for; ``(NO_LOCUS,)`` when the block is not per-locus."""
        if self.loci is None:
            return LOCI
        return self.loci or (NO_LOCUS,)

    def columns(self, tier: str = "full") -> list[str]:
        """Column names this block contributes at ``tier``, in emitted order."""
        rank = _TIER_RANK[tier]
        sel = [f for f, (t, _) in self.features.items() if _TIER_RANK[t] <= rank]
        return [f"{self.sig}:{self.name}:{loc}:{f}"
                for loc in self.emitted_loci for f in sel]


def _pcs(lo: int, hi: int, tier: str) -> dict[str, tuple[str, str]]:
    """``{PClo … PChi: (tier, "none")}`` — a contiguous run of principal-component features.

    Principal components carry no transform of their own: they are already a linear map of an
    embedding that was standardised against the frozen reference before rotation.
    """
    return feats(tier, "none", *(f"PC{i:02d}" for i in range(lo, hi + 1)))


# Retained dimensions per prototype slot. These are the *contract*, so they live in code rather
# than in the fitted artifact: a collaborator's column list must not move when the reference
# coefficients are re-fit.
#
# Set by measurement, not by taste (gate B1a, benchmark_signature_rotation.py). The rotation is
# the PCA of the bundled prototype cloud — fitted to no samples at all — so the question is how
# much *sample-level* variance it retains against a rotation fitted to the samples themselves.
# Each width below is the smallest clearing 0.98 of that ceiling on two independent cohorts with
# room to spare:
#
#   slot    k   cohort A  cohort B            k   cohort A  cohort B
#   phiv    8      0.981     0.969  fail     16      0.996     0.997  PASS
#   phij    4      0.950     0.881  fail      6      0.996     0.994  PASS
#   phic   16      0.959     0.974  fail     32      0.986     0.994  PASS
#
# The first-pass guesses (8 / 6 / 16) passed on one cohort and failed on the other, which is
# exactly why the gate runs on two. A width that falls short is cut and the signature version
# bumped — never silently widened, since that would redefine an already-shipped column.
PC_DIMS: dict[str, dict[str, int]] = {
    "phiv": {"standard": 16, "full": 24},
    "phij": {"standard": 6, "full": 12},
    "phic": {"standard": 32, "full": 48},
}


def _phi_block(name: str) -> Block:
    """One identity block — the norm at ``core``, then principal components by tier."""
    std, full = PC_DIMS[name]["standard"], PC_DIMS[name]["full"]
    return Block("rsig", name,
                 {**feats("core", "log1p", "norm"),
                  **_pcs(1, std, "standard"), **_pcs(std + 1, full, "full")},
                 attributable=True)


#: The registry, in emitted block order. ``rsig`` blocks are declared here beside the ``vsig``
#: ones so the contract is readable in one place, but they are *computed* by ``mir.signature``.
_BLOCKS: list[Block] = [
    # ---------------------------------------------------------------- provenance and coverage
    Block("vsig", "mask", feats("core", "none", "present", "estimable"), exempt=True),
    Block("vsig", "mask", feats("core", "none", "c_call", "shm"),
          loci=("IGH",), exempt=True),
    # An unrecognised V or J call does not raise: the germline lookup falls back to the maximum
    # observed distance, so a cohort using Adaptive nomenclature or an older IMGT release yields
    # a fully populated, entirely plausible, systematically wrong geometry block. The fallback
    # fraction is therefore a reported column, not a log line — it is the one number that tells a
    # collaborator their vector is not comparable to ours.
    Block("vsig", "qc", feats("core", "logit",
                              "v_fallback_frac", "j_fallback_frac", "nonstd_aa_frac")),
    Block("vsig", "qc", feats("core", "none", "n_loci_present"), loci=()),

    # ---------------------------------------------------------------------- count statistics
    Block("vsig", "depth", {**feats("core", "log10", "reads", "richness"),
                            **feats("standard", "log1p", "S_unseen")}),
    # Hill numbers standardised to a frozen coverage level, so 100 and 100,000 clonotypes are
    # comparable. estimate_d returns the linear qD; the log is applied here, not there.
    Block("vsig", "div", {**feats("core", "log10", "1D_c"),
                          **feats("standard", "log10", "0D_c", "2D_c"),
                          **feats("standard", "logit", "clonality"),
                          **feats("full", "log10", "0D_chao"),
                          **feats("full", "logit", "d50")}),
    Block("vsig", "clon", {**feats("core", "clr", "f1", "f2"),
                           **feats("core", "logit", "top")}),
    Block("vsig", "len", {**feats("core", "none", "mean", "sd"),
                          **feats("standard", "none", "skew")}),
    Block("vsig", "pair", feats("core", "none", "log_TRA_TRB", "log_TRG_TRB", "log_TRD_TRB",
                                "log_IGK_IGL", "log_IGH_TRB"), loci=()),

    # ------------------------------------------------------------------------ B-cell readouts
    # segment_usage drops a null c_call, and roughly two fifths of IGH reads carry none, so the
    # uncalled share is recovered against the sample total and closes the composition here.
    Block("vsig", "iso", {**feats("core", "clr", "IgM", "IgD", "IgG", "IgA"),
                          **feats("full", "clr", "IgE")}, loci=("IGH",)),
    Block("vsig", "shm", feats("standard", "logit", "mean_v_identity"), loci=("IGH",)),

    # ------------------------------------------------------------- recombination and sharing
    Block("vsig", "pgen", {**feats("standard", "none", "mean_log10", "sd_log10"),
                           **feats("standard", "logit", "frac_atypical")}),
    # The public-clonotype burden block is NOT here. It needs a frozen public-clonotype panel
    # built on a reference corpus, which does not exist yet; declared in a tier it contributed 28
    # permanently-nan columns to a 4,080-sample emission, indistinguishable to anyone downstream
    # from 28 columns their own samples were too shallow to support. A contract with dead columns
    # in it teaches people to ignore holes. `register()` is how it arrives once the panel ships.

    # ------------------------------------------------------------------ composition (full only)
    Block("vsig", "aa", feats("full", "arcsine", *"ACDEFGHIKLMNPQRSTVWY")),
    # Two regions only. The legacy five-region x four-quantile expansion emitted 100 numbers per
    # sample whose quantiles are unusable at the shallow depths this signature has to work at;
    # 'all' and 'center' are the two that survive there. Properties are vdjtools' own default
    # table, imported rather than re-listed so the two cannot drift apart.
    Block("vsig", "pchem", feats("full", "none", *(f"{region}_{prop}"
                                                   for region in ("all", "center")
                                                   for prop in _PCHEM_PROPERTIES))),

    # ================================================================= geometry, computed by mir
    # Embedding-native depth and diversity: n_eff is the Hill number of the clone-weight measure
    # the geometry actually uses, and mass is how much of the repertoire was ever drawn. Both
    # are properties of Phi, not of the count vector, which is why they sit beside — not instead
    # of — their vsig counterparts.
    Block("rsig", "depth", {**feats("core", "log10", "n_eff"),
                            **feats("core", "logit", "mass")}),
    # Rao's quadratic entropy of the clone-weight measure in embedding coordinates: sequence-
    # aware dispersion, which no Hill number can express because a Hill number cannot see that
    # two clonotypes are one substitution apart.
    Block("rsig", "div", feats("core", "log1p", "rao")),
    # Compartment shares recovered as mixture coefficients of Phi. Phi is linear in the
    # clone-weight measure, so Phi(S) = sum_c pi_c Phi(c) holds exactly for a partition; the
    # shares are measured, not assumed. Isotype appears here as a share of the geometry, which is
    # a different quantity from the vsig read fraction.
    Block("rsig", "band", feats("standard", "clr", "singleton", "top")),
    Block("rsig", "band", feats("standard", "clr", "IgM", "IgG", "IgA"), loci=("IGH",)),
    # Psi = mass * (Phi - naive): the signed deviation from unselected V(D)J output, scaled by
    # how much of the repertoire was actually observed. Magnitude-scaled, never centred — an
    # immune desert must land at the origin rather than at minus-the-median.
    Block("rsig", "contrast", {**feats("core", "log1p", "norm"),
                               **_pcs(1, 12, "standard"), **_pcs(13, 32, "full")},
          attributable=True, magnitude=True),
    _phi_block("phiv"),
    _phi_block("phij"),
    _phi_block("phic"),
]


def register(block: Block) -> None:
    """Add a block to the registry (used by ``mir.signature`` for late-bound vocabularies).

    Raises:
        ValueError: If a block with the same ``sig`` and ``name`` already declares any of the
            same features — a duplicate column would silently shift every index after it.
    """
    clash = {b.name for b in _BLOCKS
             if b.sig == block.sig and b.name == block.name
             and set(b.features) & set(block.features)
             and b.emitted_loci == block.emitted_loci}
    if clash:
        raise ValueError(f"{block.sig}:{block.name} already declares these features")
    _BLOCKS.append(block)


def registry(sig: str | None = None) -> list[Block]:
    """Registered blocks, optionally filtered to one signature."""
    return [b for b in _BLOCKS if sig is None or b.sig == sig]


def columns(tier: str = "standard", sig: str | None = None,
            blocks_: "tuple[str, ...] | None" = None) -> list[str]:
    """The column names of a signature, in emitted order.

    Args:
        tier: ``"core"``, ``"standard"`` or ``"full"``. Tiers are nested and each is an exact
            index subset of ``"full"``.
        sig: ``"vsig"``, ``"rsig"``, or ``None`` for both concatenated.
        blocks_: Restrict to these block names (both signatures unless ``sig`` is given).

    Returns:
        Column names as ``<sig>:<block>:<locus>:<feature>``.

    Raises:
        ValueError: If ``tier`` is unknown, or ``blocks_`` names a block that is not registered.
    """
    if tier not in _TIER_RANK:
        raise ValueError(f"tier must be one of {TIERS}; got {tier!r}")
    sel = registry(sig)
    if blocks_ is not None:
        known = {b.name for b in sel}
        unknown = [n for n in blocks_ if n not in known]
        if unknown:
            raise ValueError(f"unknown block(s) {unknown}; registered: {sorted(known)}")
        sel = [b for b in sel if b.name in blocks_]
    return [c for b in sel for c in b.columns(tier)]


def index(tier: str = "core", sig: str | None = None) -> list[int]:
    """Positions of ``tier``'s columns within the full-width layout of the same ``sig``.

    This is what makes the tiers cheap: emit ``full`` once, then slice.
    """
    wanted = set(columns(tier, sig))
    return [i for i, c in enumerate(columns("full", sig)) if c in wanted]


def describe(tier: str = "standard", sig: str | None = None):
    """The column dictionary — one row per column, so a collaborator can read the contract.

    Returns:
        A ``pl.DataFrame`` with ``column, sig, block, locus, feature, tier, transform,
        attributable, exempt, magnitude``.
    """
    import polars as pl

    rows = []
    for b in registry(sig):
        for loc in b.emitted_loci:
            for f, (t, tf) in b.features.items():
                if _TIER_RANK[t] > _TIER_RANK[tier]:
                    continue
                rows.append({"column": f"{b.sig}:{b.name}:{loc}:{f}", "sig": b.sig,
                             "block": b.name, "locus": loc, "feature": f, "tier": t,
                             "transform": tf, "attributable": b.attributable,
                             "exempt": b.exempt, "magnitude": b.magnitude})
    return pl.DataFrame(rows)


def parse(column: str) -> tuple[str, str, str, str]:
    """Split a column name into ``(sig, block, locus, feature)``.

    Raises:
        ValueError: If the name is not four colon-separated parts.
    """
    parts = column.split(":")
    if len(parts) != 4:
        raise ValueError(f"malformed signature column {column!r}: expected "
                         "'<sig>:<block>:<locus>:<feature>'")
    return parts[0], parts[1], parts[2], parts[3]


def _demo() -> None:
    """Self-check: tiers nest, tiers are index subsets, names round-trip, no duplicates."""
    for sig in ("vsig", "rsig", None):
        full = columns("full", sig)
        assert len(full) == len(set(full)), "duplicate column names in the layout"
        prev: list[str] = []
        for tier in TIERS:
            cols = columns(tier, sig)
            assert set(prev) <= set(cols), f"{tier} does not contain the previous tier"
            assert [full[i] for i in index(tier, sig)] == cols, f"{tier} is not an index subset"
            prev = cols
        for c in full:
            s, b, loc, f = parse(c)
            assert s == (sig or s) and loc in (*LOCI, NO_LOCUS), c

    v, r = columns("full", "vsig"), columns("full", "rsig")
    assert not set(v) & set(r), "vsig and rsig column names collide"
    assert columns("full") == v + r
    print(f"layout OK — vsig {len(columns('core', 'vsig'))}/{len(columns('standard', 'vsig'))}"
          f"/{len(v)}  rsig {len(columns('core', 'rsig'))}/{len(columns('standard', 'rsig'))}"
          f"/{len(r)}  (core/standard/full)")


if __name__ == "__main__":
    _demo()
