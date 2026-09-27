"""The signature column contract: raw feature groups, supports, pass-through channels.

A *signature* is a fixed-order, name-addressed feature vector for one sample. It is built in
three steps, and this module owns the vocabulary of all three while computing none of them::

    raw features  ->  rotation (per locus, from a corpus artifact)  ->  PC columns
                                                                    +  pass-through channels

Names are four colon-separated parts throughout, so one :func:`parse` reads every kind::

    vsig:div:TRB:1D_c        a raw feature      (an input to the rotation)
    vsig:pc:TRB:PC01         a rotated column   (an output)
    vsig:cov:TRB:cstar       a channel          (carried through untouched)
    vsig:qc:-:n_loci_present a channel that is not per-locus

Two things here exist because getting them wrong produced silent wrong answers before.

**``support`` is declared, never inferred.** Winsorization trims the tail a metric can run away
into, and which tail that is depends on the metric's support, *not* on its transform. Three raw
features declare ``transform="none"`` and have three different supports -- a log-probability is
bounded above by 0, a standard deviation is bounded below by 0, and a log-ratio is bounded
neither way. Deriving the trimming side from the transform silently trims the wrong end of two of
them, so the side comes from :data:`SUPPORTS` and nothing else.

**Channels are not rotated.** A provenance number that went through a rotation is no longer a
provenance number: it has been mixed with the measurements it was supposed to qualify. So the
fallback fractions, the coverage the sample actually attained, and the masks pass through in
their own units, and the rotation never sees them.

``vsig`` (statistics) is declared here. ``rsig`` (geometry) is computed in ``mir.signature`` and
registers its groups into this same registry with :func:`register_raw` / :func:`register_channel`;
nothing here imports ``mir``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

#: The seven human receptor loci, in canonical signature order.
LOCI: tuple[str, ...] = ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")

#: Placeholder in the locus slot of a column that is not per-locus.
NO_LOCUS = "-"

#: Block name reserved for rotated output. ``<sig>:pc:<locus>:PC01``.
PC_BLOCK = "pc"

#: Variance-stabilising transforms a raw feature may declare, applied where the feature is
#: computed, while its denominator is still in scope. See :mod:`vdjtools.signature.transform`.
TRANSFORMS: tuple[str, ...] = ("none", "log10", "log1p", "logit", "clr", "arcsine")

#: Supports a raw feature may declare, and the winsorization side each implies. This is the
#: closed vocabulary the trimming side is read from.
#:
#: ============ ================== ==================================================
#: support      trimmed            because
#: ============ ================== ==================================================
#: ``nonneg``   top only           ``[0, inf)`` -- cannot run away downward
#: ``nonpos``   bottom only        ``(-inf, 0]`` -- a log-probability, bounded by 0 above
#: ``real``     both               ``(-inf, inf)`` -- log-ratios, clr, logit, PC scores
#: ``unit``     neither            ``[a, b]`` -- already bounded on both sides
#: ============ ================== ==================================================
SUPPORTS: dict[str, tuple[bool, bool]] = {
    # support -> (trim_bottom, trim_top)
    "nonneg": (False, True),
    "nonpos": (True, False),
    "real": (True, True),
    "unit": (False, False),
}

#: Amino acids, anchored order. Also the k-mer alphabet.
AMINO_ACIDS: str = "ACDEFGHIKLMNPQRSTVWY"

#: Junction lengths the spectratype resolves, in residues. Fixed rather than data-driven so the
#: raw column set is the same for every sample and every corpus; anything outside falls in the
#: two saturating end bins.
SPECTRATYPE_LENGTHS: tuple[int, ...] = tuple(range(6, 32))

#: k for the k-mer group. 2, giving 400 columns per locus over the plain 20-letter alphabet.
#: k=3 is 8,000 cells against roughly 1,200 junction tokens at the corpus median depth -- over
#: 85% structural zeros, whose log-ratio coordinates read as a depth measurement wearing a motif
#: label. k=2 is ~1,200 observations over 400 parts, which is thin but genuinely estimated.
KMER_K: int = 2


def feats(transform: str, support: str, *names: str) -> dict[str, tuple[str, str]]:
    """``{name: (transform, support)}`` for a run of features sharing both.

    Groups are heterogeneous -- a diversity group carries log10 Hill numbers beside a logit
    evenness -- so both properties belong to the feature, not the group. This keeps the
    homogeneous runs terse anyway.
    """
    if transform not in TRANSFORMS:
        raise ValueError(f"unknown transform {transform!r}; known: {TRANSFORMS}")
    if support not in SUPPORTS:
        raise ValueError(f"unknown support {support!r}; known: {tuple(SUPPORTS)}")
    return {n: (transform, support) for n in names}


@dataclass(frozen=True)
class RawGroup:
    """One named family of raw features -- an input to the rotation, never an output.

    Args:
        sig: Owning signature, ``"vsig"`` or ``"rsig"``.
        name: Group name, unique within a ``sig``.
        features: ``{name: (transform, support)}``, most easily built with :func:`feats`.
            Empty for a group whose column names are not knowable without the germline; such a
            group declares ``dynamic`` instead.
        loci: Loci the group is emitted for. ``None`` means all of :data:`LOCI`; an empty tuple
            means the group is not per-locus and uses :data:`NO_LOCUS`.
        dynamic: For a group whose width depends on the germline vocabulary (``vus``, ``jus``,
            ``spec``): one ``(transform, support)`` pair shared by every column, with the names
            resolved at fit time and then **recorded in the corpus artifact**. Keeping the
            vocabulary out of here is what lets ``import vdjtools`` stay free of arda.
    """

    sig: str
    name: str
    features: dict[str, tuple[str, str]] = field(default_factory=dict)
    loci: tuple[str, ...] | None = None
    dynamic: tuple[str, str] | None = None

    def __post_init__(self) -> None:
        if self.sig not in ("vsig", "rsig"):
            raise ValueError(f"sig must be 'vsig' or 'rsig'; got {self.sig!r}")
        if bool(self.features) == bool(self.dynamic):
            raise ValueError(f"{self.sig}:{self.name} must declare exactly one of features= "
                             "(names known here) or dynamic= (names from the germline)")
        bad = {f: v for f, v in self.features.items()
               if v[0] not in TRANSFORMS or v[1] not in SUPPORTS}
        if bad:
            raise ValueError(f"{self.sig}:{self.name} has an unknown transform or support: {bad}")
        if self.dynamic and (self.dynamic[0] not in TRANSFORMS or self.dynamic[1] not in SUPPORTS):
            raise ValueError(f"{self.sig}:{self.name} dynamic={self.dynamic!r} is not a valid "
                             "(transform, support) pair")

    @property
    def emitted_loci(self) -> tuple[str, ...]:
        """Loci this group emits for; ``(NO_LOCUS,)`` when it is not per-locus."""
        if self.loci is None:
            return LOCI
        return self.loci or (NO_LOCUS,)

    def spec(self, feature: str) -> tuple[str, str]:
        """``(transform, support)`` for one of this group's features."""
        return self.dynamic if self.dynamic else self.features[feature]

    def columns(self, locus: str, names: "list[str] | None" = None) -> list[str]:
        """Raw column names for one locus.

        Args:
            locus: A member of :attr:`emitted_loci`.
            names: Required for a ``dynamic`` group -- the germline-resolved column names, in the
                order the artifact recorded them.
        """
        if self.dynamic:
            if names is None:
                raise ValueError(f"{self.sig}:{self.name} is dynamic; pass names= from the "
                                 "corpus artifact (or from the germline at fit time)")
            src = names
        else:
            src = list(self.features)
        return [f"{self.sig}:{self.name}:{locus}:{f}" for f in src]


@dataclass(frozen=True)
class Channel:
    """A named family carried through untouched -- never winsorized, never rotated.

    Args:
        sig: Owning signature.
        name: Channel name.
        features: ``{name: support}``. The support is recorded for documentation and for a
            range assertion; no trimming is applied either way.
        loci: As :class:`RawGroup`.
    """

    sig: str
    name: str
    features: dict[str, str]
    loci: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        if self.sig not in ("vsig", "rsig"):
            raise ValueError(f"sig must be 'vsig' or 'rsig'; got {self.sig!r}")
        bad = {f: s for f, s in self.features.items() if s not in SUPPORTS}
        if bad:
            raise ValueError(f"{self.sig}:{self.name} has an unknown support: {bad}")

    @property
    def emitted_loci(self) -> tuple[str, ...]:
        if self.loci is None:
            return LOCI
        return self.loci or (NO_LOCUS,)

    def columns(self, locus: str) -> list[str]:
        return [f"{self.sig}:{self.name}:{locus}:{f}" for f in self.features]


def _aa_pairs() -> tuple[str, ...]:
    return tuple(a + b for a in AMINO_ACIDS for b in AMINO_ACIDS)


# --------------------------------------------------------------------- the vsig declaration

_RAW: list[RawGroup] = [
    # -- clone-size distribution, coverage-standardised where it has to be -------------------
    RawGroup("vsig", "div", {
        **feats("log10", "nonneg", "1D_c", "0D_c", "2D_c", "0D_chao"),
        **feats("logit", "real", "clonality", "d50"),
    }),
    RawGroup("vsig", "depth", {
        **feats("log10", "nonneg", "reads", "richness"),
        **feats("log1p", "nonneg", "S_unseen"),
    }),
    RawGroup("vsig", "clon", {
        **feats("clr", "real", "f1", "f2"),
        **feats("logit", "real", "top"),
    }),
    # -- junction shape ----------------------------------------------------------------------
    RawGroup("vsig", "len", {
        **feats("none", "nonneg", "mean", "sd"),
        **feats("none", "real", "skew"),
    }),
    RawGroup("vsig", "aa", feats("arcsine", "nonneg", *AMINO_ACIDS)),
    RawGroup("vsig", "kmer", feats("arcsine", "nonneg", *_aa_pairs())),
    RawGroup("vsig", "pchem", feats("none", "real", *(
        f"{region}_{prop}" for region in ("all", "center")
        for prop in ("hydropathy", "charge", "polarity", "volume", "strength",
                     *(f"kf{i}" for i in range(1, 11)))))),
    # -- germline vocabulary; widths come from arda at fit time and the artifact thereafter --
    RawGroup("vsig", "vus", dynamic=("clr", "real")),
    RawGroup("vsig", "jus", dynamic=("clr", "real")),
    RawGroup("vsig", "spec", dynamic=("clr", "real")),
    # -- IGH only ----------------------------------------------------------------------------
    RawGroup("vsig", "iso", feats("clr", "real", "IgM", "IgD", "IgG", "IgA", "IgE"),
             loci=("IGH",)),
    RawGroup("vsig", "shm", feats("logit", "real", "mean_v_identity"), loci=("IGH",)),
    # -- cross-locus -------------------------------------------------------------------------
    RawGroup("vsig", "pair", feats("none", "real", "log_TRA_TRB", "log_TRG_TRB", "log_TRD_TRB",
                                   "log_IGK_IGL", "log_IGH_TRB"), loci=()),
]

_CHANNELS: list[Channel] = [
    # The sample's own attained Chao coverage, ALWAYS emitted -- including when the
    # coverage-standardised diversity features turn out not to be estimable, which is the case
    # it is most needed in. Previously computed, used to decide whether 28 columns were holes,
    # and then discarded: the only vsig quantity that was measured and thrown away.
    Channel("vsig", "cov", {"cstar": "unit"}),
    # Why a column is a hole. "This locus is absent", "this sample is too shallow to reach the
    # target" and "this field was never in the input" are three different facts that otherwise
    # all render as one nan.
    Channel("vsig", "mask", {"present": "unit", "estimable": "unit"}),
    Channel("vsig", "mask", {"c_call": "unit", "shm": "unit"}, loci=("IGH",)),
    # An unrecognised V or J call raises nowhere: the germline lookup falls back to the maximum
    # observed distance and the block comes out fully populated, plausible, and systematically
    # wrong. These fractions are the only thing that says so.
    Channel("vsig", "qc", {"v_fallback_frac": "unit", "j_fallback_frac": "unit",
                           "nonstd_aa_frac": "unit"}),
    Channel("vsig", "qc", {"n_loci_present": "nonneg", "winsor_frac": "unit"}, loci=()),
]


def register_raw(*groups: RawGroup) -> None:
    """Add raw groups to the registry. Used by ``mir.signature`` for the geometry half."""
    have = {(g.sig, g.name) for g in _RAW}
    for g in groups:
        if (g.sig, g.name) in have:
            raise ValueError(f"raw group {g.sig}:{g.name} is already registered")
        _RAW.append(g)


def register_channel(*channels: Channel) -> None:
    """Add pass-through channels to the registry."""
    have = {(c.sig, c.name, f) for c in _CHANNELS for f in c.features}
    for c in channels:
        clash = [f for f in c.features if (c.sig, c.name, f) in have]
        if clash:
            raise ValueError(f"channel {c.sig}:{c.name} already declares {clash}")
        _CHANNELS.append(c)


def raw_groups(sig: "str | None" = None) -> list[RawGroup]:
    """Registered raw groups, in declaration order."""
    return [g for g in _RAW if sig is None or g.sig == sig]


def channels(sig: "str | None" = None) -> list[Channel]:
    """Registered pass-through channels, in declaration order."""
    return [c for c in _CHANNELS if sig is None or c.sig == sig]


def raw_columns(sig: str, locus: str, vocab: "dict[str, list[str]] | None" = None) -> list[str]:
    """Every raw feature name for one ``(sig, locus)``, in rotation-input order.

    This is the column order of one row of the corpus matrix ``M_L``, so it is also the row order
    of the rotation. It must be reproduced exactly at apply time, which is why the corpus artifact
    stores it verbatim rather than recomputing it.

    Args:
        sig: ``"vsig"`` or ``"rsig"``.
        locus: A locus, or :data:`NO_LOCUS` for the cross-locus group.
        vocab: ``{group_name: [column names]}`` for every ``dynamic`` group emitted at this
            locus. Required unless no dynamic group applies.
    """
    out: list[str] = []
    for g in raw_groups(sig):
        if locus not in g.emitted_loci:
            continue
        names = None
        if g.dynamic:
            if vocab is None or g.name not in vocab:
                raise ValueError(f"raw_columns({sig!r}, {locus!r}) needs vocab[{g.name!r}]")
            names = vocab[g.name]
        out.extend(g.columns(locus, names))
    return out


def channel_columns(sig: str) -> list[str]:
    """Every pass-through channel column for one ``sig``, in emitted order."""
    return [c for ch in channels(sig) for loc in ch.emitted_loci for c in ch.columns(loc)]


def pc_columns(sig: str, k: dict[str, int]) -> list[str]:
    """Rotated column names, given the component count per locus.

    Args:
        sig: ``"vsig"`` or ``"rsig"``.
        k: ``{locus: n_components}``, from the corpus artifact. A locus absent from ``k`` emits
            nothing, which is how a corpus that could not fit a locus declares it.
    """
    return [f"{sig}:{PC_BLOCK}:{loc}:PC{i:02d}"
            for loc in (*LOCI, NO_LOCUS) if loc in k
            for i in range(1, k[loc] + 1)]


def signature_columns(sig: str, k: dict[str, int]) -> list[str]:
    """The emitted signature for one half: rotated columns then channels, in that order."""
    return [*pc_columns(sig, k), *channel_columns(sig)]


def parse(column: str) -> tuple[str, str, str, str]:
    """Split a column name into ``(sig, block, locus, feature)``.

    Also the definition of "is this a signature column", used as an allow-list wherever a frame
    may carry a caller's own joined columns -- so an ``age`` or ``n_reads`` column can never be
    silently winsorized or rotated.

    Raises:
        ValueError: If the name is not four colon-separated parts.
    """
    parts = column.split(":")
    if len(parts) != 4:
        raise ValueError(f"not a signature column: {column!r} (want sig:block:locus:feature)")
    return parts[0], parts[1], parts[2], parts[3]


def support_of(column: str, vocab: "dict[str, dict[str, list[str]]] | None" = None) -> str:
    """The declared support of a raw feature or channel column.

    Args:
        column: A raw or channel column name.
        vocab: Unused for lookup -- a dynamic group's support is declared on the group, not per
            column -- and accepted so callers can pass their artifact's vocab uniformly.

    Raises:
        ValueError: If no registered group or channel declares the column.
    """
    sig, block, _locus, feature = parse(column)
    if block == PC_BLOCK:
        return "real"
    for g in raw_groups(sig):
        if g.name != block:
            continue
        if g.dynamic:
            return g.dynamic[1]
        if feature in g.features:
            return g.features[feature][1]
    for c in channels(sig):
        if c.name == block and feature in c.features:
            return c.features[feature]
    raise ValueError(f"no registered group or channel declares {column!r}")


def _demo() -> None:
    """Self-check: the properties this contract exists to hold."""
    # every raw group declares exactly one of features= / dynamic=, enforced in __post_init__
    for g in raw_groups():
        assert bool(g.features) != bool(g.dynamic), g.name

    # the four supports map to four distinct trimming sides, and only 'real' trims both
    assert SUPPORTS["real"] == (True, True) and SUPPORTS["unit"] == (False, False)
    assert SUPPORTS["nonneg"] == (False, True) and SUPPORTS["nonpos"] == (True, False)

    # THE trap: transform does not determine the side. 'none' spans three supports, and the
    # three must trim differently.
    none_supports = {v[1] for g in raw_groups() for v in g.features.values() if v[0] == "none"}
    assert len(none_supports) > 1, "transform='none' must span more than one support"
    assert {SUPPORTS[s] for s in none_supports} != {SUPPORTS["real"]}

    # raw columns are per-locus and four-part; dynamic groups demand their vocabulary
    vocab = {"vus": ["TRBV1"], "jus": ["TRBJ1-1"], "spec": ["TRBV1_14"]}
    cols = raw_columns("vsig", "TRB", vocab)
    assert all(len(parse(c)) == 4 for c in cols)
    assert "vsig:div:TRB:1D_c" in cols and "vsig:vus:TRB:TRBV1" in cols
    assert "vsig:pair:TRB:log_TRA_TRB" not in cols, "cross-locus group leaked into a locus"
    assert raw_columns("vsig", NO_LOCUS) == [
        f"vsig:pair:-:{f}" for f in ("log_TRA_TRB", "log_TRG_TRB", "log_TRD_TRB",
                                     "log_IGK_IGL", "log_IGH_TRB")]
    try:
        raw_columns("vsig", "TRB")
    except ValueError as e:
        assert "vocab" in str(e)
    else:
        raise AssertionError("a dynamic group must refuse to guess its vocabulary")

    # IGH-only groups appear at IGH and nowhere else
    assert any(":iso:" in c for c in raw_columns("vsig", "IGH", vocab))
    assert not any(":iso:" in c for c in raw_columns("vsig", "TRB", vocab))

    # channels are never rotated, so they are disjoint from the pc block
    ch = channel_columns("vsig")
    assert "vsig:cov:TRB:cstar" in ch and "vsig:qc:-:winsor_frac" in ch
    assert not any(f":{PC_BLOCK}:" in c for c in ch)
    assert support_of("vsig:cov:TRB:cstar") == "unit"
    assert support_of("vsig:pc:TRB:PC01") == "real"
    assert support_of("vsig:len:TRB:sd") == "nonneg"
    assert support_of("vsig:len:TRB:skew") == "real"

    # a locus absent from k emits nothing for that locus -- how a corpus declares a failed fit
    pcs = pc_columns("vsig", {"TRB": 2, NO_LOCUS: 1})
    assert pcs == ["vsig:pc:TRB:PC01", "vsig:pc:TRB:PC02", "vsig:pc:-:PC01"], pcs
    assert signature_columns("vsig", {"TRB": 1})[0] == "vsig:pc:TRB:PC01"

    # parse is the allow-list: a caller's own column is not a signature column
    for bad in ("age", "vsig:div:TRB", "a:b:c:d:e"):
        try:
            parse(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"parse accepted {bad!r}")

    widths = {loc: len(raw_columns("vsig", loc, {"vus": [], "jus": [], "spec": []}))
              for loc in LOCI}
    print("layout OK", widths)


if __name__ == "__main__":
    _demo()
