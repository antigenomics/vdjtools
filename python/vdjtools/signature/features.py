"""Raw vsig features and pass-through channels, for one sample.

Everything here is a **pure function of one sample plus frozen vocabulary**. Nothing is fitted, no
corpus is consulted, and no value is standardised, winsorized or rotated -- those need a corpus and
live in :mod:`vdjtools.signature.corpus`. That split is the point: two people who never share a
cohort compute the same raw features, and only the artifact they rotate through can differ.

Two axes this module keeps apart, because collapsing them is how wrong answers got shipped:

* **Raw features** feed the rotation. Each carries its declared transform, applied here while its
  denominator is still in scope, because a proportion separated from the count it was observed on
  cannot be stabilised afterwards.
* **Channels** never reach the rotation. A provenance number that has been mixed with the
  measurements it was meant to qualify is no longer provenance, so the fallback fractions, the
  masks, and the coverage the sample actually attained pass through in their own units.

The coverage target is a **runtime argument**, not a frozen constant. Carrying a per-sample quantity
in a corpus artifact is how one shipped reference ended up standardising tissue samples to a level
measured on blood, digit for digit across all seven loci.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..io.schema import (
    C_CALL,
    COUNT,
    FREQ,
    J_CALL,
    JUNCTION_AA,
    V_CALL,
    resolve_duplicates,
    strip_allele,
)
from . import transform as T
from .layout import AMINO_ACIDS, KMER_K, LOCI, SPECTRATYPE_LENGTHS

#: The 20 proteinogenic amino acids, anchored. A junction containing anything else -- a stop codon,
#: an ambiguity code, the legacy out-of-frame marker -- is dropped before anything is computed. Not
#: a crash guard: most downstream code accepts a malformed junction silently and returns a finite,
#: meaningless number, so this filter is the only thing between a contaminated feature and a
#: plausible-looking wrong answer. The dropped fraction is reported as ``qc:*:nonstd_aa_frac``.
VALID_AA = r"^[ACDEFGHIKLMNPQRSTVWY]+$"

#: The 20 amino acids plus the two characters a WELL-FORMED AIRR table can legitimately carry:
#: ``*`` for a stop codon, ``_`` for the legacy out-of-frame marker. Anything outside this is a
#: broken file, not a kind of receptor -- see :func:`sanitise`.
PARSEABLE_AA = r"^[ACDEFGHIKLMNPQRSTVWY*_]+$"

#: Isotype classes and the constant-gene names that map onto them. ``IGHGP`` is a pseudogene and
#: ``IGHC`` is ambiguous, so neither is called.
ISOTYPES: dict[str, tuple[str, ...]] = {
    "IgM": ("IGHM",), "IgD": ("IGHD",),
    "IgG": ("IGHG1", "IGHG2", "IGHG3", "IGHG4"),
    "IgA": ("IGHA1", "IGHA2"), "IgE": ("IGHE",),
}

#: Clone-size weights ``g``, the **one** definition both halves of the signature use. Concave by
#: default: raw read weighting lets a single dominant clone be the entire profile, presence
#: weighting throws the expansion signal away, and ``log2(1+a)`` sits between.
#:
#: It lives in the shared dependency because ``vsig`` and ``rsig`` must weight the same repertoire
#: identically or they are not describing one measure -- and a second copy that drifted would not
#: raise anywhere, it would quietly make the two halves disagree.
WEIGHTS = {
    "log2p1": lambda a: np.log2(1.0 + a),
    "log1p": np.log1p,
    "anscombe": lambda a: np.sqrt(a + 0.375),
    "duplicate_count": lambda a: np.asarray(a, dtype=float),
    "distinct": lambda a: np.ones(np.shape(a), dtype=float),
}

#: Locus pairs whose read-count ratio is a compartment read-out: T-vs-B balance, the gamma-delta
#: share, the light-chain balance. A ratio rather than two counts, because the sequencing depth
#: that drives both cancels.
PAIRS: tuple[tuple[str, str], ...] = (
    ("TRA", "TRB"), ("TRG", "TRB"), ("TRD", "TRB"), ("IGK", "IGL"), ("IGH", "TRB"),
)


# ------------------------------------------------------------------------------ input hygiene


def assert_parseable(df: pl.DataFrame) -> None:
    """Raise if any ``junction_aa`` is not a parseable amino-acid string.

    A **data-integrity gate, not a filter**, and the first of three axes to keep apart:
    *parseable* (can we read the string at all -- here), *productive* (does the rearrangement
    encode a chain -- AIRR), *functional* (is the germline gene real -- IMGT F/ORF/P).

    Measured across 6,047,716 rows of a clinical AIRR store: zero violations. The gate therefore
    costs nothing on well-formed data and is only ever reached by a real problem.

    Raises:
        ValueError: If any junction carries an unparseable character.
    """
    if df.height == 0:
        return
    bad = df.filter(pl.col(JUNCTION_AA).is_not_null()
                    & ~pl.col(JUNCTION_AA).str.contains(PARSEABLE_AA))
    if bad.height:
        ex = bad[JUNCTION_AA].unique().head(5).to_list()
        raise ValueError(
            f"{JUNCTION_AA} has {bad.height} of {df.height} unparseable value(s): outside the 20 "
            f"amino acids and outside '*' and '_', e.g. {ex}. A well-formed AIRR table carries "
            "nothing else, so this is a data problem rather than a filtering decision. Fix the "
            "source, or call sanitise(strict=False) to drop them.")


def sanitise(df: pl.DataFrame, *, strict: bool = True,
             on_duplicate: str = "error") -> tuple[pl.DataFrame, float]:
    """Drop non-productive clonotypes; return the frame and the dropped **weight** fraction.

    Dropped by weight, not by row: losing one dominant clone matters more than losing fifty
    singletons, and a row fraction would hide that.

    Args:
        df: A clonotype frame.
        strict: Raise on unparseable characters (see :func:`assert_parseable`). ``False`` drops
            them, for a corpus known to carry ambiguity codes.
        on_duplicate: What to do when the frame has no ``junction_nt`` and repeats an amino-acid
            clonotype key -- ``"error"`` (default) or ``"sum"``. The default refuses because the
            two readings give different richness and clonality and the frame cannot say which is
            meant.

    Returns:
        ``(kept_frame, dropped_weight_fraction)``.
    """
    if df.height == 0:
        return df, 0.0
    if strict:
        assert_parseable(df)
    total = float(df[COUNT].sum())
    keep = df.filter(
        (pl.col(COUNT) > 0)
        & pl.col(JUNCTION_AA).is_not_null()
        & pl.col(JUNCTION_AA).str.contains(VALID_AA)      # anchored: a match, not a search
    )
    # AFTER the filter, not before. "Are these two rows one clonotype or two?" is only worth
    # asking about rows that reach the estimators. A frame whose duplicates are all
    # non-productive has junk to drop, not ambiguity to resolve, and raising on it would turn a
    # locus that is supposed to mask out into an error.
    keep = resolve_duplicates(keep, on_duplicate)
    kept = float(keep[COUNT].sum()) if keep.height else 0.0
    return keep, (1.0 - kept / total) if total > 0 else 0.0


def work_frame(df: pl.DataFrame, weight: str = "log2p1") -> pl.DataFrame:
    """Overwrite ``frequency`` with the normalised clone weight, so every profiler agrees.

    **Order matters.** Call this *after* filtering, and never call ``filter_functional``,
    ``downsample`` or ``select_top`` afterwards -- each recomputes ``frequency`` from the counts
    and would silently restore read weighting.
    """
    if weight not in WEIGHTS:
        raise ValueError(f"unknown weight {weight!r}; known: {sorted(WEIGHTS)}")
    a = df[COUNT].to_numpy().astype(float)
    g = WEIGHTS[weight](a)
    s = g.sum()
    return df.with_columns(pl.Series(FREQ, g / s if s > 0 else g))


def _w(df: pl.DataFrame) -> np.ndarray:
    return df[FREQ].to_numpy().astype(float)


#: ASCII byte -> index into :data:`AMINO_ACIDS`. Built once, at import. This table is what makes
#: the composition groups one numpy gather instead of a Python loop over every residue of every
#: clonotype: at the corpus depth that loop ran 1.2 million times per sample.
_AA_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _a in enumerate(AMINO_ACIDS):
    _AA_LUT[ord(_a)] = _i


def _residues(df: pl.DataFrame) -> "tuple[np.ndarray, np.ndarray, np.ndarray]":
    """Every junction residue in the frame, flattened: ``(code, weight, length-per-clonotype)``.

    ``sanitise`` has already anchored the alphabet to the 20 proteinogenic residues, so every code
    is in ``0..19`` and no validity mask is needed here. The per-residue weight is the clonotype's
    weight repeated across its residues, which is exactly what a weighted composition sums.
    """
    seqs = df[JUNCTION_AA].to_list()
    lens = np.fromiter(map(len, seqs), dtype=np.int64, count=len(seqs))
    codes = _AA_LUT[np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8)]
    return codes, np.repeat(_w(df), lens), lens


def _property_matrix() -> np.ndarray:
    """The legacy amino-acid property table as ``(20, len(DEFAULT_PROPERTIES))``, residue-indexed.

    Row order is :data:`AMINO_ACIDS`, so it is indexed directly by the codes :func:`_residues`
    returns. Reads the shipped table through its own cached loader rather than re-parsing it.
    """
    from ..features.physchem import DEFAULT_PROPERTIES, load_property_table

    tbl = load_property_table()
    by_aa = {r["amino_acid"]: r for r in tbl.iter_rows(named=True)}
    return np.array([[float(by_aa[a][p]) for p in DEFAULT_PROPERTIES] for a in AMINO_ACIDS])


def _codes(series: pl.Series, mapping: dict, default: int) -> np.ndarray:
    """Map a gene-call series onto integer indices in one polars pass, ``default`` for unknown."""
    return series.replace_strict(mapping, default=default,
                                 return_dtype=pl.Int64).to_numpy()


# ------------------------------------------------------------------------ germline vocabulary


def gene_vocab(locus: str, organism: str = "human") -> dict[str, list[str]]:
    """Column names for the three germline-width raw groups at one locus.

    arda's germline is the single source of truth for V/J vocabulary across this ecosystem, so the
    usage and spectratype widths are a property of the germline release rather than of any cohort.
    The corpus artifact records what came back here, and apply time reads it from the artifact --
    never from arda again, because a germline update would silently re-index the rotation.

    Returns:
        ``{"vus": [...], "jus": [...], "spec": [...]}``. ``spec`` is ``<V gene>_<length>`` over
        :data:`~vdjtools.signature.layout.SPECTRATYPE_LENGTHS`, V-major.
    """
    from ..model.reference import load_germline

    germ = load_germline(locus, organism=organism)
    v = sorted(set(germ.filter(pl.col("segment") == "V")["gene"].to_list()))
    j = sorted(set(germ.filter(pl.col("segment") == "J")["gene"].to_list()))
    return {"vus": v, "jus": j,
            "spec": [f"{g}_{n}" for g in v for n in SPECTRATYPE_LENGTHS]}


# ---------------------------------------------------------------------------- raw feature groups


def div_group(df: pl.DataFrame, level: float) -> dict[str, float]:
    """Hill numbers standardised to a coverage ``level``, plus distribution shape.

    Coverage standardisation is what makes a hundred clonotypes and a hundred thousand comparable:
    both are read at the same completeness rather than at their own depth. It only works where the
    sample actually reaches ``level``; beyond that the estimator extrapolates, and extrapolation is
    not a mild approximation here. Measured on one real repertoire subsampled across a 200x depth
    range: at a level every subsample attained, Shannon diversity agreed within 4%; at a level the
    shallow ones extrapolated to, the same statistic inflated roughly tenfold.

    Whether *this* sample reached it is reported by ``mask:*:estimable``, and the level it actually
    attained by ``cov:*:cstar``. When it did not reach the target the features are ``nan`` -- a hole
    a downstream model can see, rather than a confident wrong number.
    """
    from ..stats.inext import estimate_d

    keys = ("1D_c", "0D_c", "2D_c", "clonality", "0D_chao", "d50")
    nan = dict.fromkeys(keys, np.nan)
    counts = df[COUNT].to_numpy().astype(np.int64)
    if counts.size < 2 or counts.sum() < 2:
        return nan
    try:
        est = estimate_d(counts, base="coverage", level=level, q=(0, 1, 2), se=False)
    except (ValueError, ZeroDivisionError):
        return nan
    if not estimable(df, level, est):
        return nan

    qd = {int(r["order_q"]): float(r["qD"]) for r in est.to_dicts()}
    # Clonality from the *standardised* Hill numbers, not from observed evenness. Observed
    # evenness is a ratio of two depth-dependent quantities and inherits their depth dependence in
    # full: measured on one repertoire across a 667x depth range it drifted by nearly a factor of
    # two, while the same quantity built from coverage-standardised numbers is flat.
    evenness = np.log(qd[1]) / np.log(qd[0]) if qd[0] > 1 else 0.0
    f1, f2 = float((counts == 1).sum()), float((counts == 2).sum())
    chao = counts.size + (f1 * f1 / (2.0 * f2) if f2 > 0 else f1 * (f1 - 1.0) / 2.0)
    srt = np.sort(counts)[::-1]
    k = int(np.searchsorted(np.cumsum(srt), 0.5 * counts.sum()) + 1)
    return {"1D_c": T.log10(qd[1]), "0D_c": T.log10(qd[0]), "2D_c": T.log10(qd[2]),
            "clonality": T.logit(np.clip(1.0 - evenness, 0.0, 1.0), counts.size),
            "0D_chao": T.log10(chao), "d50": T.logit(k / counts.size, counts.size)}


def estimable(df: pl.DataFrame, level: float, est: "pl.DataFrame | None" = None) -> bool:
    """Whether a coverage-standardised diversity is a measurement rather than an extrapolation.

    Two necessary conditions: the estimator must not have extrapolated, and the target depth must
    be within twice the observed one. The second catches the case where it interpolates nominally
    but is leaning on almost no data.
    """
    from ..stats.inext import estimate_d

    counts = df[COUNT].to_numpy().astype(np.int64)
    if counts.size < 2 or counts.sum() < 2:
        return False
    if est is None:
        try:
            est = estimate_d(counts, base="coverage", level=level, q=(1,), se=False)
        except (ValueError, ZeroDivisionError):
            return False
    rows = est.to_dicts()
    return (all(r["method"] != "extrapolation" for r in rows)
            and max(float(r["m"]) for r in rows) <= 2.0 * counts.sum())


def depth_group(df: pl.DataFrame, n_reads: float, richness: float) -> dict[str, float]:
    """How much was seen, and how much was not.

    ``S_unseen`` is Chao's estimate of the clonotypes that exist but were never drawn. Carried
    explicitly rather than folded into a diversity estimate because it is the honest statement of
    what the sample could not resolve.
    """
    a = df[COUNT].to_numpy()
    f1, f2 = float((a == 1).sum()), float((a == 2).sum())
    s_unseen = f1 * f1 / (2.0 * f2) if f2 > 0 else f1 * (f1 - 1.0) / 2.0
    return {"reads": T.log10(n_reads), "richness": T.log10(richness),
            "S_unseen": T.log1p(max(s_unseen, 0.0))}


def clon_group(df: pl.DataFrame, n_reads: float, richness: float) -> dict[str, float]:
    """Clone-size distribution shape as a composition, plus the top clone's share.

    ``f1``/``f2``/``f3plus`` -- seen once, twice, more -- is a three-part composition read in
    log-ratio coordinates. Two of the three ship because the third is exactly determined by them.
    """
    a = df[COUNT].to_numpy()
    n = max(richness, 1.0)
    parts = {"f1": float((a == 1).sum()), "f2": float((a == 2).sum()),
             "f3plus": float((a >= 3).sum())}
    c = T.clr(parts, m=n)
    top = float(a.max()) / n_reads if n_reads > 0 and a.size else 0.0
    return {"f1": c["f1"], "f2": c["f2"], "top": T.logit(top, n_reads)}


def len_group(df: pl.DataFrame) -> dict[str, float]:
    """Weighted moments of the junction length distribution.

    Moments rather than a pooled histogram: the per-V-gene histogram is the ``spec`` group's job,
    and a locus-pooled 26-bin histogram is neither of the two useful things.
    """
    if df.height == 0:
        return dict.fromkeys(("mean", "sd", "skew"), np.nan)
    ln = df[JUNCTION_AA].str.len_chars().cast(pl.Float64).to_numpy()
    w = _w(df)
    mean = float((w * ln).sum())
    sd = float(np.sqrt((w * (ln - mean) ** 2).sum()))
    return {"mean": mean, "sd": sd,
            "skew": float((w * (ln - mean) ** 3).sum() / sd ** 3) if sd > 0 else 0.0}


def aa_group(df: pl.DataFrame) -> dict[str, float]:
    """Weighted single-residue composition of the junctions, arcsine-stabilised."""
    keys = list(AMINO_ACIDS)
    if df.height == 0:
        return dict.fromkeys(keys, np.nan)
    codes, wr, lens = _residues(df)
    acc = np.bincount(codes, weights=wr, minlength=20)
    m = float(lens.sum())                                 # residues actually observed
    tot = acc.sum()
    if tot <= 0:
        return dict.fromkeys(keys, np.nan)
    return dict(zip(keys, map(float, T.arcsine(acc / tot, m))))


def kmer_group(df: pl.DataFrame) -> dict[str, float]:
    """Weighted k-mer composition of the junctions, arcsine-stabilised.

    ``k=2`` over the plain 20-letter alphabet: 400 parts against roughly 1,200 junction tokens at
    the corpus median depth. Thin, but genuinely estimated -- ``k=3`` would be 8,000 cells, over
    85% of them structural zeros, whose coordinates read as a depth measurement wearing a motif
    label.
    """
    pairs = [a + b for a in AMINO_ACIDS for b in AMINO_ACIDS]
    if df.height == 0:
        return dict.fromkeys(pairs, np.nan)
    n = len(AMINO_ACIDS)
    codes, wr, lens = _residues(df)
    if codes.size < KMER_K:
        return dict.fromkeys(pairs, np.nan)
    # A k-mer has to sit inside ONE junction. The flattened residue array concatenates them, so the
    # window starting at any clonotype's last residue would straddle into the next clonotype: mask
    # exactly those starts out. (``cumsum(lens) - 1`` is each clonotype's last residue; the final
    # one has no successor in the window array at all.)
    ok = np.ones(codes.size - KMER_K + 1, dtype=bool)
    ok[np.cumsum(lens)[:-1] - 1] = False
    idx = (codes[:-1].astype(np.int64) * n + codes[1:])[ok]
    acc = np.bincount(idx, weights=wr[:-1][ok], minlength=n * n)
    n_tok = float(ok.sum())
    tot = acc.sum()
    if tot <= 0:
        return dict.fromkeys(pairs, np.nan)
    return dict(zip(pairs, map(float, T.arcsine(acc / tot, max(n_tok, 1.0)))))


def pchem_group(df: pl.DataFrame, regions=("all", "center")) -> dict[str, float]:
    """Weighted mean physicochemistry of the junction, over two regions.

    Two regions rather than five, and means rather than four quantiles apiece: at a hundred
    clonotypes the quantiles of a per-clonotype property are noise, while the weighted mean is a
    well-behaved average over every residue seen.
    """
    from ..features.physchem import DEFAULT_PROPERTIES

    keys = [f"{r}_{p}" for r in regions for p in DEFAULT_PROPERTIES]
    out = dict.fromkeys(keys, np.nan)
    if df.height == 0:
        return out
    props = _property_matrix()
    codes, _wres, lens = _residues(df)
    w = _w(df)
    starts = np.concatenate(([0], np.cumsum(lens)[:-1]))
    for region in regions:
        if region == "all":
            # Per-clonotype mean over its own residues, then the weighted mean over clonotypes.
            # ``reduceat`` on the flattened residue codes is the whole computation; the polars path
            # this replaces exploded one row per residue (150,000 rows on a 10,000-clonotype
            # locus) and grouped twice, at 9 ms per call against 0.3 ms here.
            per = np.add.reduceat(props[codes], starts, axis=0) / lens[:, None]
            wi = w
        elif region == "center":
            # The middle five residues; a junction shorter than five has no centre and is skipped,
            # exactly as ``_region_expr`` returns null for it.
            ok = lens >= 5
            if not ok.any():
                continue
            take = (starts[ok] + lens[ok] // 2 - 2)[:, None] + np.arange(5)
            per = props[codes[take]].mean(axis=1)
            wi = w[ok]
        else:
            raise ValueError(f"region must be 'all' or 'center'; got {region!r}")
        tot = wi.sum()
        if tot <= 0:
            continue
        vals = (per * wi[:, None]).sum(axis=0) / tot
        out |= {f"{region}_{k}": float(v) for k, v in zip(DEFAULT_PROPERTIES, vals)}
    return out


def usage_group(df: pl.DataFrame, col: str, genes: list[str]) -> dict[str, float]:
    """Weighted gene-usage composition over a fixed germline vocabulary, in clr coordinates.

    The vocabulary is fixed by the artifact, not by the sample, so a gene this sample never used is
    a structural zero of a known composition rather than a missing column. An **unrecognised** call
    is different and is not silently folded in here -- it lands in the closing residual part and is
    reported by ``qc:*:{v,j}_fallback_frac``.

    The clr is taken over the vocabulary **plus that residual**, then the residual is dropped: a clr
    of a sub-composition is a different number from the corresponding coordinate of the full one.
    """
    if df.height == 0:
        return dict.fromkeys(genes, np.nan)
    pos = {g: i for i, g in enumerate(genes)}
    calls = df.select(strip_allele(pl.col(col)).alias("g"))["g"]
    gi = _codes(calls, pos, len(genes))                 # last cell closes the composition
    acc = np.bincount(gi, weights=_w(df), minlength=len(genes) + 1)
    coords = T.clr(acc, m=df.height)
    return {g: float(coords[i]) for g, i in pos.items()}


def spec_group(df: pl.DataFrame, names: list[str], v_genes: list[str]) -> dict[str, float]:
    """Spectratype: the junction-length distribution resolved per V gene, in clr coordinates.

    A (V gene x length) table, closed by one residual cell for calls or lengths outside the
    vocabulary. This is the widest group by a wide margin -- 121 IGH V genes x 26 lengths is 3,146
    cells against a median 464 IGH clonotypes per real sample, so a single sample's table is ~0.15
    observations per cell and the great majority of it is structural zeros.

    That sparsity is real and is not hidden: a **corpus** of 10,000 repertoires puts ~4.6M
    observations behind the same cells, so the rotation's directions are well determined even where
    one sample's scores along them are noisy. How noisy is readable per sample from
    ``cov:*:cstar``, which is why that channel is always emitted.
    """
    if df.height == 0:
        return dict.fromkeys(names, np.nan)
    nl = len(SPECTRATYPE_LENGTHS)
    lo, hi = SPECTRATYPE_LENGTHS[0], SPECTRATYPE_LENGTHS[-1]
    vpos = {g: i for i, g in enumerate(v_genes)}
    resid = len(v_genes) * nl                           # last cell closes the composition
    vi = _codes(df.select(strip_allele(pl.col(V_CALL)).alias("g"))["g"], vpos, -1)
    # int64, not the UInt32 ``len_chars`` returns: the clip below subtracts ``lo`` and an unsigned
    # length underflows on a junction shorter than the shortest spectratype bin.
    ln = df[JUNCTION_AA].str.len_chars().to_numpy().astype(np.int64)
    cell = np.where(vi >= 0, vi * nl + np.clip(ln, lo, hi) - lo, resid)
    acc = np.bincount(cell, weights=_w(df), minlength=resid + 1)
    coords = T.clr(acc, m=df.height)
    return {n: float(coords[i]) for i, n in enumerate(names)}


def iso_group(df: pl.DataFrame) -> dict[str, float]:
    """Isotype composition of an IGH repertoire, in log-ratio coordinates.

    The uncalled share is a real part of the composition, not a rounding error -- roughly two
    fifths of IGH reads carry no constant-gene call -- so it closes the composition rather than
    being silently dropped the way a usage profile would drop it.
    """
    keys = list(ISOTYPES)
    if df.height == 0:
        return dict.fromkeys(keys, np.nan)
    w = _w(df)
    # Strip the allele first. The classes are gene names matched by equality, so a frame calling
    # ``IGHG1*01`` would match nothing and come back 100% uncalled -- a plausible composition
    # (some cohorts really are mostly uncalled) rather than an error anybody sees.
    calls = (df.select(strip_allele(pl.col(C_CALL).cast(pl.Utf8)).alias("c"))["c"].to_list()
             if C_CALL in df.columns else [None] * df.height)
    gene_iso = {g: i for i, (_iso, names) in enumerate(ISOTYPES.items()) for g in names}
    ci = _codes(pl.Series("c", calls, dtype=pl.Utf8), gene_iso, len(ISOTYPES))
    got = np.bincount(ci, weights=w, minlength=len(ISOTYPES) + 1)
    parts = {iso: float(got[i]) for i, iso in enumerate(ISOTYPES)}
    parts["_uncalled"] = max(1.0 - sum(parts.values()), 0.0)
    coords = T.clr(parts, m=df.height)
    return {k: coords[k] for k in keys}


def shm_group(df: pl.DataFrame) -> dict[str, float]:
    """Mean somatic hypermutation load, via ``v_identity``.

    Absent from the canonical schema -- vdjtools' readers narrow to eight columns and
    ``v_identity`` is not one of them -- so this masks out unless the caller kept it explicitly.
    """
    if df.height == 0 or "v_identity" not in df.columns:
        return {"mean_v_identity": np.nan}
    v = df["v_identity"].to_numpy().astype(float)
    ok = np.isfinite(v)
    if not ok.any():
        return {"mean_v_identity": np.nan}
    w = _w(df)[ok]
    mean = float((w * v[ok]).sum() / w.sum()) if w.sum() > 0 else np.nan
    return {"mean_v_identity": T.logit(np.clip(mean, 0.0, 1.0), int(ok.sum()))}


def pair_group(reads: dict[str, float]) -> dict[str, float]:
    """Log read-count ratios between loci -- compartment balance, with depth divided out."""
    return {f"log_{a}_{b}": float(np.log10((reads.get(a, 0.0) + 1.0) / (reads.get(b, 0.0) + 1.0)))
            for a, b in PAIRS}


# -------------------------------------------------------------------------- channel functions


def coverage_of(df: pl.DataFrame) -> float:
    """The sample's own attained Chao coverage at its observed depth, in ``[0, 1]``."""
    from ..stats.inext import sample_coverage

    counts = df[COUNT].to_numpy().astype(np.int64)
    if counts.size == 0 or counts.sum() <= 0:
        return float("nan")
    return float(sample_coverage(counts))


def qc_channel(raw: pl.DataFrame, clean: pl.DataFrame, locus: str,
               nonstd_frac: float) -> dict[str, float]:
    """Whether this sample's gene calls are in a vocabulary we recognise.

    An unrecognised V or J call raises nowhere downstream: the germline distance lookup falls back
    to the maximum observed distance, so a cohort on Adaptive nomenclature or an older IMGT release
    yields a fully populated, entirely plausible, systematically wrong result. These fractions are
    the one number telling a collaborator their vector is not comparable to ours, which is why they
    are columns and not a warning nobody reads.

    Reported as raw fractions in ``[0, 1]``, not logit-transformed: a channel is read by a human
    deciding whether to trust the row, and a logit is not.
    """
    from ..model.reference import load_germline

    out = {"nonstd_aa_frac": float(nonstd_frac)}
    try:
        germ = load_germline(locus)
    except Exception:                       # a locus with no bundled germline: unmeasurable
        return {**out, "v_fallback_frac": np.nan, "j_fallback_frac": np.nan}
    known = {seg: set(germ.filter(pl.col("segment") == seg)["gene"].to_list())
             for seg in ("V", "J")}
    for seg, col in (("V", V_CALL), ("J", J_CALL)):
        if not clean.height:
            out[f"{seg.lower()}_fallback_frac"] = np.nan
            continue
        w = _w(clean)
        genes = clean.select(strip_allele(pl.col(col)).alias("g"))["g"]
        # ``is_in`` leaves a null call null; an uncalled gene IS a fallback, so it counts as missing.
        miss = (~genes.is_in(list(known[seg]))).fill_null(True).to_numpy()
        out[f"{seg.lower()}_fallback_frac"] = float(w[miss].sum())
    return out


# --------------------------------------------------------------------------------- assembler


def raw_and_channels(frames: dict[str, pl.DataFrame], vocab: dict[str, dict[str, list[str]]], *,
                     cstar_target: "float | dict[str, float] | None" = None,
                     weight: str = "log2p1", prefiltered: bool = False,
                     on_duplicate: str = "error") -> tuple[dict[str, float], dict[str, float]]:
    """Every raw feature and every channel for one sample.

    Args:
        frames: ``{locus: clonotype frame}``, raw (this sanitises).
        vocab: ``{locus: {"vus": [...], "jus": [...], "spec": [...]}}`` from the corpus artifact,
            or from :func:`gene_vocab` at fit time.
        cstar_target: Coverage level the Hill numbers are read at. A float applies to every locus,
            a dict is per-locus, and ``None`` uses **each locus's own attained coverage** -- which
            makes the diversity features observed rather than standardised, and is the only default
            that does not smuggle in a constant from somewhere the caller cannot see.
        weight: Clone-size weight, a key of :data:`WEIGHTS`.
        prefiltered: When the caller already removed non-productive rows, report
            ``nonstd_aa_frac`` as ``nan`` rather than a confident floor of zero.
        on_duplicate: Forwarded to :func:`sanitise`.

    Returns:
        ``(raw, channels)`` -- both ``{column: value}``, holes as ``nan``.
    """
    from .layout import NO_LOCUS

    raw: dict[str, float] = {}
    chan: dict[str, float] = {}
    reads: dict[str, float] = {}
    n_present = 0

    for locus in LOCI:
        # Channels are declared for all seven loci and always emitted -- "this locus is absent"
        # is a fact worth carrying even for a corpus that does not model it. Raw features exist
        # only where the corpus has a vocabulary, since that is what the rotation is indexed by.
        modelled = locus in vocab
        df0 = frames.get(locus)
        present = df0 is not None and df0.height > 0
        clean, dropped = (sanitise(df0, on_duplicate=on_duplicate) if present
                          else (pl.DataFrame(), 0.0))
        has = clean.height > 0
        n_present += int(has)
        # The weighted frame first: qc reports fallback as a share of clone WEIGHT, not of rows,
        # so it needs the same weights every other group uses.
        wf = work_frame(clean, weight) if has else clean
        chan[f"vsig:mask:{locus}:present"] = float(has)
        chan |= {f"vsig:qc:{locus}:{k}": v for k, v in
                 qc_channel(df0 if present else pl.DataFrame(), wf, locus,
                            np.nan if prefiltered else dropped).items()}
        if not has:
            chan[f"vsig:cov:{locus}:cstar"] = np.nan
            chan[f"vsig:mask:{locus}:estimable"] = 0.0
            if locus == "IGH":
                chan["vsig:mask:IGH:c_call"] = 0.0
                chan["vsig:mask:IGH:shm"] = 0.0
            if modelled:
                for c in _locus_raw_columns(locus, vocab):
                    raw[c] = np.nan
            continue

        n_reads = float(clean[COUNT].sum())
        richness = float(clean.height)
        reads[locus] = n_reads

        attained = coverage_of(clean)
        chan[f"vsig:cov:{locus}:cstar"] = attained
        if not modelled:
            chan[f"vsig:mask:{locus}:estimable"] = 0.0
            if locus == "IGH":
                chan["vsig:mask:IGH:c_call"] = 0.0
                chan["vsig:mask:IGH:shm"] = 0.0
            continue
        level = (cstar_target.get(locus) if isinstance(cstar_target, dict)
                 else cstar_target) if cstar_target is not None else attained
        div = div_group(clean, float(level)) if level is not None and np.isfinite(level) \
            else dict.fromkeys(("1D_c", "0D_c", "2D_c", "clonality", "0D_chao", "d50"), np.nan)
        chan[f"vsig:mask:{locus}:estimable"] = float(np.isfinite(div["1D_c"]))

        groups: dict[str, dict[str, float]] = {
            "div": div,
            "depth": depth_group(clean, n_reads, richness),
            "clon": clon_group(clean, n_reads, richness),
            "len": len_group(wf),
            "aa": aa_group(wf),
            "kmer": kmer_group(wf),
            "pchem": pchem_group(wf),
            "vus": usage_group(wf, V_CALL, vocab[locus]["vus"]),
            "jus": usage_group(wf, J_CALL, vocab[locus]["jus"]),
            "spec": spec_group(wf, vocab[locus]["spec"], vocab[locus]["vus"]),
        }
        if locus == "IGH":
            groups["iso"] = iso_group(wf)
            groups["shm"] = shm_group(wf)
            chan["vsig:mask:IGH:c_call"] = float(C_CALL in clean.columns
                                                 and clean[C_CALL].is_not_null().any())
            chan["vsig:mask:IGH:shm"] = float(np.isfinite(groups["shm"]["mean_v_identity"]))
        for g, vals in groups.items():
            raw |= {f"vsig:{g}:{locus}:{k}": float(v) for k, v in vals.items()}

    raw |= {f"vsig:pair:{NO_LOCUS}:{k}": v for k, v in pair_group(reads).items()}
    chan["vsig:qc:-:n_loci_present"] = float(n_present)
    chan["vsig:qc:-:winsor_frac"] = 0.0       # corpus.apply overwrites; 0.0 = nothing clamped
    return raw, chan


def _locus_raw_columns(locus: str, vocab: dict[str, dict[str, list[str]]]) -> list[str]:
    from .layout import raw_columns

    return raw_columns("vsig", locus, vocab.get(locus))


def _demo() -> None:
    """Self-check: a synthetic two-locus sample, and the holes that must stay holes."""
    from .layout import raw_columns

    rng = np.random.default_rng(0)
    aa = np.array(list(AMINO_ACIDS))

    def frame(n, v, j, c=None):
        return pl.DataFrame({
            "v_call": [v] * n, "j_call": [j] * n, "c_call": [c] * n,
            "junction_aa": ["C" + "".join(rng.choice(aa, 12)) + "F" for _ in range(n)],
            COUNT: np.ceil(rng.zipf(1.5, n).clip(1, 900)).astype(np.int64).tolist(),
        })

    vocab = {loc: gene_vocab(loc) for loc in ("TRB", "IGH")}
    vocab |= {loc: {"vus": [], "jus": [], "spec": []} for loc in LOCI if loc not in vocab}
    frames = {"TRB": frame(600, "TRBV20-1", "TRBJ2-2"),
              "IGH": frame(300, "IGHV1-2", "IGHJ4", "IGHM")}
    raw, chan = raw_and_channels(frames, vocab)

    # every declared raw column is present, exactly once, for every locus
    want = [c for loc in (*LOCI, "-") for c in raw_columns("vsig", loc, vocab.get(loc))]
    assert set(raw) == set(want), (len(raw), len(want), set(want) ^ set(raw))

    # an absent locus is a hole everywhere, and says so in the mask -- never a zero
    assert chan["vsig:mask:TRA:present"] == 0.0 and chan["vsig:mask:TRB:present"] == 1.0
    assert all(np.isnan(raw[c]) for c in raw_columns("vsig", "TRA", vocab["TRA"]))
    assert np.isnan(chan["vsig:cov:TRA:cstar"])
    assert chan["vsig:qc:-:n_loci_present"] == 2.0

    # coverage is ALWAYS emitted for a present locus, estimable or not
    assert 0.0 <= chan["vsig:cov:TRB:cstar"] <= 1.0
    assert np.isfinite(chan["vsig:cov:IGH:cstar"])

    # the default target is the sample's own coverage, so diversity is estimable by construction
    assert chan["vsig:mask:TRB:estimable"] == 1.0
    assert np.isfinite(raw["vsig:div:TRB:1D_c"])
    # and an unreachable target is a declared hole, not a number
    _, _ = raw_and_channels(frames, vocab, cstar_target=0.999999)
    hi, hc = raw_and_channels(frames, vocab, cstar_target=0.999999)
    assert np.isnan(hi["vsig:div:TRB:1D_c"]) and hc["vsig:mask:TRB:estimable"] == 0.0

    # compositions close: clr coordinates over the FULL composition sum to zero with the residual
    v = np.array([raw[c] for c in raw_columns("vsig", "TRB", vocab["TRB"]) if ":vus:" in c])
    assert np.isfinite(v).all() and v.size == len(vocab["TRB"]["vus"])

    # IGH-only groups exist at IGH and nowhere else
    assert np.isfinite(raw["vsig:iso:IGH:IgM"])
    assert "vsig:iso:TRB:IgM" not in raw
    assert np.isnan(raw["vsig:shm:IGH:mean_v_identity"])     # no v_identity column
    assert chan["vsig:mask:IGH:shm"] == 0.0
    assert chan["vsig:mask:IGH:c_call"] == 1.0

    # the cross-locus ratios see both loci and divide depth out
    assert np.isfinite(raw["vsig:pair:-:log_IGH_TRB"])

    print(f"features OK  raw={len(raw)} channels={len(chan)}")


if __name__ == "__main__":
    _demo()
