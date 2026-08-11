"""The statistics half of the signature: one named block at a time.

Each function here takes one sample's clonotype frame (already restricted to a single locus,
already sanitised) and returns ``{feature_name: value}`` matching the features that block
declares in :mod:`vdjtools.signature.layout`. Values come out **already transformed** — the
declared transform is applied here, where the denominator is still in scope, because a
proportion separated from the count it was observed on cannot be stabilised afterwards.

Nothing here fits anything. Every number is a function of this one sample plus frozen
constants, which is what lets two people who never share a cohort produce comparable vectors.
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
    strip_allele,
)
from . import transform as T

#: The 20 proteinogenic amino acids, anchored. A junction containing anything else — a stop
#: codon, an ambiguity code, the legacy out-of-frame marker — is dropped before anything is
#: computed. This is not a crash guard: the embedder accepts most malformed junctions silently
#: and returns a finite, meaningless distance, so the filter is the only thing standing between
#: a contaminated geometry block and a plausible-looking wrong answer. The fraction dropped is
#: reported as ``qc:*:nonstd_aa_frac`` rather than logged, because a collaborator needs to see it.
VALID_AA = r"^[ACDEFGHIKLMNPQRSTVWY]+$"

#: Isotype classes, and the constant-gene prefixes that map onto them. ``IGHGP`` is a
#: pseudogene and ``IGHC`` is ambiguous, so neither is called.
ISOTYPES: dict[str, tuple[str, ...]] = {
    "IgM": ("IGHM",), "IgD": ("IGHD",),
    "IgG": ("IGHG1", "IGHG2", "IGHG3", "IGHG4"),
    "IgA": ("IGHA1", "IGHA2"), "IgE": ("IGHE",),
}

#: Locus pairs whose read-count ratio is a compartment read-out: T-vs-B balance, the
#: gamma-delta share, the light-chain balance. A ratio rather than two counts, because the
#: sequencing depth that drives both cancels.
PAIRS: tuple[tuple[str, str], ...] = (
    ("TRA", "TRB"), ("TRG", "TRB"), ("TRD", "TRB"), ("IGK", "IGL"), ("IGH", "TRB"),
)


def sanitise(df: pl.DataFrame) -> tuple[pl.DataFrame, float]:
    """Drop unusable clonotypes; return the frame and the dropped **weight** fraction.

    Dropped by weight, not by row: losing one dominant clone matters more than losing fifty
    singletons, and the row fraction would hide that.
    """
    if df.height == 0:
        return df, 0.0
    total = float(df[COUNT].sum())
    keep = df.filter(
        (pl.col(COUNT) > 0)
        & pl.col(JUNCTION_AA).is_not_null()
        & pl.col(JUNCTION_AA).str.contains(VALID_AA)      # anchored: a match, not a search
    )
    kept = float(keep[COUNT].sum()) if keep.height else 0.0
    return keep, (1.0 - kept / total) if total > 0 else 0.0


def work_frame(df: pl.DataFrame, weight: str = "log2p1") -> pl.DataFrame:
    """Overwrite ``frequency`` with the normalised clone weight, so every profiler agrees.

    The signature uses one clone-weight measure throughout — ``w = log2(1+count)/Σ`` by default
    — and every vdjtools profiler is then called with ``weight="freq"``. Weighting by raw reads
    instead would let one 30,000-read clone be the entire profile; weighting every clonotype
    equally would throw away the expansion signal. The concave weight sits between.

    **Order matters.** Call this *after* filtering and never call ``filter_functional``,
    ``downsample`` or ``select_top`` afterwards — each recomputes ``frequency`` from the counts
    and would silently restore read weighting.
    """
    a = df[COUNT].to_numpy().astype(float)
    g = {"log2p1": lambda x: np.log2(1.0 + x),
         "log1p": np.log1p,
         "anscombe": lambda x: np.sqrt(x + 0.375),
         "duplicate_count": lambda x: x,
         "distinct": np.ones_like}[weight](a)
    s = g.sum()
    return df.with_columns(pl.Series(FREQ, g / s if s > 0 else g))


def _weights(df: pl.DataFrame) -> np.ndarray:
    return df[FREQ].to_numpy().astype(float)


# --------------------------------------------------------------------------------- provenance


def qc_block(raw: pl.DataFrame, clean: pl.DataFrame, locus: str,
             nonstd_frac: float) -> dict[str, float]:
    """Whether this sample's gene calls are in a vocabulary we recognise.

    An unrecognised V or J call does **not** raise anywhere downstream: the germline distance
    lookup falls back to the maximum observed distance, so a cohort using Adaptive nomenclature
    or an older IMGT release yields a fully populated, entirely plausible, systematically wrong
    geometry block. These fractions are the one number that tells a collaborator their vector is
    not comparable to ours, so they are columns rather than a warning nobody reads.
    """
    from ..model.reference import load_germline

    out = {"nonstd_aa_frac": T.logit(nonstd_frac, max(raw.height, 1))}
    try:
        germ = load_germline(locus)
    except Exception:                       # a locus with no bundled germline: unmeasurable
        return {**out, "v_fallback_frac": np.nan, "j_fallback_frac": np.nan}

    known = {seg: set(germ.filter(pl.col("segment") == seg)["gene"].to_list())
             for seg in ("V", "J")}
    w = _weights(clean) if clean.height else np.zeros(0)
    for seg, col in (("V", V_CALL), ("J", J_CALL)):
        if not clean.height:
            out[f"{seg.lower()}_fallback_frac"] = np.nan
            continue
        genes = clean.select(strip_allele(pl.col(col)).alias("g"))["g"].to_list()
        miss = np.array([g is None or g not in known[seg] for g in genes])
        out[f"{seg.lower()}_fallback_frac"] = T.logit(float(w[miss].sum()), clean.height)
    return out


# ------------------------------------------------------------------------------------- depth


def depth_block(df: pl.DataFrame, stats: dict) -> dict[str, float]:
    """How much was seen, and how much was not.

    ``S_unseen`` is Chao's estimate of the clonotypes that exist but were never drawn. It is
    carried explicitly rather than folded into a diversity estimate because it is the honest
    statement of what the sample could not resolve.
    """
    a = df[COUNT].to_numpy()
    f1 = float((a == 1).sum())
    f2 = float((a == 2).sum())
    s_unseen = f1 * f1 / (2.0 * f2) if f2 > 0 else f1 * (f1 - 1.0) / 2.0
    return {"reads": T.log10(stats["n_reads"]),
            "richness": T.log10(stats["richness"]),
            "S_unseen": T.log1p(max(s_unseen, 0.0))}


# --------------------------------------------------------------------------------- diversity


def div_block(df: pl.DataFrame, cstar: float, tier_full: bool = False) -> dict[str, float]:
    """Hill numbers standardised to a **frozen** coverage level, plus shape statistics.

    Coverage standardisation is what makes a hundred clonotypes and a hundred thousand
    comparable: both are evaluated at the same completeness rather than at their own depth. It
    only works when the sample actually reaches ``cstar`` — beyond that the estimator
    extrapolates, and extrapolation is not a mild approximation here. Measured on one real
    repertoire subsampled across a 200x depth range: at a coverage level every subsample
    attained, the Shannon diversity agreed to within 4%; at a level the shallow ones had to
    extrapolate to, the same statistic inflated roughly tenfold.

    So ``cstar`` is frozen from the reference corpus's own attained-coverage distribution — a
    low quantile of what samples actually reach, not a textbook 0.95 — and whether *this* sample
    reached it is recorded by ``mask:*:estimable``. When it did not, the columns are ``nan``:
    a hole a downstream model can see, rather than a confident wrong number.

    Returns:
        Transformed features, or ``nan`` for every one when the sample cannot support them.
    """
    from ..stats.inext import estimate_d

    nan = {k: np.nan for k in ("1D_c", "0D_c", "2D_c", "clonality")}
    if tier_full:
        nan |= {"0D_chao": np.nan, "d50": np.nan}
    counts = df[COUNT].to_numpy().astype(np.int64)
    if counts.size < 2 or counts.sum() < 2:
        return nan

    try:
        est = estimate_d(counts, base="coverage", level=cstar, q=(0, 1, 2), se=False)
    except (ValueError, ZeroDivisionError):
        return nan
    if not estimable(df, cstar, est):
        return nan

    qd = {int(r["order_q"]): float(r["qD"]) for r in est.to_dicts()}
    # Clonality from the *standardised* Hill numbers, not from observed evenness. Observed
    # evenness is a ratio of two depth-dependent quantities and inherits their depth dependence
    # in full: measured on one repertoire across a 667x depth range it drifted by a factor of
    # nearly two, while the same quantity built from the coverage-standardised numbers is flat.
    # 1 - ln(1D)/ln(0D) is Pielou's evenness with both terms evaluated at the same completeness.
    evenness = np.log(qd[1]) / np.log(qd[0]) if qd[0] > 1 else 0.0
    out = {"1D_c": T.log10(qd[1]), "0D_c": T.log10(qd[0]), "2D_c": T.log10(qd[2]),
           "clonality": T.logit(np.clip(1.0 - evenness, 0.0, 1.0), counts.size)}
    if tier_full:
        f1, f2 = float((counts == 1).sum()), float((counts == 2).sum())
        chao = counts.size + (f1 * f1 / (2.0 * f2) if f2 > 0 else f1 * (f1 - 1.0) / 2.0)
        srt = np.sort(counts)[::-1]
        k = int(np.searchsorted(np.cumsum(srt), 0.5 * counts.sum()) + 1)
        out |= {"0D_chao": T.log10(chao), "d50": T.logit(k / counts.size, counts.size)}
    return out


def estimable(df: pl.DataFrame, cstar: float, est: "pl.DataFrame | None" = None) -> bool:
    """Whether a coverage-standardised diversity is a measurement rather than an extrapolation.

    Two conditions, both necessary: the estimator must not have had to extrapolate, and the
    target depth must be within twice the observed one. The second catches the case where the
    estimator interpolates nominally but is leaning on almost no data.
    """
    from ..stats.inext import estimate_d

    counts = df[COUNT].to_numpy().astype(np.int64)
    if counts.size < 2 or counts.sum() < 2:
        return False
    if est is None:
        try:
            est = estimate_d(counts, base="coverage", level=cstar, q=(1,), se=False)
        except (ValueError, ZeroDivisionError):
            return False
    rows = est.to_dicts()
    return (all(r["method"] != "extrapolation" for r in rows)
            and max(float(r["m"]) for r in rows) <= 2.0 * counts.sum())


def clon_block(stats: dict) -> dict[str, float]:
    """The shape of the clone-size distribution as a composition, plus the top clone's share.

    ``f1``/``f2``/``f3plus`` — how many clonotypes were seen once, twice, more — is a
    three-part composition, so it is read in log-ratio coordinates: only two of the three are
    shipped because the third is exactly determined by them.
    """
    n = max(stats["richness"], 1.0)
    c = T.clr({"f1": stats["f1"], "f2": stats["f2"], "f3plus": stats["f3plus"]}, m=n)
    return {"f1": c["f1"], "f2": c["f2"],
            "top": T.logit(stats["top_clone_fraction"], stats["n_reads"])}


# ---------------------------------------------------------------------------------- geometry-free


def len_block(df: pl.DataFrame, tier_standard: bool = True) -> dict[str, float]:
    """Weighted moments of the junction length distribution.

    Moments rather than the histogram the legacy spectratype emits: at a hundred clonotypes a
    17-bin histogram is mostly zeros, and its bins are not independent anyway.
    """
    if df.height == 0:
        return {"mean": np.nan, "sd": np.nan, **({"skew": np.nan} if tier_standard else {})}
    ln = df[JUNCTION_AA].str.len_chars().to_numpy().astype(float)
    w = _weights(df)
    mean = float((w * ln).sum())
    var = float((w * (ln - mean) ** 2).sum())
    sd = np.sqrt(var)
    out = {"mean": mean, "sd": sd}
    if tier_standard:
        out["skew"] = float((w * (ln - mean) ** 3).sum() / sd ** 3) if sd > 0 else 0.0
    return out


def iso_block(df: pl.DataFrame, tier_full: bool = False) -> dict[str, float]:
    """Isotype composition of an IGH repertoire, in log-ratio coordinates.

    The uncalled share is a real part of the composition, not a rounding error — roughly two
    fifths of IGH reads carry no constant-gene call — so it closes the composition here rather
    than being silently dropped the way a usage profile would drop it. Coordinates are taken
    over the whole composition and then selected, so a tier shipping four of six parts still
    divides by the six-part geometric mean.
    """
    keys = list(ISOTYPES) if tier_full else [k for k in ISOTYPES if k != "IgE"]
    if df.height == 0:
        return dict.fromkeys(keys, np.nan)
    w = _weights(df)
    calls = df[C_CALL].to_list() if C_CALL in df.columns else [None] * df.height
    parts = dict.fromkeys(ISOTYPES, 0.0)
    for c, wi in zip(calls, w):
        for iso, prefixes in ISOTYPES.items():
            if c in prefixes:
                parts[iso] += float(wi)
                break
    parts["_uncalled"] = max(1.0 - sum(parts.values()), 0.0)
    coords = T.clr(parts, m=df.height)
    return {k: coords[k] for k in keys}


def shm_block(df: pl.DataFrame) -> dict[str, float]:
    """Mean somatic hypermutation load, as ``1 − v_identity``.

    Absent from the canonical schema — vdjtools' readers narrow to eight columns and
    ``v_identity`` is not one of them — so this masks out unless the caller kept it explicitly.
    """
    if df.height == 0 or "v_identity" not in df.columns:
        return {"mean_v_identity": np.nan}
    v = df["v_identity"].to_numpy().astype(float)
    ok = np.isfinite(v)
    if not ok.any():
        return {"mean_v_identity": np.nan}
    w = _weights(df)[ok]
    mean_identity = float((w * v[ok]).sum() / w.sum()) if w.sum() > 0 else np.nan
    return {"mean_v_identity": T.logit(np.clip(mean_identity, 0.0, 1.0), int(ok.sum()))}


def pair_block(reads: dict[str, float]) -> dict[str, float]:
    """Log read-count ratios between loci — compartment balance, with depth divided out."""
    return {f"log_{a}_{b}": float(np.log10((reads.get(a, 0.0) + 1.0) / (reads.get(b, 0.0) + 1.0)))
            for a, b in PAIRS}


# ------------------------------------------------------------------------- composition (full)


def aa_block(df: pl.DataFrame) -> dict[str, float]:
    """Weighted residue composition of the junctions — 20 parts, arcsine-stabilised.

    Only ``k=1``. A 3-mer spectrum is 8,000 cells against roughly 1,200 tokens at the corpus
    median depth, so it is >85% structural zeros and its log-ratio coordinates end up a
    depth read-out wearing a motif label. Single residues are ~1,400 observations over 20 parts
    — genuinely estimated. Motif structure is the geometry half's job, where it is dense.
    """
    from ..features.kmer import kmer_profile

    keys = list("ACDEFGHIKLMNPQRSTVWY")
    if df.height == 0:
        return dict.fromkeys(keys, np.nan)
    prof = kmer_profile(df, k=1, weight="freq", by_locus=False)
    got = dict(zip(prof["kmer"].to_list(), prof["weight"].to_list()))
    total = sum(got.values())
    m = float(df[JUNCTION_AA].str.len_chars().sum())      # residues actually observed
    if total <= 0:
        return dict.fromkeys(keys, np.nan)
    return {k: float(T.arcsine(got.get(k, 0.0) / total, m)) for k in keys}


def pchem_block(df: pl.DataFrame, regions=("all", "center")) -> dict[str, float]:
    """Weighted mean physicochemistry of the junction, over two regions.

    Two regions rather than the legacy five, and means rather than four quantiles apiece: at a
    hundred clonotypes the quantiles of a per-clonotype property are noise, while the weighted
    mean is a well-behaved average over every residue seen.
    """
    from ..features.physchem import DEFAULT_PROPERTIES, physchem_profile

    keys = [f"{r}_{p}" for r in regions for p in DEFAULT_PROPERTIES]
    if df.height == 0:
        return dict.fromkeys(keys, np.nan)
    out = dict.fromkeys(keys, np.nan)
    for region in regions:
        prof = physchem_profile(df, group_by="locus", region=region, weight="freq")
        for row in prof.iter_rows(named=True):
            key = f"{region}_{row['property']}"
            if key in out and row["mean_value"] is not None:
                out[key] = float(row["mean_value"])
    return out


def pgen_block(df: pl.DataFrame, locus: str, *, q05: float | None = None,
               n_max: int = 2000, threads: int = 0) -> dict[str, float]:
    """Generation probability of the junctions under the bundled recombination model.

    How *surprising* this repertoire's receptors are: a clone that recombination produces
    readily is weak evidence of anything, while a low-Pgen clone that reached detectable size
    had help. ``frac_atypical`` is the share sitting below a frozen reference quantile.

    V and J are **marginalised** (``v=None, j=None``) rather than conditioned on the observed
    calls. The batch API wants allele-level names and raises on gene-level ones, and AIRR frames
    are gene-level after allele stripping — so conditioning here would either crash or, worse,
    silently condition on whichever alleles happened to survive.

    Subsampled deterministically to ``n_max`` junctions, seeded from the locus name via CRC32
    rather than :func:`hash`, whose string hashing is randomised per process and would make the
    column irreproducible across runs.
    """
    import zlib

    from ..model.bundled import load_bundled
    from ..model.native import pgen_aa_batch

    nan = {"mean_log10": np.nan, "sd_log10": np.nan, "frac_atypical": np.nan}
    if df.height == 0:
        return nan
    juncs = df[JUNCTION_AA].to_list()
    if len(juncs) > n_max:
        rng = np.random.default_rng(zlib.crc32(locus.encode()))
        juncs = [juncs[i] for i in rng.choice(len(juncs), n_max, replace=False)]
    try:
        p = np.asarray(pgen_aa_batch(load_bundled(locus), juncs, v=None, j=None, threads=threads),
                       dtype=float)
    except Exception:                       # no bundled model for this locus/species
        return nan
    p = p[np.isfinite(p) & (p > 0)]
    if p.size == 0:
        return nan
    lp = np.log10(p)
    out = {"mean_log10": float(lp.mean()), "sd_log10": float(lp.std())}
    out["frac_atypical"] = (float(T.logit(float((lp < q05).mean()), p.size))
                            if q05 is not None else np.nan)
    return out


def _demo() -> None:
    """Self-check on a synthetic sample: shapes, holes, and the depth-honesty property."""
    rng = np.random.default_rng(0)
    n = 200
    df = pl.DataFrame({
        V_CALL: ["TRBV20-1"] * (n - 5) + ["TRBV999"] * 5,
        J_CALL: ["TRBJ2-2"] * n,
        C_CALL: [None] * n,
        JUNCTION_AA: ["CASS" + "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 8)) + "F"
                      for _ in range(n)],
        COUNT: rng.integers(1, 50, n).tolist(),
        FREQ: [0.0] * n,
    })
    clean, dropped = sanitise(df)
    assert clean.height == n and dropped == 0.0
    work = work_frame(clean)
    assert abs(work[FREQ].sum() - 1.0) < 1e-12, "weights do not close"

    stats = {"n_reads": float(clean[COUNT].sum()), "richness": float(clean.height),
             "f1": float((clean[COUNT].to_numpy() == 1).sum()),
             "f2": float((clean[COUNT].to_numpy() == 2).sum()),
             "f3plus": float((clean[COUNT].to_numpy() >= 3).sum()),
             "top_clone_fraction": float(clean[COUNT].max() / clean[COUNT].sum())}

    d = depth_block(work, stats)
    assert d["reads"] > d["richness"] > 0

    q = qc_block(df, work, "TRB", dropped)
    assert q["v_fallback_frac"] > T.logit(0.0, n), "unknown V gene went unreported"

    c = clon_block(stats)
    assert np.isfinite(list(c.values())).all()

    ln = len_block(work)
    assert 10 < ln["mean"] < 16 and ln["sd"] >= 0

    iso = iso_block(work)
    assert set(iso) == {"IgM", "IgD", "IgG", "IgA"} and np.isfinite(list(iso.values())).all()

    # an unreachable coverage level must produce holes, never a confident wrong number
    assert all(np.isnan(v) for v in div_block(work, cstar=0.999).values())

    p = pair_block({"TRA": 100.0, "TRB": 100.0})
    assert abs(p["log_TRA_TRB"]) < 1e-9

    assert all(np.isnan(v) for v in shm_block(work).values()), "shm should mask out"
    print("blocks OK")


if __name__ == "__main__":
    _demo()
