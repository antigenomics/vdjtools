"""The corpus: build one, winsorize it, fit the rotation and the scaling, apply it to a sample.

This is the only module here that needs more than one sample, and everything it fits comes out of
**one pass over one matrix of repertoires**. That is the whole design, and it is a correction: the
system this replaces fitted its rotation on 10,000 individual clonotypes from a prototype panel
while every column that rotation produced was a repertoire statistic, and its centre and scale
arrived from a separate corpus of real samples. Two independent fits, stitched -- which is how one
shipped reference came to pair a centre of exactly ``0.0`` with a scale plainly fitted from data,
putting a corpus-typical sample 81 robust deviations out.

So, per locus ``L``::

    M_L (N, p_L)  one row per repertoire, columns in layout order
      -> bounds       per column, BY PERCENTILE, side from the declared support
      -> centre/scale median and 1.4826*MAD, from OBSERVED entries only
      -> rotation     R_L (p_L, k_L)
      -> pc centre/scale/quantiles   from the corpus's own PC scores

Nothing per-sample enters the artifact. The coverage a sample attained is a per-sample statistic
and ships as a channel; the coverage *target* is a runtime argument. Carrying a per-sample quantity
in a corpus artifact is how one shipped reference ended up standardising tissue to a level measured
on blood, digit for digit across all seven loci.

**Winsorization is by percentile and one-sided wherever the metric is bounded on that side.** A
robust-SD bound is not a fixed amount of probability: it depends on the distribution's shape, so an
8-SD bound trims nothing at all on a right-skewed count and a third of the corpus on a
near-degenerate column. And the side comes from the declared support, never from the transform --
``transform="none"`` spans three different supports.
"""
from __future__ import annotations

import hashlib
import json
import platform
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import polars as pl

from . import layout as L

#: ``1.4826 * MAD`` estimates the standard deviation of a normal distribution.
MAD_TO_SD = 1.4826

#: Winsorization percentiles the artifact stores, so a caller picks one at apply time without a
#: refit. ``0.01`` trims the 1st/99th percentile, ``0.05`` the 5th/95th.
WINSOR_PS: tuple[float, ...] = (0.01, 0.05)

#: A column whose corpus spread is at or below this is treated as never having varied: it is
#: centred, not scaled, and contributes nothing to the rotation. Not a tolerance on a measurement
#: -- a structural zero of a composition really is constant, and dividing by its MAD would turn
#: float noise into the corpus's dominant direction.
MIN_SPREAD = 1e-12

#: Default components per locus. A count rather than a variance fraction, because the artifact
#: stores the rotation only up to the count it was fitted at and a fraction on a wide sparse group
#: resolves to enough components to make the artifact hundreds of megabytes. The full eigenvalue
#: spectrum is stored either way, so the variance a given count reaches is always readable.
DEFAULT_COMPONENTS: int = 128

#: Modes for :func:`apply`. ``features`` clamps raw features to the corpus bounds before rotating;
#: ``pcs`` rotates first and clamps the PC scores; ``none`` clamps nothing.
MODES: tuple[str, ...] = ("features", "pcs", "none")


# ------------------------------------------------------------------------------ winsorization


def bounds_for(x: np.ndarray, supports: list[str], p: float) -> tuple[np.ndarray, np.ndarray]:
    """Per-column winsorization bounds at percentile ``p``, side chosen by declared support.

    Args:
        x: ``(N, p_L)`` corpus matrix; non-finite entries are ignored.
        supports: One support name per column, from
            :func:`~vdjtools.signature.layout.support_of`.
        p: Tail fraction, e.g. ``0.01`` for the 1st/99th percentile.

    Returns:
        ``(lo, hi)``, each ``(p_L,)``. A side that is not trimmed is ``-inf`` / ``+inf``, so the
        clamp is a no-op there rather than a bound that happens not to bite.
    """
    if not 0.0 < p < 0.5:
        raise ValueError(f"winsor p must be in (0, 0.5); got {p}")
    sides = np.array([L.SUPPORTS[s] for s in supports])          # (p_L, 2) bottom, top
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)     # all-nan column: see robust_loc_scale
        finite = np.where(np.isfinite(x), x, np.nan)
        qlo = np.nanquantile(finite, p, axis=0)
        qhi = np.nanquantile(finite, 1.0 - p, axis=0)
    lo = np.where(sides[:, 0], qlo, -np.inf)
    hi = np.where(sides[:, 1], qhi, np.inf)
    # An all-nan column yields nan quantiles, which would clamp everything to nan. A column the
    # corpus never observed has no bound to offer, so it offers none.
    lo = np.where(np.isfinite(lo), lo, -np.inf)
    hi = np.where(np.isfinite(hi), hi, np.inf)
    return lo, hi


def clamp(x: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Clamp to ``[lo, hi]``, and report which entries moved.

    The second return value is what keeps the bound from being silent. Clamping is many-to-one: it
    overwrites a real measurement of a real sample, and the sample cannot tell afterwards. A
    reference whose centre was wrong by 57 robust deviations went unnoticed for months precisely
    because every affected sample was clamped to the bound and looked like an ordinary tail case.
    """
    out = np.clip(x, lo, hi)
    moved = np.isfinite(x) & (out != x)
    return out, moved


def robust_loc_scale(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-column ``(median, 1.4826*MAD)`` from observed entries only.

    **Observed entries only, before any imputation.** Filling holes first and measuring afterwards
    deflates the scale in proportion to how sparse a column is, so the least-observed locus ends up
    with the largest apparent values and dominates every distance and every principal component.
    """
    # An all-nan column is an expected state, not an anomaly: a feature this corpus never observed
    # -- ``shm`` on synthetic data with no ``v_identity``, a locus a real cohort lacks. It comes out
    # loc=0, scale=0 and is then excluded from the rotation by MIN_SPREAD, which is the right
    # outcome, so numpy's warning about it is noise rather than information.
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        finite = np.where(np.isfinite(x), x, np.nan)
        loc = np.nanmedian(finite, axis=0)
        mad = np.nanmedian(np.abs(finite - loc), axis=0)
    return np.nan_to_num(loc), np.nan_to_num(mad) * MAD_TO_SD


def _standardize(x: np.ndarray, loc: np.ndarray, scale: np.ndarray) -> np.ndarray:
    """``(x - loc)/scale`` with holes and never-varying columns sent to the centre.

    A hole becomes 0, i.e. the corpus centre -- the least-committal imputation available once the
    location has already been estimated without it. A column whose corpus spread is below
    :data:`MIN_SPREAD` also becomes 0: it carries no information and dividing by its MAD would
    promote float noise to a principal direction.
    """
    safe = np.where(scale > MIN_SPREAD, scale, 1.0)
    z = (x - loc) / safe
    z = np.where(scale > MIN_SPREAD, z, 0.0)
    return np.where(np.isfinite(z), z, 0.0)


# -------------------------------------------------------------------------------- the artifact


@dataclass
class LocusFit:
    """Everything fitted for one locus. Arrays are column-aligned to :attr:`columns`."""

    columns: list[str]                       # raw feature names, rotation-input order
    supports: list[str]
    loc: np.ndarray                          # (p,) median
    scale: np.ndarray                        # (p,) 1.4826*MAD
    rotation: np.ndarray                     # (p, k)
    eigenvalues: np.ndarray                  # (min(N, p),) full spectrum, descending
    pc_loc: np.ndarray                       # (k,)
    pc_scale: np.ndarray                     # (k,)
    bounds: dict[float, tuple[np.ndarray, np.ndarray]]      # p -> raw (lo, hi)
    pc_bounds: dict[float, tuple[np.ndarray, np.ndarray]]   # p -> pc  (lo, hi)
    n_obs: np.ndarray                        # (p,) rows that observed each column
    n_rows: int                              # rows this locus was fitted on

    @property
    def k(self) -> int:
        return int(self.rotation.shape[1])

    def variance_at(self, k: int) -> float:
        """Cumulative variance fraction the first ``k`` components reach."""
        ev = self.eigenvalues
        tot = float(ev.sum())
        return float(ev[:k].sum() / tot) if tot > 0 else float("nan")


@dataclass
class Corpus:
    """A fitted corpus: one :class:`LocusFit` per locus, plus the manifest.

    A load **verifies** rather than trusts: the raw column order is re-derived from the stored
    vocabulary through the current layout and compared, so a layout change that would silently
    re-index the rotation raises instead. The bundled-model version is checked the same way,
    because the synthetic corpora are drawn from those models and retraining them moves the result
    by 0.1-1.6% per locus.
    """

    sig: str
    name: str
    vocab: dict[str, dict[str, list[str]]]
    fits: dict[str, LocusFit]
    meta: dict = field(default_factory=dict)

    @property
    def k(self) -> dict[str, int]:
        """``{locus: components}``. A locus absent here emits no columns -- how a corpus that
        could not fit a locus declares it, rather than shipping a rotation of zeros."""
        return {loc: f.k for loc, f in self.fits.items()}

    def columns(self, n_components: "int | float | None" = None) -> list[str]:
        """The signature columns this corpus emits: rotated columns then channels."""
        return L.signature_columns(self.sig, self.resolve_k(n_components))

    def resolve_k(self, n_components: "int | float | None" = None) -> dict[str, int]:
        """``{locus: k}`` after apply-time truncation.

        Args:
            n_components: ``None`` keeps what was fitted. An ``int`` truncates every locus to that
                count. A ``float`` in ``(0, 1)`` truncates each locus to the fewest components
                reaching that cumulative variance.

        Raises:
            ValueError: If an ``int`` exceeds what a locus was fitted with, or a ``float`` target
                is not reached by the stored spectrum. Both refuse rather than pad: the artifact
                stores the rotation only up to the count it was fitted at, so the components asked
                for genuinely do not exist and returning fewer under the requested name would make
                two matrices with the same column names carry different things.
        """
        if n_components is None:
            return self.k
        out: dict[str, int] = {}
        for loc, f in self.fits.items():
            if isinstance(n_components, int) and not isinstance(n_components, bool):
                if n_components > f.k:
                    raise ValueError(
                        f"{self.name}:{loc} was fitted with {f.k} components; {n_components} were "
                        f"asked for. Increasing the count needs a refit -- the rotation is stored "
                        f"only up to what it was fitted at. The stored spectrum says the first "
                        f"{f.k} reach {f.variance_at(f.k):.4f} of the variance.")
                out[loc] = int(n_components)
            else:
                frac = float(n_components)
                if not 0.0 < frac < 1.0:
                    raise ValueError(f"n_components as a fraction must be in (0, 1); got {frac}")
                need = _k_for_variance(f.eigenvalues, frac, cap=f.k)
                if need is None:
                    raise ValueError(
                        f"{self.name}:{loc} reaches only {f.variance_at(f.k):.4f} of the variance "
                        f"in its {f.k} stored components; {frac} was asked for. Refit with a "
                        f"larger n_components.")
                out[loc] = need
        return out

    # ------------------------------------------------------------------ persistence

    def save(self, path: "str | Path") -> Path:
        """Write ``<path>.npz`` (arrays) and ``<path>.json`` (manifest and vocabulary)."""
        path = Path(path).with_suffix(".npz")
        arrays: dict[str, np.ndarray] = {}
        for loc, f in self.fits.items():
            arrays[f"{loc}/columns"] = np.array(f.columns)
            arrays[f"{loc}/loc"] = f.loc
            arrays[f"{loc}/scale"] = f.scale
            arrays[f"{loc}/rotation"] = f.rotation.astype(np.float32)
            arrays[f"{loc}/eigenvalues"] = f.eigenvalues
            arrays[f"{loc}/pc_loc"] = f.pc_loc
            arrays[f"{loc}/pc_scale"] = f.pc_scale
            arrays[f"{loc}/n_obs"] = f.n_obs
            for p, (lo, hi) in f.bounds.items():
                arrays[f"{loc}/bound_lo_{p}"] = lo
                arrays[f"{loc}/bound_hi_{p}"] = hi
            for p, (lo, hi) in f.pc_bounds.items():
                arrays[f"{loc}/pcbound_lo_{p}"] = lo
                arrays[f"{loc}/pcbound_hi_{p}"] = hi
        # savez_compressed, then hash the BYTES we wrote. Hashing the arrays would not detect a
        # writer change, and the point of the hash is that a file in someone else's environment
        # can be proved to be the one fitted here.
        np.savez_compressed(path, **arrays)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        meta = {**self.meta, "sig": self.sig, "name": self.name,
                "n_rows": {loc: f.n_rows for loc, f in self.fits.items()},
                "k": self.k,
                "variance_at_k": {loc: f.variance_at(f.k) for loc, f in self.fits.items()},
                "winsor_ps": list(WINSOR_PS), "content_sha256": digest,
                "vocab": self.vocab}
        path.with_suffix(".json").write_text(json.dumps(meta, indent=1, sort_keys=True))
        return path

    @classmethod
    def load(cls, path: "str | Path", *, verify: bool = True) -> "Corpus":
        """Read an artifact, and check it is the one the caller thinks it is.

        Args:
            path: The ``.npz`` path; the ``.json`` sidecar must sit beside it.
            verify: Re-derive the raw column order from the stored vocabulary through the current
                layout and compare, and check the recorded bundled-model version against the
                installed one. Turn it off only to inspect a deliberately stale artifact.

        Raises:
            ValueError: On a column-order or model-version mismatch. Both would otherwise produce
                numbers in a different coordinate system that look perfectly reasonable.
        """
        path = Path(path).with_suffix(".npz")
        meta = json.loads(path.with_suffix(".json").read_text())
        z = np.load(path, allow_pickle=False)
        sig, name, vocab = meta["sig"], meta["name"], meta["vocab"]
        loci = sorted({k.split("/")[0] for k in z.files}, key=_locus_order)
        fits: dict[str, LocusFit] = {}
        for loc in loci:
            cols = [str(c) for c in z[f"{loc}/columns"]]
            fits[loc] = LocusFit(
                columns=cols,
                supports=[L.support_of(c) for c in cols],
                loc=z[f"{loc}/loc"], scale=z[f"{loc}/scale"],
                rotation=z[f"{loc}/rotation"].astype(np.float64),
                eigenvalues=z[f"{loc}/eigenvalues"],
                pc_loc=z[f"{loc}/pc_loc"], pc_scale=z[f"{loc}/pc_scale"],
                bounds={p: (z[f"{loc}/bound_lo_{p}"], z[f"{loc}/bound_hi_{p}"])
                        for p in meta["winsor_ps"] if f"{loc}/bound_lo_{p}" in z.files},
                pc_bounds={p: (z[f"{loc}/pcbound_lo_{p}"], z[f"{loc}/pcbound_hi_{p}"])
                           for p in meta["winsor_ps"] if f"{loc}/pcbound_lo_{p}" in z.files},
                n_obs=z[f"{loc}/n_obs"], n_rows=int(meta["n_rows"][loc]))
        out = cls(sig=sig, name=name, vocab=vocab, fits=fits, meta=meta)
        if verify:
            out.verify()
        return out

    def verify(self) -> None:
        """Raise unless the stored column order and model version match this installation."""
        for loc, f in self.fits.items():
            want = L.raw_columns(self.sig, loc, self.vocab.get(loc))
            if want != f.columns:
                extra, missing = set(f.columns) - set(want), set(want) - set(f.columns)
                raise ValueError(
                    f"corpus {self.name!r} locus {loc}: the raw column order it was fitted with no "
                    f"longer matches this layout ({len(f.columns)} stored, {len(want)} derived; "
                    f"{len(extra)} only in the artifact, {len(missing)} only in the layout). The "
                    f"rotation is indexed by position, so applying it now would mix features. "
                    f"Refit the corpus against this version.")
        recorded = self.meta.get("model_version")
        if recorded is not None:
            from .. import __version__ as installed
            if recorded != installed:
                raise ValueError(
                    f"corpus {self.name!r} was fitted against vdjtools {recorded} and this is "
                    f"{installed}. The synthetic corpora are drawn from that release's bundled "
                    f"recombination models, and retraining those moves the result by 0.1-1.6% per "
                    f"locus, so the mismatch is not cosmetic. Rebuild the corpus, or load with "
                    f"verify=False if you have established the models did not change.")


def _locus_order(loc: str) -> int:
    order = [*L.LOCI, L.NO_LOCUS]
    return order.index(loc) if loc in order else len(order)


def _k_for_variance(ev: np.ndarray, frac: float, cap: int) -> "int | None":
    tot = float(ev.sum())
    if tot <= 0:
        return None
    cum = np.cumsum(ev) / tot
    hit = np.nonzero(cum >= frac)[0]
    if hit.size == 0 or int(hit[0]) + 1 > cap:
        return None
    return int(hit[0]) + 1


# --------------------------------------------------------------------------------- fitting


def fit_locus(x: np.ndarray, columns: list[str], *, mode: str = "features",
              n_components: "int | float" = DEFAULT_COMPONENTS,
              winsor_p: float = 0.01) -> "LocusFit | None":
    """Fit bounds, centre, scale, rotation and PC scaling for one locus, in one pass.

    Args:
        x: ``(N, p)`` raw feature matrix. Rows that observed nothing at this locus must already be
            dropped; individual holes are fine and are handled without imputing before measuring.
        columns: Raw column names, defining the rotation's row order.
        mode: ``"features"`` winsorizes the corpus and then takes a plain SVD -- the covariance the
            SVD sees has no tails left to be dragged by, so the robustness sits in a step you can
            inspect, and the result is deterministic and bit-reproducible. ``"pcs"`` leaves the
            corpus untrimmed and uses a robust scaler with sklearn's PCA, because there the tails
            are still present at fit time and the estimator has to resist them itself.
        n_components: Count, or a variance fraction in ``(0, 1)``.
        winsor_p: The percentile whose bounds are used when ``mode="features"``. All of
            :data:`WINSOR_PS` are stored regardless.

    Returns:
        The fit, or ``None`` when the locus cannot support one (fewer rows than 2, or no column
        that ever varied). ``None`` is how a corpus declares a locus it could not fit; it is not an
        error, and it emits no columns.
    """
    if mode not in MODES:
        raise ValueError(f"unknown mode {mode!r}; known: {MODES}")
    n, p = x.shape
    if n < 2 or p == 0:
        return None
    supports = [L.support_of(c) for c in columns]
    n_obs = np.isfinite(x).sum(axis=0).astype(np.int64)

    bounds = {q: bounds_for(x, supports, q) for q in WINSOR_PS}
    work = x
    if mode == "features":
        work, _ = clamp(x, *bounds[winsor_p])

    loc, scale = robust_loc_scale(work)
    if not (scale > MIN_SPREAD).any():
        return None
    z = _standardize(work, loc, scale)

    if mode == "pcs":
        # The corpus still has its tails here, so the scaler has to resist them: median/IQR, then
        # a full-solver PCA. The randomized solver is the default above 500 features and is not
        # bit-reproducible across BLAS builds, which would make the artifact unverifiable.
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import RobustScaler
        z = RobustScaler().fit_transform(z)
        pca = PCA(n_components=min(n - 1, p), svd_solver="full").fit(z)
        ev_all = pca.explained_variance_
        comps = pca.components_.T                       # (p, min(n-1, p))
    else:
        u, s, vt = np.linalg.svd(z - z.mean(0), full_matrices=False)
        ev_all = (s ** 2) / max(n - 1, 1)
        comps = vt.T

    k = _resolve_k(ev_all, n_components, p, n)
    if k is None:
        return None
    comps = _fix_signs(comps[:, :k])
    scores = z @ comps
    pc_loc, pc_scale = robust_loc_scale(scores)
    pc_supports = ["real"] * k
    pc_bounds = {q: bounds_for((scores - pc_loc) / np.where(pc_scale > MIN_SPREAD, pc_scale, 1.0),
                               pc_supports, q) for q in WINSOR_PS}
    return LocusFit(columns=list(columns), supports=supports, loc=loc, scale=scale,
                    rotation=comps, eigenvalues=ev_all, pc_loc=pc_loc, pc_scale=pc_scale,
                    bounds=bounds, pc_bounds=pc_bounds, n_obs=n_obs, n_rows=int(n))


def _resolve_k(ev: np.ndarray, n_components: "int | float", p: int, n: int) -> "int | None":
    avail = int(min(ev.size, p, max(n - 1, 1)))
    if avail < 1:
        return None
    if isinstance(n_components, int) and not isinstance(n_components, bool):
        # At FIT time a count is a cap, not a demand. A block narrower than the request -- the
        # five cross-locus ratios against a request of 128 -- gets its full rotation, which loses
        # nothing and is recorded per locus in the artifact. The case that must refuse is the other
        # one: asking APPLY for more components than were fitted, where they genuinely do not exist
        # (see Corpus.resolve_k).
        return int(min(n_components, avail))
    frac = float(n_components)
    if not 0.0 < frac < 1.0:
        raise ValueError(f"n_components as a fraction must be in (0, 1); got {frac}")
    got = _k_for_variance(ev, frac, cap=avail)
    return avail if got is None else got


def _fix_signs(comps: np.ndarray) -> np.ndarray:
    """Flip each component so its largest-magnitude coordinate is positive.

    Eigenvector signs are arbitrary and not reproducible across BLAS builds, so they have to be
    pinned or the artifact is not comparable to itself. Largest-magnitude rather than a
    coordinate-sum rule: a sum flips whenever it happens to sit near zero, which is exactly the
    non-determinism being removed.
    """
    lead = np.argmax(np.abs(comps), axis=0)
    sign = np.sign(comps[lead, np.arange(comps.shape[1])])
    sign[sign == 0] = 1.0
    return comps * sign


def fit(rows: list[dict[str, float]], vocab: dict[str, dict[str, list[str]]], *, sig: str,
        name: str, mode: str = "features", n_components: "int | float" = DEFAULT_COMPONENTS,
        winsor_p: float = 0.01, meta: "dict | None" = None) -> Corpus:
    """Fit a corpus from raw feature rows -- one dict per repertoire.

    Each locus is fitted only on the rows that observed it, so a corpus where half the samples have
    no TRD does not learn TRD from imputed values.
    """
    loci = [loc for loc in (*L.LOCI, L.NO_LOCUS) if loc in vocab or loc == L.NO_LOCUS]
    fits: dict[str, LocusFit] = {}
    trimmed: dict[str, float] = {}
    for locus in loci:
        cols = L.raw_columns(sig, locus, vocab.get(locus))
        if not cols:
            continue
        x = np.array([[r.get(c, np.nan) for c in cols] for r in rows], dtype=float)
        keep = np.isfinite(x).any(axis=1)
        if keep.sum() < 2:
            continue
        f = fit_locus(x[keep], cols, mode=mode, n_components=n_components, winsor_p=winsor_p)
        if f is None:
            continue
        fits[locus] = f
        _, moved = clamp(x[keep], *f.bounds[winsor_p])
        obs = np.isfinite(x[keep]).sum()
        trimmed[locus] = float(moved.sum() / obs) if obs else 0.0
    return Corpus(sig=sig, name=name, vocab=vocab, fits=fits,
                  meta={**(meta or {}), "mode": mode, "winsor_p": winsor_p,
                        "n_components_requested": n_components,
                        "trimmed_frac": trimmed,
                        "platform": platform.platform(), "numpy": np.__version__})


# ----------------------------------------------------------------------------------- applying


def apply(raw: dict[str, float], chan: dict[str, float], corpus: Corpus, *,
          mode: "str | None" = None, winsor_p: "float | None" = None,
          n_components: "int | float | None" = None) -> dict[str, float]:
    """Rotate one sample's raw features through a corpus, and carry its channels through.

    Args:
        raw: Raw features for one sample, as :func:`~vdjtools.signature.features.raw_and_channels`
            returns them.
        chan: That sample's channels. Copied out untouched apart from ``qc:-:winsor_frac``.
        corpus: A fitted :class:`Corpus`.
        mode: Override the artifact's winsorization mode. ``"pcs"`` on an artifact fitted with
            ``"features"`` is allowed and is a legitimate thing to want, but the bounds it applies
            were then measured on an already-trimmed corpus; the manifest records which was fitted.
        winsor_p: Which stored percentile to clamp at. Defaults to the fitted one.
        n_components: Truncate to this many components (or this variance fraction). Exact -- the
            components are ordered, so the first ``n`` of a longer rotation are the same vectors.

    Returns:
        ``{column: value}`` in :func:`~vdjtools.signature.layout.signature_columns` order.
    """
    mode = mode or corpus.meta.get("mode", "features")
    if mode not in MODES:
        raise ValueError(f"unknown mode {mode!r}; known: {MODES}")
    p = corpus.meta.get("winsor_p", 0.01) if winsor_p is None else winsor_p
    ks = corpus.resolve_k(n_components)

    out: dict[str, float] = {}
    n_moved = n_seen = 0
    for locus, k in ks.items():
        f = corpus.fits[locus]
        x = np.array([raw.get(c, np.nan) for c in f.columns], dtype=float)
        n_seen += int(np.isfinite(x).sum())
        if mode == "features":
            if p not in f.bounds:
                raise ValueError(f"corpus {corpus.name!r} stores no bounds at p={p}; "
                                 f"available: {sorted(f.bounds)}")
            x, moved = clamp(x, *f.bounds[p])
            n_moved += int(moved.sum())
        z = _standardize(x, f.loc, f.scale)
        scores = z @ f.rotation[:, :k]
        safe = np.where(f.pc_scale[:k] > MIN_SPREAD, f.pc_scale[:k], 1.0)
        vals = (scores - f.pc_loc[:k]) / safe
        if mode == "pcs":
            if p not in f.pc_bounds:
                raise ValueError(f"corpus {corpus.name!r} stores no PC bounds at p={p}; "
                                 f"available: {sorted(f.pc_bounds)}")
            lo, hi = f.pc_bounds[p]
            vals, moved = clamp(vals, lo[:k], hi[:k])
            n_moved += int(moved.sum())
            n_seen += k
        for i in range(k):
            out[f"{corpus.sig}:{L.PC_BLOCK}:{locus}:PC{i + 1:02d}"] = float(vals[i])

    for c in L.channel_columns(corpus.sig):
        out[c] = float(chan.get(c, np.nan))
    wf = f"{corpus.sig}:qc:{L.NO_LOCUS}:winsor_frac"
    if wf in out:
        out[wf] = float(n_moved / n_seen) if n_seen else 0.0
    return out


# --------------------------------------------------------------------------- bundled corpora

#: Where the shipped artifacts live. A corpus is a rotation, a set of medians and MADs, and a set
#: of percentile bounds -- no sample is recoverable from any of it, which is why these ship publicly
#: while the per-sample matrices they were fitted on do not.
_RES = Path(__file__).resolve().parent.parent / "resources" / "signature"


def bundled_names() -> list[str]:
    """Corpus names with an installed artifact."""
    if not _RES.is_dir():
        return []
    return sorted(p.stem.removeprefix("vsig_") for p in _RES.glob("vsig_*.npz"))


def bundled_path(name: "str | Path") -> "Path | None":
    """Resolve a bundled corpus name or a filesystem path to an artifact, else ``None``."""
    direct = Path(name)
    if direct.suffix == ".npz" and direct.exists():
        return direct
    if direct.with_suffix(".npz").exists():
        return direct.with_suffix(".npz")
    cand = _RES / f"vsig_{name}.npz"
    return cand if cand.exists() else None


# ------------------------------------------------------------------- the synthetic corpora

#: Seed for the shipped synthetic corpora. Every draw is ``SEED + offset``, recorded in the
#: manifest, so a rebuild is reproducible by anyone who installs the library.
SEED: int = 20260927

#: Rank-abundance exponent for the ``memory`` corpus: clone ``i`` of ``n`` gets frequency
#: proportional to ``i ** -ZIPF_A``. Recorded in the artifact rather than assumed, because it is one
#: of the two free parameters of the selected-repertoire regime.
#:
#: NOTE this is the **Zipf law over ranks**, not ``numpy.random.Generator.zipf``, which samples
#: integers *from* a Zipf distribution. At ``a = 1.5`` that distribution has infinite mean, so
#: normalising a draw of it gives one clone almost all the mass: measured here, repertoires
#: collapsed to as few as 1 surviving clonotype, which then made every coverage-standardised
#: diversity feature a hole. The rank form is a well-behaved rank-abundance curve.
ZIPF_A: float = 1.5

#: Reads per clone in the ``memory`` multinomial. Sets how much of the rank-abundance tail survives
#: sampling, which is the realistic part: at 20 reads per clone a 500-clone repertoire keeps 336 of
#: them with 129 singletons, which is the shape a real library has. The zeros are dropped, so a
#: ``memory`` repertoire has fewer distinct clones than a ``naive`` one at the same nominal size --
#: exactly as a sampled selected repertoire does.
READS_PER_CLONE: int = 20

#: Receptors per synthetic repertoire, and the pool each is drawn from. A pool much larger than a
#: sample is what makes two synthetic repertoires nearly disjoint, the way two donors are.
DEFAULT_SIZE: int = 10_000
POOL_FACTOR: int = 100

#: Multiplicative p05-to-p95 spread of repertoire depth, per locus, from :data:`N_EFF`. Each
#: synthetic repertoire's size is drawn log-uniformly across this range around the nominal size.
#:
#: **A corpus at one fixed depth is not usable, and it fails silently.** Measured here: with every
#: ``naive`` repertoire at exactly 10,000 receptors, ``depth:reads`` and ``depth:richness`` are
#: identical in all N samples and all five ``pair:`` log-ratios are constant to the last bit -- so
#: their corpus spread is 0, they contribute nothing to the rotation, the cross-locus block cannot
#: be fitted at all, and a real sample's depth features get standardised against a degenerate
#: reference. Real repertoires span 4.1x (TRB) to 11.0x (IGH) between their 5th and 95th
#: percentiles, so drawing the depth is not an embellishment: it is what makes the corpus describe
#: the axis it is going to be asked about.
DEPTH_SPREAD: dict[str, float] = {loc: hi / lo for loc, (lo, _m, hi) in {
    "TRB": (195, 401, 790), "IGH": (111, 376, 1219), "TRA": (159, 276, 464),
    "IGK": (54, 203, 562), "IGL": (36, 123, 361), "TRG": (22, 49, 94), "TRD": (14, 34, 34),
}.items()}

#: Per-locus median effective clone count of real bulk blood repertoires, and the p05/p95 of the
#: same, measured on 1,168 samples as ``n_eff = 1/sum(w^2)`` with ``w = log2(1+count)/sum``. Used
#: by ``size="n_eff"`` and by the depth sweep, so the artifact can report how its bounds and its
#: scaling move with repertoire depth instead of leaving that unmeasured.
N_EFF: dict[str, tuple[int, int, int]] = {
    "TRB": (195, 401, 790), "IGH": (111, 376, 1219), "TRA": (159, 276, 464),
    "IGK": (54, 203, 562), "IGL": (36, 123, 361), "TRG": (22, 49, 94), "TRD": (14, 34, 34),
}


def draw_pool(locus: str, n: int, *, seed: int, source: str = "olga") -> pl.DataFrame:
    """``n`` productive rearrangements from the bundled model for one locus.

    Measured on a 16-core M-series laptop: 20,372 sequences/s on IGH, 23,113 on TRB, 34,506 on TRA,
    so a 10^7 pool is 8-14 minutes on one core and the seven loci are about ten minutes at one
    process per locus.
    """
    from ..model import load_bundled
    from ..model.generate import generate

    return generate(load_bundled(locus, source), n, seed=seed, productive_only=True)


def draw_sizes(nominal: int, n: int, spread: float, rng: np.random.Generator) -> np.ndarray:
    """``n`` repertoire sizes, log-uniform over a ``spread``-fold range around ``nominal``.

    Log-uniform rather than normal because depth is a scale, and the geometric centre is the
    nominal size so the corpus is centred where the caller asked for it. See :data:`DEPTH_SPREAD`
    for why a fixed depth is not an option.
    """
    half = np.sqrt(max(spread, 1.0))
    lo, hi = np.log(nominal / half), np.log(nominal * half)
    return np.maximum(np.exp(rng.uniform(lo, hi, size=n)).astype(np.int64), 2)


def draw_sample(pool: pl.DataFrame, size: int, regime: str, rng: np.random.Generator,
                ) -> pl.DataFrame:
    """One synthetic repertoire: ``size`` receptors drawn from ``pool``, with clone sizes.

    ``naive`` gives every clone a count of 1, which is what an unselected repertoire looks like and
    what the generator produces. ``memory`` draws Zipf frequencies, takes a multinomial sample at
    the same total, and **drops the zeros** -- so a memory repertoire has fewer distinct clones than
    a naive one at the same nominal size, exactly as a selected repertoire does.
    """
    if regime not in ("naive", "memory"):
        raise ValueError(f"regime must be 'naive' or 'memory'; got {regime!r}")
    idx = rng.choice(pool.height, size=min(size, pool.height), replace=False)
    df = pool[idx]
    if regime == "naive":
        return df.with_columns(pl.lit(1, dtype=pl.Int64).alias("duplicate_count"))
    # Zipf over RANKS: f_i ~ i^-a, normalised. Ranks are shuffled so clone size is independent of
    # where a sequence sat in the pool -- otherwise the abundance and the receptor are correlated
    # and the corpus learns an artefact of the draw order.
    rank = rng.permutation(df.height) + 1
    freq = rank.astype(float) ** -ZIPF_A
    freq /= freq.sum()
    counts = rng.multinomial(df.height * READS_PER_CLONE, freq)
    keep = counts > 0
    return df.filter(pl.Series(keep)).with_columns(
        pl.Series("duplicate_count", counts[keep].astype(np.int64)))


def synthesize(regime: str, *, loci: "tuple[str, ...]" = L.LOCI, n_samples: int = 10_000,
               size: "int | str" = DEFAULT_SIZE, seed: int = SEED,
               n_components: "int | float" = DEFAULT_COMPONENTS, mode: str = "features",
               winsor_p: float = 0.01, source: str = "olga", fit_corpus: bool = True,
               progress=None):
    """Build a synthetic corpus, and fit it.

    Args:
        regime: ``"naive"`` (every clone size 1) or ``"memory"`` (Zipf clone sizes).
        loci: Loci to build. All seven by default.
        n_samples: Repertoires in the corpus.
        size: Receptors per repertoire, or ``"n_eff"`` for the per-locus real medians in
            :data:`N_EFF`, or ``"p05"`` / ``"p95"`` for the sweep endpoints.
        seed: Base seed; every draw is a recorded offset from it.
        n_components: Components per locus, or a variance fraction.
        mode: Winsorization mode to fit at.
        winsor_p: Percentile for the fitted bounds. All of :data:`WINSOR_PS` are stored.
        source: Bundled model set -- ``"olga"``, ``"learned"`` or ``"arda"``.
        fit_corpus: ``False`` returns the raw samples instead of fitting, for scoring a held-out
            draw against a corpus fitted on another.
        progress: Optional ``callable(locus, done, total)``.

    Returns:
        ``(corpus, rows)`` when fitting, else a list of ``{locus: frame}`` samples.

    The whole build is a deterministic function of ``(regime, loci, n_samples, size, seed, source)``
    and the bundled models, so two machines at different thread counts must produce byte-identical
    artifacts. That is an acceptance criterion, not a hope: ``generate`` was once irreproducible
    *across processes* while being deterministic within one, because a marginal-table aggregation
    left its group order unspecified.
    """
    from . import features as FE

    def size_for(locus: str) -> int:
        if isinstance(size, int):
            return size
        if locus not in N_EFF:
            raise ValueError(f"size={size!r} has no measured depth for {locus}")
        return N_EFF[locus][{"p05": 0, "n_eff": 1, "p95": 2}[size]]

    vocab = {loc: FE.gene_vocab(loc) for loc in loci}
    pools, per_locus = {}, {}
    for i, locus in enumerate(sorted(loci, key=lambda x: L.LOCI.index(x))):
        n = size_for(locus)
        # The pool must cover the LARGEST repertoire the sweep can ask for, since a sample is
        # drawn without replacement.
        top = int(n * np.sqrt(DEPTH_SPREAD.get(locus, 4.0))) + 1
        pools[locus] = draw_pool(locus, min(top * POOL_FACTOR, max(top * n_samples, top)),
                                 seed=seed + i, source=source)
        rng = np.random.default_rng(seed + 1000 + i)
        sizes = draw_sizes(n, n_samples, DEPTH_SPREAD.get(locus, 4.0), rng)
        per_locus[locus] = [draw_sample(pools[locus], int(m), regime, rng) for m in sizes]
        if progress:
            progress(locus, i + 1, len(loci))
    samples = [{loc: per_locus[loc][j] for loc in per_locus} for j in range(n_samples)]
    if not fit_corpus:
        return samples

    rows = [FE.raw_and_channels(s, vocab)[0] for s in samples]
    from .. import __version__
    corpus = fit(rows, vocab, sig="vsig", name=regime, mode=mode, n_components=n_components,
                 winsor_p=winsor_p, meta={
                     "regime": regime, "n_samples": n_samples, "seed": seed, "source": source,
                     "size": size, "size_per_locus": {k: size_for(k) for k in loci},
                     "zipf_a": ZIPF_A if regime == "memory" else None,
                     "reads_per_clone": READS_PER_CLONE if regime == "memory" else 1,
                     "pool_factor": POOL_FACTOR, "model_version": __version__,
                     "depth_spread": {k: DEPTH_SPREAD.get(k) for k in loci},
                     "loci": list(loci)})
    return corpus, rows


def _demo() -> None:
    """Self-check: fit a tiny corpus, apply it, and hold the properties that matter."""
    import tempfile

    rng = np.random.default_rng(0)
    vocab = {"TRB": {"vus": ["TRBV1", "TRBV2"], "jus": ["TRBJ1-1"],
                     "spec": ["TRBV1_14", "TRBV2_14"]}}
    cols = L.raw_columns("vsig", "TRB", vocab["TRB"])
    rows = []
    for _ in range(200):
        r = {c: float(v) for c, v in zip(cols, rng.normal(size=len(cols)))}
        r |= {c: float(v) for c, v in zip(L.raw_columns("vsig", L.NO_LOCUS),
                                         rng.normal(size=5))}
        rows.append(r)
    rows[0][cols[0]] = 1e6                       # one pathological value, to be trimmed

    c = fit(rows, vocab, sig="vsig", name="demo", n_components=4)
    assert set(c.fits) == {"TRB", L.NO_LOCUS}, sorted(c.fits)
    assert c.k == {"TRB": 4, L.NO_LOCUS: 4}

    # winsorization trimmed the outlier, and only on the side the support allows
    lo, hi = c.fits["TRB"].bounds[0.01]
    i = cols.index(cols[0])
    assert hi[i] < 1e6, "a nonneg/real column left its top tail untrimmed"
    for j, col in enumerate(cols):
        bottom, top = L.SUPPORTS[L.support_of(col)]
        assert np.isfinite(lo[j]) == bottom, f"{col}: bottom bound against its support"
        assert np.isfinite(hi[j]) == top, f"{col}: top bound against its support"

    raw, chan = rows[1], {ch: 0.5 for ch in L.channel_columns("vsig")}
    sig = apply(raw, chan, c)
    assert set(sig) == set(c.columns()), "apply and columns() disagree"
    assert np.isfinite(list(sig.values())).all()

    # truncation is EXACT: the first n of a longer rotation are the same vectors
    short = apply(raw, chan, c, n_components=2)
    for loc in ("TRB", L.NO_LOCUS):
        for i in (1, 2):
            a, b = sig[f"vsig:pc:{loc}:PC{i:02d}"], short[f"vsig:pc:{loc}:PC{i:02d}"]
            assert a == b, f"truncation moved {loc} PC{i}: {a} != {b}"
    assert "vsig:pc:TRB:PC03" not in short

    # and asking for more than was fitted refuses rather than padding
    try:
        apply(raw, chan, c, n_components=99)
    except ValueError as e:
        assert "refit" in str(e)
    else:
        raise AssertionError("truncation padded instead of refusing")

    # the clamp is reported, never silent
    wild = dict(raw); wild[cols[0]] = 1e9
    assert apply(wild, chan, c)["vsig:qc:-:winsor_frac"] > 0.0
    assert apply(raw, chan, c, mode="none")["vsig:qc:-:winsor_frac"] == 0.0

    # channels pass through untouched
    assert sig["vsig:cov:TRB:cstar"] == 0.5

    # round trip, including the load-time column-order check
    with tempfile.TemporaryDirectory() as d:
        path = c.save(Path(d) / "demo")
        back = Corpus.load(path)
        assert back.k == c.k and back.vocab == c.vocab
        b = apply(raw, chan, back)
        # float32 rotation on disk, so compare at the precision that survives the round trip
        for kk in sig:
            assert abs(sig[kk] - b[kk]) < 1e-4, kk
        assert len(back.meta["content_sha256"]) == 64
        # a layout that no longer matches must raise, not silently mix features
        back.vocab["TRB"]["vus"] = ["TRBV1", "TRBV999"]
        try:
            back.verify()
        except ValueError as e:
            assert "indexed by position" in str(e)
        else:
            raise AssertionError("a stale column order loaded silently")

    print(f"corpus OK  k={c.k}  variance@4={c.fits['TRB'].variance_at(4):.3f}")


if __name__ == "__main__":
    _demo()
