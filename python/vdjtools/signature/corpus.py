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
                is not reached by the stored spectrum. Neither ever pads: a count is capped at
                each locus's stored ``k`` (the same cap the fit applied, so a 5-feature cross-locus
                block yields 5 whether you asked for 5 or 128), and a variance fraction the stored
                spectrum cannot reach raises rather than silently returning the whole rotation
                under the requested name.
        """
        if n_components is None:
            return self.k
        out: dict[str, int] = {}
        for loc, f in self.fits.items():
            if isinstance(n_components, int) and not isinstance(n_components, bool):
                # A count is an upper bound per locus, exactly as it was at fit time: a block with
                # fewer non-degenerate directions than asked for contributes the ones it has. This
                # has to mirror the fit or the flag is unusable on anything the fit produced -- the
                # cross-locus block has 5 raw features, so `--components 128` fitted it at 5 and
                # then refused to apply itself. Nothing is padded: no column is invented under a
                # name it does not hold, which is the property the refusal existed for.
                out[loc] = min(int(n_components), f.k)
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
        recorded = self.meta.get("models") or {}
        if recorded:
            have = model_fingerprint(sorted(recorded), self.meta.get("source") or "olga")
            moved = {loc: (was, have.get(loc)) for loc, was in recorded.items()
                     if have.get(loc) != was}
            if moved:
                detail = "; ".join(f"{loc}: fitted on {was}, installed {now}"
                                   for loc, (was, now) in sorted(moved.items()))
                raise ValueError(
                    f"corpus {self.name!r} was fitted on bundled recombination models that are no "
                    f"longer the installed ones ({detail}). Every synthetic repertoire in it was "
                    f"drawn from those models, and retraining one moves its locus by 0.1-1.6%, so "
                    f"the mismatch is not cosmetic. Rebuild the corpus, or load with verify=False "
                    f"if you have established the models did not change.")


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


def locus_matrix(vocab: dict[str, dict[str, list[str]]], sig: str, locus: str,
                 n: int) -> "tuple[np.ndarray, list[str]] | None":
    """A preallocated ``(n, p_L)`` buffer and its column order, or ``None`` if the locus is empty.

    Preallocated and filled in place rather than built from a list of per-sample dicts. At the
    shipped corpus size a dict of 13,483 Python floats per sample is 4-5 GB of interpreter objects
    for a matrix that is 1.08 GB as float64 -- and the dicts have to coexist with the array while it
    is assembled. See :func:`fill_row`.
    """
    cols = L.raw_columns(sig, locus, vocab.get(locus))
    if not cols:
        return None
    return np.full((n, len(cols)), np.nan), cols


def fill_row(buf: np.ndarray, cols: list[str], i: int, raw: dict[str, float]) -> None:
    """Write one sample's values for one locus into row ``i`` of a preallocated buffer."""
    buf[i] = [raw.get(c, np.nan) for c in cols]


def fit_matrices(mats: dict[str, tuple[np.ndarray, list[str]]],
                 vocab: dict[str, dict[str, list[str]]], *, sig: str, name: str,
                 mode: str = "features", n_components: "int | float" = DEFAULT_COMPONENTS,
                 winsor_p: float = 0.01, meta: "dict | None" = None) -> Corpus:
    """Fit a corpus from per-locus matrices -- the memory-bounded path.

    Each locus is fitted only on the rows that observed it, so a corpus where half the samples have
    no TRD does not learn TRD from imputed values.
    """
    fits: dict[str, LocusFit] = {}
    trimmed: dict[str, float] = {}
    for locus, (x, cols) in mats.items():
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


def fit(rows: list[dict[str, float]], vocab: dict[str, dict[str, list[str]]], *, sig: str,
        **kw) -> Corpus:
    """Fit a corpus from raw feature rows -- one dict per repertoire.

    Convenient for a small corpus and for tests. For a full-size build use
    :func:`locus_matrix` / :func:`fill_row` / :func:`fit_matrices`, which never holds more than one
    sample's dict at a time.
    """
    order = [loc for loc in (*L.LOCI, L.NO_LOCUS) if loc in vocab or loc == L.NO_LOCUS]
    mats: dict[str, tuple[np.ndarray, list[str]]] = {}
    for locus in order:
        got = locus_matrix(vocab, sig, locus, len(rows))
        if got is None:
            continue
        buf, cols = got
        for i, r in enumerate(rows):
            fill_row(buf, cols, i, r)
        mats[locus] = (buf, cols)
    return fit_matrices(mats, vocab, sig=sig, **kw)


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
        # Rotate through the WHOLE rotation, then slice. `z @ rotation[:, :k]` is the same
        # arithmetic in a different GEMM shape, and OpenBLAS accumulates it in a different
        # order: the truncated and untruncated TRB PC01 differed at 2e-15 on Linux while
        # agreeing bit-for-bit on Accelerate, which is how the exactness test passed here and
        # failed in CI. Slicing the full product makes truncation exact by construction --
        # both paths evaluate the identical expression -- at a few microseconds per locus.
        scores = (z @ f.rotation)[:k]
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


def model_fingerprint(loci, source: str = "olga") -> dict[str, str]:
    """``{locus: "<model source>@<model version>"}`` for the models a synthetic pool is drawn from.

    The version is a hash of the **germline the generator will actually draw from** (every allele's
    CDR3-region cut segment of the collapsed model), not the manifest's declared version string. A
    declared version does not move when a germline is repaired: correcting ``TRBV4-3*02``'s anchor
    changed what every TRB pool contains while leaving ``olga:human_T_beta@2.0.0`` identical, so a
    corpus fitted before the repair would have loaded silently against models that no longer produce
    it. The hash cannot miss that.

    This, not the library version, is what a synthetic corpus's numbers rest on: every pool comes
    out of these models, and retraining one moves that locus by 0.1-1.6%. Gating on the library
    version instead was both too strict -- a patch release invalidated every corpus for no reason --
    and **wrong**: ``vdjtools.__version__`` is read from installed distribution metadata, so a build
    driven by ``PYTHONPATH`` against a different installed version records *that* version rather
    than the code that ran. The first four cluster artifacts recorded ``3.6.0`` for exactly that
    reason, from a 4.0.0 tree.

    A real corpus is fitted on real repertoires and no model enters it, so its manifest carries no
    fingerprint and :meth:`Corpus.verify` has nothing to check.
    """
    import hashlib

    from ..model import load_bundled

    out = {}
    for loc in loci:
        m = load_bundled(loc, source)
        h = hashlib.sha256()
        for name in sorted(m.genomic):
            col = f"{name.split('_')[1]}_allele"
            for a, c in sorted(zip(m.genomic[name][col], m.genomic[name]["cut_segment"])):
                h.update(f"{a}\t{c}\n".encode())
        out[loc] = f"{m.manifest.source}@{h.hexdigest()[:12]}"
    return out


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


#: What a REAL repertoire looks like per locus, in the two compartments the mixture corpora are named
#: after: ``{cohort: {locus: (samples, richness, reads per expanded clone, singleton fraction, rank
#: correlations)}}``, each of the three quantities being its measured
#: p05/p25/p50/p75/p95 :data:`ladder <COHORT_QS>`.
#:
#: Measured 2026-09-27 in one streaming pass over the harmonized AIRR store, every sample with at
#: least 100 reads in the locus, grouped by ``sample_id``. **Richness** is distinct clonotypes;
#: **f1** is the fraction of them seen exactly once; **reads per expanded clone** is
#: ``(reads - singletons) / (richness - singletons)``. The sample count leads each row, so a reader
#: can see which loci the table is thin on (tissue TRD, 1,180 samples) and which it is not (tissue
#: IGK, 58,706).
#:
#: **These three ladders are the entire parameterisation of a mixture corpus**, and the middle one is
#: deliberately not the obvious quantity. Reads per *clonotype* cannot be drawn independently of the
#: singleton fraction: a repertoire with ``f1`` singletons whose other clones all carry at least 2
#: reads has at least ``2 - f1`` reads per clonotype, so a pair drawn from those two marginals lands
#: in a forbidden region about half the time on blood TRB and has to be clipped back out of it --
#: which silently rewrites the marginal it was drawn from. Drawing the read count itself is worse
#: still: measured here, a copula over (richness, reads, f1) asks for a repertoire below that floor
#: for **21.1%** of 20,000 blood TRB draws. Reads per *expanded* clone is >= 2 whatever ``f1`` is, so
#: the trio has no forbidden region at all.
#:
#: **Five quantiles, not two, and it is the resolution that makes the corpus land.** The read count is
#: a product of all three drawn quantities, so it is the sharpest check on the draw: against blood
#: TRB's measured median of 763 reads (n = 34,365 samples), a 40,000-draw simulation gives **762** at
#: five quantiles and 846 at three; blood IGH's 1,298 comes out 1,276 against 1,490 (n = 33,245).
#: Tissue is the harder case and is stated rather than smoothed over: tissue TRB's median 308 comes
#: out 425 at five quantiles, against 756 at three (n = 22,298). A product of three heavy-tailed
#: factors has a median above the product of their medians unless their joint tails are matched
#: exactly, and no marginals-plus-copula draw does that.
#:
#: The last field is the Spearman rank correlations among the three, in the order ``(richness vs
#: expanded count, richness vs f1, expanded count vs f1)``, drawn through a Gaussian copula.
#: **That is for the rotation, not for the read count**: a corpus's rotation IS its covariance
#: structure, and ``f1`` against expansion size is -0.503 on blood TRB and -0.846 on blood IGK
#: (n = 43,678) -- fewer singletons, bigger expansions. A corpus that drew the two independently
#: would hand the PCA a correlation the cohort does not have, which no amount of correct marginals
#: repairs.
#:
#: Against all this, :data:`DEPTH_SPREAD` -- 2.4x on TRD to 11.0x on IGH, from 1,168 deep blood
#: samples, and what the ``naive``/``memory`` corpora draw across -- covers about a tenth of the real
#: axis: blood TRB richness spans 43x (74 to 3,162 clonotypes) and tissue IGH 259x (40 to 10,352). A
#: corpus estimates its bounds, its centre and its per-PC scaling from its own draw, and not one of
#: the three extrapolates beyond it.
COHORT: dict[str, dict[str, tuple]] = {
    "blood": {
        "TRA": (32732, (67, 148, 276, 554, 1724), (2.08, 2.38, 2.75, 3.50, 6.98),
                (0.197, 0.566, 0.745, 0.851, 0.942), (-0.045, 0.142, -0.547)),
        "TRB": (34365, (74, 208, 457, 985, 3162), (2.10, 2.40, 2.83, 3.74, 7.91),
                (0.215, 0.580, 0.759, 0.860, 0.949), (-0.036, 0.128, -0.503)),
        "TRG": (8688, (21, 48, 72, 108, 220), (2.64, 3.51, 4.74, 7.00, 15.38),
                (0.167, 0.438, 0.583, 0.686, 0.806), (-0.421, 0.463, -0.394)),
        "TRD": (5620, (19, 48, 80, 121, 280), (2.25, 2.96, 4.17, 6.75, 17.36),
                (0.167, 0.418, 0.591, 0.714, 0.856), (-0.547, 0.484, -0.482)),
        "IGH": (33245, (84, 271, 697, 1711, 6200), (2.11, 2.42, 2.88, 3.93, 12.83),
                (0.255, 0.528, 0.697, 0.819, 0.927), (0.204, -0.132, -0.717)),
        "IGK": (43678, (80, 233, 549, 1172, 3432), (2.29, 2.83, 3.67, 5.41, 15.97),
                (0.143, 0.384, 0.555, 0.704, 0.846), (0.399, -0.301, -0.759)),
        "IGL": (37601, (61, 160, 361, 790, 2468), (2.33, 2.86, 3.63, 5.37, 15.89),
                (0.122, 0.346, 0.514, 0.671, 0.817), (0.349, -0.259, -0.714)),
    },
    "tissue": {
        "TRA": (16799, (18, 84, 148, 417, 1839), (2.10, 2.42, 2.95, 4.68, 120.24),
                (0.111, 0.476, 0.668, 0.806, 0.930), (-0.320, 0.340, -0.688)),
        "TRB": (22298, (30, 89, 158, 498, 2977), (2.11, 2.48, 3.00, 4.60, 94.74),
                (0.194, 0.517, 0.688, 0.817, 0.928), (-0.233, 0.265, -0.598)),
        "TRG": (2235, (7, 31, 57, 93, 172), (2.85, 4.42, 6.82, 13.29, 102.00),
                (0.077, 0.309, 0.500, 0.661, 0.820), (-0.540, 0.462, -0.377)),
        "TRD": (1180, (6, 21, 54, 96, 253), (2.55, 4.05, 9.29, 53.40, 373.45),
                (0.120, 0.400, 0.594, 0.728, 0.876), (-0.661, 0.373, -0.219)),
        "IGH": (47031, (40, 164, 602, 2458, 10352), (2.39, 3.30, 5.27, 10.00, 40.36),
                (0.055, 0.395, 0.561, 0.690, 0.855), (0.017, 0.079, -0.509)),
        "IGK": (58706, (34, 150, 519, 1813, 6592), (2.66, 4.09, 7.05, 14.90, 69.34),
                (0.147, 0.290, 0.410, 0.546, 0.750), (0.134, -0.147, -0.582)),
        "IGL": (51010, (29, 113, 343, 1188, 4152), (2.64, 3.88, 6.49, 13.09, 54.70),
                (0.125, 0.260, 0.377, 0.510, 0.717), (0.075, -0.094, -0.565)),
    },
}

#: The quantiles every :data:`COHORT` ladder is recorded at.
COHORT_QS: tuple[float, ...] = (0.05, 0.25, 0.50, 0.75, 0.95)

#: The named synthetic corpora: ``{name: (regime, cohort)}``. There are four, and the name is the
#: only thing a caller has to know.
#:
#: ``naive`` and ``memory`` are the two pure regimes at one nominal size. Neither is representative
#: of a real bulk sample and neither can be: ``naive`` has no clone-size structure at all, and
#: ``memory`` has no naive background and, at a fixed size, no read-depth variance either -- its read
#: count is exactly ``size * READS_PER_CLONE`` in every sample. The two ``synthetic-*`` corpora are
#: the MIXTURE, drawn across the measured per-locus richness, read-depth and singleton-fraction
#: ranges of the cohort they are named after, so they belong in the same set as the real ``blood`` /
#: ``tissue`` / ``deep-tcr`` corpora and mean what their names say.
SYNTHETIC: dict[str, tuple[str, "str | None"]] = {
    "naive": ("naive", None),
    "memory": ("memory", None),
    "synthetic-blood": ("mixed", "blood"),
    "synthetic-tissue": ("mixed", "tissue"),
}


def corpus_plan(name: str, *, size: "int | str | None" = None,
                depth_spread: "float | str | None" = None) -> tuple:
    """``(regime, cohort, size, depth_spread)`` for a named corpus -- the ONE place a name becomes a
    draw, so the two halves of a corpus cannot resolve the same name differently.

    A cohort name stands in for both ``size`` and ``depth_spread``, because both come out of the same
    measured band: the nominal size is its geometric centre and the spread is its width. An explicit
    value for either wins, which is what makes a deliberately wider or narrower variant of a named
    corpus a one-flag change rather than a new table.
    """
    if name not in SYNTHETIC:
        raise ValueError(
            f"unknown synthetic corpus {name!r}; one of {', '.join(SYNTHETIC)}, or a path to an "
            "artifact. 'mixed' is the regime, not a corpus: a mixture has no unparameterised form, "
            "and the cohort it is drawn across is what makes it one.")
    regime, cohort = SYNTHETIC[name]
    return (regime, cohort,
            (cohort or DEFAULT_SIZE) if size is None else size,
            cohort if depth_spread is None else depth_spread)


def cohort_bands(cohort: str, locus: str) -> tuple[tuple, tuple, tuple, tuple]:
    """One locus of one cohort: the richness, expanded-count and ``f1`` ladders, and their copula."""
    _samples, rich, mexp, frac, corr = COHORT[cohort][locus]
    return rich, mexp, frac, corr


def ladder_draw(q: tuple, u: np.ndarray, log: bool = True) -> np.ndarray:
    """Inverse-CDF draw from a measured :data:`COHORT_QS` ladder, at probabilities ``u``.

    Linear interpolation of the quantile function -- in log space for a scale quantity -- so the drawn
    marginal carries **all five** of the cohort's recorded quantiles and spans **exactly** its p05 to
    p95.

    That is the difference from drawing uniformly across the band, and it is not cosmetic: real
    singleton fractions sit near the top of theirs (blood TRB median 0.759 of a 0.215-0.949 band,
    n = 34,365 samples), so a uniform draw would centre the corpus at 0.582 and leave a
    cohort-median sample 0.65 robust SD off centre on every singleton-sensitive feature.
    """
    y = np.interp(u, COHORT_QS, np.log(q) if log else np.asarray(q, dtype=float))
    return np.exp(y) if log else y


def copula_uniforms(rank_corr: tuple, n: int, rng: np.random.Generator) -> np.ndarray:
    """``(n, 3)`` probabilities carrying the measured rank correlations -- a Gaussian copula.

    Spearman's rho converts to the Gaussian correlation it implies by ``2 * sin(pi * rho / 6)``, the
    three of them are made a matrix, and correlated normals are pushed through the normal CDF. The
    result is uniform in every column whatever the correlations are, so :func:`ladder_draw` still
    reproduces each measured marginal exactly -- a copula only decides how they move together, which
    is what the rotation is fitted on.

    Measured over all 14 (cohort, locus) entries of :data:`COHORT`: the realised rank correlations
    match their targets to within **0.015** on 20,000 draws, and the implied Gaussian matrix is
    comfortably positive semidefinite everywhere -- smallest eigenvalue 0.2176, on blood IGK. The
    eigenvalue clip therefore never binds on the shipped table and is defensive, for a caller passing
    three correlations that cannot come from one joint distribution; the row renormalisation is what
    restores unit marginal variance when it does bind, without which a clipped matrix would quietly
    narrow the drawn bands.
    """
    from scipy.special import ndtr

    r_nm, r_nf, r_mf = (2.0 * np.sin(np.pi * np.asarray(rank_corr, dtype=float) / 6.0))
    m = np.array([[1.0, r_nm, r_nf], [r_nm, 1.0, r_mf], [r_nf, r_mf, 1.0]])
    w, v = np.linalg.eigh(m)
    root = v * np.sqrt(np.clip(w, 0.0, None))
    z = rng.standard_normal((n, 3)) @ root.T
    z /= np.sqrt(np.maximum((root ** 2).sum(axis=1), 1e-12))
    # Onto the ladder's own support: a monotone map, so the copula's ranks survive it, and without
    # it the 10% of draws outside p05-p95 would pile up as point masses on the two band edges.
    return COHORT_QS[0] + (COHORT_QS[-1] - COHORT_QS[0]) * ndtr(z)


def depth_spread_of(locus: str, depth_spread: "float | str | None" = None) -> float:
    """The multiplicative depth range for one locus: measured, named, or a caller's override.

    ``None`` means :data:`DEPTH_SPREAD`, the p05-p95 spread of that locus among 1,168 deep blood
    samples, which is what ``naive`` and ``memory`` draw across -- 2.4x on TRD to 11.0x on IGH, so
    under one decade everywhere. A cohort name from :data:`COHORT` means that cohort's own measured
    richness band, which is 43x on blood TRB and 259x on tissue IGH. A number widens or narrows it
    deliberately: a corpus meant to describe samples across decades of sequencing depth has to be
    *drawn* across those decades, because the bounds, the centre and the per-PC scaling are all
    estimated from the draw and none of them extrapolates.
    """
    if depth_spread is None:
        return DEPTH_SPREAD.get(locus, 4.0)
    if isinstance(depth_spread, str):
        rich = cohort_bands(depth_spread, locus)[0]
        return rich[-1] / rich[0]
    return float(depth_spread)


def draw_pool(locus: str, n: int, *, seed: int, source: str = "olga") -> pl.DataFrame:
    """``n`` productive rearrangements from the bundled model for one locus.

    Measured 2026-09-27, ``generate(load_bundled("IGH"), 50_000, productive_only=True)`` on one
    core of a 16-core M-series laptop: **4,200 sequences/s**, and identical at 1 and 16 polars
    threads -- the sampler is a per-sequence Python/numpy loop, so it does not thread. Pool
    generation is therefore the dominant serial cost of a synthetic build (about 15M sequences at
    ``size=10,000`` across seven loci) and is why :func:`build_pools` fans out over processes.

    NOTE an earlier docstring here claimed 20,372 seq/s on IGH. That figure is 4.8x the measured
    rate and it is what made a 77 s stage look like a 16 s one when the build was planned.
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
                frac: "float | None" = None, mexp: "float | None" = None) -> pl.DataFrame:
    """One synthetic repertoire: ``size`` receptors drawn from ``pool``, with clone sizes.

    ``naive`` gives every clone a count of 1, which is what an unselected repertoire looks like and
    what the generator produces. ``memory`` draws Zipf frequencies, takes a multinomial sample at
    the same total, and **drops the zeros** -- so a memory repertoire has fewer distinct clones than
    a naive one at the same nominal size, exactly as a selected repertoire does.

    ``mixed`` is what a real bulk sample is -- a singleton naive background plus an expanded memory
    component -- and it needs ``frac`` (the singleton fraction) and ``mexp`` (reads per expanded
    clone), both drawn per sample by :func:`draw_plan` from the cohort's measured ladders. Neither
    pure regime can stand in for it: ``naive`` has no clone-size structure at all, and ``memory`` at
    a fixed size has a read count of exactly ``size * READS_PER_CLONE`` in every sample, so its
    depth does not vary either.
    """
    if regime not in ("naive", "memory", "mixed"):
        raise ValueError(f"regime must be 'naive', 'memory' or 'mixed'; got {regime!r}")
    idx = rng.choice(pool.height, size=min(size, pool.height), replace=False)
    df = pool[idx]
    if regime == "naive":
        return df.with_columns(pl.lit(1, dtype=pl.Int64).alias("duplicate_count"))
    if regime == "mixed":
        if frac is None or mexp is None:
            raise ValueError("regime='mixed' needs frac (the singleton fraction) and mexp (reads per "
                             "expanded clone); a mixture has no unparameterised form. Build one of "
                             f"{', '.join(k for k, v in SYNTHETIC.items() if v[0] == 'mixed')}, "
                             "whose per-locus ladders come from COHORT.")
        return _mixed_counts(df, rng, float(frac), float(mexp))
    return _zipf_counts(df, rng)


def _zipf_freq(n: int, rng: np.random.Generator) -> np.ndarray:
    """Zipf rank-abundance frequencies for ``n`` clones, over SHUFFLED ranks.

    Zipf over RANKS: ``f_i ~ i**-a``, normalised. The ranks are shuffled so clone size is
    independent of where a sequence sat in the pool -- otherwise the abundance and the receptor are
    correlated and the corpus learns an artefact of the draw order. One definition, shared by
    ``memory`` and the memory half of ``mixed``: two copies would let the two drift apart silently.
    """
    rank = rng.permutation(n) + 1
    freq = rank.astype(float) ** -ZIPF_A
    return freq / freq.sum()


def _zipf_counts(df: pl.DataFrame, rng: np.random.Generator) -> pl.DataFrame:
    """``memory`` clone sizes: a multinomial over Zipf frequencies, zeros dropped."""
    freq = _zipf_freq(df.height, rng)
    counts = rng.multinomial(df.height * READS_PER_CLONE, freq)
    keep = counts > 0
    return df.filter(pl.Series(keep)).with_columns(
        pl.Series("duplicate_count", counts[keep].astype(np.int64)))


def _mixed_counts(df: pl.DataFrame, rng: np.random.Generator, frac: float,
                  mexp: float) -> pl.DataFrame:
    """``mixed`` clone sizes: ``frac`` of the clones at one read, the rest sharing ``mexp`` each.

    Hits the drawn richness, singleton fraction and read count **exactly**, which is the whole point
    of constructing the mixture rather than sampling it: ``round(frac * n)`` clonotypes get a count
    of 1 -- that is what the singleton fraction means operationally, and it is the quantity our aging
    cohorts report as the naive share -- and the remaining ``ne`` clones split ``mexp * ne`` reads
    with Zipf rank-abundance frequencies on a floor of 2, so no expanded clone can fall back into the
    singleton class and turn the realised ``f1`` into a consequence of the sampling.

    Contrast ``memory``, where the multinomial's zeros are dropped: there both the richness and the
    singleton fraction come out of the draw rather than being asked for, which is why a memory corpus
    cannot be pointed at a measured cohort.
    """
    n = df.height
    n1 = min(int(round(frac * n)), n)
    ne = n - n1
    if ne == 0:
        return df.with_columns(pl.lit(1, dtype=pl.Int64).alias("duplicate_count"))
    counts = np.empty(n, dtype=np.int64)
    counts[:n1] = 1
    # `mexp` is reads per EXPANDED clone and is >= 2 by construction, so the floor cannot bind for a
    # drawn plan; it is here for a caller passing its own number.
    extra = max(int(round(mexp * ne)), 2 * ne) - 2 * ne
    counts[n1:] = 2 + rng.multinomial(extra, _zipf_freq(ne, rng))
    return df.with_columns(pl.Series("duplicate_count", counts))


def resolved_size(locus: str, size: "int | str") -> int:
    """The per-locus receptor count a ``size`` asks for -- an integer, or a measured depth name.

    Separate from :func:`sample_stream` because the manifest has to record what each locus was
    actually drawn at: ``size="n_eff"`` means seven different depths, and a manifest saying only
    ``"n_eff"`` cannot tell a reader which.
    """
    if isinstance(size, int):
        return size
    if size in COHORT:
        # The cohort's MEDIAN richness for that locus. `draw_plan` rescales the whole measured ladder
        # by `resolved_size / median`, so this is the nominal depth the ladder is centred on.
        return cohort_bands(size, locus)[0][COHORT_QS.index(0.50)]
    if locus not in N_EFF:
        raise ValueError(f"size={size!r} has no measured depth for {locus}")
    return N_EFF[locus][{"p05": 0, "n_eff": 1, "p95": 2}[size]]


def draw_plan(loci, *, n_samples: int, size: "int | str", seed: int,
              depth_spread: "float | str | None" = None, cohort: "str | None" = None) -> dict:
    """Every sample's drawn size, singleton fraction and expanded count: ``{locus: {name: array}}``.

    Drawn **in the parent**, once, from one generator per locus seeded at ``seed + 1000 + i``, so
    every worker is handed the same plan rather than reproducing it -- and so the plan itself can be
    inspected, which is how the realised corpus gets checked against the cohort it is named after.

    ``cohort=None`` leaves the two mixture knobs unset, which is correct for ``naive`` and ``memory``
    and keeps their draw byte-identical to what the shipped artifacts were built from.
    """
    # Compared against the cohort itself, not merely type-checked: `depth_spread="tissue"` on a blood
    # corpus would otherwise pass, change nothing, and be recorded in the manifest as the tissue
    # spread. A knob that silently does nothing is the failure this whole subsystem was rewritten to
    # end.
    if cohort is not None and depth_spread != cohort:
        raise ValueError(
            f"corpus cohort {cohort!r} draws depth from its measured richness ladder "
            f"({COHORT_QS[0]:.2f}-{COHORT_QS[-1]:.2f} of it), so depth_spread={depth_spread!r} does "
            "not apply. Move the whole ladder with `size`, which rescales it about its median, or "
            "build naive/memory, which is where an explicit spread belongs.")
    plan = {}
    for i, locus in enumerate(loci):
        rng = np.random.default_rng(seed + 1000 + i)
        if cohort is None:
            plan[locus] = {"size": draw_sizes(resolved_size(locus, size), n_samples,
                                              depth_spread_of(locus, depth_spread), rng),
                           "frac": None, "mexp": None}
            continue
        rich, mexp, frac, corr = cohort_bands(cohort, locus)
        # `size` rescales the whole ladder about its median rather than being ignored, so a smoke
        # build (--size 1000) is the same cohort SHAPE at a different nominal depth. It is 1.0 for
        # `size=<cohort>`, which is what the corpus's own default resolves to.
        scale = resolved_size(locus, size) / rich[COHORT_QS.index(0.50)]
        u = copula_uniforms(corr, n_samples, rng)
        plan[locus] = {
            "size": np.maximum((ladder_draw(rich, u[:, 0]) * scale).astype(np.int64), 2),
            "frac": ladder_draw(frac, u[:, 2], log=False),
            "mexp": ladder_draw(mexp, u[:, 1])}
    return plan


def draw_one(pools: dict, plan: dict, regime: str, j: int, seed: int) -> dict:
    """Sample ``j`` of a synthetic corpus -- a **pure function of** ``j``.

    Each sample gets its own generator, seeded from ``(seed, locus index, j)``, so any process can
    draw any sample in any order and get the same repertoire. That is what lets the build run across
    processes while staying bit-identical at every ``n_jobs``: the drawing is the only part of the
    pipeline carrying RNG state, and here it carries none between samples.
    """
    def at(v, j):
        return None if v is None else float(v[j])

    return {loc: draw_sample(pools[loc], int(plan[loc]["size"][j]), regime,
                             np.random.default_rng([seed, i, j]),
                             at(plan[loc]["frac"], j), at(plan[loc]["mexp"], j))
            for i, loc in enumerate(pools)}


def pool_target(locus: str, size: "int | str", n_samples: int,
                depth_spread: "float | str | None" = None, cohort: "str | None" = None) -> int:
    """How many receptors a locus's pool needs: the largest repertoire, times :data:`POOL_FACTOR`.

    The pool has to cover the LARGEST depth the draw can ask for, since a sample is drawn from it
    without replacement. For a cohort corpus that is the top of the measured richness ladder, scaled
    the same way :func:`draw_plan` scales it; for a pure regime it is the nominal size times the
    half-spread.
    """
    if cohort is not None:
        rich = cohort_bands(cohort, locus)[0]
        top = int(rich[-1] * resolved_size(locus, size) / rich[COHORT_QS.index(0.50)]) + 1
    else:
        top = int(resolved_size(locus, size) * np.sqrt(depth_spread_of(locus, depth_spread))) + 1
    return min(top * POOL_FACTOR, max(top * n_samples, top))


def _pool_chunk(args) -> tuple:
    """Generate one contiguous chunk of one locus's pool; write it as uncompressed Arrow IPC."""
    import os

    from .cohort import WORKER_ENV

    os.environ[WORKER_ENV] = "1"
    locus, k, n, seed, source, tmp = args
    path = Path(tmp) / f"pool_{locus}_{k:03d}.arrow"
    draw_pool(locus, n, seed=seed, source=source).write_ipc(path, compression="uncompressed")
    return locus, k, path


def build_pools(loci, *, size, n_samples: int, seed: int, source: str, tmp: Path,
                n_jobs: int = 1, progress=None, depth_spread: "float | str | None" = None,
                cohort: "str | None" = None) -> dict:
    """Generate every locus's pool across processes; return ``{locus: [chunk path, ...]}``.

    Generation is the single most expensive stage of a synthetic build and it is embarrassingly
    parallel -- measured on Aldan-3, seven pools at ``size=10,000`` took 4,850 s in one process,
    which is 81 minutes of a build whose featurisation is the part anyone cares about.

    Written as **uncompressed** IPC and read back memory-mapped, so N workers share one physical
    copy of a multi-gigabyte pool through the page cache instead of each pickling its own.
    """
    single_threaded_children()
    ordered = sorted(loci, key=lambda x: L.LOCI.index(x))
    tasks = []
    for i, locus in enumerate(ordered):
        total = pool_target(locus, size, n_samples, depth_spread, cohort)
        # Chunk sizes differ by at most one and the seeds are per chunk, so the concatenated pool
        # is a deterministic function of (locus, total, seed) and of the chunk count -- which is
        # why the chunk count is the fixed :data:`POOL_CHUNKS` and never the machine's core count.
        # Seven loci would otherwise cap the fan-out at seven, whatever the machine: measured, that
        # left 5.2M sequences taking 296 s on a 16-core box.
        for k, (a, b) in enumerate(_ranges(total, POOL_CHUNKS)):
            tasks.append((locus, k, b - a, seed + i * 1_000 + k, source, tmp))
    # One task at a time to the pool, NOT a contiguous slice per worker. The opposite of the
    # featurise stage, and for a measured reason: there are few tasks and they are large and
    # unequal (IGH's pool is 1.6x TRB's), so contiguous slicing leaves workers idle behind the big
    # ones -- it cost 122 s against a 77 s floor. Per-task submission overhead is irrelevant next
    # to a chunk that takes tens of seconds.
    done = _map_unordered(_pool_chunk, tasks, n_jobs)
    out: dict = {loc: [] for loc in ordered}
    for locus, k, path in sorted(done, key=lambda t: (t[0], t[1])):
        out[locus].append(path)
    # Per LOCUS, not per chunk. ``map`` returns only once every chunk is finished, so a per-chunk
    # callback fires 448 times in one instant -- noise that says nothing about progress.
    if progress:
        for i, locus in enumerate(ordered):
            progress(locus, i + 1, len(ordered))
    return out


def _map_unordered(fn, tasks: list, n_jobs: int) -> list:
    """``[fn(t) for t in tasks]``, load-balanced one task at a time. Raises if the pool cannot start."""
    from concurrent.futures import BrokenExecutor, ProcessPoolExecutor
    from multiprocessing import get_context

    workers = _workers(n_jobs, len(tasks))
    if workers < 2:
        return [fn(t) for t in tasks]
    try:
        with ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn")) as ex:
            return list(ex.map(fn, tasks, chunksize=1))
    except (BrokenExecutor, RuntimeError) as e:
        raise RuntimeError(
            f"could not start {workers} worker processes ({type(e).__name__}). Workers are spawned, "
            "so the calling module is re-imported: that works from an importable module and fails "
            "from `python -c` or a heredoc. Guard the call with `if __name__ == '__main__':`, or "
            "pass n_jobs=1. It is NOT falling back to one process.") from e


def read_pools(paths: dict) -> dict:
    """Memory-map the pool chunks and present one frame per locus, without copying them.

    ``rechunk=False`` is load-bearing: a rechunk would materialise the whole pool in this process,
    which is the per-worker copy this design exists to avoid.
    """
    return {loc: pl.concat([pl.read_ipc(p, memory_map=True) for p in ps], rechunk=False)
            for loc, ps in paths.items()}


#: Chunks each locus's pool is generated in. Fixed, so the pool is a function of the build and not
#: of the machine: a chunk count taken from the core count would make an artifact built on 64 cores
#: differ from the same build on 16. Sixty-four because ONE CHUNK IS THE FLOOR on wall time -- it
#: cannot be split further, however many cores the machine has. At 16 chunks IGH's 3.3M-sequence
#: pool is 207k per chunk, which is 49 s at the measured 4,200 seq/s; at 64 it is 12 s. Changing
#: this value changes the pools and so the artifact, which is why it is a recorded constant.
POOL_CHUNKS: int = 64


#: Per-worker state, set once by the pool initializer. One copy per worker **process**, never per
#: task: a worker that re-mapped the pools per task would pay that cost once per batch instead of
#: once per core, which is the mistake the initializer exists to prevent.
_W: dict = {}


def _init_worker(paths, plan, regime, seed, cols, featurise) -> None:
    import os

    from .cohort import WORKER_ENV

    os.environ[WORKER_ENV] = "1"
    _W.update(pools=read_pools(paths), plan=plan, regime=regime, seed=seed, cols=cols,
              featurise=featurise)


def _feat_batch(span) -> tuple:
    """Draw and featurise samples ``[a, b)`` in this worker; return the filled rows per locus."""
    a, b = span
    out = {loc: np.full((b - a, len(cols)), np.nan) for loc, cols in _W["cols"].items()}
    for j in range(a, b):
        raw = _W["featurise"](draw_one(_W["pools"], _W["plan"], _W["regime"], j, _W["seed"]))[0]
        for loc, cols in _W["cols"].items():
            out[loc][j - a] = [raw.get(c, np.nan) for c in cols]
    return a, out


def _ranges(n: int, parts: int) -> list:
    from .cohort import slices

    return [(a, b) for a, b in slices(n, parts) if b > a]


def _workers(n_jobs: int, cap: int) -> int:
    from ..cores import available_cores

    return max(1, min(n_jobs if n_jobs > 0 else available_cores(), cap))


def build_matrices(regime: str, *, sig: str, featurise, vocab: dict,
                   loci: "tuple[str, ...]" = L.LOCI, n_samples: int = 10_000,
                   size: "int | str" = DEFAULT_SIZE, seed: int = SEED, source: str = "olga",
                   n_jobs: int = 1, progress=None, tmp: "Path | None" = None,
                   depth_spread: "float | str | None" = None,
                   cohort: "str | None" = None) -> dict:
    """The corpus matrix, ``{locus: (array, columns)}``, built across ``n_jobs`` processes.

    Shared by both halves of the signature: ``vsig`` and ``rsig`` differ only in ``featurise`` and
    ``sig``, and the samples they see are identical -- same pools, same per-sample seeds, same
    depths -- which is what makes the two artifacts joinable on ``sample_id``.

    Parallel **across samples, never inside one**. A sample's featurisation is a pure function of
    that sample, so a worker needs no coordination; and because every worker draws its own samples
    from memory-mapped pools, the only thing that crosses a process boundary is one row of numbers
    per sample. Shipping drawn repertoires instead would be ~38 GB of pickling at the shipped size.

    Args:
        featurise: A **picklable** callable mapping ``{locus: frame}`` to ``(raw, channels)`` --
            a module-level function or a ``functools.partial`` over one, never a lambda.
        n_jobs: Worker processes; ``0`` means every core available. ``1`` runs in-process and is
            the default, so that a doctest or a 12-sample test does not spawn a pool.
        tmp: Scratch directory for the pool files. A private temporary directory by default,
            removed when the build finishes.

    Returns:
        ``{locus: (matrix, columns)}``, the matrix being ``(n_samples, p_locus)`` float64.
    """
    import shutil
    import tempfile

    ordered = sorted(loci, key=lambda x: L.LOCI.index(x))
    cols = {}
    for locus in (*ordered, L.NO_LOCUS):
        got = locus_matrix(vocab, sig, locus, n_samples)
        if got is not None:
            cols[locus] = got[1]
    mats = {loc: (np.full((n_samples, len(c)), np.nan), c) for loc, c in cols.items()}

    own_tmp = tmp is None
    tmp = Path(tempfile.mkdtemp(prefix="vdjtools-corpus-")) if own_tmp else Path(tmp)
    try:
        paths = build_pools(ordered, size=size, n_samples=n_samples, seed=seed, source=source,
                            tmp=tmp, n_jobs=n_jobs, progress=progress,
                            depth_spread=depth_spread, cohort=cohort)
        plan = draw_plan(ordered, n_samples=n_samples, size=size, seed=seed,
                         depth_spread=depth_spread, cohort=cohort)
        workers = _workers(n_jobs, n_samples)
        single_threaded_children()
        # Four spans per worker rather than one. The expensive per-worker state -- the mapped pools
        # and the germline vocabulary -- is loaded once by the initializer, not once per span, so
        # extra spans cost nothing and buy load balance plus progress that moves.
        spans = _ranges(n_samples, workers * 4 if workers > 1 else 1)
        args = (paths, plan, regime, seed, cols, featurise)
        if workers < 2:
            _init_worker(*args)
            results = [_feat_batch(sp) for sp in _progress(spans, n_samples, progress)]
        else:
            results = _pooled(spans, args, workers, n_samples, progress)
        for a, part in results:
            for loc, block in part.items():
                mats[loc][0][a:a + block.shape[0]] = block
    finally:
        if own_tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    return mats


def _progress(spans, n_samples, progress):
    for sp in spans:
        yield sp
        if progress:
            progress("featurise", sp[1], n_samples)


def single_threaded_children() -> None:
    """Make spawned workers take ONE kernel thread each, by setting the env they inherit.

    ``n_jobs`` processes each starting a kernel sized off the core count is ``cores x cores``
    threads. Measured here: seven pool workers on a 16-core box, each with a 16-thread polars,
    turned 15 s of receptor generation into 296 s. The variables have to be set in the parent,
    before any child exists -- a spawned child reads them while importing polars, which is strictly
    earlier than any initializer of ours can run. The parent's own kernel is already up, so this
    does not throttle it.
    """
    import os

    for var in ("POLARS_MAX_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "RAYON_NUM_THREADS"):
        os.environ.setdefault(var, "1")


def _pooled(spans, args, workers: int, n_samples: int, progress):
    """Run the spans in a spawned pool. A pool that cannot start **raises**, it does not degrade."""
    from concurrent.futures import BrokenExecutor, ProcessPoolExecutor
    from multiprocessing import get_context

    out, done, step = [], 0, max(n_samples // 20, 1)
    try:
        with ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn"),
                                 initializer=_init_worker, initargs=args) as ex:
            for a, part in ex.map(_feat_batch, spans):
                out.append((a, part))
                was, done = done, done + next(b - x for x, b in spans if x == a)
                # At most twenty lines, whatever the span count. One span per 4 workers means 288
                # spans on a 72-core box, and 288 progress lines in a log is not progress.
                if progress and done // step > was // step:
                    progress("featurise", done, n_samples)
    except (BrokenExecutor, RuntimeError) as e:
        raise RuntimeError(
            f"could not start {workers} worker processes ({type(e).__name__}). Workers are "
            "spawned, not forked -- polars cannot be combined with fork -- and a spawned worker "
            "re-imports the module that called this, which works from an importable module and "
            "fails from `python -c` or a heredoc. Guard the call with `if __name__ == "
            "'__main__':`, or pass n_jobs=1. It is NOT falling back to one process: that is what "
            "hid a 20x slowdown here before.") from e
    return out


def sample_stream(regime: str, *, loci: "tuple[str, ...]" = L.LOCI, n_samples: int = 10_000,
                  size: "int | str" = DEFAULT_SIZE, seed: int = SEED, source: str = "olga",
                  progress=None, tmp: "Path | None" = None, n_jobs: int = 1,
                  depth_spread: "float | str | None" = None, cohort: "str | None" = None):
    """Set up the pools in this process and return ``(ordered_loci, draw)``; ``draw(j)`` is sample j.

    The in-process path, for callers that want the repertoires themselves rather than a corpus
    matrix. ``draw`` may be called in any order -- see :func:`draw_one`.
    """
    import tempfile
    from functools import partial

    ordered = sorted(loci, key=lambda x: L.LOCI.index(x))
    tmp = Path(tempfile.mkdtemp(prefix="vdjtools-pools-")) if tmp is None else Path(tmp)
    paths = build_pools(ordered, size=size, n_samples=n_samples, seed=seed, source=source,
                        tmp=tmp, n_jobs=n_jobs, progress=progress, depth_spread=depth_spread,
                        cohort=cohort)
    pools = read_pools(paths)
    plan = draw_plan(ordered, n_samples=n_samples, size=size, seed=seed,
                     depth_spread=depth_spread, cohort=cohort)
    return ordered, partial(draw_one, pools, plan, regime, seed=seed)


def synthesize(corpus_name: str, *, loci: "tuple[str, ...]" = L.LOCI, n_samples: int = 10_000,
               size: "int | str | None" = None, seed: int = SEED,
               n_components: "int | float" = DEFAULT_COMPONENTS, mode: str = "features",
               winsor_p: float = 0.01, source: str = "olga", fit_corpus: bool = True,
               progress=None, n_jobs: int = 1, depth_spread: "float | str | None" = None):
    """Build a synthetic corpus, and fit it.

    Args:
        corpus_name: One of :data:`SYNTHETIC` -- ``"naive"`` and ``"memory"`` are the pure regimes
            at a nominal size; ``"synthetic-blood"`` and ``"synthetic-tissue"`` are the naive/memory
            mixture drawn across the named cohort's measured per-locus richness, read-depth and
            singleton-fraction bands (:data:`COHORT`), and are the two that describe a real bulk
            cohort rather than a regime.
        loci: Loci to build. All seven by default.
        n_samples: Repertoires in the corpus.
        size: Receptors per repertoire. ``None``, the default, is the corpus's own -- 10,000 for the
            pure regimes, the geometric centre of the cohort's measured richness band for a
            ``synthetic-*`` one. Also takes ``"n_eff"`` for the per-locus real medians in
            :data:`N_EFF`, ``"p05"`` / ``"p95"`` for the sweep endpoints, or a cohort name.
        seed: Base seed; every draw is a recorded offset from it.
        n_components: Components per locus, or a variance fraction.
        mode: Winsorization mode to fit at.
        winsor_p: Percentile for the fitted bounds. All of :data:`WINSOR_PS` are stored.
        source: Bundled model set -- ``"olga"``, ``"learned"`` or ``"arda"``.
        fit_corpus: ``False`` returns the raw samples instead of fitting, for scoring a held-out
            draw against a corpus fitted on another.
        progress: Optional ``callable(locus, done, total)``.
        n_jobs: Worker **processes** (not kernel threads). ``0`` means every available core; ``1``,
            the default, runs in-process. The result is bit-identical at every value.
        depth_spread: Multiplicative depth range each repertoire's size is drawn log-uniformly
            across, around ``size``. ``None`` is the corpus's own -- :data:`DEPTH_SPREAD` for a pure
            regime (2.4x-11.0x), the cohort's measured richness band for a ``synthetic-*`` one (43x
            on blood TRB, 259x on tissue IGH). Pass a number to override it: ``1000`` with
            ``size=3162`` spans 100 to 100,000 receptors per locus. The bounds, centre and per-PC
            scaling are all estimated from the draw, so a corpus describes only the depths it was
            drawn across.

    Returns:
        ``(corpus, mats)`` when fitting -- ``mats`` being ``{locus: (matrix, columns)}``, the corpus
        matrix itself -- else a list of ``{locus: frame}`` samples.

    The whole build is a deterministic function of ``(corpus_name, loci, n_samples, size, seed,
    source)``
    and the bundled models, so two machines at different thread counts must produce byte-identical
    artifacts. That is an acceptance criterion, not a hope: ``generate`` was once irreproducible
    *across processes* while being deterministic within one, because a marginal-table aggregation
    left its group order unspecified.
    """
    from functools import partial

    from . import features as FE

    regime, cohort, size, depth_spread = corpus_plan(corpus_name, size=size,
                                                    depth_spread=depth_spread)
    vocab = {loc: FE.gene_vocab(loc) for loc in loci}
    if not fit_corpus:
        _ordered, draw = sample_stream(regime, loci=loci, n_samples=n_samples, size=size,
                                       seed=seed, source=source, progress=progress, n_jobs=n_jobs,
                                       depth_spread=depth_spread, cohort=cohort)
        return [draw(j) for j in range(n_samples)]

    mats = build_matrices(regime, sig="vsig", featurise=partial(FE.raw_and_channels, vocab=vocab),
                          vocab=vocab, loci=loci, n_samples=n_samples, size=size, seed=seed,
                          source=source, n_jobs=n_jobs, progress=progress,
                          depth_spread=depth_spread, cohort=cohort)

    from .. import __version__
    corpus = fit_matrices(mats, vocab, sig="vsig", name=corpus_name, mode=mode,
                          n_components=n_components, winsor_p=winsor_p,
                          meta=corpus_meta(regime, cohort, loci, n_samples=n_samples,
                                           size=size, seed=seed, source=source,
                                           depth_spread=depth_spread,
                                           vdjtools_version=__version__))
    return corpus, mats


def corpus_meta(regime: str, cohort: "str | None", loci, *, n_samples: int,
                size: "int | str", seed: int, source: str, depth_spread, **extra) -> dict:
    """The manifest of a synthetic build -- everything needed to reproduce the draw, per locus.

    Shared by both halves, because a manifest that disagrees between ``vsig_<name>`` and
    ``rsig_<name>`` about the depths, the singleton fractions or the germline they were drawn from is
    a joinability claim nobody can check.
    """
    mix = {} if cohort is None else {loc: cohort_bands(cohort, loc) for loc in loci}
    return {
        # The corpus NAME is not here: `Corpus.name` is written beside this by `save`, and two
        # places for it is one place for them to disagree.
        "regime": regime, "cohort": cohort,
        "n_samples": n_samples, "seed": seed, "source": source, "size": size,
        "size_per_locus": {k: resolved_size(k, size) for k in loci},
        "zipf_a": ZIPF_A if regime in ("memory", "mixed") else None,
        "reads_per_clone": {"naive": 1, "memory": READS_PER_CLONE, "mixed": None}[regime],
        # The three measured p05/p50/p95 ladders a mixture corpus drew each sample from, per locus:
        # richness, reads per expanded clone, and the singleton fraction (the naive share). `None`
        # for a pure regime, which has none of them.
        "richness_band": {k: list(v[0]) for k, v in mix.items()} or None,
        "expanded_count_band": {k: list(v[1]) for k, v in mix.items()} or None,
        "singleton_frac_band": {k: list(v[2]) for k, v in mix.items()} or None,
        # ...and how they move together, which is what the rotation is fitted on.
        "rank_corr": {k: list(v[3]) for k, v in mix.items()} or None,
        "pool_factor": POOL_FACTOR,
        # The pools come out of these models; retraining one moves its locus. Gated on load.
        "models": model_fingerprint(loci, source),
        # The depths were DRAWN across this range; a corpus describes no other depths.
        "depth_spread": {k: depth_spread_of(k, depth_spread) for k in loci},
        "depth_spread_requested": depth_spread,
        "loci": list(loci),
        **extra,
    }


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
