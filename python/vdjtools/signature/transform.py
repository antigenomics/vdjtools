"""Variance-stabilising transforms, and the frozen reference rescaling.

Every signature column is a different kind of number — a read count, a Hill number, a share of
a composition, a coordinate of an embedding — and a downstream model should not have to know
which. Each feature therefore declares one transform (:mod:`vdjtools.signature.layout`), applied
where the feature is computed, and every column then passes through the same frozen reference
rescaling. What comes out is dimensionless, roughly symmetric, and on a common scale.

**Small denominators are the whole problem.** At the depths this signature has to work at — a
median of order a hundred clonotypes, a quarter of samples below ten — a "fraction" is very
often ``0/3``. A transform that sees only the ratio cannot tell ``0/3`` from ``0/500``, and maps
both to the same place, asserting a precision the data does not have. So every transform of a
proportion takes the **denominator** as well as the value, and shrinks toward the middle in
proportion to how little was counted. That is the Haldane–Anscombe correction for a logit, the
Anscombe correction for an arcsine, and a count-scaled multiplicative replacement for a CLR.

None of these are free parameters: ``1/2`` and ``3/8`` are the standard bias-minimising choices.
"""
from __future__ import annotations

import numpy as np

#: Values further than this many robust standard deviations from the reference are clipped, so a
#: single pathological sample cannot dominate a downstream model's scaling. Wide enough that a
#: genuine outlier stays an outlier.
DEFAULT_CLIP = 8.0

#: ``1.4826 * MAD`` estimates the standard deviation of a normal distribution.
_MAD_TO_SD = 1.4826


def log10(x, floor: float = 1.0):
    """``log10`` of a positive quantity, floored so an empty locus maps to 0 rather than -inf.

    Used for counts and for Hill numbers. The floor is 1 because both are counts of *things*:
    one clonotype, one read, one effective species. Zero of them is the same as the floor for
    every downstream purpose, and the presence mask already records that the locus was empty.
    """
    return np.log10(np.maximum(np.asarray(x, dtype=float), floor))


def log1p(x):
    """``log(1+x)`` for a non-negative quantity that genuinely reaches zero.

    Distinct from :func:`log10` in intent: this is for norms, dispersions and hit counts, where
    zero is a real, attainable value rather than an empty measurement.
    """
    return np.log1p(np.maximum(np.asarray(x, dtype=float), 0.0))


def logit(x, m):
    """Haldane–Anscombe logit of a proportion observed on a denominator ``m``.

    ``log((x·m + 1/2) / ((1−x)·m + 1/2))``. Adding half an observation to each side is the
    standard remedy for an empty cell; its side effect is exactly the behaviour wanted here —
    the transform of ``0`` depends on how many chances there were to see something::

        logit(0, m=3)   ->  -1.95      a fifth of the repertoire could hide here
        logit(0, m=500) ->  -6.91      it is really absent

    Args:
        x: Proportion(s) in ``[0, 1]``.
        m: Denominator(s) the proportion was observed on. Broadcasts against ``x``.

    Returns:
        The transformed value, finite for every input including exactly 0 and exactly 1.
    """
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    m = np.maximum(np.asarray(m, dtype=float), 0.0)
    return np.log((x * m + 0.5) / ((1.0 - x) * m + 0.5))


def arcsine(x, m):
    """Anscombe's variance-stabilising arcsine transform, ``asin(sqrt((x·m + 3/8)/(m + 3/4)))``.

    The right transform for a *sparse* composition — residue and gene-usage profiles, where most
    cells are structurally zero at shallow depth. Unlike a CLR it is defined at zero without any
    replacement step, and unlike a raw proportion its variance does not collapse near the
    boundary. Bounded in ``[0, π/2]``, so it cannot produce the heavy tail a log-ratio would.

    Args:
        x: Proportion(s) in ``[0, 1]``.
        m: Denominator(s) the proportion was observed on.
    """
    x = np.clip(np.asarray(x, dtype=float), 0.0, 1.0)
    m = np.maximum(np.asarray(m, dtype=float), 0.0)
    return np.arcsin(np.sqrt((x * m + 0.375) / (m + 0.75)))


def clr(parts, m=None, *, keys=None):
    """Centred log-ratio of a composition, with multiplicative zero replacement.

    A CLR is the natural coordinate for a composition whose *ratios* carry the meaning: it is
    the log of each part over the geometric mean of all of them, so it is invariant to the total
    and a difference between two coordinates is a log-ratio of two parts.

    Zeros are replaced **multiplicatively**, not additively: each zero part is set to
    ``delta = 0.5/m`` and the non-zero parts are scaled down by ``1 − n_zero·delta`` so the
    composition still closes. Adding a constant to every part instead — the common shortcut —
    distorts the ratios among the parts that *were* observed, which are the only ratios the
    coordinate system is about.

    **Compute the CLR over the whole composition, then select coordinates.** A CLR of a
    sub-composition is a different number from the corresponding coordinate of the full one, so
    a tier that ships four of six isotype parts must still divide by the six-part geometric
    mean. Doing it the other way would make the narrower tier stop being a slice of the wider
    one, which is the contract the layout exists to guarantee.

    Args:
        parts: Non-negative part values as a mapping ``{name: value}`` or an array. They need
            not sum to 1; they are closed here.
        m: Total count the composition was observed on, setting the replacement scale. Defaults
            to the sum of ``parts`` when they are counts.
        keys: Part order when ``parts`` is an array. Ignored for a mapping.

    Returns:
        ``{name: clr}`` when ``parts`` is a mapping (or ``keys`` is given), else an array. The
        coordinates sum to zero by construction, which is why the layout ships all but one of
        them: the last is exactly determined by the rest and would make any unregularised
        design matrix singular.

    Raises:
        ValueError: If fewer than two parts are given, or any part is negative.
    """
    names = list(parts) if isinstance(parts, dict) else keys
    v = np.asarray(list(parts.values()) if isinstance(parts, dict) else parts, dtype=float)
    if v.size < 2:
        raise ValueError(f"a composition needs at least 2 parts; got {v.size}")
    if np.any(v < 0):
        raise ValueError("composition parts must be non-negative")

    total = v.sum()
    if m is None:
        m = total
    m = max(float(m), 1.0)
    if total <= 0:                       # nothing observed: a flat composition, all ratios 1
        return dict.fromkeys(names, 0.0) if names else np.zeros_like(v)

    p = v / total
    zero = p <= 0
    if zero.any():
        delta = 0.5 / m
        p = np.where(zero, delta, p * (1.0 - zero.sum() * delta))

    out = np.log(p) - np.log(p).mean()
    return dict(zip(names, out)) if names else out


def reference_z(x, loc, scale, clip: float = DEFAULT_CLIP):
    """Rescale against the frozen reference: ``(x − loc) / scale``, clipped.

    ``loc`` and ``scale`` are the reference median and ``1.4826·MAD`` — robust, so a handful of
    pathological samples in the reference corpus cannot set the scale for everyone. They are
    **frozen**: a collaborator does not fit them, which is what makes their vector comparable to
    ours rather than merely internally consistent.

    A zero ``scale`` means the reference never saw this column vary; the value is passed through
    centred but unscaled rather than divided by zero.
    """
    x = np.asarray(x, dtype=float)
    scale = np.asarray(scale, dtype=float)
    safe = np.where(scale > 0, scale, 1.0)
    return np.clip((x - np.asarray(loc, dtype=float)) / safe, -clip, clip)


def robust_loc_scale(x, axis=0):
    """Reference ``(median, 1.4826·MAD)`` from observed values only.

    Non-finite entries are ignored rather than imputed. Computing the statistics *before* any
    imputation matters: filling first and then measuring deflates the scale in proportion to how
    sparse the column is, so the least-observed locus ends up with the largest apparent values
    and dominates every distance and every principal component.
    """
    x = np.asarray(x, dtype=float)
    with np.errstate(invalid="ignore"):
        loc = np.nanmedian(np.where(np.isfinite(x), x, np.nan), axis=axis)
        mad = np.nanmedian(np.abs(np.where(np.isfinite(x), x, np.nan) - loc), axis=axis)
    return np.nan_to_num(loc), np.nan_to_num(mad) * _MAD_TO_SD


def magnitude_scale(block, rms):
    """Rescale a whole block by one frozen scalar, with no centring.

    For a block whose *magnitude* is its meaning — the signed contrast to an unselected
    reference, where a sample with no immune deviation should land at the origin. Standardising
    such a block per column would force every coordinate to unit variance across samples, which
    makes a near-zero sample look exactly like a typical one and deletes the deficiency the
    block exists to carry.
    """
    return np.asarray(block, dtype=float) / (rms if rms > 0 else 1.0)


#: Emit-time transforms by the code a feature declares in the layout. Those taking a denominator
#: are called as ``f(x, m)``; the rest as ``f(x)``. ``clr`` is not here because it consumes a
#: whole composition at once rather than one value.
_UNARY = {"none": lambda x: np.asarray(x, dtype=float), "log10": log10, "log1p": log1p}
_BINARY = {"logit": logit, "arcsine": arcsine}


def apply(code: str, x, m=None):
    """Apply the transform named by ``code`` to one value or array.

    Args:
        code: A transform code from :data:`vdjtools.signature.layout.TRANSFORMS`.
        x: The value(s).
        m: The denominator, required by ``logit`` and ``arcsine``.

    Raises:
        ValueError: If ``code`` is unknown, if a denominator is missing for a transform that
            needs one, or if ``code`` is ``"clr"`` (call :func:`clr` with the whole composition).
    """
    if code in _UNARY:
        return _UNARY[code](x)
    if code in _BINARY:
        if m is None:
            raise ValueError(f"transform {code!r} needs a denominator m — a proportion without "
                             "its denominator cannot be stabilised (0/3 is not 0/500)")
        return _BINARY[code](x, m)
    if code == "clr":
        raise ValueError("clr consumes a whole composition; call clr(parts, m) directly")
    raise ValueError(f"unknown transform {code!r}")


def _demo() -> None:
    """Self-check: the properties each transform is here for."""
    # a proportion of nothing is not the same as a proportion of a lot
    assert logit(0.0, 3) > logit(0.0, 500), "logit ignores its denominator"
    assert np.isfinite([logit(0.0, 3), logit(1.0, 3)]).all(), "logit blew up at the boundary"
    assert logit(0.5, 10) == 0.0

    # monotone, and bounded where it claims to be
    p = np.linspace(0, 1, 50)
    for f, hi in ((lambda v: logit(v, 50), np.inf), (lambda v: arcsine(v, 50), np.pi / 2)):
        y = f(p)
        assert np.all(np.diff(y) > 0), "transform is not monotone"
        assert np.all(np.abs(y) <= hi + 1e-9)
    assert np.isfinite(arcsine(0.0, 0)), "arcsine undefined on an empty sample"

    # clr: closes, sums to zero, survives a structural zero, and respects ratios
    c = clr({"a": 4.0, "b": 2.0, "c": 0.0}, m=6)
    assert abs(sum(c.values())) < 1e-12, "clr coordinates do not sum to zero"
    assert abs((c["a"] - c["b"]) - np.log(2)) < 1e-9, "clr broke the a:b ratio"
    assert c["c"] < c["b"] < c["a"]
    assert all(v == 0.0 for v in clr({"a": 0.0, "b": 0.0}, m=10).values())
    # multiplicative replacement leaves the observed ratio alone; additive would not
    sparse = clr({"a": 3.0, "b": 1.0, "c": 0.0}, m=4)
    assert abs((sparse["a"] - sparse["b"]) - np.log(3)) < 1e-9

    # reference rescaling: clips, and does not divide by a zero scale
    assert reference_z(1e6, 0.0, 1.0) == DEFAULT_CLIP
    assert reference_z(3.0, 1.0, 0.0) == 2.0

    # robust statistics ignore holes rather than filling them
    loc, scale = robust_loc_scale(np.array([[1.0], [2.0], [3.0], [np.nan]]))
    assert loc[0] == 2.0 and scale[0] > 0

    # magnitude scaling leaves the origin at the origin
    assert np.allclose(magnitude_scale(np.zeros(4), 2.0), 0.0)

    assert apply("log10", 0.0) == 0.0 and apply("log1p", 0.0) == 0.0
    print("transform OK")


if __name__ == "__main__":
    _demo()
