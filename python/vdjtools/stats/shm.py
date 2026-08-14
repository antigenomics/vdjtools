"""Somatic hypermutation summary statistics for B-cell repertoires.

SHM is what makes a BCR repertoire a different object from a TCR repertoire. A T cell's V gene is
germline for life; a B cell's is rewritten in the germinal centre, and the *pattern* of that
rewriting is the readout — how much, on which isotypes, and how evenly spread.

The motivating application is **tertiary lymphoid structure** detection in tumours. A TLS is an
ectopic germinal centre, and a germinal centre leaves three joint marks on the local BCR
repertoire, none of which is sufficient alone:

1. **mutated V genes** — GC-experienced B cells carry SHM; naive infiltrate does not;
2. **class switching** — IgG/IgA over IgM/IgD, which happens in the GC;
3. **intraclonal diversification** — one lineage present at several mutation levels at once,
   because the GC is actively diversifying it rather than having imported a finished clone.

A tumour with passively infiltrating naive B cells scores high on none. A tumour with a
*resident* plasma-cell clone scores high on 1 and 2 but not 3. Only an active GC scores on all
three, which is why the switched-versus-unswitched SHM gap and the diversification breadth are
reported separately rather than folded into one index — the composite is the caller's to make,
and which of the three is doing the work is the interesting part.

What this module does **not** do is call clonal lineages. That needs junction clustering
(:mod:`vdjtools.overlap`) and is a different, heavier computation; ``shm_spectrum`` gives the
mutation-level distribution that a lineage-aware statistic would refine.

**Input requirement.** ``v_identity`` — the fraction of the V region matching germline — is not a
canonical column, because most repertoire formats do not carry it. Read it through the ``keep=``
argument that :func:`vdjtools.io.read.read_airr` and friends accept::

    df = read_airr("sample.tsv", keep=("v_identity",))
    shm_summary(df)

Without it every SHM field is ``nan`` rather than 0: a repertoire whose aligner did not report
identity is not an unmutated repertoire, and the two must not produce the same number.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..io.schema import C_CALL, weight_expr

V_IDENTITY = "v_identity"

#: Below this V-region identity a clonotype counts as mutated. 0.98 rather than a round 1.0
#: because sequencing and alignment error put a floor of a few tenths of a percent under any
#: real dataset, and a threshold at germline would call that floor somatic hypermutation.
MUTATED_MAX_IDENTITY = 0.98

#: Class-switched isotypes. IgD is *not* switched -- IgM/IgD co-expression is the naive state.
SWITCHED = ("IGHG", "IGHA", "IGHE")
UNSWITCHED = ("IGHM", "IGHD")


def _isotype(expr: pl.Expr) -> pl.Expr:
    """Constant-region call to its isotype class, e.g. ``IGHG1*01`` -> ``IGHG``.

    Subclass is dropped deliberately: IgG1 vs IgG3 is a real distinction but a much noisier one at
    the depths these statistics are computed from, and the GC signal is carried by the switch
    itself.
    """
    return expr.str.replace(r"^(IGH[MDGAE]).*$", "${1}").str.to_uppercase()


def _weighted(x: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    """Weighted mean and standard deviation, ``(nan, nan)`` when nothing is left to average."""
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not ok.any():
        return float("nan"), float("nan")
    x, w = x[ok], w[ok]
    w = w / w.sum()
    m = float(np.sum(w * x))
    return m, float(np.sqrt(max(np.sum(w * (x - m) ** 2), 0.0)))


def shm_summary(df: pl.DataFrame, *, weight: str = "freq",
                mutated_max_identity: float = MUTATED_MAX_IDENTITY) -> dict[str, float]:
    """Summarise somatic hypermutation over one repertoire.

    Args:
        df: Canonical clonotype frame carrying ``v_identity`` (see the module docstring) and,
            for the isotype fields, ``c_call``.
        weight: Clone weight, as :func:`vdjtools.io.schema.weight_expr` understands it. ``freq``
            weights by clone size, so the answer is "what fraction of the *cells* are mutated";
            ``unique`` weights every clonotype equally, answering "what fraction of the
            *lineages*". They differ a lot in a repertoire with one dominant plasma-cell clone,
            and which you want depends on the question.
        mutated_max_identity: Identity below which a clonotype counts as mutated.

    Returns:
        A flat ``{name: value}`` dict. Fields that the input cannot support are ``nan``, never 0:

        ``mean_identity``, ``sd_identity``
            Weighted V-region identity to germline.
        ``mean_shm``
            ``1 - mean_identity``, the mutation load, for readers who prefer it that way.
        ``frac_mutated``
            Weight fraction below ``mutated_max_identity``.
        ``frac_switched``
            Weight fraction of class-switched isotypes, over clonotypes with *any* isotype call.
        ``mean_identity_switched``, ``mean_identity_unswitched``
            V identity within each class.
        ``switch_shm_gap``
            ``mean_identity_unswitched - mean_identity_switched``: how much more mutated the
            switched compartment is. Positive under a working GC. This is the field that
            separates an active germinal centre from a resident plasma-cell clone, because a
            clone that switched and mutated elsewhere arrives with both arms already mutated.
        ``shm_entropy``
            Shannon entropy (nats) of the mutation-level distribution over 20 bins of identity.
            High when one repertoire holds several mutation levels at once -- intraclonal
            diversification -- and low when it is uniformly germline *or* uniformly mutated.
        ``n_identity``
            Clonotypes with a usable identity, so a caller can see what the rest rests on.
    """
    out: dict[str, float] = {
        "mean_identity": float("nan"), "sd_identity": float("nan"), "mean_shm": float("nan"),
        "frac_mutated": float("nan"), "frac_switched": float("nan"),
        "mean_identity_switched": float("nan"), "mean_identity_unswitched": float("nan"),
        "switch_shm_gap": float("nan"), "shm_entropy": float("nan"), "n_identity": 0.0,
    }
    if df is None or df.height == 0:
        return out

    # with_columns, not select: weight_expr("unique") is a scalar literal, and select
    # yields a ONE-row frame from it rather than broadcasting to the clonotypes.
    w_all = df.with_columns(weight_expr(weight).alias("w"))["w"].to_numpy().astype(float)

    if C_CALL in df.columns:
        iso = df.select(_isotype(pl.col(C_CALL)).alias("iso"))["iso"].to_numpy()
        called = np.array([i is not None and str(i).startswith("IGH") for i in iso])
        if called.any():
            sw = np.array([str(i) in SWITCHED for i in iso])
            # Denominator is the CALLED weight, not the total: `segment_usage` drops null c_call,
            # and dividing by the sample total instead would report a repertoire with 40%
            # uncalled isotype as 40% less switched than it is.
            out["frac_switched"] = float(w_all[sw & called].sum() / w_all[called].sum())
    else:
        iso = np.array([None] * df.height)
        called = np.zeros(df.height, bool)

    if V_IDENTITY not in df.columns:
        return out

    ident = df[V_IDENTITY].to_numpy().astype(float)
    ok = np.isfinite(ident)
    out["n_identity"] = float(ok.sum())
    if not ok.any():
        return out

    m, sd = _weighted(ident, w_all)
    out["mean_identity"], out["sd_identity"], out["mean_shm"] = m, sd, 1.0 - m
    wok = w_all[ok]
    out["frac_mutated"] = float(wok[ident[ok] < mutated_max_identity].sum() / wok.sum())

    sw = np.array([str(i) in SWITCHED for i in iso]) & called & ok
    un = np.array([str(i) in UNSWITCHED for i in iso]) & called & ok
    if sw.any():
        out["mean_identity_switched"] = _weighted(ident[sw], w_all[sw])[0]
    if un.any():
        out["mean_identity_unswitched"] = _weighted(ident[un], w_all[un])[0]
    if sw.any() and un.any():
        out["switch_shm_gap"] = out["mean_identity_unswitched"] - out["mean_identity_switched"]

    out["shm_entropy"] = shm_entropy(ident[ok], wok)
    return out


def shm_entropy(identity: np.ndarray, weight: np.ndarray, *, bins: int = 20) -> float:
    """Shannon entropy (nats) of the weighted mutation-level distribution.

    Binned on identity in ``[0, 1]`` rather than computed on the raw values: entropy of a
    continuous variable is not defined without a bin width, and stating the width is the honest
    version of choosing one. 20 bins puts each at 5% identity, which is coarser than the
    measurement noise and finer than the germline/mutated distinction.
    """
    h, _ = np.histogram(identity, bins=bins, range=(0.0, 1.0), weights=weight)
    p = h[h > 0]
    if p.size == 0:
        return float("nan")
    p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


def shm_spectrum(df: pl.DataFrame, *, weight: str = "freq", bins: int = 20) -> pl.DataFrame:
    """The mutation-level distribution, as a table.

    The distribution is the thing to look at before any summary of it: a bimodal spectrum (a
    germline-like mode plus a mutated mode) is a mixture of naive infiltrate and GC output, and
    reports the same ``mean_identity`` as a unimodal repertoire sitting between the two.

    Returns:
        One row per bin: ``identity_low``, ``identity_high``, ``weight``, ``fraction``.
    """
    lo = np.linspace(0.0, 1.0, bins + 1)[:-1]
    hi = np.linspace(0.0, 1.0, bins + 1)[1:]
    empty = pl.DataFrame({"identity_low": lo, "identity_high": hi,
                          "weight": np.zeros(bins), "fraction": np.full(bins, np.nan)})
    if df is None or df.height == 0 or V_IDENTITY not in df.columns:
        return empty
    ident = df[V_IDENTITY].to_numpy().astype(float)
    w = df.with_columns(weight_expr(weight).alias("w"))["w"].to_numpy().astype(float)
    ok = np.isfinite(ident) & np.isfinite(w)
    if not ok.any():
        return empty
    h, _ = np.histogram(ident[ok], bins=bins, range=(0.0, 1.0), weights=w[ok])
    total = h.sum()
    return pl.DataFrame({"identity_low": lo, "identity_high": hi, "weight": h,
                         "fraction": h / total if total > 0 else np.full(bins, np.nan)})


def _demo() -> None:
    """Two synthetic repertoires that a mean identity alone cannot tell apart."""
    # A: an active GC -- switched cells mutated, unswitched germline, several levels at once.
    gc = pl.DataFrame({
        "junction_aa": [f"CAR{i}W" for i in range(6)],
        "v_call": ["IGHV3-23*01"] * 6, "j_call": ["IGHJ4*02"] * 6,
        "c_call": ["IGHM*01", "IGHD*01", "IGHG1*01", "IGHG1*01", "IGHA1*01", "IGHG3*01"],
        "duplicate_count": [10, 10, 10, 10, 10, 10],
        "frequency": [1 / 6] * 6,
        V_IDENTITY: [1.00, 0.99, 0.92, 0.87, 0.83, 0.95],
    })
    # B: naive infiltrate -- same mean-ish identity story is NOT reproduced; everything germline.
    naive = gc.with_columns(pl.Series(V_IDENTITY, [1.0, 1.0, 0.995, 0.99, 1.0, 0.995]))

    a, b = shm_summary(gc), shm_summary(naive)
    assert a["frac_mutated"] > 0.6 > b["frac_mutated"], (a["frac_mutated"], b["frac_mutated"])
    assert a["switch_shm_gap"] > 0.05, a["switch_shm_gap"]      # switched arm is more mutated
    assert abs(b["switch_shm_gap"]) < 0.02, b["switch_shm_gap"]
    assert a["shm_entropy"] > b["shm_entropy"], (a["shm_entropy"], b["shm_entropy"])
    assert abs(a["frac_switched"] - 4 / 6) < 1e-9, a["frac_switched"]

    # No v_identity at all is nan, NOT an unmutated repertoire.
    bare = gc.drop(V_IDENTITY)
    c = shm_summary(bare)
    assert np.isnan(c["mean_identity"]) and np.isnan(c["frac_mutated"])
    assert abs(c["frac_switched"] - 4 / 6) < 1e-9        # isotype still works without identity

    sp = shm_spectrum(gc)
    assert sp.height == 20 and abs(sp["fraction"].sum() - 1.0) < 1e-9
    print("shm: GC vs naive separated on frac_mutated, switch gap and entropy; "
          "missing identity stays nan")


if __name__ == "__main__":
    _demo()
