"""Turn one sample into the ``vsig`` half of the signature — the named vector itself.

:mod:`vdjtools.signature.blocks` computes families of features; this puts them in the order
:mod:`vdjtools.signature.layout` promises, under the names it promises, with a hole wherever a
number could not honestly be produced. The result is positional: column *i* means the same thing
for every sample anyone ever computes, which is the whole point.

The geometry half (``rsig``) is assembled by ``mir.signature``, which depends on this package
rather than the other way round, and the two concatenate on ``sample_id``.

**Holes are load-bearing.** An absent locus, a locus too shallow to support a coverage-
standardised estimate, an isotype column on a sample with no constant-gene calls — each yields
``nan`` plus a ``mask`` column saying so, never a zero. Zero is a measurement; a hole is not, and
a model that cannot tell them apart will read "this donor has no IgG" and "we did not sequence
IgG here" as the same statement.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..io.schema import COUNT, LOCUS, add_locus, column_names
from . import blocks as B
from . import layout as L

#: Coverage level the Hill numbers are standardised to, per locus, until the reference artifact
#: supplies measured ones. Deliberately low: attained Good-Turing coverage on real repertoires
#: runs 0.24-0.58, so a textbook 0.95 would put every sample into extrapolation, where the
#: estimator inflates diversity roughly tenfold. See SIGNATURE.md.
DEFAULT_CSTAR = 0.20


def _locus_frames(sample) -> dict[str, pl.DataFrame]:
    """Accept either ``{locus: frame}`` or one frame carrying a ``locus`` column."""
    if isinstance(sample, dict):
        return {k: v for k, v in sample.items() if v is not None and v.height}
    df = sample if LOCUS in column_names(sample) else add_locus(sample)
    return {k[0] if isinstance(k, tuple) else k: v
            for k, v in df.partition_by(LOCUS, as_dict=True).items()}


def vsig(sample, *, tier: str = "standard", cstar: float | dict[str, float] = DEFAULT_CSTAR,
         weight: str = "log2p1", pgen_q05: dict[str, float] | None = None,
         threads: int = 0) -> dict[str, float]:
    """The ``vsig`` half of one sample's signature, as ``{column_name: value}``.

    Args:
        sample: ``{locus: clonotype frame}``, or a single frame with a ``locus`` column.
        tier: ``"core"``, ``"standard"`` or ``"full"``.
        cstar: Coverage level for the standardised Hill numbers; a scalar or ``{locus: level}``.
        weight: Clone-size weight ``g`` (see :func:`~vdjtools.signature.blocks.work_frame`).
        pgen_q05: Per-locus frozen 5th-percentile ``log10 Pgen`` for ``pgen:*:frac_atypical``;
            that column stays ``nan`` without it, since "atypical" is meaningless without a
            reference to be atypical against.
        threads: Worker threads for the Pgen batch; 0 = auto.
            Off by default so ``tier="full"`` still runs; on when a caller needs the guarantee
            that every declared column was actually computed.

    Returns:
        Every column :func:`vdjtools.signature.layout.columns` lists for ``tier`` and ``"vsig"``,
        in that order, with ``nan`` where the sample could not support one.

    Raises:
        ValueError: If ``tier`` is unknown.
    """
    want = L.columns(tier, "vsig")
    out = dict.fromkeys(want, np.nan)
    frames = _locus_frames(sample)
    full = tier == "full"
    std = tier in ("standard", "full")
    reads: dict[str, float] = {}

    for locus in L.LOCI:
        raw = frames.get(locus)
        present = raw is not None and raw.height > 0
        out[f"vsig:mask:{locus}:present"] = float(present)
        if not present:
            out[f"vsig:mask:{locus}:estimable"] = 0.0
            continue

        clean, nonstd = B.sanitise(raw)
        if clean.height == 0:
            out[f"vsig:mask:{locus}:estimable"] = 0.0
            out[f"vsig:qc:{locus}:nonstd_aa_frac"] = B.qc_block(
                raw, clean, locus, nonstd)["nonstd_aa_frac"]
            continue

        work = B.work_frame(clean, weight)
        reads[locus] = float(clean[COUNT].sum())
        stats = _stats(clean)
        level = cstar[locus] if isinstance(cstar, dict) else cstar

        _put(out, f"vsig:qc:{locus}", B.qc_block(raw, work, locus, nonstd))
        _put(out, f"vsig:depth:{locus}", B.depth_block(clean, stats))
        _put(out, f"vsig:clon:{locus}", B.clon_block(stats))
        _put(out, f"vsig:len:{locus}", B.len_block(work, tier_standard=std))

        div = B.div_block(clean, level, tier_full=full)
        _put(out, f"vsig:div:{locus}", div)
        out[f"vsig:mask:{locus}:estimable"] = float(np.isfinite(div.get("1D_c", np.nan)))

        if std:
            _put(out, f"vsig:pgen:{locus}",
                 B.pgen_block(work, locus, q05=(pgen_q05 or {}).get(locus), threads=threads))
        if full:
            _put(out, f"vsig:aa:{locus}", B.aa_block(work))
            _put(out, f"vsig:pchem:{locus}", B.pchem_block(work))

        if locus == "IGH":
            _put(out, "vsig:iso:IGH", B.iso_block(work, tier_full=full))
            _put(out, "vsig:shm:IGH", B.shm_block(work))
            out["vsig:mask:IGH:c_call"] = float(
                "c_call" in work.columns and work["c_call"].null_count() < work.height)
            out["vsig:mask:IGH:shm"] = float("v_identity" in work.columns)

    _put(out, "vsig:pair:-", B.pair_block(reads))
    out["vsig:qc:-:n_loci_present"] = float(len(reads))
    return {k: out[k] for k in want}


def _put(out: dict, prefix: str, values: dict) -> None:
    """Write ``{feature: value}`` under ``prefix``, ignoring anything the tier does not want."""
    for k, v in values.items():
        key = f"{prefix}:{k}"
        if key in out:
            out[key] = float(v)


def _stats(df: pl.DataFrame) -> dict[str, float]:
    """The clone-size summary the blocks consume.

    Deliberately computed here rather than imported from ``mir.repertoire.sample_statistics``:
    vdjtools cannot depend on mirpy (the dependency runs the other way), and these six numbers
    are a one-line reduction of the count vector.
    """
    a = df[COUNT].to_numpy()
    return {"n_reads": float(a.sum()), "richness": float(a.size),
            "f1": float((a == 1).sum()), "f2": float((a == 2).sum()),
            "f3plus": float((a >= 3).sum()),
            "top_clone_fraction": float(a.max() / a.sum()) if a.sum() else 0.0}


def vsig_cohort(samples, *, tier: str = "standard", **kw):
    """Assemble a whole cohort into one frame: ``sample_id`` plus the ``vsig`` columns.

    Args:
        samples: ``{sample_id: sample}`` or an iterable of ``(sample_id, sample)``.
        tier: Passed to :func:`vsig`.
        **kw: Passed to :func:`vsig`.

    Returns:
        A ``pl.DataFrame``, one row per sample, columns in layout order.
    """
    items = samples.items() if isinstance(samples, dict) else samples
    rows = [{"sample_id": sid, **vsig(s, tier=tier, **kw)} for sid, s in items]
    if not rows:
        return pl.DataFrame(schema={"sample_id": pl.Utf8})
    return pl.DataFrame(rows).select(["sample_id", *L.columns(tier, "vsig")])


def _demo() -> None:
    """Self-check: the contract holds, and an absent locus is a hole rather than a zero."""
    rng = np.random.default_rng(0)
    aa = list("ACDEFGHIKLMNPQRSTVWY")

    def frame(n, v, j, c=None):
        return pl.DataFrame({
            "v_call": [v] * n, "j_call": [j] * n, "c_call": [c] * n,
            "junction_aa": ["C" + "".join(rng.choice(aa, 12)) + "F" for _ in range(n)],
            "duplicate_count": np.ceil(rng.zipf(1.5, n).clip(1, 900)).astype(int).tolist(),
        })

    sample = {"TRB": frame(2000, "TRBV20-1", "TRBJ2-2"),
              "IGH": frame(800, "IGHV1-2", "IGHJ4", "IGHM")}

    for tier in L.TIERS:
        v = vsig(sample, tier=tier)
        assert list(v) == L.columns(tier, "vsig"), f"{tier} is not in layout order"

    v = vsig(sample)
    assert v["vsig:mask:TRB:present"] == 1.0 and v["vsig:mask:IGH:present"] == 1.0
    assert v["vsig:mask:TRA:present"] == 0.0
    assert np.isnan(v["vsig:div:TRA:1D_c"]), "an absent locus must be a hole, not a zero"
    assert np.isfinite(v["vsig:div:TRB:1D_c"]), "TRB should have supported a diversity estimate"
    assert v["vsig:qc:-:n_loci_present"] == 2.0
    assert np.isfinite(v["vsig:pair:-:log_IGH_TRB"])
    assert np.isnan(v["vsig:shm:IGH:mean_v_identity"]), "shm masks out without v_identity"

    # tiers are prefixes of one another, valued identically
    core, std = vsig(sample, tier="core"), vsig(sample, tier="standard")
    for k in core:
        assert (core[k] == std[k]) or (np.isnan(core[k]) and np.isnan(std[k])), k

    # a TRB-only collaborator still gets a valid, full-width vector
    solo = vsig({"TRB": sample["TRB"]})
    assert len(solo) == len(L.columns("standard", "vsig"))
    assert solo["vsig:mask:IGH:present"] == 0.0

    frame_out = vsig_cohort({"a": sample, "b": solo and {"TRB": sample["TRB"]}})
    assert frame_out.height == 2 and frame_out.columns[0] == "sample_id"
    n_finite = int(np.isfinite(np.array(list(v.values()))).sum())
    print(f"assemble OK — {len(v)} vsig columns at standard, {n_finite} finite on a 2-locus "
          f"sample, {len(v) - n_finite} holes")


if __name__ == "__main__":
    _demo()
