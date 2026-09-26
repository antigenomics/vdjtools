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

import warnings

import numpy as np
import polars as pl

from ..io.schema import COUNT, LOCUS, add_locus, column_names
from . import blocks as B
from . import layout as L

#: Coverage level the Hill numbers are standardised to, per locus, until the reference artifact
#: supplies measured ones. Deliberately low: attained Good-Turing coverage on real repertoires
#: runs 0.24-0.58, so a textbook 0.95 would put every sample into extrapolation, where the
#: estimator inflates diversity roughly tenfold.
#:
#: This is a **fallback, not a recommendation** -- it is one number for all seven loci, and the
#: measured values are neither uniform nor assay-independent (TRB attains 0.408 in amplicon data
#: and 0.126 in bulk blood RNA-seq). Pass measured constants via ``cstar=`` whenever a reference
#: artifact supplies them; ``mir.signature.signature()`` does this for you. See ``docs/signature.rst``.
#:
#: And because it is a fallback, it is **declared**: every sample carries
#: ``vsig:qc:-:cstar_fallback_frac``, the share of its present loci whose coverage level came
#: from here rather than from a measurement. A ``vsig:div`` block computed at a level nobody
#: established is fully populated and entirely plausible, so the only way a collaborator can
#: tell is if the vector says so. ``cstar=None`` refuses the fallback outright and emits holes.
DEFAULT_CSTAR = 0.20


def _locus_frames(sample) -> dict[str, pl.DataFrame]:
    """Accept either ``{locus: frame}`` or one frame carrying a ``locus`` column."""
    if isinstance(sample, dict):
        return {k: v for k, v in sample.items() if v is not None and v.height}
    df = sample if LOCUS in column_names(sample) else add_locus(sample)
    return {k[0] if isinstance(k, tuple) else k: v
            for k, v in df.partition_by(LOCUS, as_dict=True).items()}


def vsig(sample, *, tier: str = "standard",
         cstar: float | dict[str, float] | None = DEFAULT_CSTAR,
         weight: str = "log2p1", pgen_q05: dict[str, float] | None = None,
         pgen_n_max: int = 2000, columns: list[str] | None = None,
         kmer_spaces: dict | None = None, threads: int = 0,
         prefiltered: bool = False, on_duplicate: str = "error") -> dict[str, float]:
    """The ``vsig`` half of one sample's signature, as ``{column_name: value}``.

    Args:
        sample: ``{locus: clonotype frame}``, or a single frame with a ``locus`` column.
        tier: ``"core"``, ``"standard"`` or ``"full"``.
        cstar: Coverage level the Hill numbers are standardised to. ``{locus: level}`` is the
            measured form and the one to use -- a reference artifact supplies it, and
            ``mir.signature.signature()`` passes it for you. A **scalar** applies one number to
            all seven loci; that is a fallback, not a measurement, and the emitted
            ``vsig:qc:-:cstar_fallback_frac`` says what share of the present loci got one.
            ``None``, or a dict missing a locus, establishes no level for that locus: its
            ``vsig:div:*`` columns are holes and ``vsig:mask:*:estimable`` is 0, rather than a
            confident number at a level nobody chose. A plausible number in the wrong units
            cannot be detected by the caller; a hole and a fraction can.
        weight: Clone-size weight ``g`` (see :func:`~vdjtools.signature.blocks.work_frame`).
        pgen_q05: Per-locus frozen 5th-percentile ``log10 Pgen`` for ``pgen:*:frac_atypical``;
            that column stays ``nan`` without it, since "atypical" is meaningless without a
            reference to be atypical against.
        pgen_n_max: How many junctions per locus the Pgen block measures, subsampled
            deterministically by :func:`~vdjtools.signature.blocks.pgen_junctions`. Cost is
            linear in it and the Pgen block dominates this half, so this is the one knob that
            trades runtime against the noise on three columns. **Lower it only together with the
            ``pgen_q05`` reference's own draw** if you care about ``frac_atypical`` -- the mean
            and sd are unbiased at any ``n``, but a reference drawn at a different ``n_max`` is
            still a different draw. Default matches the shipped references.
        columns: Compute only these columns (intersected with the tier, kept in layout order);
            unlisted ones are dropped from the result **and their block is not run**. This is
            how you decline the Pgen block, which is ~96% of this half's cost on a seven-locus
            sample -- see :func:`~vdjtools.signature.blocks.pgen_block` for why it is expensive
            and why that is not a bug. ``mask`` and the per-locus ``div`` feeding
            ``mask:*:estimable`` are always computed; they are cheap and everything else reads
            them.
        kmer_spaces: Per-locus frozen :class:`~vdjtools.features.kmer_space.KmerSpace`. The
            ``kmer`` block is emitted only for loci present here AND only if the block has been
            registered (see :func:`vdjtools.signature.kmer.register_kmer`) -- the columns do not
            exist in the layout otherwise, so passing spaces without registering them is a no-op
            rather than a silent width change.
        threads: Worker threads for the Pgen batch; 0 = auto.
            Off by default so ``tier="full"`` still runs; on when a caller needs the guarantee
            that every declared column was actually computed.
        prefiltered: Say ``True`` when the caller already removed non-functional rearrangements
            upstream. ``vsig:qc:*:nonstd_aa_frac`` then reports ``nan`` -- a hole -- instead of a
            confident floor value. This matters because the column is measured by comparing the
            frame handed in against what survives ``sanitise``: on a pre-filtered input that
            comparison yields exactly 0, which ``logit`` turns into a large negative number
            meaning "this repertoire is exceptionally clean", when the truth is "we cannot tell".
            Every other column is unaffected either way -- they are all computed on the surviving
            rows, and pre-filtering does not change which rows survive.
        on_duplicate: What to do with a locus frame that has no ``junction_nt`` and repeats an
            amino-acid clonotype key -- ``"error"`` (default) or ``"sum"``. See
            :func:`vdjtools.io.schema.assert_resolvable`. The check runs *after* non-productive
            rows are dropped, so a locus of junk still masks out rather than raising.

    Returns:
        Every column :func:`vdjtools.signature.layout.columns` lists for ``tier`` and ``"vsig"``,
        in that order, with ``nan`` where the sample could not support one.

    Raises:
        ValueError: If ``tier`` is unknown, or a locus repeats an unresolvable amino-acid
            clonotype key and ``on_duplicate="error"``.
    """
    want = L.columns(tier, "vsig")
    if columns is not None:
        keep = set(columns)
        want = [c for c in want if c in keep]
    out = dict.fromkeys(want, np.nan)
    frames = _locus_frames(sample)
    full = tier == "full"
    std = tier in ("standard", "full")
    reads: dict[str, float] = {}
    measured_nonstd: dict[str, float] = {}
    fallback: dict[str, bool] = {}
    # A block is worth running only if some wanted column sits under it. Keying on the
    # "<sig>:<channel>:<locus>" prefix is exact -- every column name is that plus one feature.
    need = {c.rsplit(":", 1)[0] for c in want}

    for locus in L.LOCI:
        raw = frames.get(locus)
        present = raw is not None and raw.height > 0
        out[f"vsig:mask:{locus}:present"] = float(present)
        if not present:
            out[f"vsig:mask:{locus}:estimable"] = 0.0
            continue

        clean, nonstd = B.sanitise(raw, on_duplicate=on_duplicate)
        if clean.height == 0:
            out[f"vsig:mask:{locus}:estimable"] = 0.0
            out[f"vsig:qc:{locus}:nonstd_aa_frac"] = B.qc_block(
                raw, clean, locus, nonstd)["nonstd_aa_frac"]
            continue

        work = B.work_frame(clean, weight)
        reads[locus] = float(clean[COUNT].sum())
        stats = _stats(clean)
        # A dict entry is a level established FOR THIS LOCUS; a scalar is one number stretched
        # over seven loci that attain very different coverage (TRB reaches 0.408 in amplicon
        # data and 0.126 in bulk RNA-seq). Absent either, there is no level -- and no level
        # means a hole, not a guess.
        level = cstar.get(locus) if isinstance(cstar, dict) else cstar
        # "Got a level, and it was a borrowed one." A locus with no level at all is not a
        # fallback -- it is a hole, and mask:*:estimable already says so. Keeping the two
        # disjoint is what lets a reader add them up.
        fallback[locus] = level is not None and not isinstance(cstar, dict)

        measured_nonstd[locus] = float(nonstd)
        if f"vsig:qc:{locus}" in need:
            _put(out, f"vsig:qc:{locus}", B.qc_block(raw, work, locus, nonstd))
            if prefiltered:
                out[f"vsig:qc:{locus}:nonstd_aa_frac"] = np.nan
        if f"vsig:depth:{locus}" in need:
            _put(out, f"vsig:depth:{locus}", B.depth_block(clean, stats))
        if f"vsig:clon:{locus}" in need:
            _put(out, f"vsig:clon:{locus}", B.clon_block(stats))
        if f"vsig:len:{locus}" in need:
            _put(out, f"vsig:len:{locus}", B.len_block(work, tier_standard=std))

        if level is None:
            out[f"vsig:mask:{locus}:estimable"] = 0.0
        else:
            div = B.div_block(clean, level, tier_full=full)
            _put(out, f"vsig:div:{locus}", div)
            out[f"vsig:mask:{locus}:estimable"] = float(np.isfinite(div.get("1D_c", np.nan)))

        if std and f"vsig:pgen:{locus}" in need:
            _put(out, f"vsig:pgen:{locus}",
                 B.pgen_block(work, locus, q05=(pgen_q05 or {}).get(locus),
                              n_max=pgen_n_max, threads=threads))
        if full and f"vsig:aa:{locus}" in need:
            _put(out, f"vsig:aa:{locus}", B.aa_block(work))
        if full and f"vsig:pchem:{locus}" in need:
            _put(out, f"vsig:pchem:{locus}", B.pchem_block(work))

        if kmer_spaces and f"vsig:kmer:{locus}:PC01" in out:
            from .kmer import kmer_block

            _put(out, f"vsig:kmer:{locus}",
                 kmer_block(work, locus, kmer_spaces.get(locus), weight="freq"))

        if locus == "IGH":
            if "vsig:iso:IGH" in need:
                _put(out, "vsig:iso:IGH", B.iso_block(work, tier_full=full))
            if "vsig:shm:IGH" in need:
                _put(out, "vsig:shm:IGH", B.shm_block(work))
            out["vsig:mask:IGH:c_call"] = float(
                "c_call" in work.columns and work["c_call"].null_count() < work.height)
            out["vsig:mask:IGH:shm"] = float("v_identity" in work.columns)

    _put(out, "vsig:pair:-", B.pair_block(reads))
    out["vsig:qc:-:n_loci_present"] = float(len(reads))
    out["vsig:qc:-:cstar_fallback_frac"] = (
        float(sum(fallback.values())) / len(fallback) if fallback else np.nan)
    _warn_if_prefiltered(measured_nonstd, prefiltered)
    _warn_if_partial_cstar(cstar, reads)
    return {k: out[k] for k in want}


def _warn_if_partial_cstar(cstar, reads: dict[str, float]) -> None:
    """Warn when measured constants were supplied but do not cover every locus in the sample.

    A scalar fallback is a stated choice and travels in ``vsig:qc:-:cstar_fallback_frac``, so it
    does not warn. A *partial dict* is different: the caller asked for measured levels, and the
    loci the dict misses silently lose their whole diversity block. That is the right outcome --
    a hole beats a number at a level nobody established -- but it should be said once, because
    the usual cause is a reference fitted on another assay (the amplicon reference carries TRA
    and TRB only, and a B-cell sample scored against it loses five loci of ``vsig:div``).
    """
    if not isinstance(cstar, dict) or not reads:
        return
    missing = sorted(set(reads) - {k for k, v in cstar.items() if v is not None})
    if missing:
        warnings.warn(
            f"no coverage level was supplied for {', '.join(missing)}, so vsig:div:* is nan and "
            f"vsig:mask:*:estimable is 0 there. The reference covers "
            f"{', '.join(sorted(k for k in cstar if cstar[k] is not None)) or 'no locus'}. "
            "This is a hole on purpose -- standardising a Hill number to a level measured on a "
            "different locus is not a measurement -- but if you expected those loci, you are "
            "probably using a reference fitted on a different assay.", UserWarning, stacklevel=3)


def _warn_if_prefiltered(measured: dict[str, float], declared: bool) -> None:
    """Warn when every locus reports exactly zero non-functional weight.

    A single locus can genuinely reach zero -- light chains do it routinely (measured on 1,168
    blood samples from a clinical AIRR cohort: IGK 84 samples at exactly zero, 54 of them with 100+ clonotypes). Across
    EVERY locus at once it does not happen: 0 of those 1,168 samples, the minimum being 5 of 7
    loci carrying some. So the all-locus case is a reliable signal that somebody filtered
    upstream, and a per-locus threshold is not -- which is why this warns on the conjunction and
    why ``prefiltered`` is a parameter rather than something inferred.

    A warning, not an error: a synthetic or already-curated frame is a legitimate input, and the
    caller who knows that should pass ``prefiltered=True`` and get a hole instead of a floor.
    """
    if declared or len(measured) < 2:
        return
    if all(v == 0.0 for v in measured.values()):
        warnings.warn(
            f"every locus ({', '.join(sorted(measured))}) reports exactly zero non-functional "
            "weight, which does not occur in unfiltered repertoire data. The input was almost "
            "certainly filtered upstream. vsig:qc:*:nonstd_aa_frac is then measured against a "
            "frame that has nothing left to drop and reports a confident floor value rather than "
            "a hole -- pass prefiltered=True to get nan instead, or drop the qc block with a "
            "preset (classify, transfer, compact all do). Every other column is unaffected.",
            UserWarning, stacklevel=3)


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


def _one_vsig(item, tier, kw):
    """One sample's vsig row. Module-level so a worker process can unpickle it."""
    from .cohort import resolve_sample

    sid, s = item
    return {"sample_id": sid, **vsig(resolve_sample(s), tier=tier, **kw)}


def vsig_cohort(samples, *, tier: str = "standard", n_jobs: int = 1,
                columns: list[str] | None = None, **kw):
    """Assemble a whole cohort into one frame: ``sample_id`` plus the ``vsig`` columns.

    Args:
        samples: ``{sample_id: sample}`` or an iterable of ``(sample_id, sample)``. A value may
            be a zero-argument callable returning the sample, which defers the read into the
            worker and keeps peak memory at ``O(n_jobs)`` samples rather than the whole cohort.
        tier: Passed to :func:`vsig`.
        n_jobs: Worker processes; ``1`` (default) stays in-process, ``0`` uses every core this
            process may use. ``threads`` (the Pgen kernel's own threads) is the knob that
            actually scales this workload -- see :func:`vsig` and the cohort guide. The two are
            independent claims on the same cores, so do not max both.
        columns: Explicit subset, intersected with the ``vsig`` half. Unlike a select-and-drop,
            this **skips the work**: passing a list with no ``vsig:pgen:`` entry does not run
            the Pgen block at all.
        **kw: Passed to :func:`vsig`.

    Returns:
        A ``pl.DataFrame``, one row per sample, columns in layout order.
    """
    import functools

    from .cohort import parallel_rows

    want = L.columns(tier, "vsig")
    if columns is not None:
        keep = set(columns)
        want = [c for c in want if c in keep]
        kw = {**kw, "columns": want}
    items = list(samples.items() if isinstance(samples, dict) else samples)
    rows = parallel_rows(items, functools.partial(_one_vsig, tier=tier, kw=kw), n_jobs)
    if not rows:
        return pl.DataFrame(schema={"sample_id": pl.Utf8})
    return pl.DataFrame(rows).select(["sample_id", *want])


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
