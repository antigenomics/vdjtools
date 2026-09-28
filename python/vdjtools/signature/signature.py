"""``vsig`` and ``vsig_cohort``: raw features, a corpus rotation, and the emitted vector.

Three steps, each owned elsewhere and joined here::

    features.raw_and_channels(sample, vocab)   ->  raw features + channels
    corpus.apply(raw, chan, corpus)            ->  rotated columns + channels carried through

A corpus is **required**. There is no default, deliberately: a silently chosen rotation is the
mixed-units failure this rewrite exists to end, and two matrices rotated through different corpora
are not comparable no matter how alike their column names look.
"""
from __future__ import annotations

from collections.abc import Sequence

import functools
import warnings

import numpy as np
import polars as pl

from . import corpus as C
from . import features as F
from . import layout as L
from .cohort import parallel_rows, resolve_sample


def _locus_frames(sample) -> dict[str, pl.DataFrame]:
    """``{locus: frame}`` from a dict, or from one frame carrying a ``locus`` column."""
    from ..io.schema import LOCUS

    sample = resolve_sample(sample)
    if isinstance(sample, dict):
        return {k: v for k, v in sample.items() if k in L.LOCI}
    if LOCUS not in sample.columns:
        raise ValueError(
            f"a single frame needs a {LOCUS!r} column to be split by locus; pass "
            "{locus: frame} instead, or add one with vdjtools.io.schema.add_locus")
    return {loc: part.drop(LOCUS) for (loc,), part in sample.group_by([LOCUS])
            if loc in L.LOCI}


def vsig(sample, corpus: C.Corpus, *, mode: "str | None" = None,
         winsor_p: "float | None" = None, n_components: "int | float | None" = None,
         cstar_target: "float | dict[str, float] | None" = None, weight: str = "log2p1",
         prefiltered: bool = False, on_duplicate: str = "error",
         named: "bool | Sequence[str]" = (),
         columns: "list[str] | None" = None) -> dict[str, float]:
    """The statistics half of the signature for one sample.

    Args:
        sample: ``{locus: frame}``, one frame with a ``locus`` column, or a zero-argument callable
            returning either.
        corpus: A fitted :class:`~vdjtools.signature.corpus.Corpus` for ``sig="vsig"``.
        mode: Winsorization mode override -- ``"features"``, ``"pcs"`` or ``"none"``. Defaults to
            whichever the corpus was fitted with.
        winsor_p: Which stored percentile to clamp at. Defaults to the fitted one.
        n_components: Truncate to this many components, or this cumulative variance fraction.
        cstar_target: Coverage level the Hill numbers are read at. ``None`` uses each locus's own
            attained coverage, which makes them observed rather than standardised -- honest for one
            sample, and not comparable across samples. Use :func:`vsig_cohort` for a cohort, which
            resolves a shared level.
        weight: Clone-size weight; a key of :data:`~vdjtools.signature.features.WEIGHTS`.
        prefiltered: Report ``qc:*:nonstd_aa_frac`` as ``nan`` because the caller already filtered.
        on_duplicate: ``"error"`` or ``"sum"``, for a frame repeating an amino-acid clonotype key.
        named: Also return the reportable raw blocks -- ``True`` for all of them, or a sequence
            such as ``("div", "depth", "clon")``. ``()``, the default, emits exactly the rotated
            columns and channels. Values carry their declared transform (``log10``, ``clr``,
            ``logit``, ...), which :func:`~vdjtools.signature.layout.channel_table` reports.

            These are the diversity, depth, clone-size, junction-length, isotype, SHM and
            cross-locus yield numbers. They are computed either way, because the rotation is
            fitted on them; without this argument there is no supported route to reading them
            back, and a study that cannot report its own diversity floor has no floor. See
            :func:`~vdjtools.signature.layout.channel_table` for the full list and its units.
        columns: Restrict the **output** to these columns, in layout order.

            Unlike the system this replaces, this does not skip work: the rotation for a locus is
            fitted jointly over every raw group at that locus, so every group must be computed
            before any of that locus's components exist. Declining a locus entirely does skip it.

    Returns:
        ``{column: value}`` in layout order, holes as ``nan``.
    """
    if corpus.sig != "vsig":
        raise ValueError(f"corpus {corpus.name!r} is for {corpus.sig!r}, not 'vsig'")
    raw, chan = F.raw_and_channels(_locus_frames(sample), corpus.vocab, cstar_target=cstar_target,
                                  weight=weight, prefiltered=prefiltered,
                                  on_duplicate=on_duplicate)
    out = C.apply(raw, chan, corpus, mode=mode, winsor_p=winsor_p, n_components=n_components,
                  named=named)
    if columns is None:
        return out
    want = [c for c in corpus.columns(n_components, named=named) if c in set(columns)]
    return {c: out[c] for c in want}


def attained_coverage(sample, *, on_duplicate: str = "error") -> dict[str, float]:
    """Per-locus attained Chao coverage for one sample -- the cheap first pass.

    Only the count column is touched, so this is a small fraction of the cost of the full feature
    set and is what lets a cohort agree on one coverage level without anybody choosing a constant.
    """
    out: dict[str, float] = {}
    for locus, df in _locus_frames(sample).items():
        clean, _ = F.sanitise(df, on_duplicate=on_duplicate)
        if clean.height:
            out[locus] = F.coverage_of(clean)
    return out


def _one(item, corpus, kw):
    """One sample's row. Module-level so a spawned worker can unpickle it."""
    sid, sample = item
    return {"sample_id": sid, **vsig(resolve_sample(sample), corpus, **kw)}


def _cov(item, on_duplicate):
    sid, sample = item
    return {"sample_id": sid, **attained_coverage(resolve_sample(sample),
                                                  on_duplicate=on_duplicate)}


def vsig_cohort(samples, corpus: C.Corpus, *, n_jobs: int = 1,
                cstar_target: "float | dict[str, float] | str | None" = "min",
                columns: "list[str] | None" = None, **kw) -> pl.DataFrame:
    """One row per sample, ``sample_id`` first.

    Args:
        samples: ``{sample_id: sample}`` or an iterable of pairs. A sample may be a zero-argument
            **picklable** callable, which defers the read into the worker and keeps peak memory at
            ``O(n_jobs)`` samples rather than the whole cohort.
        corpus: A fitted corpus.
        n_jobs: Worker processes. ``1`` stays in-process, ``0`` uses every available core.
        cstar_target: ``"min"`` (default) reads every sample's attained coverage first and
            standardises the whole cohort to the per-locus minimum, so the diversity columns are
            comparable across samples and the level is a property of *this* cohort rather than of
            somebody's corpus. That is a cheap extra pass over the count column -- but it does
            resolve deferred samples twice, so pass an explicit float or dict for a one-pass run.
            ``None`` reads each sample at its own coverage: one pass, not comparable across samples.
        columns: Restrict the output columns.
        **kw: Forwarded to :func:`vsig` -- including ``named=``, which adds the reportable raw
            blocks in their own units.

    Returns:
        A frame whose columns are ``sample_id`` then the corpus's signature columns.
    """
    items = list(samples.items() if isinstance(samples, dict) else samples)
    if cstar_target == "min":
        cov = parallel_rows(items, functools.partial(_cov, on_duplicate=kw.get("on_duplicate",
                                                                              "error")), n_jobs)
        per_locus: dict[str, float] = {}
        for row in cov:
            for loc, v in row.items():
                if loc != "sample_id" and np.isfinite(v):
                    per_locus[loc] = min(per_locus.get(loc, v), float(v))
        cstar_target = per_locus or None
        if cstar_target is None:
            warnings.warn("no sample attained a measurable coverage at any locus, so no shared "
                          "level could be established; the diversity columns will be holes.",
                          stacklevel=2)

    rows = parallel_rows(items, functools.partial(_one, corpus=corpus,
                                                 kw={**kw, "cstar_target": cstar_target,
                                                     "columns": columns}), n_jobs)
    want = ["sample_id", *(columns if columns is not None
                           else corpus.columns(kw.get("n_components"),
                                               named=kw.get("named", ())))]
    have = [c for c in want if any(c in r for r in rows)]
    return pl.DataFrame(rows).select(have)


def _demo() -> None:
    """Self-check: fit a small corpus from synthetic samples, then score them through it."""
    from .corpus import synthesize

    corpus, _ = synthesize("naive", loci=("TRB", "TRG"), n_samples=40, size=120, seed=7,
                           n_components=4)
    samples = {f"S{i}": s for i, s in enumerate(
        synthesize("naive", loci=("TRB", "TRG"), n_samples=3, size=120, seed=99,
                   n_components=4, fit_corpus=False))}
    d = vsig_cohort(samples, corpus)
    assert d.height == 3 and d.columns[0] == "sample_id"
    assert set(d.columns[1:]) == set(corpus.columns())

    # coverage is emitted for every present locus, and is a real number in [0, 1]
    for loc in ("TRB", "TRG"):
        v = d[f"vsig:cov:{loc}:cstar"].to_numpy()
        assert np.isfinite(v).all() and ((v >= 0) & (v <= 1)).all(), loc
    # an absent locus is a declared hole, not a zero
    assert (d["vsig:mask:IGH:present"].to_numpy() == 0).all()

    # a column subset is a subset, in layout order, and changes nothing else
    sub = ["vsig:pc:TRB:PC01", "vsig:cov:TRB:cstar"]
    d2 = vsig_cohort(samples, corpus, columns=sub)
    assert d2.columns == ["sample_id", *sub], d2.columns
    for c in sub:
        assert np.allclose(d[c].to_numpy(), d2[c].to_numpy()), c

    print(f"signature OK  {d.width - 1} columns, k={corpus.k}")


if __name__ == "__main__":
    _demo()
