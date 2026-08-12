"""Consistency audit of a :class:`~vdjtools.model.model.Model` against its own germline.

:func:`~vdjtools.model.schema.validate_tables` answers a narrow question — does every marginal
normalize? — and raises on the first offender. This module answers the wider one a model builder
actually has: **is this model internally coherent, and does it agree with the germline it claims
to be built on?** It returns a tidy issue frame instead of raising, so every problem in a model is
visible at once and the result is sortable, filterable and writable like any other table.

The checks exist because each one has silently produced a wrong answer at some point:

- a functional gene left at ``P=0`` makes Pgen exactly 0 for every clonotype using it, with no error;
- deletion mass past a germline's length is unreachable, so its probability is quietly lost;
- an allele in a marginal but not in the germline frame crashes the native packer, far from the cause;
- a model whose germline drifted from arda's scores a different sequence than the one you annotated.

Severity is the contract: ``error`` means the model is broken (it will crash, or score wrongly),
``warn`` means it is suspicious but usable, ``info`` is a note.
"""
from __future__ import annotations

import polars as pl

from .events import EventKind
from .model import Model
from .reference import _ISSUE_SCHEMA, _issue
from .schema import _allele_col, normalization_keys

#: Allele column -> the genomic frame it must resolve in. Both D events share one germline set.
_ALLELE_GENOMIC = {"v_allele": "genes_v", "j_allele": "genes_j",
                   "d_allele": "genes_d", "d2_allele": "genes_d"}

#: Deletion event -> the ``palindrome_max`` key(s) bounding its most-negative (palindromic) ndel.
_DELETION_ENDS = {"v_3_del": ("v_3",), "j_5_del": ("j_5",),
                  "d_del": ("d_5", "d_3"), "d2_del": ("d_5", "d_3")}

_REQUIRED_EVENTS = {
    "VJ": {"v_choice", "j_choice", "v_3_del", "j_5_del", "vj_ins", "vj_dinucl"},
    "VDJ": {"v_choice", "j_choice", "d_gene", "n_d", "v_3_del", "j_5_del", "d_del",
            "vd_ins", "dj_ins", "vd_dinucl", "dj_dinucl"},
}
#: Events that only make sense on a VDJ model (a VJ model carrying one is malformed).
_D_EVENTS = {"d_gene", "d2_gene", "n_d", "d_del", "d2_del", "vd_ins", "dj_ins", "dd_ins",
             "vd_dinucl", "dj_dinucl", "dd_dinucl"}
#: The extra events a tandem-D (``n_D = 2``) model must carry.
_TANDEM_EVENTS = {"d2_gene", "d2_del", "dd_ins", "dd_dinucl"}


def _group_sums(df: pl.DataFrame, keys: list[str]) -> pl.DataFrame:
    """``Σ p`` per normalization group; a single ``s`` row when the event is unconditioned."""
    if keys:
        return df.group_by(keys).agg(pl.col("p").sum().alias("s"))
    return df.select(pl.col("p").sum().alias("s"))


def _cut_lengths(model: Model, frame: str) -> dict[str, int]:
    """Allele -> palindrome-extended germline length, from a genomic frame."""
    g = model.genomic.get(frame)
    if g is None:
        return {}
    col = f"{frame.split('_')[1]}_allele"   # genes_d stores 'd_allele' even for the d2 events
    return {r[col]: len(r["cut_segment"] or "") for r in g.iter_rows(named=True)}


def _choice_mass(model: Model, event: str, allele_col: str) -> dict[str, float]:
    """Total mass each allele carries in a gene-choice table (summed over any parent)."""
    t = model.tables.get(event)
    if t is None or allele_col not in t.columns:
        return {}
    agg = t.group_by(allele_col).agg(pl.col("p").sum().alias("p"))
    return {r[allele_col]: r["p"] for r in agg.iter_rows(named=True)}


def check_model(model: Model, *, germline: str | pl.DataFrame = "auto", tol: float = 1e-5,
                raise_on: str | None = None) -> pl.DataFrame:
    """Audit a model's marginals against its manifest, its germline frames, and a reference library.

    Args:
        model: The model to check.
        germline: External germline to reconcile against. ``"auto"`` (default) tries
            :func:`~vdjtools.model.reference.load_germline` for the model's own locus/organism and
            **skips silently** if arda has none (a custom library, or a non-arda organism);
            ``"none"`` disables the reconciliation; a ``pl.DataFrame`` uses that library directly.
        tol: Absolute tolerance for the "sums to 1 (or 0)" normalization check.
        raise_on: If ``"error"`` (or ``"warn"``), raise ``ValueError`` when any issue at that
            severity or worse is found, instead of returning it.

    Returns:
        Tidy issue frame ``severity, check, event, segment, allele, detail, value`` — empty when
        the model is clean. ``severity`` is one of ``error``, ``warn``, ``info``.

    Raises:
        ValueError: Only when ``raise_on`` is set and a matching issue is present.

    Example:
        >>> issues = check_model(load_bundled("TRB", "learned"))
        >>> issues.filter(pl.col("severity") == "error").height
        0
    """
    rows: list[dict] = []
    rows += _check_event_set(model)
    rows += _check_tables(model, tol)
    rows += _check_genomic_cross_refs(model)
    rows += _check_germline(model, germline)

    out = pl.DataFrame(rows, schema=_ISSUE_SCHEMA)
    if raise_on:
        levels = ("error",) if raise_on == "error" else ("error", "warn")
        bad = out.filter(pl.col("severity").is_in(levels))
        if bad.height:
            detail = "; ".join(bad["detail"].to_list()[:5])
            raise ValueError(f"model has {bad.height} {raise_on}-level issue(s): {detail}")
    return out


def _check_event_set(model: Model) -> list[dict]:
    """The declared graph matches the chain type, and tandem-D is declared consistently."""
    rows: list[dict] = []
    events = set(model.manifest.events)
    chain = model.chain_type

    for missing in sorted(_REQUIRED_EVENTS[chain] - events):
        rows.append(_issue("error", "event_set",
                           f"{chain} model is missing the required event {missing!r}", event=missing))
    if chain == "VJ":
        for extra in sorted(events & _D_EVENTS):
            rows.append(_issue("error", "event_set",
                               f"VJ model declares the D-locus event {extra!r}", event=extra))
    elif "genes_d" not in model.genomic:
        rows.append(_issue("error", "event_set", "VDJ model has no 'genes_d' germline frame"))

    # A tandem model must declare the whole n_D=2 machinery, else the Pgen DP has no D-D path to
    # spend that mass on -- `pgen.prepare` raises on this; report it as a row rather than crash.
    nd = model.tables.get("n_d")
    if nd is not None and "n_d" in nd.columns:
        p_two = nd.filter(pl.col("n_d") == 2)["p"].sum()
        if p_two > 0 and not _TANDEM_EVENTS.issubset(events):
            missing = sorted(_TANDEM_EVENTS - events)
            rows.append(_issue("error", "event_set",
                               f"P(n_D=2) = {p_two:.4g} but the tandem events {missing} are absent",
                               event="n_d", value=float(p_two)))
    return rows


def _check_tables(model: Model, tol: float) -> list[dict]:
    """Per-table numeric health: normalization, probability range, support vs germline geometry."""
    rows: list[dict] = []
    pal = model.manifest.palindrome_max

    for name, event in model.manifest.events.items():
        df = model.tables.get(name)
        if df is None:
            rows.append(_issue("error", "missing_table", f"no marginal table for event {name!r}",
                               event=name))
            continue

        bad_p = df.filter(pl.col("p").is_null() | pl.col("p").is_nan()
                          | (pl.col("p") < 0) | (pl.col("p") > 1 + tol))
        if bad_p.height:
            rows.append(_issue("error", "probability_range",
                               f"{bad_p.height} row(s) with p null/NaN or outside [0, 1]",
                               event=name, value=float(bad_p.height)))

        # All offenders at once -- validate_tables raises on the first, which hides the shape of
        # the damage (one broken gene vs a whole table that never normalized).
        keys = normalization_keys(event)
        sums = _group_sums(df, keys)
        bad = sums.filter(((pl.col("s") - 1.0).abs() > tol) & (pl.col("s").abs() > tol))
        for r in bad.head(20).iter_rows(named=True):
            where = ", ".join(f"{k}={r[k]!r}" for k in keys) or "(whole table)"
            rows.append(_issue("error", "normalization",
                               f"group [{where}] sums to {r['s']:.6g}, not 1 or 0",
                               event=name, value=float(r["s"])))
        if bad.height > 20:
            rows.append(_issue("error", "normalization",
                               f"...and {bad.height - 20} further mis-normalized group(s)",
                               event=name, value=float(bad.height - 20)))
        zeros = sums.filter(pl.col("s").abs() <= tol).height
        if zeros:
            # Legal: an undefined conditional for a gene that is never used. Worth surfacing --
            # it is also what an absorbing EM state looks like.
            rows.append(_issue("info", "zero_conditional",
                               f"{zeros} conditional group(s) sum to 0 (undefined conditional)",
                               event=name, value=float(zeros)))

        if event.kind is EventKind.DINUCLEOTIDE:
            rows += _check_dinucl(name, df)
        elif event.kind in (EventKind.DELETION, EventKind.DELETION_2D):
            rows += _check_deletion(model, name, event, df, pal)
        elif event.kind is EventKind.INS_LENGTH:
            rows += _check_insertion(name, df)
    return rows


def _check_dinucl(name: str, df: pl.DataFrame) -> list[dict]:
    rows: list[dict] = []
    cells = df.select(["from_nt", "to_nt"]).unique().height
    if cells != 16 or df.height != 16:
        rows.append(_issue("error", "dinucleotide_complete",
                           f"{df.height} row(s) / {cells} distinct (from_nt, to_nt) cell(s); "
                           f"a Markov table needs exactly 16",
                           event=name, value=float(cells)))
    outside = df.filter((pl.col("from_nt") > 3) | (pl.col("to_nt") > 3)).height
    if outside:
        rows.append(_issue("error", "dinucleotide_complete",
                           f"{outside} row(s) with a nucleotide code outside 0-3 (A,C,G,T)",
                           event=name, value=float(outside)))
    return rows


def _check_deletion(model: Model, name, event, df: pl.DataFrame, pal: dict) -> list[dict]:
    """Deletion mass must land on a trim the Pgen DP can actually reach.

    The geometry, taken from ``pgen._v_options`` / ``pgen._d_middle`` rather than assumed: a trim
    consumes the **palindrome-extended** ``cut_segment``, and ``ndel = len(cut) - contributed -
    max_pal``. So for V/J the reachable band is ``-max_pal <= ndel <= len(cut) - 1 - max_pal``
    (the ``-1`` is the invariant that V and J each contribute at least one nt), and for D it is
    ``ndel5 + ndel3 <= len(cut) - max_pal_5 - max_pal_3`` — a D may legally be deleted away
    entirely. Mass outside those bands is unreachable and silently lost from every Pgen.
    """
    rows: list[dict] = []
    if not event.given:
        return rows
    parent_allele_col = _allele_col(model.manifest.events[event.given[0]])
    frame = _ALLELE_GENOMIC.get(parent_allele_col)
    cuts = _cut_lengths(model, frame) if frame else {}
    if not cuts or parent_allele_col not in df.columns:
        return rows

    ends = _DELETION_ENDS.get(name, ())
    missing_end = [e for e in ends if e not in pal]
    if missing_end:
        rows.append(_issue("error", "palindrome_max",
                           f"event {name!r} trims end(s) {missing_end} with no declared "
                           f"palindrome_max entry", event=name))
        return rows
    pal_total = sum(pal[e] for e in ends)

    # Alleles with NO germline at all are reported once each by `unscoreable_gene_mass`; letting
    # them through here would bury a real over-long deletion under a hundred rows of one known cause.
    total = (pl.col("ndel") if event.kind is EventKind.DELETION
             else pl.col("ndel5") + pl.col("ndel3"))
    used = df.with_columns(
        _total=total,
        _cut=pl.col(parent_allele_col).replace_strict(cuts, default=-1, return_dtype=pl.Int64),
    ).filter(pl.col("_cut") > 0)
    # Largest reachable trim: V/J must leave >=1 nt, a D need not.
    slack = 1 if event.kind is EventKind.DELETION else 0
    used = used.with_columns(_max=pl.col("_cut") - pal_total - slack)

    # Report the FRACTION of each allele's deletion mass that is unreachable, not one row per cell.
    # A shared deletion-bin grid across alleles of different lengths (OLGA's array layout, and the
    # placeholder tables from `from_germline`) always strands a little mass on the short alleles;
    # what matters is how much. A few percent is a format artifact, a quarter of the distribution
    # is a real modelling error that quietly rescales every Pgen through that gene.
    per_allele = used.group_by(parent_allele_col).agg(
        tot=pl.col("p").sum(),
        lost=pl.col("p").filter(pl.col("_total") > pl.col("_max")).sum(),
    ).filter(pl.col("tot") > 0).with_columns(
        frac=pl.col("lost") / pl.col("tot")
    ).filter(pl.col("frac") > _LOST_INFO).sort("frac", descending=True)

    for r in per_allele.head(10).iter_rows(named=True):
        frac = r["frac"]
        sev = ("error" if frac > _LOST_ERROR else "warn" if frac > _LOST_WARN else "info")
        rows.append(_issue(sev, "deletion_unreachable",
                           f"{r[parent_allele_col]}: {frac:.1%} of its deletion mass lands on "
                           f"trims longer than its germline allows, so that probability is lost "
                           f"from every Pgen",
                           event=name, allele=r[parent_allele_col], value=float(frac)))
    if per_allele.height > 10:
        rows.append(_issue("info", "deletion_unreachable",
                           f"...and {per_allele.height - 10} further allele(s) with unreachable "
                           f"deletion mass (max {per_allele['frac'][10]:.1%})",
                           event=name, value=float(per_allele.height - 10)))

    for end, col in zip(ends, ("ndel",) if event.kind is EventKind.DELETION else ("ndel5", "ndel3")):
        mn = df.filter(pl.col("p") > 0)[col].min()
        if mn is not None and mn < -pal[end]:
            rows.append(_issue("error", "palindrome_max",
                               f"{col} reaches {mn} with p>0, beyond the declared "
                               f"palindrome_max[{end!r}] = {pal[end]}",
                               event=name, value=float(mn)))
    return rows


#: Mass at the largest insertion-length bin above which the support counts as clipped. Every real
#: model leaves a numerically-tiny tail there; only a meaningful amount indicates a hard cut-off.
_INS_TAIL_TOL = 1e-6

#: Thresholds on the fraction of an allele's deletion mass that is unreachable (see
#: :func:`_check_deletion`). Below ``_LOST_INFO`` it is not reported at all.
_LOST_INFO, _LOST_WARN, _LOST_ERROR = 1e-6, 0.01, 0.10


def _check_insertion(name: str, df: pl.DataFrame) -> list[dict]:
    tail = df.sort("length").tail(1)
    if tail.height and tail["p"][0] > _INS_TAIL_TOL:
        # The support was clipped before the distribution decayed -- from_arda's ins_max=40
        # placeholder does this, and EM cannot move mass past the edge it inherited.
        return [_issue("warn", "insertion_truncated",
                       f"the largest insertion length ({tail['length'][0]}) still carries "
                       f"p={tail['p'][0]:.3g}; the support is clipped, not exhausted",
                       event=name, value=float(tail["p"][0]))]
    return []


def _check_genomic_cross_refs(model: Model) -> list[dict]:
    """Every allele named in a marginal exists in the germline, and vice versa."""
    rows: list[dict] = []
    # d2_allele resolves against genes_d, whose own column is named 'd_allele'.
    known = {col: set(model.genomic[frame][f"{frame.split('_')[1]}_allele"].to_list())
             for col, frame in _ALLELE_GENOMIC.items() if frame in model.genomic}

    for name, df in model.tables.items():
        for col in df.columns:
            if col not in _ALLELE_GENOMIC or col not in known:
                continue
            unknown = sorted(set(df[col].to_list()) - known[col])
            for allele in unknown[:10]:
                rows.append(_issue("error", "allele_not_in_genomic",
                                   f"{allele!r} appears in {name!r} but not in "
                                   f"{_ALLELE_GENOMIC[col]}",
                                   event=name, segment=col.split("_")[0].upper(), allele=allele))
            if len(unknown) > 10:
                rows.append(_issue("error", "allele_not_in_genomic",
                                   f"...and {len(unknown) - 10} further unknown allele(s) in {name!r}",
                                   event=name, value=float(len(unknown) - 10)))

    for seg, event in (("v", "v_choice"), ("j", "j_choice"), ("d", "d_gene"), ("d2", "d2_gene")):
        frame = f"genes_{'d' if seg == 'd2' else seg}"
        if frame not in model.genomic or event not in model.tables:
            continue
        allele_col = f"{seg}_allele"
        gcol = "d_allele" if seg == "d2" else allele_col
        mass = _choice_mass(model, event, allele_col)
        for r in model.genomic[frame].iter_rows(named=True):
            allele, p = r[gcol], mass.get(r[gcol], None)
            has_germline = bool(r["cut_segment"])
            if p is None:
                rows.append(_issue("warn", "genomic_not_in_tables",
                                   f"{allele!r} is in {frame} but has no row in {event!r}; "
                                   f"it can never be generated",
                                   event=event, segment=seg.upper(), allele=allele))
            elif not has_germline and p > 0:
                # OLGA ships a handful of ORF alleles with an EMPTY CDR3-region germline while
                # still giving them usage; vdjtools reproduces that for exact-Pgen fidelity. The
                # consequence is real -- Pgen is exactly 0 for every clonotype using the gene --
                # so it is reported, but it is a property of the source model, not a defect here.
                # `from_olga(derive_orf=True)` reconstructs these.
                rows.append(_issue("warn", "unscoreable_gene_mass",
                                   f"{allele!r} carries p={p:.3g} in {event!r} but has no "
                                   f"CDR3-region germline, so its Pgen is always 0",
                                   event=event, segment=seg.upper(), allele=allele, value=float(p)))
            elif r["functional"] and p is not None and p <= 0:
                # The absorbing-state failure: a real gene EM (or a bad import) pinned to zero
                # makes Pgen exactly 0 for every clonotype using it, silently.
                rows.append(_issue("warn", "functional_zero_mass",
                                   f"{allele!r} is functional but carries zero mass in {event!r}",
                                   event=event, segment=seg.upper(), allele=allele, value=0.0))
    return rows


def _check_germline(model: Model, germline) -> list[dict]:
    """Reconcile the model's own germline against an external reference library."""
    if isinstance(germline, str) and germline == "none":
        return []
    if isinstance(germline, pl.DataFrame):
        ref_df = germline
    else:
        from .reference import load_germline
        try:
            ref_df = load_germline(model.locus, model.organism)
        except (ImportError, ValueError):
            # A custom library or a non-arda organism has no reference to reconcile against.
            return []

    ref = {(r["segment"], r["allele"]): r["sequence"]
           for r in ref_df.iter_rows(named=True)}
    rows: list[dict] = []
    for seg in ("v", "j", "d"):
        frame = model.genomic.get(f"genes_{seg}")
        if frame is None:
            continue
        n_absent = n_diff = 0
        for r in frame.iter_rows(named=True):
            allele = r[f"{seg}_allele"]
            want = ref.get((seg.upper(), allele))
            if want is None:
                n_absent += 1
            elif want != r["cdr3_segment"]:
                n_diff += 1
        if n_absent:
            rows.append(_issue("warn", "germline_source",
                               f"{n_absent} {seg.upper()} allele(s) are absent from the reference "
                               f"germline library",
                               segment=seg.upper(), value=float(n_absent)))
        if n_diff:
            rows.append(_issue("warn", "germline_source",
                               f"{n_diff} {seg.upper()} allele(s) have a CDR3-region germline "
                               f"differing from the reference library",
                               segment=seg.upper(), value=float(n_diff)))
    return rows
