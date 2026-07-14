"""TCRnet neighbourhood-enrichment (convergence) test — delegated to vdjmatch/seqtree.

The classic TCRnet / ``CalcDegreeStats`` analysis asks, for every clonotype in a
sample, whether it has *more* close CDR3 neighbours than a background/generative
process would produce — the signature of antigen-driven convergent selection. The
legacy tool counted a clonotype's within-sample degree and compared it to a control
sample's degree under a grouping (V/VJ/VJL) that bounded the comparison scope.

This reimplements it on a **finite-sample, control-calibrated footing** by delegating
to :func:`vdjmatch.evalue.query_evalues` (the ``seqtree`` E-value engine):

- **target** = a fuzzy ``seqtree.Index`` over the *sample's own* unique CDR3s — so a
  clonotype's neighbour count is its within-sample degree;
- **control** = a matched background repertoire (:func:`vdjmatch.evalue.background`);
- for each query CDR3, ``E = (n_target / n_control_index_size)·n_control`` is the
  neighbour count expected from the background, and ``p_enrichment`` is the Poisson
  tail probability of seeing at least the observed within-sample degree by chance.

**Locus handling.** A TCRnet neighbourhood is only meaningful within one locus (a TRA
sequence has no true neighbours in a TRB background). When neither ``control`` nor
``locus`` is given, the sample is split by the locus of its ``v_call`` and each locus
is scored against its own matched background — mirroring the legacy V-grouping that
kept loci disjoint. Clonotypes whose locus cannot be resolved (null ``v_call``) are
dropped with a warning. Passing ``control`` or ``locus`` explicitly overrides this and
scores every clonotype against that single background.

Scope note: the legacy default was ``s,id,t = 1,0,1`` (one substitution, no indels,
one total edit), i.e. vdjmatch ``"1,0,0,1"`` — the default here. ``exclude_exact`` is
``True`` by default so a clonotype's own identical copy in the target/control is not
counted as a neighbour.
"""
from __future__ import annotations

import warnings

import polars as pl

from ..io.schema import JUNCTION_AA, COUNT, J_CALL, LOCUS, V_CALL, add_locus

_VDJMATCH_HINT = (
    "vdjmatch is required for vdjtools.overlap.tcrnet; install the extra with "
    "vdjmatch is a base dependency of vdjtools -- reinstall with `pip install --force-reinstall vdjtools`."
)

_COLS = [JUNCTION_AA, V_CALL, J_CALL, COUNT, "n_neighbors", "n_control", "E", "p_enrichment", "p_any", LOCUS]


def _require_vdjmatch():
    """Import the vdjmatch pieces used by :func:`tcrnet`; raise a helpful error if missing."""
    try:
        import vdjmatch.evalue as evalue  # noqa: F401
        from vdjmatch.match.scope import search_params  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without vdjmatch
        raise ImportError(_VDJMATCH_HINT) from exc
    return evalue, search_params


def _collapse(sample: pl.DataFrame) -> pl.DataFrame:
    """Collapse to unique CDR3s, carrying the most-abundant V/J and the summed count."""
    return (sample.sort(COUNT, descending=True)
                  .group_by(JUNCTION_AA, maintain_order=True)
                  .agg(pl.col(V_CALL).first(), pl.col(J_CALL).first(),
                       pl.col(COUNT).sum().alias(COUNT)))


def _score(agg, control, params, exclude_exact, threads, evalue, seqtree, locus):
    """Score one already-collapsed (single-locus) frame against ``control``."""
    cdr3s = agg[JUNCTION_AA].to_list()
    target = seqtree.Index.build(cdr3s, alphabet="aa")
    ev = evalue.query_evalues(target, control, cdr3s, params,
                              threads=threads, exclude_exact=exclude_exact)
    ev = ev.rename({"query_cdr3": JUNCTION_AA, "n_target": "n_neighbors"})
    return (agg.join(ev, on=JUNCTION_AA, how="left")
               .with_columns(pl.lit(locus, dtype=pl.Utf8).alias(LOCUS)))


def tcrnet(sample: pl.DataFrame, control=None, scope: str = "1,0,0,1",
           locus: str | None = None, species: str = "human",
           exclude_exact: bool = True, threads: int = 0) -> pl.DataFrame:
    """Per-clonotype neighbourhood-enrichment (TCRnet) test for one sample.

    Collapses the sample to unique CDR3s, builds a fuzzy ``seqtree.Index`` over them
    (the within-sample neighbourhood target), and queries each CDR3 against that
    target and a background control via :func:`vdjmatch.evalue.query_evalues`. When
    neither ``control`` nor ``locus`` is given the sample is scored **per locus**, each
    against its own matched background (see the module docstring).

    Args:
        sample: Clonotype frame (canonical schema).
        control: A prebuilt ``seqtree.Index`` background, or ``None`` to load matched
            background(s) via :func:`vdjmatch.evalue.background`. Passing one overrides
            the per-locus split (every clonotype is scored against it).
        scope: vdjmatch edit-distance scope ``"subs,ins,dels,total"`` defining the
            neighbourhood ball (default one substitution).
        locus: Force a single background locus (e.g. ``"TRB"``); overrides the
            per-locus split. When ``None`` (and ``control`` is ``None``) the sample is
            split by the locus of its ``v_call``.
        species: Species for the background control (default ``"human"``).
        exclude_exact: Drop distance-0 (identical / self) hits from both target and
            control counts, so a clonotype is not its own neighbour (default ``True``).
        threads: Worker threads for the native search (``0`` = all cores).

    Returns:
        One row per unique clonotype (per locus) with columns ``junction_aa, v_call,
        j_call, duplicate_count, n_neighbors`` (within-sample degree), ``n_control``
        (background degree), ``E`` (expected neighbours), ``p_enrichment`` (Poisson
        tail), ``p_any``, and ``locus`` — sorted by ascending ``p_enrichment``.

    Raises:
        ImportError: If vdjmatch is not importable (it is a base dependency).
        ValueError: If ``control`` and ``locus`` are ``None`` and no clonotype has a
            resolvable locus.
    """
    evalue, search_params = _require_vdjmatch()
    import seqtree
    params = search_params(scope)

    # Explicit override: a single caller-provided (or forced-locus) background scores
    # every clonotype, regardless of its own locus.
    if control is not None or locus is not None:
        if control is None:
            control = evalue.background(locus=locus, species=species)
        out = _score(_collapse(sample), control, params, exclude_exact, threads,
                     evalue, seqtree, locus)
        return out.select(_COLS).sort("p_enrichment")

    # Auto: split by v_call locus and score each locus against its matched background.
    withloc = add_locus(sample)
    n_null = withloc.filter(pl.col(LOCUS).is_null()).height
    withloc = withloc.filter(pl.col(LOCUS).is_not_null())
    if withloc.is_empty():
        raise ValueError("cannot infer locus from sample v_call; pass locus= or control=")
    if n_null:
        warnings.warn(f"tcrnet: dropped {n_null} clonotype(s) with unresolvable locus "
                      "(null v_call); pass locus= or control= to score them.")
    parts = [
        _score(_collapse(withloc.filter(pl.col(LOCUS) == loc)),
               evalue.background(locus=loc, species=species),
               params, exclude_exact, threads, evalue, seqtree, loc)
        for loc in sorted(withloc[LOCUS].unique().to_list())
    ]
    return pl.concat(parts).select(_COLS).sort("p_enrichment")
