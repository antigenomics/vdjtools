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

Scope note: the legacy default was ``s,id,t = 1,0,1`` (one substitution, no indels,
one total edit), i.e. vdjmatch ``"1,0,0,1"`` — the default here. ``exclude_exact`` is
``True`` by default so a clonotype's own identical copy in the target/control is not
counted as a neighbour.
"""
from __future__ import annotations

import polars as pl

from ..io.schema import CDR3_AA, COUNT, J_CALL, LOCUS, V_CALL, add_locus

_VDJMATCH_HINT = (
    "vdjmatch is required for vdjtools.overlap.tcrnet; install the extra with "
    "`pip install 'vdjtools[overlap]'` (or `pip install vdjmatch>=0.0.1`)."
)


def _require_vdjmatch():
    """Import the vdjmatch pieces used by :func:`tcrnet`; raise a helpful error if missing."""
    try:
        import vdjmatch.evalue as evalue  # noqa: F401
        from vdjmatch.match.scope import search_params  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without vdjmatch
        raise ImportError(_VDJMATCH_HINT) from exc
    return evalue, search_params


def _infer_locus(sample: pl.DataFrame) -> str | None:
    """Most common three-letter locus over the sample's ``v_call`` (e.g. ``"TRB"``)."""
    loci = add_locus(sample)[LOCUS].drop_nulls()
    if loci.is_empty():
        return None
    return loci.mode().sort()[0]


def tcrnet(sample: pl.DataFrame, control=None, scope: str = "1,0,0,1",
           locus: str | None = None, species: str = "human",
           exclude_exact: bool = True, threads: int = 0) -> pl.DataFrame:
    """Per-clonotype neighbourhood-enrichment (TCRnet) test for one sample.

    Collapses the sample to unique CDR3s, builds a fuzzy ``seqtree.Index`` over them
    (the within-sample neighbourhood target), and queries each CDR3 against that
    target and a background control via :func:`vdjmatch.evalue.query_evalues`.

    Args:
        sample: Clonotype frame (canonical schema).
        control: A prebuilt ``seqtree.Index`` background, or ``None`` to load a
            matched background via :func:`vdjmatch.evalue.background`.
        scope: vdjmatch edit-distance scope ``"subs,ins,dels,total"`` defining the
            neighbourhood ball (default one substitution).
        locus: Locus for the background (e.g. ``"TRB"``); inferred from the sample's
            ``v_call`` when ``None``.
        species: Species for the background control (default ``"human"``).
        exclude_exact: Drop distance-0 (identical / self) hits from both target and
            control counts, so a clonotype is not its own neighbour (default ``True``).
        threads: Worker threads for the native search (``0`` = all cores).

    Returns:
        One row per unique clonotype with columns ``cdr3_aa, v_call, j_call,
        duplicate_count, n_neighbors`` (within-sample degree), ``n_control``
        (background degree), ``E`` (expected neighbours), ``p_enrichment`` (Poisson
        tail), and ``p_any`` — sorted by ascending ``p_enrichment``.

    Raises:
        ImportError: If vdjmatch is not installed (see the ``overlap`` extra).
        ValueError: If ``control`` is ``None`` and the locus cannot be determined.
    """
    evalue, search_params = _require_vdjmatch()
    import seqtree

    # Collapse to unique CDR3s, carrying the most-abundant clonotype's V/J and the
    # summed count (the degree is computed over the unique CDR3 footprint).
    agg = (sample.sort(COUNT, descending=True)
                 .group_by(CDR3_AA, maintain_order=True)
                 .agg(pl.col(V_CALL).first(), pl.col(J_CALL).first(),
                      pl.col(COUNT).sum().alias(COUNT)))
    cdr3s = agg[CDR3_AA].to_list()

    if control is None:
        locus = locus or _infer_locus(sample)
        if locus is None:
            raise ValueError("cannot infer locus from sample v_call; pass locus= or control=")
        control = evalue.background(locus=locus, species=species)

    target = seqtree.Index.build(cdr3s, alphabet="aa")
    params = search_params(scope)
    ev = evalue.query_evalues(target, control, cdr3s, params,
                              threads=threads, exclude_exact=exclude_exact)

    ev = ev.rename({"query_cdr3": CDR3_AA, "n_target": "n_neighbors"})
    return (agg.join(ev, on=CDR3_AA, how="left")
               .select(CDR3_AA, V_CALL, J_CALL, COUNT,
                       "n_neighbors", "n_control", "E", "p_enrichment", "p_any")
               .sort("p_enrichment"))
