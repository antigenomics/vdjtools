"""Reconstruct full-length nt contigs from (V, J, CDR3) and run them through arda.

OLGA generation (and OLGA-synthetic bootstrap data) emits only ``(V call, J call, CDR3 nt)``,
but arda — the aligner that drives scenario enumeration for *real* reads — consumes full nt
reads. ``stitch_contig`` rebuilds the read by prepending the V germline 5' of the CDR3 anchor
and appending the J germline 3' of it, so synthetic and real reads flow through the identical
``arda.annotate_sequences → scenarios → EM`` path (with known ground truth on the synthetic side).
"""
from __future__ import annotations

import polars as pl

from .model import Model


def _germline_lookup(model: Model, seg: str) -> dict[str, tuple[str, int]]:
    g = model.genomic[f"genes_{seg}"]
    return {r[0]: (r[1], r[2]) for r in g.select([f"{seg}_allele", "full_germline", "anchor"]).iter_rows()}


def stitch_contig(model: Model, v: str, j: str, cdr3_nt: str) -> str | None:
    """Rebuild a full nt contig: V framework 5' of Cys104 + CDR3 + J framework 3' of [FW]118.

    Args:
        model: Model providing full germline + anchors (``genes_v``/``genes_j``).
        v, j: V and J allele names.
        cdr3_nt: The junction-space CDR3 nt (Cys104 → [FW]118 inclusive, i.e. AIRR ``junction``).

    Returns:
        The contig, or ``None`` if either gene lacks a usable anchor / full germline.
    """
    vg = _germline_lookup(model, "v")
    jg = _germline_lookup(model, "j")
    if v not in vg or j not in jg:
        return None
    fv, av = vg[v]
    fj, aj = jg[j]
    if av < 0 or aj < 0 or not fv or not fj:
        return None
    return fv[:av] + cdr3_nt + fj[aj + 3:]


def stitch_frame(model: Model, gen: pl.DataFrame, *, cdr3_col: str = "junction_nt") -> pl.DataFrame:
    """Add a ``contig`` column to a generated frame (rows with no anchor drop to null)."""
    contigs = [stitch_contig(model, r["v_call"], r["j_call"], r[cdr3_col]) for r in gen.to_dicts()]
    return gen.with_columns(contig=pl.Series("contig", contigs, dtype=pl.Utf8))


def annotate(contigs: list[str], *, organism: str = "human") -> pl.DataFrame:
    """Annotate nt contigs with arda → a polars frame of the calls the scenario/EM path needs.

    Returns columns ``v_call, d_call, j_call, junction, junction_aa, productive`` (arda's best
    alignment). arda is the ``[model]`` extra.
    """
    import arda

    recs = arda.annotate_sequences(contigs, seqtype="nt", organism=organism)
    return pl.DataFrame(
        {
            "v_call": [r.get("v_call") for r in recs],
            "d_call": [r.get("d_call") for r in recs],
            "j_call": [r.get("j_call") for r in recs],
            "junction": [r.get("junction") for r in recs],
            "junction_aa": [r.get("junction_aa") for r in recs],
            "productive": [r.get("productive") for r in recs],
        }
    )
