"""Tandem-D (D-D) model extension: upgrade a single-D VDJ model to ``n_D ∈ {1, 2}``.

OLGA (and every bootstrap model) assumes exactly one D per rearrangement. Real IGH and TRD
repertoires contain *tandem* D-D joins — two D segments in the junction (``V-D1-D2-J``) — which
no OLGA/IGoR model can represent. :func:`to_dd` adds the four events that a tandem rearrangement
needs on top of an existing single-D :class:`~vdjtools.model.model.Model`:

- ``n_d``     — ``P(n_D)`` over ``{1, 2}`` (the D-count prior; ``2`` = a tandem join).
- ``d2_gene`` — ``P(D2 | D1)`` (the genomic-order mask lives in this table's zeros).
- ``d2_del``  — ``P(delD2_5', delD2_3' | D2)`` joint trimming of the second D.
- ``dd_ins`` / ``dd_dinucl`` — the N-region between the two Ds.

The upgrade is an *initialisation*: it seeds ``P(n_D=2) = p_nd2`` and copies the single-D
deletion/insertion profiles for the new events. EM on real tandem-D reads then learns the true
values (the enumeration in :func:`~vdjtools.model.pgen.pgen_nt` already sums the n_D=2 scenarios).
With ``p_nd2 = 0`` the model is generatively identical to its single-D input.
"""
from __future__ import annotations

import polars as pl

from .events import Event, EventKind
from .model import Model
from .schema import Manifest


def has_tandem(model: Model) -> bool:
    """True if the model places positive probability on a tandem (``n_D=2``) rearrangement.

    Guards the paths that do not yet sum/sample the ``n_D=2`` scenarios (amino-acid Pgen, the
    native ``_core`` Pgen/EM, generation, EM inference) so that a tandem model raises rather than
    silently returning a single-D result. The predicate is *positive ``P(n_D=2)``* alone — so a
    ``to_dd(..., p_nd2=0)`` model (generatively single-D, ``d2_gene`` present but unused) correctly
    returns ``False`` and flows through the single-D fast paths byte-identically.
    """
    t = model.tables.get("n_d")
    if t is None:
        return False
    row = t.filter(pl.col("n_d") == 2)
    return row.height > 0 and float(row["p"].sum()) > 0.0


def to_dd(model: Model, *, p_nd2: float = 0.05) -> Model:
    """Return a copy of ``model`` upgraded to a tandem-D (``n_D ∈ {1, 2}``) model.

    Args:
        model: A single-D VDJ :class:`~vdjtools.model.model.Model` (raises for VJ or already-DD).
        p_nd2: Initial ``P(n_D = 2)`` prior mass (``0`` reproduces the single-D model exactly).

    Returns:
        A validated D-D :class:`~vdjtools.model.model.Model` with the four tandem events added.
    """
    if model.chain_type != "VDJ":
        raise ValueError(f"D-D extension needs a VDJ model, got {model.chain_type}")
    if "d2_gene" in model.tables:
        raise ValueError("model already has a tandem-D (d2_gene) event")
    if not 0.0 <= p_nd2 < 1.0:
        raise ValueError("p_nd2 must be in [0, 1)")

    t = dict(model.tables)
    d_alleles = model.genomic["genes_d"]["d_allele"].to_list()

    # P(n_D): {1, 2}. (n_D=0 is absorbed into n_D=1 via a fully-trimmed D — see pgen._vdj_middle.)
    t["n_d"] = pl.DataFrame({"n_d": [1, 2], "p": [1.0 - p_nd2, p_nd2]}).with_columns(
        pl.col("n_d").cast(pl.UInt8)
    )

    # P(D2 | D1): initialise independent of D1 from the marginal P(D) (Σ_J P(J) P(D|J)).
    from .analyze import gene_marginal

    pd = gene_marginal(model, "d")
    d2_rows = [
        {"d_allele": d1, "d2_allele": d2, "p": pd.get(d2, 0.0)}
        for d1 in d_alleles
        for d2 in d_alleles
    ]
    t["d2_gene"] = pl.DataFrame(d2_rows)

    # D2 deletion + DD insertion: copy the single-D / DJ profiles as the initialisation.
    t["d2_del"] = model.tables["d_del"].rename({"d_allele": "d2_allele"})
    t["dd_ins"] = model.tables["dj_ins"]
    t["dd_dinucl"] = model.tables["dj_dinucl"]

    events = dict(model.manifest.events)
    events["n_d"] = Event("n_d", EventKind.N_D)
    events["d2_gene"] = Event("d2_gene", EventKind.GENE_CHOICE, ("d_gene",))
    events["d2_del"] = Event("d2_del", EventKind.DELETION_2D, ("d2_gene",))
    events["dd_ins"] = Event("dd_ins", EventKind.INS_LENGTH)
    events["dd_dinucl"] = Event("dd_dinucl", EventKind.DINUCLEOTIDE)

    manifest = Manifest(
        locus=model.manifest.locus,
        organism=model.manifest.organism,
        chain_type="VDJ",
        events=events,
        # D2 draws from the same D germline as D1, so it reuses the D palindrome maxima (d_5/d_3);
        # the enumeration reads those for both Ds — no separate d2_5/d2_3 keys.
        palindrome_max=dict(model.manifest.palindrome_max),
        source=f"{model.manifest.source}+dd",
    )
    return Model(manifest=manifest, tables=t, genomic=dict(model.genomic)).validate()
