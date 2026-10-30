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

from dataclasses import replace

import polars as pl

from .events import Event, EventKind
from .model import Model


def d_count_probabilities(model: Model) -> dict[int, float]:
    """Marginal D-history count probabilities, including the J-dependent history prior."""
    if "d_history" in model.tables:
        rows = (model.tables["d_history"].join(model.tables["j_choice"].rename({"p": "_pj"}), on="j_allele")
                .group_by("n_d").agg((pl.col("p") * pl.col("_pj")).sum().alias("p")))
        return dict(rows.iter_rows())
    table = model.tables.get("n_d")
    return dict(table.select("n_d", "p").iter_rows()) if table is not None else {}


def has_tandem(model: Model) -> bool:
    """True if the model places positive probability on a tandem (``n_D=2``) rearrangement.

    Guards the paths that do not yet sum/sample the ``n_D=2`` scenarios (amino-acid Pgen, the
    native ``_core`` Pgen/EM, generation, EM inference) so that a tandem model raises rather than
    silently returning a single-D result. The predicate is *positive ``P(n_D=2)``* alone — so a
    ``to_dd(..., p_nd2=0)`` model (generatively single-D, ``d2_gene`` present but unused) correctly
    returns ``False`` and flows through the single-D fast paths byte-identically.
    """
    return d_count_probabilities(model).get(2, 0.0) > 0.0


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
    from .reference import forbidden_tandem_pairs
    forbidden = forbidden_tandem_pairs(d_alleles,model.locus,model.organism)
    d2_rows = [
        {"d_allele": d1, "d2_allele": d2, "p": 0. if (d1,d2) in forbidden else pd.get(d2, 0.0)}
        for d1 in d_alleles
        for d2 in d_alleles
    ]
    t["d2_gene"] = pl.DataFrame(d2_rows)
    total = pl.col('p').sum().over('d_allele')
    t['d2_gene'] = t['d2_gene'].with_columns(p=pl.when(total>0).then(pl.col('p')/total).otherwise(0.))

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
    if "dj_first_nt" in model.tables:
        t["dd_first_nt"] = model.tables["dj_first_nt"]
        events["dd_first_nt"] = Event("dd_first_nt", EventKind.NUCLEOTIDE)

    manifest = replace(model.manifest,
        events=events,
        # D2 draws from the same D germline as D1, so it reuses the D palindrome maxima (d_5/d_3);
        # the enumeration reads those for both Ds — no separate d2_5/d2_3 keys.
        palindrome_max=dict(model.manifest.palindrome_max),
        source=f"{model.manifest.source}+dd",
    )
    return Model(manifest=manifest, tables=t, genomic=dict(model.genomic), training=model.training).validate()


def to_history(model: Model, *, p_nd2: float = 0.05) -> Model:
    """Build J-conditioned joint single/tandem histories on genuine physical support.

    ``d_history`` stores ``P(n_D,D1[,D2]|J)``. Single-D rows use ``d2_allele=""``.
    A J with no valid tandem receives only single-D mass. Tandem D1/D2 deletion
    tables normalize on states retaining at least one nt; single D may vanish.
    This explicit new schema leaves the legacy :func:`to_dd` representation intact.
    """
    from .reference import forbidden_dj_pairs, forbidden_tandem_pairs

    if model.chain_type != "VDJ" or "d_history" in model.tables:
        raise ValueError("to_history requires a legacy single/tandem VDJ model")
    if not 0 <= p_nd2 < 1:
        raise ValueError("p_nd2 must be in [0, 1)")
    legacy = model if "d2_gene" in model.tables else to_dd(model, p_nd2=p_nd2)
    d = model.genomic["genes_d"].filter(pl.col("functional"))["d_allele"].to_list()
    js = model.genomic["genes_j"]["j_allele"].to_list()
    bad_dj = forbidden_dj_pairs(d, js, model.locus)
    bad_dd = forbidden_tandem_pairs(d, model.locus, model.organism)

    def norm(frame, keys):
        total = pl.col("p").sum().over(keys)
        return frame.with_columns(p=pl.when(total > 0).then(pl.col("p") / total).otherwise(0.0))

    single = model.tables["d_gene"].filter(pl.col("d_allele").is_in(d))
    if bad_dj:
        single = single.join(pl.DataFrame(sorted(bad_dj), schema=["d_allele", "j_allele"], orient="row"),
                             on=["d_allele", "j_allele"], how="anti")
    single = norm(single, ["j_allele"])
    positive_j = model.tables["j_choice"].filter(pl.col("p") > 0)["j_allele"]
    supported = single.group_by("j_allele").agg(pl.col("p").sum()).filter(pl.col("p") > 0)["j_allele"]
    if set(positive_j) - set(supported):
        raise ValueError("positive J choice has no supported single-D history")
    tandem = (single.rename({"p": "_p1"}).join(legacy.tables["d2_gene"], on="d_allele")
              .filter(pl.col("d2_allele").is_in(d)).with_columns(p=pl.col("p") * pl.col("_p1"))
              .drop("_p1"))
    if bad_dd:
        tandem = tandem.join(pl.DataFrame(sorted(bad_dd), schema=["d_allele", "d2_allele"], orient="row"),
                             on=["d_allele", "d2_allele"], how="anti")
    if bad_dj:
        tandem = tandem.join(pl.DataFrame(sorted(bad_dj), schema=["d2_allele", "j_allele"], orient="row"),
                             on=["d2_allele", "j_allele"], how="anti")
    tandem = norm(tandem, ["j_allele"])
    tandem_js = tandem.filter(pl.col("p") > 0)["j_allele"].unique()
    single = single.with_columns(
        p=pl.col("p") * pl.when(pl.col("j_allele").is_in(tandem_js.implode())).then(1-p_nd2).otherwise(1.0),
        n_d=pl.lit(1, dtype=pl.UInt8), d2_allele=pl.lit(""))
    tandem = tandem.with_columns(p=pl.col("p") * p_nd2, n_d=pl.lit(2, dtype=pl.UInt8))
    history = pl.concat([single.select("j_allele", "n_d", "d_allele", "d2_allele", "p"),
                         tandem.select("j_allele", "n_d", "d_allele", "d2_allele", "p")])
    tables = {name: table for name, table in legacy.tables.items() if name not in ("n_d", "d_gene", "d2_gene")}
    tables["d_history"] = history
    pal = model.manifest.palindrome_max
    lengths = model.genomic["genes_d"].select("d_allele", pl.col("cut_segment").str.len_chars().alias("_length"))
    for source, target, allele in (("d_del", "d_tandem_del", "d_allele"),
                                   ("d2_del", "d2_del", "d2_allele")):
        nonempty = (legacy.tables[source].join(lengths.rename({"d_allele": allele}), on=allele)
                    .filter(pl.col("ndel5") + pl.col("ndel3") + pal["d_5"] + pal["d_3"] < pl.col("_length"))
                    .drop("_length"))
        nonempty = norm(nonempty, [allele])
        used = set(history.filter((pl.col("n_d") == 2) & (pl.col("p") > 0))[allele])
        if used - set(nonempty.filter(pl.col("p") > 0)[allele]):
            raise ValueError("positive tandem history has no nonempty D deletion support")
        tables[target] = nonempty
    events = {name: event for name, event in legacy.manifest.events.items()
              if name not in ("n_d", "d_gene", "d2_gene")}
    events["d_history"] = Event("d_history", EventKind.D_HISTORY, ("j_choice",))
    events["d_del"] = Event("d_del", EventKind.DELETION_2D, ("d_history",))
    events["d_tandem_del"] = Event("d_tandem_del", EventKind.DELETION_2D, ("d_history",))
    events["d2_del"] = Event("d2_del", EventKind.DELETION_2D, ("d_history",))
    return Model(manifest=replace(legacy.manifest, events=events, model_version="3.0.0",
                                  source=f"{model.manifest.source}+history"),
                 tables=tables, genomic=model.genomic, training=model.training).validate()


def validate_history(model: Model, *, tol: float = 1e-8) -> None:
    """Reject incomplete or impossible joint-history probability before any mass is dropped."""
    from .reference import forbidden_dj_pairs, forbidden_tandem_pairs

    if model.chain_type != "VDJ":
        raise ValueError("joint histories require a VDJ model")
    h = model.tables.get("d_history")
    required = {"j_allele", "n_d", "d_allele", "d2_allele", "p"}
    if h is None or not required <= set(h.columns):
        raise ValueError("joint history table has missing columns")
    genes = model.genomic.get("genes_d")
    if genes is None or "genes_j" not in model.genomic or "j_choice" not in model.tables:
        raise ValueError("joint histories require D/J germline and J choice")
    ds = genes.filter(pl.col("functional"))["d_allele"].to_list()
    js = model.genomic["genes_j"]["j_allele"].to_list()
    invalid = (~pl.col("n_d").is_in([1, 2]) | ~pl.col("d_allele").is_in(ds) |
               ~pl.col("j_allele").is_in(js) |
               ((pl.col("n_d") == 1) & (pl.col("d2_allele") != "")) |
               ((pl.col("n_d") == 2) & ~pl.col("d2_allele").is_in(ds)) |
               pl.col("p").is_null() | pl.col("p").is_nan() | (pl.col("p") < 0) | (pl.col("p") > 1 + tol))
    if h.filter(invalid.fill_null(True)).height or h.select(sorted(required - {"p"})).n_unique() != h.height:
        raise ValueError("invalid or duplicate joint-history realization")
    sums = h.group_by("j_allele").agg(pl.col("p").sum().alias("_sum"))
    needed = model.tables["j_choice"].filter(pl.col("p") > 0).select("j_allele")
    if needed.join(sums, on="j_allele", how="left").filter(
            pl.col("_sum").is_null() | ((pl.col("_sum") - 1).abs() > tol)).height:
        raise ValueError("joint histories must sum to one for every positive J")
    dj = forbidden_dj_pairs(ds, js, model.locus)
    dd = forbidden_tandem_pairs(ds, model.locus, model.organism)
    tandem = h.filter(pl.col("n_d") == 2)
    for frame, pairs, columns in ((h, dj, ["d_allele", "j_allele"]),
                                  (tandem, dj, ["d2_allele", "j_allele"]),
                                  (tandem, dd, ["d_allele", "d2_allele"])):
        if pairs and frame.join(pl.DataFrame(sorted(pairs), schema=columns, orient="row"), on=columns).height:
            raise ValueError("joint history contains a physically forbidden path")
    lengths = genes.select("d_allele", pl.col("cut_segment").str.len_chars().alias("_length"))
    positive = h.filter(pl.col("p") > 0)
    for name, column, rows, minimum in (("d_del", "d_allele", positive.filter(pl.col("n_d") == 1), 0),
                                        ("d_tandem_del", "d_allele", positive.filter(pl.col("n_d") == 2), 1),
                                        ("d2_del", "d2_allele", positive.filter(pl.col("n_d") == 2), 1)):
        table = model.tables.get(name)
        if table is None or not {column, "ndel5", "ndel3", "p"} <= set(table.columns):
            raise ValueError(f"joint histories require {name} deletion support")
        support = table.group_by(column).agg(pl.col("p").sum().alias("_sum"))
        if rows.select(column).unique().join(support, on=column, how="left").filter(
                pl.col("_sum").is_null() | ((pl.col("_sum") - 1).abs() > tol)).height:
            raise ValueError(f"positive joint history lacks normalized {name} support")
        checked = table.join(lengths.rename({"d_allele": column}), on=column, how="left")
        i5 = pl.col("ndel5") + model.manifest.palindrome_max["d_5"]
        i3 = pl.col("ndel3") + model.manifest.palindrome_max["d_3"]
        if checked.filter((pl.col("p") > 0) & (pl.col("_length").is_null() | (i5 < 0) | (i3 < 0) |
                          (pl.col("_length") - i5 - i3 < minimum))).height:
            raise ValueError(f"joint history {name} contains impossible positive deletion mass")
