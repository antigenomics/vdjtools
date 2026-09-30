"""The junction pipeline: cdr3 fix -> V/J markup -> nucleotide guess -> D markup.

One call for a bare ``(junction_aa, V, J, species)`` record, which is all a VDJdb-style database
has. Specified by M.S. on 2026-09-30; the design note is ``docs/junction_pipeline.md``.

Four stages, each of them arda's or vdjtools' existing entry point, in the only order that works --
each answer is the next one's input:

1. **arda** ``cdr3fix.markup_batch`` confirms or replaces the V/J call, repairs the junction, and
   places the boundaries (``v_end`` / ``j_start`` in residues, ``v_end_nt`` / ``j_start_nt`` in
   nucleotides).
2. **vdjtools** :func:`~vdjtools.model.viterbi.infer_nt_batch` infers the most plausible
   NUCLEOTIDE junction for that amino-acid junction under the recombination model, conditioned on
   the confirmed calls. Natively batched and threaded, one call per ``(organism, locus)``. The
   templated flanks are then written from the called allele's own germline: only the N region is
   unknowable from amino acids, and the reconstruction differed from the observed V span on 23.93 %
   of real rearrangements before this.
3. **vdjtools** :func:`~vdjtools.model.native.best_aa_scenarios_batch` names the D GENE: the
   model's own scenario weights, grouped by D and normalised. Same native call as stage 2, same
   model, already loaded -- no second model, no fitted prior table, no per-locus tempering constant.
4. **arda** ``_markup.d_local_align`` PLACES that gene: greedy gapless local alignment of its
   germline inside the ``[v_end_nt, j_start_nt)`` interior, ungated, so a row gets coordinates
   whether or not the alignment is confident. The nucleotide span is then folded onto the residues
   whose codons it touches.

**Why the nucleotide detour is the point.** A D that contributes one amino acid is invisible in a
translated junction, but it is not invisible in nucleotides: the residues at each end of the D are
part germline-D and part N region, so their codons carry D bases. Searching the translated junction
throws that away; searching the inferred nucleotide junction keeps it, and the D is then allowed to
OVERLAP ``v_end`` / ``j_start`` for exactly the same reason -- an exonuclease does not cut on a
codon boundary.

**Naming the D and placing it are separate questions, and only one estimator answers each.**
Letting the alignment choose the gene as well as the position was measured and dropped. Against the
nucleotide D calls of 4,000 real human TRB rearrangements (`isalgo/airr_control`, never shown to the
pipeline):

====================================================  ==========  ================
route                                                 gene right  has coordinates
====================================================  ==========  ================
E-value-gated alignment chooses and places              47.93 %       55.75 %
gated alignment, model posterior where it declines      71.40 %       55.75 %
**model names the gene, greedy alignment places it**  **74.30 %**   **99.70 %**
====================================================  ==========  ================

⛔ The gate is not what held the alignment back -- opening it reaches 67.10 % -- and the alignment's
gene call is not better where it IS confident: on the 2,230 rows it speaks for, the model's own
posterior is right 86.32 % against its 85.96 %, and a hybrid of the two scores *below* the model
alone (74.12 % against 74.33 %). So there is one gene estimator, and the aligner does what only it
can do: say where.

⚠ Position is the one axis the gate bought something on -- ``d_start_nt`` is exact on 60.88 % of
correctly-called rows here against 66.67 % under the gate -- but it is exact on 1,807 rows rather
than 1,278, because it answers 3,988 rows rather than 2,230. For a database drawing V/N/D/N/J that
is the trade to take.
"""
from __future__ import annotations

from typing import Iterable, Sequence

import polars as pl

__all__ = ["annotate_junctions", "JUNCTION_COLUMNS"]

#: Columns :func:`annotate_junctions` returns, in order. Stage 1 is arda's markup, stage 2 the
#: nucleotide guess, stage 3/4 the D.
JUNCTION_COLUMNS = (
    # stage 1 -- arda cdr3fix
    "cdr3_aa", "cdr3_repaired", "v_call", "j_call", "locus", "species",
    "v_end", "j_start", "v_end_nt", "j_start_nt", "v_flags", "j_flags", "good", "fix_needed",
    # ``proposed`` names the side the submission left BLANK and the junction supplied, which is a
    # different fact from the `allele` flag (the submission named another allele of the same gene).
    # It is the column that lets a consumer drop its own segment proposer, so it has to reach the
    # frame -- it was in arda's markup and filtered out here, which no test noticed.
    "proposed", "n_errors", "errors",
    # stage 2 -- the most plausible nucleotide junction
    "cdr3_nt", "pgen", "scenario_p", "runner_up_pgen", "v_alts", "j_alts",
    "v_call_nt", "j_call_nt",
    # stage 3/4 -- the D gene from the model, placed by alignment, in nt and in aa
    "d_call", "d_posterior", "d_start_nt", "d_end_nt", "d_start_aa", "d_end_aa", "np1", "np2",
)

#: The D block of a row nothing could be said about. Dtypes are pinned here too: an all-null
#: column has no dtype of its own, and a consumer joining two corpora needs one that does not move.
_NULL_D: dict[str, object] = {"d_call": None, "d_posterior": None, "d_start_nt": None,
                              "d_end_nt": None, "d_start_aa": None, "d_end_aa": None,
                              "np1": None, "np2": None}
_D_DTYPES = {"d_call": pl.String, "d_posterior": pl.Float64, "d_start_nt": pl.Int64,
             "d_end_nt": pl.Int64, "d_start_aa": pl.Int64, "d_end_aa": pl.Int64,
             "np1": pl.String, "np2": pl.String}


#: Tables whose zeros are floored per ROW. Both are a choice of allele, and this stage pins the V
#: and the J, so their contribution to a row is a constant factor.
_FLOOR_ROWWISE = ("v_choice", "j_choice")
#: Tables whose zeros are floored only where the WHOLE conditional group is zero -- an allele the fit
#: never saw. A fitted profile is left exactly alone.
_FLOOR_GROUPWISE = ("v_3_del", "j_5_del")


def _reachable(model, eps: float = 1e-9):
    """A copy of ``model`` in which every shipped allele can actually be used.

    ⛔ Without this the pipeline answers for less than half of human TRB. The bundled models are EM
    fits, and an allele the training cohort never showed comes out at **p = 0** in two places at
    once: ``v_choice`` (36 of 66 TRB V alleles, including `TRBV2`, `TRBV13`, `TRBV17`, `TRBV27`) and
    ``v_3_del`` (34 alleles with an all-zero deletion profile). Conditioning on one of those returns
    no scenario at all, so every curated record naming it would get no nucleotide junction, and
    therefore no D.

    A zero there is a statement about the cohort that was fitted, not about what a rearrangement can
    be. Two different repairs, because the two zeros are not the same kind:

    * **Allele choice** is floored row by row. This stage is asked "given THIS V and THIS J, what is
      the most plausible nucleotide reading", so the usage term is the same constant for every
      candidate of a row: flooring and renormalising **cannot change which reading wins**, only
      whether one exists. (Which also means ``pgen`` from this model is not an absolute; read
      generation probabilities off the unfloored model with ``pgen_nt``.)
    * **V/J deletion profiles** are floored only where the allele's whole profile is zero, which
      makes it uniform -- the non-informative default, leaving the alignment to decide. An allele that
      HAS a profile keeps it untouched, because there the deletion count varies between candidates of
      one row and flooring an individual zero could change the winner.

    Nothing about the **D** is touched, for the same reason and one more: the D is not pinned, so
    flooring `d_gene` would let the reconstruction invent a D the fit never saw -- measured, it moved
    the nucleotide answer on records whose V and J were both well fitted. Stage 3 calls the D by
    alignment anyway, which is the whole point of inferring nucleotides.
    """
    import copy

    out = copy.deepcopy(model)
    for ev, t in list(out.tables.items()):
        if t is None or "p" not in t.columns:
            continue
        group = [c for c in t.columns if c != "p"][:-1]      # parents; the event's own key is last
        if ev in _FLOOR_ROWWISE:
            t = t.with_columns(pl.when(pl.col("p") <= 0.0).then(pl.lit(eps))
                               .otherwise(pl.col("p")).alias("p"))
        elif ev in _FLOOR_GROUPWISE:
            total = pl.col("p").sum().over(group) if group else pl.col("p").sum()
            t = t.with_columns(pl.when(total <= 0.0).then(pl.lit(eps))
                               .otherwise(pl.col("p")).alias("p"))
        else:
            continue
        norm = pl.col("p").sum().over(group) if group else pl.col("p").sum()
        out.tables[ev] = t.with_columns((pl.col("p") / norm).alias("p"))
    return out


def _aa_span(start_nt: int, end_nt: int) -> tuple[int | None, int | None]:
    """Residues whose codons the nucleotide span ``[start_nt, end_nt]`` (1-based closed) touches.

    A nucleotide at 1-based position ``p`` sits in residue ``(p - 1) // 3``, 0-based. So a D of one
    residue's worth of nucleotides straddling a codon boundary correctly reports a TWO-residue span:
    it is part of both codons, and that is the honest amino-acid answer to a nucleotide question.
    """
    if start_nt < 1 or end_nt < 1:
        return None, None
    return (start_nt - 1) // 3, (end_nt - 1) // 3


def annotate_junctions(junction_aas: Sequence[str], v_calls: Sequence[str],
                       j_calls: Sequence[str], *, species: str | Iterable[str] = "human",
                       model_source: str = "auto", n_best: int = 4, threads: int = 0,
                       d_k: int = 4) -> pl.DataFrame:
    """Run the whole junction pipeline. **One row out per row in, in input order.**

    Args:
        junction_aas: Junction amino acids, Cys104 through Phe/Trp118 **inclusive** -- VDJdb's
            ``cdr3`` convention, not arda's ``cdr3_aa``.
        v_calls, j_calls: Per-row V and J calls. Any spelling ``arda.cdr3fix.resolve_allele``
            accepts; a comma- or ``+``-joined list takes its leading entry. **A blank is allowed**
            -- arda proposes that side from the junction (the locus comes from the side that is
            named) and the ``proposed`` column says so. A call that is present but unresolvable is
            still refused, because naming something wrong is a defect and naming nothing is a gap.
        species: One name for every row, or one per row.
        model_source: Which bundled model set supplies the nucleotide guess and names the D.
            ``"auto"`` (default) is a **chain**, not a preference: OLGA's fit first, then arda's on
            whatever OLGA left unexplained. Measured on 4,000 real human rearrangements it beats
            either set alone on every axis -- TRB nucleotide-exact 14.40 / 17.32 / **17.32 %** and D
            gene 72.58 / 74.08 / **74.35 %** for arda / OLGA / the chain, and TRA keeps 4,000 of
            4,000 nucleotide junctions where OLGA alone declines 143. Name a set explicitly
            (``"olga"``, ``"arda"``, ``"learned"``) to pin one; only arda's covers mouse.
        n_best: Candidates re-scored per row in stage 2; see :func:`infer_nt_batch`, whose own
            default is 8. **4 here, measured**: stage 2 is 69 % of the pipeline's cost and this is
            the knob that moves it, and what the nucleotides are FOR is finding the D -- they only
            have to be right where the D is. On 4,000 real human TRB rearrangements, 4 runs at 249
            us/junction against 303 and calls the D at least as accurately (72.92 % against 72.08 %,
            within the sample's own noise); 8 buys nothing and 2 starts losing rows.
        threads: Kernel threads inside the native calls. ``0`` = auto. Do **not** wrap this in a
            pool of your own -- stage 2 already parallelises across the batch.
        d_k: Scenarios kept per row when naming the D (:func:`best_aa_scenarios_batch`). **4,
            measured**: accuracy FALLS as this rises -- 74.33 % at 4, 74.22 % at 8, 73.20 % at 16,
            72.95 % at 64 -- because the truncation is doing the regularising, and the low-weight
            scenarios it admits only dilute the winner. So the cheap setting is also the good one
            (30.1 us/junction at 4 against 50.1 at 64).

    Returns:
        A frame with :data:`JUNCTION_COLUMNS`. A row the model cannot explain carries nulls in the
        stages that failed and keeps stage 1 -- never dropped, never an exception, because a corpus
        legitimately contains species and loci with no shipped model.
    """
    from arda._markup import d_local_align
    from arda.annotate.dmap import _d_germlines
    from arda.cdr3fix import load_anchors as _anchors, markup_batch, resolve_species

    from .bundled import load_bundled
    from .native import best_aa_scenarios_batch, gene_to_allele
    from .viterbi import infer_nt_batch

    n = len(junction_aas)
    if not (len(v_calls) == len(j_calls) == n):
        raise ValueError(f"ragged input: {n} junctions, {len(v_calls)} V, {len(j_calls)} J")
    sp = [species] * n if isinstance(species, str) else list(species)
    if len(sp) != n:
        raise ValueError(f"ragged input: {n} junctions, {len(sp)} species")
    if not n:
        return pl.DataFrame(schema={c: pl.Null for c in JUNCTION_COLUMNS})

    # ---- stage 1: arda's cdr3fix, one call for the whole table.
    mk = markup_batch(pl.DataFrame({"cdr3": list(junction_aas), "v": list(v_calls),
                                    "j": list(j_calls), "species": sp}),
                      cdr3="cdr3", v="v", j="j", species="species")
    mk = mk.rename({"cdr3": "cdr3_aa"})

    # ---- stage 2: the most plausible nucleotide junction, one native call per (organism, locus).
    nt: list[str | None] = [None] * n
    pg: list[float | None] = [None] * n
    sc: list[float | None] = [None] * n
    ru: list[float | None] = [None] * n
    v_nt_call: list[str | None] = [None] * n
    j_nt_call: list[str | None] = [None] * n
    d_gene: list[str | None] = [None] * n
    d_post: list[float | None] = [None] * n
    org = [resolve_species(s) for s in sp]
    loci = mk["locus"].to_list()
    repaired = mk["cdr3_repaired"].to_list()
    v_res, j_res = mk["v_call"].to_list(), mk["j_call"].to_list()
    v_end_nt, j_start_nt = mk["v_end_nt"].to_list(), mk["j_start_nt"].to_list()
    v_alts = [(a.split(",") if a else [v_res[i]]) for i, a in enumerate(mk["v_alts"].to_list())]
    j_alts = [(a.split(",") if a else [j_res[i]]) for i, a in enumerate(mk["j_alts"].to_list())]
    groups: dict[tuple[str, str], list[int]] = {}
    for i, (o, loc) in enumerate(zip(org, loci)):
        if loc and repaired[i]:
            groups.setdefault((o, loc), []).append(i)
    # ⛔ The model SET is not a free choice, and ``auto`` is a measured chain rather than a
    # preference. On 4,000 real human rearrangements, same rows, only the set changed:
    #
    #                 TRB nt exact  TRB D right | TRA nt inferred  TRA nt exact
    #   arda               14.40 %      72.58 % |    4,000/4,000        29.45 %
    #   olga               17.32 %      74.08 % |    3,857/4,000        37.88 %
    #   olga -> arda       17.32 %      74.35 % |    4,000/4,000        38.90 %
    #
    # OLGA's fit is better calibrated and ~1.9x faster, but it DECLINES rows outright (143 of 4,000
    # on TRA) and has no mouse; arda's answers everything and is thinner -- 36 of its 66 human TRB V
    # alleles sit at p = 0, which is what `_reachable` exists to floor. ⚠ The flooring is NOT the
    # gap: unfloored, arda's set names a D on 2,669 of 4,000 rows at 48.25 % correct, so it is
    # load-bearing. Running one set and then the other on what it left empty beats both on every
    # axis, which is why the default is the chain and not either name.
    sources = ("olga", "arda") if model_source == "auto" else (model_source,)
    for (o, loc), group_idx in groups.items():
        todo = list(group_idx)
        for source in sources:
            if not todo:
                break
            idx = todo
            try:
                model = _reachable(load_bundled(loc, source=source, organism=o))
            except (FileNotFoundError, ValueError, KeyError):
                continue                 # this set has no model for this (organism, locus)
            # ⚠ arda's namespace is IMGT-complete and a model's is whatever it was fitted on, so a call
            # arda resolved -- `TRBV10-3*02`, after step one re-called it -- can name an allele the model
            # has never heard of, and `infer_nt_batch` refuses a name it does not know rather than
            # silently marginalising. Fall back to the model's representative allele for that GENE, and
            # only marginalise when the gene itself is absent, which keeps the curator's call in the
            # answer wherever the model can express it at all.
            rep = gene_to_allele(model)
            known = set(model.tables["v_choice"][model.tables["v_choice"].columns[0]].to_list())
            known |= set(model.tables["j_choice"][model.tables["j_choice"].columns[-2]].to_list())

            def _in_model(calls: Sequence[str]) -> list[str] | None:
                """Every one of arda's equally-good alleles the model can express, else its gene rep."""
                out: list[str] = []
                for call in calls:
                    hit = call if call in known else rep.get(call.split("*")[0])
                    if hit and hit not in out:
                        out.append(hit)
                return out or None            # None = marginalise; the gene itself is absent

            # ⛔ The whole tie set goes in, not one allele. An amino-acid junction frequently cannot
            # separate the alleles of a gene -- `CAISE` is TRBV10-3*01, *02 and *03 alike -- and arda
            # reports that as `v_alts` rather than resolving it by name order. `infer_nt_batch` scores a
            # LIST per row, so the stage that CAN separate them (codon plausibility, and the model's own
            # usage) is the stage that does.
            got = infer_nt_batch(model, [repaired[i] for i in idx],
                                 v=[_in_model(v_alts[i]) for i in idx],
                                 j=[_in_model(j_alts[i]) for i in idx], n_best=n_best, threads=threads)
            for k, i in enumerate(idx):
                row = got.row(k, named=True)
                nt[i], pg[i] = row.get("cdr3_nt"), row.get("pgen")
                sc[i], ru[i] = row.get("scenario_p"), row.get("runner_up_pgen")
                # Which allele of the tie set the nucleotides settled on -- reported, never written back
                # over arda's call: this is model evidence about an ambiguity, not a re-annotation.
                v_nt_call[i], j_nt_call[i] = row.get("v_call"), row.get("j_call")

            # ---- stage 3: WHICH D, from the same model, in the same group, one native call.
            #
            # `best_aa_scenarios_batch` returns the top-k scenarios per row with a weight and a D
            # allele on each -- the D axis of the same Pi_L*Pi_R transfer as `pgen_aa`. Summing the
            # weights by gene and normalising is the posterior, and it needs no nucleotides, no fitted
            # prior table and no per-locus tempering constant. It is here rather than in its own pass
            # because the model is loaded and the batch is already grouped.
            try:
                scen = best_aa_scenarios_batch(
                    model, [repaired[i] for i in idx],
                    v=[(_in_model(v_alts[i]) or [None])[0] for i in idx],
                    j=[(_in_model(j_alts[i]) or [None])[0] for i in idx],
                    k=d_k, threads=threads)
            except (KeyError, ValueError):
                continue                     # no D axis on a VJ locus; the rows keep their nulls
            if "d_call" not in scen.columns or not scen.height:
                continue
            # A VJ locus scores every scenario with no D at all, and polars gives an all-null column
            # of dtype Object there -- `.str` on which is a SchemaError, not an empty result.
            top = (scen.with_columns(pl.col("d_call").cast(pl.String, strict=False))
                       .with_columns(pl.col("d_call").str.split("*").list.first().alias("_g"))
                       .filter(pl.col("_g").is_not_null())
                       .group_by("row", "_g").agg(pl.col("w").sum())
                       .with_columns((pl.col("w") / pl.col("w").sum().over("row")).alias("_p"))
                       # Sort then take the head of each group, rather than `unique(keep="first")`:
                       # a hash-unique is not order-stable, and two D genes CAN carry the same summed
                       # weight -- `_g` is in the sort so that tie breaks by name and not by thread
                       # scheduling. Determinism is a requirement here, not a nicety.
                       .sort(["row", "_p", "_g"], descending=[False, True, False])
                       .group_by("row", maintain_order=True).first())
            for k_, g_, p_ in zip(top["row"].to_list(), top["_g"].to_list(), top["_p"].to_list()):
                d_gene[idx[k_]], d_post[idx[k_]] = g_, p_

            # Rows this set could not explain fall through to the next one in the chain.
            todo = [i for i in idx if nt[i] is None]

    # ---- stage 2b: the templated flanks are GERMLINE, not a guess.
    #
    # Only the N region is unknowable from amino acids. The reconstruction is free to differ from the
    # germline inside the V and J spans -- it is scored on codon plausibility, and it is conditioned
    # on whichever allele of the gene the MODEL carries, which need not be the allele arda called --
    # and measured against 3,000 real human TRB rearrangements it did: the V span disagreed with the
    # observed sequence on 23.93 % of rows and the J span on 21.83 %. Those nucleotides are not in
    # question, so they are written from the called allele's own germline. What that buys is not
    # tidiness: the residues at each boundary are part germline and part N region, and their codons
    # are exactly the evidence the nucleotide detour exists to recover, so the D search must see the
    # real bases there.
    anchors = None
    for i in range(n):
        if not nt[i] or v_end_nt[i] < 0 or j_start_nt[i] < 0:
            continue
        anchors = _anchors(org[i])
        va, ja = anchors.get(("V", v_res[i])), anchors.get(("J", j_res[i]))
        seq, L = nt[i], len(nt[i])
        ve, js = min(v_end_nt[i], L), min(j_start_nt[i], L)
        if va and va.germline_nt and len(va.germline_nt) >= ve:
            seq = va.germline_nt[:ve].upper() + seq[ve:]
        if ja and ja.germline_nt and len(ja.germline_nt) >= L - js:
            tail = L - js
            seq = seq[:js] + (ja.germline_nt[len(ja.germline_nt) - tail:].upper() if tail else "")
        nt[i] = seq

    # ---- stage 4: PLACE that gene, greedily, on the nucleotides.
    #
    # Naming and placing are different questions and this stage only answers the second one, so it
    # is ungated: the gene is already chosen, and refusing to say where it sits does not make the
    # name any better -- it just leaves a row with nothing to draw. Every allele of the gene is
    # tried and the best-scoring placement wins, which is what "greedy" means here.
    #
    # ⚠ The interior is `[v_end_nt, j_start_nt)` and those bounds are read off the nucleotide
    # junction itself -- stage 2b wrote the templated flanks from germline, so the germline prefix
    # and suffix of `cdr3_nt` ARE `v_end_nt` and `j_start_nt`. Re-deriving them by matching germline
    # against the inferred sequence was measured and returns the same numbers (v_end_nt exact on
    # 1.85 % of rows against 1.85 %, j_start_nt 78.12 % against 78.12 %), which is the check that
    # says stage 2b already did it.
    germ_by_gene: dict[tuple[str, str], dict[str, list[str]]] = {}

    def _germlines(organism: str, locus: str) -> dict[str, list[str]]:
        """``{gene: [allele sequences]}`` for one locus. Built once, not per row."""
        hit = germ_by_gene.get((organism, locus))
        if hit is None:
            hit = {}
            for allele, seq in _d_germlines(organism).get(locus, ()):  # type: ignore[misc]
                hit.setdefault(allele.split("*")[0], []).append(seq.upper())
            germ_by_gene[(organism, locus)] = hit
        return hit

    d_rows: list[dict] = []
    for i in range(n):
        seq, gene = nt[i], d_gene[i]
        if not seq or not gene or v_end_nt[i] < 0 or j_start_nt[i] < 0:
            d_rows.append(dict(_NULL_D) | ({"d_call": gene, "d_posterior": d_post[i]}
                                           if gene else {}))
            continue
        lo, hi = max(int(v_end_nt[i]), 0), min(int(j_start_nt[i]), len(seq))
        mid = seq[lo:hi] if hi > lo else ""
        alleles = _germlines(org[i], loci[i]).get(gene, ())
        if not mid or not alleles:
            d_rows.append(dict(_NULL_D) | {"d_call": gene, "d_posterior": d_post[i]})
            continue
        best = max((d_local_align(mid, g) for g in alleles), key=lambda r: r[0])
        s_nt, e_nt = lo + best[1] + 1, lo + best[2] + 1      # 1-based closed, as everywhere
        s_aa, e_aa = _aa_span(s_nt, e_nt)
        d_rows.append({
            "d_call": gene, "d_posterior": d_post[i],
            "d_start_nt": s_nt, "d_end_nt": e_nt, "d_start_aa": s_aa, "d_end_aa": e_aa,
            # The N regions fall out of the bounds, so they are sliced rather than re-derived --
            # one definition of where the D is, not two that can disagree in a drawing.
            "np1": seq[lo:s_nt - 1] or None, "np2": seq[e_nt:hi] or None,
        })

    out = mk.with_columns(
        pl.Series("cdr3_nt", nt, dtype=pl.String),
        pl.Series("pgen", pg, dtype=pl.Float64),
        pl.Series("scenario_p", sc, dtype=pl.Float64),
        pl.Series("runner_up_pgen", ru, dtype=pl.Float64),
        pl.Series("v_call_nt", v_nt_call, dtype=pl.String),
        pl.Series("j_call_nt", j_nt_call, dtype=pl.String),
    ).hstack(pl.DataFrame(d_rows, schema=_D_DTYPES))
    return out.select([c for c in JUNCTION_COLUMNS if c in out.columns])
