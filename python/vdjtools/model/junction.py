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
3. **arda** ``annotate.dmap.map_d_junction`` finds the D -- and a tandem D2 on IGH/TRD -- by gapless
   local alignment against the V..J interior of that nucleotide junction, under the E-value gate and
   the genomic-order mask.
4. here: the nucleotide D coordinates are folded back onto amino-acid positions.

**Why the nucleotide detour is the point.** A D that contributes one amino acid is invisible in a
translated junction, but it is not invisible in nucleotides: the residues at each end of the D are
part germline-D and part N region, so their codons carry D bases. Searching the translated junction
throws that away; searching the inferred nucleotide junction keeps it, and the D is then allowed to
OVERLAP ``v_end`` / ``j_start`` for exactly the same reason -- an exonuclease does not cut on a
codon boundary.

The alignment call is reported beside :func:`~vdjtools.model.dpost.posterior_d_batch`, the
model-only answer that needs no nucleotides at all, and ``d_best`` is the two combined -- the
alignment where it speaks, the posterior where it declines. Measured against real nucleotide D calls
that is **73.05 %** correct where the alignment alone reaches 48.77 % (it answers 57.2 % of rows)
and the posterior alone 71.15 %, so neither route replaces the other and both ship.
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
    "n_errors", "errors",
    # stage 2 -- the most plausible nucleotide junction
    "cdr3_nt", "pgen", "scenario_p", "runner_up_pgen", "v_alts", "j_alts",
    "v_call_nt", "j_call_nt",
    # stage 3/4 -- D by nucleotide alignment, in nt and in aa
    "d_call", "d_start_nt", "d_end_nt", "d_start_aa", "d_end_aa", "d_support",
    "d2_call", "d2_start_nt", "d2_end_nt", "np1", "np2", "np3",
    # the model-only second opinion, and the two combined
    "d_posterior_call", "d_posterior", "d_entropy", "d_best", "d_best_source",
)

_NULL_D = {"d_call": None, "d_start_nt": None, "d_end_nt": None, "d_start_aa": None,
           "d_end_aa": None, "d_support": None, "d2_call": None, "d2_start_nt": None,
           "d2_end_nt": None, "np1": None, "np2": None, "np3": None}


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
                       model_source: str = "arda", n_best: int = 4, threads: int = 0,
                       d_max_evalue: float | None = None,
                       posterior: bool = True) -> pl.DataFrame:
    """Run the whole junction pipeline. **One row out per row in, in input order.**

    Args:
        junction_aas: Junction amino acids, Cys104 through Phe/Trp118 **inclusive** -- VDJdb's
            ``cdr3`` convention, not arda's ``cdr3_aa``.
        v_calls, j_calls: Per-row V and J calls. Any spelling ``arda.cdr3fix.resolve_allele``
            accepts; a comma- or ``+``-joined list takes its leading entry.
        species: One name for every row, or one per row.
        model_source: Which bundled model set supplies the nucleotide guess. ``"arda"`` is the
            default because it is the only set covering mouse and it shares arda's allele namespace.
        n_best: Candidates re-scored per row in stage 2; see :func:`infer_nt_batch`, whose own
            default is 8. **4 here, measured**: stage 2 is 69 % of the pipeline's cost and this is
            the knob that moves it, and what the nucleotides are FOR is finding the D -- they only
            have to be right where the D is. On 4,000 real human TRB rearrangements, 4 runs at 249
            us/junction against 303 and calls the D at least as accurately (72.92 % against 72.08 %,
            within the sample's own noise); 8 buys nothing and 2 starts losing rows.
        threads: Kernel threads inside the native calls. ``0`` = auto. Do **not** wrap this in a
            pool of your own -- stage 2 already parallelises across the batch.
        d_max_evalue: Override arda's shipped E-value gate on the D call.
        posterior: Also compute the model-only D posterior (:func:`posterior_d_batch`). It is a
            second opinion on ``d_call`` and costs about as much as the rest of the pipeline.

    Returns:
        A frame with :data:`JUNCTION_COLUMNS`. A row the model cannot explain carries nulls in the
        stages that failed and keeps stage 1 -- never dropped, never an exception, because a corpus
        legitimately contains species and loci with no shipped model.
    """
    from arda.annotate.dmap import map_d_junction
    from arda.cdr3fix import load_anchors as _anchors, markup_batch, resolve_species

    from .bundled import load_bundled
    from .native import gene_to_allele
    from .dpost import posterior_d_batch
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
    for (o, loc), idx in groups.items():
        try:
            model = _reachable(load_bundled(loc, source=model_source, organism=o))
        except (FileNotFoundError, ValueError, KeyError):
            continue                     # no bundled model for this (organism, locus)
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
    anchors = load_anchors_for = None
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

    # ---- stage 3 + 4: D by nucleotide alignment, then folded back onto residues.
    #
    # Per row, because the D search is arda's C++ gapless local aligner behind a per-record entry
    # point and the loop around it is not the cost -- the alignment is, and it is already native.
    d_rows: list[dict] = []
    for i in range(n):
        if not nt[i] or v_end_nt[i] < 0 or j_start_nt[i] < 0:
            d_rows.append(dict(_NULL_D))
            continue
        # The interior comes from STAGE ONE, not from re-matching germline against an inferred
        # sequence: the nucleotides were reconstructed under whichever allele of the gene the model
        # carries, so an exact prefix match against arda's called allele breaks at the first
        # synonymous difference and opens the window inside the V.
        call = map_d_junction(nt[i], v_res[i], j_res[i], org[i], d_max_evalue=d_max_evalue,
                              v_end=v_end_nt[i], j_start=j_start_nt[i])
        s_aa, e_aa = _aa_span(call.d_sequence_start, call.d_sequence_end)
        d_rows.append({
            "d_call": call.d_call or None,
            "d_start_nt": call.d_sequence_start if call.d_call else None,
            "d_end_nt": call.d_sequence_end if call.d_call else None,
            "d_start_aa": s_aa if call.d_call else None,
            "d_end_aa": e_aa if call.d_call else None,
            "d_support": call.d_support or None,
            "d2_call": call.d2_call or None,
            "d2_start_nt": call.d2_sequence_start if call.d2_call else None,
            "d2_end_nt": call.d2_sequence_end if call.d2_call else None,
            "np1": call.np1 or None, "np2": call.np2 or None, "np3": call.np3 or None,
        })

    # ---- the model-only second opinion on which D.
    if posterior:
        post = posterior_d_batch(repaired, v_res, j_res, species=sp)
    else:
        post = [None] * n

    # ---- the answer to use: the alignment where it speaks, the posterior where it declines.
    #
    # Measured on 4,000 real human TRB rearrangements from `isalgo/airr_control`, scored against
    # their own nucleotide D call (`appendix/bench_junction_pipeline.py`):
    #
    #   alignment on the inferred nt  called 57.2 %, correct 85.31 % of called, 48.77 % of all
    #   length-and-prior posterior    called  100 %, correct 71.17 % of called, 71.15 % of all
    #   this combination              called  100 %,                            73.05 % of all
    #
    # So neither replaces the other and the order is not arbitrary. Where both answer the alignment
    # is the better one (1,951 correct against 1,875 of the same 2,287), and on the 1,712 it declines
    # the posterior is still right 56.72 % of the time -- which is the whole reason the probabilistic
    # route survives now that the nucleotide route exists.
    d_best, d_src = [], []
    for row, p in zip(d_rows, post):
        if row["d_call"]:
            d_best.append(row["d_call"].split("*")[0])
            d_src.append("alignment")
        elif p is not None:
            d_best.append(p.d_call)
            d_src.append("posterior")
        else:
            d_best.append(None)
            d_src.append(None)

    out = mk.with_columns(
        pl.Series("cdr3_nt", nt, dtype=pl.String),
        pl.Series("pgen", pg, dtype=pl.Float64),
        pl.Series("scenario_p", sc, dtype=pl.Float64),
        pl.Series("runner_up_pgen", ru, dtype=pl.Float64),
        pl.Series("v_call_nt", v_nt_call, dtype=pl.String),
        pl.Series("j_call_nt", j_nt_call, dtype=pl.String),
        pl.Series("d_posterior_call", [p.d_call if p else None for p in post], dtype=pl.String),
        pl.Series("d_posterior", [p.posterior if p else None for p in post], dtype=pl.Float64),
        pl.Series("d_entropy", [p.entropy if p else None for p in post], dtype=pl.Float64),
        pl.Series("d_best", d_best, dtype=pl.String),
        pl.Series("d_best_source", d_src, dtype=pl.String),
    ).hstack(pl.DataFrame(d_rows, schema={k: (pl.String if k in
                                              ("d_call", "d_support", "d2_call", "np1", "np2",
                                               "np3") else pl.Int64) for k in _NULL_D}))
    return out.select([c for c in JUNCTION_COLUMNS if c in out.columns])
