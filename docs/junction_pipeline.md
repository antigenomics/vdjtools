# The junction pipeline: cdr3 fix → V/J markup → nucleotide guess → D markup

**Status: implemented in vdjtools 4.8.0 as `vdjtools.model.annotate_junctions`.**
Specified by M.S. on 2026-09-30. Every stage is arda's or vdjtools'; a consumer
(`antigenomics/vdjdb-db`) calls **one** function and keeps no annotation code of its own.

## Why it is one pipeline and not four

A VDJdb-style record is `(junction_aa, V call, J call, species)` and nothing else. Every downstream
question — is the junction well formed, where does the V stop, which D, where does the D sit — is
answered off that, and each answer is an *input to the next*:

- the D can only be looked for **between** the V and J boundaries, so the boundaries come first;
- the boundaries are only meaningful once the junction and its calls are repaired, so the repair
  comes before them;
- and a D that contributes **one amino acid** is invisible in the translated junction, so the search
  has to happen in **nucleotides** — which a record without nucleotides does not have, so they have
  to be inferred in between.

That last step is the whole reason the pipeline exists. A 1-residue D shows up in the *flanking
codons*: the residue before and after the D are each part germline-D and part N region, and only
nucleotides can say which bases those are. Translating first throws exactly that away.

## The four stages

| # | stage | library | entry point |
|---|---|---|---|
| 1 | confirm-or-replace the V/J call, repair the junction, place the aa boundaries | **arda** | `arda.cdr3fix.markup_batch` |
| 2 | most plausible nucleotide junction for that amino-acid junction | **vdjtools** | `model.infer_nt_batch` |
| 3 | D (and tandem D2) by gapless local alignment **in nucleotide space** | **arda** | `arda.annotate.dmap.map_d_junction` |
| 4 | fold the nt D coordinates back onto amino-acid positions | **vdjtools** | `model.junction` |

### Stage 1 — arda, `cdr3fix`

Re-calls the V or J when another allele of the locus explains **2 more residues contiguously from
its own anchor**, then trims framework outside the anchors, restores germline the submission was
cut inside of, and substitutes the first/last residue to the anchor on a solid match. Nothing else
is rewritten: a residue that disagrees with germline inside the templated run is flagged
(`mismatch`) and left alone, because a curation error and an allele IMGT does not record are
indistinguishable from one junction.

Out: `cdr3_repaired`, `v_call`, `j_call` (possibly re-called), `v_end` / `j_start` in residues,
`v_end_nt` / `j_start_nt` in nucleotides, `v_flags` / `j_flags`, `good`.

### Stage 2 — vdjtools, `infer_nt_batch`

The most plausible nucleotide reading of the repaired junction under the recombination model,
conditioned on the confirmed V and J. Natively batched and threaded — one call for the whole table,
never a row loop. Out: `cdr3_nt` (length exactly `3 × len(cdr3_repaired)`), plus the scenario's own
`pgen`, `scenario_p` and `runner_up_pgen`, which is what "how plausible" means here.

The model is `load_bundled(source="arda", organism=…, locus=…)`, so the nucleotide guess and the
germline namespace agree by construction. Rows are grouped by `(organism, locus)` because a model
is per locus.

### Stage 3 — arda, `map_d_junction`

Gapless local alignment of every D germline of the locus against the **V..J interior of the
inferred nucleotide junction**, with the Karlin–Altschul E-value gate and the genomic-order mask
(TRBD2 cannot reach the TRBJ1 cluster). Tandem D-D is searched on IGH and TRD. No mmseqs pass: the
boundaries from stage 1 give the interior directly.

Out: `d_call`, `d_sequence_start` / `d_sequence_end` (1-based closed, junction space), `d2_*`,
`np1` / `np2` / `np3`.

**The D is allowed to overlap the V end and the J start**, and that is deliberate: an exonuclease
does not cut on a codon boundary, so the residue at each boundary is part germline and part
junction, and forcing the D strictly inside `[v_end, j_start]` would discard the flanking-codon
evidence that made the nucleotide detour worth taking.

### Stage 4 — vdjtools, fold back to amino acids

A nucleotide position `p` (1-based) sits in residue `(p - 1) // 3` (0-based). So

    d_start_aa = (d_sequence_start - 1) // 3
    d_end_aa   = (d_sequence_end   - 1) // 3

which is the **span of residues whose codons the D touches** — the honest amino-acid answer to a
nucleotide question, and the reason a 1-residue D can still report a 2-residue span.

`vdjtools.model.posterior_d_batch` stays as the independent, model-only answer to "which D": it
marginalises the insertion-length and D-trimming distributions and needs no nucleotides at all. The
pipeline reports it beside the alignment call, so a consumer can see when the two disagree rather
than having to pick one.

## Contract

`annotate_junctions(junction_aas, v_calls, j_calls, species=…)` returns a `polars.DataFrame` with
**one row per input row, in input order**, so it joins positionally. A record the model cannot
explain is present with nulls, never dropped and never an exception: a real corpus contains species
and loci with no shipped model.

## What the consumer keeps

Nothing. `vdjdb-db`'s `annotate/cdr3fix.py`, `annotate/dgene.py`, `annotate/junction.py` and
`curate/anchors.py` collapse to one call plus a rename of the output columns into VDJdb's own
names. The build keeps the *curation* decisions (which records to flag, what to do about a
contradicted call) and no *annotation* code.
