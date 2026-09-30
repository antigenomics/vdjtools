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
`v_end_nt` / `j_start_nt` in nucleotides, `v_flags` / `j_flags`, `good`, and **`proposed`**.

**A blank V or J is answered, not refused.** A submission may leave one side out — 3,130 of VDJdb's
192,726 curation keys do — and arda proposes that side from the junction, the locus coming from the
side that *is* named. `proposed` says which side was never curated, which is a different fact from
the `allele` flag (the submission named a *different allele of the same gene*). It is the column
that lets a consumer delete its own segment proposer, so it is part of the contract.

⚠ A call that is present but **unresolvable** is still refused. `TRBVnope*01` keeps its
`FailedBadSegment` and `good = false`: naming something wrong is a defect a curator must see, naming
nothing is a gap the junction can fill.

### Stage 2 — vdjtools, `infer_nt_batch`

The most plausible nucleotide reading of the repaired junction under the recombination model,
conditioned on the confirmed V and J. Natively batched and threaded — one call for the whole table,
never a row loop. Out: `cdr3_nt` (length exactly `3 × len(cdr3_repaired)`), plus the scenario's own
`pgen`, `scenario_p` and `runner_up_pgen`, which is what "how plausible" means here.

The model is `load_bundled(source="arda", organism=…, locus=…)`, so the nucleotide guess and the
germline namespace agree by construction. Rows are grouped by `(organism, locus)` because a model
is per locus.

### Stage 3 — vdjtools, `best_aa_scenarios_batch`: WHICH D

The model's own top-`k` scenarios per row carry a weight and a D allele each — the D axis of the same
`Pi_L·Pi_R` transfer as `pgen_aa`. Summing those weights by gene and normalising **is** the posterior
over D genes. Same model as stage 2, already loaded, same `(organism, locus)` grouping, one native
threaded call. No fitted prior table, no per-locus tempering constant, no nucleotides.

`k = 4`, and that is not a speed compromise — accuracy **falls** as `k` rises (74.33 % at 4, 74.22 %
at 8, 73.20 % at 16, 72.95 % at 64, measured on one set). The truncation is doing the regularising; the low-weight
scenarios it admits only dilute the winner.

Out: `d_call` (a **gene** — an amino-acid junction does not identify a D allele) and `d_posterior`.

#### Why the aligner does not choose the gene

It was measured on 4,000 real human TRB rearrangements from `isalgo/airr_control`, whose D and
`DStart`/`DEnd` come from the nucleotides and are never shown to the pipeline:

| route | gene right, of all rows | has coordinates |
|---|---:|---:|
| E-value-gated alignment chooses and places | 47.93 % | 55.75 % |
| gated alignment, model posterior where it declines (the old `d_best`) | 71.40 % | 55.75 % |
| ungated alignment chooses and places | 67.10 % | 88.15 % |
| **model names, greedy alignment places** | **74.35 %** | **99.80 %** |

The gate is not what held the alignment back, and the alignment is not better where it *is*
confident: on the 2,230 rows it speaks for, the model's posterior is right 86.32 % against its
85.96 %, and they agree on 92.91 %. A hybrid — alignment where gated, model elsewhere — scores
**74.12 %**, below the model alone. So there is one gene estimator.

⛔ This is why `vdjtools.model.dpost` (the port of `arda.dpost`) does **not** ship. It is dominated
on both axes by a group-by over a call the pipeline already makes: 69.67 % against 74.33 % correct,
56.3 µs/junction against 30.1, and it needs a fitted `d_prior.tsv` and a per-locus `beta` that this
does not. `d_best` / `d_best_source` are gone with it — they existed only to reconcile two gene
estimators.

### Stage 3b — arda, `d_local_align`: WHERE

Greedy gapless local alignment of every allele of the **chosen gene** against the
`[v_end_nt, j_start_nt)` interior of the inferred nucleotide junction; best score wins. **Ungated**,
because the gene is already named and refusing to say where it sits does not improve the name — it
just leaves a row with nothing to draw.

⚠ The trade, stated: `d_start_nt` is exact on 64.71 % of correctly-called rows here against 66.67 %
under the E-value gate. But it is exact on **1,922 rows rather than 1,278**, because it answers 3,992
rather than 2,230. For a database drawing V/N/D/N/J that is the trade to take.

Out: `d_start_nt` / `d_end_nt` (1-based closed, junction space) and `np1` / `np2`, sliced from the
bounds rather than re-derived — one definition of where the D is, not two that can disagree in a
drawing.

**The D is allowed to overlap the V end and the J start**, and that is deliberate: an exonuclease
does not cut on a codon boundary, so the residue at each boundary is part germline and part
junction, and forcing the D strictly inside `[v_end, j_start]` would discard the flanking-codon
evidence that made the nucleotide detour worth taking.

Full D-D markup on IGH and TRD, with the E-value gate and the genomic-order mask, is still
`arda.annotate.dmap.map_d_junction`; this stage answers the single-D question a database draws.

#### ⚠ B cells: covered, unvalidated, and missing an SHM term

The bundled `arda` model set covers human **IGH**, IGK, IGL, TRA, TRB, TRD, TRG and mouse TRA/TRB,
so the D call reaches IGH — wider coverage than the `d_prior.tsv` route it replaces, which had
tables for human IGH/TRB/TRD and mouse TRB only. An IGH junction goes through end to end and comes
back with a gene, a posterior and coordinates.

⛔ **But every accuracy number on this page is human TRB.** `isalgo/airr_control` carries human and
mouse TRA/TRB and no immunoglobulin at all, so there is no B-cell truth set here and **74.35 % must
not be quoted for IGH.**

⛔ **And this pipeline has no somatic-hypermutation term anywhere in it.** Stage 2 reconstructs
nucleotides under a germline recombination model and stage 2b writes the templated flanks from
germline, so a hypermutated V tail is priced as insertion — which moves `v_end_nt`, which moves the
interior, which moves where the D can be placed. For a memory B cell that is the normal case, not a
corner case.

✅ The SHM-aware route exists and is deliberately kept: **`arda.hmm`** threads an
`arda.shmmodel.ShmModel` into the same semi-Markov recursion, so the templated V length is *priced*
rather than cut at the first mismatch. Measured on `IGHV3-30*18` / `IGHJ4*02`, one substitution in
the V tail gives `del_v >= len(v_nt) - 3` without a model and `del_v == 0` with one. It needs real
nucleotides, which is exactly the input this pipeline does not have — so the two are complements,
not alternatives, and neither replaces the other.

### Stage 4 — vdjtools, fold back to amino acids

A nucleotide position `p` (1-based) sits in residue `(p - 1) // 3` (0-based). So

    d_start_aa = (d_sequence_start - 1) // 3
    d_end_aa   = (d_sequence_end   - 1) // 3

which is the **span of residues whose codons the D touches** — the honest amino-acid answer to a
nucleotide question, and the reason a 1-residue D can still report a 2-residue span.

The **nucleotides are the authority** and every amino-acid bound is recomputed from them, so the two
alphabets cannot disagree in a drawing that shows both. That also settles what "recompute the bounds
from the most plausible nucleotide sequence" means for V and J: stage 2b writes the templated flanks
from germline, so the germline prefix and suffix of `cdr3_nt` **are** `v_end_nt` and `j_start_nt`.
Re-deriving them by matching germline against the inferred sequence was measured and returns the same
answer (`v_end_nt` exact on 1.85 % of rows either way, `j_start_nt` 78.12 % either way) — which is the
check that says stage 2b already did it.

## Cost, measured on the corpus it exists for

VDJdb's whole curation corpus — **192,726 distinct `(species, cdr3, V, J)` keys** — annotates in
**23.6 s in one process, 122 µs per key**. 192,423 keys get a nucleotide junction, **126,007 of the
126,231 D-bearing keys get a D gene**, and 122,717 of those get coordinates to draw; the 66,495 TRA
keys have no D to find.

That is 1.57× faster than pinning arda's set (35.1 s, 182 µs) as well as more accurate, because
OLGA's fit answers most rows and is the cheaper one to score. On the human TRB benchmark alone the
figure is 168 µs/junction.

| stage | µs/junction |
|---|---:|
| 1 — arda `markup_batch` | ~54 |
| 2 — `infer_nt_batch` + the germline flank splice | the remainder, and the bulk of it |
| 3 — `best_aa_scenarios_batch`, naming the D | ~30 |
| 3b — `d_local_align`, placing it | ~1.4 |

⛔ **Do not wrap this in a pool.** Every stage is already batched: one `markup_batch` for the whole
table, then one native threaded `infer_nt_batch` and one `best_aa_scenarios_batch` per
`(organism, locus)` on a model that is loaded once, then a thin Python loop over arda's C++ aligner.
A pool around it would re-import both libraries and re-load every model per worker, which is the
failure this project keeps writing down.

## Which model set — `model_source="auto"` is a chain, not a preference

Both bundled fits were scored on the same 4,000 real human rearrangements, changing nothing but the
set:

| | TRB nt exact | TRB D gene right | TRA nt inferred | TRA nt exact | µs/junction |
|---|---:|---:|---:|---:|---:|
| `"arda"` | 14.40 % | 72.58 % | 4,000/4,000 | 29.45 % | 272 |
| `"olga"` | 17.32 % | 74.08 % | 3,857/4,000 | 37.88 % | 141 |
| **`"auto"` — OLGA, then arda on what it left** | **17.32 %** | **74.35 %** | **4,000/4,000** | **38.90 %** | **168** |

OLGA's fit is better calibrated and cheaper to score; arda's answers everything and is the only set
**vendored here** for mouse. Neither dominates, so the default runs one and then the other on the
rows the first left empty — which does dominate, on every column above.

⚠ **OLGA itself is not human-only; this bundle is.** OLGA ships `mouse_T_alpha`, `mouse_T_beta`,
`mouse_B_heavy`, `mouse_B_kappa` and `mouse_B_lambda`, and all five import cleanly through
`from_olga` — they are simply not vendored yet. On mouse TRB the vendored `arda` fit wins anyway
(D gene right on **53.97 %** of 4,000 real mouse rearrangements against OLGA's **52.33 %**), so what
vendoring them would buy is mouse IGH/IGK/IGL, which no bundled set covers today.

**And the chain does not end at a bundled set.** A fitted model exists for human and mouse and for
nothing else, so every other organism arda ships germline for used to come back empty — 1,457 rhesus
keys in VDJdb answered zero times. `"auto"`'s last rung is a germline **scaffold**
(`from_arda(locus, organism=)`): the templated flanks come from that organism's own germline and
only the N region is a default. Coverage over VDJdb's 192,726 curation keys, per species:

| species | locus | keys | nucleotide junction | D gene | D coordinates |
|---|---|---:|---:|---:|---:|
| human | TRB | 115,829 | **115,819** | 115,819 | 113,232 |
| human | TRA | 58,281 | **58,203** | — | — |
| mouse | TRB | 9,015 | **8,805** | 8,805 | 8,694 |
| mouse | TRA | 8,140 | **8,140** | — | — |
| rhesus | TRB | 1,383 | **1,379** | 1,379 | 787 |
| rhesus | TRA | 74 | **73** | — | — |
| human | TRD | 4 | **4** | 4 | 4 |

**192,423 of 192,726 (99.84 %)**, against 124,022 before, with no species or locus at zero. TRA and
TRD are VJ loci here — there is no D to find. ⚠ The `keys` column is as of this release: the 461
keys naming neither V nor J had no locus before and now count under the one they resolve to, so a
per-locus denominator is not comparable to an older run's.

⚠ **The flooring is not what makes arda's set weaker, and it is load-bearing.** 36 of its 66 human
TRB V alleles sit at probability 0, so `_reachable` has to floor them before anything conditions;
*unfloored*, that set names a D on 2,669 of 4,000 rows at 48.25 % correct against the floored
72.58 %. Improving the bundled fit is the real fix and it belongs in whatever produces that table,
not here.

Pin a set by name (`"olga"`, `"arda"`, `"learned"`) when you need one model's numbers and not the
best available answer — a reproducibility run, or a comparison against a published fit.

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
