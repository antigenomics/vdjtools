# Changelog

What changed for you, per release. Anything not listed is internal.
Full release notes: <https://github.com/antigenomics/vdjtools/releases>.

## 4.8.0

**New — the junction pipeline.** `model.annotate_junctions(cdr3_aas, v, j, species=)` takes a bare
amino-acid junction with its V and J calls — all a curated database has — and returns a repaired
junction, confirmed or corrected gene calls, V/J boundaries in residues and nucleotides, the most
plausible nucleotide junction, and a D with coordinates in both alphabets. One row out per input
row, in input order.

Four stages, each an existing entry point: arda repairs the junction and places the boundaries,
`infer_nt_batch` infers the nucleotides (the templated flanks are written from germline, so only the
N region is inferred), the model names the D gene, and the aligner places it on those nucleotides —
whose coordinates then fold onto the residues their codons touch. That last step is the point: a D
contributing one amino acid is invisible in a translated junction but not in nucleotides, because
the residues at each end of it are part germline and part N region.

**Naming the D and placing it are separate questions, and one estimator answers each.** The model
names the gene (`d_call`, with `d_posterior`) from its own scenario weights — no second model, no
fitted prior table. The aligner then places that gene greedily and **ungated** (`d_start_nt` /
`d_end_nt`, with `np1` / `np2` either side), because the gene is already chosen and refusing to say
where it sits only leaves a row with nothing to draw.

On 4,000 real human TRB rearrangements whose D comes from the nucleotides: the gene is right on
**74.35 %** of all rows and **99.80 %** have coordinates. Letting the aligner choose the gene as well
reaches 47.93 % gated, 67.10 % ungated, and it is no better even where it is confident — 85.96 %
against the model's 86.32 % on those rows. So there is no `d_best`, no `d_posterior_call`, no
`d_entropy` and no `arda.dpost` port: a second D estimator was measured and is not needed.
**Requires `arda-mapper>=2.34.0`.**

**The nucleotides are the authority.** `v_end_nt`, `j_start_nt`, `d_start_nt` and `d_end_nt` are all
read off the inferred nucleotide junction, and every amino-acid bound is recomputed from them, so a
view showing both alphabets cannot draw them disagreeing.

**New — `model_source="auto"`, the default, is a measured chain rather than a preference.** OLGA's
bundled fit is better calibrated and cheaper to score; arda's answers every row and is the only one
covering mouse. Running OLGA's first and arda's on what it left empty beats both on every axis:
human TRB nucleotide-exact 14.40 / 17.32 / **17.32 %** and D gene 72.58 / 74.08 / **74.35 %** for
arda / OLGA / the chain, and human TRA keeps **4,000 of 4,000** nucleotide junctions where OLGA
alone declines 143. Pin a set by name for a reproducibility run.

**116 µs per key on the real corpus** — VDJdb's 192,726 distinct `(species, cdr3, V, J)` keys in
**22.3 s in one process**: 190,199 nucleotide junctions, a D gene on 124,022 of the 124,489 TRB keys,
and coordinates on 121,336 of those. The 65,107 TRA keys have no D to find. Every stage
is batched (one `markup_batch`, then one native threaded `infer_nt_batch` and one
`best_aa_scenarios_batch` per organism, locus and model set, on a model loaded once), so do **not**
wrap the call in a pool of your own.

**A blank V or J is handled**, not refused: arda proposes that side from the junction, and the
`proposed` column says which side was never curated. 3,130 of VDJdb's 192,726 curation keys leave a
side out, and they now get a nucleotide junction and a D like any other row.

**Fixed — `infer_nt_batch` raised on a per-row list of alleles**, one of the three call forms it
documents, when the boundary columns were computed.

**Fixed — a conditioned V could return nothing.** A bundled model is an EM fit, and an allele its
training cohort never showed comes out at probability zero — 36 of 66 human TRB V alleles — so
conditioning on one produced no scenario at all. The pipeline floors those before inferring;
verified across every allele the fit did see, no answer changes.

**⚠ What is not measured, for B cells.** Every figure above is human TRB: the truth set this is
scored against carries TRA and TRB and no immunoglobulin, so **none of them may be quoted for IGH**.
The D call does reach IGH, on a wider model set than the prior tables it replaces — but this pipeline
has **no somatic-hypermutation term**, so a hypermutated V tail is priced as insertion, which moves
`v_end_nt` and therefore where the D can sit. For B cells with real nucleotides that is
`arda.hmm(..., shm=ShmModel)`, which prices the mutated tail instead of cutting at the first
mismatch — a complement to this, not an alternative.

## 4.7.0

`model.germline_boundary(model, aas, v=, j=)` — where the V and J germlines stop, in CDR3 nucleotide
space, carried by `infer_nt_batch` as `v_end_germline` / `j_start_germline`. Read this for a
boundary and `Scenario.v_end` for a recombination history: maximising `P(sequence)` explains N-region
nucleotides as templated, so a history walks the boundary outward. `v.end` exact goes 71.8 % →
**92.9 %**.

## 4.6.0 – 4.6.1

Batched translation, and allele suffixes are no longer stripped per row. 4.6.1 is the one that
reached PyPI.

## 4.5.0

`infer_nt_batch` — `infer_nt` over a whole table, natively batched instead of a per-row Python
wrapper, which was 86 % of a consuming build at ~200,000 junctions.

## Earlier

3.x and 4.0 – 4.4 are described in their git tags (`git tag --sort=-v:refname`) and in the GitHub
release notes; v2.x is in the `legacy-1.x` history.
