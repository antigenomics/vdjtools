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
**74.30 %** of all rows and **99.70 %** have coordinates. Letting the aligner choose the gene as well
reaches 47.93 % gated, 67.10 % ungated, and it is no better even where it is confident — 85.96 %
against the model's 86.32 % on those rows. So there is no `d_best`, no `d_posterior_call`, no
`d_entropy` and no `arda.dpost` port: a second D estimator was measured and is not needed.
**Requires `arda-mapper>=2.34.0`.**

**The nucleotides are the authority.** `v_end_nt`, `j_start_nt`, `d_start_nt` and `d_end_nt` are all
read off the inferred nucleotide junction, and every amino-acid bound is recomputed from them, so a
view showing both alphabets cannot draw them disagreeing.

**A blank V or J is handled**, not refused: arda proposes that side from the junction, and the
`proposed` column says which side was never curated. 3,130 of VDJdb's 192,726 curation keys leave a
side out, and they now get a nucleotide junction and a D like any other row.

**Fixed — `infer_nt_batch` raised on a per-row list of alleles**, one of the three call forms it
documents, when the boundary columns were computed.

**Fixed — a conditioned V could return nothing.** A bundled model is an EM fit, and an allele its
training cohort never showed comes out at probability zero — 36 of 66 human TRB V alleles — so
conditioning on one produced no scenario at all. The pipeline floors those before inferring;
verified across every allele the fit did see, no answer changes.

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
