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
N region is inferred), arda aligns the D against those nucleotides, and the nucleotide D coordinates
fold onto the residues whose codons they touch. That last step is the point: a D contributing one
amino acid is invisible in a translated junction but not in nucleotides, because the residues at
each end of it are part germline and part N region.

Read **`d_best`** for the D. Measured on 4,000 real human TRB rearrangements whose D is called from
the sequence: 71.5 % correct over all rows, against 48.3 % for the alignment alone (which answers
56 % of rows, at 85.98 % accuracy) and 69.7 % for the posterior alone. 307 µs per junction.

**New — `model.posterior_d`, `posterior_d_batch`, `load_d_prior`**, moved here from `arda.dpost`
(arda 2.33.0 removes them). Ported, not rewritten: the answers are unchanged. The batch form is new
and 1.7× faster per key than the loop it replaces. **Requires `arda-mapper>=2.34.0`.**

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
