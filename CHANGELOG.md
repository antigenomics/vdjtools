# Changelog

What changed for you, per release. Anything not listed is internal.
Full release notes: <https://github.com/antigenomics/vdjtools/releases>.

## 4.8.1

**Fixed — `annotate_junctions` no longer alters the junction it was asked about** (#186, #187).
`cdr3_nt` is the amino-acid junction's nucleotides, always: `translate(cdr3_nt) == cdr3_repaired` on
every row. The templated flanks were written from germline up to the nucleotide boundary, and where
that boundary sits a residue past the amino-acid one the germline overwrote a residue. On 3,000 TRB
junctions drawn from the bundled model, 20 (0.67 %) came back with a substituted residue; now 0 of
3,000. A flank is taken from germline only where the translation survives; the V/J boundaries and the
D markup stay coordinates and never edit the sequence.

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
**Requires `arda-mapper>=2.36.0`.**

**The nucleotides are the authority.** `v_end_nt`, `j_start_nt`, `d_start_nt` and `d_end_nt` are all
read off the inferred nucleotide junction, and every amino-acid bound is recomputed from them, so a
view showing both alphabets cannot draw them disagreeing.

**New — `model_source="auto"`, the default, is a measured chain rather than a preference.** OLGA's
bundled fit is better calibrated and cheaper to score; arda's answers every row and is the only one
covering mouse. Running OLGA's first and arda's on what it left empty beats both on every axis:
human TRB nucleotide-exact 14.40 / 17.32 / **17.32 %** and D gene 72.58 / 74.08 / **74.35 %** for
arda / OLGA / the chain, and human TRA keeps **4,000 of 4,000** nucleotide junctions where OLGA
alone declines 143. Pin a set by name for a reproducibility run.

**122 µs per key on the real corpus** — VDJdb's 192,726 distinct `(species, cdr3, V, J)` keys in
**23.6 s in one process**: 192,423 nucleotide junctions, a D gene on 126,007 of the 126,231 D-bearing
keys, and coordinates on 122,717 of those. The 66,495 TRA keys have no D to find. Every stage
is batched (one `markup_batch`, then one native threaded `infer_nt_batch` and one
`best_aa_scenarios_batch` per organism, locus and model set, on a model loaded once), so do **not**
wrap the call in a pool of your own.

**A blank V or J is handled**, not refused: arda proposes that side from the junction, and the
`proposed` column says which side was never curated. 3,130 of VDJdb's 192,726 curation keys leave a
side out — 2,669 name one side, **461 name neither** — and all of them now get a nucleotide junction
and a D like any other row. For the 461 the locus is proposed too, and it agrees with the
`cdr3.alpha` / `cdr3.beta` column the record was filed under on **457 of 461 (99.13 %)**.

**Fixed — `infer_nt_batch` raised on a per-row list of alleles**, one of the three call forms it
documents, when the boundary columns were computed.

**Fixed — a conditioned V could return nothing.** A bundled model is an EM fit, and an allele its
training cohort never showed comes out at probability zero — 36 of 66 human TRB V alleles — so
conditioning on one produced no scenario at all. The pipeline floors those before inferring;
verified across every allele the fit did see, no answer changes.

**Every species in the corpus answers now, not just the fitted ones.** A fitted model exists for
human (seven loci) and mouse (TRA/TRB) and for nothing else, so every other organism arda ships
germline for came back with no nucleotides, no D and no bounds — **1,457 rhesus keys in VDJdb
answered zero times**, invisible inside a single corpus-wide total. `model_source="auto"` now ends
in a germline **scaffold** built by `from_arda(locus, organism=)`: the templated flanks come from
that organism's own germline and only the N region is a default, which is exactly the part a fitted
model would have improved. Coverage over VDJdb's own 192,726 curation keys, per species:

| species | locus | keys | nucleotide junction | D gene |
|---|---|---:|---:|---:|
| human | TRB | 115,829 | **115,819** | 115,819 |
| human | TRA | 58,281 | **58,203** | — |
| mouse | TRB | 9,015 | **8,805** | 8,805 |
| mouse | TRA | 8,140 | **8,140** | — |
| rhesus | TRB | 1,383 | **1,379** | 1,379 |
| rhesus | TRA | 74 | **73** | — |
| human | TRD | 4 | **4** | 4 |

**192,423 of 192,726 keys (99.84 %)**, against 124,022 before, and no species or locus at zero. TRA
and TRD are VJ loci for the D column's purposes — there is no D to find. ⚠ Read the `keys` column as
of this release: the 461 keys that named neither V nor J had no locus at all before and are counted
here under the locus they resolve to, so a per-locus denominator is not comparable to an older run's.

**Fixed — `infer_nt_batch` raised on a batch whose FIRST row had no D.** The native call returns its
string columns as numpy **object** arrays and polars types one from its first element alone, so
`[None, "TRBD1*01"]` came back `Object` and the declared `Utf8` cast raised `cannot cast 'Object'
type` — while `["TRBD1*01", None]` was fine. The failure therefore depended on **row order**: the
same corpus sorted differently either worked or raised.

**⚠ OLGA is not human-only — this bundle is.** `load_bundled` used to say the `olga` set covers
human only. OLGA itself ships five mouse models (`mouse_T_alpha`, `mouse_T_beta`, `mouse_B_heavy`,
`mouse_B_kappa`, `mouse_B_lambda`) and all five import cleanly through `from_olga`; they are simply
not vendored here yet. On mouse TRB the vendored `arda` fit is the better one anyway — the D gene is
right on **53.97 %** of 4,000 real mouse rearrangements against OLGA's **52.33 %** — so what
vendoring them would buy is mouse IGH/IGK/IGL, which no bundled set covers today.

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
