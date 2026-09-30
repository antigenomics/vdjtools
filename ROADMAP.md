# Roadmap

What vdjtools does today, and what is planned. The maintainer's working backlog — measured dead
ends, per-item baselines, open questions — is kept out of the published tree.

## What works

**A generative model of V(D)J recombination.** Native C++ Pgen on nucleotide and amino-acid
sequences, exact against OLGA on all seven human and mouse loci, with D-D fusions kept exact. Model
fitting by EM, and a generation sampler that reaches 2.3 M productive sequences/s on IGH.

**The junction pipeline** (`model.annotate_junctions`, `docs/junction_pipeline.md`). From a bare
amino-acid junction with its V and J calls: a repaired junction, confirmed or corrected gene calls,
V/J boundaries in residues and nucleotides, the most plausible nucleotide junction, and a D with
coordinates in both alphabets. **182 µs per junction** — VDJdb's whole 192,726-key
corpus in 35 s, one process.

**Repertoire statistics.** Diversity (including iNEXT and rarefaction), spectratype, V/J usage,
functional summaries, and hypermutation; pairwise and cohort overlap; decontamination; a
hive-partitioned cohort store.

**Single-cell interop.** One flat AIRR Rearrangement table with `cell_id`, plus adapters for scirpy,
dandelion and scRepertoire.

## Where the boundaries are

vdjtools is **the recombination model and everything probabilistic**. The germline reference, the
alignment against it, and junction repair belong to
[arda](https://github.com/antigenomics/arda), which vdjtools takes as a base dependency and calls
for markup and for the D alignment. `model.annotate_junctions` is the seam: one call that runs both
libraries in the order that works.

## Planned

- **A better nucleotide guess.** The inferred junction is exact on 14.30 % of real human TRB
  rearrangements, with a median of 2 mismatching nucleotides in 42 — all in the N region, now that
  the templated flanks are written from germline. The insertion models are Markov and the
  reconstruction takes the argmax per position, so there is room there.
- **A swept D-alignment gate.** The alignment declines 43.8 % of rows and is right on 85.98 % where
  it speaks; the gate has never been tuned on an *inferred* sequence, which is noisier than a read.
- **The same benchmark beyond human TRB** — mouse TRB first, where the bundled model is the only
  non-human fit.
- **An unfloored `pgen`** beside the scenario, so a caller reading generation probabilities does not
  have to know that `annotate_junctions` floors zero-probability alleles to make them reachable.
