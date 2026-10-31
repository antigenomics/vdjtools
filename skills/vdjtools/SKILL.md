---
name: vdjtools
description: Use when loading AIRR repertoires, computing diversity or overlap, extracting cohort features or signatures, scoring receptor generation probabilities, or handling paired single-cell receptors with vdjtools.
---

# vdjtools

Use the installed public APIs and batched operations. In a source checkout, read
its repository instructions before changing code. Documentation links below work
outside the checkout; corresponding sources are under `docs/` in this repository.
Check the installed version and CLI help before using version-dependent options.

## Choose the workflow

| Request | API entry point | Guide |
|---|---|---|
| Read, convert, filter or pool bulk repertoires | `vdjtools.io`, `vdjtools.preprocess` | [Preprocessing](https://docs.isalgo.dev/vdjtools/preprocessing.html) |
| Diversity, gene usage, overlap or TCRnet | `vdjtools.stats`, `vdjtools.overlap` | [Analysis](https://docs.isalgo.dev/vdjtools/usage.html) |
| Unrotated cohort measurements or sample PCA | `raw_features_cohort`, `TablePCA` | [Features](https://docs.isalgo.dev/vdjtools/features.html) |
| Project onto a fitted vsig corpus | `vsig`, `vsig_cohort`, `Corpus` | [Signatures](https://docs.isalgo.dev/vdjtools/signature.html) |
| Pgen, generation, model fitting or selection mass | `vdjtools.model` | [Models](https://docs.isalgo.dev/vdjtools/model.html) |
| Paired receptors, 10x, AIRR cells or interop | `vdjtools.sc` | [Single cells](https://docs.isalgo.dev/vdjtools/singlecell.html) |

Use the [API reference](https://docs.isalgo.dev/vdjtools/api.html) for import paths
and signatures, and [CLI reference](https://docs.isalgo.dev/vdjtools/cli.html) for flags.

## Start with the data unit

```python
from vdjtools import io, stats  # import before numerical runtimes

sample = io.read("clones.tsv")
diversity = stats.diversity_stats(sample)
v_usage = stats.segment_usage(sample, "v")
```

Bulk APIs consume AIRR-style Polars frames. Use single-cell readers when cell
identity, multiple chains or pairing matter; bulk collapse loses that information.
AIRR junctions include conserved anchors; CDR3s exclude them. Check which field a
scoring API expects rather than exchanging the two. Observed DNA IUPAC ambiguity
symbols normalize to A at ingress; unrelated invalid symbols raise. Preserve raw
inputs and do not enumerate ambiguous bases as alternative observations.

For large files, `io.iter_table_batches` yields ordered raw batches and optional
column projection. It does not collapse clonotypes or compute whole-sample
frequencies. Preserve whole-sample counting when combining batches.

## Features and signatures

- `vdjtools features` exports unrotated measurements and channels. Declared
  log/CLR/logit/arcsine transforms still apply; raw does not mean untransformed.
- `vdjtools signature` applies a fitted vsig corpus. mirpy rsig is a distinct
  representation; join by sample ID rather than treating coordinates as interchangeable.
- Reuse `vdjtools.signature.corpus.TablePCA` for feature-table PCA. Select measurement
  columns explicitly and preserve IDs, masks and coverage channels.
- Separate cohorts need the same feature definition, vocabulary and coverage target.
  Equal table width does not establish compatibility. Missing values are not zeros.
- Configurable `basic`/`joint` feature definitions are opt-in. A definition's
  availability does not imply that a fitted corpus exists for it. Read the feature
  and signature guides before fitting or applying a configured corpus.

## Models and Pgen

Working model APIs default to `isalgo`: seven human loci and mouse
TRA/TRB/IGH/IGK/IGL. `learned` and `arda` alias this family; `olga-legacy` selects
frozen OLGA models (`olga` is its compatibility alias). Do not substitute families
when a requested model is missing.

Use native batch APIs for many queries. For joint junction/V/J probabilities with
gene-level calls, `model.native.pgen_aa_gene_batch` sums the appropriate model
alleles in one batch. Unknown calls raise; do not catch errors and silently
marginalize over genes. Explicit multi-allele models retain allele-specific scoring.

Pgen is exact and error-free by default. The explicit AA option
`approximate_tandem=True` (CLI `--approximate-tandem`) uses rank-one tandem
coefficients for supported dense joint-history models, including human isalgo IGH.
Frame/file outputs add `pgen_method=rank1`. Unsupported structures, nucleotide
inputs and Hamming-1 requests raise. Read the model guide's accuracy and scope
before selecting it. Model parameters remain unchanged; never cache computed Pgen.

The `generate` CLI defaults to native sampling for bundled isalgo models and
reference sampling for legacy sources/custom paths. The Python
`vdjtools.model.generate.generate` primitive defaults to reference sampling.
For large new draws select `engine="native", threads=0`. Native draws are
CPU-count invariant, but reference and native engines use different seeded streams.

For fitting, ARDA calls are preliminary proposals rather than fixed gene or
boundary assignments. Preserve gene ambiguity, alternative trimming/N/P histories,
retained-germline mismatches and genomic order. Missing preliminary D2 calls do
not exclude tandem histories. Training observation likelihood is distinct from
error-free Pgen. Native inference defaults to `prune_ratio=1e-5`; zero is exact.
The cutoff does not bound total omitted mass; inspect reported omitted mass and
actual positive-likelihood observation counts. Consult the model guide before fitting.

`model.selection_factor` validates model-file identity before returning a bundled
calibration. `selection_mass` applies Q once and is not a binding probability.
Empirical control frequencies already contain selection; do not multiply them
again. Unsupported or modified models require explicit calibration choices.

## Resource use and verification

Import vdjtools before numerical runtimes. Zero threads selects the detected CPU
allocation; positive budgets must fit it. Use one native batch rather than a pool
of per-sequence calls. Sample-map workers and numerical fitting run in separate
phases; consult the owning API before interpreting `threads` or `n_jobs`.

File-backed summaries can use `io.map_samples(..., backend="process")` with
importable callbacks and a guarded script entry point. The default thread backend
supports shared models and closures. Do not silently fall back if workers fail.

For source changes, use the repository environment and focused tests. For docs,
install the docs extra in that environment and run `env -C docs make html` with
zero warnings. Check examples against real signatures and CLI help. Do not rerun
completed analyses or refit published corpora merely to update documentation.
