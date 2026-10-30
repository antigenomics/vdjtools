<p align="center">
  <!-- GitHub honours <source> and gets the SVG in both colour schemes; PyPI strips <picture>/<source>
       and cannot render SVG, so it falls through to the absolute-URL PNG in the <img>. -->
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/antigenomics/vdjtools/master/assets/vdjtools_dark.svg">
    <source srcset="https://raw.githubusercontent.com/antigenomics/vdjtools/master/assets/vdjtools_light.svg">
    <img alt="vdjtools" src="https://raw.githubusercontent.com/antigenomics/vdjtools/master/assets/vdjtools_light.png" width="340">
  </picture>
</p>

<h1 align="center">vdjtools — immune-repertoire analysis</h1>

<p align="center">
  <a href="https://pypi.org/project/vdjtools/"><img alt="PyPI" src="https://img.shields.io/pypi/v/vdjtools"></a>
  <a href="https://github.com/antigenomics/vdjtools/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/antigenomics/vdjtools/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://docs.isalgo.dev/vdjtools/"><img alt="docs" src="https://github.com/antigenomics/vdjtools/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <img alt="license" src="https://img.shields.io/badge/license-Academic%20%26%20Non--Commercial-green">
</p>

<p align="center"><b>
  <a href="https://docs.isalgo.dev/vdjtools/getting-started.html">Get started</a> ·
  <a href="https://docs.isalgo.dev/vdjtools/usage.html">Guides</a> ·
  <a href="https://docs.isalgo.dev/vdjtools/cli.html">Commands</a> ·
  <a href="https://docs.isalgo.dev/vdjtools/api.html">API</a> ·
  <a href="https://docs.isalgo.dev/vdjtools/notebooks.html">Examples</a> ·
  <a href="https://docs.isalgo.dev/vdjtools/glossary.html">Glossary</a>
</b></p>

Read receptor tables, measure diversity and overlap, fit recombination models,
and export statistical feature tables or published corpus signatures. vdjtools is
a Python/C++ rewrite of the original Groovy/Java tool, using AIRR column names and
Polars frames. arda, seqtree and vdjmatch supply the shared native/reference engines.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/antigenomics/vdjtools/master/assets/analysis-atlas-dark.svg">
  <source srcset="https://raw.githubusercontent.com/antigenomics/vdjtools/master/assets/analysis-atlas-light.svg">
  <img alt="vdjtools analysis atlas: prepare inputs, summarize repertoires, measure diversity, describe sequences, compare samples, model recombination, resolve paired cells and build cohort features." src="https://raw.githubusercontent.com/antigenomics/vdjtools/master/assets/analysis-atlas.png" width="100%">
</picture>

## Choose an analysis

| Analysis type | What you can do | Explore |
|---|---|---|
| **Input & preparation** | Read AIRR, MiXCR, MiGEC, MiTCR, immunoSEQ/Adaptive, IMGT, Vidjil, RTCR, TRUST4, arda, legacy vdjtools and Parquet; filter, correct sequencing errors, downsample, pool and recompute frequencies. | [Formats & preprocessing](https://docs.isalgo.dev/vdjtools/preprocessing.html) |
| **Summary & basic statistics** | Measure V/D/J/C and VJ usage, nucleotide/amino-acid junction-length spectratypes and isotype composition. Weight by reads or unique clonotypes. | [Repertoire statistics](https://docs.isalgo.dev/vdjtools/usage.html) |
| **Diversity & sampling** | Summarize clone-size distributions and estimate richness, Shannon and Simpson diversity; calculate Hill-number rarefaction and extrapolation, sample coverage and bootstrap intervals. | [Diversity guide](https://docs.isalgo.dev/vdjtools/usage.html) |
| **Sequence composition & motifs** | Profile k-mers, V-anchored motifs and amino-acid physicochemical properties, including hydropathy; summarize somatic hypermutation when identity annotations are available. | [Features](https://docs.isalgo.dev/vdjtools/features.html) |
| **Comparison & analysis models** | Compare exact/fuzzy repertoire overlap and sample distances, test TCRnet neighbourhood enrichment against control repertoires, analyse biomarkers and dynamics, and reduce cross-batch V/J-usage bias. | [Analysis guide](https://docs.isalgo.dev/vdjtools/usage.html) |
| **Recombination & Pgen** | Score nucleotide and amino-acid generation probabilities, generate synthetic receptors and fit recombination models with latent gene, trimming and insertion histories. Use current `isalgo` models or explicit frozen `olga-legacy` models. | [Model workshop](https://docs.isalgo.dev/vdjtools/model.html) |
| **Paired data & single cells** | Import 10x/CellRanger, AIRR cell tables and arda cells; inspect chain multiplicity, resolve and pair chains, flag mispairing, calculate paired Pgen and bridge single-cell ecosystems. | [Single-cell guide](https://docs.isalgo.dev/vdjtools/singlecell.html) |
| **Cohort features & signatures** | Export statistical feature tables, fit and reuse sample PCA, or project onto published vsig corpora; join complementary mirpy rsig coordinates by sample. | [Features & PCA](https://docs.isalgo.dev/vdjtools/features.html) · [Signatures](https://docs.isalgo.dev/vdjtools/signature.html) |
| **Downstream tools** | Query VDJdb for receptor specificity with VDJmatch; build receptor and repertoire embeddings with mirpy. | [VDJmatch](https://github.com/antigenomics/vdjmatch) · [mirpy](https://github.com/antigenomics/mirpy) |

**Choose the unit first:** bulk analyses operate on clonotypes; single-cell readers
preserve cell identity and chain pairing. AIRR junctions include the conserved
anchors; CDR3s exclude them. Pgen describes recombination probability, not antigen
specificity. Model fitting accounts for observation errors separately from generation.

## Install

```bash
pip install vdjtools
```

Prebuilt CPython wheels include the native extension. A source install needs a C++
toolchain. Optional extras include `[overlap]` for sample MDS, `[sc]` for single-cell
bridges and `[examples]` for marimo notebooks. See [getting started](https://docs.isalgo.dev/vdjtools/getting-started.html).

## Quickstart

```python
from vdjtools import io, stats   # import before numerical libraries

sample = io.read("clones.tsv")
diversity = stats.diversity_stats(sample)
v_usage = stats.segment_usage(sample, "v")
```

```bash
vdjtools diversity samples/*.tsv -o diversity.tsv
vdjtools features samples/*.tsv --coverage 0.8 --threads 0 -o features.tsv
vdjtools signature --corpus blood samples/*.tsv -o vsig.tsv
```

`features` exports unrotated measurements and channels without a corpus.
Declared log/CLR/logit/arcsine transforms remain; “raw” means before corpus scaling,
winsorization and PCA. `signature` applies the published vsig corpus transform.
These are different outputs. mirpy's rsig geometry half joins vsig on `sample_id`;
mirpy's sample RFF vectors are a separate representation.

For large new sequence draws, select the native sampler and its allocated CPU budget:

```bash
vdjtools generate -m TRB -n 10000 --engine native --threads 0 --seed 7 -o generated.tsv
```

Native draws are identical across thread counts. The default reference sampler preserves
historical seeded draws; the two engines use different random streams.

Fit your own sample PCA on selected measurements, compatible control rows, or reuse
a saved rotation:

```bash
vdjtools pca features.tsv --columns feature_columns.txt --fit rotation.npz --components 32 -o pcs.tsv
vdjtools pca new_features.tsv --rotation rotation.npz -o new_pcs.tsv
```

The [raw feature/PCA guide](https://docs.isalgo.dev/vdjtools/features.html) explains
feature selection, coverage, imputation, feature-space identities and saved rotations.

## Find the right page

| Task | Guide |
|---|---|
| Load, measure and compare your first samples | [Getting started](https://docs.isalgo.dev/vdjtools/getting-started.html) |
| Filter, correct, downsample or pool | [Pre-processing](https://docs.isalgo.dev/vdjtools/preprocessing.html) |
| Diversity, overlap, gene usage, biomarkers and dynamics | [User guide](https://docs.isalgo.dev/vdjtools/usage.html) |
| Unrotated features and optional sample PCA | [Raw features](https://docs.isalgo.dev/vdjtools/features.html) |
| Published vsig/rsig corpus coordinates | [Signatures](https://docs.isalgo.dev/vdjtools/signature.html) · [Methods](https://docs.isalgo.dev/vdjtools/signature-methods.html) · [Channels](https://docs.isalgo.dev/vdjtools/channels.html) |
| Pgen, generation and model fitting | [Model workshop](https://docs.isalgo.dev/vdjtools/model.html) |
| Single-cell data and interoperability | [Single cell](https://docs.isalgo.dev/vdjtools/singlecell.html) |
| Functions, flags, examples and terms | [API](https://docs.isalgo.dev/vdjtools/api.html) · [Commands](https://docs.isalgo.dev/vdjtools/cli.html) · [Examples](https://docs.isalgo.dev/vdjtools/notebooks.html) · [Glossary](https://docs.isalgo.dev/vdjtools/glossary.html) |

Recorded [performance measurements](https://docs.isalgo.dev/vdjtools/index.html#performance)
and [model validation](https://docs.isalgo.dev/vdjtools/model.html) describe the
native Pgen, generation and inference engines.

## Resource budgets

Import `vdjtools` before numerical libraries in a fresh interpreter. Native batch
APIs use allocation-aware CPU budgets; sample workers use serial inner kernels.
Signature `--jobs` and analysis `--threads` select bounded workers, with `0` for the
detected allocation. Excessive or nested budgets raise. The final feature table
still grows with sample count and feature width. Configured randomized corpus fits
and `TablePCA.fit_features` use sparse rows and bounded statistic chunks; their
remaining workspace grows with nonzero entries and the requested component count. See
[CPU/container limits](https://docs.isalgo.dev/vdjtools/usage.html#cpu-container-limits).
For file-backed Python-heavy summaries, select `--backend process` (or
`io.map_samples(..., backend="process")`) to use coarse sample batches. The default
thread backend shares resources and supports closures; process startup can dominate
small cohorts.

## Development

```bash
bash setup.sh --dev-parents --tests
source .venv/bin/activate
uv pip install -e ".[docs]"
pytest tests/python -q -m "not slow and not heavy and not benchmark"
env -C docs make html
```

Uses a repo-local uv environment and a C++ toolchain. See the API reference for
public functions and the model workshop for native model workflows.

## Citing

For the v2 rewrite, cite this repository. For the original tool, cite
[Shugay et al., *PLoS Computational Biology* 2015](https://doi.org/10.1371/journal.pcbi.1004503).
The VDJtrack recapture model in `vdjtools.dynamics` is
[Pavlova, Zvyagin and Shugay (2024)](https://doi.org/10.3389/fimmu.2024.1321603).

## License

ISALGO LAB Academic and Non-Commercial Research License. The legacy Groovy/Java v1.x tool lives on the
[`legacy-1.x`](https://github.com/antigenomics/vdjtools/tree/legacy-1.x) branch, with its releases
under the repository tags `v0.0.1` … `1.2.1`.

### Configurable statistical features

Opt-in `basic` uses diversity, V–J usage, junction-length spectratype and separate
isotype composition; `joint` uses diversity, V–J–length usage and separate isotype.
Their extended variants add IMGT11 or frozen BLOSUM62 motifs (`XXX`, `XX_X`, `X_XX`),
with optional gene or clustered V context. Existing feature/signature defaults stay
legacy. See [features](https://docs.isalgo.dev/vdjtools/features.html) and
[corpus fitting](https://docs.isalgo.dev/vdjtools/signature.html) for fit/save/apply
and the distinction between an available definition and a published rotation.

Public recombination models, corpus bundles and example sources are mirrored on
[Hugging Face](https://huggingface.co/datasets/isalgo/airr_feature_signatures).
The dataset card shows selective downloads; downloading the whole mirror is unnecessary.

Working model APIs and CLIs default to `isalgo`, targeting seven human loci and
mouse TRA/TRB/IGH/IGK/IGL. All twelve fitted models ship in the package; see the
[validation table](https://docs.isalgo.dev/vdjtools/model.html#release-5-4-validation)
for exact Pgen support and fixed-model observation likelihood. `olga-legacy` selects frozen OLGA models
(`olga` is its alias). `learned` and `arda` alias fresh `isalgo`; their previous
parameter families require a pinned previous release.

For a sample too large to hold in memory, `vdjtools.io.iter_table_batches` yields
ordered raw CSV/TSV/gzip/Parquet batches with optional column projection. It does
not collapse clonotypes or calculate chunk-local frequencies; consumers retain
whole-sample counting semantics. Prefer `iter_samples` when one sample fits.

For joint junction/V/J probabilities with gene-level calls, use
`vdjtools.model.native.pgen_aa_gene_batch`: one native batch sums the model alleles
and supports exact or closed Hamming-1 mass. ALICE uses the same helper.
