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
