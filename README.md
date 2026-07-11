<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/vdjtools_dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="assets/vdjtools_light.svg">
    <!-- Absolute PNG fallback: PyPI strips <picture>/<source> and cannot render a relative or
         raw-served SVG, so the logo must be an absolute-URL raster here. GitHub uses the SVG sources. -->
    <img alt="vdjtools" src="https://raw.githubusercontent.com/antigenomics/vdjtools/master/assets/vdjtools_dark.png" width="320">
  </picture>
</p>

<h1 align="center">vdjtools — immune-repertoire analysis</h1>

<p align="center">
  <a href="https://pypi.org/project/vdjtools/"><img alt="PyPI" src="https://img.shields.io/pypi/v/vdjtools"></a>
  <a href="https://github.com/antigenomics/vdjtools/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/antigenomics/vdjtools/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://docs.isalgo.dev/vdjtools/"><img alt="docs" src="https://github.com/antigenomics/vdjtools/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <img alt="license" src="https://img.shields.io/badge/license-GPLv3-green">
</p>

TCR/BCR immune-repertoire analysis — a clean-room **Python + C++** rewrite of the legacy
Groovy/Java [vdjtools](https://doi.org/10.1371/journal.pcbi.1004503), standardised on the
**AIRR schema** and **polars** DataFrames with minimal object-orientation.

Built on the antigenomics ecosystem:
[seqtree](https://github.com/antigenomics/seqtree) (fuzzy search / e-value engine),
[vdjmatch](https://github.com/antigenomics/vdjmatch) (overlap + TCRnet),
[arda](https://github.com/antigenomics/arda) (AIRR annotation + markup repair).

> **Status: v2.0.0 under active development.** The legacy v1.x tool lives on the
> [`legacy-1.x`](https://github.com/antigenomics/vdjtools/tree/legacy-1.x) branch and its
> releases remain available under the repository tags (`v0.0.1` … `1.2.1`).

## Install

```bash
pip install vdjtools
```

Prebuilt wheels ship for CPython 3.10–3.13 on Linux, macOS (Apple Silicon), and Windows; the
native `_core` C++ extension is bundled (the source distribution compiles it on install). The
pure-analytics paths (diversity / spectratype / usage / overlap) work out of the box; the model
and annotation paths additionally pull in [arda](https://github.com/antigenomics/arda) (MMseqs2):

```bash
pip install "vdjtools[model]"
```

### Development

```bash
conda env create -f environment.yml   # python + mmseqs2 (arda backend) + C++ toolchain
conda activate vdjtools
pip install -e ".[dev,test]"          # builds the _core C++ extension
```

Or run the bootstrap script: `bash setup.sh --dev-parents --tests`.

## Capabilities (rolling out by phase — see [ROADMAP.md](ROADMAP.md))

- **Model** — native V(D)J recombination model: generation probability (Pgen), sequence
  generation, and EM model inference. Supersedes OLGA and IGoR: arda-driven scenario
  enumeration, polars marginal tables, D-D tandem support, and a native (pybind11) Pgen/EM core.
- **Stats** — diversity (Chao1/Shannon/Simpson/…), spectratype, V/J/VJ usage.
- **Features** — CDR physicochemical profiles, k-mer / V+k-mer summaries.
- **Overlap** — sample overlap and TCRnet (via vdjmatch/seqtree).
- **Preprocess** — downsampling, error-correction, VJ-usage batch-effect correction, pooling/joining.
- **Biomarker** — incidence-based association (Fisher) vs HLA / condition / chain-pairing; metaclonotypes.
- **Single-cell** — AIRR Cell / 10x interoperability.

## License

GPL-3.0-or-later.
