<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/antigenomics/vdjtools/master/assets/vdjtools_dark.png">
    <img alt="vdjtools" src="https://raw.githubusercontent.com/antigenomics/vdjtools/master/assets/vdjtools_light.png" width="440">
  </picture>
</p>

<p align="center">
  <a href="https://pypi.org/project/vdjtools/"><img alt="PyPI" src="https://img.shields.io/pypi/v/vdjtools"></a>
  <a href="https://github.com/antigenomics/vdjtools/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/antigenomics/vdjtools/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://antigenomics.github.io/vdjtools/"><img alt="docs" src="https://github.com/antigenomics/vdjtools/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-GPLv3-green"></a>
</p>

TCR/BCR immune-repertoire analysis — a clean-room **Python + C++** rewrite of the legacy
Groovy/Java [vdjtools](https://doi.org/10.1371/journal.pcbi.1004503), standardised on the
**AIRR schema** and **polars** DataFrames with minimal object-orientation.

Built on the antigenomics ecosystem:
[seqtree](https://github.com/antigenomics/seqtree) (fuzzy search / e-value engine),
[vdjmatch](https://github.com/antigenomics/vdjmatch) (overlap + TCRnet),
[arda](https://github.com/antigenomics/arda) (AIRR annotation + markup repair).

> **Status: v2.0.0 pre-release under active development** (latest: `v2.0.0-alpha.2` — the native
> V(D)J model engine). Install a pre-release with `pip install --pre vdjtools`. The legacy v1.x tool
> lives on the [`legacy-1.x`](https://github.com/antigenomics/vdjtools/tree/legacy-1.x) branch and
> its releases remain available under the repository tags (`v0.0.1` … `1.2.1`).

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

## Quickstart — recombination model engine

Precomputed models for all **7 human loci** ship in the wheel — no OLGA or download needed:

```python
from vdjtools.model import load_bundled, native
from vdjtools.model.generate import generate

model = load_bundled("TRB", source="olga")     # or source="learned" (fit to real repertoires)

native.pgen_nt(model, "TGTGCCAGCAGC...")        # nucleotide generation probability (native C++)
native.pgen_aa(model, "CASSLAPGATNEKLFF")       # amino-acid Pgen (codon-marginalised)
native.pgen_aa(model, "CASSLAPGATNEKLFF", mismatches=1)   # + the whole Hamming-1 ball
generate(model, 1000)                           # sample a repertoire -> polars DataFrame
```

Matches OLGA's Pgen to machine precision across all 7 loci, and adds tandem-D (D-D) support that
OLGA/IGoR lack. Learn a model from your own out-of-frame reads with `model.infer.infer_native`.

Explore any model's recombination **Bayes net** interactively (entropy, mutual information, marginals):

```bash
pip install "vdjtools[examples]"
marimo edit notebooks/model_explorer.py
```

## Capabilities (rolling out by phase — see [ROADMAP.md](ROADMAP.md))

- **Model** — native V(D)J recombination model: generation probability (Pgen — nt, aa,
  1-mismatch, V/J-agnostic), sequence generation, and EM inference, all in a native (pybind11)
  core. Supersedes OLGA and IGoR: arda-driven scenario enumeration, polars marginal tables,
  read-parallelised EM, and **tandem-D (D-D)** support. Concordant with OLGA across all 7 loci;
  precomputed OLGA + real-data-learned models bundled ([`load_bundled`](python/vdjtools/model/bundled.py)).
- **Stats** — diversity (Chao1/Shannon/Simpson/…), spectratype, V/J/VJ usage.
- **Features** — CDR physicochemical profiles, k-mer / V+k-mer summaries.
- **Overlap** — sample overlap and TCRnet (via vdjmatch/seqtree).
- **Preprocess** — downsampling, error-correction, VJ-usage batch-effect correction, pooling/joining.
- **Biomarker** — incidence-based association (Fisher) vs HLA / condition / chain-pairing; metaclonotypes.
- **Single-cell** — AIRR Cell / 10x interoperability.

## License

GPL-3.0-or-later.
