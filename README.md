# vdjtools 2.0

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

## Install

```bash
conda env create -f environment.yml   # python + mmseqs2 (arda backend) + C++ toolchain
conda activate vdjtools
pip install -e ".[dev,test]"          # builds the _core C++ extension
```

Or run the bootstrap script: `bash setup.sh --dev-parents --tests`.

Pure-analytics use (diversity / spectratype / usage) does not require MMseqs2; the model and
annotation paths do (via arda).

## License

GPL-3.0-or-later.
