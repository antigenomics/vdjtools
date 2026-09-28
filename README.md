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

<p align="center"><b>
  <a href="https://docs.isalgo.dev/vdjtools/getting-started.html">Get started</a> ·
  <a href="https://docs.isalgo.dev/vdjtools/usage.html">Guides</a> ·
  <a href="https://docs.isalgo.dev/vdjtools/cli.html">Commands</a> ·
  <a href="https://docs.isalgo.dev/vdjtools/api.html">API</a> ·
  <a href="https://docs.isalgo.dev/vdjtools/glossary.html">Glossary</a>
</b></p>

Analysis of T- and B-cell receptor repertoire sequencing: read any format, measure diversity and
overlap, correct batch effects, score generation probability against a native V(D)J model, and
reduce a whole repertoire to one comparable feature vector.

vdjtools reads the output of every common repertoire pipeline into one canonical
[AIRR](https://docs.airr-community.org/) table backed by [polars](https://pola.rs), and every
analysis takes and returns that same table — so results chain together instead of needing a
converter between each step.

A clean-room Python + C++ rewrite of the Groovy/Java
[vdjtools](https://doi.org/10.1371/journal.pcbi.1004503), built on
[arda](https://github.com/antigenomics/arda) (germline reference and annotation),
[seqtree](https://github.com/antigenomics/seqtree) (fuzzy search and e-values) and
[vdjmatch](https://github.com/antigenomics/vdjmatch) (overlap and TCRnet).

## Install

```bash
pip install vdjtools
```

That is the whole installation — every capability below works from it, with no extras to opt into.
Wheels are prebuilt for CPython 3.10–3.13 on Linux, macOS (Apple Silicon) and Windows and include
the native C++ extension; a source install compiles it. The three engines above are base
dependencies, imported lazily so `import vdjtools` stays light.

Two optional extras exist because each has a working alternative:

```bash
pip install "vdjtools[overlap]"   # scikit-learn, for cluster_samples(method="mds")
pip install "vdjtools[sc]"        # single-cell bridges: anndata, awkward, mudata, pyyaml
```

## Quickstart

```python
from vdjtools import io as vio, stats, overlap

sample = vio.read("clones.tsv")          # MiXcr, immunoSEQ, AIRR, Parquet, native — detected
stats.diversity_stats(sample)            # richness, Chao1, Shannon, inverse Simpson, d50
stats.segment_usage(sample, "v")         # V-gene usage

cohort = vio.read_samples(vio.read_metadata("metadata.txt"), base_dir="samples/")
overlap.overlap_matrix(cohort)           # pairwise repertoire overlap
```

The same thing with no Python:

```bash
vdjtools diversity     sampleA.tsv sampleB.tsv -o diversity.tsv
vdjtools segment-usage *.tsv --segment v -o usage.tsv
vdjtools overlap       *.tsv -o overlap.tsv
```

Every command writes to `-o` — TSV, or Parquet when the path ends in `.parquet` / `.pq` — or to
stdout, so commands pipe. Analysis commands also take a metadata table (`-m` plus `--base-dir`) and
parallelise over samples with `-t/--threads`.

**[Getting started](https://docs.isalgo.dev/vdjtools/getting-started.html)** walks through this in
about ten minutes, with no data to download.

## Where to start

| You want to | Go to |
|---|---|
| Get a first result from a file you already have | [Getting started](https://docs.isalgo.dev/vdjtools/getting-started.html) |
| Clean, filter, downsample or batch-correct a cohort | [Pre-processing](https://docs.isalgo.dev/vdjtools/preprocessing.html) |
| Measure diversity, overlap, gene usage, spectratype | [Guides](https://docs.isalgo.dev/vdjtools/usage.html) |
| Score Pgen, or fit a model to your own reads | [Model workshop](https://docs.isalgo.dev/vdjtools/model.html) |
| Reduce each sample to one feature vector for a classifier | [Signatures](https://docs.isalgo.dev/vdjtools/signature.html) |
| Work with 10x data, or hand off to scirpy / dandelion | [Single cell](https://docs.isalgo.dev/vdjtools/singlecell.html) |
| Look up a command or a function | [Commands](https://docs.isalgo.dev/vdjtools/cli.html) · [API](https://docs.isalgo.dev/vdjtools/api.html) |
| Check what a term means | [Glossary](https://docs.isalgo.dev/vdjtools/glossary.html) |

## What is in the box

| Module | What it does |
|---|---|
| `vdjtools.io` | Readers for native vdjtools, AIRR Rearrangement TSV and Parquet; format-detecting converters for MiXcr, MiGec, Adaptive immunoSEQ, IMGT/HighV-QUEST, Vidjil, RTCR, TRUST4 and arda; metadata-driven batches and hive-partitioned cohorts |
| `vdjtools.stats` | Diversity (richness, Chao1, Efron-Thisted, Shannon, Simpson, d50), coverage-based Hill-number rarefaction with bootstrap intervals, spectratype, V/D/J/C usage |
| `vdjtools.model` | Native V(D)J recombination model: Pgen over nucleotides and amino acids, the Hamming-1 ball, sequence generation, EM inference, tandem-D support, and a 13-command model workshop |
| `vdjtools.overlap` | Exact and similarity-aware repertoire overlap, TCRnet and ALICE neighbourhood enrichment, sample clustering |
| `vdjtools.preprocess` | Format conversion, the three filtering axes, frequency handling, downsampling, error correction, V/J-usage batch correction, pooling and joining |
| `vdjtools.features` | CDR3 physicochemical profiles and k-mer summaries |
| `vdjtools.biomarker` | Incidence association against binary, per-HLA-allele or stratified conditions; co-occurrence pairing; metaclonotypes |
| `vdjtools.dynamics` | Longitudinal clonotype tracking: paired within-donor expansion testing, the VDJtrack recapture model, an edgeR NB-exact caller |
| `vdjtools.signature` | One fixed, named, already-standardised feature vector per repertoire, rotated through a published corpus |
| `vdjtools.sc` | Single-cell ingestion, chain pairing and QC, paired Pgen, round-trip interop with scirpy, dandelion and scRepertoire |

## Recombination models, built in

Models for all **seven human loci** ship in the wheel — no OLGA install and nothing to download:

```python
from vdjtools.model import load_bundled, native

model = load_bundled("TRB", source="olga")            # or source="learned"
native.pgen_aa(model, "CASSLAPGATNEKLFF")             # 3.595e-08
native.pgen_aa(model, "CASSLAPGATNEKLFF", mismatches=1)   # and its whole Hamming-1 ball
```

Pgen agrees with OLGA to 1e-15 on every locus while being several times faster, and adds tandem-D
support that OLGA and IGoR lack. You can also fit a model to your own non-productive reads, audit it
against its germline, compare two models parameter by parameter, and render its recombination Bayes
net — see the [model workshop](https://docs.isalgo.dev/vdjtools/model.html).

## Signatures

One repertoire in, one fixed-width named feature vector out, standardised against a published
reference corpus so your matrix and a collaborator's are in the same coordinate system without
either of you fitting a scaler:

```bash
vdjtools signature --corpus blood samples/*.tsv.gz -o vsig.tsv
mir       signature --corpus blood samples/*.tsv.gz -o rsig.tsv   # the geometry half, from mirpy
```

Nine corpora are published, all at 256 components per locus. `blood` (11,117 real samples, 947 study
groups), `tissue` (21,131 / 1,934) and `deep-tcr` (3,936 amplicon samples, TRA+TRB) are fitted on
real repertoires; the synthetic ones need no cohort at all and rebuild bit-identically anywhere.
Artifacts are fetched and digest-verified on first use rather than bundled.

`--corpus` is required: a signature is comparable to another one only if both were rotated through
the same corpus, and nothing about the numbers would say otherwise.
[Signatures](https://docs.isalgo.dev/vdjtools/signature.html) ·
[How they were fitted](https://docs.isalgo.dev/vdjtools/signature-methods.html) ·
[Channels](https://docs.isalgo.dev/vdjtools/channels.html)

## Performance

The Pgen, generation, EM and diversity hot paths are a native C++ core reached through pybind11;
everything else is polars. Single thread, Apple M3, bundled human TRB model:

| Operation | Throughput | vs OLGA |
|---|--:|--:|
| Nucleotide Pgen, single-D VDJ | 0.5 ms/seq | **9×** |
| Amino-acid Pgen | 0.6–0.9 ms/seq | **8.6×** |
| Amino-acid Pgen over the Hamming-1 ball | 15 ms/seq | **8.7×** |
| Sequence generation, raw draws | 32,100 seq/s | — |
| Sequence generation, `productive_only=True` | 8,960 seq/s | — |

The productive filter costs 3.6× on TRB and 5.1× on IGH (19,900 raw draws/s against 3,930
productive), because out-of-frame and stop-codon draws are discarded and redrawn. Batched Pgen
parallelises over sequences (11× on 16 cores, bit-identical to serial); the EM E-step parallelises
over reads (6.7× on 8 threads). Memory: 63 MB resident for `import vdjtools` plus one model, 123 MB
with all seven loaded.

Reproduce the in-repo half with `RUN_BENCHMARK=1 pytest tests/python -k benchmark`.

## Development

Uses [uv](https://docs.astral.sh/uv/) — one repo-local `.venv`, no conda:

```bash
uv venv && source .venv/bin/activate
uv pip install -e ".[dev,test]"       # builds the _core C++ extension
pytest tests/python -q
```

`bash setup.sh --dev-parents --tests` does the same and editable-installs sibling checkouts of
seqtree, arda and vdjmatch if present. You need a C++ toolchain (Xcode Command Line Tools on macOS,
`build-essential` on Linux). MMseqs2 is needed only for arda's annotation path.

## Citing

For the v2 rewrite, cite this repository. For the original tool, cite
[Shugay et al., *PLoS Computational Biology* 2015](https://doi.org/10.1371/journal.pcbi.1004503).
The VDJtrack recapture model in `vdjtools.dynamics` is Pavlova, Zvyagin and Shugay (2024).

## License

GPL-3.0-or-later. The legacy Groovy/Java v1.x tool lives on the
[`legacy-1.x`](https://github.com/antigenomics/vdjtools/tree/legacy-1.x) branch, with its releases
under the repository tags `v0.0.1` … `1.2.1`.
