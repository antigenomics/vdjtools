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

> **Status: `v2.2.0`** — the native V(D)J model engine plus the full analytics suite (diversity,
> overlap/TCRnet, preprocessing, biomarkers, single-cell), CDR features, and legacy-format ingestion
> (MiXcr, MiGec, immunoSEQ, IMGT/HighV-QUEST, Vidjil, RTCR, TRUST4, arda). Clonotype columns follow the AIRR
> **junction** convention (`junction_nt` / `junction_aa`). The legacy v1.x tool lives on the
> [`legacy-1.x`](https://github.com/antigenomics/vdjtools/tree/legacy-1.x) branch and its releases
> remain available under the repository tags (`v0.0.1` … `1.2.1`).

## Install

```bash
pip install vdjtools
```

Prebuilt wheels ship for CPython 3.10–3.13 on Linux, macOS (Apple Silicon), and Windows; the
native `_core` C++ extension is bundled (the source distribution compiles it on install).

That one command gives you **everything vdjtools advertises** — no extras to opt into. The three
antigenomics engines it delegates to are **base dependencies**:

| Engine | Powers |
|---|---|
| [arda](https://github.com/antigenomics/arda) | the germline reference (V/D/J + CDR3 anchors), model engine, annotation |
| [seqtree](https://github.com/antigenomics/seqtree) | fuzzy search / e-values → error correction, similarity overlap, TCRnet |
| [vdjmatch](https://github.com/antigenomics/vdjmatch) | sample overlap, TCRnet, metaclonotypes |

So `preprocess.correct()`, `overlap.tcrnet()`, `biomarker.metaclonotypes()` and
`model.reference.load_germline()` all just work from a plain install — and downstream libraries
can depend on plain `vdjtools` and rely on them being there. All three are imported **lazily**, so
`import vdjtools` stays light; between them they add only `requests` on top of what vdjtools
already needs. arda fetches its IMGT reference once, on first use.

Two things are still optional, because each has a working alternative or is a side integration:

```bash
pip install "vdjtools[overlap]"   # scikit-learn — only for cluster_samples(method="mds");
                                  # method="hclust" works out of the box (scipy)
pip install "vdjtools[sc]"        # single-cell extras: anndata (scverse bridge), pyyaml
```

MMseqs2 is needed **only** for arda's alignment/annotation path (`model.stitch.annotate`) — never
for germline lookup, Pgen, generation, or the analytics. Install it via conda/brew if you need it.

### Development

```bash
conda env create -f environment.yml   # python + C++ toolchain + mmseqs2 (arda's aligner)
conda activate vdjtools
pip install -e ".[dev,test]"          # builds the _core C++ extension
```

The conda env is a convenience, not a requirement — `pip install -e ".[dev,test]"` in any venv
works. It supplies MMseqs2 so the slow arda annotation round-trips in the test suite run too.

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
native.pgen_aa_batch(model, seqs, mismatches=1, threads=0)  # Pgen over many CDR3s, thread-parallel (~11×)
generate(model, 1000)                           # sample a repertoire -> polars DataFrame
```

Matches OLGA's Pgen to machine precision across all 7 loci, and adds tandem-D (D-D) support that
OLGA/IGoR lack. Learn a model from your own **non-functional** reads (out-of-frame *or* stop-codon — both escaped
selection, which is all a generative model needs) with `model.infer.infer_native`.

Explore any model's recombination **Bayes net** interactively (entropy, mutual information, marginals):

```bash
pip install "vdjtools[examples]"
marimo edit notebooks/model_explorer.py
```

## Command line

`pip install vdjtools` installs the `vdjtools` command — the model engine (OLGA/IGoR-style) and the
repertoire analytics (over sample files or a metadata table, like the legacy tool):

```bash
# recombination model engine — built-in models for all 7 loci (no download)
vdjtools models                                # list the bundled models
vdjtools generate -m TRB -n 1000 -o gen.tsv    # sample sequences   (cf. olga-generate_sequences)
vdjtools pgen seqs.tsv -m TRB -o pgen.tsv      # Pgen per CDR3       (cf. olga-compute_pgen)
vdjtools pgen seqs.tsv -m TRB --mismatches 1   # + the Hamming-1 ball; --v-col/--j-col to condition

# repertoire analytics — sample files, or a cohort via -m/--metadata + --base-dir
vdjtools diversity      sampleA.tsv sampleB.tsv -o diversity.tsv
vdjtools overlap        *.tsv -o overlap.tsv
vdjtools segment-usage  *.tsv --segment v -o usage.tsv
vdjtools spectratype    *.tsv -o spectra.tsv
```

Native vdjtools and AIRR Rearrangement inputs are auto-detected; every command writes TSV to `-o`
(or stdout, so it pipes). Run `vdjtools <command> --help` for options.

## Analytics (Python API)

Every reader returns one canonical `polars` clonotype frame (AIRR **junction** columns), and every
analysis function takes and returns such frames — so results chain together and drop straight into
plotting. A tour of the analysis modules (full runnable walkthrough in the
[**User guide**](https://docs.isalgo.dev/vdjtools/usage.html)):

```python
from vdjtools import io as vio, stats, features, overlap, preprocess

# load (auto-detects MiXcr / immunoSEQ / AIRR / native / … and converts), or a whole cohort:
sample = vio.read("clones.tsv")
cohort = vio.read_samples(vio.read_metadata("metadata.txt"), base_dir="samples/")

# diversity, rarefaction, segment usage, spectratype
stats.diversity_stats(sample)                 # observed, Chao1, Shannon, inverse-Simpson, d50, …
stats.inext(sample, q=(0, 1, 2))              # Hill-number rarefaction/extrapolation + bootstrap CIs
stats.segment_usage(sample, "v")              # V (or "j") usage;  stats.spectratype(sample)

# CDR3 physicochemistry & k-mers
features.physchem_profile(sample, region="all")

# repertoire overlap & TCRnet (fuzzy/similarity/TCRnet via the [overlap] engine)
overlap.overlap_metrics(sampleA, sampleB)     # F / D / Jaccard / Morisita-Horn …
overlap.tcrnet(sample)                         # per-clonotype neighbourhood enrichment

# preprocessing: downsample to a common depth, error-correct, filter, pool
preprocess.downsample(sample, 100_000)
preprocess.correct(preprocess.filter_functional(sample))

# cross-batch V/J-usage bias: batch-correct usage, then resample the clonotype table
usage = preprocess.correct_vj_usage(cohort, batch_col="batch", transform="sigmoid")  # Vlasova 2026
fixed = preprocess.apply_vj_correction(sampleA, usage, sample_id="A0")
```

Incidence-based clonotype association (Emerson 2017 / Howie 2015 / De Witt 2018 / Vlasova 2026)
— a choice of test, condition, and co-occurrence — and single-cell paired-chain Pgen:

```python
from vdjtools import biomarker, sc
from vdjtools.biomarker import association, condition

# feature vs condition: Fisher / chi2 / Bayesian / permutation; binary, per-HLA-allele, or CMH-stratified
association(cohort, condition.binary(meta, "cmv"), test=["fisher", "bayes_bf"])
association(cohort, condition.stratified(meta, "cmv", "hla"), stratum_col="_stratum")  # CMV | HLA (CMH)

# feature vs feature: in-silico α-β pairing / same-chain co-specificity (θ lift + Fisher + FDR)
biomarker.cooccurrence(cohort, chain_a="TRA", chain_b="TRB", evalue=True)

sc.paired_pgen(sc.pair_chains(sc.read_10x("filtered_contig_annotations.csv")))  # pgen_alpha·pgen_beta
```

## Performance

The Pgen / generation / EM / diversity hot paths are a native C++ (pybind11) core; everything else is
polars. Amino-acid Pgen matches OLGA to machine precision (1e-15) across all 7 loci while being
several times faster, and the built-in models keep the resident set small. Single thread, Apple M3
(arm64), bundled human TRB model:

| operation | throughput | vs OLGA |
|---|---|---|
| nucleotide Pgen (single-D VDJ) | **~0.5 ms/seq** | **9×** |
| amino-acid Pgen | **~0.6–0.9 ms/seq** | **8.6×** |
| Pgen + Hamming-1 ball (1 substitution) | **~15 ms/seq** | **8.7×** |
| sequence generation | **~32 000 seq/s** | — |

Nucleotide Pgen (via the same transfer-matrix DP as the aa path — an in-frame CDR3 is an aa query with
one codon fixed per position) is exact vs OLGA across all loci. Batched Pgen / 1-mismatch over many
CDR3s parallelises over sequences (`native.pgen_aa_batch`, **~11× on 16 cores**, bitwise-identical to
the serial result); the EM E-step parallelises over reads (~6.7× on 8 threads); diversity/rarefaction
run on a native iNEXT kernel (bootstrap + parallel batch). Memory
stays light — **~63 MB** resident for `import vdjtools` plus one loaded model, **~123 MB** with all
seven bundled models resident. Reproduce with `~/vcs/projects/2026-vdjtools-benchmark/bench/bench_pgen.py` and the `test_*_benchmark.py`
suites (`RUN_BENCHMARK=1`).

## Capabilities (see the [User guide](https://docs.isalgo.dev/vdjtools/usage.html), the [API reference](https://docs.isalgo.dev/vdjtools/), and [ROADMAP.md](ROADMAP.md))

- **IO** — canonical clonotype frame on AIRR **junction** columns (`junction_nt` / `junction_aa`);
  readers for native vdjtools, AIRR Rearrangement TSV, and Parquet, plus format-detecting converters
  for MiXcr (v1/2 + v3/4, incl. C-gene / BCR isotype), MiGec, Adaptive immunoSEQ (v1/v2),
  IMGT/HighV-QUEST, Vidjil, RTCR, TRUST4, and arda AIRR output
  ([`vdjtools.io.convert`](python/vdjtools/io/convert.py)); metadata-driven batch + hive-partitioned cohorts.
- **Model** — native V(D)J recombination model: generation probability (Pgen — nt, aa,
  1-mismatch, V/J-agnostic, **thread-parallel batch**), sequence generation, and EM inference, all in a
  native (pybind11) core. Supersedes OLGA and IGoR: arda-driven scenario enumeration, polars marginal
  tables, read-parallelised EM, and **tandem-D (D-D)** support. Concordant with OLGA across all 7 loci;
  precomputed OLGA + real-data-learned models bundled ([`load_bundled`](python/vdjtools/model/bundled.py)).
- **Stats** — diversity (Chao1/Shannon/Simpson/…), spectratype, V/J/VJ usage.
- **Features** — CDR physicochemical profiles, k-mer / V+k-mer summaries.
- **Overlap** — sample overlap and TCRnet (via vdjmatch/seqtree), similarity-aware overlap, clustering.
- **Preprocess** — downsampling, error-correction, VJ-usage batch-effect correction, pooling/joining.
- **Biomarker** — incidence association (Fisher / χ² / Bayesian / permutation) vs binary / HLA-allele / CMH-stratified conditions; α-β & same-chain co-occurrence pairing; metaclonotypes.
- **Single-cell** — AIRR Cell / 10x interoperability, chain pairing + QC, and paired α/β Pgen.

## License

GPL-3.0-or-later.
