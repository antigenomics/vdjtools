# vdjtools examples

Every example is an interactive [marimo](https://marimo.io) notebook — `pip install -e ".[examples]"`
(add `,overlap` / `,sc` where noted), then `marimo edit examples/<name>.py`.

**Data** auto-loads from HuggingFace, but each notebook first looks in a gitignored **`./data_dump/`**
directory at the repo root — drop the datasets there, or symlink your local copies, and nothing is
re-downloaded:

```bash
mkdir -p data_dump
ln -s /path/to/airr_benchmark data_dump/airr_benchmark        # etc. per dataset
# VDJdb: notebooks fetch the latest antigenomics/vdjdb-db release automatically,
# or drop vdjdb.slim.txt(.gz) in data_dump/ to skip the download.
```

## `aging.py` — TCR repertoire aging (streaming + iNEXT + overlap)

A [marimo](https://marimo.io) notebook reproducing the classic aging signals of the human
TCR-beta repertoire on the full-depth Britanova **"Cord Blood to Centenarians"** cohort (78
donors, ages 0–103), reading each signal three complementary ways from a single cohort load:

- **Cohort-streaming stats** — `stats.diversity_cohort`, the singleton→hyperexpanded clone-size
  distribution, and the CDR3 `spectratype`, each a single streamed `group_by` over a
  `io.scan_cohort` LazyFrame (peak memory independent of cohort size).
- **Diversity declines with age** — coverage-standardized Hill-number diversity via iNEXT
  (`sample_coverage`, `estimate_d`, `inext_batch`) that removes the sequencing-depth confound,
  plus the classic rarefaction/extrapolation curves (`rarefaction`).
- **Repertoires diverge with age** — pairwise exact-match overlap (`vdjtools.overlap`, the `F`
  metric on CDR3aa+V+J) after equal-depth downsampling, embedded with metric MDS: young samples
  cluster centrally and older donors scatter to the periphery.
- **Repertoires become more clonal with age** — top-clone read share.

```bash
pip install -e ".[examples,overlap]"          # overlap = vdjmatch/seqtree + scikit-learn (MDS)
examples/run.sh                               # interactive marimo editor
# or directly:
marimo edit examples/aging.py
marimo run  examples/aging.py                 # read-only served app
```

### Data

Auto-loads from the HuggingFace dataset
[`isalgo/airr_benchmark`](https://huggingface.co/datasets/isalgo/airr_benchmark) (folder
`vdjtools/`, full sequencing depth — **~0.5 GB total**), preferring a local **`./data_dump/`** copy
(gitignored — symlink your data there, e.g. `ln -s /path/to/airr_benchmark data_dump/airr_benchmark`).
The selected samples are ingested once into the **gitignored** `examples/.data/aging_nb/`
hive-partitioned Parquet cohort; HuggingFace verifies integrity and caches, so a re-run fetches
nothing. The `samples` slider trades coverage of the age range against runtime (overlap is O(n²)).

## `single_cell.py` — paired-chain single-cell TCR (10x dCODE)

A [marimo](https://marimo.io) notebook running the `vdjtools.sc` single-cell path on the
public **dCODE donor 4** dataset: `read_10x` (ingest 10x contigs), `chain_multiplicity`
(TRA/TRB presence-quadrant QC), `resolve_chains`/`pair_chains` (α/β receptors with
doublet handling), then an unsupervised **1-substitution β-CDR3 clustering** graded
against the dextramer antigen labels with `cluster_eval` (high purity/q-measure; a
shuffled labelling collapses the scores). Needs the `[sc]` + `[overlap]` extras:

```bash
pip install -e ".[examples,sc,overlap]"
marimo edit examples/single_cell.py
```

## `cdr_features.py` — CDR3 physicochemistry & k-mer features

A [marimo](https://marimo.io) notebook computing CDR3 amino-acid features with
`vdjtools.features` on the aging cohort: `physchem_profile` (hydropathy, charge,
volume, the 10 Kidera factors) per sample, correlated with donor age (the strongest is
CDR3 **hydropathy**, Spearman r ≈ −0.45), and `kmer_profile` 3-mer spectra embedded by
PCA and coloured by age. Pure polars:

```bash
pip install -e ".[examples]"
marimo edit examples/cdr_features.py
```

## `preprocess.py` — the repertoire preprocessing pipeline

A [marimo](https://marimo.io) notebook walking real Britanova samples (three sequencing
batches) through `vdjtools.preprocess`: `filter_functional` (drop non-coding),
`correct` (collapse PCR/sequencing-error variants), `downsample` (equalise depth),
`filter_frequency`/`filter_segment`, `decontaminate` (cross-sample bleed),
`pool_samples`/`join_samples`, and `correct_vj_usage` (VJ-usage batch-effect correction
— a before/after PCA where the batch separation collapses). Pure polars:

```bash
pip install -e ".[examples]"
marimo edit examples/preprocess.py
```

## `overlap_similarity.py` — exact / fuzzy / similarity-aware overlap

A [marimo](https://marimo.io) notebook contrasting three notions of repertoire
overlap with `vdjtools.overlap`, on a slice of the Britanova aging cohort: **exact**
(`Z = I`), **fuzzy** (`Z = 1[≤1 substitution]`), and **similarity-weighted** (the
TINA / Leinster-Cobbold form `pᵀZq` with a BLOSUM62 kernel `Z = exp(−P/τ)`) — the
continuous kernel that neither legacy vdjtools nor mirpy has. It builds all-pairs
distance matrices (`pairwise_distances`), embeds the cohort (`cluster_samples`, metric
MDS), shows where the similarity kernel finds graded overlap between repertoires that
share **no identical clonotype**, and runs a convergence test (`tcrnet`). Needs the
`[overlap]` extra (`vdjmatch`, `seqtree`, `scikit-learn`):

```bash
pip install -e ".[examples,overlap]"
marimo edit examples/overlap_similarity.py
```

## `emerson_biomarker.py` — CMV / HLA biomarker discovery (interactive)

A [marimo](https://marimo.io) notebook reproducing the core of Emerson et al.
(*Nat Genet* 2017) on the Emerson **HIP** cohort: an incidence-based **Fisher's
exact** screen (`vdjtools.biomarker.fisher_association`) for public TCRβ chains
associated with **CMV serostatus** or **HLA-A\*02**, validated live against a local
**VDJdb** dump by CMV epitope + HLA allele. The two options of the method are
interactive dropdowns — the **V/J-match requirement** (CDR3 / +V / +V+J) and
**exact vs 1-mismatch** CDR3 matching (metaclonotypes) — plus phenotype,
min-incidence, and the significance threshold. It rediscovers known CMV clones
(e.g. `CASSLAPGATNEKLFF` ↔ pp65 `NLVPMVATV` / HLA-A\*02:01) from raw repertoires.

```bash
pip install -e ".[examples,overlap]"          # overlap = vdjmatch, for the 1-mismatch option
marimo edit examples/emerson_biomarker.py
```

Data: a **balanced 400-subject subset** of [`isalgo/airr_hip`](https://huggingface.co/datasets/isalgo/airr_hip)
(the Emerson HIP cohort) auto-downloads into the gitignored `examples/.data/emerson_nb/` cache
(HuggingFace verifies integrity; a re-run fetches nothing). VDJdb validation resolves a
`./data_dump/` copy, else fetches the latest `antigenomics/vdjdb-db` release (both cached), and is
skipped gracefully if unavailable.

## `emerson_cmv_hla.py` — the same screen at cohort scale (streaming Fisher)

A [marimo](https://marimo.io) notebook running the CMV / HLA-A\*02 screen at scale: the cohort is
streamed into a hive-partitioned Parquet dataset one sample at a time (`ingest_cohort`), analysed
as one out-of-core LazyFrame (`scan_cohort`), and the per-feature Fisher tests are vectorised
through the hypergeometric tail — the cohort never fully in RAM (peak RSS reported per step). A
**subjects** slider scales toward the full 786; **CMV** is one-tailed, **HLA-A\*02** two-tailed,
each with an inline volcano and a VDJdb overlay, plus an optional 1-mismatch metaclonotype screen.
The interactive, condition×test×scope companion is [`emerson_biomarker.py`](emerson_biomarker.py).

```bash
pip install -e ".[examples,overlap]"
marimo edit examples/emerson_cmv_hla.py
```

## `scale_cohort.py` — large-cohort analytics (Parquet + streaming, flat memory)

A [marimo](https://marimo.io) notebook demonstrating the out-of-core cohort pattern on a
**synthetic** cohort you size with a slider: `ingest_cohort` streams every sample into a
hive-partitioned Parquet dataset one at a time (peak RSS ≈ one sample, reported), then every
statistic — the V-usage matrix, per-sample richness, an age-filtered count — is a `group_by`
collected with `engine="streaming"` over a single `scan_cohort` LazyFrame. Slide the sample count
up (toward thousands) and watch peak RSS stay flat. No data download.

```bash
pip install -e ".[examples]"
marimo edit examples/scale_cohort.py
```

## `vaccination_tracking.py` — longitudinal clonotype tracking + recapture model

A [marimo](https://marimo.io) notebook tracking clonotypes across **vaccination time
courses** (yellow-fever [`isalgo/airr_yfv19`](https://huggingface.co/datasets/isalgo/airr_yfv19),
influenza `isalgo/airr_flu_vac`, TBE `isalgo/airr_tbev_vac`) with `vdjtools.dynamics`: the paired
within-donor expansion test (`test_pair` → emergent / expanded / persistent / contracted /
vanishing) as sunken/alluvial and trajectory plots, metaclonotype-grouped testing, and the
**VDJtrack recapture model** (`capture_rates` / `capture_test`, Beta credible bands).

```bash
pip install -e ".[examples]"
marimo edit examples/vaccination_tracking.py
```

## `ankspond_motif.py` — the ankylosing-spondylitis "AS27" motif

A [marimo](https://marimo.io) notebook reproducing the Komech 2018 **TRBV9 / TRBJ2-3** CDR3β motif
in ankylosing spondylitis on [`isalgo/airr_ankspond`](https://huggingface.co/datasets/isalgo/airr_ankspond)
(60 donors), with the disease-vs-HLA-B27-carriage contrast — B27 is 26/27 confounded with AS, so
only the B27-matched comparison separates disease from carriage (AS/B27+ 16/26 vs HD/B27+ 1/12,
OR ≈ 17.6) — plus a metaclonotype family view.

```bash
pip install -e ".[examples,overlap]"
marimo edit examples/ankspond_motif.py
```

## `biomarker_explorer.py` — public-TCR association + co-occurrence (interactive)

A [marimo](https://marimo.io) notebook over the Emerson HIP cohort: `biomarker.association`
(condition × test × match-scope, with a live VDJdb overlay) plus a `biomarker.cooccurrence` panel
— the interactive superset of [`emerson_biomarker.py`](emerson_biomarker.py).

```bash
pip install -e ".[examples,overlap]"
marimo edit examples/biomarker_explorer.py
```

## `model_explorer.py` — recombination Bayes-net explorer

A [marimo](https://marimo.io) notebook exploring any bundled recombination model (OLGA vs learned):
the Bayes-net graph, per-event entropy, mutual information, and the marginal tables
(`vdjtools.model.analyze`). No download — uses the models shipped in the wheel.

```bash
pip install -e ".[examples]"
marimo edit examples/model_explorer.py
```

## `model_workshop.py` — build, fit, check, compare and score a model

The full [model workshop](https://docs.isalgo.dev/vdjtools/model.html) end to end: a custom germline
library → a model scaffold → EM with its training log → `check_model` → comparison against the
bundled TRB models (per-event Jensen-Shannon, gene usage, the comparison graph) → Pgen
distributions, log-likelihood and BIC on a held-out set → entropy and total diversity →
`extend_alleles` → `rescale_usage`.

Runs offline in seconds on a three-gene toy locus. A checkbox swaps in a small pre-annotated real
TRB subset from HuggingFace (needs access; falls back to the simulated set if unavailable).
Graphviz `dot` is optional — without it the comparison graph is skipped, nothing else changes.

```bash
pip install -e ".[examples]"
marimo edit examples/model_workshop.py
```

## `signature_features.py` — signature transforms, V-call resolution, the V+k-mer space

A [marimo](https://marimo.io) notebook on the three feature choices that are easy to get wrong,
each measurable in seconds on repertoires sampled from the bundled models (no download):
why the amino-acid block is **Anscombe arcsine** and not `log1p` of counts (`log1p` moves ~7.5×
across a depth change that leaves arcsine at 1.04×), why an ambiguous V call must go through
`resolve_gene` rather than `strip_allele`, and what `fit_kmer_space` builds — gapped/ungapped
patterns, BLOSUM62-clustered alphabets, TF-IDF, truncated SVD, and why its components must not be
selected by explained variance.

```bash
pip install -e ".[examples]"
marimo edit examples/signature_features.py
```

## `signature_features.py` — see also: feature presets

`vdjtools presets` lists the named, ranked feature sets and `vdjtools presets <name>` explains one
in full (what it contains, how it is computed, when to use it, and its caveats). To produce a table
over a whole dataset in parallel:

```bash
vdjtools presets                                                    # the ranked table
vdjtools signature *.tsv --preset statistics --threads 0 --out vsig.parquet
```

For the full vector — statistics **and** embedding geometry — use mirpy's `mir signature
--preset ...`, and see its `examples/feature_vectors.py` notebook.
