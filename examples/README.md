# vdjtools v2 examples

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
`vdjtools/`, full sequencing depth — **~0.5 GB total**), preferring a local `~/hf/` or `./` copy.
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
(the Emerson HIP cohort) auto-downloads into the gitignored
`examples/.data/emerson_nb/` cache (HuggingFace verifies integrity; a re-run
fetches nothing). VDJdb validation needs a local `vdjdb-db` slim dump and is
skipped gracefully if absent. The **full 786-subject, non-interactive** version is
[`emerson_cmv_hla.py`](emerson_cmv_hla.py) — run it with
`python examples/emerson_cmv_hla.py` (writes volcano plots + a vdjdb-validated hit
list; peak ~22 GB RAM at full scale).

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
