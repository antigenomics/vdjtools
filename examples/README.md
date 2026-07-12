# vdjtools v2 examples

## `aging_airr_benchmark.py` — TCR repertoire aging

A [marimo](https://marimo.io) notebook that reproduces the classic aging signals
of the human TCR-beta repertoire on the full-depth Britanova **"Cord Blood to
Centenarians"** cohort (the legacy vdjtools aging example: 79 donors, ages 0–103,
including 8 cord-blood samples) using the basic-analytics layer of vdjtools v2. It
loads the native `.txt.gz` files with `vdjtools.io` (`read_metadata` / `vio.read`)
and shows three things:

- **Diversity declines with age** — coverage-standardized Hill-number diversity
  via iNEXT (`inext_batch`, `estimate_d`, `sample_coverage`), Spearman r ≈ −0.71.
- **Repertoires diverge with age** — pairwise repertoire overlap
  (`vdjtools.overlap`, the exact-match `F` metric) after equal-depth downsampling,
  embedded with metric MDS: cord-blood samples cluster centrally and older donors
  scatter to the periphery (distance-from-centroid vs age r ≈ +0.64).
- **Repertoires become more clonal with age** — top-clone read share, r ≈ +0.70 —
  plus the classic rarefaction/extrapolation curves (`rarefaction`).

### Run it

```bash
pip install -e ".[examples]"
examples/run.sh                              # interactive marimo editor
# or directly:
marimo edit examples/aging_airr_benchmark.py
marimo run  examples/aging_airr_benchmark.py # read-only served app
```

### Data

The notebook's first cell auto-downloads its inputs from the HuggingFace dataset
[`isalgo/airr_benchmark`](https://huggingface.co/datasets/isalgo/airr_benchmark)
(folder `vdjtools/`, full sequencing depth — **~0.5 GB total**) into the
**gitignored** `examples/.data/aging/` directory. Every file is verified against
the committed `aging_manifest.json` (`{filename: md5}`): a file already present
with the right md5 is skipped with no network call, so a second run downloads
nothing. The cache directory is never committed.

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
