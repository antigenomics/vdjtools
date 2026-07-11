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
