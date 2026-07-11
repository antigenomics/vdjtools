# vdjtools v2 examples

## `aging_airr_benchmark.py` — TCR repertoire aging

A [marimo](https://marimo.io) notebook that reproduces the classic aging signal in
the human TCR-beta repertoire (the legacy vdjtools *aging_lite* example: the
Britanova cohort, 41 donors, ages 6–90) using the basic-analytics layer of
vdjtools v2. It loads the native `.txt.gz` files with `vdjtools.io`
(`read_metadata` / `read_samples`), computes **coverage-standardized Hill-number
diversity** vs age with iNEXT (`inext_batch`, `estimate_d`, `sample_coverage`) —
showing that diversity declines with age (Spearman r ≈ −0.64) — draws
rarefaction/extrapolation curves for young vs old donors (`rarefaction`), and
summarizes age-associated **clonal expansion** (top-clone read share) straight
from the canonical frame.

### Run it

```bash
pip install -e ".[examples]"
marimo edit examples/aging_airr_benchmark.py     # interactive
# or headless:
marimo run examples/aging_airr_benchmark.py
```

### Data

The notebook's first cell auto-downloads its inputs from the HuggingFace dataset
[`isalgo/airr_benchmark`](https://huggingface.co/datasets/isalgo/airr_benchmark)
(folder `vdjtools_lite/`) into the **gitignored** `examples/.data/aging/`
directory. Every file is verified against the committed `aging_manifest.json`
(`{filename: md5}`): a file already present with the right md5 is skipped with no
network call, so a second run downloads nothing. The cache directory is never
committed.
