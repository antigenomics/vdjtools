# CLAUDE.md — vdjtools v2

## What this is
Clean-room **Python + C++** rewrite (v2.0.0, GPL-3.0) of the legacy Groovy/Java vdjtools.
TCR/BCR repertoire analysis on the **AIRR schema + polars**, minimal OO, built on the
antigenomics ecosystem: **seqtree** (fuzzy search / e-value engine), **vdjmatch** (overlap +
TCRnet), **arda** (AIRR annotation + markup repair; brings conda/mmseqs2).

## Layout
- `python/vdjtools/` — package (src-layout via `wheel.packages`). Subpackages: `io model stats
  features overlap preprocess biomarker sc cli` (lazy-loaded; `import vdjtools` pulls no heavy deps).
- `src/`, `include/vdjtools/` — C++ core; `src/_bindings.cpp` → the `vdjtools._core` pybind11 ext.
  **Native only for Pgen DP, generation sampler, EM E-step** (Python-first everywhere else).
- `tests/{cpp,python}/`, `docs/` (Sphinx + pydata, gh-pages), `.github/workflows/{ci,docs,publish}.yml`.
- `CMakeLists.txt` (scikit-build-core + pybind11, C++20), `environment.yml` (conda), `setup.sh`.

## Build / test / run
```bash
conda env create -f environment.yml && conda activate vdjtools   # or reuse .venv
pip install -e ".[dev,test]"                                      # builds _core
pytest tests/python -q
cmake -S . -B build -DVDJTOOLS_TESTS=ON && cmake --build build && ctest --test-dir build
sphinx-build -W --keep-going -b html docs docs/_build/html        # docs gate (zero warnings)
```
Co-developed parents are early-alpha: `bash setup.sh --dev-parents` editable-installs
`../seqtree ../arda ../vdjmatch` if present (else PyPI: `seqtree`, `arda-mapper`, `vdjmatch`).
Deps trace to real imports — `arda-mapper`/`vdjmatch` live in the `[model]`/`[overlap]` extras
until the phase that imports them promotes them to base deps.

## Git model
`master` = v2 (tagged releases) ← `dev` (integration) ← `feature/*` (one per phase).
**Legacy v1.x is on the `legacy-1.x` branch and under tags `v0.0.1`..`1.2.1`** — do not disturb.
The v2 history is an orphan root (no shared ancestry with legacy). Carried-over legacy resource
files (`aa_property_table.txt`, `cdr3contact.txt`, `vj_families.txt`) and format-conversion test
fixtures live on `legacy-1.x` (`src/test/resources/samples/`), pull them over when a phase needs them.

## Conventions
- AIRR Rearrangement/Cell + polars `pl.DataFrame` in and out; minimal OO (thin index classes only).
- Delegate rather than reimplement: overlap/TCRnet → vdjmatch (`cluster.overlap`,
  `evalue.query_evalues`); annotation/markup/scenarios → arda; search/e-value → seqtree.
- Native code goes through the single `_core` ext. Flip `editable.rebuild=true` in pyproject
  during C++-heavy work for recompile-on-import (needs the `build-dir` already set).

## Open loops / next steps
Phase 0 (scaffold, git surgery, CI) is **done and pushed**. See `ROADMAP.md` for all phases and
`SOURCES.md` for data provenance. The full approved design lives in the session plan file
(`~/.claude/plans/i-want-to-complenely-gleaming-snail.md`).

**Phase 1 — `feature/model-engine`** (native V(D)J model, supersedes OLGA+IGoR). Model = directory
of **long-format polars parquet** marginal tables + `manifest.json` (Bayes net declared as data via
each event's `given`). VJ loci degrade cleanly (no D tables). Bootstrap data: mirpy's OLGA models
(7 loci) + OLGA-synthetic out-of-frame seqs (`mirpy/mir/resources/olga/default_models/`, ref only;
**no tandem D in bootstrap**). Progress on branch:

- **DONE 1a** `model/{events,schema,model,io}.py` — OLGA→polars loader, `manifest.json`, parquet
  round-trip. Lossless vs OLGA's arrays across all 7 loci (`tests/python/test_model_loader.py`).
- **DONE 1b (nt)** `model/pgen.py` — reference **nucleotide** Pgen (direct scenario sum over the
  tables, no OLGA at runtime). Matches OLGA `compute_nt_CDR3_pgen` exactly, all 7 loci
  (`tests/python/test_pgen_nt.py`; exhaustive check is `-m slow`). This is the quantity EM needs.
- **TODO aa Pgen** — needs the transfer-matrix codon-marginalizing DP; do it in the native `_core`
  port (P1f) which serves nt+aa in one DP. Reference: OLGA `generation_probability.py` (spec was
  extracted; grid is (4,3L), split-point dot product, `Tvd/Svd/Dvd/lTvd/lDvd` insertion matrices).
- **TODO 1c** `scenario.py` + `stitch.py` — `stitch_contig(V,J,CDR3)` rebuilds full nt reads (OLGA
  emits only V+J+CDR3); **arda** best alignment → plausible scenario set (slide V3'/J5' window; seed
  D via `arda.map_d_junction`; emit n_D∈{0,1,2}). Needs the **D-D schema extension**: add `n_d`,
  `d2_gene`, `d2_del`, `dd_ins`/`dd_dinucl` events (the loader already emits `n_d`=δ(1) for OLGA VDJ).
  Watch: arda↔OLGA coordinate convention (1-based junction vs 0-based CDR3 grid) — one authoritative
  converter, fixtures both ways.
- **TODO 1d** `infer.py`/`validate.py` — EM (E-step soft counts over scenarios; M-step polars
  group-normalize). Closed-loop oracle: EM must recover OLGA marginals on the synthetic bootstrap.
- **TODO 1e** `generate.py` — ancestral sampler → polars DataFrame.
- **TODO 1f** native `_core`: fast transfer-matrix Pgen (nt+aa, +D-D term), sampler, EM E-step;
  assert C++ == reference-Python == OLGA, benchmark faster than OLGA. `arda` is a hard dep here.

Model schema notes: `ndel` is **biological** (neg = palindromic P-nt); dinucleotide row
`(from_nt,to_nt,p)=P(next|prev)` (OLGA's col-stochastic `R[next,prev]`); validation allows a group
to sum to 1 **or 0** (undefined conditional for an unused gene, kept for gene-index alignment).
