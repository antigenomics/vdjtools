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

**Next: Phase 1 — `feature/model-engine`** (native V(D)J recombination model, supersedes OLGA+IGoR):
- `python/vdjtools/model/`: `events schema reference io model scenario stitch pgen generate infer validate`.
- Model = directory of **long-format polars parquet** marginal tables + `manifest.json` (Bayes net
  declared as data via `events{}.given`). **D-D capable**: `p_nd(n_d∈{0,1,2})`, `p_d1`, `p_d2`, `DD`
  insertion junction. VJ loci degrade cleanly (no D tables).
- Scenario enumeration from **arda** best alignment → plausible set (slide V3'/J5' over a window;
  seed D via `arda.map_d_junction`; emit n_D∈{0,1,2}). `stitch_contig` rebuilds full nt reads from
  OLGA's (V,J,CDR3) output so synthetic + real reads share the arda→scenarios→EM path.
- **C++ (`_core`)**: Pgen DP (OLGA transfer-matrix factorization + D-D term; one DP for nt & aa),
  generation sampler, EM E-step. Python/polars: I/O, marginals, EM driver, validation.
- Bootstrap data: mirpy's OLGA models (7 loci) + ~100k OLGA-synthetic out-of-frame seqs/chain
  (`mirpy/mir/resources/olga/default_models/`, reference only). **No tandem D in bootstrap** — build
  D-D but unit-test it structurally; closed-loop oracle: EM must recover OLGA marginals, Pgen must
  match OLGA. Real out-of-frame data (all 7 loci) ships later from the owner.
- Watch: arda↔OLGA coordinate convention (1-based junction vs 0-based CDR3 grid) — one authoritative
  converter in `reference.py` with fixtures both ways; Pgen float64 underflow on long aa-CDR3.
