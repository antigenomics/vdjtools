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
- **arda germline = single source of truth**: all V/D/J germline + CDR3 anchors resolve from arda
  by allele name via `model.reference.load_germline` (arda's anchor convention is byte-identical to
  OLGA's: 0-based Cys104/[FW]118 offset into full germline). Never mix germline sources within a
  model — OLGA bootstrap models keep OLGA germline (exact-Pgen fidelity); arda-native (EM) models use
  arda. Raw anchor *indices* can differ by whole framework codons (IMGT drift) though the CDR3-region
  germline is identical — harmless as long as sources aren't mixed.
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
- **DONE 1a′** `model/reference.py` — **arda germline as source of truth**: `load_germline(locus,
  organism)` from arda (`cdr3fix.load_anchors` + `d_germlines.fasta`), `cut_segment` palindrome
  derivation (reproduces OLGA's cut segs), `reconcile_olga` audit. Shared frame verified vs OLGA
  (`tests/python/test_reference.py`). `from_olga` deliberately keeps OLGA germline (exact-Pgen
  invariant); arda is canonical for arda-native models + scenarios + stitching. arda is the `[model]`
  extra (`pip install -e ../arda` for dev). Gap: arda ships no full-length V/J germline (needs a
  helper) — a **P1c/stitching prerequisite**.
- **DONE aa Pgen** `model/pgen.py::pgen_aa` — codon-marginalizing left-to-right DP. VJ is fast;
  VDJ enumerates D placements + a multi-block DP handling the DJ insertion 3'→5' (reference impl:
  correct but 1–30 s/seq — the prime native-port target). Matches OLGA exactly (VJ 200 seqs; VDJ
  beta oracle 1.2036e-10 + shorts); `aa == Σ nt over synonymous` (`tests/python/test_pgen_aa.py`).
  **nt Pgen exactness fix**: V and J must each contribute ≥1 nt (never fully deleted) — OLGA excludes
  `len==0`; we now match nt+aa to ratio 1.0 (was ≤0.34% high on heavily-deleted seqs).
- **DONE 1e** `model/generate.py` — ancestral sampler → polars DataFrame; every draw scoreable
  (functional genes, ≥1 nt V/J); usage/length dists match OLGA (`test_generate.py`).
- **DONE 1d** `model/infer.py` — EM (E-step scenario soft counts; M-step polars normalize; align-init
  seeds gene usage). Closed-loop recovery of OLGA marginals (insertion 0.99, dinucl 0.999, Jmarg 0.94,
  aggP(delV) 0.97, germline-group V 0.93; per-allele V is germline-ambiguity-limited). `test_infer.py`.
- **DONE 1c** `model/stitch.py` — `stitch_contig(model,v,j,cdr3)` rebuilds full nt reads; `annotate`
  wraps `arda.annotate_sequences`. arda round-trips stitched synthetic contigs (junction + V/J gene),
  `test_stitch.py` (slow, needs arda+mmseqs). The plausible scenario *set* is what pgen/EM enumerate
  for the arda-called (V,J); arda supplies gene identification for real reads.
- **TODO 1f (native) — REQUIRED for VDJ, not just nice-to-have.** The **VDJ** paths (aa Pgen and the
  EM E-step `_accum_vdj`) are ~tens of s/seq in pure Python — the D × delD5 × delD3 × position
  enumeration × the many V candidates that all share the conserved Cys prefix. VJ is fine (~ms/seq).
  Port the validated hot loops to `src/` behind `_core` (PackedModel + pgen enumeration/DP, aa VDJ
  split-DP, sampler, EM E-step); assert C++ == Python == OLGA, benchmark faster than OLGA. Alternative
  that also fixes VDJ EM speed: **arda-masked E-step** — annotate reads with arda, restrict each read's
  scenario enumeration to arda's called (V,J,D) instead of enumerating all genes. Algorithms are all
  done + OLGA-validated in Python; this is the performance layer.
- **TODO D-D extension** (not in OLGA bootstrap): add `n_d`∈{0,1,2}, `d2_gene`, `d2_del`,
  `dd_ins`/`dd_dinucl` events + enumeration; the loader already emits `n_d`=δ(1). Ships with real
  tandem-D data (owner). arda full-length V/J germline helper needed for arda-native stitching.

Model schema notes: `ndel` is **biological** (neg = palindromic P-nt); dinucleotide row
`(from_nt,to_nt,p)=P(next|prev)` (OLGA's col-stochastic `R[next,prev]`); validation allows a group
to sum to 1 **or 0** (undefined conditional for an unused gene, kept for gene-index alignment).
Pgen/gen/EM invariant: **V and J each contribute ≥1 nt** to the CDR3 (OLGA-compatible).
