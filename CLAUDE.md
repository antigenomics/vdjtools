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
- **Never modify non-dependency libraries.** The only dependencies are `arda`, `vdjmatch`,
  `seqtree`. Everything else under `~/vcs/code/` (mirpy, IGoR, OLGA, pygor3, …) is
  **reference/oracle only — read-only**. If you find a bug in one, surface it (note it here or
  tell the owner); never edit it. Cross-validate against them; don't touch them.
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
- **DONE 1f (native `_core`)** — `include/vdjtools/model.hpp` (`PackedModel`, `Counts`) + `src/pgen.cpp`
  + `python/vdjtools/model/native.py` (`pack`, `pgen_nt`, `pgen_aa`) and `infer.py::infer_native`.
  `pack` reconstructs the polars model into dense C++ arrays; the hot loops are ported behind the
  pybind11 `_core` module. Verified exact (`tests/python/test_native.py`): native **nt Pgen** ==
  Python/OLGA (machine-eps), **89x** faster for VDJ (210→2.4 ms/seq); native **EM E-step** soft counts
  == Python (4e-16), **~100x** faster (TRA 13→0.1 s/it), VDJ **masked** EM now practical (~12 ms/seq).
  **arda-masked E-step**: `gene_masks`/`arda_masks` + `infer(masks=)`/`infer_native(masks=)` restrict
  enumeration to aligned genes (15x in Python; combines with native).
- **DONE 1g (native aa transfer matrix — VDJ *and* VJ)** `src/pgen.cpp` — replaced both enumerations
  with the Murugan/OLGA `Pi_L·Pi_R` split-DP: build left (V+VD/VJ-ins) and right partial sums once,
  stitch at the D placement (VDJ) or thread the J germline per J (VJ; J plays D's role). Key insight:
  cross-block coupling is codon-only (never the insertion Markov), so left/right factor cleanly. Exact
  vs OLGA machine-precision; **8.6x** faster on TRB, **1.9x** on TRA (`test_native_aa_vdj_matches_olga`,
  `appendix/bench_pgen.py`). The VJ enumeration had been ~1000x slower than OLGA — the TM fixes it.
- **DONE 1h (v/j-agnostic + 1-mismatch aa Pgen)** — the codon check is now a 64-bit **allowed-codon
  mask** per position (`ok_codon`), so wildcards/motifs run in one pass. `native.pgen_aa(m, aa, v, j,
  mismatches=)`: `v`/`j`=None marginalizes that gene (V/J-agnostic, always supported); `mismatches=1`
  sums the Hamming-1 ball via OLGA's inclusion-exclusion identity `Σ_k Pgen(a_{k→*}) − (L−1)Pgen(a)`
  but each term is one fast TM pass. Matches `compute_hamming_dist_1_pgen` to ~1e-15; **8.7x** faster on
  TRB, **2.5x** on TRA (`test_native_aa_hamming1_matches_olga`). Documented in `murugan_model.tex` §M.6.
- **TODO native perf gaps**: (a) **VJ / Hamming-1 codon-boundary sweep** — the 1-mm ball does L+1 TM
  passes; a forward/backward codon-boundary sweep would collapse them to ~1 pass for VJ loci (VDJ's
  D-placement sum couples positions, so L+1 is retained there). (b) native **generation sampler** (Python
  is already fast — low priority). (c) parallelize `estep_batch` over reads (GIL released) for Nx on multicore.
- **DONE model diagnostics** `model/analyze.py` — Bayes-net→graphviz DOT (nodes=marginal entropy H,
  edges=mutual information I; bnlearn-style, rendered via the `dot` CLI — no python-graphviz dep),
  `entropy_table`/`mutual_information`/`compare_entropy`, works on any Model. Cross-locus H table +
  I(V;J)=0 (VDJ independence made visible), I(delD5;delD3|D)≈1.18 bit (within-D conditional coupling,
  averaged over D — not the D-marginal). Single-parent factorizations only (raises on ≥2 parents). `test_analyze.py`.
- **DONE D-D consistency guards** (audit-driven): the not-yet-tandem paths — native `_core`
  (`native.pack`), amino-acid `pgen_aa`, `generate`, `infer`/`infer_native` — **raise
  `NotImplementedError` on a model with P(n_D=2)>0** (`dd.has_tandem`) instead of silently returning
  single-D. Only Python `pgen_nt` sums tandems. `prepare` rejects a malformed model (n_D=2 mass but no
  d2 tables). D-D model with p_nd2=0 stays byte-identical single-D (native included). `test_dd.py` guards.
- **DONE D-D Python reference** `model/pgen.py::_dd_middle` + `model/dd.py::to_dd` — n_D∈{1,2}
  enumeration (0-D folds into 1-D via a fully-trimmed D; tandem requires each D ≥1 nt → disjoint
  partition, resolves tandem-vs-long-insertion identifiability). `prepare` reads `n_d`/`d2_gene`/
  `d2_del`/`dd_ins`/`dd_dinucl` when present; `pgen_nt` weights P(n_D=1)·single + P(n_D=2)·tandem.
  Backward-compatible (p_nd2=0 == single-D, machine-eps on real TRD). Reference is correct-but-slow on
  TRD (the native port is the speed job). `test_dd.py` (tiny hand-checked model + TRD backward-compat).
  Real signal: **TRD 4.15% nonfunc / 3.42% func tandem-D** (28-bucket survey, new HF revision).
- **DONE appendix** `appendix/murugan_model.tex` — §M.4 tandem-D (Prop: disjoint n_D partition) + §M.9
  diagnostics (entropy/MI, Bayes-net figure `bn_trb.pdf`/`bn_trd_dd.pdf`, cross-locus H table). 9 pages.
- **TODO native D-D port** (task): add `p_nd`, `d2`/`dd_ins`/`dd_dinucl` to `PackedModel`; extend
  `pgen_nt`/`pgen_aa`/`estep_batch` over n_D∈{1,2} (one extra D block in Πᵣ, weighted by P(n_D=2)).
  Match `_dd_middle`. Then EM on real TRD learns P(n_D=2)>0. arda full-length V/J germline helper still
  needed for arda-native stitching.
- **TODO real-data EM comparison** (task): `infer_native` on real nonfunc reads (TRB, TRD) → compare
  inferred vs legacy-OLGA marginals via `analyze` (the "ours vs OLGA on the same data" deliverable).

Model schema notes: `ndel` is **biological** (neg = palindromic P-nt); dinucleotide row
`(from_nt,to_nt,p)=P(next|prev)` (OLGA's col-stochastic `R[next,prev]`); validation allows a group
to sum to 1 **or 0** (undefined conditional for an unused gene, kept for gene-index alignment).
Pgen/gen/EM invariant: **V and J each contribute ≥1 nt** to the CDR3 (OLGA-compatible).
