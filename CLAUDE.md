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
- `appendix/` — **library-only**: the LaTeX theory appendix (`murugan_model.tex` + `refs.bib` +
  `bn_*.pdf`), `build_bundled_models.py` (builds the models *shipped in the wheel*) and
  `concordance.py` (validates them against the OLGA oracle). Nothing that merely *uses* the library.

## Repo split (2026-07-17) — benchmarks live elsewhere
**`~/vcs/projects/2026-vdjtools-benchmark`** now holds the benchmark/campaign work: `bench/` (the
`bm_*`/`bench_*` scripts, the confound gates — `accept_gate`, `depth_gate`, `theta_ceiling`,
`cooccurrence_fpr`, `assoc_depth` — and `validate_batch_covid.py`) plus `scripts/`
(`biomarker_bench.sbatch`). They are analyses that *use* vdjtools, not part of it — the same split
mirpy made between `~/vcs/code/mirpy` and `~/vcs/projects/2026-mirpy-analysis`. This repo stays
library + tests + docs; the **numbers of record stay in `SOURCES.md` here** (the Phase-6b tables),
which now cite the benchmark repo's paths.

⚠ **As of the move it is not a git repo** — nothing there is under version control yet.
⚠ The scripts hardcode cluster paths (`/projects/biomarkers/{raw,results}`, `/projects/fmba_covid`)
and are run ad hoc; they were never importable from here, so nothing in the library broke.

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
Deps trace to real imports — a phase's parent is promoted from an extra to a base dep once it
lands. **All three antigenomics engines are BASE deps**: `arda-mapper` (v2.4.0), `seqtree` +
`vdjmatch` (v2.5.0). Rule: *if the README advertises it, a plain `pip install vdjtools` must
deliver it.* Every advertised capability delegates to one of the three and none has a fallback —
arda → germline reference + model engine (mirpy depends on vdjtools *for* the germline); seqtree →
fuzzy search/e-values (`preprocess.correct`, similarity overlap, TCRnet); vdjmatch → overlap,
TCRnet, metaclonotypes. All are imported lazily, so `import vdjtools` stays light; between them
they add only `requests`. mmseqs2 is needed only for arda's annotate path. `[model]`/`[preprocess]`
are kept-but-empty aliases so existing pins still resolve; `[overlap]` now carries only
scikit-learn (MDS; `hclust` works without it). `test_smoke.py` pins this contract.

**Heavy tests run on Aldan-3, not the Mac.** Unit tests stay local + fast; anything heavy —
benchmarks (`RUN_BENCHMARK=1`), full-locus EM/concordance, large-data (100k+ read) runs — goes to
the Aldan-3 HPC cluster via the `aldan3` CLI (repo `../aldan3-client`) instead of the 32 GB laptop.
Submit + monitor deterministically (every subcommand takes `--json`):
`aldan3 slurm submit <script.sh> [-- ARGS…] [--env <e>] [--cpus/--mem/--time/--gpus …]`,
then `aldan3 slurm queue` · `log <id>` · `hist <id>` (sacct usage) · `cancel <id>`;
`aldan3 slurm template cpu|gpu|array -o job.sbatch` scaffolds a starter script. The runnable
`scripts/*.sbatch` now live in the **benchmark repo** (see the Repo split above), not here;
`aldan3` is just the external driver.

## Git model
`master` = v2 (tagged releases) ← `dev` (integration) ← `feature/*` (one per phase).
**Legacy v1.x is on the `legacy-1.x` branch and under tags `v0.0.1`..`1.2.1`** — do not disturb.
The v2 history is an orphan root (no shared ancestry with legacy). Carried-over legacy resource
files (`aa_property_table.txt`, `cdr3contact.txt`, `vj_families.txt`) and format-conversion test
fixtures live on `legacy-1.x` (`src/test/resources/samples/`), pull them over when a phase needs them.

**Worktrees convention**: develop each phase/feature in its own git worktree, not in the main
checkout — `git worktree add .claude/worktrees/<name> -b feature/<name>` (one worktree ↔ one
`feature/*` branch). This isolates parallel agents from `master` and from each other; never run
two features in one worktree. `.claude/` (including `.claude/worktrees/`) is gitignored — never
commit it. Consolidate a finished phase by merging its `feature/*` branch into `dev`, then
`git worktree remove .claude/worktrees/<name>`.

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
  invariant); arda is canonical for arda-native models + scenarios + stitching. arda is a **base
  dep** (`pip install -e ../arda` for dev). Gap: arda ships no full-length V/J germline (needs a
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
  Python/OLGA (machine-eps) — VDJ nt was originally a per-scenario D enumeration (~28 ms/seq TRB, i.e.
  ~6× *slower* than OLGA; a docs figure once wrongly claimed 2.4 ms/89×) and is now routed through the
  aa transfer matrix — see the nt-TM bullet below; native **EM E-step** soft counts
  == Python (4e-16), **~100x** faster (TRA 13→0.1 s/it), VDJ **masked** EM now practical (~12 ms/seq).
  **arda-masked E-step**: `gene_masks`/`arda_masks` + `infer(masks=)`/`infer_native(masks=)` restrict
  enumeration to aligned genes (15x in Python; combines with native).
- **DONE 1g (native aa transfer matrix — VDJ *and* VJ)** `src/pgen.cpp` — replaced both enumerations
  with the Murugan/OLGA `Pi_L·Pi_R` split-DP: build left (V+VD/VJ-ins) and right partial sums once,
  stitch at the D placement (VDJ) or thread the J germline per J (VJ; J plays D's role). Key insight:
  cross-block coupling is codon-only (never the insertion Markov), so left/right factor cleanly. Exact
  vs OLGA machine-precision; **8.6x** faster on TRB, **1.9x** on TRA (`test_native_aa_vdj_matches_olga`,
  `~/vcs/projects/2026-vdjtools-benchmark/bench/bench_pgen.py`). The VJ enumeration had been ~1000x slower than OLGA — the TM fixes it.
- **DONE 1h (v/j-agnostic + 1-mismatch aa Pgen)** — the codon check is now a 64-bit **allowed-codon
  mask** per position (`ok_codon`), so wildcards/motifs run in one pass. `native.pgen_aa(m, aa, v, j,
  mismatches=)`: `v`/`j`=None marginalizes that gene (V/J-agnostic, always supported); `mismatches=1`
  sums the Hamming-1 ball via OLGA's inclusion-exclusion identity `Σ_k Pgen(a_{k→*}) − (L−1)Pgen(a)`
  but each term is one fast TM pass. Matches `compute_hamming_dist_1_pgen` to ~1e-15; **8.7x** faster on
  TRB, **2.5x** on TRA (`test_native_aa_hamming1_matches_olga`). Documented in `murugan_model.tex` §M.6.
- **DONE 1i (native nt Pgen via the aa transfer matrix — single-D *and* D-D)** `src/pgen.cpp::pgen_nt`
  — an in-frame nt CDR3 (length always ×3, Cys→Phe/Trp) is exactly an aa query with a **singleton
  allowed-codon mask** per position, so the same `pgen_aa_masked` Pi_L·Pi_R DP gives the identical value
  far faster than the per-`(V,J,delV,delJ)` D enumeration it used before. `pgen_aa_masked` mixes the
  D-count prior itself (`p_nd1`·single-D + `p_nd2`·`pgen_aa_vdj_dd`), so **both single-D and D-D nt route
  through the TM**. Gated only on `m.vdj && N%3==0`; the enumeration (`pgen_nt` outer loop → `d_middle`/
  `dd_middle`) is retained solely for non-in-frame nt (`N%3≠0`, never a real CDR3) and as the oracle.
  Speed: single-D VDJ nt **0.53 ms/seq TRB** (was ~28 ms → **53×**, **9×** faster than OLGA's 4.8 ms);
  D-D nt **~15 ms/seq TRD** (was ~350 ms via `dd_middle` → **24×**). Exact vs OLGA on all 7 loci
  (`test_native.py::test_native_matches_python` rtol 1e-9; slow exhaustive + `appendix/concordance.py`
  nt r(log10)=1.0) and vs the D-D enumeration/`_dd_middle` reference (real TRD max-rel 3.8e-15;
  `test_dd.py::test_aa_dd_equals_nt_sum` now pins native D-D nt per synonymous codon incl. tandem-only
  `CHHF`). NB: the earlier "89× / 2.4 ms" nt figure was never real — see the 1f note.
- **DONE native pgen batch-parallelization** `src/pgen.cpp::pgen_aa_batch` + `native.pgen_aa_batch(model,
  seqs, v=, j=, mismatches=, threads=)` — Pgen / 1-mm ball over many CDR3s, partitioned across worker
  threads (GIL released, disjoint writes → **bitwise-identical** to per-sequence, thread-count-invariant).
  **11.3× exact / 11.6× 1-mm on 16 cores**. GIL also released on the single `pgen_aa`/`pgen_aa_hamming1`
  bindings. This is the exact real-workload speedup (Pgen over many clonotypes). `test_native_pgen_batch.py`.
- **TODO native perf gaps**: (a) **VJ / Hamming-1 codon-boundary sweep** (collapse the 1-mm L+1 TM
  passes to ~1) — **set aside**: `pgen_aa_vj`'s V/J combine boundary migrates with the delJ sum, so a
  wildcarded codon has no clean O(1)-per-codon leave-one-out; a forced rewrite risks the exact-Pgen
  invariant for a non-bottleneck. Batch parallelization (above) is the exact win instead. (b) native
  **generation sampler** (Python is already fast — low priority). (c) `estep_batch` read-parallelization
  is done. (d) **DONE `overlap.pairwise_distances` O(n²) re-hashing** (`overlap/cluster.py`): the exact
  metrics (F/F2/D/R/jaccard) now aggregate each sample's clonotype key→freq ONCE (`metrics._overlap_from_agg`
  factored out of `overlap_pair`), then every pair is a join over the pre-aggregated frames — previously each
  pair re-aggregated both frames, so each sample was collapsed (n-1)×. Distances **bitwise-identical** (same
  join + numpy; `array_equal` across all 5 metrics), **3.2× on n=40×4k** (grows with cohort depth/raw rows).
  Fuzzy (`scope=`) / `similarity_*` still delegate to vdjmatch/seqtree per pair (engine-indexed, out of scope).
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
- **DONE native D-D nt Pgen** `src/pgen.cpp::dd_middle` + `PackedModel` (`p_nd1/p_nd2`, `pd2_given_d1`,
  `del_d2`, `ins_dd`/`R_dd`/`bias_dd`, `dd` flag). Factorized per-D1 left/right partial sums + O(N²)
  DD-insertion sweep (O(nD²L²N + nD·N²), not naive O(nD²L⁴N²)). Exact vs Python reference on the tiny
  model; ~255 ms/seq on real TRD (naive was seconds/seq → timed out). `native.pgen_nt` supports tandems
  (`pack` no longer guards); OLGA cannot compute D-D Pgen at all. `test_native_dd_matches_python_reference`.
  **Superseded for in-frame nt by the aa transfer matrix (1i, ~15 ms/seq, 24×)** — `dd_middle` now runs
  only for non-in-frame nt and as the correctness oracle; it stays exact and is the reference the TM checks against.
- **DONE real-data EM comparison** `~/vcs/projects/2026-vdjtools-benchmark/bench/bench_em.py` — `infer_native` on real nonfunc TRB+TRD reads
  (single-D, arda-masked) vs legacy OLGA via `analyze`. Finding: **real repertoires have broader
  trim/insertion entropy than OLGA's synthetic model** (TRB d_del 6.4→7.6 bit, vd_ins 3.8→4.5); within-D
  coupling I(delD5;delD3|D)≈1.1 bit robust across both. Held-out loglik improves.
- **DONE native D-D E-step** `src/pgen.cpp::accum_dd` + `Counts` (`n_d, d2_gene, d2_del, ins_dd,
  dinucl_dd`) + `_mstep_native`/`infer_native` (guard removed). Factorized forward/backward: same per-D1
  A[e1]/B[s2] partial sums as `dd_middle` (so the soft-count normalizer == Pgen by construction), plus a
  backward message C[e1]=Σ insDD·B and forward message Dmsg[s2]=Σ A·insDD from the combine sweep;
  re-enumerating each block once attributes every per-realization soft count. Soft counts == Python
  `_accum_dd` **exactly** (2e-16, all 15 events); closed-loop native EM recovers P(n_D=2) (0.0565 vs
  0.058 empirical on 500 synthetic TRD reads); **~370x** faster than naive (68 ms/read masked vs
  >25 s/read timeout). Single-D unchanged (p_nd1=1 → n_d renorm to δ(1)). `test_native_dd_estep_matches_
  python_reference`, `test_native_dd_em_recovers_p_nd2`. **Bug found+fixed** (pre-existing): `native.pack`
  cached by `id(model)` returned a stale `PackedModel` after CPython id-reuse (TRB→TRD in one process
  crashed the M-step, 89 vs 18 V) — now stores+verifies the model ref.
- **DONE real TRD D-D EM** (`appendix/real_trd_dd` pattern) — native D-D EM on real out-of-frame TRD
  reads (arda-mapped 11.7k unique clonotypes; arda's own d2_call rate **4.15%**) learns generative
  **P(n_D=2)=0.028** (1k reads, arda V+D-masked, 26 ms/read, held-out LL −40289→−35431). Below arda's
  hard-call rate as expected — EM marginalizes the tandem-vs-long-insertion ambiguity, so the generative
  tandem prob is more conservative than the alignment flag. (Convergence gated on V-usage TV → stops at
  2 iters under masking; an n_d-focused stop would refine the number.)
- **D-D skip shortcut — benchmarked, no exact skip exists.** D-D Pgen is 2.5× single-D (`dd_middle` ≈60%
  of the compute) but **0%** of reads have a zero D-D contribution (two 1-nt Ds + insertions tile almost
  any mid), so there is no exact per-read skip. Dropping `dd_middle` is <1% median Pgen error but up to
  100% on genuinely-tandem reads. n_D=2 mass is concentrated (top 20% of reads = 85%); a length-gate that
  skips the shortest 25% biases learned P(n_D=2) by −2.5% (skipping 50% → −11%). Conclusion: keep D-D
  exact by default (correctness ethos); the clean **exact** speedup is read-parallelization, not gating.
- **DONE native aa D-D + Hamming-1 + v/j-agnostic** `src/pgen.cpp::pgen_aa_vdj_dd` — tandem aa Pgen
  transfer matrix: Lf (V+insVD) reused; per-J D-less right DP (`mk_right_tm` D=−1); D1 middle-left `Mf`
  (thread D1 germline out of Lf + `extend_ins_into` through insDD); D2 threaded + `combine_tm`. J looped
  explicitly (P(D1|J) couples D1↔J). `pgen_aa_masked` mixes p_nd1·single + p_nd2·tandem. Hamming-1 and
  v/j-agnostic fall out of the same masked call. **No Python fallback** (removed). native == Python
  `_dd_aa_middle` == Σnt exact on the codon-aligned tiny model; Hamming-1 == brute ball; fast on real TRD
  (fixed 0ms/agnostic 10ms/ham1 60ms) where Python is intractable. Also added `pgen._dd_aa_middle` (oracle).
- **DONE native estep threading** — `estep_batch(…, threads=0)` partitions reads over workers (private
  Counts, fixed-order reduce, GIL released); <64-read batches stay single-threaded (bitwise-exact tests).
  **6.7× on 8 threads**, exact (Δcount 2e-13). This is the exact speedup that replaces a biased read-gate.
- **DONE D-D default for the D-bearing loci** — `infer`/`infer_native(single_d=False, p_nd2_init=0.02)`
  promote a single-D template to D-D for `DD_DEFAULT_LOCI={TRB,TRD,IGH}` (via `_maybe_promote_dd`→`to_dd`);
  EM learns P(n_D=2) out of the box. `single_d=True`, VJ loci, already-tandem all unchanged. `from_olga`
  stays single-D (exact OLGA fidelity — the default applies at inference).
- **DONE 7-chain concordance** `appendix/concordance.py` — native nt & aa Pgen vs the OLGA oracle across
  all 7 loci: **r(log10 Pgen)=1.00000 everywhere**, max-rel ~1e-14–1e-16 (aa==Σnt confirmed). The two
  larger outliers (TRG nt 1.6e-2, TRD nt 5.7e-4) are deep-tail sequences (Pgen ~1e-28) where FP
  summation order dominates — absolute agreement ~1e-30; `test_pgen_nt` proves exactness on all 7 loci.
- **DONE bundled models + loader** `model/bundled.py` (`load_bundled(locus, source)`, `list_bundled`) —
  ship all 7 loci × {`olga`, `learned`} in the wheel (`model/_bundled/`, ~1 MB; scikit-build-core packs
  them automatically). `learned` = native EM on the FULL real HF **non-functional** read set —
  out-of-frame AND stop-codon, since both escaped selection and keeping only out-of-frame conditions
  the training set on junction length mod 3 (`appendix/build_bundled_models.py`; **no cap, no
  subsampling** — every clonotype surviving the germline filter, and the printed `n_used == n_clono`
  is the check). ⚠ The old *"2k clonotypes/locus, held-out LL improves on every locus"* was wrong
  twice: the cap was real, and the LL was the EM's **own training objective**, which EM increases
  monotonically by construction — it validated nothing. Real held-out + oracle comparison:
  `appendix/compare_models.py`. The bundled **`olga`** models come from the repo's own
  `tests/python/fixtures/olga/default_models` (all 7 human loci; pip olga ships only 5 and no
  TRG/TRD — those two trace to mirpy `legacy-v2` commit aeccd75, verified byte-identical).
- **DONE arda-anchored D-D learning** — unregularized D-D EM over-attributes tandems on real data
  (identifiability; TRB→0.28). Two regularizers, both native==Python exact: **`dd_allowed`** per-read gate
  (a read may be n_D=2 only where arda called a `d2_call`) and **`nd_prior`** Dirichlet single-D pseudocount.
  `infer/infer_native(..., dd_allowed=, nd_prior=)`; native `estep_batch(..., dd_allowed)`. Anchored learned
  D-loci: TRB **0.000**, TRD **0.006**, IGH **0.009** (plausible; arda hard-call ~4%). `test_dd_anchor_and_prior…`.
- **DONE marimo explorer** `notebooks/model_explorer.py` — reactive Bayes-net/entropy/MI/marginal explorer
  over any bundled model (OLGA vs learned); `[examples]` extra. README/docs/SOURCES updated.
- **DONE V-zeroing fix (2 root causes)** — the shipped learned models zeroed 68/89 TRB V alleles; NOT a
  masking/soft-realign problem (soft realign avalanches mass onto the most permissive germline,
  IGKV3-20 0.10→0.74 — removed). (1) `infer._align_init` collapsed germline-identical paralogs: `max(...)`
  takes the first of an exact tie (TRBV6-2/6-5/6-6, IGKV2-28/2D-28 tie identically), seeding the rest at
  P(V)=0 which the E-step's `if pv==0: continue` makes absorbing — fixed by splitting each read's vote
  across all tied genes. (2) `io.from_olga(derive_orf=True)` (opt-in, builder-only) reconstructs the CDR3
  germline OLGA leaves empty for ORF alleles (TRBV23-1, 8.6% of real TRB) from `full[anchor:]`; the oracle
  default keeps it off (exact-OLGA-Pgen invariant). Result: **0 functional genes zeroed**, Pearson(arda,
  learned) 0.97 (TRB). `test_infer.py::test_align_init_splits_germline_identical_ties`,
  `test_model_loader.py::test_derive_orf_is_opt_in_and_preserves_the_oracle`. Learned models rebuilt (masked, arda-anchored).
- **TODO** arda full-length V/J germline helper still needed for arda-native stitching (P1c residual); the
  `derive_orf` reconstruction covers the ORF-usage case but not full-length stitching.

**AS/B27 motif campaign — `feature/as-b27-motif`** (`~/vcs/projects/2026-vdjtools-benchmark/bench/bm_ankspond.py`, runs locally in
~27 s; HF `isalgo/airr_ankspond`, 60 donors). Reproduces Komech 2018's TRBV9/TRBJ2-3 motif and
fixed two real bugs on the way:

- **`model/native` Pgen allele guard (was a silent wrong answer)** — `vi.get(v, -1)` mapped any
  unrecognised V/J to `-1` = *marginalize over all V/J*. The model is keyed by **allele**, and real
  repertoires carry **gene-level** `v_call` (`TRBV9`), so `pgen_aa(m, cdr3, "TRBV9", "TRBJ2-3")`
  returned the V/J-agnostic value — **2.38× too high** — and raised nothing. `_gene_idx` now raises
  and names the alleles to pass. Exact values unchanged. If you pin an older vdjtools, pass `*01`.
- **`io.read_mitcr`** — the MiTCR/tcR dotted dialect (`Read.count`, `CDR3.nucleotide.sequence`,
  `V.gene`). ankspond `old/` (the actual 2018 cohort, and the only part with nucleotide CDR3 +
  V/D/J markup) previously **raised** in `sniff_format`. `_lower_map` is exact-lowercase, so
  MiGEC's space-separated picks never match dotted headers — it needs its own reader.
- **`features.kmer` is no longer descriptive-only**: `flank` drops the conserved anchors (verified
  == `seqtree.seeds.core_kmers` over 840 comparisons), and `kmer_cohort` → `association(key=
  ("v_call","kmer"), match="exact")` is the V+k-mer test. `_feature_frame` no longer requires
  `junction_aa` for `match="exact"` (fuzzy/1mm still *search* on it and still do); a key of
  germline calls **alone** still raises — that is `stats.segment_usage`, not a biomarker.
  ⚠ `str.len_chars()` is **UInt32** — `len - 2*flank` underflows on short junctions; cast first.

Findings worth not re-deriving: **B27 is 26/27 confounded with AS** in this cohort, so only the
**B27-matched** contrast separates disease from carriage (AS/B27+ **16/26** vs HD/B27+ **1/12**,
OR=17.6, p=0.0023; batch-matched 26 vs 7 → OR=9.6, p=0.035). B27 carriage among *healthy* is
**null** (p=0.60) ⇒ disease, not carriage. **V-pinning is load-bearing** (unpinned, the healthy arm
gains 4 wrong-V convergents). Not depth (MWU p=0.94). **38 of 40 `old/` donors reappear in `new/`**
⇒ no independent replication exists in this dataset; never pool. At 26-vs-12 **BH cannot clear
0.05** over 273 pinned features (min attainable Fisher p = 3.6e-3) — the covid19 lesson again; the
*ranking* is the result (motif at ranks 1,2,5,7,13; V+4mer `VGLY` rank 1, OR=25.0, beating the best
single clonotype). **VDJdb release matters**: the 2024-06 checkout has **zero** records for this
motif; 2025-12-29 has 7 (Yang 2022 *Nature*, B\*27:05, self + *E. coli* epitopes) — but that oracle
is **partly circular** (same group, plausibly these donors). Handoff plan for mirpy: `~/vcs/projects/2026-mirpy-analysis` branch `as-b27-embedding`.
The campaign script moved out with the rest of the benchmarks (see the Repo split above).

Model schema notes: `ndel` is **biological** (neg = palindromic P-nt); dinucleotide row
`(from_nt,to_nt,p)=P(next|prev)` (OLGA's col-stochastic `R[next,prev]`); validation allows a group
to sum to 1 **or 0** (undefined conditional for an unused gene, kept for gene-index alignment).
Pgen/gen/EM invariant: **V and J each contribute ≥1 nt** to the CDR3 (OLGA-compatible).
