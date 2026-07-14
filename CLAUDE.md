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
Deps trace to real imports — a phase's parent is promoted from an extra to a base dep once it
lands. **`arda-mapper` is a BASE dep** (since v2.4.0): a plain `pip install vdjtools` must ship
the model engine *and* the germline reference, because downstream libs (mirpy) depend on vdjtools
*for* that reference. It's imported lazily, so `import vdjtools` stays light; mmseqs2 is needed
only for arda's annotate path. `vdjmatch`/`seqtree` remain the `[overlap]` extra. `[model]` is a
kept-but-empty alias so existing `vdjtools[model]` pins still resolve.

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
  `appendix/bench_pgen.py`). The VJ enumeration had been ~1000x slower than OLGA — the TM fixes it.
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
- **DONE real-data EM comparison** `appendix/bench_em.py` — `infer_native` on real nonfunc TRB+TRD reads
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
  them automatically). `learned` = native EM on real HF out-of-frame reads (`appendix/build_bundled_models.py`,
  2k clonotypes/locus, held-out LL improves on every locus).
- **DONE arda-anchored D-D learning** — unregularized D-D EM over-attributes tandems on real data
  (identifiability; TRB→0.28). Two regularizers, both native==Python exact: **`dd_allowed`** per-read gate
  (a read may be n_D=2 only where arda called a `d2_call`) and **`nd_prior`** Dirichlet single-D pseudocount.
  `infer/infer_native(..., dd_allowed=, nd_prior=)`; native `estep_batch(..., dd_allowed)`. Anchored learned
  D-loci: TRB **0.000**, TRD **0.006**, IGH **0.009** (plausible; arda hard-call ~4%). `test_dd_anchor_and_prior…`.
- **DONE marimo explorer** `notebooks/model_explorer.py` — reactive Bayes-net/entropy/MI/marginal explorer
  over any bundled model (OLGA vs learned); `[examples]` extra. README/docs/SOURCES updated.
- **TODO** arda full-length V/J germline helper still needed for arda-native stitching (P1c residual).

Model schema notes: `ndel` is **biological** (neg = palindromic P-nt); dinucleotide row
`(from_nt,to_nt,p)=P(next|prev)` (OLGA's col-stochastic `R[next,prev]`); validation allows a group
to sum to 1 **or 0** (undefined conditional for an unused gene, kept for gene-index alignment).
Pgen/gen/EM invariant: **V and J each contribute ≥1 nt** to the CDR3 (OLGA-compatible).
