# CLAUDE.md — vdjtools v2

## What this is
Clean-room **Python + C++** rewrite (v2, GPL-3.0) of the legacy Groovy/Java vdjtools. TCR/BCR
repertoire analysis on the **AIRR schema + polars**, minimal OO, built on the antigenomics ecosystem:
**seqtree** (fuzzy search / e-values), **vdjmatch** (overlap + TCRnet), **arda** (AIRR annotation +
markup repair; brings conda/mmseqs2).

API surface is in [`skills/vdjtools/SKILL.md`](skills/vdjtools/SKILL.md); the release-by-release
narrative in `CHANGELOG.md`; data provenance and the numbers of record in `SOURCES.md`.

## Layout
- `python/vdjtools/` — package (src-layout). Subpackages `io model stats features overlap preprocess
  biomarker sc cli`, **lazy-loaded** — `import vdjtools` pulls no heavy deps.
- `src/`, `include/vdjtools/` — C++ core; `src/_bindings.cpp` → the `vdjtools._core` pybind11 ext.
  **Native only for Pgen DP, generation sampler, EM E-step** — Python-first everywhere else.
- `tests/{cpp,python}/`, `docs/` (Sphinx + pydata, gh-pages), `.github/workflows/`,
  `CMakeLists.txt` (scikit-build-core + pybind11, C++20), `environment.yml`, `setup.sh`.
- `appendix/` — **library-only**: the LaTeX theory appendix, `build_bundled_models.py` (builds the
  models shipped in the wheel), `concordance.py` (validates them against the OLGA oracle). Nothing
  that merely *uses* the library.

**Repo split (2026-07-17)** — benchmarks live in `~/vcs/projects/2026-vdjtools-benchmark`
(`bench/` scripts, the confound gates, `scripts/*.sbatch`). They *use* vdjtools, they aren't part of
it. NOTE: That directory **is not a git repo yet**, and its scripts hardcode cluster paths.

## Build / test / run
```bash
uv venv && source .venv/bin/activate && uv pip install -e ".[dev,test]"   # builds _core
pytest tests/python -q
cmake -S . -B build -DVDJTOOLS_TESTS=ON && cmake --build build && ctest --test-dir build
sphinx-build -W --keep-going -b html docs docs/_build/html                # zero-warning gate
```
`bash setup.sh --dev-parents` editable-installs `../seqtree ../arda ../vdjmatch` if present. Flip
`editable.rebuild=true` in pyproject during C++-heavy work for recompile-on-import.

**All three antigenomics engines are BASE deps** (`arda-mapper`, `seqtree`, `vdjmatch`). Rule: *if
the README advertises it, a plain `pip install vdjtools` must deliver it* — every advertised
capability delegates to one of the three and none has a fallback. All are imported lazily, so
`import vdjtools` stays light; between them they add only `requests`. mmseqs2 is needed only for
arda's annotate path. `[model]`/`[preprocess]` are kept-but-empty aliases so old pins resolve.
`test_smoke.py` pins this contract.

**Heavy tests run on Aldan-3, not the Mac** — benchmarks (`RUN_BENCHMARK=1`), full-locus
EM/concordance, 100k+ read runs. Drive it with the `aldan3` CLI (`../aldan3-client`); every
subcommand takes `--json`: `aldan3 slurm submit <script.sh> [-- ARGS…] [--env/--cpus/--mem/--time]`,
then `queue` · `log <id>` · `hist <id>` · `cancel <id>`; `slurm template cpu|gpu|array` scaffolds.

## Git model
`master` = v2 (tagged releases) ← `dev` (integration) ← `feature/*` (one per phase). **Legacy v1.x is
on `legacy-1.x` under tags `v0.0.1`..`1.2.1` — do not disturb**; the v2 history is an orphan root.
Carried-over legacy resources (`aa_property_table.txt`, `cdr3contact.txt`, `vj_families.txt`) and
format-conversion fixtures live there — pull them over when a phase needs them.

**Worktrees**: one worktree ↔ one `feature/*` branch —
`git worktree add .claude/worktrees/<name> -b feature/<name>`. Never two features in one worktree;
`.claude/` is gitignored, never commit it. Merge the finished branch into `dev`, then
`git worktree remove`.

## Conventions & invariants
- AIRR Rearrangement/Cell + polars `pl.DataFrame` in and out; minimal OO (thin index classes only).
- **arda germline = single source of truth**: all V/D/J germline + CDR3 anchors resolve from arda by
  allele name via `model.reference.load_germline` (arda's anchor convention is byte-identical to
  OLGA's). **Never mix germline sources within a model** — OLGA bootstrap models keep OLGA germline
  (exact-Pgen fidelity); arda-native EM models use arda. Raw anchor *indices* can differ by whole
  framework codons (IMGT drift) while the CDR3-region germline is identical — harmless unless mixed.
- **Delegate rather than reimplement**: overlap/TCRnet → vdjmatch; annotation/markup/scenarios →
  arda; search/e-value → seqtree.
- **Never modify non-dependency libraries.** Only `arda`, `vdjmatch`, `seqtree` are dependencies.
  Everything else under `~/vcs/code/` (mirpy, IGoR, OLGA, pygor3, …) is **reference/oracle only,
  read-only** — cross-validate against them, surface bugs to the owner, never edit.
- Native code goes through the single `_core` ext.
- **Model schema**: `ndel` is *biological* (negative = palindromic P-nt); dinucleotide row
  `(from_nt,to_nt,p) = P(next|prev)` (OLGA's col-stochastic `R[next,prev]`); validation allows a
  group to sum to 1 **or 0** (undefined conditional for an unused gene, kept for index alignment).
- **Pgen/gen/EM invariant: V and J each contribute ≥1 nt to the CDR3** (OLGA-compatible). Getting
  this wrong made nt Pgen up to 0.34% high on heavily-deleted sequences.

## WARNING: Traps that produced silent wrong answers
- **`native` Pgen allele guard.** `vi.get(v, -1)` mapped an unrecognised V/J to `-1` = *marginalize
  over all V/J*. The model is keyed by **allele** and real repertoires carry **gene-level** `v_call`
  (`TRBV9`), so `pgen_aa(m, cdr3, "TRBV9", "TRBJ2-3")` silently returned the V/J-agnostic value —
  **2.38× too high**, no error. `_gene_idx` now raises; pass `*01` on older pins.
- **`native.pack` cached by `id(model)`** — CPython id-reuse returned a stale `PackedModel` after a
  TRB→TRD switch in one process (89 vs 18 V, M-step crash). It now stores and verifies the model ref.
- **`_align_init` collapsed germline-identical paralogs**: `max(...)` takes the first of an exact tie
  (TRBV6-2/6-5/6-6 tie identically), seeding the rest at P(V)=0, which the E-step's `if pv==0:
  continue` makes **absorbing** — 68/89 TRB V alleles zeroed in shipped models. Votes now split
  across tied genes. Soft realign was *not* the cause and made it worse (mass avalanches onto the
  most permissive germline) — removed.
- **A held-out-LL claim that validated nothing**: it was the EM's own training objective, which EM
  increases monotonically by construction. Use `appendix/compare_models.py` for real held-out +
  oracle comparison. The same note's "2k clonotypes/locus" cap was also real and wrong.
- **`sc.paired_pgen` was 100% null on real CellRanger data** — and silently. CellRanger reports
  **gene**-level V/J (`TRBV10-3`); the model is keyed by allele; `native.pgen_aa` raises on a gene
  name *on purpose* (see the `-1` trap above) — and `_chain_pgen` swallowed it with a bare
  `except Exception: return None`. 27,268/27,268 dCODE receptors scored null with no signal.
  `paired_pgen` now resolves gene → representative allele (`*01`) explicitly (`resolve_genes=False`
  to opt out), warns when a whole locus is null, and catches only `(KeyError, ValueError)`.
  **A bare `except` over a call that raises deliberately re-creates the very bug the raise prevents.**
- **`collapse_alleles` relabelled an ORF allele's germline `*01`** (fixed 3.9.1). The representative
  key was `(len(cut_segment), usage, name)` and **usage was only passed for V** — so on J and D the
  key degenerated to the allele *name* after the length tie and `max` took the lexicographically
  last one. `TRBJ2-7*01` (F, `SYEQYF`) and `*02` (**ORF**, `SYEQYV`) are both 19 nt, so `*02` won
  and shipped under `*01`'s label: `Pgen("CASSIRSSYEQYF"|TRBJ2-7*01)` was **exactly 0** on the
  default `collapse=True` path, no error, in every version that shipped these models — **864/864**
  TRBJ2-7 junctions drawn from the model itself scored `pgen_nt == 0`. Key is now
  `(length, IMGT functionality, usage, prefer *01, name)`; **length still leads** — leading with
  functionality installs an empty germline on `TRBV23/OR9-2` (its only non-P allele has none), the
  same silent zero from the other side. 65 genes' germlines moved; `collapse=False` is the escape
  hatch. **arda's `cdr3_anchors.tsv` (`functionality`/`status`/`templated_aa`) is the reference for
  anchor questions** — vdjtools' own `anchor` column is an nt offset into `full_germline`, and
  reading it as a codon index into `cut_segment` flags 822 correct entries.
- **`from_olga(derive_orf=True)` sliced J germlines on the V side of the anchor** (fixed 3.9.1):
  V's CDR3 region is `full[anchor:]`, J's is `full[:anchor+3]`, and one line used the V form for
  both. 11 J alleles in the bundled `learned` human TRA model still carry the framework *downstream*
  of Phe118 — **a model rebuild is needed to clear them**; `test_collapse.py` pins them by name
  until then. Size it against the unaffected `arda` TRA model, not against the `learned` model's own
  usage: **`TRAJ35` (functional) holds 2.93e-03 in `arda` and 3.79e-06 in `learned` — 774× lower**,
  while unaffected controls (`TRAJ33`, `TRAJ42`) agree between the two within 25%.
- **A germline defect in an EM-fit model is not a fixed-size error — the model learns around it.**
  Wrong germline → reads stop scoring against the gene → EM drives its usage toward 0 → the gene's
  mass, measured *after* the fit, looks negligible. The 11 alleles above hold 0.33% of `learned`
  TRA's J mass, which reads as "ignore it"; against an unaffected reference the real suppression is
  three orders of magnitude on a functional gene. **Never size a germline/annotation bug from the
  affected model's own posterior.** Compare against a reference fit that does not share the defect
  (here `arda` vs `learned`), or against pre-EM read counts.
- **A barcoded AIRR table sniffed as bulk `"airr"`** and `read_airr` pooled reads across cells,
  dropping `cell_id` with no error. `sniff_format` now returns `"airr_cell"` and `io.read` refuses
  it, pointing at `sc.read_airr_cell`; `fmt="airr"` still pools deliberately.
- **`str.len_chars()` is UInt32** — `len - 2*flank` underflows on short junctions; cast first.
- **`_lower_map` is exact-lowercase**, so MiGEC's space-separated column picks never match the
  MiTCR/tcR dotted dialect (`Read.count`, `CDR3.nucleotide.sequence`) — it needs its own reader.
- A docs figure once claimed nt Pgen "2.4 ms / 89×" — **never real**. Quote benchmark numbers from
  `2026-vdjtools-benchmark/bench/`, not from prose.

## Open loops / next steps
- **Single-cell interop (`feature/single-cell-interop`) landed for 3.8.0.** Everything routes
  through ONE flat AIRR Rearrangement table (`sequence_id` + `cell_id`) in `sc/airr.py` — that is
  what scirpy, dandelion AND scRepertoire all read, so it is one emitter + thin adapters, not four
  bridges. **None of them consumes AIRR `Cell` objects**; cell state lives in `adata.obs` /
  `Dandelion.metadata` / Seurat `meta.data`. `write_airr_cell` stays a spec-faithful export, not an
  interop path. Asymmetry on purpose: *writing* a container delegates to the library that owns it
  (no schema copy to drift), *reading* is ours (`from_scirpy` needs only `awkward`, `read_h5ddl`
  only `h5py`). TODO: nothing blocking — possible next steps are a `Dandelion` polars-backend
  fast path (`ddl.set_backend("polars")` takes polars frames directly) and reconciling
  `resolve_chains` with scirpy's `chain_qc` `receptor_subtype` vocabulary (report alongside, never
  overwrite — that is how a QC call gets lost).
- **Dev-env note**: the worktree needs its OWN venv — the editable install's meta-path finder wins
  over `PYTHONPATH`, so a symlinked `_core` will NOT redirect `import vdjtools` to a worktree.
  Also: `cd` inside a backgrounded/`/tmp` command resets the shell cwd back to the MAIN repo, so
  relative-path writes silently land there. Use absolute paths for edits.
- **Phase 1 (`feature/model-engine`) is functionally complete** — native nt/aa Pgen via the
  Murugan/OLGA `Pi_L·Pi_R` transfer matrix (single-D and D-D), batch-parallel Pgen, threaded EM
  E-step, D-D learning with arda anchoring, 7-locus concordance `r(log10 Pgen)=1.00000`, bundled
  `olga` + `learned` models shipped in **v2.9.0**.
- **TODO — arda full-length V/J germline helper** for arda-native stitching (the P1c residual).
  `derive_orf` covers the ORF-usage case but not full-length stitching.
- **TODO — native perf gaps**: (a) VJ / Hamming-1 codon-boundary sweep is **set aside** — the V/J
  combine boundary migrates with the delJ sum, so there's no clean O(1)-per-codon leave-one-out and a
  forced rewrite risks the exact-Pgen invariant for a non-bottleneck; batch parallelization is the
  exact win instead. (b) native generation sampler — low priority, Python is already fast.
- **D-D has no exact skip.** It's 2.5× single-D, but **0%** of reads have a zero D-D contribution
  (two 1-nt Ds plus insertions tile almost any mid). A length gate biases learned P(n_D=2) by −2.5%
  (−11% if it skips half). Keep D-D exact by default; the clean exact speedup is read-parallelization.
- **AS/B27 campaign findings not to re-derive**: B27 is 26/27 confounded with AS in that cohort, so
  only the **B27-matched** contrast separates disease from carriage; carriage among healthy is null.
  **V-pinning is load-bearing.** 38 of 40 `old/` donors reappear in `new/` — **no independent
  replication exists there, never pool.** At that n, BH cannot clear 0.05 over 273 features
  (min attainable p = 3.6e-3) — the *ranking* is the result. The VDJdb oracle for this motif is
  **partly circular** (same group). Handoff: `2026-mirpy-analysis` branch `as-b27-embedding`.
- **Phase 15 (`feature/model-workshop`) landed in v3.3.0** — `from_germline` (custom V(D)J
  libraries; `from_arda` is now a wrapper, output byte-identical), `check_model`, `model/score.py`
  (nt-Pgen log-likelihood + AIC/BIC, Pgen-distribution comparison, scenario/sequence entropy and
  Hill q=1/q=2 diversity), `compare_models`/`compare_net_dot`, a persisted EM training log,
  `infer_frame`, `extend_alleles`, `data.build_all`, the `vdjtools model` CLI sub-app,
  `docs/model.rst` and `examples/model_workshop.py`.
- **arda CLI drift bit us once**: `annotate_reads` shelled out to `arda rnaseq map -o …`, but arda
  2.19 turned `rnaseq` into the full map→assemble→correct preset with no stage positional and
  `-p/--out-prefix` instead of `-o`, so the real invocation exited 2. It is `arda map` now, and the
  pin is `arda-mapper>=2.19.0`. Note a `--help` smoke test would NOT have caught it (typer
  short-circuits `--help` before argument parsing) — exercise the real argv.
- **Unreachable deletion mass — diagnosed, and it was two separate things.**
  1. *Ours, fixed*: `collapse_alleles` chose each gene's representative germline by **usage**, and
     IMGT ships some truncated alleles — `IGKV3-20*02` is 11 nt against `*01`'s 30 and had the
     higher learned usage, so the collapsed gene got an 11-nt germline (relabelled `*01`) and
     stranded 25% of its own deletion distribution. Now: longest germline first, usage second, then
     the conditionals are projected onto that germline's reachable support. `collapse=True` (the
     default) is clean on every bundled model; `test_collapse.py` pins it.
  2. *OLGA's, deliberately kept*: the **uncollapsed** `olga` models carry it because OLGA's own
     `model_marginals.txt` does — verified against OLGA's raw arrays, same fractions to 4 dp
     (`IGHV4-30-4*01` 100%, `IGKJ4*02` 80.9%, `TRAV20*03` 54.7%), and our Pgen matches olga-pip
     exactly (ratio 1.000000) on every sequence OLGA will score. `IGHV4-30-4*01` has Pgen ≡ 0 in
     OLGA too. **Do not "fix" this** — it would break the exact-OLGA-Pgen invariant. `check_model`
     reports it as `warn` (not `error`) when `manifest.source` starts with `olga`.
- **`infer_nt` (aa→nt) landed in v3.6.0, natively.** Stage 1 is `native.best_aa_scenarios` — the
  same `Pi_L·Pi_R` transfer matrix `pgen_aa` sums over, with `max` and the winning `(V,delV)`/
  `(J,delJ)` carried in the state; stage 2 re-scores the survivors with the exact `pgen_nt`.
  2.5 ms/TRB, 0.5 ms/TRA (all of VDJdb in ~3 min); reproduces the brute-force oracle 25/25 TRG,
  19/19 TRA. **The Python scenario enumeration in `viterbi.py` is the reference, not the product**
  — it is ~600× slower on TRB and is selected only when a `prepare()`-d model is passed;
  `test_native_search_agrees_with_the_python_reference` is what makes it worth keeping.
  NOTE: Two things measured and rejected along the way: fixing the germline trim first and then picking
  the best codon per residue agrees with the oracle only 9/25 and 4/19 (a trim chosen before the
  codons pins a codon the optimum would have trimmed), and bounding the codon DP by 1.0 in
  branch-and-bound prunes nothing — the insertion chain costs ~0.4/nt, so the cutoff sits six orders
  of magnitude too high.
- No emoji anywhere in the repo (v3.6.0); the old `⛔`/`⚠` markers are `WARNING:` / `NOTE:`.
- `ruff check .` is green as of v3.6.1, with `[tool.ruff.lint]` ignoring only `E702`/`E741`/`E731`
  (house style). Keep it green so a real finding is not buried in style hits. `ruff format` is NOT
  the house style — it would rewrite 171 files; do not run it.
  (The old `rescale.py` `empty_as_null` note is resolved: it is set explicitly at `model/rescale.py:62`.)
